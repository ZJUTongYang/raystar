#include <raystar/raystar_core.h>

#include <geometry_msgs/msg/pose_stamped.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <nav_msgs/msg/path.hpp>

#include <path_collision_oracle.h>

#include <sys/resource.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using raystar::GridMap;
using raystar::PathSolution;
using raystar::PlanningLimitReached;
using raystar::PlanningLimits;
using raystar::PlanningOutcome;
using raystar::PlanResult;
using raystar::Point2d;
using Clock = std::chrono::steady_clock;

constexpr int kDefaultWarmups = 3;
constexpr int kDefaultIterations = 20;
constexpr int kDefaultTimeoutMs = 5000;
constexpr std::size_t kDefaultMaxNodes = 10000;
constexpr std::uint64_t kFnvOffset = 1469598103934665603ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

struct Scenario {
  std::string name;
  GridMap map;
  Point2d start;
  Point2d goal;
  std::function<std::size_t(int)> expected_paths;
};

struct ScenarioDescriptor {
  const char* name;
  const char* description;
};

constexpr std::array<ScenarioDescriptor, 6> kScenarioCatalog{{
  {"open_256", "256x256 open grid; isolates map construction and the direct path"},
  {"single_obstacle_256", "One central 64x64 obstacle with two non-homotopic routes"},
  {"narrow_gate_256", "Two long obstacles leave a four-cell central gate and outer detours"},
  {"dense_lattice_192", "A deterministic 7x7 lattice of separated obstacle islands"},
  {"large_open_1024", "1024x1024 open grid; stresses large-map geometry construction"},
  {"bundled_testmap_50", "Repository testmap.pgm irregular-obstacle regression scene"},
}};

struct Options {
  std::string scenario = "all";
  std::vector<int> k_values = {1, 3, 10, 50};
  int warmups = kDefaultWarmups;
  int iterations = kDefaultIterations;
  std::string testmap_path;
  bool list_scenarios = false;
};

struct TrialRecord {
  bool accepted = false;
  bool success = false;
  bool request_satisfied = false;
  bool search_complete = false;
  bool path_validation_pass = false;
  bool terminal_consistency_pass = false;
  bool scenario_contract_pass = false;
  bool determinism_pass = true;
  bool rss_available = false;
  std::size_t found_paths = 0;
  std::size_t expected_paths = 0;
  std::size_t expanded_nodes = 0;
  std::size_t total_path_points = 0;
  double map_time_ms = 0.0;
  double plan_time_ms = 0.0;
  double wall_time_ms = 0.0;
  double shortest_cost = 0.0;
  double longest_cost = 0.0;
  long process_hwm_kib_after_plan = 0;
  PlanningOutcome outcome = PlanningOutcome::failed;
  PlanningLimitReached limit = PlanningLimitReached::none;
  std::uint64_t path_digest = kFnvOffset;
  std::string verdict = "FAIL_INTERNAL";
  std::string diagnostic;
};

std::string defaultTestmapPath() {
  std::error_code error;
  const std::filesystem::path executable = std::filesystem::read_symlink("/proc/self/exe", error);
  if (!error) {
    const std::filesystem::path installed =
      executable.parent_path().parent_path().parent_path() / "share/raystar/test/testmap.pgm";
    if (std::filesystem::is_regular_file(installed, error))
      return installed.string();
  }
  return RAYSTAR_PROFILE_TESTMAP_PATH;
}

GridMap makeOpenMap(unsigned int width, unsigned int height) {
  GridMap map;
  map.width = width;
  map.height = height;
  map.resolution = 1.0F;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.assign(static_cast<std::size_t>(width) * height, 0);
  return map;
}

void fillRectangle(
  GridMap& map, unsigned int left, unsigned int bottom, unsigned int right, unsigned int top) {
  if (left >= right || bottom >= top || right > map.width || top > map.height)
    throw std::invalid_argument("Invalid profiling obstacle rectangle");
  for (unsigned int y = bottom; y < top; ++y) {
    for (unsigned int x = left; x < right; ++x)
      map.data[static_cast<std::size_t>(y) * map.width + x] = 1;
  }
}

std::string readPgmToken(std::istream& stream) {
  while (stream) {
    stream >> std::ws;
    if (stream.peek() == '#') {
      stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }
    std::string token;
    stream >> token;
    return token;
  }
  throw std::runtime_error("Unexpected end of PGM header");
}

unsigned int parsePgmUnsigned(const std::string& token, const char* field) {
  std::size_t consumed = 0;
  unsigned long value = 0;
  try {
    value = std::stoul(token, &consumed);
  } catch (const std::exception&) {
    throw std::runtime_error(std::string("Invalid PGM ") + field);
  }
  if (consumed != token.size() || value > std::numeric_limits<unsigned int>::max())
    throw std::runtime_error(std::string("Invalid PGM ") + field);
  return static_cast<unsigned int>(value);
}

GridMap loadBundledPgm(const std::string& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream)
    throw std::runtime_error("Cannot open bundled profiling map: " + path);

  const std::string magic = readPgmToken(stream);
  const unsigned int width = parsePgmUnsigned(readPgmToken(stream), "width");
  const unsigned int height = parsePgmUnsigned(readPgmToken(stream), "height");
  const unsigned int maximum = parsePgmUnsigned(readPgmToken(stream), "maximum value");
  if (!stream || magic != "P5" || width == 0 || height == 0 || maximum != 255)
    throw std::runtime_error("Bundled profiling map must be an 8-bit binary PGM");
  char separator = '\0';
  stream.get(separator);
  if (!stream || !std::isspace(static_cast<unsigned char>(separator)))
    throw std::runtime_error("Bundled profiling map has no raster separator");
  if (separator == '\r' && stream.peek() == '\n')
    stream.get();

  if (height > std::numeric_limits<std::size_t>::max() / width)
    throw std::runtime_error("Bundled profiling map dimensions overflow size_t");
  const std::size_t cell_count = static_cast<std::size_t>(width) * height;
  if (cell_count > PlanningLimits::kDefaultMaxMapCells)
    throw std::runtime_error("Bundled profiling map exceeds the default map-cell limit");
  if (cell_count > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()))
    throw std::runtime_error("Bundled profiling map exceeds the stream-size limit");
  std::vector<unsigned char> pixels(cell_count);
  stream.read(reinterpret_cast<char*>(pixels.data()), static_cast<std::streamsize>(pixels.size()));
  if (stream.gcount() != static_cast<std::streamsize>(pixels.size()))
    throw std::runtime_error("Bundled profiling map has truncated pixel data");

  GridMap map = makeOpenMap(width, height);
  map.resolution = 0.1F;
  map.origin_x = 2.0;
  map.origin_y = 3.0;
  for (std::size_t index = 0; index < pixels.size(); ++index)
    map.data[index] = pixels[index] > 200 ? 0 : 1;
  return map;
}

Scenario makeScenario(const std::string& name, const std::string& testmap_path) {
  if (name == "open_256") {
    return Scenario{
      name, makeOpenMap(256, 256), {4.5, 4.5}, {250.5, 250.5}, [](int) { return 1U; }};
  }
  if (name == "single_obstacle_256") {
    GridMap map = makeOpenMap(256, 256);
    fillRectangle(map, 96, 96, 160, 160);
    return Scenario{name, std::move(map), {16.5, 128.5}, {239.5, 128.5}, [](int k) {
                      return static_cast<std::size_t>(std::min(k, 2));
                    }};
  }
  if (name == "narrow_gate_256") {
    GridMap map = makeOpenMap(256, 256);
    fillRectangle(map, 8, 120, 124, 136);
    fillRectangle(map, 128, 120, 248, 136);
    return Scenario{name, std::move(map), {126.5, 16.5}, {126.5, 239.5}, [](int k) {
                      return static_cast<std::size_t>(std::min(k, 5));
                    }};
  }
  if (name == "dense_lattice_192") {
    GridMap map = makeOpenMap(192, 192);
    for (unsigned int row = 0; row < 7; ++row) {
      for (unsigned int column = 0; column < 7; ++column) {
        const unsigned int left = 20 + 22 * column;
        const unsigned int bottom = 20 + 22 * row;
        fillRectangle(map, left, bottom, left + 10, bottom + 10);
      }
    }
    return Scenario{name, std::move(map), {4.5, 95.5}, {187.5, 95.5}, [](int k) {
                      return static_cast<std::size_t>(k);
                    }};
  }
  if (name == "large_open_1024") {
    return Scenario{
      name, makeOpenMap(1024, 1024), {4.5, 4.5}, {1018.5, 1018.5}, [](int) { return 1U; }};
  }
  if (name == "bundled_testmap_50") {
    return Scenario{name, loadBundledPgm(testmap_path), {2.5, 2.5}, {46.5, 46.5}, [](int k) {
                      return static_cast<std::size_t>(k);
                    }};
  }
  throw std::invalid_argument("Unknown scenario: " + name);
}

const char* outcomeName(PlanningOutcome outcome) {
  switch (outcome) {
  case PlanningOutcome::invalid_request:
    return "invalid_request";
  case PlanningOutcome::complete:
    return "complete";
  case PlanningOutcome::no_path:
    return "no_path";
  case PlanningOutcome::limit_reached:
    return "limit_reached";
  case PlanningOutcome::failed:
    return "failed";
  }
  return "unknown";
}

const char* limitName(PlanningLimitReached limit) {
  switch (limit) {
  case PlanningLimitReached::none:
    return "none";
  case PlanningLimitReached::max_nodes:
    return "max_nodes";
  case PlanningLimitReached::timeout:
    return "timeout";
  case PlanningLimitReached::cancelled:
    return "cancelled";
  case PlanningLimitReached::max_path_points:
    return "max_path_points";
  }
  return "unknown";
}

void hashBytes(std::uint64_t& digest, const void* data, std::size_t size) {
  const auto* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t index = 0; index < size; ++index) {
    digest ^= bytes[index];
    digest *= kFnvPrime;
  }
}

void hashDouble(std::uint64_t& digest, double value) {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  hashBytes(digest, &bits, sizeof(bits));
}

std::uint64_t mapDigest(const GridMap& map) {
  std::uint64_t digest = kFnvOffset;
  hashBytes(digest, &map.width, sizeof(map.width));
  hashBytes(digest, &map.height, sizeof(map.height));
  hashBytes(digest, &map.resolution, sizeof(map.resolution));
  hashDouble(digest, map.origin_x);
  hashDouble(digest, map.origin_y);
  hashBytes(digest, map.data.data(), map.data.size() * sizeof(map.data.front()));
  return digest;
}

nav_msgs::msg::OccupancyGrid makeOracleMap(const GridMap& map) {
  nav_msgs::msg::OccupancyGrid result;
  result.header.frame_id = "map";
  result.info.width = map.width;
  result.info.height = map.height;
  result.info.resolution = map.resolution;
  result.info.origin.position.x = map.origin_x;
  result.info.origin.position.y = map.origin_y;
  result.info.origin.orientation.w = 1.0;
  result.data.resize(map.data.size());
  for (unsigned int y = 0; y < map.height; ++y) {
    for (unsigned int x = 0; x < map.width; ++x) {
      const std::size_t index = static_cast<std::size_t>(y) * map.width + x;
      const bool forced_border = x == 0 || y == 0 || x + 1 == map.width || y + 1 == map.height;
      result.data[index] = map.data[index] != 0 || forced_border ? 100 : 0;
    }
  }
  return result;
}

nav_msgs::msg::Path makeOraclePath(const GridMap& map, const PathSolution& solution) {
  nav_msgs::msg::Path path;
  path.header.frame_id = "map";
  for (const auto& point : solution.projectedPath()) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = path.header.frame_id;
    pose.pose.position.x = map.origin_x + point.first * static_cast<double>(map.resolution);
    pose.pose.position.y = map.origin_y + point.second * static_cast<double>(map.resolution);
    pose.pose.orientation.w = 1.0;
    path.poses.emplace_back(std::move(pose));
  }
  return path;
}

bool validateAncestry(const raystar::RaystarCore& core,
                      const PathSolution& solution,
                      std::string& diagnostic) {
  const auto& nodes = core.getNodes();
  const auto& chain = solution.path_node_index_;
  if (chain.empty() || chain.front() != 0) {
    diagnostic = "path ancestry does not begin at root node zero";
    return false;
  }
  if (chain.size() != solution.turning_points_.size() + 1) {
    diagnostic = "path ancestry length does not match its turning-point count";
    return false;
  }
  std::set<int> seen;
  for (std::size_t index = 0; index < chain.size(); ++index) {
    const int node_index = chain[index];
    if (node_index < 0 || static_cast<std::size_t>(node_index) >= nodes.size()) {
      diagnostic = "path ancestry contains an out-of-range node";
      return false;
    }
    if (!seen.insert(node_index).second) {
      diagnostic = "path ancestry contains a cycle";
      return false;
    }
    const int expected_parent = index == 0 ? -1 : chain[index - 1];
    if (nodes[static_cast<std::size_t>(node_index)].parentIndex() != expected_parent) {
      diagnostic = "path ancestry parent link is inconsistent";
      return false;
    }
  }
  if (nodes[static_cast<std::size_t>(chain.back())].localShortestPath() !=
      solution.turning_points_) {
    diagnostic = "path ancestry leaf geometry differs from the returned turning points";
    return false;
  }
  return true;
}

bool validateReturnedPaths(const Scenario& scenario,
                           std::uint64_t original_map_digest,
                           const raystar::RaystarCore& core,
                           const PlanResult& result,
                           TrialRecord& record) {
  if (mapDigest(scenario.map) != original_map_digest) {
    record.diagnostic = "planner mutated its input map";
    return false;
  }

  const auto oracle_map = makeOracleMap(scenario.map);
  double previous_cost = -std::numeric_limits<double>::infinity();
  std::set<std::uint64_t> individual_digests;
  for (const auto& solution : result.path_solutions) {
    const auto points = solution.projectedPath();
    if (points.size() < 2 || points.front() != scenario.start || points.back() != scenario.goal) {
      record.diagnostic = "path does not preserve the continuous request endpoints";
      return false;
    }

    double measured_cost = 0.0;
    std::uint64_t individual_digest = kFnvOffset;
    for (std::size_t index = 0; index < points.size(); ++index) {
      if (!std::isfinite(points[index].first) || !std::isfinite(points[index].second)) {
        record.diagnostic = "path contains a non-finite coordinate";
        return false;
      }
      hashDouble(individual_digest, points[index].first);
      hashDouble(individual_digest, points[index].second);
      if (index > 0) {
        measured_cost += std::hypot(points[index].first - points[index - 1].first,
                                    points[index].second - points[index - 1].second);
      }
    }
    if (!std::isfinite(solution.path_cost_) ||
        std::abs(measured_cost - solution.path_cost_) > 1e-9 * std::max(1.0, measured_cost)) {
      record.diagnostic = "path cost differs from its projected polyline length";
      return false;
    }
    if (std::isfinite(previous_cost)) {
      const double cost_order_tolerance =
        1e-9 * std::max({1.0, std::abs(previous_cost), std::abs(solution.path_cost_)});
      if (solution.path_cost_ + cost_order_tolerance < previous_cost) {
        record.diagnostic = "paths are not ordered by nondecreasing cost";
        return false;
      }
    }
    previous_cost = solution.path_cost_;
    if (!individual_digests.insert(individual_digest).second) {
      record.diagnostic = "planner returned a duplicate projected path";
      return false;
    }
    std::uint64_t solution_digest = individual_digest;
    hashDouble(solution_digest, solution.path_cost_);
    hashBytes(record.path_digest, &solution_digest, sizeof(solution_digest));

    const auto oracle_path = makeOraclePath(scenario.map, solution);
    const auto collision = raystar::test_oracle::validatePath(oracle_map, oracle_path);
    if (!collision.collisionFree()) {
      record.diagnostic =
        "independent collision oracle rejected a returned path: " + collision.diagnostic;
      return false;
    }
    const auto self_intersection = raystar::test_oracle::validateNoSelfIntersection(oracle_path);
    if (!self_intersection.intersectionFree()) {
      record.diagnostic = "independent self-intersection oracle rejected a returned path: " +
                          self_intersection.diagnostic;
      return false;
    }
    if (!validateAncestry(core, solution, record.diagnostic))
      return false;

    record.total_path_points += points.size();
  }

  if (!result.path_solutions.empty()) {
    record.shortest_cost = result.path_solutions.front().path_cost_;
    record.longest_cost = result.path_solutions.back().path_cost_;
  }
  return true;
}

bool validateTerminalConsistency(const raystar::RaystarCore& core,
                                 const PlanResult& result,
                                 TrialRecord& record) {
  const auto fail = [&](const std::string& message) {
    if (record.diagnostic.empty())
      record.diagnostic = message;
    return false;
  };
  if (!std::isfinite(record.map_time_ms) || record.map_time_ms < 0.0 ||
      !std::isfinite(record.plan_time_ms) || record.plan_time_ms < 0.0 ||
      !std::isfinite(record.wall_time_ms) || record.wall_time_ms < 0.0) {
    return fail("planner returned a non-finite or negative timing value");
  }
  if (!record.rss_available)
    return fail("getrusage could not provide a process high-water RSS sample");
  if (result.success != !result.path_solutions.empty())
    return fail("success does not match whether at least one path was returned");
  if (result.expanded_nodes != core.getNodes().size())
    return fail("expanded_nodes differs from the retained search tree size");

  if (result.outcome == PlanningOutcome::complete) {
    if (result.limit_reached != PlanningLimitReached::none)
      return fail("complete outcome reports a reached resource limit");
    if (result.path_solutions.empty())
      return fail("complete outcome returned no path");
  } else if (result.outcome == PlanningOutcome::no_path) {
    if (result.limit_reached != PlanningLimitReached::none || result.success ||
        !result.path_solutions.empty()) {
      return fail("no_path outcome has inconsistent success, path, or limit fields");
    }
  } else if (result.outcome == PlanningOutcome::limit_reached) {
    if (result.limit_reached == PlanningLimitReached::none)
      return fail("limit_reached outcome does not identify a resource limit");
  } else if (result.limit_reached != PlanningLimitReached::none) {
    return fail("non-limit outcome reports a reached resource limit");
  }
  return true;
}

bool validateScenarioContract(const PlanResult& result, TrialRecord& record) {
  const auto fail = [&](const std::string& message) {
    if (record.diagnostic.empty())
      record.diagnostic = message;
    return false;
  };
  if (result.outcome != PlanningOutcome::complete ||
      result.limit_reached != PlanningLimitReached::none) {
    return fail("planner did not complete the scenario search without a resource limit");
  }
  if (result.path_solutions.size() != record.expected_paths)
    return fail("returned path count differs from the scenario contract");
  return true;
}

std::string classifyVerdict(const TrialRecord& record) {
  if (!record.determinism_pass)
    return "FAIL_DETERMINISM";
  if (!record.path_validation_pass)
    return "FAIL_CORRECTNESS";
  if (!record.terminal_consistency_pass)
    return "FAIL_TERMINAL";
  if (record.scenario_contract_pass)
    return "PASS";
  if (record.outcome == PlanningOutcome::limit_reached)
    return std::string("INCOMPLETE_") + limitName(record.limit);
  if (record.outcome == PlanningOutcome::no_path)
    return "FAIL_NO_PATH";
  if (record.outcome == PlanningOutcome::invalid_request)
    return "FAIL_INVALID_REQUEST";
  if (record.outcome == PlanningOutcome::failed)
    return "FAIL_INTERNAL";
  return "FAIL_CONTRACT";
}

void enforceDeterminism(const TrialRecord& baseline, TrialRecord& record) {
  const bool matches = record.success == baseline.success && record.outcome == baseline.outcome &&
                       record.limit == baseline.limit &&
                       record.request_satisfied == baseline.request_satisfied &&
                       record.search_complete == baseline.search_complete &&
                       record.found_paths == baseline.found_paths &&
                       record.expected_paths == baseline.expected_paths &&
                       record.expanded_nodes == baseline.expanded_nodes &&
                       record.total_path_points == baseline.total_path_points &&
                       record.path_digest == baseline.path_digest;
  if (matches)
    return;
  record.determinism_pass = false;
  record.accepted = false;
  record.verdict = classifyVerdict(record);
  record.diagnostic = "repeat differs from the first invocation's result or expanded node count";
}

std::pair<bool, long> processHwmKib() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return {false, 0};
  return {true, usage.ru_maxrss};
}

TrialRecord runTrial(const Scenario& scenario,
                     int k,
                     raystar::RaystarCore& core,
                     const PlanningLimits& limits) {
  TrialRecord record;
  record.expected_paths = scenario.expected_paths(k);
  const std::uint64_t original_map_digest = mapDigest(scenario.map);
  const auto wall_start = Clock::now();
  const PlanResult result =
    core.plan(scenario.map, scenario.start, scenario.goal, k, false, limits);
  const auto wall_end = Clock::now();

  record.success = result.success;
  record.found_paths = result.path_solutions.size();
  record.expanded_nodes = result.expanded_nodes;
  record.map_time_ms = result.map_time_ms;
  record.plan_time_ms = result.plan_time_ms;
  record.wall_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(wall_end - wall_start).count() / 1000.0;
  record.outcome = result.outcome;
  record.limit = result.limit_reached;
  record.request_satisfied = result.outcome == PlanningOutcome::complete &&
                             record.found_paths == static_cast<std::size_t>(k);
  record.search_complete =
    result.outcome == PlanningOutcome::complete || result.outcome == PlanningOutcome::no_path;
  const auto [rss_available, process_hwm_kib] = processHwmKib();
  record.rss_available = rss_available;
  record.process_hwm_kib_after_plan = process_hwm_kib;
  record.path_validation_pass =
    validateReturnedPaths(scenario, original_map_digest, core, result, record);
  record.terminal_consistency_pass = validateTerminalConsistency(core, result, record);
  record.scenario_contract_pass = validateScenarioContract(result, record);
  record.accepted = record.path_validation_pass && record.terminal_consistency_pass &&
                    record.scenario_contract_pass && record.determinism_pass;
  record.verdict = classifyVerdict(record);
  return record;
}

std::string digestString(std::uint64_t digest) {
  std::ostringstream stream;
  stream << std::hex << std::setw(16) << std::setfill('0') << digest;
  return stream.str();
}

std::size_t occupiedCellCount(const GridMap& map) {
  return static_cast<std::size_t>(
    std::count_if(map.data.begin(), map.data.end(), [](std::uint8_t value) { return value != 0; }));
}

void printHeader() {
  std::cout
    << "schema_version,scenario,width,height,resolution,input_occupied_cells,start_x,start_y,"
       "goal_x,goal_y,k,phase,iteration,allow_self_crossing,max_nodes,timeout_ms,success,"
       "outcome,limit_reached,request_satisfied,search_complete,found_paths,"
       "expected_paths,expanded_nodes,map_time_ms,plan_time_ms,wall_time_ms,"
       "total_path_points,shortest_cost_cells,longest_cost_cells,path_digest,"
       "path_validation_pass,terminal_consistency_pass,scenario_contract_pass,"
       "determinism_pass,rss_available,process_hwm_kib_after_plan,verdict,acceptance\n";
}

void printRecord(const Scenario& scenario,
                 int k,
                 const std::string& phase,
                 int iteration,
                 const PlanningLimits& limits,
                 const TrialRecord& record) {
  std::cout << std::fixed << std::setprecision(3) << "2," << scenario.name << ','
            << scenario.map.width << ',' << scenario.map.height << ',' << scenario.map.resolution
            << ',' << occupiedCellCount(scenario.map) << ',' << scenario.start.first << ','
            << scenario.start.second << ',' << scenario.goal.first << ',' << scenario.goal.second
            << ',' << k << ',' << phase << ',' << iteration << ",false," << limits.max_nodes << ','
            << limits.planning_timeout.count() << ',' << (record.success ? "true" : "false") << ','
            << outcomeName(record.outcome) << ',' << limitName(record.limit) << ','
            << (record.request_satisfied ? "true" : "false") << ','
            << (record.search_complete ? "true" : "false") << ',' << record.found_paths << ','
            << record.expected_paths << ',' << record.expanded_nodes << ',' << record.map_time_ms
            << ',' << record.plan_time_ms << ',' << record.wall_time_ms << ','
            << record.total_path_points << ',' << record.shortest_cost << ',' << record.longest_cost
            << ',' << digestString(record.path_digest) << ','
            << (record.path_validation_pass ? "true" : "false") << ','
            << (record.terminal_consistency_pass ? "true" : "false") << ','
            << (record.scenario_contract_pass ? "true" : "false") << ','
            << (record.determinism_pass ? "true" : "false") << ','
            << (record.rss_available ? "true" : "false") << ',' << record.process_hwm_kib_after_plan
            << ',' << record.verdict << ',' << (record.accepted ? "PASS" : "FAIL") << '\n';
  std::cout.flush();
}

int parseNonnegative(const std::string& value, const std::string& option) {
  std::size_t consumed = 0;
  int parsed = 0;
  try {
    parsed = std::stoi(value, &consumed);
  } catch (const std::exception&) {
    throw std::invalid_argument(option + " requires an integer");
  }
  if (consumed != value.size() || parsed < 0)
    throw std::invalid_argument(option + " requires a nonnegative integer");
  return parsed;
}

std::vector<int> parseKValues(const std::string& value) {
  std::vector<int> result;
  std::stringstream stream(value);
  std::string item;
  while (std::getline(stream, item, ',')) {
    const int parsed = parseNonnegative(item, "--k-values");
    if (parsed <= 0 || parsed > 100)
      throw std::invalid_argument("--k-values entries must be between 1 and 100");
    result.push_back(parsed);
  }
  if (result.empty())
    throw std::invalid_argument("--k-values must contain at least one value");
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

void printUsage(const char* program) {
  std::cout << "Usage: " << program << " [options]\n"
            << "  --scenario NAME       Run one scenario (default: all)\n"
            << "  --k-values LIST       Comma-separated K values (default: 1,3,10,50)\n"
            << "  --warmups N           Unreported warm-up runs per case (default: 3)\n"
            << "  --iterations N        Measured repeats per case (default: 20)\n"
            << "  --testmap PATH        Override bundled testmap.pgm path\n"
            << "  --list-scenarios      Print scenario names and exit\n"
            << "  --help                Show this message\n";
}

Options parseOptions(int argc, char** argv) {
  Options options;
  options.testmap_path = defaultTestmapPath();
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    const auto require_value = [&](const std::string& option) {
      if (index + 1 >= argc)
        throw std::invalid_argument(option + " requires a value");
      return std::string(argv[++index]);
    };
    if (argument == "--scenario")
      options.scenario = require_value(argument);
    else if (argument == "--k-values")
      options.k_values = parseKValues(require_value(argument));
    else if (argument == "--warmups")
      options.warmups = parseNonnegative(require_value(argument), argument);
    else if (argument == "--iterations")
      options.iterations = parseNonnegative(require_value(argument), argument);
    else if (argument == "--testmap")
      options.testmap_path = require_value(argument);
    else if (argument == "--list-scenarios")
      options.list_scenarios = true;
    else if (argument == "--help") {
      printUsage(argv[0]);
      std::exit(0);
    } else {
      throw std::invalid_argument("Unknown option: " + argument);
    }
  }
  if (options.iterations <= 0)
    throw std::invalid_argument("--iterations must be positive");
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parseOptions(argc, argv);

    if (options.list_scenarios) {
      for (const auto& descriptor : kScenarioCatalog)
        std::cout << descriptor.name << "\t" << descriptor.description << '\n';
      return 0;
    }

    const auto selected =
      std::find_if(kScenarioCatalog.begin(), kScenarioCatalog.end(), [&](const auto& item) {
        return item.name == options.scenario;
      });
    if (options.scenario != "all" && selected == kScenarioCatalog.end())
      throw std::invalid_argument("Unknown scenario: " + options.scenario);

    PlanningLimits limits;
    limits.max_k = 100;
    limits.max_nodes = kDefaultMaxNodes;
    limits.planning_timeout = std::chrono::milliseconds(kDefaultTimeoutMs);

    bool all_passed = true;
    printHeader();
    for (const auto& descriptor : kScenarioCatalog) {
      if (options.scenario != "all" && descriptor.name != options.scenario)
        continue;
      const Scenario scenario = makeScenario(descriptor.name, options.testmap_path);
      for (const int k : options.k_values) {
        raystar::RaystarCore core;

        const TrialRecord first = runTrial(scenario, k, core, limits);
        printRecord(scenario, k, "first", 0, limits, first);
        all_passed = all_passed && first.accepted;
        if (!first.accepted)
          std::cerr << scenario.name << " K=" << k << " first: " << first.diagnostic << '\n';

        for (int warmup = 0; warmup < options.warmups; ++warmup) {
          TrialRecord record = runTrial(scenario, k, core, limits);
          enforceDeterminism(first, record);
          all_passed = all_passed && record.accepted;
          if (!record.accepted) {
            std::cerr << scenario.name << " K=" << k << " warmup " << warmup << ": "
                      << record.diagnostic << '\n';
          }
        }

        for (int iteration = 0; iteration < options.iterations; ++iteration) {
          TrialRecord record = runTrial(scenario, k, core, limits);
          enforceDeterminism(first, record);
          printRecord(scenario, k, "measured", iteration + 1, limits, record);
          all_passed = all_passed && record.accepted;
          if (!record.accepted) {
            std::cerr << scenario.name << " K=" << k << " measured " << iteration + 1 << ": "
                      << record.diagnostic << '\n';
          }
        }
      }
    }
    return all_passed ? 0 : 1;
  } catch (const std::exception& exception) {
    std::cerr << "raystar_profile: " << exception.what() << '\n';
    return 2;
  }
}
