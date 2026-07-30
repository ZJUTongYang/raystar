#include <algorithm>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

#include <QApplication>
#include <pluginlib/class_loader.hpp>
#include <rviz_common/panel.hpp>

namespace {

constexpr char kPluginId[] = "raystar_rviz_plugins/RaystarPanel";
constexpr char kPluginType[] = "raystar_rviz_plugins::RaystarPanel";

bool isWithinPrefix(const std::filesystem::path& path, const std::filesystem::path& prefix) {
  const std::filesystem::path relative = path.lexically_relative(prefix);
  return !relative.empty() && !relative.is_absolute() && relative.begin() != relative.end() &&
         *relative.begin() != "..";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: plugin_loader_smoke <expected-install-prefix>\n";
    return 2;
  }

  try {
    const std::filesystem::path expected_prefix = std::filesystem::canonical(argv[1]);
    QApplication application(argc, argv);
    pluginlib::ClassLoader<rviz_common::Panel> loader("rviz_common", "rviz_common::Panel");
    const auto declared = loader.getDeclaredClasses();
    if (std::count(declared.begin(), declared.end(), kPluginId) != 1 ||
        !loader.isClassAvailable(kPluginId) || loader.getClassType(kPluginId) != kPluginType) {
      std::cerr << "canonical Raystar Panel registration is unavailable\n";
      return 3;
    }

    const std::filesystem::path library_path =
      std::filesystem::canonical(loader.getClassLibraryPath(kPluginId));
    if (!isWithinPrefix(library_path, expected_prefix)) {
      std::cerr << "plugin resolved outside relocated prefix: " << library_path << '\n';
      return 4;
    }

    auto panel = loader.createSharedInstance(kPluginId);
    if (!panel) {
      std::cerr << "pluginlib returned a null Raystar Panel\n";
      return 5;
    }
    std::cout << "loaded=" << kPluginId << " library=" << library_path << '\n';
  } catch (const std::exception& exception) {
    std::cerr << "plugin smoke failed: " << exception.what() << '\n';
    return 6;
  }
  return 0;
}
