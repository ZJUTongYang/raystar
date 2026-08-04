# Raystar Cost-Bounded Non-Homotopic Path Enumeration 实现方案

> 目标仓库：`src/raystar`（相对于 ROS 2 workspace 根目录）
> 当前基线：branch `ros2`，commit `955158baeaf466f17da5ae21334d461b576d474f`
> ROS 2 workspace：包含上述仓库的工作区根目录
> 目标版本：`raystar` / `raystar_interfaces` / `raystar_rviz_plugins` `0.2.0`

## 实现快照（2026-08-03）

历史整理前的 `0.2.0` working-tree 快照已经包含：single-goal cost-bounded enumeration、
shared-tree multi-goal enumeration、基于 free-triangle portal sleeve 与 funnel 的
UPS、content/policy/semantics environment identity、协作式 timeout/cancellation、
完整 radical-sum 长度证书，以及按最终公开 binary64 cost 的稳定排序。原有
top-K overload 与零初始化后的请求源码语义保留；但 rosidl schema 已变化，
server、Panel 和所有客户端仍必须作为同一 `0.2.0` overlay 重新编译和部署。

第二个全新 ROS 2 Humble clean overlay 使用 `Release`、`-Werror` 和单 worker
完成九包构建。原始测试汇总为 962 条记录，可分解为 919 个 detailed testcases、
27 个关联 CTest wrappers 和 16 个 CTest-only checks；结果为 0 error、0 failure、
2 skipped。其中 Raystar 三包的 detailed testcases 分别为 `raystar` 277、
`raystar_interfaces` 9、`raystar_rviz_plugins` 10。真实 ROS 进程测试还覆盖了
single/multi-goal、UPS、timeout 与精确 UUID cancellation。该证据只验证当前组合
快照；下文历史拆分产生的每个中间 commit 仍须重新 clean-build，不能追溯性地把
最终快照的结果宣称为每个中间 commit 的独立测试结果。

该次验证仍是 dirty-tree development evidence，并非 accepted release artifact。
Jazzy、sanitizer、relocated-install 和发布来源 pin 等门槛必须由各自证据满足，
不能由上述 Humble clean build 替代。

## Multi-goal 原生扩展（2026-08-02）

在 single-goal bounded 功能之上，Raystar Core 已增加
`planToGoalsWithinCosts()`：输入同一个 start、按序排列的 goals，以及每个 goal
各自的 inclusive cost bound。它不是在 ROS 层循环调用 single-goal planner，而是
复用一棵 topology tree。共享候选优先级为：

```text
min_i(G(v) + distance(v, goal_i) - bound_i)
```

当共享 open frontier 的最小值严格大于零时，所有仍 active 的 goal 都获得 bound
穷尽证书。一个 goal 达到 `max_cost_bounded_paths` 时只退出该 goal；其它 goal
继续使用已生成的 tree。不可达 goal 在一次 flood-fill 后单独返回 `no_path`，不会
使同批次中的可达 goal 失败。single-goal bounded overload 已反向复用该实现，
top-K 仍保持原搜索路径。

ROS 2 新增独立的 `PlanRaystarGoalSet` Action，是因为其数组输入和逐 goal 输出与
single-goal wire schema 本质不同；它仍与 `PlanRaystarPaths` 和兼容 Service 共享
map cache、worker、cancel state、capacity-one admission 和 Core instance。因此，
下文“扩展现有接口，不新增平行 Action”的决策仅指为 single-goal 增加 bounded
mode 时不另建 Action，并不否定 multi-goal 的原生数组接口。

## UPS 与公开数值契约（2026-08-03）

UPS 不把输出限制在 triangular mesh edges 上。Raystar 先把 topology reference
trace 为保留重复 occurrence 的有向 portal sleeve，再在 lifted sleeve 内运行 funnel；
结果可以穿过 triangle interiors。`PathResult.topology_path` 是 UPS reference，dense
`PathResult.path` 仅用于执行与显示。相同 configuration 的有向 pair 是合法的
zero-length identity transition。

公开路径 cost 不是逐段 `hypot` ceiling 的普通浮点和。实现对 binary64 endpoints
构造 exact-dyadic squared segment lengths，对完整 radical sum 做自适应区间比较，
并 fail closed 于无法证明的边界。最终结果按该公开 cost stable-sort，同时对 source
solution indices 与 visualization association 应用同一 permutation；相等公开 cost
保留 Core topology order。

## 1. 功能定义

在保留现有 top-K 行为的同时，为 Raystar 增加第二种一等搜索模式：

```text
给定 start、goal 和最大路径长度 B，
返回所有长度 <= B、拓扑互异的 locally shortest paths，
并明确报告该长度阈值内的搜索是否已经完备穷尽。
```

这里的“所有 paths”严格指：

- Raystar 当前地图与几何模型中的 locally shortest paths；
- non-homotopic 的定义沿用当前 Raystar tree representation；
- 是否允许 self-crossing 由已有 `allow_self_crossing` 控制；
- 路径长度使用当前 `PathResult.cost` 的欧氏 polyline length；
- ROS API 中单位为米，Core API 中单位为 grid coordinate。

它不表示枚举同一同伦类中的任意非最短曲线。

## 2. API 决策

### 2.1 扩展现有接口，不新增平行 Action

扩展：

- `raystar_interfaces/action/PlanRaystarPaths.action`
- `raystar_interfaces/srv/GetRaystarPaths.srv`
- `raystar_interfaces/msg/PlanningResultInfo.msg`

理由：现有 Service 和 Action 已共用 `executePlanning()`，并具有成熟的 map identity、capacity-one admission、cancellation、资源限制、结果序列化和 RViz 缓存语义。新增另一套 Action 会复制 worker 和取消状态机，并增加两类 Action 竞争同一个 stateful Core 的复杂度。

### 2.2 Request schema

在 Action goal 和 Service request 中加入相同字段：

```text
uint8 SEARCH_MODE_TOP_K=0
uint8 SEARCH_MODE_ALL_WITHIN_LENGTH=1

uint8   search_mode
float64 max_path_length  # ROS API 中为 metres；仅 bounded mode 使用
```

现有 `k` 字段保留。请求合法组合为：

| mode | `k` | `max_path_length` | 语义 |
|---|---:|---:|---|
| `TOP_K` | `> 0` | `0.0` | 完全保持当前行为 |
| `ALL_WITHIN_LENGTH` | `0` | finite 且 `> 0` | 返回阈值内全部 paths |

其它组合直接返回 `STATUS_INVALID_REQUEST`。不采用“忽略无关字段”的宽松策略，避免客户端同时给 K 和 bound 时产生歧义。

新增字段放在现有 goal/request 字段之后。旧客户端源码在针对新接口重新编译后，由于新字段默认零初始化，仍进入 `TOP_K`；但 rosidl schema 已改变，旧二进制与新 server 不保证 wire compatibility。

#### Zero-bound identity 的边界

公开 ROS single-goal、goal-set Action 和兼容 Service 继续要求 bounded-mode
`max_path_length` 为 finite 且严格大于零。这里的零同时是 top-K request metadata 的
sentinel，因此 ROS wire contract 不用零表达 bounded request。

Core contract 则有意允许 finite、non-negative 的 `max_path_cost`。当 bound 为零时，
只有 start 与 goal 完全相同产生的 zero-length identity path 可以被接纳；不同端点的
任何正长度路径都在界外，搜索可以返回经过证明的空集合。该规则使 shared-tree Core
能够精确表达 identity configuration，但它不放宽 ROS 请求校验，也不应与 UPS 的
same-configuration zero-length transition 混为一谈。

### 2.3 Result schema

在 `PlanningResultInfo.msg` 追加：

```text
uint8   search_mode
float64 requested_max_path_length  # metres；top-K 时为 0
bool    cost_bound_exhausted
```

语义：

- `cost_bound_exhausted=true`：Core 已证明不存在尚未返回、长度 `<= bound` 的 path。
- 该字段描述搜索是否穷尽，不描述 ROS output 是否完整。
- 若搜索已穷尽但 ROS response budget 省略了路径，可以出现：

```text
cost_bound_exhausted = true
output_complete = false
request_satisfied = false
status = STATUS_PARTIAL_OUTPUT
```

`request_satisfied` 的新定义：

- top-K：保持当前定义，即实际返回数等于请求 K；
- bounded：`cost_bound_exhausted && output_complete`；允许零路径时为 true；
- `success` 继续表示至少返回一条路径，所以“阈值内没有路径”可以是 `success=false`、`request_satisfied=true`、`STATUS_NO_PATH`。

### 2.4 Status 复用

不新增高层 status：

| bounded outcome | status |
|---|---|
| 阈值穷尽，返回至少一条 | `STATUS_COMPLETE` |
| 阈值穷尽，但零条 | `STATUS_NO_PATH` |
| timeout / max nodes / max paths / cancel | `STATUS_PARTIAL_SEARCH` 或 `STATUS_CANCELLED` |
| 搜索完整、ROS 输出被截断 | `STATUS_PARTIAL_OUTPUT` |

新增 limit bit：

```text
uint16 LIMIT_MAX_PATHS=32
```

## 3. Core API 设计

### 3.1 搜索目标类型

在 `raystar/include/raystar/raystar_core.h` 增加：

```cpp
enum class SearchMode {
  top_k,
  all_within_cost,
};

struct SearchObjective {
  SearchMode mode{SearchMode::top_k};
  int k{1};
  double max_path_cost{0.0};  // Core/grid units；bounded mode 允许 identity bound 0

  static SearchObjective topK(int requested_k);
  static SearchObjective allWithinCost(double maximum_cost);
};

enum class PlanningCompletion {
  none,
  requested_k_reached,
  frontier_exhausted,
  cost_bound_exhausted,
};
```

`PlanResult` 增加：

```cpp
PlanningCompletion completion{PlanningCompletion::none};
```

保留全部现有 `plan(..., int K, ...)` overload，并让它们包装：

```cpp
SearchObjective::topK(K)
```

另增加接收 `SearchObjective` 的 overload。这样 Core 的现有源码调用者不需要修改，新的 tether GCP 可以直接使用 bounded objective。

### 3.2 资源限制

在 `PlanningLimits` 增加：

```cpp
size_t max_cost_bounded_paths = 1000;
```

在 ROS node 增加同名动态参数，并给出 `[1, INT_MAX]` descriptor。它只限制 bounded mode 的已发现 paths 数量：

- 达到该值时，如果 frontier lower bound 已经大于 B 或 frontier 已空，则仍可成功完成；
- 否则以 `PlanningLimitReached::max_paths` 停止，保留已找到前缀并令 `cost_bound_exhausted=false`。

该限制不能被当成正常搜索终止条件。它只是一项内存/响应保护措施。

## 4. 完备停止条件

### 4.1 当前队列为什么能提供证书

当前候选优先级为：

```text
F = G + H
H = EuclideanDistance(candidate_endpoint, goal)
```

欧氏启发式对欧氏 polyline cost 是 admissible 且 consistent。只要继续保持以下不变量：

```text
child.F >= parent.F
```

priority queue 的最小 `F` 就是所有尚未展开分支可达到路径长度的下界。因此：

```text
queue.min_F > B
```

意味着任何尚未发现的 path 都严格长于 B。

### 4.2 搜索循环修改

修改 `raystar/src/raystar_core.cpp` 主循环。逻辑顺序应为：

```cpp
while (!queue.empty()) {
  if (stop_token.poll())
    return stop_for_request();

  const double frontier_lower_bound = queue.front().Fcost_;

  if (objective.mode == SearchMode::all_within_cost &&
      frontier_lower_bound > objective.max_path_cost) {
    result.completion = PlanningCompletion::cost_bound_exhausted;
    break;
  }

  if (objective.mode == SearchMode::all_within_cost &&
      result.path_solutions.size() >= limits.max_cost_bounded_paths) {
    return stop_for_limit(PlanningLimitReached::max_paths);
  }

  if (N_.size() >= limits.max_nodes)
    return stop_for_limit(PlanningLimitReached::max_nodes);

  Candidate best_candidate = pop_min(queue);
  expand(best_candidate);

  if (goal_is_visible) {
    const double path_cost = ...;
    if (objective.mode == SearchMode::top_k ||
        path_cost <= objective.max_path_cost) {
      append_solution();
    }

    if (objective.mode == SearchMode::top_k &&
        result.path_solutions.size() >= objective.k) {
      result.completion = PlanningCompletion::requested_k_reached;
      break;
    }
  }
}

if (queue.empty() && result.completion == PlanningCompletion::none)
  result.completion = PlanningCompletion::frontier_exhausted;
```

判断必须使用 `>`，因为 bound 是 inclusive，`F == B` 的 candidate 仍需展开。

同时对最终 `new_path_length` 再做一次 `<= B` 检查。理论上，在 consistent lower-bound invariant 下，超过 B 的 solution 不应被加入；显式检查可防止将未来算法改动中的 heuristic/rounding regression 变成错误结果。

### 4.3 完成语义

bounded mode 中以下两种 completion 都令 `cost_bound_exhausted=true`：

- `cost_bound_exhausted`：frontier 尚有候选，但其最小下界已超过 B；
- `frontier_exhausted`：不存在任何剩余候选，因此任意有限 B 都自然穷尽。

以下情况必须为 false：

- timeout；
- cancellation；
- max nodes；
- max cost-bounded paths；
- Core `max_path_points`；
- 几何计算失败或无效输入。

### 4.4 数值边界

- ROS 输入 `max_path_length` 必须 finite 且 `> 0`；Core 的 zero-bound identity 是上文
  明确限定的内部 API 语义。
- Node 先构造能覆盖 world/grid rounding contraction 的 Core search-superset bound，
  不能把一次 rounded division 当成完备搜索边界。
- metric/grid product 使用 exact-dyadic comparison 校准；serialized topology 与 dense
  path 的最终长度则由完整 radical-sum certificate 独立证明不超过原始 metric bound。
- equality、`std::nextafter(bound, ±inf)`、subnormal、signed zero、不同 rounding mode
  以及非常大/非常小 resolution 都必须由测试覆盖。
- 不引入任意、用户不可见的 epsilon。无法在有界 refinement precision 内证明的
  inclusive comparison 必须 fail closed，并且不得宣称 `cost_bound_exhausted`。
- 最终 `PathResult.cost` 是同时覆盖 topology/dense serialized polylines 的最小可证明
  binary64 ceiling；发布前按该值稳定排序，而不是假定 Core grid cost 顺序永远相同。

## 5. ROS Node 修改

### 5.1 请求解析

在 `raystar/src/raystar_node.cpp` 增加统一 helper：

```cpp
template <typename RequestT>
bool parseSearchObjective(const RequestT& request,
                          double resolution,
                          const PlanningLimits& limits,
                          SearchObjective& objective,
                          std::string& error);
```

`executePlanning()` 的变化：

1. 先验证 mode 及其 scalar 字段；
2. map 转换完成并取得 resolution 后生成 Core objective；
3. 调用新的 `core_.plan(..., objective, ...)`；
4. 根据 `PlanResult.completion` 填充 result info；
5. 用 mode-aware helper 计算 `request_satisfied`；
6. 日志打印 mode、K 或 bound，而不是始终打印 K。

### 5.2 初始化与错误路径

当前 `initializePlanningResponse()`、exception/cancel branches 都只接收 K。重构为接收小型 metadata：

```cpp
struct RequestResultMetadata {
  uint8_t search_mode;
  uint32_t requested_path_count;
  double requested_max_path_length;
  bool debug_requested;
};
```

确保 invalid request、busy、cancel、exception、stale map 等所有早退路径也能保持一致的 mode metadata；没有可信请求时使用零初始化。

### 5.3 capacity 与 cancellation

不改变现有：

- 一个 `RaystarCore`；
- Service 与 Action 共享 `planning_busy_`；
- capacity one；
- Action worker cooperative cancellation；
- cached-map `map_id` 校验。

bounded enumeration 可能比 top-K 更长，因此 cancellation 和 timeout 测试是发布阻断项。

## 6. RViz Panel

修改：

- `raystar_rviz_plugins/include/.../raystar_panel.h`
- `raystar_rviz_plugins/src/raystar_panel.cpp`
- `raystar_rviz_plugins/test/test_raystar_panel_plugin.cpp`

UI 增加：

- `QComboBox`：`Top K` / `All within length`；
- `QDoubleSpinBox`：最大路径长度（m）；
- mode 切换时只启用相关输入：top-K 启用 K，bounded 启用 length；
- bounded 完成时显示 `Bound exhausted: yes/no`；
- 保存/恢复 RViz config 中的 mode 和 bound。

Panel 不根据人类可读 `message` 判断完成状态，只读取 `PlanningResultInfo`。

初次实现不强制增加 Action feedback；如果 paper-scale GCP 的单次请求明显超过当前交互容忍时间，再增加节流的 `found_path_count / expanded_nodes / frontier_lower_bound` feedback。

## 7. 文件级修改清单

| 文件 | 修改 |
|---|---|
| `raystar_interfaces/action/PlanRaystarPaths.action` | 新增 mode、bound 及注释 |
| `raystar_interfaces/srv/GetRaystarPaths.srv` | 与 Action 保持相同搜索输入 |
| `raystar_interfaces/msg/PlanningResultInfo.msg` | mode metadata、bound certificate、`LIMIT_MAX_PATHS` |
| `raystar_interfaces/test/test_interface_schema.py` | 锁定新 schema 和常量 |
| `raystar/include/raystar/planning_limits.h` | `max_cost_bounded_paths` |
| `raystar/include/raystar/raystar_core.h` | `SearchObjective`、`PlanningCompletion`、overloads |
| `raystar/src/raystar_core.cpp` | objective validation、bound stopping、completion reason、max paths limit |
| `raystar/src/raystar_node.h/.cpp` | request parsing、单位转换、结果映射、参数、日志 |
| `raystar/test/test_raystar_core.cpp` | Core 边界、完备性与限制测试 |
| `raystar/test/test_integration.cpp` | Action/Service bounded workflows、取消和动态参数测试 |
| `raystar_rviz_plugins/.../raystar_panel.*` | mode/bound UI 与结果显示 |
| `raystar_rviz_plugins/test/test_raystar_panel_plugin.cpp` | UI、配置持久化、goal schema 测试 |
| `raystar/launch/raystar_demo.launch.py` | 暴露 `max_cost_bounded_paths` |
| `raystar/README.md`、根 `README.md`、`LAUNCH.md` | API、示例、限制、迁移文档 |
| 三个 `package.xml` | 版本同步升级至 `0.2.0` |

## 8. 测试矩阵

### 8.1 Core 单元测试

1. **Open map / bound below direct path**：零结果，`cost_bound_exhausted`，`no_path`。
2. **Bound exactly equal to direct path**：包含该路径。
3. **One obstacle / bound between first and second classes**：只返回第一类且证书为 true。
4. **Bound covering both classes**：返回两类，排序不变。
5. **Bound above all existing classes**：以 `frontier_exhausted` 完成。
6. **Bounded result equals high-K prefix**：相同路径 digest、cost 和顺序。
7. **NaN / ±inf / negative bound**：Core invalid；ROS 另外拒绝 zero bound。
8. **Core zero bound**：相同端点只返回 identity path，不同端点返回 certified empty set。
9. **max paths reached**：保留前缀，`max_paths`，证书 false。
10. **max nodes / timeout / cancel / max path points**：证书 false，后续正常请求不受污染。
11. **self-crossing true/false**：分别与同策略的 high-K baseline 一致。
12. **lower-bound monotonicity**：所有新入队 child 的 F 不小于被展开 parent F（允许明确定义的 ULP 容差）。
13. **determinism**：重复运行 path digest、expanded nodes、completion reason 一致。

### 8.2 ROS 集成测试

- Action 和 legacy Service 的 bounded 结果一致。
- ROS bound 单位为米；在非 1.0 resolution 下验证转换。
- zero-initialized `search_mode` 的现有 top-K tests 全部不变。
- stale map ID、busy、invalid mode/field combination。
- bounded Action cancel 后 terminal state 为 CANCELED，证书 false。
- output budget 截断时搜索证书保留，但 `request_satisfied=false`。
- 动态修改 `max_cost_bounded_paths`：新请求读取新 snapshot，进行中的请求不变。
- marker cache 从多路径 bounded 结果切回 top-K 时正确清理旧 namespace。

### 8.3 性能与回归

- 现有 24 个 `(scenario,K)` profiling cases 结果 digest 和 expanded nodes 不变。
- profiling runner 增加可选 `--max-path-length`，输出：bound、completion reason、bound exhausted、found count。
- 在 `single_obstacle_256`、`narrow_gate_256`、`dense_lattice_192` 上选择多个 bound，确认结果数随 bound 单调不减。
- bounded mode 达到相同最终 path set 时，结果必须等于足够大 K 的现有模式；运行时间允许不同。

## 9. 五提交历史整理方案（2026-08-03）

这不是重新设计功能，而是把已经组合验证的大型 dirty tree 整理为最少且可审阅的
本地历史。整理应在不推送的本地 feature branch 上进行；每个 implementation commit
必须同时包含自己的 tests，并在形成后独立 clean-build。最终 tree 必须与整理前快照
一致，除非审计明确记录了额外修正。

### Commit 1：公开 interfaces contract（breaking）

- bounded single-goal Action/Service fields 与完成证书；
- shared-tree goal-set Action 及逐 goal result messages；
- UPS transition Action/messages、`topology_path` 与 environment identity fields；
- rosidl generation、schema contract tests 和 `raystar_interfaces` `0.2.0` metadata。

这是 source-additive 但 wire-incompatible 的 schema commit。commit message 必须包含
`BREAKING CHANGE`，并说明 `0.1.x` binaries 不得与 `0.2.0` server 混用。

### Commit 2：Core bounded/shared-tree planning 与 portal-sleeve UPS

- `SearchObjective`、per-goal completion、resource limits 与旧 top-K wrappers；
- exact metric-bound/radical-sum certificates；
- multi-endpoint Polymap construction、free-triangle trace、portal witness 与 funnel；
- Core/geometry/numerical unit tests，以及对应 Boost build dependency。

bounded、multi-goal 与 UPS 已交织在 `raystar_core.cpp`、`polymap.cpp` 和公共测试夹具中；
为减少构造未经验证中间状态的风险，不再强行拆成多个 Core commits。

### Commit 3：ROS node、launch 与 RViz integration

- 三类 Action、兼容 Service、shared capacity-one worker 与 cancellation/timeout mapping；
- environment cache/identity、metric admission、public-cost ordering 与结果净化；
- launch parameters、RViz bounded-mode controls、integration/plugin tests；
- `raystar` / `raystar_rviz_plugins` interface constraints 与 `0.2.0` metadata。

### Commit 4：Profiling、acceptance 与 relocatable resource lookup

- top-K、bounded 与 multi-goal profiling modes 及 schema-aware summarizer；
- failure/smoke/summary contracts；
- 通过 `ament_index_cpp` 查找 installed bundled map，禁止 source-tree path fallback。

resource relocation 与算法语义正交，但它和新 profiling schema 已在同一工具文件中
高度交织；保留为一个 tooling commit 比再次进行脆弱 hunk 拆分更安全。

### Commit 5：Documentation 与 migration record

- 根 README、package README、launch/compatibility 文档与本实现记录；
- 明确 0.2.0 wire migration、zero-bound identity 边界和历史 profiling evidence boundary；
- 不把 dirty-tree 或最终组合测试结果伪装成中间 commits 的独立证据。

所有提交只保留在本地，未经用户单独授权不得 push。五个提交完成后应逐项运行
`git diff --check`、schema/Core/ROS/RViz/profiling tests，并对比最终 tree 与整理前
snapshot，确认没有遗失未跟踪的新接口、helper 或 regression tests。

## 10. 验收标准

功能完成必须同时满足：

1. top-K 模式所有既有行为、结果 digest 和测试保持不变；
2. bounded mode 只返回 `cost <= max_path_length` 的路径；
3. 返回路径仍然 collision-free、按 cost 非递减且 topology-distinct；
4. `cost_bound_exhausted=true` 只在数学停止条件或 frontier exhaustion 下出现；
5. 任何资源限制、取消或 Core failure 都不会伪造完备证书；
6. bounded 结果与足够大 K 结果的阈值内前缀完全一致；
7. ROS metre/grid conversion 在非单位 resolution 下通过边界测试；
8. Action、Service、RViz、launch 和文档使用相同语义；
9. Humble/Jazzy CI、sanitizers 和 relocated-install tests 全部通过。

满足以上条件后，tether workspace convexity repository 可以把：

```text
search_mode = ALL_WITHIN_LENGTH
max_path_length = maximum_tether_length
allow_self_crossing = false
```

作为 GCP 请求，并把 `cost_bound_exhausted && output_complete` 作为 goal configuration 集合完整性的证书。
