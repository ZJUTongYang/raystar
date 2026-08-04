# Raystar dependency and platform contract

This file defines the supported source-build and binary-runtime combinations
for `raystar`, `raystar_interfaces`, and `raystar_rviz_plugins`. Patch releases
from the ROS and Ubuntu repositories are intentionally not pinned; CI records
the exact package versions used by every run.

## Supported matrix

| Support level | OS | ROS 2 | Toolchain baseline | Geometry | RViz plugin |
|---|---|---|---|---|---|
| Minimum | Ubuntu 22.04 | Humble | GCC 11, CMake 3.22, C++17 | CGAL 5.4, GMP 6.2.1, MPFR 4.1.0 | Qt 5.15, `rviz_common` 11.2, pluginlib 5.1 |
| Current LTS | Ubuntu 24.04 | Jazzy | GCC 13, CMake 3.28, C++17 | CGAL 5.6.1 with distribution GMP/MPFR | Qt 5.15 with Jazzy `rviz_common` and pluginlib |

Both rows are required lanes in `.github/workflows/ci.yml` and run the complete
build and test suite. Humble/Jammy is the minimum supported platform; Jazzy/Noble
is the newer compatibility lane.

CGAL versions newer than 5.6.1, including CGAL 6.x, are not rejected by CMake
when they satisfy the `>= 5.4` source contract, but they are not currently a
CI-backed release claim. Qt 6, ROS 2 Kilted/Rolling, non-Linux platforms, and
32-bit builds are unverified.

## 0.2.0 interface migration

Raystar `0.2.0` changes the rosidl wire schema. Existing
`PlanRaystarPaths.action`, `GetRaystarPaths.srv`, `PlanningResultInfo.msg`,
`PathResult.msg`, and `MapStatus.msg` records gained fields, and `0.2.0` adds
the goal-set and transition-graph Actions with their supporting messages.
Appending a field is still a DDS/rosidl type change: a `0.1.x` client binary
must not communicate with a `0.2.0` server binary.

Rebuild and deploy `raystar_interfaces`, `raystar`, `raystar_rviz_plugins`, and
every downstream client in one consistent overlay. Source code that
zero-initializes the new `search_mode` retains top-K request semantics after it
is recompiled; this source-level default is not binary compatibility.
Downstream consumers must also observe these new contracts:

- `PathResult.path` is the dense execution/display path, while
  `PathResult.topology_path` is the unsampled topology reference required by
  UPS and other homotopy consumers.
- `PlanningResultInfo.environment_id` binds map content, occupancy policy, and
  explicit semantic versions. It complements rather than replaces the cached
  message `map_id`.
- `environment_identity.hpp` is installed beside `map_identity.hpp`; clients
  should use its helper instead of independently reproducing the hash layout.
- Public path sequences are ordered by finalized conservative binary64 cost,
  which can differ by one ULP from Core grid-cost order.

## Dependency ownership

| Consumer | Build dependency | Export dependency | Runtime dependency |
|---|---|---|---|
| `raystar_node` / internal `raystar_core` | CGAL >= 5.4, GMP, MPFR, header-only Boost.Multiprecision, and ROS message/client libraries | None; `raystar_core` and its headers are package-private | ROS libraries, GMP and MPFR; CGAL and Boost.Multiprecision are header-only here |
| optional `raystar_profile` | `raystar_core`, `ament_index_cpp`, and the profiling oracle dependencies | None | `ament_index_cpp` and an installed `raystar` package resource index entry are required to locate `share/raystar/test/testmap.pgm` unless `--testmap` is supplied |
| `raystar_interfaces` | rosidl generators and the messages embedded in its IDL | Generated interface targets plus the public `map_identity.hpp` and `environment_identity.hpp` dependencies | rosidl runtime and message type-support libraries |
| `raystar_rviz_plugins` | Qt 5.15 Widgets, RViz, pluginlib and its public ROS dependencies | Installed Panel header, imported library target, Qt Widgets and public ROS dependencies | Qt Core/Widgets, RViz/pluginlib and generated Raystar interface libraries |

The portable rosdep keys include `boost`, `cgal`, `libgmp`, `mpfr`,
`ament_index_cpp`, `qtbase5-dev`, `libqt5-core`, and `libqt5-widgets`. Names
such as `libboost-dev`, `libgmp-dev`, and `libmpfr-dev` are Debian package
names, not portable rosdep keys.

The exact-constructions kernel makes GMP and MPFR real runtime dependencies:
on Jammy, `raystar_node` resolves `libgmp.so.10` and `libmpfr.so.6`. They are
therefore declared for both source compilation and deployment. CGAL remains a
build-only dependency because no Raystar CGAL library or header is installed.
Boost.Multiprecision supplies the exact-dyadic and radical-sum integer
arithmetic and is likewise a header-only build dependency; it introduces no
Boost shared-library runtime requirement. The optional installed profiling
runner uses `ament_index_cpp` at runtime so relocating a complete install
prefix does not leave a compiled-in source-tree map path.

## Enforced checks

- `raystar/CMakeLists.txt` rejects CGAL below 5.4 and CGAL configurations
  without GMP/MPFR, requires Boost for the exact-length implementation, and
  deliberately allows a later CGAL major release.
- `raystar_interfaces/test/test_interface_schema.py` locks the complete 0.2.0
  Action/Service/message field and status-constant contract.
- `raystar_rviz_plugins/CMakeLists.txt` requires Qt 5.15, `rviz_common` 11.2,
  and pluginlib 5.1. Its package config explicitly loads the Qt Widgets
  component before exposing the imported Panel target.
- The profiling relocation contract rejects `RAYSTAR_PROFILE_TESTMAP_PATH`,
  `/proc/self/exe`, and any embedded source-tree map path; normal installed-map
  lookup goes through `ament_index_cpp`, while `--testmap` remains an explicit
  override.
- The normal CI matrix builds and runs all tests on Humble/Jammy and
  Jazzy/Noble.
- The relocated-install lane performs a fresh non-symlink build, moves the
  complete install prefix, and then uses only the moved setup to verify package
  discovery, runtime linking, an out-of-tree exported-target consumer, actual
  pluginlib loading, direct node startup, and a real headless launch.

A developer `--symlink-install` tree is intentionally not relocatable. Release
or deployment checks must omit `--symlink-install` and move the complete install
tree rather than an individual package prefix.
