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

## Dependency ownership

| Consumer | Build dependency | Export dependency | Runtime dependency |
|---|---|---|---|
| `raystar_node` | CGAL >= 5.4, GMP, MPFR and ROS message/client libraries | None; `raystar_core` and its headers are package-private | ROS libraries, GMP and MPFR; CGAL itself is header-only in the supported distributions |
| `raystar_interfaces` | rosidl generators and the messages embedded in its IDL | Generated interface targets plus the public `map_identity.hpp` dependencies | rosidl runtime and message type-support libraries |
| `raystar_rviz_plugins` | Qt 5.15 Widgets, RViz, pluginlib and its public ROS dependencies | Installed Panel header, imported library target, Qt Widgets and public ROS dependencies | Qt Core/Widgets, RViz/pluginlib and generated Raystar interface libraries |

The portable rosdep keys are `cgal`, `libgmp`, `mpfr`, `qtbase5-dev`,
`libqt5-core`, and `libqt5-widgets`. Names such as `libgmp-dev` and
`libmpfr-dev` are Debian package names, not portable rosdep keys.

The exact-constructions kernel makes GMP and MPFR real runtime dependencies:
on Jammy, `raystar_node` resolves `libgmp.so.10` and `libmpfr.so.6`. They are
therefore declared for both source compilation and deployment. CGAL remains a
build-only dependency because no Raystar CGAL library or header is installed.

## Enforced checks

- `raystar/CMakeLists.txt` rejects CGAL below 5.4 and CGAL configurations
  without GMP/MPFR, while deliberately allowing a later CGAL major release.
- `raystar_rviz_plugins/CMakeLists.txt` requires Qt 5.15, `rviz_common` 11.2,
  and pluginlib 5.1. Its package config explicitly loads the Qt Widgets
  component before exposing the imported Panel target.
- The normal CI matrix builds and runs all tests on Humble/Jammy and
  Jazzy/Noble.
- The relocated-install lane performs a fresh non-symlink build, moves the
  complete install prefix, and then uses only the moved setup to verify package
  discovery, runtime linking, an out-of-tree exported-target consumer, actual
  pluginlib loading, direct node startup, and a real headless launch.

A developer `--symlink-install` tree is intentionally not relocatable. Release
or deployment checks must omit `--symlink-install` and move the complete install
tree rather than an individual package prefix.
