# The exported Panel target has Qt5::Widgets in its link interface.  A plain
# ament_export_dependencies(Qt5) cannot carry CMake components, so load the
# required component explicitly before downstream consumers import the target.
find_package(Qt5 5.15 QUIET REQUIRED COMPONENTS Widgets)
