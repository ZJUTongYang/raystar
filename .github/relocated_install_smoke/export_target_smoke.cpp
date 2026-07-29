#include <type_traits>

#include <QApplication>
#include <rviz_common/panel.hpp>

#include "raystar_rviz_plugins/raystar_panel.h"

static_assert(std::is_base_of_v<rviz_common::Panel, raystar_rviz_plugins::RaystarPanel>);

int main(int argc, char** argv) {
  QApplication application(argc, argv);
  raystar_rviz_plugins::RaystarPanel panel;
  return panel.metaObject() == nullptr ? 1 : 0;
}
