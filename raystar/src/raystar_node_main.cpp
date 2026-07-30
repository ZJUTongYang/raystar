#include "raystar_node.h"
#include <rclcpp/rclcpp.hpp>

#include <exception>

int main(int argc, char** argv) {
  rclcpp::init(argc, argv);
  int exit_code = 0;
  try {
    auto node = std::make_shared<raystar::RaystarNode>();
    rclcpp::spin(node);
  } catch (const std::exception& exception) {
    RCLCPP_FATAL(rclcpp::get_logger("raystar"), "Failed to initialize: %s", exception.what());
    exit_code = 1;
  }
  rclcpp::shutdown();
  return exit_code;
}
