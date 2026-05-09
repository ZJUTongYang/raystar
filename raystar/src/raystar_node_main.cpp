#include "raystar_node.h"
#include <rclcpp/rclcpp.hpp>

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  auto node = std::make_shared<raystar::RaystarNode>();
  rclcpp::spin(node);
  rclcpp::shutdown();
  return 0;
}
