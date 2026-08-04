# Identity and provenance of the server's currently admitted cached map.
std_msgs/Header               header
unique_identifier_msgs/UUID  map_id
uint64                        generation
uint32                        width
uint32                        height
float32                       resolution
# Startup-only occupancy policy and compile-time planning semantics. These
# fields explain how environment_id values returned by planning Actions were
# derived; allow_unknown remains a per-request policy.
uint32                        occupied_threshold
uint32                        environment_identity_version
uint32                        occupancy_semantics_version
uint32                        geometry_semantics_version
uint32                        topology_semantics_version
