# Machine-readable outcome of one planning request. `message` in the enclosing
# Service/Action result is diagnostic text only and must not be parsed by
# clients to determine success, completeness, or resource limits.

uint8 STATUS_UNKNOWN=0
uint8 STATUS_COMPLETE=1               # Top-K satisfied, or a bounded search was exhausted
uint8 STATUS_FEWER_PATHS=2            # Search completed; fewer than K paths exist
uint8 STATUS_NO_PATH=3                # Search completed without finding a path
uint8 STATUS_PARTIAL_SEARCH=4          # Search stopped before its normal completion
uint8 STATUS_PARTIAL_OUTPUT=5          # Some metric-eligible paths could not be returned
uint8 STATUS_INVALID_REQUEST=6
uint8 STATUS_INVALID_CONFIGURATION=7
uint8 STATUS_BUSY=8
uint8 STATUS_CANCELLED=9
uint8 STATUS_FAILED=10

# `limits_reached` is a bitmask because a partial search can subsequently also
# encounter an output limit. LIMIT_NONE is represented by zero.
uint16 LIMIT_NONE=0
uint16 LIMIT_TIMEOUT=1
uint16 LIMIT_MAX_NODES=2
uint16 LIMIT_MAX_PATH_POINTS=4
uint16 LIMIT_MAX_RESPONSE_BYTES=8
uint16 LIMIT_CANCELLED=16
uint16 LIMIT_MAX_PATHS=32

unique_identifier_msgs/UUID map_id # Exact cached/request map used by this result
unique_identifier_msgs/UUID environment_id # Content/policy/semantics identity used by this result
uint8  status
uint16 limits_reached
bool   request_satisfied     # Top-K: returned K; bounded: cost_bound_exhausted && output_complete
bool   search_complete       # K found, frontier exhausted, or the requested cost bound exhausted
bool   output_complete       # Every authoritative metric-eligible path was returned
bool   debug_requested       # Request explicitly opted in to debug-node output
bool   debug_output_complete # No debug requested, or every expanded node was included
uint32 requested_path_count
uint32 found_path_count      # Serialized-world-topology-certified eligible paths before resource filtering
uint32 returned_path_count   # Must equal path_results.size()
uint64 expanded_nodes
float64 map_time_ms
float64 plan_time_ms
uint8   search_mode              # Echoes the request mode; zero when no request metadata was available
float64 requested_max_path_length # Metres; zero for TOP_K
bool    cost_bound_exhausted     # Search proof only; combine with output_complete to know every bounded path was delivered
