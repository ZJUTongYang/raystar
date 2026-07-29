# One bounded search-tree debug record. Keeping the related values together
# prevents clients from having to align parallel arrays.
int32   index
float64 g_cost  # Accumulated cost in metres
float64 f_cost  # G + H cost in metres
