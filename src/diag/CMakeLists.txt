# Diagnostics static library
add_library(neort_diaglib STATIC
  diag_bounce_nonlin.f90
  diag_atten_map.f90
  diag_contrib_map.f90
  diag_dnorm_breakdown.f90
)

target_link_libraries(neort_diaglib PUBLIC
  neo_rt
  fortplot::fortplot
)

target_include_directories(neort_diaglib PUBLIC
  ${PROJECT_BINARY_DIR}
)

set_target_properties(neort_diaglib PROPERTIES
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/modules"
)
