module neort_gc_wall_context
    !! Configuration-only seam for the production wall path.  The cylindrical
    !! backend copies these values into its data-oriented context; orbit RHS
    !! algebra does not depend on this module.
    implicit none

    character(len=1024), public, save :: configured_wall_file = ''
    character(len=16), public, save :: configured_wall_units = 'm'

end module neort_gc_wall_context
