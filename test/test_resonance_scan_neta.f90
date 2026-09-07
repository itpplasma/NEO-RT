program test_resonance_scan_neta
    ! The optional eta resolution is diagnostic-only.  Keep its default
    ! backward-compatible while allowing a reproducible refinement census.
    use diag_resonance_scan, only: default_neta, resolve_neta
    implicit none

    if (default_neta /= 180) error stop "resonance scan default NETA changed"

    if (resolve_neta() /= default_neta) error stop "missing NETA did not use default"
    if (resolve_neta(720) /= 720) error stop "explicit NETA was not preserved"

    print *, "test_resonance_scan_neta PASSED"
end program test_resonance_scan_neta
