module neort_gc_physical_return_contract
    !! Typed boundary between numerical return evidence and an independent
    !! exactly-two theorem/provider.
    implicit none
    private

    type, public :: gc_cylindrical_physical_return_certificate_t
        integer :: certificate_id = 0
        integer :: crossing_count = 0
        logical :: exactly_two_proved = .false.
    end type gc_cylindrical_physical_return_certificate_t

end module neort_gc_physical_return_contract
