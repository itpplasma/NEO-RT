module neort_gc_operational_class_adapter_bridge
    !! Convert one certified operational class assembly into adapter classes.
    !!
    !! The operational assembly owns the ordered class topology and its
    !! canonical total variations.  This bridge only validates the partition,
    !! transfers that data to the cylindrical adapter representation, and
    !! retains the candidate's root evidence at the two outer boundaries.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use neort_gc_certified_interval_roots, only: gc_interval_t
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t, operator(+)
    use neort_gc_cylindrical_class_adapter, only: &
        gc_cylindrical_class_interval_t
    use neort_gc_operational_class_assembly, only: &
        GC_CLASS_BOUNDARY_REGULAR, GC_CLASS_BOUNDARY_TURNING, &
        GC_CLASS_BOUNDARY_SEPARATRIX, GC_CLASS_BOUNDARY_X, &
        GC_CLASS_ASSEMBLY_SUCCESS, &
        gc_operational_class_assembly_result_t, &
        gc_operational_class_interval_t
    implicit none
    private

    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_SUCCESS = 0
    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_INVALID_INPUT = 1
    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_INVALID_ASSEMBLY = 2
    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR = 3
    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR = 4
    integer, parameter, public :: GC_CLASS_ADAPTER_BRIDGE_ID_OVERFLOW = 5

    ! Component IDs are packed as parent_id*stride+class_ordinal.  The
    ! fixed stride makes ID blocks disjoint for different parent components;
    ! int64 arithmetic is used before conversion to the adapter's integer.
    integer(int64), parameter, public :: &
        GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE = 1000000_int64

    public :: bridge_gc_operational_class_assembly

contains

    subroutine bridge_gc_operational_class_assembly(candidate, assembly, &
            split_classes, status)
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        type(gc_operational_class_assembly_result_t), intent(in) :: assembly
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            split_classes(:)
        integer, intent(out) :: status

        integer :: i, nclasses, allocation_status
        integer(int64) :: packed_id

        allocate(split_classes(0))
        status = GC_CLASS_ADAPTER_BRIDGE_INVALID_INPUT

        if (.not. valid_candidate(candidate)) return
        status = validate_assembly(candidate, assembly)
        if (status /= GC_CLASS_ADAPTER_BRIDGE_SUCCESS) then
            return
        end if
        nclasses = assembly%nclasses
        if (.not. valid_component_id_range(candidate%component_id, nclasses)) then
            status = GC_CLASS_ADAPTER_BRIDGE_ID_OVERFLOW
            return
        end if

        deallocate(split_classes)
        allocate(split_classes(nclasses), stat=allocation_status)
        if (allocation_status /= 0) then
            allocate(split_classes(0))
            status = GC_CLASS_ADAPTER_BRIDGE_INVALID_ASSEMBLY
            return
        end if

        do i = 1, nclasses
            call packed_component_id(candidate%component_id, i, packed_id)
            split_classes(i) = gc_cylindrical_class_interval_t()
            split_classes(i)%component_id = int(packed_id)
            split_classes(i)%sigma = assembly%classes(i)%sigma
            split_classes(i)%ifuntype = assembly%classes(i)%ifuntype
            split_classes(i)%left_boundary_id = &
                assembly%classes(i)%left_boundary_id
            split_classes(i)%right_boundary_id = &
                assembly%classes(i)%right_boundary_id
            split_classes(i)%rc_min = assembly%classes(i)%x_lo
            split_classes(i)%rc_max = assembly%classes(i)%x_hi
            split_classes(i)%psi_star_min = &
                assembly%classes(i)%canonical_lo
            split_classes(i)%psi_star_max = &
                assembly%classes(i)%canonical_hi
            split_classes(i)%canonical_measure = &
                assembly%classes(i)%canonical_total_variation
            split_classes(i)%rc_min_enclosure = gc_interval_t( &
                split_classes(i)%rc_min, split_classes(i)%rc_min)
            split_classes(i)%rc_max_enclosure = gc_interval_t( &
                split_classes(i)%rc_max, split_classes(i)%rc_max)
            split_classes(i)%psi_star_min_enclosure = gc_interval_t( &
                split_classes(i)%psi_star_min, split_classes(i)%psi_star_min)
            split_classes(i)%psi_star_max_enclosure = gc_interval_t( &
                split_classes(i)%psi_star_max, split_classes(i)%psi_star_max)
            split_classes(i)%canonical_measure_enclosure = gc_interval_t( &
                assembly%classes(i)%canonical_measure_enclosure%lo, &
                assembly%classes(i)%canonical_measure_enclosure%hi)
            split_classes(i)%allowed_interval = .true.
            split_classes(i)%topology_certified = .true.
            split_classes(i)%root_isolation_certified = .true.
            split_classes(i)%orbit_return_certified = .false.
            split_classes(i)%lower_boundary_kind = &
                boundary_kind_name(assembly%classes(i)%left_kind)
            split_classes(i)%upper_boundary_kind = &
                boundary_kind_name(assembly%classes(i)%right_kind)
            split_classes(i)%limiting_chart = candidate%limiting_chart
        end do

        ! Root evidence belongs to the candidate's outer endpoints.  Internal
        ! assembly boundaries are separatrix/class metadata, not new root
        ! certificates.
        if (candidate%lower_root) then
            split_classes(1)%lower_root = .true.
            split_classes(1)%lower_tangent = candidate%lower_tangent
            split_classes(1)%lower_root_certificate = &
                candidate%lower_root_certificate
            split_classes(1)%lower_root_coordinate_map = &
                candidate%lower_root_coordinate_map
            split_classes(1)%rc_min_enclosure = candidate%rc_min_enclosure
            split_classes(1)%psi_star_min_enclosure = &
                candidate%psi_star_min_enclosure
            split_classes(1)%lower_boundary_kind = candidate%lower_boundary_kind
        end if
        if (candidate%upper_root) then
            split_classes(nclasses)%upper_root = .true.
            split_classes(nclasses)%upper_tangent = candidate%upper_tangent
            split_classes(nclasses)%upper_root_certificate = &
                candidate%upper_root_certificate
            split_classes(nclasses)%upper_root_coordinate_map = &
                candidate%upper_root_coordinate_map
            split_classes(nclasses)%rc_max_enclosure = candidate%rc_max_enclosure
            split_classes(nclasses)%psi_star_max_enclosure = &
                candidate%psi_star_max_enclosure
            split_classes(nclasses)%upper_boundary_kind = candidate%upper_boundary_kind
        end if
        status = GC_CLASS_ADAPTER_BRIDGE_SUCCESS
    end subroutine bridge_gc_operational_class_assembly

    logical function valid_candidate(candidate)
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate

        valid_candidate = .false.
        if (candidate%component_id <= 0) return
        if (abs(candidate%sigma) /= 1) return
        if (.not. candidate%allowed_interval) return
        if (.not. candidate%topology_certified) return
        if (.not. candidate%root_isolation_certified) return
        if (.not. all(ieee_is_finite([candidate%rc_min, candidate%rc_max, &
                candidate%psi_star_min, candidate%psi_star_max, &
                candidate%canonical_measure]))) return
        if (candidate%rc_max <= candidate%rc_min) return
        if (candidate%canonical_measure <= 0.0_dp) return
        if (.not. valid_interval(candidate%rc_min_enclosure)) return
        if (.not. valid_interval(candidate%rc_max_enclosure)) return
        if (.not. valid_interval(candidate%psi_star_min_enclosure)) return
        if (.not. valid_interval(candidate%psi_star_max_enclosure)) return
        if (.not. valid_interval(candidate%canonical_measure_enclosure)) return
        if (candidate%canonical_measure_enclosure%lo < 0.0_dp) return
        if (candidate%canonical_measure < &
                candidate%canonical_measure_enclosure%lo) return
        if (candidate%canonical_measure > &
                candidate%canonical_measure_enclosure%hi) return
        valid_candidate = .true.
    end function valid_candidate

    integer function validate_assembly(candidate, assembly)
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        type(gc_operational_class_assembly_result_t), intent(in) :: assembly

        integer :: i, j
        real(dp) :: coordinate_tolerance, canonical_tolerance
        type(gc_outward_interval_t) :: measure_sum
        type(gc_outward_interval_t) :: class_measure

        validate_assembly = GC_CLASS_ADAPTER_BRIDGE_INVALID_ASSEMBLY
        if (assembly%status /= GC_CLASS_ASSEMBLY_SUCCESS) return
        if (.not. assembly%complete) return
        if (assembly%nclasses < 1) return
        if (.not. allocated(assembly%classes)) return
        if (size(assembly%classes) /= assembly%nclasses) return

        coordinate_tolerance = 512.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            abs(candidate%rc_min), abs(candidate%rc_max))
        canonical_tolerance = 512.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            abs(candidate%psi_star_min), abs(candidate%psi_star_max), &
            abs(candidate%canonical_measure))
        measure_sum = gc_outward_interval(0.0_dp, 0.0_dp)
        do i = 1, assembly%nclasses
            if (assembly%classes(i)%class_id /= i) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (assembly%classes(i)%component_id /= candidate%component_id) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (assembly%classes(i)%sigma /= candidate%sigma) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (.not. all(ieee_is_finite([assembly%classes(i)%x_lo, &
                    assembly%classes(i)%x_hi, &
                    assembly%classes(i)%canonical_lo, &
                    assembly%classes(i)%canonical_hi, &
                    assembly%classes(i)%canonical_total_variation]))) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (assembly%classes(i)%x_hi <= assembly%classes(i)%x_lo) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (assembly%classes(i)%canonical_total_variation <= 0.0_dp) then
                if (assembly%classes(i)%canonical_measure_enclosure%hi > &
                        0.0_dp) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
                    return
                end if
            end if
            if (.not. assembly%classes(i)%canonical_measure_certified) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
                return
            end if
            if (assembly%classes(i)%canonical_measure_certificate_id <= 0) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
                return
            end if
            if (.not. valid_measure_interval(assembly%classes(i)% &
                    canonical_measure_enclosure)) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
                return
            end if
            if (assembly%classes(i)%canonical_total_variation < &
                    assembly%classes(i)%canonical_measure_enclosure%lo .or. &
                    assembly%classes(i)%canonical_total_variation > &
                    assembly%classes(i)%canonical_measure_enclosure%hi) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
                return
            end if
            if (.not. valid_boundary_kind(assembly%classes(i)%left_kind)) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (.not. valid_boundary_kind(assembly%classes(i)%right_kind)) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (assembly%classes(i)%ifuntype /= &
                    10*assembly%classes(i)%left_kind + &
                    assembly%classes(i)%right_kind) then
                validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                return
            end if
            if (i == 1) then
                if (assembly%classes(i)%left_boundary_id /= 0) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%x_lo-candidate%rc_min) > &
                        coordinate_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%canonical_lo- &
                        candidate%psi_star_min) > canonical_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
            else
                if (assembly%classes(i)%left_boundary_id <= 0) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (assembly%classes(i)%left_boundary_id /= &
                        assembly%classes(i-1)%right_boundary_id) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%x_lo- &
                        assembly%classes(i-1)%x_hi) > coordinate_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%canonical_lo- &
                        assembly%classes(i-1)%canonical_hi) > &
                        canonical_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
            end if
            if (i == assembly%nclasses) then
                if (assembly%classes(i)%right_boundary_id /= 0) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%x_hi-candidate%rc_max) > &
                        coordinate_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
                if (abs(assembly%classes(i)%canonical_hi- &
                        candidate%psi_star_max) > canonical_tolerance) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
            else
                if (assembly%classes(i)%right_boundary_id <= 0) then
                    validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                    return
                end if
            end if
            do j = 1, i-1
                if (i < assembly%nclasses) then
                    if (assembly%classes(i)%right_boundary_id == &
                            assembly%classes(j)%right_boundary_id .or. &
                            assembly%classes(i)%right_boundary_id == &
                            assembly%classes(j)%left_boundary_id) then
                        validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                        return
                    end if
                end if
                if (assembly%classes(i)%left_boundary_id > 0) then
                    if (assembly%classes(i)%left_boundary_id == &
                            assembly%classes(j)%left_boundary_id .and. &
                            j /= i-1) then
                        validate_assembly = GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR
                        return
                    end if
                end if
            end do
            class_measure = gc_outward_interval( &
                assembly%classes(i)%canonical_measure_enclosure%lo, &
                assembly%classes(i)%canonical_measure_enclosure%hi)
            measure_sum = measure_sum + class_measure
        end do
        if (.not. gc_outward_interval_is_valid(measure_sum)) then
            validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
            return
        end if
        if (measure_sum%lo < candidate%canonical_measure_enclosure%lo) then
            validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
            return
        end if
        if (measure_sum%hi > candidate%canonical_measure_enclosure%hi) then
            validate_assembly = GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR
            return
        end if
        validate_assembly = GC_CLASS_ADAPTER_BRIDGE_SUCCESS
    end function validate_assembly

    pure logical function valid_boundary_kind(kind)
        integer, intent(in) :: kind

        valid_boundary_kind = kind >= GC_CLASS_BOUNDARY_REGULAR .and. &
            kind <= GC_CLASS_BOUNDARY_X
    end function valid_boundary_kind

    pure logical function valid_interval(value)
        type(gc_interval_t), intent(in) :: value

        valid_interval = all(ieee_is_finite([value%lo, value%hi])) .and. &
            value%lo <= value%hi
    end function valid_interval

    pure logical function valid_measure_interval(value)
        use neort_gc_certified_interval_roots, only: gc_interval_t
        type(gc_interval_t), intent(in) :: value

        valid_measure_interval = valid_interval(value)
        if (.not. valid_measure_interval) return
        valid_measure_interval = value%lo >= 0.0_dp
    end function valid_measure_interval

    logical function valid_component_id_range(component_id, nclasses)
        integer, intent(in) :: component_id, nclasses
        integer(int64) :: upper, class_count

        valid_component_id_range = .false.
        if (component_id <= 0 .or. nclasses < 1) return
        class_count = int(nclasses, int64)
        if (class_count >= GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE) return
        upper = int(huge(0), int64)
        valid_component_id_range = int(component_id, int64) <= &
            (upper-class_count)/GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE
    end function valid_component_id_range

    subroutine packed_component_id(parent_id, class_ordinal, packed_id)
        integer, intent(in) :: parent_id, class_ordinal
        integer(int64), intent(out) :: packed_id

        packed_id = int(parent_id, int64)* &
            GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE + &
            int(class_ordinal, int64)
    end subroutine packed_component_id

    pure function boundary_kind_name(kind) result(name)
        integer, intent(in) :: kind
        character(len=24) :: name

        name = 'unresolved'
        select case (kind)
            case (GC_CLASS_BOUNDARY_REGULAR)
                name = 'operational-regular'
            case (GC_CLASS_BOUNDARY_TURNING)
                name = 'operational-turning'
            case (GC_CLASS_BOUNDARY_SEPARATRIX)
                name = 'operational-separatrix'
            case (GC_CLASS_BOUNDARY_X)
                name = 'operational-x-point'
        end select
    end function boundary_kind_name

end module neort_gc_operational_class_adapter_bridge
