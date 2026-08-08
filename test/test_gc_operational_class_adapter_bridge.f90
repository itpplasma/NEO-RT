program test_gc_operational_class_adapter_bridge
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use neort_gc_certified_interval_roots, only: &
        gc_interval_root_box_t, gc_interval_t
    use neort_gc_cylindrical_class_adapter, only: &
        gc_cylindrical_class_interval_t
    use neort_gc_cylindrical_topology, only: &
        gc_cylindrical_root_coordinate_map_t
    use neort_gc_operational_class_adapter_bridge, only: &
        GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE, &
        GC_CLASS_ADAPTER_BRIDGE_ID_OVERFLOW, &
        GC_CLASS_ADAPTER_BRIDGE_INVALID_ASSEMBLY, &
        GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR, &
        GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR, &
        GC_CLASS_ADAPTER_BRIDGE_SUCCESS, &
        bridge_gc_operational_class_assembly
    use neort_gc_operational_class_assembly, only: &
        GC_CLASS_ASSEMBLY_INVALID_INPUT, GC_CLASS_ASSEMBLY_SUCCESS, &
        GC_CLASS_BOUNDARY_REGULAR, GC_CLASS_BOUNDARY_SEPARATRIX, &
        GC_CLASS_BOUNDARY_X, gc_operational_class_assembly_result_t, &
        gc_operational_class_interval_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_cylindrical_class_interval_t) :: candidate
    type(gc_operational_class_assembly_result_t) :: assembly
    type(gc_cylindrical_class_interval_t), allocatable :: split(:)
    integer :: status

    call make_fixture(candidate, assembly)
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_SUCCESS, &
        'validated operational assembly was rejected')
    call require(allocated(split) .and. size(split) == 3, &
        'bridge did not create the expected partition')
    call require(maxval(abs([split%rc_min] - [0.0_dp, 3.0_dp, 7.0_dp])) &
        < 1.0e-12_dp .and. maxval(abs([split%rc_max] - &
        [3.0_dp, 7.0_dp, 10.0_dp])) < 1.0e-12_dp, &
        'physical class endpoints were not mapped')
    call require(maxval(abs([split%psi_star_min] - &
        [2.0_dp, 5.0_dp, 1.0_dp])) < 1.0e-12_dp .and. &
        maxval(abs([split%psi_star_max] - &
        [5.0_dp, 1.0_dp, 8.0_dp])) < 1.0e-12_dp, &
        'canonical class endpoints were not mapped')
    call require(maxval(abs([split%canonical_measure] - &
        [4.0_dp, 7.0_dp, 1.0_dp])) < 1.0e-12_dp, &
        'canonical class measures were not mapped')
    call require(abs(sum(split%canonical_measure) - &
        candidate%canonical_measure) < 1.0e-12_dp, &
        'canonical measure does not partition the candidate')
    call require(all(split%ifuntype == [14, 43, 31]), &
        'POTATO ifuntype metadata was not preserved')
    call require(all(split%left_boundary_id == [0, 11, 12]) .and. &
        all(split%right_boundary_id == [11, 12, 0]), &
        'boundary ownership metadata was not preserved')
    call require(all(split%topology_certified) .and. &
        all(split%root_isolation_certified) .and. all(split%allowed_interval), &
        'bridge did not certify the converted intervals')
    call require(all(.not. split%orbit_return_certified), &
        'bridge claimed an unproved physical orbit return')
    call require(all([split%component_id] == &
        [int(7_int64*GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE+1_int64), &
        int(7_int64*GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE+2_int64), &
        int(7_int64*GC_CLASS_ADAPTER_BRIDGE_COMPONENT_STRIDE+3_int64)]), &
        'component identifiers are not deterministic and unique')

    call require(split(1)%lower_root .and. .not. split(1)%upper_root .and. &
        .not. split(2)%lower_root .and. .not. split(2)%upper_root .and. &
        .not. split(3)%lower_root .and. split(3)%upper_root, &
        'outer root ownership was not restricted to first/last classes')
    call require(split(1)%lower_root_certificate%cut_id == 101 .and. &
        split(3)%upper_root_certificate%cut_id == 202 .and. &
        split(2)%lower_root_certificate%cut_id == 0 .and. &
        split(2)%upper_root_certificate%cut_id == 0, &
        'root certificates were copied to the wrong classes')
    call require(split(1)%lower_root_coordinate_map%map_certificate_id == 301 .and. &
        split(3)%upper_root_coordinate_map%map_certificate_id == 302 .and. &
        split(2)%lower_root_coordinate_map%map_certificate_id == 0 .and. &
        split(2)%upper_root_coordinate_map%map_certificate_id == 0, &
        'root coordinate maps were copied to the wrong classes')

    assembly%status = GC_CLASS_ASSEMBLY_INVALID_INPUT
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_INVALID_ASSEMBLY, &
        'invalid assembly status was accepted')

    call make_fixture(candidate, assembly)
    assembly%classes(2)%x_lo = 3.5_dp
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR, &
        'non-contiguous physical partition was accepted')

    call make_fixture(candidate, assembly)
    assembly%classes(3)%canonical_total_variation = 2.0_dp
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_MEASURE_ERROR, &
        'incorrect total measure was accepted')

    call make_fixture(candidate, assembly)
    assembly%classes(2)%left_boundary_id = 12
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_PARTITION_ERROR, &
        'duplicate boundary ownership was accepted')

    call make_fixture(candidate, assembly)
    candidate%component_id = huge(candidate%component_id)
    call bridge_gc_operational_class_assembly(candidate, assembly, split, status)
    call require(status == GC_CLASS_ADAPTER_BRIDGE_ID_OVERFLOW, &
        'component-ID overflow was not rejected')

    call pass_test

contains

    subroutine make_fixture(candidate, assembly)
        type(gc_cylindrical_class_interval_t), intent(out) :: candidate
        type(gc_operational_class_assembly_result_t), intent(out) :: assembly

        candidate = gc_cylindrical_class_interval_t()
        candidate%component_id = 7
        candidate%sigma = 1
        candidate%rc_min = 0.0_dp
        candidate%rc_max = 10.0_dp
        candidate%psi_star_min = 2.0_dp
        candidate%psi_star_max = 8.0_dp
        candidate%canonical_measure = 12.0_dp
        candidate%rc_min_enclosure = gc_interval_t(0.0_dp, 0.0_dp)
        candidate%rc_max_enclosure = gc_interval_t(10.0_dp, 10.0_dp)
        candidate%psi_star_min_enclosure = gc_interval_t(2.0_dp, 2.0_dp)
        candidate%psi_star_max_enclosure = gc_interval_t(8.0_dp, 8.0_dp)
        candidate%canonical_measure_enclosure = gc_interval_t(12.0_dp, 12.0_dp)
        candidate%lower_root = .true.
        candidate%upper_root = .true.
        candidate%lower_tangent = .true.
        candidate%lower_boundary_kind = 'candidate-lower-root'
        candidate%upper_boundary_kind = 'candidate-upper-root'
        candidate%limiting_chart = 'manufactured-certified-cut'
        candidate%topology_certified = .true.
        candidate%root_isolation_certified = .true.
        candidate%allowed_interval = .true.
        candidate%lower_root_certificate = root_certificate(101, 0.0_dp)
        candidate%upper_root_certificate = root_certificate(202, 10.0_dp)
        candidate%lower_root_coordinate_map = root_map(301, 0.0_dp)
        candidate%upper_root_coordinate_map = root_map(302, 10.0_dp)

        assembly = gc_operational_class_assembly_result_t()
        assembly%status = GC_CLASS_ASSEMBLY_SUCCESS
        assembly%complete = .true.
        assembly%nclasses = 3
        allocate(assembly%classes(3))
        call set_class(assembly%classes(1), 1, 0.0_dp, 3.0_dp, 2.0_dp, &
            5.0_dp, GC_CLASS_BOUNDARY_REGULAR, GC_CLASS_BOUNDARY_X, 14, &
            0, 11, 4.0_dp)
        call set_class(assembly%classes(2), 2, 3.0_dp, 7.0_dp, 5.0_dp, &
            1.0_dp, GC_CLASS_BOUNDARY_X, GC_CLASS_BOUNDARY_SEPARATRIX, 43, &
            11, 12, 7.0_dp)
        call set_class(assembly%classes(3), 3, 7.0_dp, 10.0_dp, 1.0_dp, &
            8.0_dp, GC_CLASS_BOUNDARY_SEPARATRIX, GC_CLASS_BOUNDARY_REGULAR, &
            31, 12, 0, 1.0_dp)
    end subroutine make_fixture

    subroutine set_class(class, class_id, x_lo, x_hi, canonical_lo, &
            canonical_hi, left_kind, right_kind, ifuntype, left_id, right_id, &
            variation)
        type(gc_operational_class_interval_t), intent(out) :: class
        integer, intent(in) :: class_id, left_kind, right_kind, ifuntype
        integer, intent(in) :: left_id, right_id
        real(dp), intent(in) :: x_lo, x_hi, canonical_lo, canonical_hi, variation

        class = gc_operational_class_interval_t()
        class%class_id = class_id
        class%component_id = 7
        class%sigma = 1
        class%x_lo = x_lo
        class%x_hi = x_hi
        class%canonical_lo = canonical_lo
        class%canonical_hi = canonical_hi
        class%left_kind = left_kind
        class%right_kind = right_kind
        class%ifuntype = ifuntype
        class%left_boundary_id = left_id
        class%right_boundary_id = right_id
        class%canonical_total_variation = variation
    end subroutine set_class

    function root_certificate(cut_id, coordinate) result(root)
        integer, intent(in) :: cut_id
        real(dp), intent(in) :: coordinate
        type(gc_interval_root_box_t) :: root

        root = gc_interval_root_box_t()
        root%lo = coordinate
        root%hi = coordinate
        root%cut_id = cut_id
    end function root_certificate

    function root_map(map_certificate_id, coordinate) result(mapping)
        integer, intent(in) :: map_certificate_id
        real(dp), intent(in) :: coordinate
        type(gc_cylindrical_root_coordinate_map_t) :: mapping

        mapping = gc_cylindrical_root_coordinate_map_t()
        mapping%map_certificate_id = map_certificate_id
        mapping%mapped_class_enclosure = gc_interval_t(coordinate, coordinate)
    end function root_map

end program test_gc_operational_class_adapter_bridge
