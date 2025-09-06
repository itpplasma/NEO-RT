module neo2_reader
    use hdf5
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    type :: neo2_data_t
        real(dp), allocatable :: stor(:)
        real(dp), allocatable :: spol(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: n_spec(:,:)
        real(dp), allocatable :: T_spec(:,:)
        real(dp), allocatable :: Er(:)
        real(dp), allocatable :: MtOvR(:,:)
    end type neo2_data_t

    type profiles_t
        real(dp), allocatable :: n(:)
        real(dp), allocatable :: T(:)
        real(dp), allocatable :: phi_e(:)
    end type profiles_t

    interface h5read
        module procedure h5_read_1d, h5_read_2d
    end interface h5read
    
    contains

    subroutine read_neo2_data(neo2_file, neo2_data)
        character(len=*), intent(in) :: neo2_file
        type(neo2_data_t), intent(out) :: neo2_data

        integer(hid_t) :: file_id
        integer :: hdferr
        
        call h5open_f(hdferr)
        call h5fopen_f(neo2_file, H5F_ACC_RDONLY_F, file_id, hdferr)
        
        call h5read(file_id, "boozer_s", neo2_data%stor)
        call h5read(file_id, "aiota", neo2_data%q)
        neo2_data%q = 1.0_dp / neo2_data%q

        allocate(neo2_data%spol(size(neo2_data%stor)))
        call stor_to_spol(neo2_data%stor, neo2_data%q, neo2_data%spol)

        call h5read(file_id, "n_spec", neo2_data%n_spec)
        call h5read(file_id, "T_spec", neo2_data%T_spec)
        call h5read(file_id, "Er", neo2_data%Er)
        call h5read(file_id, "MtOvR", neo2_data%MtOvR)
        
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)
        
    end subroutine read_neo2_data
    
    subroutine stor_to_spol(stor, q, spol)
        ! Convert toroidal flux (stor) to poloidal flux (spol)
        ! Uses trapezoidal integration: dpsi_pol/dpsi_tor = 1/q
        ! So psi_pol = integral(1/q * dpsi_tor)
        
        real(dp), intent(in) :: stor(:), q(:)
        real(dp), intent(out) :: spol(:)
        
        integer :: i, n
        real(dp) :: dpsi_tor, dpsi_pol
        
        n = size(stor)
        
        ! Start at axis
        spol(1) = 0.0_dp
        
        ! Trapezoidal integration
        do i = 2, n
            dpsi_tor = stor(i) - stor(i-1)
            dpsi_pol = 0.5_dp * (1.0_dp/q(i) + 1.0_dp/q(i-1)) * dpsi_tor
            spol(i) = spol(i-1) + dpsi_pol
        enddo
        
        ! Normalize to [0,1]
        spol = (spol - spol(1)) / (spol(n) - spol(1))
        
    end subroutine stor_to_spol

    subroutine polyfit(x, y, order, coef)
        ! Fits polynomial of given order to (x,y) data using least squares
        ! Returns coefficients in POTATO convention: coef(0) is highest power coefficient
        ! Polynomial: p(x) = coef(0)*x^order + coef(1)*x^(order-1) + ... + coef(order)
        
        real(dp), intent(in) :: x(:), y(:)
        integer, intent(in) :: order
        real(dp), intent(out) :: coef(0:)
        
        integer :: n, nrhs, lda, ldb, info, i, j
        real(dp), allocatable :: A(:,:), b(:), work(:)
        
        n = order + 1
        nrhs = 1
        
        allocate(A(size(x), n))
        allocate(b(size(x)))
        allocate(work(max(1, min(size(x), n) + max(size(x), n))))
        
        ! Build Vandermonde matrix
        do i = 1, size(x)
            do j = 1, n
                A(i, j) = x(i)**(order - j + 1)
            enddo
        enddo
        
        b = y
        
        ! Solve least squares problem
        lda = size(x)
        ldb = size(x)
        call dgels('N', size(x), n, nrhs, A, lda, b, ldb, work, size(work), info)
        
        ! Extract coefficients in POTATO convention
        do i = 0, order
            coef(i) = b(i + 1)
        enddo
        
        deallocate(A, b, work)
        
    end subroutine polyfit

    subroutine polyval(x, coef, y)
        ! Evaluate polynomial with given coefficients at points x
        ! Coefficients in POTATO convention: coef(0) is highest power coefficient
        
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: coef(0:)
        real(dp), intent(out) :: y(:)
        
        integer :: i, j, order
        order = size(coef) - 1
        
        y = 0.0_dp
        do i = 1, size(x)
            do j = 0, order
                y(i) = y(i) + coef(j) * x(i)**(order - j)
            enddo
        enddo
        
    end subroutine polyval
    
    subroutine polyint(coef_in, coef_out, const)
        ! Integrate polynomial coefficients
        ! Input: coef_in(0:n) represents polynomial of degree n
        ! Output: coef_out(0:n+1) represents integrated polynomial of degree n+1
        ! const: integration constant (value at x=0)
        
        real(dp), intent(in) :: coef_in(0:)
        real(dp), intent(out) :: coef_out(0:)
        real(dp), intent(in), optional :: const
        
        integer :: i, n
        n = size(coef_in) - 1
        
        ! Integrate each term: x^k -> x^(k+1)/(k+1)
        do i = 0, n
            coef_out(i) = coef_in(i) / real(n - i + 1, dp)
        enddo
        
        ! Add integration constant as the lowest order term
        if (present(const)) then
            coef_out(n+1) = const
        else
            coef_out(n+1) = 0.0_dp
        endif
        
    end subroutine polyint
    
    subroutine h5_read_1d(file_id, name, array)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(out) :: array(:)
        
        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(1)
        integer :: hdferr
        
        call h5dopen_f(file_id, name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, hdferr)
        allocate(array(dims(1)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine h5_read_1d
    
    subroutine h5_read_2d(file_id, name, array)
        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: name
        real(dp), allocatable, intent(out) :: array(:,:)
        
        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(2)
        integer :: hdferr
        
        call h5dopen_f(file_id, name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_dims_f(space_id, dims, dims, hdferr)
        allocate(array(dims(1), dims(2)))
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine h5_read_2d

end module neo2_reader