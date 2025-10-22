module neort_datatypes

    implicit none

    type :: magfie_params_t
        real(8) :: s
        real(8) :: R0
        real(8) :: a
        real(8) :: eps
        real(8) :: psi_pr
        real(8) :: B0
        real(8) :: Bthcov
        real(8) :: Bphcov
        real(8) :: dBthcovds
        real(8) :: dBphcovds
        real(8) :: q
        real(8) :: iota
        real(8) :: dVds
        real(8) :: M_t
        real(8) :: Om_tE
        real(8) :: Om_tBref
        real(8) :: vth
        real(8) :: T_in_eV
        real(8) :: m0
        real(8) :: n0
        real(8) :: Dp
        real(8) :: Drp_over_Dp
        real(8) :: etatp
        real(8) :: etadt
        logical :: pertfile
        logical :: nonlin
        real(8) :: dpp
        real(8) :: dhh
        real(8) :: fpeff
    end type magfie_params_t

    type :: magfie_tensors_t
        real(8) :: theta
        real(8) :: bmod
        real(8) :: sqrtg
        real(8) :: hder(3)
        real(8) :: hcovar(3)
        real(8) :: hctrvr(3)
        real(8) :: hcurl(3)
        complex(8) :: bn
        complex(8) :: eps_exp  ! TODO: more meaningful name
    end type magfie_tensors_t

    type :: magfie_data_t
        type(magfie_params_t) :: params
        integer :: n_points = 0
        type(magfie_tensors_t), allocatable :: tensors(:)
    end type magfie_data_t

end module neort_datatypes
