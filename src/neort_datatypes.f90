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

    type :: transport_summary_t
        real(8) :: M_t  ! toroidal Mach number of electric precession
        real(8) :: Dco(2)  ! D11, D12 transport coefficients for co-passing particles, non-axisymmetric
        real(8) :: Dctr(2)  ! D11, D12 transport coefficients for counter-passing particles, non-axisymmetric
        real(8) :: Dt(2)  ! D11, D12 transport coefficients for trapped particles, non-axisymmetric
    end type transport_summary_t

    type :: torque_summary_t
        logical :: has_torque
        real(8) :: s
        real(8) :: dVds  ! derivative of volume inside flux surface by s
        real(8) :: M_t  ! toroidal Mach number of electric precession
        real(8) :: Tco  ! torque by co-passing particles
        real(8) :: Tctr  ! torque by counter-passing particles
        real(8) :: Tt  ! torque by trapped particles
    end type torque_summary_t

    type :: transport_harmonic_t
        integer :: mth  ! canonical poloidal harmonic number
        real(8) :: Dresco(2)  ! D11, D12 transport coefficients for co-passing particles, resonant, per harmonic
        real(8) :: Dresctr(2)  ! D11, D12 transport coefficients for counter-passing particles, resonant, per harmonic
        real(8) :: Drest(2)  ! D11, D12 transport coefficients for trapped particles, resonant, per harmonic
        real(8) :: Tresco  ! torque by co-passing particles, resonant, per harmonic
        real(8) :: Tresctr  ! torque by counter-passing particles, resonant, per harmonic
        real(8) :: Trest  ! torque by trapped particles, resonant, per harmonic
        real(8) :: vminp_over_vth  ! lower integration limit of normalized velocity for passing particles
        real(8) :: vmaxp_over_vth  ! upper integration limit of normalized velocity for passing particles
        real(8) :: vmint_over_vth  ! lower integration limit of normalized velocity for trapped particles
        real(8) :: vmaxt_over_vth  ! upper integration limit of normalized velocity for trapped particles
    end type transport_harmonic_t

    type :: transport_data_t
        type(transport_summary_t) :: summary
        type(torque_summary_t) :: torque
        type(transport_harmonic_t), allocatable :: harmonics(:)
    end type transport_data_t

end module neort_datatypes
