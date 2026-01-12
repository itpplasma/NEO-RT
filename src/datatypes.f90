module neort_datatypes
    use iso_fortran_env, only: dp => real64

    implicit none

    type :: magfie_params_t
        real(dp) :: s  ! Normalised toroidal flux
        real(dp) :: R0  ! Major radius [cgs]
        real(dp) :: a  ! Minor radius [cgs]
        real(dp) :: eps  ! Inverse aspect ratio (a/R0)
        real(dp) :: psi_pr  ! Radial derivative of poloidal flux
        real(dp) :: B0  ! Reference magnetic field magnitude [cgs]
        real(dp) :: Bthcov  ! Covariant poloidal field component
        real(dp) :: Bphcov  ! Covariant toroidal field component
        real(dp) :: dBthcovds  ! Derivative of Bthcov with respect to s
        real(dp) :: dBphcovds  ! Derivative of Bphcov with respect to s
        real(dp) :: q  ! Safety factor
        real(dp) :: iota  ! Rotational transform (1/q)
        real(dp) :: dVds  ! Derivative of flux-surface volume with respect to s
        real(dp) :: M_t  ! Toroidal Mach number
        real(dp) :: Om_tE  ! Toroidal rotation frequency from electric field
        real(dp) :: Om_tBref  ! Reference toroidal rotation frequency
        real(dp) :: vth  ! Thermal velocity [cgs]
        real(dp) :: T_in_eV  ! Temperature [eV]
        real(dp) :: m0  ! Poloidal mode number of perturbation
        real(dp) :: n0  ! Toroidal mode number of perturbation
        real(dp) :: Dp  ! Thermodynamic force parameter
        real(dp) :: Drp_over_Dp  ! Ratio of rotation force to thermodynamic force
        real(dp) :: etatp  ! Normalized temperature gradient
        real(dp) :: etadt  ! Normalized density gradient
        logical :: pertfile  ! Flag: use perturbation from file (not analytic)
        logical :: nonlin  ! Flag: enable nonlinear physics calculations
        real(dp) :: dpp  ! Nonlinear diagnostic parameter
        real(dp) :: dhh  ! Nonlinear diagnostic parameter
        real(dp) :: fpeff  ! Effective passing fraction parameter
    end type magfie_params_t

    type :: magfie_tensors_t
        real(dp) :: theta  ! Boozer poloidal angle [rad]
        real(dp) :: bmod  ! Magnetic field magnitude |B| [cgs]
        real(dp) :: sqrtg  ! Jacobian determinant sqrt(g)
        real(dp) :: hder(3)  ! Contravariant basis vectors
        real(dp) :: hcovar(3)  ! Covariant basis vectors
        real(dp) :: hctrvr(3)  ! Contravariant components of field-line curvature
        real(dp) :: hcurl(3)  ! Components of curl of basis vectors
        complex(dp) :: bn  ! Perturbation amplitude normalised by |B|
        complex(dp) :: eps_exp  ! Reference analytic perturbation epsmn*exp(i*m0*theta)
    end type magfie_tensors_t

    type :: magfie_data_t
        type(magfie_params_t) :: params  ! Magnetic field parameters and configuration
        integer :: n_points = 0  ! Number of poloidal angle sampling points
        type(magfie_tensors_t), allocatable :: tensors(:)  ! Magnetic field data
                                                           ! sampled along poloidal angle
    end type magfie_data_t

    type :: transport_summary_t
        real(dp) :: M_t  ! Toroidal Mach number
        real(dp) :: Dco(2)  ! Transport coefficients (D11, D12) for co-passing particles
        real(dp) :: Dctr(2)  ! Transport coefficients (D11, D12) for counter-passing particles
        real(dp) :: Dt(2)  ! Transport coefficients (D11, D12) for trapped particles
        ! D11: particle flux coefficient, D12: momentum flux coefficient
    end type transport_summary_t

    type :: torque_summary_t
        logical :: has_torque  ! Flag: torque calculation enabled
        real(dp) :: s  ! Normalised toroidal flux
        real(dp) :: dVds  ! Derivative of flux-surface volume with respect to s
        real(dp) :: M_t  ! Toroidal Mach number
        real(dp) :: Tco  ! Torque density from co-passing particles
        real(dp) :: Tctr  ! Torque density from counter-passing particles
        real(dp) :: Tt  ! Torque density from trapped particles
    end type torque_summary_t

    type :: transport_harmonic_t
        integer :: mth  ! Poloidal mode number associated with resonance
        real(dp) :: Dresco(2)
        ! Transport coefficients (D11, D12) for co-passing particles, per harmonic
        real(dp) :: Dresctr(2)
        ! Transport coefficients (D11, D12) for counter-passing particles, per harmonic
        real(dp) :: Drest(2)  ! Transport coefficients (D11, D12) for trapped particles, per harmonic
        real(dp) :: Tresco  ! Torque from co-passing particles, per harmonic
        real(dp) :: Tresctr  ! Torque from counter-passing particles, per harmonic
        real(dp) :: Trest  ! Torque from trapped particles, per harmonic
        real(dp) :: vminp_over_vth  ! Lower bound of passing velocity grid normalised
        real(dp) :: vmaxp_over_vth  ! Upper bound of passing velocity grid normalised
        real(dp) :: vmint_over_vth  ! Lower bound of trapped velocity grid normalised
        real(dp) :: vmaxt_over_vth  ! Upper bound of trapped velocity grid normalised
    end type transport_harmonic_t

    type :: transport_data_t
        type(transport_summary_t) :: summary  ! Summed transport coefficients across all harmonics
        type(torque_summary_t) :: torque  ! Summed torque density across all harmonics
        type(transport_harmonic_t), allocatable :: harmonics(:)  ! Per-harmonic transport and torque
    end type transport_data_t

end module neort_datatypes
