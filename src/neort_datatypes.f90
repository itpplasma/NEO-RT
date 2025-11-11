module neort_datatypes

    implicit none

    type :: magfie_params_t
        real(8) :: s  ! Normalised toroidal flux
        real(8) :: R0  ! Major radius [cgs]
        real(8) :: a  ! Minor radius [cgs]
        real(8) :: eps  ! Inverse aspect ratio (a/R0)
        real(8) :: psi_pr  ! Radial derivative of poloidal flux
        real(8) :: B0  ! Reference magnetic field magnitude [cgs]
        real(8) :: Bthcov  ! Covariant poloidal field component
        real(8) :: Bphcov  ! Covariant toroidal field component
        real(8) :: dBthcovds  ! Derivative of Bthcov with respect to s
        real(8) :: dBphcovds  ! Derivative of Bphcov with respect to s
        real(8) :: q  ! Safety factor
        real(8) :: iota  ! Rotational transform (1/q)
        real(8) :: dVds  ! Derivative of flux-surface volume with respect to s
        real(8) :: M_t  ! Toroidal Mach number
        real(8) :: Om_tE  ! Toroidal rotation frequency from electric field
        real(8) :: Om_tBref  ! Reference toroidal rotation frequency
        real(8) :: vth  ! Thermal velocity [cgs]
        real(8) :: T_in_eV  ! Temperature [eV]
        real(8) :: m0  ! Poloidal mode number of perturbation
        real(8) :: n0  ! Toroidal mode number of perturbation
        real(8) :: Dp  ! Thermodynamic force parameter
        real(8) :: Drp_over_Dp  ! Ratio of rotation force to thermodynamic force
        real(8) :: etatp  ! Normalized temperature gradient
        real(8) :: etadt  ! Normalized density gradient
        logical :: pertfile  ! Flag: use perturbation from file (not analytic)
        logical :: nonlin  ! Flag: enable nonlinear physics calculations
        real(8) :: dpp  ! Nonlinear diagnostic parameter
        real(8) :: dhh  ! Nonlinear diagnostic parameter
        real(8) :: fpeff  ! Effective passing fraction parameter
    end type magfie_params_t

    type :: magfie_tensors_t
        real(8) :: theta  ! Boozer poloidal angle [rad]
        real(8) :: bmod  ! Magnetic field magnitude |B| [cgs]
        real(8) :: sqrtg  ! Jacobian determinant sqrt(g)
        real(8) :: hder(3)  ! Contravariant basis vectors
        real(8) :: hcovar(3)  ! Covariant basis vectors
        real(8) :: hctrvr(3)  ! Contravariant components of field-line curvature
        real(8) :: hcurl(3)  ! Components of curl of basis vectors
        complex(8) :: bn  ! Perturbation amplitude normalised by |B|
        complex(8) :: eps_exp  ! Reference analytic perturbation epsmn*exp(i*m0*theta)
    end type magfie_tensors_t

    type :: magfie_data_t
        type(magfie_params_t) :: params  ! Magnetic field parameters and configuration
        integer :: n_points = 0  ! Number of poloidal angle sampling points
        type(magfie_tensors_t), allocatable :: tensors(:)  ! Magnetic field data
                                                           ! sampled along poloidal angle
    end type magfie_data_t

    type :: transport_summary_t
        real(8) :: M_t  ! Toroidal Mach number
        real(8) :: Dco(2)  ! Transport coefficients (D11, D12) for co-passing particles
        real(8) :: Dctr(2)  ! Transport coefficients (D11, D12) for counter-passing particles
        real(8) :: Dt(2)  ! Transport coefficients (D11, D12) for trapped particles
        ! D11: particle flux coefficient, D12: momentum flux coefficient
    end type transport_summary_t

    type :: torque_summary_t
        logical :: has_torque  ! Flag: torque calculation enabled
        real(8) :: s  ! Normalised toroidal flux
        real(8) :: dVds  ! Derivative of flux-surface volume with respect to s
        real(8) :: M_t  ! Toroidal Mach number
        real(8) :: Tco  ! Torque density from co-passing particles
        real(8) :: Tctr  ! Torque density from counter-passing particles
        real(8) :: Tt  ! Torque density from trapped particles
    end type torque_summary_t

    type :: transport_harmonic_t
        integer :: mth  ! Poloidal mode number associated with resonance
        real(8) :: Dresco(2)
        ! Transport coefficients (D11, D12) for co-passing particles, per harmonic
        real(8) :: Dresctr(2)
        ! Transport coefficients (D11, D12) for counter-passing particles, per harmonic
        real(8) :: Drest(2)  ! Transport coefficients (D11, D12) for trapped particles, per harmonic
        real(8) :: Tresco  ! Torque from co-passing particles, per harmonic
        real(8) :: Tresctr  ! Torque from counter-passing particles, per harmonic
        real(8) :: Trest  ! Torque from trapped particles, per harmonic
        real(8) :: vminp_over_vth  ! Lower bound of passing velocity grid normalised
        real(8) :: vmaxp_over_vth  ! Upper bound of passing velocity grid normalised
        real(8) :: vmint_over_vth  ! Lower bound of trapped velocity grid normalised
        real(8) :: vmaxt_over_vth  ! Upper bound of trapped velocity grid normalised
    end type transport_harmonic_t

    type :: transport_data_t
        type(transport_summary_t) :: summary  ! Summed transport coefficients across all harmonics
        type(torque_summary_t) :: torque  ! Summed torque density across all harmonics
        type(transport_harmonic_t), allocatable :: harmonics(:)  ! Per-harmonic transport and torque
    end type transport_data_t

end module neort_datatypes
