module vd_constants
    implicit none

    ! DEFINE CONSTANTS USED IN CALCULATIONS 
    real, parameter  :: vk = 0.4         ! Add comments and units to these
    real, parameter :: rmin = 1.0        ! Add comments and units to these
    real, parameter :: rmax = 1.e5       ! Add comments and units to these
    real, parameter :: xmfp = 6.5e-8     ! Add comments and units to these
    real, parameter :: g = 9.8           ! Add comments and units to these
    real, parameter :: vabs = 1.81e-2    ! dynamic viscosity of air
    real, parameter :: boltz = 1.38e-20  ! Add comments and units to these
    real, parameter :: pi = 3.1415927    ! Add comments and units to these
    real, parameter :: vair = 1.5e-5     ! kinetic viscosity of air
    real, parameter :: diffh2o = 2.3e-5  ! Add comments and units to these

end module vd_constants
