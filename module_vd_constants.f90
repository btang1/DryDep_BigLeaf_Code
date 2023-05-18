module vd_constants
    implicit none

    ! DEFINE CONSTANTS USED IN CALCULATIONS 
    real, parameter  :: vk = 0.4         ! Von Karman constant                       (1)
    real, parameter :: rmin = 1.0        ! minimum resistance                        (?)
    real, parameter :: rmax = 1.e5       ! maximum resistance                        (?)
    real, parameter :: xmfp = 6.5e-8     ! mean free path                            (?)
    real, parameter :: g = 9.8           ! gravitational acceleration constant       (1)
    real, parameter :: vabs = 1.81e-2    ! dynamic viscosity of air                  (g/m-s)
    real, parameter :: boltz = 1.38e-20  ! Boltzmann constant.                       (m^2*g/(s^2*K))
    real, parameter :: pi = 3.1415927    ! Circle constant.                          (1)
    real, parameter :: vair = 1.5e-5     ! kinematic viscosity of air                (m^2/s)   
    real, parameter :: diffh2o = 2.3e-5  ! molecular diffusivity of water in air     (m^2/s)

end module vd_constants
