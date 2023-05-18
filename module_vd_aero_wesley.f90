module vd_aero_wesely
   
use vd_constants, only :: xmfp, g, vabs, boltz, vair, vk, rmin, pi

implicit none

contains

   subroutine aero_wesely(z0,deltaz,psih,ustar,diam,rhop,ts,vd)

      !developed by Dr.Beiming Tang based on CAMx code
      !NOAA ARL UFS project
      !05/06/2023

      !Modifications:

      !Input arguments:
      !    z0        surface roughness length                (m)
      !    deltaz    Layer 1 midpoint height                 (m)
      !    psih      similarity stability correction term    (-)
      !    ustar     friction velocity                       (m/s)
      !    diam      log-mean sectional aerosol diameter     (m)
      !    rhop      aerosol density                         (g/m^3)
      !    ts        surface temperature                     (K)

      !Output arguments:
      !    vd        deposition velocity                     (m/s)

      !Routines called:
      !    vd_constants

      !Called by:
      !    future routine in UFS dry deposition

      !Step 0. Definal variables type
      !!Step 0-1.Input and output variables
      real, intent(in) :: z0       ! surface roughness length                (m)
      real, intent(in) :: deltaz   ! Layer 1 midpoint height                 (m)
      real, intent(in) :: psih     ! similarity stability correction term    (-)
      real, intent(in) :: ustar    ! friction velocity                       (m/s)
      real, intent(in) :: diam     ! log-mean sectional aerosol diameter     (m)
      real, intent(in) :: rhop     ! aerosol density                         (g/m^3)
      real, intent(in) :: ts       ! surface temperature                     (K)
      real, intent(out) :: vd      ! deposition velocity                     (m/s)
      
      !Step 0-2. Local varibles
      real :: power1               ! 
      real :: scf                  ! Cunningham correction for small particles
      real :: vsed.                ! gravitational settling velocity. (m/s)
      real :: difbrwn              ! Brownian diffusivity
      real :: schmidt              ! Schmidt number
      real :: stokes               ! Stokes number
      real :: ra                   ! Atmospheric resistance
      real :: sc23                 ! Schmidt number to power of (-2/3)
      real :: power2               ! -3/schmidt number      
      real :: xinert               ! 10 ^(-3/sstokes)
      real :: rd                   ! boundary layer resistance
      
      !Step 1. Calculate speed correction factor and sedimendation velocity
      power1 = amin1(7.6,0.55*diam/xmfp)
      scf   = 1. + (2.514 + 0.8 * exp(-power1))*xmfp/diam
      vsed  = rhop*g*(diam*diam)*scf/(18.*vabs)

      !Step 2. Calcualte Brownian diffusivity, Schmidt number, and Stokes number
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
      stokes  = vsed*(ustar*ustar)/(vair*g)

      !Step 3. Calculate atmospheric resistance, Ra
      ra = (log(deltaz/z0)-psih)/(vk*ustar)
      ra = amax1(ra,rmin)

      !Step 4. Calcualte deposition boundary layer resistance, Rd
      sc23 = schmidt **(-2./3.)
      power2 = -3./stokes
      if (power2 < -37) then
         xinert = 10. **(-37.)
      else
         xinert = 10. **(power2)
      endif

      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd, rmin)

      !Step 5. Calculate aerosol depostion velocity for this cell, land use, and aerosol size
      vd = vsed + 1./(ra + rd + ra* rd * vsed)

      return

   end subroutine aerosol_wesely

end module  vd_aero_wesely
























