module vd_aero_wesley
   
use vd_constants

implicit none

contains

   subroutine aero_wesley(z0,deltaz,psih,ustar,diam,rhop,ts,vd)

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
      !    none

      !Called by:
      !    future routine

      !Step 0. Definal variables type & constant values
      real, intent(in) :: z0       ! Add comments and units to these
      real, intent(in) :: deltaz   ! Add comments and units to these
      real, intent(in) :: psih     ! Add comments and units to these
      real, intent(in) :: ustar    ! Add comments and units to these
      real, intent(in) :: diam     ! Add comments and units to these
      real, intent(in) :: rhop     ! Add comments and units to these
      real, intent(in) :: ts       ! Add comments and units to these
      real, intent(out) :: vd       ! Add comments and units to these

      !Step 1. Calculate speed correction factor and sedimendation velocity
      power = amin1(7.6,0.55*diam/xmfp)
      scf   = 1. + (2.514 + 0.8 * exp(-power))*xmfp/diam
      vsed  = rhop*g*(diam*diam)*scf/(18.*vabs)

      !Step 2. Calcualte Brownian diffusivity and Schmidt number and Stokes number
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
      stokes  = vsed*(ustar*ustar)/(vair*g)

      !Step 3. Calculate atmospheric resistance, Ra
      ra = (log(deltaz/z0)-psih)/(vk*ustar)
      ra = amax1(ra,rmin)

      !Step 4. Calcualte deposition layer resistance, Rd
      sc23 = schmidt **(-2./3.)
      power = -3./stokes
      if (power <= -37) then
         xinert = 10. **(-37.)
      else
         xinert = 10. **(power)
      endif

      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd, rmin)

      !Step 5. Calculate aerosol depostion velocity for this cell, land use, and aerosol size
      vd = vsed + 1./(ra + rd + ra* rd * vsed)

      return

   end subroutine Vd_Aerosol_Wesely

end module  vd_aero_wesley
























