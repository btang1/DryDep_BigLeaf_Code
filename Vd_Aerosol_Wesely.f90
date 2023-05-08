subroutine Vd_Aerosol_Wesely(z0,deltaz,psih,ustar,diam,rhop,ts,vd)

!developed by Dr.Beiming Tang based on CAMx code
!NOAA ARL UFS project
!05/06/2023

!Modifications:

!Input arguments:
!    z0        surface roughness length                (unit = m)
!    deltaz    Layer 1 midpoint height                 (unit = m)
!    psih      similarity stability correction term    (unit = 1)
!    ustar     friction velocity                       (unit = m/s)
!    diam      log-mean sectional aerosol diameter     (unit = m)
!    rhop      aerosol density                         (unit = g/m^3)
!    ts        surface temperature                     (unit = K)

!Output arguments:
!    vd        deposition velocity                     (unit = m/s)

!Routines called:
!    none

!Called by:
!    should be part of dry deposition code

!Step 0. Definal variables type & constant values
    real :: z0
    real :: deltaz
    real :: psih
    real :: ustar
    real :: diam
    real :: rhop
    real :: ts
    real :: vd

    real :: vk = 0.4
    real :: rmin = 1.0
    real :: xmfp = 6.5e-8
    real :: g = 9.8
    real :: vabs = 1.81e-2               !dynamic viscosity of air
    real :: boltz = 1.38e-20
    real :: pi = 3.1415927
    real :: vair = 1.5e-5                !kinetic viscosity of air

!
!Main Program Start
!

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
    end
























