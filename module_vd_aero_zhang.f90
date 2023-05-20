module vd_aerosol_zhang
    !Developed by Dr.Beiming Tang based on CAMx code
    !NOAA ARL UFS project
    !05/09/2023
    
    !LUC No.        Vegetation type
    !    1          water
    !    2          ice
    !    3          inland lake
    !    4          evengreen needleleaf trees
    !    5          evengreen broadleaf  trees
    !    6          deciduous needlelead trees
    !    7          deciduous broadleaf  trees
    !    8          tropical broadleaf   trees
    !    9          drought deciduous    trees
    !   10          evergreen broadlead  shrubs
    !   11          deciduous            shrubs
    !   12          thorn                shrubs
    !   13          short grass and forbs
    !   14          long grass
    !   15          crops
    !   16          rice
    !   17          sugar
    !   18          maize
    !   19          cotton
    !   20          irrigated crops
    !   21          urban
    !   22          tundra
    !   23          swamp
    !   24          desert
    !   25          mixed wood forests
    !   26          transitional forests
    
    use vd_constants, only :: ???
    
    implicit none

    contains
    
    subroutine Aerosol_Zhang(z2,zl,z0_f,utar,t2,ts,lai_f,mlu,vdp,diam,rhoprt)
    !Input arguments
    !    z2         meteorology reference height          (m)
    !    zl         Z/L stability parameter               (1)
    !    z0_f       surface roughness                     (m)
    !    ustar      friction velocity                     (m/s)
    !    t2         temperature at Z2                     (K)
    !    ts         surface temperature                   (K)
    !    lai_f      leaf area index                       (1)
    !    mlu        land use type                         (1)
    !    vdp        particle dry deposition velocity      (m/s)
    !    diam       log-mean sectional aerosol diameter   (m)
    !    rhoprt     aerosol density                       (g/m^3) 

    !    ra         aerosol resistance                    (s/m)
    !    aest       parameter for calculating EIM         (?)
    !    binsize    radius of a size bin                  (?)
    !    eb         Brownian collection efficiency        (?)
    !    eim        impaction collection efficiency       (?)
    !    ein        interception collection efficiency    (?)
    !    gama       parameter for calculating EB          (?)
    !    pllp       leaf dimension for calculating EIN    (?)

    !Output argument:
    !    vdp        deposition velocity                   (unit = m/s)

    !Step 0. Define input variables type & value constant
    !!Step 0-1. Define global variables
        real    :: z2
        real    :: zl
        real    :: z0_f
        real    :: ustar
        real    :: t2
        real    :: ts
        real    :: lai_f
        integer :: mlu
        real    :: vdp
        real    :: diam
        real    :: rhopt

    !!Step 0-2. Define local varible
        integer :: i                             !land category
        real    :: rhop
        real    :: ra
        real    :: pllp
        real    :: binsize
        real    :: tave
        real    :: amu
        real    :: anu
        real    :: prii
        real    :: priiv
        real    :: vphil
        real    :: cfac
        real    :: taurel                                     !Relaxation time
        real    :: amob
        real    :: diff                                       !Brownian diffusivity
        real    :: schm
        real    :: pdepv
        real    :: st
        real    :: eb
        real    :: eim
        real    :: ein
        real    :: r1
        real    :: rs
        real    :: vdsize
    
    !!Step 0-3. Define constant values
        real,parameter    :: aa1    = 1.257
        real,parameter    :: aa2    = 0.4
        real,parameter    :: aa3    = 1.1
        real,parameter    :: amfp   = 6.53e-8                  !Lambda, mean free path for air molecules
        real,parameter    :: roarow = 1.19                     ! air density at 20C, (unit = kg/m^3)
        real,parameter    :: boltzk = 1.3806044503487214e-23
        real,parameter    :: dair = 0.369*29.+6.29
        real,parameter    :: dh2o = 0.369*18.+6.29

        real,dimension(26)    :: aest  
        real,dimension(26)    :: pllp1 !mim leaf dimension for each LUC
        real,dimension(26)    :: pllp2 !max leaf dimension for each LUC
        real,dimension(26)    :: gama

        aset  = (/100., 50.,100., 1.0, 0.8, 1.1, 0.8, 0.6, 1.0, 1.1,  &
                   1.1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1, 1.2, 1.2,  &
                   1.5, 50., 2.0, 50., 0.8, 0.8/)
        pllp1 = (/-0.9,-0.9,-0.9, 2.0, 5.0, 2.0, 5.0, 5.0, 5.0, 5.0,  &
                   5.0, 2.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 2.0,  &
                  10.0,-0.9,10.0,-0.9, 5.0, 5.0/)
        pllp2 = (/-0.9,-0.9,-0.9, 2.0, 5.0, 5.0,10.0, 5.0,10.0, 5.0,  &
                  10.0, 5.0, 5.0, 5.0, 5.0, 5.0,10.0,10.0,10.0, 5.0,  &
                  10.0,-0.9,10.0,-0.9, 5.0, 5.0/)
        gama  = (/ 0.5,0.54, 0.5,0.56,0.56,0.56,0.56,0.58,0.56,0.55,  &
                  0.55,0.54,0.54,0.55,0.54,0.54,0.54,0.55,0.54,0.54,  &
                  0.56,0.54,0.54,0.54,0.56,0.56/)
    
        !particle density          
        rhop  = rhoprt*1.0e-3           !Covert a unit from g/m^3 to kg/m^3

        vdp   = 0.0    
        i     = mlu     

    !
    !Main Program Start
    !


    !Step 1. Calculate aerodynamic resistance, Ra 
        if (zl >= 0) then
            ra = (0.74*log(z2/z0_f) + 4.7 *zl)/(0.4*ustar)
        else
            ra = 0.74*(log(z2/z0_f) - 2*log((1+sqrt(1-9.*zl))*0.5))/(0.4*ustar)
        endif
    
        ra = amax1(ra,5.0)
    
        if (i == 1 .or. i == 3) then   !for water or inland lake
            ra = amin1(ra,2000.)
        else
            ra = amin1(ra,1000.)
        endif

    !Step 2. Calculate settling velocity, Vg 
    !!Loop for partical species ???
        if (lai_ref(i,15) == lat_ref(i,14)) then    !lai_ref IS NOT DEFINED
            pllp = pllp2(i)
        else
            pllp = pllp2(i) - (lai_f-lai_ref(i,14))/(lai_ref(i,15)-lai_ref(i,14)+1.e-10)*(pllp2(i)-pllp1(i))
        endif

        !calculate radius of particle
        binsize  = diam/2

        tave     = 0.5*(t2 + ts)

        !calculate air dynamic viscosity
        amu      = 145.8 * 1.e-8 * tave ** 1.5/ (tave + 110.4)

        !calculate air kinetic viscosity
        anu      = amu/roarow
  
        !calculate cunningham slip correction factor (for samll particles) & relaxation time (=vg/gravity)
    
        prii     = 2.* 9.81 /(9. * amu)
        priiv    = prii * (rhop -roarow)
        vphil    = 0.
        cfac     = 1. + amfp/binsize*(aa1+aa2*exp(-aa3*binsize/amfp))            !Follow Zhang et al., (2001) eqn (3)
        taurel   = amax1(priiv*binsize**2*cfac/9.81,0)                  

        !calculate stokes friction and diffusion coefficients                    !Same as Wesely method for aerosol
        amob     = 6. * 3.14 * amu *binsize /cfac
        diff     = boltzk *tave /amob                                  
        schm     = anu/diff                                                      !Follow Zhang et al., (2001)

        !calculate gravitational settling velocity
        pdepv    = taurel * 9.81                                                 !Follow Zhang et al., (2001) eqn (2)


    !Step 3. Calculate surface resistance, Rs
        !calculate efficiency by diffusion, impaction, interception and particle rebound

        !!smooth surface
        if (pllp <= 0) then
            st   = taurel*g*ustar*ustar/anu                                      !Follow Slinn(1982), listed in Zhang et al.,(2001).CAMx original code is wrong, show *g??
        !!vegetated surface
        else
            st   = taurel*ustar/pllp*1000                                        !Follow Giorgi(1988),listed in Zhang et al.,(2001). where A = pllp/1000
        endif

        eb       = schm **(-gama(i))                                             !Follow Zhang et al.,(2001) eqn (6)
        eim      = (st/(st+aest(i)))**2                                          !Follow Zhang et al.,(2001) eqn (7c)
        ein      = 0.0
        if (pllp > 0.00) then
            ein  = (1000.*2.*binsize/pllp)**2*0.5                                !Follow Zhang et al.,(2001) eqn (8)
        end if

        r1       = exp(-st**0.5)                                                 !Follow Zhang et al.,(2001) eqn (9)
    
        if (r1 < 0.5) then
            r1   = 0.5
        rs       = 1./(3.*ustar*r1*(eb + eim + ein))                             !Follow Zhang et al.,(2001) eqn (5),where epsilon == 3

        !Step 4. Calculate deposition velocity
            vdsize   = pdepv + 1./(ra + rs)                                          !Follow Zhang et al.,(2001) eqn (1)
            vdp      = vdsize

        return
    
    end subroutine Aerosol_Zhang
    
 end module vd_aerosol_zhang






















































































