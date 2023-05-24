module vd_aerosol_zhang
    
    Use vd_constants, only ::xmfp,roarow,aa1,aa2,aa3,boltzk,dair,dh2o,aest,pllp1,pllp2,gama,lai_ref

    Implicit none
    
    contains

    subroutine aerosol_zhang(z2,zl,z0_f,utar,t2,ts,lai_f,mlu,diam,rhoprt,vd)

        !Developed by Dr. Beiming Tang based on CAMx code
        !NOAA ARL UFS project
        !05/20/2023
        
        !Modifications:
        
        !Input arguments
        !    z2         meteorology reference height          (m)
        !    zl         Z/L stability parameter               (1)
        !    z0_f       surface roughness                     (m)
        !    ustar      friction velocity                     (m/s)
        !    t2         temperature at Z2                     (K)
        !    ts         surface temperature                   (K)
        !    lai_f      leaf area index                       (1)
        !    mlu        land use type                         (1)
        !    diam       log-mean sectional aerosol diameter   (m)
        !    rhoprt     aerosol density                       (g/m^3) 

        !Output argument:
        !    vd         deposition velocity                   (m/s)
        
        !Routine called:
        !    vd_constants

        !Called by:
        !    future routine in UFS dry deposition
        
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
        !   26          transitional forest



        !Step 0. Define input variables type & value constant
        !!Step 0-1. Input and output variables
            real    :: z2
            real    :: zl
            real    :: z0_f
            real    :: ustar
            real    :: t2
            real    :: ts
            real    :: lai_f
            integer :: mlu
            real    :: diam
            real    :: rhoprt
            real    :: vd

        !!Step 0-2. Local variables
            integer :: i
            real    :: rhop
            real    :: ra
            real    :: pllp
            real    :: binsize                                      !radius of aerosol.  (m)
            real    :: tave
            real    :: amu                                          !air dynamic viscosity
            real    :: anu                                          !air kinematic viscosity
            real    :: prii
            real    :: priiv
            real    :: vphil
            real    :: cfac
            real    :: taurel                                       !Relaxation time
            real    :: amob
            real    :: diff                                         !Brownian diffusivity
            real    :: schm
            real    :: vg
            real    :: st
            real    :: eb
            real    :: eim
            real    :: ein
            real    :: r1
            real    :: rs
            real    :: vdsize
   

        !!Step 0-3. Initialization       
            rhop  = rhoprt*1.0e-3                                  !Covert a unit from g/m^3 to kg/m^3
            vd   = 0.0    
            i    = mlu     

        !Step 1. Calculate aerodynamic resistance, Ra 
            
            if (zl >= 0) then
                ra = (0.74*log(z2/z0_f) + 4.7 *zl)/(0.4*ustar)                          !source of this equation unclear
            else
                ra = 0.74*(log(z2/z0_f) - 2*log((1+sqrt(1-9.*zl))*0.5))/(0.4*ustar)     !source of this equation unclear
            endif
    
            ra = amax1(ra,5.0)
    
            if (i == 1 .or. i == 3) then                                                !for water or inland lake
                ra = amin1(ra,2000.)
            else
                ra = amin1(ra,1000.)
            endif

        !Step 2. Calculate settling velocity, Vg 

            !!Step 2-1. Calculate leaf dimension 
            if (lai_ref(i,15) == lat_ref(i,14)) then    
                pllp = pllp2(i)
            else
                pllp = pllp2(i) - (lai_f-lai_ref(i,14))/(lai_ref(i,15)-lai_ref(i,14)+1.e-10)*(pllp2(i)-pllp1(i))
            endif

            !!Step 2-2. calculate radius of particle
            binsize  = diam/2

            !!Step 2-3. calculate air dynamic viscosity
            tave     = 0.5*(t2 + ts)
            amu      = 145.8 * 1.e-8 * tave ** 1.5/ (tave + 110.4)

            !!Step 2-4. calculate air kinetic viscosity
            anu      = amu/roarow
  
            !!Step 2-5. calculate cunningham slip correction factor (for samll particles) & relaxation time (=vg/gravity)
    
            prii     = 2.* 9.81 /(9. * amu)
            priiv    = prii * (rhop -roarow)
            vphil    = 0.
            cfac     = 1. + xmfp/binsize*(aa1+aa2*exp(-aa3*binsize/xmfp))            !Follow Zhang et al., (2001) eqn (3)
            taurel   = xmax1(priiv*binsize**2*cfac/9.81,0)                  

            !!Step 2-6. calculate stokes friction and diffusion coefficients         !Same as Wesely method for aerosol
            amob     = 6. * 3.14 * amu *binsize /cfac
            diff     = boltzk *tave /amob                                  
            schm     = anu/diff                                                      !Follow Zhang et al., (2001)

            !!Step 2-7. calculate gravitational settling velocity
            vg    = taurel * 9.81                                                    !Follow Zhang et al., (2001) eqn (2)


        !Step 3. Calculate surface resistance, Rs

            !!Step 3-1. Calculate Stokes number for different surface type
            
            !!!Step 3-1-1. smooth surface
            if (pllp <= 0) then
                st   = vg*ustar*ustar/anu                                            !Giorgi(1998), used in Zhang et al., (2001)            
            
            !!!Step 3-1-2. vegetated surface
            else
                st   = vg*ustar/(9.81* pllp/1000)                                    !Slinn(1982), used in Zhang et al.,(2001). where A = pllp/1000         
            endif

            !!Step 3-2. Calculate efficiency by diffusion(Eb), impaction(Eim), interception(Win) and particle rebound(R1)

            eb       = schm **(-gama(i))                                             !Follow Zhang et al.,(2001) eqn (6)
            
            eim      = (st/(st+aest(i)))**2                                          !Follow Zhang et al.,(2001) eqn (7c)
            
            ein      = 0.0
            if (pllp > 0.0) then
                ein  = (2.*binsize/(pllp/1000.))**2*0.5                              !Follow Zhang et al.,(2001) eqn (8)
            end if

            r1       = exp(-st**0.5)                                                 !Follow Zhang et al.,(2001) eqn (9)
            if (r1 < 0.5) then
                r1   = 0.5

            !!Step 3-3. Calculate surface resistance, Rs
            rs       = 1./(3.*ustar*r1*(eb + eim + ein))                             !Follow Zhang et al.,(2001) eqn (5). where epsilon = 3

        !Step 4. Calculate deposition velocity
            
            vd  = vg + 1./(ra + rs)                                                  !Follow Zhang et al.,(2001) eqn (1)
 

            return

            end subroutine aerosol_zhang

    end module vd_aerosol_zhang





















































































