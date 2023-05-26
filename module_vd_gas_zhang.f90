module vd_gas_zhang
    use vd_constants, only :: dair,dh2o,rhoh2o,diffh2o,rmin,rmax_zhang,CP_A,CP_B,lai_ref,gamst,gmag,Rac1,Rac2,  &
                              RcutdO,RcutwO,RgO,RcutdS,RgS,rsmin,brs,tmin,tmax,topt,bvpd,psi1,psi2,sdmax
    implicit none

    contains
  
    subroutine Gas_Zhang(z2,zl,z0_f,ustar,t2,ts,srad,rh,fcld,pp,coszen,mlu,snow,  &
                         io3,inh3,conc_nh3,henry,henso2,f0,diffrat,rscale,lat_f,vd,f_net)

    !developed by Dr.Beiming Tang based on CAMx code
    !NOAA ARL UFS project
    !05/24/2023

    !Modifications:

    !Input arguments:
    !    z2        met reference height                     (m)
    !    zl        Z/L stability parameter                  (1)
    !    z0_f      surface roughnesss                       (m)
    !    ustar     friction velocity                        (m/s)
    !    t2        temperature at z2                        (K)
    !    ts        surface temperature                      (K)
    !    srad      solar radiance                           (W/m^2)
    !    rh        relative humidity                        (0-1,NOT percentage)
    !    fcld      cloud fraction                           (0-1)
    !    pp        precipitation water content              (g/m^3)
    !    coszen    cosine of solar zenith angle             (1)
    !    mlu       land use index                           (category 1-26)
    !    snow      snow water equivalent                    (m)
    !    io3       ozone flag                               (1)
    !    inh3      ammonia flag                             (1)
    !    conc_nh3  ambient NH3 concentration                (ug/m^3)
    !    henry     Henry's Law constant                     (M/atm)
    !    henso2    Henry's Law constant of SO2              (M/atm)
    !    f0        Reactivity parameter                     (??)
    !    diffrat   Ratio of diffusivity                     (1=H20:Species)
    !    rscale    Rs scaling for extremely soluble gases   (1)
    !    lat_f     leaf area index                          (1)

    !Output arguments:
    !    vd        deposition velocity                      (m/s)
    !    f_net     net flux of NH3 exchange                 (ug/m^2-s)

    !Routines called:
    !    vd_constants

    !Called by:
    !    future code of UFS dry deposition

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


    !Step 0. Define variables of dry deposition code
    !!Step 0-1. Input and output variables
        real*8, intent(in)    :: z2
        real*8, intent(in)    :: zl
        real*8, intent(in)    :: z0_f
        real*8, intent(in)    :: ustar
        real*8, intent(in)    :: t2
        real*8 ,intent(in)    :: ts
        real*8, intent(in)    :: srad
        real*8, intent(in)    :: rh
        real*8, intent(in)    :: fcld
        real*8, intent(in)    :: pp
        real*8, intent(in)    :: coszen
        integer,intent(in)    :: mlu
        real*8, intent(in)    :: snow
        integer,intent(in)    :: io3
        integer,intent(in)    :: inh3
        real*8, intent(in)    :: conc_nh3
        real*8, intent(in)    :: henry
        real*8, intent(in)    :: henso2
        real*8, intent(in)    :: f0
        real*8, intent(in)    :: diffrat
        real*8, intent(in)    :: rscale
        real*8, intent(in)    :: lai_f
        real*8, intent(out)   :: vd
        real*8, intent(out)   :: f_net

    !!Step 0-2. Define local varibales
        real*8    :: sd                           !snow depth                     (cm) 
        real      :: volrat
        real*8    :: prec                         !hourly precipitation           (mm/hour)
        real*8    :: alpha                        !scaling factor based on SO2    (no unit)
        real*8    :: beta                         !scaling factor based on O3     (no unit)
        real*8    :: rm                           !mesophyll resistance           (s/m)
        real*8    :: di                           !gas diffusivity                (cm^2/s)
        integer   :: i                            !land use type index            (1)
        logical   :: is_rain                      !rain flag                      (1)
        logical   :: is_dew                       !dew flag                       (1)
        real*8    :: es                           !                               (mb)        
        real*8    :: temp                                                 
        real*8    :: ra                                                  
        real*8    :: rst                                                   
        real*8    :: rdu                                                  
        real*8    :: rdv                                                  
        real*8    :: ww                           
        real*8    :: rdm                           
        real*8    :: rdn
        real*8    :: rv
        real*8    :: rn
        real*8    :: ratio
        real*8    :: sv
        real*8    :: fv
        real*8    :: pardir
        real*8    :: pardif
        real*8    :: pshad
        real*8    :: psun
        real*8    :: rshad
        real*8    :: rsun
        real*8    :: gshad
        real*8    :: gsun
        real*8    :: fsun
        real*8    :: fshd
        real*8    :: gspar
        real*8    :: t
        real*8    :: bt
        real*8    :: gt
        real*8    :: d0             !              (kPa)
        real*8    :: gd
        real*8    :: psi
        real*8    :: gw
        real*8    :: coedew
        real*8    :: dg              !              (g/kg)
        real*8    :: usmin
        real*8    :: wst
        real*8    :: racz
        real*8    :: rgo_f
        real*8    :: rgs_f
        real*8    :: rcuto_f
        real*8    :: rcuts_f
        real*8    :: fsnow
        real*8    :: rsnows
        real*8    :: vi
        real*8    :: rb
        real*8    :: dvh2o
        real*8    :: rs
        real*8    :: rcut
        real*8    :: rg
        real*8    :: rc
        real*8    :: x_st                                  !Compensation point conc (ug/m^3)
        real*8    :: x_g                                   !Compensation point conc (ug/m^3)        
        real*8    :: x_c                                   !Compensation point conc (ug/m^3)
        real*8    :: Rabi                                  !Inverse resistance terms
        real*8    :: Rsi                                   !Inverse resistance terms
        real*8    :: Racgi                                 !Inverse resistance terms
        real*8    :: Rcuti                                 !Inverse resistance terms
        
        
    !Step 1. Initialization 
    !!Step 1-1. Calculate snow depth 
        sd          = snow*1000.                           !*100 for m->cm, *10 for SWE->depth

    !!Step 1-2. Calcualte precipitation rate 
        volrat      = pp/rhoh2o
        prec        = (vlorat/1.0e-7)**1.27

    !!Step 1-3. Calculate mesophyll resistance, Rm
        alpha       = henry/henso2
        beta        = f0
        rm          = 1./(henry/3000/ + 100.*f0)           !mesophyll resistance follow Wesely et al.,(1989) eqn (6)
        rm          = amax1(rmin, rm)
        rm          = amin1(rmax_zhang, rm)

    !!Step 1-4. Calculate diffusivity of gas species
        di          = diffh2o/diffrat                      !unit = m^2/s
        di          = di*1.e4                              !unit = cm^2/s
        
    !!Step 1-5. Define function for saturation vaport pressure 
        es(temp)    = 6.108*EXP(17.27 * (temp - 273.16)/(temp - 35.86))

    !!Step 1-6. Initialization land use type & assign initial vd value
        vd          = 0.0
        i           = mlu


    !Step 2. Calculate aerodynamic resistance above canopy, Ra
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
        end if


    !Step 3. Calculate stomatal resistance for water vapor only, all other species in step 9, Rst(water)
    !!Step 3-0. set big value for stomatal resistance when stoma are closed.
        rst = 99999.9

    !!Step 3-1. Calcualte direct and diffuse PAR from solar radiation and solar zenith angle    
        !only calculate stomatal resistance if there is solar radiation, leaf area index is NOT zero, and within reasnable temperature range
        if (srad >= 0.1 .and. ts < (tmax(i) + 273.15) .and. ts > (tmin(i) + 273.15) .and. lai_f >0.001 .and. coszen > 0.001) then
            rdu    = 600.*exp(-0.185/coszen)*coszen
            rdv    = 0.4*(600.-rdu)*coszen
            ww     = -log(coszen)/2.302585
            ww     = -1.195+0.4459*ww -0.0345*ww**2
            ww     = 1320*10**ww
            rdm    = (720.*exp(-0.06/coszen)-ww)*coszen
            rdn    = 0.6*(720 - rdm -ww)*coszen
            rv     = amax1(0.1,rdu+rdv)
            rn     = amax1(0.01,rdm + rdn)
            ratio  = amin1(0.9,srad/(rv+rn))
            sv     = ratio*rv
            fv     = amin1(0.99,(0.9-ratio)/0.7)
            fv     = amax1(0.01,rdu/rv*(1.0-fv**0.6667)                  ! fraction of PAR in the direct beam
            pardir = fv * sv                                             ! par from direct radiation
            pardif = sv - pardir                                         ! Par from diffuse radiation            

     !!Step 3-2. Calculate sublit and shaded leaf area, PAR for sublit and shaded leaves
         if (lai_f > 2.5 .and. srad > 200.) then
             pshad = pardif * exp(-0.5*lai_f**0.8) + 0.07*pardir *(1.1-0.1*lai_f)*exp(-coszen)
             psun  = pardir**0.8*0.5/coszen + pshad
         else
             pshad = pardif*exp(-0.5*lai_f**0.7)+0.07*pardir*(1.1-0.1*lai_f)*exp(-coszen)
             psun  = pardir *0.5/coszen + pshad
         end if

         rshad     = rsmin(i) + brs(i) *rsmin(i)/pshad
         rsun      = rsmin(i) + brs(i) *rsmin(i)/psun
         gshad     = 1./rshad
         gsun      = 1./rsun
         fsun      = 2.*coszen*(1.-exp(-0.5*lat_f/coszen)                  ! sublit leaf area
         fshd      = lai_f -fsun                                           ! shaded leaf area 

     !!Step 3-3. Stomatal conductance before including effects of temperature, vapor pressure deficit, and water stress.
         gspar     = fsun*gsun + fshd*gshad
         
         !!!Step 3-3-1. function for temparature effect
         t  = ts -273.15
         bt = (tmax(i) -topt(i))/(tmax(i) - tmin(i))
         gt = (tmax(i) - t)/(tmax(i) -topt(i))
         gt = gt **bt
         gt = gt*(T-tmin(i))/(topt(i)-tmin(i))

         !!!Step 3-3-2. function for vapor pressure deficit
         d0 = es(ts)*(1.-rh)/10.     
         gd = 1.- bvpd(i)*d0

         !!!Step 3-3-3. functiion for water stress effect
         psi = -0.72-0.0013*srad
         gw  = (psi-psi2(i))/(psi1(i)-psi2(i))
         if (gw > 1.) then
             gw = 1.0
         end if
         if (gw < 0.1) then
             gw = 0.1
         end if
         if (gd > 1.) then
             gd = 1.0
         end if
         if (gd < 0.1) then
             gd = 0.1
         end if

      !!Step 3-4. Stomatal resistance for water vapor
          rst = 1.0/(GSPAR*gt*gd*gw)
      end if
      
      
      !Step 4. Calculate fraction of stomatal blocking under wet conditions
      !!Step 4-1. Decide if dew or rain occurs
      if (fcld < 0.25) then
          coedew = 0.3
      else if (fcld >= 0.25 .and. fcld < 0.75) then
          coedew = 0.2
      else
          coedew = 0.1
      end if

      dq = 0.622/1000.*es(ts)*(1.-rh)*1000.     
      dq = amax1(0.0001,dq)
      usmin = 1.5/dq*coedew

      is_rain = .false.
      is_dew  = .false.
      if (ts > 273.15 .and. prec > 0.20) then
          is_rain = .true.
      else if (ts > 273.15 .and. ustar < usmin) then
          is_dew = .true.
      end if
      
      !!Step 4-2.Decide fraction of stomatal blocking due to wet conditions
      Wst = 0
      if ((is_dew .or. is_rain) .and. srad > 200.) then
          Wst = (srad -200.)/800.
          Wst = amin1(wst, 0.5)
      end if 


    !Step 5. Calculate in-canopy aerodynamic resistance, Rac
        Racz = Rac1(i) + (lai_f - lai_ref(i,14))/(lai_ref(i,15)-lai_ref(i,14)+1.D-10)*(rac2(i)-rac1(i))    !Zhang et al., (2003) eqn (7a)
        Racz = Racz*lai_f**0.25/(ustar*ustar)                                                              !Zhang et al.,(2003) eqn (7)


    !Step 6. Calculate ground resistace,Rg
    !!Step 6-1. Calcualte ground resistance for O3
        if (i >= 4 .and. ts < 272.15) then
            rgo_f = amin1(rgo(i)*2., rgo(i)*exp(0.2*(272.15-ts)))
        else
            rgo_f = rgo(i)
        end if

    !!Step 6-2. Calcualte ground resistance for SO2
        if (i == 2) then
            rgs_f = amin1(rgs(i)*(273.15-ts),500.))
            rgs_f = amax1(rgs(i),100.)
        else if (i >= 4 .and. is_rain) then
            rgs_f = 50.
        else if (i >= 4 .and. is_dew) then
            rgs_f = 100.
        else if (i >= 4 .and. ts< 272.15) then
            rgs_f = amin1(rgs(i)*2.,rgs(i)*exp(0.2*(272.15-ts)))
        else
            rgs_f = rgs(i)
        end if


    !Step 7. Calculate cuticle resistance,Rcut for o3 and so2
        if (rcutdo(i) <= -1) then
            rcuto_f = 1.e25
            rcuts_f = 1.e25
        else if (is_rain) then
            rcuto_f = rcutwo(i)/(lai_f**0.5*ustar)                    !Zhang et al., (2003) eqn (9b)
            rcuts_f = 50./(lai_f**0.5*ustar)
            rcuts_f = amax1(rcuts_f,20.)
        else if (is_dew) then
            rcuto_f = rcutwo(i)/(lai_f**0.5*ustar)                    !Zhang et al.,(2003) eqn (9b)
            rcuts_f = 100./(lai_f**0.5*ustar)
            rcuts_f = amax1(rcuts_f,20.)
        else if (ts < 272.15) then
            rcuto_f = rcutdo(i)/(exp(3*rh)*lai_f**0.25*ustar)         !Zhang et al., (2003) eqn (9a),Notice in paper use RH precentage, here use RH 0-1
            rcuts_f = rcutds(i)/(exp(3*rh)*lai_f**0.25*ustar)
            rcuto_f = amin1(rcuto_f*2.,rcuto_f*exp(0.2*(272.15-ts)))  !Zhang et al., (2003) eqn (10b)
            rcuts_f = amin1(rcuts_f*2.,rcuts_f*exp(0.2*(272.15-ts))) 
            rcuto_f = amax1(rcuto_f,100.)
            rcuts_f = amax1(rcuts_f,100.)
        else
            rcuto_f = rcutdo(i)/(exp(3*rh)*lai_f**0.25*ustar)
            rcuts_f = rcutds(i)/(exp(3*rh)*lai_f**0.25*ustar)
            rcuto_f = amax1(rcuto_f,100.)
            rcuts_f = amax1(rcust_f,100.)
        end if

     !Step 6&7 supplemental.if snow occurs, Rg and Rcut are adjusted by snow cover fraction
        fsnow = sd/sdmax(i)
        fsnow = amin1(1.0,fsnow)                                !snow cover fraction for leaves
        if (fsnow > 0.0001 .and. i >= 4) then
            rsnows = amin1(70.*(275.15-ts),500.)                !zhang et al., (2003) eqn (8a)
            rsnows = amax1(rsnows,100.)
            rcuts_f = 1.0/((1.-fsnow)/rcuts_f +fsnow/rsnows)    !zhang et al., (2003) eqn (10d)
            rcuto_f = 1.0/((1.-fsnow)/rcuto_f +fsnow/2000.)
            fsnow = amin1(1.0,fsnow*2.)                         !snow cover fraction for ground
            rgs_f = 1.0/((1.-fsnow)/rgs_f + fsnow/rsnows)
            rgo_f = 1.0/((1.-fsnow)/rgo_f + fsnow/2000.)
         end if
 
 
      !Step 8. Calculate quasi-lamina resistance, Rb
         vi = 145.8 * 1.e-4*(ts*0.5 + t2 *0.5)**1.5/(ts*0.5 + t2*0.5 +110.4)    !diffusivity for each gas species,source unclear
         rb = 5.*(v1/di)**0.666667/ustar


      !Step 9. Calculate stamatal resistance for each species from the ratio of diffusivity of water vapor to the gas species
         dvh2o = 0.001*ts**1.75 * sqrt((29.+18.)/29./18.)
         dvh2o = dvh2o/(dair**0.3333 + dh2o **0.3333)**2
         rs    = rst*dvh2o/di + rm
         
         
      !Step 10. rescale cuticle and ground resistances for each species
         rcut = 1./(alpha/rcuts_f + beta/rcuto_f)
         rg   = 1./(alpha/rgs_f + beta/rgo_f)


      !Step 11. bi-directional NH3 flux calculation
          if (inh3 == 1) then
              if (lai_f < 0.5) then       !set a minimum lai for stomatal exchange
                  x_st = 0.D0
              else
                  x_st = (cp_a/t2)*exp(-cp_b/t2)*gamst(i)
              end if

              if (sd > 0) then             !if there is snow cover
                  x_g = 0.D0
              else
                  x_g = (cp_a/ts)*exp(-cp_b/ts)*gamg(i)
              end if
             
              rabi  = 1.D0/(ra+rb)
              rsi   = (1.D0 - Wst)/Rs        ! need verify source??
              racgi = 1.D0/(racz + rg)
              rcuti = 1.D0/rcut

              if (i >= 13 .and. i <= 20) then
                  rc = rsi + racgi +rcuti
                  rc = amax1(10.0,1./rc)                      !set minimum surface resistance as 10 (s/m)
                  vd = 1.D0/(ra + rb + rc)
                  f_net = 0.D0
              else
                  x_c   = ((c_nh3*rabi) + (x_st*rsi) + (x_g*racgi))/(rabi +rsi +racgi + rcuti)
                  f_net = -(c_nh3 -x_c) * rabi                !net flux (ug/m^2-s), positive flux --> emissions
                  vd    = max(-f_net/c_nh3, 0.D0)
                  f_net = max(f_net,0.D0)
              end if
              return
           end if

      !Step 12. Calculate total surface resistance
         rc = (1.- Wst)/Rs + 1./(Racz + Rg) + 1./Rcut
         rc = amax1(10.0,1/rc)
         if (io3 == 1 .and. i == 1) then
             rc = 1./(1.e-4 + 5.e-6*henry *ustar *(ts-273.15)**3.)
             rc = amax1(rc, 1500.)
         end if
         rc = rc*rscalsp                                         !scale for extremely soluble gas

       !Step 13. calculate dry deposition velocity
           vd = 1./(ra+rb+rc)

    return
    
    end subroutine Gas_Zhang

end module vd_gas_zhang






























