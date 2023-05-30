module vd_constants
    implicit none

    ! DEFINE CONSTANTS USED IN CALCULATIONS 
    real, parameter :: vk = 0.4              ! Von Karman constant                       (1)
    real, parameter :: rmin = 1.0            ! minimum resistance                        (s/m)
    real, parameter :: rmax = 1.e5           ! maximum resistance                        (s/m)
    real, parameter :: xmfp = 6.5e-8         ! mean free path                            (??)
    real, parameter :: g = 9.8               ! gravitational acceleration constant       (1)
    real, parameter :: vabs = 1.81e-2        ! dynamic viscosity of air                  (g/m-s)
    real, parameter :: boltz = 1.38e-20      ! Boltzmann constant.                       (m^2*g/(s^2*K))
    real, parameter :: boltzk = 1.38e-23     ! Boltzmann constant                        (m^2*kg/(s^2*K))
    real, parameter :: pi = 3.1415927        ! Circle constant.                          (1)
    real, parameter :: vair = 1.5e-5         ! kinematic viscosity of air                (m^2/s)   
    real, parameter :: diffh2o = 2.3e-5      ! molecular diffusivity of water in air     (m^2/s)
    real, parameter :: roarow = 1.19         ! air density at 20C                        (kg/m^3)
    real, parameter :: aa1    = 1.257        !??
    real, parameter :: aa2    = 0.4          !??
    real, parameter :: aa3    = 1.1          !??
    real, parameter :: dair = 0.369*29.+6.29 !??
    real, parameter :: dh2o = 0.369*18.+6.29 !??
    real, parameter :: rhoh2o = 1.e6         ! water density                             (g/m^3)
    real, parameter :: rmax_zhang = 100.0    ! max resistance in zhang                   (s/m)
    real*8,parameter :: CP_A = 2.7457D15     ! Parameters for Bi-Di NH3 drydep, Whaley et al., (2018), &
                                             ! constant A for compensation point Eq (K-ug/m^3)
    real*8,parameter :: CP_B = 10378.D0      ! Parameters for Bi-Di NH3 drydep, Whaley et al., (2018), &
                                             ! constant B for compensation point Eq (K)
  
    ! Baseline resistances used in Wesely model  (11 land time * 5 seasons)
    real :: rj(55,1)
    real :: rlu(55,1)
    real :: rac(55,1)
    real :: rgss(55,1)
    real :: rgso(55,1)
    real :: rlcs(55,1)
    real :: rlco(55,1)
    
    data rj  /9999.,  60., 120.,  70., 130., 100.,9999.,9999., 80., 100., 150.,   &
              9999.,9999.,9999.,9999., 250., 500.,9999.,9999.,9999.,9999.,9999.,  &
              9999.,9999.,9999.,9999., 250., 500.,9999.,9999.,9999.,9999.,9999.,  &
              9999.,9999.,9999.,9999., 400., 800.,9999.,9999.,9999.,9999.,9999.,  &
              9999., 120.,240., 140., 250., 190.,9999.,9999., 160., 200., 300./ 
      
    data rlu /9999.,2000.,2000.,2000.,2000.,2000.,9999.,9999.,2500.,2000.,4000.,  &
              9999.,9000.,9000.,9000.,4000.,8000.,9999.,9999.,9000.,9000.,9000.,  &
              9999.,9999.,9000.,9000.,4000.,8000.,9999.,9999.,9000.,9000.,9000.,  &
              9999.,9999.,9999.,9999.,6000.,9000.,9999.,9999.,9000.,9000.,9000.,  &
              9999.,4000.,4000.,4000.,2000.,3000.,9999.,9999.,4000.,4000.,8000./ 
    
    data rac / 100., 200., 100.,2000.,2000.,2000.,0.001,0.001, 300., 150., 200.,  &
               100., 150., 100.,1500.,2000.,1700.,0.001,0.001, 200., 120., 140.,  &
               100.,  10., 100.,1000.,2000.,1500.,0.001,0.001, 100.,  50., 120.,  &
               100.,  10.,  10.,1000.,2000.,1500.,0.001,0.001,  50.,  10.,  50.,  &
               100.,  50.,  80.,1200.,2000.,1500.,0.001,0.001, 200.,  60., 120./ 
    
    data rgss /400., 150., 350., 500., 500., 100.,0.001,1000.,0.001, 220., 400.,  &
               400., 200., 350., 500., 500., 100.,0.001,1000.,0.001, 300., 400.,  &
               400., 150., 350., 500., 500., 200.,0.001,1000.,0.001, 200., 400.,  &
               100., 100., 100., 100., 100., 100.,0.001,1000., 100., 100.,  50.,  &
               500., 150., 350., 500., 500., 200.,0.001,1000.,0.001, 250., 400./ 
      
    data rgso /300., 150., 200., 200., 200., 300.,2000., 400.,1000., 180., 200.,  &
               300., 150., 200., 200., 200., 300.,2000., 400., 800., 180., 200.,  &
               300., 150., 200., 200., 200., 300.,2000., 400.,1000., 180., 200.,  &
               600.,3500.,3500.,3500.,3500.,3500.,2000., 400.,3500.,3500.,3500.,  &
               300., 150., 200., 200., 200., 300.,2000., 400.,1000., 180., 200./ 
    
    data rlcs /9999.,2000.,2000.,2000.,2000.,2000.,9999.,9999.,2500.,2000.,4000.,  &
               9999.,9000.,9000.,9000.,2000.,4000.,9999.,9999.,9000.,9000.,9000.,  &
               9999.,9999.,9000.,9000.,3000.,6000.,9999.,9999.,9000.,9000.,9000.,  &
               9999.,9999.,9999.,9000., 200., 400.,9999.,9999.,9000.,9999.,9000.,  &
               9999.,4000.,4000.,4000.,2000.,3000.,9999.,9999.,4000.,4000.,8000./ 
    
    data rlco /9999.,1000.,1000.,1000.,1000.,1000.,9999.,9999.,1000.,1000.,1000.,  &
               9999., 400., 400., 400.,1000., 600.,9999.,9999., 400., 400., 400.,  &
               9999.,1000., 400., 400.,1000., 600.,9999.,9999., 800., 600., 600.,  &
               9999.,1000.,1000., 400.,1500., 600.,9999.,9999., 800.,1000., 800.,  &
               9999.,1000., 500., 500.,1500., 700.,9999.,9999., 600., 800., 800./ 
    
    ! Define reference lai used in zhang model (gas&aerosol)          
    data lai_ref /0.0, 0.0, 0.0, 5.0, 6.0, 0.1, 0.1, 6.0, 4.0, 3.0,  &
                  0.5, 3.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 1.0, 4.0, 0.0, 3.0, 3.0,                      &!Jan
                  0.0, 0.0, 0.0, 5.0, 6.0, 0.1, 0.1, 6.0, 4.0, 3.0,  &
                  0.5, 3.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 1.0, 4.0, 0.0, 3.0, 3.0,                      &!Feb
                  0.0, 0.0, 0.0, 5.0, 6.0, 0.5, 0.5, 6.0, 4.0, 3.0,  &
                  1.0, 3.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 0.5, 4.0, 0.0, 3.0, 3.0,                      & !Mar
                  0.0, 0.0, 0.0, 5.0, 6.0, 1.0, 1.0, 6.0, 4.0, 3.0,  &
                  1.0, 3.0, 1.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.0,  &
                  0.1, 0.1, 4.0, 0.0, 4.0, 4.0,                      & !Apr
                  0.0, 0.0, 0.0, 5.0, 6.0, 2.0, 2.0, 6.0, 4.0, 3.0,  &
                  1.5, 3.0, 1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  &
                  0.5, 0.1, 4.0, 0.0, 4.5, 4.5,                      & !May
                  0.0, 0.0, 0.0, 5.0, 6.0, 4.0, 4.0, 6.0, 4.0, 3.0,  &
                  2.0, 3.0, 1.0, 0.5, 2.0, 2.5, 3.0, 2.0, 3.0, 1.0,  &
                  1.0, 0.1, 4.0, 0.0, 5.0, 5.0,                      & !Jun
                  0.0, 0.0, 0.0, 5.0, 6.0, 5.0, 5.0, 6.0, 4.0, 3.0,  &
                  3.0, 3.0, 1.0, 1.0, 3.0, 4.0, 4.0, 3.0, 4.0, 1.0,  &
                  1.0, 0.1, 4.0, 0.0, 5.0, 5.0,                      & !Jul
                  0.0, 0.0, 0.0, 5.0, 6.0, 5.0, 5.0, 6.0, 4.0, 3.0,  &
                  3.0, 3.0, 1.0, 2.0, 3.5, 5.0, 4.5, 3.5, 4.5, 1.0,  &
                  1.0, 1.0, 4.0, 0.0, 5.0, 5.0,                      & !Aug
                  0.0, 0.0, 0.0, 5.0, 6.0, 4.0, 4.0, 6.0, 4.0, 3.0,  &
                  2.0, 3.0, 1.0, 2.0, 4.0, 6.0, 5.0, 4.0, 5.0, 1.0,  &
                  1.0, 2.0, 4.0, 0.0, 4.0, 4.0,                      & !Sep
                  0.0, 0.0, 0.0, 5.0, 6.0, 2.0, 2.0, 6.0, 4.0, 3.0,  &
                  1.5, 3.0, 1.0, 1.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  1.0, 1.5, 4.0, 0.0, 3.0, 3.0,                      & !Oct
                  0.0, 0.0, 0.0, 5.0, 6.0, 1.0, 1.0, 6.0, 4.0, 3.0,  &
                  1.0, 3.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.4, 1.5, 4.0, 0.0, 3.0, 3.0,                      & !Nov
                  0.0, 0.0, 0.0, 5.0, 6.0, 0.1, 0.1, 6.0, 4.0, 3.0,  &
                  0.5, 3.0, 1.0, 1.0, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 1.0, 4.0, 0.0, 3.0, 3.0,                      & !Dec
                  0.0, 0.0, 0.0, 5.0, 6.0, 0.1, 0.1, 6.0, 4.0, 3.0,  &
                  0.5, 3.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 1.0, 4.0, 0.0, 3.0, 3.0,                      & !Jan
                  0.0, 0.0, 0.0, 5.0, 6.0, 0.1, 0.1, 6.0, 4.0, 3.0,  &
                  0.5, 3.0, 1.0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0,  &
                  0.1, 0.1, 4.0, 0.0, 3.0, 3.0,                      & !MIN
                  0.0, 0.0, 0.0, 5.0, 6.0, 5.0, 5.0, 6.0, 4.0, 3.0,  &
                  3.0, 3.0, 1.0, 2.0, 4.0, 6.0, 5.0, 4.0, 5.0, 1.0,  &
                  1.0, 2.0, 4.0, 0.0, 5.0, 5.0/                      & !MAX
                  
    ! Define array of constants used in Zhang aerosol model
    real,dimension(26)    :: aest                          !??
    real,dimension(26)    :: pllp1                         !mim leaf dimension for each LUC
    real,dimension(26)    :: pllp2                         !max leaf dimension for each LUC
    real,dimension(26)    :: gama                          !??

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
                      
                      
    ! Define array of constants used in Zhang gas model
    real*8,dimension(26)   :: gamst                      ! Stomatal emission potential for NH3 [(mol/L)/(mol/L)]
    real*8,dimension(26)   :: gamg                       ! Ground emission potential for NH3 [(mol/L)/(mol/L)]
    real,  dimension(26)   :: Rac1                       ! In-canopy aerodynamic resistance minimum for Rac0 (s/m)
    real,  dimension(26)   :: Rac2                       ! In-canopy aerodynamic resistance maximum for Rac0 (s/m)
    real,  dimension(26)   :: RcutdO                     ! dry cuticle resistance for o3 (s/m)
    real,  dimension(26)   :: RcutwO                     ! wet cuticle resistance for o3 (s/m)
    real,  dimension(26)   :: RgO                        ! ground resistance for o3 (s/m)
    real,  dimension(26)   :: RcutdS                     ! Dry cuticle resistance of so2 (s/m)    
    real,  dimension(26)   :: RgS                        ! Ground resistance for SO2 (s/m)
    real,  dimension(26)   :: rsmin                      ! minimum stomatal resistance 
    real,  dimension(26)   :: brs                        ! empirical light response coefficient
    real,  dimension(26)   :: tmin                       ! minimum temperature below which complete closure occurs (degree C)
    real,  dimension(26)   :: tmax                       ! maximum temperature above which complete closure occurs (degree C)
    real,  dimension(26)   :: topt                       ! optimum temperature of maximum stomatal opening (degree C)
    real,  dimension(26)   :: bvpd                       ! water-vapor-pressure deficit constant     (k/Pa)
    real,  dimension(26)   :: psi1                       ! parameter specify leaf-water potential dependency
    real,  dimension(26)   :: psi2                       ! parameter specify leaf-water potential dependency
    real,  dimension(26)   :: sdmax                      ! Maximum snow depth over which snow fraction for leaves is 1.0,  &
                                                         ! Snow fraction for ground is treated 2 times of that for leaves
    
    gamst  = (/  0.D0,    0.D0,    0.D0, 3000.D0, 3000.D0, 300.D0,  300.D0,  300.D0,  300.D0,  300.D0,    & !   LU# 01 - 10
               300.D0,  300.D0,  300.D0,  300.D0,  800.D0, 800.D0,  800.D0,  800.D0,  800.D0,  800.D0,    & !   LU# 11 - 20
                 0.D0,   20.D0,  100.D0,    0.D0, 1000.D0, 300.D0   /)                                      !   LU# 21 - 26           
  
    gamg   = (/  0.D0,    0.D0,    0.D0, 1000.D0, 1000.D0, 1000.D0,  200.D0,   20.D0,  500.D0,   20.D0,   &    
               200.D0,   20.D0, 2000.D0, 2000.D0,  800.D0,  800.D0,  800.D0,  800.D0,  800.D0,  800.D0,   &    
                 0.D0,   20.D0,   20.D0,    0.D0, 1000.D0,   20.D0   /)                                          
               
    Rac1   = (/  0  ,  0   ,  0   ,  100 ,  250 ,  60  ,  100 ,  300 ,  100 ,  60  ,    &
                20  ,  40  ,  20  ,  10  ,  10  ,  10  ,  10  ,  10  ,  10  ,  20  ,    &
                40  ,  0   ,  20  ,  0   ,  100 ,  100    /)
             
    Rac2   = (/ 0   ,  0   ,  0   ,  100 ,  250 ,  100 ,  250 ,  300 ,  100 ,  60  ,    &
                60  ,  40  ,  20  ,  40  ,  40  ,  40  ,  40  ,  50  ,  40  ,  20  ,    &
                40  ,  0   ,  20  ,  0   ,  100 ,  100    /)
    
    RcutdO = (/-999 , -999 , -999 , 4000 , 6000 , 4000 , 6000 , 6000 , 8000 , 6000 ,    &
               5000 , 5000 , 4000 , 4000 , 4000 , 4000 , 4000 , 5000 , 5000 , 4000 ,    &
               6000 , 8000 , 5000 , -999 , 4000 , 4000   /)
               
    RcutwO = (/-999 , -999 , -999 ,  200 ,  400 ,  200 ,  400 ,  400 ,  400 ,  400 ,    &
                300 ,  300 ,  200 ,  200 ,  200 ,  200 ,  200 ,  300 ,  300 ,  200 ,    &
                400 ,  400 ,  300 , -999 ,  200 ,  200    /)
                
    RgO    = (/2000 , 2000 , 2000 ,  200 ,  200 , 200 ,  200 ,  200 ,  200 ,  200 ,     &
                200 ,  200 ,  200 ,  200 ,  200 , 200 ,  200 ,  200 ,  200 ,  500 ,     &
                500 ,  500 ,  500 ,  500 ,  200 , 200    /)     
                
    RcutdS = (/-999 , -999 , -999 , 2000 , 2500 , 2000 , 2500 , 2500 , 6000 , 2000 ,    &
               2000 , 2000 , 1000 , 1000 , 1500 , 1500 , 2000 , 2000 , 2000 , 2000 ,    &
               4000 , 2000 , 1500 , -999 , 2500 , 2500   /)       
               
    RgS    = (/  20 ,   70 ,  20  ,  200 ,  100 ,  200 ,  200 ,  100 ,  300 ,  200 ,    &
                200 ,  200 ,  200 ,  200 ,  200 ,   50 ,  200 ,  200 ,  200 ,   50 ,    &
                300 ,  300 ,   50 ,  700 ,  200 ,  200    /)    
                
    rsmin  = (/-999 , -999 , -999 ,  250 ,  150 ,  250 ,  150 ,  150 ,  250 ,  150 ,    &
                150 ,  250 ,  150 ,  100 ,  120 ,  120 ,  120 ,  250 ,  125 ,  150 ,    &
                200 ,  150 ,  150 , -999 ,  150 ,  150    /)
                
    brs    = (/-999 , -999 , -999 ,   44 ,   40 ,   44 ,   43 ,   40 ,   44 ,   40 ,    &
                 44 ,   44 ,   50 ,   20 ,   40 ,   40 ,   50 ,   65 ,   65 ,   40 ,    & 
                 42 ,   25 ,   40 , -999 ,   44 ,   43     /)
                 
    tmin   = (/-999 , -999 , -999 ,   -5 ,    0 ,   -5 ,    0 ,    0 ,    0 ,    0 ,    & 
                 -5 ,    0 ,    5 ,    5 ,    5 ,    5 ,    5 ,    5 ,   10 ,    5 ,    & 
                  0 ,   -5 ,    0 , -999 ,   -3 ,    0      /)
                  
    tmax   = (/-999 , -999 , -999 ,   40 ,   45 ,   40 ,   45 ,   45 ,   45 ,   45 ,    & 
                 40 ,   45 ,   40 ,   45 ,   45 ,   45 ,   45 ,   45 ,   45 ,   45 ,    & 
                 45 ,   40 ,   45 , -999 ,   42 ,   45     /)
                 
    topt   = (/-999 , -999 , -999 ,   15 ,   30 ,   15 ,   27 ,   30 ,   25 ,   30 ,    & 
                 15 ,   25 ,   30 ,   25 ,   27 ,   27 ,   25 ,   25 ,   30 ,   25 ,    & 
                 22 ,   20 ,   20 , -999 ,   21 ,   25     /)
                 
    bvpd   = (/-999 , -999 , -999 ,  0.31,  0.27,  0.31,  0.36,  0.27,  0.31,  0.27,    &  
                0.27,  0.27,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,    &  
                0.31,  0.24,  0.27, -999 ,  0.34,  0.31   /)
                
    psi1   = (/-999 , -999 , -999 , -2.0 , -1.0 , -2.0 , -1.9 , -1.0 , -1.0 , -2.0 ,    &
               -2.0 , -2.0 , -1.5 , -1.5 , -1.5 , -1.5 , -1.5 , -1.5 , -1.5 , -1.5 ,    &
               -1.5 ,    0 , -1.5 , -999 , -2.0 , -2.0    /)
               
    psi2  = (/ -999 , -999 , -999 , -2.5 , -5.0 , -2.5 , -2.5 , -5.0 , -4.0 , -4.0 ,    &
               -4.0 , -3.5 , -2.5 , -2.5 , -2.5 ,-2.5  , -2.5 , -2.5 , -2.5 , -2.5 ,    &
               -3.0 , -1.5 , -2.5 , -999 , -2.5 ,-3.0     /)        
               
    SDMAX  = (/9999.,  10. , 9999.,  200.,  200.,  200.,  200.,  200.,  200.,   30.,    &
                 20.,  30. ,   10.,   20.,   10.,   10.,   10.,   10.,   10.,   10.,    &
                200.,  10. ,   10.,   10.,  200.,  200.    /)
                          
end module vd_constants
