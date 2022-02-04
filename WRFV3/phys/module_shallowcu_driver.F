!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_shallowcu_driver
CONTAINS
   SUBROUTINE shallowcu_driver(                                       &
                 ! memory reset (YLH)
                     reset_memory_flag,reset_memory_timestep          &  
                 ! record vars (YLH memory experiment)
                     ,record_vars_flag,record_timestep                &  
                 ! fix vars (YLH memory experiment)
                     ,fix_vars_flag                                   & 
                 ! Order dependent args for domain, mem, and tile dims
                     ,ids,ide, jds,jde, kds,kde                       &
                     ,ims,ime, jms,jme, kms,kme                       &
                     ,ips,ipe, jps,jpe, kps,kpe                       &
                     ,i_start,i_end,j_start,j_end,kts,kte,num_tiles   &
                 ! Order independent args (use VAR= in call)
                 ! --Prognostic
                     ,u,v,th,t                                        &
                     ,p,pi,rho,moist                                  &
                 ! --Other arguments
                     ,num_moist                                       &
                     ,itimestep,dt,dx,cudt,curr_secs,adapt_step_flag  &
                     ,rainsh,pratesh,nca,rainshv                      &
                     ,z,z_at_w,dz8w,mavail,pblh,p8w                   &
                     ,tke_pbl                                         &
                     ,cldfra,cldfra_old,cldfra_old_mp,cldfra_conv     &
                     ,cldfrash                                        &
                     ,htop,hbot                                       &
                 ! Package selection variables
                     ,shcu_physics                                    &
                 ! Optional moisture tracers 
                     ,qv_curr, qc_curr, qr_curr                       &
                     ,qi_curr, qs_curr, qg_curr                       & 
                     ,qnc_curr,qni_curr                               &
                 !+ YLH (org)
                     ,org                                             &
                     ,org_adv_curr                                    &
                     ,rorgshten                                       &
                     ,org2thpert,org2qpert,org2wpert,org2cbmf1        &
                     ,org2cbmf2,org2rkm,org2overlap                   &
                     ,cbmf1,cbmf2                                     &
                 !- YLH (org)

#if ( WRF_CHEM == 1 )
                     ,chem, chem_opt                                  &
#endif
                 ! Optional output arguments for CAMZM scheme
                     ,dlf, rliq, rliq2,dlf2  &
                     ,cmfmc, cmfmc2       &
                 ! Optional output arguments for CAMUW scheme
                     ,cush, snowsh, icwmrsh, rprdsh, cbmf, cmfsl      &
                     ,cmflq, evapcsh                                  &
                 ! Optional moisture and other tendencies
                     ,rqvshten,rqcshten,rqrshten                      &
                     ,rqishten,rqsshten,rqgshten                      &
                     ,rqcnshten,rqinshten                             &
                     ,rqvblten,rqvften                                &
                     ,rushten,rvshten                                 &
                     ,rthshten,rthraten,rthblten,rthften              &
                 ! Optional moisture tracer flags 
                     ,f_qv,f_qc,f_qr                                  &
                     ,f_qi,f_qs,f_qg                                  &
                     ,ht,shfrc3d,is_CAMMGMP_used                      &
                 ! for grims shallow convection with ysupbl
                     ,wstar,delta,kpbl,znu,raincv                   &
                 ! for nscv shallow convection
                     ,w,xland,hfx,qfx,mp_physics,pgcon                                           &
                                                                      )
!----------------------------------------------------------------------
   USE module_model_constants

!+ YLH   
   !USE module_state_description, ONLY: CAMUWSHCUSCHEME    &
   USE module_state_description, ONLY: CAMUWSHCUSCHEME,CAMUWORGSCHEME    &
!- YLH

# if (EM_CORE == 1)
                                       , CAMMGMPSCHEME    & 
# endif
                                       , G3SHCUSCHEME     & 
                                       , GRIMSSHCUSCHEME  &
                                       , NSCVSHCUSCHEME 
! *** add new modules of schemes here

   USE module_shcu_camuwshcu_driver, ONLY : camuwshcu_driver

!+ YLH   
   USE module_shcu_camuworg_driver, ONLY : camuworg_driver
!- YLH

   USE module_shcu_grims           , ONLY : grims
   USE module_shcu_nscv            , ONLY : shcu_nscv
   USE module_dm
   USE module_domain, ONLY: domain
#if ( WRF_CHEM == 1 )
      USE module_state_description, ONLY: num_chem
#endif

   !  This driver calls subroutines for the shallow cumulus
   !  parameterizations.
   !
   !  1. G3 shallow cumulus
   !  2. UW shallow cumulus from CAM
   !  3. GRIMs shallow cumulus from GRIMs (available only with ysupbl)
   !  4. NCEP shallow cumulus scheme (separated from NSAS cumulus scheme) 
   !                                 (Han and Pan 2011, YSU)
   !

!+ YLH
   ! Additionally, this module is also used to call to UW ORG scheme
   ! (Mapes and Neale 2011)
!- YLH

!----------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Rho_d      dry density (kg/m^3)
! Theta_m    moist potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- DT            time step (second)
!-- CUDT          cumulus time step (minute)
!-- curr_secs     current forecast time (seconds)
!-- itimestep     number of time step (integer)   
!-- DX            horizontal space interval (m)
!-- rr            dry air density (kg/m^3)
!
!-- RUSHTEN       Zonal wind tendency due to shallow
!                 cumulus scheme precipitation (m/s/s)
!-- RVSHTEN       Meridional wind tendency due to 
!                 cumulus scheme precipitation (m/s/s)
!-- RTHSHTEN      Theta tendency due to shallow
!                 cumulus scheme precipitation (K/s)
!-- RQVSHTEN      Qv tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQRSHTEN      Qr tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQCSHTEN      Qc tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQSSHTEN      Qs tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQISHTEN      Qi tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!-- RQGSHTEN      Qg tendency due to shallow
!                 cumulus scheme precipitation (kg/kg/s)
!
!-- RAINSH        accumulated total shallow cumulus scheme precipitation (mm)
!-- RAINSHV       time-step shallow cumulus scheme precipitation (mm)
!-- PRATESH       precipitiation rate from shallow cumulus scheme (mm/s)
!-- NCA           counter of the cloud relaxation 
!                 time in KF cumulus scheme (integer)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- t_phy         temperature (K)
!-- tke_pbl       turbulent kinetic energy from PBL scheme (m2/s2)
!-- w             vertical velocity (m/s)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- z             height above sea level at middle of layers (m)
!-- z_at_w        height above sea level at layer interfaces (m)
!-- dz8w          dz between full levels (m)
!-- pblh          planetary boundary layer height (m)
!-- mavail        soil moisture availability
!-- p8w           pressure at full levels (Pa)
!-- p_phy         pressure (Pa)
!-- pi_phy        the exner function, (p/p0)**(R/Cp) (dimensionless)
!                 points (dimensionless)
!-- hfx           upward heat flux at surface (W/m2)
!-- RTHRATEN      radiative temp forcing for Grell-Devenyi scheme
!-- RTHBLTEN      PBL temp forcing for Grell-Devenyi scheme
!-- RQVBLTEN      PBL moisture forcing for Grell-Devenyi scheme
!-- RTHFTEN
!-- RQVFTEN
!-- cldfra        cloud fraction
!-- cldfra_old    cloud fraction from previous time step
!-- cldfrash      cloud fraction from shallow Cu
!-- cldfra_old_mp cloud fraction from previous time step if CAMMGMP microphysics is used
!-- cldfra_conv   convective cloud fraction
!-- rho           density (kg/m^3)
!-- XLV0          latent heat of vaporization constant 
!                 used in temperature dependent formula (J/kg)
!-- XLV1          latent heat of vaporization constant 
!                 used in temperature dependent formula (J/kg/K)
!-- XLS0          latent heat of sublimation constant 
!                 used in temperature dependent formula (J/kg)
!-- XLS1          latent heat of sublimation constant
!                 used in temperature dependent formula (J/kg/K)
!-- R_d           gas constant for dry air ( 287. J/kg/K)
!-- R_v           gas constant for water vapor (461 J/k/kg)
!-- Cp            specific heat at constant pressure (1004 J/k/kg)
!-- rvovrd        R_v divided by R_d (dimensionless)
!-- G             acceleration due to gravity (m/s^2)
!-- EP_1          constant for virtual temperature 
!                 (R_v/R_d - 1) (dimensionless)
!--shfrc3d        Shallow cloud fraction
!-- cmfmc         Deep + Shallow Convective m
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- HBOT          index of lowest model layer with convection
!-- HTOP          index of highest model layer with convection
!-- LBOT          index of lowest model layer with convection
!-- LTOP          index of highest model layer with convection
!-- periodic_x    T/F this is using periodic lateral boundaries in the X direction
!-- periodic_y    T/F this is using periodic lateral boundaries in the Y-direction
!
!======================================================================
   LOGICAL,      INTENT(IN   )    :: is_CAMMGMP_used !BSINGH:01/31/2013: Added for CAMUWSHCU
   INTEGER,      INTENT(IN   )    ::                             &
                                      ids,ide, jds,jde, kds,kde, &
                                      ims,ime, jms,jme, kms,kme, &
                                                        kts,kte, &
                                      itimestep, num_tiles
#if ( WRF_CHEM == 1 )
   INTEGER,      INTENT(IN   )    :: chem_opt
#endif

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   ! YLH change to inout
   INTEGER,      INTENT(INOUT   )    ::                             &
                  num_moist

   INTEGER,      INTENT(IN   )    ::               shcu_physics

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),      &
         INTENT(INOUT)  ::                                       &
                                                          moist
#if ( WRF_CHEM == 1 )
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),       &
         INTENT(INOUT)  ::                                    &
                                                          chem
#endif
   

   ! YLH change to inout
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(INOUT ) ::                                       &
                                                         cldfra  &
                                                    ,cldfra_old  &
                                                  ,cldfra_old_mp &
                                                    ,cldfra_conv &
                                                      ,       z  &
                                                      ,  z_at_w  &
                                                      ,    dz8w  &
                                                      ,     p8w  &
                                                      ,       p  &
                                                      ,      pi  &
                                                      ,       u  &
                                                      ,       v  &
                                                      ,      th  &
                                                      ,       t  &
                                                      , tke_pbl  &
                                                      ,     rho


   ! YLH change to inout
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) ::        &
                  MAVAIL,PBLH,ht

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
          INTENT(INOUT) ::                               RAINSH  &
                                                    ,       NCA  & 
                                                    ,      HTOP  & 
                                                    ,      HBOT
 

   REAL, DIMENSION( ims:ime , jms:jme ),INTENT(INOUT),OPTIONAL :: &
        PRATESH, RAINSHV
   REAL, DIMENSION( ims:ime , jms:jme ) :: tmppratesh
                                                    
   REAL,  INTENT(IN   ) :: DT, DX
   INTEGER,      INTENT(IN   ),OPTIONAL    ::                             &
                                      ips,ipe, jps,jpe, kps,kpe
   REAL,  INTENT(IN   ),OPTIONAL :: CUDT
   REAL,  INTENT(IN   ),OPTIONAL :: CURR_SECS
   LOGICAL,INTENT(IN   ),OPTIONAL    ::     adapt_step_flag
   REAL   :: cudt_pass, curr_secs_pass
   LOGICAL :: adapt_step_flag_pass

!
! optional arguments
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr                  & 
                     ! optional scalar tracers !BSINGH
                     ,qnc_curr,qni_curr                          &
                     ! YLH
                     ,org_adv_curr                               &
                      ! optional moisture and other tendencies
                     ,rqvshten,rqcshten,rqrshten                 &
                     ,rqishten,rqsshten,rqgshten                 &
                     ,rqcnshten,rqinshten                        &
                     ,rqvblten,rqvften                           &
                     ,rthraten,rthblten                          &
                     ,rthften,rushten,rvshten,rthshten

!+ YLH
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                               &
                      rorgshten

    REAL, DIMENSION( ims:ime , jms:jme ),                        &
         OPTIONAL, INTENT(INOUT) ::                              &
                            org

    REAL, INTENT(IN) :: org2thpert,org2qpert,org2wpert,org2cbmf1,  &
                        org2cbmf2,org2rkm,org2overlap,cbmf1,cbmf2

!- YLH

   REAL, DIMENSION( ims:ime , jms:jme ),                         &
                    OPTIONAL, INTENT(INOUT) ::                   &
                rliq, rliq2 &
               ,cbmf, cush, snowsh
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                  cldfrash, cmfsl, cmflq, icwmrsh,               &
                  dlf, evapcsh,                                  &
                  cmfmc, cmfmc2, rprdsh
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
        INTENT(OUT) ::                                 &
                  dlf2                                             ! Required by CAMMGMP Microphysics
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
        INTENT(OUT) ::                                 &
                 shfrc3d                                           ! Required by wet scavenging code in WRF_CHEM
! for grims shallow convection with ysupbl
!
   REAL, DIMENSION( ims:ime, jms:jme )                         , &
         OPTIONAL, INTENT(IN   )                 ::      wstar
   REAL, DIMENSION( ims:ime, jms:jme )                         , &
         OPTIONAL, INTENT(IN   )                 ::        delta
   REAL, DIMENSION( ims:ime, jms:jme )                         , &
         OPTIONAL, INTENT(IN   )                 ::       raincv
   REAL, DIMENSION( kms:kme )                                  , &
         OPTIONAL, INTENT(IN   )       ::                    znu
   INTEGER, DIMENSION( ims:ime , jms:jme )                     , &
         OPTIONAL, INTENT(IN)                    ::         kpbl
!
! for nscv shallow convection
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                , &
         INTENT(IN   ) ::                                      w
   REAL, DIMENSION( ims:ime, jms:jme )                         , &
         INTENT(IN   ) ::                                  xland
   REAL, DIMENSION( ims:ime, jms:jme )                         , &
         OPTIONAL, INTENT(IN   ) ::                      hfx,qfx
   INTEGER, INTENT(IN), OPTIONAL ::                   mp_physics
   REAL,    INTENT(IN), OPTIONAL ::                        pgcon
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg


! LOCAL  VAR

   INTEGER :: i,j,k,its,ite,jts,jte,ij
   CHARACTER(len=200) :: message


!-----------------------------------------------------------------

!#################################################################
! Vars for recording variables during run-time (YLH wrf-python)

 ! derived data types 
 type :: var_dim_0 ! single value (integer or real)
         integer :: arr_dim_0_int
         real    :: arr_dim_0_real
 end type

 type :: var_dim_1 ! only k (vertical dimension)
         real,dimension(75) :: arr_dim_1_real ! need to specify dimension, cannot use kms:kme
 end type

 type :: var_dim_2 ! only i and j (horizontal dimensions)
         integer,dimension(11,11) :: arr_dim_2_int
         real,dimension(11,11) :: arr_dim_2_real
 end type

 type :: var_dim_3 ! i,j and k (3 dimensions)
         integer,dimension(11,75,11) :: arr_dim_3_int
         real,dimension(11,75,11) :: arr_dim_3_real
 end type

 type :: var_str
     character(:), allocatable :: arr_str
 end type

 ! config vars for recording
 logical,intent(in)           :: record_vars_flag,fix_vars_flag
 logical,intent(in)           :: reset_memory_flag
 integer,intent(in)           :: record_timestep ! timestep to record vars
 integer,intent(in)           :: reset_memory_timestep 

 ! vars
 type(var_dim_0), allocatable :: var_array_dim_0_int(:)
 type(var_dim_0), allocatable :: var_array_dim_0_real(:)
 type(var_dim_1), allocatable :: var_array_dim_1_real(:)
 type(var_dim_2), allocatable :: var_array_dim_2_int(:)
 type(var_dim_2), allocatable :: var_array_dim_2_real(:)
 type(var_dim_3), allocatable :: var_array_dim_3_int(:)
 type(var_dim_3), allocatable :: var_array_dim_3_real(:)

 integer                                     :: target_var_dim_0_int
 real                                        :: target_var_dim_0_real
 real,dimension(kms:kme)                     :: target_var_dim_1_real
 integer,dimension(ims:ime,jms:jme)          :: target_var_dim_2_int
 real,dimension(ims:ime,jms:jme)             :: target_var_dim_2_real
 integer,dimension(ims:ime,kms:kme,jms:jme)  :: target_var_dim_3_int
 real,dimension(ims:ime,kms:kme,jms:jme)     :: target_var_dim_3_real
 character(:), allocatable                   :: target_var_str

 real,dimension(ims:ime,kms:kme,jms:jme) :: target_var1,target_var2,target_var3
 integer :: var_array_length
 integer :: iTarget,aa

 ! vars for writing to file
 type(var_str), allocatable :: var_array_names_dim_0_int(:)      &
                              ,var_array_names_dim_0_real(:)     &
                              ,var_array_names_dim_1_real(:)     &
                              ,var_array_names_dim_2_int(:)      &
                              ,var_array_names_dim_2_real(:)     &
                              ,var_array_names_dim_3_int(:)      &
                              ,var_array_names_dim_3_real(:) 

!#################################################################

   if (.not. PRESENT(CURR_SECS)) then
      curr_secs_pass = -1
   else 
      curr_secs_pass = curr_secs
   endif

   if (.not. PRESENT(CUDT)) then
      cudt_pass = -1
   else
      cudt_pass = cudt
   endif

   if (.not. PRESENT(adapt_step_flag)) then
      adapt_step_flag_pass = .false.
   else
      adapt_step_flag_pass = adapt_step_flag
   endif

   ! Initialize tmppratesh to pratesh

   if ( PRESENT ( pratesh ) ) then
      tmppratesh(:,:) = pratesh(:,:)
   else
      tmppratesh(:,:) = 0.
   end if
   

   IF (shcu_physics .eq. 0) return
   
   ! YLH
   !print * ,"YLH num_moist:", num_moist 
   !print * ,"YLH itimestep:", itimestep
   
   
! DON'T JUDGE TIME STEP HERE, SINCE KF NEEDS ACCUMULATED W FIELD.
! DO IT INSIDE THE INDIVIDUAL CUMULUS SCHEME

! SET START AND END POINTS FOR TILES
!$OMP PARALLEL DO   &
!$OMP PRIVATE ( ij ,its,ite,jts,jte, i,j,k)
   DO ij = 1 , num_tiles
      its = i_start(ij)
      ite = i_end(ij)
      jts = j_start(ij)
      jte = j_end(ij)

!*********************************** YLH (MEMORY EXPERIMENTS) ***********************************!
!
! RECORD CONFIG VARS AND CONSTANTS 

if (record_vars_flag .eqv. .true.) then

    if (itimestep .eq. record_timestep) then

        print *,"#--------------START RECORDING-------------#"    
        print *,"YLH, itimestep = ",itimestep
    
        ! record config vars
        print *,"###START###,config_vars"
    
        print *,"dt = ",dt
        print *,"itimestep = ",itimestep
        print *,"dx = ",dx
        print *,"ids = ",ids
        print *,"ide = ",ide
        print *,"jds = ",jds
        print *,"jde = ",jde
        print *,"kds = ",kds
        print *,"kde = ",kde
        print *,"ims = ",ims
        print *,"ime = ",ime
        print *,"jms = ",jms
        print *,"jme = ",jme
        print *,"kms = ",kms
        print *,"kme = ",kme
        print *,"its = ",its ! local var
        print *,"ite = ",ite ! local var
        print *,"jts = ",jts ! local var
        print *,"jte = ",jte ! local var
        print *,"kts = ",kts
        print *,"kts = ",kte
        print *,"f_qv = ",f_qv
        print *,"f_qc = ",f_qc
        print *,"f_qr = ",f_qr
        print *,"f_qi = ",f_qi
        print *,"f_qs = ",f_qs
        print *,"adapt_step_flag_pass = ",adapt_step_flag_pass
    
        print *,"###END###,config_vars"
    
 
!#################################################################
! RECORD STATE VARIABLES 
    
    
        !################## single value vars (INTEGER) ###################!  

        var_array_length = 1
        allocate(var_array_dim_0_int(var_array_length))
        allocate(var_array_names_dim_0_int(var_array_length))

        !------------ DEFINE VARS HERE (single values integer) ------------!

        var_array_dim_0_int(1)%arr_dim_0_int        = num_moist
        var_array_names_dim_0_int(1)%arr_str        = "num_moist"

        !------------------------------------------------------------------!

        do iTarget = 1,var_array_length

           target_var_dim_0_int = var_array_dim_0_int(iTarget)%arr_dim_0_int
           target_var_str = var_array_names_dim_0_int(iTarget)%arr_str

           print *,"###START###,",target_var_str

           write(*,*) target_var_dim_0_int

           print *,"###END###,",target_var_str

        end do

        !################## single value vars (REAL) ###################!

        !var_array_length = 1
        !allocate(var_array_dim_0_real(var_array_length))
        !allocate(var_array_names_dim_0_real(var_array_length))

        !!-------------- DEFINE VARS HERE (single values real) -------------!

        !var_array_dim_0_real(1)%arr_dim_0_real        = thBinSize
        !var_array_names_dim_0_real(1)%arr_str         = "thBinSize"

        !!------------------------------------------------------------------!

        !do iTarget = 1,var_array_length

        !   target_var_dim_0_real = var_array_dim_0_real(iTarget)%arr_dim_0_real
        !   target_var_str = var_array_names_dim_0_real(iTarget)%arr_str

        !   print *,"###START###,",target_var_str

        !   write(*,*) target_var_dim_0_real
        !   !write(*,100) target_var_dim_0_real

        !   print *,"###END###,",target_var_str

        !end do       

        !################### 1-dimensional variables (REAL) #####################!

        !var_array_length = 1
        !allocate(var_array_dim_1_real(var_array_length))
        !allocate(var_array_names_dim_1_real(var_array_length))

        !!-------------- DEFINE VARS HERE (1-dimensional vars REAL) ------------!

        !var_array_dim_1_real(1)%arr_dim_1_real          = znu
        !var_array_names_dim_1_real(1)%arr_str           = "znu"

        !!------------------------------------------------------------------!

        !do iTarget = 1,var_array_length

        !   target_var_dim_1_real = var_array_dim_1_real(iTarget)%arr_dim_1_real
        !   target_var_str = var_array_names_dim_1_real(iTarget)%arr_str

        !   print *,"###START###,",target_var_str

        !   do k = kms,kme

        !       write(*,*) target_var_dim_1_real(k)
        !       !write(*,100) target_var_dim_1_real(k)

        !   end do

        !   print *,"###END###,",target_var_str

        !end do

        !################### 2-dimensional variables (INTEGER)  #####################!

        !var_array_length = 1
        !allocate(var_array_dim_2_int(var_array_length))
        !allocate(var_array_names_dim_2_int(var_array_length))        
        !
        !!-------------- DEFINE VARS HERE (2-dimensional vars INT) ------------!

        !var_array_dim_2_int(1)%arr_dim_2_int           = k22_shallow
        !var_array_names_dim_2_int(1)%arr_str           = "k22_shallow"

        !!------------------------------------------------------------------!

        !do iTarget = 1,var_array_length

        !   target_var_dim_2_int = var_array_dim_2_int(iTarget)%arr_dim_2_int
        !   target_var_str = var_array_names_dim_2_int(iTarget)%arr_str

        !   print *,"###START###,",target_var_str

        !   !do j = jms,jme
        !   !    do i = ims,ime

        !   do j = 1,2
        !       do i = 1,2

        !           write(*,100,advance="no") target_var_dim_2_int(i,j)

        !       end do
        !       write (*,*) ''
        !   end do

        !   print *,"###END###,",target_var_str

        !end do


        !################### 2-dimensional variables (REAL)  #####################!

        var_array_length = 2
        allocate(var_array_dim_2_real(var_array_length))
        allocate(var_array_names_dim_2_real(var_array_length))

        !-------------- DEFINE VARS HERE (2-dimensional vars REAL) ------------!

        var_array_dim_2_real(1)%arr_dim_2_real          = pblh
        var_array_names_dim_2_real(1)%arr_str           = "pblh"
        var_array_dim_2_real(2)%arr_dim_2_real          = ht
        var_array_names_dim_2_real(2)%arr_str           = "ht"

        !------------------------------------------------------------------!

        do iTarget = 1,var_array_length

           target_var_dim_2_real = var_array_dim_2_real(iTarget)%arr_dim_2_real
           target_var_str = var_array_names_dim_2_real(iTarget)%arr_str

           print *,"###START###,",target_var_str

           !do j = jms,jme
           !    do i = ims,ime

           do j = 1,2
               do i = 1,2

                   write(*,100,advance="no") target_var_dim_2_real(i,j)

               end do
               write (*,*) ''
           end do

           print *,"###END###,",target_var_str

        end do


        !################### 3-dimensional variables (REAL) #####################!

        var_array_length = 19
        allocate(var_array_dim_3_real(var_array_length))
        allocate(var_array_names_dim_3_real(var_array_length))
    
        !-------------- DEFINE VARS HERE (3-dimensional vars REAL) ------------!

        var_array_dim_3_real(1)%arr_dim_3_real          = t
        var_array_names_dim_3_real(1)%arr_str           = "t"
        var_array_dim_3_real(2)%arr_dim_3_real          = u
        var_array_names_dim_3_real(2)%arr_str           = "u"
        var_array_dim_3_real(3)%arr_dim_3_real          = v
        var_array_names_dim_3_real(3)%arr_str           = "v"
        var_array_dim_3_real(4)%arr_dim_3_real          = p
        var_array_names_dim_3_real(4)%arr_str           = "p"
        var_array_dim_3_real(5)%arr_dim_3_real          = p8w
        var_array_names_dim_3_real(5)%arr_str           = "p8w"
        var_array_dim_3_real(6)%arr_dim_3_real          = pi
        var_array_names_dim_3_real(6)%arr_str           = "pi"        
        var_array_dim_3_real(7)%arr_dim_3_real          = dz8w
        var_array_names_dim_3_real(7)%arr_str           = "dz8w"
        var_array_dim_3_real(8)%arr_dim_3_real          = z
        var_array_names_dim_3_real(8)%arr_str           = "z"
        var_array_dim_3_real(9)%arr_dim_3_real          = z_at_w
        var_array_names_dim_3_real(9)%arr_str           = "z_at_w"
        var_array_dim_3_real(10)%arr_dim_3_real         = tke_pbl
        var_array_names_dim_3_real(10)%arr_str          = "tke_pbl"
        var_array_dim_3_real(11)%arr_dim_3_real         = cldfra
        var_array_names_dim_3_real(11)%arr_str          = "cldfra"
        var_array_dim_3_real(12)%arr_dim_3_real         = cldfra_old
        var_array_names_dim_3_real(12)%arr_str          = "cldfra_old" 
        var_array_dim_3_real(13)%arr_dim_3_real         = cldfra_old_mp
        var_array_names_dim_3_real(13)%arr_str          = "cldfra_old_mp"
        var_array_dim_3_real(14)%arr_dim_3_real         = cldfra_conv
        var_array_names_dim_3_real(14)%arr_str          = "cldfra_conv"
        var_array_dim_3_real(15)%arr_dim_3_real         = qv_curr
        var_array_names_dim_3_real(15)%arr_str          = "qv_curr"
        var_array_dim_3_real(16)%arr_dim_3_real         = qc_curr
        var_array_names_dim_3_real(16)%arr_str          = "qc_curr"
        var_array_dim_3_real(17)%arr_dim_3_real         = qi_curr
        var_array_names_dim_3_real(17)%arr_str          = "qi_curr"
        var_array_dim_3_real(18)%arr_dim_3_real         = qnc_curr
        var_array_names_dim_3_real(18)%arr_str          = "qnc_curr"
        var_array_dim_3_real(19)%arr_dim_3_real         = qni_curr
        var_array_names_dim_3_real(19)%arr_str          = "qni_curr"

        !------------------------------------------------------------------!

        do iTarget = 1,var_array_length
    
           target_var_dim_3_real = var_array_dim_3_real(iTarget)%arr_dim_3_real
           target_var_str = var_array_names_dim_3_real(iTarget)%arr_str
    
           print *,"###START###,",target_var_str
    
           !do j = jms,jme
           !    do k = kms,kme
           !       do i = ims,ime
    
           do j = 1,2
               do k = kms,kme
                   do i = 1,2

                       write(*,100,advance="no") target_var_dim_3_real(i,k,j)
    
                   end do
                   write (*,*) ''
               end do
           end do
    
           print *,"###END###,",target_var_str
    
        end do
    
        ! formats
        100 format(E15.6)

        ! Check if recorded variables are correct              
        print *,"T", t(1,2,1)
        print *, "th", th(1,2,1)
    
        target_var1 = var_array_dim_3_real(1)%arr_dim_3_real
        target_var2 = var_array_dim_3_real(2)%arr_dim_3_real
    
        print *,"CHECK"
        print *,target_var1(1,2,1)
        print *,target_var2(1,2,1)
    

        print *,"#--------------END RECORDING-------------#"

    endif     


endif

!*********************************** YLH (MEMORY EXPERIMENTS) ***********************************!

   scps_select: SELECT CASE(shcu_physics)

   CASE (G3SHCUSCHEME)
      ! This setting takes the place of ishallow in v3.1.1+

   CASE (CAMUWSHCUSCHEME)
      CALL wrf_debug(100,'in camuw_scps')
      IF(.not.f_qi)THEN
         WRITE( message , * ) 'This shallow cumulus option requires ice microphysics option: f_qi = ', f_qi
         CALL wrf_error_fatal ( message )
      ENDIF
      CALL camuwshcu_driver(                                             &
           !+ YLH (memory experiment)
            fix_vars_flag=fix_vars_flag                                  &
           ,record_timestep=record_timestep                              &
           !- YLH
           ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde              &
           ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme              &
           ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte              &
           ,NUM_MOIST=num_moist, DT=dt                                   &
           ,P=p, P8W=p8w, PI_PHY=pi                                      &
           ,Z=z, Z_AT_W=z_at_w, DZ8W=dz8w                                &
           ,T_PHY=t, U_PHY=u, V_PHY=v                                    &
           ,MOIST=moist, QV=qv_curr, QC=qc_curr, QI=qi_curr              &
           ,QNC=qnc_curr, QNI=qni_curr                                   &
#if ( WRF_CHEM == 1 )
           ,CHEM=chem, CHEM_OPT=chem_opt                                 &
#endif
           ,PBLH_IN=pblh, TKE_PBL=tke_pbl                                &
           ,CLDFRA=cldfra, CLDFRA_OLD=cldfra_old                         &
           ,CLDFRA_OLD_MP=cldfra_old_mp                                  &
           ,CLDFRA_CONV=cldfra_conv,IS_CAMMGMP_USED=is_CAMMGMP_used      &
           ,CLDFRASH=cldfrash                                            &
           ,CUSH_INOUT=cush, PRATESH=tmppratesh                          &
           ,SNOWSH=snowsh                                                &
           ,ICWMRSH=icwmrsh, CMFMC=cmfmc, CMFMC2_INOUT=cmfmc2            &
           ,RPRDSH_INOUT=rprdsh, CBMF_INOUT=cbmf                         &
           ,CMFSL=cmfsl, CMFLQ=cmflq, DLF=dlf,DLF2=dlf2                  & !DLF2 is required by CAMMGMP microphysics
           ,EVAPCSH_INOUT=evapcsh                                        &
           ,RLIQ=rliq, RLIQ2_INOUT=rliq2, CUBOT=hbot, CUTOP=htop         &
           ,RUSHTEN=rushten, RVSHTEN=rvshten, RTHSHTEN=rthshten          &
           ,RQVSHTEN=rqvshten, RQCSHTEN=rqcshten, RQRSHTEN=rqrshten      &
           ,RQISHTEN=rqishten, RQSSHTEN=rqsshten, RQGSHTEN=rqgshten      &
           ,RQCNSHTEN=rqcnshten,RQINSHTEN=rqinshten                      &
           ,HT=ht,SHFRC3D=shfrc3d,ITIMESTEP=itimestep                    &
                                                                         )


!# YLH (ORG CASE)      
   CASE (CAMUWORGSCHEME)
      CALL wrf_debug(100,'in camuw_org_shcu')
      IF(.not.f_qi)THEN
         WRITE( message , * ) 'This shallow cumulus option requires ice microphysics option: f_qi = ', f_qi
         CALL wrf_error_fatal ( message )
      ENDIF
      CALL camuworg_driver(                                             &
           !+ YLH (reset memory experiment)
           reset_memory_flag=reset_memory_flag                           &
           ,reset_memory_timestep=reset_memory_timestep                  &
           !+ YLH (memory experiment)
           ,fix_vars_flag=fix_vars_flag                                  &
           ,record_timestep=record_timestep                              & 
           !- YLH
           ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde              &
           ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme              &
           ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte              &
           ,NUM_MOIST=num_moist, DT=dt                                   &
           ,P=p, P8W=p8w, PI_PHY=pi                                      &
           ,Z=z, Z_AT_W=z_at_w, DZ8W=dz8w                                &
           ,T_PHY=t, U_PHY=u, V_PHY=v                                    &
           ,MOIST=moist, QV=qv_curr, QC=qc_curr, QI=qi_curr              &
           ,QNC=qnc_curr, QNI=qni_curr                                   & 
           !+ YLH (org)
           ,org=org                                                      &
           ,org_adv=org_adv_curr                                         &
           ,rorgshten=rorgshten                                          &
           ,org2thpert=org2thpert,org2qpert=org2qpert                    &
           ,org2wpert=org2wpert,org2cbmf1=org2cbmf1                      &
           ,org2cbmf2=org2cbmf2,org2rkm=org2rkm,org2overlap=org2overlap  &
           ,cbmf1=cbmf1,cbmf2=cbmf2                                      &
           !- YLH
#if ( WRF_CHEM == 1 )
           ,CHEM=chem, CHEM_OPT=chem_opt                                 &
#endif
           ,PBLH_IN=pblh, TKE_PBL=tke_pbl                                &
           ,CLDFRA=cldfra, CLDFRA_OLD=cldfra_old                         &
           ,CLDFRA_OLD_MP=cldfra_old_mp                                  &
           ,CLDFRA_CONV=cldfra_conv,IS_CAMMGMP_USED=is_CAMMGMP_used      &
           ,CLDFRASH=cldfrash                                            &
           ,CUSH_INOUT=cush, PRATESH=tmppratesh                          &
           ,SNOWSH=snowsh                                                &
           ,ICWMRSH=icwmrsh, CMFMC=cmfmc, CMFMC2_INOUT=cmfmc2            &
           ,RPRDSH_INOUT=rprdsh, CBMF_INOUT=cbmf                         &
           ,CMFSL=cmfsl, CMFLQ=cmflq, DLF=dlf,DLF2=dlf2                  & !DLF2 is required by CAMMGMP microphysics
           ,EVAPCSH_INOUT=evapcsh                                        &
           ,RLIQ=rliq, RLIQ2_INOUT=rliq2, CUBOT=hbot, CUTOP=htop         &
           ,RUSHTEN=rushten, RVSHTEN=rvshten, RTHSHTEN=rthshten          &
           ,RQVSHTEN=rqvshten, RQCSHTEN=rqcshten, RQRSHTEN=rqrshten      &
           ,RQISHTEN=rqishten, RQSSHTEN=rqsshten, RQGSHTEN=rqgshten      &
           ,RQCNSHTEN=rqcnshten,RQINSHTEN=rqinshten                      & 
           ,HT=ht,SHFRC3D=shfrc3d,ITIMESTEP=itimestep                    &
                                                                         )
!- YLH
                                                                 
   CASE (GRIMSSHCUSCHEME)
      CALL wrf_debug(100,'in grims_scps')
      IF ( PRESENT( wstar ) ) THEN
      CALL grims(                                                        &
            QV3D=qv_curr,T3D=t                                           &
           ,P3DI=p8w,P3D=p,PI3D=pi,Z3DI=Z_AT_W                           &
           ,WSTAR=wstar,HPBL=pblh,DELTA=delta                        &
           ,RTHSHTEN=rthshten,RQVSHTEN=rqvshten                          &
           ,DT=dt,G=g,XLV=xlv,RD=r_d,RV=r_v                        &
           ,RCP=rcp,P1000MB=p1000mb                                      &
           ,KPBL2D=kpbl,ZNU=znu,RAINCV=raincv                            &
           ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde              &
           ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme              &
           ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte              &
                                                                         )
      ENDIF

   CASE (NSCVSHCUSCHEME)
      CALL wrf_debug(100,'in nscv_scps')
      IF ( PRESENT ( QFX ) .AND. PRESENT( HFX ) ) THEN
         CALL SHCU_NSCV(                                        &
                DT=dt,P3DI=p8w,P3D=p,PI3D=pi                    &
               ,QC3D=QC_CURR,QI3D=QI_CURR,RHO3D=rho             &
               ,QV3D=QV_CURR,T3D=t                              &
               ,RAINCV=RAINCV                                   &
               ,XLAND=XLAND,DZ8W=dz8w,W=w,U3D=u,V3D=v           &
               ,HPBL=pblh,HFX=hfx,QFX=qfx                       &
               ,MP_PHYSICS=mp_physics                           &
               ,pgcon=pgcon                                     &
               ,CP=cp,CLIQ=cliq,CPV=cpv,G=g,XLV=xlv,R_D=r_d     &
               ,R_V=r_v,EP_1=ep_1,EP_2=EP_2                     &
               ,CICE=cice,XLS=xls,PSAT=psat                     &
               ,F_QI=f_qi,F_QC=f_qc                             &
               ,RTHSHTEN=RTHSHTEN,RQVSHTEN=RQVSHTEN             &
               ,RQCSHTEN=RQCSHTEN,RQISHTEN=RQISHTEN             &
               ,RUSHTEN=RUSHTEN,RVSHTEN=RVSHTEN                 &
               ,PRATESH=tmppratesh                              &
               ,HBOT=HBOT,HTOP=HTOP                             &
               ,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde &
               ,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme &
               ,ITS=its,ITE=ite,JTS=jts,JTE=jte,KTS=kts,KTE=kte &
                                                                )
      ELSE
         CALL wrf_error_fatal('Lacking arguments for SHCU_NSCV in shallow cumulus driver')
      ENDIF

   CASE DEFAULT 
      WRITE( message , * ) 'The shallow cumulus option does not exist: shcu_physics = ', shcu_physics
      CALL wrf_error_fatal ( message )

   END SELECT scps_select

   ENDDO
   !$OMP END PARALLEL DO

   !
   ! Copy pratesh back to output array, if necessary.
   !
   if (PRESENT(PRATESH)) then
      pratesh(:,:) = tmppratesh(:,:)
      if (PRESENT(RAINSHV)) then
         rainshv(:,:) = pratesh(:,:)*dt
      endif
   endif

   END SUBROUTINE shallowcu_driver

END MODULE module_shallowcu_driver
