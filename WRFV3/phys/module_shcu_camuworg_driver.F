MODULE module_shcu_camuworg_driver
  USE shr_kind_mod,    only: r8 => shr_kind_r8

!+ YLH
! ORG scheme (Mapes and Neale 2011) modified from original uwshcu_driver file
!- YLH

! Roughly based on convect_shallow_tend in convect_shallow.F90 from CAM
! but tailored for the UW shallow cumulus scheme.

  !-------------------------------------------
  !Future modifications and important warnings (BSINGH:02/01/2013- Notes from WIG):
  !===========================================
  !1. UWShCu is hard-wired for certain moisture variables that could cause trouble 
  !   depending on which MP is used 
  !2. Mixing for rain, snow, and graupel moist vars is commented out. 
  !   For other MP routines that treat this prognostically, we need to implement 
  !   ShCu mixing of these (and other possible) moist species.
  !3. Fractional occurrence of shallow convection is currently not calculated.
  !-------------------------------------------




  IMPLICIT NONE

  PRIVATE                  !Default to private
  PUBLIC :: &              !Public entities
       camuworg_driver

  CHARACTER (LEN=256) :: message ! YLH

CONTAINS

!------------------------------------------------------------------------
SUBROUTINE camuworg_driver(                                  &
              ! YLH
              reset_memory_flag,reset_memory_timestep         &
             ,fix_vars_flag,record_timestep                   & 
             ,ids,ide, jds,jde, kds,kde                       &
             ,ims,ime, jms,jme, kms,kme                       &
             ,its,ite, jts,jte, kts,kte                       &
             ,num_moist, dt                                   &
             ,p, p8w, pi_phy, z, z_at_w, dz8w                 &
             ,t_phy, u_phy, v_phy                             &
             ,moist, qv, qc, qi, qnc, qni                     &
             !+ YLH
             ,org                                             &
             ,org_adv                                         &
             ,rorgshten                                       &
             !- YLH
#if ( WRF_CHEM == 1 )
             ,chem, chem_opt                                  &
#endif
             ,pblh_in, tke_pbl, cldfra, cldfra_old            &
             ,cldfra_old_mp,cldfra_conv, is_CAMMGMP_used      &
             ,cldfrash                                        &
             ,cush_inout, pratesh, snowsh, icwmrsh    &
             ,cmfmc, cmfmc2_inout, rprdsh_inout, cbmf_inout   &
             ,cmfsl, cmflq, dlf, dlf2, evapcsh_inout          &
             ,rliq, rliq2_inout, cubot, cutop                 &
             ,rushten, rvshten, rthshten                      &
             ,rqvshten, rqcshten, rqrshten                    &
             ,rqishten, rqsshten, rqgshten                    &
             ,rqcnshten,rqinshten                             &
             ,ht, shfrc3d,itimestep                           &
!+ YLH             
             ,org2thpert,org2qpert,org2wpert,org2cbmf1        &
             ,org2cbmf2,org2rkm,org2overlap,cbmf1,cbmf2       &
                                                              )
!- YLH

! This routine is based on convect_shallow_tend in CAM. It handles the
! mapping of variables from the WRF to the CAM framework for the UW
! shallow convective parameterization.
!
! Author: William.Gustafson@pnl.gov, Jan. 2010
!------------------------------------------------------------------------
  USE module_state_description, only: param_first_scalar, &
                                      p_qc, p_qr, p_qi, p_qs, p_qg, p_qnc, p_qni, &
!+ YLH
                                      p_org_adv
!- YLH                                      
  USE module_cam_support,       only: pcols, pver, pcnst =>pcnst_runtime
#if ( WRF_CHEM == 1 )
  USE module_cam_support,       only:  cam_mam_aerosols
#endif
  USE constituents,             only: cnst_get_ind
  USE physconst,                only: latice,cpair, gravit, latvap

!+ YLH  
  !USE uwshcu,                   only: compute_uwshcu_inv
  USE uworg,                   only: compute_uworg_inv
!- YLH

  USE wv_saturation,            only: fqsatd

!+ YLH: comment out if statement  
!#if ( WRF_CHEM == 1 )
   use module_state_description,  only: num_chem, param_first_scalar
   use module_mp_cammgmp_driver,  only: physics_update, physics_ptend_init

  !use module_state_description,  only: num_chem, param_first_scalar,CBMZ_CAM_MAM3_NOAQ, &
  !     CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_NOAQ,CBMZ_CAM_MAM7_AQ
  !use module_data_cam_mam_asect, only: lptr_chem_to_q, factconv_chem_to_q
  !use module_mp_cammgmp_driver,  only: physics_update, physics_ptend_init
!#endif
!- YLH

! Subroutine arguments...
  LOGICAL, INTENT(IN) :: is_CAMMGMP_used
  INTEGER, INTENT(IN   ) ::    ids,ide, jds,jde, kds,kde,  &
                               ims,ime, jms,jme, kms,kme,  &
                               its,ite, jts,jte, kts,kte,  &
                               itimestep

  ! YLH change to inout
  INTEGER, INTENT(INOUT   ) ::    num_moist

#if ( WRF_CHEM == 1 )
  INTEGER, INTENT(IN   ) ::    chem_opt
#endif

 ! YLH change to inout
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ), INTENT(INOUT) :: &
                             moist    !moist tracer array

#if ( WRF_CHEM == 1 )
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ), INTENT(INOUT) :: &
                              chem    !moist tracer array
#endif

!+ YLH
  REAL, INTENT(IN) :: org2thpert,org2qpert,org2wpert,org2cbmf1,  &
                      org2cbmf2,org2rkm,org2overlap,cbmf1,cbmf2 
!- YLH

  ! YLH change to inout
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
                             cldfra, & !cloud fraction
                         cldfra_old, & !previous time step cloud fraction
                      cldfra_old_mp, &
                        cldfra_conv, &
                               dz8w, & !height between layer interface (m)
                                  p, & !pressure at mid-level (Pa)
                                p8w, & !pressure at level interface (Pa)
                             pi_phy, & !exner function, (p0/p)^(R/cpair) (none)
                                 qv, & !water vapor mixing ratio (kg/kg-dry air)
                              t_phy, & !temperature (K)
                            tke_pbl, & !turbulent kinetic energy from PBL (m2/s2)
                              u_phy, & !zonal wind component on T points (m/s)
                              v_phy, & !meridional wind component on T points (m/s)
                                  z, & !height above sea level at mid-level (m)
                             z_at_w    !height above sea level at interface (m)

  ! YLH (org)
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &  
                               org_adv ! Advected 3D org
 
  ! YLH (org)                     
  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
                                org    ! Org 2D
  
  ! YLH change to inout
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT), OPTIONAL :: &
                                 qc, & !cloud droplet mixing ratio (kg/kg-dry air)
                                 qi, & !cloud ice crystal mixing ratio (kg/kg-dry air)
                                qnc, & !cloud water  number concentration (#/kg)
                                qni    !cloud ice number concentration (#/kg)

  ! YLH change to inout
  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
                            pblh_in, & !height of PBL (m)
                            ht         !Terrain height (m)

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: &
                           cldfrash, & !shallow convective cloud fraction
                              cmfmc, & !deep+shalow cloud fraction (already contains deep part from ZM)
                       cmfmc2_inout, & !shallow cloud fraction
                              cmflq, & !convective flux of total water in energy unit (~units)
                              cmfsl, & !convective flux of liquid water static energy (~units)
                                dlf, & !dq/dt due to export of cloud water (input=deep from ZM, output=deep+shallow)
                      evapcsh_inout, & !output array for evaporation of shallow convection precipitation (kg/kg/s)
                            icwmrsh, & !shallow cumulus in-cloud water mixing ratio (kg/m2)
                       rprdsh_inout, & !dq/dt due to deep(~?) & shallow convective rainout (~units?)
                            rushten, & !UNcoupled zonal wind tend from shallow Cu scheme (m/s2)
                            rvshten, & !UNcoupled meridional wind tend from shallow Cu scheme (m/s2)
                           rthshten, & !UNcoupled potential temperature tendendcy from shallow cu scheme (K/s)
                           rqvshten, & !UNcoupled water vapor mixing ratio tend from shallow Cu scheme (kg/kg/s)
                           rqcshten, & !UNcoupled clod droplet mixing ratio tend from shallow Cu scheme (kg/kg/s)
                           rqrshten, & !UNcoupled raindrop mixing ratio tend from shallow Cu scheme (kg/kg/s)
                           rqishten, & !UNcoupled ice crystal mixing ratio tend from shallow Cu scheme (kg/kg/s)
                           rqsshten, & !UNcoupled snow mixing ratio tend from shallow Cu scheme (kg/kg/s)
                           rqgshten, & !UNcoupled graupel mixing ratio tend from shallow Cu scheme (kg/kg/s)
                          rqcnshten, & !PMA
                          rqinshten

  !+ YLH 
  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) :: rorgshten ! org_adv tendency
  !- YLH

  REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: &
                         cbmf_inout, & !cloud base mass flux (kg/m2/s)
                              cubot, & !level number of base of convection
                              cutop, & !level number of top of convection
                         cush_inout, & !convective scale height (~units?)
                            pratesh, & !time-step shallow cumulus precip rate at surface (mm/s)
                               rliq, & !vertically-integrated reserved cloud condensate (m/s)
                        rliq2_inout, & !vertically-integrated reserved cloud condensate for shallow (m/s)
                             snowsh    !accumulated convective snow rate at surface for shallow Cu (m/s) ~are these the units we should use for WRF?

  REAL, INTENT(IN) :: &
                                 dt    !time step (s)

  REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT) ::   &
                        dlf2           ! Required by CAMMGMP Microphysics                            
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT) ::   &
                           shfrc3d     !Shallow cloud fraction

   ! YLH (memory experiment)
   LOGICAL, INTENT(IN) :: fix_vars_flag
   LOGICAL, INTENT(IN) :: reset_memory_flag
   INTEGER, INTENT(IN) :: record_timestep 
   INTEGER, INTENT(IN) :: reset_memory_timestep 

! Local variables...
  !Variables dimensioned for input to CAM routines
  REAL(r8), DIMENSION(pcols, kte, pcnst) ::  &
                             moist8, & !tracer array for CAM routines
                         tnd_tracer    !tracer tendency

!+ YLH: Placeholder array to hold original tnd_tracer values (tnd_tracer is reset in physics_update)
  REAL(r8), DIMENSION(pcols, kte, pcnst) :: tnd_tracer_ph
!- YLH

  REAL(r8), DIMENSION(pcols, kte+1) ::  &
                              pint8, & !pressure at layer interface (Pa)
                                zi8, & !height above the ground at interfaces (m)
                               tke8, & !turbulent kinetic energy at level interfaces (m2/s2)
                              slflx, & !convective liquid water static energy flux (~units?)
                              qtflx, & !convective total water flux (~units?)
                            flxprec, & ! Shallow convective-scale flux of precip (rain+snow) at interfaces [ kg/m2/s ]
                            flxsnow, & ! Shallow convective-scale flux of snow at interfaces [ kg/m2/s ]
                            cmfmc2     !cloud fraction

  ! YLH
  REAL(r8), DIMENSION(pcols) :: &
                               org8

  REAL(r8), DIMENSION(pcols, kte) ::  &
                               cld8, & !cloud fraction
                            cldold8, & !previous time step cloud fraction ~should this be just the convective part?
                             cmfdqs, & !convective snow production (~units?)
                            evapcsh, & !evaporation of shallow convection precipitation >= 0. (kg/kg/s)
                           iccmr_uw, & !in-cloud cumulus LWC+IWC (kg/m2)
                           icwmr_uw, & !in-cloud cumulus LWC (kg/m2)
                           icimr_uw, & !in-cloud cumulus IWC (kg/m2)
                              pdel8, & !pressure difference between layer interfaces (Pa)
                           pdeldry8, & !pressure difference between layer interfaces for dry atm (Pa)
                              pmid8, & !pressure at layer middle (Pa)
                                qc2, & !dq/dt due to export of cloud water
                                qh8, & !specific humidity (kg/kg-moist air)
                                qc8, & !cloud liquid water (~units?)
                                qi8, & !cloud ice (~units?)
                              qhtnd, & !specific humidity tendency (kg/kg/s)
                              qctnd, & !cloud water mixing ratio tendency
                              qitnd, & !cloud ice mixing ratio tendency
                             rprdsh, & !dq/dt due to deep(~?) & shallow convective rainout (~units?)
                                 s8, & !dry static energy (J/kg)
                              shfrc, & !shallow cloud fraction
                               stnd, & !heating rate (dry static energy tendency, W/kg)
                                 t8, & !temperature (K)
                                 u8, & !environment zonal wind (m/s)
                               utnd, & !zonal wind tendency (m/s2)
                                 v8, & !environment meridional wind (m/s)
                               vtnd, & !meridional wind tendency (m/s2)
                                zm8    !height between interfaces (m)

  REAL(r8), DIMENSION(pcols, kte) ::  &
                           qcten_det, &
                           qiten_det, &
                          qcnten_det, &
                          qinten_det, &
                           qsten_det

  REAL(r8), DIMENSION(pcols) ::  &
                               cbmf, & !cloud base mass flux (kg/m2/s)
                               cnb2, & !bottom level of convective activity
                               cnt2, & !top level of convective activity
                               cush, & !convective scale height (~units?)
                               pblh, & !pblh height (m)
                              precc, & !convective precip (rain+snow) at surface for shallow Cu (m/s)
                              rliq2, & !vertically-integrated reserved cloud condensate for shallow (m/s)
                               snow    !convective snow rate at surface (m/s)

  !Other local vars
  REAL(r8) :: ztodt,dum1
  INTEGER :: i, j, k, kflip, m, mp1
  INTEGER :: cnb, cnt      !index of cloud base and top in CAM world (indices decrease with height)
  INTEGER :: lchnk         !chunk identifier, used to map 2-D to 1-D arrays in WRF
  INTEGER :: ncnst         !number of tracers
  INTEGER :: ncol          !number of atmospheric columns in chunk
  CHARACTER(LEN=1000) :: msg

  character*24 :: ptend_name            !ptend%name in CAM5 - used in parameterization
  logical      :: ptend_ls              !ptend%ls   in CAM5 - used for calling physics_update subroutine
  logical      :: ptend_lq(pcnst)       !ptend%lq   in CAM5

  integer :: l, l2
  real(r8) :: state_s(pcols,kte)
  real(r8) :: ptend_s(pcols,kte)                   !Dummy arguments for physics_update call
  
!*********************** YLH ***********************!     

  LOGICAL :: end_of_file 
  real, dimension(kte) :: t_k,qv_k,qc_k,qi_k
  real, dimension(kte) :: qnc_k, qni_k
  real, dimension(kte) :: dz8w_k,p_k,p8w_k,pi_k 
  real, dimension(kte) :: z_k,z_at_w_k,tke_pbl_k 
  real, dimension(kte) :: u_k,v_k
  real, dimension(kte) :: cldfra_k,cldfra_old_k,   &
                          cldfra_old_mp_k,cldfra_conv_k
   
  real :: pblh_in_k,ht_k

  REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: t_ori,qv_ori,qc_ori,qi_ori,  &
                                                qnc_ori, qni_ori,dz8w_ori,p_ori,p8w_ori,pi_ori, &
                                                z_ori,z_at_w_ori,tke_pbl_ori,u_ori,v_ori, &
                                                cldfra_ori,cldfra_old_ori,   &
                                                cldfra_old_mp_ori,cldfra_conv_ori

  REAL, DIMENSION(ims:ime, jms:jme) ::  pblh_in_ori,ht_ori       


!*********************** YLH MEMORY RESET EXPERIMENT ***********************!

  !print *,"YLH org:",org(1,1),org(1,2),org(2,1),org(2,2)
 
  !print *,"YLH org_adv:"

  !do k = kts,kte
  !     print *,org_adv(1,k,1)
  !enddo

  if (reset_memory_flag .eqv. .true.) then

      if (itimestep .eq. reset_memory_timestep) then 

        !call wrf_debug ( 100, " -------------- YLH reset memory")
        !call wrf_debug ( 100, "itimestep:")
        !write( message, *) itimestep
        !call wrf_debug (100, message)

        print *,"YLH reset memory"
        print *,"timestep:",itimestep

        do j = jts,jte
           do i = its,ite
               
               org(i,j) = 0.0

           enddo
        enddo

        !print *,"org:",org(1,1),org(1,2),org(2,1),org(2,2)

        do j = jts,jte
            do k = kts,kte
                do i = its,ite

                    org_adv(i,k,j) = 0.0
                    !print *,"i,k,j,org_adv",i,k,j,org_adv(i,k,j)

                enddo
            enddo
        enddo

      endif
  endif


!*********************** YLH MEMORY RESET EXPERIMENT ***********************!


!*********************** YLH ***********************!

!--------------
! YLH FIX input state vars to test memory
!--------------      

if (fix_vars_flag .eqv. .true.) then

   if (itimestep .ge. record_timestep) then

       call wrf_debug ( 100, ' -------- YLH FIXING STARTS' )
       call wrf_debug ( 100, ' YLH BEFORE VALUES' )
       call wrf_debug(100,'k T QV')

       do j=jts,jte
          do k=kts,kte
             do i=its,ite

                write( message , * ) k,t_phy(i,k,j),qv(i,k,j)
                call wrf_debug ( 100, message )

                t_ori(i,k,j)=t_phy(i,k,j)
                qv_ori(i,k,j)=qv(i,k,j)
                qc_ori(i,k,j)=qc(i,k,j)
                qi_ori(i,k,j)=qi(i,k,j)
                qnc_ori(i,k,j)=qnc(i,k,j)
                qni_ori(i,k,j)=qni(i,k,j)
                dz8w_ori(i,k,j)=dz8w(i,k,j)
                p_ori(i,k,j)=p(i,k,j)
                p8w_ori(i,k,j)=p8w(i,k,j)
                pi_ori(i,k,j)=pi_phy(i,k,j)
                z_ori(i,k,j)=z(i,k,j)
                z_at_w_ori(i,k,j)=z_at_w(i,k,j)
                tke_pbl_ori(i,k,j)=tke_pbl(i,k,j)
                u_ori(i,k,j)=u_phy(i,k,j)
                v_ori(i,k,j)=v_phy(i,k,j)
                cldfra_ori(i,k,j)=cldfra(i,k,j)
                cldfra_old_ori(i,k,j)=cldfra_old(i,k,j)
                cldfra_old_mp_ori(i,k,j)=cldfra_old_mp(i,k,j)
                cldfra_conv_ori(i,k,j)=cldfra_conv(i,k,j)

             enddo
          enddo
       enddo

       do j=jts,jte
          do i=its,ite

                ht_ori(i,j)=ht(i,j)
                pblh_in_ori(i,j)=pblh_in(i,j)

           enddo
       enddo
          

       !-------------------- GRID 0 0
       call read_input00(kte,t_k,qv_k,qc_k,qi_k,qnc_k,qni_k, &
                         dz8w_k,p_k,p8w_k,pi_k,z_k,z_at_w_k, &
                         tke_pbl_k,u_k,v_k,cldfra_k,cldfra_old_k, &
                         cldfra_old_mp_k,cldfra_conv_k, &
                         ht_k,pblh_in_k)

       !do j=jts,jte
       !   do i=its,ite
             do k=kts,kte

                t_phy(1,k,1) = t_k(k)
                qv(1,k,1) = qv_k(k)
                qc(1,k,1) = qc_k(k)
                qi(1,k,1) = qi_k(k)
                qnc(1,k,1) = qnc_k(k)
                qni(1,k,1) = qni_k(k)          
                dz8w(1,k,1) = dz8w_k(k)
                p(1,k,1) = p_k(k)
                p8w(1,k,1) = p8w_k(k)
                z(1,k,1) = z_k(k)
                z_at_w(1,k,1) = z_at_w_k(k)
                pi_phy(1,k,1) = pi_k(k)
                tke_pbl(1,k,1) = tke_pbl_k(k)
                u_phy(1,k,1) = u_k(k)
                v_phy(1,k,1) = v_k(k)
                cldfra(1,k,1) = cldfra_k(k)
                cldfra_old(1,k,1) = cldfra_old_k(k)
                cldfra_old_mp(1,k,1) = cldfra_old_mp_k(k)
                cldfra_conv(1,k,1) = cldfra_conv_k(k)

             enddo
       !   enddo
       !enddo

       !do j=jts,jte
       !   do i=its,ite

                pblh_in(1,1) = pblh_in_k
                ht(1,1) = ht_k

       !   enddo
       !enddo

       !-------------------- GRID 0 1
       call read_input01(kte,t_k,qv_k,qc_k,qi_k,qnc_k,qni_k, &
                         dz8w_k,p_k,p8w_k,pi_k,z_k,z_at_w_k, &
                         tke_pbl_k,u_k,v_k,cldfra_k,cldfra_old_k, &
                         cldfra_old_mp_k,cldfra_conv_k, &
                         ht_k,pblh_in_k)

       !do j=jts,jte
       !   do i=its,ite
             do k=kts,kte

                t_phy(1,k,2) = t_k(k)
                qv(1,k,2) = qv_k(k)
                qc(1,k,2) = qc_k(k)
                qi(1,k,2) = qi_k(k)
                qnc(1,k,2) = qnc_k(k)
                qni(1,k,2) = qni_k(k)          
                dz8w(1,k,2) = dz8w_k(k)
                p(1,k,2) = p_k(k)
                p8w(1,k,2) = p8w_k(k)
                z(1,k,2) = z_k(k)
                z_at_w(1,k,2) = z_at_w_k(k)
                pi_phy(1,k,2) = pi_k(k)
                tke_pbl(1,k,2) = tke_pbl_k(k)
                u_phy(1,k,2) = u_k(k)
                v_phy(1,k,2) = v_k(k)
                cldfra(1,k,2) = cldfra_k(k)
                cldfra_old(1,k,2) = cldfra_old_k(k)
                cldfra_old_mp(1,k,2) = cldfra_old_mp_k(k)
                cldfra_conv(1,k,2) = cldfra_conv_k(k)

             enddo
       !   enddo
       !enddo

       !do j=jts,jte
       !   do i=its,ite

                pblh_in(1,2) = pblh_in_k
                ht(1,2) = ht_k

       !   enddo
       !enddo

       !-------------------- GRID 1 0
       call read_input10(kte,t_k,qv_k,qc_k,qi_k,qnc_k,qni_k, &
                         dz8w_k,p_k,p8w_k,pi_k,z_k,z_at_w_k, &
                         tke_pbl_k,u_k,v_k,cldfra_k,cldfra_old_k, &
                         cldfra_old_mp_k,cldfra_conv_k, &
                         ht_k,pblh_in_k)

       !do j=jts,jte
       !   do i=its,ite
             do k=kts,kte

                t_phy(2,k,1) = t_k(k)
                qv(2,k,1) = qv_k(k)
                qc(2,k,1) = qc_k(k)
                qi(2,k,1) = qi_k(k)
                qnc(2,k,1) = qnc_k(k)
                qni(2,k,1) = qni_k(k)          
                dz8w(2,k,1) = dz8w_k(k)
                p(2,k,1) = p_k(k)
                p8w(2,k,1) = p8w_k(k)
                z(2,k,1) = z_k(k)
                z_at_w(2,k,1) = z_at_w_k(k)
                pi_phy(2,k,1) = pi_k(k)
                tke_pbl(2,k,1) = tke_pbl_k(k)
                u_phy(2,k,1) = u_k(k)
                v_phy(2,k,1) = v_k(k)
                cldfra(2,k,1) = cldfra_k(k)
                cldfra_old(2,k,1) = cldfra_old_k(k)
                cldfra_old_mp(2,k,1) = cldfra_old_mp_k(k)
                cldfra_conv(2,k,1) = cldfra_conv_k(k)

             enddo
       !   enddo
       !enddo

       !do j=jts,jte
       !   do i=its,ite

                pblh_in(2,1) = pblh_in_k
                ht(2,1) = ht_k

       !   enddo
       !enddo

       !-------------------- GRID 1 1
       call read_input11(kte,t_k,qv_k,qc_k,qi_k,qnc_k,qni_k, &
                         dz8w_k,p_k,p8w_k,pi_k,z_k,z_at_w_k, &
                         tke_pbl_k,u_k,v_k,cldfra_k,cldfra_old_k, &
                         cldfra_old_mp_k,cldfra_conv_k, &
                         ht_k,pblh_in_k)

       !do j=jts,jte
       !   do i=its,ite
             do k=kts,kte

                t_phy(2,k,2) = t_k(k)
                qv(2,k,2) = qv_k(k)
                qc(2,k,2) = qc_k(k)
                qi(2,k,2) = qi_k(k)
                qnc(2,k,2) = qnc_k(k)
                qni(2,k,2) = qni_k(k)          
                dz8w(2,k,2) = dz8w_k(k)
                p(2,k,2) = p_k(k)
                p8w(2,k,2) = p8w_k(k)
                z(2,k,2) = z_k(k)
                z_at_w(2,k,2) = z_at_w_k(k)
                pi_phy(2,k,2) = pi_k(k)
                tke_pbl(2,k,2) = tke_pbl_k(k)
                u_phy(2,k,2) = u_k(k)
                v_phy(2,k,2) = v_k(k)
                cldfra(2,k,2) = cldfra_k(k)
                cldfra_old(2,k,2) = cldfra_old_k(k)
                cldfra_old_mp(2,k,2) = cldfra_old_mp_k(k)
                cldfra_conv(2,k,2) = cldfra_conv_k(k)

             enddo
       !   enddo
       !enddo

       !do j=jts,jte
       !   do i=its,ite

                pblh_in(2,2) = pblh_in_k
                ht(2,2) = ht_k

       !   enddo
       !enddo

       !------------------------

       call wrf_debug ( 50, ' YLH AFTER VALUE' )
       call wrf_debug(50,'k T QV DZ8W')

       do k=kts,kte

         write( message , * ) k,t_phy(1,k,1),qv(1,k,1),dz8w(1,k,1)
         call wrf_debug ( 50, message )
         write( message , * ) k,t_phy(1,k,2),qv(1,k,2),dz8w(1,k,2)
         call wrf_debug ( 50, message )

       enddo

       call wrf_debug ( 50, ' YLH AFTER VALUE' )
       call wrf_debug(50,'k U V cldfra tke_pbl')

       do k=kts,kte

         write( message , * ) k,u_phy(1,k,1),v_phy(1,k,1),cldfra(1,k,1),tke_pbl(1,k,1)
         call wrf_debug ( 50, message )
         write( message , * ) k,u_phy(1,k,2),v_phy(1,k,2),cldfra(1,k,2),tke_pbl(1,k,2)
         call wrf_debug ( 50, message )
       enddo

       call wrf_debug ( 100, ' YLH AFTER VALUE' )
       call wrf_debug(100,'k HT PBLH')


       write( message , * ) ht(1,1),pblh_in(1,1)
       call wrf_debug ( 100, message )
       write( message , * ) ht(1,2),pblh_in(1,2)
       call wrf_debug ( 100, message )

   endif

endif   

!-----------------------------




#if ( WRF_CHEM == 1 )
  !BSINGH:02/01/2013: Sanity check for Non-MAM simulations
  if(.NOT.cam_mam_aerosols .AND. chem_opt .NE. 0) then
     write(msg,*)'CAMUWSHACU DRIVER - camuworg_driver is valid for only MAM aerosols ', &
          '(chem_opts:',CBMZ_CAM_MAM3_NOAQ,CBMZ_CAM_MAM3_AQ,CBMZ_CAM_MAM7_NOAQ,CBMZ_CAM_MAM7_AQ ,')'
     call wrf_error_fatal( msg )
  endif
#endif


!
! Initialize...
!
  ncol  = 1     !chunk size in WRF is 1 since we loop over all columns in a tile
  ncnst = pcnst !Balwinder.Singh@pnnl.gov
  ztodt = dt

  CALL wrf_debug(100,'in camuw_org_shcu driver')

!
! Map variables to inputs for zm_convr and call it...
! Loop over the points in the tile and treat them each as a CAM chunk.
!

!+ YLH:
!print *, "YLH pcnst, ncnst:", pcnst, ncnst
!print *, "YLH num_moist:",num_moist

!+ YLH: test if moist(1) is qvapor
  !print *,"YLH moist(1,2,1,1)",moist(1,2,1,1)
  !print *,"YLH its,ite,jts,jte",its,ite,jts,jte
  !print *,"YLH kts,kte,kms,kme",kts,kte,kms,kme
  !print *,"YLH ims,ime,jms,jme",ims,ime,jms,jme
  !print *,"YLH qv(1,2,1)",qv(1,2,1)
  !print *,"YLH moist(1,12,1,1)",moist(1,12,1,1)
  !print *,"YLH qv(1,12,1)",qv(1,12,1)

  !do j = jms,jme
  !   do i = ims,ime

  !      print *,"YLH i,j,moist(i,1,j,1)",i,j,moist(i,1,j,1)
  !      print *,"YLH i,j,moist(i,11,j,1)",i,j,moist(i,11,j,1)     
  !      print *,"YLH i,j,moist(i,2,j,1)",i,j,moist(i,2,j,1)
  !      print *,"YLH i,j,moist(i,12,j,1)",i,j,moist(i,12,j,1)
  !      print *,"YLH i,j,moist(i,3,j,1)",i,j,moist(i,3,j,1)
  !      print *,"YLH i,j,moist(i,13,j,1)",i,j,moist(i,13,j,1)        

  !      print *,"YLH i,j,moist(i,13,j,2)",i,j,moist(i,13,j,2)
  !      print *,"YLH i,j,moist(i,13,j,3)",i,j,moist(i,13,j,3)
  !      print *,"YLH i,j,moist(i,13,j,4)",i,j,moist(i,13,j,4)
  !      print *,"YLH i,j,moist(i,13,j,6)",i,j,moist(i,13,j,6)
  !   enddo
  !enddo
!- YLH

  ij_loops : do j = jts,jte
     do i = its,ite
     
!+ YLH
        !print *,"YLH LOOP i,j",i,j
!- YLH

        lchnk = (j-jts)*(ite-its+1) + (i-its+1) !1-D index location from the 2-D tile

        !Flip variables on the layer interfaces
        do k = kts,kte+1
           kflip = kte-k+2

           pint8(1,kflip) = p8w(i,k,j)
           zi8(1,kflip)   = z_at_w(i,k,j) - ht(i,j) ! height above the ground at interfaces
        end do

        !Flip variables on the layer middles
        do k = kts,kte
           kflip = kte-k+1
           if(is_CAMMGMP_used) then
              cld8(1,kflip)    = cldfra_old_mp(i,k,j)
              cldold8(1,kflip) = cldfra_conv(i,k,j)
           else
              cld8(1,kflip)    = cldfra(i,k,j)
              cldold8(1,kflip) = cldfra_old(i,k,j)
           endif
           if (itimestep .eq. 1) then
             cld8(1,kflip)    = 0._r8
             cldold8(1,kflip) = 0._r8
           end if
           cld8(1,kflip) = min(max((cld8(1,kflip) + cldold8(1,kflip)),0._r8),1._r8)

           pdel8(1,kflip)   = p8w(i,k,j) - p8w(i,k+1,j)
           pmid8(1,kflip)   = p(i,k,j)
           qh8(1,kflip)     = max( qv(i,k,j)/(1. + qv(i,k,j)), 1e-30 ) !values of 0 cause a crash in entropy
           if( present(qc) ) then
              qc8(1,kflip)  = qc(i,k,j)/(1. + qv(i,k,j)) !Convert to moist mix ratio
           else
              qc8(1,kflip)  = 0.
           end if
           if( present(qi) ) then
              qi8(1,kflip)  = qi(i,k,j)/(1. + qv(i,k,j)) !Used in convtran, ditto for conversion
           else
              qi8(1,kflip)  = 0.
           end if
           pdeldry8(1,kflip)= pdel8(1,kflip)*(1._r8 - qh8(1,kflip))
           t8(1,kflip)      = t_phy(i,k,j)
           s8(1,kflip)      = cpair*t8(1,kflip) + gravit*(z(i,k,j)-ht(i,j))
           u8(1,kflip)      = u_phy(i,k,j)
           v8(1,kflip)      = v_phy(i,k,j)
           zm8(1,kflip)     = z(i,k,j)-ht(i,j)
        end do

        !BSINGH - TKE at the interfaces
        do k = kts, kte+1
           kflip = kte - k + 2
           
           tke8(1,kflip) = tke_pbl(i,k,j)  !Turbulent kinetic energy            
        end do

        !Flip the tracer array -
        !shift tracer dimension down one to remove "blank" index and
        !convert to wet instead of dry mixing ratios.

!+ YLH: print out moist value for org_adv before org scheme
!        print *, "YLH moist(1,30,1,6) BEFORE org scheme:",moist(1,30,1,6)
!- YLH
        do k = kts,kte
           kflip = kte-k+1

           moist8(1,kflip,1:ncnst) = 0.

           moist8(1,kflip,1) = max(0.0,qv(i,k,j)/(1. + qv(i,k,j)))

           call cnst_get_ind( 'CLDLIQ', m )
           moist8(1,kflip,m) = max(0.0,qc(i,k,j)/(1. + qv(i,k,j)))
           
           !+ YLH:
           !print *, "YLH CLDLIQ m:",m

           call cnst_get_ind( 'CLDICE', m )
           moist8(1,kflip,m) = max(0.0,qi(i,k,j)/(1. + qv(i,k,j)))
           
           !+ YLH:
           !print *, "YLH CLDICE m:",m

           call cnst_get_ind( 'NUMLIQ', m )
           moist8(1,kflip,m) = max(0.0,qnc(i,k,j)/(1. + qv(i,k,j)))
           
           !+ YLH:
           !print *, "YLH NUMLIQ m:",m        

           call cnst_get_ind( 'NUMICE', m )
           moist8(1,kflip,m) = max(0.0,qni(i,k,j)/(1. + qv(i,k,j)))

           !+ YLH:
           !print *, "YLH NUMICE m:",m
           

!+ YLH: get org_adv and org
           call cnst_get_ind( 'ORG_ADV', m )
           !moist(i,k,j,m) = 4e-006
           !moist8(1,kflip,m) = moist(i,k,j,m)
           moist8(1,kflip,m) = org_adv(i,k,j)
           !print *, "YLH1 i,k,j,moist before org scheme:",i,k,j,moist(i,k,j,m)
           !print *, "YLH1 k,moist8 before org scheme:",kflip,moist8(1,kflip,m)

           org8(1) = org(i,j)
           
!- YLH           


#if ( WRF_CHEM == 1 )
           !Following Do-Loop is obtained from chem/module_cam_mam_aerchem_driver.F 
           do l = param_first_scalar, num_chem
              l2 = lptr_chem_to_q(l)
              if ((l2 >= 1) .and. (l2 <= pcnst)) then
                 moist8(1,kflip,l2) = max(0.0,chem(i,k,j,l)*factconv_chem_to_q(l))
              end if
           end do ! l             
#endif           

        end do

        !Some remapping to get arrays to pass into the routine
        pblh(1) = pblh_in(i,j)
        cush(1) = cush_inout(i,j)
     
!+ YLH
        org8(1) = org(i,j)
!+ YLH        
        !+ YLH:
        !print *, "YLH shape of moist8:",shape(moist8)
        !print *, "YLH shape of moist:",shape(moist)
!
! Main guts of the routine...
! This is a bit inefficient because we are flippling the arrays and they
! will then get flipped back again by compute_uwshcu_inv. We are doing
! this to preserve the CAM code as much as possible for maintenance.
!
        call compute_uworg_inv(                        &
             pcols, pver, ncol, ncnst, ztodt,           &
             pint8, zi8, pmid8, zm8, pdel8,             &
             u8, v8, qh8, qc8, qi8,                     &
             t8, s8, moist8,                            &
             tke8, cld8, cldold8, pblh, cush,           &
             cmfmc2, slflx, qtflx,                      &
             flxprec, flxsnow,                          &
             qhtnd, qctnd, qitnd,                       &
             stnd, utnd, vtnd, tnd_tracer,              &
             rprdsh, cmfdqs, precc, snow,               &
             evapcsh, shfrc, iccmr_UW, icwmr_UW,        &
             icimr_UW, cbmf, qc2, rliq2,                &
!+ YLH            
             cnt2, cnb2, fqsatd, lchnk, pdeldry8,       &
             org8,                                      &
             org2thpert,org2qpert,org2wpert,org2cbmf1,  &
             org2cbmf2,org2rkm,org2overlap,             & 
             cbmf1,cbmf2                                )
             !cnt2, cnb2, fqsatd, lchnk, pdeldry8        )
!- YLH

!
! Map output into WRF-dimensioned arrays...
!
        cush_inout(i,j) = cush(1)
 
!+ YLH
        org(i,j) = org8(1)
!+ YLH

!+ YLH: set tnd_tracer placeholder to hold original values of tnd_tracer
       tnd_tracer_ph = tnd_tracer
       !print *,"YLH tnd_tracer(1,30,6):",tnd_tracer(1,30,6)
       !do k = kts,kte
       !    kflip = kte-k+1
       !    print *, "YLH2 k,moist8 AFTER org scheme:",kflip,moist8(1,kflip,m)
       !enddo
!- YLH

!PMA>
       do k = kts,kte
           kflip = kte-k+1
           qc2(1,kflip)=max(0._r8,min(1.e-6_r8,qc2(1,kflip)))
            if( t8(1,kflip) > 268.15_r8 ) then
              dum1 = 0.0_r8
            elseif( t8(1,kflip) < 238.15_r8 ) then
              dum1 = 1.0_r8
            else
              dum1 = ( 268.15_r8 - t8(1,kflip) ) / 30._r8
            endif
           qcten_det(1,kflip) = qc2(1,kflip) * ( 1._r8 - dum1 )
           qiten_det(1,kflip) = qc2(1,kflip) * dum1 
           qcnten_det(1,kflip) = 3._r8 * (qc2(1,kflip)    * ( 1._r8 - dum1 ) ) / (4._r8*3.14159_r8*(10.e-6_r8**3)*997._r8) 
           qinten_det(1,kflip) = 3._r8 * (qc2(1,kflip)    *  dum1 ) / (4._r8*3.14159_r8*(50.e-6_r8**3)*500._r8)  
           qsten_det(1,kflip)      =  qc2(1,kflip) * dum1 * latice                    ! liquid to ice heating
        end do
        do k = kts,kte
           kflip = kte-k+1
           dlf2(i,k,j)         = qc2(1,kflip)
           shfrc3d(i,k,j)      = shfrc(1,kflip)   ! Required by CAM's wet scavenging - Balwinder.Singh@pnnl.gov

           !Add shallow reserved cloud condensate to deep reserved cloud condensate
           ! dlf (kg/kg/s, qc in CAM),  rliq done below
           dlf(i,k,j)          = dlf(i,k,j) + qc2(1,kflip)

           evapcsh_inout(i,k,j)= evapcsh(1,kflip)
           icwmrsh(i,k,j)      = icwmr_uw(1,kflip)

           rprdsh(1,kflip)     = rprdsh(1,kflip) + cmfdqs(1,kflip)
           rprdsh_inout(i,k,j) = rprdsh(1,kflip)
           !Not doing rprdtot for now since not yet used by other CAM routines in WRF

           !Tendencies of winds, potential temperature, and moisture
           !fields treated specifically by UW scheme
           rushten(i,k,j)  = utnd(1,kflip)
           rvshten(i,k,j)  = vtnd(1,kflip)
           rthshten(i,k,j) = (stnd(1,kflip)+qsten_det(1,kflip))/cpair/pi_phy(i,k,j)
           rqvshten(i,k,j) = qhtnd(1,kflip)*(1. + qv(i,k,j))**2
           if( p_qc >= param_first_scalar ) &
                rqcshten(i,k,j) = (qctnd(1,kflip)+qcten_det(1,kflip))*(1. + qv(i,k,j))
           if( p_qi >= param_first_scalar ) &
                rqishten(i,k,j) = (qitnd(1,kflip)+qiten_det(1,kflip))*(1. + qv(i,k,j))

           if( p_qnc >= param_first_scalar ) then
              call cnst_get_ind( 'NUMLIQ', m )
              rqcnshten(i,k,j) = (tnd_tracer(1,kflip,m)+qcnten_det(1,kflip))*(1. + qv(i,k,j))
           endif
           if( p_qni >= param_first_scalar ) then
              call cnst_get_ind( 'NUMICE', m )
              rqinshten(i,k,j) = (tnd_tracer(1,kflip,m)+qinten_det(1,kflip))*(1. + qv(i,k,j))
           endif

!+ YLH
           if( p_org_adv >= param_first_scalar ) then
              call cnst_get_ind( 'ORG_ADV', m )
              rorgshten(i,k,j) = tnd_tracer(1,kflip,m)
              !print *,"---YLH,i,k,kflip,j,p_org_adv",i,k,kflip,j,p_org_adv
              !print *,"---YLH rorgshten(i,k,j)",rorgshten(i,k,j)
              !print *,"---YLH tnd_tracer(1,kflip,m)",tnd_tracer(1,kflip,m)
           endif
!- YLH           

        end do 
           
!+ YLH: Apply ORG_TEND to ORG_ADV to update org3D tracer  

        call cnst_get_ind( 'ORG_ADV', m ) 
        
        !print *,"--- YLH UPDATE ORG_ADV, index:",m

        state_s(:,:)  = 0.0_r8
        ptend_s(:,:)  = 0.0_r8
        ptend_ls      = .FALSE.
        ptend_lq(:)   = .False.
        ptend_lq(m)   = .TRUE.
        ptend_name    = 'org_adv_3d'
        
        !do k = kts,kte
        !   kflip = kte-k+1

        !   print *, "YLH3 k,moist8 before physics_update:",kflip,moist8(1,kflip,m)

        !enddo

        !print *, "----YLH,org_adv(1,10,1) before physics_update:",org_adv(1,10,1)
        !print *, "----YLH,tnd_tracer(1,65,m) before physics_update:",tnd_tracer(1,65,m)
        !print *, "----YLH,tnd_tracer(1,10,m) before physics_update:",tnd_tracer(1,10,m)
        !print *, "----YLH,rorgshten(1,10,1) before physics_update:",rorgshten(1,10,1)

        call physics_update(lchnk,ztodt,moist8,tnd_tracer,state_s,ptend_s,ptend_name,ptend_lq,ptend_ls,pcnst)

        tnd_tracer = tnd_tracer_ph
           
        do k = kts,kte
           kflip = kte-k+1

           !print *,"YLH kflip,tnd_tracer(1,k,iorg):",kflip,tnd_tracer(1,kflip,m)
           !print *,"YLH2 tnd_tracer(1,2,6):",tnd_tracer(1,2,m)
           !print *,"YLH2 tnd_tracer(1,73,6):",tnd_tracer(1,73,m)

           !moist(i,k,j,m) = moist8(1,kflip,m)
           org_adv(i,k,j) = moist8(1,kflip,m)
           !print *, "YLH4 i,k,j,moist AFTER org scheme:",i,k,j,moist(i,k,j,m)
           !print *, "YLH4 k,moist8 AFTER physics_update:",kflip,moist8(1,kflip,m)
                  
        end do !k-loop to kte              

        !print *,"YLH tnd_tracer(1,30,6):",tnd_tracer(1,30,6)
        !print *, "----YLH,org_adv(1,10,1) after physics_update:",org_adv(1,10,1)
        !print *, "----YLH,org(1,1) after physics_update:",org(1,1)
        !print *, "----YLH,tnd_tracer(1,65,m) after physics_update:",tnd_tracer(1,65,m)
        !print *, "----YLH,tnd_tracer(1,10,m) after physics_update:",tnd_tracer(1,10,m)        

!- YLH

!PMA<

       
#if ( WRF_CHEM == 1 )
           !BSINGH - update moist8 by physics update call
           !Update chem array and state constituents
           !populate state_s, ptend_s, ptend_ls with dummy values (zeros) for physics update call
           state_s(:,:)  = 0.0_r8
           ptend_s(:,:)  = 0.0_r8
           ptend_ls      = .FALSE.
           ptend_lq(:)   = .TRUE.
           ptend_lq(1:5) = .FALSE.
           ptend_name    = 'convect_shallow'
           
           call physics_update(lchnk,ztodt,moist8,tnd_tracer,state_s,ptend_s,ptend_name,ptend_lq,ptend_ls,pcnst)
           do k = kts,kte
              kflip = kte-k+1
              do l = param_first_scalar, num_chem
                 l2 = lptr_chem_to_q(l)
                 if ((l2 >= 1) .and. (l2 <= pcnst)) then
                    chem(i,k,j,l) = moist8(1,kflip,l2)/factconv_chem_to_q(l)
                 end if
              end do ! l
           end do !k-loop to kte              
#endif
        


        do k = kts,kte+1
           kflip = kte-k+2

           !Convective fluxes of 'sl' and 'qt' in energy unit
           cmfsl(i,k,j) = slflx(1,kflip)
           cmflq(i,k,j) = qtflx(1,kflip)*latvap
           !BSINGH - Storing CMFMC and CMFMC2 at the interfaces
           cmfmc2_inout(i,k,j) = cmfmc2(1,kflip)
           cmfmc(i,k,j)        = cmfmc(i,k,j) + cmfmc2(1,kflip)
        end do !k-loop to kte+1

        !Calculate fractional occurance of shallow convection
        !~Not doing this since it would require adding time averaging ability across output times

        !Rain rate for shallow convection
        pratesh(i,j) = precc(1)*1e3/dt !~this will need changing for adaptive time steps and cudt

        !Get indices of convection top and bottom based on deep+shallow
        !Note: cnt2 and cnb2 have indices decreasing with height, but
        !      cutop and cubot have indicies increasing with height
        kflip = kte - cutop(i,j) + 1
        cnt = kflip
        if( cnt2(1) < kflip ) cnt = cnt2(1)
        cutop(i,j) = kte - cnt + 1

        kflip = kte - cubot(i,j) + 1
        cnb = kflip
        if( cnb2(1) > kflip ) cnb = cnb2(1)
        cubot(i,j) = kte - cnb + 1

        !Add shallow reserved cloud condensate to deep reserved cloud condensate
        !dlf done above, rliq (m/s)
        rliq2_inout(i,j) = rliq2(1)
        rliq(i,j)        = rliq(i,j) + rliq2(1)

     end do
  end do ij_loops
  
   !------------------
   ! YLH return original value
   !------------------

if (fix_vars_flag .eqv. .true.) then

   if (itimestep .ge. record_timestep) then


       call wrf_debug ( 100, " YLH TEST RETURNED VALUES" )
       call wrf_debug(100,'k T t_ori QV')

       do j=jts,jte
          do k=kts,kte
             do i=its,ite

                t_phy(i,k,j)=t_ori(i,k,j)
                qv(i,k,j)=qv_ori(i,k,j)
                qc(i,k,j)=qc_ori(i,k,j)
                qi(i,k,j)=qi_ori(i,k,j)
                qnc(i,k,j)=qnc_ori(i,k,j)
                qni(i,k,j)=qni_ori(i,k,j)
                dz8w(i,k,j)=dz8w_ori(i,k,j)
                p(i,k,j)=p_ori(i,k,j)
                p8w(i,k,j)=p8w_ori(i,k,j)
                pi_phy(i,k,j)=pi_ori(i,k,j)
                z(i,k,j)=z_ori(i,k,j)
                z_at_w(i,k,j)=z_at_w_ori(i,k,j)
                tke_pbl(i,k,j)=tke_pbl_ori(i,k,j)
                u_phy(i,k,j)=u_ori(i,k,j)
                v_phy(i,k,j)=v_ori(i,k,j)
                cldfra(i,k,j)=cldfra_ori(i,k,j)
                cldfra_old(i,k,j)=cldfra_old_ori(i,k,j)
                cldfra_old_mp(i,k,j)=cldfra_old_mp_ori(i,k,j)
                cldfra_conv(i,k,j)=cldfra_conv_ori(i,k,j)

                WRITE( message , * ) &
                k,t_phy(i,k,j),t_ori(i,k,j),qv(i,k,j)
                CALL wrf_debug ( 100, message )                

             enddo
          enddo
       enddo

       do j=jts,jte
          do i=its,ite

                ht(i,j)=ht_ori(i,j)
                pblh_in(i,j)=pblh_in_ori(i,j)
                
           enddo
       enddo

   endif
endif   


   !------------------
   
  
END SUBROUTINE camuworg_driver

! ****************************************************************************

   ! YLH subroutine to read temperature and moisture input file

   !-------------------- GRID 0 0
   subroutine read_input00(kte,t,qv,qc,qi,qnc,qni, &
                         dz8w,p,p8w,pi,z,z_at_w, &
                         tke_pbl,u,v,cldfra,cldfra_old, &
                         cldfra_old_mp,cldfra_conv, &
                         ht,pblh_in)
                     
      implicit none
      integer, intent(in)  :: kte
      real, dimension(kte) :: t,qv,qc,qi,qnc,qni, &
                              dz8w,p,p8w,pi,z,z_at_w, &
                              tke_pbl,u,v,cldfra,cldfra_old, &
                              cldfra_old_mp,cldfra_conv
      real ::  ht,pblh_in
      integer :: k
      logical :: end_of_file

      ! YLH 3-dim vars
      open(unit=10,file='uw_profiles_3dim_00',form='formatted',status='old')
      rewind(10)

      end_of_file = .false.
      k = 1

      do while(.not. end_of_file)

          if (k .gt. kte) then

              end_of_file = .true.

          else if (k .le. kte) then

              read(10,*) t(k),qv(k),qc(k),qi(k),qnc(k),qni(k), &
                         dz8w(k),p(k),p8w(k),pi(k),z(k),z_at_w(k), &
                         tke_pbl(k),u(k),v(k),cldfra(k),cldfra_old(k), &
                         cldfra_old_mp(k),cldfra_conv(k) 
              !print*, t(k), qv(k)

              k = k+1


          end if

      enddo

      close(unit=10,status='keep')

      ! YLH 2-dim vars

      open(unit=10,file='uw_profiles_2dim_00',form='formatted',status='old')

      rewind(10)

      end_of_file = .false.

      read(10,*) ht,pblh_in

      end_of_file = .true.

      close(unit=10,status='keep')

   end subroutine read_input00

   
   !-------------------- GRID 0 1
   subroutine read_input01(kte,t,qv,qc,qi,qnc,qni, &
                         dz8w,p,p8w,pi,z,z_at_w, &
                         tke_pbl,u,v,cldfra,cldfra_old, &
                         cldfra_old_mp,cldfra_conv, &
                         ht,pblh_in)
                     
      implicit none
      integer, intent(in)  :: kte
      real, dimension(kte) :: t,qv,qc,qi,qnc,qni, &
                              dz8w,p,p8w,pi,z,z_at_w, &
                              tke_pbl,u,v,cldfra,cldfra_old, &
                              cldfra_old_mp,cldfra_conv
      real ::  ht,pblh_in
      integer :: k
      logical :: end_of_file

      ! YLH 3-dim vars
      open(unit=10,file='uw_profiles_3dim_01',form='formatted',status='old')
      rewind(10)

      end_of_file = .false.
      k = 1

      do while(.not. end_of_file)

          if (k .gt. kte) then

              end_of_file = .true.

          else if (k .le. kte) then

              read(10,*) t(k),qv(k),qc(k),qi(k),qnc(k),qni(k), &
                         dz8w(k),p(k),p8w(k),pi(k),z(k),z_at_w(k), &
                         tke_pbl(k),u(k),v(k),cldfra(k),cldfra_old(k), &
                         cldfra_old_mp(k),cldfra_conv(k) 
              !print*, t(k), qv(k)

              k = k+1


          end if

      enddo

      close(unit=10,status='keep')

      ! YLH 2-dim vars

      open(unit=10,file='uw_profiles_2dim_01',form='formatted',status='old')

      rewind(10)

      end_of_file = .false.

      read(10,*) ht,pblh_in

      end_of_file = .true.

      close(unit=10,status='keep')

   end subroutine read_input01

   
   !-------------------- GRID 1 0
      subroutine read_input10(kte,t,qv,qc,qi,qnc,qni, &
                         dz8w,p,p8w,pi,z,z_at_w, &
                         tke_pbl,u,v,cldfra,cldfra_old, &
                         cldfra_old_mp,cldfra_conv, &
                         ht,pblh_in)
                     
      implicit none
      integer, intent(in)  :: kte
      real, dimension(kte) :: t,qv,qc,qi,qnc,qni, &
                              dz8w,p,p8w,pi,z,z_at_w, &
                              tke_pbl,u,v,cldfra,cldfra_old, &
                              cldfra_old_mp,cldfra_conv
      real ::  ht,pblh_in
      integer :: k
      logical :: end_of_file

      ! YLH 3-dim vars
      open(unit=10,file='uw_profiles_3dim_10',form='formatted',status='old')
      rewind(10)

      end_of_file = .false.
      k = 1

      do while(.not. end_of_file)

          if (k .gt. kte) then

              end_of_file = .true.

          else if (k .le. kte) then

              read(10,*) t(k),qv(k),qc(k),qi(k),qnc(k),qni(k), &
                         dz8w(k),p(k),p8w(k),pi(k),z(k),z_at_w(k), &
                         tke_pbl(k),u(k),v(k),cldfra(k),cldfra_old(k), &
                         cldfra_old_mp(k),cldfra_conv(k) 
              !print*, t(k), qv(k)

              k = k+1


          end if

      enddo

      close(unit=10,status='keep')

      ! YLH 2-dim vars

      open(unit=10,file='uw_profiles_2dim_10',form='formatted',status='old')

      rewind(10)

      end_of_file = .false.

      read(10,*) ht,pblh_in

      end_of_file = .true.

      close(unit=10,status='keep')

   end subroutine read_input10


   !-------------------- GRID 1 1
      subroutine read_input11(kte,t,qv,qc,qi,qnc,qni, &
                         dz8w,p,p8w,pi,z,z_at_w, &
                         tke_pbl,u,v,cldfra,cldfra_old, &
                         cldfra_old_mp,cldfra_conv, &
                         ht,pblh_in)
                     
      implicit none
      integer, intent(in)  :: kte
      real, dimension(kte) :: t,qv,qc,qi,qnc,qni, &
                              dz8w,p,p8w,pi,z,z_at_w, &
                              tke_pbl,u,v,cldfra,cldfra_old, &
                              cldfra_old_mp,cldfra_conv
      real ::  ht,pblh_in
      integer :: k
      logical :: end_of_file

      ! YLH 3-dim vars
      open(unit=10,file='uw_profiles_3dim_11',form='formatted',status='old')
      rewind(10)

      end_of_file = .false.
      k = 1

      do while(.not. end_of_file)

          if (k .gt. kte) then

              end_of_file = .true.

          else if (k .le. kte) then

              read(10,*) t(k),qv(k),qc(k),qi(k),qnc(k),qni(k), &
                         dz8w(k),p(k),p8w(k),pi(k),z(k),z_at_w(k), &
                         tke_pbl(k),u(k),v(k),cldfra(k),cldfra_old(k), &
                         cldfra_old_mp(k),cldfra_conv(k) 
              !print*, t(k), qv(k)

              k = k+1


          end if

      enddo

      close(unit=10,status='keep')

      ! YLH 2-dim vars

      open(unit=10,file='uw_profiles_2dim_11',form='formatted',status='old')

      rewind(10)

      end_of_file = .false.

      read(10,*) ht,pblh_in

      end_of_file = .true.

      close(unit=10,status='keep')

   end subroutine read_input11

! ****************************************************************************

END MODULE module_shcu_camuworg_driver
