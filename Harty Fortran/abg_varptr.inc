C*****************************************************************
C*---------------------------MDI---------------------------------*
C*---------------------AUTHOR : pkumar---------------------------*
C*-----------------------DATE : 04/19/99-------------------------*
C*****************************************************************
C*                                                               *
C* who          when          what                               *
C*----          ----          ----                               *
C* pkumar       04/01/99      Routine written                    *
C* joost        04/05/04      Extended for A/Aircraft            *
C*****************************************************************
C*                                                               *
C* This routine is an include file for tire models that defines  *
C* pointer locations for tire outputs into VARINF array according*
C* to TYDEX standard                                             *
C*                                                               *
C* Input                                                         *
C* -----                                                         *
C*      None                                                     *
C*                                                               *
C* Output                                                        *
C* ------                                                        *
C*      None                                                     *
C*                                                               *
C*****************************************************************

C--- Tire forces in W-frame (TYDEX standard)

      integer fx_w_ptr, fy_w_ptr, fz_w_ptr
      integer mx_w_ptr, my_w_ptr, mz_w_ptr

      parameter (fx_w_ptr = 1, fy_w_ptr = 2, fz_w_ptr = 3)
      parameter (mx_w_ptr = 4, my_w_ptr = 5, mz_w_ptr = 6)


C--- Tire slip an camber (with relaxation effects)- In ISO Frame
C---  slipi_i_ptr - Lateral slip (rad)
C---  slipx_i_ptr - Longitudinal slip (%)
C---  camb_ptr    - Camber (rad)

      integer    slipi_ptr,  slipx_ptr, camb_ptr
      parameter (slipi_ptr = 7,  slipx_ptr = 8, camb_ptr = 9)


C--- Tire forces in SAE-frame

      integer fx_sae_ptr, fy_sae_ptr, fz_sae_ptr
      integer mx_sae_ptr, my_sae_ptr, mz_sae_ptr

      parameter (fx_sae_ptr = 26, fy_sae_ptr = 27, fz_sae_ptr = 28)
      parameter (mx_sae_ptr = 29, my_sae_ptr = 30, mz_sae_ptr = 31)


C--- Tire forces at Hub center 

      integer fx_hub_ptr, fy_hub_ptr, fz_hub_ptr
      integer mx_hub_ptr, my_hub_ptr, mz_hub_ptr

      parameter (fx_hub_ptr = 32, fy_hub_ptr = 33, fz_hub_ptr = 34)
      parameter (mx_hub_ptr = 35, my_hub_ptr = 36, mz_hub_ptr = 37)


C--- Tire forces in ISO-frame

      integer fx_iso_ptr, fy_iso_ptr, fz_iso_ptr
      integer mx_iso_ptr, my_iso_ptr, mz_iso_ptr

      parameter (fx_iso_ptr = 38, fy_iso_ptr = 39, fz_iso_ptr = 40)
      parameter (mx_iso_ptr = 41, my_iso_ptr = 42, mz_iso_ptr = 43)


C--- Tire kinematics
C---  effpen_ptr - Effective penetration
C---  vcpvrt_ptr - Tire closing velocity 
C---  vcplon_ptr - Longitudinal slip velocity at CP
C---  vcplat_ptr - Lateral slip velocity at CP
C---  vlon_ptr   - Longitudinal velocity at Wheel Center
C---  radius_ptr - Effective rolling radius
C---  omegar_ptr - Tire rotational velocity

      integer    effpen_ptr, vcpvrt_ptr, vcplon_ptr
      integer    vcplat_ptr, vlon_ptr, radius_ptr
      integer    omegar_ptr

      parameter (effpen_ptr = 44, vcpvrt_ptr = 45, vcplon_ptr = 46)
      parameter (vcplat_ptr = 47, vlon_ptr = 48, radius_ptr = 49)
      parameter (omegar_ptr = 50)


C--- Tire slip definition (without relaxation effects)
C---  slipi_d_ptr - lateral slip (rad)
C---  slipx_d_ptr - longitudinal slip (%)

      integer    slipi_d_ptr, slipx_d_ptr
      parameter (slipi_d_ptr = 51, slipx_d_ptr = 52)


C--- Derivatives of state variables
C---  dudt_ptr - Longitudinal direction
C---  dvdt_ptr - Lateral direction

      integer    dudt_ptr, dvdt_ptr
      parameter (dudt_ptr = 53, dvdt_ptr = 54)


C--- Friction coefficients
C---  muxtyr_ptr - Longitudinal direction
C---  muytyr_ptr - Lateral direction

      integer    muxtyr_ptr, muytyr_ptr
      parameter (muxtyr_ptr = 55, muytyr_ptr = 56)


C--- Tire characteristics
C---  pt_ptr     - pneumatic trail
C---  mzr_ptr    - residual torque
C---  s_ptr      - moment arm of Fx
C---  sigkp0_ptr - longitudinal relaxation length
C---  sigal0_ptr - lateral relaxation length
C---  mgyr_ptr   - gyroscopic moment
C---  svykap_ptr - kappa induced side force
C---  svx_ptr    - 
C---  svy_ptr    - 

      integer pt_ptr, mzr_ptr, s_ptr, sigkp0_ptr, sigal0_ptr
      integer mgyr_ptr, svykap_ptr, svx_ptr, svy_ptr

      parameter (pt_ptr = 57, mzr_ptr = 58, s_ptr = 59)
      parameter (sigkp0_ptr = 60, sigal0_ptr = 61)
      parameter (mgyr_ptr = 62, svykap_ptr = 63, svx_ptr = 64)
      parameter (svy_ptr = 65)


C--- Contact point coordinates

      integer    rcp1_ptr, rcp2_ptr, rcp3_ptr
      parameter (rcp1_ptr = 66, rcp2_ptr = 67, rcp3_ptr = 68)


C--- Road normal

      integer    rnorm1_ptr, rnorm2_ptr, rnorm3_ptr
      parameter (rnorm1_ptr = 69, rnorm2_ptr = 70, rnorm3_ptr = 71)


C--- Road surface parameters

      integer surfac1_ptr, surfac2_ptr, surfac3_ptr
      parameter (surfac1_ptr = 72, surfac2_ptr = 73, surfac3_ptr = 74)

C--- SWIFT Contact Parameters:

      integer dis_tra_ptr
      parameter ( dis_tra_ptr = 75 )

      integer eff_pla_hei_ptr, eff_pla_ang_ptr, eff_pla_cur_ptr
      parameter ( eff_pla_hei_ptr = 76 )
      parameter ( eff_pla_ang_ptr = 77 )
      parameter ( eff_pla_cur_ptr = 78 )

      integer con_len_ptr
      parameter ( con_len_ptr = 79 )

C--- Additional A/Air parameters

      integer flength_ptr
      parameter (flength_ptr = 80)

C--- Loaded Radius Pointer
      integer lrad_ptr
      parameter (lrad_ptr = 95)

