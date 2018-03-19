MODULE SAO_OMIL2_allocation_module

  USE  SAO_OMIL2_ReadLib_Basic_module, ONLY: r8, r4, i4, i2
  IMPLICIT NONE

  ! -----------------------
  ! Gridded Array Variables
  ! -----------------------
  REAL    (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: gcol_reg, gcol_cor, gcol_err, &
                                                     gcol_aer, gcol_alb, gcol_amf, &
                                                     gslt_reg, gslt_cor, gslt_err, &
                                                     gslt_aer, gcol_rms, gsrf_alt
  REAL    (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: gcld_cfr, gcld_ctp
  INTEGER (KIND=i2), DIMENSION (:,:), ALLOCATABLE :: good_idx, gcol_qfl
  INTEGER (KIND=i4), DIMENSION (:,:), ALLOCATABLE :: gcol_num

  ! ------------------
  ! OMI Data Variables
  ! ------------------
  INTEGER (KIND=i2), DIMENSION (:,:), ALLOCATABLE :: qflag, xtqfe, xtqf, corflag, &
                                                     amfflag, srf_alt
  REAL    (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: lat, lon, sza, cld_cfr, cld_ctp, &
                                                     vza, amf, amf_epr
  REAL    (KIND=r4), DIMENSION (:,:), ALLOCATABLE :: corlat, corlon
  REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: col_reg, col_cor, col_err, &
                                                     slt_reg, slt_cor, slt_err, &
                                                     col_rms, srf_alb

END MODULE SAO_OMIL2_allocation_module
