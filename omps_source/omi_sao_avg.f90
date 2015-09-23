PROGRAM omi_sao_avg

  ! =====================================================
  ! Program for tessellation-based gridding of OMI data
  ! =====================================================

  USE SAO_OMIL2_ReadLib_basic_module
  IMPLICIT NONE

  ! ---------------------------
  ! Various CHARACTER variables
  ! ---------------------------
  CHARACTER (LEN=maxchlen) :: input_control_file
  CHARACTER (LEN=maxchlen) :: listfile, outfile, pge_esdt, output_format, swathname

  ! ------------------------------
  ! Collection of LOGICALs
  ! ------------------------------
  LOGICAL :: yn_use_rbszoom  ! Include rebinned spatial zoom granules
  LOGICAL :: yn_norm_output  ! Normalization of output
  LOGICAL :: yn_remove_bg    ! Remove background values (HCHO only)
  LOGICAL :: yn_gpix_weight  ! Ground-pixel size weighted gridding
  LOGICAL :: yn_ucert_weight ! Uncertainty-weighted gridding
  LOGICAL :: yn_amf_geo      ! Geometric air mass factor

  ! -----------------------------------------------------
  ! Normazliation value for uncertainty-weighted gridding
  ! -----------------------------------------------------
  REAL (KIND=r8) :: errwght

  ! -----------------------------------
  ! The Grid and all related quantities
  ! -----------------------------------
  INTEGER (KIND=i4) :: nlongr, nlatgr
  REAL    (KIND=r4) :: dlongr, dlatgr
  REAL    (KIND=r4) :: latmin, latmax, lonmin, lonmax

  ! ------------------------------------------------
  ! Maxima for cloud fraction and solar zenith angle
  ! ------------------------------------------------
  REAL (KIND=r4) :: cld_frc_min, cld_frc_max, szamax

  ! ---------------------------------------
  ! Maximum and minimum for xtrack position
  ! ---------------------------------------
  INTEGER (KIND=i2) :: xtrackmin, xtrackmax

  ! -----------------------------------------------
  ! Maximum quality flag value to include (0, 1, 2)
  ! -----------------------------------------------
  INTEGER (KIND=i2) :: qflg_max


  ! --------------------------------
  ! Make PGE write STD/IO unbuffered
  ! --------------------------------
  CALL unbufferSTDout()

  ! ----------------------------------------------
  ! Read name of process control file from STD I/O
  ! ----------------------------------------------
  READ (*,'(A)') input_control_file

  ! -------------------------------------------------
  ! Read averaging control file with process settings
  ! -------------------------------------------------
  CALL omi_avg_read_input (                                           &
     TRIM(ADJUSTL(input_control_file)),                               &
     pge_esdt, listfile, outfile, output_format, swathname,           &
     yn_norm_output, yn_gpix_weight, yn_ucert_weight, yn_use_rbszoom, &
     yn_amf_geo, yn_remove_bg,                                        &
     qflg_max, szamax, cld_frc_min, cld_frc_max, errwght,             &
     lonmin, lonmax, dlongr, latmin, latmax, dlatgr,                  &
     xtrackmin, xtrackmax                                             )
  
  ! -----------------------
  ! Set up the lon/lat grid
  ! -----------------------------------------------------------------------------------
  ! Note that we add one entry to the number of grid points. This will provide an upper
  ! boundary at the dateline, and it will save us some headache when figuring out the 
  ! upper and lower grid cell indices that cover the satellite pixel.
  ! -----------------------------------------------------------------------------------
  nlongr = NINT( (lonmax-lonmin) / dlongr ) 
  nlatgr = NINT( (latmax-latmin) / dlatgr )

  CALL gridding_process (                                            &
       TRIM(ADJUSTL(pge_esdt)), TRIM(ADJUSTL(listfile)),             &
       TRIM(ADJUSTL(outfile)), nlongr, nlatgr, dlongr, dlatgr,       &
       lonmin, lonmax, latmin, latmax, szamax, cld_frc_min,          &
       cld_frc_max, yn_gpix_weight, yn_ucert_weight, errwght,        &
       yn_use_rbszoom, yn_amf_geo, qflg_max,                         &
       yn_remove_bg, yn_norm_output, TRIM(ADJUSTL(output_format)),   &
       TRIM(ADJUSTL(swathname)), xtrackmin, xtrackmax                )

  STOP 
END PROGRAM omi_sao_avg


SUBROUTINE gridding_process (                                        &
       pge_esdt, listfile, outfile, nlongr, nlatgr, dlongr, dlatgr,  &
       lonmin, lonmax, latmin, latmax, szamax, cld_frc_min,          &
       cld_frc_max, yn_gpix_weight, yn_ucert_weight, errwght,        &
       yn_use_rbszoom, yn_amf_geo, qflg_max,                         &
       yn_remove_bg, yn_norm_output, output_format, swathname,       &
       xtrackmin, xtrackmax                                          )

  USE SAO_OMIL2_ReadLib_basic_module
  USE SAO_OMIL2_ReadLib_he5_module
  USE SAO_OMIL2_CreateHDF5_module
  USE SAO_OMIL2_allocation_module, ONLY: &
       gcol_reg, gcol_cor, gcol_err, gcol_aer, gcol_alb, good_idx, &
       gcol_qfl, gcol_num, gcol_amf, gslt_reg, gslt_cor, gslt_err, &
       gslt_aer, gcol_rms, gsrf_alt, gcld_cfr, gcld_ctp

  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  CHARACTER (LEN=*), INTENT (IN) :: pge_esdt, listfile, outfile, output_format, swathname
  LOGICAL,           INTENT (IN) :: yn_gpix_weight, yn_ucert_weight, yn_use_rbszoom
  LOGICAL,           INTENT (IN) :: yn_norm_output, yn_remove_bg, yn_amf_geo
  INTEGER (KIND=i4), INTENT (IN) :: nlongr, nlatgr
  REAL    (KIND=r8), INTENT (IN) :: errwght
  REAL    (KIND=r4), INTENT (IN) :: dlongr, dlatgr
  REAL    (KIND=r4), INTENT (IN) :: lonmin, lonmax, latmin, latmax, szamax, cld_frc_min, cld_frc_max
  INTEGER (KIND=i2), INTENT (IN) :: qflg_max, xtrackmin, xtrackmax


  ! ------------------------------
  ! Local Variables and Parameters
  ! ------------------------------
  INTEGER (KIND=i4) :: i, ilon, ilat, ios

  ! -----------
  ! OMI L2 data
  ! -----------
  CHARACTER (LEN=maxchlen) :: l2_fname, l2_swath
  INTEGER   (KIND=i4)      :: he5stat, l2_swathfile_id, l2_swath_id
  INTEGER   (KIND=i4)      :: nTimes, nXtrack
  REAL      (KIND=r4)      :: max_area

  ! --------------
  ! Grid variables
  ! --------------
  INTEGER (KIND=i4), PARAMETER :: idx_col = 1, idx_cor = 2, idx_err = 3, idx_aer = 4,  &
                                  idx_col_slt = 5, idx_cor_slt = 6, idx_err_slt = 7,   &
                                  idx_num = 8, idx_amf = 9, idx_flg = 10, idx_alb = 11,&
                                  idx_rms = 12
  REAL    (KIND=r4)                                :: low_no_good, frac, deno, frac_slt
  REAL    (KIND=r4), DIMENSION (1:nlongr)          :: grid_lon
  REAL    (KIND=r4), DIMENSION (1:nlatgr)          :: grid_lat
  REAL    (KIND=r4), DIMENSION (idx_rms)           :: good_norm

  ! -------------------------------------------------------------------------
  ! Allocatable arrays. For whatever reason, declaring the following arrays
  ! with DIMENSION(1:nlongr,1:nlatgr) leads to segmentation faults (Intel/
  ! ifort 32bit compilation) at run-time. Hence the, admittedly less elegant,
  ! dynamic memory allocation.
  !
  ! NOTE: THIS HAS NOW BEEN MOVED TO A SEPARATE MODULE
  !
  ! Also note that "WHERE" statements have been mostly replaced by loops in
  ! order to reduce memory strain. In this manner the code runs even for 
  ! fairly large grid sizes () on a laptop of average memory.
  ! -------------------------------------------------------------------------
  INTEGER (KIND=i4) :: estat


  ALLOCATE ( gcol_alb(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_alb'
  ALLOCATE ( gcol_reg(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_reg'
  ALLOCATE ( gcol_cor(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_cor'
  ALLOCATE ( gcol_err(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_err'
  ALLOCATE ( good_idx(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'good_idx'
  ALLOCATE ( gcol_aer(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_aer'
  ALLOCATE ( gcol_qfl(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_qfl'
  ALLOCATE ( gcol_num(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_num'
  ALLOCATE ( gcol_amf(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_amf'
  ALLOCATE ( gslt_reg(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gslt_reg'
  ALLOCATE ( gslt_cor(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gslt_cor'
  ALLOCATE ( gslt_err(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gslt_err'
  ALLOCATE ( gslt_aer(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gslt_aer'
  ALLOCATE ( gcol_rms(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcol_rms'
  ALLOCATE ( gsrf_alt(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gsrf_alt'
  ALLOCATE ( gcld_cfr(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcld_cfr'
  ALLOCATE ( gcld_ctp(1:nlongr,1:nlatgr), STAT=estat ) ; IF ( estat /= 0 ) STOP 'gcld_ctp'

  ! ------------------------
  ! Initialize output arrays
  ! ------------------------
  good_idx = 0_i2   ;  gcol_qfl =   0_i2
  gcol_reg = 0.0_r4 ;  gcol_err = 0.0_r4
  gcol_cor = 0.0_r4 ;  gcol_aer = 0.0_r4
  gcol_alb = 0.0_r4 ;  gcol_num =   0_i4
  gcol_amf = 0.0_r4 ;  gslt_reg = 0.0_r4
  gslt_cor = 0.0_r4 ;  gslt_err = 0.0_r4
  gslt_aer = 0.0_r4 ;  gcol_rms = 0.0_r4
  gsrf_alt = 0.0_r4 ;  gcld_cfr = 0.0_r4
  gcld_ctp = 0.0_r4

  ! ------------------------------------------------------------------
  ! Calculate the maximally possible area a tessellated pixel may have
  ! ------------------------------------------------------------------
  max_area = (lonmax-lonmin)*(latmax-latmin)

  ! ------------------------------------------
  ! Compose longitude and latitude grid values
  ! ------------------------------------------
  grid_lon(1:nlongr) = lonmin + (/ (REAL(i,KIND=r4)*dlongr, i = 0, nlongr-1) /)
  grid_lat(1:nlatgr) = latmin + (/ (REAL(i,KIND=r4)*dlatgr, i = 0, nlatgr-1) /)

  ! -------------------------------
  ! Open list with input file names
  ! -------------------------------
  OPEN ( UNIT=listunit, FILE=TRIM(ADJUSTL(listfile)), STATUS='OLD', ACTION='READ', IOSTAT=ios )
  IF ( ios /= 0 ) THEN
     WRITE (*, '(3A)') 'ERROR: Failed to open <', TRIM(ADJUSTL(listfile)), '> for reading'
     STOP 1
  END IF
 
  ios = 0
  datfiles: DO WHILE ( ios == 0 )

     READ (UNIT=listunit, FMT='(A)', IOSTAT=ios) l2_fname

     ! -------------------------------------------------------
     ! Check READ status and L2 FILENAME; exit/skip if no good
     ! -------------------------------------------------------
     IF ( ios                <  0 ) EXIT datfiles
     IF ( LEN_TRIM(l2_fname) == 0 ) CYCLE

     WRITE (*,'(A)') TRIM(ADJUSTL(l2_fname))

     ! ----------------
     ! Read OMI L2 data
     ! ----------------
     l2_swath=''
     CALL saopge_l2_open_attach ( &
          he5stat, TRIM(ADJUSTL(l2_fname)), l2_swathfile_id, l2_swath_id, l2_swath )
     IF ( he5stat /= 0 .OR. l2_swathfile_id < 0 ) CYCLE

     CALL saopge_l2_read_dimensions ( he5stat, l2_swath_id, nTimes_k=nTimes, nXtrack_k=nXtrack )
     IF ( he5stat /= 0 .OR. l2_swathfile_id < 0 ) THEN
        he5stat = HE5_SWdetach ( l2_swath_id )
        he5stat = HE5_SWclose  ( l2_swathfile_id )
        CYCLE
     END IF

     CALL datafile_loop (                                                         &
          TRIM(ADJUSTL(pge_esdt)),                                                &
          nXtrack, nTimes, l2_swathfile_id, l2_swath_id,                          &
          nlongr, nlatgr, dlongr, dlatgr, grid_lon(1:nlongr), grid_lat(1:nlatgr), &
          latmin, latmax, lonmin, lonmax, szamax, cld_frc_min, cld_frc_max,       &
          max_area, yn_gpix_weight, yn_ucert_weight, errwght, yn_use_rbszoom,     &
          yn_amf_geo, yn_remove_bg, qflg_max, xtrackmin, xtrackmax  )

     ! --------------------------------
     ! Detach from swath and close file
     ! --------------------------------
     he5stat = HE5_SWdetach ( l2_swath_id )
     he5stat = HE5_SWclose  ( l2_swathfile_id )

  END DO datfiles
  CLOSE ( listunit )

  ! ----------------------------------------------------------------------------
  ! The final normalization to the total area covered in a gived grid cell pixel.
  ! -----------------------------------------------------------------------------
  DO ilon = 1, nlongr
     DO ilat = 1, nlatgr
        IF ( gcol_aer(ilon,ilat) /= 0.0_r4 .AND. &
             gcol_err(ilon,ilat) >  0.0_r4 .AND. &
             gslt_aer(ilon,ilat) /= 0.0_r4 .AND. &
             gslt_err(ilon,ilat) >  0.0_r4 .AND. &
             gcol_num(ilon,ilat) >  0_i4          ) THEN
           frac                = 1.0_r4 / gcol_aer(ilon,ilat) 
           frac_slt            = 1.0_r4 / gslt_aer(ilon,ilat)
           good_idx(ilon,ilat) = 1_i2
           gcol_qfl(ilon,ilat) = 1_i2
           gcol_reg(ilon,ilat) = gcol_reg(ilon,ilat) * frac
           gcol_cor(ilon,ilat) = gcol_cor(ilon,ilat) * frac
           gslt_reg(ilon,ilat) = gslt_reg(ilon,ilat) * frac_slt
           gslt_cor(ilon,ilat) = gslt_cor(ilon,ilat) * frac_slt
           gcol_alb(ilon,ilat) = gcol_alb(ilon,ilat) * frac
           gcol_amf(ilon,ilat) = gcol_amf(ilon,ilat) * frac
           gsrf_alt(ilon,ilat) = gsrf_alt(ilon,ilat) * frac
           gcld_cfr(ilon,ilat) = gcld_cfr(ilon,ilat) * frac
           gcld_ctp(ilon,ilat) = gcld_ctp(ilon,ilat) * frac
           gcol_rms(ilon,ilat) = gcol_rms(ilon,ilat) * frac
           ! ------------------------------------------------------------------
           ! GCOL_ERR holds the sum of the inverse squares of the unweighted
           ! fitting uncertainties. The final uncertainties are the SQRT of
           ! the inverse of that, divided by the SQRT of the number of samples.
           ! See for example <http://en.wikipedia.org/wiki/Weighted_mean>.
           ! Multiplication with ERRWGHT comes in because we have normalized
           ! the colunm errors to it in order to avoid numerical overflows
           ! or underflows when computing the squares of them.
           ! ------------------------------------------------------------------
           gcol_err(ilon,ilat) = REAL ( errwght, KIND=r4) / &
                ( SQRT(gcol_err(ilon,ilat)*REAL(gcol_num(ilon,ilat), KIND=r4)) )
           gslt_err(ilon,ilat) = REAL ( errwght, KIND=r4) / &
                ( SQRT(gslt_err(ilon,ilat)*REAL(gcol_num(ilon,ilat), KIND=r4)) )
        ELSE
           good_idx(ilon,ilat) = 0_i2
           gcol_qfl(ilon,ilat) = 0_i2
           gcol_num(ilon,ilat) = 0_i4
           gcol_reg(ilon,ilat) = r4missval
           gcol_cor(ilon,ilat) = r4missval
           gcol_err(ilon,ilat) = r4missval
           gslt_reg(ilon,ilat) = r4missval
           gslt_cor(ilon,ilat) = r4missval
           gslt_err(ilon,ilat) = r4missval
           gcol_alb(ilon,ilat) = r4missval
           gcol_amf(ilon,ilat) = r4missval
           gsrf_alt(ilon,ilat) = r4missval
           gcld_cfr(ilon,ilat) = r4missval
           gcld_ctp(ilon,ilat) = r4missval
           gcol_rms(ilon,ilat) = r4missval
        END IF
     END DO
  END DO
  ! ----------------------------
  ! Normalize output quantities?
  ! ----------------------------
  good_norm = 1.0_r4
  IF ( yn_norm_output ) THEN
     low_no_good = r4missval
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gcol_reg(1:nlongr,1:nlatgr), good_norm(idx_col)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gcol_cor(1:nlongr,1:nlatgr), good_norm(idx_cor)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gcol_err(1:nlongr,1:nlatgr), good_norm(idx_err)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gslt_reg(1:nlongr,1:nlatgr), good_norm(idx_col_slt)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gslt_cor(1:nlongr,1:nlatgr), good_norm(idx_cor_slt)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gslt_err(1:nlongr,1:nlatgr), good_norm(idx_err_slt)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gcol_aer(1:nlongr,1:nlatgr), good_norm(idx_aer)   )
     CALL normalize_gridded_values (                     &
          nlongr, nlatgr, low_no_good,                   &
          gcol_rms(1:nlongr,1:nlatgr), good_norm(idx_rms)   )

  END IF

  ! --------------------------
  ! Write output to file
  ! --------------------------
  IF ( INDEX(TRIM(ADJUSTL(output_format)), 'ascii') > 0 ) THEN
     CALL gridded_ascii_output (                                    &
          TRIM(ADJUSTL(outfile)),                                   &
          nlongr, nlatgr, dlongr, dlatgr, good_norm(1:4),           &
          grid_lon(1:nlongr), grid_lat(1:nlatgr),                   &
          gcol_reg(1:nlongr,1:nlatgr), gcol_err(1:nlongr,1:nlatgr), &
          gcol_cor(1:nlongr,1:nlatgr), gcol_aer(1:nlongr,1:nlatgr)   )
  END IF

  IF ( INDEX(TRIM(ADJUSTL(output_format)), 'he5') > 0 ) THEN

     swath_name = TRIM(ADJUSTL(swathname))

     CALL create_he5_file (                               &
          TRIM(ADJUSTL(outfile)), nlongr, nlatgr,         &
          grid_lon(1:nlongr), grid_lat(1:nlatgr), dlongr, &
          dlatgr, cld_frc_min, cld_frc_max, qflg_max,     &
          xtrackmin, xtrackmax                            )

     CALL write_he5_data_r4 ( 'SurfaceAlbedo',     nlongr, nlatgr,       &
          gcol_alb(1:nlongr,1:nlatgr), good_norm(idx_alb)           )
     CALL write_he5_data_r4 ( 'AMF',          nlongr, nlatgr,            &
          gcol_amf(1:nlongr,1:nlatgr), good_norm(idx_amf)           )
     CALL write_he5_data_r4 ( 'GridArea',          nlongr, nlatgr,       &
          gcol_aer(1:nlongr,1:nlatgr), good_norm(idx_aer)           )
     CALL write_he5_data_i2 ( 'QualityFlag',       nlongr, nlatgr,       &
          gcol_qfl(1:nlongr,1:nlatgr), good_norm(idx_flg)           )
     CALL write_he5_data_i4 ( 'NumberOfSamples',   nlongr, nlatgr,       &
          gcol_num(1:nlongr,1:nlatgr), good_norm(idx_num)           )
     CALL write_he5_data_r4 ( 'SurfaceAltitude', nlongr, nlatgr,         &
          gsrf_alt(1:nlongr,1:nlatgr), good_norm(idx_err)           )
     CALL write_he5_data_r4 ( 'CloudFraction', nlongr, nlatgr,           &
          gcld_cfr(1:nlongr,1:nlatgr), good_norm(idx_err)           )
     CALL write_he5_data_r4 ( 'CloudPressure', nlongr, nlatgr,           &
          gcld_ctp(1:nlongr,1:nlatgr), good_norm(idx_err)           )

     CALL write_he5_data_r4 ( 'ColumnFitted',      nlongr, nlatgr,       &
          gcol_reg(1:nlongr,1:nlatgr), good_norm(idx_col)           )
     CALL write_he5_data_r4 ( 'ColumnFittedCorrected',   nlongr, nlatgr, &
          gcol_cor(1:nlongr,1:nlatgr), good_norm(idx_cor)           )
     CALL write_he5_data_r4 ( 'ColumnUncertainty', nlongr, nlatgr,       &
          gcol_err(1:nlongr,1:nlatgr), good_norm(idx_err)           )
     CALL write_he5_data_r4 ( 'RMS', nlongr, nlatgr,                     &
          gcol_rms(1:nlongr,1:nlatgr), good_norm(idx_err)           )

     CALL write_he5_data_r4 ( 'SlantFitted',          nlongr, nlatgr,      &
          gslt_reg(1:nlongr,1:nlatgr), good_norm(idx_col_slt)       )
     CALL write_he5_data_r4 ( 'SlantFittedCorrected', nlongr, nlatgr,      &
          gslt_cor(1:nlongr,1:nlatgr), good_norm(idx_cor_slt)       )
     CALL write_he5_data_r4 ( 'SlantUncertainty',          nlongr, nlatgr, &
          gslt_err(1:nlongr,1:nlatgr), good_norm(idx_err_slt)       )

     CALL close_he5_output_file ( )

  END IF

  IF ( ALLOCATED ( gcol_alb ) )  DEALLOCATE ( gcol_alb )
  IF ( ALLOCATED ( gcol_amf ) )  DEALLOCATE ( gcol_amf )
  IF ( ALLOCATED ( gcol_reg ) )  DEALLOCATE ( gcol_reg )
  IF ( ALLOCATED ( gcol_cor ) )  DEALLOCATE ( gcol_cor )
  IF ( ALLOCATED ( gcol_err ) )  DEALLOCATE ( gcol_err )
  IF ( ALLOCATED ( gcol_reg ) )  DEALLOCATE ( gslt_reg )
  IF ( ALLOCATED ( gcol_cor ) )  DEALLOCATE ( gslt_cor )
  IF ( ALLOCATED ( gcol_err ) )  DEALLOCATE ( gslt_err )
  IF ( ALLOCATED ( good_idx ) )  DEALLOCATE ( good_idx )
  IF ( ALLOCATED ( gcol_aer ) )  DEALLOCATE ( gcol_aer )
  IF ( ALLOCATED ( gcol_qfl ) )  DEALLOCATE ( gcol_qfl )
  IF ( ALLOCATED ( gcol_num ) )  DEALLOCATE ( gcol_num )
  IF ( ALLOCATED ( gcol_rms ) )  DEALLOCATE ( gcol_rms )
  IF ( ALLOCATED ( gsrf_alt ) )  DEALLOCATE ( gsrf_alt )
  IF ( ALLOCATED ( gcld_cfr ) )  DEALLOCATE ( gcld_cfr )
  IF ( ALLOCATED ( gcld_ctp ) )  DEALLOCATE ( gcld_ctp )

  RETURN
END SUBROUTINE gridding_process

SUBROUTINE datafile_loop (                                             &
     pge_esdt,                                                         &
     nXtrack, nTimes, l2_swathfile_id, l2_swath_id,                    &
     nlongr, nlatgr, dlongr, dlatgr, grid_lon, grid_lat,               &
     latmin, latmax, lonmin, lonmax, szamax, cld_frc_min, cld_frc_max, &
     max_area, yn_gpix_weight, yn_ucert_weight, errwght,               &
     yn_use_rbszoom, yn_amf_geo, yn_remove_bg, qflg_max,               &
     xtrack_min, xtrack_max)

  USE SAO_OMIL2_ReadLib_basic_module
  USE SAO_OMIL2_ReadLib_he5_module
  USE SAO_OMIL2_ReadLib_Tessel_module
  USE SAO_OMIL2_median_module
  USE SAO_OMIL2_allocation_module, ONLY: &
       qflag, lat, lon, sza, srf_alb, srf_alt, amf, vza,   &
       corlat, corlon, col_rms, col_reg, col_cor, col_err, &
       slt_reg, slt_cor, slt_err, xtqf, amf, col_rms,     &
       srf_alt, cld_cfr, cld_ctp, corflag, amfflag
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*), INTENT (IN) :: pge_esdt
  INTEGER (KIND=i4), INTENT (IN) :: nXtrack, nTimes, l2_swathfile_id, l2_swath_id, nlongr, nlatgr
  REAL    (KIND=r4), INTENT (IN) :: latmin, latmax, lonmin, lonmax, szamax, cld_frc_min, cld_frc_max
  REAL    (KIND=r4), INTENT (IN) :: max_area, dlongr, dlatgr
  REAL    (KIND=r8), INTENT (IN) :: errwght
  LOGICAL,           INTENT (IN) :: yn_gpix_weight, yn_ucert_weight, yn_use_rbszoom
  LOGICAL,           INTENT (IN) :: yn_remove_bg, yn_amf_geo
  REAL    (KIND=r4), DIMENSION (1:nlongr), INTENT (IN) :: grid_lon
  REAL    (KIND=r4), DIMENSION (1:nlatgr), INTENT (IN) :: grid_lat
  INTEGER (KIND=i2),                       INTENT (IN) :: qflg_max, xtrack_min, xtrack_max

  ! ----------------
  ! Output variables
  ! ----------------
  ! <<<ALL PASSED BY MODULE NOW>>>
  ! -------------------------------

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4), DIMENSION (0:nTimes-1) :: ascending_part
  REAL    (KIND=r8), DIMENSION (4)          :: pix_lat, pix_lon

  ! ------------------------------------------------
  ! Variables related to cross-track pixel positions
  ! ------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: xtrshift = 15
  INTEGER (KIND=i4)            :: spix, epix
  LOGICAL                      :: yn_reb_spat_zoom

  ! ---------
  ! Yes or No
  ! ---------
  LOGICAL, PARAMETER :: yes = .TRUE., no = .FALSE.

  INTEGER (KIND=i4) :: it, ix, iup, ido, he5stat, j1, j2, estat

  ! --------------------------------------------
  ! Variables related to the actual tessellation
  ! --------------------------------------------
  LOGICAL                                 :: yn_polar_exception, yn_good_pix
  INTEGER (KIND=i4)                       :: dateline_offset, tess_orient
  INTEGER (KIND=i4)                       :: nlon_tess, nlat_tess
  INTEGER (KIND=i4), DIMENSION (5)        :: tess_idx
  INTEGER (KIND=i4), DIMENSION (2)        :: tess_idx_abs
  REAL    (KIND=r8), DIMENSION (6,4)      :: tess_pars
  REAL    (KIND=r8), DIMENSION (4,2)      :: tess_satpix
  REAL    (KIND=r8), DIMENSION (1:nlongr) :: tess_lonpts, longr_r8
  REAL    (KIND=r8), DIMENSION (1:nlatgr) :: tess_latpts, latgr_r8

  ! --------
  ! OMI data
  ! --------
  INTEGER (KIND=i4), DIMENSION (3) :: ymd


  ! ------------------
  ! Auxilary variables
  ! ------------------
  REAL (KIND=r4), PARAMETER :: dp5 = 0.5_r4
  REAL (KIND=r4)            :: lonmin_p5, lonmax_p5, latmin_p5, latmax_p5, bgval
  REAL (KIND=r8)            :: dlongr_r8, dlatgr_r8, amfr8

  ! ----------------------------------------------------------------
  ! Variables for eliminating along-track stripes (OMCHOCHO problem)
  ! ----------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (nxtrack+1) :: xtravg, n_xtravg
  REAL (KIND=r8)                        :: mean_xtr, median_xtr, thresh_xtr


  ! -----------------------------
  ! ALLOCATE memory for OMI data
  ! -----------------------------
  ALLOCATE ( qflag  (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'qflag  '
  ALLOCATE ( lat    (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'lat    '
  ALLOCATE ( lon    (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'lon    '
  ALLOCATE ( sza    (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'sza    '
  ALLOCATE ( vza    (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'vza    '
  ALLOCATE ( cld_cfr(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'cld_cfr'
  ALLOCATE ( cld_ctp(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'cld_ctp'
  ALLOCATE ( srf_alb(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'srf_alb'
  ALLOCATE ( srf_alt(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'srf_alt'
  ALLOCATE ( corlat (0:nxtrack,0:ntimes  ), STAT=estat ) ; IF ( estat /= 0 ) STOP 'corlat '
  ALLOCATE ( corlon (0:nxtrack,0:ntimes  ), STAT=estat ) ; IF ( estat /= 0 ) STOP 'corlon '
  ALLOCATE ( col_reg(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'col_reg'
  ALLOCATE ( col_cor(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'col_cor'
  ALLOCATE ( col_err(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'col_err'
  ALLOCATE ( col_rms(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'col_rms'
  ALLOCATE ( xtqf   (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'xtqf   '
  ALLOCATE ( corflag(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'corflag'
  ALLOCATE ( amfflag(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'amfflag'
  ALLOCATE ( amf    (1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'amf    '
  ALLOCATE ( slt_reg(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'slt_reg'
  ALLOCATE ( slt_cor(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'slt_cor'
  ALLOCATE ( slt_err(1:nXtrack,0:nTimes-1), STAT=estat ) ; IF ( estat /= 0 ) STOP 'slt_err'


  ! ------------------------------------------------------------------
  ! Enlarge longitude and latitude range variables by 0.5 deg. Forcing
  ! EXACT limits can lead to ragged image edges when plotted.
  ! ------------------------------------------------------------------
  lonmin_p5 = lonmin - dp5 ; lonmax_p5 = lonmax + dp5
  latmin_p5 = latmin - dp5 ; latmax_p5 = latmax + dp5

  he5stat = HE5_EHrdglatt ( l2_swathfile_id, 'GranuleYear',  YMD(1) )
  he5stat = HE5_EHrdglatt ( l2_swathfile_id, 'GranuleMonth', YMD(2) )
  he5stat = HE5_EHrdglatt ( l2_swathfile_id, 'GranuleDay',   YMD(3) )
 
  CALL saopge_l2_read_geoloc (                        &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, &
       Lat_k=lat, Lon_k=lon, SZA_k=sza, VZA_k=vza,    &
       TerHgt_k=srf_alt, XTQF_k = xtqf)

  CALL saopge_l2_read_datafields (                                         &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0,             &
       fitcol_k=col_reg, fiterr_k=col_err, fitrms_k=col_rms,               &
       qaflg_k=qflag, pclon_k=corlon, pclat_k=corlat, fitcolcor_k=col_cor)
  
  ! ------------------
  ! Screen for clouds. 
  ! ------------------
  CALL saopge_l2_read_datafields (                              &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0,  &
       amfcfr_k=cld_cfr                                           )
  CALL saopge_l2_read_datafields (                              &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0,  &
       amfprs_k=cld_ctp                                           )

  ! ------------------------------------------------
  ! Read molecular or geometric AMF & AMF diagnostic
  ! ------------------------------------------------
  IF (.NOT. yn_amf_geo) THEN
     CALL saopge_l2_read_datafields (                                       &
          he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0, amf_k=amf, amfdiag_k=amfflag)
     CALL saopge_l2_read_datafields (                                       &
          he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0, amfdiag_k=amfflag)
  ELSE IF (yn_amf_geo) THEN
     CALL saopge_l2_read_datafields (                                       &
          he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0, amfgeo_k=amf )
  ENDIF

  ! ---------------------------------------------------------------
  ! Apply the AMF to regular and reference sector corrected columns
  ! ---------------------------------------------------------------
  DO ix = 1, nXtrack
     DO it = 0, nTimes-1
        ! ---------------------------------
        ! Skip current pixel if AMF is <= 0 
        ! ---------------------------------    
        IF ( amf(ix,it) <= 0.0_r4) CYCLE
        ! ------------------------------------------
        ! Undo AMF, but only for non-missing columns
        ! ------------------------------------------
        amfr8 = REAL(amf(ix,it), KIND=r8)
        IF ( col_reg(ix,it) > -1.0E+30_r8 ) slt_reg(ix,it) = col_reg(ix,it) * amfr8
        IF ( col_err(ix,it) > -1.0E+30_r8 ) slt_cor(ix,it) = col_cor(ix,it) * amfr8
        IF ( col_err(ix,it) > -1.0E+30_r8 ) slt_err(ix,it) = col_err(ix,it) * amfr8
     END DO
  END DO

  ! --------------------
  ! Read surface albedo.
  ! --------------------
  CALL saopge_l2_read_datafields (                              &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0, amfalb_k=srf_alb )

  ! ------------------------------
  ! Read reference correction flag
  ! ------------------------------
  CALL saopge_l2_read_datafields (                              &
       he5stat, l2_swath_id, 1, nXtrack, 0, nTimes-1, 0, 0, 0, corflag_k=corflag )

  ! -------------------------------
  ! Check for background adjustment
  ! -------------------------------
  IF ( (TRIM(ADJUSTL(pge_esdt)) == 'OMHCHO') .AND. yn_remove_bg ) THEN
     CALL hcho_get_background_value ( ymd(1), ymd(2), bgval )
  ELSE
     bgval = 0.0_r4
  END IF

  ! -----------------------------------
  ! Set range of cross-track positions
  ! ----------------------------------
  spix = 1  ;  epix = nXtrack ; yn_reb_spat_zoom = .FALSE.
  IF ( ALL ( lon(                 1:xtrshift,0:nTimes-1) < -180.0_r4 ) .AND. &
       ALL ( lon(nxtrack/2+xtrshift:nxtrack, 0:nTimes-1) < -180.0_r4 )        ) THEN
     spix = xtrshift+1 ; epix = nXtrack/2+xtrshift ; yn_reb_spat_zoom = .TRUE.
  END IF

  ascending_part = 1
  DO it = 1, ntimes/2
     iup = MIN ( ntimes/2+it, ntimes-1 )
     ido = MAX ( ntimes/2-it, 0        )
     IF ( lat(nxtrack/2,ido)  >= -90.0_r4    .AND. &
          ascending_part(ido) == 1           .AND. &
          lat(nxtrack/2,ido)  >  lat(nxtrack/2,ido+1) ) ascending_part(0:ido) = 0
     IF ( lat(nxtrack/2,iup)  <= +90.0_r4    .AND. &
          ascending_part(iup) == 1           .AND. &
          lat(nxtrack/2,iup)  <  lat(nxtrack/2,iup-1) ) ascending_part(iup:ntimes-1) = 0
  END DO

  swathlines: DO it = 0, ntimes - 1             

     IF ( ascending_part(it) == 0 ) CYCLE

     xtrack: DO ix = 1, nxtrack

        ! ------------------------------------------------------------
        ! Various reasons for skipping the current pixels:
        !   (*) Geolocation, with margin, outside the gridding domain
        !   (*) Quality Flag shows datum isn't good enough
        !   (*) Solar Zenith Angle too large
        !   (*) Pixel cloud coverage too large
        !   (*) We have a rebinned spatial zoom but don't want
        !   (*) Don't use pixels affected by row anomaly
        !   (*) Only pixels with AMF calculation
        !   (*) Only pixels inbetween xtrack_min and xtrack_max.
        !       Posible conflict with zoom files but I'm not planning
        !       on using them
        ! ------------------------------------------------------------
        IF ( ( lon(ix,it) > lonmax_p5 .OR. lon(ix,it) < lonmin_p5 ) .OR. &
             ( lat(ix,it) > latmax_p5 .OR. lat(ix,it) < latmin_p5 ) .OR. &
             ( ABS(qflag(ix,it))  >  qflg_max                     ) .OR. &
             ( col_err(ix,it)     <= 0.0_r8                       ) .OR. &
             ( ABS(   sza(ix,it)) >  szamax                       ) .OR. &
             ( cld_cfr(ix,it)     <  cld_frc_min                  ) .OR. &
             ( cld_cfr(ix,it)     >  cld_frc_max                  ) .OR. &
             ( yn_reb_spat_zoom .AND. (.NOT. yn_use_rbszoom)      ) .OR. &
             ( xtqf(ix,it) > 0_i2                                 ) .OR. &             
             ( amfflag(ix,it) < 0_i2                              ) .OR. &
             ( corflag(ix,it) < 0_i2                              ) .OR. &
             ( ix .LT. xtrack_min                                 ) .OR. &
             ( ix .GT. xtrack_max                                 ) ) CYCLE
        
        DO j1 = ix-1, ix
           DO j2 = it, it+1
              IF ( ABS(corlon(j1,j2)) < 1.0e-10_r4 ) corlon(j1,j2) = 0.0_r4
              IF ( ABS(corlat(j1,j2)) < 1.0e-10_r4 ) corlat(j1,j2) = 0.0_r4
           END DO
        END DO

        ! -------------------------------------------------------------------
        ! Assign pixel coordinates clockwise starting from bottom left corner
        ! -------------------------------------------------------------------
        pix_lon(1:4) = REAL ( (/ corlon(ix-1,it), corlon(ix,it), corlon(ix,it+1), corlon(ix-1,it+1) /), KIND=r8 )
        pix_lat(1:4) = REAL ( (/ corlat(ix-1,it), corlat(ix,it), corlat(ix,it+1), corlat(ix-1,it+1) /), KIND=r8 )
        
        ! ------------------------------------------------------------------
        ! Cycle if we don't have physical boundaries (edges of Zoom orbits!)
        ! ------------------------------------------------------------------
        IF ( ANY( ABS(pix_lon) > 180.0_r8 ) .OR. ANY( ABS(pix_lat) > 90.0_r8 ) ) CYCLE

        longr_r8 (1:nlongr) = REAL ( grid_lon(1:nlongr), KIND=r8 )
        latgr_r8 (1:nlatgr) = REAL ( grid_lat(1:nlatgr), KIND=r8 )
        dlongr_r8           = REAL(dlongr, KIND=r8)
        dlatgr_r8           = REAL(dlatgr, KIND=r8)
        CALL prepare_omi_tessell (                                                    &
             pix_lat, pix_lon, dlongr_r8, dlatgr_r8, nlongr, nlatgr,                  &
             longr_r8(1:nlongr), latgr_r8(1:nlatgr),                                  &
             tess_polar_lat, nlon_tess, tess_lonpts, nlat_tess, tess_latpts, tess_idx,&
             tess_satpix, tess_pars, tess_orient, tess_idx_abs, dateline_offset,      &
             yn_polar_exception, yn_good_pix )

        IF ( yn_polar_exception .OR. (.NOT. yn_good_pix) ) CYCLE

        ! ----------------------------
        ! Correct for background value
        ! ---------------------------- 
        col_reg(ix,it) = col_reg(ix,it) - bgval

        CALL gridding_loop (                                                 &
             nlongr, nlatgr, nlon_tess, nlat_tess,                           &
             tess_lonpts(1:nlon_tess), tess_latpts(1:nlat_tess),             &
             tess_idx, tess_satpix, tess_pars, tess_orient, tess_idx_abs,    &
             dateline_offset, max_area,                                      &
             col_reg(ix,it), col_cor(ix,it), col_err(ix,it), col_rms(ix,it), &
             srf_alb(ix,it), slt_reg(ix,it), slt_cor(ix,it), slt_err(ix,it), &
             amf(ix,it), srf_alt(ix,it), cld_cfr(ix,it), cld_ctp(ix,it),     &
             yn_gpix_weight, yn_ucert_weight, errwght        )
           
     END DO xtrack
  END DO swathlines

  ! ------------------------------
  ! DEALLOCATE memory for OMI data
  ! ------------------------------
  IF ( ALLOCATED(qflag  ) ) DEALLOCATE (qflag  )
  IF ( ALLOCATED(lat    ) ) DEALLOCATE (lat    )
  IF ( ALLOCATED(lon    ) ) DEALLOCATE (lon    )
  IF ( ALLOCATED(sza    ) ) DEALLOCATE (sza    )
  IF ( ALLOCATED(vza    ) ) DEALLOCATE (vza    )
  IF ( ALLOCATED(srf_alb) ) DEALLOCATE (srf_alb)
  IF ( ALLOCATED(corlat ) ) DEALLOCATE (corlat )
  IF ( ALLOCATED(corlon ) ) DEALLOCATE (corlon )
  IF ( ALLOCATED(col_reg) ) DEALLOCATE (col_reg)
  IF ( ALLOCATED(col_cor) ) DEALLOCATE (col_cor)
  IF ( ALLOCATED(col_err) ) DEALLOCATE (col_err)
  IF ( ALLOCATED(col_rms) ) DEALLOCATE (col_rms)
  IF ( ALLOCATED(xtqf   ) ) DEALLOCATE (xtqf   )
  IF ( ALLOCATED(corflag) ) DEALLOCATE (corflag)
  IF ( ALLOCATED(amf    ) ) DEALLOCATE (amf    )
  IF ( ALLOCATED(amfflag) ) DEALLOCATE (amfflag)
  IF ( ALLOCATED(srf_alt) ) DEALLOCATE (srf_alt)
  IF ( ALLOCATED(cld_cfr) ) DEALLOCATE (cld_cfr)
  IF ( ALLOCATED(cld_ctp) ) DEALLOCATE (cld_ctp)
  IF ( ALLOCATED(slt_reg) ) DEALLOCATE (slt_reg)
  IF ( ALLOCATED(slt_cor) ) DEALLOCATE (slt_cor)
  IF ( ALLOCATED(slt_err) ) DEALLOCATE (slt_err)

  RETURN
END SUBROUTINE datafile_loop

SUBROUTINE gridding_loop (                                           &
     nlongr, nlatgr, nlon_tess, nlat_tess, tess_lonpts, tess_latpts, &
     tess_idx, tess_satpix, tess_pars, tess_orient, tess_idx_abs,    &
     dateline_offset, max_area, col_reg, col_cor, col_err, col_rms,  &
     srf_alb, slt_reg, slt_cor, slt_err, amf, srf_alt, cld_cfr,      &
     cld_ctp, yn_gpix_weight, yn_ucert_weight, errwght               )

  USE SAO_OMIL2_ReadLib_basic_module
  USE SAO_OMIL2_ReadLib_he5_module
  USE SAO_OMIL2_ReadLib_Tessel_module
  USE SAO_OMIL2_allocation_module, ONLY: gcol_aer, gcol_alb, gcol_num, &
                                         gcol_reg, gcol_cor, gcol_err, &
                                         gslt_reg, gslt_cor, gslt_err, &
                                         gcol_amf, gslt_aer, gcol_rms, &
                                         gsrf_alt, gcld_cfr, gcld_ctp
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                          INTENT (IN) :: nlongr, nlatgr, nlon_tess, nlat_tess
  INTEGER (KIND=i4),                          INTENT (IN) :: dateline_offset, tess_orient
  REAL    (KIND=r8), DIMENSION (1:nlon_tess), INTENT (IN) :: tess_lonpts
  REAL    (KIND=r8), DIMENSION (1:nlat_tess), INTENT (IN) :: tess_latpts
  REAL    (KIND=r4),                          INTENT (IN) :: max_area
  REAL    (KIND=r8),                          INTENT (IN) :: errwght
  INTEGER (KIND=i4), DIMENSION (5),           INTENT (IN) :: tess_idx
  INTEGER (KIND=i4), DIMENSION (2),           INTENT (IN) :: tess_idx_abs
  REAL    (KIND=r8), DIMENSION (6,4),         INTENT (IN) :: tess_pars
  REAL    (KIND=r8), DIMENSION (4,2),         INTENT (IN) :: tess_satpix
  REAL    (KIND=r8),                          INTENT (IN) :: col_reg, col_cor, col_err, col_rms
  REAL    (KIND=r8),                          INTENT (IN) :: slt_reg, slt_cor, slt_err, srf_alb
  REAL    (KIND=r4),                          INTENT (IN) :: amf, cld_cfr, cld_ctp
  INTEGER (KIND=i2),                          INTENT (IN) :: srf_alt
  LOGICAL,                                    INTENT (IN) :: yn_gpix_weight, yn_ucert_weight

  ! ------------------
  ! Modified variables
  ! ------------------
  ! <<<ALL PASSED BY MODULE NOW>>>
  ! ------------------------------

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4), DIMENSION (1:nlon_tess)             :: ylimit_upp, ylimit_low
  REAL    (KIND=r8), DIMENSION (1:nlon_tess,1:nlat_tess) :: tess_area

  INTEGER (KIND=i4) :: i, j, ilat, ilon
  REAL    (KIND=r8) :: tess_sum, max_area_r8, frac, errtmp, col_err_norm, &
                       tess_sum_slt, frac_slt, errtmp_slt,  slt_err_norm

  ! --------------------------------------------------------------------------
  ! Here is a good point to make sure MAX_AREA is consistent with the weighting
  ! scheme we use. If we are using ground-pixel size based weighting, then the
  ! R8 variable holds the MAX_AREA as passed from the calling routine. If not,
  ! the value is simply "1.0", because that is what TESS_SUM below gets set to
  ! in that case, and we are checking for "TESS_SUM <= MAX_AREA".
  ! --------------------------------------------------------------------------
  max_area_r8 = 1.0_r8
  IF ( yn_gpix_weight ) max_area_r8 = REAL (max_area, KIND=r8)

  ! --------------------
  ! Tessellate the pixel
  ! --------------------
  tess_area = 0.0_r8 ; tess_sum = 0.0_r8 ; ylimit_low = 0 ; ylimit_upp = 0
  CALL tesselate_areamaster (  &
       tess_lonpts(1:nlon_tess), tess_latpts(1:nlat_tess), nlon_tess, nlat_tess, & !i
       tess_satpix, tess_pars, tess_orient,                                      & !i
       tess_idx(1), tess_idx(2), tess_idx(3), tess_idx(4), tess_idx(5),          & !i
       tess_area(1:nlon_tess,1:nlat_tess), tess_sum, ylimit_low(1:nlon_tess),    & !o
       ylimit_upp(1:nlon_tess) )                                                   !o

  ! ---------------------------------------------------------------------------------
  ! Use tessellated pixel areas as weights to distribute column values over the grid.
  ! ---------------------------------------------------------------------------------
  ! The weight for each grid cell is the area TESS_AREA of overlap of the ground
  ! pixel with the grid cell, divided by the size TESS_SUM of the ground pixel.
  ! We also have to keep track of the total area that goes into one grid cell; 
  ! this value will be used to normalize the grid cell value after all contributions
  ! have been added.
  ! ---------------------------------------------------------------------------------
  DO i = 1, nlon_tess
     DO j = 1, nlat_tess
        IF ( tess_area(i,j) < 0.0_r8 ) tess_area(i,j) = 0.0_r8
     END DO
  END DO
  ! --------------------------------------------------------------
  ! Assign TESS_SUM depending on whether or not we are using
  ! ground-pixel size based weighting. If not, it is simply "1.0".
  ! --------------------------------------------------------------
  tess_sum = 1.0_r8 
  tess_sum_slt = 1.0_r8

  IF ( yn_gpix_weight ) tess_sum = SUM ( tess_area, MASK = tess_area >= 0.0_r8 )

  ! -------------------------------------------------
  ! Proceed only if sum of parts of tessellated pixel
  ! is within physically sensible bounds.
  ! -------------------------------------------------
  IF ( tess_sum > 0.0_r8 .AND. tess_sum <= max_area_r8 ) THEN

     ! ------------------------------------------------------------
     ! Include fitting errors in weighting. COL_ERR_NORM holds the
     ! normalized fitting uncertainties, which are added as inverse
     ! squares. If they are not normalized, we are running the risk
     ! of numerical over- or underflows.
     ! -------------------------------------------------------------
     col_err_norm = col_err
     slt_err_norm = slt_err
     IF ( errwght /= 0.0_r8 ) THEN
        col_err_norm = col_err_norm / errwght
        slt_err_norm = slt_err_norm / errwght
     END IF
     IF ( yn_ucert_weight ) THEN
        tess_sum = tess_sum * (col_err_norm*col_err_norm)
        tess_sum_slt = tess_sum * (slt_err_norm*slt_err_norm)
     END IF

     ! -------------------------
     ! Short-cut for "1/TESS_SUM
     ! -------------------------
     frac     = 1.0_r8 / tess_sum
     frac_slt = 1.0_r8 / tess_sum_slt
     ! ------------------------------
     ! Loop over tessellation indices
     ! ------------------------------
     DO i = 1, nlon_tess

        ilon = tess_idx_abs(1) + dateline_offset + i - 1
        IF ( ilon > nlongr ) ilon = ilon - nlongr

        DO j = 1, nlat_tess
           ilat = tess_idx_abs(2) + j - 1
           
           IF ( tess_area(i,j) <= 0.0_r8 .OR. tess_area(i,j) > max_area_r8 ) CYCLE

           gcol_num(ilon,ilat) = gcol_num(ilon,ilat) + 1
           gcol_aer(ilon,ilat) = gcol_aer(ilon,ilat) + REAL(tess_area(i,j)         * frac,     KIND=r4)
           gslt_aer(ilon,ilat) = gslt_aer(ilon,ilat) + REAL(tess_area(i,j)         * frac_slt, KIND=r4)
           gcol_reg(ilon,ilat) = gcol_reg(ilon,ilat) + REAL(tess_area(i,j)*col_reg * frac,     KIND=r4)
           gcol_cor(ilon,ilat) = gcol_cor(ilon,ilat) + REAL(tess_area(i,j)*col_cor * frac,     KIND=r4)
           gcol_rms(ilon,ilat) = gcol_rms(ilon,ilat) + REAL(tess_area(i,j)*col_rms * frac,     KIND=r4)
           gcol_alb(ilon,ilat) = gcol_alb(ilon,ilat) + REAL(tess_area(i,j)*srf_alb * frac,     KIND=r4)
           gcol_amf(ilon,ilat) = gcol_amf(ilon,ilat) + REAL(tess_area(i,j)*amf     * frac,     KIND=r4)
           gsrf_alt(ilon,ilat) = gsrf_alt(ilon,ilat) + REAL(tess_area(i,j)*REAL(srf_alt, KIND=r4) * frac, KIND=r4)
           gcld_cfr(ilon,ilat) = gcld_cfr(ilon,ilat) + REAL(tess_area(i,j)*cld_cfr * frac,     KIND=r4)
           gcld_ctp(ilon,ilat) = gcld_ctp(ilon,ilat) + REAL(tess_area(i,j)*cld_ctp * frac,     KIND=r4)
           gslt_reg(ilon,ilat) = gslt_reg(ilon,ilat) + REAL(tess_area(i,j)*slt_reg * frac_slt, KIND=r4)
           gslt_cor(ilon,ilat) = gslt_cor(ilon,ilat) + REAL(tess_area(i,j)*slt_cor * frac_slt, KIND=r4)

           ! ----------------------------------------------------------------------
           ! To compute the weighted error, we have to add the inverse square of
           ! the individual errors and take the square root of the inverse at the
           ! end. See for example <http://en.wikipedia.org/wiki/Weighted_mean>.
           !
           ! NOTE that we do NOT use inverse pixel size (or error) weighting here.
           ! Somehow it doesn't sound right to derive an error-weighted uncertainty
           ! and then reduce that by the SQRT of the number of samples.
           ! ----------------------------------------------------------------------
           errtmp     = REAL ( col_err_norm, KIND=r4)
           errtmp_slt = REAL ( slt_err_norm, KIND=r4)
           ! -----------------------------------
           ! Sanity check to avoid division by 0
           ! -----------------------------------
           IF ( ABS(errtmp)     == 0.0_r4 ) errtmp     = 1.0_r4
           IF ( ABS(errtmp_slt) == 0.0_r4 ) errtmp_slt = 1.0_r4
           ! ------------------
           ! Compute the square
           ! ------------------
           errtmp     = errtmp     * errtmp
           errtmp_slt = errtmp_slt * errtmp_slt
           ! ----------------------
           ! Add to composite error
           ! ----------------------
           gcol_err(ilon,ilat) = gcol_err(ilon,ilat) + 1.0_r4 / errtmp
           gslt_err(ilon,ilat) = gslt_err(ilon,ilat) + 1.0_r4 / errtmp_slt

        END DO

     END DO

  END IF

END SUBROUTINE gridding_loop

SUBROUTINE gridded_ascii_output (                              &
     outfile, nlongr, nlatgr, dlongr, dlatgr, grid_norm,       &
     grid_lon, grid_lat, gcol_reg, gcol_err, gcol_dst, gcol_aer )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i4, r4, r8, ounit
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                             INTENT (IN) :: nlongr, nlatgr
  CHARACTER (LEN=*),                             INTENT (IN) :: outfile
  REAL (KIND=r4),                                INTENT (IN) :: dlongr, dlatgr
  REAL (KIND=r4), DIMENSION (1:4),               INTENT (IN) :: grid_norm
  REAL (KIND=r4), DIMENSION (1:nlongr),          INTENT (IN) :: grid_lon
  REAL (KIND=r4), DIMENSION (1:nlatgr),          INTENT (IN) :: grid_lat
  REAL (KIND=r4), DIMENSION (1:nlongr,1:nlatgr), INTENT (IN) :: gcol_reg, gcol_err, gcol_dst, gcol_aer

  ! --------------
  ! Local variable
  ! --------------
  INTEGER (KIND=i4) :: ios

  ! ---------------------------------------------------------------
  ! Writing the results to file. Exit with error if unable to open.
  ! ---------------------------------------------------------------
  ios = 0
  OPEN (UNIT=ounit, FILE=TRIM(ADJUSTL(outfile)), STATUS='UNKNOWN', ACTION='WRITE', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE (*,'(A)') 'ERROR: Unable to open output file <'//TRIM(ADJUSTL(outfile))//'>'
     STOP 1
  END IF
  
  WRITE (UNIT=ounit,FMT='(2I6)')                 nlongr, nlatgr
  WRITE (UNIT=ounit,FMT='(1F12.6)')              dlongr
  WRITE (UNIT=ounit,FMT='(5000000(7F13.4:/))')   grid_lon(1:nlongr)
  WRITE (UNIT=ounit,FMT='(1F12.6)')              dlatgr
  WRITE (UNIT=ounit,FMT='(5000000(7F13.4:/))')   grid_lat(1:nlatgr)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') grid_norm(1)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') gcol_reg(1:nlongr,1:nlatgr)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') grid_norm(2)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') gcol_dst(1:nlongr,1:nlatgr)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') grid_norm(3)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') gcol_err(1:nlongr,1:nlatgr)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') grid_norm(4)
  WRITE (UNIT=ounit,FMT='(5000000(1P7E13.4:/))') gcol_aer(1:nlongr,1:nlatgr)

  CLOSE (ounit)

  RETURN
END SUBROUTINE gridded_ascii_output

SUBROUTINE normalize_gridded_values ( m, n, low_no_good, grid_2d, valn )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i4, r4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                  INTENT (IN) :: m, n
  REAL    (KIND=r4),                  INTENT (IN) :: low_no_good

  ! -------------------------
  ! Output/modified variables
  ! -------------------------
  REAL    (KIND=r4),                  INTENT (OUT)   :: valn
  REAL    (KIND=r4), DIMENSION (m,n), INTENT (INOUT) :: grid_2d

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: ngrid, i, j
  REAL    (KIND=r4) :: meav

  ! ----------------------------------------------------
  ! NOTE: MEDIAN was replaced by MEAN because the MEDIAN
  !       calculation, in its current implementation, is
  !       a rather time-consuming affair.
  ! ----------------------------------------------------

  ! ---------------------------------------------------------------
  ! Find number of grid entries that are > 0; this setting has been
  ! taken care of before the call of this normalization routine.
  ! ---------------------------------------------------------------
  ngrid = COUNT ( grid_2d(1:m,1:n) > low_no_good )

  IF ( ngrid > 0 ) THEN
     ! ----------------------------------
     ! Calculate MEAN of good grid points
     ! ----------------------------------
     meav = SUM ( grid_2d(1:m,1:n), MASK=(grid_2d(1:m,1:n) > low_no_good) ) / REAL(ngrid, KIND=r4)

     ! -------------------------------------------------------
     ! Find a nice-looking value for normalization - something
     ! like the closest INTEGER base of the MEAN.
     ! -------------------------------------------------------
     CALL get_normval ( meav, valn )

     ! ----------------------------------------------------
     ! Normalize. Set the normalization factor to something
     ! sensible, since it will be written to file.
     ! ----------------------------------------------------
     IF ( valn == 0.0_r8 ) valn = 1.0_r8
     DO i = 1, m
        DO j = 1, n
           IF ( grid_2d(i,j) > low_no_good ) grid_2d(i,j) = grid_2d(i,j) / valn
        END DO
     END DO
  END IF

  RETURN
END SUBROUTINE normalize_gridded_values


SUBROUTINE get_normval ( val, valn )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i4, r4, r8
  IMPLICIT NONE

  REAL (KIND=r4), INTENT (IN)  :: val
  REAL (KIND=r4), INTENT (OUT) :: valn

  ! ---------------------------------------------------
  ! Find an integral value for normalization, based on
  ! a given input value.
  ! ---------------------------------------------------

  ! ---------------
  ! Local variables
  ! ---------------
  REAL    (KIND=r4) :: bas
  INTEGER (KIND=i4) :: pow

  pow    = INT(log10(val))
  bas    = REAL(ANINT(val/10.0_r4**pow))
  valn   = bas*(10.0_r4**pow)

  RETURN
END SUBROUTINE get_normval

SUBROUTINE hcho_get_background_value ( year, month, bgval )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ------------
  ! Input values
  ! ------------
  INTEGER (KIND=i4), INTENT (IN) :: year, month

  ! ------------
  ! Output value
  ! ------------
  REAL (KIND=r4), INTENT (OUT) :: bgval


  ! ----------------------------------------------------------
  ! The background values are based on HCHO monthly averages 
  ! derived over the remote Pacific, in the region 
  !            +/-50 deg latitude, 180+/-15 deg longitude
  ! for the time period of 10/2004 through 09/2009. The MEDIAN
  ! of the averaged data was taken, and a 4th order polynomial
  ! was fitted to the time series. Due to lack of data, values
  ! for 08+09/2004 and 10-12/2009 were extrapolated from the 
  ! fitted curve.
  ! ----------------------------------------------------------

  ! ----------------------------------------------------------
  ! We have 65 months: 5 in 2004 and 12 each in 2005-2009
  ! ----------------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: nmonths = 65
  REAL    (KIND=r4), DIMENSION (nmonths), PARAMETER :: hcho_background = (/ &
       1.86348e+15_r4, 2.10683e+15_r4, 2.33003e+15_r4, 2.53413e+15_r4, 2.72016e+15_r4, &
       2.88913e+15_r4, 3.04204e+15_r4, 3.17985e+15_r4, 3.30352e+15_r4, 3.41397e+15_r4, &
       3.51213e+15_r4, 3.59889e+15_r4, 3.67511e+15_r4, 3.74166e+15_r4, 3.79936e+15_r4, &
       3.84903e+15_r4, 3.89147e+15_r4, 3.92745e+15_r4, 3.95772e+15_r4, 3.98302e+15_r4, &
       4.00407e+15_r4, 4.02155e+15_r4, 4.03616e+15_r4, 4.04854e+15_r4, 4.05932e+15_r4, &
       4.06914e+15_r4, 4.07857e+15_r4, 4.08821e+15_r4, 4.09861e+15_r4, 4.11030e+15_r4, &
       4.12381e+15_r4, 4.13962e+15_r4, 4.15824e+15_r4, 4.18010e+15_r4, 4.20565e+15_r4, &
       4.23531e+15_r4, 4.26948e+15_r4, 4.30854e+15_r4, 4.35284e+15_r4, 4.40274e+15_r4, &
       4.45855e+15_r4, 4.52057e+15_r4, 4.58908e+15_r4, 4.66435e+15_r4, 4.74662e+15_r4, &
       4.83610e+15_r4, 4.93301e+15_r4, 5.03751e+15_r4, 5.14979e+15_r4, 5.26998e+15_r4, &
       5.39820e+15_r4, 5.53455e+15_r4, 5.67913e+15_r4, 5.83200e+15_r4, 5.99320e+15_r4, &
       6.16275e+15_r4, 6.34066e+15_r4, 6.52691e+15_r4, 6.72148e+15_r4, 6.92431e+15_r4, &
       7.13532e+15_r4, 7.35442e+15_r4, 7.58150e+15_r4, 7.81642e+15_r4, 8.05903e+15_r4   /)

  INTEGER (KIND=i4) :: imonth

  bgval = 0.0_r4

  ! --------------------------------------------------------------
  ! Set IMONTH based on input of YEAR and MONTH. The data start in
  ! 08/2004, so YEAR=2004, MONTH=8 has to translate to IMONTH=1.
  ! --------------------------------------------------------------
  imonth = ( year - 2004_i4 ) * 12_i4 + month - 7

  IF ( imonth < 1 .OR. imonth > nmonths ) RETURN

  bgval = hcho_background ( imonth )

  RETURN
END SUBROUTINE hcho_get_background_value

SUBROUTINE omi_avg_read_input (                                       &
     input_control_file,                                              &
     pge_esdt, listfile, outfile, output_format, swathname,           &
     yn_norm_output, yn_gpix_weight, yn_ucert_weight, yn_use_rbszoom, &
     yn_amf_geo,  yn_remove_bg,                                       &
     qflg_max, szamax, cld_frc_min, cld_frc_max, errwght,             &
     lonmin, lonmax, dlongr, latmin, latmax, dlatgr,                  &
     xtrackmin, xtrackmax                                             )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i2, i4, r4, r8
  IMPLICIT NONE

  ! --------------
  ! Input variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: input_control_file

  ! ----------------
  ! Output variables
  ! ----------------
  CHARACTER (LEN=*),   INTENT (OUT) :: pge_esdt, listfile, outfile, output_format, swathname
  LOGICAL,             INTENT (OUT) :: yn_norm_output, yn_gpix_weight, yn_ucert_weight, yn_use_rbszoom
  LOGICAL,             INTENT (OUT) :: yn_amf_geo, yn_remove_bg
  INTEGER   (KIND=i2), INTENT (OUT) :: qflg_max, xtrackmin, xtrackmax
  REAL      (KIND=r8), INTENT (OUT) :: errwght
  REAL      (KIND=r4), INTENT (OUT) :: szamax, cld_frc_min, cld_frc_max
  REAL      (KIND=r4), INTENT (OUT) :: lonmin, lonmax, dlongr, latmin, latmax, dlatgr

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4), PARAMETER :: ipunit = 66
  INTEGER (KIND=i4)            :: ios
  LOGICAL                      :: yn_fail

  ! -----------------------------------------------------
  ! File Mark Strings to search for in input control file
  ! -----------------------------------------------------
  CHARACTER (LEN=19), PARAMETER :: fm_esdt = 'OMI ESDT to average'
  CHARACTER (LEN=34), PARAMETER :: fm_list = 'File name with list of input files'
  CHARACTER (LEN=16), PARAMETER :: fm_outf = 'Output file name'
  CHARACTER (LEN=13), PARAMETER :: fm_ofmt = 'Output format'
  CHARACTER (LEN=16), PARAMETER :: fm_norm = 'Normalize output'
  CHARACTER (LEN=31), PARAMETER :: fm_qflg = 'Maximum quality flag to include'
  CHARACTER (LEN=26), PARAMETER :: fm_szen = 'Maximum solar zenith angle'
  CHARACTER (LEN=29), PARAMETER :: fm_glim = 'Longitude and latitude limits'
  CHARACTER (LEN=22), PARAMETER :: fm_cfrc = 'Maximum cloud fraction'
  CHARACTER (LEN=27), PARAMETER :: fm_pixw = 'Ground pixel size weighting'
  CHARACTER (LEN=21), PARAMETER :: fm_errw = 'Uncertainty weighting'
  CHARACTER (LEN=34), PARAMETER :: fm_rszg = 'Use rebinned spatial zoom granules'
  CHARACTER (LEN=25), PARAMETER :: fm_amfg = 'Geometric air mass factor'
  CHARACTER (LEN=22), PARAMETER :: fm_bgrm = 'HCHO remove background'
  CHARACTER (LEN=18), PARAMETER :: fm_xtra = 'xtrack pixel range'

  ! -----------------------
  ! Open input control file
  ! -----------------------
  OPEN ( UNIT=ipunit, FILE=TRIM(ADJUSTL(input_control_file)), STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE (*,'(A)') 'ERROR: Unable to open file <'//TRIM(ADJUSTL(input_control_file))//'>'
     STOP 1
  END IF

  ! ---------------------------------------------------------
  ! Name of the OMI product to grid. E.g., OMBRO, OMHCHO, etc
  ! ---------------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_esdt, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, '(A)') pge_esdt

  ! -----------------------------------------------------
  ! Name of file with list of input files (absolute path)
  ! -----------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_list, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, '(A)') listfile

  ! -----------------------------------
  ! Name of output file (absolute path)
  ! -----------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_outf, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, '(A)') outfile

  ! ---------------------------------------------------------
  ! Output format ("ascii", "he5"; can be more than one) and,
  ! if "he5", HE5 Swath Name
  ! ---------------------------------------------------------
  REWIND (ipunit)
  swathname = '?????'
  CALL skip_to_filemark ( ipunit, fm_ofmt, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, '(A)') output_format
  IF ( INDEX(TRIM(ADJUSTL(output_format)), 'he5') > 0 ) READ  ( ipunit, '(A)' ) swathname
  
  ! -----------------------------------------------------------------
  ! Output normalization ("T" or "F"); a scale factor is determined
  ! automatically, based on the rounded MEDIAN value of the averages.
  ! -----------------------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_norm, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) yn_norm_output

  ! -------------------------------------------------------
  ! Maximum quality flag to include (currently  0, 1, or 2)
  ! -------------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_qflg, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) qflg_max

  ! --------------------------
  ! Maximum solar zenith angle
  ! --------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_szen, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) szamax

  ! ------------------------------
  ! Longitude and latitude limits;
  ! (1) Lon_min Lon_max dLon
  ! (2) Lat_min Lat_max dLat
  ! ------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_glim, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) lonmin, lonmax, dlongr
  READ (ipunit, *) latmin, latmax, dlatgr

  ! -------------------------------------------------------------
  ! Maximum cloud fraction (ignored for some ESDTs, e.g., OMOCLO)
  ! Input range: 0.0 <= cfr <= 1.0
  ! -------------------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_cfrc, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) cld_frc_min, cld_frc_max

  ! -----------------------------------------------------------------
  ! Weighting of averages inversely proportional to ground pixel size
  ! "T" or "F"
  ! -----------------------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_pixw, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) yn_gpix_weight

  ! -------------------------------------------------------------------
  ! Weighting of averages inversely proportional to fitting uncertainty
  ! "T" or "F" followed by REAL number to use as scale factor to bring
  ! uncertainties to ~1 (e.g., 1.0E+13 for OMBRO, 1.0e15 for OMHCHO)
  ! ------------------------------------------------------------------
  REWIND (ipunit)
  errwght = 1.0_r8
  CALL skip_to_filemark ( ipunit, fm_errw, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) yn_ucert_weight, errwght
  IF ( errwght == 0.0_r8 ) errwght = 1.0_r8

  ! ----------------------------------------------------
  ! Include rebinned spatial zoom granules in averaging?
  ! ----------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_rszg, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) yn_use_rbszoom

  ! ---------------------------------------------------
  ! Use geometric air mass factor ("T" or "F")
  ! For OMOCLO, "F" will product slant column averages.
  ! ---------------------------------------------------
  REWIND (ipunit)
  CALL skip_to_filemark ( ipunit, fm_amfg, yn_fail ) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) yn_amf_geo

  ! -----------------------------------
  ! HCHO only: remove background values
  ! -----------------------------------
  yn_remove_bg = .FALSE.
  IF ( TRIM(ADJUSTL(pge_esdt)) == 'OMHCHO' ) THEN
     REWIND (ipunit)
     CALL skip_to_filemark ( ipunit, fm_bgrm, yn_fail ) ; IF ( yn_fail ) STOP 1
     READ (ipunit, *) yn_remove_bg
  END IF

  ! ------------------
  ! Xtrack pixel range
  ! ------------------
  REWIND(ipunit)
  CALL skip_to_filemark (ipunit, fm_xtra, yn_fail) ; IF ( yn_fail ) STOP 1
  READ (ipunit, *) xtrackmin, xtrackmax

  CLOSE ( ipunit )

  RETURN
END SUBROUTINE omi_avg_read_input


SUBROUTINE skip_to_filemark ( funit, lm_string, yn_fail )

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i4, maxchlen
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4), INTENT (IN) :: funit
  CHARACTER (LEN=*),   INTENT (IN) :: lm_string

  ! ---------------
  ! Output variable
  ! ---------------
  LOGICAL, INTENT (OUT) :: yn_fail

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: lmlen, ios
  CHARACTER (LEN=maxchlen) :: tmpline

  yn_fail = .FALSE.

  ! -------------------------------------------
  ! Determine the length of the string landmark
  ! -------------------------------------------
  lmlen = LEN(TRIM(ADJUSTL(lm_string)))

  ! -------------------------------------------------------
  ! Read lines in the file until we either find the string,
  ! reach the end of the file, or reading fails otherwise.
  ! ----------------------------------------------------
  ios = 0
  getlm: DO WHILE ( ios == 0 )
     READ (UNIT=funit, FMT='(A)', IOSTAT=ios) tmpline
     tmpline = TRIM(ADJUSTL(tmpline))
     IF ( ios /= 0 .OR. tmpline(1:lmlen) == lm_string ) EXIT getlm
  END DO getlm
  
  ! -------------------------
  ! Check for successful READ
  ! -------------------------
  yn_fail = .FALSE.
  IF ( ios /= 0 ) THEN
     yn_fail = .TRUE.
     WRITE (*,'(A)') 'ERROR: Unable to locate string <'//TRIM(ADJUSTL(lm_string))//'> in input file.'
  END IF

  RETURN
END SUBROUTINE skip_to_filemark
