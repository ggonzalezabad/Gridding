MODULE SAO_OMIL2_ReadLib_he5_module

  USE  SAO_OMIL2_ReadLib_Basic_module
  IMPLICIT NONE

  INCLUDE 'hdfeos5.inc'
  INTEGER (KIND = i4), EXTERNAL :: &
       he5_ehrdglatt , he5_ehwrglatt, he5_swattach,  he5_swclose,   he5_swcreate, &
       he5_swdefdfld , he5_swdefdim,  he5_swdefgfld, he5_swdetach,  he5_swopen,   &
       he5_swrdattr,   he5_swrdfld,   he5_swrdgattr, he5_swrdlattr, he5_swwrattr, &
       he5_swwrfld,    he5_swwrgattr, he5_swwrlattr

  INTEGER (KIND = i48), EXTERNAL :: &
       HE5_SWinqswath, HE5_SWinqdims, HE5_SWinqdfldalias

CONTAINS

  SUBROUTINE saopge_l2_open_attach ( he5stat, l2_fname, l2_swathfile_id, l2_swath_id, l2_swath )

    ! ------------------------------------------------------
    ! Opens the L2 output file READ ONLY and attach to swath
    ! ------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: l2_fname

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER   (KIND=i4),      INTENT (OUT) :: he5stat, l2_swathfile_id, l2_swath_id
    CHARACTER (LEN=maxchlen), INTENT (OUT) :: l2_swath

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)  :: locestat
    INTEGER (KIND=i48) :: swlen48, locestat48


    he5stat = estat_ok

    ! -----------------------------------------------------------
    ! Open HE5 output file and check SWATH_FILE_ID ( -1 if error)
    ! -----------------------------------------------------------
    l2_swathfile_id = HE5_SWopen ( TRIM(ADJUSTL(l2_fname)), HE5F_ACC_RDonly )
    IF ( l2_swathfile_id == he5_stat_fail ) THEN
       WRITE (*,'(A)') 'FATAL ERROR: HE5_SWOpen failed.'
       he5stat = estat_fatal; RETURN
    END IF

    ! ---------------------------------------------
    ! Check for existing HE5 swath and attach to it
    ! ---------------------------------------------
    l2_swath = ' '
    locestat48  = HE5_SWinqswath  ( TRIM(ADJUSTL(l2_fname)), l2_swath, swlen48 )
    l2_swath_id = HE5_SWattach ( l2_swathfile_id, TRIM(ADJUSTL(l2_swath)) )
    IF ( l2_swath_id == he5_stat_fail ) THEN
       IF ( l2_swath_id == he5_stat_fail )  WRITE (*,'(A)') 'FATAL ERROR: HE5_SWattach failed.'
       he5stat = estat_fatal
    END IF

    RETURN
  END SUBROUTINE saopge_l2_open_attach


  SUBROUTINE saopge_l2_read_dimensions ( &
       he5stat, l2_swath_id, nTimes_k, nXtrack_k, nFitElem_k, nClenFitElem_k )

    ! -----------------------
    ! Return swath dimensions
    ! -----------------------
 
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: l2_swath_id

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),           INTENT (OUT) :: he5stat
    INTEGER (KIND=i4), OPTIONAL, INTENT (OUT) :: &
         nTimes_k, nXtrack_k, nFitElem_k, nClenFitElem_k

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4), PARAMETER             :: maxdim = 100
    INTEGER   (KIND=i4)                        :: ndim, nt, nx, nf, nc, nmin
    INTEGER   (KIND=i4)                        :: i, j, locestat, swlen, iend, istart
    INTEGER   (KIND=i48)                       :: ndim48
    INTEGER   (KIND=i4),  DIMENSION(0:maxdim)  :: dim_array
    INTEGER   (KIND=i48), DIMENSION(0:maxdim)  :: dim_array48
    CHARACTER (LEN=maxchlen)                   :: dim_chars

    he5stat = estat_ok

    ! ------------------------------
    ! Inquire about swath dimensions
    ! ------------------------------
    ndim48 = HE5_SWinqdims  ( l2_swath_id, dim_chars, dim_array48(0:maxdim) )
    ndim   = INT(ndim48, KIND=i48)
    dim_array(0:ndim) = INT ( dim_array48(0:ndim), KIND=i4 )

!    ! ----------------------------------------
!    ! Since we don't know the total number of dimensions, we keep increasing
!    ! them until the first ZERO appears. ZERO indicates unset dimensions.
!    ! ----------------------------------------------------------------------
!    nmin = 1 ; j = 0
!    readdim: DO WHILE ( nmin > 0 )
!       locestat48  = HE5_SWinqdims  ( l2_swath_id, dim_chars, dim_array(0:j) )
!       nmin = MINVAL(dim_array(0:j))
!       IF ( nmin == 0 ) THEN
!          ndim = j - 1
!          EXIT readdim
!       ELSE
!          j = j + 1
!       END IF
!    END DO readdim

    dim_chars = TRIM(ADJUSTL(dim_chars))
    swlen = LEN_TRIM(ADJUSTL(dim_chars))

    ! -------------------------------------
    ! Extract the dimensions from the swath
    ! -------------------------------------
    nt = -1_i4 ; nx = -1_i4 ; nf = -1_i4 ; nc = -1_i4

    ! ----------------------------------------------------------------------
    ! Find the positions of separators (commas, ",") between the dimensions.
    ! Add a "pseudo separator" at the end to fully automate the consecutive
    ! check for nTimes and nXtrack.
    ! ----------------------------------------------------------------------
    istart = 1  ;  iend = 1  ;  j = 0
    getdim: DO i = 1, swlen 
       IF ( dim_chars(i:i) == ',' ) THEN
          iend = i-1
          IF ( dim_chars(istart:iend) == "nTimes"  )             nt = dim_array(j)
          IF ( dim_chars(istart:iend) == "nXtrack" )             nx = dim_array(j)
          IF ( dim_chars(istart:iend) == "nFitElements"  )       nf = dim_array(j)
          IF ( dim_chars(istart:iend) == "nCharLenFitElements" ) nc = dim_array(j)
          istart = i + 1 ; j = j + 1
          IF ( j > ndim ) EXIT getdim
       END IF
    END DO getdim

    ! ----------------------------------------------------
    ! Assign whatever output variables have been requested
    ! ----------------------------------------------------
    IF ( PRESENT (nTimes_k)       ) nTimes_k       = nt
    IF ( PRESENT (nXtrack_k)      ) nXtrack_k      = nx
    IF ( PRESENT (nFitElem_k)     ) nFitElem_k     = nf
    IF ( PRESENT (nClenFitElem_k) ) nClenFitElem_k = nc

    RETURN
  END SUBROUTINE saopge_l2_read_dimensions

  SUBROUTINE saopge_l2_read_geoloc ( &
       he5stat,                                                       & ! Output
       l2_swath_id, nxs, nxe, nts, nte,                               &
       Lat_k, Lon_k, SZA_k, VZA_k, SAA_k, VAA_k, SpcrAlt_k, TerHgt_k, &
       Time_k, UTC_k, XTQFE_k                                  )

    ! -------------------------------------------------
    ! Read geolocation fields from SAO PGE L2 data file
    ! -------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: l2_swath_id, nxs, nxe, nts, nte

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                                        INTENT (OUT) :: he5stat
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: Lat_k
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: Lon_k
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: SZA_k
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: VZA_k
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: SAA_k
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: VAA_k
    REAL    (KIND=r4), DIMENSION (nts:nte        ), OPTIONAL, INTENT (OUT) :: SpcrAlt_k
    INTEGER (KIND=i2), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: TerHgt_k
    REAL    (KIND=r8), DIMENSION (nts:nte        ), OPTIONAL, INTENT (OUT) :: Time_k  
    INTEGER (KIND=i2), DIMENSION (6,nts:nte      ), OPTIONAL, INTENT (OUT) :: UTC_k
    INTEGER (KIND=i2), DIMENSION (nxs:nxe,nts:nte), OPTIONAL, INTENT (OUT) :: XTQFE_k

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                :: locestat
    INTEGER (KIND=i48)                               :: nxs48, nxe48, nts48, nte48
    INTEGER (KIND=i48)                               :: start_1d, stride_1d, edge_1d
    INTEGER (KIND=i48), DIMENSION (2)                :: start_2d, stride_2d, edge_2d
    REAL    (KIND=r4),  DIMENSION (:,:), ALLOCATABLE :: ESZA_loc
    INTEGER (KIND=i2),  DIMENSION (:,:), ALLOCATABLE :: TerHgt_loc, XTQFE_loc

    he5stat = estat_ok

    ! ------------------------------------------
    ! Convert dimension arguments to INTEGER 4/8
    ! ------------------------------------------
    nxs48 = INT(nxs,KIND=i48) ; nxe48 = INT(nxe,KIND=i48)
    nts48 = INT(nts,KIND=i48) ; nte48 = INT(nte,KIND=i48)


    ! --------
    ! Latitude
    ! --------
    IF ( PRESENT ( Lat_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'Latitude', &
            start_2d, stride_2d, edge_2d, Lat_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------
    ! Longitude
    ! ---------
    IF ( PRESENT ( Lon_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'Longitude', &
            start_2d, stride_2d, edge_2d, Lon_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ------------------
    ! Solar Zenith Angle
    ! ------------------
    IF ( PRESENT ( SZA_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SolarZenithAngle', &
            start_2d, stride_2d, edge_2d, SZA_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! --------------------
    ! Viewing Zenith Angle
    ! --------------------
    IF ( PRESENT ( VZA_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ViewingZenithAngle', &
            start_2d, stride_2d, edge_2d, VZA_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ----------------------
    ! Solar Azimuth Angle
    ! ----------------------
    IF ( PRESENT ( SAA_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SolarAzimuthAngle', &
            start_2d, stride_2d, edge_2d, SAA_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ----------------------
    ! Viewing Azimuth Angle
    ! ----------------------
    IF ( PRESENT ( VAA_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ViewingAzimuthAngle', &
            start_2d, stride_2d, edge_2d, VAA_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -------------------
    ! Spacecraft Altitude
    ! -------------------
    IF ( PRESENT ( SpcrAlt_k ) ) THEN
       start_1d = nts48 ; stride_1d = 1_i48 ; edge_1d = INT(nte,KIND=i48)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SpacecraftAltitude', &
            start_1d, stride_1d, edge_1d, SpcrAlt_k(nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! --------------
    ! Terrain Height
    ! --------------
    IF ( PRESENT ( TerHgt_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       ALLOCATE ( TerHgt_loc(nxs:nxe,nts:nte) )
       locestat = HE5_SWrdfld ( l2_swath_id, 'TerrainHeight', &
            start_2d, stride_2d, edge_2d, TerHgt_loc(nxs:nxe,nts:nte) )
       TerHgt_k(nxs:nxe,nts:nte) = TerHgt_loc(nxs:nxe,nts:nte)
       IF ( ALLOCATED ( TerHgt_loc ) ) DEALLOCATE ( TerHgt_loc ) 
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ----
    ! Time
    ! ----
    IF ( PRESENT ( Time_k ) ) THEN
       start_1d = nts48 ; stride_1d = 1_i48 ; edge_1d = INT(nte,KIND=i48)
       locestat = HE5_SWrdfld ( l2_swath_id, 'Time', &
            start_1d, stride_1d, edge_1d, Time_k(nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------------------------
    ! Xtrack Quality Flags Expanded
    ! -----------------------------
    IF ( PRESENT ( XTQFE_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       ALLOCATE ( XTQFE_loc(nxs:nxe,nts:nte) )
       locestat = HE5_SWrdfld ( l2_swath_id, 'XTrackQualityFlagsExpanded', &
            start_2d, stride_2d, edge_2d, XTQFE_loc(nxs:nxe,nts:nte) )
       XTQFE_k(nxs:nxe,nts:nte) = XTQFE_loc(nxs:nxe,nts:nte)
       IF ( ALLOCATED ( XTQFE_loc ) ) DEALLOCATE ( XTQFE_loc ) 
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    IF ( PRESENT ( UTC_k ) ) THEN
       start_2d  = (/ 0_i48,   nts48 /)
       stride_2d = (/ 1_i48,   1_i48 /)
       edge_2d   = (/ 6_i48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'TimeUTC', &
            start_2d, stride_2d, edge_2d, UTC_k(1:6,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    RETURN
  END SUBROUTINE saopge_l2_read_geoloc

  SUBROUTINE saopge_l2_read_datafields ( &
       he5stat,                                                           & ! Output
       l2_swath_id, nxs, nxe, nts, nte, nfs, nfe, nfecl,                  &
       amfgeo_k, fitcorr_k, fitcorrstr_k, amf_k, amfdiag_k,               &
       amfcfr_k, amfprs_k, fitcol_k, amfalb_k,                            &
       fitcoldstr_k, corrfac_k, fiterr_k, fitrms_k, convflg_k, qaflg_k,   &
       itercnt_k, avgcol_k, avgerr_k, avgrms_k,                           &
       parea_k, pclon_k, pclat_k,                                         &
       o3colt1_k, o3colt2_k, o3colt3_k, o3errt1_k, o3errt2_k, o3errt3_k,  &
       fitcolrs_k)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: l2_swath_id, nxs, nxe, nts, nte
    INTEGER (KIND=i4), INTENT (IN) :: nfs, nfe, nfecl

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                                              INTENT (OUT) :: he5stat
    REAL    (KIND=r4), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amfgeo_k
    REAL    (KIND=r8), DIMENSION(nfs:nfe,nxs:nxe,nts:nte),OPTIONAL, INTENT (OUT) :: fitcorr_k
    INTEGER (KIND=i2), DIMENSION(nfs:nfe,nfecl),          OPTIONAL, INTENT (OUT) :: fitcorrstr_k
    REAL    (KIND=r4), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amf_k
    INTEGER (KIND=i2), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amfdiag_k
    REAL    (KIND=r4), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amfcfr_k
    REAL    (KIND=r4), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amfprs_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: amfalb_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe),                OPTIONAL, INTENT (OUT) :: corrfac_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: fitcol_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: fitcoldstr_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: fiterr_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: fitrms_k
    INTEGER (KIND=i2), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: convflg_k
    INTEGER (KIND=i2), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: qaflg_k
    INTEGER (KIND=i2), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: itercnt_k
    REAL    (KIND=r8),                                    OPTIONAL, INTENT (OUT) :: avgcol_k, avgerr_k, avgrms_k
    REAL    (KIND=r4), DIMENSION(nxs:nxe),                OPTIONAL, INTENT (OUT) :: parea_k
    ! Pixel Corners have dimensions "0:nXtrack,0:nTimes"
    REAL    (KIND=r4), DIMENSION(nxs-1:nxe,nts:nte+1),    OPTIONAL, INTENT (OUT) :: pclon_k, pclat_k
    ! Optional arguments: Present only for O3 prefits
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3colt1_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3colt2_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3colt3_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3errt1_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3errt2_k
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: o3errt3_k

    ! gga to be able to read Reference sector corrected total columns
    REAL    (KIND=r8), DIMENSION(nxs:nxe,nts:nte),        OPTIONAL, INTENT (OUT) :: fitcolrs_k

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                 :: locestat
    INTEGER (KIND=i48)                :: iline, nxs48, nxe48, nts48, nte48, nfs48, nfe48, nfecl48
    INTEGER (KIND=i48)                :: start_1d, stride_1d, edge_1d
    INTEGER (KIND=i48), DIMENSION (2) :: start_2d, stride_2d, edge_2d
    INTEGER (KIND=i48), DIMENSION (3) :: start_3d, stride_3d, edge_3d
    LOGICAL                           :: yn_have_amf, yn_have_alb

    he5stat = estat_ok

    ! ------------------------------------------
    ! Convert dimension arguments to INTEGER 4/8
    ! ------------------------------------------
    nxs48 = INT(nxs,KIND=i48) ; nxe48 = INT(nxe,KIND=i48)
    nts48 = INT(nts,KIND=i48) ; nte48 = INT(nte,KIND=i48)
    nfs48 = INT(nfs,KIND=i48) ; nfe48 = INT(nfe,KIND=i48)
    nfecl48 = INT(nfecl,KIND=i48)

    ! -------------------------
    ! Geometric air mass factor
    ! -------------------------
    IF ( PRESENT ( amfgeo_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'AirMassFactorGeometric', &
            start_2d, stride_2d, edge_2d, amfgeo_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! --------------------
    ! Fitting correlations
    ! --------------------
    IF ( PRESENT (fitcorr_k) ) THEN
       start_3d  = (/ nfs48-1_i48, nxs48-1_i48, nts48 /)
       stride_3d = (/       1_i48,       1_i48, 1_i48 /)
       edge_3d   = (/       nfe48,       nxe48, nte48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'FittingCorrelations', &
            start_3d, stride_3d, edge_3d, fitcorr_k(nfs:nfe,nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -------------------------------------
    ! Elements in fitting correlation array
    ! -------------------------------------
    IF ( PRESENT (fitcorrstr_k) ) THEN
       start_2d  = (/ nfs48-1_i48,   0_i48 /)
       stride_2d = (/       1_i48,   1_i48 /)
       edge_2d   = (/       nfe48, nfecl48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'FittingCorrelationElements', &
            start_2d, stride_2d, edge_2d, fitcorrstr_k(nfs:nfe,1:nfecl) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -------------------------
    ! Molecular air mass factor
    ! -------------------------
    IF ( PRESENT ( amf_k )  ) THEN
       CALL locate_datafield ( l2_swath_id, 'AirMassFactor',  yn_have_amf )
       IF ( yn_have_amf ) THEN
          start_2d  = (/ nxs48-1_i48,       nts48 /)
          stride_2d = (/       1_i48,       1_i48 /)
          edge_2d   = (/       nxe48, nte48+1_i48 /)
          locestat = HE5_SWrdfld ( l2_swath_id, 'AirMassFactor', &
               start_2d, stride_2d, edge_2d, amf_k(nxs:nxe,nts:nte) )
          IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
       ELSE
          amf_k = 0.0_r8
       END IF
    END IF

    ! -------------------------------------
    ! Molecular air mass factor diagnostics
    ! -------------------------------------
    IF ( PRESENT ( amfdiag_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'AirMassFactorDiagnosticFlag', &
            start_2d, stride_2d, edge_2d, amfdiag_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------
    ! Adjusted Scene Albedo
    ! ---------------------
    IF ( PRESENT ( amfalb_k )  ) THEN
       CALL locate_datafield ( l2_swath_id, 'Albedo',  yn_have_alb )
       IF ( yn_have_alb ) THEN
          start_2d  = (/ nxs48-1_i48,       nts48 /)
          stride_2d = (/       1_i48,       1_i48 /)
          edge_2d   = (/       nxe48, nte48+1_i48 /)
          locestat = HE5_SWrdfld ( l2_swath_id, 'Albedo', &
               start_2d, stride_2d, edge_2d, amfalb_k(nxs:nxe,nts:nte) )
          IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
       ELSE
          amfalb_k = 0.0_r8
       END IF
    END IF

    ! ----------------
    ! AMF Cloud fields
    ! ----------------
    IF ( PRESENT ( amfcfr_k )  ) THEN
       CALL locate_datafield ( l2_swath_id, 'AMFCloudFraction',  yn_have_amf )
       IF ( yn_have_amf ) THEN
          start_2d  = (/ nxs48-1_i48,       nts48 /)
          stride_2d = (/       1_i48,       1_i48 /)
          edge_2d   = (/       nxe48, nte48+1_i48 /)
          locestat = HE5_SWrdfld ( l2_swath_id, 'AMFCloudFraction', &
               start_2d, stride_2d, edge_2d, amfcfr_k(nxs:nxe,nts:nte) )
          IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
       ELSE
          amfcfr_k = 0.0_r4
       END IF
    END IF
    IF ( PRESENT ( amfprs_k )  ) THEN
       CALL locate_datafield ( l2_swath_id, 'AMFCloudPressure',  yn_have_amf )
       IF ( yn_have_amf ) THEN
          start_2d  = (/ nxs48-1_i48,       nts48 /)
          stride_2d = (/       1_i48,       1_i48 /)
          edge_2d   = (/       nxe48, nte48+1_i48 /)
          locestat = HE5_SWrdfld ( l2_swath_id, 'AMFCloudPressure', &
               start_2d, stride_2d, edge_2d, amfprs_k(nxs:nxe,nts:nte) )
          IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
       ELSE
          amfprs_k = 0.0_r4
       END IF
    END IF

    ! -------------
    ! Column amount
    ! -------------
    IF ( PRESENT ( fitcol_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ColumnAmount', &
            start_2d, stride_2d, edge_2d, fitcol_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------------------
    ! Column amount destriped
    ! -----------------------
    IF ( PRESENT ( fitcoldstr_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ColumnAmountDestriped', &
            start_2d, stride_2d, edge_2d, fitcoldstr_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ----------------------------
    ! Destriping correction factor
    ! ----------------------------
    IF ( PRESENT ( corrfac_k ) ) THEN
       start_2d  = (/ nxs48-1_i48, 0_i48 /)
       stride_2d = (/       1_i48, 0_i48 /)
       edge_2d   = (/       nxe48, 0_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'CrossTrackStripeCorrection', &
            start_2d(1), stride_2d(1), edge_2d(1), corrfac_k(nxs:nxe) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -------------------------
    ! Column amount uncertainty
    ! -------------------------
    IF ( PRESENT ( fiterr_k ) ) THEN
       start_2d  = (/ nxs48-1_i48, nts48   /)
       stride_2d = (/ 1_i48,               1_i48               /)
       edge_2d   = (/ nxe48,   nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ColumnUncertainty', &
            start_2d, stride_2d, edge_2d, fiterr_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------
    ! Fitting RMS
    ! -----------
    IF ( PRESENT ( fitrms_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'FittingRMS', &
            start_2d, stride_2d, edge_2d, fitrms_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------------------
    ! Fitting convergence flag (ELSUNC)
    ! ---------------------------------
    IF ( PRESENT ( convflg_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'FitConvergenceFlag', &
            start_2d, stride_2d, edge_2d, convflg_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------------------
    ! Main Data Quality Flag
    ! ---------------------------------
    IF ( PRESENT ( qaflg_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'MainDataQualityFlag', &
            start_2d, stride_2d, edge_2d, qaflg_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------------
    ! Convergence iteration count
    ! ---------------------------
    IF ( PRESENT ( itercnt_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'IterationCount', &
            start_2d, stride_2d, edge_2d, itercnt_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------
    ! Average column amount
    ! ---------------------
    IF ( PRESENT ( avgcol_k ) ) THEN
       start_1d = 0_i48 ; stride_1d = 1_i48 ; edge_1d = 1_i48
       locestat = HE5_SWrdfld ( l2_swath_id, 'AverageColumnAmount', &
            start_1d, stride_1d, edge_1d, avgcol_k )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! --------------------------
    ! Average column uncertainty
    ! --------------------------
    IF ( PRESENT ( avgerr_k ) ) THEN
       start_1d = 0_i48 ; stride_1d = 1_i48 ; edge_1d = 1_i48
       locestat = HE5_SWrdfld ( l2_swath_id, 'AverageColumnUncertainty', &
            start_1d, stride_1d, edge_1d, avgerr_k )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -------------------
    ! Average fitting RMS
    ! -------------------
    IF ( PRESENT ( avgrms_k ) ) THEN
       start_1d = 0_i48 ; stride_1d = 1_i48 ; edge_1d = 1_i48
       locestat = HE5_SWrdfld ( l2_swath_id, 'AverageFittingRMS', &
            start_1d, stride_1d, edge_1d, avgrms_k )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------------------
    ! Pixel Area (Tiled)
    ! -----------------------
    IF ( PRESENT ( parea_k ) ) THEN
       start_2d  = (/ nxs48-1_i48, 0_i48 /)
       stride_2d = (/       1_i48, 0_i48 /)
       edge_2d   = (/       nxe48, 0_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'PixelArea', &
            start_2d(1), stride_2d(1), edge_2d(1), parea_k(nxs:nxe) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ----------------------
    ! Pixel Corner Latitudes
    ! ----------------------
    IF ( PRESENT ( pclat_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/ nxe48+1_i48, nte48+2_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'PixelCornerLatitudes', &
            start_2d, stride_2d, edge_2d, pclat_k(nxs-1:nxe,nts:nte+1) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------------------
    ! Pixel Corner Longitudes
    ! -----------------------
    IF ( PRESENT ( pclon_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/ nxe48+1_i48, nte48+2_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'PixelCornerLongitudes', &
            start_2d, stride_2d, edge_2d, pclon_k(nxs-1:nxe,nts:nte+1) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! ---------------------------------------------
    ! Column amount reference sector corrected !gga
    ! ---------------------------------------------
    IF ( PRESENT ( fitcolrs_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'ReferenceSectorCorrectedVerticalColumn', &
            start_2d, stride_2d, edge_2d, fitcolrs_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF


    ! -----------------------------------------
    ! O3 prefit slant columns at 3 temperatures
    ! -----------------------------------------
    IF ( PRESENT ( o3colt1_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnAmountTemperatureT1', &
            start_2d, stride_2d, edge_2d, o3colt1_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF
    IF ( PRESENT ( o3colt2_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnAmountTemperatureT2', &
            start_2d, stride_2d, edge_2d, o3colt2_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF
    IF ( PRESENT ( o3colt3_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnAmountTemperatureT3', &
            start_2d, stride_2d, edge_2d, o3colt3_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    ! -----------------------------------------
    ! O3 prefit slant columns at 3 temperatures
    ! -----------------------------------------
    IF ( PRESENT ( o3errt1_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnUncertaintyTemperatureT1', &
            start_2d, stride_2d, edge_2d, o3errt1_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF
    IF ( PRESENT ( o3errt2_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnUncertaintyTemperatureT2', &
            start_2d, stride_2d, edge_2d, o3errt2_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF
    IF ( PRESENT ( o3errt3_k ) ) THEN
       start_2d  = (/ nxs48-1_i48,       nts48 /)
       stride_2d = (/       1_i48,       1_i48 /)
       edge_2d   = (/       nxe48, nte48+1_i48 /)
       locestat = HE5_SWrdfld ( l2_swath_id, 'SlantColumnUncertaintyTemperatureT3', &
            start_2d, stride_2d, edge_2d, o3errt3_k(nxs:nxe,nts:nte) )
       IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )
    END IF

    RETURN
  END SUBROUTINE saopge_l2_read_datafields

  SUBROUTINE saopge_l2_read_2dr4 ( &
       he5stat, l2_swath_id, fieldname, nxs, nxe, nts, nte, data )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER   (KIND=i4), INTENT (IN) :: l2_swath_id, nxs, nxe, nts, nte
    CHARACTER (LEN=*),   INTENT (IN) :: fieldname

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                             INTENT (OUT) :: he5stat
    REAL    (KIND=r4), DIMENSION(nxs:nxe,nts:nte), INTENT (OUT) :: data

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                 :: locestat
    INTEGER (KIND=i48)                :: nxs48, nxe48, nts48, nte48
    INTEGER (KIND=i48), DIMENSION (2) :: start_2d, stride_2d, edge_2d

    he5stat = estat_ok

    ! ------------------------------------------
    ! Convert dimension arguments to INTEGER 4/8
    ! ------------------------------------------
    nxs48 = INT(nxs,KIND=i48) ; nxe48 = INT(nxe,KIND=i48)
    nts48 = INT(nts,KIND=i48) ; nte48 = INT(nte,KIND=i48)

    ! -------------------------
    ! Geometric air mass factor
    ! -------------------------
    start_2d  = (/ nxs48-1_i48,       nts48 /)
    stride_2d = (/       1_i48,       1_i48 /)
    edge_2d   = (/       nxe48, nte48+1_i48 /)
    locestat = HE5_SWrdfld ( l2_swath_id, TRIM(ADJUSTL(fieldname)), &
         start_2d, stride_2d, edge_2d, data(nxs:nxe,nts:nte) )
    IF ( locestat /= he5_stat_ok ) he5stat = MAX ( he5stat, estat_error )

    RETURN
  END SUBROUTINE saopge_l2_read_2dr4



  SUBROUTINE read_omiclouds ( cloudfile, o3str, nxs, nxe, nts, nte, cfr, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nxs, nxe, nts, nte
    CHARACTER (LEN=*), INTENT (IN) :: cloudfile, o3str

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                              INTENT (OUT) :: errstat
    REAL    (KIND=r4), DIMENSION (nxs:nxe,nts:nte), INTENT (OUT) :: cfr

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER   (KIND=i4)               :: locestat, swath_id, swath_file_id
    REAL      (KIND=r8)               :: scale_cfr, offset_cfr, missval_cfr, scale_ctp, offset_ctp
    INTEGER   (KIND=i2)               :: missval_ctp
    CHARACTER (LEN=maxchlen)          :: swath_name
    INTEGER (KIND=i48)                :: locestat48, swlen48
    INTEGER (KIND=i48)                :: nxs48, nxe48, nts48, nte48
    INTEGER (KIND=i48), DIMENSION (2) :: start_2d, stride_2d, edge_2d

    errstat = estat_ok

    ! ------------------------------------------
    ! Convert dimension arguments to INTEGER 4/8
    ! ------------------------------------------
    nxs48 = INT(nxs,KIND=i48) ; nxe48 = INT(nxe,KIND=i48)
    nts48 = INT(nts,KIND=i48) ; nte48 = INT(nte,KIND=i48)

    ! -------------------------
    ! Attach to OMI Cloud Swath
    ! -------------------------
    swath_file_id = HE5_SWopen     ( TRIM(ADJUSTL(cloudfile)), he5f_acc_rdonly )
    IF ( swath_file_id < 0 ) THEN
       errstat = estat_fatal; RETURN
    END IF
    swath_name = ' '
    locestat48    = HE5_SWinqswath ( TRIM(ADJUSTL(cloudfile)), swath_name, swlen48 )
    swath_id      = HE5_SWattach   ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id < 0 ) THEN
       errstat = estat_fatal
       locestat = HE5_SWclose  ( swath_file_id )
       RETURN
    END IF

    ! -------------------------------------------------------------------------
    ! The first time around we read the scale factors since we don't need to do
    ! this every time we get more scan lines.
    ! -------------------------------------------------------------------------
    scale_cfr = 1.0_r4 ; offset_cfr = 0.0_r4 ; missval_cfr = -9999.0_r4
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudFraction', 'ScaleFactor',  scale_cfr )
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudFraction', 'Offset',       offset_cfr )
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudFraction', 'MissingValue', missval_cfr )

    !scale_ctp = 0.0_r4 ; offset_ctp = 0.0_r4 ; missval_ctp = 0
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudPressure', 'ScaleFactor',  scale_ctp )
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudPressure', 'Offset',       offset_ctp )
    !locestat = HE5_SWrdlattr ( swath_id, 'CloudPressure', 'MissingValue', missval_ctp )

    ! -----------------------
    ! Read current data block
    ! -----------------------
    start_2d  = (/ nxs48-1_i48,        nts48 /)
    stride_2d = (/       1_i48,        1_i48 /)
    edge_2d   = (/       nxe48,  nte48+1_i48 /)

    ! -------------------------------------------
    ! Read cloud fraction and cloud top pressure.
    ! -------------------------------------------
    locestat = HE5_SWrdfld ( swath_id, 'CloudFraction'//TRIM(ADJUSTL(o3str)), start_2d, stride_2d, edge_2d, cfr(nxs:nxe,nts:nte) )
    !locestat = HE5_SWrdfld ( swath_id, 'CloudPressure'//TRIM(ADJUSTL(o3str)), start_2d, stride_2d, edge_2d, ctp(nxs:nxe,nts:nte) )

    ! -------------------
    ! Apply scale factors
    ! -------------------
    cfr = cfr * scale_cfr

    locestat = HE5_SWdetach ( swath_id )
    locestat = HE5_SWclose  ( swath_file_id )

    RETURN
  END SUBROUTINE read_omiclouds

  SUBROUTINE locate_datafield ( l2_swath_id, field_name, yn_have_field )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: l2_swath_id
    CHARACTER (LEN=*), INTENT (IN) :: field_name

    ! ----------------
    ! Output variables
    ! ----------------
    LOGICAL, INTENT (OUT) :: yn_have_field

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i48)           :: nflds48, nstrsize48
    INTEGER (KIND=i4)            :: nflds, nstrsize, sstr, estr, i, tlen
    CHARACTER (LEN=100*maxchlen) :: alias_name, tmpstr

    ! --------------------------
    ! Initialize return variable
    ! --------------------------
    yn_have_field = .FALSE.

    ! ---------------------------------------------------------
    ! Find number of datafields in swath, and return if obvious
    ! ---------------------------------------------------------
    nflds48 = HE5_SWinqdfldalias ( l2_swath_id, alias_name, nstrsize48 )
    nflds    = INT ( nflds48,    KIND=i4 )
    nstrsize = INT ( nstrsize48, KIND=i4 )
    
    IF ( ( nflds < 1 ) .OR. ( nstrsize < LEN_TRIM(ADJUSTL(field_name)) ) ) RETURN

    sstr = 1 ; estr = nstrsize
    get_name: DO i = 1, nflds-1
       tmpstr = ""
       tmpstr = TRIM(ADJUSTL(alias_name(sstr:nstrsize))) ; tlen = nstrsize - sstr + 1
       estr = INDEX ( tmpstr(1:tlen), ',' )       
       IF ( TRIM(ADJUSTL(tmpstr(1:estr-1))) == TRIM(ADJUSTL(field_name)) ) THEN
          yn_have_field = .TRUE.
          EXIT get_name
       END IF
       sstr = sstr + estr
    END DO get_name

    RETURN
  END SUBROUTINE locate_datafield

  SUBROUTINE read_2dr4_field (  swid, fname, n1s, n1e, n2s, n2e, arr, estat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: swid, n1s, n1e, n2s, n2e
    CHARACTER (LEN=*), INTENT (IN) :: fname

    ! -----------------------------
    ! Output and Modified Variables
    ! -----------------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: estat
    REAL (KIND=r4), DIMENSION (n1s:n1e,n2s:n2e), INTENT (OUT) :: arr

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4)                 :: locestat
    INTEGER (KIND=i48)                :: j1, j2, k1, k2
    INTEGER (KIND=i48), DIMENSION (2) :: start_2d, stride_2d, edge_2d


    ! ------------------------------------------
    ! Convert dimension arguments to INTEGER 4/8
    ! ------------------------------------------
    j1 = INT(n1s,KIND=i48) ; j2 = INT(n1e,KIND=i48)
    k1 = INT(n2s,KIND=i48) ; k2 = INT(n2e,KIND=i48)

    start_2d  = (/ j1-1_i48,       j2 /)
    stride_2d = (/    1_i48,    1_i48 /)
    edge_2d   = (/       k1, k2+1_i48 /)
    locestat = HE5_SWrdfld ( swid, TRIM(ADJUSTL(fname)), &
         start_2d, stride_2d, edge_2d, arr(j1:j2,k1:k2) )
    IF ( locestat /= he5_stat_ok ) estat = MAX ( estat, estat_error )

    RETURN
  END SUBROUTINE read_2dr4_field

END MODULE SAO_OMIL2_ReadLib_he5_module
