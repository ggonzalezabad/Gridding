MODULE SAO_OMIL2_CreateHDF5_module

  USE SAO_OMIL2_ReadLib_basic_module, ONLY: i2, i4, i48, r4, r8, maxchlen

  IMPLICIT NONE
  INCLUDE 'hdfeos5.inc'

  INTEGER (KIND=i4), EXTERNAL, PRIVATE ::                                    &
       he5_swattach, he5_swclose, he5_swcreate, he5_swdefdfld, he5_swdefdim, &
       he5_swdetach,  he5_swopen, he5_swwrfld,  he5_swwrlattr, he5_swdefcomch, &
       he5_swdefgfld, he5_swwrattr


  ! -------------------------------------
  ! Names with dimensions and data fields
  ! -------------------------------------
  CHARACTER (LEN=10), PARAMETER :: nlondim  = 'nLongitude', gridlons = 'Longitudes'
  CHARACTER (LEN= 9), PARAMETER :: nlatdim  = 'nLatitude',  gridlats = 'Latitudes'
  CHARACTER (LEN=20), PARAMETER :: nlatlon  = nlondim//','//nlatdim


  ! ------------------------------------------
  ! Variables that hold Swath File information
  ! ------------------------------------------
  CHARACTER (LEN=maxchlen)  :: swath_name     ! The Swath Name
  INTEGER                   :: swath_file_id  ! ID number for HE5 output file
  INTEGER                   :: swath_id       ! ID number for swath
  INTEGER, PRIVATE          :: he5stat        ! General status variable

  INTEGER (KIND=i48)               :: he5_start_1, he5_stride_1, he5_edge_1
  INTEGER (KIND=i48), DIMENSION(2) :: he5_start_2, he5_stride_2, he5_edge_2

  ! --------------------------------------------------------
  ! Variables associated with field compression and chunking
  ! --------------------------------------------------------
  INTEGER, PARAMETER                         :: ncchunkdim = 2
  INTEGER, PARAMETER                         :: comp_par = 9, comp_type = HE5_HDFE_COMP_SHUF_DEFLATE
  INTEGER (KIND=i48), DIMENSION (ncchunkdim) :: cchunk_dim


  ! ------------------------------------
  ! Number of significant digits to keep
  ! ------------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_roff_dig = 6

  ! -------------
  ! Missing Value
  ! -------------
  REAL (KIND=r4), PARAMETER :: r8MissVal_r4 = -1.0E+30_r4
  REAL (KIND=r8), PARAMETER :: r8MissVal_r8 = -1.0E+30_r8

CONTAINS

  SUBROUTINE create_he5_file ( &
       he5file, nlon, nlat, lons, lats, dlon, dlat, cld_frc_min, cld_frc_max, qflg_max, &
        xtrack_min, xtrack_max, cld_prs_min_k, cld_prs_max_k                            )

    !------------------------------------------------------------------------------
    ! This subroutine creates an HE5 file with gridded averages of SAO PGEs
    !
    ! Input:
    !   he5file     - name of HE5 output file
    !   nlon        - number of longitudes
    !   nlat        - number of latitudes
    !   lons        - grid longitudes
    !   lats        - grid latitudes
    !   dlon        - longitude grid spacing
    !   dlat        - latitude grid spacing
    !   cld_frc_min - minimum pixel cloud fraction
    !   cld_frc_max - maximum pixel cloud fraction
    !   qflg_max    - maximum data quality flag
    !   cld_prs_min - minimum pixel cloud fraction, OPTIONAL
    !   cld_prs_max - maximum pixel cloud fraction, OPTIONAL
    !   xtrack_min  - minimum xtrack pixel number
    !   xtrack_max  - maximum xtrack pixel number
    !------------------------------------------------------------------------------

    IMPLICIT NONE


    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                  INTENT (IN) :: he5file
    INTEGER (KIND=i4),                  INTENT (IN) :: nlon, nlat
    INTEGER (KIND=i2),                  INTENT (IN) :: qflg_max, xtrack_min, xtrack_max
    REAL    (KIND=r4),                  INTENT (IN) :: dlon, dlat, cld_frc_min, cld_frc_max
    REAL    (KIND=r4), DIMENSION(nlon), INTENT (IN) :: lons
    REAL    (KIND=r4), DIMENSION(nlat), INTENT (IN) :: lats
    REAL    (KIND=r4), OPTIONAL,        INTENT (IN) :: cld_prs_min_k, cld_prs_max_k
    
    he5stat = 0

    ! ---------------------------------------------------------------
    ! Open HE5 output file and check AMF_SWATH_FILE_ID ( -1 if error)
    ! ---------------------------------------------------------------
    swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(he5file)), he5f_acc_trunc )
    IF ( swath_file_id == -1 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWopen failed!'; STOP 1
    END IF

    ! ------------------------------------------------------
    ! Create HE5 swath and check SWATH_ID ( -1 if error)
    ! ------------------------------------------------------
    swath_id = HE5_SWcreate ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == -1 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWcreate failed!'; STOP 1
    END IF

    ! ----------------------------------
    ! Define new dimensions in HE5 swath
    ! ----------------------------------
    he5stat = HE5_SWdefdim  ( swath_id, nlondim, INT(nlon, KIND=i48) )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdefdim failed for swath dimensions!'; STOP 1
    END IF
    he5stat = HE5_SWdefdim  ( swath_id, nlatdim, INT(nlat, KIND=i48) )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdefdim failed for swath dimensions!'; STOP 1
    END IF

    ! -----------------------------------------------
    ! Define array value data fields in HE5 swath
    ! -----------------------------------------------
    he5stat = HE5_SWdefgfld ( swath_id, gridlons, nlondim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )
    he5stat = HE5_SWdefgfld ( swath_id, gridlats, nlatdim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdefgfld failed for Array Data Fields!'; STOP 1
    END IF


    ! -------------------------------------------------------------------
    ! Set the chunking dimension for higher-dimensional array compression
    ! -------------------------------------------------------------------
    cchunk_dim(1:ncchunkdim) = (/ INT(nlon,KIND=i48), INT(nlat,KIND=i48) /)
    he5stat = HE5_SWdefcomch ( &
         swath_id, comp_type, comp_par, ncchunkdim, cchunk_dim(1:ncchunkdim) )


    ! ----------------------------------------------
    ! Compose the data field names and define fields
    ! ----------------------------------------------
    he5stat = HE5_SWdefdfld (swath_id,'AMF',                   nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'CloudFraction',         nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'CloudPressure',         nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'ColumnFitted',          nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'ColumnFittedCorrected', nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'SlantFitted',           nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'SlantFittedCorrected',  nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'ColumnUncertainty',     nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'SlantUncertainty',      nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'GridArea',              nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'NumberOfSamples',       nlatlon," ",HE5T_Native_INT,  he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'QualityFlag',           nlatlon," ",HE5T_Native_INT16,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'SurfaceAlbedo',         nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )
    he5stat = HE5_SWdefdfld (swath_id,'SurfaceAltitude',       nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )   
    he5stat = HE5_SWdefdfld (swath_id,'RMS          ',         nlatlon," ",HE5T_Native_FLOAT,he5_hdfe_nomerge )

    ! -------------------------------------------------------------------------------
    ! Detach from and re-attach to created swath (recommended before adding to swath)
    ! -------------------------------------------------------------------------------
    he5stat  = HE5_SWdetach ( swath_id )
    swath_id = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
    IF ( swath_id == -1 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWattach failed!'; STOP 1
    END IF

    ! ---------------------------------
    ! Write Global Swath Attributes:
    !   (1) Minimum Cloud Fraction
    !   (2) Maximum Cloud Fraction
    !   (3) Maximum Data Quality Flag
    ! ---------------------------------
    he5stat = he5_swwrattr ( swath_id, "Minimum Cloud Fraction",    HE5T_NATIVE_FLOAT, 1_i48, cld_frc_min )
    he5stat = he5_swwrattr ( swath_id, "Maximum Cloud Fraction",    HE5T_NATIVE_FLOAT, 1_i48, cld_frc_max )
    he5stat = he5_swwrattr ( swath_id, "Maximum Data Quality Flag", HE5T_NATIVE_INT16, 1_i48, qflg_max    )
    he5stat = he5_swwrattr ( swath_id, "Minimum xtrack pixel",      HE5T_NATIVE_INT16, 1_i48, xtrack_min  )
    he5stat = he5_swwrattr ( swath_id, "Maximum xtrack pixel",      HE5T_NATIVE_INT16, 1_i48, xtrack_max  )
    IF ( PRESENT ( cld_prs_min_k ) ) &
         he5stat = he5_swwrattr ( swath_id, "Minimum Cloud Top Pressure", HE5T_NATIVE_FLOAT, 1_i48, cld_prs_min_k )
    IF ( PRESENT ( cld_prs_max_k ) ) &
         he5stat = he5_swwrattr ( swath_id, "Maximum Cloud Top Pressure", HE5T_NATIVE_FLOAT, 1_i48, cld_prs_max_k )
         
    ! -------------------------------
    ! Write dimension-defining arrays
    ! -------------------------------
    he5_start_1 = 0_i48 ; he5_stride_1 = 1_i48 ; he5_edge_1 = INT(nlon,KIND=i48)
    he5stat = HE5_SWwrfld   (swath_id,gridlons,he5_start_1,he5_stride_1,he5_edge_1,lons(1:nlon))
    he5stat = HE5_SWwrlattr (swath_id,gridlons, "GridSpacing", HE5T_NATIVE_FLOAT, 1_i48, dlon )

    he5_start_1 = 0_i48 ; he5_stride_1 = 1_i48 ; he5_edge_1 = INT(nlat,KIND=i48)
    he5stat = HE5_SWwrfld   (swath_id,gridlats,he5_start_1,he5_stride_1,he5_edge_1,lats(1:nlat))
    he5stat = HE5_SWwrlattr (swath_id,gridlats, "GridSpacing", HE5T_NATIVE_FLOAT, 1_i48, dlat )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdwrfld failed for grid latitudes and longitudes!'; STOP 1
    END IF


    RETURN
  END SUBROUTINE create_he5_file

  SUBROUTINE write_he5_data_r4 ( fieldname, nlon, nlat, grid_2d, scale_fac )

    !------------------------------------------------------------------------------
    ! This subroutine write r4 values to the HE5 file
    !
    ! Input:
    !   fieldname  - string with the field name to write to
    !   nlon       - number of longitudes
    !   nlat       - number of latitudes
    !   grid_2d    - gridded data to write to file
    !   scale_fac  - scaling factor for GRID_2D
    !------------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                         INTENT (IN)    :: fieldname
    INTEGER,                                   INTENT (IN)    :: nlon, nlat
    REAL     (KIND=r4),                        INTENT (IN)    :: scale_fac
    REAL     (KIND=r4), DIMENSION (nlon,nlat), INTENT (INOUT) :: grid_2d

    ! ------------------------------------------------
    ! Skip writing if we don't have a valid Field Name
    ! ------------------------------------------------
    IF ( LEN_TRIM(ADJUSTL(fieldname)) <= 0 ) RETURN

    ! -------------------------
    ! Round-Off of gridded data
    ! -------------------------
    CALL roundoff_2darr_r4 ( n_roff_dig, nlon, nlat, r8MissVal_r4, grid_2d(1:nlon,1:nlat) )

    ! ------------------------
    ! Write the gridded arrays
    ! ------------------------
    he5_start_2  = (/ 0_i48,              0_i48 /)
    he5_stride_2 = (/ 1_i48,              1_i48 /)
    he5_edge_2   = (/ INT(nlon,KIND=i48), INT(nlat,KIND=i48) /)

    he5stat = HE5_SWwrfld   ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         he5_start_2, he5_stride_2, he5_edge_2, grid_2d(1:nlon,1:nlat)  )
    he5stat = HE5_SWwrlattr ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         "ScaleFactor", HE5T_NATIVE_FLOAT, 1_I48, scale_fac                 )
    he5stat = HE5_SWwrlattr ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         "MissingValue", HE5T_NATIVE_FLOAT, 1_i48, r8MissVal_r4             )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdwrfld failed for '//TRIM(ADJUSTL(fieldname))
       STOP 1
    END IF


    RETURN
  END SUBROUTINE write_he5_data_r4

  SUBROUTINE write_he5_data_i4 ( fieldname, nlon, nlat, grid_2d, scale_fac )

    !------------------------------------------------------------------------------
    ! This subroutine write i4 values to the HE5 file
    !
    ! Input:
    !   fieldname  - string with the field name to write to
    !   nlon       - number of longitudes
    !   nlat       - number of latitudes
    !   grid_2d    - gridded data to write to file
    !   scale_fac  - scaling factor for GRID_2D
    !------------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                          INTENT (IN) :: fieldname
    INTEGER,                                    INTENT (IN) :: nlon, nlat
    INTEGER  (KIND=i4),  DIMENSION (nlon,nlat), INTENT (IN) :: grid_2d
    REAL     (KIND=r4),                         INTENT (IN) :: scale_fac

    ! ------------------------------------------------
    ! Skip writing if we don't have a valid Field Name
    ! ------------------------------------------------
    IF ( LEN_TRIM(ADJUSTL(fieldname)) <= 0 ) RETURN

    ! ------------------------
    ! Write the gridded arrays
    ! ------------------------
    he5_start_2  = (/ 0_I48,              0_i48              /)
    he5_stride_2 = (/ 1_I48,              1_i48              /)
    he5_edge_2   = (/ INT(nlon,KIND=i48), INT(nlat,KIND=i48) /)

    he5stat = HE5_SWwrfld   ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         he5_start_2, he5_stride_2, he5_edge_2, grid_2d(1:nlon,1:nlat)  )
    he5stat = HE5_SWwrlattr ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         "ScaleFactor", HE5T_NATIVE_FLOAT, 1_i48, scale_fac                 )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdwrfld failed for '//TRIM(ADJUSTL(fieldname))
       STOP 1
    END IF

    RETURN
  END SUBROUTINE write_he5_data_i4

  SUBROUTINE write_he5_data_i2 ( fieldname, nlon, nlat, grid_2d, scale_fac )

    !------------------------------------------------------------------------------
    ! This subroutine write i2 values to the HE5 file
    !
    ! Input:
    !   fieldname  - string with the field name to write to
    !   nlon       - number of longitudes
    !   nlat       - number of latitudes
    !   grid_2d    - gridded data to write to file
    !   scale_fac  - scaling factor for GRID_2D
    !------------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*),                          INTENT (IN) :: fieldname
    INTEGER,                                    INTENT (IN) :: nlon, nlat
    INTEGER  (KIND=i2),  DIMENSION (nlon,nlat), INTENT (IN) :: grid_2d
    REAL     (KIND=r4),                         INTENT (IN) :: scale_fac

    ! ------------------------------------------------
    ! Skip writing if we don't have a valid Field Name
    ! ------------------------------------------------
    IF ( LEN_TRIM(ADJUSTL(fieldname)) <= 0 ) RETURN

    ! ------------------------
    ! Write the gridded arrays
    ! ------------------------
    he5_start_2  = (/ 0_I48,              0_i48              /)
    he5_stride_2 = (/ 1_I48,              1_i48              /)
    he5_edge_2   = (/ INT(nlon,KIND=i48), INT(nlat,KIND=i48) /)

    he5stat = HE5_SWwrfld   ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         he5_start_2, he5_stride_2, he5_edge_2, grid_2d(1:nlon,1:nlat)  )
    he5stat = HE5_SWwrlattr ( swath_id, TRIM(ADJUSTL(fieldname)),     &
         "ScaleFactor", HE5T_NATIVE_FLOAT, 1_i48, scale_fac                 )
    IF ( he5stat /= 0 ) THEN
       WRITE (*,*) 'ERROR: HE5_SWdwrfld failed for '//TRIM(ADJUSTL(fieldname))
       STOP 1
    END IF

    RETURN
  END SUBROUTINE write_he5_data_i2

  SUBROUTINE close_he5_output_file ( )

    IMPLICIT NONE

    ! -----------------------------------------------
    ! Detach from HE5 swath and close HE5 output file
    ! -----------------------------------------------
    he5stat = HE5_SWdetach ( swath_id )
    he5stat = HE5_SWclose  ( swath_file_id )
    IF ( he5stat /= 0 ) WRITE (*,*) 'WARNING: HE5_SWdetach/HE5_SWclose failed!'

    RETURN
  END SUBROUTINE close_he5_output_file

  SUBROUTINE roundoff_r8 ( ndecim, r8value )

    IMPLICIT NONE

    ! ---------------------------------------------------------------------
    ! Explanation of SUBROUTINE arguments:
    !
    !    ndecim ........... number of significant digits to keep
    !    r8value .......... DOUBLE PRECISION number to be truncated
    ! ---------------------------------------------------------------------

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN)    :: ndecim
    REAL    (KIND=r8), INTENT (INOUT) :: r8value

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: pow
    REAL    (KIND=r8) :: tmpval


    ! ------------------------------------------
    ! Nothing to be done if we have a Zero value
    ! ------------------------------------------
    IF (  r8value == 0.0_r8 ) RETURN

    ! ---------------------------------
    ! Power of 10 of the original value
    ! ---------------------------------
    pow = INT ( LOG10 (ABS(r8value)), KIND=i4 )

    ! ------------------------------------------------------
    ! Remove original power of 10 and shift by NDECIM powers
    ! ------------------------------------------------------
    tmpval = r8value * 10.0_r8**(ndecim-pow)

    ! ---------------------------------------------
    ! Find the nearest INTEGER and undo power-shift
    ! ---------------------------------------------
    r8value = ANINT ( tmpval )  * 10.0_r8**(pow-ndecim)

    RETURN
  END SUBROUTINE roundoff_r8

  SUBROUTINE roundoff_2darr_r8 ( ndecim, n1, n2, missval, r8value )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
    REAL    (KIND=r8),                    INTENT (IN)    :: missval
    REAL    (KIND=r8), DIMENSION (n1,n2), INTENT (INOUT) :: r8value
  
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: pow, i, j
    REAL    (KIND=r8) :: tmpval, r8val
  
    DO i = 1, n1
       DO j = 1, n2
          r8val = r8value(i,j)
          IF ( r8val == 0.0_r8 .OR. r8val <= missval ) CYCLE

          ! ---------------------------------
          ! Power of 10 of the original value
          ! ---------------------------------
          pow = INT ( LOG10 (ABS(r8val)), KIND=i4 )

          ! ------------------------------------------------------
          ! Remove original power of 10 and shift by NDECIM powers
          ! ------------------------------------------------------
          tmpval = r8val * 10.0_r8**(ndecim-pow)

          ! ---------------------------------------------
          ! Find the nearest INTEGER and undo power-shift
          ! ---------------------------------------------
          r8value(i,j) = ANINT ( tmpval )  * 10.0_r8**(pow-ndecim)
  
       END DO
    END DO
  
    RETURN
  END SUBROUTINE roundoff_2darr_r8
  

  SUBROUTINE roundoff_2darr_r4 ( ndecim, n1, n2, missval, r4value )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
    REAL    (KIND=r4),                    INTENT (IN)    :: missval
    REAL    (KIND=r4), DIMENSION (n1,n2), INTENT (INOUT) :: r4value
  
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: pow, i, j
    REAL    (KIND=r4) :: tmpval, r4val
  
    DO i = 1, n1
       DO j = 1, n2

          r4val = r4value(i,j)
          IF ( r4val == 0.0_r4 .OR. r4val <= missval ) CYCLE

          ! ---------------------------------
          ! Power of 10 of the original value
          ! ---------------------------------
          pow = INT ( LOG10 (ABS(r4val)), KIND=i4 )

          ! ------------------------------------------------------
          ! Remove original power of 10 and shift by NDECIM powers
          ! ------------------------------------------------------
          tmpval = r4val * 10.0_r4**(ndecim-pow)

          ! ---------------------------------------------
          ! Find the nearest INTEGER and undo power-shift
          ! ---------------------------------------------
          r4value(i,j) = ANINT ( tmpval )  * 10.0_r4**(pow-ndecim)
  
       END DO
    END DO
  
    RETURN
  END SUBROUTINE roundoff_2darr_r4


END MODULE SAO_OMIL2_CreateHDF5_module
