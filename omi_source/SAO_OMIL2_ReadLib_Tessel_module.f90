MODULE SAO_OMIL2_ReadLib_Tessel_module

  USE SAO_OMIL2_ReadLib_Basic_module
  USE SAO_OMIL2_ReadLib_SphGeo_module

  IMPLICIT NONE

  ! ----------------------------
  ! Parameter values and indices
  ! ----------------------------
  ! Small adjustment to make pixel corners conform to tessellation requirements
  REAL    (KIND=r8), PARAMETER     :: small_dpix = 1.0E-5_r8
  ! Absolute latitude value above which to perform polar and dateline checking
  REAL    (KIND=r8), PARAMETER     :: tess_polar_lat = 0.0_r8



CONTAINS

  SUBROUTINE pmt_master ( &
       pix_lon, pix_lat, nlon_grid, nlat_grid, lon_grid, lat_grid, dlon_grid, dlat_grid,    &
       idx_lon, idx_lat, ilon_min, ilon_max, ilat_min, ilat_max, nfac_fg,                   &
       tessarea, tess_sum, yn_good_pix )

    IMPLICIT NONE

    
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlon_grid, nlat_grid, nfac_fg
    INTEGER (KIND=i4),                        INTENT (IN) :: ilon_min, ilon_max, ilat_min, ilat_max
    REAL    (KIND=r8),                        INTENT (IN) :: dlon_grid, dlat_grid
    REAL    (KIND=r8), DIMENSION (4),         INTENT (IN) :: pix_lon, pix_lat
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (IN) :: lon_grid
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid
    INTEGER (KIND=i4), DIMENSION (4),         INTENT (IN) :: idx_lat, idx_lon


    ! ----------------
    ! Output Variables
    ! ----------------
    LOGICAL,                                                            INTENT (OUT) :: yn_good_pix
    REAL    (KIND=r8),                                                  INTENT (OUT) :: tess_sum
    REAL    (KIND=r8), DIMENSION (ilon_min:ilon_max,ilat_min:ilat_max), INTENT (OUT) :: tessarea

    ! ---------------
    ! Local Variables
    ! ---------------
    INTEGER (KIND=i4)                :: estat, i, ilat, ilon, jlat, jlon
    INTEGER (KIND=i4)                :: nlon_tess_fg, nlat_tess_fg, nlon_fg, nlat_fg
    INTEGER (KIND=i4)                :: ilon_min_fg, ilon_max_fg, ilat_min_fg, ilat_max_fg
    REAL    (KIND=r8)                :: box_area
    INTEGER (KIND=i4), DIMENSION (4) :: idx_lat_fg, idx_lon_fg

    REAL    (KIND=r8), DIMENSION (:),   ALLOCATABLE :: lon_fg, lat_fg
    REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: latlon_isec_fg
    INTEGER (KIND=i4), DIMENSION (:,:), ALLOCATABLE :: tess_latlon_idx_fg

    REAL    (KIND=r8),  PARAMETER  :: r8missval = -1.0E+30_r8

    yn_good_pix     = .FALSE.

    ! ----------------------
    ! The fine sampling grid
    ! ----------------------
    nlon_fg = nfac_fg*(ilon_max-ilon_min+1)+1  ;  nlat_fg = nfac_fg*(ilat_max-ilat_min+1)+1
    ALLOCATE ( lon_fg(nlon_fg), STAT = estat )  ;  IF ( estat /= 0 ) GO TO 666
    ALLOCATE ( lat_fg(nlat_fg), STAT = estat )  ;  IF ( estat /= 0 ) GO TO 666

    lon_fg(1:nlon_fg) = lon_grid(ilon_min) + (/ (REAL(i,KIND=r8)*dlon_grid/REAL(nfac_fg,KIND=r8), i = 0, nlon_fg-1) /)
    lat_fg(1:nlat_fg) = lat_grid(ilat_min) + (/ (REAL(i,KIND=r8)*dlat_grid/REAL(nfac_fg,KIND=r8), i = 0, nlat_fg-1) /)
    

    CALL pmt_find_indices_for_tessellation (                    &
         pix_lon, pix_lat, nlon_fg, lon_fg(1:nlon_fg), nlat_fg, lat_fg(1:nlat_fg),  &
         idx_lon_fg, idx_lat_fg, ilon_min_fg, ilon_max_fg, ilat_min_fg, ilat_max_fg, nlon_tess_fg, nlat_tess_fg   )
     
    ALLOCATE ( tess_latlon_idx_fg (1:nlon_fg,1:nlat_fg), STAT=estat )  ;  IF ( estat /= 0 ) GO TO 666
    ALLOCATE (     latlon_isec_fg (1:nlat_fg,1:3),       STAT=estat )  ;  IF ( estat /= 0 ) GO TO 666

    tess_latlon_idx_fg = 0 ; latlon_isec_fg = r8missval

    CALL pmt_find_grid_intersections (                     &
         nlat_fg, lat_fg(1:nlat_fg), pix_lon, pix_lat, idx_lon_fg, idx_lat_fg, ilat_min_fg, ilat_max_fg, &
         latlon_isec_fg(ilat_min_fg:ilat_max_fg,1:3) )

    DO ilat = ilat_min_fg, ilat_max_fg
       WHERE ( &
            lon_fg(1:nlon_fg) >= latlon_isec_fg(ilat,2) .AND. &
            lon_fg(1:nlon_fg) <  latlon_isec_fg(ilat,3)         )
          tess_latlon_idx_fg(1:nlon_fg,ilat) = 1
       ENDWHERE
    END DO

    DO ilat = ilat_min, ilat_max
       jlat = nfac_fg*(ilat-ilat_min) + 1
       DO ilon = ilon_min, ilon_max
          jlon = nfac_fg*(ilon-ilon_min) + 1
          tessarea       (ilon,ilat) = REAL ( SUM(tess_latlon_idx_fg(jlon:jlon+nfac_fg-1,jlat:jlat+nfac_fg-1)), KIND=r8 )
       END DO
    END DO

    box_area = dlon_grid * dlat_grid / REAL(nfac_fg*nfac_fg, KIND=r8)

    tess_sum = SUM ( tess_latlon_idx_fg ) * box_area
    tessarea = tessarea * box_area

    yn_good_pix = .TRUE.

666 CONTINUE

    IF ( ALLOCATED ( lon_fg             ) ) DEALLOCATE ( lon_fg             )
    IF ( ALLOCATED ( lat_fg             ) ) DEALLOCATE ( lat_fg             )
    IF ( ALLOCATED ( tess_latlon_idx_fg ) ) DEALLOCATE ( tess_latlon_idx_fg )
    IF ( ALLOCATED ( latlon_isec_fg     ) ) DEALLOCATE ( latlon_isec_fg     )


    RETURN
  END SUBROUTINE pmt_master




  SUBROUTINE pmt_find_indices_for_tessellation (                     &
       pix_lon, pix_lat, nlon_grid, lon_grid, nlat_grid, lat_grid,   &
       idx_lon, idx_lat, ilon_min, ilon_max, ilat_min, ilat_max, nlon_tess, nlat_tess   )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlon_grid, nlat_grid
    REAL    (KIND=r8), DIMENSION (4),         INTENT (IN) :: pix_lon, pix_lat
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (IN) :: lon_grid
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: ilon_min, ilon_max, ilat_min, ilat_max
    INTEGER (KIND=i4), INTENT (OUT) :: nlon_tess, nlat_tess
    INTEGER (KIND=i4), DIMENSION (4), INTENT (OUT) :: idx_lon, idx_lat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                :: i

    DO i = 1, 4
       IF ( pix_lon(i) >= MAXVAL(lon_grid(1:nlon_grid)) ) THEN
          idx_lon(i) = nlon_grid
       ELSE IF ( pix_lon(i) <= MINVAL(lon_grid(1:nlon_grid)) ) THEN
          idx_lon(i) = 1
       ELSE
          idx_lon(i) = &
               MAXVAL( MAXLOC( lon_grid(1:nlon_grid), MASK=(lon_grid(1:nlon_grid) <= pix_lon(i) ) ) )
       END IF
       IF ( pix_lat(i) >= MAXVAL(lat_grid(1:nlat_grid)) ) THEN
          idx_lat(i) = nlat_grid
       ELSE IF ( pix_lat(i) <= MINVAL(lat_grid(1:nlat_grid)) ) THEN
          idx_lat(i) = 1
       ELSE
          idx_lat(i) = &
               MAXVAL( MAXLOC( lat_grid(1:nlat_grid), MASK=(lat_grid(1:nlat_grid) <= pix_lat(i) ) ) )
       END IF
    END DO

    IF ( ANY ( idx_lon == 0 ) ) THEN
       PRINT *, 'BAD Lons'
       PRINT *, idx_lon
       PRINT *, pix_lon
       PRINT *, nlon_grid, lon_grid(1), lon_grid(nlon_grid)
    END IF
    IF ( ANY ( idx_lat == 0 ) ) THEN
       PRINT *, 'BAD Lats'
       PRINT *, idx_lat
       PRINT *, pix_lat
       PRINT *, nlat_grid, lat_grid(1), lat_grid(nlat_grid)
    END IF

    ilon_min = MINVAL(idx_lon) ; ilon_max = MAXVAL(idx_lon)
    ilat_min = MINVAL(idx_lat) ; ilat_max = MAXVAL(idx_lat)

    nlon_tess = ilon_max - ilon_min + 1
    nlat_tess = ilat_max - ilat_min + 1

    RETURN
  END SUBROUTINE pmt_find_indices_for_tessellation


  SUBROUTINE pmt_find_grid_intersections (                     &
       nlat_grid, lat_grid, pix_lon, pix_lat, idx_lon, idx_lat, ilat_min, ilat_max, &
       latlon_isec )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlat_grid, ilat_min, ilat_max
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid
    REAL    (KIND=r8), DIMENSION (4),         INTENT (IN) :: pix_lon, pix_lat
    INTEGER (KIND=i4), DIMENSION (4),         INTENT (IN) :: idx_lon, idx_lat

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), DIMENSION (ilat_min:ilat_max,3), INTENT (OUT) :: latlon_isec


    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: j, j1, j2, pA, pB, ilons, ilone, ilats, ilate, ilat, ilon
    REAL    (KIND=r8) :: a,  b,  c,  alp,  bet,  gam
    REAL    (KIND=r8) :: a0, b0, c0, alp0, bet0, gam0
    REAL      (KIND=r8),  PARAMETER                   :: r8missval = -1.0E+30_r8


    latlon_isec = r8missval
    
    DO j = 1, 4

       j1 = j ; j2 = j + 1 ; IF ( j2 > 4 ) j2 = 1
        
       IF ( pix_lat(j1) > pix_lat(j2) ) THEN
          pA = j1 ; pB = j2 
       ELSE
          pA = j2 ; pB = j1
       END IF

       ilons = MIN(idx_lon(pA),idx_lon(pB)) ; ilone = MAX(idx_lon(pA),idx_lon(pB))
       ilats = MIN(idx_lat(pA),idx_lat(pB)) ; ilate = MAX(idx_lat(pA),idx_lat(pB))

       ! --------------------------------------------------------------------
       ! The base quantities for the triangle defined by the first two points
       ! --------------------------------------------------------------------
       b0   = pihalf - pix_lat(pA)*deg2rad 
       a0   = pihalf - pix_lat(pB)*deg2rad
       gam0 = ABS(pix_lon(pA) - pix_lon(pB))*deg2rad
       CALL sphergeom_sas ( a0, b0, gam0, c0, alp0, bet0 )

       DO ilat = ilats+1, ilate
          a   = pihalf - lat_grid(ilat)*deg2rad
          CALL sphergeom_ssa ( a, b0, alp0, c, bet, gam )
          IF ( gam < 0.0_r8 ) gam = gam + pi
          gam = angle_minus_twopi ( gam, pi ) * rad2deg
          gam = pix_lon(pA)-SIGN(gam,pix_lon(pA)-pix_lon(pB))

          ! -----------------------------------------------------------------
          ! Save the results in an array. There can only be two intersections
          ! of each latitude line by pixel boundaries.
          ! -----------------------------------------------------------------
          latlon_isec(ilat,1) = lat_grid(ilat)
          IF ( latlon_isec(ilat,2) == r8missval ) THEN
             latlon_isec(ilat,2) = gam
          ELSE
             ! --------------------------------------------
             ! Take care of the correct ordering right here
             ! --------------------------------------------
             IF ( latlon_isec(ilat,2) <= gam ) THEN
                latlon_isec(ilat,3) = gam
             ELSE
                latlon_isec(ilat,3) = latlon_isec(ilat,2)
                latlon_isec(ilat,2) = gam
             END IF
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE pmt_find_grid_intersections


  SUBROUTINE pmt_prepare_omi_tessell (                                         &
       pix_lon, pix_lat, nlon_grid, nlat_grid, lon_grid, lat_grid, dlon_grid,  &
       nlon_tess, nlat_tess, ilon_min, ilon_max, ilat_min, ilat_max,           &
       idx_lon, idx_lat, dateline_offset, yn_polar_exception )

    ! -------------------------------------------------------------------------
    ! Preparation for PMT (Poor Man's Tessellation). Similar to regular set-up,
    ! but doesn't require tweaking of "irregular" pixels that can't be handled
    ! by R. Spurr's true tessellation code.
    ! -------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlon_grid, nlat_grid
    REAL    (KIND=r8),                        INTENT (IN) :: dlon_grid
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (IN) :: lon_grid
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (4),         INTENT (INOUT) :: pix_lat, pix_lon

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                        INTENT (OUT) :: dateline_offset, nlon_tess, nlat_tess
    INTEGER (KIND=i4),                        INTENT (OUT) :: ilon_min, ilon_max, ilat_min, ilat_max
    LOGICAL,                                  INTENT (OUT) :: yn_polar_exception
    INTEGER (KIND=i4), DIMENSION (4),         INTENT (OUT) :: idx_lon, idx_lat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                :: i

    ! ---------------------
    ! Check for polar pixel
    ! ---------------------
    yn_polar_exception = .FALSE.
    CALL check_for_polar_pixel ( &
         pix_lat, pix_lon, 0.0_r8, dlon_grid, dateline_offset, yn_polar_exception )
    IF ( yn_polar_exception ) RETURN


    CALL pmt_find_indices_for_tessellation (                    &
         pix_lon, pix_lat, nlon_grid, lon_grid(1:nlon_grid), nlat_grid, lat_grid(1:nlat_grid),  &
         idx_lon, idx_lat, ilon_min, ilon_max, ilat_min, ilat_max, nlon_tess, nlat_tess   )

    RETURN
  END SUBROUTINE pmt_prepare_omi_tessell


  SUBROUTINE prepare_omi_tessell ( &
       pix_lat, pix_lon, dlon_grid, dlat_grid, nlon_grid, nlat_grid, lon_grid, lat_grid, &
       polar_lat, nlon_tess, tess_lonpts, nlat_tess, tess_latpts, tess_idx, tess_satpix, &
       tess_pars, tess_orient, tess_idx_abs, dateline_offset, yn_polar_exception, yn_good_pix )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlon_grid, nlat_grid
    REAL    (KIND=r8),                        INTENT (IN) :: dlon_grid, dlat_grid, polar_lat
    REAL    (KIND=r8), DIMENSION (4),         INTENT (IN) :: pix_lat, pix_lon
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (IN) :: lon_grid
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                        INTENT (OUT) :: dateline_offset, tess_orient
    INTEGER (KIND=i4),                        INTENT (OUT) :: nlon_tess, nlat_tess
    LOGICAL,                                  INTENT (OUT) :: yn_polar_exception, yn_good_pix
    INTEGER (KIND=i4), DIMENSION (5),         INTENT (OUT) :: tess_idx
    INTEGER (KIND=i4), DIMENSION (2),         INTENT (OUT) :: tess_idx_abs
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (OUT) :: tess_lonpts
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (OUT) :: tess_latpts
    REAL    (KIND=r8), DIMENSION (4,2),       INTENT (OUT) :: tess_satpix
    REAL    (KIND=r8), DIMENSION (6,4),       INTENT (OUT) :: tess_pars

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                :: i
    REAL    (KIND=r8), DIMENSION (4) :: locpix_lat, locpix_lon, tmplat, tmplon
    LOGICAL                          :: yn_side1, yn_side3, yn_irregular, yn_modified, yn_problem

    ! -------------------------------------------------------
    ! Save pixel coordinates in local arrays for manipulation
    ! -------------------------------------------------------
    tmplat = pix_lat ; tmplon = pix_lon

    ! ----------------------------------------------------------------------
    ! We can't have pixel coordinates that coincide with a grid point.
    ! For those cases we have to make a small adjustment, which is applied
    ! in the negative: subtracted from positive values and added to negative
    ! ones. This assures that we don't adjust beyond the maximum values of
    ! allowed geolocations.
    ! Also, for OMI, this adjustment can't be too large due to the small
    ! pixel sizes.
    ! ----------------------------------------------------------------------
    WHERE ( MOD(pix_lon, dlon_grid) == 0.0_r8 ) &
         tmplon = pix_lon + SIGN(small_dpix, pix_lon)
    WHERE ( MOD(pix_lat, dlat_grid) == 0.0_r8 ) &
         tmplat = pix_lat + SIGN(small_dpix, pix_lat)

    ! ---------------------
    ! Check for polar pixel
    ! ---------------------
    dateline_offset = 0_i4  ;  yn_polar_exception = .FALSE.
    CALL check_for_polar_pixel ( &
         tmplon, tmplat, polar_lat, dlon_grid, dateline_offset, yn_polar_exception )
    IF ( yn_polar_exception ) RETURN

    ! -------------------------------------------------------------------------
    ! Yet another layer of backup: In order to adjust the shape of pixels that
    ! don't conform to the tessellation requirements, we want to save and grid-
    ! and date-line adjusted pixels in their original ordering.
    ! -------------------------------------------------------------------------
    locpix_lon = tmplon  ;  locpix_lat = tmplat

    ! ---------------------------------------------------
    ! Find the correct order of the pixels as required by the tessellation
    ! package. See the following routine for details.
    ! ---------------------------------------------------
    CALL order_pixel_for_tessellation ( tmplon, tmplat )

    ! ---------------------------------------------------------------
    ! Now we have the check whether we can actually process the pixel,
    ! and whether we may have to make adjustments.
    ! ---------------------------------------------------------------
    yn_good_pix        = .TRUE.
    tess_satpix(1:4,1) = tmplon(1:4)
    tess_satpix(1:4,2) = tmplat(1:4)

    yn_problem = .FALSE.

    CALL parsetup ( tess_satpix, tess_pars, tess_orient, yn_problem )
    IF ( yn_problem ) THEN
       yn_good_pix = .FALSE.
       RETURN
    END IF

    CALL check_goodness_of_pixel ( &
         tess_satpix, tess_pars, yn_side1, yn_side3, yn_irregular )

    ! ----------------------------------------------------------------------
    ! If any of the above three logicals are .TRUE. then we have encountered
    ! problem.But we may be able to tweak the pixel into conforming shape.
    ! ----------------------------------------------------------------------
    IF ( ANY( (/ yn_side1, yn_side3, yn_irregular /) ) ) THEN
       yn_good_pix = .FALSE. ;  yn_modified = .FALSE.

       ! ----------------------------------------------------------------
       ! Restore pixel corners before reordering. This makes it easier to
       ! manipulate, since we don't have to worry about the orientation
       ! of the pixel.
       ! ----------------------------------------------------------------
       tmplon = locpix_lon  ;  tmplat = locpix_lat

       IF ( yn_irregular ) CALL tweak_irregular_pixel ( tmplon, tmplat, yn_modified )

       IF ( yn_modified ) THEN
          ! -------------------------------------------
          ! Modification requires reordering of pixels.
          ! -------------------------------------------
          CALL order_pixel_for_tessellation ( tmplon, tmplat )
          tess_satpix(1:4,1) = tmplon(1:4)
          tess_satpix(1:4,2) = tmplat(1:4)
          ! ----------------------------------------------------------
          ! If we have modified the pixel, it better conform to the
          ! required tessellation shape by now. We check this as above
          ! and if it still doesn't conform, we quit on it.
          ! ----------------------------------------------------------
          CALL parsetup ( tess_satpix, tess_pars, tess_orient, yn_problem )
          IF ( yn_problem ) THEN
             yn_good_pix = .FALSE.
          ELSE

             CALL check_goodness_of_pixel ( &
                  tess_satpix, tess_pars, yn_side1, yn_side3, yn_irregular )
             IF ( .NOT. ANY( (/ yn_side1, yn_side3, yn_irregular /) ) )  yn_good_pix = .TRUE.
          END IF

       END IF
    END IF
    ! -----------------------------------------------------------------
    ! If after all of this we still don't have a good pixel, we return.
    ! -----------------------------------------------------------------
    IF ( .NOT. yn_good_pix ) RETURN

    ! -------------------------------------------------------------------
    ! Next we find the lower- and upper-bounding indices of the satellite
    ! footprint in the surface grid.
    ! -------------------------------------------------------------------
    tmplon(1:4) = tess_satpix(1:4,1)
    tmplat(1:4) = tess_satpix(1:4,2)
    CALL find_indices_for_tessellation (                                         &
         tmplon, tmplat, nlon_grid, lon_grid, nlat_grid, lat_grid,               &
         nlon_tess, tess_lonpts, nlat_tess, tess_latpts, tess_idx, tess_idx_abs, &
         yn_good_pix )

    RETURN
  END SUBROUTINE prepare_omi_tessell

  SUBROUTINE tweak_irregular_pixel ( lon, lat, yn_modified )        

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,                       INTENT (INOUT) :: yn_modified  
    REAL (KIND=r8), DIMENSION (4), INTENT (INOUT) :: lon, lat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                :: j0, k0, j1, k1
    REAL    (KIND=r8), DIMENSION (4) :: tmplon, tmplat

    REAL (KIND=r8) :: adjlon_14, adjlon_23, adjlat_12, adjlat_34


    tmplon = lon ; tmplat = lat

    j0 = MAXVAL(MINLOC(tmplon)) ; j1 = MAXVAL(MAXLOC(tmplon))
    k0 = MAXVAL(MINLOC(tmplat)) ; k1 = MAXVAL(MAXLOC(tmplat))

    ! ------------------------------------------------------------
    ! The pedestrian approach to getting the pixels into shape:
    ! Check for all the eventualities and correct them one by one.
    ! For simplicity, we only adjust longitudes.
    ! -----------------------------------------------------------
    adjlon_14 = ABS(tmplon(4) - tmplon(1)) + small_dpix
    adjlon_23 = ABS(tmplon(3) - tmplon(2)) + small_dpix
    adjlat_12 = ABS(tmplat(1) - tmplat(2)) + small_dpix
    adjlat_34 = ABS(tmplat(3) - tmplat(4)) + small_dpix

    ! -------------------------------------------------
    ! The problem lies with the minimum latitude corner
    ! -------------------------------------------------
    IF ( k0 == j0 ) THEN
       ! ----------------------------------------
       ! Case 1: MIN(lon) coincides with MIN(lat)
       ! ----------------------------------------
       SELECT CASE ( k0 )
       CASE (2:3)
          tmplon(k0) = tmplon(k0) + adjlon_23
       CASE DEFAULT
          tmplon(k0) = tmplon(k0) + adjlon_14
       END SELECT
    ELSE IF ( k0 == j1 ) THEN
       ! ----------------------------------------
       ! Case 2: MAX(lon) coincides with MIN(lat)
       ! ----------------------------------------
       SELECT CASE ( k0 )
       CASE (2:3)
          tmplon(k0) = tmplon(k0) - adjlon_23
       CASE DEFAULT
          tmplon(k0) = tmplon(k0) - adjlon_14
       END SELECT
    END IF

    ! -------------------------------------------------
    ! The problem lies with the maximum latitude corner
    ! -------------------------------------------------
    IF ( k1 == j0 ) THEN
       ! ----------------------------------------
       ! Case 3: MIN(lon) coincides with MAX(lat)
       ! ----------------------------------------
       SELECT CASE ( k1 )
       CASE (2:3)
          tmplon(k1) = tmplon(k1) + adjlon_23
       CASE DEFAULT
          tmplon(k1) = tmplon(k1) + adjlon_14
       END SELECT
    ELSE IF ( k1 == j1 ) THEN
       ! ----------------------------------------
       ! Case 3: MAX(lon) coincides with MAX(lat)
       ! ----------------------------------------
       SELECT CASE ( k1 )
       CASE ( 1 )
          tmplon(k1) = tmplon(k1) - adjlon_14
          !tmplat(k1) = tmplat(k1) - SIGN(adjlat_12,tmplat(k1))
       CASE ( 2 )
          tmplon(k1) = tmplon(k1) - adjlon_23
          !tmplat(k1) = tmplat(k1) - SIGN(adjlat_12,tmplat(k1))
       CASE ( 3 )
          tmplon(k1) = tmplon(k1) - adjlon_23
          !tmplat(k1) = tmplat(k1) - SIGN(adjlat_34,tmplat(k1))
       CASE ( 4 )
          tmplon(k1) = tmplon(k1) - adjlon_14
          !tmplat(k1) = tmplat(k1) - SIGN(adjlat_34,tmplat(k1))
       END SELECT
    END IF

    ! -----------------------------
    ! Check if anything has changed
    ! -----------------------------
    IF ( ANY ( tmplon /= lon ) .OR. ANY ( lat /= tmplat ) ) THEN
       lon = tmplon ; lat = tmplat ; yn_modified = .TRUE.
    END IF

    RETURN
  END SUBROUTINE tweak_irregular_pixel

  SUBROUTINE check_goodness_of_pixel ( &
       tess_satpix, tess_pars, yn_side1, yn_side3, yn_irregular )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), DIMENSION (6,4), INTENT (IN) :: tess_pars
    REAL (KIND=r8), DIMENSION (4,2), INTENT (IN) :: tess_satpix

    ! ----------------
    ! Output variables
    ! ----------------
    LOGICAL, INTENT (OUT) :: yn_side1, yn_side3, yn_irregular

    ! ---------------
    ! Local variables
    ! ---------------
    LOGICAL                       :: yn_good_pix
    REAL (KIND=r8), DIMENSION (2) :: lonvals
    REAL (KIND=r8)                :: tmp1, tmp2

    ! -----------------
    ! Initialize output
    ! -----------------
    yn_side1 = .FALSE. ; yn_side3 = .FALSE. ;  yn_irregular = .FALSE.

    yn_good_pix = .TRUE.
    lonvals(1:2) = tess_satpix(1:2,1)
    CALL check_pixel_side ( lonvals, tess_pars(2,1), yn_good_pix )
    IF ( .NOT. yn_good_pix ) yn_side1 = .TRUE.

    yn_good_pix = .TRUE.
    lonvals(1:2) = tess_satpix(3:4,1)
    CALL check_pixel_side ( lonvals, tess_pars(2,1), yn_good_pix )
    IF ( .NOT. yn_good_pix ) yn_side3 = .TRUE.

    yn_good_pix = .TRUE.
    CALL check_pixel_irregular ( tess_satpix(1:4,2), yn_good_pix )
    IF ( .NOT. yn_good_pix ) yn_irregular = .TRUE.

    RETURN
  END SUBROUTINE check_goodness_of_pixel

  SUBROUTINE check_pixel_side ( lon_satpix, tess_par, yn_good_pix )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),                INTENT (IN) :: tess_par
    REAL (KIND=r8), DIMENSION (2), INTENT (IN) :: lon_satpix

    ! ---------------
    ! Output variable
    ! ---------------
    LOGICAL, INTENT (INOUT) :: yn_good_pix

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2

    tmp1 = tess_par + pihalf_deg ; tmp2 = tess_par - pihalf_deg
    IF ( ( tmp1 > lon_satpix(1) .AND. tmp1 < lon_satpix(2) )   .OR. &
         ( tmp2 > lon_satpix(1) .AND. tmp2 < lon_satpix(2) ) ) THEN
       yn_good_pix = .FALSE.
    ENDIF

    RETURN
  END SUBROUTINE check_pixel_side

  SUBROUTINE check_pixel_irregular ( lat_satpix, yn_good_pix )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), DIMENSION (4), INTENT (IN) :: lat_satpix

    ! -----------------
    ! Modified variable
    ! -----------------
    LOGICAL, INTENT (INOUT) :: yn_good_pix

    ! check for corner 2 and 4 to see if top most and bottom most
    IF ( MAXVAL(MAXLOC(lat_satpix(1:4))) /= 2   .OR. &
         MAXVAL(MINLOC(lat_satpix(1:4))) /= 4  ) THEN
       yn_good_pix = .FALSE.
    ENDIF

    RETURN
  END SUBROUTINE check_pixel_irregular


  SUBROUTINE check_for_polar_pixel ( &
       locpix_lon, locpix_lat, polar_lat, dlon_grid, dateline_offset, yn_polar_pix )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: polar_lat, dlon_grid

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (4), INTENT (INOUT) :: locpix_lon, locpix_lat

    ! ---------------
    ! Output variable
    ! ---------------
    LOGICAL,           INTENT (OUT) :: yn_polar_pix
    INTEGER (KIND=i4), INTENT (OUT) :: dateline_offset

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: adj_lon


    ! -----------------------------------------------------------------------
    ! Before we delve into setting up any indices, it is the time to deal
    ! with special cases. Those should come in basically two flavors:
    !
    !   (1) Longitudes run across the dateline;
    !   (2) Latitudes run across the pole
    !
    ! where case "(2)" contains "(1)" as an added benefit. We won't attempt
    ! to solve the polar problem at this point, but the first one should be
    ! fixed without too much trouble by adding an offset that pulls the pixel
    ! away from the dateline. Note that this does NOT involve pixels with an
    ! upper longitude in the last grid row, since we have set up our grid 
    ! just so that we always run from [-180,+180].
    ! -----------------------------------------------------------------------
    dateline_offset = 0_i4  ;  yn_polar_pix = .FALSE.
    IF ( ABS(MAXVAL(locpix_lat)) > polar_lat .AND. &
         ABS ( MAXVAL(locpix_lon) - MINVAL(locpix_lon) ) > pi_deg ) THEN

       ! -----------------------------------------------------------
       ! Find the offset of the eastern-most negative value from the
       ! dateline boundary. The offset is the number of grid cells
       ! by which the pixel will have to be pulled to the west in 
       ! order to bring it clear of the dateline. We add an extra
       ! cell for safety.
       ! -----------------------------------------------------------
       adj_lon = pi_deg - ABS(MAXVAL( locpix_lon, MASK=locpix_lon <= 0.0_r8 ))
       dateline_offset = INT( ANINT ( adj_lon / dlon_grid ) ) + 1

       ! ----------------------------------------------------------------
       ! The minimum must lie in the negative half. We figure out the 
       ! distance from -180.0 and subtract the appropriate number of grid
       ! cell entries from all longitudes. Anything < 0 is adjusted by 
       ! 2*Pi. If after this we still have longitudes differing by more 
       ! than Pi, this has to be a polar pixel and we can't proceed.
       !
       ! After the tessellation,  we have to adjust the output by the
       ! offset determined here.
       !
       ! ....., of course, if the grid is not regular, we are SOL.
       ! ----------------------------------------------------------------
       locpix_lon = locpix_lon - REAL(dateline_offset, KIND=r8) * dlon_grid
       WHERE ( locpix_lon < -pi_deg ) locpix_lon = twopi_deg + locpix_lon

       ! --------------------------------------------------------------------
       ! If we still have differences > Pi after the adjustment, we must
       ! have a truly polar pixel. And we don't have any remedy for this yet.
       ! --------------------------------------------------------------------
       IF ( ABS ( MAXVAL(locpix_lon) - MINVAL(locpix_lon) ) > pi_deg ) THEN
          ! -------------------------------
          ! This must be a polar pixel ...
          ! -------------------------------
          yn_polar_pix = .TRUE.
          RETURN
       END IF
    END IF


    ! ------------------------------------------------------------------------
    ! And yet another adjustment: Occasionally we encounter pixels with too
    ! wide a spread in longitude, usually going along with at least some of
    ! the latitude (and some longitude?) values being close to 0.0. This leads
    ! to ragged East-West striping around the Lon=0 line. Since we don't know
    ! what causes those exceptions, we treat them for now just like polar
    ! pixels, will say, we skip them.
    ! ------------------------------------------------------------------------
    IF ( MINVAL(ABS(locpix_lat)) < 0.0001_r8 ) THEN
       IF ( (MAXVAL(locpix_lon)-MINVAL(locpix_lon) > 1.0_r8) .OR. &
            (MAXVAL(locpix_lat)-MINVAL(locpix_lat) > 1.0_r8) ) THEN
          yn_polar_pix = .TRUE. 
          !WRITE (*,'(4F10.4, A, 4F10.4, 3X, L1)') locpix_lat, ' --- ', locpix_lon, yn_polar_pix
       ELSE
          !WRITE (*,'(4F10.4, A, 4F10.4, 3X, L1)') locpix_lat, ' ||| ', locpix_lon, yn_polar_pix
       END IF
    END IF

    RETURN
  END SUBROUTINE check_for_polar_pixel

  SUBROUTINE order_pixel_for_tessellation ( locpix_lon, locpix_lat )

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (4), INTENT (INOUT) :: locpix_lon, locpix_lat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                :: i, ipos
    INTEGER (KIND=i4), DIMENSION (4) :: pix_order
    REAL    (KIND=r8), DIMENSION (4) :: tmplon
    INTEGER (KIND=i4), PARAMETER     :: huge_latlon = 999.9_r8

    ! ---------------------------------------------------
    ! Find the correct order of the pixels:
    !   * Lon(pix(1)) is smallest (western-most)
    !   * Lon(pix(3)) is largest  (eastern-most)
    !   * Lat(pix(2)) > Lat(pix(1))
    ! ---------------------------------------------------
    tmplon = locpix_lon
    DO i = 1, 4
       ipos = MINVAL(MINLOC(tmplon))
       pix_order(i) = ipos
       tmplon(ipos) = huge_latlon
    END DO

    ! -----------------------------------------------------
    ! At this point PIX_ORDER(1) and PIX_ORDER(4) contain
    ! the smallest and largest longitude respectively. So
    ! now let's make sure that PIX_ORDER(2) indeed contains
    ! the larger of the two latitudes in the middle.
    ! -----------------------------------------------------
    ! First, check whether Lat(2) > Lat(3)
    ! ------------------------------------
    IF ( locpix_lat(pix_order(2)) < locpix_lat(pix_order(3)) ) &
         pix_order(2:3) = (/ pix_order(3), pix_order(2) /)

    ! ----------------------------------------------------
    ! The maximum longitude will end up in position 3, so
    ! we flip the order in PIX_ORDER for easier assignment
    ! of the properly ordered arrays.
    ! ----------------------------------------------------
    pix_order(3:4) = (/ pix_order(4), pix_order(3) /)

    ! ---------------------------------------
    ! The reassignment is now a piece of cake
    ! ---------------------------------------
    locpix_lon = locpix_lon(pix_order) ; locpix_lat = locpix_lat(pix_order)

    RETURN
  END SUBROUTINE order_pixel_for_tessellation

  SUBROUTINE find_indices_for_tessellation (                                   &
       locpix_lon, locpix_lat, nlon_grid, lon_grid, nlat_grid, lat_grid,       &
       nlon_tess, tess_lonpts, nlat_tess, tess_latpts, tess_idx, tess_idx_abs, &
       yn_good_pix )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nlon_grid, nlat_grid
    REAL    (KIND=r8), DIMENSION (4),         INTENT (IN) :: locpix_lon, locpix_lat
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (IN) :: lon_grid
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (IN) :: lat_grid

    ! -----------------
    ! Modified variable
    ! -----------------
    LOGICAL, INTENT (INOUT) :: yn_good_pix

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                        INTENT (OUT) :: nlon_tess, nlat_tess
    INTEGER (KIND=i4), DIMENSION (5),         INTENT (OUT) :: tess_idx
    INTEGER (KIND=i4), DIMENSION (2),         INTENT (OUT) :: tess_idx_abs
    REAL    (KIND=r8), DIMENSION (nlon_grid), INTENT (OUT) :: tess_lonpts
    REAL    (KIND=r8), DIMENSION (nlat_grid), INTENT (OUT) :: tess_latpts

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4), PARAMETER     :: int_missing = -9999_i4
    INTEGER (KIND=i4), DIMENSION (4) :: istart, jstart
    INTEGER (KIND=i4)                :: i
    !INTEGER (KIND=i4), EXTERNAL      :: offset
    
    ! -------------------------------------------------------------------
    ! Next we find the lower- and upper-bounding indices of the satellite
    ! footprint in the surface grid.
    ! -------------------------------------------------------------------
    istart = int_missing ; jstart = int_missing
    DO i = 1, 4
       istart(i) = offset ( 1, lon_grid(1:nlon_grid), nlon_grid, 1, locpix_lon(i) )
       jstart(i) = offset ( 1, lat_grid(1:nlat_grid), nlat_grid, 1, locpix_lat(i) )
    END DO
    IF ( ANY ( istart(1:4) == int_missing ) .OR. &
         ANY ( jstart(1:4) == int_missing ) )  THEN
       yn_good_pix = .FALSE.
       RETURN
    END IF

    IF ( ANY ( istart == 0 ) ) THEN
       WHERE ( istart == 0 )
          istart = 1
       END WHERE
    END IF
    ! TPK: Added the following on 09/09 to prevent jstart(4)=0
    IF ( ANY ( jstart == 0 ) ) THEN
       WHERE ( jstart == 0 )
          jstart = 1
       END WHERE
    END IF
    ! ---------------------------------------------------------------
    ! Save the minimum values of ISTART and JSTART, since they define
    ! the lower left bounding grid box, relative to which the other
    ! tessellation indices are defined, and which have to be used in
    ! the assignment of the final output arrays.
    ! ---------------------------------------------------------------
    tess_idx_abs(1:2) = (/ MINVAL(istart), MINVAL(jstart) /)

    ! -----------------------------------------------------------
    ! Now we are able to determine the number of grid points used
    ! in the tessellation.
    ! -----------------------------------------------------------
    nlon_tess = istart(3)-istart(1) + 2
    nlat_tess = jstart(2)-jstart(4) + 2

    ! ---------------------------------------
    ! Minimum/Maximum check for array indices
    ! ---------------------------------------
    IF ( istart(3)+nlon_tess-1 > nlon_grid ) THEN
       istart(3) = istart(3) - 1
       nlon_tess = nlon_tess - 1
    END IF
    IF ( jstart(4)+nlat_tess-1 > nlat_grid ) THEN
       jstart(4) = jstart(4) - 1
       nlat_tess = nlat_tess - 1
    END IF

    ! -----------------------------------------------------------------------
    ! Finally we asign the longitudes and latitudes used in the tessellation,
    ! and we save all start/stop indices in a single array for easier use.
    ! -----------------------------------------------------------------------
    tess_lonpts (1:nlon_tess) = lon_grid(istart(1):istart(3)+nlon_tess-1)
    tess_latpts (1:nlat_tess) = lat_grid(jstart(4):jstart(4)+nlat_tess-1)

    tess_idx = (/ 1, istart(2)-istart(1)+1, nlon_tess-1, istart(4)-istart(1)+1, nlat_tess-1 /)

    RETURN
  END SUBROUTINE find_indices_for_tessellation

  ! ==========================
  ! Start of R. Spurr Routines
  ! ==========================

  SUBROUTINE tesselate_areamaster  &
       ( xgrid, ygrid, nxdim, nydim,  &
       corner_coords, side_params, orient,  &
       istart1, istart2, istart3, istart4, jfinis,  &
       area, sum, ylimit_lower, ylimit_upper )

    !    ********************************************************************
    !    *   robert spurr, november 1998					*
    !    *   sao, 60 garden street, cambridge, ma 02138, usa		*
    !    *   +1 (617) 496 7819; email rspurr@cfa.harvard.edu		*
    !    *									*
    !    *   algorithm may be subject of licensing agreement		*
    !    ********************************************************************

    !  master module inputs
    !  ====================

    !  corenr coordinates and side curves
    !  ----------------------------------

    !  we have 4 corner coordinates corner_coords, and parameters side_params 
    !  describing the 4 curves joining these corners, and an orientation (orient) 
    !  for the footprint as follows :

    !        2
    !        /\
    !       /  \
    !      /    \
    !     1\     \3		orientation = -1 (corner 2 before corner 4)
    !       \    /
    !        \  /
    !         \/4

    !          2
    !         /\
    !        /  \
    !       /    \3
    !      /     / 
    !     1\    /		orientation = +1 (corner 4 before corner 2)
    !       \  /
    !        \/4

    !  note that the actual area routines do not need to know the curve
    !  equations - all information is contained in the side_params array.

    !  offset and box limits
    !  ---------------------

    !  the following inputs are all offsets in the x-y grid which contains
    !  the footprint defined by corner_coords. the user should ensure that any 
    !  original grid of values be reduced to cover the footprint (avoids
    !  unnecessary computation inside the tessellation module)

    !		istart1, istart2, istart3, istart4, jfinis

    !   jfinis is the offset in y-space of corner 2, (the offset for corner 4
    !   should be set to 1).
    !   istart1, istart2, istart3, istart4 are the offsets in the x-direction
    !   for corners 1 through 4 respectively. istart1 = 1 upon input.

    !  grid and data input
    !  -------------------

    !  integers nxdim and nydim are just dimensioning values for the x and y
    !  grid values xgrid and ygrid. the array of data bins albused covers
    !  only the box containing the footprint, and this should be prepared
    !  before entry to the tessellation master.

    !  master module output
    !  ====================

    !  the outputs include the tessellated area bins area which cover the pixel
    !  and the total area of the pixel (sum), the area-weighted albedo (albedo)

    !  also output (for reference) are the ranges over which grid bins in the
    !  y direction are required. these are dependent on the x-slice.

    !  subroutine declarations
    !  =======================

    !  input

    INTEGER		nxdim, nydim, orient
    REAL*8		xgrid(nxdim), ygrid(nydim)
    INTEGER		istart1, istart2, istart3, istart4, jfinis
    REAL*8		corner_coords(4,2), side_params(6,4)

    !  output

    REAL*8		area(nydim, nxdim), sum
    INTEGER		ylimit_lower(nxdim)
    INTEGER		ylimit_upper(nxdim)

    !  local variables
    !  ===============

    !  options for pathway throught the tessellation

    LOGICAL		double_first, double_last
    LOGICAL		do_group1, do_group2
    LOGICAL		triple_first, triple_last, double_middle

    !  counters

    INTEGER		i, j

    !  selection module, chooses path through the x-slicing
    !  ====================================================

    CALL tesselate_options_chooser  &
         ( xgrid, ygrid, nxdim, nydim,  &
         istart1, istart2, istart3, istart4, jfinis, orient,  &
         do_group1, double_first, double_last,  &
         do_group2, triple_first, triple_last, double_middle )

    !  initialise output
    !  =================

    DO i = istart1, istart3
       DO j = 1, jfinis
          area(j,i) = 0.0
       ENDDO
    ENDDO

    !  tessellate for group 1 options
    !  ==============================

    IF ( nxdim > 0 .AND. nydim > 0 ) THEN  !tpk

       IF ( do_group1 ) THEN

          CALL tesselate_options_1  &
               ( xgrid, ygrid, nxdim, nydim, corner_coords, orient,   &
               double_first, double_last, side_params,  &
               istart1, istart2, istart3, istart4, jfinis,  &
               area, ylimit_lower, ylimit_upper )

          !  tessellate for group 2 options
          !  ==============================

       ELSE IF ( do_group2 ) THEN

          CALL tesselate_options_2  &
               ( xgrid, ygrid, nxdim, nydim, corner_coords, orient,   &
               triple_first, triple_last, double_middle, side_params,    &
               istart1, istart2, istart3, istart4, jfinis,  &
               area, ylimit_lower, ylimit_upper )

       ENDIF

    END IF !tpk

    !  area sum and weighted albedo computation
    !  ========================================

    sum = 0.0
    DO i = istart1, istart3
       DO j = ylimit_lower(i), ylimit_upper(i)
          IF ( j > 0 ) sum = sum + area(j,i)   !tpk
       ENDDO
    ENDDO

    !  finish

  END SUBROUTINE tesselate_areamaster

  !

  SUBROUTINE tesselate_options_chooser  &
       ( xgrid, ygrid, nxdim, nydim,  &
       istart1, istart2, istart3, istart4,  &
       jfinis, orient,  &
       do_group1, double_first, double_last,  &
       do_group2, triple_first, triple_last, double_middle )

    !  input

    INTEGER		nxdim, nydim
    REAL*8		xgrid(nxdim)
    REAL*8		ygrid(nydim)
    INTEGER		istart1, istart2, istart3, istart4
    INTEGER		jfinis, orient

    !  output (pathway options

    LOGICAL		double_first, double_last
    LOGICAL		do_group1, do_group2
    LOGICAL		triple_first, triple_last, double_middle

    !   initialise options

    double_first = .FALSE.
    double_last = .FALSE.
    triple_first = .FALSE.
    triple_last = .FALSE.
    double_middle = .FALSE.
    do_group1 = .FALSE.
    do_group2 = .FALSE.

    !  group 1
    !  =======

    !  determining factor: corners 2 and 4 are not in the same x-slice

    IF ( istart2.NE.istart4 ) THEN
       do_group1 = .TRUE.
       double_first = ((istart1.EQ.istart2).OR.(istart1.EQ.istart4))
       double_last  = ((istart3.EQ.istart2).OR.(istart3.EQ.istart4))
       RETURN
    ENDIF

    !  group 2
    !  =======

    !  determining factor: corners 2 and 4 are in the same x-slice

    IF ( istart2.EQ.istart4 ) THEN
       do_group2 = .TRUE.
       triple_first = ((istart1.EQ.istart2).AND.do_group2)
       triple_last  = ((istart3.EQ.istart2).AND.do_group2)
       double_middle= ((.NOT.triple_first).AND.(.NOT.triple_last))
       RETURN
    ENDIF

    !  finish

  END SUBROUTINE tesselate_options_chooser
  !  list of modules in tessellation package :

  !	tesselate_options_1
  !	tesselate_options_2
  !	singleside_slicer
  !	doubleside_slicer
  !	single_left_slicer
  !	single_right_slicer
  !	single_topbot_slicer
  !	double_left_slicer
  !	double_right_slicer
  !	double_middle_slicer
  !	triple_left_slicer
  !	triple_right_slicer
  !	quadruple_slicer
  !	offset

  !    ********************************************************************
  !    *   robert spurr, november 1998					*
  !    *   sao, 60 garden street, cambridge, ma 02138, usa		*
  !    *   +1 (617) 496 7819; email rspurr@cfa.harvard.edu		*
  !    *									*
  !    *   algorithm may be subject of licensing agreement		*
  !    ********************************************************************

  !  function parameters are stored in array pars which should
  !  be pre-computed. the modules do not know what function is being
  !  used in the area and line intercept determinations.

  SUBROUTINE tesselate_options_1  &
       ( xgrid, ygrid, nxdim, nydim,  &
       cc, orient, double_first, double_last, pars,  &
       istart1, istart2, istart3, istart4, jfinis,  &
       area, ylimit_lower, ylimit_upper )

    !  input 
    !  =====

    !  grid values

    INTEGER		nxdim, nydim
    REAL*8		xgrid(nxdim), ygrid(nydim)

    !  corner coordinates

    REAL*8		cc(4,2)

    !  offsets

    INTEGER		istart1, istart2, istart4, istart3, jfinis

    !  function parameters

    REAL*8		pars(6,4)

    !  control of options

    INTEGER		orient
    LOGICAL		double_first, double_last

    !  output
    !  ======

    REAL*8		area(nydim,nxdim)
    INTEGER		ylimit_lower(nxdim)
    INTEGER		ylimit_upper(nxdim)

    !  local variables
    !  ===============

    REAL*8		xold, xnew, yold_l, yold_u, ynew_l, ynew_u
    INTEGER		is, is_next, is_first, is_last, is_beg, is_end
    INTEGER		i_next, jc1, jc2, jold_l, jold_u, jnew_l, jnew_u
    INTEGER		side_l, side_u, side_c, side_a, istart_dummy

    !  parameter values and indices

    INTEGER		side_1, side_2, side_3, side_4
    INTEGER		asct, parallel, perpen
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )
    PARAMETER	( asct = 1, parallel = 1, perpen = -1 )

    !  external functions

    !REAL*8		sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset
    !EXTERNAL	offset, sidefunction

    !  double first option
    !  ###################

    IF ( double_first ) THEN

       is_first = 1
       i_next = istart1 + 1
       is_next = is_first
       xnew = xgrid(i_next)

       !  orientation = 1  (corners 4 and 1; upper/lower sides are 1 and 3)
       !  orientation = -1 (corners 2 and 1; upper/lower sides are 2 and 4)

       IF ( orient .EQ. 1 ) THEN
          side_c = side_4
          side_u = side_1
          side_l = side_3
       ELSE IF ( orient .EQ. -1 ) THEN
          side_c = side_2
          side_u = side_2
          side_l = side_4
       ENDIF

       !  find where upper/lower sides cut next x-level

       ynew_u = sidefunction ( pars, xnew, side_u )
       ynew_l = sidefunction ( pars, xnew, side_l )

       ! tpk additions
       ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
       ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )

       jnew_l  = offset ( 1, ygrid, nydim, asct, ynew_l )
       jc1     = offset ( 1, ygrid, nydim, asct, cc(side_1,2) )
       jc2     = offset ( 1, ygrid, nydim, asct, cc(side_c,2) )
       jnew_u  = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  get double corner areas

       CALL double_left_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jnew_l, jc1, jc2, jnew_u, xnew, ynew_l, ynew_u,  &
            cc(side_1,1), cc(side_1,2),  &
            cc(side_c,1), cc(side_c,2),  &
            area(1,is_first) )

       !  assign limits

       IF ( orient .EQ. 1 ) THEN
          ylimit_lower(is_first) = jc2
          ylimit_upper(is_first) = jnew_u
       ELSE IF ( orient .EQ. -1 ) THEN
          ylimit_lower(is_first) = jnew_l
          ylimit_upper(is_first) = jc2
       ENDIF

    ENDIF

    !  single first option
    !  ###################

    IF ( .NOT. double_first ) THEN

       !  first slice (includes corner 1)
       !  -------------------------------

       is_first = 1
       i_next = istart1 + 1
       is_next = is_first

       !  find where sides 1 and 4 cut next x-level, 

       xnew   = xgrid(i_next)
       ynew_u = sidefunction ( pars, xnew, side_1 )
       ynew_l = sidefunction ( pars, xnew, side_4 )
       jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
       jc1    = offset ( 1, ygrid, nydim, asct, cc(side_1,2) )
       jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  assign corner area

       CALL single_left_slicer  &
            ( ygrid, nydim, nydim, pars,  &
            jnew_l, jc1, jnew_u,  &
            xnew, ynew_l, ynew_u,  &
            cc(side_1,1), cc(side_1,2),  &
            area(1,is_first) )

       !  assign limits

       ylimit_lower(is_first) = jnew_l
       ylimit_upper(is_first) = jnew_u

       !  slices between corners 1 and 2/4 (depends on orientation)
       !  --------------------------------

       !  how many x-slices 

       is_beg = is_first + 1
       IF ( orient.EQ.1 ) THEN
          istart_dummy = istart4
       ELSE IF ( orient.EQ.-1 ) THEN
          istart_dummy = istart2
       ENDIF
       is_end = istart_dummy - istart1

       !  start loop over x-slices (not done if is_beg > is_end)

       DO is = is_beg, is_end

          !  update left slice boundary values

          jold_l = jnew_l
          jold_u = jnew_u
          yold_l = ynew_l
          yold_u = ynew_u
          xold   = xnew

          !  find where sides borders cut next x-level

          i_next = i_next + 1
          IF ( i_next > 0 .AND. i_next <= nxdim ) THEN  !tpk
             xnew   = xgrid(i_next)
             ynew_u = sidefunction ( pars, xnew, side_1 )
             ynew_l = sidefunction ( pars, xnew, side_4 )
             jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
             jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

             !  assign areas for each x-slice

             CALL doubleside_slicer  &
                  ( ygrid, nydim, nydim, side_4, side_1, pars,  &
                  jold_l, jold_u, jnew_l, jnew_u, perpen, orient,  &
                  xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
                  area(1,is) )

             !  assign limits

             ylimit_lower(is) = jnew_l
             ylimit_upper(is) = jnew_u
          END IF
       ENDDO

       !  slice including corner 2/4 (depends on orientation)
       !  --------------------------

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold   = xnew

       is_next = is_end + 1
       i_next = i_next + 1

       IF ( i_next > 0 .AND. i_next <= nxdim ) THEN  !tpk
          xnew   = xgrid(i_next)

          !  orientation = 1  (corner 4; upper/lower sides are 1 and 3)
          !  orientation = -1 (corner 2; upper/lower sides are 2 and 4)

          IF ( orient .EQ. 1 ) THEN
             side_c = side_4
             side_u = side_1
             side_l = side_3
             side_a = side_u
             jc1 = 1
          ELSE IF ( orient .EQ. -1 ) THEN
             side_c = side_2
             side_u = side_2
             side_l = side_4
             side_a = side_l
             jc1 = jfinis
          ENDIF

          !  find where upper/lower sides cut next x-level

          ynew_u = sidefunction ( pars, xnew, side_u )
          ynew_l = sidefunction ( pars, xnew, side_l )

          ! tpk additions
          ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
          ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )
          
          jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
          jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )
          
          !  assign area

          CALL single_topbot_slicer  &
               ( ygrid, nydim, nydim, side_a, pars, -orient, orient,  &
               jold_l, jold_u, jc1, jnew_l, jnew_u,  &
               xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
               cc(side_c,1), cc(side_c,2),   &
               area(1,is_next) )

          !  assign limits

          IF ( orient .EQ. 1 ) THEN
             ylimit_lower(is_next) = jc1
             ylimit_upper(is_next) = jnew_u
          ELSE IF ( orient .EQ. -1 ) THEN
             ylimit_lower(is_next) = jnew_l
             ylimit_upper(is_next) = jc1
          ENDIF
       END IF !tpk

    ENDIF

    !  central section
    !  ###############

    !  slices between corner 2/4 and corner 4/2  (depends on orientation)
    !  ----------------------------------------

    !  how many slices

    is_beg = is_next + 1
    IF ( orient.EQ.1 ) THEN
       istart_dummy = istart2
       side_u = side_1
       side_l = side_3
    ELSE IF ( orient.EQ.-1 ) THEN
       istart_dummy = istart4
       side_u = side_2
       side_l = side_4
    ENDIF
    is_end = istart_dummy - istart1

    !  start loop over x-slices (not done if is_beg > is_end)

    DO is = is_beg, is_end

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold  = xnew

       !  find where upper and lower sides cut next x-level

       i_next = i_next + 1

       ! tpk addition
       IF ( i_next > nxdim ) CYCLE

       xnew = xgrid(i_next)
       ynew_u = sidefunction ( pars, xnew, side_u )
       ynew_l = sidefunction ( pars, xnew, side_l )

       ! tpk additions
       ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
       ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )

       jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
       jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  assign areas for each x-slice

       CALL doubleside_slicer                                  &
            ( ygrid, nydim, nydim, side_l, side_u, pars,       &
            jold_l, jold_u, jnew_l, jnew_u, parallel, orient,  &
            xold, xnew, yold_l, yold_u, ynew_l, ynew_u,        &
            area(1,is) )

       !  assign limits

       IF ( orient .EQ. 1 ) THEN
          ylimit_lower(is) = jold_l
          ylimit_upper(is) = jnew_u
       ELSE IF ( orient .EQ. -1 ) THEN
          ylimit_lower(is) = jnew_l
          ylimit_upper(is) = jold_u
       ENDIF
    ENDDO

    !  double last option
    !  ##################

    IF ( double_last .AND. is_end < nxdim ) THEN  !tpk addition

       !  final slice (includes corners 2/4 and 3)
       !  ----------------------------------------

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold  = xnew
       is_last = is_end + 1

       !  orientation = 1  (corners 2 and 3; upper/lower sides are 1 and 3)
       !  orientation = -1 (corners 4 and 3; upper/lower sides are 2 and 4)

       IF ( orient .EQ. 1 ) THEN
          side_c = side_2
       ELSE IF ( orient .EQ. -1 ) THEN
          side_c = side_4
       ENDIF

       !  corner limits

       jc1 = offset ( 1, ygrid, nydim, asct, cc(side_c,2) )
       jc2 = offset ( 1, ygrid, nydim, asct, cc(side_3,2) )

       !  get double corner

       CALL double_right_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jold_l, jc1, jc2, jold_u, xold, yold_l, yold_u,  &
            cc(side_c,1), cc(side_c,2), cc(side_3,1), cc(side_3,2),  &
            area(1,is_last) )

       !  assign limits

       IF ( orient .EQ. 1 ) THEN
          ylimit_lower(is_last) = jold_l
          ylimit_upper(is_last) = jc1
       ELSE IF ( orient .EQ. -1 ) THEN
          ylimit_lower(is_last) = jc1
          ylimit_upper(is_last) = jold_u
       ENDIF

    ENDIF

    !  single last option
    !  ##################

    IF ( .NOT. double_last ) THEN

       !  slice including corner 2/4 (depends on orientation)
       !  --------------------------

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold   = xnew
       is_next = is_end + 1
       i_next = i_next + 1

       ! *******************************************************
       ! the following is a fix to i_next going out of bounds.
       ! it is a departure from the original program and may be
       ! introducing an error. it did not do so for the case it
       ! was discovered, but that doesn't mean it won't do so 
       ! at other places. of course, indices should not be left
       ! to their own devices in such a cavalrier manner in the
       ! first place.  tpk, 19 feb 2005.
       ! *******************************************************
       i_next = MIN ( i_next, nxdim ) !tpk

       xnew   = xgrid(i_next)

       !  orientation = 1  (corner 2; upper/lower sides are 2 and 3)
       !  orientation = -1 (corner 4; upper/lower sides are 2 and 3)

       IF ( orient .EQ. 1 ) THEN
          side_c = side_2
          side_u = side_2
          side_l = side_3
          side_a = side_3
          jc1 = jfinis
       ELSE IF ( orient .EQ. -1 ) THEN
          side_c = side_4
          side_u = side_2
          side_l = side_3
          side_a = side_2
          jc1 = 1
       ENDIF

       !  find where upper/lower sides cut next x-level

       ynew_u = sidefunction ( pars, xnew, side_u )
       ynew_l = sidefunction ( pars, xnew, side_l )

       ! tpk additions
       ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
       ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )

       jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
       jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  assign area (orient controls top or bottom)

       CALL single_topbot_slicer  &
            ( ygrid, nydim, nydim, side_a, pars, orient, orient,  &
            jold_l, jold_u, jc1, jnew_l, jnew_u,  &
            xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
            cc(side_c,1), cc(side_c,2),   &
            area(1,is_next) )

       !  assign limits

       IF ( orient .EQ. 1 ) THEN
          ylimit_lower(is_next) = jold_l
          ylimit_upper(is_next) = jc1
       ELSE IF ( orient .EQ. -1 ) THEN
          ylimit_lower(is_next) = jc1
          ylimit_upper(is_next) = jold_u
       ENDIF

       !  slices between corner 2 and corner 3
       !  ------------------------------------

       is_beg = is_next + 1
       is_end = istart3 - istart1

       DO is = is_beg, is_end

          !  update left slice boundary values

          jold_l = jnew_l
          jold_u = jnew_u
          yold_l = ynew_l
          yold_u = ynew_u
          xold  = xnew

          !  find where sides 2 and 3 cut next x-level

          i_next = i_next + 1
          xnew = xgrid(i_next)
          ynew_u = sidefunction ( pars, xnew, side_2 )
          ynew_l = sidefunction ( pars, xnew, side_3 )
          jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
          jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

          !  assign area

          CALL doubleside_slicer  &
               ( ygrid, nydim, nydim, side_3, side_2, pars,  &
               jold_l, jold_u, jnew_l, jnew_u, perpen, orient,  &
               xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
               area(1,is) )

          !  assign limits

          ylimit_lower(is) = jold_l
          ylimit_upper(is) = jold_u

       ENDDO

       !  last slice (including corner 3)
       !  -------------------------------

       !  update left slice boundary values; find y-limit for this corner

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold   = xnew
       is_last = is_end + 1
       jc2 = offset ( 1, ygrid, nydim, asct, cc(side_3,2) )

       CALL single_right_slicer  &
            ( ygrid, nydim, nydim, pars,  &
            jold_l, jc2, jold_u,  &
            xold, yold_l, yold_u, cc(side_3,1), cc(side_3,2),  &
            area(1,is_last) )

       !  assign limits

       ylimit_lower(is_last) = jold_l
       ylimit_upper(is_last) = jold_u

    ENDIF

    !  finish

  END SUBROUTINE tesselate_options_1

  !


  !  function parameters are stored in array pars which should
  !  be pre-computed. the modules do not know what function is being
  !  used in the area and line intercept determinations.

  SUBROUTINE tesselate_options_2  &
       ( xgrid, ygrid, nxdim, nydim,  &
       cc, orient, triple_first, triple_last, double_middle, pars,    &
       istart1, istart2, istart3, istart4, jfinis,  &
       area, ylimit_lower, ylimit_upper )

    !  input 
    !  =====

    !  grid values

    INTEGER		nxdim, nydim
    REAL*8		xgrid(nxdim), ygrid(nydim)

    !  corner coordinates

    REAL*8		cc(4,2)

    !  offsets

    INTEGER		istart1, istart2, istart4, istart3, jfinis

    !  function parametaers

    REAL*8		pars(6,4)

    !  control of options

    INTEGER		orient
    LOGICAL		triple_first, triple_last, double_middle

    !  output
    !  ======

    REAL*8		area(nydim,nxdim)
    INTEGER		ylimit_lower(nxdim)
    INTEGER		ylimit_upper(nxdim)

    !  local variables
    !  ===============

    REAL*8		xold, xnew, yold_l, yold_u, ynew_l, ynew_u
    INTEGER		is, is_next, is_first, is_last, is_beg, is_end
    INTEGER		jc1, jc2, jc3, jc4, jold_l, jold_u, jnew_l, jnew_u
    INTEGER		side_f, side_d, i_next, is_single, istart_dummy

    !  parameter values and indices

    INTEGER		side_1, side_2, side_3, side_4
    INTEGER		asct, parallel, perpen
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )
    PARAMETER	( asct = 1, parallel = 1, perpen = -1 )

    !  external functions

    !REAL*8		sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset
    !INTEGER		offset
    !EXTERNAL	offset, sidefunction

    !  quadruple option (all 4 corners in one slice)
    !  ################

    IF ( triple_first .AND. triple_last ) THEN

       !   orientation = +1, in order 1-4-2-3
       !   orientation = -1, in order 1-2-4-3

       IF ( orient .EQ. 1 ) THEN
          side_f = side_4
          side_d = side_2
          jc2 = 1
          jc3 = jfinis
       ELSE IF ( orient .EQ. -1 ) THEN
          side_f = side_2
          side_d = side_4
          jc2 = jfinis
          jc3 = 1
       ENDIF

       !  get area

       is_single = 1
       jc1 = offset ( 1, ygrid, nydim, asct, cc(side_1,2) )
       jc4 = offset ( 1, ygrid, nydim, asct, cc(side_3,2) )
       CALL quadruple_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jc1, jc2, jc3, jc4,  &
            cc(side_1,1), cc(side_1,2),  &
            cc(side_f,1), cc(side_f,2),  &
            cc(side_d,1), cc(side_d,2),  &
            cc(side_3,1), cc(side_3,2),  &
            area(1:nydim,is_single) )       ! tpk

       !  assign limits

       ylimit_lower(is_single) = 1
       ylimit_upper(is_single) = jfinis

       !  return after this option - no more to do

       RETURN

    ENDIF

    !  triple first option
    !  ###################

    IF ( triple_first ) THEN

       is_first = 1
       i_next = istart1 + 1
       is_next = is_first

       ! *******************************************************
       ! the following is a fix to i_next going out of bounds.
       ! it is a departure from the original program and may be
       ! introducing an error. it did not do so for the case it
       ! was discovered, but that doesn't mean it won't do so 
       ! at other places. of course, indices should not be left
       ! to their own devices in such a cavalrier manner in the
       ! first place.  tpk, 19 feb 2005.
       ! *******************************************************
       i_next = MIN ( i_next, nxdim ) !tpk

       xnew = xgrid(i_next)

       !  orientation = 1  (corners 1-4-2; upper/lower sides are 2 and 3)
       !  orientation = -1 (corners 1-2-4; upper/lower sides are 2 and 3)

       IF ( orient .EQ. 1 ) THEN
          side_f = side_4
          side_d = side_2
          jc2 = 1
          jc3 = jfinis
       ELSE IF ( orient .EQ. -1 ) THEN
          side_f = side_2
          side_d = side_4
          jc2 = jfinis
          jc3 = 1
       ENDIF

       !  find where upper/lower sides cut next x-level

       ynew_u = sidefunction ( pars, xnew, side_2 )
       ynew_l = sidefunction ( pars, xnew, side_3 )

       ! tpk additions
       ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
       ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )

       jnew_l  = offset ( 1, ygrid, nydim, asct, ynew_l )
       jc1 = offset ( 1, ygrid, nydim, asct, cc(side_1,2) )
       jnew_u  = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  get triple corner

       CALL triple_left_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jc1, jc2, jc3, jnew_l, jnew_u,  &
            xnew, ynew_l, ynew_u,  &
            cc(side_1,1), cc(side_1,2),  &
            cc(side_f,1), cc(side_f,2),  &
            cc(side_d,1), cc(side_d,2),  &
            area(1,is_first) )

       !  assign limits

       ylimit_lower(is_first) = 1
       ylimit_upper(is_first) = jfinis

    ENDIF

    !  single first option
    !  ###################

    IF ( .NOT. triple_first ) THEN

       !  first slice (includes corner 1)
       !  -------------------------------

       is_first = 1
       i_next = istart1 + 1
       is_next = is_first

       !  find where sides 1 and 4 cut next x-level, assign area and check

       xnew   = xgrid(i_next)
       ynew_u = sidefunction ( pars, xnew, side_1 )
       ynew_l = sidefunction ( pars, xnew, side_4 )
       jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
       jc1    = offset ( 1, ygrid, nydim, asct, cc(side_1,2) )
       jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )
       CALL single_left_slicer  &
            ( ygrid, nydim, nydim, pars,  &
            jnew_l, jc1, jnew_u,  &
            xnew, ynew_l, ynew_u,  &
            cc(side_1,1), cc(side_1,2),  &
            area(1,is_first) )

       !  assign limits

       ylimit_lower(is_first) = jnew_l
       ylimit_upper(is_first) = jnew_u

       !  slices between corners 1 and 2/4 (depends on orientation)
       !  --------------------------------

       !  how many x-slices 

       is_beg = is_first + 1
       IF ( orient.EQ.1 ) THEN
          istart_dummy = istart4
       ELSE IF ( orient.EQ.-1 ) THEN
          istart_dummy = istart2
       ENDIF
       is_end = istart_dummy - istart1

       !  start loop over x-slices (not done if is_beg > is_end)

       DO is = is_beg, is_end

          !  update left slice boundary values

          jold_l = jnew_l
          jold_u = jnew_u
          yold_l = ynew_l
          yold_u = ynew_u
          xold   = xnew

          !  find where sides borders cut next x-level

          i_next = i_next + 1
          IF ( i_next > nxdim ) CYCLE   !tpk
          xnew   = xgrid(i_next)
          ynew_u = sidefunction ( pars, xnew, side_1 )
          ynew_l = sidefunction ( pars, xnew, side_4 )
          jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
          jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

          !  assign area

          CALL doubleside_slicer  &
               ( ygrid, nydim, nydim, side_4, side_1, pars,  &
               jold_l, jold_u, jnew_l, jnew_u, perpen, orient,  &
               xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
               area(1,is) )

          !  assign limits

          ylimit_lower(is) = jnew_l
          ylimit_upper(is) = jnew_u

       ENDDO

    ENDIF

    !  double middle option
    !  ####################

    IF ( double_middle ) THEN

       !  slice including corners 2 and 4 (together)
       !  -------------------------------

       is_next = is_end + 1

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold  = xnew

       !  orientation = 1  (corners 4-2; upper/lower sides are 2 and 3)
       !  orientation = -1 (corners 2-4; upper/lower sides are 2 and 3)

       IF ( orient .EQ. 1 ) THEN
          side_f = side_4
          side_d = side_2
          jc1 = 1
          jc2 = jfinis
       ELSE IF ( orient .EQ. -1 ) THEN
          side_f = side_2
          side_d = side_4
          jc1 = jfinis
          jc2 = 1
       ENDIF

       !  find where upper/lower sides cut next x-level

       i_next = i_next + 1
       
       ! *******************************************************
       ! the following is a fix to i_next going out of bounds.
       ! it is a departure from the original program and may be
       ! introducing an error. it did not do so for the case it
       ! was discovered, but that doesn't mean it won't do so 
       ! at other places. of course, indices should not be left
       ! to their own devices in such a cavalrier manner in the
       ! first place.  tpk, 19 feb 2005 and 9 jan 2008.
       ! *******************************************************
       i_next  = MIN ( i_next,  nxdim ) !tpk
       is_next = MIN ( is_next, nxdim ) !tpk

       xnew   = xgrid(i_next)
       ynew_u = sidefunction ( pars, xnew, side_2 )
       ynew_l = sidefunction ( pars, xnew, side_3 )

       ! tpk additions
       ynew_l = MAX ( ynew_l, ygrid(1)+1.0d-05 )
       ynew_u = MAX ( ynew_u, ygrid(1)+1.0d-05 )

       jnew_l  = offset ( 1, ygrid, nydim, asct, ynew_l )
       jnew_u  = offset ( 1, ygrid, nydim, asct, ynew_u )

       !  assign area

       CALL double_middle_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jold_l, jold_u, jc1, jc2, jnew_l, jnew_u,  &
            xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
            cc(side_f,1), cc(side_f,2),  &
            cc(side_d,1), cc(side_d,2),  &
            area(1, is_next) )

       !  assign limits

       ylimit_lower(is_next) = 1
       ylimit_upper(is_next) = jfinis

    ENDIF

    !  triple last option
    !  ##################

    IF ( triple_last ) THEN

       !  final slice (includes corners 2/4 and 3)
       !  ----------------------------------------

       !  update left slice boundary values

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold  = xnew
       is_last = is_end + 1

       !  orientation = 1  (corners 4,2,3; upper/lower sides are 1 and 4)
       !  orientation = -1 (corners 2,4,3; upper/lower sides are 1 and 4)

       IF ( orient .EQ. 1 ) THEN
          side_f = side_4
          side_d = side_2
          jc1 = 1
          jc2 = jfinis
       ELSE IF ( orient .EQ. -1 ) THEN
          side_f = side_2
          side_d = side_4
          jc1 = jfinis
          jc2 = 1
       ENDIF

       !  corner limits

       jc3 = offset ( 1, ygrid, nydim, asct, cc(side_3,2) )

       CALL triple_right_slicer  &
            ( ygrid, nydim, nydim, orient, pars,  &
            jc1, jc2, jc3, jold_l, jold_u,  &
            xold, yold_l, yold_u,  &
            cc(side_f,1), cc(side_f,2),  &
            cc(side_d,1), cc(side_d,2),  &
            cc(side_3,1), cc(side_3,2),  &
            area(1, is_last) )

       !  assign limits

       ylimit_lower(is_last) = 1
       ylimit_upper(is_last) = jfinis

    ENDIF

    !  single last option
    !  ##################

    IF ( .NOT. triple_last ) THEN

       !  slices between corner 2 and corner 3
       !  ------------------------------------

       is_beg = is_next + 1
       is_end = istart3 - istart1

       DO is = is_beg, is_end

          !  update left slice boundary values

          jold_l = jnew_l
          jold_u = jnew_u
          yold_l = ynew_l
          yold_u = ynew_u
          xold  = xnew

          !  find where sides 2 and 3 cut next x-level

          i_next = i_next + 1
          xnew = xgrid(i_next)
          ynew_u = sidefunction ( pars, xnew, side_2 )
          ynew_l = sidefunction ( pars, xnew, side_3 )
          jnew_l = offset ( 1, ygrid, nydim, asct, ynew_l )
          jnew_u = offset ( 1, ygrid, nydim, asct, ynew_u )

          !  assign area

          CALL doubleside_slicer  &
               ( ygrid, nydim, nydim, side_3, side_2, pars,  &
               jold_l, jold_u, jnew_l, jnew_u, perpen, orient,  &
               xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
               area(1,is) )

          !  assign limits

          ylimit_lower(is) = jold_l
          ylimit_upper(is) = jold_u

       ENDDO

       !  last slice (including corner 3)
       !  -------------------------------

       !  update left slice boundary values; find y-limit for this corner

       jold_l = jnew_l
       jold_u = jnew_u
       yold_l = ynew_l
       yold_u = ynew_u
       xold   = xnew
       is_last = is_end + 1
       jc2 = offset ( 1, ygrid, nydim, asct, cc(side_3,2) )

       ! tpk additions
       IF ( is_last > 0 ) THEN
          CALL single_right_slicer                              &
               ( ygrid, nydim, nydim, pars,                          &
               jold_l, jc2, jold_u,                                  &
               xold, yold_l, yold_u, cc(side_3,1), cc(side_3,2),     &
               area(1,is_last) )

          !       assign limits

          ylimit_lower(is_last) = jold_l
          ylimit_upper(is_last) = jold_u
       END IF
    ENDIF

    !  finish

  END SUBROUTINE tesselate_options_2

  SUBROUTINE singleside_slicer  &
       ( yg, local_nydim, local_nybox, side, pars,  &
       jold, jnew, xold, yold, xnew, ynew,  &
       area )

    !  input/output

    INTEGER		local_nydim, local_nybox, side
    REAL*8		yg(local_nydim), pars(6,4)
    INTEGER		jold, jnew
    REAL*8		xold, xnew, yold, ynew
    REAL*8		area(local_nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4, j
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables (including external functions)

    REAL*8		x1, y1, x2, y2, quad, corner
    !REAL*8		inv_sidefunction, grid_area, corner_area
    !EXTERNAL	inv_sidefunction, grid_area, corner_area

    !  trivial case when j_old = j_new
    !  ===============================

    IF ( (jold.EQ.jnew .AND. jold > 0) .OR. jold == 1 ) THEN !tpk
       area(jold) = corner_area(pars,xold,xnew,yold,ynew,side)
       RETURN
    ENDIF
    IF ( jold == 0 ) RETURN !tpk

    !  for side 1
    !  ==========

    IF ( side .EQ. side_1 ) THEN
       !  first corner situation
       x1 = xold
       y1 = yold

       ! tpk additions
       IF ( jold+1 <= local_nydim ) THEN

          y2 = yg(jold+1)
          x2 = inv_sidefunction ( pars, y2, side )
          quad   = grid_area (pars,x2,xnew,y1,y2)
          corner = corner_area(pars,x1,x2,y1,y2,side)
          IF ( jold > 0 ) area(jold) = corner + quad  !tpk

          !       write(*, *) 'side1: ', corner, quad
          !       line sides situation
          DO j = jold + 1, jnew - 1
             IF ( j <= 0 ) CYCLE  !tpk
             y1 = y2
             x1 = x2
             y2 = yg(j+1)
             x2 = inv_sidefunction ( pars, y2, side )
             quad   = grid_area (pars,x2,xnew,y1,y2)
             corner = corner_area(pars,x1,x2,y1,y2,side)
             area(j) = corner + quad
          ENDDO
          !       last corner situation
          y1 = y2
          x1 = x2
          y2 = ynew
          x2 = xnew
          IF ( jnew > 0 ) area(jnew) = corner_area(pars,x1,x2,y1,y2,side) !tpk

       END IF ! tpk

       !  for side 2
       !  ==========

    ELSE IF ( side .EQ. side_2 ) THEN
       !  first corner situation
       x1 = xold
       y1 = yold
       y2 = yg(jold)
       x2 = inv_sidefunction ( pars, y2, side )
       area(jold) = corner_area(pars,x1,x2,y1,y2,side)
       !	  write(*, *) 'side2: ', area(jold)

       !  line sides situation
       DO j = jold - 1, jnew + 1, -1
          IF ( j <= 0 ) CYCLE  !tpk
          y1 = y2
          x1 = x2
          y2 = yg(j)
          x2 = inv_sidefunction ( pars, y2, side )
          quad   = grid_area (pars,xold,x1,y2,y1)
          corner = corner_area(pars,x1,x2,y1,y2,side)
          area(j) = corner + quad
       ENDDO
       !  last corner situation
       y1 = y2
       x1 = x2
       y2 = ynew
       x2 = xnew
       quad   = grid_area (pars,xold,x1,y2,y1)
       corner = corner_area(pars,x1,x2,y1,y2,side)
       IF ( jnew >= 1 .AND. jnew <= local_nybox ) &  !tpk
            area(jnew) = corner + quad

       !  for side 3
       !  ==========

    ELSE IF ( side .EQ. side_3 ) THEN
       !  first corner situation

       ! tpk additions
       IF ( jold > 0 .AND. jold+1 <= local_nydim ) THEN
          x1 = xold
          y1 = yold
          y2 = yg(jold+1)
          x2 = inv_sidefunction ( pars, y2, side )
          area(jold) = corner_area(pars,x1,x2,y1,y2,side)
          !	  write(*, *) 'side3: ', area(jold)
          !  line sides situation
          DO j = jold + 1, jnew - 1
             IF ( j <= 0 ) CYCLE  !tpk
             y1 = y2
             x1 = x2
             y2 = yg(j+1)
             x2 = inv_sidefunction ( pars, y2, side )
             quad   = grid_area (pars,xold,x1,y1,y2)
             corner = corner_area(pars,x1,x2,y1,y2,side)
             area(j) = corner + quad
          ENDDO
          !  last corner situation
          y1 = y2
          x1 = x2
          y2 = ynew
          x2 = xnew
          quad   = grid_area (pars,xold,x1,y1,y2)
          corner = corner_area(pars,x1,x2,y1,y2,side)
          area(jnew) = corner + quad
       END IF
       !  for side 4
       !  ==========

    ELSE IF ( side .EQ. side_4 .AND. jold > 0 ) THEN !tpk
       !  first corner situation
       x1 = xold
       y1 = yold
       y2 = yg(jold)
       x2 = inv_sidefunction ( pars, y2, side )
       quad   = grid_area (pars,x2,xnew,y2,y1)
       corner = corner_area(pars,x1,x2,y1,y2,side)
       area(jold) = corner + quad
       !	  write(*, *) 'side4: ', corner, quad
       !  line sides situation
       DO j = jold - 1, jnew + 1, -1
          IF ( j <= 0 ) CYCLE   !tpk
          y1 = y2
          x1 = x2
          y2 = yg(j)
          x2 = inv_sidefunction ( pars, y2, side )
          quad   = grid_area (pars,x2,xnew,y2,y1)
          corner = corner_area(pars,x1,x2,y1,y2,side)
          area(j) = corner + quad
       ENDDO
       !  last corner situation
       y1 = y2
       x1 = x2
       y2 = ynew
       x2 = xnew
       IF ( jnew > 0 ) &      ! tpk
            area(jnew) = corner_area(pars,x1,x2,y1,y2,side)

    ENDIF

    !  finish

  END SUBROUTINE singleside_slicer

  !

  SUBROUTINE doubleside_slicer  &
       ( yg, nydim, nybox, side_l, side_u, pars,  &
       jold_l, jold_u, jnew_l, jnew_u, parallel, ornt,  &
       xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
       area )

    !  input/output

    INTEGER		nybox, nydim, side_l, side_u, ornt, parallel
    REAL*8		yg(nydim), pars(6,4)
    INTEGER		jold_l, jold_u, jnew_l, jnew_u
    REAL*8		xold, xnew, yold_l, yold_u, ynew_l, ynew_u
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables (including external function)

    INTEGER		j, jmin, jmax, local_nybox
    PARAMETER	( local_nybox = 5000 )
    REAL*8		area_l(local_nybox), area_u(local_nybox)
    REAL*8		grid !, grid_area
    !EXTERNAL	grid_area

    !  initialise

    jmin = MAX(MIN ( jold_l, jnew_l ),          1) ! tpk
    jmax = MIN(MAX ( jold_u, jnew_u ),local_nybox) ! tpk
    DO j = jmin, jmax
       IF ( j <= 0 ) CYCLE  !tpk
       area_l(j) = 0.0
       area_u(j) = 0.0
    ENDDO

    !  side slicing
    !  ============

    !  slicing lower side (labelled l = either 3 or 4)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_l, pars,  &
         jold_l, jnew_l, xold, yold_l, xnew, ynew_l,  &
         area_l )
    !	  write(*, *) 'area_l = ', area_l(1:10)

    !  slicing upper side (labelled u = either 1 or 2)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_u, pars,  &
         jold_u, jnew_u, xold, yold_u, xnew, ynew_u,  &
         area_u )

    !	  write(*, *) 'area_u = ', area_u(1:10)
    !	  write(*, *) 'parallel = ', parallel

    !  parallel lines (or at least running together)
    !  =============================================

    IF ( parallel .EQ. 1 ) THEN

       !  orientation 1
       !  #############

       IF ( ornt .EQ. 1 ) THEN

          !  overlapping case

          IF ( jnew_l .GT. jold_u ) THEN

             DO j = jold_l, jold_u-1
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_l(j)
             ENDDO
             DO j = jnew_l+1, jnew_u
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_u(j)
             ENDDO

             ! tpk additions
             IF ( jnew_l < nydim .AND. MIN(jnew_l,jold_u) > 0 ) THEN
                grid = grid_area(pars,xold,xnew,ynew_l,yg(jnew_l+1))
                area_l(jnew_l) = area_l(jnew_l) + grid
                grid = grid_area(pars,xold,xnew,yg(jold_u),yold_u)
                area_u(jold_u) = area_u(jold_u) + grid
                grid = grid_area(pars,xold,xnew,yg(jnew_l),yg(jnew_l+1))
                area(jnew_l) = area_l(jnew_l) + area_u(jnew_l) - grid
                grid = grid_area(pars,xold,xnew,yg(jold_u),yg(jold_u+1))
                area(jold_u) = area_l(jold_u) + area_u(jold_u) - grid
                DO j = jold_u+1, jnew_l-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
                ENDDO
             END IF

             !  same grid cases

          ELSE IF ( jnew_l .EQ. jold_u ) THEN

             DO j = jold_l, jnew_l-1
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_l(j)
             ENDDO
             DO j = jold_u+1, jnew_u
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_u(j)
             ENDDO
             grid = grid_area(pars,xold,xnew,ynew_l,yold_u)

             IF ( jold_u > 0 .AND. jnew_l > 0 ) &  !tpk
                  area(jold_u) = area_l(jnew_l) + area_u(jold_u) + grid

             ! non-overlapping case :
             !    complete lower and upper, assigning directly and fill in grids

          ELSE

             IF ( jold_l > 0 .AND. jnew_l > 0 ) THEN !tpk
                DO j = jold_l, jnew_l-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = area_l(j)
                ENDDO
                DO j = jold_u+1, jnew_u
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = area_u(j)
                ENDDO
                grid = grid_area(pars,xold,xnew,ynew_l,yg(jnew_l+1))
                area(jnew_l) = area_l(jnew_l) + grid
                grid = grid_area(pars,xold,xnew,yg(jold_u),yold_u)
                area(jold_u) = area_u(jold_u) + grid
                DO j = jnew_l+1, jold_u-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
                ENDDO
             END IF !tpk

          ENDIF

          !  orientation -1
          !  ##############

       ELSE IF ( ornt .EQ. -1 ) THEN

          !  overlapping case

          !tpk additions
          !tpk if ( jnew_u .lt. jold_l ) then
          !tpk IF ( jnew_u .LT. jold_l .AND. jold_l < nydim ) THEN
          IF ( jnew_u .LT. jold_l .AND. jold_l < nydim &
               .AND. jnew_l > 0 .AND. jnew_u > 0 ) THEN   !tpk

             DO j = jold_l+1, jold_u
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_u(j)
             ENDDO
             DO j = jnew_l, jnew_u-1
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_l(j)
             ENDDO
             grid = grid_area(pars,xold,xnew,yold_l,yg(jold_l+1))
             area_l(jold_l) = area_l(jold_l) + grid
             grid = grid_area(pars,xold,xnew,yg(jnew_u),ynew_u)
             area_u(jnew_u) = area_u(jnew_u) + grid
             grid = grid_area(pars,xold,xnew,yg(jold_l),yg(jold_l+1))
             area(jold_l) = area_l(jold_l) + area_u(jold_l) - grid
             grid = grid_area(pars,xold,xnew,yg(jnew_u),yg(jnew_u+1))
             area(jnew_u) = area_l(jnew_u) + area_u(jnew_u) - grid
             DO j = jnew_u+1, jold_l-1
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
             ENDDO

             !  same grid cases

          ELSE IF ( jnew_u .EQ. jold_l ) THEN

             DO j = jnew_l, jold_l-1
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_l(j)
             ENDDO
             DO j = jnew_u+1, jold_u
                IF ( j <= 0 ) CYCLE  !tpk
                area(j) = area_u(j)
             ENDDO
             grid = grid_area(pars,xold,xnew,yold_l,ynew_u)
             area(jold_l) = area_l(jold_l) + area_u(jold_l) + grid

             ! non-overlapping case :
             !    complete lower and upper, assigning directly and fill in grids

          ELSE


             ! tpk addtions
             IF ( jold_l < nydim .AND. jnew_u > 0 .AND. jold_l > 0 ) THEN
                DO j = jnew_l, jold_l-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = area_l(j)
                ENDDO
                DO j = jnew_u+1, jold_u
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = area_u(j)
                ENDDO
                grid = grid_area(pars,xold,xnew,yold_l,yg(jold_l+1))
                area(jold_l) = area_l(jold_l) + grid
                grid = grid_area(pars,xold,xnew,yg(jnew_u),ynew_u)
                area(jnew_u) = area_u(jnew_u) + grid
                DO j = jold_l+1, jnew_u-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
                ENDDO
             END IF

          ENDIF

       ENDIF

       !  non-parallel cases
       !  ==================

    ELSE IF ( parallel .EQ. -1 ) THEN

       !  lower 4 upper 1
       !  ###############

       IF ( side_l .EQ. side_4 ) THEN

          !  assignation

          DO j = jnew_l, jold_l-1 
             IF ( j <= 0 ) CYCLE   !tpk
             area(j) = area_l(j)
          ENDDO
          DO j = jold_u+1, jnew_u
             IF ( j <= 0 ) CYCLE  !tpk
             area(j) = area_u(j)
          ENDDO

          !  add grids as required

          IF ( jold_l .EQ. jold_u .AND. jold_l > 0 ) THEN

             grid = grid_area(pars,xold,xnew,yold_l,yold_u)
             area(jold_l) = area_l(jold_l) + area_u(jold_l) + grid

          ELSE
             IF ( jold_l > 0 .AND. jold_l < nydim .AND. jold_u > 0 ) THEN  !tpk
                grid = grid_area(pars,xold,xnew,yold_l,yg(jold_l+1))
                area(jold_l) = area_l(jold_l) + grid
                grid = grid_area(pars,xold,xnew,yg(jold_u),yold_u)
                area(jold_u) = area_u(jold_u) + grid
                DO j = jold_l+1, jold_u-1
                   IF ( j <= 0 ) CYCLE   !tpk
                   area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
                ENDDO
             END IF

          ENDIF

          !  lower 3 upper 2
          !  ###############

       ELSE IF ( side_l .EQ. side_3 ) THEN

          !  assignation

          DO j = jold_l, jnew_l-1
             IF ( j <= 0 ) CYCLE  !tpk
             area(j) = area_l(j)
          ENDDO
          DO j = jnew_u+1, jold_u
          IF ( j <= 0 ) CYCLE  !tpk
             area(j) = area_u(j)
          ENDDO

          !  add grids as required

          IF ( jnew_l .EQ. jnew_u ) THEN

             IF ( jnew_l > 0 .AND. jnew_l <= nydim ) THEN    !tpk
                grid = grid_area(pars,xold,xnew,ynew_l,ynew_u)
                area(jnew_l) = area_l(jnew_l) + area_u(jnew_l) + grid
             END IF

          ELSE

             IF ( jnew_l > 0 .AND. jnew_u > 0 .AND. jnew_l < nydim ) THEN   !tpk
                grid = grid_area(pars,xold,xnew,ynew_l,yg(jnew_l+1))
                area(jnew_l) = area_l(jnew_l) + grid
                grid = grid_area(pars,xold,xnew,yg(jnew_u),ynew_u)
                area(jnew_u) = area_u(jnew_u) + grid
                DO j = jnew_l+1, jnew_u-1
                   IF ( j <= 0 ) CYCLE  !tpk
                   area(j) = grid_area(pars,xold,xnew,yg(j),yg(j+1))
                ENDDO
             END IF  !tpk

          ENDIF

       ENDIF

    ENDIF

    !  finish

  END SUBROUTINE doubleside_slicer

  !	  

  SUBROUTINE single_left_slicer  &
       ( yg, nydim, nybox, pars,  &
       jnew_l, jcnr, jnew_u,  &
       xnew, ynew_l, ynew_u, xcnr, ycnr,  &
       area )

    !  input/output

    INTEGER		nydim, nybox
    REAL*8		yg(nydim), pars(6,4)
    INTEGER		jcnr, jnew_l, jnew_u
    REAL*8		xnew,  ynew_l, ynew_u, xcnr, ycnr
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox
    PARAMETER	( local_nybox = 5000 )
    REAL*8		area_l(local_nybox), area_u(local_nybox)

    !  side slicing from corner to next lower boundary (side = 4)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_4, pars,  &
         jcnr, jnew_l, xcnr, ycnr, xnew, ynew_l,  &
         area_l )
    DO j = jnew_l, jcnr - 1
       IF ( j <= 0 ) CYCLE  !tpk
       area(j) = area_l(j)
    ENDDO

    !  side slicing from corner to next upper boundary (side = 1)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_1, pars,  &
         jcnr, jnew_u, xcnr, ycnr, xnew, ynew_u,  &
         area_u )
    DO j = jcnr+1, jnew_u
       IF ( j <= 0 ) CYCLE  !tpk
       area(j) = area_u(j)
    ENDDO

    !  add contributions in corner square

    IF ( jcnr > 0 ) area(jcnr) = area_l(jcnr) + area_u(jcnr) !tpk

    !  finish

  END SUBROUTINE single_left_slicer

  !	  

  SUBROUTINE single_right_slicer  &
       ( yg, nydim, nybox, pars,  &
       jold_l, jcnr, jold_u,  &
       xold, yold_l, yold_u, xcnr, ycnr,  &
       area )

    !  input/output

    INTEGER		nydim, nybox
    REAL*8		yg(nydim), pars(6,4)
    INTEGER		jcnr, jold_l, jold_u
    REAL*8		xold,  yold_l, yold_u, xcnr, ycnr
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox
    PARAMETER	( local_nybox = 5000 )
    REAL*8		area_l(local_nybox), area_u(local_nybox)

    !  side slicing from lower boundary to corner (side = 3)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_3, pars,  &
         jold_l, jcnr, xold, yold_l, xcnr, ycnr,  &
         area_l )
    DO j = jold_l, jcnr - 1
       IF ( j <= 0 ) CYCLE  !tpk
       area(j) = area_l(j)
    ENDDO

    !  side slicing from upper boundary to corner (side = 2)

    CALL singleside_slicer  &
         ( yg, nydim, local_nybox, side_2, pars,  &
         jold_u, jcnr, xold, yold_u, xcnr, ycnr,  &
         area_u )
    DO j = jcnr+1, jold_u
       IF ( j <= 0 ) CYCLE  !tpk
       area(j) = area_u(j)
    ENDDO

    !  add contributions in corner square
    IF ( jcnr > 0 .AND. jcnr <= nybox )  &      !tpk
         area(jcnr) = area_l(jcnr) + area_u(jcnr)

    !  finish

  END SUBROUTINE single_right_slicer

  !

  SUBROUTINE single_topbot_slicer  &
       ( yg, nydim, nybox, side_a, pars, do_top, ornt,  &
       jold_l, jold_u, jcnr, jnew_l, jnew_u,  &
       xold, xnew, yold_l, yold_u, ynew_l, ynew_u, xcnr, ycnr,  &
       area )

    !  input/output

    INTEGER		nydim, nybox, side_a, do_top, ornt
    REAL*8		yg(nydim), pars(6,4)
    INTEGER		jold_l, jold_u, jcnr, jnew_l, jnew_u
    REAL*8		xold, xnew, yold_l, yold_u, ynew_l, ynew_u
    REAL*8		xcnr, ycnr
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables (including external functions)

    INTEGER		j, jmax, jmin, local_nybox, asct
    PARAMETER	( local_nybox = 5000, asct = +1 )
    INTEGER		jnew_l_f, jnew_u_f, jold_l_d, jold_u_d
    REAL*8		area_f(local_nybox), area_d(local_nybox)
    REAL*8		xnew_f, ynew_l_f, ynew_u_f
    REAL*8		xold_d, yold_l_d, yold_u_d !, sidefunction
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    !  bottom slicer
    !  #############

    IF ( do_top .EQ. -1 ) THEN

       !c  initialise

       jmax = MAX(jold_u,jnew_u)
       DO j = jcnr, jmax
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice ( lower side = 4 )

       xnew_f   = xcnr
       ynew_l_f = ycnr
       ynew_u_f = sidefunction ( pars, xnew_f, side_a )
       jnew_l_f = jcnr
       jnew_u_f = offset ( 1, yg, nydim, asct, ynew_u_f )

       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_a, pars,  &
            jold_l, jold_u, jnew_l_f, jnew_u_f, -ornt, ornt,  &
            xold, xnew_f, yold_l, yold_u, ynew_l_f, ynew_u_f,  &
            area_f )

       !  second slice ( lower side = 3 )

       jold_l_d = jcnr
       jold_u_d = jnew_u_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_a, pars,  &
            jold_l_d, jold_u_d, jnew_l, jnew_u, ornt, ornt,  &
            xold_d, xnew, yold_l_d, yold_u_d, ynew_l, ynew_u,  &
            area_d )

       !  put the areas together (straightforward)

       jmax = MAX(jold_u,jnew_u)
       DO j = jcnr, jmax
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

    ENDIF

    !  top slicer
    !  ##########

    IF ( do_top .EQ. +1 ) THEN

       !c  initialise

       jmin = MIN(jold_l,jnew_l)
       DO j = jmin, jcnr
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice ( upper side = 1 )

       xnew_f   = xcnr
       ynew_u_f = ycnr
       ynew_l_f = sidefunction ( pars, xnew_f, side_a )
       jnew_u_f = jcnr
       jnew_l_f = offset ( 1, yg, nydim, asct, ynew_l_f ) 
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_a, side_1, pars,  &
            jold_l, jold_u, jnew_l_f, jnew_u_f, ornt, ornt,  &
            xold, xnew_f, yold_l, yold_u, ynew_l_f, ynew_u_f,  &
            area_f )

       !  second slice ( upper side = 2 )

       jold_u_d = jcnr
       jold_l_d = jnew_l_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_a, side_2, pars,  &
            jold_l_d, jold_u_d, jnew_l, jnew_u, -ornt, ornt,  &
            xold_d, xnew, yold_l_d, yold_u_d, ynew_l, ynew_u,  &
            area_d )

       !  put the areas together (straightforward)

       jmin = MIN(jold_l,jnew_l)
       DO j = jmin, jcnr
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE single_topbot_slicer

  !	  

  SUBROUTINE double_left_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jnew_l, jcnr_f, jcnr_d, jnew_u, xnew, ynew_l, ynew_u,  &
       xcnr_f, ycnr_f, xcnr_d, ycnr_d,  &
       area )

    !  input/output

    INTEGER		nydim, nybox, ornt
    REAL*8		yg(nydim)
    INTEGER		jnew_l, jcnr_f, jcnr_d, jnew_u
    REAL*8		xnew, ynew_l, ynew_u, pars(6,4)
    REAL*8		xcnr_f, ycnr_f, xcnr_d, ycnr_d
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    INTEGER		jnew_l_f, jnew_u_f, jold_l_d, jold_u_d
    REAL*8		area_f(local_nybox), area_d(local_nybox)
    REAL*8		xnew_f, ynew_l_f, ynew_u_f !, sidefunction
    REAL*8		xold_d, yold_l_d, yold_u_d
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jcnr_d, jnew_u
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice containing corner

       xnew_f   = xcnr_d
       ynew_l_f = ycnr_d
       ynew_u_f = sidefunction ( pars, xnew_f, side_1 )
       jnew_l_f = jcnr_d
       jnew_u_f = offset(1,yg,nydim,ascending,ynew_u_f)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_f, jcnr_f, jnew_u_f,  &
            xnew_f, ynew_l_f, ynew_u_f, xcnr_f, ycnr_f,   &
            area_f )

       !  second  slice between lines

       jold_l_d = jcnr_d
       jold_u_d = jnew_u_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f

       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l_d, jold_u_d, jnew_l, jnew_u, ornt, ornt,  &
            xold_d, xnew, yold_l_d, yold_u_d, ynew_l, ynew_u,  &
            area_d )

       !  put the areas together (straightforward)

       DO j = jcnr_d, jnew_u
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jnew_l, jcnr_d
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice containing corner

       xnew_f   = xcnr_d
       ynew_u_f = ycnr_d
       ynew_l_f = sidefunction ( pars, xnew_f, side_4 )
       jnew_u_f = jcnr_d
       jnew_l_f = offset(1,yg,nydim,ascending,ynew_l_f)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_f, jcnr_f, jnew_u_f,  &
            xnew_f, ynew_l_f, ynew_u_f, xcnr_f, ycnr_f,  &
            area_f )

       !  second  slice between lines

       jold_u_d = jcnr_d
       jold_l_d = jnew_l_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l_d, jold_u_d, jnew_l, jnew_u, -ornt, ornt,  &
            xold_d, xnew, yold_l_d, yold_u_d, ynew_l, ynew_u,  &
            area_d )

       !  put the areas together (straightforward)

       DO j = jnew_l, jcnr_d
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE double_left_slicer

  !	  

  SUBROUTINE double_right_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jold_l, jcnr_f, jcnr_d, jold_u,  &
       xold, yold_l, yold_u, xcnr_f, ycnr_f, xcnr_d, ycnr_d,  &
       area )

    !  input/output

    INTEGER		nydim, nybox, ornt
    REAL*8		yg(nydim)
    INTEGER		jold_l, jcnr_f, jcnr_d, jold_u
    REAL*8		xold, yold_l, yold_u, pars(6,4)
    REAL*8		xcnr_f, ycnr_f, xcnr_d, ycnr_d
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    INTEGER		jnew_l_f, jnew_u_f, jold_l_d, jold_u_d
    REAL*8		area_f(local_nybox), area_d(local_nybox)
    REAL*8		xnew_f, ynew_l_f, ynew_u_f !, sidefunction
    REAL*8		xold_d, yold_l_d, yold_u_d
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jold_l, jcnr_f
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice for lines 1 and 3

       xnew_f   = xcnr_f
       ynew_u_f = ycnr_f
       ynew_l_f = sidefunction ( pars, xnew_f, side_3 )
       jnew_u_f = jcnr_f
       jnew_l_f = offset(1,yg,nydim,ascending,ynew_l_f)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l, jold_u, jnew_l_f, jnew_u_f, ornt, ornt,  &
            xold, xnew_f, yold_l, yold_u, ynew_l_f, ynew_u_f,  &
            area_f )

       !  second  slice to final corner (lines 2 and 3)

       jold_l_d = jnew_l_f
       jold_u_d = jnew_u_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_d, jcnr_d, jold_u_d,  &
            xold_d, yold_l_d, yold_u_d, xcnr_d, ycnr_d,  &
            area_d )

       !  put the areas together (straightforward)

       DO j = jold_l, jcnr_f
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jcnr_f, jold_u
          IF ( j <= 0 ) CYCLE  !tpk
          area_f(j) = 0.0
          area_d(j) = 0.0
       ENDDO

       !  first slice for lines 2 and 4

       xnew_f   = xcnr_f
       ynew_l_f = ycnr_f
       ynew_u_f = sidefunction ( pars, xnew_f, side_2 )
       jnew_l_f = jcnr_f
       jnew_u_f = offset(1,yg,nydim,ascending,ynew_u_f)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l, jold_u, jnew_l_f, jnew_u_f, -ornt, ornt,  &
            xold, xnew_f, yold_l, yold_u, ynew_l_f, ynew_u_f,  &
            area_f )

       !  second  slice to final corner (lines 2 and 3)

       jold_l_d = jnew_l_f
       jold_u_d = jnew_u_f
       xold_d   = xnew_f
       yold_l_d = ynew_l_f
       yold_u_d = ynew_u_f
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_d, jcnr_d, jold_u_d,  &
            xold_d, yold_l_d, yold_u_d, xcnr_d, ycnr_d,  &
            area_d )

       !  put the areas together (straightforward)

       DO j = jcnr_f, jold_u
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_f(j) + area_d(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE double_right_slicer

  !	  

  SUBROUTINE double_middle_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jold_l, jold_u, jcnr_f, jcnr_d, jnew_l, jnew_u,   &
       xold, xnew, yold_l, yold_u, ynew_l, ynew_u,  &
       xcnr_f, ycnr_f, xcnr_d, ycnr_d,  &
       area )

    !  input/output

    INTEGER		nydim, nybox, ornt
    REAL*8		yg(nydim)
    INTEGER		jcnr_f, jcnr_d
    INTEGER		jold_l, jnew_l, jold_u, jnew_u
    REAL*8		xold, yold_l, yold_u
    REAL*8		xcnr_f, ycnr_f, xcnr_d, ycnr_d
    REAL*8		xnew, ynew_l, ynew_u, pars(6,4)
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    REAL*8		area_s1(local_nybox), area_s2(local_nybox)
    REAL*8		area_s3(local_nybox) !, sidefunction
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    INTEGER		jnew_l_s1, jnew_u_s1, jold_l_s1, jold_u_s1
    INTEGER		jnew_l_s2, jnew_u_s2, jold_l_s2, jold_u_s2
    INTEGER		jnew_l_s3, jnew_u_s3, jold_l_s3, jold_u_s3

    REAL*8		xold_s1, yold_l_s1, yold_u_s1
    REAL*8		xnew_s1, ynew_l_s1, ynew_u_s1
    REAL*8		xold_s2, yold_l_s2, yold_u_s2
    REAL*8		xnew_s2, ynew_l_s2, ynew_u_s2
    REAL*8		xold_s3, yold_l_s3, yold_u_s3
    REAL*8		xnew_s3, ynew_l_s3, ynew_u_s3

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jcnr_f, jcnr_d
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (left edge to first corner)

       xold_s1   = xold
       yold_l_s1 = yold_l
       yold_u_s1 = yold_u
       jold_l_s1 = jold_l
       jold_u_s1 = jold_u
       xnew_s1   = xcnr_f
       ynew_l_s1 = ycnr_f
       ynew_u_s1 = sidefunction ( pars, xnew_s1, side_1 )
       jnew_l_s1 = jcnr_f
       jnew_u_s1 = offset(1,yg,nydim,ascending,ynew_u_s1)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_1, pars,  &
            jold_l_s1, jold_u_s1, jnew_l_s1, jnew_u_s1, -ornt, ornt,  &
            xold_s1, xnew_s1, yold_l_s1,  &
            yold_u_s1, ynew_l_s1, ynew_u_s1,  &
            area_s1 )

       !  second slice (first corner to second corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xcnr_d
       ynew_u_s2 = ycnr_d
       ynew_l_s2 = sidefunction ( pars, xnew_s2, side_3 )
       jnew_u_s2 = jcnr_d
       jnew_l_s2 = offset(1,yg,nydim,ascending,ynew_l_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (second corner to right edge)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       xnew_s3   = xnew
       ynew_u_s3 = ynew_u
       ynew_l_s3 = ynew_l
       jnew_u_s3 = jnew_u
       jnew_l_s3 = jnew_l
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_2, pars,  &
            jold_l_s3, jold_u_s3, jnew_l_s3, jnew_u_s3, -ornt, ornt,  &
            xold_s3, xnew_s3, yold_l_s3,  &
            yold_u_s3, ynew_l_s3, ynew_u_s3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jcnr_f, jcnr_d
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jcnr_d, jcnr_f
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (left edge to first corner)

       xold_s1   = xold
       yold_l_s1 = yold_l
       yold_u_s1 = yold_u
       jold_l_s1 = jold_l
       jold_u_s1 = jold_u
       xnew_s1   = xcnr_f
       ynew_u_s1 = ycnr_f
       ynew_l_s1 = sidefunction ( pars, xnew_s1, side_4 )
       jnew_u_s1 = jcnr_f
       jnew_l_s1 = offset(1,yg,nydim,ascending,ynew_l_s1)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_1, pars,  &
            jold_l_s1, jold_u_s1, jnew_l_s1, jnew_u_s1, ornt, ornt,  &
            xold_s1, xnew_s1, yold_l_s1,  &
            yold_u_s1, ynew_l_s1, ynew_u_s1,  &
            area_s1 )

       !  second slice (first corner to second corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xcnr_d
       ynew_l_s2 = ycnr_d
       ynew_u_s2 = sidefunction ( pars, xnew_s2, side_2 )
       jnew_l_s2 = jcnr_d
       jnew_u_s2 = offset(1,yg,nydim,ascending,ynew_u_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, -ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (second corner to right edge)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       xnew_s3   = xnew
       ynew_u_s3 = ynew_u
       ynew_l_s3 = ynew_l
       jnew_u_s3 = jnew_u
       jnew_l_s3 = jnew_l
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_2, pars,  &
            jold_l_s3, jold_u_s3, jnew_l_s3, jnew_u_s3, ornt, ornt,  &
            xold_s3, xnew_s3, yold_l_s3,  &
            yold_u_s3, ynew_l_s3, ynew_u_s3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jcnr_d, jcnr_f
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE double_middle_slicer

  !	  

  SUBROUTINE triple_left_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jc1, jc2, jc3, jnew_l, jnew_u,  &
       xnew, ynew_l, ynew_u, xc1, yc1, xc2, yc2, xc3, yc3,  &
       area )

    !  order of corners is 1-4-2 (orient=1), 1-2-4 (orient=-1)

    !  input/output

    INTEGER		nydim, nybox
    REAL*8		yg(nydim)
    INTEGER		jc1, jc2, jc3, ornt
    INTEGER		jnew_l, jnew_u
    REAL*8		xnew, ynew_l, ynew_u, pars(6,4)
    REAL*8		xc1, yc1, xc2, yc2, xc3, yc3
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    REAL*8		area_s1(local_nybox), area_s2(local_nybox)
    REAL*8		area_s3(local_nybox) !, sidefunction
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    INTEGER		jnew_l_s1, jnew_u_s1
    INTEGER		jnew_l_s2, jnew_u_s2, jold_l_s2, jold_u_s2
    INTEGER		jnew_l_s3, jnew_u_s3, jold_l_s3, jold_u_s3
    REAL*8		xnew_s1, ynew_l_s1, ynew_u_s1
    REAL*8		xold_s2, yold_l_s2, yold_u_s2
    REAL*8		xnew_s2, ynew_l_s2, ynew_u_s2
    REAL*8		xold_s3, yold_l_s3, yold_u_s3
    REAL*8		xnew_s3, ynew_l_s3, ynew_u_s3

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jc2, jc3
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (first corner to second corner)

       xnew_s1   = xc2
       ynew_l_s1 = yc2
       ynew_u_s1 = sidefunction ( pars, xnew_s1, side_1 )
       jnew_l_s1 = jc2
       jnew_u_s1 = offset(1,yg,nydim,ascending,ynew_u_s1)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_s1, jc1, jnew_u_s1,  &
            xnew_s1, ynew_l_s1, ynew_u_s1, xc1, yc1,  &
            area_s1 )

       !  second slice (second corner to third corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc3
       ynew_u_s2 = yc3
       ynew_l_s2 = sidefunction ( pars, xnew_s2, side_3 )
       jnew_u_s2 = jc3
       jnew_l_s2 = offset(1,yg,nydim,ascending,ynew_l_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (third corner to right edge)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       xnew_s3   = xnew
       ynew_u_s3 = ynew_u
       ynew_l_s3 = ynew_l
       jnew_u_s3 = jnew_u
       jnew_l_s3 = jnew_l
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_2, pars,  &
            jold_l_s3, jold_u_s3, jnew_l_s3, jnew_u_s3, -ornt, ornt,  &
            xold_s3, xnew_s3, yold_l_s3,  &
            yold_u_s3, ynew_l_s3, ynew_u_s3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc2, jc3
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jc3, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (first corner to second corner)

       xnew_s1   = xc2
       ynew_u_s1 = yc2
       ynew_l_s1 = sidefunction ( pars, xnew_s1, side_4 )
       jnew_u_s1 = jc2
       jnew_l_s1 = offset(1,yg,nydim,ascending,ynew_l_s1)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_s1, jc1, jnew_u_s1,  &
            xnew_s1, ynew_l_s1, ynew_u_s1, xc1, yc1,  &
            area_s1 )

       !  second slice (second corner to third corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc3
       ynew_l_s2 = yc3
       ynew_u_s2 = sidefunction ( pars, xnew_s2, side_2 )
       jnew_l_s2 = jc3
       jnew_u_s2 = offset(1,yg,nydim,ascending,ynew_u_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, -ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (third corner to right edge)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       xnew_s3   = xnew
       ynew_u_s3 = ynew_u
       ynew_l_s3 = ynew_l
       jnew_u_s3 = jnew_u
       jnew_l_s3 = jnew_l
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_2, pars,  &
            jold_l_s3, jold_u_s3, jnew_l_s3, jnew_u_s3, ornt, ornt,  &
            xold_s3, xnew_s3, yold_l_s3,  &
            yold_u_s3, ynew_l_s3, ynew_u_s3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc3, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE triple_left_slicer

  !  

  SUBROUTINE triple_right_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jc1, jc2, jc3, jold_l, jold_u,  &
       xold, yold_l, yold_u, xc1, yc1, xc2, yc2, xc3, yc3,  &
       area )

    !  order of corners is 4-2-3 (orient=1), 2-4-3 (orient=-1)

    !  input/output

    INTEGER		nydim, nybox
    REAL*8		yg(nydim)
    INTEGER		jc1, jc2, jc3, ornt
    INTEGER		jold_l, jold_u
    REAL*8		xold, yold_l, yold_u, pars(6,4)
    REAL*8		xc1, yc1, xc2, yc2, xc3, yc3
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    REAL*8		area_s1(local_nybox), area_s2(local_nybox)
    REAL*8		area_s3(local_nybox) !, sidefunction
    !EXTERNAL	offset, sidefunction
    !INTEGER (KIND=i4), EXTERNAL      :: offset

    INTEGER		jnew_l_s1, jnew_u_s1, jold_l_s1, jold_u_s1
    INTEGER		jnew_l_s2, jnew_u_s2, jold_l_s2, jold_u_s2
    INTEGER		jold_l_s3, jold_u_s3
    REAL*8		xold_s1, yold_l_s1, yold_u_s1
    REAL*8		xnew_s1, ynew_l_s1, ynew_u_s1
    REAL*8		xold_s2, yold_l_s2, yold_u_s2
    REAL*8		xnew_s2, ynew_l_s2, ynew_u_s2
    REAL*8		xold_s3, yold_l_s3, yold_u_s3

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jc1, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (left edge to first corner)

       xold_s1   = xold
       yold_l_s1 = yold_l
       yold_u_s1 = yold_u
       jold_l_s1 = jold_l
       jold_u_s1 = jold_u
       xnew_s1   = xc1
       ynew_l_s1 = yc1
       ynew_u_s1 = sidefunction ( pars, xnew_s1, side_1 )
       jnew_l_s1 = jc1
       jnew_u_s1 = offset(1,yg,nydim,ascending,ynew_u_s1)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_1, pars,  &
            jold_l_s1, jold_u_s1, jnew_l_s1, jnew_u_s1, -ornt, ornt,  &
            xold_s1, xnew_s1, yold_l_s1,  &
            yold_u_s1, ynew_l_s1, ynew_u_s1,  &
            area_s1 )

       !  second slice (first corner to second corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc2
       ynew_u_s2 = yc2
       ynew_l_s2 = sidefunction ( pars, xnew_s2, side_3 )
       jnew_u_s2 = jc2
       jnew_l_s2 = offset(1,yg,nydim,ascending,ynew_l_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (second corner to third corner)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_s3, jc3, jold_u_s3,  &
            xold_s3, yold_l_s3, yold_u_s3, xc3, yc3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc1, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jc2, jc1
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (left edge to first corner)

       xold_s1   = xold
       yold_l_s1 = yold_l
       yold_u_s1 = yold_u
       jold_l_s1 = jold_l
       jold_u_s1 = jold_u
       xnew_s1   = xc1
       ynew_u_s1 = yc1
       ynew_l_s1 = sidefunction ( pars, xnew_s1, side_4 )
       jnew_u_s1 = jc1
       jnew_l_s1 = offset(1,yg,nydim,ascending,ynew_l_s1)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_1, pars,  &
            jold_l_s1, jold_u_s1, jnew_l_s1, jnew_u_s1, ornt, ornt,  &
            xold_s1, xnew_s1, yold_l_s1,  &
            yold_u_s1, ynew_l_s1, ynew_u_s1,  &
            area_s1 )

       !  second slice (first corner to second corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc2
       ynew_l_s2 = yc2
       ynew_u_s2 = sidefunction ( pars, xnew_s2, side_2 )
       jnew_l_s2 = jc2
       jnew_u_s2 = offset(1,yg,nydim,ascending,ynew_u_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, -ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (second corner to third corner)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_s3, jc3, jold_u_s3,  &
            xold_s3, yold_l_s3, yold_u_s3, xc3, yc3,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc2, jc1
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE triple_right_slicer

  !	  

  SUBROUTINE quadruple_slicer  &
       ( yg, nydim, nybox, ornt, pars,  &
       jc1, jc2, jc3, jc4,  &
       xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4,  &
       area )

    !  order of corners is 1-4-2-3 (orient=1), 1-2-4-3 (orient=-1)

    !  input/output

    INTEGER		nydim, nybox
    REAL*8		yg(nydim), pars(6,4)
    INTEGER		jc1, jc2, jc3, jc4, ornt
    REAL*8		xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4
    REAL*8		area(nybox)

    !  side parameters

    INTEGER		side_1, side_2, side_3, side_4
    PARAMETER	( side_1 = 1, side_2 = 2 )
    PARAMETER	( side_3 = 3, side_4 = 4 )

    !  local variables

    INTEGER		j, local_nybox, ascending
    PARAMETER	( local_nybox = 5000, ascending = +1 )
    REAL*8		area_s1(local_nybox), area_s2(local_nybox)
    REAL*8		area_s3(local_nybox)
    !EXTERNAL	offset, sidefunction

    INTEGER		jnew_l_s1, jnew_u_s1
    INTEGER		jnew_l_s2, jnew_u_s2, jold_l_s2, jold_u_s2
    INTEGER		jold_l_s3, jold_u_s3
    REAL*8		xnew_s1, ynew_l_s1, ynew_u_s1
    REAL*8		xold_s2, yold_l_s2, yold_u_s2
    REAL*8		xnew_s2, ynew_l_s2, ynew_u_s2
    REAL*8		xold_s3, yold_l_s3, yold_u_s3

    !  orientation = 1
    !  ===============

    IF ( ornt .EQ. 1 ) THEN

       !  initialise

       DO j = jc2, jc3
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (first corner to second corner)

       xnew_s1   = xc2
       ynew_l_s1 = yc2
       ynew_u_s1 = sidefunction ( pars, xnew_s1, side_1 )
       jnew_l_s1 = jc2
       jnew_u_s1 = offset(1,yg,nydim,ascending,ynew_u_s1)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_s1, jc1, jnew_u_s1,  &
            xnew_s1, ynew_l_s1, ynew_u_s1, xc1, yc1,  &
            area_s1 )

       !  second slice (second corner to third corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc3
       ynew_u_s2 = yc3
       ynew_l_s2 = sidefunction ( pars, xnew_s2, side_3 )
       jnew_u_s2 = jc3
       jnew_l_s2 = offset(1,yg,nydim,ascending,ynew_l_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_3, side_1, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (third corner to fourth corner)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_s3, jc4, jold_u_s3,  &
            xold_s3, yold_l_s3, yold_u_s3, xc4, yc4,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc2, jc3
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

       !  orientation = -1
       !  ================

    ELSE IF ( ornt .EQ. -1 ) THEN

       !  initialise

       DO j = jc3, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area_s1(j) = 0.0
          area_s2(j) = 0.0
          area_s3(j) = 0.0
       ENDDO

       !  first slice (first corner to second corner)

       xnew_s1   = xc2
       ynew_u_s1 = yc2
       ynew_l_s1 = sidefunction ( pars, xnew_s1, side_4 )
       jnew_u_s1 = jc2
       jnew_l_s1 = offset(1,yg,nydim,ascending,ynew_l_s1)
       CALL single_left_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jnew_l_s1, jc1, jnew_u_s1,  &
            xnew_s1, ynew_l_s1, ynew_u_s1, xc1, yc1,  &
            area_s1 )

       !  second slice (second corner to third corner)

       xold_s2   = xnew_s1
       yold_l_s2 = ynew_l_s1
       yold_u_s2 = ynew_u_s1
       jold_l_s2 = jnew_l_s1
       jold_u_s2 = jnew_u_s1
       xnew_s2   = xc3
       ynew_l_s2 = yc3
       ynew_u_s2 = sidefunction ( pars, xnew_s2, side_2 )
       jnew_l_s2 = jc3
       jnew_u_s2 = offset(1,yg,nydim,ascending,ynew_u_s2)
       CALL doubleside_slicer  &
            ( yg, nydim, local_nybox, side_4, side_2, pars,  &
            jold_l_s2, jold_u_s2, jnew_l_s2, jnew_u_s2, -ornt, ornt,  &
            xold_s2, xnew_s2, yold_l_s2,  &
            yold_u_s2, ynew_l_s2, ynew_u_s2,  &
            area_s2 )

       !  third slice (third corner to fourth corner)

       xold_s3   = xnew_s2
       yold_l_s3 = ynew_l_s2
       yold_u_s3 = ynew_u_s2
       jold_l_s3 = jnew_l_s2
       jold_u_s3 = jnew_u_s2
       CALL single_right_slicer  &
            ( yg, nydim, local_nybox, pars,  &
            jold_l_s3, jc4, jold_u_s3,  &
            xold_s3, yold_l_s3, yold_u_s3, xc4, yc4,  &
            area_s3 )

       !  put the areas together (straightforward)

       DO j = jc3, jc2
          IF ( j <= 0 ) CYCLE  !tpk
          area(j) = area_s1(j) + area_s2(j) + area_s3(j)
       ENDDO

    ENDIF

    !  finish

  END SUBROUTINE quadruple_slicer

  !  tessellation_sphericalfunc.f
  !  ============================

  !  compilation of functions and modules required to carry out the
  !  tessellation algorithm in spherical geometry. comprises

  !	sidefunction
  !	great circle function expressing latitude on gc, given longitude

  !	inv_sidefunction
  !	great circle inverse function expressing longitude on gc, given
  !       latitude, allowing for phase differences, and taking care that
  !       longitude value lies between appropriate corner limits.

  !	corner_area
  !	returns the area of spherical triangle bound by one pixel side,
  !       one latitude ordinate and one longitude ordinate.

  !	grid_area
  !	returns the area of single latitude/longitude bin on a sphere.

  !       parameter_setup
  !       calculates the parameters for each of the 4 sides of the footprint
  !          pars(1,side) = k (see equation in text)
  !          pars(2,side) = phi_s (see equation in text)
  !          pars(3,side) = q (see text), only required for corner_area
  !          pars(4,side) = degrees to radains conversion (same for all sides)
  !          pars(5,side) = lower longitude limit (inv_sidefunction)
  !          pars(6,side) = upper longitude limit (inv_sidefunction)

  !        also returns the orientation of the pixel

  !    ********************************************************************
  !    *   robert spurr, november 1998					*
  !    *   sao, 60 garden street, cambridge, ma 02138, usa		*
  !    *   +1 (617) 496 7819; email rspurr@cfa.harvard.edu		*
  !    *									*
  !    *   algorithm may be subject of licensing agreement		*
  !    ********************************************************************

  !  parameter setup for great circles

  SUBROUTINE parsetup ( cc, pars, orient, yn_problem )

    !  input/output
    REAL*8		cc(4,2),pars(6,4)
    INTEGER		orient

    LOGICAL, INTENT (INOUT) :: yn_problem    !tpk

    !  local variables

    REAL*8		l,k,l1(4),l2(4),p1(4),p2(4),dtr
    REAL*8		a,b,c,q1,q2,ps,q,k0,l0
    INTEGER		side

    !  degrees to radians conversion

    dtr = datan(1.0d0)/45.0d0

    !  limiting latitude and longitude values at corners

    side = 1
    p2(side) = cc(2,1)
    p1(side) = cc(1,1)
    l2(side) = cc(2,2)
    l1(side) = cc(1,2)
    side = 2
    p2(side) = cc(3,1)
    p1(side) = cc(2,1)
    l2(side) = cc(3,2)
    l1(side) = cc(2,2)
    side = 3
    p2(side) = cc(4,1)
    p1(side) = cc(3,1)
    l2(side) = cc(4,2)
    l1(side) = cc(3,2)
    side = 4
    p2(side) = cc(1,1)
    p1(side) = cc(4,1)
    l2(side) = cc(1,2)
    l1(side) = cc(4,2)

    !  for each side

    DO side = 1, 4

       !  constants involved in great circle equation and corner area results

       a = dtan(l1(side)*dtr)
       b = dtan(l2(side)*dtr)
       c = dsin((p2(side)-p1(side))*dtr)
       IF ( c == 0.0_r8 ) THEN
          yn_problem = .TRUE.
          RETURN
       END IF
       q1 = -a*dcos(p2(side)*dtr) + b*dcos(p1(side)*dtr)
       q2 = -a*dsin(p2(side)*dtr) + b*dsin(p1(side)*dtr)
       k0 = dsqrt(q1*q1+q2*q2)/c
       k = k0
       ps = datan(q2/q1)/dtr
       l0 = datan(k*dsin((p1(side)-ps)*dtr))/dtr
       IF ( dabs((l1(side)/l0)-1).GT.1.0d0)k=-k0
       l = dsqrt(1.0d0+k*k)
       q = k/l

       !  assign parameters as listed above

       pars(1,side) = k
       pars(2,side) = ps
       pars(3,side) = q
       pars(4,side) = dtr
       pars(5,side) = MIN(p2(side),p1(side))
       pars(6,side) = MAX(p2(side),p1(side))

    ENDDO

    IF ( ANY(pars == 0.0_r8) ) yn_problem = .TRUE.

    !  orientation

    IF ( cc(4,1).LT.cc(2,1) ) THEN
       orient = +1
    ELSE
       orient = -1
    ENDIF

    !  finish

  END SUBROUTINE parsetup


  !  offset function
  INTEGER (KIND=i4) FUNCTION offset ( start, grid, ngrid, ascending, value )

    INTEGER		start, ngrid, n, ascending
    REAL*8		grid(ngrid), value
    LOGICAL		loop

    offset = 0

    !if ( ascending. eq.1 ) then
    !     offset = maxval( maxloc( grid(1:ngrid), mask=(grid(1:ngrid) <= value ) ) )
    !else
    !     offset = maxval( maxloc( grid(1:start), mask=(grid(1:start) <= value ) ) )
    !end if
    !return

    IF ( ascending .eq. 1 ) THEN
       offset = MAXVAL( MAXLOC( grid(1:ngrid), MASK=(grid(1:ngrid) < value ) ) )
    ELSE
       loop = .TRUE.
       n = start
       DO WHILE (loop)
          n = n - 1
          IF (value.GT.grid(n)) THEN
             loop = .FALSE.
             offset = n + 1
          ENDIF
       ENDDO
    ENDIF

  END FUNCTION offset

  !  great circle side function

  REAL*8 FUNCTION sidefunction ( pars, phi, side )
    REAL*8		pars(6,4), phi, dtr, help
    INTEGER		side
    dtr  = pars(4,side)
    help = dsin ( dtr * ( phi - pars(2,side) ) )
    sidefunction = datan ( pars(1,side) * help ) / dtr
  END FUNCTION sidefunction

  !  great circle inverse side function

  REAL*8 FUNCTION inv_sidefunction ( pars, lam, side )
    REAL*8		pars(6,4), lam, dtr, help, ps, x, xr
    INTEGER		side
    dtr  = pars(4,side)
    ps = pars(2,side)
    help = MIN(1.0_r8,MAX(-1.0_r8,dtan ( dtr * lam ) / pars(1,side))) ! tpk
    help = dasin(help)/dtr
    x = ps + help
    xr = x + 360.0d0

    !	write(*, *) '*** ', x, xr, side, pars(5, side), pars(6, side),
    !     * ps, help, lam, dtr, pars(1, side)

    IF ( x.GE.pars(5,side).AND.x.LE.pars(6,side) ) THEN
       inv_sidefunction = x
    ELSE IF ( xr.GE.pars(5,side).AND.xr.LE.pars(6,side) ) THEN
       inv_sidefunction = xr
       !  
    ELSE IF (x .LE. pars(5, side)) THEN
       inv_sidefunction = ps + 180.0 - help
    ELSE
       inv_sidefunction = ps - 180.0 - help
    ENDIF
  END FUNCTION inv_sidefunction

  !  corner area function

  REAL*8 FUNCTION corner_area (pars,x1,x2,y1,y2,side)
    REAL*8 		pars(6,4),x1,x2,y1,y2
    INTEGER		side
    REAL*8		dtr, q1, q2, astar, dx, corner
    dtr = pars(4,side)
    q1 = pars(3,side)*dcos((x1-pars(2,side))*dtr)
    q2 = pars(3,side)*dcos((x2-pars(2,side))*dtr)
    astar = dasin(q1)-dasin(q2)
    dx = ( x2 - x1 ) * dtr
    IF ( side.EQ.1 ) THEN
       corner = astar - dx*dsin(y1*dtr)
    ELSE IF  ( side.EQ.2 ) THEN
       corner = astar - dx*dsin(y2*dtr)
    ELSE IF  ( side.EQ.3 ) THEN
       corner = dx*dsin(y2*dtr) - astar
    ELSE IF  ( side.EQ.4 ) THEN
       corner = dx*dsin(y1*dtr) - astar
    ENDIF
    corner_area = dabs(corner) / dtr / dtr

    !	if (corner_area .gt. 1000.) then
    !		 print *, 'bad case', side
    !		 print *,  x1, x2, y1, y2, dtr, q1, q2, astar, dx, corner, corner_area
    !	endif
  END FUNCTION corner_area

  !  grid area function ( convention x2 > x1, y2 > y1 )

  REAL*8 FUNCTION grid_area (pars,x1,x2,y1,y2)
    REAL*8 		pars(6,4),x1,x2,y1,y2,dtr
    dtr = pars(4,1)
    grid_area = (x2-x1)*(dsin(y2*dtr)-dsin(y1*dtr))/dtr
  END FUNCTION grid_area



!  SUBROUTINE tesselations_patch1  &
!       ( maxruns, max_ydim, max_xdim, max_adim,  &
!       nruns, foot_coords, global_albedo_1,   &
!       etop05_climdata, oxgrid, oygrid,  &
!       foot_albedos, foot_heights, foot_areas )
!
!    !  arguments
!    !  =========
!
!    !  input numbers
!
!    INTEGER		max_adim, max_xdim, max_ydim
!    INTEGER		maxruns, nruns
!
!    !  latitude and longitudes + albedo set
!
!    REAL*8		oxgrid(max_xdim), oygrid(max_ydim)
!    REAL		global_albedo_1(max_ydim,max_adim)
!    INTEGER		etop05_climdata(max_ydim,max_adim)
!
!    !  input coordinates
!
!    DOUBLE PRECISION foot_coords ( 5, 2, maxruns )
!
!    !  output albedos
!
!    DOUBLE PRECISION foot_albedos ( maxruns )
!    DOUBLE PRECISION foot_heights ( maxruns )
!    DOUBLE PRECISION foot_areas   ( maxruns )
!
!    !  local variables
!    !  ===============
!
!    INTEGER		max_xbox, max_ybox
!    PARAMETER	( max_xbox = 50,   max_ybox = 50 )
!
!    !  geolocation input
!
!    REAL*8		geo_cc(5,2)
!
!    !  corner coordinates and side equation parameters
!
!    REAL*8		pars(6,4), cc(4,2)
!
!    !  limits output
!
!    INTEGER		istart1, istart2, istart3, istart4
!    INTEGER		jstart, jfinis, orient
!
!    !  box-grid offsets and limits
!
!    INTEGER		istart1_t, istart2_t, istart3_t, istart4_t
!    INTEGER		jstart_t, jfinis_t
!    INTEGER		nxused, nyused
!
!    !  box grids, tesselated areas and albedos in box
!
!    REAL*8		xgrid_t(max_xbox), ygrid_t(max_ybox)
!    REAL*8		area(max_ybox,max_xbox)
!    REAL*8		albused(max_ybox,max_xbox)
!    REAL*8		tpgused(max_ybox,max_xbox)
!
!    !  circularity flag
!
!    LOGICAL		circularity
!
!    !  special case flags
!
!    LOGICAL		special_side1
!    LOGICAL		special_side3
!
!    !  major output results (area sum and weighted albedo)
!
!    REAL*8		sum, albedo, height
!    INTEGER		ylimit_lower(max_xbox)
!    INTEGER		ylimit_upper(max_xbox)
!
!    !  other local variables
!
!    INTEGER		i, j, i1, j1, ig, jg, iz, rank(4), k, m,n
!    REAL*8		lima, minl, lim1, lim2, ep
!    LOGICAL		todo(4)
!    INTEGER		offset
!    !EXTERNAL	offset
!
!    !  nruns loop
!
!    DO n = 1, nruns
!
!       !  initialise output
!
!       foot_albedos(n) = 0.0d0
!       foot_areas(n)   = 0.0d0
!
!       !  set local variables
!
!       DO i = 1, 5
!          DO j = 1, 2
!             geo_cc(i,j) = foot_coords(i,j,n)
!          ENDDO
!       ENDDO
!
!       !  test for circularity (corners not yet ranked so must do generally)
!
!       lima = 180.0
!       circularity = .FALSE.
!       DO k = 1, 4
!          IF ( geo_cc(k,2) .GT. lima ) THEN
!             DO i = 1, 4
!                IF ( i.NE.k ) THEN
!                   IF ( dabs(geo_cc(i,2)-geo_cc(k,2)).GT.lima ) THEN
!                      geo_cc(i,2) = 360.0 + geo_cc(i,2)
!                      circularity = .TRUE.
!                   ENDIF
!                ENDIF
!             ENDDO
!          ENDIF
!       ENDDO
!
!       !  rank corners by increasing longitude, 1 = smallest.
!
!       DO k = 1, 4
!          todo(k) = .TRUE.
!       ENDDO
!       DO k = 1, 4
!          minl = +1000.0
!          DO i = 1, 4
!             IF ( todo(i) ) THEN
!                minl = MIN(minl,geo_cc(i,2))
!                IF ( minl .EQ. geo_cc(i,2) ) rank(k) = i
!             ENDIF
!          ENDDO
!          todo(rank(k)) = .FALSE.
!       ENDDO
!
!       ! assign ranked corner coordinates for use in tessellation algorithm
!
!       cc(1,1) = geo_cc(rank(1),2)
!       cc(1,2) = geo_cc(rank(1),1)
!       cc(3,1) = geo_cc(rank(4),2)
!       cc(3,2) = geo_cc(rank(4),1)
!       IF ( geo_cc(rank(2),1) .GT. cc(1,2) ) THEN
!          cc(2,1) = geo_cc(rank(2),2)
!          cc(2,2) = geo_cc(rank(2),1)
!          cc(4,1) = geo_cc(rank(3),2)
!          cc(4,2) = geo_cc(rank(3),1)
!       ELSE
!          cc(2,1) = geo_cc(rank(3),2)
!          cc(2,2) = geo_cc(rank(3),1)
!          cc(4,1) = geo_cc(rank(2),2)
!          cc(4,2) = geo_cc(rank(2),1)
!       ENDIF
!
!       !  exact values are not allowed, make small displacements
!
!       ep = 1.0d-05
!       IF (cc(1,1)-dfloat(INT(cc(1,1))).EQ.0.0 )cc(1,1)=cc(1,1)+ep
!       IF (cc(1,2)-dfloat(INT(cc(1,2))).EQ.0.0 )cc(1,2)=cc(1,2)+ep
!       IF (cc(2,1)-dfloat(INT(cc(2,1))).EQ.0.0 )cc(2,1)=cc(2,1)-ep
!       IF (cc(2,2)-dfloat(INT(cc(2,2))).EQ.0.0 )cc(2,2)=cc(2,2)-ep
!       IF (cc(3,1)-dfloat(INT(cc(3,1))).EQ.0.0 )cc(3,1)=cc(3,1)-ep
!       IF (cc(3,2)-dfloat(INT(cc(3,2))).EQ.0.0 )cc(3,2)=cc(3,2)-ep
!       IF (cc(4,1)-dfloat(INT(cc(4,1))).EQ.0.0 )cc(4,1)=cc(4,1)+ep
!       IF (cc(4,2)-dfloat(INT(cc(4,2))).EQ.0.0 )cc(4,2)=cc(4,2)+ep
!
!       !  tessellation setup
!       !  =============+====
!
!       !  find parameters for the footprint sides
!
!       CALL parsetup ( cc, pars, orient )
!
!       !  special cases (only flagged in this version, not computed)
!
!       special_side1 = .FALSE.
!       lim1 = pars(2,1) + 90.0
!       lim2 = pars(2,1) - 90.0
!       IF ( ( lim1.GT.cc(1,1).AND.lim1.LT.cc(2,1) )   .OR.  &
!            ( lim2.GT.cc(1,1).AND.lim2.LT.cc(2,1) ) ) THEN
!          special_side1 = .TRUE.
!       ENDIF
!       special_side3 = .FALSE.
!       lim1 = pars(2,3) + 90.0
!       lim2 = pars(2,3) - 90.0
!       IF ( ( lim1.GT.cc(4,1).AND.lim1.LT.cc(3,1) )   .OR.  &
!            ( lim2.GT.cc(4,1).AND.lim2.LT.cc(3,1) ) ) THEN
!          special_side3 = .TRUE.
!       ENDIF
!       IF ( special_side1 .OR. special_side3 ) THEN
!          WRITE(*,'(i5,l2,a)')m,circularity,' special cases'
!          GOTO 999
!       ENDIF
!
!       !  find corner offsets in the original x-y grids
!
!       istart1 = offset ( 1,       oxgrid, max_xdim, 1, cc(1,1) )
!       istart3 = offset ( istart1, oxgrid, max_xdim, 1, cc(3,1) )
!       jstart  = offset ( 1,       oygrid, max_ydim, 1, cc(4,2) )
!       jfinis  = offset ( 1,       oygrid, max_ydim, 1, cc(2,2) )
!       istart2 = offset ( istart1, oxgrid, max_xdim, 1, cc(2,1) )
!       istart4 = offset ( istart1, oxgrid, max_xdim, 1, cc(4,1) )
!
!       !  prepare albedo data in box around footprint
!       !  (  take account of the circularity at meridian crossing)
!
!       ig = 0
!       IF ( circularity ) THEN
!          DO i = istart1, max_adim
!             ig = ig + 1
!             jg = 0
!             DO j = jstart, jfinis
!                jg = jg + 1
!                albused(jg,ig) = global_albedo_1(j,i)
!                tpgused(jg,ig) = dfloat(etop05_climdata(j,i))/1000.0d0
!             ENDDO
!          ENDDO
!          DO i = max_adim+1, istart3
!             ig = ig + 1
!             jg = 0
!             DO j = jstart, jfinis
!                jg = jg + 1
!                iz = i-max_adim
!                albused(jg,ig) = global_albedo_1(j,iz)
!                tpgused(jg,ig) = dfloat(etop05_climdata(j,iz))/1000.0d0
!             ENDDO
!          ENDDO
!       ELSE
!          DO i = istart1, istart3
!             ig = ig + 1
!             jg = 0
!             DO j = jstart, jfinis
!                jg = jg + 1
!                albused(jg,ig) = global_albedo_1(j,i)
!                tpgused(jg,ig) = dfloat(etop05_climdata(j,i))/1000.0d0
!             ENDDO
!          ENDDO
!       ENDIF
!
!       !  define offsets in narrow box covering only the pixel
!
!       nxused = istart3 - istart1 + 2
!       nyused = jfinis  - jstart  + 2
!       istart3_t = nxused - 1
!       istart1_t = 1
!       istart2_t = istart2 - istart1 + 1
!       istart4_t = istart4 - istart1 + 1
!       jfinis_t  = nyused - 1
!       jstart_t  = 1
!
!       !  define tesselation grid box round footprint
!
!       DO j = 1, nyused
!          j1 = j + jstart - 1
!          ygrid_t(j) = oygrid(j1)
!       ENDDO
!       DO i = 1, nxused
!          i1 = i + istart1 - 1
!          xgrid_t(i) = oxgrid(i1)
!       ENDDO
!
!       !  single call to tessellation master routine
!
!       CALL tesselate_areamaster  &
!            ( xgrid_t, ygrid_t, max_xbox, max_ybox,  &
!            cc, pars, orient,   &
!            istart1_t, istart2_t, istart3_t, istart4_t, jfinis_t,  &
!            area, sum, ylimit_lower, ylimit_upper )
!
!       !  assign output
!
!       albedo = 0.0d0
!       height = 0.0d0
!       DO i = istart1_t, istart3_t
!          DO j = ylimit_lower(i), ylimit_upper(i)
!             albedo = albedo + area(j,i)*albused(j,i)
!             height = height + area(j,i)*tpgused(j,i)
!          ENDDO
!       ENDDO
!       foot_albedos(n) = albedo/sum
!       foot_heights(n) = height/sum
!       foot_areas(n)   = sum
!
!       !  control point for special cases
!
!999    CONTINUE
!
!       !  end loop
!
!    ENDDO
!
!    !  finish
!
!    RETURN
!  END SUBROUTINE tesselations_patch1

END MODULE SAO_OMIL2_ReadLib_Tessel_module
