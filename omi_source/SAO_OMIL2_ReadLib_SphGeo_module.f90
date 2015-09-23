MODULE SAO_OMIL2_ReadLib_SphGeo_module

  USE  SAO_OMIL2_ReadLib_Basic_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE sphergeom_asa ( c, alp, bet, a, b )

    ! -------------------------------------------------------------------
    ! ASA: Angle-Side-Angle 
    ! 
    ! Computation of quantities in a spherical triangle: Given two angles
    ! ALP and BET, and the side C inbetween them, computes the sides A
    ! and B.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: c, alp, bet

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: a, b

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2, gam

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    a = 0.0_r8 ; b = 0.0_r8

    ! -------------------------------------
    ! Compute the angle BET opposite side B
    ! -------------------------------------
    tmp1 = ATAN(TAN(c/2.0_r8) * COS((alp-bet)/2.0_r8) / COS((alp+bet)/2.0_r8))
    tmp2 = ATAN(TAN(c/2.0_r8) * SIN((alp-bet)/2.0_r8) / SIN((alp+bet)/2.0_r8))
    
    a = tmp1 + tmp2
    b = tmp1 - tmp2

    gam = ACOS(-COS(alp)*COS(bet) + SIN(alp)*SIN(bet)*COS(c))
    a   = ASIN(SIN(c)*SIN(alp)/SIN(gam))
    b   = ASIN(SIN(c)*SIN(bet)/SIN(gam))

    RETURN
  END SUBROUTINE sphergeom_asa

  SUBROUTINE sphergeom_ssa ( a, b, alp, c, bet, gam )

    ! -------------------------------------------------------------------
    ! SSA: Side-Side-Angle 
    ! 
    ! Computation of quantities in a spherical triangle: Given two sides 
    ! A, B and an angle ALP opposite to one of the sides, computes the
    ! third side and the remaining two angles in the triangle.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a, b, alp

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: c, bet, gam

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c = 0.0_r8 ; bet = 0.0_r8 ; gam = 0.0_r8

    ! -------------------------------------
    ! Compute the angle BET opposite side B
    ! -------------------------------------
    tmp1 = SIN(b) * SIN(alp) ; tmp2 = SIN(a)
    IF ( ABS(ABS(tmp1)-ABS(tmp2)) < eps ) THEN
       bet = pihalf
    ELSE IF ( ABS(tmp1) < eps ) THEN
       bet = 0.0_r8
    ELSE
       bet = ASIN(tmp1/tmp2)
    END IF

    ! -----------------------------------
    ! Now the length of C, the third side
    ! -----------------------------------
    tmp1 = TAN(b) * COS(alp) ; tmp2 = TAN(a) * COS(bet)
    c = ATAN(tmp1) + ATAN(tmp2)

    ! -------------------------------------------
    ! Finally the angle GAM between sides A and B
    ! -------------------------------------------
    ! Note: COT(x) = -TAN(x+pi/2)
    ! -----------------------------------------------------------
    ! COT(phi) = COS(b) * TAN(alp) ; COT(psi) = COS(a) * TAN(bet)
    ! gam = phi + psi
    ! -----------------------------------------------------------
    tmp1 = COS(b) * TAN(alp) ; tmp2 = COS(a) * TAN(bet)
    gam = ATAN(-tmp1) + ATAN(-tmp2) - pi
    IF ( ABS(ABS(gam)-pi) < eps ) gam = 0.0_r8

    RETURN
  END SUBROUTINE sphergeom_ssa

  SUBROUTINE sphergeom_sas ( a, b, gam, c, alp, bet )

    ! -------------------------------------------------------------------
    ! SAS: Side-Angle-Side 
    ! 
    ! Computation of quantities in a spherical triangle: Given two sides 
    ! A, B and the angle GAM inbetween them, computes the third side and
    ! one of the remaining two angles in the triangle.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a, b, gam

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: c, alp, bet

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c = 0.0_r8 ; alp = 0.0_r8

    ! -----------------------------------------------------------
    ! Compute length of baseline c between the two points A and B
    ! -----------------------------------------------------------
    tmp1 = COS(a) * COS(b) + SIN(a) * SIN(b) * COS(gam)
    IF ( ABS(tmp1) < eps ) THEN
       c = pihalf
    ELSE
       c = ACOS(tmp1)
    END IF

    ! ----------------------------------------------------
    ! Now the angle ALP at point A (between sides a and c)
    ! ----------------------------------------------------
    tmp1 = COS(a) - COS(b) * COS(c) ; tmp2 = SIN(b) * SIN(c)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          alp = ACOS(tmp1/tmp2)
       ELSE
          alp = 0.0_r8
       END IF
    ELSE
       alp = 0.0_r8
    END IF

    ! ----------------------------------------------------
    ! Now the angle BET at point B (between sides b and c)
    ! ----------------------------------------------------
    tmp1 = COS(b) - COS(c) * COS(a) ; tmp2 = SIN(c) * SIN(a)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          bet = ACOS(tmp1/tmp2)
       ELSE
          bet = 0.0_r8
       END IF
    ELSE
       bet = 0.0_r8
    END IF

    RETURN
  END SUBROUTINE sphergeom_sas

  SUBROUTINE sphergeom_coordinates_at_point ( a0, b0, gam0, c_in, abs_or_frac, a, gam )

    ! -----------------------------------------------------------
    ! Finds the co-ordinates of the baseline extended from two
    ! lon/lat points on a sphere given the hypotenuse C_IN.
    ! -----------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a0, b0, gam0, c_in
    CHARACTER (LEN=3), INTENT (IN) :: abs_or_frac

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: a, gam

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2, c0, c, alp0, gamsig

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    a = 0.0_r8 ; gam = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp1 = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp1) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp1)
    END IF

    tmp1 = COS(a0) - COS(b0) * COS(c0) ; tmp2 = SIN(b0) * SIN(c0)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          alp0 = ACOS(tmp1/tmp2)
       ELSE
          alp0 = 0.0_r8
       END IF
    ELSE
       alp0 = 0.0_r8
    END IF

    ! -------------------------------------------------------------
    ! Now that we have the basic parameters of the triangle, we can
    ! compute the lat/lon at the desired point.
    ! -------------------------------------------------------------
    SELECT CASE ( abs_or_frac )
    CASE ( 'abs' )
       c = c_in
    CASE ( 'frc' )
       c = c0 * c_in
    CASE DEFAULT
       ! ------------------------------------------------------------
       ! Temporary fix for unknown option: Complain and set c to c0/2
       ! ------------------------------------------------------------
       WRITE (*,'(A,A)') 'ERROR: Unknown option -- ', abs_or_frac
       c = c0 / 2.0_r8
    END SELECT

    ! -----------------------------
    ! The "latitude" of the mid-point
    ! -----------------------------
    tmp1 = COS(b0) * COS(c) + SIN(b0) * SIN(c) * COS(alp0)
    IF ( ABS(tmp1) < eps ) THEN
       a = pihalf
    ELSE
       a = ACOS(tmp1)
    END IF

    ! --------------------------------
    ! The "longitude" of the mid-point
    ! --------------------------------
    tmp1 = COS(c) - COS(a) * COS(b0)
    tmp2 = SIN(a) * SIN(b0)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( gam0 /= 0.0_r8 .AND. ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          gam = ACOS(tmp1/tmp2)
       ELSE
          gam = 0.0_r8
       END IF
    ELSE
       gam = 0.0_r8
    END IF

    IF ( ((gam0 > 0.0_r8) .AND. (gam0 < pi)) .OR. (gam0 < -pi) ) THEN
       gamsig = 1.0_r8
    ELSE
       gamsig = -1.0_r8
    END IF

    ! --------------------------------------
    ! Return the angle with the correct sign
    ! --------------------------------------
    gam = ABS(gam) * gamsig

    RETURN
  END SUBROUTINE sphergeom_coordinates_at_point

  SUBROUTINE sphergeom_baseline_comp ( a0, b0, gam0, c0 )

    ! -------------------------------------------------------
    ! Finds the lengh of the baseline of a spherical triangle
    ! -------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: a0, b0, gam0

    ! ---------------
    ! Output variable
    ! ---------------
    REAL (KIND=r8), INTENT (OUT) :: c0

    ! --------------
    ! Local variable
    ! --------------
    REAL (KIND=r8) :: tmp

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c0 = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp)
    END IF

    RETURN
  END SUBROUTINE sphergeom_baseline_comp

  SUBROUTINE lonlat_to_pi ( lon, lat )

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), INTENT (INOUT) :: lon, lat

    ! ------------------------------------
    ! Adjust longitude values to [-pi,+pi]
    ! ------------------------------------
    IF ( ABS(lon) > twopi ) lon = MOD(lon, twopi)
    IF ( lon >  pi ) lon = lon - twopi
    IF ( lon < -pi ) lon = lon + twopi
    ! ---------------------------------------
    ! Adjust latitude values to [-pi/2,+pi/2]
    ! ---------------------------------------
    IF ( ABS(lat) > pihalf ) lat = MOD(lat, pihalf)
    IF ( lat >  pihalf ) lat =   pi - lat
    IF ( lat < -pihalf ) lat = -(pi + lat)
    
    RETURN
  END SUBROUTINE lonlat_to_pi

  REAL (KIND=r8) FUNCTION angle_minus_twopi ( gamma0, pival ) RESULT ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: gamma0, pival

    IF ( ABS(gamma0) > pival ) THEN
       gamma = gamma0 - 2.0_r8 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    ELSE
       gamma = gamma0
    END IF

    RETURN
  END FUNCTION angle_minus_twopi


  SUBROUTINE sphgeo_comp_pixel_corners ( &
       nxtrack, ntimes, lon, lat, clon, clat, estat, yn_omi_pixel_adjust_k )

    ! ------------------------------------------------------------------------
    ! Compute corner coordinates of ground pixels given only the pixel centers
    ! ------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                    INTENT (IN) :: nxtrack, ntimes
    REAL    (KIND=r4), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT (IN) :: lon, lat
    LOGICAL, OPTIONAL,                                    INTENT (IN) :: yn_omi_pixel_adjust_k

    ! ----------------
    ! Output variables
    ! ---------------------------------------------------
    ! NOTE: The corner coordinates are R8 rather that R4!
    ! ---------------------------------------------------
    INTEGER (KIND=i4),                                  INTENT (OUT) :: estat
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes), INTENT (OUT) :: clon, clat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: fstat, i, j
    REAL    (KIND=r8)                                   :: a0, b0, c0, gam0, alp0, a, gam
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: lonrad, latrad
    REAL    (KIND=r8), DIMENSION (0:nxtrack,0:ntimes)   :: tmplat, tmplon
    LOGICAL                                             :: yn_pixel_adjust_crosstrack

    estat = estat_ok

    ! --------------------------------------------------------------------------------
    ! Check whether we want to perform some pixel adjustment. This is an issue for
    ! OMI, where across-track pixels are of different size (larger towards the edges of
    ! the swath). Using the "center" method between adjacent across-track pixels will
    ! introduce a small error in the pixel boundaries, and the pixels at the sides of
    ! the swath will be underestimated. 
    !
    ! Note that this adjustment is still experimental until the distortions around the
    ! poles have been solved. Hence the default is not to make this adjustment.
    ! --------------------------------------------------------------------------------
    yn_pixel_adjust_crosstrack = .FALSE.
    IF ( PRESENT (yn_omi_pixel_adjust_k) ) yn_pixel_adjust_crosstrack = yn_omi_pixel_adjust_k

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------
    lonrad = REAL ( lon, KIND=r8 ) * deg2rad
    latrad = REAL ( lat, KIND=r8 ) * deg2rad

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    clon   = -999.9_r4 ; clat   = -999.9_r4
    tmplon = -999.9_r8 ; tmplat = -999.9_r8

    a0   = 0.0_r8 ; b0 = 0.0_r8 ; c0  = 0.0_r8
    gam0 = 0.0_r8 ; a  = 0.0_r8 ; gam = 0.0_r8

    ! ---------------------------------------------------
    ! First interpolate between center points ALONG-TRACK
    ! ---------------------------------------------------
    ! Array points filled: [1:nxtrack, 0:ntimes]
    ! ---------------------------------------------------
    DO i = 1, nxtrack
       DO j = 0, ntimes -2
          a0   = pihalf - latrad(i,j  )
          b0   = pihalf - latrad(i,j+1) 
          gam0 = angle_minus_twopi ( lonrad(i,j) - lonrad(i,j+1), pi )

          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 0.5_r8, 'frc', a, gam )
          tmplat(i,j+1) = pihalf - a
          tmplon(i,j+1) = lonrad(i,j) - gam
          CALL lonlat_to_pi (tmplon(i,j+1), tmplat(i,j+1) )
       END DO

       j = 0
       a0 = pihalf - tmplat(i,j+1)
       b0 = pihalf - tmplat(i,j+2)
       gam0 = angle_minus_twopi ( tmplon(i,j+1) - tmplon(i,j+2), pi )
       CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
       tmplat(i,j) = pihalf - a
       tmplon(i,j) = tmplon(i,j+2) + gam
       CALL lonlat_to_pi (tmplon(i,j), tmplat(i,j) )

       j = ntimes
       a0 = pihalf - tmplat(i,j-1)
       b0 = pihalf - tmplat(i,j-2)
       gam0 = angle_minus_twopi ( tmplon(i,j-1) - tmplon(i,j-2), pi )
       CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
       tmplat(i,j) = pihalf - a
       tmplon(i,j) = tmplon(i,j-2) + gam
       CALL lonlat_to_pi (tmplon(i,j), tmplat(i,j) )
    END DO

    !OPEN (UNIT=66, FILE='temporary_coordinates.dat', STATUS='UNKNOWN', ACTION='WRITE')
    !WRITE (66,'(5I10)') ntimes+1, nxtrack+1
    !DO j = 0, ntimes
    !   DO i = 0, nxtrack
    !      WRITE (66,'(2I5, 8F18.10)') j, i, tmplon(i,j), tmplat(i,j)
    !   END DO
    !END DO
    !CLOSE (66)

    ! ----------------------------------------------------------
    ! Now interpolate ACROSS-TRACK whatever points are available
    ! ----------------------------------------------------------  
    DO j = 0, ntimes
       ! ----------------------------------
       ! First the center line of the swath
       ! ----------------------------------------------------------------------
       ! To be save, we compute ALL possible pixels using the center cross
       ! method. Only in favorable circumstances do we use a more sophisticated
       ! approach, and having all pixels precomputed saves us some IF THEN ELSE
       ! logic and several lines of code.
       ! ----------------------------------------------------------------------
       DO i = 1, nxtrack-1
          a0   = pihalf - tmplat(i+1,j)
          b0   = pihalf - tmplat(i,  j) 
          gam0 = angle_minus_twopi ( tmplon(i+1,j) - tmplon(i,j), pi )
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 0.5_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = tmplon(i+1,j) - gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )

          IF ( i == 1 ) THEN
             a0 = pihalf - tmplat(i,j)
             b0 = pihalf -   clat(i,j)
             gam0 = angle_minus_twopi ( tmplon(i,j) - clon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
             clat(i-1,j) = pihalf - a
             clon(i-1,j) = clon(i,j) + gam
             CALL lonlat_to_pi (clon(i-1,j), clat(i-1,j) )
          END IF
          IF ( i == nxtrack-1 ) THEN
             a0 = pihalf - tmplat(i+1,j)
             b0 = pihalf -   clat(i,j)
             gam0 = angle_minus_twopi ( tmplon(i+1,j) - clon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
             clat(i+1,j) = pihalf - a
             clon(i+1,j) = clon(i,j) + gam
             CALL lonlat_to_pi (clon(i+1,j), clat(i+1,j) )
          END IF
       END DO
    END DO

    !!! ===============================================================================!!!
    !!! The following code is experimental and doesn't quite work yet. It is supposed  !!!
    !!! to adjust the corner coordinates according to the true distance between the    !!!
    !!! pixel centers, and it does so by working outwards from the center of the swath !!!
    !!! to the across-track edges. Known issues are:                                   !!!
    !!!   * At the along-track edges, we seem to be mis-aligning the pixels and are    !!!
    !!!     losing the proper last line of corners                                     !!!
    !!!   * At the poles, the distortion of pixels progresses outwards, resulting in   !!!
    !!!     heavy distortions past the poles (going outward). This can be ameliorated  !!!
    !!!     somewhat by excluding pole-most latitudes (within 3 deg of pole), but this !!!
    !!!     still leaves some undesirable distortions.                                 !!!
    !!! ===============================================================================!!!
    IF ( yn_pixel_adjust_crosstrack ) THEN
       DO j = 0, ntimes
          ! -------------------------------------------------------------------
          ! From the center of the swath to lower cross-track pixel numbers
          ! -------------------------------------------------------------------
          DO i = nxtrack/2-1, 1, -1
             a0 = pihalf -   clat(i+1,j)
             b0 = pihalf - tmplat(i+1,j) 
             gam0 = clon(i+1,j) - tmplon(i+1,j)
             gam0 = angle_minus_twopi ( gam0, pi )
             CALL sphergeom_baseline_comp ( a0, b0, gam0, c0 )
             a0 = pihalf - tmplat(i+1,j)
             b0 = pihalf - tmplat(i,  j) 
             gam0 = angle_minus_twopi ( tmplon(i+1,j) - tmplon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, c0, 'abs', a, gam )
             clat(i,j) = pihalf - a
             clon(i,j) = tmplon(i+1,j) - gam
             CALL lonlat_to_pi (clon(i,j), clat(i,j) )
          END DO
          i = 0
          a0 = pihalf - tmplat(i+1,j)
          b0 = pihalf -   clat(i+1,j)
          gam0 = tmplon(i+1,j) - clon(i+1,j)
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = clon(i+1,j) + gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )
       
       ! -----------------------------------------------------------------------
       ! From the center of the swath to higer cross-track pixel numbers
       ! -----------------------------------------------------------------------
          DO i = nxtrack/2+1, nxtrack-1
             a0 = pihalf - tmplat(i,  j)
             b0 = pihalf -   clat(i-1,j)
             gam0 = tmplon(i,j) - clon(i-1,j)
             gam0 = angle_minus_twopi ( gam0, pi )
             CALL sphergeom_baseline_comp ( a0, b0, gam0, c0 )
             IF ( j == ntimes .AND. i == nxtrack ) WRITE (*,'(4F10.6)') a0,b0,gam0,c0
             a0 = pihalf - tmplat(i,  j)
             b0 = pihalf - tmplat(i+1,j) 
             gam0 = angle_minus_twopi ( tmplon(i,j) - tmplon(i+1,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, c0, 'abs', a, gam )
             IF ( j == ntimes .AND. i == nxtrack ) WRITE (*,'(6F10.6)') a0,b0,gam0,c0, a, gam
             clat(i,j) = pihalf - a
             clon(i,j) = tmplon(i,j) - gam
             CALL lonlat_to_pi (clon(i,j), clat(i,j) )
          END DO
          i = nxtrack
          a0 = pihalf - tmplat(i,  j)
          b0 = pihalf - clat  (i-1,j)
          gam0 = tmplon(i,j) - clon(i-1,j)
          gam0 = angle_minus_twopi ( gam0, pi )
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = clon(i-1,j) + gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )
       END DO
    END IF

    ! ------------------------------------
    ! Convert lons and lats back to Degree
    ! ------------------------------------
    clon = clon * rad2deg
    clat = clat * rad2deg

    RETURN
  END SUBROUTINE sphgeo_comp_pixel_corners

END MODULE SAO_OMIL2_ReadLib_SphGeo_module
