MODULE SAO_OMIL2_ReadLib_Misc_module

  USE  SAO_OMIL2_ReadLib_Basic_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE fourpt_center ( p1, p2, p3, p4, pc, estat )

    ! ----------------------------------------------------------
    ! Finds the corner point between 4 center points; the corner
    ! is defined as the intersection of the two diagonal lines
    ! between P1-P3 and P2-P4. A counter-clockwise orientation
    ! is assumed for the points, and the points themselves are
    ! 2-dim arrays of [x,y] <==> [lon,lat]
    ! ----------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r4), DIMENSION (2), INTENT (IN) :: p1, p2, p3, p4

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r4), DIMENSION (2), INTENT (OUT) :: pc
    INTEGER (KIND=i4),                INTENT (OUT) :: estat

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r4) :: x1, x2, x3, x4, xc, y1, y2, y3, y4, yc, frac13, frac24

    estat = estat_ok

    ! Some short-hand
    x1 = p1(1) ; y1 = p1(2) ; x2 = p2(1) ; y2 = p2(2)
    x3 = p3(1) ; y3 = p3(2) ; x4 = p4(1) ; y4 = p4(2)

    ! More short-hand
    frac13 = 0.0_r4 ; frac24 = 0.0_r4
    IF ( x1 /= x3 .AND. x2 /= x4 ) THEN
       frac13 = ( y3 - y1 ) / ( x3 - x1 )
       frac24 = ( y4 - y2 ) / ( x4 - x2 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    xc = 0.0_r4
    IF ( frac13-frac24 /= 0.0_r4 ) THEN
       xc = ( y2 - y1 - x2*frac24 + x1*frac13 ) / ( frac13 - frac24 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    yc = y1 + (xc - x1) * frac13

    pc(1:2) = (/ xc, yc /)

    RETURN
  END SUBROUTINE fourpt_center


  SUBROUTINE twoline_intersec ( p1, p2, p3, p4, pc, estat )

    ! ----------------------------------------------------------
    ! Finds the intersection of two lines defined by two points
    ! each. Line 1: P1-P2, Line 2: P3-P4. The points themselves
    ! are 2-dim arrays of [x,y] <==> [lon,lat]
    ! ----------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r4), DIMENSION (2), INTENT (IN) :: p1, p2, p3, p4

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r4), DIMENSION (2), INTENT (OUT) :: pc
    INTEGER (KIND=i4),                INTENT (OUT) :: estat

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r4) :: x1, x2, x3, x4, xc, y1, y2, y3, y4, yc, frac12, frac34

    estat = estat_ok

    ! Some short-hand
    x1 = p1(1) ; y1 = p1(2) ; x2 = p2(1) ; y2 = p2(2)
    x3 = p3(1) ; y3 = p3(2) ; x4 = p4(1) ; y4 = p4(2)

    ! More short-hand
    frac12 = 0.0_r4 ; frac34 = 0.0_r4
    IF ( x1 /= x2 .AND. x3 /= x4 ) THEN
       frac12 = ( y2 - y1 ) / ( x2 - x1 )
       frac34 = ( y4 - y3 ) / ( x4 - x3 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    xc = 0.0_r4
    IF ( frac12-frac34 /= 0.0_r4 ) THEN
       xc = ( y3 - y1 - x3*frac34 + x1*frac12 ) / ( frac12 - frac34 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    yc = y1 + (xc - x1) * frac12

    pc(1:2) = (/ xc, yc /)

    RETURN
  END SUBROUTINE twoline_intersec


  SUBROUTINE twoline_intersec_r8 ( p1, p2, p3, p4, pc, estat )

    ! ----------------------------------------------------------
    ! Finds the intersection of two lines defined by two points
    ! each. Line 1: P1-P2, Line 2: P3-P4. The points themselves
    ! are 2-dim arrays of [x,y] <==> [lon,lat]
    ! ----------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), DIMENSION (2), INTENT (IN) :: p1, p2, p3, p4

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8), DIMENSION (2), INTENT (OUT) :: pc
    INTEGER (KIND=i4),                INTENT (OUT) :: estat

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: x1, x2, x3, x4, xc, y1, y2, y3, y4, yc, frac12, frac34

    estat = estat_ok

    ! Some short-hand
    x1 = p1(1) ; y1 = p1(2) ; x2 = p2(1) ; y2 = p2(2)
    x3 = p3(1) ; y3 = p3(2) ; x4 = p4(1) ; y4 = p4(2)

    ! More short-hand
    frac12 = 0.0_r8 ; frac34 = 0.0_r8
    IF ( x1 /= x2 .AND. x3 /= x4 ) THEN
       frac12 = ( y2 - y1 ) / ( x2 - x1 )
       frac34 = ( y4 - y3 ) / ( x4 - x3 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    xc = 0.0_r8
    IF ( frac12-frac34 /= 0.0_r8 ) THEN
       xc = ( y3 - y1 - x3*frac34 + x1*frac12 ) / ( frac12 - frac34 )
    ELSE
       estat = estat_error
       RETURN
    END IF

    yc = y1 + (xc - x1) * frac12

    pc(1:2) = (/ xc, yc /)

    RETURN
  END SUBROUTINE twoline_intersec_r8


END MODULE SAO_OMIL2_ReadLib_Misc_module
