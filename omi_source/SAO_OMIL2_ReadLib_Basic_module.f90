MODULE SAO_OMIL2_ReadLib_Basic_module

  USE ISO_C_BINDING, ONLY: C_LONG
  IMPLICIT NONE

  ! =====================================================
  ! Define KIND variables for single and double precision
  ! =====================================================

  INTEGER, PARAMETER :: i1 = 1, i2 = 2, i3 = 3, i4 = 4, i8 = 8
  INTEGER, PARAMETER :: r4 = KIND(1.0), r8 = KIND(1.0D0)

  ! ---------------------------------------------------------------------
  ! HDF-EOS5 requires long integers (i.e., KIND=8) for some of the
  ! variables in the Fortran wrappers for the C interface when compiling
  ! with 64bit compilers. 32bit compilers require KIND=4 for those
  ! variables and will fail if they are declared as long integers. 
  !
  ! With help from the HDF/HDF-EOS developers, the ISO_C_BINDING C_LONG
  ! is now being used to define a KIND parameter "i48", which has the
  ! correct value based on the compiler being used. Works like a charm.
  ! ---------------------------------------------------------------------
  INTEGER, PARAMETER :: i48 = C_LONG

  ! ----------------
  ! Some error stati
  ! ----------------
  INTEGER (KIND=i4), PARAMETER :: &
       estat_ok = 0, estat_warning = 1, estat_error = 2, estat_fatal = 3, &
       he5_stat_ok = 0, he5_stat_fail = -1

  ! ------------------------
  ! Maximum character length
  ! ------------------------
  INTEGER (KIND=i4), PARAMETER :: maxchlen = 256

  ! -------------------------------------------------------------------
  ! Values for Pi (rad, deg) and Conversions between Degree and Radians
  ! -------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER :: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER :: pihalf     = 0.5_r8  * pi
  REAL (KIND=r8), PARAMETER :: twopi      = 2.0_r8  * pi
  REAL (KIND=r8), PARAMETER :: pi_deg     = 180.0_r8
  REAL (KIND=r8), PARAMETER :: pihalf_deg =  90.0_r8
  REAL (KIND=r8), PARAMETER :: twopi_deg  = 360.0_r8
  REAL (KIND=r8), PARAMETER :: deg2rad    = pi       / 180.0_r8
  REAL (KIND=r8), PARAMETER :: rad2deg    = 180.0_r8 / pi

  REAL (KIND=r8), PARAMETER :: earth_radius = 6378.0_r8

  ! ------------------------------------------------------------
  ! Precison for DEG <-> RAD conversion - anything less than EPS
  ! is effectively ZERO.
  ! ------------------------------------------------------------
  REAL (KIND=r8), PARAMETER :: eps = 1.0E-10_r8

  ! ---------------------------------------------------
  ! Dimension parameters for ELSUNC auxiliary variables
  ! ---------------------------------------------------
  INTEGER (KIND=i4), PARAMETER :: elsunc_np = 11, elsunc_nw = 6

  ! --------------
  ! I/O file units
  ! --------------
  INTEGER (KIND=i4), PARAMETER :: listunit = 66, ounit = 88

  ! ----------------------
  ! OMI Missing Data value
  ! ----------------------
  REAL (KIND=r8),  PARAMETER :: r8missval = -1.0E+30_r8
  REAL (KIND=r4),  PARAMETER :: r4missval = -1.0E+30_r4

END MODULE SAO_OMIL2_ReadLib_Basic_module

