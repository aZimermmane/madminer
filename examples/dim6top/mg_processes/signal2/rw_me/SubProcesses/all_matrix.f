
C     PY ((21, 21), (-12, -11, -5, 5, 11, 12)) : (21, 21, 5, -11, 12,
C      -5, 11, -12) # M0_ 1
C     PY ((21, 21), (-13, -12, -5, 5, 11, 14)) : (21, 21, 5, -13, 14,
C      -5, 11, -12) # M1_ 1
C     PY ((21, 21), (-14, -11, -5, 5, 12, 13)) : (21, 21, 5, -11, 12,
C      -5, 13, -14) # M2_ 1
C     PY ((21, 21), (-14, -13, -5, 5, 13, 14)) : (21, 21, 5, -13, 14,
C      -5, 13, -14) # M3_ 1
C     PY ((2, -2), (-12, -11, -5, 5, 11, 12)) : (2, -2, 5, -11, 12,
C      -5, 11, -12) # M4_ 1
C     PY ((2, -2), (-13, -12, -5, 5, 11, 14)) : (2, -2, 5, -13, 14,
C      -5, 11, -12) # M5_ 1
C     PY ((2, -2), (-14, -11, -5, 5, 12, 13)) : (2, -2, 5, -11, 12,
C      -5, 13, -14) # M6_ 1
C     PY ((2, -2), (-14, -13, -5, 5, 13, 14)) : (2, -2, 5, -13, 14,
C      -5, 13, -14) # M7_ 1
C     PY ((4, -4), (-12, -11, -5, 5, 11, 12)) : (4, -4, 5, -11, 12,
C      -5, 11, -12) # M8_ 1
C     PY ((4, -4), (-13, -12, -5, 5, 11, 14)) : (4, -4, 5, -13, 14,
C      -5, 11, -12) # M9_ 1
C     PY ((4, -4), (-14, -11, -5, 5, 12, 13)) : (4, -4, 5, -11, 12,
C      -5, 13, -14) # M10_ 1
C     PY ((4, -4), (-14, -13, -5, 5, 13, 14)) : (4, -4, 5, -13, 14,
C      -5, 13, -14) # M11_ 1
C     PY ((1, -1), (-12, -11, -5, 5, 11, 12)) : (1, -1, 5, -11, 12,
C      -5, 11, -12) # M12_ 1
C     PY ((3, -3), (-12, -11, -5, 5, 11, 12)) : (3, -3, 5, -11, 12,
C      -5, 11, -12) # M12_ 1
C     PY ((1, -1), (-13, -12, -5, 5, 11, 14)) : (1, -1, 5, -13, 14,
C      -5, 11, -12) # M13_ 1
C     PY ((3, -3), (-13, -12, -5, 5, 11, 14)) : (3, -3, 5, -13, 14,
C      -5, 11, -12) # M13_ 1
C     PY ((1, -1), (-14, -11, -5, 5, 12, 13)) : (1, -1, 5, -11, 12,
C      -5, 13, -14) # M14_ 1
C     PY ((3, -3), (-14, -11, -5, 5, 12, 13)) : (3, -3, 5, -11, 12,
C      -5, 13, -14) # M14_ 1
C     PY ((1, -1), (-14, -13, -5, 5, 13, 14)) : (1, -1, 5, -13, 14,
C      -5, 13, -14) # M15_ 1
C     PY ((3, -3), (-14, -13, -5, 5, 13, 14)) : (3, -3, 5, -13, 14,
C      -5, 13, -14) # M15_ 1
C     PY ((5, -5), (-12, -11, -5, 5, 11, 12)) : (5, -5, 5, -11, 12,
C      -5, 11, -12) # M28_ 1
C     PY ((5, -5), (-13, -12, -5, 5, 11, 14)) : (5, -5, 5, -13, 14,
C      -5, 11, -12) # M29_ 1
C     PY ((5, -5), (-14, -11, -5, 5, 12, 13)) : (5, -5, 5, -11, 12,
C      -5, 13, -14) # M30_ 1
C     PY ((5, -5), (-14, -13, -5, 5, 13, 14)) : (5, -5, 5, -13, 14,
C      -5, 13, -14) # M31_ 1
      SUBROUTINE SMATRIXHEL(PDGS, PROCID, NPDG, P, ALPHAS, SCALE2,
     $  NHEL, ANS)
      IMPLICIT NONE
C     ALPHAS is given at scale2 (SHOULD be different of 0 for loop
C      induced, ignore for LO)  

CF2PY double precision, intent(in), dimension(0:3,npdg) :: p
CF2PY integer, intent(in), dimension(npdg) :: pdgs
CF2PY integer, intent(in):: procid
CF2PY integer, intent(in) :: npdg
CF2PY double precision, intent(out) :: ANS
CF2PY double precision, intent(in) :: ALPHAS
CF2PY double precision, intent(in) :: SCALE2
      INTEGER PDGS(*)
      INTEGER NPDG, NHEL, PROCID
      DOUBLE PRECISION P(*)
      DOUBLE PRECISION ANS, ALPHAS, PI,SCALE2
      INCLUDE 'coupl.inc'

      PI = 3.141592653589793D0
      G = 2* DSQRT(ALPHAS*PI)
      CALL UPDATE_AS_PARAM()
C     if (scale2.ne.0d0) stop 1

      IF(21.EQ.PDGS(1).AND.21.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 0
        CALL M0_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(21.EQ.PDGS(1).AND.21.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 1
        CALL M1_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(21.EQ.PDGS(1).AND.21.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 2
        CALL M2_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(21.EQ.PDGS(1).AND.21.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 3
        CALL M3_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 4
        CALL M4_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 5
        CALL M5_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 6
        CALL M6_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(2.EQ.PDGS(1).AND.-2.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 7
        CALL M7_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 8
        CALL M8_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 9
        CALL M9_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 10
        CALL M10_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(4.EQ.PDGS(1).AND.-4.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 11
        CALL M11_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 12
        CALL M12_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 13
        CALL M12_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 14
        CALL M13_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 15
        CALL M13_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 16
        CALL M14_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 17
        CALL M14_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(1.EQ.PDGS(1).AND.-1.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 18
        CALL M15_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(3.EQ.PDGS(1).AND.-3.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 19
        CALL M15_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 20
        CALL M28_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.11.EQ.PDGS(7).AND.-12.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 21
        CALL M29_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -11.EQ.PDGS(4).AND.12.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 22
        CALL M30_SMATRIXHEL(P, NHEL, ANS)
      ELSE IF(5.EQ.PDGS(1).AND.-5.EQ.PDGS(2).AND.5.EQ.PDGS(3).AND.
     $ -13.EQ.PDGS(4).AND.14.EQ.PDGS(5).AND.-5.EQ.PDGS(6)
     $ .AND.13.EQ.PDGS(7).AND.-14.EQ.PDGS(8)
     $ .AND.(PROCID.LE.0.OR.PROCID.EQ.1)) THEN  ! 23
        CALL M31_SMATRIXHEL(P, NHEL, ANS)
      ENDIF

      RETURN
      END

      SUBROUTINE INITIALISE(PATH)
C     ROUTINE FOR F2PY to read the benchmark point.
      IMPLICIT NONE
      CHARACTER*512 PATH
CF2PY INTENT(IN) :: PATH
      CALL SETPARA(PATH)  !first call to setup the paramaters
      RETURN
      END


      SUBROUTINE CHANGE_PARA(NAME, VALUE)
      IMPLICIT NONE
CF2PY intent(in) :: name
CF2PY intent(in) :: value

      CHARACTER*512 NAME
      DOUBLE PRECISION VALUE

      LOGICAL M14_HELRESET
      COMMON /M14_HELRESET/ M14_HELRESET
      LOGICAL M7_HELRESET
      COMMON /M7_HELRESET/ M7_HELRESET
      LOGICAL M13_HELRESET
      COMMON /M13_HELRESET/ M13_HELRESET
      LOGICAL M12_HELRESET
      COMMON /M12_HELRESET/ M12_HELRESET
      LOGICAL M10_HELRESET
      COMMON /M10_HELRESET/ M10_HELRESET
      LOGICAL M8_HELRESET
      COMMON /M8_HELRESET/ M8_HELRESET
      LOGICAL M1_HELRESET
      COMMON /M1_HELRESET/ M1_HELRESET
      LOGICAL M6_HELRESET
      COMMON /M6_HELRESET/ M6_HELRESET
      LOGICAL M29_HELRESET
      COMMON /M29_HELRESET/ M29_HELRESET
      LOGICAL M4_HELRESET
      COMMON /M4_HELRESET/ M4_HELRESET
      LOGICAL M28_HELRESET
      COMMON /M28_HELRESET/ M28_HELRESET
      LOGICAL M2_HELRESET
      COMMON /M2_HELRESET/ M2_HELRESET
      LOGICAL M11_HELRESET
      COMMON /M11_HELRESET/ M11_HELRESET
      LOGICAL M15_HELRESET
      COMMON /M15_HELRESET/ M15_HELRESET
      LOGICAL M5_HELRESET
      COMMON /M5_HELRESET/ M5_HELRESET
      LOGICAL M31_HELRESET
      COMMON /M31_HELRESET/ M31_HELRESET
      LOGICAL M30_HELRESET
      COMMON /M30_HELRESET/ M30_HELRESET
      LOGICAL M9_HELRESET
      COMMON /M9_HELRESET/ M9_HELRESET
      LOGICAL M0_HELRESET
      COMMON /M0_HELRESET/ M0_HELRESET
      LOGICAL M3_HELRESET
      COMMON /M3_HELRESET/ M3_HELRESET

      INCLUDE '../Source/MODEL/input.inc'
      INCLUDE '../Source/MODEL/coupl.inc'

      M14_HELRESET = .TRUE.
      M7_HELRESET = .TRUE.
      M13_HELRESET = .TRUE.
      M12_HELRESET = .TRUE.
      M10_HELRESET = .TRUE.
      M8_HELRESET = .TRUE.
      M1_HELRESET = .TRUE.
      M6_HELRESET = .TRUE.
      M29_HELRESET = .TRUE.
      M4_HELRESET = .TRUE.
      M28_HELRESET = .TRUE.
      M2_HELRESET = .TRUE.
      M11_HELRESET = .TRUE.
      M15_HELRESET = .TRUE.
      M5_HELRESET = .TRUE.
      M31_HELRESET = .TRUE.
      M30_HELRESET = .TRUE.
      M9_HELRESET = .TRUE.
      M0_HELRESET = .TRUE.
      M3_HELRESET = .TRUE.

      SELECT CASE (NAME)
      CASE ('Lambda')
      MDL_LAMBDA = VALUE
      CASE ('DIM6_1')
      MDL_LAMBDA = VALUE
      CASE ('ctG')
      MDL_CTG = VALUE
      CASE ('DIM6_16')
      MDL_CTG = VALUE
      CASE ('ctq8')
      MDL_CTQ8 = VALUE
      CASE ('DIM6_55')
      MDL_CTQ8 = VALUE
      CASE ('MB')
      MDL_MB = VALUE
      CASE ('MASS_5')
      MDL_MB = VALUE
      CASE ('MT')
      MDL_MT = VALUE
      CASE ('MASS_6')
      MDL_MT = VALUE
      CASE ('MTA')
      MDL_MTA = VALUE
      CASE ('MASS_15')
      MDL_MTA = VALUE
      CASE ('MZ')
      MDL_MZ = VALUE
      CASE ('MASS_23')
      MDL_MZ = VALUE
      CASE ('MH')
      MDL_MH = VALUE
      CASE ('MASS_25')
      MDL_MH = VALUE
      CASE ('aEWM1')
      AEWM1 = VALUE
      CASE ('SMINPUTS_1')
      AEWM1 = VALUE
      CASE ('Gf')
      MDL_GF = VALUE
      CASE ('SMINPUTS_2')
      MDL_GF = VALUE
      CASE ('aS')
      AS = VALUE
      CASE ('SMINPUTS_3')
      AS = VALUE
      CASE ('ymb')
      MDL_YMB = VALUE
      CASE ('YUKAWA_5')
      MDL_YMB = VALUE
      CASE ('ymt')
      MDL_YMT = VALUE
      CASE ('YUKAWA_6')
      MDL_YMT = VALUE
      CASE ('ymtau')
      MDL_YMTAU = VALUE
      CASE ('YUKAWA_15')
      MDL_YMTAU = VALUE
      CASE ('WT')
      MDL_WT = VALUE
      CASE ('DECAY_6')
      MDL_WT = VALUE
      CASE ('WZ')
      MDL_WZ = VALUE
      CASE ('DECAY_23')
      MDL_WZ = VALUE
      CASE ('WW')
      MDL_WW = VALUE
      CASE ('DECAY_24')
      MDL_WW = VALUE
      CASE ('WH')
      MDL_WH = VALUE
      CASE ('DECAY_25')
      MDL_WH = VALUE
      CASE DEFAULT
      WRITE(*,*) 'no parameter matching', NAME, VALUE
      END SELECT

      RETURN
      END

      SUBROUTINE UPDATE_ALL_COUP()
      IMPLICIT NONE
      CALL COUP()
      RETURN
      END


      SUBROUTINE GET_PDG_ORDER(PDG, ALLPROC)
      IMPLICIT NONE
CF2PY INTEGER, intent(out) :: PDG(24,8)
CF2PY INTEGER, intent(out) :: ALLPROC(24)
      INTEGER PDG(24,8), PDGS(24,8)
      INTEGER ALLPROC(24),PIDS(24)
      DATA PDGS/ 21,21,21,21,2,2,2,2,4,4,4,4,1,3,1,3,1,3,1,3,5,5,5,5
     $ ,21,21,21,21,-2,-2,-2,-2,-4,-4,-4,-4,-1,-3,-1,-3,-1,-3,-1,-3,-5
     $ ,-5,-5,-5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,-11,
     $ -13,-11,-13,-11,-13,-11,-13,-11,-13,-11,-13,-11,-11,-13,-13,-11
     $ ,-11,-13,-13,-11,-13,-11,-13,12,14,12,14,12,14,12,14,12,14,12
     $ ,14,12,12,14,14,12,12,14,14,12,14,12,14,-5,-5,-5,-5,-5,-5,-5,-5
     $ ,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,11,11,13,13,11
     $ ,11,13,13,11,11,13,13,11,11,11,11,13,13,13,13,11,11,13,13,-12,
     $ -12,-14,-14,-12,-12,-14,-14,-12,-12,-14,-14,-12,-12,-12,-12,-14
     $ ,-14,-14,-14,-12,-12,-14,-14 /
      DATA PIDS/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /
      PDG = PDGS
      ALLPROC = PIDS
      RETURN
      END

      SUBROUTINE GET_PREFIX(PREFIX)
      IMPLICIT NONE
CF2PY CHARACTER*20, intent(out) :: PREFIX(24)
      CHARACTER*20 PREFIX(24),PREF(24)
      DATA PREF / 'M0_','M1_','M2_','M3_','M4_','M5_','M6_','M7_'
     $ ,'M8_','M9_','M10_','M11_','M12_','M12_','M13_','M13_','M14_'
     $ ,'M14_','M15_','M15_','M28_','M29_','M30_','M31_'/
      PREFIX = PREF
      RETURN
      END



