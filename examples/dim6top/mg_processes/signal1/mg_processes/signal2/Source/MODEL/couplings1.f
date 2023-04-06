ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1 = -(MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_2 = (2.000000D+00*MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_649 = (MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_650 = -(MDL_CW*MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_SW)
      GC_660 = -(MDL_EE*MDL_COMPLEXI*MDL_SW)/(6.000000D+00*MDL_CW)
      GC_590 = (MDL_CTQ8*MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_725 = (MDL_CTG*MDL_COMPLEXI*MDL_GSTRONG*MDL_VEV)
     $ /(MDL_LAMBDA__EXP__2*MDL_SQRT__2)
      GC_953 = -((MDL_COMPLEXI*MDL_YB)/MDL_SQRT__2)
      GC_954 = -((MDL_COMPLEXI*MDL_YT)/MDL_SQRT__2)
      END
