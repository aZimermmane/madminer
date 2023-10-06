ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP5()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_210 = -(MDL_CQTQD8TI/MDL_LAMBDA__EXP__2)+(MDL_CQTQD8T
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_211 = MDL_CQTQD8TI/MDL_LAMBDA__EXP__2-(MDL_CQTQD8T
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_212 = MDL_CQTQD8TI/MDL_LAMBDA__EXP__2+(MDL_CQTQD8T
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_213 = -(MDL_CQTQD8I/MDL_LAMBDA__EXP__2)+(4.000000D+00
     $ *MDL_CQTQD8TI)/MDL_LAMBDA__EXP__2+(MDL_CQTQD8*MDL_COMPLEXI)
     $ /MDL_LAMBDA__EXP__2-(4.000000D+00*MDL_CQTQD8T*MDL_COMPLEXI)
     $ /MDL_LAMBDA__EXP__2
      GC_215 = -(MDL_CQU1IX3133/MDL_LAMBDA__EXP__2)+(MDL_CQU1X3133
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_216 = MDL_CQU1IX3133/MDL_LAMBDA__EXP__2+(MDL_CQU1X3133
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_219 = -(MDL_CQU1IX3233/MDL_LAMBDA__EXP__2)+(MDL_CQU1X3233
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_220 = MDL_CQU1IX3233/MDL_LAMBDA__EXP__2+(MDL_CQU1X3233
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_231 = -(MDL_CQU8IX3133/MDL_LAMBDA__EXP__2)+(MDL_CQU8X3133
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_232 = MDL_CQU8IX3133/MDL_LAMBDA__EXP__2+(MDL_CQU8X3133
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_235 = -(MDL_CQU8IX3233/MDL_LAMBDA__EXP__2)+(MDL_CQU8X3233
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_236 = MDL_CQU8IX3233/MDL_LAMBDA__EXP__2+(MDL_CQU8X3233
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_247 = -(MDL_CQUQD1IX1333/MDL_LAMBDA__EXP__2)-(MDL_CQUQD1X1333
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_248 = MDL_CQUQD1IX1333/MDL_LAMBDA__EXP__2-(MDL_CQUQD1X1333
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_251 = -(MDL_CQUQD1IX2333/MDL_LAMBDA__EXP__2)-(MDL_CQUQD1X2333
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_252 = MDL_CQUQD1IX2333/MDL_LAMBDA__EXP__2-(MDL_CQUQD1X2333
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_265 = -(MDL_CQUQD1IX3313/MDL_LAMBDA__EXP__2)+(MDL_CQUQD1X3313
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_266 = MDL_CQUQD1IX3313/MDL_LAMBDA__EXP__2+(MDL_CQUQD1X3313
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_269 = -(MDL_CQUQD1IX3323/MDL_LAMBDA__EXP__2)+(MDL_CQUQD1X3323
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_270 = MDL_CQUQD1IX3323/MDL_LAMBDA__EXP__2+(MDL_CQUQD1X3323
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_271 = -(MDL_CQUQD1IX3331/MDL_LAMBDA__EXP__2)-(MDL_CQUQD1X3331
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_272 = MDL_CQUQD1IX3331/MDL_LAMBDA__EXP__2-(MDL_CQUQD1X3331
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_273 = -(MDL_CQUQD1IX3331/MDL_LAMBDA__EXP__2)+(MDL_CQUQD1X3331
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_274 = MDL_CQUQD1IX3331/MDL_LAMBDA__EXP__2+(MDL_CQUQD1X3331
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      GC_275 = -(MDL_CQUQD1IX3332/MDL_LAMBDA__EXP__2)-(MDL_CQUQD1X3332
     $ *MDL_COMPLEXI)/MDL_LAMBDA__EXP__2
      END