C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(-1,2,-2)*Gamma(-1,4,-3)*ProjM(-2,3)*ProjP(-3,1)
C     
      SUBROUTINE FFFF15P1N_4(F1, F2, F3, COUP,F4)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 F3(*)
      COMPLEX*16 F4(6)
      F4(3)= COUP*(-2D0 )* CI * F3(3)*(F2(5)*F1(5)+F2(6)*F1(6))
      F4(4)= COUP*(-2D0 )* CI * F3(4)*(F2(5)*F1(5)+F2(6)*F1(6))
      F4(5)= COUP*0D0
      F4(6)= COUP*0D0
      END

