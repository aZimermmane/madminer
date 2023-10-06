C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjM(-1,1) -
C      Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjM(-1,1) +
C      Gamma(3,2,-2)*Gamma(4,-2,-1)*ProjP(-1,1) -
C      Gamma(3,-2,-1)*Gamma(4,2,-2)*ProjP(-1,1)
C     
      SUBROUTINE FFVV4_0(F1, F2, V3, V4, COUP,VERTEX)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      COMPLEX*16 COUP
      COMPLEX*16 F1(*)
      COMPLEX*16 F2(*)
      COMPLEX*16 TMP16
      COMPLEX*16 TMP17
      COMPLEX*16 TMP18
      COMPLEX*16 TMP19
      COMPLEX*16 V3(*)
      COMPLEX*16 V4(*)
      COMPLEX*16 VERTEX
      TMP16 = (F1(3)*(F2(3)*(V3(3)*(V4(3)+V4(6))+(V3(4)*(-1D0)*(V4(4)
     $ +CI*(V4(5)))+(V3(5)*(+CI*(V4(4))-V4(5))-V3(6)*(V4(3)+V4(6)))))
     $ +F2(4)*(V3(3)*(V4(4)+CI*(V4(5)))+(V3(4)*(-1D0)*(V4(3)+V4(6))
     $ +(V3(5)*(-1D0)*(+CI*(V4(3)+V4(6)))+V3(6)*(V4(4)+CI*(V4(5)))))))
     $ +F1(4)*(F2(3)*(V3(3)*(V4(4)-CI*(V4(5)))+(V3(4)*(-V4(3)+V4(6))
     $ +(V3(5)*(+CI*(V4(3))-CI*(V4(6)))+V3(6)*(-V4(4)+CI*(V4(5))))))
     $ +F2(4)*(V3(3)*(V4(3)-V4(6))+(V3(4)*(-V4(4)+CI*(V4(5)))+(V3(5)*(
     $ -1D0)*(+CI*(V4(4))+V4(5))+V3(6)*(V4(3)-V4(6)))))))
      TMP17 = (F1(3)*(F2(3)*(V3(3)*(V4(3)-V4(6))+(V3(4)*(-V4(4)+CI
     $ *(V4(5)))+(V3(5)*(-1D0)*(+CI*(V4(4))+V4(5))+V3(6)*(V4(3)-V4(6)))
     $ ))+F2(4)*(V3(3)*(-1D0)*(V4(4)+CI*(V4(5)))+(V3(4)*(V4(3)+V4(6))
     $ +(V3(5)*(+CI*(V4(3)+V4(6)))-V3(6)*(V4(4)+CI*(V4(5)))))))+F1(4)
     $ *(F2(3)*(V3(3)*(-V4(4)+CI*(V4(5)))+(V3(4)*(V4(3)-V4(6))+(V3(5)
     $ *(-CI*(V4(3))+CI*(V4(6)))+V3(6)*(V4(4)-CI*(V4(5))))))+F2(4)
     $ *(V3(3)*(V4(3)+V4(6))+(V3(4)*(-1D0)*(V4(4)+CI*(V4(5)))+(V3(5)*(
     $ +CI*(V4(4))-V4(5))-V3(6)*(V4(3)+V4(6)))))))
      TMP18 = (F1(5)*(F2(5)*(V3(3)*(V4(3)-V4(6))+(V3(4)*(-1D0)*(V4(4)
     $ +CI*(V4(5)))+(V3(5)*(+CI*(V4(4))-V4(5))+V3(6)*(V4(3)-V4(6)))))
     $ +F2(6)*(V3(3)*(-1D0)*(V4(4)+CI*(V4(5)))+(V3(4)*(V4(3)-V4(6))
     $ +(V3(5)*(+CI*(V4(3))-CI*(V4(6)))+V3(6)*(V4(4)+CI*(V4(5)))))))
     $ +F1(6)*(F2(5)*(V3(3)*(-V4(4)+CI*(V4(5)))+(V3(4)*(V4(3)+V4(6))
     $ +(V3(5)*(-1D0)*(+CI*(V4(3)+V4(6)))+V3(6)*(-V4(4)+CI*(V4(5))))))
     $ +F2(6)*(V3(3)*(V4(3)+V4(6))+(V3(4)*(-V4(4)+CI*(V4(5)))+(V3(5)*(
     $ -1D0)*(+CI*(V4(4))+V4(5))-V3(6)*(V4(3)+V4(6)))))))
      TMP19 = (F1(5)*(F2(5)*(V3(3)*(V4(3)+V4(6))+(V3(4)*(-V4(4)+CI
     $ *(V4(5)))+(V3(5)*(-1D0)*(+CI*(V4(4))+V4(5))-V3(6)*(V4(3)+V4(6)))
     $ ))+F2(6)*(V3(3)*(V4(4)+CI*(V4(5)))+(V3(4)*(-V4(3)+V4(6))+(V3(5)
     $ *(-CI*(V4(3))+CI*(V4(6)))-V3(6)*(V4(4)+CI*(V4(5)))))))+F1(6)
     $ *(F2(5)*(V3(3)*(V4(4)-CI*(V4(5)))+(V3(4)*(-1D0)*(V4(3)+V4(6))
     $ +(V3(5)*(+CI*(V4(3)+V4(6)))+V3(6)*(V4(4)-CI*(V4(5))))))+F2(6)
     $ *(V3(3)*(V4(3)-V4(6))+(V3(4)*(-1D0)*(V4(4)+CI*(V4(5)))+(V3(5)*(
     $ +CI*(V4(4))-V4(5))+V3(6)*(V4(3)-V4(6)))))))
      VERTEX = COUP*(-CI*(TMP16+TMP18)+CI*(TMP17+TMP19))
      END

