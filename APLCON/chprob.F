

      DOUBLE PRECISION FUNCTION CHPROB(CHISQ,N) ! prob from chisquare
*     chi square probability for N degrees of freedom at CHISQ
*     =integral from 0 ... CHISQ   
      IMPLICIT NONE
      DOUBLE PRECISION CHISQ,DGAMIN
      INTEGER N
*     ...
      IF(CHISQ.LE.0.0) THEN
         CHPROB=1.0D0
      ELSE
         CHPROB=1.0D0-DGAMIN(0.5D0*N,0.5D0*CHISQ)
      END IF 
      END 

      DOUBLE PRECISION FUNCTION DGAMIN(A,X)
*     incomplete gamma function P(a,x)
*     returns -1.0 for x < 0 or a =< 0
      IMPLICIT NONE
      DOUBLE PRECISION A,X,EPS,FPMIN
      INTEGER ITMAX,I,N
      PARAMETER (ITMAX=300,EPS=3.0E-7,FPMIN=1.0E-30)
      DOUBLE PRECISION AN,AP,B,C,D,DEL,H,SUM,DGAMML,GLN
*     ...
      IF(X.LT.0.0D0.OR.A.LE.0.0D0) THEN
         DGAMIN=-1.0D0  ! error
      ELSE IF(X.EQ.0.0) THEN
         DGAMIN=0.0D0
      ELSE IF(X.LT.A+1.0D0) THEN ! series representation
         GLN=DGAMML(A)         ! ln[Gamma(a)]
         AP=A  
         SUM=1.0D0/A
         DEL=SUM
         DO N=1,ITMAX
          AP=AP+1.
          DEL=DEL*X/AP
          SUM=SUM+DEL
          IF(ABS(DEL).LT.ABS(SUM)*EPS) GOTO 10
         END DO
 10      DGAMIN=SUM*EXP(-X+A*LOG(X)-GLN)
      ELSE                     ! continued fraction representation
         GLN=DGAMML(A)         ! ln[Gamma(a)]
         B=X+1.0D0-A
         C=1.0D0/FPMIN
         D=1.0D0/B
         H=D
         DO I=1,ITMAX
          AN=-I*(I-A)
          B=B+2.0D0
          D=AN*D+B
          IF(ABS(D).LT.FPMIN) D=FPMIN
          C=B+AN/C
          IF(ABS(C).LT.FPMIN) C=FPMIN
          D=1.0D0/D
          DEL=D*C
          H=H*DEL
          IF(ABS(DEL-1.0D0).LT.EPS) GOTO 20
         END DO
         RETURN
 20      DGAMIN=1.0D0-EXP(-X+A*LOG(X)-GLN)*H
      END IF
      END

      DOUBLE PRECISION FUNCTION DGAMML(X)       ! ln[Gamma(x)]
      IMPLICIT NONE
      DOUBLE PRECISION X
      INTEGER J
      DOUBLE PRECISION SER,STP,TMP,XX,YY,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,-86.50532032941677D0,
     +24.01409824083091D0,-1.231739572450155D0,0.1208650973866179D-2,
     +-.05395239384953D-5,2.5066282746310005D0/
*     ...
      XX=X
      YY=XX
      TMP=(XX+0.5D0)*LOG(XX+5.5D0)-XX-5.5D0
      SER=1.000000000190015D0
      DO J=1,6
        YY=YY+1.0D0
        SER=SER+COF(J)/YY
      END DO
      DGAMML=TMP+LOG(STP*SER/XX)
      END
