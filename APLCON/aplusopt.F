
      SUBROUTINE ADUMMY
*     __________________________________________________________________
      IMPLICIT NONE
      DOUBLE PRECISION STEP,XLOW,XHIG,POW,ARG
      INTEGER IT,LUNP,JPR,NBINOM,I1,I2
#include "comcfit.inc"
#include "nauxfit.inc"
#include "declarefl.inc"
*     __________________________________________________________________
*
*     parameters for fit method
*     __________________________________________________________________
      ENTRY APRINT(LUNP,JPR)
      LUNSIM=LUNP       ! print unit
      IF(LUNSIM.LE.0)     LUNSIM=6
      IPR=JPR           ! print flag
      RETURN
*     __________________________________________________________________
      ENTRY APDEPS(ARG)                 ! constraint accuracy
      EPSF  =ARG        ! |F| accuracy
      RETURN
*     __________________________________________________________________
      ENTRY APEPSCHI(ARG)               ! chi2 accuracy
      EPSCHI=ARG        
      RETURN      
*     __________________________________________________________________
      ENTRY APDERF(ARG)                 ! factor for step definition
      DERFAC=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APDERU(ARG)                 ! factor for step definition
      DERUFC=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APDLOW(ARG)                 ! factor for step definition
      DERLOW=ARG
      RETURN
*     __________________________________________________________________
      ENTRY APITER(IT)                  ! iteration limit
      ITERMX=MAX(3,IT)  ! max number of iterations
      RETURN
    
*     __________________________________________________________________      
      ENTRY APSTEP(I,STEP)              ! step size for numdif
      IF(I.LT.1.OR.I.GT.NX) RETURN 
      IPAK=I
#include "unpackfl.inc"     
      AUX(INDST+I)=ABS(STEP)            ! ST(I)= ...
      IF(STEP.NE.0.0D0) THEN
         NTINE=0 ! variable
      ELSE
         NTINE=1 ! fixed by user
      END IF
      GOTO 100
*     __________________________________________________________________
      ENTRY APFIX(I)                    ! fixed parameter
      IF(I.LT.1.OR.I.GT.NX) RETURN 
      IPAK=I
#include "unpackfl.inc"
      NTINE=1    ! fixed by user 
      GOTO 100
*     __________________________________________________________________
      ENTRY APLIMT(I,XLOW,XHIG)        ! range of variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      AUX(INDLM+2*(I-1)+1)=MIN(XLOW,XHIG) ! lower limit XL(1,I)
      AUX(INDLM+2*(I-1)+2)=MAX(XLOW,XHIG) ! upper limit XL(2,I)
      NTLIM=4 
      GOTO 100
*     __________________________________________________________________
      ENTRY APTRIN(I)                   ! inverse value
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=1    ! transformation to inverse
      GOTO 100
*     __________________________________________________________________
      ENTRY APOISS(I)                   ! Poisson distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=2    ! Poisson distributed variable
      NTLIM=1    
c      WRITE(*,*) 'APOISS I,IPAK,NTVAR,NTLIM ',I,IPAK,NTVAR,NTLIM    
      GOTO 100
*     __________________________________________________________________
      ENTRY ABINOM(I,NBINOM)           ! Binomial distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=3    ! Binomial distributed variable
      AUX(INDLM+2*I)=NBINOM
      NTLIM=1  
      GOTO 100
*     __________________________________________________________________
      ENTRY APLOGN(I)                  ! Lognormal distributed variable
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=4    ! Lognormal distributed variable
      NTLIM=1
      GOTO 100
*     __________________________________________________________________
      ENTRY APSQRT(I)                  ! SQRT transformation
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=5    ! SQRT transformation
      NTLIM=1 
      GOTO 100
*     __________________________________________________________________
      ENTRY APOWER(I,POW)              ! x^power transformation
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTVAR=6    ! x^power transformation
      AUX(INDLM+2*I)=POW
      NTLIM=1 
      GOTO 100
*     __________________________________________________________________
      ENTRY APOSIT(I)                  ! positive 
      IF(I.LT.1.OR.I.GT.NX) RETURN
      IPAK=I
#include "unpackfl.inc"
      NTLIM=1    ! positive
      GOTO 100
*     __________________________________________________________________
 100  CONTINUE
#include "packfl.inc"
      END 

      SUBROUTINE APSTAT(FOPT,NFUN,NITER)          ! return Fopt and Nfun      
*     __________________________________________________________________
*     return information after the fit
*     __________________________________________________________________
      IMPLICIT NONE
#include "comcfit.inc"
      INTEGER NFUN,NITER
      DOUBLE PRECISION FOPT
*     ...
      FOPT =CHISQ
      NFUN =NCALLS
      NITER=ITER 
      END 




