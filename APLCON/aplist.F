      SUBROUTINE APPULL(PULLS)
      DOUBLE PRECISION PULLS(*)
#include "comcfit.inc"
#include "nauxfit.inc" 
*     ...
      DO I=1,NX
       PULLS(I)=AUX(INDPU+I)
      END DO
      END


      
