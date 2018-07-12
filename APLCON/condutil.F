
      INTEGER FUNCTION IJSYM(I,J) ! index (I,J)=(J,I) in symmetric matrix
      IMPLICIT NONE
      INTEGER I,J !,NCOUNT
      IF(I.LE.J) THEN
         IJSYM=(J*J-J)/2+I
      ELSE
         IJSYM=(I*I-I)/2+J
      END IF
      END

      SUBROUTINE DUMINV(A, W,B,NX,NF,MB, NRANK, AUX,QNEXT) ! matrix inversion
*     Obtain solution of a system of linear equations V *  X  =  B  with
*     symmetric matrix V and inverse (for M =  1)  or  matrix  inversion
*     only (for M = 0).
*
*                   - - - -
*        CALL SMINV(W,B,N,M,NRANK,AUX,QNEXT)
*                   - -     -----
*
*           W = symmetric N-by-N matrix in symmetric storage mode
*               W(1) = W11, W(2) = W12, W(3) = W22, W(4) = W13, . . .
*               replaced by inverse matrix
*           B = N-vector   (for M = 0 use a dummy argument)
*               replaced by solution vector
*           M = see above
*
*
*     Method of solution is by elimination selecting the  pivot  on  the
*     diagonal each stage. The rank of the matrix is returned in  NRANK.
*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
*     matrix V and the corresponding elements of  B  are  set  to  zero.
*     SMINV can be used for a dimension up to 100 (see INVCDR).
*
      IMPLICIT NONE
      DOUBLE PRECISION VKK 
      DOUBLE PRECISION W(*),B(*),A(*),AUX(*),QNEXT(*)
      DOUBLE PRECISION SUM,EPS,VJK
      INTEGER NX,NF
      INTEGER I,J,IJ,IA,JFIRST,JLAST,MB,NRANK,
     +        K,LAST,JK,NMEAS,N,JJ,L,KK,M,JL,LK  

      LOGICAL SOLVE
       
      DATA EPS/1.0D-6/

*     special entry for partial inversion ******************************

*     ...  
      NRANK=0 
      N=NX+NF                   ! dimension parameter 

*     make sure AUX is zero, prevents uninit access
*     in continued execution of DBMINV
      DO I=1,N
         AUX(I)=0.0D0
      END DO

*     -VX(NX-sym) is already inserted in W(NX+NF-sym) ------------------ 


      IJ=(NX*NX+NX)/2           ! number of elements of V(.) 

      DO I=1,N
       QNEXT(I)=0.0D0           ! reset pointer
      END DO

      IA=0
      DO J=1,NF
       DO I=1,NX
        W(IJ+I)=A(IA+I)         ! copy A(.) into W_12
       END DO
       DO I=1,J
        W(IJ+NX+I)=0.0          ! reset last submatrix W_22 of W(.)
       END DO
       IJ=IJ+NX+J 
       IA=IA+NX
      END DO

*     distinguish between measured and unmeasured variables ------------

      JFIRST=0                  ! first index of measured variable 
      NMEAS=0                   ! number of measured variables
      DO I=1,NX
       IF(W((I*I+I)/2).LT.0.0) THEN  ! measured variable
          IF(JFIRST.EQ.0) THEN
             JFIRST=I           ! first index of measured variable
          ELSE
             QNEXT(JLAST)=I     ! insert index at previous index
          END IF
          JLAST=I               ! save index
          NMEAS=NMEAS+1
       END IF
      END DO
      IF(JLAST.EQ.0) GOTO 10    ! nothing to do
      QNEXT(JLAST)=-1           ! stop index for last measured variable 


*     apply exchange algorithm to sub-matrices -------------------------

      DO I=NX+1,N               ! loop I over constraint equations

       J=JFIRST                 ! first index of unmeasured variable
       DO M=1,NMEAS             ! already inverted element index J
        SUM=0.0D0
        JK=(J*J-J)/2            ! index of diagonal element before
        DO K=1,NX
         IF(K.LE.J) JK=JK+1     ! index in j column
         IF(QNEXT(K).NE.0.0D0) SUM=SUM+W(JK)*W((I*I-I)/2+K)
         IF(K.GE.J) JK=JK+K     ! index in j row
        END DO
        AUX(J)=SUM              ! = A-row * VX-row/col
        J=QNEXT(J)              ! next index of unmeasured variable
       END DO

       DO K=I,N 
        SUM=0.0D0
        J=JFIRST                ! first index of unmeasured variable 
        DO M=1,NMEAS            ! already inverted element index J
         SUM=SUM+W((K*K-K)/2+J)*AUX(J) ! = A-row * H
         J=QNEXT(J)             ! next index of unmeasured variable
        END DO
        W((K*K-K)/2+I)=W((K*K-K)/2+I)+SUM   ! add to diagonal W_22
       END DO

       J=JFIRST                 ! first index of unmeasured variable
       DO M=1,NMEAS
        W((I*I-I)/2+J)=-AUX(J)  ! add to off-diagonal W_22
        J=QNEXT(J)              ! next index of unmeasured variable
       END DO

      END DO

*     set pointer for unmeasured variables ----------------------------- 

      JFIRST=0
      JLAST=0
      DO I=1,N                  
       IF(QNEXT(I).EQ.0.0D0) THEN ! unmeasured variable
          IF(JFIRST.EQ.0) THEN  
             JFIRST=I           ! first index of unmeasured variable
          ELSE
             QNEXT(JLAST)=I     ! next index of unmeasured variable
          END IF
          JLAST=I
       ELSE
          QNEXT(I)=0.0D0        ! reset index for measured variable
       END IF 
      END DO
      IF(JLAST.EQ.0) GOTO 10    ! no unmeasured variable
      QNEXT(JLAST)=-1           ! end flag

      
*     common code for inversion and (M=1) solution of matrix equation  

 10   SOLVE=.TRUE. 
      IF(MB.EQ.0) SOLVE=.FALSE.  ! solution flag

*     loop begin (loop on all remaining rows/cols)

      DO I=1,N                  ! loop on all remaining elements
       VKK=0.0D0                ! search for pivot element
       K=0                      ! pivot index
       J=JFIRST                 ! first candidate index
       LAST=0
 20    IF(J.GT.0) THEN          ! test for linearity and zero matrix 
          JJ=(J*J+J)/2          ! diagonal index
          IF(ABS(W(JJ)).GT.MAX(ABS(VKK),EPS*AUX(J))) THEN
             VKK=W(JJ)          ! largest pivot candidate so far
             K=J                ! index of largest
             L=LAST
          END IF
          LAST=J
          J=QNEXT(J)             ! index of next candidate
          GOTO 20
       END IF
       IF(K.NE.0) THEN          ! pivot element found - proceed
          NRANK=NRANK+1         ! increase rank counter
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             JFIRST=QNEXT(K)    ! new first index
          ELSE
             QNEXT(L)=QNEXT(K)  ! bridge used index
          END IF
          QNEXT(K)=0.0D0        ! reset used index
          NRANK=NRANK+1         ! increase rank

          VKK    =1.0/VKK       ! invert pivot 
          W(KK)  =-VKK          
          IF(SOLVE) B(K)=B(K)*VKK
          JK     =KK-K
          JL     =0
          DO J=1,N              ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =W(JK)
              W(JK)=VKK*VJK
              IF(SOLVE) B(J) =B(J)-B(K)*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  W(JL)=W(JL)-W(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE                     ! no pivot candadate found - reset
          DO K=1,N
           IF(QNEXT(K).NE.0.0D0) THEN   ! undefined variable
              IF(SOLVE) B(K)=0.0D0 ! clear undefined vector element
              DO J=1,K
               IF(QNEXT(J).NE.0.0D0) W((K*K-K)/2+J)=0.0D0  ! clear matrix row/col
              END DO
           END IF
          END DO
          GOTO 30
       END IF

      END DO                    ! end of inversion loop

 30   DO I=1,(N*N+N)/2
       W(I)=-W(I)               ! finally reverse sign
      END DO
      END 


      DOUBLE PRECISION FUNCTION SCALXY(X,Y,N) ! scalar vector product
*     Scalar product of two vectors
*                - - -                 T
*        S = VXY(X,Y,N)           S = X  * Y (scalar product)
*
      IMPLICIT NONE
      INTEGER N,J ! I,M
      DOUBLE PRECISION X(*),Y(*)
      DOUBLE PRECISION SUM
*     ...
      SUM=0.0D0
      DO J=1,N
       SUM=SUM+X(J)*Y(J)
      END DO
      SCALXY=SUM
      END
