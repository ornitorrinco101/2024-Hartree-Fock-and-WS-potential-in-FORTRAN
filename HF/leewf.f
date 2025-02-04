      implicit real*8 (a-h,o-z)
      dimension nlj(100,6),bkw(100,10),ener(100),ws(100,100)
      do isp=1,6
      read (5,1000,end=10) (nlj(isp,j),j=1,5),nlj(isp,6),ener(isp)
      read (5,1001) (bkw(isp,j),j=1,nlj(isp,6))
      enddo
 1000 FORMAT (6X,5I1,I5,f10.4)
 1001 FORMAT (1x,5D14.7)
c
c
 10   nsps=6
      b=1.6d0
c
c
C
C#################
      CALL OAINIT
C#################
c
      DO 30 J=1,100
      X=float(j-1)*.1
C
      DO 40 I=1,NSPS
      LL=NLJ(I,4)
      WF=0.D0
      DO 50 K=1,NLJ(I,6)
      WF=WF+bkw(i,k)*RAD(K-1,LL,X,B)
   50 CONTINUE
C
      WS(I,J)=WF
   40 CONTINUE
   30 CONTINUE
C
      do i=1,100
      X=float(i-1)*.1
      write (6,1004) x,(ws(k,i),k=1,nsps)
 1004 format(1x,f6.2,2x,6(e12.5,1x))
      enddo
      stop
      end

c-----------------------------------------------------
      SUBROUTINE OAINIT
C
C-----------------------------------------------
C      SUBRUTINA DE INICIALIZACION PARA EL USO
C      DE LAS FUNCIONES RADIALES DE OSCILADOR
C      Y SUS DERIVADAS
C
C      VERSION EN DOBLE PRECISION
C
C      CALCULA:
C      F(K)= (K-1)!
C      G(K)= GAMMA((2*K-1)/2)
C      K=1,31
C-----------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(31),G(31)
C
      COMMON /FACGAM/ F,G
C
      F(1)=1.D0
      G(1)=DSQRT(4.D0*DATAN(1.D0))
      DO 10 I=1,30
      XI=I
      F(I+1)=F(I)*XI
      G(I+1)=G(I)*(XI-.5D0)
   10 CONTINUE
C
      RETURN
      END
      FUNCTION RAD (N,L,X,B)
C
C------------------------------------------------
C      CALCULA LA FUNCION RADIAL DE OSCILADOR
C      ARMONICO DE PARAMETRO DE OSCILADOR B
C      EN EL PUNTO X [FERMI].
C      N=0,1,...,30
C      L=0,1,...,30
C
C      VERSION DOBLE PRECISION
C
C      !! LA SUBRUTINA .OAINIT_D. DEBE LLAMARSE
C      ANTES DE UTILIZAR ESTA FUNCION !!
C------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(31),G(31)
C
      COMMON /FACGAM/ F,G
C
      RAD=0.D0
      AN=DSQRT(2.D0*F(N+1)*G(L+N+2)/(B*B*B))
      IF (X-1.D-8) 10,20,20
   10 IF (L) 30,40,30
   40 RAD=AN/(G(2)*F(N+1))
      GOTO 30
   20 Y=X/B
      NN=N+1
      C=0.D0
      DO 50 K=1,NN
      C=-C+(Y**(K+K-2))/(F(K)*F(N-K+2)*G(L+K+1))
   50 CONTINUE
      P=1+2*(N/2*2-N)
      RAD=AN*P*C*DEXP(-.5D0*Y*Y)*(Y**L)
C
   30 RETURN
      END
