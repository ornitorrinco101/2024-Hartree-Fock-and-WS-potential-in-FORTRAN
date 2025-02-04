C
C     THIS PROGRAM INTEGRATES THE SCHRODINGER EQUATION
C     FOR THE WOODS-SAXON POTENTIAL USING NUMEROV METHOD
C     THE CENTRIFUGAL KERNEL IS USED TO DETERMINE THE
C     EINGENVALUES AND SINGLE-PARTICLE WAVE FUNCTIONS.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      CALL INPUT
      CALL CONST
      CALL SCHRO
      CALL ESCRI
c      CALL TEST
C
      STOP
      END
C
C--------------------------------------------------------
      SUBROUTINE INPUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE       
      COMMON /COLA/ ICOLA
C
      READ (5,*)
      READ (5,*) IA,IZ,ICOLA
      WRITE (6,*) IA,IZ,ICOLA
      READ (5,*)
      READ (5,*) (VPARAM(2,J),J=1,7)
      WRITE (6,*) (VPARAM(2,J),J=1,7)
      READ (5,*) (VPARAM(1,J),J=1,7)
      WRITE (6,*) (VPARAM(1,J),J=1,7)
      READ (5,*)
      READ (5,*) RMAX,NPOIN,NPWF
      WRITE (6,*) RMAX,NPOIN,NPWF
      READ (5,*)
      READ (5,*) STEP,EMAX,EERR
      WRITE (6,*) STEP,EMAX,EERR
      READ (5,*)
      READ (5,*) IWRITE
      WRITE (6,*) IWRITE
C
      H=RMAX/DFLOAT(NPOIN)
C
      RETURN
      END
C
C--------------------------------------------------------
      SUBROUTINE CONST
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2) 
      COMMON /MASAS/ XNMAS(2),PIMAS
C
C     NUMERICAL CONSTANTS
C
      HBARC=197.3285151D0
      ALFA=1.D0/137.0359821D0
      XNMAS(1)=939.5728D0
      XNMAS(2)=938.2592D0      
      PIMAS=139.5669D0
C
C     EVALUATION OF THE POTENTIAL CONSTANTS
C
      DO ISP=1,2
        CENER(ISP)=2.D0*XNMAS(ISP)/HBARC/HBARC
        CCENT(ISP)=CENER(ISP)*VPARAM(ISP,1)
        CSO(ISP)=2.D0*XNMAS(ISP)*VPARAM(ISP,4)/VPARAM(ISP,6)/PIMAS/PIMAS
      ENDDO
      CCOUL(1)=0.D0              
      CCOUL(2)=2.D0*DFLOAT(IZ-1)*ALFA*XNMAS(2)/HBARC
C
      RETURN
      END
C
C--------------------------------------------------------
      SUBROUTINE SCHRO
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      DIMENSION V(3000),IZERO(100),WF(3000)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /DATANG/ ISP,LL,JJ,WW 
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2)
      COMMON /EWAVE/  NLJ(100,5),NRST(2),WAVE(100,3000),
     *                VPOT(100,3000),EE(100)
C
      DO 10 ISP=2,1,-1
        NRST(ISP)=0
C
        LL=0
   20   WW=2.D0*DBLE(FLOAT(LL))+1.D0
C
        DO 30 JJ=LL+1,LL,-1
          IF (LL.EQ.0.AND.JJ.EQ.LL) GOTO 30

          CALL WSPOT (V)
C
          NN=1
          E0=VPARAM(ISP,1)
          CALL WSEING (E0,WF)
          Y00=WF(1)
C
          DO IENE=1,10000
            IF (E0.GE.EMAX) GOTO 40
            E1=E0+STEP
            IF (E1.GT.0.0) GOTO 40
            CALL WSEING (E1,WF)
            Y01=WF(1)
C            IF(JJ.EQ.2.AND.LL.EQ.2.AND.NN.EQ.1.AND.ISP.EQ.1) THEN
C             WRITE(6,*) E0,E1
C             WRITE(6,*) Y00,Y01
C            ENDIF
            IF (Y00/Y01.LT.0.D0) THEN
              IZERO(NN)=IENE-1
              NN=NN+1
            ENDIF
            E0=E1
            Y00=Y01
          ENDDO
C
   40     IF (NN.EQ.1.AND.JJ.EQ.LL+1) GOTO 10
C
          DO IND=1,NN-1
            E0=VPARAM(ISP,1)+DFLOAT(IZERO(IND))*STEP
            E1=E0+STEP
   50       IF (DABS(E1-E0).LT.EERR) THEN
              IF (E1.GT.EMAX) GOTO 30
              NRST(ISP)=NRST(ISP)+1
              INDIC=NRST(ISP)
              IF (ISP.EQ.1) INDIC=INDIC+NRST(2)
              NLJ(INDIC,2)=ISP-1
              NLJ(INDIC,3)=IND
              NLJ(INDIC,4)=LL
              NLJ(INDIC,5)=JJ
 111          FORMAT (1X,3I4,2D13.5)
              EE(INDIC)=E1
              WRITE (6,111) IND,LL,JJ,E1
              CALL XNORM(INDIC,E1,WF)
              DO K=1,NPOIN+1
                WAVE(INDIC,K)=WF(K)
                VPOT(INDIC,K)=V(K)
              END DO
C
            ELSE
C-- ENERGY CORRECTION -------------------------
C
              CALL WSEING(E0,WF)
              Y00=WF(1)
              CALL WSEING(E1,WF)
              Y01=WF(1)
              ECOR=(E0*Y01-E1*Y00)/(Y01-Y00)
              E0=E1
              Y00=Y01
              E1=ECOR
              GOTO 50
            ENDIF
          ENDDO
C
   30     CONTINUE
C
          LL=LL+1
          GOTO 20
C
   10   CONTINUE            
C
          CALL ORDEN
C
        RETURN
        END
C
C--------------------------------------------------------
      SUBROUTINE WSPOT (V)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      DIMENSION V(3000),U1(3),U2(3)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /DATANG/ ISP,LL,JJ,WW 
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2)
      COMMON /DIAGON/ EP(3000),EC(3000),EMN,
     *                CYP(3000),CYC(3000),CYM(3000)
C        
      HT=H*H/12.D0
C
      DO 10 I=2,NPOIN+1
        X=FLOAT(I-1)*H
        VCENT=CCENT(ISP)/(1.D0+DEXP((X-VPARAM(ISP,2))/VPARAM(ISP,3)))
        ARGSO=DEXP((X-VPARAM(ISP,5))/VPARAM(ISP,6)/2.D0)
        ARGSO=ARGSO+1.D0/ARGSO
        VSO=CSO(ISP)*DFLOAT(JJ*JJ-LL*LL-LL-1)/X/ARGSO/ARGSO
        VCOUL=CCOUL(ISP)/X
        IF (X.GE.VPARAM(ISP,7)) GOTO 20
        XRCH=X/VPARAM(ISP,7)
        VCOUL=VCOUL/2.D0*(3.D0*XRCH-XRCH*XRCH*XRCH)
C
   20   V(I)=VCENT+VSO+VCOUL
        EP(I)=V(I)
C
   10   CONTINUE
C
C       ------------------------------------------
C       EVALUATION OF THE COEFFICIENTS CYP,CYC,CYM
C
        DO 30 I=2,3
        X=FLOAT(I-1)*H
        U1(I)=X**(LL+1)
        U2(I)=(X/RMAX)**(LL+1)/RMAX**LL-1.D0/X**LL
   30   CONTINUE
C
        CYP(2)=-U1(2)
        CYC(2)=U1(3)
        CYM(2)=U1(2)*U2(3)-U1(3)*U2(2)
C
        EC(2)=CYC(2)*HT+H*CYM(2)*CYP(2)/WW
        CYC(2)=H*CYM(2)*CYP(2)*EP(2)/WW-CYC(2)*(1.D0-HT*EP(2))
        EP2=CYP(2)*HT
        CYP(2)=CYP(2)*(1.D0-HT*EP(3))
C
        DO 40 I=3,NPOIN
        X=FLOAT(I)*H
C
        DO 50 J=1,2
        U1(J)=U1(J+1)
        U2(J)=U2(J+1)
   50   CONTINUE
C
        U1(3)=X**(LL+1)
        U2(3)=(X/RMAX)**(LL+1)/RMAX**LL-1.D0/X**LL
C
        CYP(I)=U1(1)*U2(2)-U1(2)*U2(1)
        CYC(I)=U1(3)*U2(1)-U1(1)*U2(3)
C
   40   CONTINUE
C
        DO 60 I=3,NPOIN
        EC(I)=H*CYP(I+1)*CYP(I)/WW+HT*CYC(I)
        CYC(I)=H*CYP(I+1)*CYP(I)*EP(I)/WW-
     *         CYC(I)*(1.D0-HT*EP(I))
        EP3=HT*CYP(I)
        CYP(I)=CYP(I)*(1.D0-HT*EP(I+1))
        CYM(I)=CYP(I+1)*(1.D0-HT*EP(I-1))
        EP(I-1)=EP2
        EP2=EP3
   60   CONTINUE
C
        U1(1)=U1(2)
        U1(2)=U1(3)
        U2(1)=U2(2)
        U2(2)=U2(3)
        X=FLOAT(NPOIN+1)*H
        U1(3)=X**(LL+1)
        U2(3)=(X/RMAX)**(LL+1)/RMAX**LL-1.D0/X**LL
C
        CYMN=U1(1)*U2(2)-U1(2)*U2(1)
        EMN=CYMN*HT
        EC(NPOIN)=H*CYMN*CYP(NPOIN)/WW+CYC(NPOIN)*HT
        EP3=CYP(NPOIN)*HT
        CYC(NPOIN)=H*CYMN*CYP(NPOIN)*EP(NPOIN)/WW-
     *             CYC(NPOIN)*(1.D0-HT*EP(NPOIN-1))
        CYP(NPOIN)=CYP(NPOIN)*(1.D0-HT*EP(NPOIN+1))
        CYM(NPOIN)=CYMN*(1.D0-HT*EP(NPOIN-1))
        EP(NPOIN-1)=EP2
        EP(NPOIN)=EP3
C
        RETURN
        END
C
C--------------------------------------------------------
        SUBROUTINE WSEING (EE,WF)
C
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
        DIMENSION WF(3000)
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /DATANG/ ISP,LL,JJ,WW 
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2) 
      COMMON /DIAGON/ EP(3000),EC(3000),EMN,
     *                CYP(3000),CYC(3000),CYM(3000)
      COMMON /COLA/ ICOLA
C
        VENER=CENER(ISP)*EE
        DO I=1,NPOIN+1
          WF(I)=0.D0
        ENDDO
C
        IF (ICOLA.EQ.0) THEN 
          YP=0.D0
          YC=1.D0
        ENDIF
        IF (ICOLA.EQ.1) THEN 
          GAMMA=DSQRT(DABS(VENER))
          YP=DEXP(-GAMMA*RMAX)
          YC=YP*DEXP(H*GAMMA)
C          YC=1.D0       
        ENDIF
C
        YM=(CYC(NPOIN)-EC(NPOIN)*VENER)*YC/(CYM(NPOIN)+EMN*VENER)-
     *     (CYP(NPOIN)+EP(NPOIN)*VENER)*YP/(CYM(NPOIN)+EMN*VENER)
        IF (ICOLA.EQ.0) THEN
          WF(NPOIN+1)=0.D0
          WF(NPOIN)=1.D0
        ENDIF
        IF (ICOLA.EQ.1) THEN
          WF(NPOIN+1)=DEXP(-GAMMA*RMAX)
          WF(NPOIN)=WF(NPOIN+1)*DEXP(H*GAMMA)
C          WF(NPOIN)=1.D0
        ENDIF
        WF(NPOIN-1)=YM
C
        DO 10 I=NPOIN-1,3,-1
C
        YP=YC
        YC=YM
        YM=(CYC(I)-EC(I)*VENER)*YC/(CYM(I)+EP(I+1)*VENER)-
     *     (CYP(I)+EP(I)*VENER)*YP/(CYM(I)+EP(I+1)*VENER)
        WF(I-1)=YM
C
   10   CONTINUE
C
        Y0=(CYC(2)-EC(2)*VENER)*YM-(CYP(2)+EP(2)*VENER)*YC
        WF(1)=Y0
C
        RETURN
        END
C
C--------------------------------------------------------
        SUBROUTINE XNORM (INDI,EE,WF)
C
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
        DIMENSION WF(3000),DF(3000)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2) 
      COMMON /DATANG/ ISP,LL,JJ,WW 
      COMMON /COLA/ ICOLA
      COMMON /NORM/ UNORM(100)
C
        DO I=0,NPOIN
        DF(I)=WF(I+1)*WF(I+1)
        END DO
C 
        IF (ICOLA.EQ.0) TAIL=0.D0
        IF (ICOLA.EQ.1) THEN 
          GAMMA=DSQRT(DABS(EE*CENER(ISP)))
          XGAM=DEXP(GAMMA*RMAX)
          TAIL=1/2/GAMMA/XGAM/XGAM
        ENDIF
        UNORM(INDI)=SMP(0.D0,NPOIN,H,DF)+TAIL
C
        SSI=1.D0
        IF (WF(2).LT.0.D0) SSI=-1.D0
        DO I=2,NPOIN+1
        R=H*DFLOAT(I-1)
        WF(I)=SSI*WF(I)/DSQRT(UNORM(INDI))/R
        END DO
C
        RETURN
        END
C
C--------------------------------------------------------
      FUNCTION SMP(F0,IO,DR,WF)
C
      IMPLICIT REAL*8 (A-H,O-Z) 
C
      DIMENSION WF(1)
C
      S=1/3*(F0+WF(IO))
      TT=2.0/3.0
      DO 1 J1=1,IO-1
      W=1+J1-2*(J1/2)
      WEI=TT*W*DR
      S=S + WF(J1)*WEI
 1    CONTINUE
      SMP=S
      RETURN
      END
C
C----------------------------------------------------------
      SUBROUTINE ORDEN
C
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
C
      DIMENSION V(3000),IZERO(100),WF(3000)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE        
      COMMON /EWAVE/  NLJ(100,5),NRST(2),WAVE(100,3000),
     *                VPOT(100,3000),EE(100)
C
      DO ISP=2,1,-1
C
        IND0=1-NRST(2)*(ISP-2)
        IND1=NRST(ISP)-1-NRST(2)*(ISP-2)
C
        DO 100 INDIC=IND0,IND1
C
        DO 110 INDIZ=INDIC+1,IND1+1
        IF (EE(INDIZ).GT.EE(INDIC)) GOTO 110
        AEN=EE(INDIZ)
        EE(INDIZ)=EE(INDIC)
        EE(INDIC)=AEN
C
        DO J=2,5
        IAUX=NLJ(INDIZ,J)
        NLJ(INDIZ,J)=NLJ(INDIC,J)
        NLJ(INDIC,J)=IAUX
        ENDDO
C
        DO J=1,NPOIN+1
        WAUX=WAVE(INDIZ,J)
        WAVE(INDIZ,J)=WAVE(INDIC,J)
        WAVE(INDIC,J)=WAUX
        VAUX=VPOT(INDIZ,J)
        VPOT(INDIZ,J)=VPOT(INDIC,J)
        VPOT(INDIC,J)=VAUX
        ENDDO
C
  110   CONTINUE
C
  100   CONTINUE
C
        IEM=0
        IPMAX=-IA*(ISP-2)+(-1)**ISP*IZ
        DO 130 I=IND0, IND1+1
        IEM=IEM+2*NLJ(I,5)
        NLJ(I,1)=1
        IF (IEM.LE.IPMAX) NLJ(I,1)=0
  130   CONTINUE
C
      ENDDO
C
        RETURN
        END
C
C--------------------------------------------------------
        SUBROUTINE ESCRI
C
        IMPLICIT REAL*8 (A-H,O-Z)
        IMPLICIT INTEGER*4 (I-N)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE   
      COMMON /EWAVE/  NLJ(100,5),NRST(2),WAVE(100,3000),
     *                VPOT(100,3000),EE(100)
C
        WRITE (*,1004) IZ,IA,((VPARAM(I,J),J=1,7),I=2,1,-1)
        WRITE (*,1005) RMAX,NPOIN,STEP,EMAX,EERR
        WRITE (*,1002) (NRST(I),I=2,1,-1)
        DO  I=1,NRST(2)
        WRITE (*,1003) (NLJ(I,J),J=1,5),EE(I)
        END DO
        WRITE (*,*)
        DO  I=NRST(2)+1,NRST(2)+NRST(1)
       WRITE (*,1003) (NLJ(I,J),J=1,5),EE(I)
       END DO
       WRITE (*,*)
       WRITE (*,*)
C
       IF(IWRITE.EQ.1) THEN
C
        OPEN(2,FILE='FDO-PROT.WS',STATUS='UNKNOWN')
        WRITE(2,1111) ((NLJ(I,J),J=1,5),I=1,NRST(2))
 1111   FORMAT('#######',2X,80(4(1X,I1),1X,I2,1X))
        DO I=2,NPOIN+1,10
         R=H*DFLOAT(I-1)
         WRITE (2,1001) R,(WAVE(INDIC,I),INDIC=1,NRST(2))
       END DO
c
        CLOSE(2)
C
        OPEN(2,FILE='FDO-NEUT.WS',STATUS='UNKNOWN')
        WRITE(2,1111) ((NLJ(I,J),J=1,5),I=NRST(2)+1,NRST(2)+NRST(1))
        DO I=2,NPOIN+1,10
         R=H*DFLOAT(I-1)
         WRITE (2,1001) R,(WAVE(INDIC,I),
     >                  INDIC=NRST(2)+1,NRST(2)+NRST(1))
        END DO
        CLOSE(2)
C
C        DO IND=1,NRST(2),6
C         INDEND=IND+5
C         IF (NRST(2).LT.IND+5) INDEND=NRST(2)
C          DO I=2,NPOIN+1,10
C         R=H*DFLOAT(I-1)
C         WRITE (2,1001) R,(WAVE(INDIC,I),INDIC=IND,INDEND)
C        END DO
C        WRITE (2,*)
C       END DO

C        DO IND=NRST(2)+1,NRST(2)+NRST(1),6
C        INDEND=IND+5
C        IF (NRST(2)+NRST(1).LT.IND+5) INDEND=NRST(2)+NRST(1)
C        DO I=2,NPOIN+1,10
C        R=H*DFLOAT(I-1)
C        WRITE (2,1001) R,(WAVE(INDIC,I),INDIC=IND,INDEND)
C        END DO
C      WRITE (2,*)
C      END DO

      ENDIF       
C
    1   FORMAT (A20)
C 1001   FORMAT (1X,F6.3,2X,6(D10.3,1X))
 1001   FORMAT (1X,F6.3,2X,80(e11.4,1X))
 1002   FORMAT (1X,'NR. PROTON STATES=',I3,2X,
     *          'NR. NEUTRON STATES=',I3,//)
 1003   FORMAT (1X,5(I2,1X),F7.3)
 1004   FORMAT (//,1X,'Z=',I2,2X,'A=',I3,//,7F9.4,/,7F9.4)
 1005   FORMAT (/,1X,'RMAX=',F6.3,'  NPOIN=',I4,2X,'STEP=',F4.1,
     *          2X,'EMAX=',F9.4,2X,'EERR=',D7.1,//)
C
        RETURN
        END

      SUBROUTINE TEST
C
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT INTEGER*4 (I-N)
C
        DIMENSION FF(100,3000),DF(3000),ERED(100),F(3000)
        DIMENSION FP(3000)
        DIMENSION GAMMA(100),XGAM(100)
C
      COMMON /DATPOT/ IA,IZ,VPARAM(2,7),
     *                RMAX,NPOIN,H,NPWF,
     *                STEP,EMAX,EERR,IWRITE 
      COMMON /CTEPOT/ CENER(2),CCENT(2),CSO(2),CCOUL(2) 
      COMMON /EWAVE/  NLJ(100,5),NRST(2),WAVE(100,3000),
     *                VPOT(100,3000),EE(100)
      COMMON /COLA/ ICOLA
      COMMON /NORM/ UNORM(100)
      COMMON /MASAS/ XNMAS(2),PIMAS
C
C
        DO 10 INDI=1,NRST(1)+NRST(2)
        ISP=NLJ(INDI,2)+1
C
        ERED(INDI)=CENER(ISP)*EE(INDI)
        GAMMA(INDI)=DSQRT(DABS(ERED(INDI)))
        XGAM(INDI)=DEXP(-GAMMA(INDI)*RMAX)
        LL=NLJ(INDI,4)
        JJ=NLJ(INDI,5)
C
        FF(INDI,1)=0.D0
        DO K=2,NPOIN+1
        X=FLOAT(K-1)*H
C        
        CCC=X*X*(ERED(INDI)-VPOT(INDI,K))-DFLOAT(LL*LL+LL)
        FF(INDI,K)=CCC*WAVE(INDI,K)
        ENDDO
C
   10   CONTINUE
C
        DO 140 I=1,NRST(2)+NRST(1)-1
        DO 150 J=I+1,NRST(2)+NRST(1)
C
c        I=7
c        J=8
        DO K=1,NPOIN
c          X=K*H
          DF(K)=FF(I,K+1)*WAVE(J,K+1)-FF(J,K+1)*WAVE(I,K+1)
c          XARG=DEXP(-(X-VPARAM(1,5))/VPARAM(1,6))
c          XARG=XARG/X/(1+XARG)**2
c          PPP=-6*VPARAM(1,4)*XNMAS(1)*XARG/(VPARAM(1,6)*PIMAS**2)
c     *        +ERED(7)-ERED(8)        
c          FP(K)=PPP*WAVE(7,K+1)*WAVE(8,K+1)*X**2
C          F(K)=(ERED(6)-ERED(10))*WAVE(6,K+1)*WAVE(10,K+1)*X**2
c        PRINT *, K,DF(K),FP(K)
        END DO
C 
        IF (ICOLA.EQ.0) TAIL=0.D0
        IF (ICOLA.EQ.1) THEN 
          TAIL=(GAMMA(J)-GAMMA(I))*XGAM(I)*XGAM(J)/
     *        (DSQRT(UNORM(J))*DSQRT(UNORM(I)))
        ENDIF
C        TSAL=TAIL+SMP(0.D0,NPOIN,H,DF)
        TSAL=SMP(0.D0,NPOIN,H,DF)
C        TSALP=TAIL+SMP(0.D0,NPOIN,H,F)
c        TTT=TAIL+SMP(0.D0,NPOIN,H,FP)
        WRITE(2,100) (NLJ(I,K),K=1,5),(NLJ(J,K),K=1,5),TSAL
  100   FORMAT (2X,4I1,I2,4X,4I1,I2,4X,D16.10)
C
  150   CONTINUE
  140   CONTINUE 
C
        RETURN
        END



