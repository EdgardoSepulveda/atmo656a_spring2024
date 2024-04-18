
c     PROGRAM GETPHASE
  
      COMMON / SIZDIS / NEWSD,         RGV,           SIGMAG,
     2                  RFRS,          RFIS,          RFRC,  
     3                  RFIC,          ALAMB,         IPHASE
 
      COMMON / STUFF /  OMEGA,         ASY,           EXT,   
     2                  SCAT,          QEXT,          QSCAT, 
     3                  QBS
      DIMENSION         COSPHI(181),       SCTPHS(181),  radius(251)
      PIE    = ACOS( -1.0 )
      TWOPI  = 2.0 * PIE
      
c      do 999 rr=0.01,10. ,0.1  ! for question 4
      radius(1) = 0.01
      do 999 my_int = 2, 251
      rr = radius(my_int-1)
c      do 999 rr=0.01,2.5,0.01
      if ( rr .eq. 0.01 ) then
        NEWSD = 1
      else
        NEWSD = 0
      end if
c      NEWSD  = 1    ! = 1 IF SIZE DIS VALUES HAVE NOT PREVIOUSLY BEEN USED. New Size Distribution?
c     RGV = 4.0       ! GEOMETRIC MEAN RADIUS FOR THE VOLUME DISTRIBUTION
      SIGMAG = 1.8    ! GEOMETRIC STANDARD DEVIATION
      RFRS   = 1.3873 ! REAL INDEX OF REFRACTION OF PARTICLES
      RFIS   = 0.0066 ! IMAGINARY INDEX OF REFRACTION OF PARTICLES
      RFRC   = RFRS
      RFIC   = RFIS
      ALAMB  = 0.5    ! WAVELENGTH AT WHICH CALCULATIONS ARE TO BE PERFORMED

      WNUM   = TWOPI / ALAMB
      RGV    = rr
      xx     = WNUM*rr
      radius(my_int) = rr + 0.01 
      CALL GETPHS( NSCATA, COSPHI, SCTPHS )
      WRITE(*,*) xx,QEXT
 999  continue

      STOP 99999
      END

C=====================================================
      SUBROUTINE GETPHS( NSCATA, COSPHI, SCTPHS )
C
C THIS SUBROUTINE COMPUTES THE PHASE FUNCTION SCTPHS(I) AT NSCATA
C ANGLES BETWEEN 0.0 AND 180.0 DEGREES SPECIFIED BY COSPHI(I) WHICH
C CONTAINS THE ANGLE COSINES. THESE VALUES ARE RETURNED TO FUNCTION
C PHASFN FOR AZIMUTHAL AVERAGING.
C INPUT DATA FOR THIS ROUTINE IS PASSED THROUGH COMMON /SIZDIS/
C AND CONSISTS OF THE FOLLOWING PARAMETERS
C NEWSD  = 1 IF SIZE DIS VALUES HAVE NOT PREVIOUSLY BEEN USED IN
C          THIS ROUTINE,  = 0 OTHERWISE.
C RGV    = GEOMETRIC MEAN RADIUS FOR THE VOLUME DISTRIBUTION OF THE
C          SPECIFIED PARTICLES
C SIGMAG = GEOMETRIC STANDARD DEVIATION
C RFR,I  = REAL AND IMAGINARY INDEX OF REFRACTION OF PARTICLES
C ALAMB  = WAVELENGTH AT WHICH CALCULATIONS ARE TO BE PERFORMED
C
  
      COMMON / SIZDIS / NEWSD,         RGV,           SIGMAG,
     2                  RFRS,          RFIS,          RFRC,  
     3                  RFIC,          ALAMB,         IPHASE
 
      COMMON / STUFF /  OMEGA,         ASY,           EXT,   
     2                  SCAT,          QEXT,          QSCAT, 
     3                  QBS
 
      Real*8             THETA(91),        ELTRMX(4,91,2),    
     2                   PI(3,91),         TAU(3,91),         
     3                   CSTHT(91),        SI2THT(91)

      Real*8   rout, rfro, rfio, dqext, dqscat, ctbrqs, dqbs, 
     2         rin,  rfri, rfii, wnum
 
      COMPLEX*16           ACAP(15000) ! Original is 15000 
 
      DIMENSION         COSPHI(181),       SCTPHS(181), 
     2                   RV(31),           RS(31),       
     3                   RN(31),           RC(31),       
     4                   ALGRV(31),        ALGRS(31),    
     5                   ALGRN(31),        DV(31),       
     6                   DS(31),           DN(31),       
     7                   S1(181),           S2(181)
c    8                   LABEL( 9 )
 
      DATA INITIAL / 1 /
 
      IF ( INITIAL .eq. 1 )  then
         IT     = 91
         LL     = 15000 ! Original is 15000
         IT2    = 181
         NSDI   = 31
         NSDIH  = ( NSDI - 1 ) / 2
         NSDI   = 2 * NSDIH + 1
 
         PIE    = ACOS( -1.0 )
         TWOPI  = 2.0 * PIE
         SQRT2  = SQRT( 2.0 )
         SQRT2P = SQRT( TWOPI )
 
         TRUNC  = 4.25
 
C TRUNC IS THE VALUE WHICH DETERMINES THE MAX AND MIN RADIUS IN THE SIZE
C DISTRIBUTION INTEGRATION ACCORDING TO THE FORMULAE
C RMIN = RG*EXP(-TRUNC**2)  AND  RMAX = RG*EXP(TRUNC**2)
 
         nscath = 5
         NSCATA = 2.0 * NSCATH - 1
 
         IF ( IPHASE  .le.  0 )  then
              NSCATH  =  0.0
              NSCATA  =  0.0
         endif
 
C NSCATH IS HALF THE NUMBER OF ANGLES AT WHICH THE PHASE FUNCTION IS TO
C BE COMPUTED
 
         IF ( NSCATA .gt. IT2  .OR.  NSCATH .gt. IT)  then
              WRITE( 6,105 )  NSCATA, NSCATH, IT2, IT
              STOP 11
         endif
 
         DANG = ( PIE / 2.0 ) / FLOAT( NSCATH - 1 )
         do  20  i=1,nscath
              ANG       = ( I-1 ) * DANG
              THETA(I)  = 180.0 * ANG / PIE
              COSPHI(I) = COS( ANG )
              J         = NSCATA - I + 1
              COSPHI(J) = -COSPHI(I)
   20 continue
 
         INITIAL = 0
      endif

c
c add bypass around size distribution calculations -- must redo this
c  calculation if new size distribution is specified
c
      if ( newsd .eq. 1 )  then 
         ALGSIG = ALOG( SIGMAG )
         ALGSSQ = ALGSIG * ALGSIG
         TLGSSQ = 2.0 * ALGSSQ
         ALGRGV = ALOG( RGV )
         ALGRGS = ALGRGV -       ALGSSQ
         ALGRGN = ALGRGV - 3.0 * ALGSSQ
         RGS    = EXP( ALGRGS )
         RGN    = EXP( ALGRGN )
         CN     = 1.0 / ( SQRT2P * ALGSIG )
         CS     = PIE * CN * EXP( 2.0*ALGRGN + 2.0*ALGSSQ )
         CV     = (4.0/3.0) * PIE * CN * EXP( 3.0*ALGRGN + 4.5*ALGSSQ )
         DLOGR  = TRUNC * SQRT2 * ALGSIG / FLOAT( NSDIH )
         DRM    = EXP( -DLOGR )
         DRP    = EXP(  DLOGR )
         K      = NSDIH + 1
         ALGRV(K) = ALGRGV
         ALGRS(K) = ALGRGS
         ALGRN(K) = ALGRGN
         RV(K) = RGV
         RS(K) = RGS
         RN(K) = RGN

         DO 25  K=1,NSDIH
              KM = NSDIH + 1 - K
              KP = NSDIH + 1 + K
              ALGRV(KM) = ALGRGV  -  K * DLOGR
              ALGRS(KM) = ALGRGS  -  K * DLOGR
              ALGRN(KM) = ALGRGN  -  K * DLOGR
              ALGRV(KP) = ALGRGV  +  K * DLOGR
              ALGRS(KP) = ALGRGS  +  K * DLOGR
              ALGRN(KP) = ALGRGN  +  K * DLOGR
              RV(KM) = RV( KM+1 ) * DRM
              RS(KM) = RS( KM+1 ) * DRM
              RN(KM) = RN( KM+1 ) * DRM
              RV(KP) = RV( KP-1 ) * DRP
              RS(KP) = RS( KP-1 ) * DRP
              RN(KP) = RN( KP-1 ) * DRP
   25 continue
 
         VOL = 0.0
         VLC = 0.0
         SFC = 0.0
         XND = 0.0
 
         DO 30  K=1,NSDI
              X     = ALGRV(K)
              DV(K) = CV * EXP( -(X - ALGRGV)**2 / TLGSSQ )
              VOL   = VOL + DV(K)
              X     = ALGRS(K)
              DS(K) = CS * EXP( -(X - ALGRGS)**2 / TLGSSQ )
              SFC   = SFC + DS(K)
              X     = ALGRN(K)
              DN(K) = CN * EXP( -(X - ALGRGN)**2 / TLGSSQ )
              XND   = XND + DN(K)
   30 continue

         DO 35  K=1,NSDI
              RC(K) = 1.0e-12 * rn(k)
              VLC   = VLC + DN(K) * RC(K) ** 3
   35 continue
 
         VOL = VOL * DLOGR
         VLC = ( 4.0 / 3.0 ) * PIE * VLC * DLOGR
         VLC = VLC / VOL
         VLS = 1.0 - VLC
         SFC = SFC * DLOGR
         XND = XND * DLOGR
c        WRITE( 6,120 )  VOL, VLC, VLS, SFC, XND
c        WRITE( 6,125 )  LABEL

c        write( 9,120 )  vol, vlc, vls, sfc, xnd
         newsd  =  0
      endif

 
C REPEAT CALCULATIONS FROM THIS POINT FOR EACH WAVELENGTH
C REPEAT FROM THIS POINT FOR EACH NEW REFRACTIVE INDEX
 
      WNUM = TWOPI / ALAMB
      RFRO = RFRS
      RFIO = RFIS
      RFRI = RFRC
      RFII = RFIC
      EXT  = 0.0
      SCAT = 0.0
      ASY  = 0.0
      BSCT = 0.0
 
      do 40  j=1,nscata
         S1(J) = 0.0
         S2(J) = 0.0
   40 continue

      do 60  k=1,nsdi 
         ROUT = RN(K)
         RIN  = RC(K)
 
         IF ( NSCATH  .eq.  0.0 )  then
              JX  =  1
         ELSE
              JX   = NSCATH
         endif
c         if (k .le. 29 )  go to 60
c          print *, 'll = ', ll

         CALL DMIESS(  ROUT,    RFRO,    RFIO,    THETA,   JX,  
     2                 DQEXT,   DQSCAT,   CTBRQS,  ELTRMX,  PI, 
     3                 TAU,     CSTHT,   SI2THT,  ACAP,    DQBS,  IT,
     4                 LL,      RIN,     RFRI,    RFII,    WNUM   )
 
c         print *, ' '
c         print *, 'Completed k = ', k
c         print *, '    qext, qscat = ', dqext, dqscat

         S1SV = S1(1)
         do 50  j=1,jx
              S1(J) = S1(J) + ELTRMX(1,J,1) * DN(K)
              S2(J) = S2(J) + ELTRMX(2,J,1) * DN(K)
   50    continue

         do 52  j=jx+1,nscata 
              JJ    = NSCATA - J + 1
              S1(J) = S1(J) + ELTRMX(1,JJ,2) * DN(K)
              S2(J) = S2(J) + ELTRMX(2,JJ,2) * DN(K)
   52    continue
 
         X    =  PIE * RN(K) * RN(K) * DN(K)
         EXT  =  EXT  +   DQEXT * X
         SCAT = SCAT  +  DQSCAT * X
         ASY  =  ASY  +  CTBRQS * X
         BSCT = BSCT  +    DQBS * X
   60 continue
 
      ASY   =  ASY / SCAT
      SCAT  = SCAT * DLOGR
      EXT   =  EXT * DLOGR
      BSCT  = BSCT * DLOGR
      OMEGA = SCAT / EXT
      QEXT  =  EXT / SFC
      QSCAT = SCAT / SFC
      QBS   = BSCT / SFC
 
C TEST FOR CONVERGENCE OF PHASE FUNCTION INTEGRATION
      IF ( IPHASE  .gt.  0 )  then
         TEST = ABS( S1(1) - S1SV ) / S1(1)
         IF ( TEST .gt. .01 )
     2         WRITE( 6,150 )  S1SV, S1(1), XLAST
 
         X     = TWOPI / ( WNUM**2  * SCAT )

         do 65  j=1,nscata
              SCTPHS(J) =  X * DLOGR * ( S1(J) + S2(J) )
   65    continue
      endif
 
      RETURN
 
  100 FORMAT ( 7X, I3 )
  105 FORMAT ( ///, 1X,'NUMBER OF ANGLES SPECIFIED =',2I6, / 
     2             10X,'EXCEEDS ARRAY DIMENSIONS =',2I6 )
 
  120 FORMAT (/10X,'INTEGRATED VOLUME',           T40,'=',1PE14.5,/
     2         15X,'PERCENT VOLUME IN CORE',      T40,'=',0PF10.5,/
     3         15X,'PERCENT VOLUME IN SHELL',     T40,'=',0PF10.5,/
     4         10X,'INTEGRATED SURFACE AREA',     T40,'=',1PE14.5,/
     5         10X,'INTEGRATED NUMBER DENSITY',   T40,'=',1PE14.5 )
  125 FORMAT ( 10X,'CORE RADIUS COMPUTED FROM :', /, 20X, 9A8, /  )
 
  150 FORMAT ( ///,1X,'* * * WARNING * * *', / 
     2         10X,'PHASE FUNCTION CALCULATION MAY NOT HAVE CONVERGED'/
     3         10X,'VALUES OF S1 AT NSDI-1 AND NSDI ARE :', 2E14.6, /  
     4         10X,'VALUE OF X AT NSDI =', E14.6 )
 
      END
C====================================================================
      SUBROUTINE DMIESS(  RO,      RFR,     RFI,     THETD,     JX,
     2                    QEXT,    QSCAT,   CTBRQS,  ELTRMX,    PI,
     3                    TAU,     CSTHT,   SI2THT,  ACAP, QBS, IT,
     4                    LL,      R,       RE2,     TMAG2,     WVNO  )
c
      implicit real*8  (a-h,o-z)
C
      COMMON / WARRAY / W(3,30000)
C
      COMPLEX*16    FNAP,      FNBP,      ACAP(LL),     W,
     2              FNA,       FNB,       RF,           RRF,
     3              RRFX,      WM1,       FN1,          FN2,
     4              TC1,       TC2,       WFN(2),       Z(4),
     5              K1,        K2,        K3,
     6              RC,        U(8),         DH1,
     7              DH2,       DH4,       P24H24,       P24H21,
     8              PSTORE,    HSTORE,    DUMMY,        DUMSQ
C
      real*8        T(5),      TA(4),     TB(2),        TC(2),
     2              TD(2),     TE(2),     PI( 3,IT ),   TAU( 3,IT ),
     3              CSTHT(IT), THETD(IT), SI2THT(IT),   ELTRMX( 4,IT,2 )
C
      EQUIVALENCE   (  FNA,TB(1) ),   ( FNB,TC(1) ),
     2              ( FNAP,TD(1) ),   (FNBP,TE(1) )
C
      IFLAG = 1
      IF ( R/RO .LT. 1.0E-06 )   IFLAG = 2
C
      IF ( JX .LE. IT )   GO TO 20
         WRITE( 6,7 )
         WRITE( 6,6 )
         STOP 30
   20 RF =  DCMPLX( RFR,  -RFI )
      RC =  DCMPLX( RE2,-TMAG2 )
      X  =  RO * WVNO
      K1 =  RC * WVNO
      K2 =  RF * WVNO
      K3 =  DCMPLX( WVNO, 0.0d+00 )
      Z(1) =  K2 * RO
      Z(2) =  K3 * RO
      Z(3) =  K1 * R
      Z(4) =  K2 * R
C
      X1   =  DREAL( Z(1) )
      Y1   =  DIMAG( Z(1) )
      X4   =  DREAL( Z(4) )
      Y4   =  DIMAG( Z(4) )
      RRF  =  1.0 / RF
      RX   =  1.0 / X
      RRFX =  RRF * RX
      T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
      T(1) =  DSQRT( T(1) )
      NMX1 =  1.30 * T(1)
C
      IF ( NMX1 .LE. LL-1 )   GO TO 21
         WRITE(6,8)
         STOP 32
   21 NMX2 = 1.1 * T(1)
      IF ( NMX1 .GT.  200 )   GO TO 22
         NMX1 = 200
         NMX2 = 190
C
   22 ACAP( NMX1+1 )  =  ( 0.0,0.0 )
      IF ( IFLAG .EQ. 2 )   GO TO 26
         DO 29   N = 1,3
   29    W( N,NMX1+1 )  =  ( 0.0,0.0 )
   26 CONTINUE
      DO 23   N = 1,NMX1
         NN = NMX1 - N + 1
         ACAP(NN) = (NN+1) * RRFX - 1.0 / ( (NN+1) * RRFX + ACAP(NN+1) )
C
         IF ( IFLAG .EQ. 2 )   GO TO 23
            DO 31   M = 1,3
   31       W( M,NN ) = (NN+1) / Z(M+1)  -
     1                       1.0 / (  (NN+1) / Z(M+1)  +  W( M,NN+1 )  )
   23 CONTINUE
C
      DO 30   J = 1,JX
      IF ( THETD(J) .LT. 0.0 )  THETD(J) =  ABS( THETD(J) )
      IF ( THETD(J) .GT. 0.0 )  GO TO 24
      CSTHT(J)  = 1.0
      SI2THT(J) = 0.0
      GO TO 30
   24 IF ( THETD(J) .GE. 90.0 )  GO TO 25
      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0
      CSTHT(J)  =  DCOS( T(1) )
      SI2THT(J) =  1.0 - CSTHT(J)**2
      GO TO 30
   25 IF ( THETD(J) .GT. 90.0 )  GO TO 28
      CSTHT(J)  =  0.0
      SI2THT(J) =  1.0
      GO TO 30
   28 WRITE( 6,5 )  THETD(J)
      WRITE( 6,6 )
      STOP 34
   30 CONTINUE
C
      DO 35  J = 1,JX
      PI(1,J)  =  0.0
      PI(2,J)  =  1.0
      TAU(1,J) =  0.0
      TAU(2,J) =  CSTHT(J)
   35 CONTINUE
C
C   INITIALIZATION OF HOMOGENEOUS SPHERE
      T(1)   =  DCOS(X)
      T(2)   =  DSIN(X)
      WM1    =  DCMPLX( T(1),-T(2) )
      WFN(1) =  DCMPLX( T(2), T(1) )
      TA(1)  =  T(2)
      TA(2)  =  T(1)
      WFN(2) =  RX * WFN(1) - WM1
      TA(3)  =   DREAL(WFN(2))
      TA(4)  =   DIMAG(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 560
C
      N = 1
C
C INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
C
      SINX1   =  DSIN( X1 )
      SINX4   =  DSIN( X4 )
      COSX1   =  DCOS( X1 )
      COSX4   =  DCOS( X4 )
      EY1     =  DEXP( Y1 )
      E2Y1    =  EY1 * EY1
      EY4     =  DEXP( Y4 )
      EY1MY4  =  DEXP( Y1 - Y4 )
      EY1PY4  =  EY1 * EY4
      EY1MY4  =  DEXP( Y1 - Y4 )
      AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
      BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
      CC  =  SINX1 * ( E2Y1 + 1.0 )
      DD  =  COSX1 * ( E2Y1 - 1.0 )
      DENOM   =  1.0  +  E2Y1 * ( 4.0 * SINX1 * SINX1 - 2.0 + E2Y1 )
      REALP   =  ( AA * CC  +  BB * DD ) / DENOM
      AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
      DUMMY   =  DCMPLX( REALP, AMAGP )
      AA  =  SINX4 * SINX4 - 0.5
      BB  =  COSX4 * SINX4
      P24H24  =  0.5 + DCMPLX( AA,BB ) * EY4 * EY4
      AA  =  SINX1 * SINX4  -  COSX1 * COSX4
      BB  =  SINX1 * COSX4  +  COSX1 * SINX4
      CC  =  SINX1 * SINX4  +  COSX1 * COSX4
      DD  = -SINX1 * COSX4  +  COSX1 * SINX4
      P24H21  =  0.5 * DCMPLX( AA,BB ) * EY1 * EY4  +
     2           0.5 * DCMPLX( CC,DD ) * EY1MY4
      DH4  =  Z(4) / ( 1.0 + ( 0.0,1.0 ) * Z(4) )  -  1.0 / Z(4)
      DH1  =  Z(1) / ( 1.0 + ( 0.0,1.0 ) * Z(1) )  -  1.0 / Z(1)
      DH2  =  Z(2) / ( 1.0 + ( 0.0,1.0 ) * Z(2) )  -  1.0 / Z(2)
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
c 
c temporary print
c
      print 742,  n
      print 744,  u(1), u(2)
      print 745,  u(3), u(4) 
      print 746,  u(5), u(6) 
      print 748,  u(7), u(8)
      print 752,  fna, fnb

c
      GO TO 561
  560 TC1  =  ACAP(1) * RRF  +  RX
      TC2  =  ACAP(1) * RF   +  RX
      FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
      FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
C
  561 CONTINUE
      FNAP = FNA
      FNBP = FNB
      T(1) = 1.50
C
      TB(1) = T(1) * TB(1)
      TB(2) = T(1) * TB(2)
      TC(1) = T(1) * TC(1)
      TC(2) = T(1) * TC(2)
      DO 60 J = 1,JX
          ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
          ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
          ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
          ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
          ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
          ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
          ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
          ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
   60 CONTINUE
C
      QEXT   = 2.0 * ( TB(1) + TC(1))
      QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75
C
      QBSR   = -2.0*(TC(1) - TB(1))
      QBSI   = -2.0*(TC(2) - TB(2))
      RMM    = -1.0
C
      CTBRQS = 0.0
      N = 2
   65 T(1) = 2*N - 1
      T(2) =   N - 1
      T(3) = 2*N + 1
      DO 70  J = 1,JX
          PI(3,J)  = ( T(1) * PI(2,J) * CSTHT(J) - N * PI(1,J) ) / T(2)
          TAU(3,J) = CSTHT(J) * ( PI(3,J) - PI(1,J) )  -
     1                          T(1) * SI2THT(J) * PI(2,J)  +  TAU(1,J)
   70 CONTINUE
C
C   HERE SET UP HOMOGENEOUS SPHERE
      WM1    =  WFN(1)
      WFN(1) =  WFN(2)
      TA(1)  =  DREAL(WFN(1))
      TA(2)  =  DIMAG(WFN(1))
      WFN(2) =  T(1) * RX * WFN(1)  -  WM1
      TA(3)  =  DREAL(WFN(2))
      TA(4)  =  DIMAG(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 1000
C
C   HERE SET UP STRATIFIED SPHERE
C
      DH2  =  - N / Z(2)  +  1.0 / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0 / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0 / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
c
c temporary output to get u functions and a,b coefficients
c
      print 742,  n
      print 744,  u(1), u(2)
      print 745,  u(3), u(4) 
      print 746,  u(5), u(6) 
      print 748,  u(7), u(8)
      print 752,  fna, fnb

  742 format( /, 1x, 'Output for N = ',  i6 )
  744 format(    5x, 'U1,2:   ',  3( 1p2e14.6, 5x) )
  745 format(    5x, 'U3,4:   ',  3( 1p2e14.6, 5x) )
  746 format(    5x, 'U5,6:   ',  3( 1p2e14.6, 5x) )
  748 format(    5x, 'U7,8:   ',  3( 1p2e14.6, 5x) )
  752 format(    5x, 'An, Bn: ',  2( 1p2e14.6, 5x) )

C
 1000 CONTINUE
      TC1  =  ACAP(N) * RRF  +  N * RX
      TC2  =  ACAP(N) * RF   +  N * RX
      FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
      FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
      M    =  WVNO * R
      IF ( N .LT. M )   GO TO 1002
      IF ( IFLAG .EQ. 2 )   GO TO 1001
      IF ( CDABS(  ( FN1-FNA ) / FN1  )  .LT.  1.0E-09   .AND.
     1     CDABS(  ( FN2-FNB ) / FN2  )  .LT . 1.0E-09  )    IFLAG = 2
      IF ( IFLAG .EQ. 1 )   GO TO 1002
 1001 FNA  =  FN1
      FNB  =  FN2
C
 1002 CONTINUE
      T(5)  =  N
      T(4)  =  T(1) / ( T(5) * T(2) )
      T(2)  =  (  T(2) * ( T(5) + 1.0 )  ) / T(5)
C
      CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2)  +
     1                               TE(1) * TC(1)  +  TE(2) * TC(2) )
     2                   +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
      QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
      T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT   =  QSCAT  +  T(3) * T(4)
C
      RMM     =  -RMM
      QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
      QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
C
      T(2)    =  N * (N+1)
      T(1)    =  T(3) / T(2)
      K = (N/2)*2
      DO 80 J = 1,JX
       ELTRMX(1,J,1) = ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1) = ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1) = ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1) = ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      IF ( K .EQ. N )   GO TO 75
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
      GO TO 80
   75  ELTRMX(1,J,2) =ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2) =ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2) =ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2) =ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
   80 continue
C
      IF ( T(4) .LT. 1.0E-10 )   GO TO 100
      N = N + 1
      DO 90 J = 1,JX
         PI(1,J)   =   PI(2,J)
         PI(2,J)   =   PI(3,J)
         TAU(1,J)  =  TAU(2,J)
         TAU(2,J)  =  TAU(3,J)
   90 CONTINUE
      FNAP  =  FNA
      FNBP  =  FNB
      IF ( N .LE. NMX2 )   GO TO 65
         WRITE( 6,8 )

         st1 =  ( X**2 ) * ( RFR**2 + RFI**2 )
         write( 6,888 ) x, rfr, rfi, st1, n, nmx2
         write( 6,889 ) t(4)

         IF ( T(4) .LT. 1.0E-10 )   GO TO 100
         STOP 36
  100 DO 120 J = 1,JX
      DO 120 K = 1,2
         DO  115  I= 1,4
         T(I)  =  ELTRMX(I,J,K)
  115    CONTINUE
         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
  120 CONTINUE
      T(1)    =    2.0 * RX**2
C
C QBS IS THE BACK SCATTER CROSS SECTION
C
      PIG   = ACOS(-1.0)
      RXP4  = RX*RX/(4.0*PIG)
      QBS   = RXP4*(QBSR**2 + QBSI**2)
C
      QEXT    =  QEXT * T(1)
      QSCAT   =  QSCAT * T(1)
      CTBRQS  =  2.0 * CTBRQS * T(1)
C
      RETURN
C
    5 FORMAT( 10X, 'THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN
     1 90.0 DEGREES. IT IS ', E15.4 )
    6 FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
    7 FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
    8 FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST
     1 GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )
c
  888 format( 1h0, 'Output values of x, rfr, rfi, st1, n, nmx2:',  /,
     2         1x,  1p4e12.4, 3x, 2i10  )
  889 format( 1h0, 'Convergence criterion = ', 1pe14.5 )
c
      END



