C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C 'dlight.f'
C Copyright (C) 2017 Jessica Kirkby-Kent & Pierre Maxted
C Contact: j.kirkbykent@gmail.com
C Ver. 1.000 (For python 2.7)
C --Part of pyebop.py--
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      SUBROUTINE DLIGHT (V,PH,NDAT,MAG)
C Double precision version, p.maxted@keele.ac.uk, 26 Jan 2015
C
C  Must be compiled before it will run. To compile use
C
C      f2py -c -m dlight dlight.f
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
cf2py intent(in) v, ph
cf2py intent(out) mag
cf2py integer intent(hide), depend(ph) ::  ndat = len(ph)
      integer NDAT
      double precision V(18),MAG(NDAT)
      double precision PH(NDAT),TP,TS
      double precision LP,LS,LECL,LE
C      COMMON /VARIB/ BS,RP,RATIO,UP,US,FI,ECOSW,ESINW,YP,YS,SP,SS,Q,
C     $       TANGL,EL,DPH,SFACT,DGAM,EXTRA(13)
C     REF.  THE ASTROPHYSICAL JOURNAL, 174 617-628,1972 JUNE 15
C     ECLIPSING-BINARY SOLUTIONS BY SEQUENTIAL OPTIMIZATION
C     OF THE PARAMETERS
C
C
      DATA PI,TWOPI,RAD/3.141592653589793D0,6.28318530717959D0,
     +                  0.01745329251994D0/ 
C
C        DETERMINE PRIMARY AND SECONDARY BIAXIAL DIMENSIONS
C        USE SPERICAL RADII FOR THE COMPUTATION OF ECLIPSE FUNCTIONS
C        USE OBLATENESSES FOR THE COMPUTATION OF THE OUTSIDE ECLIPSE
C        PHOTOMETRIC VARIATIONS WITH LIMB AND GRAVITY DARKENING
C

      BS     = V(1)
      RP     = V(2)/(1.0D0+V(3))
      RATIO  = V( 3)
      UP     = V( 4)
      US     = V( 5)
      FI     = V( 6)
      ECOSW  = real(V( 7))
      ESINW  = real(V( 8))

      YP     = real(V( 9))
      YS     = real(V(10))
      SP     = real(V(11))
      SS     = real(V(12))
      Q      = real(V(13))
      TANGL  = real(V(14))
      EL     = real(V(15))
      DPH    = V(16)
      SFACT  = V(17)
      DGAM   = V(18)
 
      if ( Q <= 0.0D0 ) then
        CALL BIAX (RP,0.0D0,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,0.0D0,RSA,RSB,ES)
      else
        CALL BIAX (RP,Q,RPA,RPB,EP)
        RS=RP*RATIO
        CALL BIAX (RS,1.0D0/Q,RSA,RSB,ES)
      end if
c      
      do i=1,ndat
      phase=ph(i)    
C
C        CORRECT THE OBSERVED PHASE FOR ANY EPOCH ERROR IN EPHEMERIS
C
      THETA=PHASE+DPH
C
      SINI  = SIN(FI*RAD)
      SINI2 = SINI*SINI
      COSI2 = 1.0D0  - SINI2
C
C        TRANSLATE TIDAL LEAD/LAG ANGLE TO RADIANS
      TANGR=TANGL*RAD
C
C     EQUATION 9
C        CONVERT PHASE TO RADIANS
      FMN=THETA*TWOPI
C
C        GET CURRENT VALUES OF E, AND W
      CALL GETEW (ECOSW,ESINW,E,W)
C
C        TEST FOR CIRCULAR ORBIT
      IF (E)   17,20,17
   20 COSVW=COS(FMN)
      SINVW=SIN(FMN)
      RV=1.0D0
      GO TO 25
C
C        SOLUTION OF KEPLER'S EQUATION BY DIFFERENTIAL CORRECTIONS
C        (NON-ZERO ECCENTRICITY ONLY . . . )
C
C     EQUATION 6
C
   17 OMEGA = 450.0D0  - W
   23 IF (OMEGA - 360.0D0)         22,21,21
   21 OMEGA = OMEGA - 360.0D0
      GO TO 23
   22 OMEGA = OMEGA*RAD
C        SINE AND COSINE OF OMEGA
      COSW=COS(OMEGA)
      SINW=SIN(OMEGA)
C
C        COMPUTE MEAN ANOMALY CORRECTION TO PHASE
C        CORRESPONDING TO V=OMEGA=90-W
C        AT WHICH PHASE COS(V-OMEGA)=1
      E0=ATAN2(SQRT(1.0D0-E*E)*SINW,COSW+E)
C
C        MEAN ANOMALY OF MID-PRIMARY ECLIPSE
      FMA0=E0-E*SIN(E0)
C
C        MEAN ANOMALY
      FMA=FMN+FMA0
C     FIRST APPROXIMATION OF ECCENTRIC ANOMALY
      EA=FMA+E*SIN(FMA)
C
      DO 10 J=1,15
C        EVALUATE SINE AND COSINE OF ECCENTRIC ANOMALY
      SINE=SIN(EA)
      COSE=COS(EA)
      DENOM=1.0D0-E*COSE
      DISC=FMA-EA+E*SINE
      EA=EA+DISC/DENOM
C        TEST FOR CONVERGENCE
      IF (ABS(DISC) - 2.0D-05)     15,15,10
   10 CONTINUE
C
C
C        EVALUATE SINE AND COSINE OF TRUE ANOMALY
   15 COSV=(COSE-E)/DENOM
      SINV=SINE*SQRT(1.0D0-E*E)/DENOM
C
C        RADIUS VECTOR
      RV = (1.0E0-E*E)/(1.0D0+E*COSV)
C
C        THE PHOTOMETRIC PHASE ARGUMENT IN TERMS OF ORBIT PARAMETERS
C        VW = V-OMEGA
      COSVW=COSV*COSW+SINV*SINW
      SINVW=SINV*COSW-COSV*SINW
C
   25 COS2=COSVW*COSVW
      SIN2=1.0D0-COS2
C
      CSVWT=COS(TANGR)*COSVW-SIN(TANGR)*SINVW
C
C
C        PHOTOMETRIC EFFECTS
C
C
C        TEST FOR SIMPLE CASE OF TWO SPHERICAL STARS
      IF (EP .EQ. 0.0D0  .AND.  ES .EQ. 0.0D0) GO TO 26
C
C        EITHER OR BOTH STARS ARE OBLATE
C
      FMAXP=((1.0D0-UP)+0.666666667D0*UP*(1.0D0+0.2D0*EP))
     1      *(1.0D0+3.0D0*YP*EP)/(1.0D0-EP)
      FMAXS=((1.0D0-US)+0.666666667D0*US*(1.0D0+0.2D0*ES))
     1      *(1.0D0+3.0D0*YS*ES)/(1.0D0-ES)
C        CHANGE IN INTENSITY RATIO DUE TO OBLATENESS RELATED VARIABLES
C        FROM QUADRATURE TO MINIMUM
C        FACE ON TO END ON
      DELTP=(15.0D0+UP)/(15.0D0-5.0D0*UP)*(1.0D0+YP)*EP
      DELTS=(15.0D0+US)/(15.0D0-5.0D0*US)*(1.0D0+YS)*ES
C        FORE-SHORTENING FUNCTION OF OBLATENESS
      SHORT=SINI2*CSVWT*CSVWT
      GO TO 27
C
C        BOTH STARS ARE SPHERICAL
C
   26 FMAXP=1.0D0-UP/3.0D0
      FMAXS=1.0D0-US/3.0D0
      DELTP=0.0D0
      DELTS=0.0D0
      SHORT=0.0
C
C        UN-NORMALIZED BRIGHTNESS OF STELLAR COMPONENTS AT QUADRATURE
   27 OP=PI*RPB*RPB*FMAXP
      OS=PI*RSB*RSB*FMAXS*BS
C        THE NORMALIZING FACTOR
      OTOT=OP+OS
C        BRIGHTNESS CONTRIBUTION FROM EACH COMPONENT
      LP=OP/OTOT*(1.0D0-DELTP*SHORT)
      LS=OS/OTOT*(1.0D0-DELTS*SHORT)
C
C        REFLECTION AND RERADIATION EQUATION
      IF (SP .EQ. 0.0D0  .AND.  SS .EQ. 0.0D0)   GO TO 28
      HEAT=SINI*COSVW
      HEAT2=0.5D0+0.5D0*HEAT*HEAT
      DLP=SP*(HEAT2+HEAT)
      DLS=SS*(HEAT2-HEAT)
      GO TO 29
   28 DLP=0.0D0
      DLS=0.0D0
C
C        WHICH ECLIPSE COULD THIS BE
   29 IF (COSVW)         40,40,30
C
C     PRIMARY ECLIPSE
C
   30 R1 = RP
      R2 = RS
      UU = UP
      LE=LP
      DLE=DLP
      TP = 1.0D0
      TS = 0.0D0
      GO TO 60
C
C
C     SECONDARY ECLIPSE
C
   40 R1 = RS
      R2 = RP
      UU = US
      LE=LS
      TP = 0.0D0
      TS = 1.0D0
      DLE=DLS
C
   60 SUM = 0.0D0
      ALAST = 0.0D0
      AREA=0.0D0
C
C     EQUATION  5
C
      DD = SINVW*SINVW + COSVW*COSVW*COSI2
      IF (DD .LE. 1.0D-06)  DD=0.0D0
      DD = DD*RV*RV
      D = SQRT(ABS(DD))
      R22 = R2*R2
C
C     EQUATION 17
C
      GAMN = 90.01D0*RAD
      DGAMA = DGAM*RAD
      DGM = DGAMA/2.0D0
      RK = 0.0D0
      GAM = 0.0D0
   50 GAM = GAM + DGAMA
C        HAS LIMIT OF INTEGRATION BEEN REACHED
      IF (GAM - GAMN)              48,48,49
C
   48 RR = R1*SIN(GAM)
      R12 = RR*RR
C
      AA = 0.0D0
C        ARE THE PROJECTED DISKS CONCENTRIC
      IF (D)                       405,406,405
  406 IF (RR - R2)                 230,230,403
  403 IF (RK - R2)                 404, 49, 49
  404 AA = PI*R22
      GO TO 215
C        TEST FOR NO ECLIPSE
  405 IF (D-R1-R2)                 240,216,216
  216 SUM = 0.0D0
      GO TO 49
C        DECIDE WHICH AREA EQUATIONS FOR NON-CONCENTRIC ECLIPSE
  240 IF (D-RR-R2)                 245,215,215
  245 IF (D-R2+RR)                 230,230,250
  250 IF (R1-R2)                   255,255,280
  255 IF (DD-R22+R12)              205,210,210
  280 IF (D-RR+R2)                 290,260,260
  260 IF (RR-R2)                   255,255,265
  265 IF (DD-R12+R22)              270,210,210
C
C     EQUATION 12
C
  270 S1 = ABS((R12 - R22 - DD)*0.5D0/D)
      A1 = ABS(R2-S1)
      B2 = ABS(RR-S1-D  )
      AA=PI*R22-(R22*ACOS((R2-A1)/R2)
     1   - (R2-A1)*SQRT(2.0D0*R2*A1-A1*A1))
     2   +R12*ACOS((RR-B2)/RR)-(RR-B2)*SQRT(2.0D0*RR*B2-B2*B2)
      GO TO 215
C
  290 IF (R1 - R2 - D)             260,260,295
  295 IF (RK - R2 - D)             300,215,215
  300 RR = R2 + D
      R12 = RR*RR
      GAMN = 0.0D0
      GO TO 260
C
  230 AA = PI*R12
      GO TO 215
C
C     EQUATION 10
C
  205 S = ABS((R12 - R22 + DD)*0.5D0/D)
      A = ABS(RR-S)
      B1 = ABS(R2-S-D)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0D0*RR*A - A*A)
      AB1 = R22*ACOS((R2-B1)/R2) - (R2-B1)*SQRT(2.0D0*R2*B1-B1*B1)
      AA = PI*R12 - A1 + AB1
      GO TO 215
C
C     EQUATION 1
C
  210 S = ABS((R12 - R22 + DD)*0.5D0/D)
      A = ABS(RR-S)
      B = ABS(S-D+R2)
      A1 = R12*ACOS((RR-A)/RR) - (RR-A)*SQRT(2.0D0*RR*A - A*A)
      AA1 = R22*ACOS((R2-B)/R2) - (R2-B)*SQRT(2.0D0*R2*B - B*B)
      AA = A1 + AA1
C
  215 DAREA = AA - ALAST
      SUM = SUM + DAREA*(1.0D0  - UU + UU*COS(GAM-DGM))
      ALAST = AA
      AREA = AREA + DAREA
C
      RK = RR
      GO TO 50
C
C        LIGHT LOSS FROM ECLIPSE
C
   49 ADISK = PI*R1*R1
      ALPHA = SUM/(ADISK*(1.0D0-UU/3.0D0))
      LECL = ALPHA*LE
      AREA = AREA/ADISK
      REFL=DLP+DLS-AREA*DLE
C
C        THEORETICAL INTENSITY WITH THIRD LIGHT AND QUADRATURE
C        SCALE FACTOR APPLIED
C
      FLITE = ((LP+LS-LECL+REFL)*(1.0D0-EL)+EL)

      MAG(i)=SFACT-2.5D0*DLOG10(FLITE)
C
      ENDDO
      RETURN
      END
      
      SUBROUTINE BIAX (R,Q,A,B,EPS)
            ! EBOP subroutine to calculate biaxial ellipsoid dimensions
            ! and oblateness for each star after Chandrasekhar (1933).
      DOUBLE PRECISION R,Q,A,B,EPS
      if ( Q <= 0.0D0 )  then
        A = R
        B = R
        EPS = 0.0D0
      else
        A = R * ( 1.0D0 + (1.0D0 + 7.0D0*Q)/6.0D0 * R**3)
        B = R * ( 1.0D0 + (1.0D0 - 2.0D0*Q)/6.0D0 * R**3)
        EPS = (A - B) / A
        B=( (1.0D0 - EPS) * R**3) ** (1.0D0/3.0D0)
        A = B / (1.0D0 - EPS)
      end if

      END SUBROUTINE BIAX      


      SUBROUTINE GETEW (ECOSW,ESINW,E,W)
            ! EBOP subroutine to calculate e and w from e(cos)w e(sin)w
      DOUBLE PRECISION ECOSW,ESINW,E,W
      if ( ECOSW == 0.0D0  .and.  ESINW == 0.0D0 ) then
        E = 0.0D0
        W = 0.0D0
      else
        W = atan2( ESINW,ECOSW )
        E = sqrt( ESINW*ESINW + ECOSW*ECOSW )
        W = W * 180.0D0 / 3.141592653589793D0
      end if

      END SUBROUTINE GETEW
