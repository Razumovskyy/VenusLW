! *** ONLY for 1 term *** 18 August, 2022 (8.11.2002).
SUBROUTINE Planck_separate_H(V1,V2,N_K)
        IMPLICIT REAL*8 (A-H,O-Z)
CHARACTER KISHKANAME*15,N_K*3 
      PARAMETER (N_T=8001,Tmin=100.d0,Tmax=900.d0,N_SEP=20,H=2.D0/2048.D0)
        DIMENSION T(N_T),PLA(N_T),PLA_INT_SEP(N_SEP,N_T),PLK(N_T)
	COMMON/KSHK/KISHKANAME
PLANK2(V1,ttt1)=3.7403E-8*V1**3/(DEXP(1.43868*V1/ttt1)-1.)*2.
	    DATA WPL_OLD/-100./
	    VSTART=V1
	    VFINISH=V2
             DT=(Tmax-Tmin)/(N_T-1)
              DO L =1,N_T
              T(L)=Tmin+DT*(L-1)
              DO M =1,N_SEP
              PLA_INT_SEP(M,L)=0.D0
              END DO 
              END DO
         OPEN(10,FILE=KISHKANAME) 
            READ(10,*)WO
            IF(WO.LT.0.0)WO=-WO
              DO I=1,90000000
              READ(10,*,END=1)WN,KAN
               IF(WN.GT.0.0D0) THEN
                  W=(WO+WN)/2.D0
!                  H=(WN-WO) 
                  WO=WN
        IF(W.GT.VSTART.AND.W.LT.VFINISH)THEN ! ***********
           IF(DABS(WPL_OLD-W).ge.0.1D0)THEN
           WPL_OLD = W + 0.05D0
                  DO L=1,N_T
                  TT=T(L)
                   PLK(L)=PLANK2(WPL_OLD,TT)*H ! **************
                  END DO
            END IF
 
                  DO L=1,N_T
	PLA_INT_SEP(KAN,L)=PLA_INT_SEP(KAN,L)+PLK(L)
                  END DO
        END IF          
                               ELSE
                WO=-WN
                               END IF 

                END DO
  1       CONTINUE
! ***********************************************************
!    *** Only 1-term *** !
M=2 ! 1-term ONLY.
OPEN(11,FILE='PLANCKseparate.'//N_K) 
	WRITE(11,*)N_T
	 DO L=1,N_T
       WRITE(11,*)T(L),PLA_INT_SEP(M,L)
	 END DO
	CLOSE(11)
      END 
       
