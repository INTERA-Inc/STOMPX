!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LN_LN( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Linear Isotherm for Aqueous-NAPL Solute Equilibrium
!     Linear Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LN_LN'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Linear Kd sorption  ---
!
      IF( IPCSL(IZN,NSL).EQ.1 ) THEN
        PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!---  Linear, Freundlich, or Langmuir sorption  ---
!
      ELSE
        PCSLX = PCSL(1,IZN,NSL)
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX/(MCLX + PCSLX)
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSLX.LE.EPSL ) THEN
        CLX = CTX/(MCLX + MCNX*PCLN(1,NSL))
        CNX = PCLN(1,NSL)*CLX
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX/(MCLX + MCNX*PCLN(1,NSL) + PCSLX)
        CNX = PCLN(1,NSL)*CLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LN_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FR_LN( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Freundlich Isotherm for Aqueous-NAPL Solute Equilibrium
!     Linear Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FR_LN'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Linear Kd sorption  ---
!
      IF( IPCSL(IZN,NSL).EQ.1 ) THEN
        PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!---  Linear, Freundlich, or Langmuir sorption  ---
!
      ELSE
        PCSLX = PCSL(1,IZN,NSL)
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX/(MCLX + PCSLX)
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSLX.LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL)
        DF = MCLX + (CLX**VAR(4))*MCNX*PCLN(1,NSL)*
     &    PCLN(2,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL) +
     &    CLX*PCSLX
        DF = MCLX + (CLX**VAR(4))*MCNX*PCLN(1,NSL)*
     &    PCLN(2,NSL) + PCSLX
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FR_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LG_LN( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Langmuir Isotherm for Aqueous-NAPL Solute Equilibrium
!     Linear Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE GRID
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LG_LN'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Linear Kd sorption  ---
!
      IF( IPCSL(IZN,NSL).EQ.1 ) THEN
        PCSLX = RHOS(IZN)*PCSL(1,IZN,NSL)*(1.D+0-PORT(2,N))
!
!---  Linear, Freundlich, or Langmuir sorption  ---
!
      ELSE
        PCSLX = PCSL(1,IZN,NSL)
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX/(MCLX + PCSLX)
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSLX.LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1)
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1) +
     &    CLX*PCSLX
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1) + PCSLX
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Linear'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LG_LN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LN_FR( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Linear Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LN_FR'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCSL(2,IZN,NSL))*PCSL(1,IZN,NSL)
        DF = MCLX + (CLX**VAR(3))*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX/(MCLX + MCNX*PCLN(1,NSL))
        CNX = PCLN(1,NSL)*CLX
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        F = -CTX + CLX*(MCLX + MCNX*PCLN(1,NSL)) +
     &    (CLX**PCSL(2,IZN,NSL))*PCSL(1,IZN,NSL)
        DF = MCLX + MCNX*PCLN(1,NSL) + (CLX**VAR(3))*
     &    PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LN_FR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FR_FR( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Freundlich Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FR_FR'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCSL(2,IZN,NSL))*PCSL(1,IZN,NSL)
        DF = MCLX + (CLX**VAR(3))*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = CTX + CLX*MCLX - (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL)
        DF = MCLX - (CLX**VAR(4))*MCNX*PCLN(1,NSL)*PCLN(2,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL) +
     &    (CLX**PCSL(2,IZN,NSL))*PCSL(1,IZN,NSL)
        DF = MCLX + (CLX**VAR(4))*MCNX*PCLN(1,NSL)*
     &    PCLN(2,NSL) + (CLX**VAR(3))*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FR_FR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LG_FR( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Langmuir Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LG_FR'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCSL(2,IZN,NSL))*
     &    PCSL(1,IZN,NSL)
        DF = MCLX + (CLX**VAR(3))* PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1)
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        VAR(3) = -1.D+0 + PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1) +
     &    (CLX**PCSL(2,IZN,NSL))*PCSL(1,IZN,NSL)
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1) + (CLX**VAR(3))*
     *     PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Freundlich'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LG_FR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LN_LG( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Langmuir Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LN_LG'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX - (CLX*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL))/(VAR(2)**2) +
     &     PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX/(MCLX + MCNX*PCLN(1,NSL))
        CNX = PCLN(1,NSL)*CLX
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + CLX*MCNX*PCLN(1,NSL) +
     &    (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX + MCNX*PCLN(1,NSL) - (CLX*PCSL(1,IZN,NSL)*
     &    PCSL(2,IZN,NSL))/(VAR(2)**2) + PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Linear/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LN_LG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FR_LG( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Langmuir Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FR_LG'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX - (CLX*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL))/(VAR(2)**2) +
     &     PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL)
        DF = MCLX + (CLX**VAR(4))*MCNX*PCLN(1,NSL)*
     &    PCLN(2,NSL)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL )  THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        VAR(4) = -1.D+0 + PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX**PCLN(2,NSL))*MCNX*PCLN(1,NSL) +
     &     (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX + (CLX**VAR(4))*MCNX*PCLN(1,NSL)*
     &    PCLN(2,NSL) - (CLX*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL))/
     &    (VAR(2)**2) + PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Freundlich/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*(CLX**PCLN(2,NSL))
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FR_LG group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LG_LG( CTX,CNX,CLX,MCLX,MCNX,N,IZN,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Langmuir Isotherm for Aqueous-NAPL Solute Equilibrium
!     Freundlich Isotherm for Aqueous-Solid Solute Equilibrium
!
!     Compute aqueous and NAPL solute concentrations from total
!     solute concentration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 1, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 VAR(9),MCNX,MCLX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LG_LG'
!
!---  Zero total concentration  ---
!
      IF( CTX/EPSL.LT.EPSL ) THEN
        CNX = 0.D+0
        CLX = 0.D+0
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Zero nonaqueous concentration  ---
!
      IF( MCNX*PCLN(1,NSL).LE.EPSL ) THEN
        CNX = 0.D+0
        CLX = CTX
        NC = 0
   10   CONTINUE
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX - (CLX*PCSL(1,IZN,NSL)*PCSL(2,IZN,NSL))/(VAR(2)**2) +
     &     PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 10
        CLX = MAX( CLX,0.D+0 )
!
!---  Zero sorbed concentration  ---
!
      ELSEIF( PCSL(1,IZN,NSL).LE.EPSL ) THEN
        CLX = CTX
        NC = 0
   20   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1)
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 20
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
!
!---  Solve for concentrations  ---
!
      ELSE
        CLX = CTX
        NC = 0
   30   CONTINUE
        VAR(1) = 1.D+0 + CLX*PCLN(2,NSL)
        VAR(2) = 1.D+0 + CLX*PCSL(2,IZN,NSL)
        F = -CTX + CLX*MCLX + (CLX*MCNX*PCLN(1,NSL))/VAR(1) +
     &     (CLX*PCSL(1,IZN,NSL))/VAR(2)
        DF = MCLX - (CLX*MCNX*PCLN(1,NSL)*PCLN(2,NSL))/(VAR(1)**2) +
     &     (MCNX*PCLN(1,NSL))/VAR(1) - (CLX*PCSL(1,IZN,NSL)*
     &     PCSL(2,IZN,NSL))/(VAR(2)**2) +
     &     PCSL(1,IZN,NSL)/VAR(2)
!
!-----  Newton iteration error, zero partial derivative  ---
!
        IF( ABS(DF).LT.EPSL ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        DCLX = -F/DF
        CLX = CLX + DCLX
        NC = NC+1
!
!-----  Newton iteration error, excessive number of iterations  ---
!
        IF( NC.GT.32 ) THEN
          INDX = 17
          N_DB = N
          CHMSG = 'Solute Distribution Convergence Failure: '//
     &      'Equilibrium (Total Solute): '//
     &      'L:N Langmuir/S:L Langmuir'
          CALL WRMSGS(INDX)
        ENDIF
        IF( ABS(DCLX/(CLX+EPSL)).GT.1.D-6 .AND. CLX.GT.0.D+0 ) GOTO 30
        CLX = MAX( CLX,0.D+0 )
        CNX = PCLN(1,NSL)*CLX/(1.D+0+PCLN(2,NSL)*CLX)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LG_LG group  ---
!
      RETURN
      END

