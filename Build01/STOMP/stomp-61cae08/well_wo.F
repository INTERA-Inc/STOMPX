!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCFW4
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
!     Water-Oil Mode
!     Well boundary condition fluxes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2001.
!     Last Modified by M.D. White, PNNL, 5 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE JACOB
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
      INTEGER MP(5),MB(5)
      SAVE MP,MB
!
!----------------------Data Statements---------------------------------!
!
      DATA MP /2,2,3,2,4/
      DATA MB /2,3,2,4,2/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCFW4'
!
!---  Loop over wells  ---
!
      DO 1000 NWL = 1,NWLS
        TMZ = TM
        IF( IWM(NWL).EQ.0 ) GOTO 1000
        IF( IWCC(NWL).EQ.1 ) TMZ = MOD( TM,WLVR(1,IWM(NWL),NWL) )
        IF( TMZ.LE.WLVR(1,1,NWL) ) GOTO 1000
        IF( IWM(NWL).GT.1 .AND. TMZ.GT.WLVR(1,IWM(NWL),NWL) ) GOTO 1000
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        K = IWLDM(4,NWL)
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NQZ = NSZW(NW)+1
!
!---    Zero well boundary fluxes  ---
!
        DO 10 M = 1,ISVF
          WL_W(M,NQZ) = 0.D+0
          WN_W(M,NQZ) = 0.D+0
          WDLO_W(M,NQZ) = 0.D+0
   10   CONTINUE
!
!---  Well permeabilities, equivalent to
!     a saturated hydraulic conductivity of 100 m/s  ---
!
        PRNWX = 1.D-5/VISRN
        PRLWX = 1.D-5/VISRL
!
!---    Monitoring well  ---
!
        IF( IWT(NWL).EQ.1 ) THEN
          DO 100 M = 1,ISVF
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.5 ) THEN
!
!---          Advective aqueous-phase flux out well head  ---
!
              HDLZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHOLW(MP(M),NW)
              HDLZ = MAX( HDLZ,0.D+0 )
              RKLWX = SLW(MP(M),NW)
              WL_W(M,NQZ) = 2.D+0*PRLWX*RKLWX*HDLZ/DZGF(N)
!
!---          Advective NAPL flux out well head  ---
!
              HDNZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHONW(MP(M),NW)
              HDNZ = MAX( HDNZ,0.D+0 )
              RKNWX = SNW(MP(M),NW)
              WN_W(M,NQZ) = 2.D+0*PRNWX*RKNWX*HDNZ/DZGF(N)
            ENDIF
  100     CONTINUE
!
!---    Aqueous slurping well  ---
!
        ELSEIF( IWT(NWL).EQ.2 ) THEN
          DO 200 M = 1,ISVF
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.5 ) THEN
!
!---          Advective aqueous-phase flux out well head  ---
!
              HDLZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHOLW(MP(M),NW)
              HDLZ = MAX( HDLZ,0.D+0 )
              RKLWX = SLW(MP(M),NW)
              WL_W(M,NQZ) = 2.D+0*PRLWX*RKLWX*HDLZ/DZGF(N)
!
!---          Advective NAPL flux out well head  ---
!
              HDNZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHONW(MP(M),NW)
              HDNZ = MAX( HDNZ,0.D+0 )
              RKNWX = SNW(MP(M),NW)
              WN_W(M,NQZ) = 2.D+0*PRNWX*RKNWX*HDNZ/DZGF(N)
            ENDIF
  200     CONTINUE
!
!---    NAPL slurping well  ---
!
        ELSEIF( IWT(NWL).EQ.3 ) THEN
          DO 300 M = 1,ISVF
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.5 ) THEN
!
!---          Advective aqueous-phase flux out well head  ---
!
              HDLZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHOLW(MP(M),NW)
              HDLZ = MAX( HDLZ,0.D+0 )
              RKLWX = SLW(MP(M),NW)
              WL_W(M,NQZ) = 2.D+0*PRLWX*RKLWX*HDLZ/DZGF(N)
!
!---          Advective NAPL flux out well head  ---
!
              HDNZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHONW(MP(M),NW)
              HDNZ = MAX( HDNZ,0.D+0 )
              RKNWX = SNW(MP(M),NW)
              WN_W(M,NQZ) = 2.D+0*PRNWX*RKNWX*HDNZ/DZGF(N)
            ENDIF
  300     CONTINUE
!
!---    Multifluid pumping well  ---
!
        ELSEIF( IWT(NWL).EQ.4 ) THEN
          DO 400 M = 1,ISVF
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.5 ) THEN
!
!---          Advective aqueous-phase flux out well head  ---
!
              HDLZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHOLW(MP(M),NW)
              HDLZ = MAX( HDLZ,0.D+0 )
              RKLWX = SLW(MP(M),NW)
              WL_W(M,NQZ) = 2.D+0*PRLWX*RKLWX*HDLZ/DZGF(N)
!
!---          Advective NAPL flux out well head  ---
!
              HDNZ = PW(MP(M),NW)-PWB(MB(M),NWL)-5.D-1*GRAVZ*DZGF(N)*
     &          RHONW(MP(M),NW)
              HDNZ = MAX( HDNZ,0.D+0 )
              RKNWX = SNW(MP(M),NW)
              WN_W(M,NQZ) = 2.D+0*PRNWX*RKNWX*HDNZ/DZGF(N)
            ENDIF
  400     CONTINUE
        ENDIF
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCFW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCJW4
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
!     Water-Oil Mode
!     Modify the Jacobian matrix for the well boundary conditions.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 FW(LUK+1),FO(LUK+1)
      INTEGER MP(3),MF(3)
      INTEGER MJ(2)
      SAVE MP,MF,MJ
!
!----------------------Data Statements---------------------------------!
!
      DATA MP /2,3,4/
      DATA MF /1,3,5/
      DATA MJ /2,3/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCJW4'
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
!
!---  Loop over wells  ---
!
      DO 1000 NWL = 1,NWLS
        TMZ = TM
        IF( IWM(NWL).EQ.0 ) GOTO 1000
        IF( IWCC(NWL).EQ.1 ) TMZ = MOD( TM,WLVR(1,IWM(NWL),NWL) )
        IF( TMZ.LE.WLVR(1,1,NWL) ) GOTO 1000
        IF( IWM(NWL).GT.1 .AND. TMZ.GT.WLVR(1,IWM(NWL),NWL) ) GOTO 1000
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        K = IWLDM(4,NWL)
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NQZ = NSZW(NW)+1
        RDW = WBR(NWL)
        ARWX = GPI*(RDW**2)
!
!---    Well flow (top flux)  ---
!
        DO 100 M = 1,ISVC+1
!
!---      Water equation  ---
!
          FLW = XLWW(MP(M),NW)*RHOLW(MP(M),NW)
          FW(M) = ARWX*WL_W(MF(M),NQZ)*FLW
!
!---      Oil equation  ---
!
          FLW = XLOW(MP(M),NW)*RHOLW(MP(M),NW)
          FNW = RHONW(MP(M),NW)
          FO(M) = ARWX*(WL_W(MF(M),NQZ)*FLW + WN_W(MF(M),NQZ)*FNW)
  100   CONTINUE
!
!---    Load Jacobian matrix  ---
!
        NMD = IXP(N)
        MPW = IM(IEQWW,NMD)
        MPA = IM(IEQOW,NMD)
        DO 200 M = 1,ISVC
!
!---      Banded solver  ---
!
          IF( ILES.EQ.1 ) THEN
            MCOL = IM(M+ISVC,NMD)
!
!---        Water equation  ---
!
            MROW = MPW-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(FW(MJ(M))-FW(1))/DNRW(M,NW)
!
!---        Oil equation  ---
!
            MROW = MPA-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(FO(MJ(M))-FO(1))/DNRW(M,NW)
!
!---      Conjugate gradient solver  ---
!
          ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---        Water equation  ---
!
            MCOL = KLU(MPW,M+5*ISVC)
            DLU(MCOL) = DLU(MCOL)+(FW(MJ(M))-FW(1))/DNRW(M,NW)
!
!---        Oil equation  ---
!
            MCOL = KLU(MPA,M+5*ISVC)
            DLU(MCOL) = DLU(MCOL)+(FO(MJ(M))-FO(1))/DNRW(M,NW)
          ELSE
            WRITE(ISC,'(A)') 'ERROR: Unknown Linear Equation Solver'
            WRITE(IWR,'(A)') 'ERROR: Unknown Linear Equation Solver'
            STOP
          ENDIF
  200   CONTINUE
!
!---    Water equation  ---
!
        BLU(MPW) = BLU(MPW) - FW(1)
        RSDL(IEQWW,N) = BLU(MPW)
!
!---    Oil equation  ---
!
        BLU(MPA) = BLU(MPA) - FO(1)
        RSDL(IEQOW,N) = BLU(MPA)
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCJW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BCPW4
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
!     Water-Oil Mode
!     Well boundary condition properties.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2001.
!     Last Modified by M.D. White, PNNL, 5 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE JACOB
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
      REAL*8 WLVRX(LNWV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BCPW4'
!
!---  Loop over wells  ---
!
      DO 1000 NWL = 1,NWLS
        TMZ = TM
        IF( IWM(NWL).EQ.0 ) GOTO 1000
        IF( IWCC(NWL).EQ.1 ) TMZ = MOD( TM,WLVR(1,IWM(NWL),NWL) )
        IF( TMZ.LE.WLVR(1,1,NWL) ) GOTO 1000
        IF( IWM(NWL).GT.1 .AND. TMZ.GT.WLVR(1,IWM(NWL),NWL) ) GOTO 1000
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        K = IWLDM(4,NWL)
        N = ND(I,J,K)
        NW = ABS(IXW(N))
!
!---  Assign local well boundary condition variables  ---
!
        IF( IWM(NWL).EQ.1 ) THEN
          DO 10 M = 1,LNWV
            WLVRX(M) = WLVR(M,1,NWL)
   10     CONTINUE
        ELSE
          DO 30 L = 2,IWM(NWL)
            IF( TMZ.LE.WLVR(1,L,NWL) ) THEN
              TDBC = (WLVR(1,L,NWL)-WLVR(1,L-1,NWL))
              TFBC = (TMZ-WLVR(1,L-1,NWL))/TDBC
              DO 20 M = 1,LNWV
                WLVRX(M) = WLVR(M,L-1,NWL) +
     &           TFBC*(WLVR(M,L,NWL)-WLVR(M,L-1,NWL))
   20         CONTINUE
              GOTO 40
            ENDIF
   30     CONTINUE
          GOTO 1000
        ENDIF
   40   CONTINUE
!
!---    Assign values to well boundary properties  ---
!
        TX = T(2,N)
        PWX = WLVRX(2)
!
!---    Monitoring well  ---
!
        IF( IWT(NWL).EQ.1 ) THEN
          DO 100 M = 2,ISVC+2
            PWB(M,NWL) = PWX - PATM
  100     CONTINUE
!
!---    Aqueous slurping well  ---
!
        ELSEIF( IWT(NWL).EQ.2 ) THEN
          DO 200 M = 2,ISVC+2
            PWB(M,NWL) = PWX - PATM
  200     CONTINUE
!
!---    NAPL slurping well  ---
!
        ELSEIF( IWT(NWL).EQ.3 ) THEN
          DO 300 M = 2,ISVC+2
            PWB(M,NWL) = PWX - PATM
  300     CONTINUE
!
!---    Multifluid pumping well  ---
!
        ELSEIF( IWT(NWL).EQ.4 ) THEN
          DO 400 M = 2,ISVC+2
            PWB(M,NWL) = PWX - PATM
  400     CONTINUE
        ENDIF
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BCPW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE CHKW4
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
!     Water-Oil Mode
!     Initialize well state variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2001.
!     Last Modified by M.D. White, PNNL, 5 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE GRID
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/CHKW4'
!
!---  Restart simulation  ---
!
      IF( IEO.EQ.2 ) THEN
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Initialize well primary variables  ---
!
      DO 20 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        PWB(2,NWL) = WHP(NWL)-PATM
        PWX = WHP(NWL)
!
!---    Assume equilibrium conditions, determine whether NAPL is
!       denser or lighter than water, using the well-bottom
!       temperature and well-head pressure  ---
!
        N = ND(I,J,KSX)
        CALL WATLQD( T(2,N),PWX,RHOLX )
        CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONX )
!
!---  Loop over well nodes from the top down  ---
!
        DO 10 K = KEX,KSX,-1
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          DZTX = ZE(5,N)-ZE(1,ND(I,J,KSX))
          DZPX = ZP(N)-ZE(1,ND(I,J,KSX))
          DZBX = ZE(1,N)-ZE(1,ND(I,J,KSX))
!
!---      LNAPL  ---
!
          IF( RHOLX.GT.RHONX ) THEN
            DNTWX = WWD(NWL)
            DNBWX = WWD(NWL)-WWDN(NWL)
            DWTWX = WWD(NWL)-WWDN(NWL)
            DWBWX = 0.D+0
            DIWX =  WWD(NWL)-WWDN(NWL)
!
!---      DNAPL  ---
!
          ELSEIF( RHOLX.LT.RHONX ) THEN
            DNTWX = WWDN(NWL)
            DNBWX = 0.D+0
            DWTWX = WWD(NWL)
            DWBWX = WWDN(NWL)
            DIWX = WWDN(NWL)
!
!---      Equal densities  ---
!
          ELSE
            DNTWX = WWD(NWL)
            DNBWX = WWD(NWL)
            DWTWX = WWD(NWL)
            DWBWX = WWD(NWL)
            IF( WWDN(NWL).LT.EPSL ) DWBWX = 0.D+0
            IF( WWD(NWL)-WWDN(NWL).LT.EPSL ) DNBWX = 0.D+0
            DIWX = DZPX
          ENDIF

!
!---      Total-liquid saturated well node  ---
!
          IF( WWD(NWL).GE.DZTX ) THEN
!
!---        Aqueous only well node  ---
!
            IF( DWTWX.GE.DZTX .AND. DWBWX.LE.DZBX ) THEN
              NPHAZW(2,NW) = 1
              CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
              CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
              YMLOW(2,NW) = WIDO(NWL)*PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
              STW(2,NW) = 1.D+0
              IF( YMLOW(2,NW).GT.XMLOMX ) THEN
                YLOWX = (YMLOW(2,NW)*WTMO)/
     &            (YMLOW(2,NW)*WTMO+(1.D+0-YMLOW(2,NW))*WTMW)
                XLOMX = (XMLOMX*WTMO)/(XMLOMX*WTMO+(1.D+0-XMLOMX)*WTMW)
                SNW(2,NW) = STW(2,NW)*RHOLW(2,NW)*(YLOWX-XLOMX)/
     &            (RHONW(2,NW)+RHOLW(2,NW)*(YLOWX-XLOMX))
              ELSE
                SNX = 0.D+0
              ENDIF
              SLW(2,NW) = STW(2,NW)-SNW(2,NW)
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
!
!---        NAPL only well node  ---
!
            ELSEIF( DNTWX.GE.DZTX .AND. DNBWX.LE.DZBX ) THEN
              NPHAZW(2,NW) = 5
              SLW(2,NW) = 0.D+0
              SNW(2,NW) = 1.D+0
              STW(2,NW) = 1.D+0
              CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
              CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
              YMLOW(2,NW) = PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
!
!---        Aqueous and NAPL well node  ---
!
            ELSEIF( DIWX.LT.DZTX .AND. DIWX.GT.DZBX ) THEN
              NPHAZW(2,NW) = 5
              STW(2,NW) = 1.D+0
              IF( RHOLX.GT.RHONX ) THEN
                SLW(2,NW) = (DIWX-DZBX)/(DZTX-DZBX)
                SNW(2,NW) = MAX( 1.D+0-SLW(2,NW),0.D+0 )
              ELSEIF( RHOLX.LT.RHONX ) THEN
                SNW(2,NW) = (DIWX-DZBX)/(DZTX-DZBX)
                SLW(2,NW) = MAX( 1.D+0-SNW(2,NW),0.D+0 )
              ELSE
                SNW(2,NW) = WWDN(NWL)/WWD(NWL)
                SLW(2,NW) = MAX( 1.D+0-SNW(2,NW),0.D+0 )
              ENDIF
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
              YMLOW(2,NW) = PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
            ENDIF
!
!---  Total-liquid partially saturated well node  ---
!
          ELSEIF( WWD(NWL).LT.DZTX .AND. WWD(NWL).GE.DZBX ) THEN
!
!---        Aqueous only well node  ---
!
            IF( RHOLX.LT.RHONX .AND. DWBWX.LE.DZBX ) THEN
              NPHAZW(2,NW) = 2
              STW(2,NW) = (WWD(NWL)-DZBX)/(DZTX-DZBX)
              CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
              CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
              YMLOW(2,NW) = WIDO(NWL)*PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
              IF( YMLOW(2,NW).GT.XMLOMX ) THEN
                YLOWX = (YMLOW(2,NW)*WTMO)/
     &            (YMLOW(2,NW)*WTMO+(1.D+0-YMLOW(2,NW))*WTMW)
                XLOMX = (XMLOMX*WTMO)/(XMLOMX*WTMO+(1.D+0-XMLOMX)*WTMW)
                SNW(2,NW) = STW(2,NW)*RHOLW(2,NW)*(YLOWX-XLOMX)/
     &            (RHONW(2,NW)+RHOLW(2,NW)*(YLOWX-XLOMX))
              ELSE
                SNX = 0.D+0
              ENDIF
              SLW(2,NW) = STW(2,NW)-SNW(2,NW)
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
!
!---        NAPL only well node  ---
!
            ELSEIF( RHOLX.GT.RHONX .AND. DNBWX.LE.DZBX ) THEN
              NPHAZW(2,NW) = 6
              STW(2,NW) = (WWD(NWL)-DZBX)/(DZTX-DZBX)
              SNW(2,NW) = STW(2,NW)
              SLW(2,NW) = 0.D+0
              CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
              CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
              YMLOW(2,NW) = PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
!
!---        Aqueous and NAPL well node  ---
!
            ELSEIF( DIWX.LT.DZTX .AND. DIWX.GT.DZBX ) THEN
              NPHAZW(2,NW) = 6
              STW(2,NW) = (WWD(NWL)-DZBX)/(DZTX-DZBX)
              IF( RHOLX.GT.RHONX ) THEN
                SLW(2,NW) = (DIWX-DZBX)/(WWD(NWL)-DZBX)
                SNW(2,NW) = MAX( 1.D+0-SLW(2,NW),0.D+0 )
              ELSEIF( RHOLX.LT.RHONX ) THEN
                SNW(2,NW) = (DIWX-DZBX)/(WWD(NWL)-DZBX)
                SLW(2,NW) = MAX( 1.D+0-SNW(2,NW),0.D+0 )
              ELSE
                SNW(2,NW) = STW(2,NW)*WWDN(NWL)/WWD(NWL)
                SLW(2,NW) = MAX( STW(2,NW)-SNW(2,NW),0.D+0 )
              ENDIF
              CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
              CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
              PWTX = PWX
              CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &          PWTX,PW(2,NW),PWBX,NW )
              PWX = PWBX
              YMLOW(2,NW) = PSO(2,N)/HCOW
              IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
            ENDIF
!
!---  Total-liquid unsaturated well node  ---
!
          ELSE
!
!---        DNAPL: Unsaturated  ---
!
            IF( RHONX.GT.RHOLX ) THEN
              NPHAZW(2,NW) = 3
!
!---        LNAPL: Unsaturated  ---
!
            ELSE
              NPHAZW(2,NW) = 4
            ENDIF
            STW(2,NW) = 0.D+0
            SLW(2,NW) = 0.D+0
            SNW(2,NW) = 0.D+0
            CALL WATLQD( T(2,N),PWX,RHOLW(2,NW) )
            CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONW(2,NW) )
            PW(2,NW) = PWX - PATM
            YMLOW(2,NW) = 0.D+0
            XMLOMX = PSO(2,N)/HCOW
            XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
          ENDIF
!
!---  Ratio aqueous/total liquid saturation  ---
!
          PWLW(2,NW) = SLW(2,NW)/(STW(2,NW)+EPSL)
   10   CONTINUE
   20 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of CHKW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFLUXW4
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
!     Water-Oil Mode
!     Compute dissolved-oil flux in well.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE JACOB
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
      INTEGER MP(5),MT(5),MW(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA MP /2,2,3,2,4/
      DATA MT /2,3,2,4,2/
      DATA MW /2,3,2,4,2/
      SAVE MP,MT,MW
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFLUXW4'
!
!---  Loop over well nodes  ---
!
      DO 1000 NWL = 1,NWLS
      I = IWLDM(1,NWL)
      J = IWLDM(2,NWL)
      KSX = IWLDM(3,NWL)
      KEX = IWLDM(4,NWL)
      RDW = WBR(NWL)
      DO 900 K = KSX,KEX-1
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NPZ = NSZW(NW)
        NQZ = NPZ+1
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
          RDE = XE(2,N)
        ELSE
          RDE = SQRT( AFZ(NSZ(N))/GPI )
        ENDIF
        ACWX = 2.D+0*GPI*RDW*DZGF(N)
        DRD2 = (RDE**2-RDW**2)
        NT = N+IJFLD
        NWT = ABS(IXW(NT))
        DXMLOZ = XMLOW(2,NW)-XMLOW(2,NWT)
        DXMLOX = XMLO(2,N)-XMLOW(2,NW)
        DO 100 M = 1,ISVF
!
!---      Vertical dissolved-oil diffusion within the well  ---
!
          INDX = -8
          SLWPX = SLW(MP(M),NW)
          SLWTX = SLW(MT(M),NWT)
          SLWX = DIFMN( SLWPX,SLWTX,DZGF(N),DZGF(NT),DXMLOZ,INDX )
          WTMX = XMLOW(MP(M),NW)*WTMO + (1.D+0-XMLOW(MP(M),NW))*WTMW
          RHOMLPX = RHOLW(MP(M),NW)/WTMX
          WTMX = XMLOW(MT(M),NWT)*WTMO + (1.D+0-XMLOW(MT(M),NWT))*WTMW
          RHOMLTX = RHOLW(MT(M),NWT)/WTMX
          DFP = DFLOW(MP(M),NW)*RHOMLPX
          DFT = DFLOW(MT(M),NWT)*RHOMLTX
          INDX = 12
          DFM = DIFMN( DFP,DFT,DZGF(N),DZGF(NT),DXMLOZ,INDX )
          DXMLOW = (XMLOW(MP(M),NW)-XMLOW(MT(M),NWT))
          WDLO_W(M,NQZ) = DFM/DZGP(NQZ)*DXMLOW*SLWX
!
!---      Axial dissolved-oil diffusion across the well screen
!         Positive flow defined into the well  ---
!
          IF( IXW(N).LT.0 ) THEN
            SLWWX = MAX(SLW(MW(M),NW),0.D+0)
            WTMX = XMLOW(MW(M),NW)*WTMO + (1.D+0-XMLOW(MW(M),NW))*WTMW
            RHOMLWX = RHOLW(MW(M),NW)/WTMX
            DFW = SLWWX*DFLOW(MW(M),NW)*RHOMLWX
            DFP = TORL(MP(M),N)*SL(MP(M),N)*DFLO(MP(M),N)*RHOML(MP(M),N)
            INDX = 12
            DFM = DIFMN( DFP,DFW,RDE,RDW,DXMLOX,INDX )
            GEOMFX = 2.D+0*GPI*DRD2*DZGF(N)/
     &       (((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
            UDLO_W(M,NW) = DFM*GEOMFX*(XMLO(MP(M),N)-XMLOW(MW(M),NW))
     &        /ACWX
          ENDIF
  100   CONTINUE
  900 CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFLUXW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLUXW4
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
!     Water-Oil Mode
!     Compute aqueous volumetric, NAPL volumetric flux in well.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE PORMED
      USE NAPL
      USE JACOB
      USE GRID
      USE FDVP
      USE FDVN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER MP(5),MT(5),MW(5)
!
!----------------------Data Statements---------------------------------!
!
      DATA MP /2,2,3,2,4/
      DATA MT /2,3,2,4,2/
      DATA MW /2,3,2,4,2/
      SAVE MP,MT,MW
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLUXW4'
!
!---  Well permeabilities, equivalent to
!     a saturated hydraulic conductivity of 100 m/s  ---
!
      PRNWX = 1.D-5/VISRN
      PRLWX = 1.D-5/VISRL
!
!---  Loop over well nodes  ---
!
      DO 1000 NWL = 1,NWLS
      I = IWLDM(1,NWL)
      J = IWLDM(2,NWL)
      KSX = IWLDM(3,NWL)
      KEX = IWLDM(4,NWL)
      RDW = WBR(NWL)
      QL_W(2,NWL) = 0.D+0
      QN_W(2,NWL) = 0.D+0
      QT_W(2,NWL) = 0.D+0
      DO 900 K = KSX,KEX-1
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NPZ = NSZW(NW)
        NQZ = NPZ+1
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
          RDE = XE(2,N)
        ELSE
          RDE = SQRT( AFZ(NSZ(N))/GPI )
        ENDIF
        ACWX = 2.D+0*GPI*RDW*DZGF(N)
        PERMX = SQRT( PERM(1,IZ(N))*PERM(2,IZ(N)) )
        DRD2 = (RDE**2-RDW**2)
        NT = N+IJFLD
        NWT = ABS(IXW(NT))
        DO 100 M = 1,ISVF
!
!---      Bulk fluid vertical flow within the well  ---
!
          HDBZ = PW(MP(M),NW)-PW(MT(M),NWT)-0.5D+0*GRAVZ*
     &      (RHOGW(MP(M),NW)*DZGF(NT)+RHOGW(MT(M),NWT)*DZGF(N))
          IF( M.EQ.1 ) HDB = HDBZ
          WB_WX = PRNWX*HDBZ/DZGP(NQZ)
!
!---      Bouyant NAPL vertical flow within the well  ---
!
          INDX = -1
          RHONX = DIFMN( RHONW(MP(M),NW),RHONW(MT(M),NWT),
     &      DZGF(N),DZGF(NT),HDB,INDX )
          RHOLX = DIFMN( RHOLW(MP(M),NW),RHOLW(MT(M),NWT),
     &      DZGF(N),DZGF(NT),HDB,INDX )
          WBN_WX = (RHOLX-RHONX)*GRAVZ*SQRT(PRNWX*PRLWX)
          INDX = -4
          SNX = DIFMN( SNW(MP(M),NW),SNW(MT(M),NWT),
     &      DZGF(N),DZGF(NT),WBN_WX,INDX )
          WBL_WX = -WBN_WX
          SLX = DIFMN( SLW(MP(M),NW),SLW(MT(M),NWT),
     &      DZGF(N),DZGF(NT),WBL_WX,INDX )
          WBN_WX = WBN_WX*MIN( SNX,SLX )
!
!---      NAPL vertical flow within the well  ---
!
          RKNPX = SNW(MP(M),NW)
          RKNTX = SNW(MT(M),NWT)
          INDX = -4
          RKNX = DIFMN( RKNPX,RKNTX,DZGF(N),DZGF(NT),HDB,INDX )
          WN_W(M,NQZ) = WB_WX*RKNX + WBN_WX
!
!---      Aqueous vertical flow within the well  ---
!
          RKLPX = SLW(MP(M),NW)
          RKLTX = SLW(MT(M),NWT)
          INDX = -4
          RKLX = DIFMN( RKLPX,RKLTX,DZGF(N),DZGF(NT),HDB,INDX )
          WL_W(M,NQZ) = WB_WX*RKLX - WBN_WX
!
!---      Axial flow across the well screen
!         Positive flow defined into the well  ---
!
          IF( IXW(N).LT.0 ) THEN
           DPNPW = PN(MP(M),N)-PW(MW(M),NW)
           IF( M.EQ.1 ) DPNX = DPNPW
           AFNWX = SNW(MW(M),NW)
           AFNPX = RKN(MP(M),N)
           INDX = -4
           AFNX = DIFMN( AFNPX,AFNWX,RDE,RDW,DPNX,INDX )
           INDX = 6
           VNM = DIFMN(VISN(MP(M),N),VISNW(MW(M),NW),RDE,RDW,DPNX,INDX)
           HCNX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &       (VNM*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
           UN_W(M,NW) = AFNX*HCNX*DPNPW/ACWX
           DPLPW = PL(MP(M),N)-PW(MW(M),NW)
           IF( M.EQ.1 ) DPLX = DPLPW
           AFLWX = SLW(MW(M),NW)
           AFLPX = SQRT( RKL(1,MP(M),N)*RKL(2,MP(M),N) )
           INDX = -4
           AFLX = DIFMN( AFLPX,AFLWX,RDE,RDW,DPLX,INDX )
           INDX = 5
           VLM = DIFMN(VISL(MP(M),N),VISLW(MW(M),NW),RDE,RDW,DPLX,INDX)
           HCLX = 2.D+0*GPI*PERMX*DRD2*DZGF(N)/
     &       (VLM*((RDE**2)*LOG(RDE/RDW)-5.D-1*DRD2))
           UL_W(M,NW) = AFLX*HCLX*DPLPW/ACWX
           IF( M.EQ.1 ) THEN
             QL_W(2,NWL) = QL_W(2,NWL) + UL_W(M,NW)*ACWX
             QN_W(2,NWL) = QN_W(2,NWL) + UN_W(M,NW)*ACWX
             QT_W(2,NWL) = QT_W(2,NWL) + QL_W(2,NWL) + QN_W(2,NWL)
           ENDIF
          ENDIF
  100   CONTINUE
  900 CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLUXW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBOW4
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
!     Load the Jacobian matrix for the well-oil equation with
!     aqueous-phase and NAPL contributions
!
!     RWS - residual no variance
!     RWP - residual variance in well node
!     RWB - residual variance in bottom well node
!     RWT - residual variance in top well node
!     RWF - residual variance in field node
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!     well_wo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE JACOB
      USE GRID
      USE FILES
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
      REAL*8 STWX(LUK+1),FW(2*LUK+1,3)
      REAL*8 RWP(LUK),RWB(LUK),RWT(LUK),RWF(LUK)
      INTEGER MP(5),MW(5),MB(5),MT(5)
      INTEGER MJ(2),MA(2),MS(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE MW,MP,MB,MT,MS,MJ,MA
      DATA MW /2,3,2,4,2/
      DATA MP /2,2,3,2,4/
      DATA MB /1,2,3,4,5/
      DATA MT /1,3,2,5,4/
      DATA MS /2,3,4/
      DATA MJ /2,4/
      DATA MA /3,5/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBOW4'
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
!
!---  Loop over all nodes containing cased or screened well sections ---
!
      DO 1000 NWL = 1,NWLS
      I = IWLDM(1,NWL)
      J = IWLDM(2,NWL)
      KSX = IWLDM(3,NWL)
      KEX = IWLDM(4,NWL)
      RDW = WBR(NWL)
      RDWS = WBRS(NWL)
      DO 900 K = KSX,KEX
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NPZ = NSZW(NW)
        NQZ = NPZ+1
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
          RDE = XE(2,N)
        ELSE
          RDE = SQRT( AFZ(NSZ(N))/GPI )
        ENDIF
        ACWX = 2.D+0*GPI*RDW*DZGF(N)
        ARWX = GPI*(RDW**2)
        ARWSX = GPI*(RDWS**2)
        VOLX = ARWSX*DZGF(N)
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = (XLOW(1,NW)*RHOLW(1,NW)*SLW(1,NW) +
     &    RHONW(1,NW)*SNW(1,NW))
        DO 100 M = 1,ISVC+1
          STW0 = (XLOW(MS(M),NW)*RHOLW(MS(M),NW)*SLW(MS(M),NW) +
     &      RHONW(MS(M),NW)*SNW(MS(M),NW))
          STWX(M) = (STW0-STW1)*DTI*VOLX
  100   CONTINUE
!
!---  Initialize well fluxes  ---
!
        DO 120 MD = 1,3
        DO 120 M = 1,ISVF
          FW(M,MD) = 0.D+0
  120   CONTINUE
!
!---  Screen flow into well  ---
!
        IF( IXW(N).LT.0 ) THEN
          DO 200 M = 1,ISVF
            IF( UL_W(M,NW).GT.0.D+0 ) THEN
              RHOLWX = RHOL(MP(M),N)*XLO(MP(M),N)
            ELSE
              RHOLWX = RHOLW(MW(M),NW)*XLOW(MW(M),NW)
            ENDIF
            IF( UN_W(M,NW).GT.0.D+0 ) THEN
              RHONWX = RHON(MP(M),N)
            ELSE
              RHONWX = RHONW(MW(M),NW)
            ENDIF
            FW(M,2) = -ACWX*(UL_W(M,NW)*RHOLWX + UN_W(M,NW)*RHONWX
     &        + UDLO_W(M,NW)*WTMO)
  200     CONTINUE
        ENDIF
!
!---  Well flow (bottom flux)  ---
!
        IF( K.GT.KSX ) THEN
          NB = N-IJFLD
          NWB = ABS(IXW(NB))
          DO 300 M = 1,ISVF
            FLWB = XLOW(MP(M),NWB)*RHOLW(MP(M),NWB)
            FLWP = XLOW(MW(M),NW)*RHOLW(MW(M),NW)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL_W(1,NPZ),INDX )
            FNWB = RHONW(MP(M),NWB)
            FNWP = RHONW(MW(M),NW)
            INDX = 3
            FNW = DIFMN( FNWB,FNWP,DZGF(NB),DZGF(N),WN_W(1,NPZ),INDX )
            FW(M,1) = -ARWX*(WL_W(MB(M),NPZ)*FLW + WN_W(MB(M),NPZ)*FNW
     &        + WDLO_W(MB(M),NPZ)*WTMO)
  300     CONTINUE
        ENDIF
!
!---  Well flow (top flux)  ---
!
        IF( K.LT.KEX ) THEN
          NT = N+IJFLD
          NWT = ABS(IXW(NT))
          DO 400 M = 1,ISVF
            FLWT = XLOW(MP(M),NWT)*RHOLW(MP(M),NWT)
            FLWP = XLOW(MW(M),NW)*RHOLW(MW(M),NW)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL_W(1,NQZ),INDX )
            FNWT = RHONW(MP(M),NWT)
            FNWP = RHONW(MW(M),NW)
            INDX = 3
            FNW = DIFMN( FNWP,FNWT,DZGF(N),DZGF(NT),WN_W(1,NQZ),INDX )
            FW(M,3) = ARWX*(WL_W(MT(M),NQZ)*FLW + WN_W(MT(M),NQZ)*FNW
     &        + WDLO_W(MT(M),NQZ)*WTMO)
  400     CONTINUE
        ENDIF
!
!---  Compute oil equation residuals
!     No variance  ---
!
      RWS = STWX(1) - SRCO_W(2,NW)
      DO 410 L = 1,3
        RWS = RWS + FW(1,L)
  410 CONTINUE
!
!---  Variance in well node  ---
!
      DO 420 M = 1,ISVC
        RWP(M) = STWX(M+1) - SRCO_W(M+2,NW)
  420 CONTINUE
      DO 424 M = 1,ISVC
        DO 422 L = 1,3
          RWP(M) = RWP(M) + FW(MJ(M),L)
  422   CONTINUE
  424 CONTINUE
!
!---  Variance in bottom well node  ---
!
      DO 430 M = 1,ISVC
        RWB(M) = RWS - FW(1,1) + FW(MA(M),1)
  430 CONTINUE
!
!---  Variance in field node  ---
!
      DO 440 M = 1,ISVC
        RWF(M) = RWS - FW(1,2) + FW(MA(M),2)
  440 CONTINUE
!
!---  Variance in top well node  ---
!
      DO 450 M = 1,ISVC
        RWT(M) = RWS - FW(1,3) + FW(MA(M),3)
  450 CONTINUE
!
!---  Load Jacobian matrix  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---  Well-oil equation (well node variables)   ---
!
        NMD = IXP(N)
        MJP = IM(IEQOW,NMD)
        DO 500 M = 1,ISVC
          MCOL = IM(M+ISVC,NMD)
          MROW = MJP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWP(M)-RWS)/DNRW(M,NW)
  500   CONTINUE
        BLU(MJP) = BLU(MJP) - RWS
        RSDL(IEQOW,N) = BLU(MJP)
!
!---  Well-oil equation (bottom well variables)   ---
!
        IF( K.NE.KSX ) THEN
          NB = N-IJFLD
          NMD = IXP(NB)
          DO 510 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWB(M)-RWS)/DNRW(M,NWB)
  510     CONTINUE
        ENDIF
!
!---  Well-oil equation (top well variables)   ---
!
        IF( K.NE.KEX ) THEN
          NT = N+IJFLD
          NMD = IXP(NT)
          DO 520 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWT(M)-RWS)/DNRW(M,NWT)
  520     CONTINUE
        ENDIF
!
!---  Well-oil equation (node field variables)
!     for screened well sections   ---
!
        IF( IXW(N).LT.0 ) THEN
          NMD = IXP(N)
          MJP = IM(IEQOW,NMD)
          DO 530 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWF(M)-RWS)/DNR(M,N)
  530     CONTINUE
!
!---  Field-oil equation (node well variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQO,NMD)
          DO 540 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)-(FW(MJ(M),2)-FW(1,2))/
     &        DNRW(M,NW)
  540     CONTINUE
!
!---  Field-oil equation (node field variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQO,NMD)
          DO 550 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)-(FW(MA(M),2)-FW(1,2))/
     &        DNR(M,N)
  550     CONTINUE
          BLU(MJP) = BLU(MJP) + FW(1,2)
          RSDL(IEQO,N) = BLU(MJP)
        ENDIF
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---  Well-oil equation (well node variables)   ---
!
        NMD = IXP(N)
        MJP = IM(IEQOW,NMD)
        DO 600 M = 1,ISVC
          MCOL = KLU(MJP,M+5*ISVC)
          DLU(MCOL) = DLU(MCOL)+(RWP(M)-RWS)/DNRW(M,NW)
  600   CONTINUE
        BLU(MJP) = BLU(MJP) - RWS
        RSDL(IEQOW,N) = BLU(MJP)
!
!---  Well-oil equation (bottom well variables)   ---
!
        IF( K.NE.KSX ) THEN
          NB = N-IJFLD
          NMD = IXP(NB)
          DO 610 M = 1,ISVC
            MCOL = KLU(MJP,M+ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWB(M)-RWS)/DNRW(M,NWB)
  610     CONTINUE
        ENDIF
!
!---  Well-oil equation (top well variables)   ---
!
        IF( K.NE.KEX ) THEN
          NT = N+IJFLD
          NMD = IXP(NT)
          DO 620 M = 1,ISVC
            MCOL = KLU(MJP,M+9*ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWT(M)-RWS)/DNRW(M,NWT)
  620     CONTINUE
        ENDIF
!
!---  Well-oil equation (node field variables)
!     for screened well sections   ---
!
        IF( IXW(N).LT.0 ) THEN
          NMD = IXP(N)
          MJP = IM(IEQOW,NMD)
          DO 630 M = 1,ISVC
            MCOL = KLU(MJP,M+4*ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWF(M)-RWS)/DNR(M,N)
  630     CONTINUE
!
!---  Field-oil equation (node well variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQO,NMD)
          DO 640 M = 1,ISVC
            MCOL = KLU(MJP,M+5*ISVC)
            DLU(MCOL) = DLU(MCOL)-(FW(MJ(M),2)-FW(1,2))/DNRW(M,NW)
  640     CONTINUE
!
!---  Field-oil equation (node field variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQO,NMD)
          DO 650 M = 1,ISVC
            MCOL = KLU(MJP,M+4*ISVC)
            DLU(MCOL) = DLU(MCOL)-(FW(MA(M),2)-FW(1,2))/DNR(M,N)
  650     CONTINUE
          BLU(MJP) = BLU(MJP) + FW(1,2)
          RSDL(IEQO,N) = BLU(MJP)
        ENDIF
      ELSE
        WRITE(ISC,'(A)') 'ERROR: Unknown Linear Equation Solver'
        WRITE(IWR,'(A)') 'ERROR: Unknown Linear Equation Solver'
        STOP
      ENDIF
  900 CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBOW4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBWW4
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
!     Load the Jacobian matrix for the well-water equation with
!     aqueous-phase and NAPL contributions
!
!     RWS - residual no variance
!     RWP - residual variance in well node
!     RWB - residual variance in bottom well node
!     RWT - residual variance in top well node
!     RWF - residual variance in field node
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!     well_wo.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
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
      REAL*8 STWX(LUK+1),FW(2*LUK+1,3)
      REAL*8 RWP(LUK),RWB(LUK),RWT(LUK),RWF(LUK)
      INTEGER MP(5),MW(5),MB(5),MT(5)
      INTEGER MJ(2),MA(2),MS(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE MW,MP,MB,MT,MS,MJ,MA
      DATA MW /2,3,2,4,2/
      DATA MP /2,2,3,2,4/
      DATA MB /1,2,3,4,5/
      DATA MT /1,3,2,5,4/
      DATA MS /2,3,4/
      DATA MJ /2,4/
      DATA MA /3,5/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBWW4'
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
!
!---  Loop over all nodes containing cased or screened well sections ---
!
      DO 1000 NWL = 1,NWLS
      I = IWLDM(1,NWL)
      J = IWLDM(2,NWL)
      KSX = IWLDM(3,NWL)
      KEX = IWLDM(4,NWL)
      RDW = WBR(NWL)
      RDWS = WBRS(NWL)
      DO 900 K = KSX,KEX
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        NPZ = NSZW(NW)
        NQZ = NPZ+1
        IF( (ICS.EQ.2 .OR. ICS.EQ.6) .AND. (I.EQ.1) ) THEN
          RDE = XE(2,N)
        ELSE
          RDE = SQRT( AFZ(NSZ(N))/GPI )
        ENDIF
        ACWX = 2.D+0*GPI*RDW*DZGF(N)
        ARWX = GPI*(RDW**2)
        ARWSX = GPI*(RDWS**2)
        VOLX = ARWSX*DZGF(N)
!
!---    First-order, forward-difference, time differential  ---
!
        STW1 = XLWW(1,NW)*RHOLW(1,NW)*SLW(1,NW)
        DO 100 M = 1,ISVC+1
          STW0 = XLWW(MS(M),NW)*RHOLW(MS(M),NW)*SLW(MS(M),NW)
          STWX(M) = (STW0-STW1)*DTI*VOLX
  100   CONTINUE
!
!---  Initialize well fluxes  ---
!
        DO 120 MD = 1,3
        DO 120 M = 1,ISVF
          FW(M,MD) = 0.D+0
  120   CONTINUE
!
!---  Screen flow into well  ---
!
        IF( IXW(N).LT.0 ) THEN
          DO 200 M = 1,ISVF
            IF( UL_W(M,NW).GT.0.D+0 ) THEN
              RHOLWX = XLW(MP(M),N)*RHOL(MP(M),N)
            ELSE
              RHOLWX = XLWW(MW(M),NW)*RHOLW(MW(M),NW)
            ENDIF
            FW(M,2) = -ACWX*(UL_W(M,NW)*RHOLWX - UDLO_W(M,NW)*WTMW)
  200     CONTINUE
        ENDIF
!
!---  Well flow (bottom flux)  ---
!
        IF( K.GT.KSX ) THEN
          NB = N-IJFLD
          NWB = ABS(IXW(NB))
          DO 300 M = 1,ISVF
            FLWB = XLWW(MP(M),NWB)*RHOLW(MP(M),NWB)
            FLWP = XLWW(MW(M),NW)*RHOLW(MW(M),NW)
            INDX = 2
            FLW = DIFMN( FLWB,FLWP,DZGF(NB),DZGF(N),WL_W(1,NPZ),INDX )
            FW(M,1) = -ARWX*(WL_W(MB(M),NPZ)*FLW-WDLO_W(MB(M),NPZ)*WTMW)
  300     CONTINUE
        ENDIF
!
!---  Well flow (top flux)  ---
!
        IF( K.LT.KEX ) THEN
          NT = N+IJFLD
          NWT = ABS(IXW(NT))
          DO 400 M = 1,ISVF
            FLWT = XLWW(MP(M),NWT)*RHOLW(MP(M),NWT)
            FLWP = XLWW(MW(M),NW)*RHOLW(MW(M),NW)
            INDX = 2
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(NT),WL_W(1,NQZ),INDX )
            FW(M,3) = ARWX*(WL_W(MT(M),NQZ)*FLW-WDLO_W(MT(M),NQZ)*WTMW)
  400     CONTINUE
        ENDIF
!
!---  Compute oil equation residuals
!     No variance  ---
!
      RWS = STWX(1) - SRCW_W(2,NW)
      DO 410 L = 1,3
        RWS = RWS + FW(1,L)
  410 CONTINUE
!
!---  Variance in well node  ---
!
      DO 420 M = 1,ISVC
        RWP(M) = STWX(M+1) - SRCW_W(M+2,NW)
  420 CONTINUE
      DO 424 M = 1,ISVC
        DO 422 L = 1,3
          RWP(M) = RWP(M) + FW(MJ(M),L)
  422   CONTINUE
  424 CONTINUE
!
!---  Variance in bottom well node  ---
!
      DO 430 M = 1,ISVC
        RWB(M) = RWS - FW(1,1) + FW(MA(M),1)
  430 CONTINUE
!
!---  Variance in field node  ---
!
      DO 440 M = 1,ISVC
        RWF(M) = RWS - FW(1,2) + FW(MA(M),2)
  440 CONTINUE
!
!---  Variance in top well node  ---
!
      DO 450 M = 1,ISVC
        RWT(M) = RWS - FW(1,3) + FW(MA(M),3)
  450 CONTINUE
!
!---  Load Jacobian matrix  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---  Well-water equation (well node variables)   ---
!
        NMD = IXP(N)
        MJP = IM(IEQWW,NMD)
        DO 500 M = 1,ISVC
          MCOL = IM(M+ISVC,NMD)
          MROW = MJP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWP(M)-RWS)/DNRW(M,NW)
  500   CONTINUE
        BLU(MJP) = BLU(MJP) - RWS
        RSDL(IEQWW,N) = BLU(MJP)
!
!---  Well-water equation (bottom well variables)   ---
!
        IF( K.NE.KSX ) THEN
          NB = N-IJFLD
          NMD = IXP(NB)
          DO 510 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWB(M)-RWS)/DNRW(M,NWB)
  510     CONTINUE
        ENDIF
!
!---  Well-water equation (top well variables)   ---
!
        IF( K.NE.KEX ) THEN
          NT = N+IJFLD
          NMD = IXP(NT)
          DO 520 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWT(M)-RWS)/DNRW(M,NWT)
  520     CONTINUE
        ENDIF
!
!---  Well-water equation (node field variables)
!     for screened well sections   ---
!
        IF( IXW(N).LT.0 ) THEN
          NMD = IXP(N)
          MJP = IM(IEQWW,NMD)
          DO 530 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)+(RWF(M)-RWS)/DNR(M,N)
  530     CONTINUE
!
!---  Field-water equation (node well variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQW,NMD)
          DO 540 M = 1,ISVC
            MCOL = IM(M+ISVC,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)-(FW(MJ(M),2)-FW(1,2))/
     &        DNRW(M,NW)
  540     CONTINUE
!
!---  Field-water equation (node field variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQW,NMD)
          DO 550 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MJP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL)-(FW(MA(M),2)-FW(1,2))/
     &        DNR(M,N)
  550     CONTINUE
          BLU(MJP) = BLU(MJP) + FW(1,2)
          RSDL(IEQW,N) = BLU(MJP)
        ENDIF
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---  Well-water equation (well node variables)   ---
!
        NMD = IXP(N)
        MJP = IM(IEQWW,NMD)
        DO 600 M = 1,ISVC
          MCOL = KLU(MJP,M+5*ISVC)
          DLU(MCOL) = DLU(MCOL)+(RWP(M)-RWS)/DNRW(M,NW)
  600   CONTINUE
        BLU(MJP) = BLU(MJP) - RWS
        RSDL(IEQWW,N) = BLU(MJP)
!
!---  Well-water equation (bottom well variables)   ---
!
        IF( K.NE.KSX ) THEN
          NB = N-IJFLD
          NMD = IXP(NB)
          DO 610 M = 1,ISVC
            MCOL = KLU(MJP,M+ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWB(M)-RWS)/DNRW(M,NWB)
  610     CONTINUE
        ENDIF
!
!---  Well-water equation (top well variables)   ---
!
        IF( K.NE.KEX ) THEN
          NT = N+IJFLD
          NMD = IXP(NT)
          DO 620 M = 1,ISVC
            MCOL = KLU(MJP,M+9*ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWT(M)-RWS)/DNRW(M,NWT)
  620     CONTINUE
        ENDIF
!
!---  Well-water equation (node field variables)
!     for screened well sections   ---
!
        IF( IXW(N).LT.0 ) THEN
          NMD = IXP(N)
          MJP = IM(IEQWW,NMD)
          DO 630 M = 1,ISVC
            MCOL = KLU(MJP,M+4*ISVC)
            DLU(MCOL) = DLU(MCOL)+(RWF(M)-RWS)/DNR(M,N)
  630     CONTINUE
!
!---  Field-water equation (node well variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQW,NMD)
          DO 640 M = 1,ISVC
            MCOL = KLU(MJP,M+5*ISVC)
            DLU(MCOL) = DLU(MCOL)-(FW(MJ(M),2)-FW(1,2))/DNRW(M,NW)
  640     CONTINUE
!
!---  Field-water equation (node field variables)
!     for screened well sections   ---
!
          NMD = IXP(N)
          MJP = IM(IEQW,NMD)
          DO 650 M = 1,ISVC
            MCOL = KLU(MJP,M+4*ISVC)
            DLU(MCOL) = DLU(MCOL)-(FW(MA(M),2)-FW(1,2))/DNR(M,N)
  650     CONTINUE
          BLU(MJP) = BLU(MJP) + FW(1,2)
          RSDL(IEQW,N) = BLU(MJP)
        ENDIF
      ELSE
        WRITE(ISC,'(A)') 'ERROR: Unknown Linear Equation Solver'
        WRITE(IWR,'(A)') 'ERROR: Unknown Linear Equation Solver'
        STOP
      ENDIF
  900 CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBWW4 group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE LDOW4
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
!     Water-Oil Mode
!     Load old time step well arrays.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/LDOW4'
!
!---  Loop over well nodes  ---
!
      DO 1000 NWL = 1,NWLS
      I = IWLDM(1,NWL)
      J = IWLDM(2,NWL)
      KSX = IWLDM(3,NWL)
      KEX = IWLDM(4,NWL)
      QL_W(1,NWL) = QL_W(2,NWL)
      QN_W(1,NWL) = QN_W(2,NWL)
      QT_W(1,NWL) = QT_W(2,NWL)
      DO 900 K = KSX,KEX
        N = ND(I,J,K)
        NW = ABS(IXW(N))
        IF( NPHAZW(2,NW).EQ.3.OR.NPHAZW(2,NW).EQ.4 ) YMLOW(2,NW) = 0.D+0
        PW(1,NW) = PW(2,NW)
        PWLW(1,NW) = PWLW(2,NW)
        RHOLW(1,NW) = RHOLW(2,NW)
        RHONW(1,NW) = RHONW(2,NW)
        XLWW(1,NW) = XLWW(2,NW)
        XLOW(1,NW) = XLOW(2,NW)
        SLW(1,NW) = SLW(2,NW)
        SNW(1,NW) = SNW(2,NW)
        STW(1,NW) = STW(2,NW)
        XMLOW(1,NW) = XMLOW(2,NW)
        YMLOW(1,NW) = YMLOW(2,NW)
        VISLW(1,NW) = VISLW(2,NW)
        VISNW(1,NW) = VISNW(2,NW)
        NPHAZW(1,NW) = NPHAZW(2,NW)
  900 CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of LDOW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROPW4
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
!     Water-Oil Mode
!     Compute thermo-physical properties for well fluids.
!
!     Phase Conditions
!     Number              Description            Primary Unknowns
!       1            Aqueous w/ Dissolved Oil        PW,YMLOW
!       2            Aqueous-Gas w/ Dissolved Oil    PW,YMLOW
!       3            DNAPL Unsaturated               SLW,SNW
!       4            LNAPL Unsaturated               SLW,SNW
!       5            Aqueous-NAPL                    PWLW,PW
!       6            Aqueous-NAPL-Gas                PWLW,PW
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 12 October 2001.
!     Last Modified by M.D. White, PNNL, 12 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE JACOB
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
      REAL*8 WWDX(LSV,LNW)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROPW4'
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
      SW_TOL = 1.D-4
!
!---  Loop over wells  ---
!
      DO 100 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        PWX = PWB(2,NWL) + PATM
!
!---  Loop over well nodes from the top down  ---
!
        DO 95 K = KEX,KSX,-1
          N = ND(I,J,K)
          NW = ABS(IXW(N))
!
!---    Saturated or partially saturated dissolved-oil conditions
!         Water equation: well pressure
!         Oil equation:  mole fraction of dissolved oil  ---
!
          IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
            PWPX = PW(2,NW) + PATM
            CALL WATLQD( T(2,N),PWPX,RHOLWX )
            CALL VOCLQD( T(2,N),PWPX,PSO(2,N),RHONWX )
            XMLOMX = PSO(2,N)/HCOW
            XMLOW(2,NW) = MIN( YMLOW(2,NW),XMLOMX )
            CALL SUPSW4( DZGF(N),PW(2,NW),PWB(2,NWL),RHOLWX,RHONWX,
     &        SLW(2,NW),SNW(2,NW),STW(2,NW),XMLOMX,YMLOW(2,NW) )
            PWLW(2,NW) = SLW(2,NW)/(STW(2,NW)+EPSL)
            IF( SLW(2,NW).LT.EPSL ) THEN
              IF( RHONWX.GT.RHOLWX )THEN
                NPHAZW(2,NW) = 3
              ELSE
                NPHAZW(2,NW) = 4
              ENDIF
            ELSEIF( STW(2,NW).GE.1.D+0 ) THEN
              IF( YMLOW(2,NW).GT.(1.5D+0*(PSO(2,N)/HCOW)) )THEN
                NPHAZW(2,NW) = 5
              ELSE
                NPHAZW(2,NW) = 1
              ENDIF
            ELSE
              IF( YMLOW(2,NW).GT.(1.5D+0*(PSO(2,N)/HCOW)) )THEN
                NPHAZW(2,NW) = 6
              ELSE
                NPHAZW(2,NW) = 2
              ENDIF
            ENDIF
!
!---    Unsaturated conditions
!         Water equation: aqueous saturation
!         Oil equation:  NAPL saturation  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
            CALL WATLQD( T(2,N),PWX,RHOLWX )
            CALL VOCLQD( T(2,N),PWX,PSO(2,N),RHONWX )
            STW(2,NW) = SNW(2,NW) + SLW(2,NW)
            PWLW(2,NW) = SLW(2,NW)/(STW(2,NW)+EPSL)
            PWTX = PWX
            CALL PRESSW4( SLW(2,NW),SNW(2,NW),T(2,N),
     &        PWTX,PW(2,NW),PWBX,NW )
            IF( STW(2,NW).LT.SW_TOL ) THEN
              IF( RHONWX.GT.RHOLWX )THEN
                NPHAZW(2,NW) = 3
              ELSE
                NPHAZW(2,NW) = 4
              ENDIF
            ELSEIF( STW(2,NW).GE.1.D+0 ) THEN
              IF( SNW(2,NW).GT.SW_TOL ) THEN
                NPHAZW(2,NW) = 5
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZW(2,NW) = 1
              ELSE
                NPHAZW(2,NW) = 5
              ENDIF
            ELSE
              IF( SNW(2,NW).GT.SW_TOL ) THEN
                NPHAZW(2,NW) = 6
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZW(2,NW) = 2
              ELSE
                NPHAZW(2,NW) = 6
              ENDIF
            ENDIF
!
!---    Saturated or partially saturated NAPL conditions
!         Water equation: ratio aqueous over total-liquid saturation
!         Oil equation:  well pressure  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
            PWPX = PW(2,NW) + PATM
            CALL WATLQD( T(2,N),PWPX,RHOLWX )
            CALL VOCLQD( T(2,N),PWPX,PSO(2,N),RHONWX )
            CALL STW4( DZGF(N),PW(2,NW),PWB(2,NWL),PWLW(2,NW),
     &        RHOLWX,RHONWX,SLW(2,NW),SNW(2,NW),STW(2,NW) )
            IF( STW(2,NW).LT.EPSL ) THEN
              IF( RHONWX.GT.RHOLWX )THEN
                NPHAZW(2,NW) = 3
              ELSE
                NPHAZW(2,NW) = 4
              ENDIF
            ELSEIF( STW(2,NW).GE.1.D+0 ) THEN
              IF( SNW(2,NW).GT.EPSL ) THEN
                NPHAZW(2,NW) = 5
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZW(2,NW) = 1
              ELSE
                NPHAZW(2,NW) = 5
              ENDIF
            ELSE
              IF( SNW(2,NW).GT.EPSL ) THEN
                NPHAZW(2,NW) = 6
              ELSEIF( PSO(2,N)/HCOW.GT.1.D-12 ) THEN
                NPHAZW(2,NW) = 2
              ELSE
                NPHAZW(2,NW) = 6
              ENDIF
            ENDIF
          ENDIF
          RHOWX = (SLW(2,NW)*RHOLWX + SNW(2,NW)*RHONWX)/
     &      (SLW(2,NW) + SNW(2,NW) + EPSL)
          STWX = MIN( MAX( STW(2,NW),0.D+0 ),1.D+0 )
          PWX = PWX + STWX*DZGF(N)*GRAVZ*RHOWX
!
!---    Primary variable increments  ---
!
          IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
            DNRW(IEQW,NW) = MAX( 1.D-1,1.D-6*(PW(2,NW)+PATM) )
            XMLOX = PSO(2,N)/HCOW
            DNRW(IEQO,NW) = SIGN( 1.D-6*XMLOX,5.D-1*XMLOX-YMLOW(2,NW) )
            DNRW(IEQO,NW) = 1.D-5*XMLOX
          ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
            DNRW(IEQW,NW) = 1.D-6
            DNRW(IEQO,NW) = 1.D-6
          ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
            DNRX = 1.D-6
            DNRW(IEQW,NW) = SIGN( DNRX,5.D-1-PWLW(2,NW) )
            DNRW(IEQO,NW) = MAX( 1.D-1,1.D-6*(PW(2,NW)+PATM) )
          ENDIF
!
!---      Increment primary variables  ---
!
          DO 90 M = 3,2*ISVC
            YMLOW(M,NW) = YMLOW(2,NW)
            SNW(M,NW) = SNW(2,NW)
            SLW(M,NW) = SLW(2,NW)
            PW(M,NW) = PW(2,NW)
            PWLW(M,NW) = PWLW(2,NW)
!
!---      Aqueous or aqueous-gas w/ dissolved-oil conditions
!           Water equation: well pressure
!           Oil equation:  mole fraction of dissolved oil  ---
!
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
              IF( M.EQ.IEQWW ) THEN
                PW(M,NW) = PW(M,NW) + DNRW(IEQW,NW)
              ELSEIF( M.EQ.IEQOW ) THEN
                YMLOW(M,NW) = YMLOW(M,NW) + DNRW(IEQO,NW)
              ENDIF
!
!---      Desaturated conditions
!           Water equation: aqueous saturation
!           Oil equation:  NAPL saturation  ---
!
            ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
              IF( M.EQ.IEQWW ) THEN
                SLW(M,NW) = SLW(M,NW) + DNRW(IEQW,NW)
              ELSEIF( M.EQ.IEQOW ) THEN
                SNW(M,NW) = SNW(M,NW) + DNRW(IEQO,NW)
              ENDIF
!
!---      Aqueous-NAPL conditions
!           Water equation: ratio aqueous over total-liquid saturation
!           Oil equation:  well pressure  ---
!
            ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
              IF( M.EQ.IEQWW ) THEN
                PWLW(M,NW) = PWLW(M,NW) + DNRW(IEQW,NW)
              ELSEIF( M.EQ.IEQOW ) THEN
                PW(M,NW) = PW(M,NW) + DNRW(IEQO,NW)
              ENDIF
            ENDIF
   90     CONTINUE
   95   CONTINUE
  100 CONTINUE
!
!---  Loop over wells  ---
!
      DO 300 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        PWX = PWB(2,NWL) + PATM
!
!---    Loop over well nodes from the top down  ---
!
        DO 250 K = KEX,KSX,-1
          N = ND(I,J,K)
          NW = ABS(IXW(N))
!
!---      Compute well physical properties  ---
!
          DO 200 M = 2,ISVC+2
!
!---        Aqueous or aqueous-gas w/ dissolved-oil conditions
!             Water equation: well pressure
!             Oil equation:  mole fraction of dissolved oil  ---
!
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
!
!---          Compute the aqueous and NAPL density  ---
!
              PWPX = PW(M,NW) + PATM
              CALL WATLQD( T(2,N),PWPX,RHOLW(M,NW) )
              CALL VOCLQD( T(2,N),PWPX,PSO(2,N),RHONW(M,NW) )
              IF( YMLOW(M,NW).LT.EPSL ) YMLOW(M,NW) = 0.D+0
              XMLOMX = PSO(2,N)/HCOW
              XMLOW(M,NW) = MIN( YMLOW(M,NW),XMLOMX )
              CALL SUPSW4( DZGF(N),PW(M,NW),PWB(2,NWL),RHOLW(M,NW),
     &          RHONW(M,NW),SLW(M,NW),SNW(M,NW),STW(M,NW),XMLOMX,
     &          YMLOW(M,NW) )
              PWLW(M,NW) = SLW(M,NW)/(STW(M,NW)+EPSL)
!
!---        Deaturated conditions
!             Water equation: aqueous saturation
!             Oil equation:  NAPL saturation  ---
!
            ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
              XMLOW(M,NW) = 0.D+0
              PWTX = PWX
              CALL PRESSW4( SLW(M,NW),SNW(M,NW),T(2,N),
     &          PWTX,PW(M,NW),PWBX,NW )
!
!---          Compute the aqueous and NAPL density  ---
!
              PWPX = PW(M,NW) + PATM
              CALL WATLQD( T(2,N),PWPX,RHOLW(M,NW) )
              CALL VOCLQD( T(2,N),PWPX,PSO(2,N),RHONW(M,NW) )
              STW(M,NW) = SLW(M,NW) + SNW(M,NW)
              STW(M,NW) = MIN( MAX( STW(M,NW),0.D+0 ),1.D+0 )
              PWLW(M,NW) = SLW(M,NW)/(STW(M,NW)+EPSL)
!
!---        Partially saturated conditions
!             Water equation: aqueous saturation
!             Oil equation:  well pressure  ---
!
            ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
!
!---          Compute the aqueous and NAPL density  ---
!
              PWPX = PW(M,NW) + PATM
              CALL WATLQD( T(2,N),PWPX,RHOLW(M,NW) )
              CALL VOCLQD( T(2,N),PWPX,PSO(2,N),RHONW(M,NW) )
              YMLOW(M,NW) = PSO(2,N)/HCOW
              IF( YMLOW(M,NW).LT.EPSL ) YMLOW(M,NW) = 0.D+0
              XMLOW(M,NW) = YMLOW(M,NW)
              CALL STW4( DZGF(N),PW(M,NW),PWB(2,NWL),PWLW(M,NW),
     &          RHOLW(M,NW),RHONW(M,NW),SLW(M,NW),SNW(M,NW),STW(M,NW) )
            ENDIF
!
!---        Compute the aqueous-phase mole fraction
!           and mass fraction  ---
!
            XMLWWX = MAX( 1.D+0-XMLOW(M,NW),0.D+0 )
            IF( XMLWWX.LT.EPSL ) XMLWWX = 0.D+0
            XLOW(M,NW) = XMLOW(M,NW)*WTMO/(XMLOW(M,NW)*WTMO +
     &        XMLWWX*WTMW)
            IF( XLOW(M,NW).LT.EPSL ) XLOW(M,NW) = 0.D+0
            XLWW(M,NW) = MAX( 1.D+0-XLOW(M,NW),0.D+0 )
            IF( XLWW(M,NW).LT.EPSL ) XLWW(M,NW) = 0.D+0
!
!---        Compute the liquid water molar density and viscosity  ---
!
            RHOMLX = RHOLW(M,NW)/(XMLWWX*WTMW + XMLOW(M,NW)*WTMO)
            CALL WATLQV( T(2,N),PWPX,PSW(2,N),VISLW(M,NW) )
!
!---        Compute the NAPL viscosity  ---
!
            CALL VOCLQV( T(2,N),VISNW(M,NW),VISCO,IVISO )
!
!---        Compute liquid diffusion coefficients  ---
!
            CALL VOCDFL( T(2,N),VISLW(M,NW),DFLOW(M,NW) )
!
!---        Bulk fluid density  ---
!
            RHOGW(M,NW) = ((SLW(M,NW)+EPSL)*RHOLW(M,NW)+
     &        SNW(M,NW)*RHONW(M,NW))/(SLW(M,NW)+SNW(M,NW)+EPSL)
!
!---        Compute well interfacial node pressure  ---
!
            IF( M.EQ.2 ) THEN
              DPWX = STW(M,NW)*DZGF(N)*GRAVZ*RHOGW(M,NW)
            ENDIF
  200     CONTINUE
          PWX = PWX + DPWX
  250   CONTINUE
  300 CONTINUE
!
!---  Compute the well total-liquid depth  ---
!
      DO 420 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        WWDX(2,NWL) = 0.D+0
        WWDX(IEQW+2,NWL) = 0.D+0
        WWDX(IEQO+2,NWL) = 0.D+0
        DO 410 K = KSX,KEX
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          WWDX(2,NWL) = WWDX(2,NWL) + STW(2,NW)*DZGF(N)
          WWDX(IEQW+2,NWL) = WWDX(IEQW+2,NWL) + SLW(2,NW)*DZGF(N)
          WWDX(IEQO+2,NWL) = WWDX(IEQO+2,NWL) + SNW(2,NW)*DZGF(N)
  410   CONTINUE
  420 CONTINUE
!
!---  Load well variables into field variables for output  ---
!
      DO 520 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        DO 510 K = KSX,KEX
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          PI(2,N) = WWDX(2,NWL)
          PI(IEQW+2,N) = WWDX(IEQW+2,NWL)
          PI(IEQO+2,N) = WWDX(IEQO+2,NWL)
  510   CONTINUE
  520 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROPW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PRESSW4( SLWX,SNWX,TWX,PWTX,PWX,PWBX,NW )
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
!     Water-Oil Mode
!     Newton loop to compute the well pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 18 October 2001.
!     Last Modified by M.D. White, PNNL, 18 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
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
      REAL*8 FX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PRESSW4'
!
!---  Newton loop to compute the well pressure  ---
!
      N = IWN(NW)
      K = KD(N)
      STWX = MIN( MAX( SLWX+SNWX,0.D+0 ),1.D+0 )
      CALL WATLQD( TWX,PWTX,RHOLWX )
      CALL VOCLQD( TWX,PWTX,PSO(2,N),RHONWX )
      RHOWX = ((SLWX+EPSL)*RHOLWX + SNWX*RHONWX)/(SLWX + SNWX + EPSL)
      PWGX = PWTX - PATM
      PWX = PWGX + (STWX-5.D-1)*DZGF(N)*GRAVZ*RHOWX
      IF( RHOWX/EPSL.LT.EPSL ) GOTO 30
      NC = 0
   10 CONTINUE
      DPWPX = 1.D-1
      DO 20 L = 1,2
        PWPX = PWX
        IF( L.EQ.2 ) PWPX = PWX + DPWPX
        PWAX = PWPX + PATM
        CALL WATLQD( TWX,PWAX,RHOLWX )
        CALL VOCLQD( TWX,PWAX,PSO(2,N),RHONWX )
        RHOWX = ((SLWX+EPSL)*RHOLWX + SNWX*RHONWX)/
     &    (SLWX + SNWX + EPSL)
        FX(L) = STWX - 5.D-1 - (PWPX-PWGX)/
     &    (DZGF(N)*GRAVZ*RHOWX+EPSL)
   20 CONTINUE
      DFX = (FX(2)-FX(1))/DPWPX
      DPWPX = -FX(1)/(DFX+EPSL)
      PWX = PWX + DPWPX
      NC = NC + 1
      IF( NC.GT.32 ) THEN
        INDX = 15
        IMSG = NW
        CHMSG = 'Convergence Failure on Well Pressure'
     &    // ' at Well Node: '
        CALL WRMSGS( INDX )
      ENDIF
      IF( ABS(DPWPX).GT.1.D-6 ) GOTO 10
   30 CONTINUE
      PWBX = PWTX + STWX*DZGF(N)*GRAVZ*RHOWX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PRESSW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDWELL4
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
!     Water-Oil Mode
!     Reads the well card.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2001.
!     Last Modified by M.D. White, PNNL, 5 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE GRID
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
      INTEGER NCH
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDWELL4'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Well Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Check well parameter switch  ---
!
      IF( LWELL.NE.1 ) THEN
        INDX = 5
        CHMSG = 'Well Parameter Turned Off: Parameter LWELL /= 1'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Read number of wells  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Wells'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NWLS)
!
!---  Initialize well node and well surface counters  ---
!
      NWLN = 0
      NSZX = 0
      DO 1000 NWL = 1,NWLS
!
!---    Read well type  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Well Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        WRITE(IWR,'(/,2A,$)') VARB(1:IVR),': '
        IF( INDEX(ADUM(1:),'monitoring').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Monitoring Well'
          IWT(NWL) = 1
        ELSEIF( INDEX(ADUM(1:),'slurping').NE.0 .AND.
     &    INDEX(ADUM(1:),'aqueous').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Aqueous Slurping Well'
          IWT(NWL) = 2
        ELSEIF( INDEX(ADUM(1:),'slurping').NE.0 .AND.
     &    INDEX(ADUM(1:),'napl').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'NAPL Slurping Well'
          IWT(NWL) = 3
        ELSEIF( INDEX(ADUM(1:),'pumping').NE.0 ) THEN
          WRITE(IWR,'(2X,A)') 'Pumping Well'
          IWT(NWL) = 4
        ENDIF
!
!---    Read well-bore radius  ---
!
        VARB = 'Well-Bore Radius'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WBR(NWL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',WBR(NWL)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,WBR(NWL),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',WBR(NWL),', m)'
!
!---    Read initial well-aqueous depth  ---
!
        VARB = 'Initial Well-Aqueous Depth'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WWDL(NWL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',WWDL(NWL)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,WWDL(NWL),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',WWDL(NWL),', m)'
!
!---    Read initial well-NAPL depth  ---
!
        VARB = 'Initial Well-NAPL Depth'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WWDN(NWL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',WWDN(NWL)
        INDX = 0
        IUNM = 1
        CALL RDUNIT(UNTS,WWDN(NWL),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',WWDN(NWL),', m)'
!
!---    Write initial well-total liquid depth  ---
!
        WWD(NWL) = WWDL(NWL) + WWDN(NWL)
        WRITE(IWR,'(A,1PE11.4,A)') 'Initial Well-Total-Liquid Depth: ',
     &    WWD(NWL),', m'
!
!---    Read initial well-head pressure  ---
!
        VARB = 'Initial Well-Head Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WHP(NWL))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',WHP(NWL)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,WHP(NWL),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',WHP(NWL),', Pa)'
!
!---    Read initial well dissolved-oil saturation  ---
!
        VARB = 'Initial Well Dissolved-Oil Saturation'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WIDO(NWL))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',WIDO(NWL)
!
!---    Read well-bore radius for storage  ---
!
        WBRS(NWL) = WBR(NWL)
        CALL CHKDPR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          VARB = 'Well-Bore Radius for Storage'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,WBRS(NWL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',WBRS(NWL)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,WBRS(NWL),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',WBRS(NWL),', m)'
        ENDIF
!
!---    Read well domain  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Well I Index'
        CALL RDINT( ISTART,ICOMMA,CHDUM,I )
        VARB = 'Well J Index'
        CALL RDINT( ISTART,ICOMMA,CHDUM,J )
        VARB = 'Well Lower K Index'
        CALL RDINT( ISTART,ICOMMA,CHDUM,KSX )
        VARB = 'Well Upper K Index'
        CALL RDINT( ISTART,ICOMMA,CHDUM,KEX )
        I = MAX( 1,I )
        I = MIN( I,IFLD )
        J = MAX( 1,J )
        J = MIN( J,JFLD )
        KSX = MAX( 1,KSX )
        KSX = MIN( KSX,KEX,KFLD )
        KEX = MAX( 1,KSX,KEX )
        KEX = MIN( KEX,KFLD )
        IWLDM(1,NWL) = I
        IWLDM(2,NWL) = J
        IWLDM(3,NWL) = KSX
        IWLDM(4,NWL) = KEX
        WRITE(IWR,'(/,A)' ) 'Well Domain'
        WRITE (IWR,'(2X,A,I6)') 'I = ',I
        WRITE (IWR,'(2X,A,I6)') 'J = ',J
        WRITE (IWR,'(2(2X,A,I6))') 'K = ',KSX,' to ',KEX
        DO 100 K = KSX,KEX
          N = ND(I,J,K)
          NWLN = NWLN+1
          IXW(N) = NWLN
          IWL(NWLN) = NWL
          IWN(NWLN) = N
          NSZX = NSZX+1
          NSZW(NWLN) = NSZX
  100   CONTINUE
        NSZX = NSZX+1
!
!---    Read number of well screen intervals  ---
!
        VARB = 'Number of Well Screen Intervals'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NWSI)
        IF( NWSI.GT.LNWS ) THEN
          INDX = 5
          CHMSG = 'Number of Well Screen Intervals > Parameter LNWS'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read well screen intervals  ---
!
        DO 110 NWS = 1,NWSI
          VARB = 'Screen Lower K Index'
          CALL RDINT( ISTART,ICOMMA,CHDUM,KSX )
          VARB = 'Screen Upper K Index'
          CALL RDINT( ISTART,ICOMMA,CHDUM,KEX )
          KSX = MAX( KSX,KSX )
          KSX = MIN( KEX,KSX )
          KEX = MAX( KSX,KEX )
          KEX = MIN( KEX,KEX )
          IWLDM(3+NWS*2,NWL) = KSX
          IWLDM(4+NWS*2,NWL) = KEX
          WRITE(IWR,'(/,A)' ) 'Well-Screen Interval'
          WRITE (IWR,'(2(2X,A,I6))') 'K = ',KSX,' to ',KEX
          DO 110 K = KSX,KEX
            N = ND(I,J,K)
            IXW(N) = -IXW(N)
  110     CONTINUE
!
!---    Read number of well times  ---
!
        VARB = 'Number of Well Times'
        CALL RDINT(ISTART,ICOMMA,CHDUM,IWM(NWL))
        IF( IWM(NWL).LE.-3 ) THEN
          IWCC(NWL) = 1
          IWM(NWL) = ABS(IWM(NWL))
          WRITE(IWR,'(A)') 'Cyclic Well Times'
        ELSEIF( IWM(NWL).GE.1 ) THEN
          IWCC(NWL) = 0
          WRITE(IWR,'(A)') 'Noncyclic Well Times'
        ELSE
          INDX = 4
          CHMSG = 'Number of Cyclic Well Times < 3'
          CALL WRMSGS( INDX )
        ENDIF
        IF( IWM(NWL).GT.LNWT ) THEN
          INDX = 5
          CHMSG = 'Number of Well Times > Parameter LNWT'
          CALL WRMSGS( INDX )
        ENDIF
        DO 200 M = 1,IWM(NWL)
!
!---      Read new line  ---
!
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
!
!---      Read well time  ---
!
          VARB = 'Well Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(1,M,NWL))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(/,4A,1PE11.4)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',WLVR(1,M,NWL)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,WLVR(1,M,NWL),INDX)
!
!---      Monitoring well variables  ---
!
          IF( IWT(NWL).EQ.1 ) THEN
!
!---        Well-head pressure  ---
!
            VARB = 'Monitoring Well: Well-Head Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(2,M,NWL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',WLVR(2,M,NWL)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,WLVR(2,M,NWL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',WLVR(2,M,NWL),', Pa)'
!
!---      Pumping well variables  ---
!
          ELSEIF( IWT(NWL).EQ.4 ) THEN
!
!---        Well-head pressure  ---
!
            VARB = 'Pumping Well: Well-Head Pressure'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(2,M,NWL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',WLVR(2,M,NWL)
            INDX = 0
            IUNM = -1
            IUNKG = 1
            IUNS = -2
            CALL RDUNIT(UNTS,WLVR(2,M,NWL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',WLVR(2,M,NWL),', Pa)'
!
!---        Injection/withdrawl elevation above well bottom  ---
!
            VARB = 'Pumping Well: Injection/Withdrawl Elevation '
            CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(3,M,NWL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',WLVR(3,M,NWL)
            INDX = 0
            IUNM = 1
            CALL RDUNIT(UNTS,WLVR(3,M,NWL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',WLVR(3,M,NWL),', m)'
            ZBX = 5.D-1*(ZE(1,ND(I,J,KSX))+ZE(4,ND(I,J,KSX)))
            ZTX = 5.D-1*(ZE(5,ND(I,J,KEX))+ZE(8,ND(I,J,KSX)))
            ZWX = ZTX-ZBX
            IF( WLVR(3,M,NWL).LT.0.D+0 .OR.  WLVR(3,M,NWL).GT.ZWX ) THEN
              INDX = 7
              IMSG = NWL
              CHMSG = 'Injection/Withdrawl Elevation Outside Well: '
              CALL WRMSGS( INDX )
            ENDIF
!
!---        Volumetric injection/withdrawl pumping rate;
!           where, pumping rates are defined positive into the well  ---
!
            VARB = 'Pumping Well: Volumetric Pumping Rate'
            CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(4,M,NWL))
            CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
            WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',WLVR(4,M,NWL)
            INDX = 0
            IUNM = 3
            IUNS = -1
            CALL RDUNIT(UNTS,WLVR(4,M,NWL),INDX)
            WRITE(IWR,'(A,1PE11.4,A)') ' (',WLVR(4,M,NWL),', m^3/s)'
            IF( WLVR(4,M,NWL).GT.EPSL ) THEN
!
!---          Injection aqueous fraction  ---
!
              VARB = 'Pumping Well: Injection Aqueous Volumetric Frac.'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(5,M,NWL))
              WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',WLVR(5,M,NWL)
!
!---          Injection dissolved-oil saturation  ---
!
              VARB = 'Pumping Well: Injection Dissolved-Oil Saturation'
              CALL RDDPR(ISTART,ICOMMA,CHDUM,WLVR(6,M,NWL))
              WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',WLVR(6,M,NWL)
            ENDIF
          ENDIF
  200   CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDWELL4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDLW4
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
!     Water-Oil Mode
!     Compute the maximum well relative residuals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 5 October 2001.
!     Last Modified by M.D. White, PNNL, 5 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE OUTPU
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RSDLW4'
!
!---  Reset primary well variables for convergence failure
!     in field nodes  ---
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) THEN
        DO 20 NWL = 1,NWLS
          I = IWLDM(1,NWL)
          J = IWLDM(2,NWL)
          KSX = IWLDM(3,NWL)
          KEX = IWLDM(4,NWL)
          DO 10 K = KSX,KEX
            N = ND(I,J,K)
            NW = ABS(IXW(N))
            PW(2,NW) = PW(1,NW)
            PWLW(2,NW) = PWLW(1,NW)
            YMLOW(2,NW) = YMLOW(1,NW)
            SNW(2,NW) = SNW(1,NW)
            SLW(2,NW) = SLW(1,NW)
            NPHAZW(2,NW) = NPHAZW(1,NW)
   10     CONTINUE
   20   CONTINUE
        GOTO 1000
      ENDIF
!
!---  Initialize well residuals  ---
!
      DO 30 M = 1,ISVC
        RSD(M+ISVC) = 0.D+0
        NSD(M+ISVC) = 0
   30 CONTINUE
!
!---  Loop over well nodes  ---
!
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
      DO 200 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        RDW = WBR(NWL)
        RDWS = WBRS(NWL)
        ARWX = GPI*(RDW**2)
        ARWSX = GPI*(RDWS**2)
        DO 100 K = KSX,KEX
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          NMD = IXP(N)
          MPL = IM(IEQWW,NMD)
          MPN = IM(IEQOW,NMD)
          DPWX = BLU(MPL)
          DPOX = BLU(MPN)
          VOLX = ARWSX*DZGF(N)
!
!---      Aqueous or aqueous-gas w/ dissolved-oil conditions
!         Water equation: well pressure
!         Oil equation:  mole fraction of dissolved oil  ---
!
          IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLWW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPWX)/(ABS(PW(2,NW))+PATM),
     &        ABS(RSDL(IEQWW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQWW) ) THEN
              RSD(IEQWW) = RSDX
              NSD(IEQWW) = NW
            ENDIF
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLOW(2,NW) +
     &        RHONW(2,NW)*SNW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPOX)/(PSO(2,N)/HCOW),
     &        ABS(RSDL(IEQOW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQOW) ) THEN
              RSD(IEQOW) = RSDX
              NSD(IEQOW) = NW
            ENDIF
!
!---      Desaturated conditions
!         Water equation: aqueous saturation
!         Oil equation:  NAPL saturation  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLWW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPWX),
     &        ABS(RSDL(IEQWW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQWW) ) THEN
              RSD(IEQWW) = RSDX
              NSD(IEQWW) = NW
            ENDIF
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLOW(2,NW) +
     &        RHONW(2,NW)*SNW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPOX),
     &        ABS(RSDL(IEQOW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQOW) ) THEN
              RSD(IEQOW) = RSDX
              NSD(IEQOW) = NW
            ENDIF
!
!---      Aqueous-NAPL or aqueous-NAPL-gas
!         Water equation: ratio aqueous/total-liquid saturation
!         Oil equation: well pressure  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLWW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPWX),
     &        ABS(RSDL(IEQWW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQWW) ) THEN
              RSD(IEQWW) = RSDX
              NSD(IEQWW) = NW
            ENDIF
            ACP = (RHOLW(2,NW)*SLW(2,NW)*XLOW(2,NW) +
     &        RHONW(2,NW)*SNW(2,NW))*DTI*VOLX
            RSDX = MIN( ABS(DPOX)/(ABS(PW(2,NW))+PATM),
     &        ABS(RSDL(IEQOW,N)/(ACP+EPSL)) )
            IF( RSDX.GT.RSD(IEQOW) ) THEN
              RSD(IEQOW) = RSDX
              NSD(IEQOW) = NW
            ENDIF
          ENDIF
  100   CONTINUE
  200 CONTINUE
!
!---  Assign a convergence index  ---
!
      DO 300 M = 1,ISVC
        IF( RSD(M+ISVC).GT.RSDMX ) ICNV = 2
  300 CONTINUE
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        WRITE(ISC,'(10X,A)') '---  Convergence Failure  ---'
        WRITE(IWR,'(10X,A)') '---  Convergence Failure  ---'
        WRITE(ISC,'(4X,A,1PE11.4,A,I6)') 'Well Water Equation Maximum '
     &  // 'Residual = ',RSD(IEQWW),' Well Node = ',NSD(IEQWW)
        WRITE(IWR,'(4X,A,1PE11.4,A,I6)') 'Well Water Equation Maximum '
     &  // 'Residual = ',RSD(IEQWW),' Well Node = ',NSD(IEQWW)
        WRITE(ISC,'(4X,A,1PE11.4,A,I6)') 'Well Oil Equation Maximum '
     &  // 'Residual = ',RSD(IEQOW),' Well Node = ',NSD(IEQOW)
        WRITE(IWR,'(4X,A,1PE11.4,A,I6)') 'Well Oil Equation Maximum '
     &  // 'Residual = ',RSD(IEQOW),' Well Node = ',NSD(IEQOW)
!
!---  Reduce time step  ---
!
        IF( NTSR.LT.4 ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - 8.D-1*DT
          DT = 2.D-1*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---      Reset primary field variables  ---
!
          DO 400 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---      Reset primary well variables  ---
!
          DO 420 NWL = 1,NWLS
            I = IWLDM(1,NWL)
            J = IWLDM(2,NWL)
            KSX = IWLDM(3,NWL)
            KEX = IWLDM(4,NWL)
            DO 410 K = KSX,KEX
              N = ND(I,J,K)
              NW = ABS(IXW(N))
              PW(2,NW) = PW(1,NW)
              PWLW(2,NW) = PWLW(1,NW)
              YMLOW(2,NW) = YMLOW(1,NW)
              SNW(2,NW) = SNW(1,NW)
              SLW(2,NW) = SLW(1,NW)
              NPHAZW(2,NW) = NPHAZW(1,NW)
  410       CONTINUE
  420     CONTINUE
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Reset primary field variables  ---
!
          DO 500 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  500     CONTINUE
!
!---      Reset primary well variables  ---
!
          DO 520 NWL = 1,NWLS
            I = IWLDM(1,NWL)
            J = IWLDM(2,NWL)
            KSX = IWLDM(3,NWL)
            KEX = IWLDM(4,NWL)
            DO 510 K = KSX,KEX
              N = ND(I,J,K)
              NW = ABS(IXW(N))
              PW(2,NW) = PW(1,NW)
              PWLW(2,NW) = PWLW(1,NW)
              YMLOW(2,NW) = YMLOW(1,NW)
              SNW(2,NW) = SNW(1,NW)
              SLW(2,NW) = SLW(1,NW)
              NPHAZW(2,NW) = NPHAZW(1,NW)
  510       CONTINUE
  520     CONTINUE
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit ' //
     &      'Exceeded  ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit ' //
     &      'Exceeded  ---'
          ICNV = 4
        ENDIF
      ENDIF
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDLW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SFINW4
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
!     Water-Oil Mode
!     Integrate multifluid well pumping rates.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 2 November 2001.
!     Last Modified by M.D. White, PNNL, 2 November 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_CL
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SFINW4'
!
!---  Loop over wells  ---
!
      DO 10 NWL = 1,NWLS
        QL_W(4,NWL) = QL_W(4,NWL) + QL_W(3,NWL)*DT
        QN_W(4,NWL) = QN_W(4,NWL) + QN_W(3,NWL)*DT
        QT_W(4,NWL) = QT_W(4,NWL) + QT_W(3,NWL)*DT
   10 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SFINW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORCW4
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
!     Water-Oil Mode
!     Well sources.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 1 November 2001.
!     Last Modified by M.D. White, PNNL, 1 November 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE NAPL
      USE JACOB
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
      REAL*8 WLVRX(LNWV)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORCW4'
!
!---  Loop over wells  ---
!
      DO 1000 NWL = 1,NWLS
        IF( IWT(NWL).EQ.1 ) GOTO 1000
        TMZ = TM
        IF( IWM(NWL).EQ.0 ) GOTO 1000
        IF( IWCC(NWL).EQ.1 ) TMZ = MOD( TM,WLVR(1,IWM(NWL),NWL) )
        IF( TMZ.LE.WLVR(1,1,NWL) ) GOTO 1000
        IF( IWM(NWL).GT.1 .AND. TMZ.GT.WLVR(1,IWM(NWL),NWL) ) GOTO 1000
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
!
!---    Assign local well source variables  ---
!
        IF( IWM(NWL).EQ.1 ) THEN
          DO 10 M = 1,LNWV
            WLVRX(M) = WLVR(M,1,NWL)
   10     CONTINUE
        ELSE
          DO 30 L = 2,IWM(NWL)
            IF( TMZ.LE.WLVR(1,L,NWL) ) THEN
              TDBC = (WLVR(1,L,NWL)-WLVR(1,L-1,NWL))
              TFBC = (TMZ-WLVR(1,L-1,NWL))/TDBC
              DO 20 M = 1,LNWV
                WLVRX(M) = WLVR(M,L-1,NWL) +
     &           TFBC*(WLVR(M,L,NWL)-WLVR(M,L-1,NWL))
   20         CONTINUE
              GOTO 40
            ENDIF
   30     CONTINUE
          GOTO 1000
        ENDIF
   40   CONTINUE
!
!---  Loop over well nodes  ---
!
        ZBX = 0.D+0
        SRCFX = 1.D+0
        DO 900 K = KSX,KEX
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          ZTX = ZBX + DZGF(N)
!
!---      Zero sources  ---
!
          DO 50 M = 2,ISVC+2
            SRCW_W(M,NW) = 0.D+0
            SRCO_W(M,NW) = 0.D+0
   50     CONTINUE
!
!---      Multifluid pumping well  ---
!
          IF( IWT(NWL).EQ.4 ) THEN
!
!---        Injection pumping  ---
!
            IF( WLVRX(4).GT.EPSL ) THEN
!
!---          Check elevation point against node elevations  ---
!
              IF( WLVRX(3)-ZBX.GE.-EPSL .AND.
     &          WLVRX(3)-ZTX.LE.EPSL ) THEN
!
!---            Check elevation point near upper surface and
!               for elevation points within "epsilon" of
!               the inter-node surface inject into both nodes  ---
!
                IF( K.LT.KEX ) THEN
                  IF( WLVRX(3)-ZTX.GE.-EPSL .AND.
     &              WLVRX(3)-ZTX.LE.EPSL ) SRCFX = 5.D-1
                ENDIF
!
!---          Define injection sources  ---
!
                DO 400 M = 2,ISVC+2
                  XMLOWX = WLVRX(6)*PSO(2,N)/HCOW
                  XMLWWX = MAX( 1.D+0-XMLOWX,0.D+0 )
                  XLOWX = XMLOWX*WTMO/(XMLWWX*WTMW + XMLOWX*WTMO)
                  XLWWX = MAX( 1.D+0-XLOWX,0.D+0 )
                  SRCW_W(M,NW) = WLVRX(4)*SRCFX*
     &              WLVRX(5)*RHOLW(2,NW)*XLWWX
                  SRCO_W(M,NW) = WLVRX(4)*SRCFX*
     &              (WLVRX(5)*RHOLW(2,NW)*XLOWX +
     &              (1.D+0-WLVRX(5))*RHONW(2,NW))
                  IF( M.EQ.2 ) THEN
                    QL_W(3,NWL) = WLVRX(4)*WLVRX(5)
                    QN_W(3,NWL) = WLVRX(4)*(1.D+0-WLVRX(5))
                    QT_W(3,NWL) = WLVRX(4)
                  ENDIF
  400           CONTINUE
              ENDIF
!
!---        Withdrawl pumping  ---
!
            ELSEIF( WLVRX(4).LT.-EPSL ) THEN
!
!---          Check elevation point against node elevations  ---
!
              IF( WLVRX(3).GE.ZBX .AND. WLVRX(3).LT.ZTX ) THEN
!
!---            Define withdrawl sources,
!               limit withdrawl volumetric rate according to
!               the total liquid volume within the well casing
!               above the pump and the net total-liquid into
!               the well  ---
!
                DO 410 M = 2,ISVC+2
                  ZTLX = STW(M,NW)*DZGF(N) + ZBX
                  ARWSX = GPI*(WBRS(NWL)**2)
                  DZTLX = MAX( PI(2,N)-WLVRX(3),0.D+0 )
                  IF( WLVRX(3)-ZTLX.GT.EPSL ) DZTLX = 0.D+0
                  QTLX = -(ARWSX*DZTLX*DTI + QT_W(2,NWL))
                  QTLX = MAX( WLVRX(4),MIN( QTLX,0.D+0 ) )
                  SLRX = SLW(M,NW)/STW(M,NW)
                  SNRX = SNW(M,NW)/STW(M,NW)
                  SRCW_W(M,NW) = QTLX*SLRX*XLWW(M,NW)*RHOLW(M,NW)
                  SRCO_W(M,NW) = QTLX*(SLRX*XLOW(M,NW)*RHOLW(M,NW) +
     &              SNRX*RHONW(M,NW))
                  IF( M.EQ.2 ) THEN
                    QL_W(3,NWL) = SLRX*QTLX
                    QN_W(3,NWL) = SNRX*QTLX
                    QT_W(3,NWL) = QTLX
                  ENDIF
  410           CONTINUE
              ENDIF
            ENDIF
          ENDIF
          ZBX = ZTX
  900   CONTINUE
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORCW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE STW4( DZX,PWX,PWBX,PWLWX,RHOLWX,RHONWX,SLWX,SNWX,STWX )
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
!     Water-Oil Mode
!     Total-liquid, aqueous, and NAPL saturations from ratio
!     aqueous/total-liquid saturation and well pressure
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 29 October 2001.
!     Last Modified by M.D. White, PNNL, 29 October 2001.
!     Last Modified by M.D. White, PNNL, 30 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GRID
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/STW4'
!
!---  Ratio of aqueous to total-liquid saturation held in
!     variable PWLWX  ---
!
      STWX = 5.D-1 + (PWX-PWBX)/
     &  (GRAVZ*DZX*(PWLWX*(RHOLWX-RHONWX)+RHONWX))
      STWX = MIN( 1.D+0,MAX( 0.D+0,STWX ) )
      SLWX = PWLWX*STWX
      SNWX = STWX-SLWX
      ISUB_LOG = ISUB_LOG-1
!
!---  End of STW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SUPSW4( DZX,PWX,PWBX,RHOLWX,RHONWX,SLWX,SNWX,STWX,
     &  XMLOMX,YMLOWX )
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
!     Water-Oil Mode
!     Supersaturated aqueous oil mole fraction converted to NAPL
!     saturation.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 23 October 2001.
!     Last Modified by M.D. White, PNNL, 23 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE GRID
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SUPSW4'
!
!---  Convert mole fractions to mass fractions  ---
!
      YLOWX = (YMLOWX*WTMO)/(YMLOWX*WTMO+(1.D+0-YMLOWX)*WTMW)
      XLOMX = (XMLOMX*WTMO)/(XMLOMX*WTMO+(1.D+0-XMLOMX)*WTMW)
      XLOFX = RHOLWX*(YLOWX-XLOMX)/RHONWX
!
!---  Supersaturated dissolved-oil conditions  ---
!
      IF( XLOFX.GT.EPSL ) THEN
        SLWX = 1.D+0/(1.D+0 + XLOFX)
        SNWX = 1.D+0 - SLWX
        STWX = 5.D-1 + (PWX-PWBX)*(SLWX+SNWX)/
     &    (GRAVZ*DZX*(SLWX*RHOLWX+SNWX*RHONWX))
!
!---    Total-liquid saturated conditions  ---
!
        IF( STWX.GE.1.D+0 ) THEN
          STWX = 1.D+0
!
!---    Total-liquid unsaturated conditions  ---
!
        ELSE
          SLWX = 5.D-1*(GRAVZ*DZX*(RHOLWX+XLOFX*RHONWX) +
     &      2.D+0*(PWX-PWBX)*(1.D+0+XLOFX))/(GRAVZ*DZX*
     &      (1.D+0+XLOFX)*(RHOLWX+XLOFX*RHONWX))
          SNWX = SLWX*XLOFX
          STWX = SLWX + SNWX
        ENDIF
!
!---  Subsaturated dissolved-oil conditions  ---
!
      ELSE
        SNWX = 0.D+0
        SLWX = 5.D-1 + (PWX-PWBX)/(GRAVZ*DZX*RHOLWX)
        SLWX = MIN( MAX( SLWX,0.D+0 ),1.D+0 )
        STWX = SLWX
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SUPSW4 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTW4
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
!     Water-Oil Mode
!     Update well variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by M.D. White, PNNL, 10 October 2001.
!     Last Modified by M.D. White, PNNL, 10 October 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FD
      USE WELL_CL
      USE SOLTN
      USE OUTPU
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FILES
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/UPDTW4'
!
!---  Reset primary well variables for excessive primary variable
!     error in field nodes  ---
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) THEN
        DO 20 NWL = 1,NWLS
          I = IWLDM(1,NWL)
          J = IWLDM(2,NWL)
          KSX = IWLDM(3,NWL)
          KEX = IWLDM(4,NWL)
          DO 10 K = KSX,KEX
            N = ND(I,J,K)
            NW = ABS(IXW(N))
            PW(2,NW) = PW(1,NW)
            PWLW(2,NW) = PWLW(1,NW)
            YMLOW(2,NW) = YMLOW(1,NW)
            SNW(2,NW) = SNW(1,NW)
            SLW(2,NW) = SLW(1,NW)
            NPHAZW(2,NW) = NPHAZW(1,NW)
   10     CONTINUE
   20   CONTINUE
        GOTO 1000
      ENDIF
      IEQWW = IEQW+ISVC
      IEQOW = IEQO+ISVC
!
!---  Loop over well nodes  ---
!
      DO 200 NWL = 1,NWLS
        I = IWLDM(1,NWL)
        J = IWLDM(2,NWL)
        KSX = IWLDM(3,NWL)
        KEX = IWLDM(4,NWL)
        DO 100 K = KSX,KEX
          N = ND(I,J,K)
          NW = ABS(IXW(N))
          NMD = IXP(N)
          MPL = IM(IEQWW,NMD)
          MPN = IM(IEQOW,NMD)
          DPWX = BLU(MPL)
          DPOX = BLU(MPN)
!
!---      Aqueous or aqueous-gas w/ dissolved-oil
!         Water equation: well pressure
!         Oil equation:  mole fraction of dissolved oil  ---
!
          IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
!
!---        Limit changes in well pressure to 1.D+6 Pa  ---
!
            IF( ABS(DPWX).GT.1.D+6 ) IERR = 1
            PW(2,NW) = PW(2,NW) + DPWX
!
!---        Limit changes in mole fraction of dissolved oil to 25% of
!           solubility limit  ---
!
            XMLOMX = 2.5D-1*PSO(2,N)/HCOW
            DPOX = SIGN( MIN( XMLOMX,ABS(DPOX) ),DPOX )
            YMLOW(2,NW) = YMLOW(2,NW) + DPOX
            IF( YMLOW(2,NW).LT.EPSL ) YMLOW(2,NW) = 0.D+0
!
!---      Desaturated conditions
!         Water equation: aqueous saturation
!         Oil equation:  NAPL saturation  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
!
!---        Limit changes in saturation to 0.125  ---
!
            DPWX = SIGN( MIN( 1.25D-1,ABS(DPWX) ),DPWX )
            SLW(2,NW) = SLW(2,NW) + DPWX
            IF( SLW(2,NW).LT.EPSL ) SLW(2,NW) = 0.D+0
            SLW(2,NW) = MIN( SLW(2,NW),1.D+0 )
            DPOX = SIGN( MIN( 1.25D-1,ABS(DPOX) ),DPOX )
            SNW(2,NW) = SNW(2,NW) + DPOX
            IF( SNW(2,NW).LT.EPSL ) SNW(2,NW) = 0.D+0
            SNW(2,NW) = MIN( SNW(2,NW),1.D+0 )
!
!---      Aqueous-NAPL or aqueous-NAPL-gas conditions
!         Water equation: ratio aqueous over total-liquid saturation
!         Oil equation:  well pressure  ---
!
          ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
            PWLW(2,NW) = PWLW(2,NW) + DPWX
            IF( PWLW(2,NW).LT.EPSL ) PWLW(2,NW) = 0.D+0
            PWLW(2,NW) = MIN( PWLW(2,NW),1.D+0 )
!
!---        Limit changes in well pressure to 1.D+6 Pa  ---
!
            IF( ABS(DPOX).GT.1.D+6 ) IERR = 1
            PW(2,NW) = PW(2,NW) + DPOX
          ENDIF
!
!---      Reduce time step for excessive changes in
!         well primary variables   ---
!
          IF( IERR.EQ.1 ) THEN
            ICNV = 1
            WRITE(ISC,'(10X,A)') '---  Excessive Primary Variable ' //
     &        'Change  ---'
            WRITE(IWR,'(10X,A)') '---  Excessive Primary Variable ' //
     &        'Change  ---'
            IF( NPHAZW(2,NW).EQ.1 .OR. NPHAZW(2,NW).EQ.2 ) THEN
              WRITE(ISC,'(4X,A)') 'Aqueous-Saturated Well Node'
              WRITE(IWR,'(4X,A)') 'Aqueous-Saturated Well Node'
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Pressure = ',PW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Pressure = ',PW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &          'Well Dissolved-Oil Mole Fraction = ',YMLOW(2,NW),
     &          ' Well No. = ',NWL,' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &          'Well Dissolved-Oil Mole Fraction = ',YMLOW(2,NW),
     &          ' Well No. = ',NWL,' Well Node = ',NW
            ELSEIF( NPHAZW(2,NW).EQ.3 .OR. NPHAZW(2,NW).EQ.4 ) THEN
              WRITE(ISC,'(4X,A)') 'Desaturated Well Node'
              WRITE(IWR,'(4X,A)') 'Desaturated Well Node'
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Aqueous Saturation = ',SLW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Aqueous Saturation = ',SLW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &          'Well NAPL Saturation = ',SNW(2,NW),
     &          ' Well No. = ',NWL,' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &          'Well NAPL Saturation = ',SNW(2,NW),
     &          ' Well No. = ',NWL,' Well Node = ',NW
            ELSEIF( NPHAZW(2,NW).EQ.5 .OR. NPHAZW(2,NW).EQ.6 ) THEN
              WRITE(ISC,'(4X,A)') 'Aqueous-NAPL-Saturated Well Node'
              WRITE(IWR,'(4X,A)') 'Aqueous-NAPL-Saturated Well Node'
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Ratio Well Aqueous/Total-Liquid Saturation = ',
     &        PWLW(2,NW),' Well No. = ',NWL,' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Ratio Well Aqueous/Total-Liquid Saturation = ',
     &        PWLW(2,NW),' Well No. = ',NWL,' Well Node = ',NW
              WRITE(ISC,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Pressure = ',PW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
              WRITE(IWR,'(4X,A,1PE12.5,A,I6,A,I6)')
     &        'Well Pressure = ',PW(2,NW),' Well No. = ',NWL,
     &        ' Well Node = ',NW
            ENDIF
            GOTO 300
          ENDIF
  100   CONTINUE
  200 CONTINUE
!
!---      Reduce time step  ---
!
  300 CONTINUE
      IF( ICNV.EQ.1 ) THEN
        IF( NTSR.LT.4 ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - 8.D-1*DT
          DT = 2.D-1*DT
          DTO = DT
          DTI = 1.D+0/DT
          VAR = DT
          VARX = DTX
          IF( UNTM.NE.'null' ) THEN
            INDX = 1
            IUNS = 1
            CALL RDUNIT(UNTM,VAR,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,VARX,INDX)
            NCH = INDEX( UNTM,'  ')-1
          ENDIF
          WRITE(ISC,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(4X,A,1PE11.4,1X,2A,1PE11.4,1X,A)')
     &      'Time Step Reduced From ',VARX,UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---      Reset primary field variables  ---
!
          DO 400 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  400     CONTINUE
!
!---      Reset primary well variables  ---
!
          DO 420 NWL = 1,NWLS
            I = IWLDM(1,NWL)
            J = IWLDM(2,NWL)
            KSX = IWLDM(3,NWL)
            KEX = IWLDM(4,NWL)
            DO 410 K = KSX,KEX
              N = ND(I,J,K)
              NW = ABS(IXW(N))
              PW(2,NW) = PW(1,NW)
              PWLW(2,NW) = PWLW(1,NW)
              YMLOW(2,NW) = YMLOW(1,NW)
              SNW(2,NW) = SNW(1,NW)
              SLW(2,NW) = SLW(1,NW)
              NPHAZW(2,NW) = NPHAZW(1,NW)
  410       CONTINUE
  420     CONTINUE
!
!---  Number of time step reductions failure: stop simulation  ---
!
        ELSE
!
!---      Reset primary variables  ---
!
          DO 500 N = 1,NFLD
            PL(2,N) = PL(1,N)
            PN(2,N) = PN(1,N)
            XMLO(2,N) = XMLO(1,N)
            SN(2,N) = SN(1,N)
            NPHAZ(2,N) = NPHAZ(1,N)
  500     CONTINUE
!
!---      Reset primary well variables  ---
!
          DO 520 NWL = 1,NWLS
            I = IWLDM(1,NWL)
            J = IWLDM(2,NWL)
            KSX = IWLDM(3,NWL)
            KEX = IWLDM(4,NWL)
            DO 510 K = KSX,KEX
              N = ND(I,J,K)
              NW = ABS(IXW(N))
              PW(2,NW) = PW(1,NW)
              PWLW(2,NW) = PWLW(1,NW)
              YMLOW(2,NW) = YMLOW(1,NW)
              SNW(2,NW) = SNW(1,NW)
              SLW(2,NW) = SLW(1,NW)
              NPHAZW(2,NW) = NPHAZW(1,NW)
  510       CONTINUE
  520     CONTINUE
          WRITE(ISC,'(10X,A)') '---  Time Step Reduction Limit ' //
     &      'Exceeded  ---'
          WRITE(IWR,'(10X,A)') '---  Time Step Reduction Limit ' //
     &      'Exceeded  ---'
          ICNV = 4
        ENDIF
      ENDIF
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTW4 group  ---
!
      RETURN
      END

