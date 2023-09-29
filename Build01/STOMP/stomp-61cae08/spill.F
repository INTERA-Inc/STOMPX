!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLXSPL
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
!     Aqueous surface spill flux and infiltration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 07 December, 2007.
!     Last Modified by SK Wurstner PNNL, 04 November, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!
!----------------------Type Declarations-------------------------------!
!
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FLXSPL'
!
!---  X-direction aqueous surface spill flux  ----
!
      IF( IFLD.GT.1 ) THEN
        DO 200 J = 1,JFLD
        DO 200 I = 2,IFLD
          NSP = (J-1)*IFLD + I
          NSW = NSP-1
          NSPX = (J-1)*(IFLD-1) + I - 1
          K = KSPS(I,J)
          N = ND(I,J,K)
          NW = N-1
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            THET = TXSP(NSPX)+THXZ
!            HLSPN(MN,NSW) = ((ZSP(NSW)+HLSP(MN,NSW))/COS(THET)) 
!     &        - ZSP(NSW)
!            HLSPN(MP,NSP) = ((ZSP(NSP)+HLSP(MP,NSP))/COS(THET)) 
!     &        - ZSP(NSP)
!            HLSPN(MP,NSP) = ((ZSP(NSP)+HLSP(MP,NSP))/COS(THET)) 
!     &        - ZSP(NSP)
            HLSPNWX = HLSP(MN,NSW)*COS(THET) 
            HLSPNPX = HLSP(MP,NSP)*COS(THET)
            DHDX = (HLSPNWX-HLSPNPX)/PDXSP(NSPX)
            IF( M.EQ.1 ) DHL = DHDX
!
!---        Upwind averaging  ----
!
            INDX = -4
            RHOAVG = DIFMN(RHOLSP(MN,NSW),RHOLSP(MP,NSP),
     &        DXGF(NW),DXGF(N),DHL,INDX)
            VISAVG = DIFMN(VISLSP(MN,NSW),VISLSP(MP,NSP),
     &        DXGF(NW),DXGF(N),DHL,INDX)
            HLSPAVG = DIFMN(HLSPNWX,HLSPNPX,
     &        DXGF(NW),DXGF(N),DHL,INDX)
            R = GRAV*RHOAVG*COS(THET)/3.D+0/VISAVG
            ULSP(M,NSPX) = R*(HLSPAVG**3*DHDX+TAN(THET)*HLSPAVG**3)
!            
            IF( HLSPNWX.LT.EPSL ) 
     &          ULSP(M,NSPX) = MIN(ULSP(M,NSPX),0.D+0)
            IF( HLSPNPX.LT.EPSL ) 
     &          ULSP(M,NSPX) = MAX(ULSP(M,NSPX),0.D+0)
!     
  100     CONTINUE
  200   CONTINUE
      ENDIF      
!
!---  Y-direction aqueous surface spill flux  ---
!
      IF( JFLD.GT.1 ) THEN
        DO 400 J = 2,JFLD
        DO 400 I = 1,IFLD
          NSP = (J-1)*IFLD + I
          NSS = NSP-IFLD
          NSPY = (J-2)*IFLD + I
          K = KSPS(I,J)
          N = ND(I,J,K)
          NS = N-IFLD
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            THET = TYSP(NSPY)+THYZ
!            HLSPN(MN,NSS) = ((ZSP(NSS)+HLSP(MN,NSS))/COS(THET))
!     &        - ZSP(NSS)
!            HLSPN(MP,NSP) = ((ZSP(NSP)+HLSP(MP,NSP))/COS(THET))
!     &        - ZSP(NSP)
            HLSPNSX = HLSP(MN,NSS)*COS(THET) 
            HLSPNPX = HLSP(MP,NSP)*COS(THET)
            DHDY = (HLSPNSX-HLSPNPX)/PDYSP(NSPY)
            IF( M.EQ.1 ) DHL = DHDY
!
!---        Upwind averaging  ----
!
            INDX = -4
            RHOAVG = DIFMN(RHOLSP(MN,NSS),RHOLSP(MP,NSP),
     &        DYGF(NS),DYGF(N),DHL,INDX)
            VISAVG = DIFMN(VISLSP(MN,NSS),VISLSP(MP,NSP),
     &        DYGF(NS),DYGF(N),DHL,INDX)
            HLSPAVG = DIFMN(HLSPNSX,HLSPNPX,
     &        DYGF(NS),DYGF(N),DHL,INDX)
            R = GRAV*RHOAVG*COS(THET)/3.D+0/VISAVG
            VLSP(M,NSPY) =  R*(HLSPAVG**3*DHDY+TAN(THET)*HLSPAVG**3)  
            IF( HLSPNSX.LT.EPSL ) 
     &          VLSP(M,NSPY) = MIN(VLSP(M,NSPY),0.D+0)
            IF( HLSPNPX.LT.EPSL ) 
     &          VLSP(M,NSPY) = MAX(VLSP(M,NSPY),0.D+0)
 !
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Aqueous infiltration  ---
!
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        K = KSPS(I,J)
        N = ND(I,J,K)
        NQZ = NSZ(N)+IJFLD
        IF( IXP(N).GT.0 ) THEN
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDLZ = PL(MN,N) - PL(2,N)
     &        - 5.D-1*GRVZ(NQZ)*DZGF(N)*(RHOL(MN,N)-RHOL(2,N))
     &        - HLSP(MP,NSP)*RHOLSP(MP,NSP)*GRVZ(NQZ)
            IF( HLSP(MP,NSP).LT.EPSL ) HDLZ = MAX( HDLZ,0.D+0 )
            IF( PL(2,N).LT.PG(2,N) ) HDLZ = MIN( HDLZ,0.D+0 )
            IF( M.EQ.1 ) HDZ = HDLZ
            RKLSPX = 1.D+0
            INDX = 10
            RKLZ = DIFMN(RKL(3,MN,N),RKLSPX,DZGF(N),DZGF(N),HDZ,INDX)
            INDX = 7
            VLZ = DIFMN(VISL(MN,N),VISLSP(MP,NSP),DZGF(N),
     &        DZGF(N),HDZ,INDX)
            PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
            WL(M,NQZ) = 2.D+0*PERM_PX*RKLZ*HDLZ/(VLZ*DZGF(N))
  500     CONTINUE
        ENDIF
  600 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLXSPL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FLXSPN
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
!     NAPL surface spill flux and infiltration.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 07 December, 2007.
!     Last Modified by SK Wurstner PNNL, 04 November, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
      USE FDVN
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
      SUB_LOG(ISUB_LOG) = '/FLXSPN'
      HNSPMN = 1.D-4
!
!---  X-direction NAPL surface spill flux  ----
!
!---  Equations are written in terms of model datum so that 
!---  head is expressed normal to the tilted surface ---
!
      IF( IFLD.GT.1 ) THEN
        DO 200 J = 1,JFLD
        DO 200 I = 2,IFLD
          NSP = (J-1)*IFLD + I
          NSW = NSP-1
          NSPX = (J-1)*(IFLD-1) + I - 1
          K = KSPS(I,J)
          N = ND(I,J,K)
          NW = N-1
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            THET = TXSP(NSPX)+THXZ
!            HNSPN(MN,NSW) = ((ZSP(NSW)+HNSP(MN,NSW))/COS(THET)) 
!     &        - ZSP(NSW)
!            HNSPN(MP,NSP) = ((ZSP(NSP)+HNSP(MP,NSP))/COS(THET)) 
!     &        - ZSP(NSP)
            HNSPNWX = HNSP(MN,NSW)*SPNORM(NSP)
            HNSPNPX = HNSP(MP,NSP)*SPNORM(NSP)
!            HNSPN(MN,NSW) = HNSP(MN,NSW)*COS(THET)
!            HNSPN(MP,NSP) = HNSP(MP,NSP)*COS(THET)
            DHDX = (HNSPNWX-HNSPNPX)/PDXSP(NSPX)
            IF( M.EQ.1 ) DHN = DHDX
!
!---        Upwind averaging  ----
!
            INDX = -4
            RHOAVG = DIFMN(RHONSP(MN,NSW),RHONSP(MP,NSP),
     &        DXGF(NW),DXGF(N),DHN,INDX)
            VISAVG = DIFMN(VISNSP(MN,NSW),VISNSP(MP,NSP),
     &        DXGF(NW),DXGF(N),DHN,INDX)
            HNSPAVG = DIFMN(HNSPNWX,HNSPNPX,
     &        DXGF(NW),DXGF(N),DHN,INDX)
            R = GRAV*RHOAVG*COS(THET)/3.D+0/VISAVG
!            GAM = TAN(TXSP(NSPX)+THXZ)
            UNSP(M,NSPX) = R*(HNSPAVG**3*DHDX+TAN(THET)*HNSPAVG**3)
!
!---        Check flux against minimum flux based on mass  ----
! 
            IF( UNSP(M,NSPX).GT.0.D+0 ) THEN
              HNSPX = MAX(HNSP(MN,NSW)-SPHMIN(IZ(NW)),0.D+0)
              USP = HNSPX*AFZSP(NSW)*DTI/DYSP(NSPX)
              UNSP(M,NSPX) = MIN(UNSP(M,NSPX),USP)
            ELSE
              HNSPX = MAX(HNSP(MP,NSP)-SPHMIN(IZ(N)),0.D+0)
              USP = -1.D+0*HNSPX*AFZSP(NSP)*DTI/DYSP(NSPX)
              UNSP(M,NSPX) = MAX(UNSP(M,NSPX),USP)
            ENDIF                  
!
!---        Check spill height against specified minimum height  ----
! 
            IF( HNSPNWX.LT.SPHMIN(IZ(N)) .AND. 
     &        HNSPNPX.LT.SPHMIN(IZ(N)) ) UNSP(M,NSPX) = 0.D+0
!
  100     CONTINUE
  200   CONTINUE
      ENDIF      
!
!---  Y-direction NAPL surface spill flux  ---
!
!      WRITE(110,*) 'I, J, M, NSPY, VNSP(M,NSPY)'
      IF( JFLD.GT.1 ) THEN
        DO 400 J = 2,JFLD
        DO 400 I = 1,IFLD
          NSP = (J-1)*IFLD + I
          NSS = NSP-IFLD
          NSPY = (J-2)*IFLD + I
          K = KSPS(I,J)
          N = ND(I,J,K)
          NS = N-IFLD
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            THET = TYSP(NSPY)+THYZ
!            HNSPN(MN,NSS) = ((ZSP(NSS)+HNSP(MN,NSS))/COS(THET))
!     &        - ZSP(NSS)
!            HNSPN(MP,NSP) = ((ZSP(NSP)+HNSP(MP,NSP))/COS(THET))
!     &        - ZSP(NSP)
            HNSPNSX = HNSP(MN,NSS)*SPNORM(NSP)
            HNSPNPX = HNSP(MP,NSP)*SPNORM(NSP)
!            HNSPN(MN,NSS) = HNSP(MN,NSS)*COS(THET)
!           HNSPN(MP,NSP) = HNSP(MP,NSP)*COS(THET)
            DHDY = (HNSPNSX-HNSPNPX)/PDYSP(NSPY)
            IF( M.EQ.1 ) DHN = DHDY
!
!---        Upwind averaging  ----
!
            INDX = -4
            RHOAVG = DIFMN(RHONSP(MN,NSS),RHONSP(MP,NSP),
     &        DYGF(NS),DYGF(N),DHN,INDX)
            VISAVG = DIFMN(VISNSP(MN,NSS),VISNSP(MP,NSP),
     &        DYGF(NS),DYGF(N),DHN,INDX)
            HNSPAVG = DIFMN(HNSPNSX,HNSPNPX,
     &        DYGF(NS),DYGF(N),DHN,INDX)
            R = GRAV*RHOAVG*COS(THET)/3.D+0/VISAVG
!            GAM = TAN(TYSP(NSPY)+THYZ)
            VNSP(M,NSPY) = R*(HNSPAVG**3*DHDY+TAN(THET)*HNSPAVG**3)
! 
!
!---        Check flux against minimum flux based on mass  ----
!
            IF( VNSP(M,NSPY).GT.0.D+0 ) THEN
              HNSPY = MAX(HNSP(MN,NSS)-SPHMIN(IZ(NS)),0.D+0)
              VSP = HNSPY*AFZSP(NSS)*DTI/DXSP(NSPY)
              VNSP(M,NSPY) = MIN(VNSP(M,NSPY),VSP)
            ELSE
              HNSPY = MAX(HNSP(MP,NSP)-SPHMIN(IZ(N)),0.D+0)
              VSP = -1.D+0*HNSPY*AFZSP(NSP)*DTI/DXSP(NSPY)
              VNSP(M,NSPY) = MAX(VNSP(M,NSPY),VSP)
            ENDIF              
!
!  DEBUG CHECK
!
!            IF( NSPY.EQ.35 .OR. NSPY.EQ.45) then
!           IF( I.EQ.6 .AND. J.EQ.36) then

!              WRITE(111,*) 'VNSP BEFORE CONDITIOPNAL'
!              WRITE(111,*) I, J, M, NSPY, VNSP(M,NSPY)
!              WRITE(111,*) THET,HNSPN(MN,NSS),HNSPN(MP,NSP)      
!            endif
!
!            IF( HNSPN(MN,NSS).LT.EPSL ) 
!     &          VNSP(M,NSPY) = MIN(VNSP(M,NSPY),0.D+0)
!            IF( HNSPN(MP,NSP).LT.EPSL ) 
!     &          VNSP(M,NSPY) = MAX(VNSP(M,NSPY),0.D+0)
!
!---        Check spill height against specified minimum height  ----
! 
            IF( HNSPNSX.LT.SPHMIN(IZ(N)) .AND. 
     &        HNSPNPX.LT.SPHMIN(IZ(N)) ) VNSP(M,NSPY) = 0.D+0
!
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  NAPL Infiltration  ---
!
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        K = KSPS(I,J)
        N = ND(I,J,K)
        NQZ = NSZ(N)+IJFLD
        IF( IXP(N).GT.0 ) THEN
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            HDNZ = PN(MN,N) - PN(2,N)
     &        - 5.D-1*GRVZ(NQZ)*DZGF(N)*(RHON(MN,N)-RHON(2,N))
     &        - HNSP(MP,NSP)*RHONSP(MP,NSP)*GRVZ(NQZ)
            IF( HNSP(MP,NSP)/EPSL.LT.EPSL ) HDNZ = MAX( HDNZ,0.D+0 )
            IF( PN(2,N).LT.PG(2,N) ) HDNZ = MIN( HDNZ,0.D+0 )
            IF( M.EQ.1 ) HDZ = HDNZ
            RKNSPX = 1.D+0
            INDX = 10
            RKNZ = DIFMN(RKN(MN,N),RKNSPX,DZGF(N),DZGF(N),HDZ,INDX)
            INDX = 7
            VNZ = DIFMN(VISN(MN,N),VISNSP(MP,NSP),DZGF(N),
     &        DZGF(N),HDZ,INDX)
            PERM_PX = PERMRF(MN,N)*PERM(3,IZ(N))
            WN(M,NQZ) = 2.D+0*PERM_PX*RKNZ*HDNZ/(VNZ*DZGF(N))
  500     CONTINUE
!
!---      Ignore critical NAPL pressure changes with
!         aqueous pressure changes for no NAPL conditions  ---
!
          IF( SN(2,N).LT.EPSL ) WN(3,NQZ) = WN(1,NQZ)
        ENDIF
  600 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FLXSPN group  ---
!
      RETURN
      END
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INCRMSP
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
!     Increment aqueous and NAPL surface spill heads.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 22 February, 2008.
!     Last Modified by SK Wurstner PNNL, 22 February, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/INCRMSP'
!
!---  Increment aqueous and NAPL surface spill heads  ---
!
      DO 200 J = 1,JFLD
      DO 200 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        DNRSP(IEQW,NSP) = 1.D-8
        DNRSP(IEQO,NSP) = 1.D-8
        IF( HLSP(2,NSP).GT.EPSL ) DNRSP(IEQW,NSP) = 1.D-6*HLSP(2,NSP)
        IF( HNSP(2,NSP).GT.EPSL ) DNRSP(IEQO,NSP) = 1.D-6*HNSP(2,NSP)
        DO 100 M = 2,ISVC+2
          HNSP(M,NSP) = HNSP(2,NSP)
          HLSP(M,NSP) = HLSP(2,NSP)
          IF( M.EQ.IEQW+2 ) THEN
            HLSP(M,NSP) = HLSP(M,NSP) + DNRSP(IEQW,NSP)
          ELSEIF( M.EQ.IEQO+2 ) THEN
            HNSP(M,NSP) = HNSP(M,NSP) + DNRSP(IEQO,NSP)
          ENDIF
  100   CONTINUE
  200 CONTINUE     
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INCRMSP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLDFD( RSS,RSP,RSA,I,J,K,MEQ )
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
!     Load the Jacobian matrix for field nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 March 2008.
!     Last Modified by MD White, PNNL, 11 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,5)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLDFD'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Field node  ---
!
        N = ND(I,J,K)
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        DO 100 M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSP(M)/DNR(M,N)
  100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
!
!---    Surface node  ---
!
        NSP = (J-1)*IFLD + I
        NMD = IXP(NSP+NFLD)
        DO 200 M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,1)/DNRSP(M,NSP)
  200   CONTINUE
!
!---  Splib solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Field node  ---
!
        N = ND(I,J,K)
        NMD = IXP(N)
        MP = IM(MEQ,NMD)
        MA = 3*ISVC
        DO 2100 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + RSP(M)/DNR(M,N)
 2100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
!
!---    Surface node  ---
!
        NSP = (J-1)*IFLD + I
        NMD = IXP(NSP+NFLD)
        MA = 6*ISVC
        DO 2700 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + RSA(M,1)/DNRSP(M,NSP)
 2700   CONTINUE
!
!---  Unrecognized linear equation solver ---
!
      ELSE
        WRITE(ISC,'(A)') 'ERROR: Unrecognized Linear Equation Solver'
        WRITE(IWR,'(A)') 'ERROR: Unrecognized Linear Equation Solver'
        STOP
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLDFD group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBLDSP( RSS,RSP,RSA,I,J,K,MEQ )
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
!     Load the Jacobian matrix for surface spill nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 4 March 2008.
!     Last Modified by MD White, PNNL, 4 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE JACOB
      USE GRID
      USE FILES
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 RSP(LUK),RSA(LUK,5)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBLDSP'
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
!
!---    Surface node  ---
!
        NSP = (J-1)*IFLD + I
        NMD = IXP(NSP+NFLD)
        MP = IM(MEQ,NMD)
        DO 100 M = 1,ISVC
          MCOL = IM(M,NMD)
          MROW = MP-MCOL+MDC
          ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSP(M)/DNRSP(M,NSP)
  100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
!
!---    Field node  ---
!
        N = ND(I,J,K)
        IF( IXP(N).GT.0 ) THEN
          NMD = IXP(N)
          DO 200 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,1)/DNR(M,N)
  200     CONTINUE
        ENDIF
!
!---    South surface node ---
!
        IF( J.NE.1 ) THEN
          NSS = NSP-IFLD
          NMD = IXP(NSS+NFLD)
          DO 300 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,2)/DNRSP(M,NSS)
  300     CONTINUE
        ENDIF
!
!---    West surface node ---
!
        IF( I.NE.1 ) THEN
          NSW = NSP-1
          NMD = IXP(NSW+NFLD)
          DO 400 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,3)/DNRSP(M,NSW)
  400     CONTINUE
        ENDIF
!
!---    East surface node ---
!
        IF( I.NE.IFLD ) THEN
          NSE = NSP+1
          NMD = IXP(NSE+NFLD)
          DO 500 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,4)/DNRSP(M,NSE)
  500     CONTINUE
        ENDIF
!
!---    North surface node ---
!
        IF( J.NE.JFLD ) THEN
          NSN = NSP+IFLD
          NMD = IXP(NSN+NFLD)
          DO 600 M = 1,ISVC
            MCOL = IM(M,NMD)
            MROW = MP-MCOL+MDC
            ALU(MROW,MCOL) = ALU(MROW,MCOL) + RSA(M,5)/DNRSP(M,NSN)
  600     CONTINUE
        ENDIF
!
!---  Splib solver  ---
!
      ELSEIF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Surface node  ---
!
        NSP = (J-1)*IFLD + I
        NMD = IXP(NSP+NFLD)
        MP = IM(MEQ,NMD)
        MA = 3*ISVC
        DO 2100 M = 1,ISVC
          MCOL = KLU(MP,M+MA)
          DLU(MCOL) = DLU(MCOL) + RSP(M)/DNRSP(M,NSP)
 2100   CONTINUE
        BLU(MP) = BLU(MP) - RSS
!        RSDL(MEQ,N) = BLU(MP)
!
!---    Field node  ---
!
        N = ND(I,J,K)
        IF( IXP(N).GT.0 ) THEN
          NMD = IXP(N)
          DO 2200 M = 1,ISVC
            MCOL = KLU(MP,M)
            DLU(MCOL) = DLU(MCOL) + RSA(M,1)/DNR(M,N)
 2200     CONTINUE
        ENDIF
!
!---    South surface node ---
!
        IF( J.NE.1 ) THEN
          NSS = NSP-IFLD
          NMD = IXP(NSS+NFLD)
          MA = ISVC
          DO 2300 M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + RSA(M,2)/DNRSP(M,NSS)
 2300     CONTINUE
        ENDIF
!
!---    West surface node ---
!
        IF( I.NE.1 ) THEN
          NSW = NSP-1
          NMD = IXP(NSW+NFLD)
          MA = 2*ISVC
          DO 2400 M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + RSA(M,3)/DNRSP(M,NSW)
 2400     CONTINUE
        ENDIF
!
!---    East surface node ---
!
        IF( I.NE.IFLD ) THEN
          NSE = NSP+1
          NMD = IXP(NSE+NFLD)
          MA = 4*ISVC
          DO 2500 M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + RSA(M,4)/DNRSP(M,NSE)
 2500     CONTINUE
        ENDIF
!
!---    North surface node ---
!
        IF( J.NE.JFLD ) THEN
          NSN = NSP+IFLD
          NMD = IXP(NSN+NFLD)
          MA = 5*ISVC
          DO 2600 M = 1,ISVC
            MCOL = KLU(MP,M+MA)
            DLU(MCOL) = DLU(MCOL) + RSA(M,5)/DNRSP(M,NSN)
 2600     CONTINUE
        ENDIF
!
!---  Unrecognized linear equation solver ---
!
      ELSE
        WRITE(ISC,'(A)') 'ERROR: Unrecognized Linear Equation Solver'
        WRITE(IWR,'(A)') 'ERROR: Unrecognized Linear Equation Solver'
        STOP
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBLDSP group
!
      RETURN
      END
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBPSP
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
!     Configure the Jacobian matrix pointer arrays.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 March 2008.
!     Last Modified by MD White, PNNL, 6 March 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/JCBPSP'
!
!---  Jacobian matrix space  ---
!
      IF( (NFLD-NXP).GT.LAN ) THEN
        INDX = 12
        CHMSG = 'Grid Refinement: Number of Active Nodes > ' //
     &    'Parameter LAN: '
        IMSG = NFLD-NXP
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Skip to alternate algorithm for well simulations  ---
!
      IF( LWELL.EQ.0 ) THEN
!
!---  Compute the Jacobian matrix pointers  ---
!
      DO 110 K = 1,KFLD
      DO 110 J = 1,JFLD
      DO 110 I = 1,IFLD
        N = ND(I,J,K)
        IF( IXP(N).EQ.0 ) GOTO 110
        NMD = IXP(N)
        DO 100 M = 1,ISVC
          IM(M,NMD) = (NMD-1)*ISVC + M
  100   CONTINUE
  110 CONTINUE
!
!---  Compute the Jacobian matrix pointers for surface spill nodes ---
!
      DO 130 N = 1,IJFLD
        NMD = IXP(NFLD+N)
        DO 120 M = 1,ISVC
          IM(M,NMD) = (NMD-1)*ISVC + M
  120   CONTINUE
  130 CONTINUE
!
!---  Determine the matrix half-band widths ---
!
      MHBW = 0
      DO 530 I = 1,IFLD
        DO 520 J = 1,JFLD
          DO 510 K = 1,KFLD
            N = ND(I,J,K)
            IF( IXP(N).EQ.0 ) GOTO 510
!
!---        Connection to bottom node ---
!
            IF( K.GT.1 ) THEN
              IF( IXP(N-IJFLD).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N-IJFLD)-IXP(N)))
            ENDIF
!
!---        Connection to south node ---
!
            IF( J.GT.1 ) THEN
              IF( IXP(N-IFLD).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N-IFLD)-IXP(N)))
            ENDIF
!
!---        Connection to west node ---
!
            IF( I.GT.1 ) THEN
              IF( IXP(N-1).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N-1)-IXP(N)))
            ENDIF
!
!---        Connection to east node ---
!
            IF( I.LT.IFLD ) THEN
              IF( IXP(N+1).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N+1)-IXP(N)))
            ENDIF
!
!---        Connection to north node ---
!
            IF( J.LT.JFLD ) THEN
              IF( IXP(N+IFLD).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N+IFLD)-IXP(N)))
            ENDIF
!
!---        Connection to top node ---
!
            IF( K.LT.KFLD ) THEN
              IF( IXP(N+IJFLD).GT.0 )
     &          MHBW = MAX(MHBW,ABS(IXP(N+IJFLD)-IXP(N)))
            ENDIF
 510      CONTINUE
 520    CONTINUE
 530  CONTINUE
!
!---  Determine the matrix half-band widths including 
!     surface spill layer ---
!
      DO 560 I = 1,IFLD
        DO 550 J = 1,JFLD
          NX = NFLD + (J-1)*IFLD + I
          K = KSPS(I,J)
!
!---      Connection to field node (infiltration) ---
!
          IF( K.GT.0 ) THEN
            N = ND(I,J,K)
            IF( IXP(N).GT.0 ) MHBW = MAX(MHBW,ABS(IXP(N)-IXP(NX)))
          ENDIF
!
!---      Connection to southern surface node (surface flux) ---
!
          IF( J.GT.1 ) MHBW = MAX(MHBW,ABS(IXP(NX)-IXP(NX-IFLD)))
!
!---      Connection to western surface node (surface flux) ---
!
          IF( I.GT.1 ) MHBW = MAX(MHBW,ABS(IXP(NX)-IXP(NX-1)))
!
!---      Connection to eastern surface node (surface flux) ---
!
          IF( I.LT.IFLD ) MHBW = MAX(MHBW,ABS(IXP(NX)-IXP(NX+1)))
!
!---      Connection to northern surface node (surface flux) ---
!
          IF( J.LT.JFLD ) MHBW = MAX(MHBW,ABS(IXP(NX)-IXP(NX+IFLD)))
 550    CONTINUE
 560  CONTINUE
!
      MLC = ISVC*MHBW + ISVC - 1
      MLT = ISVT*MHBW + ISVT - 1
      MUC = ISVC*MHBW + ISVC - 1
      MUT = ISVT*MHBW + ISVT - 1
      MDC = MLC + MUC + 1
      MDT = MLT + MUT + 1
!
!---  SPLIB Solver  ---
!
      IF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    Load Jacobian matrix in the increment order I,J,K  ---
!
        NC = 0
        NCC = 0
        NLU(1) = NC+1
        NLUC(1) = NCC+1
        DO 690 K = 1,KFLD
        DO 690 J = 1,JFLD
        DO 690 I = 1,IFLD
          N = ND(I,J,K)
          NP = IXP(N)
          IF( NP.LE.0 ) GOTO 690
          DO 680 L = 1,ISVC
!
!---        Bottom node  ---
!
            IF( K.GT.1 ) THEN
              NMD = IM(L,NP)
              NB = N - IJFLD
              NB = IXP(NB)
              IF( NB.LE.0 ) GOTO 611
              DO 610 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NB)
                KLU(NMD,M) = NC
  610         CONTINUE
            ENDIF
  611       CONTINUE
!
!---        South node  ---
!
            IF( J.GT.1 ) THEN
              NMD = IM(L,NP)
              NS = N - IFLD
              NS = IXP(NS)
              IF( NS.LE.0 ) GOTO 621
              DO 620 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NS)
                KLU(NMD,M+ISVC) = NC
  620         CONTINUE
            ENDIF
  621       CONTINUE
!
!---        West node  ---
!
            IF( I.GT.1 ) THEN
              NMD = IM(L,NP)
              NW = N - 1
              NW = IXP(NW)
              IF( NW.LE.0 ) GOTO 631
              DO 630 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NW)
                KLU(NMD,M+2*ISVC) = NC
  630         CONTINUE
            ENDIF
  631       CONTINUE
!
!---        Node  ---
!
            DO 640 M = 1,ISVC
              NMD = IM(L,NP)
              NC = NC+1
              MLU(NC) = IM(M,NP)
              KLU(NMD,M+3*ISVC) = NC
  640       CONTINUE
!
!---        East node  ---
!
            IF( I.LT.IFLD ) THEN
              NMD = IM(L,NP)
              NE = N + 1
              NE = IXP(NE)
              IF( NE.LE.0 ) GOTO 651
              DO 650 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NE)
                KLU(NMD,M+4*ISVC) = NC
  650         CONTINUE
            ENDIF
  651       CONTINUE
!
!---        North node  ---
!
            IF( J.LT.JFLD ) THEN
              NMD = IM(L,NP)
              NN = N + IFLD
              NN = IXP(NN)
              IF( NN.LE.0 ) GOTO 661
              DO 660 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NN)
                KLU(NMD,M+5*ISVC) = NC
  660         CONTINUE
            ENDIF
  661       CONTINUE
!
!---        Top node  ---
!
            IF( K.LT.KFLD .AND. K.LT.KSPS(I,J) ) THEN
              NMD = IM(L,NP)
              NT = N + IJFLD
              NT = IXP(NT)
              IF( NT.LE.0 ) GOTO 671
              DO 670 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NT)
                KLU(NMD,M+6*ISVC) = NC
  670         CONTINUE
            ENDIF
  671       CONTINUE
!
!---        Surface node  ---
!
            IF( K.EQ.KSPS(I,J) ) THEN
              NMD = IM(L,NP)
              NT = NFLD + (J-1)*IFLD + I
              NT = IXP(NT)
              IF( NT.LE.0 ) GOTO 673
              DO 672 M = 1,ISVC
                NC = NC+1
                MLU(NC) = IM(M,NT)
                KLU(NMD,M+6*ISVC) = NC
  672         CONTINUE
            ENDIF
  673       CONTINUE
            NMD = IM(L,NP)
            NLU(NMD+1) = NC+1
  680     CONTINUE
!
!---      Transport: Node  ---
!
          MA = 1
          NCC = NCC+1
          MLUC(NCC) = NP
          KLUC(NP,MA) = NCC
          MA = MA + 1
!
!---      Transport: Bottom node  ---
!
          IF( K.GT.1 ) THEN
            NB = N - IJFLD
            NB = IXP(NB)
            IF( NB.LE.0 ) GOTO 681
            NCC = NCC+1
            MLUC(NCC) = NB
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  681     CONTINUE
!
!---      Transport: South node  ---
!
          IF( J.GT.1 ) THEN
            NS = N - IFLD
            NS = IXP(NS)
            IF( NS.LE.0 ) GOTO 682
            NCC = NCC+1
            MLUC(NCC) = NS
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  682     CONTINUE
!
!---      Transport: West node  ---
!
          IF( I.GT.1 ) THEN
            NW = N - 1
            NW = IXP(NW)
            IF( NW.LE.0 ) GOTO 683
            NCC = NCC+1
            MLUC(NCC) = NW
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  683     CONTINUE
!
!---      Transport: East node  ---
!
          IF( I.LT.IFLD ) THEN
            NE = N + 1
            NE = IXP(NE)
            IF( NE.LE.0 ) GOTO 685
            NCC = NCC+1
            MLUC(NCC) = NE
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  685     CONTINUE
!
!---      Transport: North node  ---
!
          IF( J.LT.JFLD ) THEN
            NN = N + IFLD
            NN = IXP(NN)
            IF( NN.LE.0 ) GOTO 686
            NCC = NCC+1
            MLUC(NCC) = NN
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  686     CONTINUE
!
!---      Transport: Top node  ---
!
          IF( K.LT.KFLD ) THEN
            NT = N + IJFLD
            NT = IXP(NT)
            IF( NT.LE.0 ) GOTO 687
            NCC = NCC+1
            MLUC(NCC) = NT
            KLUC(NP,MA) = NCC
            MA = MA + 1
          ENDIF
  687     CONTINUE
          NLUC(NP+1) = NCC+1
  690   CONTINUE
        MKT = NCC
        DO 990 I = 1,IFLD
        DO 990 J = 1,JFLD
          NX = NFLD + (J-1)*IFLD + I
          K = KSPS(I,J)
          NP = IXP(NX)
          DO 980 L = 1,ISVC
!
!---        Connection to field node ---
!
            IF( K.GT.0 ) THEN
              NMD = IM(L,NP)
              NB = IXP(ND(I,J,K))
              IF( NB.GT.0 ) THEN
                DO 910 M = 1,ISVC
                  NC = NC+1
                  MLU(NC) = IM(M,NB)
                  KLU(NMD,M) = NC
  910           CONTINUE
              ENDIF
            ENDIF
!
!---        Connection to southern surface node ---
!
            IF( J.GT.1 ) THEN
              NMD = IM(L,NP)
              NS = IXP(NX-IFLD)
              IF( NS.GT.0 ) THEN
                DO 920 M = 1,ISVC
                  NC = NC+1
                  MLU(NC) = IM(M,NS)
                  KLU(NMD,M+ISVC) = NC
  920           CONTINUE
              ENDIF
            ENDIF
!
!---        Connection to western surface node ---
!
            IF( I.GT.1 ) THEN
              NMD = IM(L,NP)
              NW = IXP(NX-1)
              IF( NW.GT.0 ) THEN
                DO 930 M = 1,ISVC
                  NC = NC+1
                  MLU(NC) = IM(M,NW)
                  KLU(NMD,M+2*ISVC) = NC
  930           CONTINUE
              ENDIF
            ENDIF
!
!---        Surface node ---
!
            DO 940 M = 1,ISVC
              NMD = IM(L,NP)
              NC = NC+1
              MLU(NC) = IM(M,NP)
              KLU(NMD,M+3*ISVC) = NC
  940       CONTINUE
!
!---        Connection to eastern surface node ---
!
            IF( I.LT.IFLD ) THEN
              NMD = IM(L,NP)
              NE = IXP(NX+1)
              IF( NE.GT.0 ) THEN
                DO 950 M = 1,ISVC
                  NC = NC+1
                  MLU(NC) = IM(M,NE)
                  KLU(NMD,M+4*ISVC) = NC
  950           CONTINUE
              ENDIF
            ENDIF
!
!---        Connection to northern surface node ---
!
            IF( J.LT.JFLD ) THEN
              NMD = IM(L,NP)
              NN = IXP(NX+IFLD)
              IF( NN.GT.0 ) THEN
                DO 960 M = 1,ISVC
                  NC = NC+1
                  MLU(NC) = IM(M,NN)
                  KLU(NMD,M+5*ISVC) = NC
  960           CONTINUE
              ENDIF
            ENDIF
  980     CONTINUE
  990   CONTINUE
        MKC = NC
      ENDIF
!
!---  Alternate algorithm for well simulations  ---
!
      ELSE
!      
!---    Determine maximum well length  ---
!
        KWX = 0
        DO 1010 NWL = 1,NWLS
          KWX = MAX( KWX,(IWLDM(4,NWL)-IWLDM(3,NWL)+1))
 1010   CONTINUE
        IJFLDX = IFLD*JFLD
        JKFLDX = JFLD*(KFLD+KWX)
        KIFLDX = IFLD*(KFLD+KWX)
        NC = 0
!
!---    X-Y Plane yields the lowest band width;
!       Load Jacobian matrix in the increment order I,J,K  ---
!
        IF( IJFLDX.LE.JKFLDX .AND. IJFLDX.LE.KIFLDX ) THEN
          DO 1030 K = 1,KFLD
          DO 1030 J = 1,JFLD
          DO 1030 I = 1,IFLD
            N = ND(I,J,K)
            IF( IXP(N).EQ.0 ) GOTO 1030
            NMD = IXP(N)
            ISVCX = ISVC
            IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
            DO 1020 M = 1,ISVCX
              NC = NC+1
              IM(M,NMD) = NC
 1020       CONTINUE
 1030     CONTINUE
!
!---    Y-Z Plane yields the lowest band width;
!       load Jacobian matrix in the increment order J,K,I  ---
!
        ELSEIF( JKFLDX.LE.IJFLDX .AND. JKFLDX.LE.KIFLDX ) THEN
          DO 1050 I = 1,IFLD
          DO 1050 K = 1,KFLD
          DO 1050 J = 1,JFLD
            N = ND(I,J,K)
            IF( IXP(N).EQ.0 ) GOTO 1050
            NMD = IXP(N)
            ISVCX = ISVC
            IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
            DO 1040 M = 1,ISVCX
              NC = NC+1
              IM(M,NMD) = NC
 1040       CONTINUE
 1050     CONTINUE
!
!---    Z-X Plane yields the lowest band width;
!       load Jacobian matrix in the increment order K,I,J  ---
!
        ELSEIF( KIFLDX.LE.IJFLDX .AND. KIFLDX.LE.JKFLDX ) THEN
          DO 1070 J = 1,JFLD
          DO 1070 I = 1,IFLD
          DO 1070 K = 1,KFLD
            N = ND(I,J,K)
            IF( IXP(N).EQ.0 ) GOTO 1070
            NMD = IXP(N)
            ISVCX = ISVC
            IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
            DO 1060 M = 1,ISVCX
              NC = NC+1
              IM(M,NMD) = NC
 1060       CONTINUE
 1070     CONTINUE
        ENDIF
!
!---  Determine the matrix half-band widths  ---
!
        MHBC = 0
        MHBT = 0
        DO 1100 K = 1,KFLD
        DO 1100 J = 1,JFLD
        DO 1100 I = 1,IFLD
          N = ND(I,J,K)
          IF( IXP(N).EQ.0 ) GOTO 1100
          NP = IXP(N)
          IF( K.GT.1 ) THEN
            IF( IXP(N-IJFLD).GT.0 ) THEN
              NB = IXP(N-IJFLD)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NB)),
     &          ABS(IM(ISVC,NP)-IM(1,NB)))
              MHBT = MAX(MHBT,ABS(NP-NB))
            ENDIF
          ENDIF
          IF( J.GT.1 ) THEN
            IF( IXP(N-IFLD).GT.0 ) THEN
              NS = IXP(N-IFLD)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NS)),
     &          ABS(IM(ISVC,NP)-IM(1,NS)))
              MHBT = MAX(MHBT,ABS(NP-NS))
            ENDIF
          ENDIF
          IF( I.GT.1 ) THEN
            IF( IXP(N-1).GT.0 ) THEN
              NW = IXP(N-1)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NW)),
     &          ABS(IM(ISVC,NP)-IM(1,NW)))
              MHBT = MAX(MHBT,ABS(NP-NW))
            ENDIF
          ENDIF
          IF( I.LT.IFLD ) THEN
            IF( IXP(N+1).GT.0 ) THEN
              NE = IXP(N+1)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NE)),
     &          ABS(IM(ISVC,NP)-IM(1,NE)))
              MHBT = MAX(MHBT,ABS(NP-NE))
            ENDIF
          ENDIF
          IF( J.LT.JFLD ) THEN
            IF( IXP(N+IFLD).GT.0 ) THEN
              NN = IXP(N+IFLD)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NN)),
     &           ABS(IM(ISVC,NP)-IM(1,NN)))
              MHBT = MAX(MHBT,ABS(NP-NN))
            ENDIF
          ENDIF
          IF( K.LT.KFLD ) THEN
            IF( IXP(N+IJFLD).GT.0 ) THEN
              NT = IXP(N+IJFLD)
              MHBC = MAX(MHBC,ABS(IM(1,NP)-IM(ISVC,NT)),
     &          ABS(IM(ISVC,NP)-IM(1,NT)))
              MHBT = MAX(MHBT,ABS(NP-NT))
            ENDIF
          ENDIF
 1100   CONTINUE
        MLC = MHBC
        MLT = ISVT*MHBT + ISVT - 1
        MUC = MHBC
        MUT = ISVT*MHBT + ISVT - 1
        MDC = MLC + MUC + 1
        MDT = MLT + MUT + 1
!
!---    SPLIB Solver  ---
!
        IF( ILES.EQ.3 .OR. ILES.EQ.4 ) THEN
!
!---    X-Y Plane yields the lowest band width.
!       Load Jacobian matrix in the increment order I,J,K  ---
!
          IF( IJFLDX.LE.JKFLDX .AND. IJFLDX.LE.KIFLDX ) THEN
            NC = 0
            NCC = 0
            NLU(1) = NC+1
            NLUC(1) = NCC+1
            DO 2000 K = 1,KFLD
            DO 2000 J = 1,JFLD
            DO 2000 I = 1,IFLD
              N = ND(I,J,K)
              NP = IXP(N)
              IF( NP.LE.0 ) GOTO 2000
              ISVCX = ISVC
              IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
              DO 1900 L = 1,ISVCX
                NMD = IM(L,NP)
!
!---            Bottom node  ---
!
                IF( K.GT.1 ) THEN
                  NB = N-IJFLD
                  NB = IXP(NB)
                  IF( NB.LE.0 ) GOTO 1240
!
!---              Well node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 1200 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NB)
                        KLU(NMD,M) = NC
 1200                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N-IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 1210 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NB)
                          KLU(NMD,M) = NC
 1210                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    DO 1220 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NB)
                      KLU(NMD,M) = NC
 1220               CONTINUE
                  ENDIF
                ENDIF
 1240         CONTINUE
!
!---          South node  ---
!
                IF( J.GT.1 ) THEN
                  NS = N-IFLD
                  NS = IXP(NS)
                  IF( NS.LE.0 ) GOTO 1260
                  MA = ISVC
                  IF( IXW(N).NE.0 ) MA = 2*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 1250 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NS)
                      KLU(NMD,M+MA) = NC
 1250               CONTINUE
                  ENDIF
                ENDIF
 1260           CONTINUE
!
!---            West node  ---
!
                IF( I.GT.1 ) THEN
                  NW = N-1
                  NW = IXP(NW)
                  IF( NW.LE.0 ) GOTO 1280
                  MA = 2*ISVC
                  IF( IXW(N).NE.0 ) MA = 3*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 1270 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NW)
                      KLU(NMD,M+MA) = NC
 1270               CONTINUE
                  ENDIF
                ENDIF
 1280           CONTINUE
!
!---            Node  ---
!
                MA = 3*ISVC
                IF( IXW(N).NE.0 ) MA = 4*ISVC
!
!---            Well screen node  ---
!
                IF( IXW(N).LT.0 ) THEN
!
!---              Field equations  ---
!
                  DO 1290 M = 1,2*ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 1290             CONTINUE
!
!---            Well casing node  ---
!
                ELSEIF( IXW(N).GT.0 ) THEN
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 1300 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 1300               CONTINUE
                  ENDIF
!
!---              Well equations  ---
!
                  IF( L.GT.ISVC ) THEN
                    DO 1310 M = ISVC+1,2*ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 1310               CONTINUE
                  ENDIF
!
!---            Non-well node  ---
!
                ELSE
!
!---              Field equations  ---
!
                  DO 1320 M = 1,ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 1320             CONTINUE
                ENDIF
!
!---            East node  ---
!
                IF( I.LT.IFLD ) THEN
                  NE = N+1
                  NE = IXP(NE)
                  IF( NE.LE.0 ) GOTO 1340
                  MA = 4*ISVC
                  IF( IXW(N).NE.0 ) MA = 6*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 1330 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NE)
                      KLU(NMD,M+MA) = NC
 1330               CONTINUE
                  ENDIF
                ENDIF
 1340           CONTINUE
!
!---            North node  ---
!
                IF( J.LT.JFLD ) THEN
                  NN = N+IFLD
                  NN = IXP(NN)
                  IF( NN.LE.0 ) GOTO 1360
                  MA = 5*ISVC
                  IF( IXW(N).NE.0 ) MA = 7*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 1350 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NN)
                      KLU(NMD,M+MA) = NC
 1350               CONTINUE
                  ENDIF
                ENDIF
 1360           CONTINUE
!
!---            Top node  ---
!
                IF( K.LT.KFLD ) THEN
                  NT = N+IJFLD
                  NT = IXP(NT)
                  IF( NT.LE.0 ) GOTO 1400
                  MA = 6*ISVC
                  IF( IXW(N).NE.0 ) MA = 8*ISVC
!
!---              Well screen node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 1370 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 1370                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N+IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 1380 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NT)
                          KLU(NMD,M+MA) = NC
 1380                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 1390 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 1390                 CONTINUE
                    ENDIF
                  ENDIF
                ENDIF
 1400           CONTINUE
                NLU(NMD+1) = NC+1
 1900         CONTINUE
!
!---          Solute transport equations  ---
!
              MA = 1
              NCC = NCC+1
              MLUC(NCC) = NP
              KLUC(NP,MA) = NCC
              MA = MA + 1
              IF( K.GT.1 ) THEN
                NB = N - IJFLD
                NB = IXP(NB)
                IF( NB.LE.0 ) GOTO 1910
                NCC = NCC+1
                MLUC(NCC) = NB
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1910         CONTINUE
              IF( J.GT.1 ) THEN
                NS = N - IFLD
                NS = IXP(NS)
                IF( NS.LE.0 ) GOTO 1920
                NCC = NCC+1
                MLUC(NCC) = NS
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1920         CONTINUE
              IF( I.GT.1 ) THEN
                NW = N - 1
                NW = IXP(NW)
                IF( NW.LE.0 ) GOTO 1930
                NCC = NCC+1
                MLUC(NCC) = NW
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1930         CONTINUE
              IF( I.LT.IFLD ) THEN
                NE = N + 1
                NE = IXP(NE)
                IF( NE.LE.0 ) GOTO 1940
                NCC = NCC+1
                MLUC(NCC) = NE
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1940         CONTINUE
              IF( J.LT.JFLD ) THEN
                NN = N + IFLD
                NN = IXP(NN)
                IF( NN.LE.0 ) GOTO 1950
                NCC = NCC+1
                MLUC(NCC) = NN
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1950         CONTINUE
              IF( K.LT.KFLD ) THEN
                NT = N + IJFLD
                NT = IXP(NT)
                IF( NT.LE.0 ) GOTO 1960
                NCC = NCC+1
                MLUC(NCC) = NT
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 1960         CONTINUE
              NLUC(NP+1) = NCC+1
 2000       CONTINUE
            MKC = NC
            MKT = NCC
!
!---    Y-Z Plane yields the lowest band width.
!       Load Jacobian matrix in the increment order J,K,I  ---
!
          ELSEIF( JKFLDX.LE.IJFLDX .AND. JKFLDX.LE.KIFLDX ) THEN
            NC = 0
            NCC = 0
            NLU(1) = NC+1
            NLUC(1) = NCC+1
            DO 4000 I = 1,IFLD
            DO 4000 K = 1,KFLD
            DO 4000 J = 1,JFLD
              N = ND(I,J,K)
              NP = IXP(N)
              IF( NP.LE.0 ) GOTO 4000
              ISVCX = ISVC
              IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
              DO 3900 L = 1,ISVCX
                NMD = IM(L,NP)
!
!---            Bottom node  ---
!
                IF( K.GT.1 ) THEN
                  NB = N-IJFLD
                  NB = IXP(NB)
                  IF( NB.LE.0 ) GOTO 3240
!
!---              Well node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 3200 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NB)
                        KLU(NMD,M) = NC
 3200                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N-IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 3210 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NB)
                          KLU(NMD,M) = NC
 3210                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    DO 3220 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NB)
                      KLU(NMD,M) = NC
 3220               CONTINUE
                  ENDIF
                ENDIF
 3240         CONTINUE
!
!---          South node  ---
!
                IF( J.GT.1 ) THEN
                  NS = N-IFLD
                  NS = IXP(NS)
                  IF( NS.LE.0 ) GOTO 3260
                  MA = ISVC
                  IF( IXW(N).NE.0 ) MA = 2*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 3250 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NS)
                      KLU(NMD,M+MA) = NC
 3250               CONTINUE
                  ENDIF
                ENDIF
 3260           CONTINUE
!
!---            West node  ---
!
                IF( I.GT.1 ) THEN
                  NW = N-1
                  NW = IXP(NW)
                  IF( NW.LE.0 ) GOTO 3280
                  MA = 2*ISVC
                  IF( IXW(N).NE.0 ) MA = 3*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 3270 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NW)
                      KLU(NMD,M+MA) = NC
 3270               CONTINUE
                  ENDIF
                ENDIF
 3280           CONTINUE
!
!---            Node  ---
!
                MA = 3*ISVC
                IF( IXW(N).NE.0 ) MA = 4*ISVC
!
!---            Well screen node  ---
!
                IF( IXW(N).LT.0 ) THEN
!
!---              Field equations  ---
!
                  DO 3290 M = 1,2*ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 3290             CONTINUE
!
!---            Well casing node  ---
!
                ELSEIF( IXW(N).GT.0 ) THEN
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 3300 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 3300               CONTINUE
                  ENDIF
!
!---              Well equations  ---
!
                  IF( L.GT.ISVC ) THEN
                    DO 3310 M = ISVC+1,2*ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 3310               CONTINUE
                  ENDIF
!
!---            Non-well node  ---
!
                ELSE
!
!---              Field equations  ---
!
                  DO 3320 M = 1,ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 3320             CONTINUE
                ENDIF
!
!---            East node  ---
!
                IF( I.LT.IFLD ) THEN
                  NE = N+1
                  NE = IXP(NE)
                  IF( NE.LE.0 ) GOTO 3340
                  MA = 4*ISVC
                  IF( IXW(N).NE.0 ) MA = 6*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 3330 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NE)
                      KLU(NMD,M+MA) = NC
 3330               CONTINUE
                  ENDIF
                ENDIF
 3340           CONTINUE
!
!---            North node  ---
!
                IF( J.LT.JFLD ) THEN
                  NN = N+IFLD
                  NN = IXP(NN)
                  IF( NN.LE.0 ) GOTO 3360
                  MA = 5*ISVC
                  IF( IXW(N).NE.0 ) MA = 7*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 3350 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NN)
                      KLU(NMD,M+MA) = NC
 3350               CONTINUE
                  ENDIF
                ENDIF
 3360           CONTINUE
!
!---            Top node  ---
!
                IF( K.LT.KFLD ) THEN
                  NT = N+IJFLD
                  NT = IXP(NT)
                  IF( NT.LE.0 ) GOTO 3400
                  MA = 6*ISVC
                  IF( IXW(N).NE.0 ) MA = 8*ISVC
!
!---              Well screen node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 3370 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 3370                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N+IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 3380 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NT)
                          KLU(NMD,M+MA) = NC
 3380                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 3390 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 3390                 CONTINUE
                    ENDIF
                  ENDIF
                ENDIF
 3400           CONTINUE
                NLU(NMD+1) = NC+1
 3900         CONTINUE
!
!---          Solute transport equations  ---
!
              MA = 1
              NCC = NCC+1
              MLUC(NCC) = NP
              KLUC(NP,MA) = NCC
              MA = MA + 1
              IF( K.GT.1 ) THEN
                NB = N - IJFLD
                NB = IXP(NB)
                IF( NB.LE.0 ) GOTO 3910
                NCC = NCC+1
                MLUC(NCC) = NB
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3910         CONTINUE
              IF( J.GT.1 ) THEN
                NS = N - IFLD
                NS = IXP(NS)
                IF( NS.LE.0 ) GOTO 3920
                NCC = NCC+1
                MLUC(NCC) = NS
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3920         CONTINUE
              IF( I.GT.1 ) THEN
                NW = N - 1
                NW = IXP(NW)
                IF( NW.LE.0 ) GOTO 3930
                NCC = NCC+1
                MLUC(NCC) = NW
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3930         CONTINUE
              IF( I.LT.IFLD ) THEN
                NE = N + 1
                NE = IXP(NE)
                IF( NE.LE.0 ) GOTO 3940
                NCC = NCC+1
                MLUC(NCC) = NE
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3940         CONTINUE
              IF( J.LT.JFLD ) THEN
                NN = N + IFLD
                NN = IXP(NN)
                IF( NN.LE.0 ) GOTO 3950
                NCC = NCC+1
                MLUC(NCC) = NN
                 KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3950         CONTINUE
              IF( K.LT.KFLD ) THEN
                NT = N + IJFLD
                NT = IXP(NT)
                IF( NT.LE.0 ) GOTO 3960
                NCC = NCC+1
                MLUC(NCC) = NT
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 3960         CONTINUE
              NLUC(NP+1) = NCC+1
 4000       CONTINUE
            MKC = NC
            MKT = NCC
!
!---    Z-X Plane yields the lowest band width.
!       Load Jacobian matrix in the increment order K,I,J  ---
!
          ELSEIF( KIFLDX.LE.IJFLDX .AND. KIFLDX.LE.JKFLDX ) THEN
            NC = 0
            NCC = 0
            NLU(1) = NC+1
            NLUC(1) = NCC+1
            DO 6000 J = 1,JFLD
            DO 6000 I = 1,IFLD
            DO 6000 K = 1,KFLD
              N = ND(I,J,K)
              NP = IXP(N)
              IF( NP.LE.0 ) GOTO 6000
              ISVCX = ISVC
              IF( IXW(N).NE.0 ) ISVCX = 2*ISVC
              DO 5900 L = 1,ISVCX
                NMD = IM(L,NP)
!
!---            Bottom node  ---
!
                IF( K.GT.1 ) THEN
                  NB = N-IJFLD
                  NB = IXP(NB)
                  IF( NB.LE.0 ) GOTO 5240
!
!---              Well node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 5200 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NB)
                        KLU(NMD,M) = NC
 5200                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N-IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 5210 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NB)
                          KLU(NMD,M) = NC
 5210                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    DO 5220 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NB)
                      KLU(NMD,M) = NC
 5220               CONTINUE
                  ENDIF
                ENDIF
 5240         CONTINUE
!
!---          South node  ---
!
                IF( J.GT.1 ) THEN
                  NS = N-IFLD
                  NS = IXP(NS)
                  IF( NS.LE.0 ) GOTO 5260
                  MA = ISVC
                  IF( IXW(N).NE.0 ) MA = 2*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 5250 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NS)
                      KLU(NMD,M+MA) = NC
 5250               CONTINUE
                  ENDIF
                ENDIF
 5260           CONTINUE
!
!---            West node  ---
!
                IF( I.GT.1 ) THEN
                  NW = N-1
                  NW = IXP(NW)
                  IF( NW.LE.0 ) GOTO 5280
                  MA = 2*ISVC
                  IF( IXW(N).NE.0 ) MA = 3*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 5270 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NW)
                      KLU(NMD,M+MA) = NC
 5270               CONTINUE
                  ENDIF
                ENDIF
 5280           CONTINUE
!
!---            Node  ---
!
                MA = 3*ISVC
                IF( IXW(N).NE.0 ) MA = 4*ISVC
!
!---            Well screen node  ---
!
                IF( IXW(N).LT.0 ) THEN
!
!---              Field equations  ---
!
                  DO 5290 M = 1,2*ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 5290             CONTINUE
!
!---            Well casing node  ---
!
                ELSEIF( IXW(N).GT.0 ) THEN
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 5300 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 5300               CONTINUE
                  ENDIF
!
!---              Well equations  ---
!
                  IF( L.GT.ISVC ) THEN
                    DO 5310 M = ISVC+1,2*ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NP)
                      KLU(NMD,M+MA) = NC
 5310               CONTINUE
                  ENDIF
!
!---            Non-well node  ---
!
                ELSE
!
!---              Field equations  ---
!
                  DO 5320 M = 1,ISVC
                    NC = NC+1
                    MLU(NC) = IM(M,NP)
                    KLU(NMD,M+MA) = NC
 5320             CONTINUE
                ENDIF
!
!---            East node  ---
!
                IF( I.LT.IFLD ) THEN
                  NE = N+1
                  NE = IXP(NE)
                  IF( NE.LE.0 ) GOTO 5340
                  MA = 4*ISVC
                  IF( IXW(N).NE.0 ) MA = 6*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 5330 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NE)
                      KLU(NMD,M+MA) = NC
 5330               CONTINUE
                  ENDIF
                ENDIF
 5340           CONTINUE
!
!---            North node  ---
!
                IF( J.LT.JFLD ) THEN
                  NN = N+IFLD
                  NN = IXP(NN)
                  IF( NN.LE.0 ) GOTO 5360
                  MA = 5*ISVC
                  IF( IXW(N).NE.0 ) MA = 7*ISVC
!
!---              Field equations  ---
!
                  IF( L.LE.ISVC ) THEN
                    DO 5350 M = 1,ISVC
                      NC = NC+1
                      MLU(NC) = IM(M,NN)
                      KLU(NMD,M+MA) = NC
 5350               CONTINUE
                  ENDIF
                ENDIF
 5360           CONTINUE
!
!---            Top node  ---
!
                IF( K.LT.KFLD ) THEN
                  NT = N+IJFLD
                  NT = IXP(NT)
                  IF( NT.LE.0 ) GOTO 5400
                  MA = 6*ISVC
                  IF( IXW(N).NE.0 ) MA = 8*ISVC
!
!---              Well screen node  ---
!
                  IF( IXW(N).NE.0 ) THEN
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 5370 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 5370                 CONTINUE
                    ENDIF
!
!---                Well equations  ---
!
                    IF( IXW(N+IJFLD).NE.0 ) THEN
                      IF( L.GT.ISVC ) THEN
                        DO 5380 M = ISVC+1,2*ISVC
                          NC = NC+1
                          MLU(NC) = IM(M,NT)
                          KLU(NMD,M+MA) = NC
 5380                   CONTINUE
                      ENDIF
                    ENDIF
!
!---              Non-well node  ---
!
                  ELSE
!
!---                Field equations  ---
!
                    IF( L.LE.ISVC ) THEN
                      DO 5390 M = 1,ISVC
                        NC = NC+1
                        MLU(NC) = IM(M,NT)
                        KLU(NMD,M+MA) = NC
 5390                 CONTINUE
                    ENDIF
                  ENDIF
                ENDIF
 5400           CONTINUE
                NLU(NMD+1) = NC+1
 5900         CONTINUE
!
!---          Solute transport equations  ---
!
              MA = 1
              NCC = NCC+1
              MLUC(NCC) = NP
              KLUC(NP,MA) = NCC
              MA = MA + 1
              IF( K.GT.1 ) THEN
                NB = N - IJFLD
                NB = IXP(NB)
                IF( NB.LE.0 ) GOTO 5910
                NCC = NCC+1
                MLUC(NCC) = NB
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5910         CONTINUE
              IF( J.GT.1 ) THEN
                NS = N - IFLD
                NS = IXP(NS)
                IF( NS.LE.0 ) GOTO 5920
                NCC = NCC+1
                MLUC(NCC) = NS
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5920         CONTINUE
              IF( I.GT.1 ) THEN
                NW = N - 1
                NW = IXP(NW)
                IF( NW.LE.0 ) GOTO 5930
                NCC = NCC+1
                MLUC(NCC) = NW
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5930         CONTINUE
              IF( I.LT.IFLD ) THEN
                NE = N + 1
                NE = IXP(NE)
                IF( NE.LE.0 ) GOTO 5940
                NCC = NCC+1
                MLUC(NCC) = NE
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5940         CONTINUE
              IF( J.LT.JFLD ) THEN
                NN = N + IFLD
                NN = IXP(NN)
                IF( NN.LE.0 ) GOTO 5950
                NCC = NCC+1
                MLUC(NCC) = NN
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5950         CONTINUE
              IF( K.LT.KFLD ) THEN
                NT = N + IJFLD
                NT = IXP(NT)
                IF( NT.LE.0 ) GOTO 5960
                NCC = NCC+1
                MLUC(NCC) = NT
                KLUC(NP,MA) = NCC
                MA = MA + 1
              ENDIF
 5960         CONTINUE
              NLUC(NP+1) = NCC+1
 6000       CONTINUE
            MKC = NC
            MKT = NCC
          ENDIF
        ENDIF
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBPSP group  ---
!
      RETURN
      END
!  
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCSPO
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
!     Load the Jacobian matrix for the NAPL surface spill head.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 12 December, 2007.
!     Last Modified by SK Wurstner PNNL, 04 November, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXN
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,5),FW(LSFV,5)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCSPO'
!
!---  Loop over surface nodes ---
!
      DO 2000 J = 1,JFLD
      DO 2000 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        NSPX = (J-1)*(IFLD-1) + I - 1
        NSPY = (J-2)*IFLD + I
        NSQX = NSPX + 1
        NSQY = NSPY + IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
!
!---    First-order, forward-difference, time differential  ---
!
        STWOX = RHONSP(1,NSP)*HNSP(1,NSP)
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STWX(M) = (RHONSP(MP,NSP)*HNSP(MP,NSP)-STWOX)*AFZSP(NSP)*DTI
  100   CONTINUE
!
!---    Initialize surface fluxes  ---
!
        DO 210 MD = 1,5
          DO 200 M = 1,ISVF
            FW(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---    Infiltration to field node (surface node perspective)  ---
!
        IF( IXP(N).GT.0 ) THEN
          NQZ = NSZ(N)+IJFLD
          DO 300 M = 1,ISVF
            MB = MADJ(M)
            MP = MNOD(M)
            FNOB = RHON(MB,N)
            FNOP = RHONSP(MP,NSP)
            INDX = 4
            FNO = DIFMN( FNOB,FNOP,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
            FW(M,1) = -AFZ(NQZ)*WN(M,NQZ)*FNO
  300     CONTINUE
        ENDIF
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          NSS = NSP-IFLD
          DO 400 M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FWS = RHONSP(MS,NSS)
            FWP = RHONSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWS,FWP,DYGF(NS),DYGF(N),VNSP(1,NSPY),INDX )
            FW(M,2) = -VNSP(M,NSPY)*FWM*DXSP(NSPY)
  400     CONTINUE
        ENDIF
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          NSW = NSP-1
          DO 500 M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FWW = RHONSP(MW,NSW)
            FWP = RHONSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWW,FWP,DXGF(NW),DXGF(N),UNSP(1,NSPX),INDX )
            FW(M,3) = -UNSP(M,NSPX)*FWM*DYSP(NSPX)
  500     CONTINUE
        ENDIF
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          NSE = NSP+1
          DO 600 M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FWE = RHONSP(ME,NSE)
            FWP = RHONSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWP,FWE,DXGF(N),DXGF(NE),UNSP(1,NSQX),INDX )
            FW(M,4) = UNSP(MF,NSQX)*FWM*DYSP(NSQX)
  600     CONTINUE
        ENDIF
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          NSN = NSP+IFLD
          DO 700 M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FWN = RHONSP(MN,NSN)
            FWP = RHONSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWP,FWN,DYGF(N),DYGF(NN),VNSP(1,NSQY),INDX )
            FW(M,5) = VNSP(MF,NSQY)*FWM*DXSP(NSQY)
  700     CONTINUE
        ENDIF
!
!---    Compute oil equation residuals  ---
!
        RWS = STWX(1) - SRCOSP(2,NSP)
        DO 900 MD = 1,5
          RWS = RWS + FW(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RWP(M) = STWX(M+1) - STWX(1) - SRCOSP(M+2,NSP) + SRCOSP(2,NSP)
          MM = 2*M
          DO 910 MD = 1,5
            RWP(M) = RWP(M) + FW(MM,MD) - FW(1,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,5
            RWA(M,MD) = FW(MM,MD) - FW(1,MD)
  930     CONTINUE
  940   CONTINUE
!
!---    Load Jacobian matrix for the surface node  ---
!
        CALL JCBLDSP( RWS,RWP,RWA,I,J,K,IEQO )
!
!---    Infiltration from the surface node (field node perspective)  ---
!
        IF( IXP(N).GT.0 ) THEN
          DO 1000 M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FNOP = RHON(MP,N)
            FNOT = RHONSP(MT,NSP)
            INDX = 4
            FNO = DIFMN( FNOP,FNOT,DZGF(N),DZGF(N),WN(1,NQZ),INDX )
            FW(M,1) = AFZ(NQZ)*WN(MF,NQZ)*FNO
 1000     CONTINUE
!
!---      Compute oil equation residuals  ---
!
          RWS = FW(1,1)
          DO 1020 M = 1,ISVC
            MM = 2*M
            RWP(M) = FW(MM,1) - FW(1,1)
 1020     CONTINUE
          DO 1040 M = 1,ISVC
            MM = 2*M + 1
            RWA(M,1) = FW(MM,1) - FW(1,1)
 1040     CONTINUE
!
!---      Modify the Jacobian matrix for the field node  ---
!
          CALL JCBLDFD( RWS,RWP,RWA,I,J,K,IEQO )
        ENDIF
!
!---  Continue to next surface node  ---
!
 2000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCSPO group
!
      RETURN
      END
!  
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCSPW
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
!     Load the Jacobian matrix for the aqueous surface spill head.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 12 December, 2007.
!     Last Modified by SK Wurstner PNNL, 04 November, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 STWX(LUK+1),RWP(LUK),RWA(LUK,5),FW(LSFV,5)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCSPW'
!
!---  Loop over surface nodes  ---
!
      DO 2000 J = 1,JFLD
      DO 2000 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        NSPX = (J-1)*(IFLD-1) + I - 1
        NSPY = (J-2)*IFLD + I
        NSQX = NSPX + 1
        NSQY = NSPY + IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
!
!---    First-order, forward-difference, time differential  ---
!
        STWOX = RHOLSP(1,NSP)*HLSP(1,NSP)
        DO 100 M = 1,ISVC+1
          MP = M + 1
          STWX(M) = (RHOLSP(MP,NSP)*HLSP(MP,NSP)-STWOX)*AFZSP(NSP)*DTI
  100   CONTINUE
!
!---    Initialize surface fluxes  ---
!
        DO 210 MD = 1,5
          DO 200 M = 1,ISVF
            FW(M,MD) = 0.D+0
  200     CONTINUE
  210   CONTINUE
!
!---    Infiltration to field node  ---
!
        K = KSPS(I,J)
        N = ND(I,J,K)
        NQZ = NSZ(N)+IJFLD
        DO 300 M = 1,ISVF
          MB = MADJ(M)
          MP = MNOD(M)
          FLWB = RHOL(MB,N)
          FLWP = RHOLSP(MP,NSP)
          INDX = 2
          FLW = DIFMN( FLWB,FLWP,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
          FW(M,1) = -AFZ(NQZ)*WL(M,NQZ)*FLW
  300   CONTINUE
!
!---    South ---
!
        IF( J.NE.1 ) THEN
          NS = N-IFLD
          NSS = NSP-IFLD
          DO 400 M = 1,ISVF
            MS = MADJ(M)
            MP = MNOD(M)
            FWS = RHOLSP(MS,NSS)
            FWP = RHOLSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWS,FWP,DYGF(NS),DYGF(N),VLSP(1,NSPY),INDX )
            FW(M,2) = -VLSP(M,NSPY)*FWM*DXSP(NSPY)
  400     CONTINUE
        ENDIF
!
!---    West ---
!
        IF( I.NE.1 ) THEN
          NW = N-1
          NSW = NSP-1
          DO 500 M = 1,ISVF
            MW = MADJ(M)
            MP = MNOD(M)
            FWW = RHOLSP(MW,NSW)
            FWP = RHOLSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWW,FWP,DXGF(NW),DXGF(N),ULSP(1,NSPX),INDX )
            FW(M,3) = -ULSP(M,NSPX)*FWM*DYSP(NSPX)
  500     CONTINUE
        ENDIF
!
!---    East ---
!
        IF( I.NE.IFLD ) THEN
          NE = N+1
          NSE = NSP+1
          DO 600 M = 1,ISVF
            ME = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FWE = RHOLSP(ME,NSE)
            FWP = RHOLSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWP,FWE,DXGF(N),DXGF(NE),ULSP(1,NSQX),INDX )
            FW(M,4) = ULSP(MF,NSQX)*FWM*DYSP(NSQX)
  600     CONTINUE
        ENDIF
!
!---    North ---
!
        IF( J.NE.JFLD ) THEN
          NN = N+IFLD
          NSN = NSP+IFLD
          DO 700 M = 1,ISVF
            MN = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FWN = RHOLSP(MN,NSN)
            FWP = RHOLSP(MP,NSP)
            INDX = 4
            FWM = DIFMN( FWP,FWN,DYGF(N),DYGF(NN),VLSP(1,NSQY),INDX )
            FW(M,5) = VLSP(MF,NSQY)*FWM*DXSP(NSQY)
  700     CONTINUE
        ENDIF
!
!---    Compute water equation residuals  ---
!
        RWS = STWX(1)
        DO 900 MD = 1,5
          RWS = RWS + FW(1,MD)
  900   CONTINUE
        DO 920 M = 1,ISVC
          RWP(M) = STWX(M+1) - STWX(1)
          MM = 2*M
          DO 910 MD = 1,5
            RWP(M) = RWP(M) + FW(MM,MD) - FW(1,MD)
  910     CONTINUE
  920   CONTINUE
        DO 940 M = 1,ISVC
          MM = 2*M + 1
          DO 930 MD = 1,5
            RWA(M,MD) = FW(MM,MD) - FW(1,MD) 
  930     CONTINUE
  940   CONTINUE
!
!---    Load Jacobian Matrix  ---
!
        CALL JCBLDSP( RWS,RWP,RWA,I,J,K,IEQW )
!
!---    Infiltration from the surface node (field node perspective)  ---
!
        IF( IXP(N).GT.0 ) THEN
          DO 1000 M = 1,ISVF
            MT = MADJ(M)
            MP = MNOD(M)
            MF = MFLX(M)
            FLWP = RHOL(MP,N)
            FLWT = RHOLSP(MT,NSP)
            INDX = 4
            FLW = DIFMN( FLWP,FLWT,DZGF(N),DZGF(N),WL(1,NQZ),INDX )
            FW(M,1) = AFZ(NQZ)*WL(MF,NQZ)*FLW
 1000     CONTINUE
!
!---      Compute oil equation residuals  ---
!
          RWS = FW(1,1)
          DO 1020 M = 1,ISVC
            MM = 2*M
            RWP(M) = FW(MM,1) - FW(1,1)
 1020     CONTINUE
          DO 1040 M = 1,ISVC
            MM = 2*M + 1
            RWA(M,1) = FW(MM,1) - FW(1,1)
 1040     CONTINUE
!
!---      Modify the Jacobian matrix for the field node  ---
!
          CALL JCBLDFD( RWS,RWP,RWA,I,J,K,IEQW )
        ENDIF
!
!---  Continue to next surface node  ---
!
 2000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCSPW group
!
      RETURN
      END
!  
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RSDSP
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
!     Surface spill residuals.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 10 March, 2008.
!     Last Modified by SK Wurstner PNNL, 10 March, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOURC
      USE SOLTN
      USE OUTPU
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
      SUB_LOG(ISUB_LOG) = '/RSDSP'
      NWELL = LWELL
!
!---  Reset surface spill heads for convergence failure
!     in field nodes  ---
!
      IF( ICNV.EQ.1 .OR. ICNV.EQ.4 ) THEN
        DO 10 J = 1,JFLD
        DO 10 I = 1,IFLD
          NSP = (J-1)*IFLD + I
          HNSP(2,NSP) = HNSP(1,NSP)
          HLSP(2,NSP) = HLSP(1,NSP)
   10   CONTINUE
        GOTO 100
      ENDIF
!
!---  Initialize surface spill residuals  ---
!
      DO 20 M = 1,ISVC
        RSD(M+ISVC+ISVC*NWELL) = 0.D+0
        NSD(M+ISVC+ISVC*NWELL) = 0
   20 CONTINUE
!
!---  Loop over surface spill nodes  ---
!
      DO 30 J = 1,JFLD
      DO 30 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        NX = NFLD + NSP
        NMD = IXP(NX)
        MPL = IM(IEQW,NMD)
        MPN = IM(IEQO,NMD)
!        IF( HLSP(1,NSP).GT.EPSL .OR. HLSP(2,NSP).GT.EPSL ) THEN
        IF( HLSP(2,NSP).GT.1.D-09 .AND. 
     &    (HLSP(2,NSP)+BLU(MPL)).GT.EPSL ) THEN
          RSDX = ABS(BLU(MPL))
        ELSE
          RSDX = 0.D+0
        ENDIF
        IF( RSDX.GT.RSD(IEQW+ISVC+ISVC*NWELL) ) THEN
          RSD(IEQW+ISVC+ISVC*NWELL) = RSDX
          NSD(IEQW+ISVC+ISVC*NWELL) = NSP
        ENDIF
!        IF( HNSP(1,NSP).GT.EPSL .OR. HNSP(2,NSP).GT.EPSL )  THEN
        IF( HNSP(2,NSP).GT.1.D-09 .AND. 
     &    (HNSP(2,NSP)+BLU(MPN)).GT.EPSL ) THEN
          RSDX = ABS(BLU(MPN))
        ELSE
          RSDX = 0.D+0
        ENDIF
        IF( RSDX.GT.RSD(IEQO+ISVC+ISVC*NWELL) ) THEN
          RSD(IEQO+ISVC+ISVC*NWELL) = RSDX
          NSD(IEQO+ISVC+ISVC*NWELL) = NSP
        ENDIF
   30 CONTINUE
!
!---  Assign a convergence index  ---
!
      DO 40 M = 1,ISVC
        IF( RSD(M+ISVC+ISVC*NWELL).GT.RSDMX ) ICNV = 2
   40 CONTINUE
      IF( ICNV.EQ.2 .AND. NITER.GE.NRIMX ) ICNV = 1
!
!---  Unconverged solution Newton-Raphson iteration limit exceeded  ---
!
      IF( ICNV.EQ.1 ) THEN
        WRITE(ISC,'(A)') 'Execution Note: Newton Convergence Failure'
        WRITE(IWR,'(A)') 'Execution Note: Newton Convergence Failure'
        WRITE(ISC,'(A,1PE11.4)') '  Surface Spill Water Residual: ',
     &    RSD(IEQW+ISVC+ISVC*NWELL)
        WRITE(ISC,'(A,I6)') '  Surface Spill Node: ',
     &    NSD(IEQW+ISVC+ISVC*NWELL)
        WRITE(IWR,'(A,1PE11.4)') '  Surface Spill Water Residual: ',
     &    RSD(IEQW+ISVC+ISVC*NWELL)
        WRITE(IWR,'(A,I6)') '  Surface Spill Node: ',
     &    NSD(IEQW+ISVC+ISVC*NWELL)
        WRITE(ISC,'(A,1PE11.4)') '  Surface Spill Oil Residual: ',
     &    RSD(IEQO+ISVC+ISVC*NWELL)
        WRITE(ISC,'(A,I6)') '  Surface Spill Node: ',
     &    NSD(IEQO+ISVC+ISVC*NWELL)
        WRITE(IWR,'(A,1PE11.4)') '  Surface Spill Oil Residual: ',
     &    RSD(IEQO+ISVC+ISVC*NWELL)
        WRITE(IWR,'(A,I6)') '  Surface Spill Node: ',
     &    NSD(IEQO+ISVC+ISVC*NWELL)
        WRITE(ISC,'(A,1PE11.4)') '  Surface Spill Height: ',
     &    HNSP(2,IEQO+ISVC+ISVC*NWELL)
        WRITE(IWR,'(A,1PE11.4)') '  Surface Spill Height: ',
     &    HNSP(2,IEQO+ISVC+ISVC*NWELL)
!
!---    Reset primary field variables  ---
!
        DO 50 N = 1,NFLD
          PL(2,N) = PL(1,N)
          PN(2,N) = PN(1,N)
          XMLO(2,N) = XMLO(1,N)
          SN(2,N) = SN(1,N)
          NPHAZ(2,N) = NPHAZ(1,N)
   50   CONTINUE
!
!---    Reset AWTAP well variables  ---
!
#ifdef awtap
        DO 60 N = 1,NSR
          PLWB(2,N) = PLWB(1,N)
          PGW(2,N) = PGW(1,N)
   60   CONTINUE
#endif
!
!---    Reset surface spill head variables  ---
!
        DO 70 J = 1,JFLD
        DO 70 I = 1,IFLD
          NSP = (J-1)*IFLD + I
          HNSP(2,NSP) = HNSP(1,NSP)
          HLSP(2,NSP) = HLSP(1,NSP)
   70   CONTINUE
!
!---    Reduce time step  ---
!
        IF( NTSR.LT.4 ) THEN
          NTSR = NTSR + 1
          DTX = DT
          TM = TM - (1.D+0-DTCF)*DT
          DT = DTCF*DT
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
          WRITE(ISC,'(A,1PE11.4,1X,$)') '  Time step reduced from ',
     &      VARX
          WRITE(ISC,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
          WRITE(IWR,'(A,1PE11.4,1X,$)') '  Time step reduced from ',
     &      VARX
          WRITE(IWR,'(2A,1PE11.4,1X,A)') UNTM(1:NCH),' to ',
     &      VAR,UNTM(1:NCH)
!
!---    Number of time step reductions failure: stop simulation  ---
!
        ELSE
          WRITE(ISC,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          WRITE(IWR,'(A)') 'Execution Note:  Reduction Limit Exceeded'
          ICNV = 4
        ENDIF
      ENDIF
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RSDSP group
!
      RETURN
      END
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SP_GRID
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
!     Spill grid variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 07 December, 2007.
!     Last Modified by SK Wurstner PNNL, 07 December, 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
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
      SUB_LOG(ISUB_LOG) = '/SP_GRID'
!
!---  Define upper active node for surface spills  ---
!
      KSPSX = 0
      DO 930 J = 1,JFLD
        DO 920 I = 1,IFLD
          DO 910 K = KFLD,1,-1
            N = ND(I,J,K)
            IF( IXP(N).NE.0 ) THEN
              KSPS(I,J) = K
              KSPSX = MAX( KSPSX,KFLD-K ) 
              GOTO 920
            ENDIF
  910     CONTINUE
  920   CONTINUE
  930 CONTINUE
!
!---  Check spill node to field node connection against
!     band width parameter  ---
!
!      IF( (KSPSX+1)*LFX*LFY.GT.LMNP ) THEN
!        INDX = 5
!        CHMSG = 'Spill Node-Field Node Connection > Parameter LMNP'
!        CALL WRMSGS( INDX )
!      ENDIF
!
!---  Load Jacobian matrix for surface spill nodes  ---
!
      NC = NFLD-NXP
      DO 950 J = 1,JFLD
        DO 940 I = 1,IFLD
          NX = NFLD + (J-1)*IFLD + I
          NC = NC+1
          IXP(NX) = NC
  940   CONTINUE
  950 CONTINUE
!
!---  Assign top node coordinates using top active node 
!     to define surface ---
!
      DO 1200 J = 1,JFLD
      DO 1200 I = 1,IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
        NQZ = NSZ(N)+IJFLD
        NSP = (J-1)*IFLD + I
        XSP(NSP) = 2.5D-1*(XE(5,N)+XE(6,N)+XE(7,N)+XE(8,N))
        YSP(NSP) = 2.5D-1*(YE(5,N)+YE(6,N)+YE(7,N)+YE(8,N))
        ZSP(NSP) = 2.5D-1*(ZE(5,N)+ZE(6,N)+ZE(7,N)+ZE(8,N))
        AFZSP(NSP) = AFZ(NQZ)
 1200 CONTINUE
!
!---  Calculate slope of surface between adjacent surface nodes, 
!     calculate the lengths of cell faces in original grid space and
!     the distance between adjacent nodes in the plane of the surface
!     (x direction points downhill) ---
! 
!
      DO 1300 J = 1,JFLD
      DO 1300 I = 2,IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
        NSP = (J-1)*IFLD + I
        NSW = NSP-1
        NSPX = (J-1)*(IFLD-1) + I - 1
        TXSP(NSPX) = ATAN((ZSP(NSP)-ZSP(NSW))/
     &    (SQRT((XSP(NSP)-XSP(NSW))**2 + (YSP(NSP)-YSP(NSW))**2)))
        DYSP(NSPX) = SQRT( (XE(7,N)-XE(5,N))**2
     &    + (YE(7,N)-YE(5,N))**2
     &    + (ZE(7,N)-ZE(5,N))**2 )
        PDXSP(NSPX) = (XSP(NSP)-XSP(NSW))/COS(TXSP(NSPX))
 1300 CONTINUE
!
      DO 1400 J = 2,JFLD
      DO 1400 I = 1,IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
        NSP = (J-1)*IFLD + I
        NSS = NSP-IFLD
        NSPY = (J-2)*IFLD + I
        TYSP(NSPY) = ATAN((ZSP(NSP)-ZSP(NSS))/
     &    (SQRT((XSP(NSP)-XSP(NSS))**2 + (YSP(NSP)-YSP(NSS))**2)))
        DXSP(NSPY) = SQRT( (XE(6,N)-XE(5,N))**2
     &    + (YE(6,N)-YE(5,N))**2
     &    + (ZE(6,N)-ZE(5,N))**2 )
        PDYSP(NSPY) = (YSP(NSP)-YSP(NSS))/COS(TYSP(NSPY))
 1400 CONTINUE
!
!---  Calculate cosine of angle defining normal to the surface for each node
! 
!
      DO 1500 J = 1,JFLD
      DO 1500 I = 1,IFLD
        NSP = (J-1)*IFLD + I
        SPNORM(NSP) = COS(THXZ)*COS(THYZ)
 1500 CONTINUE
!
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SP_GRID group  ---
!
      RETURN
      END
!     
!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE PROPSP
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
!     Aqueous and NAPL surface spill density and viscosity.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 14 December, 2007.
!     Last Modified by SK Wurstner PNNL, 14 December, 2007.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PROPSP'
!
!---  Calculate density and viscosity for each surface node ---
!
      DO 1000 J = 1,JFLD
      DO 1000 I = 1,IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
        NSP = (J-1)*IFLD + I
        TX = T(2,N)
        PX = PATM
        CALL WATSP( TX,PSWX )
        CALL VOCSP( TX,PSOX )
        DO 500 M = 2,ISVC+2   
          CALL WATLQD( TX,PX,RHOLSP(M,NSP) )
          CALL WATLQV( TX,PX,PSWX,VISLSP(M,NSP) )
          CALL VOCLQD( TX,PX,PSOX,RHONSP(M,NSP) )
          CALL VOCLQV( TX,VISNSP(M,NSP),VISCO,IVISO )
  500   CONTINUE  
 1000 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of PROPSP group  ---
!
      RETURN
      END
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORCSP4
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
!
!     Compute source terms for surface spill model.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, September, 1993.
!     Last Modified by Sk Wurstner, PNNL, May 22, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOURC
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
#ifdef awtap
      REAL*8 HCFL(LFZ,LSR),HCFN(LFZ,LSR)
#endif
      REAL*8 SRX(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORCSP4'
!
!---  Load oil sources for surface spill option   ---
!
!---    Zero source terms  ---
!
        IF( NSR.GT.0 ) THEN
          DO 60 N = 1,IJFLD
            DO 50 M = 2,ISVC+2
              SRCWSP(M,N) = 0.D+0
              SRCOSP(M,N) = 0.D+0
   50       CONTINUE
   60     CONTINUE
        ENDIF
!
!---  Loop over sources  ---
!
      DO 600 NS = 1,NSR
        IF( TM.LE.SRC(1,1,NS) ) GOTO 600
        SRX(1) = TM
        IF( ISRM(NS).EQ.1 ) THEN
          DO 70 N = 1,8
            SRX(N) = SRC(N,1,NS)
   70     CONTINUE
        ELSE
          DO 100 M = 2,ISRM(NS)
            IF( TM.LE.SRC(1,M,NS) ) THEN
             DTSR = MIN( SRC(1,M,NS)-TM,DT )
             TFSR = (TM-0.5D+0*DTSR-SRC(1,M-1,NS))/
     &         (SRC(1,M,NS)-SRC(1,M-1,NS))
             DO 80 N = 1,8
               SRX(N) = SRC(N,M-1,NS) + TFSR*(SRC(N,M,NS)-SRC(N,M-1,NS))
   80        CONTINUE
             GOTO 110
            ENDIF
  100     CONTINUE
          GOTO 600
        ENDIF
  110   CONTINUE
!
!---  Loop source domain and translate to spill  ---
!
        DO 500 I = ISRDM(1,NS),ISRDM(2,NS)
          DO 500 J = ISRDM(3,NS),ISRDM(4,NS)
              K = KSPS(I,J)
              N = ND(I,J,K)
              NSP = (J-1)*IFLD + I
              IF( IXP(N).EQ.0 ) GOTO 500
              DO 400 M = 2,ISVC+2
                PGX = PG(M,N) + PATM
!
!---  Aqueous Volumetric Rate w/ Component Mass Fractions ---
!
                IF( ISRT(NS).EQ.41 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                   SRCWSP(M,NSP) = SRCWSP(M,NSP) +
     &               SRX(4)*RHOLSP(M,NSP)*(1.D+0-SRX(6))
                   SRCOSP(M,NSP) = SRCOSP(M,NSP) + SRX(4)*RHOLSP(M,NSP)
     &               *SRX(6)
                  ELSE
                   SRCWSP(M,NSP) = SRCWSP(M,NSP) + SRX(4)*RHOLSP(M,NSP)
     &               *XLW(M,N)
                   SRCOSP(M,NSP) = SRCOSP(M,NSP) + SRX(4)*RHOLSP(M,NSP)
     &               *XLO(M,N)
                  ENDIF
!
!---  NAPL Volumetric Rate  ---
!
                ELSEIF( ISRT(NS).EQ.40 ) THEN
                  IF( SRX(4).GE.0.D+0 ) THEN
                    SRCOSP(M,NSP) = SRCOSP(M,NSP) + SRX(4)
     &                 *RHONSP(M,NSP)
                  ELSE
                    SRCOSP(M,NSP) = SRCOSP(M,NSP) + SRX(4)
     &                 *RHONSP(M,NSP)
                  ENDIF
                ENDIF
  400         CONTINUE
  500   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORCSP4 group  ---
!
      RETURN
      END

!
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE UPDTSP
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
!     Update aqueous and NAPL surface spill heads.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner PNNL, 22 February, 2008.
!     Last Modified by SK Wurstner PNNL, 22 February, 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE JACOB
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
      SUB_LOG(ISUB_LOG) = '/UPDTSP'
!
!---  Update aqueous and NAPL surface spill heads  ---
!
      DO 100 J = 1,JFLD
      DO 100 I = 1,IFLD
        K = KSPS(I,J)
        N = ND(I,J,K)
        NSP = (J-1)*IFLD + I
        NX = NFLD + NSP
        NMD = IXP(NX)
        MPL = IM(IEQW,NMD)
        MPN = IM(IEQO,NMD)
!        IF( SPHMIN(IZ(N)).GT.EPSL ) THEN
!          DHLSPX = 1.D-1*SPHMIN(IZ(N))
!          DHLSPX = SIGN(MIN(DHLSPX,ABS(BLU(MPL))),BLU(MPL))
!        ELSE
          DHLSPX = BLU(MPL)
!        ENDIF
        HLSP(2,NSP) = HLSP(2,NSP)+DHLSPX
!        IF( SPHMIN(IZ(N)).GT.EPSL ) THEN
!          DHNSPX = 1.D-1*SPHMIN(IZ(N))
!          DHNSPX = SIGN(MIN(DHNSPX,ABS(BLU(MPN))),BLU(MPN))
!        ELSE
          DHNSPX = BLU(MPN)
!        ENDIF
        HNSP(2,NSP) = HNSP(2,NSP)+DHNSPX
        IF( HLSP(2,NSP).LT.0.D+0 ) HLSP(2,NSP) = 0.D+0
        IF( HNSP(2,NSP).LT.0.D+0 ) HNSP(2,NSP) = 0.D+0
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of UPDTSP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDSPILL
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
!     Read input file for surface spill information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by SK Wurstner, PNNL, 07 November 2008.
!     Last Modified by SK Wurstner, PNNL, 07 November 2008.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SPILL
      USE SOLTN
      USE PORMED
      USE GRID
      USE FILES
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
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDSPILL'
!
!---  Check for undefined rock/soil types in active nodes  ---
!
      CARD = 'Rock/Soil Zonation Card'
      DO 1 N = 1,NFLD
        IF( IXP(N).EQ.0 ) GOTO 1
        IF( IZ(N).EQ.0 ) THEN
          INDX = 7
          IMSG = N
          CHMSG = 'Undefined Rock/Soil Type @ Node'
          CALL WRMSGS( INDX )
        ENDIF
    1 CONTINUE
!
!---  Write card information to ouput file  ---
!
      CARD = 'Spill Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil mechanical information lines  ---
!
      N = 0
      IJK = 0
   10 CONTINUE
        IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Rock/Soil Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---  IJK, KIJ, or JKI indexing ---
!
        IF( INDEX(ADUM(1:),'indexing').NE.0 ) THEN
          IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Indexing Option Not Declared ' // 
     &        'in Rock/Soil Zonation Card'
            CALL WRMSGS( INDX )
          ENDIF
          IF( INDEX(ADUM,'ijk').NE.0 ) THEN
            IJK = 1
          ELSEIF( INDEX(ADUM,'jki').NE.0 ) THEN
            IJK = 2
          ELSEIF( INDEX(ADUM,'kij').NE.0 ) THEN
            IJK = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          GOTO 220
        ENDIF
!
!---  Search known rock types for a matching type ---
!
        DO 100 M = 1, NROCK
          IF( ADUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 200
          ENDIF
  100   CONTINUE
!
!---  Search known scaling groups for a matching type ---
!
        IF( ISLC(19).EQ.1 ) THEN
          DO 110 M = 1,NSCALE
             IF( ADUM.EQ.SCALNM(M) ) THEN
                ISGRP = M
                IROCK = 1
                GOTO 200
             ENDIF
  110     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
          GOTO 10
        ENDIF
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type: ' // ADUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
  200   CONTINUE
!
!---  Loop over rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
          DO 202 M = IROCK,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              IROCK = M
              GOTO 204
            ENDIF
  202     CONTINUE
        ENDIF
  204   CONTINUE
!
!---    Write rock/soil name  ---
!
        WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
        N = N + 1
  220   CONTINUE
!
!---  Read minimum spill height ---
!
        VARB = 'Minimum Spill Height'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          CALL RDIJK( ISTART,IJK,CHDUM,UNTS,SPHMIN )
        ELSE
          IDFLT = 1
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SPHMIN(IROCK))
          IDFLT = 1
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SPHMIN(IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,SPHMIN(IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SPHMIN(IROCK),', m)'
        ENDIF
!
!---  Loop over remaining rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
          DO 490 M = IROCK+1,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              N = N+1
              SPHMIN(M) = SPHMIN(IROCK)
            ENDIF
  490     CONTINUE
        ENDIF
!
!---  Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDSPILL group ---
!
      RETURN
      END
