!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD_EOR
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
!     Compute the contribution to the energy flux by thermal conduction
!     for nonboundary node faces.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXT
      USE FDVT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THD_EOR'
!
!---  X-direction thermal conduction, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DTK = T(2,NW)-T(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) + 
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(1,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NW)).EQ.2 ) THEN
              TKW = MAX(1.D+0-PORD(MN,NW),0.D+0)*THKS(1,IZ(NW)) +
     &          PORD(MN,NW)*(THKL(MN,NW)*SL(MN,NW) +
     &          THKN(MN,NW)*SN(MN,NW) + THKG(MN,NW)*SG(MN,NW))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NW)).EQ.4 ) THEN
              TKW = THKS(1,IZ(NW)) + 
     &          SQRT(SL(MN,NW))*(THKS(4,IZ(NW))-THKS(1,IZ(NW)))
            ENDIF
            INDX = 1
            TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,INDX )
            UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction thermal conduction, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DTK = T(2,NS)-T(2,N)
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(2,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NS)).EQ.2 ) THEN
              TKS = MAX(1.D+0-PORD(MN,NS),0.D+0)*THKS(2,IZ(NS)) +
     &          PORD(MN,NS)*(THKL(MN,NS)*SL(MN,NS) + 
     &          THKN(MN,NS)*SN(MN,NS) + THKG(MN,NS)*SG(MN,NS))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NS)).EQ.4 ) THEN
              TKS = THKS(2,IZ(NS)) + 
     &          SQRT(SL(MN,NS))*(THKS(5,IZ(NS))-THKS(2,IZ(NS)))
            ENDIF
            INDX = 1
            TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,INDX )
            VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/(DYGP(NPY)*RP(I))
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction thermal conduction, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DTK = T(2,NB)-T(2,N)
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(3,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NB)).EQ.2 ) THEN
              TKB = MAX(1.D+0-PORD(MN,NB),0.D+0)*THKS(3,IZ(NB)) +
     &          PORD(MN,NB)*(THKL(MN,NB)*SL(MN,NB) + 
     &          THKN(MN,NB)*SN(MN,NB) + THKG(MN,NB)*SG(MN,NB))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NB)).EQ.4 ) THEN
              TKB = THKS(3,IZ(NB)) + 
     &          SQRT(SL(MN,NB))*(THKS(6,IZ(NB))-THKS(3,IZ(NB)))
            ENDIF
            INDX = 1
            TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,INDX )
            WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THD_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THG_EOR
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
!     Compute the contribution to the energy flux by gas advection
!     for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THG_EOR'
!
!---  X-direction gas advection, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HW = HG(MN,NW)*RHOG(MN,NW)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UG(M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction gas advection, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HS = HG(MN,NS)*RHOG(MN,NS)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VG(M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction gas advection, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HB = HG(MN,NB)*RHOG(MN,NB)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WG(M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THG_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THL_EOR
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
!     Compute the contribution to the energy flux by aqueous advection
!     for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THL_EOR'
!
!---  X-direction aqueous advection, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HW = HL(MN,NW)*RHOL(MN,NW)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UL(M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction aqueous advection, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HS = HL(MN,NS)*RHOL(MN,NS)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VL(M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction aqueous advection, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HB = HL(MN,NB)*RHOL(MN,NB)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WL(M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THL_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THN_EOR
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
!     Compute the contribution to the energy flux by nonaqueous-liquid
!     advection for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE GLB_PAR
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THN_EOR'
!
!---  X-direction gas advection, excluding boundaries
!
      DO 200 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 200
        DO 150 NC = 1,4
          NW = ICM(NC,3,N)
          IF( NW.EQ.0 ) EXIT
!
!---      Multiple west connections  ---
!
          IF( INBS(3,N).LT.0 ) THEN
            NPX = NSSX(NW)
          ELSE
            NPX = NSX(N)
          ENDIF
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HW = HN(MN,NW)*RHON(MN,NW)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UN(M,NPX)
  100     CONTINUE
  150   CONTINUE
  200 CONTINUE
!
!---  Y-direction gas advection, excluding boundaries
!
      DO 400 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 400
        DO 350 NC = 1,4
          NS = ICM(NC,2,N)
          IF( NS.EQ.0 ) EXIT
          I = ID(N)
!
!---      Multiple south connections  ---
!
          IF( INBS(2,N).LT.0 ) THEN
            NPY = NSSY(NS)
          ELSE
            NPY = NSY(N)
          ENDIF
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HS = HN(MN,NS)*RHON(MN,NS)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VN(M,NPY)
  300     CONTINUE
  350   CONTINUE
  400 CONTINUE
!
!---  Z-direction gas advection, excluding boundaries
!
      DO 600 N = 1,NFBN
        IF( IXP(N).EQ.0 .OR. IBR(4,N).NE.N ) GOTO 600
        DO 550 NC = 1,4
          NB = ICM(NC,1,N)
          IF( NB.EQ.0 ) EXIT
!
!---      Multiple bottom connections  ---
!
          IF( INBS(1,N).LT.0 ) THEN
            NPZ = NSSZ(NB)
          ELSE
            NPZ = NSZ(N)
          ENDIF
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HB = HN(MN,NB)*RHON(MN,NB)
!
!---        Upwind mean  ---
!
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WN(M,NPZ)
  500     CONTINUE
  550   CONTINUE
  600 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGB_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGB_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WG(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &      - HP*MAX(-WG(M,NPZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGE_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGE_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UG(M,NQX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
     &      - HB*MAX(-UG(M,NQX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGN_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGN_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VG(M,NQY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
     &      - HB*MAX(-VG(M,NQY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGS_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGS_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VG(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
     &      - HP*MAX(-VG(M,NPY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGT_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGT_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WG(M,NQZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
     &      - HB*MAX(-WG(M,NQZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGT_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGW_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THGW_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UG(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
     &      - HP*MAX(-UG(M,NPX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLB_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLB_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WL(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &      - HP*MAX(-WL(M,NPZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLE_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLE_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UL(M,NQX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
     &      - HB*MAX(-UL(M,NQX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLN_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLN_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VL(M,NQY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
     &      - HB*MAX(-VL(M,NQY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLS_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLS_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VL(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
     &      - HP*MAX(-VL(M,NPY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLT_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLT_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WL(M,NQZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
     &      - HB*MAX(-WL(M,NQZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLT_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLW_EOR( N,NB )
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
!     Compute the contribution to the energy flux by aqueous advection
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THLW_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UL(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
     &      - HP*MAX(-UL(M,NPX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNB_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNB_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WN(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
     &      - HP*MAX(-WN(M,NPZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNE_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNE_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UN(M,NQX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
     &      - HB*MAX(-UN(M,NQX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNN_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNN_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VN(M,NQY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
     &      - HB*MAX(-VN(M,NQY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNS_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNS_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VN(M,NPY),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VN(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VN(M,NPY),ZERO)
     &      - HP*MAX(-VN(M,NPY),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNT_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNT_EOR'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WN(M,NQZ),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
     &      - HB*MAX(-WN(M,NQZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNT_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNW_EOR( N,NB )
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
!     Compute the contribution to the energy flux by gas advection
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
      USE FDVP
      USE CONST
      USE BCVT
      USE BCVP
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THNW_EOR'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy inflow boundary condition  ---
!
        IF( IBCT(IEQT,NB).EQ.5 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UN(M,NPX),ZERO)
!
!---    Energy outflow boundary condition  ---
!
        ELSEIF( IBCT(IEQT,NB).EQ.6 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UN(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UN(M,NPX),ZERO)
     &      - HP*MAX(-UN(M,NPX),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNW_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDB_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a bottom boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDB_EOR'
      K = KD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPZ = NSZ(N)
          WQ(M,NPZ) = TKP*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NPZ = NSZ(N)
          WQ(M,NPZ) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DZGF(N)
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDB_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDE_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a east boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDE_EOR'
      I = ID(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N) 
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQX = NSX(N)+1
          UQ(M,NQX) = TKP*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NQX = NSX(N)+1
          UQ(M,NQX) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DXGF(N)
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDE_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDN_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a north boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDN_EOR'
      I = ID(N)
      J = JD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQY = NSY(N)+IFLD
          VQ(M,NQY) = TKP*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NQY = NSY(N)+IFLD
          VQ(M,NQY) = TKP*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDS_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a south boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDS_EOR'
      I = ID(N)
      J = JD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPY = NSY(N)
          VQ(M,NPY) = TKP*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NPY = NSY(N)
          VQ(M,NPY) = TKP*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDS_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDT_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a top boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDT_EOR'
      K = KD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N)  + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          WQ(M,NQZ) = TKP*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          WQ(M,NQZ) = TKP*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDT_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDW_EOR( N,NB )
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
!     Compute the contribution to the energy flux by thermal conduction
!     on a west boundary, using the effective thermal conductivity
!     of the adjacent node.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 29 August 2013
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE BCVP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDW_EOR'
      I = ID(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO 100 M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPX = NSX(N)
          UQ(M,NPX) = TKP*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
  100   CONTINUE
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO 200 M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NPX = NSX(N)
          UQ(M,NPX) = TKP*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
  200   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDW_EOR group  ---
!
      RETURN
      END

