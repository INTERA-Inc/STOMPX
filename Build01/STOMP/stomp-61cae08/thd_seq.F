!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THC_SEQ
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
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, PNNL, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
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
      SUB_LOG(ISUB_LOG) = '/THC_SEQ'
!
!---  X-direction thermal conduction, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DTK = T(2,NW)-T(2,N)
          DO 100 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) + 
     &        PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &        THKG(MP,N)*SG(MP,N) + THKN(MP,N)*SN(MP,N))
            TKW = MAX(1.D+0-PORD(MN,NW),0.D+0)*THKS(1,IZ(NW)) +
     &        PORD(MN,NW)*(THKL(MN,NW)*SL(MN,NW) +
     &        THKG(MN,NW)*SG(MN,NW) + THKN(MN,NW)*SN(MN,NW))
            INDX = 1
            TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,INDX )
            UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DTK = T(2,NS)-T(2,N)
          DO 300 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) +
     &        PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &        THKG(MP,N)*SG(MP,N) + THKN(MP,N)*SN(MP,N))
            TKS = MAX(1.D+0-PORD(MN,NS),0.D+0)*THKS(2,IZ(NS)) +
     &        PORD(MN,NS)*(THKL(MN,NS)*SL(MN,NS) + 
     &        THKG(MN,NS)*SG(MN,NS) + THKN(MN,NS)*SN(MN,NS))
            INDX = 1
            TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,INDX )
            VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/(DYGP(NPY)*RP(I))
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
            IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DTK = T(2,NB)-T(2,N)
          DO 500 M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
            TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) +
     &        PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &        THKG(MP,N)*SG(MP,N) + THKN(MP,N)*SN(MP,N))
            TKB = MAX(1.D+0-PORD(MN,NB),0.D+0)*THKS(3,IZ(NB)) +
     &        PORD(MN,NB)*(THKL(MN,NB)*SL(MN,NB) + 
     &        THKG(MN,NB)*SG(MN,NB) + THKN(MN,NB)*SN(MN,NB))
            INDX = 1
            TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,INDX )
            WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THC_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAG_SEQ
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
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, Battelle, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THAG_SEQ'
!
!---  X-direction gas advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HW = HG(MN,NW)*RHOG(MN,NW)
            UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UG(M,NPX),ZERO)
     &        -HP*MAX(-UG(M,NPX),ZERO)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction gas advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HS = HG(MN,NS)*RHOG(MN,NS)
            VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VG(M,NPY),ZERO)
     &        -HP*MAX(-VG(M,NPY),ZERO)
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction gas advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HB = HG(MN,NB)*RHOG(MN,NB)
            WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &        -HP*MAX(-WG(M,NPZ),ZERO)
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAG_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAN_SEQ
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, Battelle, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THAN_SEQ'
!
!---  X-direction gas advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HW = HN(MN,NW)*RHON(MN,NW)
            UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UN(M,NPX),ZERO)
     &        -HP*MAX(-UN(M,NPX),ZERO)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction gas advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HS = HN(MN,NS)*RHON(MN,NS)
            VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VN(M,NPY),ZERO)
     &        -HP*MAX(-VN(M,NPY),ZERO)
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction gas advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HB = HN(MN,NB)*RHON(MN,NB)
            WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
     &        -HP*MAX(-WN(M,NPZ),ZERO)
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAL_SEQ
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
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, PNNL, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THAL_SEQ'
!
!---  X-direction aqueous advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HW = HL(MN,NW)*RHOL(MN,NW)
            UQ(M,NPX) = UQ(M,NPX) + HW*MAX(UL(M,NPX),ZERO)
     &        -HP*MAX(-UL(M,NPX),ZERO)
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction aqueous advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HS = HL(MN,NS)*RHOL(MN,NS)
            VQ(M,NPY) = VQ(M,NPY) + HS*MAX(VL(M,NPY),ZERO)
     &        -HP*MAX(-VL(M,NPY),ZERO)
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction aqueous advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HB = HL(MN,NB)*RHOL(MN,NB)
            WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &        -HP*MAX(-WL(M,NPZ),ZERO)
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAL_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDG_SEQ
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
!     Compute the contribution to the energy flux by water vapor,
!     voc vapor, and air molecular diffusion for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, PNNL, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THDG_SEQ'
!
!---  X-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            UDGAX = -UDGW(M,NPX)
            UQ(M,NPX) = UQ(M,NPX) 
     &        + HGW(MN,NW)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &        + HGA(MN,NW)*MAX(UDGAX,ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            VDGAX = -VDGW(M,NPY)
            VQ(M,NPY) = VQ(M,NPY) 
     &        + HGW(MN,NS)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &        + HGA(MN,NS)*MAX(VDGAX,ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            WDGAX = -WDGW(M,NPZ)
            WQ(M,NPZ) = WQ(M,NPZ) 
     &        + HGW(MN,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &        + HGA(MN,NB)*MAX(WDGAX,ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDG_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDN_SEQ
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
!     Compute the contribution to the energy flux by water vapor,
!     voc vapor, and air molecular diffusion for nonboundary node faces.
!     Donor cell interfacial averaging.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 15 September 2004.
!     Last Modified by MD White, PNNL, 15 September 2004.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE FLUXT
      USE FLUXN
      USE FDVT
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
      SUB_LOG(ISUB_LOG) = '/THDN_SEQ'
!
!---  X-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( IFLD.GT.1 ) THEN
        DO 200 N = 1,NFLD+NWN_LW
          I = ID(N)
          IF( I.EQ.1 ) GOTO 200
          J = JD(N)
          K = KD(N)
          NW = ICM(1,3,N)
          IF( NW.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
          NPX = NSX(N)
          DO 100 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            UDNAX = -UDNW(M,NPX)
            UQ(M,NPX) = UQ(M,NPX) 
     &        + HGW(MN,NW)*MAX(UDNW(M,NPX),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-UDNW(M,NPX),ZERO)*WTMW
     &        + HGO(MN,NW)*MAX(UDNAX,ZERO)*WTMA
     &        - HGO(MP,N)*MAX(-UDNAX,ZERO)*WTMA
  100     CONTINUE
  200   CONTINUE
      ENDIF
!
!---  Y-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( JFLD.GT.1 ) THEN
        DO 400 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) GOTO 400
          K = KD(N)
          NS = ICM(1,2,N)
          IF( NS.EQ.0 ) CYCLE
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
          NPY = NSY(N)
          DO 300 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            VDNAX = -VDNW(M,NPY)
            VQ(M,NPY) = VQ(M,NPY) 
     &        + HGW(MN,NS)*MAX(VDNW(M,NPY),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-VDNW(M,NPY),ZERO)*WTMW
     &        + HGO(MN,NS)*MAX(VDNAX,ZERO)*WTMA
     &        - HGO(MP,N)*MAX(-VDNAX,ZERO)*WTMA
  300     CONTINUE
  400   CONTINUE
      ENDIF
!
!---  Z-direction water vapor, VOC vapor, and air diffusion,
!     excluding boundaries  ---
!
      IF( KFLD.GT.1 ) THEN
        DO 600 N = 1,NFLD+NWN_LW
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) GOTO 600
          NB = ICM(1,1,N)
          IF( NB.EQ.0 ) CYCLE
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
          NPZ = NSZ(N)
          DO 500 M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            WDNAX = -WDNW(M,NPZ)
            WQ(M,NPZ) = WQ(M,NPZ) 
     &        + HGW(MN,NB)*MAX(WDNW(M,NPZ),ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-WDNW(M,NPZ),ZERO)*WTMW
     &        + HGO(MN,NB)*MAX(WDNAX,ZERO)*WTMA
     &        - HGO(MP,N)*MAX(-WDNAX,ZERO)*WTMA
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGB_SEQ( N,NB )
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
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
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
      SUB_LOG(ISUB_LOG) = '/THAGB_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WG(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &      -HP*MAX(-WG(M,NPZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGE_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGE_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
     &      -HB*MAX(-UG(M,NQX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGN_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGN_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
     &      -HB*MAX(-VG(M,NQY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGS_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGS_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VG(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
     &      -HP*MAX(-VG(M,NPY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGT_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGT_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
     &      - HB*MAX(-WG(M,NQZ),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGW_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGW_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UG(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
     &      -HP*MAX(-UG(M,NPX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THAGW_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANB_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a bottom boundary
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
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
      SUB_LOG(ISUB_LOG) = '/THANB_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WN(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
     &      -HP*MAX(-WN(M,NPZ),ZERO)
        ENDIF 
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANE_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THANE_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
     &      -HB*MAX(-UN(M,NQX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANN_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THANN_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
     &      -HB*MAX(-VN(M,NQY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANS_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THANS_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VN(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VN(M,NPY),ZERO)
     &      -HP*MAX(-VN(M,NPY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANT_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THANT_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
     &      - HB*MAX(-WN(M,NQZ),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THANW_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by nonaqueous
!     liquid advection on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THANW_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UN(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UN(M,NPX),ZERO)
     &      -HP*MAX(-UN(M,NPX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THANW_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALB_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THALB_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WL(M,NPZ),ZERO)
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &      -HP*MAX(-WL(M,NPZ),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALE_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THALE_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
     &      - HB*MAX(-UL(M,NQX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALN_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THALN_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
     &      -HB*MAX(-VL(M,NQY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALS_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THALS_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VL(M,NPY),ZERO)
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
     &      -HP*MAX(-VL(M,NPY),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALT_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THALT_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
     &      - HB*MAX(-WL(M,NQZ),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALW_SEQ( N,NB )
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THALW_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
        IF( IBCT(IEQT,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UL(M,NPX),ZERO)
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
     &      - HP*MAX(-UL(M,NPX),ZERO)
        ENDIF
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THALW_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGB_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGB_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        WDGAX = -WDGW(M,NPZ)
        WQ(M,NPZ) = WQ(M,NPZ) 
     &    + HGWB(MP,NB)*MAX(WDGW(M,NPZ),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-WDGW(M,NPZ),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(WDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-WDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGE_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGE_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        UDGAX = -UDGW(M,NQX)
        UQ(M,NQX) = UQ(M,NQX) 
     &    - HGWB(MN,NB)*MAX(-UDGW(M,NQX),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(UDGW(M,NQX),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-UDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(UDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGN_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGN_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        VDGAX = -VDGW(M,NQY)
        VQ(M,NQY) = VQ(M,NQY) 
     &    - HGWB(MN,NB)*MAX(-VDGW(M,NQY),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(VDGW(M,NQY),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-VDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(VDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGS_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGS_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        VDGAX = -VDGW(M,NPY)
        VQ(M,NPY) = VQ(M,NPY) 
     &    + HGWB(MP,NB)*MAX(VDGW(M,NPY),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-VDGW(M,NPY),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(VDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-VDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGT_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGT_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        WDGAX = -WDGW(M,NQZ)
        WQ(M,NQZ) = WQ(M,NQZ) 
     &    - HGWB(MN,NB)*MAX(-WDGW(M,NQZ),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(WDGW(M,NQZ),ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-WDGAX,ZERO)*WTMA
     &    + HGA(MN,N)*MAX(WDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGW_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDGW_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        UDGAX = -UDGW(M,NPX)
        UQ(M,NPX) = UQ(M,NPX) 
     &    + HGWB(MP,NB)*MAX(UDGW(M,NPX),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-UDGW(M,NPX),ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(UDGAX,ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-UDGAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGW_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNB_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, April 8, 1994.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNB_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        WDNAX = -WDNW(M,NPZ)
        WQ(M,NPZ) = WQ(M,NPZ) 
     &    + HGWB(MP,NB)*MAX(WDNW(M,NPZ),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-WDNW(M,NPZ),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(WDNAX,ZERO)*WTMA
     &    - HGO(MP,N)*MAX(-WDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNE_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNE_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        UDNAX = -UDNW(M,NQX)
        UQ(M,NQX) = UQ(M,NQX) 
     &    - HGWB(MN,NB)*MAX(-UDNW(M,NQX),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(UDNW(M,NQX),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-UDNAX,ZERO)*WTMA
     &    + HGO(MN,N)*MAX(UDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNN_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNN_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        VDNAX = -VDNW(M,NQY)
        VQ(M,NQY) = VQ(M,NQY) 
     &    - HGWB(MN,NB)*MAX(-VDNW(M,NQY),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(VDNW(M,NQY),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-VDNAX,ZERO)*WTMA
     &    + HGO(MN,N)*MAX(VDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNS_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNS_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        VDNAX = -VDNW(M,NPY)
        VQ(M,NPY) = VQ(M,NPY) 
     &    + HGWB(MP,NB)*MAX(VDNW(M,NPY),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-VDNW(M,NPY),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(VDNAX,ZERO)*WTMA
     &    - HGO(MP,N)*MAX(-VDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNT_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNT_SEQ'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        WDNAX = -WDNW(M,NQZ)
        WQ(M,NQZ) = WQ(M,NQZ) 
     &    - HGWB(MN,NB)*MAX(-WDNW(M,NQZ),ZERO)*WTMW
     &    + HGW(MN,N)*MAX(WDNW(M,NQZ),ZERO)*WTMW
     &    - HGOB(MN,NB)*MAX(-WDNAX,ZERO)*WTMA
     &    + HGO(MN,N)*MAX(WDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDNW_SEQ( N,NB )
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
!     Compute the contribution to the energy flux by water vapor,
!     and air molecular diffusion on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, February 11, 1993.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE CONST
      USE BCVT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THDNW_SEQ'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        UDNAX = -UDNW(M,NPX)
        UQ(M,NPX) = UQ(M,NPX) 
     &    + HGWB(MP,NB)*MAX(UDNW(M,NPX),ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-UDNW(M,NPX),ZERO)*WTMW
     &    + HGOB(MP,NB)*MAX(UDNAX,ZERO)*WTMA
     &    - HGO(MP,N)*MAX(-UDNAX,ZERO)*WTMA
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDNW_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCB_SEQ( N,NB )
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
!     on a bottom boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, April 8, 1994.
!     thd_co2e.F 1227 2018-07-04 14:50:35Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCB_SEQ'
      K = KD(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) 
     &    + PORD(MP,N)*(THKL(MP,N)*SL(MP,N) 
     &    + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
        TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &    + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
        INDX = 1
        NPZ = NSZ(N)
        TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,INDX )
        WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCB_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCE_SEQ( N,NB )
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
!     on a east boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 8, 1995.
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCE_SEQ'
      I = ID(N)
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(1,IZ(N)) 
     &    + PORD(MN,N)*(THKL(MN,N)*SL(MN,N) 
     &    + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
        TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &    + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
        INDX = 1
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,INDX )
        UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCE_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCN_SEQ( N,NB )
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
!     on a north boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 8, 1995.
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCN_SEQ'
      I = ID(N)
      J = JD(N)
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(2,IZ(N)) 
     &    + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &    + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
        TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &    + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
        INDX = 1
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,INDX )
        VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCN_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCS_SEQ( N,NB )
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
!     on a south boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 8, 1995.
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCS_SEQ'
      I = ID(N)
      J = JD(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) 
     &    + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &    + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
        TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &    + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
        INDX = 1
        NPY = NSY(N)
        TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,INDX )
        VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCS_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCT_SEQ( N,NB )
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
!     on a top boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 8, 1995.
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCT_SEQ'
      K = KD(N)
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(3,IZ(N)) 
     &    + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &    + THKN(MN,N)*SN(MN,N)  + THKG(MN,N)*SG(MN,N))
        TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &    + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
        INDX = 1
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,INDX )
        WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCT_SEQ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THCW_SEQ( N,NB )
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
!     on a west boundary
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1993.
!     Last Modified by MD White, PNNL, December 8, 1995.
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
      USE BCVT
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
      SUB_LOG(ISUB_LOG) = '/THCW_SEQ'
      I = ID(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) 
     &    + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &    + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
        TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(3,IZ(N))
     &    + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &    + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
        INDX = 1
        NPX = NSX(N)
        TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,INDX )
        UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THCW_SEQ group  ---
!
      RETURN
      END

