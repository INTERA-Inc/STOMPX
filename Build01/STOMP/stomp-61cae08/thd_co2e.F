!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD33
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
!     thd_co2e.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THD33'
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
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) + 
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
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
     &          THKG(MN,NW)*SG(MN,NW))
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
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
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
     &          THKG(MN,NS)*SG(MN,NS))
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
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN          
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKG(MP,N)*SG(MP,N))
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
     &          THKG(MN,NB)*SG(MN,NB))
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
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THD33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAG33
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
!     thd_co2e.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THAG33'
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
!---  End of THAG33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAL33
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
!     thd_co2e.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THAL33'
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
!---  End of THAL33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDG33
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
!     thd_co2e.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/THDG33'
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
            UDGWX = -UDGA(M,NPX)
            UQ(M,NPX) = UQ(M,NPX) + HGW(MN,NW)*MAX(UDGWX,ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-UDGWX,ZERO)*WTMW
     &        + HGA(MN,NW)*MAX(UDGA(M,NPX),ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-UDGA(M,NPX),ZERO)*WTMA
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
            VDGWX = -VDGA(M,NPY)
            VQ(M,NPY) = VQ(M,NPY) + HGW(MN,NS)*MAX(VDGWX,ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-VDGWX,ZERO)*WTMW
     &        + HGA(MN,NS)*MAX(VDGA(M,NPY),ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-VDGA(M,NPY),ZERO)*WTMA
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
            WDGWX = -WDGA(M,NPZ)
            WQ(M,NPZ) = WQ(M,NPZ) + HGW(MN,NB)*MAX(WDGWX,ZERO)*WTMW
     &        - HGW(MP,N)*MAX(-WDGWX,ZERO)*WTMW
     &        + HGA(MN,NB)*MAX(WDGA(M,NPZ),ZERO)*WTMA
     &        - HGA(MP,N)*MAX(-WDGA(M,NPZ),ZERO)*WTMA
  500     CONTINUE
  600   CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDG33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGB33( N,NB )
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
!     thd_co2e.F 1336 2020-06-23 14:31:14Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
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
      SUB_LOG(ISUB_LOG) = '/THAGB33'
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
!---  End of THAGB33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGE33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THAGE33'
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
!---  End of THAGE33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGN33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THAGN33'
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
!---  End of THAGN33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGS33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THAGS33'
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
!---  End of THAGS33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGT33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THAGT33'
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
!---  End of THAGT33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THAGW33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THAGW33'
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
!---  End of THAGW33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALB33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALB33'
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
!---  End of THALB33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALE33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALE33'
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
!---  End of THALE33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALN33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALN33'
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
!---  End of THALN33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALS33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALS33'
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
!---  End of THALS33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALT33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALT33'
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
!---  End of THALT33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THALW33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THALW33'
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
!---  End of THALW33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGB33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGB33'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        WDGWX = -WDGA(M,NPZ)-WDGO(M,NPZ)
        WQ(M,NPZ) = WQ(M,NPZ) 
     &    + HGWB(MP,NB)*MAX(WDGWX,ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-WDGWX,ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(WDGA(M,NPZ),ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-WDGA(M,NPZ),ZERO)*WTMA
     &    + HGOB(MP,NB)*MAX(WDGO(M,NPZ),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-WDGO(M,NPZ),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGB33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGE33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGE33'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
        UDGWX = -UDGA(M,NQX)-UDGO(M,NQX)
        UQ(M,NQX) = UQ(M,NQX) 
     &    - HGWB(MN,NB)*MAX(-UDGWX,ZERO)*WTMW
     &    + HGW(MN,N)*MAX(UDGWX,ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-UDGA(M,NQX),ZERO)*WTMA
     &    + HGA(MN,N)*MAX(UDGA(M,NQX),ZERO)*WTMA
     &    - HGOB(MN,NB)*MAX(-UDGO(M,NQX),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(UDGO(M,NQX),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGE33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGN33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGN33'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
        VDGWX = -VDGA(M,NQY)-VDGO(M,NQY)
        VQ(M,NQY) = VQ(M,NQY) 
     &    - HGWB(MN,NB)*MAX(-VDGWX,ZERO)*WTMW
     &    + HGW(MN,N)*MAX(VDGWX,ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-VDGA(M,NQY),ZERO)*WTMA
     &    + HGA(MN,N)*MAX(VDGA(M,NQY),ZERO)*WTMA
     &    - HGOB(MN,NB)*MAX(-VDGO(M,NQY),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(VDGO(M,NQY),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGN33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGS33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGS33'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        VDGWX = -VDGA(M,NPY)-VDGO(M,NPY)
        VQ(M,NPY) = VQ(M,NPY) 
     &    + HGWB(MP,NB)*MAX(VDGWX,ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-VDGWX,ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(VDGA(M,NPY),ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-VDGA(M,NPY),ZERO)*WTMA
     &    + HGOB(MP,NB)*MAX(VDGO(M,NPY),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-VDGO(M,NPY),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGS33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGT33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGT33'
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
        WDGWX = -WDGA(M,NQZ)-WDGO(M,NQZ)
        WQ(M,NQZ) = WQ(M,NQZ) 
     &    - HGWB(MN,NB)*MAX(-WDGWX,ZERO)*WTMW
     &    + HGW(MN,N)*MAX(WDGWX,ZERO)*WTMW
     &    - HGAB(MN,NB)*MAX(-WDGA(M,NQZ),ZERO)*WTMA
     &    + HGA(MN,N)*MAX(WDGA(M,NQZ),ZERO)*WTMA
     &    - HGOB(MN,NB)*MAX(-WDGO(M,NQZ),ZERO)*WTMO
     &    + HGO(MN,N)*MAX(WDGO(M,NQZ),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGT33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDGW33( N,NB )
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
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE NAPL
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
      SUB_LOG(ISUB_LOG) = '/THDGW33'
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        UDGWX = -UDGA(M,NPX)-UDGO(M,NPX)
        UQ(M,NPX) = UQ(M,NPX) 
     &    + HGWB(MP,NB)*MAX(UDGWX,ZERO)*WTMW
     &    - HGW(MP,N)*MAX(-UDGWX,ZERO)*WTMW
     &    + HGAB(MP,NB)*MAX(UDGA(M,NPX),ZERO)*WTMA
     &    - HGA(MP,N)*MAX(-UDGA(M,NPX),ZERO)*WTMA
     &    + HGOB(MP,NB)*MAX(UDGO(M,NPX),ZERO)*WTMO
     &    - HGO(MP,N)*MAX(-UDGO(M,NPX),ZERO)*WTMO
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDGW33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDB33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDB33'
      K = KD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(3,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
          INDX = 1
          NPZ = NSZ(N)
          IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))     
          TKB = THKS(3,IZ(N))
     &      + SQRT(SLB(MP,NB))*(THKS(6,IZ(N))-THKS(3,IZ(N)))     
          INDX = 1
          NPZ = NSZ(N)
          IF( INBS(1,N).GT.0 ) NPZ = INBS(1,N)
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDB33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDE33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDE33'
      I = ID(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N) 
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(1,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQX = NSX(N)+1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          TKB = THKS(1,IZ(N))
     &      + SQRT(SLB(MN,NB))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NQX = NSX(N)+1
          IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
        ENDDO     
      ENDIF 
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDE33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDN33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDN33'
      I = ID(N)
      J = JD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(2,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQY = NSY(N)+IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(2,IZ(N)) 
     &        + SQRT(SL(MN,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          TKB = THKS(2,IZ(N))
     &        + SQRT(SLB(MN,NB))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NQY = NSY(N)+IFLD
          IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDN33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDS33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDS33'
      I = ID(N)
      J = JD(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(2,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
        INDX = 1
        NPY = NSY(N)
        TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,INDX )
        VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
      ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          TKB = THKS(2,IZ(N))
     &      + SQRT(SL(MP,NB))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NPY = NSY(N)
          TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,INDX )
          VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDS33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDT33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDT33'
      K = KD(N)
      DTK = T(2,N)-TB(2,NB)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKN(MN,N)*SN(MN,N)  + THKG(MN,N)*SG(MN,N))
          TKB = MAX(1.D+0-PORDB(MN,NB),0.D+0)*THKS(3,IZ(N))
     &      + PORDB(MN,NB)*(THKLB(MN,NB)*SLB(MN,NB)
     &      + THKNB(MN,NB)*SNB(MN,NB) + THKGB(MN,NB)*SGB(MN,NB))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          TKB = THKS(3,IZ(N))
     &      + SQRT(SL(MN,NB))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,INDX )
          WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDT33 group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDW33( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/THDW33'
      I = ID(N)
      DTK = TB(2,NB)-T(2,N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          TKB = MAX(1.D+0-PORDB(MP,NB),0.D+0)*THKS(1,IZ(N))
     &      + PORDB(MP,NB)*(THKLB(MP,NB)*SLB(MP,NB)
     &      + THKNB(MP,NB)*SNB(MP,NB) + THKGB(MP,NB)*SGB(MP,NB))
          INDX = 1
          NPX = NSX(N)
          TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          TKB = THKS(1,IZ(N))
     &      + SQRT(SL(MP,NB))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NPX = NSX(N)
          TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,INDX )
          UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
        ENDDO
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDW33 group  ---
!
      RETURN
      END

