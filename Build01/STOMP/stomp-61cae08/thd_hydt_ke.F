!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THD_HYDT_KE
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THD_HYDT_KE'
!
!---  X-direction thermal conduction, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DTK = T(2,NW)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) + 
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(1,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
              TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))
     &          *(THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &          PORD(MP,N)*THKI(MP,N)*SI(MP,N)
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NW)).EQ.2 ) THEN
              TKW = MAX(1.D+0-PORD(MN,NW),0.D+0)*THKS(1,IZ(NW)) +
     &          PORD(MN,NW)*(THKL(MN,NW)*SL(MN,NW) +
     &          THKH(MN,NW)*SH(MN,NW) + THKI(MN,NW)*SI(MN,NW) + 
     &          THKN(MN,NW)*SN(MN,NW) + THKG(MN,NW)*SG(MN,NW))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NW)).EQ.4 ) THEN
              TKW = THKS(1,IZ(NW)) + 
     &          SQRT(SL(MN,NW))*(THKS(4,IZ(NW))-THKS(1,IZ(NW)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(NW)).EQ.8 ) THEN
              TKW = THKS(1,IZ(NW)) + (SQRT(SL(MN,NW)) + SQRT(SH(MN,NW)))
     &          *(THKS(2,IZ(NW))-THKS(1,IZ(NW))) +
     &          PORD(MN,NW)*THKI(MN,NW)*SI(MN,NW)
            ENDIF
            INDX = 1
            TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,INDX )
            UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DTK = T(2,NS)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(2,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
              TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))
     &          *(THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &          PORD(MP,N)*THKI(MP,N)*SI(MP,N)
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NS)).EQ.2 ) THEN
              TKS = MAX(1.D+0-PORD(MN,NS),0.D+0)*THKS(2,IZ(NS)) +
     &          PORD(MN,NS)*(THKL(MN,NS)*SL(MN,NS) + 
     &          THKH(MN,NS)*SH(MN,NS) + THKI(MN,NS)*SI(MN,NS) + 
     &          THKN(MN,NS)*SN(MN,NS) + THKG(MN,NS)*SG(MN,NS))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NS)).EQ.4 ) THEN
              TKS = THKS(2,IZ(NS)) + 
     &          SQRT(SL(MN,NS))*(THKS(5,IZ(NS))-THKS(2,IZ(NS)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(NS)).EQ.8 ) THEN
              TKS = THKS(1,IZ(NS)) + (SQRT(SL(MN,NS)) + SQRT(SH(MN,NS)))
     &          *(THKS(2,IZ(NS))-THKS(1,IZ(NS))) +
     &          PORD(MN,NS)*THKI(MN,NS)*SI(MN,NS)
            ENDIF
            INDX = 1
            TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,INDX )
            VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/(DYGP(NPY)*RP(I))
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DTK = T(2,NB)-T(2,N)
          DO M = 1,ISVF
            MN = MNEG(M)
            MP = MPOS(M)
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(N)).EQ.2 ) THEN
              TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) +
     &          PORD(MP,N)*(THKL(MP,N)*SL(MP,N) + 
     &          THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) + 
     &          THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
              TKP = THKS(3,IZ(N)) + 
     &          SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
              TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))*
     &          (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &          PORD(MP,N)*THKI(MP,N)*SI(MP,N)
            ENDIF
!
!---        Parallel function  ---
!
            IF( ITHK(IZ(NB)).EQ.2 ) THEN
              TKB = MAX(1.D+0-PORD(MN,NB),0.D+0)*THKS(3,IZ(NB)) +
     &          PORD(MN,NB)*(THKL(MN,NB)*SL(MN,NB) + 
     &          THKH(MN,NB)*SH(MN,NB) + THKI(MN,NB)*SI(MN,NB) + 
     &          THKN(MN,NB)*SN(MN,NB) + THKG(MN,NB)*SG(MN,NB))
!
!---        Somerton function  ---
!
            ELSEIF( ITHK(IZ(NB)).EQ.4 ) THEN
              TKB = THKS(3,IZ(NB)) + 
     &          SQRT(SL(MN,NB))*(THKS(6,IZ(NB))-THKS(3,IZ(NB)))
!
!---        Composite function  ---
!
            ELSEIF( ITHK(IZ(NB)).EQ.8 ) THEN
              TKB = THKS(1,IZ(NB)) + (SQRT(SL(MN,NB)) + SQRT(SH(MN,NB)))
     &          *(THKS(2,IZ(NB))-THKS(1,IZ(NB))) +
     &          PORD(MN,NB)*THKI(MN,NB)*SI(MN,NB)
            ENDIF
            INDX = 1
            TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,INDX )
            WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THD_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THG_HYDT_KE
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
!     Written by MD White, PNNL, 24 October 2012.
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THG_HYDT_KE'
!
!---  X-direction gas advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HW = HG(MN,NW)*RHOG(MN,NW)
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UG(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UG(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction gas advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HS = HG(MN,NS)*RHOG(MN,NS)
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VG(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VG(M,NPY)
          ENDDO
        ENDDO    
      ENDIF
!
!---  Z-direction gas advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HG(MP,N)*RHOG(MP,N)
            HB = HG(MN,NB)*RHOG(MN,NB)
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WG(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WG(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THG_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THL_HYDT_KE
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
!     Written by MD White, PNNL, 24 October 2012.
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THL_HYDT_KE'
!
!---  X-direction aqueous advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HW = HL(MN,NW)*RHOL(MN,NW)
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UL(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UL(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction aqueous advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HS = HL(MN,NS)*RHOL(MN,NS)
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VL(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VL(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction aqueous advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HL(MP,N)*RHOL(MP,N)
            HB = HL(MN,NB)*RHOL(MN,NB)
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WL(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WL(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THL_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THN_HYDT_KE
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
!     Written by MD White, PNNL, 24 October 2012.
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
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THN_HYDT_KE'
!
!---  X-direction gas advection, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          IF( I.EQ.1 ) CYCLE
          J = JD(N)
          K = KD(N)
          NW = N-1
          IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &      INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) CYCLE
          NPX = NSX(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HW = HN(MN,NW)*RHON(MN,NW)
            INDX = -4
            HM = DIFMN( HW,HP,DXGF(NW),DXGF(N),UN(1,NPX),INDX )
            UQ(M,NPX) = UQ(M,NPX) + HM*UN(M,NPX)
          ENDDO
        ENDDO
      ENDIF
!
!---  Y-direction gas advection, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          IF( J.EQ.1 ) CYCLE
          K = KD(N)
          NS = N-IFLD
          IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &      INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) CYCLE
          NPY = NSY(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HS = HN(MN,NS)*RHON(MN,NS)
            INDX = -4
            HM = DIFMN( HS,HP,DYGF(NS),DYGF(N),VN(1,NPY),INDX )
            VQ(M,NPY) = VQ(M,NPY) + HM*VN(M,NPY)
          ENDDO
        ENDDO
      ENDIF
!
!---  Z-direction gas advection, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
        DO N = 1,NFLD
          I = ID(N)
          J = JD(N)
          K = KD(N)
          IF( K.EQ.1 ) CYCLE
          NB = N-IJFLD
          IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &      INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) CYCLE
          NPZ = NSZ(N)
          DO M = 1,ISVF
            MP = MPOS(M)
            MN = MNEG(M)
            HP = HN(MP,N)*RHON(MP,N)
            HB = HN(MN,NB)*RHON(MN,NB)
            INDX = -4
            HM = DIFMN( HB,HP,DZGF(NB),DZGF(N),WN(1,NPZ),INDX )
            WQ(M,NPZ) = WQ(M,NPZ) + HM*WN(M,NPZ)
          ENDDO
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THN_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGB_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGB_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WG(M,NPZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WG(M,NPZ),ZERO)
     &      - HP*MAX(-WG(M,NPZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGB_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGE_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGE_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UG(M,NQX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UG(M,NQX),ZERO)
     &      - HB*MAX(-UG(M,NQX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGE_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGN_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGN_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VG(M,NQY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VG(M,NQY),ZERO)
     &      - HB*MAX(-VG(M,NQY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGN_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGS_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGS_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VG(M,NPY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VG(M,NPY),ZERO)
     &      - HP*MAX(-VG(M,NPY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGS_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGT_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGT_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HG(MN,N)*RHOG(MN,N)
        HB = HGB(MN,NB)*RHOGB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WG(M,NQZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WG(M,NQZ),ZERO)
     &      - HB*MAX(-WG(M,NQZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGT_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THGW_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THGW_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HG(MP,N)*RHOG(MP,N)
        HB = HGB(MP,NB)*RHOGB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UG(M,NPX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UG(M,NPX),ZERO)
     &      - HP*MAX(-UG(M,NPX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THGW_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLB_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLB_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WL(M,NPZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WL(M,NPZ),ZERO)
     &      - HP*MAX(-WL(M,NPZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLB_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLE_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLE_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UL(M,NQX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UL(M,NQX),ZERO)
     &      - HB*MAX(-UL(M,NQX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLE_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLN_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLN_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VL(M,NQY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VL(M,NQY),ZERO)
     &      - HB*MAX(-VL(M,NQY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLN_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLS_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLS_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VL(M,NPY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VL(M,NPY),ZERO)
     &      - HP*MAX(-VL(M,NPY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLS_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLT_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLT_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HL(MN,N)*RHOL(MN,N)
        HB = HLB(MN,NB)*RHOLB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WL(M,NQZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WL(M,NQZ),ZERO)
     &      - HB*MAX(-WL(M,NQZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLT_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THLW_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THLB_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HL(MP,N)*RHOL(MP,N)
        HB = HLB(MP,NB)*RHOLB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UL(M,NPX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UL(M,NPX),ZERO)
     &      - HP*MAX(-UL(M,NPX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THLW_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNB_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNB_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPZ = NSZ(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) - HP*MAX(-WN(M,NPZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NPZ) = WQ(M,NPZ) + HB*MAX(WN(M,NPZ),ZERO)
     &      - HP*MAX(-WN(M,NPZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNB_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNE_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNE_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQX = NSX(N)+1
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NQX) = UQ(M,NQX) - HB*MAX(-UN(M,NQX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NQX) = UQ(M,NQX) + HP*MAX(UN(M,NQX),ZERO)
     &      - HB*MAX(-UN(M,NQX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNE_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNN_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNN_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQY = NSY(N)+IFLD
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NQY) = VQ(M,NQY) - HB*MAX(-VN(M,NQY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NQY) = VQ(M,NQY) + HP*MAX(VN(M,NQY),ZERO)
     &      - HB*MAX(-VN(M,NQY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNN_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNS_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNS_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPY = NSY(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          VQ(M,NPY) = VQ(M,NPY) - HP*MAX(-VN(M,NPY),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VN(M,NPY),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          VQ(M,NPY) = VQ(M,NPY) + HB*MAX(VN(M,NPY),ZERO)
     &      - HP*MAX(-VN(M,NPY),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNS_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNT_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNT_HYDT_KE'
      DO M = 1,ISVF
        MN = MNEG(M)
        NQZ = NSZ(N)+IJFLD
        HP = HN(MN,N)*RHON(MN,N)
        HB = HNB(MN,NB)*RHONB(MN,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          WQ(M,NQZ) = WQ(M,NQZ) - HB*MAX(-WN(M,NQZ),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          WQ(M,NQZ) = WQ(M,NQZ) + HP*MAX(WN(M,NQZ),ZERO)
     &      - HB*MAX(-WN(M,NQZ),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNT_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THNW_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      SUB_LOG(ISUB_LOG) = '/THNW_HYDT_KE'
      DO M = 1,ISVF
        MP = MPOS(M)
        NPX = NSX(N)
        HP = HN(MP,N)*RHON(MP,N)
        HB = HNB(MP,NB)*RHONB(MP,NB)
!
!---    Energy outflow boundary  ---
!
        IF( IBCT(1,NB).EQ.7 ) THEN
          UQ(M,NPX) = UQ(M,NPX) - HP*MAX(-UN(M,NPX),ZERO)
!
!---    Energy inflow boundary  ---
!
        ELSEIF( IBCT(1,NB).EQ.13 ) THEN
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UN(M,NPX),ZERO)
!
!---    Energy boundary  ---
!
        ELSE
          UQ(M,NPX) = UQ(M,NPX) + HB*MAX(UN(M,NPX),ZERO)
     &      - HP*MAX(-UN(M,NPX),ZERO)
        ENDIF
      ENDDO
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THNW_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDB_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDB_HYDT_KE'
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N) 
     &      + THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPZ = NSZ(N)
          WQ(M,NPZ) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DZGF(N)
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NPZ = NSZ(N)
          WQ(M,NPZ) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DZGF(N)
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MP,N)*THKI(MP,N)*SI(MP,N)
          INDX = 1
          NPZ = NSZ(N)
          WQ(M,NPZ) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DZGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDB_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDE_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDE_HYDT_KE'
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N) 
     &      + THKH(MN,N)*SH(MN,N) + THKI(MN,N)*SI(MN,N) 
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQX = NSX(N)+1
          UQ(M,NQX) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DXGF(N)
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(4,IZ(N))-THKS(1,IZ(N)))
          INDX = 1
          NQX = NSX(N)+1
          UQ(M,NQX) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DXGF(N)
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MN,N)) + SQRT(SH(MN,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MN,N)*THKI(MN,N)*SI(MN,N)
          INDX = 1
          NQX = NSX(N)+1
          UQ(M,NQX) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DXGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDE_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDN_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDN_HYDT_KE'
      I = ID(N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKH(MN,N)*SH(MN,N) + THKI(MN,N)*SI(MN,N) 
     &      + THKN(MN,N)*SN(MN,N) + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQY = NSY(N)+IFLD
          VQ(M,NQY) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/((DYGF(N))*RP(I))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NQY = NSY(N)+IFLD
          VQ(M,NQY) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/((DYGF(N))*RP(I))
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MN,N)) + SQRT(SH(MN,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MN,N)*THKI(MN,N)*SI(MN,N)
          INDX = 1
          INDX = 1
          NQY = NSY(N)+IFLD
          VQ(M,NQY) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/((DYGF(N))*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDN_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDS_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDS_HYDT_KE'
      I = ID(N)
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(2,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPY = NSY(N)
          VQ(M,NPY) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/((DYGF(N))*RP(I))
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(2,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(5,IZ(N))-THKS(2,IZ(N)))
          INDX = 1
          NPY = NSY(N)
          VQ(M,NPY) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/((DYGF(N))*RP(I))
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MP,N)*THKI(MP,N)*SI(MP,N)
          INDX = 1
          NPY = NSY(N)
          VQ(M,NPY) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/((DYGF(N))*RP(I))
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDS_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDT_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDT_HYDT_KE'
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = MAX(1.D+0-PORD(MN,N),0.D+0)*THKS(3,IZ(N)) 
     &      + PORD(MN,N)*(THKL(MN,N)*SL(MN,N)
     &      + THKH(MN,N)*SH(MN,N) + THKI(MN,N)*SI(MN,N) 
     &      + THKN(MN,N)*SN(MN,N)  + THKG(MN,N)*SG(MN,N))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          WQ(M,NQZ) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DZGF(N)
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MN,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          WQ(M,NQZ) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DZGF(N)
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MN = MNEG(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MN,N)) + SQRT(SH(MN,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MN,N)*THKI(MN,N)*SI(MN,N)
          INDX = 1
          NQZ = NSZ(N)+IJFLD
          WQ(M,NQZ) = 2.D+0*TKP*(T(MN,N)-TB(MN,NB))/DZGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDT_HYDT_KE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THDW_HYDT_KE( N,NB )
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
!     Written by MD White, PNNL, 24 October 2012.
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
      USE FDVI
      USE FDVH
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
      SUB_LOG(ISUB_LOG) = '/THDW_HYDT_KE'
!
!---  Parallel function  ---
!
      IF( ITHK(IZ(N)).EQ.2 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = MAX(1.D+0-PORD(MP,N),0.D+0)*THKS(1,IZ(N)) 
     &      + PORD(MP,N)*(THKL(MP,N)*SL(MP,N)
     &      + THKH(MP,N)*SH(MP,N) + THKI(MP,N)*SI(MP,N) 
     &      + THKN(MP,N)*SN(MP,N) + THKG(MP,N)*SG(MP,N))
          INDX = 1
          NPX = NSX(N)
          UQ(M,NPX) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DXGF(N)
        ENDDO
!
!---  Somerton function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.4 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(3,IZ(N)) 
     &      + SQRT(SL(MP,N))*(THKS(6,IZ(N))-THKS(3,IZ(N)))
          INDX = 1
          NPX = NSX(N)
          UQ(M,NPX) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DXGF(N)
        ENDDO
!
!---  Composite function  ---
!
      ELSEIF( ITHK(IZ(N)).EQ.8 ) THEN
        DO M = 1,ISVF
          MP = MPOS(M)
          TKP = THKS(1,IZ(N)) + (SQRT(SL(MP,N)) + SQRT(SH(MP,N)))*
     &      (THKS(2,IZ(N))-THKS(1,IZ(N))) +
     &      PORD(MP,N)*THKI(MP,N)*SI(MP,N)
          INDX = 1
          NPX = NSX(N)
          UQ(M,NPX) = 2.D+0*TKP*(TB(MP,NB)-T(MP,N))/DXGF(N)
        ENDDO
      ENDIF
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THDW_HYDT_KE group  ---
!
      RETURN
      END

