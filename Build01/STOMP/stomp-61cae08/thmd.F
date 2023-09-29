!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMD
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, March 23, 2000.
!     thmd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
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
      SUB_LOG(ISUB_LOG) = '/THMD'
!
!---  X-direction thermal conduction, excluding boundaries
!
      INDX = 1
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
      DO 200 J = 1,JFLD
      DO 200 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
        NPX = NSX(N)
        DTK = T(2,NW)-T(2,N)
        DO 100 M = 1,ISVF
           MN = MNEG(M)
           MP = MPOS(M)
           IF( ISLC(5).EQ.1 ) THEN
             TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &         THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &         PORD(MP,N),PORT(MP,N),INDX )
             TKW = THKE( IZ(NW),SL(MN,NW),SN(MN,NW),SI(MN,NW),
     &         THKL(MN,NW),THKG(MN,NW),THKN(MN,NW),THKI(MN,NW),
     &         PORD(MN,NW),PORT(MN,NW),INDX )
           ELSE
             TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,THKL(MP,N),
     &         THKG(MP,N),THKN(MP,N),ZERO,PORD(MP,N),PORT(MP,N),INDX )
             TKW = THKE( IZ(NW),SL(MN,NW),SN(MN,NW),ZERO,
     &         THKL(MN,NW),THKG(MN,NW),THKN(MN,NW),ZERO,
     &         PORD(MN,NW),PORT(MN,NW),INDX )
           ENDIF
           JNDX = 1
           TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,JNDX )
           UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
  100   CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
!
      INDX = 2
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
      DO 400 K = 1,KFLD
      DO 400 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
        NPY = NSY(N)
        DTK = T(2,NS)-T(2,N)
        DO 300 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          IF( ISLC(5).EQ.1 ) THEN
            TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &        THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &        PORD(MP,N),PORT(MP,N),INDX )
            TKS = THKE( IZ(NS),SL(MN,NS),SN(MN,NS),SI(MN,NS),
     &        THKL(MN,NS),THKG(MN,NS),THKN(MN,NS),THKI(MN,NS),
     &        PORD(MN,NS),PORT(MN,NS),INDX )
          ELSE
            TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,
     &        THKL(MP,N),THKG(MP,N),THKN(MP,N),ZERO,
     &        PORD(MP,N),PORT(MP,N),INDX )
            TKS = THKE( IZ(NS),SL(MN,NS),SN(MN,NS),ZERO,
     &        THKL(MN,NS),THKG(MN,NS),THKN(MN,NS),ZERO,
     &        PORD(MN,NS),PORT(MN,NS),INDX )
          ENDIF
          JNDX = 1
          TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,JNDX )
          VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/DYGP(NPY)/RP(I)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
!
      INDX = 3
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
      DO 600 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
        NPZ = NSZ(N)
        DTK = T(2,NB)-T(2,N)
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          IF( ISLC(5).EQ.1 ) THEN
            TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &        THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &        PORD(MP,N),PORT(MP,N),INDX )
            TKB = THKE( IZ(NB),SL(MN,NB),SN(MN,NB),SI(MN,NB),
     &        THKL(MN,NB),THKG(MN,NB),THKN(MN,NB),THKI(MN,NB),
     &        PORD(MN,NB),PORT(MN,NB),INDX )
          ELSE
            TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,THKL(MP,N),
     &        THKG(MP,N),THKN(MP,N),ZERO,PORD(MP,N),PORT(MP,N),INDX )
            TKB = THKE( IZ(NB),SL(MN,NB),SN(MN,NB),ZERO,
     &        THKL(MN,NB),THKG(MN,NB),THKN(MN,NB),ZERO,
     &        PORD(MN,NB),PORT(MN,NB),INDX )
          ENDIF
          JNDX = 1
          TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,JNDX )
          WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMD group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDL
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, March 23, 2000.
!     thmd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
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
      SUB_LOG(ISUB_LOG) = '/THMDL'
!
!---  X-direction thermal conduction, excluding boundaries
!
      INDX = 1
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
      DO 200 J = 1,JFLD
      DO 200 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
        NPX = NSX(N)
        DTK = T(2,NW)-T(2,N)
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &      PORT(MP,N),INDX )
          TKW = THKE_L( IZ(NW),SL(MN,NW),THKL(MN,NW),PORD(MN,NW),
     &      PORT(MN,NW),INDX )
          JNDX = 1
          TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,JNDX )
          UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
  100   CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
!
      INDX = 2
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
      DO 400 K = 1,KFLD
      DO 400 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
        NPY = NSY(N)
        DTK = T(2,NS)-T(2,N)
        DO 300 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &      PORT(MP,N),INDX )
          TKS = THKE_L( IZ(NS),SL(MN,NS),THKL(MN,NS),PORD(MN,NS),
     &      PORT(MN,NS),INDX )
          JNDX = 1
          TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,JNDX )
          VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/DYGP(NPY)/RP(I)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
!
      INDX = 3
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
      DO 600 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
        NPZ = NSZ(N)
        DTK = T(2,NB)-T(2,N)
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &      PORT(MP,N),INDX )
          TKB = THKE_L( IZ(NB),SL(MN,NB),THKL(MN,NB),PORD(MN,NB),
     &      PORT(MN,NB),INDX )
          JNDX = 1
          TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,JNDX )
          WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDFW
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
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, March 23, 2000.
!     thmd.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THMDFW'
!
!---  X-direction thermal conduction, excluding boundaries
!
      INDX = 1
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
      DO 200 J = 1,JFLD
      DO 200 I = 2,IFLD
        N = ND(I,J,K)
        NW = N-1
        IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &    INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
        NPX = NSX(N)
        DTK = T(2,NW)-T(2,N)
        DO 100 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),THKL(MP,N),
     &      THKG(MP,N),THKI(MP,N),PORD(MP,N),PORT(MP,N),INDX )
          TKW = THKE_FW( IZ(NW),SL(MN,NW),SI(MN,NW),THKL(MN,NW),
     &      THKG(MN,NW),THKI(MN,NW),PORD(MN,NW),PORT(MN,NW),INDX )
          JNDX = 1
          TK = DIFMN( TKW,TKP,DXGF(NW),DXGF(N),DTK,JNDX )
          UQ(M,NPX) = TK*(T(MN,NW)-T(MP,N))/DXGP(NPX)
  100   CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction thermal conduction, excluding boundaries
!
      INDX = 2
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
      DO 400 K = 1,KFLD
      DO 400 J = 2,JFLD
        N = ND(I,J,K)
        NS = N-IFLD
        IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &    INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
        NPY = NSY(N)
        DTK = T(2,NS)-T(2,N)
        DO 300 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),THKL(MP,N),
     &      THKG(MP,N),THKI(MP,N),PORD(MP,N),PORT(MP,N),INDX )
          TKS = THKE_FW( IZ(NS),SL(MN,NS),SI(MN,NS),THKL(MN,NS),
     &      THKG(MN,NS),THKI(MN,NS),PORD(MN,NS),PORT(MN,NS),INDX )
          JNDX = 1
          TK = DIFMN( TKS,TKP,DYGF(NS),DYGF(N),DTK,JNDX )
          VQ(M,NPY) = TK*(T(MN,NS)-T(MP,N))/DYGP(NPY)/RP(I)
  300   CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction thermal conduction, excluding boundaries
!
      INDX = 3
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
      DO 600 I = 1,IFLD
      DO 600 K = 2,KFLD
        N = ND(I,J,K)
        NB = N-IJFLD
        IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &    INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
        NPZ = NSZ(N)
        DTK = T(2,NB)-T(2,N)
        DO 500 M = 1,ISVF
          MN = MNEG(M)
          MP = MPOS(M)
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),THKL(MP,N),
     &      THKG(MP,N),THKI(MP,N),PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE_FW( IZ(NB),SL(MN,NB),SI(MN,NB),THKL(MN,NB),
     &      THKG(MN,NB),THKI(MN,NB),PORD(MN,NB),PORT(MN,NB),INDX )
          JNDX = 1
          TK = DIFMN( TKB,TKP,DZGF(NB),DZGF(N),DTK,JNDX )
          WQ(M,NPZ) = TK*(T(MN,NB)-T(MP,N))/DZGP(NPZ)
  500   CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDFW group  ---
!
      RETURN
      END

