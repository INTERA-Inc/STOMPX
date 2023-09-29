!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLA
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
!     Compute dissolved air molar diffusion rates through the
!     aqueous phase.
!     dffla.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, February, 1994.
!     Last Modified by MD White, PNNL, October 27, 1997.
!     dffla.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
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
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(4).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFLA'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 200 K = 1,KFLD
        DO 200 J = 1,JFLD
          DO 200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 200
            NPX = NSX(N)
            DXLA = (XMLA(2,NW)-XMLA(2,N))
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &          *RHOML(MP,N)
              DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
     &          *RHOML(MN,NW)
              INDX = 14
              DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
              UDLA(M,NPX) = DFM*(XMLA(MN,NW)-XMLA(MP,N))/DXGP(NPX)
  100       CONTINUE
  200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 400 I = 1,IFLD
        DO 400 K = 1,KFLD
          DO 400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 400
            NPY = NSY(N)
            DXLA = (XMLA(2,NS)-XMLA(2,N))
            DO 300 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &          *RHOML(MP,N)
              DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
     &          *RHOML(MN,NS)
              INDX = 14
              DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
              VDLA(M,NPY) = DFM*(XMLA(MN,NS)-XMLA(MP,N))/DYGP(NPY)/RP(I)
  300       CONTINUE
  400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 600 J = 1,JFLD
        DO 600 I = 1,IFLD
          DO 600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 600
            NPZ = NSZ(N)
            DXLA = (XMLA(2,NB)-XMLA(2,N))
            DO 500 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &          *RHOML(MP,N)
              DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
     &          *RHOML(MN,NB)
              INDX = 14
              DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
              WDLA(M,NPZ) = DFM*(XMLA(MN,NB)-XMLA(MP,N))/DZGP(NPZ)
  500       CONTINUE
  600 CONTINUE
      ENDIF
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
!
!---  X-direction molar diffusion, excluding boundaries
!
      IF( IFLD.GT.1 ) THEN
      DO 1200 K = 1,KFLD
        DO 1200 J = 1,JFLD
          DO 1200 I = 2,IFLD
            N = ND(I,J,K)
            NW = N-1
            IF( IXP(N).EQ.0 .OR. IXP(NW).EQ.0 .OR.
     &        INBS(3,N).GT.0 .OR. INBS(4,NW).GT.0 ) GOTO 1200
            NPX = NSX(N)
            DXLA = (XMLA(2,NW)*RHOML(2,NW)-XMLA(2,N)*RHOML(2,N))
            DO 1100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
              DFW = TORL(MN,NW)*PORD(MN,NW)*SL(MN,NW)*DFLA(MN,NW)
              INDX = 14
              DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXLA,INDX)
              UDLA(M,NPX) = DFM*(XMLA(MN,NW)*RHOML(MN,NW)
     &         - XMLA(MP,N)*RHOML(MP,N))/DXGP(NPX)
 1100       CONTINUE
 1200 CONTINUE
      ENDIF
!
!---  Y-direction molar diffusion, excluding boundaries
!
      IF( JFLD.GT.1 ) THEN
      DO 1400 I = 1,IFLD
        DO 1400 K = 1,KFLD
          DO 1400 J = 2,JFLD
            N = ND(I,J,K)
            NS = N-IFLD
            IF( IXP(N).EQ.0 .OR. IXP(NS).EQ.0 .OR.
     &        INBS(2,N).GT.0 .OR. INBS(5,NS).GT.0 ) GOTO 1400
            NPY = NSY(N)
            DXLA = (XMLA(2,NS)*RHOML(2,NS)-XMLA(2,N)*RHOML(2,N))
            DO 1300 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
              DFS = TORL(MN,NS)*PORD(MN,NS)*SL(MN,NS)*DFLA(MN,NS)
              INDX = 14
              DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXLA,INDX )
              VDLA(M,NPY) = DFM*(XMLA(MN,NS)*RHOML(MN,NS)
     &          - XMLA(MP,N)*RHOML(MP,N))/DYGP(NPY)/RP(I)
 1300       CONTINUE
 1400 CONTINUE
      ENDIF
!
!---  Z-direction molar diffusion, excluding boundaries
!
      IF( KFLD.GT.1 ) THEN
      DO 1600 J = 1,JFLD
        DO 1600 I = 1,IFLD
          DO 1600 K = 2,KFLD
            N = ND(I,J,K)
            NB = N-IJFLD
            IF( IXP(N).EQ.0 .OR. IXP(NB).EQ.0 .OR.
     &        INBS(1,N).GT.0 .OR. INBS(6,NB).GT.0 ) GOTO 1600
            NPZ = NSZ(N)
            DXLA = (XMLA(2,NB)*RHOML(2,NB)-XMLA(2,N)*RHOML(2,N))
            DO 1500 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
              DFB = TORL(MN,NB)*PORD(MN,NB)*SL(MN,NB)*DFLA(MN,NB)
              INDX = 14
              DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXLA,INDX)
              WDLA(M,NPZ) = DFM*(XMLA(MN,NB)*RHOML(MN,NB)
     &          - XMLA(MP,N)*RHOML(MP,N))/DZGP(NPZ)
 1500       CONTINUE
 1600 CONTINUE
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLA group  ---
!
      RETURN
      END
