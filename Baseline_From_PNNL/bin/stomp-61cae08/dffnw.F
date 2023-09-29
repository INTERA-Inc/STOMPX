!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNW
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
!     Compute dissolved water molar diffusion rates through the NAPL.
!     dffnw.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, October, 2000.
!     Last Modified by MD White, PNNL, October 4, 2000.
!     dffnw.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE POINTE
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXP
      USE FDVP
      USE FDVN
      USE FDVA
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
      SUB_LOG(ISUB_LOG) = '/DFFNW'
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
            DXNW = XMNW(2,NW)-XMNW(2,N)
            DO 100 M = 1,ISVF
              MN = MNEG(M)
              MP = MPOS(M)
              DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*
     &          DFNW(MP,N)*RHOMN(MP,N)
              DFW = TORN(MN,NW)*PORD(MN,NW)*(SN(MN,NW)-SNT(MN,NW))*
     &          DFNW(MN,NW)*RHOMN(MN,NW)
              INDX = 14
              DFM = DIFMN( DFW,DFP,DXGF(NW),DXGF(N),DXNW,INDX)
              UDGW(M,NPX) = DFM*(XMNW(MN,NW)-XMNW(MP,N))
     &          /DXGP(NPX)
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
            DXNW = XMNW(2,NS)-XMNW(2,N)
            DO 300 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*
     &          DFNW(MP,N)*RHOMN(MP,N)
              DFS = TORN(MN,NS)*PORD(MN,NS)*(SN(MN,NS)-SNT(MN,NS))*
     &          DFNW(MN,NS)*RHOMN(MN,NS)
              INDX = 14
              DFM = DIFMN( DFS,DFP,DYGF(NS),DYGF(N),DXNW,INDX )
              VDGW(M,NPY) = DFM*(XMNW(MN,NS)-XMNW(MP,N))
     &          /DYGP(NPY)/RP(I)
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
            DXNW = XMNW(2,NB)-XMNW(2,N)
            DO 500 M = 1,ISVF
              MP = MPOS(M)
              MN = MNEG(M)
              DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*
     &          DFNW(MP,N)*RHOMN(MP,N)
              DFB = TORN(MN,NB)*PORD(MN,NB)*(SN(MN,NB)-SNT(MN,NB))*
     &          DFNW(MN,NB)*RHOMN(MN,NB)
              INDX = 14
              DFM = DIFMN( DFB,DFP,DZGF(NB),DZGF(N),DXNW,INDX)
              WDGW(M,NPZ) = DFM*(XMNW(MN,NB)-XMNW(MP,N))
     &          /DZGP(NPZ)
  500       CONTINUE
  600 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNW group  ---
!
      RETURN
      END
