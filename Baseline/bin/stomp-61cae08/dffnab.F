!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAB( N,NB )
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
!     Compute water-NAPL mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
!     dffnab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAB'
      K = KD(N)
      DXMNA = XMNAB(2,NB)-XMNA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNA(MP,N)
     &   *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNAB(MP,NB)
     &   *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNA,INDX )
       WDGO(M,NSZ(N)) = 2.D+0*DFM*(XMNAB(MP,NB)-XMNA(MP,N))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAE( N,NB )
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
!     Compute water-NAPL mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAE'
      I = ID(N)
      DXMNA = XMNA(2,N)-XMNAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNA(MN,N)
     &   *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNAB(MN,NB)
     &   *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNA,INDX )
       UDGO(M,NSX(N)+1) = 2.D+0*DFM*(XMNA(MN,N)-XMNAB(MN,NB))
     &   /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAN( N,NB )
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
!     Compute water-NAPL mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAN'
      I = ID(N)
      J = JD(N)
      DXMNA = XMNA(2,N)-XMNAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNA(MN,N)
     &   *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNAB(MN,NB)
     &   *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNA,INDX )
       VDGO(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMNA(MN,N)-XMNAB(MN,NB))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAS( N,NB )
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
!     Compute water-NAPL mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAS'
      I = ID(N)
      J = JD(N)
      DXMNA = XMNAB(2,NB)-XMNA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNA(MP,N)
     &   *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNAB(MP,NB)
     &   *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNA,INDX )
       VDGO(M,NSY(N)) = 2.D+0*DFM*(XMNAB(MP,NB)-XMNA(MP,N))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAT( N,NB )
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
!     Compute water-NAPL mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAT'
      K = KD(N)
      DXMNA = XMNA(2,N)-XMNAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNA(MN,N)
     &  *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNAB(MN,NB)
     &  *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNA,INDX )
       WDGO(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMNA(MN,N)-XMNAB(MN,NB))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNAW( N,NB )
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
!     Compute water-NAPL mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNNL, October, 2000.
!     Last Modified by MD White, Battelle, PNNL, October 6, 2000.
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
      USE BCVP
      USE BCVN
      USE BCVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      IF( ISLC(2).LT.1 ) RETURN
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/DFFNAW'
      I = ID(N)
      DXMNA = XMNAB(2,NB)-XMNA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNA(MP,N)
     &  *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNAB(MP,NB)
     &  *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNA,INDX )
       UDGO(M,NSX(N)) = 2.D+0*DFM*(XMNAB(MP,NB)-XMNA(MP,N))
     &   /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNAW group  ---
!
      RETURN
      END
