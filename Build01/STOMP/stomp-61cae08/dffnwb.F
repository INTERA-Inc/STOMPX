!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWB( N,NB )
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
!     dffnwb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFNWB'
      K = KD(N)
      DXMNW = XMNWB(2,NB)-XMNW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNW(MP,N)
     &   *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &   *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMNW,INDX )
       WDGW(M,NSZ(N)) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWE( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/DFFNWE'
      I = ID(N)
      DXMNW = XMNW(2,N)-XMNWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNW(MN,N)
     &   *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNWB(MN,NB)
     &   *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMNW,INDX )
       UDGW(M,NSX(N)+1) = 2.D+0*DFM*(XMNW(MN,N)-XMNWB(MN,NB))
     &  /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWN( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/DFFNWN'
      I = ID(N)
      J = JD(N)
      DXMNW = XMNW(2,N)-XMNWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNW(MN,N)
     &  *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)* DFNWB(MN,NB)
     &   *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMNW,INDX )
       VDGW(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMNW(MN,N)- XMNWB(MN,NB))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWS( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/DFFNWS'
      I = ID(N)
      J = JD(N)
      DXMNW = XMNWB(2,NB)-XMNW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNW(MP,N)
     &   *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &   *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMNW,INDX )
       VDGW(M,NSY(N)) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWT( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/DFFNWT'
      K = KD(N)
      DXMNW = XMNW(2,N)-XMNWB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORN(MN,N)*PORD(MN,N)*(SN(MN,N)-SNT(MN,N))*DFNW(MN,N)
     &   *RHOMN(MN,N)
       DFB = TORNB(MN,NB)*PORDB(MN,NB)*SNB(MN,NB)*DFNWB(MN,NB)
     &   *RHOMNB(MN,NB)
       INDX = 12
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMNW,INDX )
       WDGW(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMNW(MN,N)-XMNWB(MN,NB))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFNWW( N,NB )
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
      SUB_LOG(ISUB_LOG) = '/DFFNWW'
      I = ID(N)
      DXMNW = XMNWB(2,NB)-XMNW(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORN(MP,N)*PORD(MP,N)*(SN(MP,N)-SNT(MP,N))*DFNW(MP,N)
     &   *RHOMN(MP,N)
       DFB = TORNB(MP,NB)*PORDB(MP,NB)*SNB(MP,NB)*DFNWB(MP,NB)
     &   *RHOMNB(MP,NB)
       INDX = 12
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMNW,INDX )
       UDGW(M,NSX(N)) = 2.D+0*DFM*(XMNWB(MP,NB)-XMNW(MP,N))
     &   /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFNWW group  ---
!
      RETURN
      END
