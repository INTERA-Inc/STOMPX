!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOB( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
!     dfflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLOB'
      K = KD(N)
      DXMLO = XMLOB(2,NB)-XMLO(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 15
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLO,INDX )
       WDLO(M,NSZ(N)) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOE( N,NB )
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
!     Compute VOC vapor mole diffusion rates on an east boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLOE'
      I = ID(N)
      DXMLO = XMLO(2,N)-XMLOB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 15
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLO,INDX )
       UDLO(M,NSX(N)+1) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))
     &   /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLON( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a north boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLON'
      I = ID(N)
      J = JD(N)
      DXMLO = XMLO(2,N)-XMLOB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 15
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLO,INDX )
       VDLO(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLON group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOS( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a south boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLOS'
      I = ID(N)
      J = JD(N)
      DXMLO = XMLOB(2,NB)-XMLO(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 15
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLO,INDX )
       VDLO(M,NSY(N)) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOT( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a top boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLOT'
      K = KD(N)
      DXMLO = XMLO(2,N)-XMLOB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLO(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLOB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 15
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLO,INDX )
       WDLO(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLO(MN,N)-XMLOB(MN,NB))
     &   /DZGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLOW( N,NB )
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
!     Compute VOC vapor mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
!     dfflob.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVP
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
      SUB_LOG(ISUB_LOG) = '/DFFLOW'
      I = ID(N)
      DXMLO = XMLOB(2,NB)-XMLO(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLO(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLOB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 15
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLO,INDX )
       UDLO(M,NSX(N)) = 2.D+0*DFM*(XMLOB(MP,NB)-XMLO(MP,N))
     &   /DXGF(N)
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLOW group  ---
!
      RETURN
      END
