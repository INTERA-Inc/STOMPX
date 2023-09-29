!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAB( N,NB )
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
!     Compute dissolved air mole diffusion rates on a bottom boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
!     dfflab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLAB'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NSZ(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))
     &  /DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NSZ(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   - XMLA(MP,N)*RHOML(MP,N))/DZGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAE( N,NB )
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
!     Compute dissolved air mole diffusion rates on an east boundary.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAE'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &  *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NSX(N)+1) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &  /DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NSX(N)+1) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/DXGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAN( N,NB )
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
!     Compute dissolved air mole diffusion rates on a north boundary.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAN'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLA(MN,N)-XMLAB(MN,NB))
     &   /(DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NSY(N)+IFLD) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAS( N,NB )
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
!     Compute dissolved air mole diffusion rates on a south boundary.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAS'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      J = JD(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NSY(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/
     &   (DYGF(N)*RP(I))
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      J = JD(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DYGF(N),DYGF(N),DXMLA,INDX )
       VDLA(M,NSY(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   -XMLA(MP,N)*RHOML(MP,N))/(DYGF(N)*RP(I))
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAT( N,NB )
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
!     Compute dissolved air mole diffusion rates on a top boundary.
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
      SUB_LOG(ISUB_LOG) = '/DFFLAT'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      K = KD(N)
      DXMLA = XMLA(2,N)-XMLAB(2,NB)
      DO 100 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
     &   *RHOML(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
     &   *RHOMLB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLA(MN,N)-
     &   XMLAB(MN,NB))/DZGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      K = KD(N)
      DXMLA = XMLA(2,N)*RHOML(2,N)-XMLAB(2,NB)*RHOMLB(2,NB)
      DO 200 M = 1,ISVF
       MN = MNEG(M)
       DFP = TORL(MN,N)*PORD(MN,N)*SL(MN,N)*DFLA(MN,N)
       DFB = TORLB(MN,NB)*PORDB(MN,NB)*SLB(MN,NB)*DFLAB(MN,NB)
       INDX = 14
       DFM = DIFMN( DFP,DFB,DZGF(N),DZGF(N),DXMLA,INDX )
       WDLA(M,NSZ(N)+IJFLD) = 2.D+0*DFM*(XMLA(MN,N)*RHOML(MN,N)
     &   - XMLAB(MN,NB)*RHOMLB(MN,NB))/DZGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE DFFLAW( N,NB )
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
!     Compute dissolved air mole diffusion rates on a west boundary.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!     Last Modified by MD White, Battelle, PNL, October 27, 1997.
!     dfflab.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      SUB_LOG(ISUB_LOG) = '/DFFLAW'
!
!---  Mole-fraction diffusion gradient option ---
!
      IF( ISLC(27).EQ.1 ) THEN
      I = ID(N)
      DXMLA = XMLAB(2,NB)-XMLA(2,N)
      DO 100 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
     &   *RHOML(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
     &   *RHOMLB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NSX(N)) = 2.D+0*DFM*(XMLAB(MP,NB)-XMLA(MP,N))/DXGF(N)
  100 CONTINUE
!
!---  Molar-density diffusion gradient option ---
!
      ELSE
      I = ID(N)
      DXMLA = XMLAB(2,NB)*RHOMLB(2,NB)-XMLA(2,N)*RHOML(2,N)
      DO 200 M = 1,ISVF
       MP = MPOS(M)
       DFP = TORL(MP,N)*PORD(MP,N)*SL(MP,N)*DFLA(MP,N)
       DFB = TORLB(MP,NB)*PORDB(MP,NB)*SLB(MP,NB)*DFLAB(MP,NB)
       INDX = 14
       DFM = DIFMN( DFB,DFP,DXGF(N),DXGF(N),DXMLA,INDX )
       UDLA(M,NSX(N)) = 2.D+0*DFM*(XMLAB(MP,NB)*RHOMLB(MP,NB)
     &   - XMLA(MP,N)*RHOML(MP,N))/DXGF(N)
  200 CONTINUE
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of DFFLAW group  ---
!
      RETURN
      END
