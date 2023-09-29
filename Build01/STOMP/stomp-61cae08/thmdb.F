!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDB( N,NB )
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
!     thmdb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDB'
      K = KD(N)
      NPZ = NSZ(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 3
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,
     &      THKL(MP,N),THKG(MP,N),THKN(MP,N),ZERO,
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),ZERO,
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),ZERO,
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,JNDX )
        ENDIF
        WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDE( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDE'
      I = ID(N)
      NQX = NSX(N)+1
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 1
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKN(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),ZERO,THKL(MN,N),
     &      THKG(MN,N),THKN(MN,N),ZERO,PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),ZERO,
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),ZERO,
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,JNDX )
        ENDIF
        UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDE group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDN( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDN'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 2
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKN(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),ZERO,THKL(MN,N),
     &      THKG(MN,N),THKN(MN,N),ZERO,PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),ZERO,
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),ZERO,
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,JNDX )
        ENDIF
        VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDS( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDS'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 2
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,THKL(MP,N),
     &      THKG(MP,N),THKN(MP,N),ZERO,PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),ZERO,
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),ZERO,
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,JNDX )
        ENDIF
        VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDS group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDT( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDT'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 3
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKN(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MN,N),SN(MN,N),ZERO,THKL(MN,N),
     &      THKG(MN,N),THKN(MN,N),ZERO,PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE( IZ(N),SLB(MN,NB),SNB(MN,NB),ZERO,
     &      THKLB(MN,NB),THKGB(MN,NB),THKNB(MN,NB),ZERO,
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,JNDX )
        ENDIF
        WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDW( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE CONST
      USE BCVT
      USE BCVP
      USE BCVI
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THMDW'
      I = ID(N)
      NPX = NSX(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 1
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKN(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
          TKP = THKE( IZ(N),SL(MP,N),SN(MP,N),ZERO,THKL(MP,N),
     &      THKG(MP,N),THKN(MP,N),ZERO,PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE( IZ(N),SLB(MP,NB),SNB(MP,NB),ZERO,
     &      THKLB(MP,NB),THKGB(MP,NB),THKNB(MP,NB),ZERO,
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ENDIF
        JNDX = 1
        TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,JNDX )
        UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDW group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDBL( N,NB )
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
!     thmdb.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDBL'
      K = KD(N)
      NPZ = NSZ(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 3
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &     PORT(MP,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MP,NB),THKLB(MP,NB),PORDB(MP,NB),
     &     PORTB(MP,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKB,TKP,DZGF(N),DZGF(N),DTK,JNDX )
        ENDIF
        WQ(M,NPZ) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDBL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDEL( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDEL'
      I = ID(N)
      NQX = NSX(N)+1
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 1
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MN,N),THKL(MN,N),PORD(MN,N),
     &     PORT(MN,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MN,NB),THKLB(MN,NB),PORDB(MN,NB),
     &     PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DXGF(N),DXGF(N),DTK,JNDX )
        ENDIF
        UQ(M,NQX) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDEL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDNL( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDNL'
      I = ID(N)
      J = JD(N)
      NQY = NSY(N)+IFLD
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 2
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MN,N),THKL(MN,N),PORD(MN,N),
     &     PORT(MN,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MN,NB),THKLB(MN,NB),PORDB(MN,NB),
     &     PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DYGF(N),DYGF(N),DTK,JNDX )
        ENDIF
        VQ(M,NQY) = TK*(T(MN,N)-TB(MN,NB))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDNL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDSL( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDSL'
      I = ID(N)
      J = JD(N)
      NPY = NSY(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 2
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &     PORT(MP,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MP,NB),THKLB(MP,NB),PORDB(MP,NB),
     &     PORTB(MP,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKB,TKP,DYGF(N),DYGF(N),DTK,JNDX )
        ENDIF
        VQ(M,NPY) = TK*(TB(MP,NB)-T(MP,N))/((5.D-1*DYGF(N))*RP(I))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDSL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDTL( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE BCVT
      USE BCVP
      USE BCVI
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
      SUB_LOG(ISUB_LOG) = '/THMDTL'
      K = KD(N)
      NQZ = NSZ(N)+IJFLD
      DTK = T(2,N)-TB(2,NB)
      DO 100 M = 1,ISVF
        MN = MNEG(M)
        INDX = 3
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MN,N),SI(MN,N),
     &      THKL(MN,N),THKG(MN,N),THKI(MN,N),
     &      PORD(MN,N),PORT(MN,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MN,NB),SIB(MN,NB),
     &      THKLB(MN,NB),THKGB(MN,NB),THKIB(MN,NB),
     &      PORDB(MN,NB),PORTB(MN,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MN,N),THKL(MN,N),PORD(MN,N),
     &     PORT(MN,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MN,NB),THKLB(MN,NB),PORDB(MN,NB),
     &     PORTB(MN,NB),INDX )
        ENDIF
        IF( IBCT(IEQT,NB).EQ.2 .OR. IBCT(IEQT,NB).EQ.18 ) THEN
          TK = TKP
        ELSE
          JNDX = 1
          TK = DIFMN( TKP,TKB,DZGF(N),DZGF(N),DTK,JNDX )
        ENDIF
        WQ(M,NQZ) = TK*(T(MN,N)-TB(MN,NB))/(5.D-1*DZGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDTL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE THMDWL( N,NB )
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
      USE POINTE
      USE JACOB
      USE GRID
      USE FLUXT
      USE FDVT
      USE FDVP
      USE FDVI
      USE BCVT
      USE BCVP
      USE BCVI
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/THMDWL'
      I = ID(N)
      NPX = NSX(N)
      DTK = TB(2,NB)-T(2,N)
      DO 100 M = 1,ISVF
        MP = MPOS(M)
        INDX = 1
        IF( ISLC(5).EQ.1 ) THEN
          TKP = THKE_FW( IZ(N),SL(MP,N),SI(MP,N),
     &      THKL(MP,N),THKG(MP,N),THKI(MP,N),
     &      PORD(MP,N),PORT(MP,N),INDX )
          TKB = THKE_FW( IZ(N),SLB(MP,NB),SIB(MP,NB),
     &      THKLB(MP,NB),THKGB(MP,NB),THKIB(MP,NB),
     &      PORDB(MP,NB),PORTB(MP,NB),INDX )
        ELSE
         TKP = THKE_L( IZ(N),SL(MP,N),THKL(MP,N),PORD(MP,N),
     &     PORT(MP,N),INDX )
         TKB = THKE_L( IZ(N),SLB(MP,NB),THKLB(MP,NB),PORDB(MP,NB),
     &     PORTB(MP,NB),INDX )
        ENDIF
        JNDX = 1
        TK = DIFMN( TKB,TKP,DXGF(N),DXGF(N),DTK,JNDX )
        UQ(M,NPX) = TK*(TB(MP,NB)-T(MP,N))/(5.D-1*DXGF(N))
  100 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of THMDWL group  ---
!
      RETURN
      END
