!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBIMTS
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
!     Correct Jacobian matrix for interphase mass transfer of
!     VOC between the NAPL and aqueous phase.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on August 19, 1996.
!     jcbimts.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE PORMED
      USE POINTE
      USE NAPL
      USE JACOB
      USE GRID
      USE FLUXP
      USE FDVS
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 XMSS(LCMP),XMOL(LCMP),RS(LUK+1)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBIMTS'
!
!---  Loop over all nodes, skipping inactive nodes  ---
!
      DO 1000 K = 1,KFLD
       DO 1000 J = 1,JFLD
        DO 1000 I = 1,IFLD
        N = ND(I,J,K)
        IZN = IZ(N)
        IF( IXP(N).EQ.0 ) GOTO 1000
        IF( SN(2,N).LT.EPSL ) GOTO 1000
        NPX = NSX(N)
        NPY = NSY(N)
        NPZ = NSZ(N)
        NQX = NPX + 1
        NQY = NPY + IFLD
        NQZ = NPZ + IJFLD
        DO 800 M = 2,ISVC+2
          MN = MNEGB(M-1)
          MP = MPOSB(M-1)
!
!---  Compute the interphase mass transfer rate  ---
!
          ULMX = 5.D-1*( ABS(UL(MP,NPX))+ABS(UL(MN,NQX)) )
          VLMX = 5.D-1*( ABS(VL(MP,NPY))+ABS(VL(MN,NQY)) )
          WLMX = 5.D-1*( ABS(WL(MP,NPZ))+ABS(WL(MN,NQZ)) )
          VLMX = SQRT( ULMX**2 + VLMX**2 + WLMX**2 )
          REYN = VLMX*RHOL(M,N)*PMDD(IZN)/
     &      (VISL(M,N)*PORD(M,N)*(1.D+0-SCHR(4,IZN)))
          SCMN = VISL(M,N)/(DFLO(M,N)*RHOL(M,N))
!
!---  Welty correlation form  ---
!
          IF( IMTC.EQ.1 ) THEN
            SHWN = CIMTC(1) + CIMTC(2)*(REYN**CIMTC(3))*(SCMN**CIMTC(4))
          ELSE
            INDX = 12
            CHMSG = 'Unknown Mass Transfer Rate Correlation: Node: '
            IMSG = N
            CALL WRMSGS( INDX )
          ENDIF
!
!---  Interphase mass transfer flux  ---
!
          CIMTX = SHWN*DFLO(M,N)/(PMDD(IZN)**2)
          PX = MAX( PL(M,N),PG(M,N),PN(M,N) ) + PATM
          CALL SFCSOL( XLS(M,N),T(M,N),PX,PSO(M,N),HCOW,XLOMX )
          RIMTX = PORD(M,N)*SL(M,N)*CIMTX*(XLOMX-XLO(M,N))*RHOL(M,N)
          IF( XLOMX-XLO(M,N).LT.ZERO ) RIMTX = 1.D+3*RIMTX
          RIMTX = MIN( PORD(M,N)*SN(M,N)*RHON(M,N)*DTI,RIMTX )
          RS(M-1) = -RIMTX*VOL(N)
  800   CONTINUE
!
!---  Modify residual for the dissolved VOC equation  ---
!
        CALL JCBBMB( RS,N,IEQD )
!
!---  Modify residual for the oil equation  ---
!
        DO 900 M = 1,ISVC+1
          RS(M) = -RS(M)
  900   CONTINUE
        CALL JCBBMB( RS,N,IEQO )
 1000 CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBIMTS group  ---
!
      RETURN
      END
