!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE BRNLQV( TX,XLSX,VISLX )
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
!     Compute aqueous-brine (sodium chloride) liquid viscosity.
!
!     Leijnse, A.  1992.  Three-Dimensional Modeling of Coupled Flow
!     and Transport in Porous Media.  Ph.D. Dissertation, Department
!     of Civil Engineering and Geological Sciences, University of
!     Notre Dame, Notre Dame, Indiana.
!
!     Mahiuddin, S., and K. Ismail.  1996.  Temperature and
!     concentration dependence of the viscosity of aqueous sodium
!     nitrate and sodium thiosulfate electrolytic solutions.
!     Fluid Phase Equilibria, 123:231-243.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, July 1995.
!     Last Modified by MD White, Battelle, 13 Sept 2001.
!     brnlqv.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 A_C(3,2),B_C(6,2),T_C(3,2)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A_C,B_C,T_C
      DATA A_C /9.91828D-4,2.98772D-3,4.22179D-3,
     &  1.00716D-3,3.65659D-3,1.76215D-3 /
      DATA B_C /642.237D+0,-689.63D+0,8273.85D+0,
     &  -51214.8D+0,125889.D+0,-104430.D+0
     &  643.915D+0,121.295D+0,-3997.75D+0,
     &  19616.9D+0,-34793.9D+0,21344.6D+0 /
      DATA T_C /134.711D+0,15.6163D+0,29.475D+0,
     &  134.045D+0,16.4398D+0,136.606D+0 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/BRNLQV'
!
!---  Brine viscosity computed from an empirical relation by Leijnse ---
!
      IF( ISLC(22).EQ.0 ) THEN
        VISLX = VISLX*(1.D+0 + 1.85D+0*XLSX + 4.1D+0*(XLSX**2)
     &    + 4.45D+1*(XLSX**3))
!
!---  Function for sodium nitrate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISLC(22).EQ.1 ) THEN
        A_CX = A_C(1,1)
        DO 100 N = 1,2
          A_CX = A_CX + A_C(N,1)*(XLSX**N)
  100   CONTINUE
        B_CX = B_C(1,1)
        DO 110 N = 1,5
          B_CX = B_CX + B_C(N,1)*(XLSX**N)
  110   CONTINUE
        T_CX = T_C(1,1)
        DO 120 N = 1,2
          T_CX = T_CX + T_C(N,1)*(XLSX**N)
  120   CONTINUE
        VISLX = A_CX*SQRT(TX+TABS)*EXP(B_CX/(TX+TABS-T_CX))
!
!---  Function for sodium thiosulfate by S. Mahiuddin and K. Ismail  ---
!
      ELSEIF( ISLC(22).EQ.2 ) THEN
        A_CX = A_C(1,2)
        DO 200 N = 1,2
          A_CX = A_CX + A_C(N,2)*(XLSX**N)
  200   CONTINUE
        B_CX = B_C(1,2)
        DO 210 N = 1,5
          B_CX = B_CX + B_C(N,2)*(XLSX**N)
  210   CONTINUE
        T_CX = T_C(1,2)
        DO 220 N = 1,2
          T_CX = T_CX + T_C(N,2)*(XLSX**N)
  220   CONTINUE
        VISLX = A_CX*SQRT(TX+TABS)*EXP(B_CX/(TX+TABS-T_CX))
!
!---  Function for sodium chloride Henry's problem  ---
!
      ELSEIF( ISLC(22).EQ.3 ) THEN
        VISLX = VISLX
!
!---  Function for sodium chloride Elder's problem  ---
!
      ELSEIF( ISLC(22).EQ.4 ) THEN
        VISLX = 1.0D-3
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of BRNLQV group  ---
!
      RETURN
      END

