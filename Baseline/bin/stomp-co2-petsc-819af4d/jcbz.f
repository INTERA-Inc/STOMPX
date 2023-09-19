!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE JCBZ( ISLV,MU,ML,MK,INDX )
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
!     Zero the Jacobian matrix elements and solution vector.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1996.
!     Last Modified by MD White on September 5, 1996.




!     jcbz.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_CL
      USE SOLTN
      USE LEAK_WELL
      USE JACOB
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/JCBZ'
      ICNV = 3
      NSPILL = LSPILL
!
!---  Number of unknowns equals number of active nodes for
!     transport equations ---
!
      IF( INDX.EQ.1 ) THEN
        NUKX = NFBN-NRFN-NXP+NWN_LW
!
!---  Geomechanics - number of unknowns equals number of active
!     finite element nodes x 3
!
      ELSEIF( INDX.EQ.2 ) THEN
        NUKX = 3*NFEN_GM
!
!---  Fully coupled matrix-fracture-borehole flow - number of unknowns 
!     equals (number of active nodes + fracture triangles 
!     + borehole nodes) x number of coupled equations
!     + number of coupled wells  ---
!
      ELSEIF( INDX.EQ.5 ) THEN
        NUKX = ISLV*(NFLD-NXP+NT_FRC-NXP_FRC+NBN_BH) + N_CW
!
!---  Fully coupled matrix-fracture-borehole transport - number of  
!     unknowns equals (number of active nodes + fracture triangles 
!     + borehole nodes)
!
      ELSEIF( INDX.EQ.6 ) THEN
        NUKX = NFLD-NXP+NT_FRC-NXP_FRC+NBN_BH
!
!---  Number of unknowns equals number of active nodes x number of
!     coupled equations + number of coupled wells  ---
!
      ELSEIF( N_CW.GT.0 ) THEN
        NUKX = ISLV*(NFBN-NRFN-NXP+NWN_LW) + N_CW
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells +
!     number of surface nodes)
!
      ELSEIF( LSPILL.GT.0 ) THEN
        NUKX = ISLV*(NFBN-NRFN-NXP+NWLN+NWN_LW+IJFLD)
!
!---  Number of unknowns equals number of number of coupled equations
!     x (number of active nodes + number of wells)
!
      ELSE
        NUKX = ISLV*(NFBN-NRFN-NXP+NWLN+NWN_LW)
      ENDIF
!
!---  Dual-porosity option  ---
!
      IF( ISLC(11).EQ.1 ) NUKX = NUKX + ISLV*(NFBN-NRFN-NXP+NWN_LW)
!
!---  Banded solver  ---
!
      IF( ILES.EQ.1 ) THEN
        NCOL = NUKX
        NROW = 2*MU + ML + 1
        DO 110 MCOL = 1,NCOL
          DO 100 MROW = 1,NROW
            ALU(MROW,MCOL) = 0.D+0
  100     CONTINUE
          BLU(MCOL) = 0.D+0
          ILU(MCOL) = MCOL
  110   CONTINUE
!
!---  SPLib Solver  ---
!
      ELSEIF( ILES.EQ.3 ) THEN
        NROW = NUKX
        DO MROW = 1,MK
          DLU(MROW) = 0.D+0
        ENDDO
        DO MROW = 1,NROW
          BLU(MROW) = 0.D+0
        ENDDO

!
!---  PETSc Solver  ---
!
      ELSEIF( ILES.EQ.5 ) THEN
        NROW = NUKX
        DO MROW = 1,MK
          DLU(MROW) = 0.D+0
        ENDDO
        DO MROW = 1,NROW
          BLU(MROW) = 0.D+0
        ENDDO

!
!---  Unknown Solver  ---
!
      ELSE






        INDX = 3
        CHMSG = 'Unknown Linear Equation Solver'
        CALL WRMSGS( INDX )
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of JCBZ group  ---
!
      RETURN
      END

