!----------------------Function----------------------------------------!
!
      FUNCTION ESNR( ASLX,ASTX,ASTMAXX,ASNRMX )
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
!     Compute effective residual-NAPL saturation.
!
!     Lenhard, R.J., M. Oostrom, and J.H. Dane.  2003.  "A constitutive
!     model for air-NAPL-water flow in the vadose zone accounting for
!     residual NAPL in strongly water-wet porous media."  Submitted to
!     Journal of Contaminant Hydrology.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 January 2003.
!     Last Modified by Mark White, PNNL, 9 January 2003.
!     rt_napl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ESNR'
!
!---  Compute residual-NAPL for positive effective maximum residual-NAPL
!     saturation  ---
!
      IF( ASNRMX.GT.EPSL ) THEN
!
!---    Linear van Geel et al. formulation  ---
!
        IF( ISNR.EQ.1 ) THEN
          ESNR = ASNRMX*(ASTMAXX-ASTX)
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
!
!---    Modified linear van Geel et al. formulation  ---
!
        ELSEIF( ISNR.EQ.2 ) THEN
          ESNR = ASNRMX*(ASTMAXX-ASLX)
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
!
!---    van Geel et al. Formulation with Land's parameter  ---
!
        ELSEIF( ISNR.EQ.3 ) THEN
          RX = (1.D+0/ASNRMX) - 1.D+0
          ESNR = (ASTMAXX-ASTX)/(1.D+0+(RX*ASTMAXX))
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
!
!---    Modified van Geel et al. Formulation with Land's parameter  ---
!
        ELSEIF( ISNR.EQ.4 ) THEN
          RX = (1.D+0/ASNRMX) - 1.D+0
          ESNR = (ASTMAXX-ASLX)/(1.D+0+(RX*ASTMAXX))
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
!
!---    Lenhard et al. formulation  ---
!
        ELSEIF( ISNR.EQ.5 ) THEN
          ESNR = ASNRMX*SQRT(ASTMAXX-ASLX)*((1.D+0-ASLX)**1.5D+0)
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
!
!---    Modified Lenhard et al. formulation  ---
!
        ELSEIF( ISNR.EQ.6 ) THEN
          ESNR = ASNRMX*SQRT(ASTMAXX-ASTX)*((1.D+0-ASTX)**1.5D+0)
          IF( ESNR.LT.EPSL ) ESNR = 0.D+0
        ENDIF
!
!---  Zero trapped-NAPL for zero effective maximum trapped-NAPL
!     saturation  ---
!
      ELSE
        ESNR = 0.D+0
      ENDIF
!
!---  End of ESNR group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION ESNT( ASLX,ASLMINX,ASNTMX )
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
!     Compute effective trapped-NAPL saturation.
!
!     Kaluarachchi, J.J. and J.C. Parker.  1992.  "Multiphase Flow with
!     a Simplified Model for Oil Entrapment."  Transport in Porous
!     Media, 7:1-14.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 January 2003.
!     Last Modified by Mark White, PNNL, 9 January 2003.
!     rt_napl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ESNT'
!
!---  Compute trapped-NAPL for positive effective maximum trapped-NAPL
!     saturation  ---
!
      IF( (ASLX-ASLMINX).GT.EPSL .AND. ASNTMX.GT.EPSL ) THEN
!
!---    Land's parameter  ---
!
        RX = (1.D+0/ASNTMX) - 1.D+0
!
!---    Formulation of Kaluarachchi and Parker  ---
!
        ESNT = (1.D+0 - ASLMINX)/(1.D+0 + RX*(1.D+0 - ASLMINX)) -
     &    (1.D+0 - ASLX)/(1.D+0 + RX*(1.D+0 - ASLX))
        IF( ESNT.LT.EPSL ) ESNT = 0.D+0
!
!---  Zero trapped-NAPL for zero effective maximum trapped-NAPL
!     saturation  ---
!
      ELSE
        ESNT = 0.D+0
      ENDIF
!
!---  End of ESNT group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION ASTMX( ASLX,ASTX,ESNRX,ASNRMX )
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
!     Compute apparent maximum total-liquid saturation.
!
!     Lenhard, R.J., M. Oostrom, and J.H. Dane.  2003.  "A constitutive
!     model for air-NAPL-water flow in the vadose zone accounting for
!     residual NAPL in strongly water-wet porous media."  Submitted to
!     Journal of Contaminant Hydrology.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 January 2003.
!     Last Modified by Mark White, PNNL, 9 January 2003.
!     rt_napl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ASTMX'
!
!---  Compute apparent maximum total-liquid saturation for positive
!     effective maximum residual-NAPL saturation  ---
!
      IF( ASNRMX.GT.EPSL ) THEN
!
!---    Linear van Geel et al. formulation  ---
!
        IF( ISNR.EQ.1 ) THEN
          ASTMX = ESNRX/ASNRMX + ASTX
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
!
!---    Modified linear van Geel et al. formulation  ---
!
        ELSEIF( ISNR.EQ.2 ) THEN
          ASTMX = ESNRX/ASNRMX + ASLX
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
!
!---    van Geel et al. Formulation with Land's parameter  ---
!
        ELSEIF( ISNR.EQ.3 ) THEN
          RX = (1.D+0/ASNRMX) - 1.D+0
          ASTMX = (ESNRX+ASTX)/(1.D+0-(RX*ESNRX))
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
!
!---    Modified van Geel et al. Formulation with Land's parameter  ---
!
        ELSEIF( ISNR.EQ.4 ) THEN
          RX = (1.D+0/ASNRMX) - 1.D+0
          ASTMX = (ESNRX+ASLX)/(1.D+0-(RX*ESNRX))
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
!
!---    Lenhard et al. formulation  ---
!
        ELSEIF( ISNR.EQ.5 ) THEN
          ASTMX = ASLX + (ESNRX**2)/((ASNRMX**2)*((1.D+0-ASLX)**3))
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
!
!---    Modified Lenhard et al. formulation  ---
!
        ELSEIF( ISNR.EQ.6 ) THEN
          ASTMX = ASTX + (ESNRX**2)/((ASNRMX**2)*((1.D+0-ASTX)**3))
          ASTMX = MIN( ASTMX,1.D+0 )
          IF( ASTMX.LT.EPSL ) ASTMX = 0.D+0
        ENDIF
!
!---  Set apparent maximum total-liquid saturation to the apparent
!     aqueous saturation for zero effective maximum residual-NAPL
!     saturation  ---
!
      ELSE
        ASTMX = ASLX
      ENDIF
!
!---  End of ASTMX group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION ASLMN( ASLX,ESNTX,ASNTMX )
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
!     Compute apparent minimum aqueous saturation..
!
!     Kaluarachchi, J.J. and J.C. Parker.  1992.  "Multiphase Flow with
!     a Simplified Model for Oil Entrapment."  Transport in Porous
!     Media, 7:1-14.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 9 January 2003.
!     Last Modified by Mark White, PNNL, 9 January 2003.
!     rt_napl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ASLMN'
!
!---  Compute apparent minimum aqueous saturation for positive
!     effective maximum trapped-NAPL saturation  ---
!
      IF( ASNTMX.GT.EPSL ) THEN
!
!---    Land's parameter  ---
!
        RX = (1.D+0/ASNTMX) - 1.D+0
!
!---    Formulation of Kaluarachchi and Parker  ---
!
        ESNRX = ESNTX*RX
        ESNRX2 = ESNRX*RX
        ASLMN = (ESNTX + 2.D+0*ESNRX - ESNRX*ASLX + ESNRX2
     &    - ESNRX2*ASLX - ASLX)/(ESNRX + ESNRX2 - ESNRX2*ASLX - 1.D+0)
        ASLMN = MIN( ASLMN,1.D+0 )
        IF( ASLMN.LT.EPSL ) ASLMN = 0.D+0
!
!---  Zero trapped-NAPL for zero effective maximum trapped-NAPL
!     saturation  ---
!
      ELSE
        ASLMN = ASLX
      ENDIF
!
!---  End of ASLMN group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

!----------------------Function----------------------------------------!
!
      FUNCTION PNC( ENPR,PGX,PLX )
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
!     Compute NAPL critical pressure.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 23 January 2003.
!     Last Modified by Mark White, PNNL, 23 January 2003.
!     rt_napl.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/PNC'
!
!---  Compute NAPL critical pressure assuming the NAPL-aqueous and
!     gas-NAPL interfacial tensions sum to the gas-aqueous interfacial
!     tension  ---
!
      PNC = MAX( ((BNL*PLX+BGN*PGX)/(BNL+BGN)),(PLX+ENPR/BNL) )
!
!---  End of PNC group
!
      ISUB_LOG = ISUB_LOG-1
      RETURN
      END

