!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE ADSORP( PCX,INDX,XLX,SLX,PORDX,PORTX,RHOLX,RHOSX,YX )
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
!     Compute the solid-liquid adsorption factor.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on January 30, 1997.
!     adsorp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 PCX(2)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/ADSORP'
!
!---  Linear adsorption isotherm  ---
!
      IF( INDX.EQ.0 ) THEN
        YX = 1.D+0 + PCX(1)/(PORDX*SLX + EPSL)
!
!---  Linear adsorption isotherm w/ Kd  ---
!
      ELSEIF( INDX.EQ.1 ) THEN
        YX = 1.D+0 + PCX(2)*RHOSX*(1.D+0-PORTX)/(PORDX*SLX + EPSL)
!
!---  Freundlich adsorption isotherm  ---
!
      ELSEIF( INDX.EQ.2 ) THEN
        YX = 1.D+0 + PCX(2)*PCX(1)*((XLX*RHOLX)**(PCX(2)-1.D+0))/
     &    (PORDX*SLX + EPSL)
!
!---  Langmuir adsorption isotherm  ---
!
      ELSEIF( INDX.EQ.3 ) THEN
        YX = 1.D+0 + PCX(1)/(PORDX*SLX*((1.D+0+PCX(2)*XLX*RHOLX)**2)
     &    + EPSL)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of ADSORP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE SORB( CLOX,RHOBX,XSOX,IZN )
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
!     Compute the solid-liquid adsorption factor.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle's Pacific Northwest Division, 1997.
!     Last Modified by MD White on January 30, 1997.
!     adsorp.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/SORB'
!
!---  Linear adsorption isotherm w/ Kd  ---
!
      IF( IPCSLD(IZN).EQ.1 ) THEN
        XSOX = PCSLD(1,IZN)*CLOX
!
!---  Linear adsorption isotherm  ---
!
      ELSEIF( IPCSLD(IZN).EQ.2 ) THEN
        XSOX = PCSLD(1,IZN)*CLOX/RHOBX
!
!---  Freundlich adsorption isotherm  ---
!
      ELSEIF( IPCSLD(IZN).EQ.3 ) THEN
        XSOX = PCSLD(1,IZN)*(CLOX**PCSLD(2,IZN))
!
!---  Langmuir adsorption isotherm  ---
!
      ELSEIF( IPCSLD(IZN).EQ.4 ) THEN
        XSOX = (PCSLD(1,IZN)*CLOX)/(1.D+0 + PCSLD(2,IZN)*CLOX)
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of SORB group  ---
!
      RETURN
      END
