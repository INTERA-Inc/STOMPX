!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATLQK( TX,PX,PSWX,THK )
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
!     Calculate the subcooled or saturated thermal conductivity, as a
!     function of temperature per the steam table equations
!     as given by the 1967 International Formulation Committee:
!     Formulation for Industrial Use.
!
!     Thermodynamic and Transport Properties of Steam.
!     1967. ASME Steam Tables.
!     The American Society of Mechanical Engineers.
!     United Engineering Center, 345 East 47th Street, New York, N.Y.
!
!     The temperature is limited in this subroutine to the following
!     values:  0.01 C < T > 364.0 !
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, January, 1992.
!     Last Modified by MD White, Battelle, PNL, January 6, 1992.
!     watlqk.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 CFA(5),CFB(4),CFC(4)
!
!----------------------Data Statements---------------------------------!
!
      SAVE CFA,CFB,CFC
      DATA CFA /-922.47D+0,2839.5D+0,-1800.7D+0,525.77D+0,-73.440D+0/
      DATA CFB /-0.94730D+0,2.5186D+0,-2.0012D+0,0.51536D+0/
      DATA CFC /1.6563D-3,-3.8929D-3,2.9323D-3,-7.1693D-4/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATLQK'
      TR = MIN( (TX+TABS)/TABS,TCRW/TABS )
      PR = (PX - PSWX)*1.D-5
      THK = CFA(1)+CFA(2)*TR+CFA(3)*TR**2+CFA(4)*TR**3+CFA(5)*TR**4
      THK = THK + PR*(CFB(1)+CFB(2)*TR+CFB(3)*TR**2+CFB(4)*TR**3)
      THK = THK + PR*PR*(CFC(1)+CFC(2)*TR+CFC(3)*TR**2+CFC(4)*TR**3)
      THK = THK*1.D-3
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATLQK group  ---
!
      RETURN
      END

