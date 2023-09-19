!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WATICK( TX,THKIX )
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
!     Calculate the thermal conductivity of ice as a function of
!     temperature using polynomial fit of data from Dickerson (1969).
!
!     Dickerson, R. W., Jr., Thermal properties of food, in The Freezing
!     Preservation of Foods, 4th ed., Vol. 2, D. K. Tressler, W. B. Van
!     Arnsdel, and M. J. Copley (Editors), AVI Publishing Co., Westport,
!     Conn.
!
!     The temperature is limited in this subroutine to the following
!     values:  -100.0 C < T <= 0.0 !
!
!----------------------Authors-----------------------------------------!
!
!     Written by WE Nichols, Battelle, PNL, September, 1994.
!     Last Modified by WE Nichols, Battelle, PNL, September, 1994.
!     watick.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 A(3)
!
!----------------------Data Statements---------------------------------!
!
      SAVE A
      DATA A / 7.39519D+0,-2.86936D-2,3.54452D-5 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WATICK'
      TK = TX +TABS
      THKIX = A(1) +A(2)*TK +A(3)*TK**2
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WATICK group  ---
!
      RETURN
      END

