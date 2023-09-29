!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDTRMX
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
!     Read input file for ternary mixture information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, July 31, 2000.
!     Last Modified by Mark White, PNNL, July 31, 2000.
!     rdtrmx.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDTRMX'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Ternary Mixture Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Ternary Phase Diagram Constants  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Binodal Curve Type Option'
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'mole') .NE. 0 ) THEN
        ITERDC = 1
        WRITE(IWR,'(2A)') VARB(1:IVR),': Mole Fraction'
      ELSEIF( INDEX(ADUM(1:),'mass') .NE. 0 ) THEN
        ITERDC = 2
        WRITE(IWR,'(2A)') VARB(1:IVR),': Mass Fraction'
      ELSEIF( INDEX(ADUM(1:),'volu') .NE. 0 ) THEN
        ITERDC = 3
        WRITE(IWR,'(2A)') VARB(1:IVR),': Volume Fraction'
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Binodal Curve Type Option: '//ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
      VARB = 'Binodal Curve Constant: Ah'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(1))
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(1)
      VARB = 'Binodal Curve Constant: Bh'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(2))
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(2)
      VARB = 'Binodal Curve Constant: Eh'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(3))
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(3)
      VARB = 'Binodal Curve Constant: Fh'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(4))
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(4)
      VARB = 'Solubility Limit of Water in Alcohol-Free Liquid Oil'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(8))
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(8)
      TERDC(9) = TERDC(8)
      VARB = 'Solubility Limit of Oil in Alcohol-Free Liquid Water'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TERDC(10))
      TERDC(11) = TERDC(10)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),': ',TERDC(10)
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDTRMX group  ---
!
      RETURN
      END
