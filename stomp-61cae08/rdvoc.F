!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDVOC
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
!     Read input file for volatile organic compound information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, December 1992.
!     Last Modified by MD White, Battelle, PNL, December 31, 1992.
!     rdvoc.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE NAPL
      USE FILES
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
      CHARACTER*64 ADUM,UNTS
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDVOC'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Oil Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read organic compound name  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Volatile Organic Component Name: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      WRITE (IWR,'(/,2A)') VARB(1:IVR),ADUM
      WRITE (IWR,'(/,A)') 'Thermophysical Constants '
!
!---  Read Molecular weight, gm/mole ---
!---  Read Freezing point, C ---
!---  Read Normal boiling point, C ---
!---  Read Critical temperature, C ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Molecular Weight, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,WTMO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',WTMO
      INDX = 0
      IUNKG = 1
      IUNMOL = -1
      CALL RDUNIT(UNTS,WTMO,INDX)
      RCO = RCU/WTMO
      WRITE (IWR,'(2X,A,1PE11.4)') 'Ideal Gas Constant, J/kg K: ',RCO
      VARB = 'Freezing Point Temperature, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TFPO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',TFPO
      INDX = 0
      IUNK = 1
      CALL RDUNIT(UNTS,TFPO,INDX)
      TFPO = TFPO + TABS
      VARB = 'Normal Boiling Point, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TBO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',TBO
      INDX = 0
      IUNK = 1
      CALL RDUNIT(UNTS,TBO,INDX)
      TBO = TBO + TABS
      VARB = 'Critical Temperature, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,TCRO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',TCRO
      INDX = 0
      IUNK = 1
      CALL RDUNIT(UNTS,TCRO,INDX)
      TCRO = TCRO + TABS
!
!---  Read Critical pressure, Pa  ---
!---  Read Critical volume, cm^3/mole  ---
!---  Read Critical compressibility  ---
!---  Read Pitzner acentric factor  ---
!---  Read Dipole moment, debyes  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Critical Pressure, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,PCRO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':
     & ',PCRO
      INDX = 0
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,PCRO,INDX)
      VARB = 'Critical Molar Volume, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,VCRO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':
     & ',VCRO
      INDX = 0
      IUNM = 3
      IUNMOL = -1
      CALL RDUNIT(UNTS,VCRO,INDX)
!
!---  Convert to units of cm^3/mole  ---
!
      VCRO = 1.D+3*VCRO
      VARB = 'Critical Compressibility: '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,ZCRO)
      WRITE (IWR,'(2X,A,1PE11.4)') VARB(1:IVR),ZCRO
      VARB = 'Pitzner acentric factor: '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,PAFO)
      WRITE (IWR,'(2X,A,1PE11.4)') VARB(1:IVR),PAFO
      VARB = 'Dipole moment, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,DPMO)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',DPMO
      INDX = 0
      CALL RDUNIT(UNTS,DPMO,INDX)
!
!---  Isobaric molar specific heat inputs  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Isobaric Molar Specific Heat Constants: '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,CPOC(1))
      CALL RDDPR(ISTART,ICOMMA,CHDUM,CPOC(2))
      CALL RDDPR(ISTART,ICOMMA,CHDUM,CPOC(3))
      CALL RDDPR(ISTART,ICOMMA,CHDUM,CPOC(4))
      WRITE (IWR,'(/,A,$)') VARB(1:IVR)
      WRITE (IWR,'(A)') 'Cp = a + bT + cT^2 + dT^3'
      WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',CPOC(1)
      WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',CPOC(2)
      WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',CPOC(3)
      WRITE (IWR,'(2X,A,1PE11.4)') 'Constant d: ',CPOC(4)
!
!---  Saturated vapor pressure inputs  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Saturated Vapor Pressure Function: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      WRITE (IWR,'(/,A,$)') 'Vapor Pressure: '
      IF( INDEX(ADUM(1:),'equation 1') .NE. 0 )  THEN
        IVAPO = 1
        VARB = 'Saturated VOC Vapor Pressure Constants (Eqn 1): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(2))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(3))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(4))
        WRITE (IWR,'(A)') 'Pvap = (ax + bx^1.5 + cx^3 + dx^6)/(1-x), whe
     &re x = 1 - T/Tc '
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',SATOC(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',SATOC(2)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',SATOC(3)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant d: ',SATOC(4)
      ELSEIF( INDEX(ADUM(1:),'equation 2') .NE. 0 )  THEN
        IVAPO = 2
        VARB = 'Saturated VOC Vapor Pressure Constants (Eqn 2): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(2))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(3))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(4))
        WRITE (IWR,'(A)') 'ln(Pvap) = a - b/T + cln(T) + dPvap/T^2'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',SATOC(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',SATOC(2)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',SATOC(3)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant d: ',SATOC(4)
      ELSEIF( INDEX(ADUM(1:),'equation 3') .NE. 0 )  THEN
        IVAPO = 3
        VARB = 'Saturated VOC Vapor Pressure Constants (Eqn 3): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(2))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(3))
        WRITE (IWR,'(A)') 'ln(Pvap) = a - b/(T + c)'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',SATOC(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',SATOC(2)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',SATOC(3)
      ELSEIF( INDEX(ADUM(1:),'constant') .NE. 0 ) THEN
        IVAPO = 4
        VARB = 'Saturated VOC Vapor Pressure'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,3A,1PE11.4)')VARB(1:IVR),UNTS(1:NCH),':',SATOC(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,SATOC(1),INDX)
      ELSEIF( INDEX(ADUM(1:),'equation 5') .NE. 0 )  THEN
        IVAPO = 5
        VARB = 'Saturated VOC Vapor Pressure Constants (Eqn 5): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(2))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SATOC(3))
        WRITE (IWR,'(A)') 'ln(Pvap) = a - b/(T + c)'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',SATOC(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',SATOC(2)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',SATOC(3)
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Oil Saturated Vapor Pressure Function: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Liquid density inputs  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      WRITE (IWR,'(/,A,$)') 'Liquid Density: '
      VARB = 'Liquid Density Function: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'hbt technique') .NE. 0 )  THEN
        IRHOO = 1
        WRITE(IWR,'(A)') 'Hankinson-Brobst-Thomson Technique'
        VARB = 'Pure Component Characteristic Volume, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PCCVO)
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':
     & ',PCCVO
        VARB = 'HBT Accentric Factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WHBTO)
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),WHBTO
      ELSEIF( INDEX(ADUM(1:),'modified rackett w/ refer') .NE. 0 )  THEN
        IRHOO = 3
        WRITE(IWR,'(A)') 'Modified Rackett Technique w/ Reference'
        VARB = 'Rackett Compressibility Factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ZRAO)
        IF( ZRAO .LT. SMALL ) ZRAO = 2.9056D-1 - 8.775D-2*PAFO
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),ZRAO
        VARB = 'HBT Accentric Factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WHBTO)
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),WHBTO
        VARB = 'Reference Liquid Density, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RHORO)
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':',RHORO
        INDX = 0
        IUNM = -3
        IUNKG = 1
        CALL RDUNIT(UNTS,RHORO,INDX)
        VARB = 'Reference Temperature, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,TDRO)
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':',TDRO
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,TDRO,INDX)
      ELSEIF( INDEX(ADUM(1:),'modified rackett') .NE. 0 )  THEN
        IRHOO = 2
        WRITE(IWR,'(A)') 'Modified Rackett Technique'
        VARB = 'Rackett Compressibility Factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,ZRAO)
        IF( ZRAO .LT. SMALL ) ZRAO = 2.9056D-1 - 8.775D-2*PAFO
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),ZRAO
        VARB = 'HBT Accentric Factor: '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,WHBTO)
        WRITE(IWR,'(2X,A,1PE11.4)') VARB(1:IVR),WHBTO
      ELSEIF( INDEX(ADUM(1:),'constant') .NE. 0 )  THEN
        IRHOO = 4
        WRITE(IWR,'(A)') 'Constant'
        VARB = 'Reference Liquid Density, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,RHORO)
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),':',RHORO
        INDX = 0
        IUNM = -3
        IUNKG = 1
        CALL RDUNIT(UNTS,RHORO,INDX)
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Oil Liquid Density Function: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Liquid density reference value  ---
!
      TX = TSPRF - TABS
      CALL VOCSP( TX,PSOX )
      CALL VOCLQD(TX,PATM,PSOX,RHORN)
!
!---  Liquid viscosity inputs  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      WRITE (IWR,'(/,A,$)') 'Liquid viscosity: '
      VARB = 'Liquid Viscosity Function: '
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
      IF( INDEX(ADUM(1:),'reference') .NE. 0 ) THEN
        IVISO = 1
        WRITE (IWR,'(A)') 'Reference Values'
        VARB = 'Reference Liquid Viscosity, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,VISCO(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -1
        CALL RDUNIT(UNTS,VISCO(1),INDX)
        VARB = 'Reference Temperature, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,VISCO(2)
        INDX = 0
        IUNK = 1
        CALL RDUNIT(UNTS,VISCO(2),INDX)
      ELSEIF( INDEX(ADUM(1:),'equation 1') .NE. 0 ) THEN
        IVISO = 2
        VARB = 'Liquid Viscosity Constants (Eqn 1): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(2))
        WRITE (IWR,'(A)') 'Semi-Empirical (viscosity = aT^b)'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',VISCO(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',VISCO(2)
      ELSEIF( INDEX(ADUM(1:),'equation 2') .NE. 0 ) THEN
        IVISO = 3
        VARB = 'Liquid Viscosity Constants (Eqn 2): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(2))
        WRITE (IWR,'(A)') 'Semi-Empirical (viscosity = exp(a + b/T))'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',VISCO(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',VISCO(2)
      ELSEIF( INDEX(ADUM(1:),'equation 3') .NE. 0 ) THEN
        IVISO = 4
        VARB = 'Liquid Viscosity Constants (Eqn 3): '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(1))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(2))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(3))
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(4))
        WRITE (IWR,'(A)') 'Semi-Empirical (Viscosity = exp(a + b/T + cT
     &+ dT^2)'
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant a: ',VISCO(1)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant b: ',VISCO(2)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant c: ',VISCO(3)
        WRITE (IWR,'(2X,A,1PE11.4)') 'Constant d: ',VISCO(4)
      ELSEIF( INDEX(ADUM(1:),'constant') .NE. 0 ) THEN
        IVISO = 5
        WRITE (IWR,'(A)') 'Constant'
        VARB = 'Reference Liquid Viscosity, '
        CALL RDDPR(ISTART,ICOMMA,CHDUM,VISCO(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE (IWR,'(2X,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': '
     &,VISCO(1)
        INDX = 0
        IUNM = -1
        IUNKG = 1
        IUNS = -1
        CALL RDUNIT(UNTS,VISCO(1),INDX)
      ELSE
        INDX = 4
        CHMSG = 'Unrecognized Oil Liquid Viscosity Function: '
     &    // ADUM(1:NCH)
        CALL WRMSGS( INDX )
      ENDIF
      TX = TSPRF - TABS
      CALL VOCLQV( TX,VISRN,VISCO,IVISO )
!
!---  Aqueous solubility  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      IDFLT = 1
      HCOW = 1.D+20
      VARB = 'Henry''s Constant for Aqueous Solubility, '
      CALL RDDPR(ISTART,ICOMMA,CHDUM,HCOW)
      CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
      WRITE (IWR,'(/,3A,1PE11.4)') VARB(1:IVR),UNTS(1:NCH),': ',HCOW
      INDX = 0
      IUNM = -1
      IUNKG = 1
      IUNS = -2
      CALL RDUNIT(UNTS,HCOW,INDX)
!
!---  Aqueous-gas surface tension  ---
!
      CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
!
!---  Single exponential function  ---
!
      IF( INDX.EQ.1 ) THEN
        VARB = 'Aqueous-Gas Surface Tension Coefficient: y0'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SFCSF(1))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',SFCSF(1)
        INDX = 0
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,SFCSF(1),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',SFCSF(1),', N/m)'
        VARB = 'Aqueous-Gas Surface Tension Coefficient: a1'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SFCSF(2))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',SFCSF(2)
        INDX = 0
        IUNKG = 1
        IUNS = -2
        CALL RDUNIT(UNTS,SFCSF(2),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',SFCSF(2),', N/m)'
        VARB = 'Aqueous-Gas Surface Tension Coefficient: t1'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,SFCSF(3))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',SFCSF(3)
        INDX = 0
        IUNM = 3
        IUNKG = -1
        CALL RDUNIT(UNTS,SFCSF(3),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',SFCSF(3),', m^3/kg)'
        CALL CHKDPR(ISTART,ICOMMA,CHDUM,INDX)
!
!---    Double exponential function  ---
!
        IF( INDX.EQ.1 ) THEN
          VARB = 'Aqueous-Gas Surface Tension Coefficient: a2'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SFCSF(4))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SFCSF(4)
          INDX = 0
          IUNKG = 1
          IUNS = -2
          CALL RDUNIT(UNTS,SFCSF(4),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SFCSF(4),', N/m)'
          VARB = 'Aqueous-Gas Surface Tension Coefficient: t2'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SFCSF(5))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SFCSF(5)
          INDX = 0
          IUNM = 3
          IUNKG = -1
          CALL RDUNIT(UNTS,SFCSF(5),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SFCSF(5),', m^3/kg)'
        ENDIF
      ENDIF
!
!---  Gas-phase reference viscosity (saturated conditions @ 20 C)  ---
!
      CALL WATSP( TX,PSWX )
      CALL VOCSP( TX,PSOX )
      PSAX = PATM-PSWX-PSOX
      IF( PSAX .LE. ZERO ) THEN
        INDX = 2
        CHMSG = 'Saturated Oil Pressure > Atmospheric Pressure'
        CALL WRMSGS( INDX )
      ENDIF
      PGX = MAX( PSWX+PSOX,PATM )
      XGMW = PSWX/PGX
      XGMO = PSOX/PGX
      XGMA = PSAX/PGX
      CALL WATGSV( TX,VISW )
      CALL AIRGSV( TX,VISA )
      CALL VOCGSV( TX,VISO )
      CALL GASVIS( XGMW,XGMO,XGMA,VISW,VISO,VISA,VISRG )
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDVOC group  ---
!
      RETURN
      END
