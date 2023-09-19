!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDPLANT
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
!     Read input file for plant property information.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Gene Freeman, PNNL, 9 April 2002.
!     Last Modified by MD White, PNNL, 9 April 2002.
!     rdplant.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!-------------------List of input parameters and units ----------------!
!
!     Maximum extent in Z (Zm), m
!     Maximum extent in X (Xm), m
!     Maximum extent in Y (Ym), m
!     Fit parameter in Z (z*), m
!     Fit parameter in X (x*), m
!     Fit parameter in Y (y*), m
!     Leaf area index (LAI), unitless
!     Plant canopy height (PCH), m
!     Water stress head 1 (h1), m
!     Water stress head 2 (h2), m
!     Water stress head 3 (h3), m
!     Water stress head 4 (h4), m
!     Water uptake reduced by 50% (h50), m
!     Crop coefficient - stage1, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage2, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage3, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Crop coefficient - stage4, unitless
!     Crop start time jan.1=1, dec.31=365, day
!     Evapotranspiration (ET), m/day
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
      USE FILES
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 UNTS,ADUM
      CHARACTER*512 CHDUM
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDPLANT'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Plant Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of plant varietals  ---
!
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      ISTART = 1
      VARB = 'Number of Plants'
      CALL RDINT( ISTART,ICOMMA,CHDUM,NPLANT )
      IF( NPLANT.GT.LPLANT ) THEN
        INDX = 5
        CHMSG = 'Number of Plant Varietals > Parameter LPLANT'
        CALL WRMSGS( INDX )
      ENDIF
      IF( NPLANT.GT.(LBCV-2) ) THEN
        INDX = 5
        CHMSG = 'Number of Plant Varietals > Parameter (LBCV-2)'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Check for multiple plant temperature option and
!     rainfall interception-condensation shedding option  ---
!
      CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
      IF( INDX.EQ.1 ) THEN
        VARB = 'Plant Temperature and Rainfall Interception Options'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!        IF( INDEX(ADUM(1:),'multiple').NE.0 ) THEN
!          ISLC(24) = ISLC(24)+2
!          WRITE(IWR,'(A)') '  Multiple Plant Temperature Option'
!        ELSE
          ISLC(24) = ISLC(24)+1
          WRITE(IWR,'(A)') '  Single Plant Temperature Option'
!        ENDIF
        IF( INDEX(ADUM(1:),'rainfall').NE.0 .OR.
     &    INDEX(ADUM(1:),'interception').NE.0 ) THEN
          ISLC(26) = 1
          WRITE(IWR,'(A)') '  Rainfall Interception'
        ELSE
          WRITE(IWR,'(A)') '  No Rainfall Interception'
        ENDIF
      ELSE
        ISLC(24) = ISLC(24)+1
        WRITE(IWR,'(A)') '  Single Plant Temperature Option'
        WRITE(IWR,'(A)') '  No Rainfall Interception'
      ENDIF
!
!---  Loop over the plants information lines  ---
!
      DO 500 IP = 1,NPLANT
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plant Name: '
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,PLANT(IP))
        DO 100 M = 1,IP-1
          IF( PLANT(M).EQ.PLANT(IP) ) THEN
            INDX = 4
            CHMSG = 'Duplicate Plant Name: ' // PLANT(IP)
            CALL WRMSGS( INDX )
          ENDIF
  100   CONTINUE
        WRITE (IWR,'(/,2A)') ' Plant Name: ',PLANT(IP)
!
!---    Check for root stress or stomatal resistance options  ---
!
        CALL CHKCHR( ISTART,ICOMMA,CHDUM,INDX )
        IF( INDX.EQ.1 ) THEN
          VARB = 'Plant Options'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          IF( INDEX(ADUM(1:),'vrugt').NE.0 ) THEN
            IRSM_P(IP) = 1
            WRITE(IWR,'(A)') '  Vrugt Root Stress Model'
          ELSEIF( INDEX(ADUM(1:),'jarvis').NE.0 ) THEN
            IRSM_P(IP) = 2
            WRITE(IWR,'(A)') '  Jarvis Root Stress Model'
          ELSEIF( INDEX(ADUM(1:),'stress').NE.0 ) THEN
            IRSM_P(IP) = 1
            WRITE(IWR,'(A)') '  Vrugt Root Stress Model'
          ELSE
            WRITE(IWR,'(A)') '  No Root Stress Model'
          ENDIF
          IF( INDEX(ADUM(1:),'hicks').NE.0 ) THEN
            ISRM_P(IP) = 1
            WRITE(IWR,'(A)') '  Hicks Stomatal Resistance Model'
          ELSE
            WRITE(IWR,'(A)') '  No Stomatal Resistance Model'
          ENDIF
        ELSE
          WRITE(IWR,'(A)') '  No Root Stress Model'
          WRITE(IWR,'(A)') '  No Stomatal Resistance Model'
        ENDIF
!
!---    Read root (Z) depth characteristics  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Max. Root Depth'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(1,IP) = 1.D-3
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(1,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(1,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(1,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(1,IP),', m)'
!
!---    Minimum maximum plant root depth of 1 mm  ---
!
        PARMS_P(1,IP) = MAX( 1.D-3,PARMS_P(1,IP) )
!
!---    Read root z* characteristics  ---
!
        VARB = 'Null Root Depth'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(4,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(4,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(4,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(4,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(4,IP),', m)'
!
!---    Read root pz characteristics  ---
!
        VARB = 'Root Depth fit parameter'
        IDFLT = 1
        PARMS_P(7,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(7,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(7,IP)
!
!---    Read plant short-wave albedo  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Plant Solar Albedo'
!
!---    Check for temporal input to be read in conjunction with crop
!       coefficient times  ---
!
        CALL RDCHR(ISTART,ICOMMA,NCHA,CHDUM,ADUM)
        IF( INDEX(ADUM,'temporal').NE.0 ) THEN
          IALB_P(IP) = 1
!
!---      Read plant solar albedo (initial stage start)  ---
!
          VARB = 'Plant Solar Albedo: Initial Stage Start'
          IDFLT = 1
          PARMS_P(5,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(5,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(5,IP)
!
!---      Read plant solar albedo (crop development start)  ---
!
          VARB = 'Plant Solar Albedo: Crop Development Start'
          IDFLT = 1
          PARMS_P(6,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(6,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(6,IP)
!
!---      Read plant solar albedo (mid-season start)  ---
!
          VARB = 'Plant Solar Albedo: Mid-Season Start'
          IDFLT = 1
          PARMS_P(8,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(8,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(8,IP)
!
!---      Read plant solar albedo (late-season start)  ---
!
          VARB = 'Plant Solar Albedo: Late-Season Start'
          IDFLT = 1
          PARMS_P(9,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(9,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(9,IP)
!
!---      Read plant solar albedo (late-season stop)  ---
!
          VARB = 'Plant Solar Albedo: Late-Season Stop'
          IDFLT = 1
          PARMS_P(10,IP) = PARMS_P(5,IP)
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(10,IP))
          WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &        ': ',PARMS_P(10,IP)
        ELSE
          ISTART = 1
          IDFLT = 1
          PARMS_P(10,IP) = 1.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(10,IP))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),
     &        ': ',PARMS_P(10,IP)
        ENDIF
!
!---    Read plant height characteristics  ---
!
        VARB = 'Plant Canopy Height'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(11,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(11,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(11,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(11,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(11,IP),', m)'
!
!---    Minimum maximum plant canopy height of 1 mm  ---
!
        PARMS_P(11,IP) = MAX( 1.D-3,PARMS_P(11,IP) )
!
!---    Maximum condensate depth  ---
!
        IF( ISLC(26).EQ.1 ) THEN
          VARB = 'Maximum Condensate Depth'
          UNTS = 'm'
          IUNM = 1
          IDFLT = 1
          PARMS_P(16,IP) = 0.2D-3
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(16,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',PARMS_P(16,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(16,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(16,IP),', m)'
        ENDIF
!
!---  Vrugt root stress model
!
      IF( IRSM_P(IP).EQ.1 ) THEN
!
!---  Read first stress point head  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Root Stress-Point 1 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(12,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(12,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(12,IP),', m)'
!
!---  Read second stress point head  ---
!
        VARB = 'Root Stress-Point 2 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(13,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(13,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(13,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(13,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(13,IP),', m)'
!
!---  Read third stress point head  ---
!
        VARB = 'Root Stress-Point 3 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(14,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(14,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(14,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(14,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(14,IP),', m)'
!
!---  Read fourth stress point head  ---
!
        VARB = 'Root Stress-Point 4 (Capillary Head)'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(15,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(15,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(15,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(15,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(15,IP),', m)'
!
!---  Jarvis root stress model
!
      ELSEIF( IRSM_P(IP).EQ.2 ) THEN
!
!---  Read first Jarvis stress point  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Wilting-Point Water Content'
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(12,IP)
!
!---  Read second Jarvis stress point  ---
!
        VARB = 'Normalized Soil Water Content: Critical Point 1'
        IDFLT = 1
        PARMS_P(13,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(13,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(13,IP)
!
!---  Read third Jarvis stress point  ---
!
        VARB = 'Normalized Soil Water Content: Critical Point 2'
        IDFLT = 1
        PARMS_P(14,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(14,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(14,IP)
!
!---  Read fourth Jarvis stress point  ---
!
        VARB = 'Saturated Water Content'
        IDFLT = 1
        PARMS_P(15,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(15,IP))
        WRITE(IWR,'(2A,1PE11.4)') VARB(1:IVR),': ',PARMS_P(15,IP)
      ELSE
!
!---  Read head when uptake reduce 50%  ---
!
       CALL RDINPL( CHDUM )
       CALL LCASE( CHDUM )
       ISTART = 1
       VARB = 'root uptake reduced 50%'
        UNTS = 'm'
        IUNM = 1
        IDFLT = 1
        PARMS_P(12,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(12,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(12,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(12,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(12,IP),', m)'
      ENDIF
!
!---  Read crop coefficient (initial stage start)  ---
!
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        ISTART = 1
        VARB = 'Crop Coefficient: Initial Stage Start'
        IDFLT = 1
        PARMS_P(17,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(17,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(17,IP)
!
!---  Read crop coefficient day of year (initial stage start)  ---
!
        VARB = 'Initial Stage Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(18,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(18,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(18,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(18,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(18,IP),', s)'
!
!---  Read crop coefficient (crop development start)  ---
!
        VARB = 'Crop Coefficient: Crop Development Start'
        IDFLT = 1
        PARMS_P(19,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(19,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(19,IP)
!
!---  Read crop coefficient day of year (crop development start)  ---
!
        VARB = 'Crop Development Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(20,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(20,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(20,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(20,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(20,IP),', s)'
!
!---  Read crop coefficient (mid-season start)  ---
!
        VARB = 'Crop Coefficient: Mid-Season Start'
        IDFLT = 1
        PARMS_P(21,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(21,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(21,IP)
!
!---  Read crop coefficient day of year (mid-season start)  ---
!
        VARB = 'Mid-Season Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(22,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(22,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(22,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(22,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(22,IP),', s)'
!
!---  Read crop coefficient (late-season start)  ---
!
        VARB = 'Crop Coefficient: Late-Season Start'
        IDFLT = 1
        PARMS_P(23,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(23,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(23,IP)
!
!---  Read crop coefficient day of year (late-season start)  ---
!
        VARB = 'Late-Season Start Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(24,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(24,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(24,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(24,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(24,IP),', s)'
!
!---  Read crop coefficient (late-season stop)  ---
!
        VARB = 'Crop Coefficient: Late-Season Stop'
        IDFLT = 1
        PARMS_P(25,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(25,IP))
        WRITE(IWR,'(2X,2A,1PE11.4,$)') VARB(1:IVR),
     &      ': ',PARMS_P(25,IP)
!
!---  Read crop coefficient day of year (late-season stop)  ---
!
        VARB = 'Late-Season Stop Day'
        UNTS = 'day'
        IUNS = 1
        IDFLT = 1
        PARMS_P(26,IP) = 1.D+0
        CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(26,IP))
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',PARMS_P(26,IP)
        INDX = 0
        CALL RDUNIT(UNTS,PARMS_P(26,IP),INDX)
        WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(26,IP),', s)'
!
!---    Hicks stomatal resistance model  ---
!
        IF( ISRM_P(IP).EQ.1 ) THEN
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
          ISTART = 1
!
!---      Read minimum stomatal resistance, s/m  ---
!
          VARB = 'Stomatal Resistance: Minimum Stomatal Resistance'
          UNTS = 's/m'
          IUNM = -1
          IUNS = 1
          IDFLT = 1
          PARMS_P(2,IP) = 100.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(2,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(2,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(2,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(2,IP),', s/m)'
!
!---      Read light response coefficient, W/m^2  ---
!
          VARB = 'Stomatal Resistance: Light Response Coefficient'
          UNTS = 'w/m^2'
          IUNKG = 1
          IUNS = -3
          IDFLT = 1
          PARMS_P(3,IP) = 20.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(3,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(3,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(3,IP),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',PARMS_P(3,IP),', W/m^2)'
!
!---      Read minimum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Minimum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(27,IP) = 5.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(27,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(27,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(27,IP),INDX)
          PARMS_P(27,IP) = PARMS_P(27,IP)+TABS
!
!---      Read maximum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Maximum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(28,IP) = 45.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(28,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(28,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(28,IP),INDX)
          PARMS_P(28,IP) = PARMS_P(28,IP)+TABS
!
!---      Read optimum temperature for stomatal opening, C  ---
!
          VARB = 'Stomatal Resistance: Optimum Temperature'
          UNTS = 'c'
          IUNK = 1
          IDFLT = 1
          PARMS_P(29,IP) = 25.D+0
          CALL RDDPR(ISTART,ICOMMA,CHDUM,PARMS_P(29,IP))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &        ': ',PARMS_P(29,IP)
          INDX = 0
          CALL RDUNIT(UNTS,PARMS_P(29,IP),INDX)
          PARMS_P(29,IP) = PARMS_P(29,IP)+TABS
        ENDIF
!
!---  Read next plant type  ---
!
        IF( IP.LT.NPLANT ) WRITE(IWR,'(/)')
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDPLNT group ---
!
      RETURN
      END
