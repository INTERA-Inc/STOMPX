!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDFRAC
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
!     Read input file for rock/soil fracture property information.
!
!     FRAC_P(1,IZN) x-direction fracture spacing, m
!     FRAC_P(2,IZN) x-direction fracture aperture, m
!     FRAC_P(3,IZN) y-direction fracture spacing, m
!     FRAC_P(4,IZN) y-direction fracture aperture, m
!     FRAC_P(5,IZN) z-direction fracture spacing, m
!     FRAC_P(6,IZN) z-direction fracture aperture, m
!     FRAC_P(7,IZN) fracture volume fraction
!     FRAC_P(8,IZN) fracture porosity
!
!     ISLC(11) = 1 Equivalent Continuum Model
!     FRAC_P(9,IZN) fracture bulk x-direction eff. permeability, m^2
!     FRAC_P(10,IZN) fracture bulk y-direction eff. permeability, m^2
!     FRAC_P(11,IZN) fracture z-direction eff. permeability, m^2
!
!     ISLC(11) = 2 Dual Porosity Model
!     FRAC_P(9,IZN) fracture x-direction permeability, m^2
!     FRAC_P(10,IZN) fracture y-direction permeability, m^2
!     FRAC_P(11,IZN) fracture z-direction permeability, m^2
!
!     FRAC_P(12,IZN) effective fracture aperture
!     FRAC_P(13,IZN) fraction of fracture surface contact area
!     FRAC_P(14,IZN) fracture van Genuchten s-P alpha parameter, 1/m
!     FRAC_P(15,IZN) fracture van Genuchten s-P n parameter
!
!     FRAC_P(16,IZN) fracture van Genuchten k-s kappa parameter
!     FRAC_P(17,IZN) fracture van Genuchten k-s lambda parameter
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
      USE FILES
      USE DUAL_POR
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
      SUB_LOG(ISUB_LOG) = '/RDFRAC'
!
!---  Write card information to ouput file  ---
!
      CARD = 'Rock/Soil Fracture Properties Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE (IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Loop over the rock/soil fracture properties information lines  ---
!
      N = 0
      IJK = 0
      ISGRP = 0
   10 CONTINUE
        IF( N.GE.NROCK .OR. IJK.GT.0 ) GOTO 500
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Rock/Soil Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    IJK, KIJ, or JKI indexing ---
!
        IF( INDEX(ADUM(1:),'indexing').NE.0 ) THEN
          IF( INDEX(ROCK(1)(1:),'indexing').EQ.0 ) THEN
            INDX = 4
            CHMSG = 'Indexing Option Not Declared ' // 
     &        'in Rock/Soil Zonation Card'
            CALL WRMSGS( INDX )
          ENDIF
          IF( INDEX(ADUM,'ijk').NE.0 ) THEN
            IJK = 1
          ELSEIF( INDEX(ADUM,'jki').NE.0 ) THEN
            IJK = 2
          ELSEIF( INDEX(ADUM,'kij').NE.0 ) THEN
            IJK = 3
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Indexing Option' // ADUM(1:NCH)
            CALL WRMSGS( INDX )
          ENDIF
          GOTO 220
        ENDIF
!
!---    Search known rock types for a matching type ---
!
        DO 100 M = 1,NROCK
          IF( ADUM .EQ. ROCK(M)) THEN
            IROCK = M
            GOTO 200
          ENDIF
  100   CONTINUE
!
!---    Search known scaling groups for a matching type ---
!
        IF( ISLC(19).EQ.1 ) THEN
          DO 110 M = 1,NSCALE
             IF( ADUM.EQ.SCALNM(M) ) THEN
                ISGRP = M
                IROCK = 1
                GOTO 200
             ENDIF
  110     CONTINUE
          INDX = 2
          CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
          GOTO 10
        ENDIF
        INDX = 2
        CHMSG = 'Unrecognized Rock/Soil Type: ' // ADUM(1:NCH)
        CALL WRMSGS( INDX )
        GOTO 10
  200   CONTINUE
!
!---    Loop over rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. ISGRP.NE.0 ) THEN
          DO 202 M = IROCK,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              IROCK = M
              GOTO 204
            ENDIF
  202     CONTINUE
        ENDIF
  204   CONTINUE
!
!---    Write rock/soil name  ---
!
        WRITE (IWR,'(/,2A)') 'Rock/Soil Name: ',ROCK(IROCK)
        N = N + 1
  220   CONTINUE
!
!---    Read x-direction fracture spacing ---
!
        VARB = 'X-Direction Fracture Spacing'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 1
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(1,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(1,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(1,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(1,IROCK),', m)'
        ENDIF
!
!---    Read x-direction fracture aperture ---
!
        VARB = 'X-Direction Fracture Aperture'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 2
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(2,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(2,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(2,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(2,IROCK),', m)'
        ENDIF
!
!---    Read y-direction fracture spacing ---
!
        VARB = 'Y-Direction Fracture Spacing'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 3
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(3,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(3,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(3,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(3,IROCK),', m)'
        ENDIF
!
!---    Read y-direction fracture aperture ---
!
        VARB = 'Y-Direction Fracture Aperture'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 4
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(4,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(4,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(4,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(4,IROCK),', m)'
        ENDIF
!
!---    Read z-direction fracture spacing ---
!
        VARB = 'Z-Direction Fracture Spacing'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 5
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(5,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(5,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(5,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(5,IROCK),', m)'
        ENDIF
!
!---    Read z-direction fracture aperture ---
!
        VARB = 'Z-Direction Fracture Aperture'
        UNTS = 'm'
        IUNM = 1
        IF( IJK.GT.0 ) THEN
          INDX = 6
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(6,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',FRAC_P(6,IROCK)
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(6,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(6,IROCK),', m)'
        ENDIF
!
!---    Fracture volume fraction ---
!
        IF( IJK.GT.0 ) THEN
          DO IROCK = 1,NFLD
            FRAC_P(7,IROCK) = 1.D+0 - ((FRAC_P(1,IROCK)-FRAC_P(2,IROCK))
     &        *(FRAC_P(3,IROCK)-FRAC_P(4,IROCK))
     &        *(FRAC_P(5,IROCK)-FRAC_P(6,IROCK)))/
     &        (FRAC_P(1,IROCK)*FRAC_P(3,IROCK)*FRAC_P(5,IROCK))
          ENDDO
        ELSE
          FRAC_P(7,IROCK) = 1.D+0 - ((FRAC_P(1,IROCK)-FRAC_P(2,IROCK))
     &      *(FRAC_P(3,IROCK)-FRAC_P(4,IROCK))
     &      *(FRAC_P(5,IROCK)-FRAC_P(6,IROCK)))/
     &      (FRAC_P(1,IROCK)*FRAC_P(3,IROCK)*FRAC_P(5,IROCK))
        ENDIF
!
!---    Read fracture porosity ---
!
        VARB = 'Fracture Porosity'
        UNTS = 'null'
        INDX = 8
        LNDX = 15
        IF( IJK.GT.0 ) THEN
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(8,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',FRAC_P(8,IROCK)
        ENDIF
!
!---    Read bulk x-direction fracture permeability for 
!       zero-permeability matrix  ---
!
        VARB = 'Bulk X-Direction Fracture Intrinsic Permeability'
        IF( IJK.GT.0 ) THEN
          ISTX = ISTART
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            UNTS = 'm^2'
            IUNM = 2
          ELSE
            UNTS = 'm/s'
            IUNM = 1
            IUNS = -1
          ENDIF
          ISTART = ISTX
          INDX = 9
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(9,IROCK))
          IF( ICOMMA.EQ.ISTART .AND. IFLD.GT.1 ) THEN
            INDX = 4
            CHMSG = 'Zero Bulk X-Direction Fracture Intrinsic Perm.'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(9,IROCK)
            IUNM = 2
          ELSE
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(9,IROCK)
            IUNM = 1
            IUNS = -1
          ENDIF
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(9,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(9,IROCK),', m^2)'
        ENDIF
!
!---    Read bulk y-direction fracture permeability for 
!       zero-permeability matrix  ---
!
        VARB = 'Bulk Y-Direction Fracture Intrinsic Permeability'
        IF( IJK.GT.0 ) THEN
          ISTX = ISTART
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            UNTS = 'm^2'
            IUNM = 2
          ELSE
            UNTS = 'm/s'
            IUNM = 1
            IUNS = -1
          ENDIF
          ISTART = ISTX
          INDX = 10
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(10,IROCK))
          IF( ICOMMA.EQ.ISTART .AND. IFLD.GT.1 ) THEN
            INDX = 4
            CHMSG = 'Zero Bulk Y-Direction Fracture Intrinsic Perm.'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(10,IROCK)
            IUNM = 2
          ELSE
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(10,IROCK)
            IUNM = 1
            IUNS = -1
          ENDIF
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(10,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(10,IROCK),', m^2)'
        ENDIF
!
!---    Read bulk z-direction fracture permeability for 
!       zero-permeability matrix  ---
!
        VARB = 'Bulk Z-Direction Fracture Intrinsic Permeability'
        IF( IJK.GT.0 ) THEN
          ISTX = ISTART
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            UNTS = 'm^2'
            IUNM = 2
          ELSE
            UNTS = 'm/s'
            IUNM = 1
            IUNS = -1
          ENDIF
          ISTART = ISTX
          INDX = 11
          LNDX = 15
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,FRAC_P,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,FRAC_P(11,IROCK))
          IF( ICOMMA.EQ.ISTART .AND. IFLD.GT.1 ) THEN
            INDX = 4
            CHMSG = 'Zero Bulk Z-Direction Fracture Intrinsic Perm.'
            CALL WRMSGS( INDX )
          ENDIF
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          IHC = INDEX( UNTS(1:),'hc' )
          IF( IHC.EQ.0 ) THEN
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(11,IROCK)
            IUNM = 2
          ELSE
            WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &        UNTS(1:NCH),': ',FRAC_P(11,IROCK)
            IUNM = 1
            IUNS = -1
          ENDIF
          INDX = 0
          CALL RDUNIT(UNTS,FRAC_P(11,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',FRAC_P(11,IROCK),', m^2)'
        ENDIF
!
!---    Read effective fracture aperture, m ---
!
        VARB = 'Effective Fracture Aperture'
        IF( IJK.GT.0 ) THEN
          INDX = 15
          LNDX = LSCHR
          UNTS = 'm'
          IUNM = 1
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(15,IROCK))
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',SCHR(15,IROCK)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,SCHR(15,IROCK),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',SCHR(15,IROCK),', m)'
        ENDIF
        VARB = 'Fraction of Fracture Surface Contact Area'
        IF( IJK.GT.0 ) THEN
          INDX = 16
          LNDX = LSCHR
          UNTS = 'null'
          CALL RDIJKD( ISTART,IJK,CHDUM,UNTS,SCHR,INDX,LNDX )
        ELSE
          CALL RDDPR(ISTART,ICOMMA,CHDUM,SCHR(16,IROCK))
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',SCHR(16,IROCK)
        ENDIF
        IF( IJK.GT.0 ) THEN
          DO N = 1,NFLD
            CALL FRAC_KSP(N,IERR)
            IF( IERR.EQ.1 ) THEN
              INDX = 4
              CHMSG = 'Fracture ksP Property Failure'
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
        ELSE
          CALL FRAC_KSP(IROCK,IERR)
          IF( IERR.EQ.1 ) THEN
            INDX = 4
            CHMSG = 'Fracture ksP Property Failure'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Fracture s-P van Genuchten ''alpha'' Parameter'
          IVR = INDEX( VARB,'  ')-1 
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),', m :',
     &      FRAC_P(14,IROCK)
          VARB = 'Fracture s-P van Genuchten ''n'' Parameter'
          IVR = INDEX( VARB,'  ')-1 
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),' :',
     &      FRAC_P(15,IROCK)
          VARB = 'Fracture k-s van Genuchten ''kappa'' Parameter'
          IVR = INDEX( VARB,'  ')-1 
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),' :',
     &      FRAC_P(16,IROCK)
          VARB = 'Fracture k-s van Genuchten ''lambda'' Parameter'
          IVR = INDEX( VARB,'  ')-1 
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),' :',
     &      FRAC_P(17,IROCK)
        ENDIF
!
!---    Loop over remaining rock/soils within scaling group  ---
!
        IF( ISLC(19).EQ.1 .AND. IROCK.LT.NROCK ) THEN
          DO 490 M = IROCK+1,NROCK
            IF( ISCALE(M).EQ.ISGRP ) THEN
              N = N+1
              DO 480 L = 1,15
                FRAC_P(L,M) = FRAC_P(L,IROCK)
  480         CONTINUE
            ENDIF
  490     CONTINUE
        ENDIF
!
!---    Read next rock/soil type or scaling group  ---
!
        IF( N.LT.NROCK ) WRITE(IWR,'(/)')
        GOTO 10
 500  CONTINUE
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDFRAC group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FRAC_KSP( IZN,IERR )
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
!     Development of van Genucthen retention parameters from
!     effective fracture aperture and fraction of fracture surface
!     area contact.
!
!     Wang, J. S. Y., and T. N. Narasimhan, “ Hydrologic Mechanisms 
!     Governing Fluid Flow in a Partially Saturated, Fractured Porous 
!     Medium.” Water Resources Research, 12(12):1861-1874.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 November 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE DUAL_POR
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 AJ(2,2),BJ(2),GX(2)
      REAL*8 SX(256),RKX(256),SIGX(256)
      REAL*8 AX(2)
      REAL*8 COVARX(2,2),WORKX(2,2)
      INTEGER LISTA(2)
      INTEGER IJ(2)
      EXTERNAL FKRVG
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/FRAC_KSP'
!
!---  Guess beta and bc  ---
!
      IERR = 0
      BETAX = 3.84D+4
      BCX = 1.55D-5
      BBX = SCHR(15,IZN)
      WX = SCHR(16,IZN)
!
!---  Two-variable Newton-Raphson loop to find beta and bc from the
!     effective fracture aperture b, and the fraction of fracture
!     surface contact area w  ---
!
   10 CONTINUE
      BJ(1) = -(WX-1.D+0+EXP(-BETAX*BCX)*(1.D+0+BETAX*BCX))
      AJ(1,1) = -(BCX**2)*EXP(-BETAX*BCX)*BETAX
      AJ(1,2) = -(BETAX**2)*EXP(-BETAX*BCX)*BCX
      BJ(2) = -((BETAX**3)-(6.D+0*(4.D+0 + BETAX*BCX)/(BBX**3))*
     &  EXP(-BETAX*BCX))
      AJ(2,1) = 3.D+0*((BETAX**2)*(BBX**3) + 
     &  6.D+0*BCX*EXP(-BETAX*BCX) + 
     &  2.D+0*(BCX**2)*EXP(-BETAX*BCX)*BETAX)/(BBX**3)
      AJ(2,2) = 6.D+0*BETAX*EXP(-BETAX*BCX)*(3.D+0+BETAX*BCX)/(BBX**3)
!
!---  Solve linear system  ---
!
      JP = 2
      KP = 2
      CALL LUDCMP( AJ,JP,KP,IJ,DJ )
      CALL LUBKSB( AJ,JP,KP,IJ,BJ )
      DBETAX = BJ(1)
      DBCX = BJ(2)
      BETAX = BETAX + DBETAX
      BCX = BCX + DBCX
      IF( ABS(DBETAX)/ABS(BETAX).GT.1.D-6 .OR. 
     &  ABS(DBCX)/ABS(BCX).GT.1.D-6 ) GOTO 10
!
!---  Assign surface tension of liquid water to 0.0781 kg/s^2  ---
!
      GAMMAX = 0.0718
      PCX = (1.D+0/3.D+0)*GAMMAX*(BETAX*BCX - 4.D+0 + 
     &  SQRT((BETAX**2)*(BCX**2) + 4.D+0*BETAX*BCX + 16.D+0))/(BCX)
!
!---  Guess van Genucthen n  ---
!
      CNX = 2.D+0
!
!---  Newton-Raphson loop to find van Genuchten n
!     by assuming slr = 0, and equating the capillary pressure, and 
!     capillary pressure times the slope of the saturation versus
!     capillary pressure function for both the Wang and Narasimhan
!     and van Genuchen functions  ---
!
   20 CONTINUE
      DCNY = 1.D-6
      DO M = 1,2
        CNY = CNX
        IF( M.EQ.2 ) CNY = CNX + DCNY
        ALPHAX = ((CNY-1.D+0)/((PCX**CNY)*CNY))**(1.D+0/CNY)
        VARX = (ALPHAX*PCX)**CNY
        GX(M) = ((1.D+0/(1.D+0 + VARX)**((CNY-1.D+0)/CNY)))*
     &    ((CNY-1.D+0)/(1.D+0 + VARX))*VARX -
     &    4.D+0*(BETAX**3)*((GAMMAX**2)/(PCX**3))*
     &    EXP(-2.D+0*BETAX*GAMMAX/PCX)*
     &    ((2.D+0*GAMMAX + BCX*PCX)/(2.D+0 + BETAX*BCX))
      ENDDO
      DCNX = -GX(1)/((GX(2)-GX(1))/DCNY)
      CNX = CNX + DCNX
      CNX = MAX( 1.D+0,CNX )
      IF( ABS(DCNX)/ABS(CNX).GT.1.D-6 ) GOTO 20
      ALPHAX = ((CNX-1.D+0)/((PCX**CNX)*CNX))**(1.D+0/CNX)
!
!---  Return van Genuchten alpha (1/m) and n parameters  ---
!
      FRAC_P(14,IZN) = ALPHAX*RHORL*GRAV
      FRAC_P(15,IZN) = CNX
!
!---  Generate 256 data points for the Wang and Narasimham relative
!     permeability model for fractures  ---
!
      NDATA = 256
      PCLX = 1.D+1
      PCHX = 1.D+7
      DO M = 1,NDATA
        PCX = EXP( REAL(M-1)*(LOG(PCHX)-LOG(PCLX))/REAL(NDATA-1) 
     &    + LOG(PCLX))
        BSX = 2.D+0*GAMMAX/PCX
        SX(M) = (1.D+0/(2.D+0 + BETAX*BCX))*((2.D+0 - EXP(-BETAX*BSX)*
     &    (2.D+0 + 2.D+0*BETAX*BSX + (BETAX*BSX)**2)) + BETAX*BCX*
     &    (1.D+0 - EXP(-BETAX*BSX)*(1.D+0 + BETAX*BSX)))
        RKX(M) = (1.D+0/(6.D+0*(4.D+0 + BETAX*BCX)))*
     &    ((24.D+0 - EXP(-BETAX*BSX)*(24.D+0 + 24.D+0*BETAX*BSX + 12.D+0
     &    *(BETAX*BSX)**2 + 4.D+0*(BETAX*BSX)**3 + (BETAX*BSX)**4))
     &    + BETAX*BCX*(6.D+0 - EXP(-BETAX*BSX)*(6.D+0 + 6.D+0*BETAX*BSX
     &    + 3.D+0*(BETAX*BSX)**2 + (BETAX*BSX)**3)))
        SIGX(M) = 1.D-2
      ENDDO
      MA = 2
      MFIT = 2
      LISTA(1) = 1
      LISTA(2) = 2
      NCA = 2
      AX(1) = -5.D-1
      AX(2) = 7.D-1
      ALAMDAX = -1.D+0
      CALL MRQMIN( SX,RKX,SIGX,NDATA,AX,MA,LISTA,MFIT,COVARX,WORKX,
     &  NCA,CHISQX,FKRVG,ALAMDAX )
      DO M = 1,32
        OCHISQX = CHISQX
        CALL MRQMIN( SX,RKX,SIGX,NDATA,AX,MA,LISTA,MFIT,COVARX,WORKX,
     &    NCA,CHISQX,FKRVG,ALAMDAX )
        IF( ABS(OCHISQX-CHISQX).LT.1.D-9 ) EXIT
      ENDDO
      IF( M.GT.32 ) IERR = 1
!
!---  Return van Genuchten kappa and lambda parameters  ---
!
      FRAC_P(16,IZN) = AX(1)
      FRAC_P(17,IZN) = AX(2)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of FRAC_KSP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE FKRVG( X,AX,Y0X,DYDAX,NA )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     van Genuchten relative permeability model
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 22 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(NA) :: AX,DYDAX
!
!----------------------Executable Lines--------------------------------!
!
      DAX = 1.D-7
      VAR1X = 1.D+0-(X**(1.D+0/AX(2)))
      VAR2X = VAR1X**AX(2)
      Y0X = (X**AX(1))*(1.D+0-VAR2X)**2
      A1X = AX(1) + DAX
      VAR1X = 1.D+0-(X**(1.D+0/AX(2)))
      VAR2X = VAR1X**AX(2)
      Y1X = (X**A1X)*(1.D+0-VAR2X)**2
      A2X = AX(2) + DAX
      VAR1X = 1.D+0-(X**(1.D+0/A2X))
      VAR2X = VAR1X**A2X
      Y2X = (X**A1X)*(1.D+0-VAR2X)**2
      DYDAX(1) = (Y1X-Y0X)/DAX
      DYDAX(2) = (Y2X-Y0X)/DAX
!
!---  End of FKRVG group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MRQMIN( X,Y,SIG,NDATA,A,MA,LISTA,MFIT,COVAR,ALPHA,
     &  NCA,CHISQ,FUNCS,ALAMDA )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Levenberg-Marquardt method
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 22 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      PARAMETER (MMAX=20)
      REAL*8, DIMENSION(NDATA) :: X,Y,SIG
      REAL*8, DIMENSION(MA) :: A
      REAL*8, DIMENSION(NCA,NCA) :: COVAR,ALPHA
      REAL*8, DIMENSION(MMAX) :: ATRY,BETAX,DA
      INTEGER, DIMENSION(NCA) :: IJ
      INTEGER, DIMENSION(MA) :: LISTA
      EXTERNAL FUNCS
!
!----------------------Executable Lines--------------------------------!
!
      IF( ALAMDA.LT.0.D+0 ) THEN
        KK = MFIT + 1
        DO J = 1,MA
          IHIT = 0
          DO K = 1,MFIT
            IF( LISTA(K).EQ.J ) IHIT = IHIT + 1
          ENDDO
          IF( IHIT.EQ.0 ) THEN
            LISTA(KK) = J
            KK = KK+1
          ELSEIF( IHIT.GT.1 ) THEN
            INDX = 26
            CHMSG = 'Improper permutation in LISTA'
            CALL WRMSGS( INDX )
          ENDIF
        ENDDO
        IF( KK.NE.(MA+1) ) THEN
          INDX = 26
          CHMSG = 'Improper permutation in LISTA'
          CALL WRMSGS( INDX )
        ENDIF
        ALAMDA = 1.D-3
        CALL MRQCOF( X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETAX,NCA,
     &    CHISQ,FUNCS )
        OCHISQ = CHISQ
        DO J = 1,MA
          ATRY(J) = A(J)
        ENDDO
      ENDIF
      DO J = 1,MFIT
        DO K = 1,MFIT
          COVAR(J,K) = ALPHA(J,K)
        ENDDO
        COVAR(J,J) = ALPHA(J,J)*(1.D+0+ALAMDA)
        DA(J) = BETAX(J)
      ENDDO
      CALL LUDCMP( COVAR,MFIT,NCA,IJ,DJ )
      CALL LUBKSB( COVAR,MFIT,NCA,IJ,DA )
      IF( ALAMDA.LT.1.D-14 ) THEN
        CALL COVSRT( COVAR,NCA,MA,LISTA,MFIT )
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      DO J = 1,MFIT
        ATRY(LISTA(J)) = A(LISTA(J))+DA(J)
      ENDDO
      OCHISQ = CHISQ
      CALL MRQCOF( X,Y,SIG,NDATA,ATRY,MA,LISTA,MFIT,COVAR,DA,NCA,
     &    CHISQ,FUNCS )
      IF( CHISQ.LT.OCHISQ ) THEN
        ALAMDA = 1.D-1*ALAMDA
        OCHISQ = CHISQ
        DO J = 1,MFIT
          DO K = 1,MFIT
            ALPHA(J,K) = COVAR(J,K)
          ENDDO
          BETAX(J) = DA(J)
          A(LISTA(J)) = ATRY(LISTA(J))
        ENDDO
      ELSE
        ALAMDA = 1.D+0*ALAMDA
        CHISQ = OCHISQ
      ENDIF
!
!---  End of MRQMIN group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE MRQCOF( X,Y,SIG,NDATA,A,MA,LISTA,MFIT,ALPHA,BETAX,
     &  NALP,CHISQ,FUNCS )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Used by MRQMIN to evaluate the linearized fitting matrix ALPHA
!     and vector BETAX
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 22 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      PARAMETER (MMAX=20)
      REAL*8, DIMENSION(NDATA) :: X,Y,SIG
      REAL*8, DIMENSION(NALP,NALP) :: ALPHA
      REAL*8, DIMENSION(MA) :: A,BETAX
      REAL*8, DIMENSION(MMAX) :: DYDA
      INTEGER, DIMENSION(MFIT) :: LISTA
      EXTERNAL FUNCS
!
!----------------------Executable Lines--------------------------------!
!
      DO J = 1,MFIT
        DO K = 1,J
          ALPHA(J,K) = 0.D+0
        ENDDO
        BETAX(J) = 0.D+0
      ENDDO
      CHISQ = 0.D+0
      DO I = 1,NDATA
        CALL FUNCS( X(I),A,YMOD,DYDA,MA )
        SIG2I = 1.D+0/(SIG(I)*SIG(I))
        DY = Y(I) - YMOD
        DO J = 1,MFIT
          WT = DYDA(LISTA(J))*SIG2I
          DO K = 1,J
            ALPHA(J,K) = ALPHA(J,K) + WT*DYDA(LISTA(K))
          ENDDO
          BETAX(J) = BETAX(J) + DY*WT
        ENDDO
        CHISQ = CHISQ + DY*DY*SIG2I
      ENDDO
      DO J = 2,MFIT
        DO K = 1,J-1
          ALPHA(K,J) = ALPHA(J,K)
        ENDDO
      ENDDO
!
!---  End of MRQCOF group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE COVSRT( COVAR,NCVM,MA,LISTA,MFIT )
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
!     EOR (Enhanced Oil Recovery and Carbon Utilization) Mode
!
!     Given the covariance matrix COVAR of a fit for MFIT of MA total
!     parameters, and their ordering LISTA(I), repack the covariance
!     matrix to the true order of the parameters.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark D White, PNNL, 22 December 2015
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(NCVM,NCVM) :: COVAR
      INTEGER, DIMENSION(MFIT) :: LISTA
!
!----------------------Executable Lines--------------------------------!
!
      DO J = 1,MA-1
        DO I = J+1,MA
          COVAR(I,J) = 0.D+0
        ENDDO
      ENDDO
      DO I = 1,MFIT-1
        DO J = I+1,MFIT
          IF( LISTA(J).GT.LISTA(I) ) THEN
            COVAR(LISTA(J),LISTA(I)) = COVAR(I,J)
          ELSE
            COVAR(LISTA(I),LISTA(J)) = COVAR(I,J)
          ENDIF
        ENDDO
      ENDDO
      SWAP = COVAR(1,1)
      DO J = 1,MA
        COVAR(1,J) = COVAR(J,J)
        COVAR(J,J) = 0.D+0
      ENDDO
      COVAR(LISTA(1),LISTA(1)) = SWAP
      DO J = 2,MFIT
        COVAR(LISTA(J),LISTA(J)) = COVAR(1,J)
      ENDDO
      DO J = 2,MA
        DO I = 1,J-1
          COVAR(I,J) = COVAR(J,I)
        ENDDO
      ENDDO
!
!---  End of COVSRT group
!
      RETURN
      END
