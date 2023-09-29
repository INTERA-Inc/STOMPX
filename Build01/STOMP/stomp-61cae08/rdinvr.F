!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDINVR
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
!     Read inverse control card.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 May 2001.
!     Last Modified by MD White, PNNL, 21 May 2001.
!     rdinvr.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PORMED
      USE GRID
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
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*10 FORM1,SIUNIT
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1
      DATA FORM1 / '(2X,2A,I2)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDINVR'
!
!---  Write card information to output file  ---
!
      CARD = 'Inverse Control Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read inverse parameter associations  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Inverse Parameters'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NINVP)
      IF( NINVP.GT.LINVP ) THEN
        INDX = 6
        CHMSG = 'Number of Inverse Parameters > Parameter LINVP'
        CALL WRMSGS( INDX )
      ENDIF
      DO 100 N = 1,NINVP
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Inverse Parameter Name'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_INVP(1,N))
        WRITE(IWR,'(/,3A)') VARB(1:IVR),': ',C_INVP(1,N)(1:NCH)
        VARB = 'Inverse Parameter Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        IF( INDEX(ADUM(1:),'total porosity').NE.0 ) THEN
          I_INVP(3,N) = 1
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),': Total Porosity'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'diffusive porosity').NE.0 ) THEN
          I_INVP(3,N) = 2
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),': Diffusive Porosity'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'porosity').NE.0 ) THEN
          I_INVP(3,N) = 3
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Total and Diffusive Porosity'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'x-dir').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'hydraulic').NE.0 .OR.
     &      INDEX(ADUM(1:),'permeability').NE.0 ) ) THEN
          I_INVP(3,N) = 4
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': X-Dir. Hydraulic Conductivity / Intrinsic Permeability'
          INVUM = 1
          INVUKG = 0
          INVUS = -1
          SIUNIT = 'm/s'
        ELSEIF( INDEX(ADUM(1:),'y-dir').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'hydraulic').NE.0 .OR.
     &      INDEX(ADUM(1:),'permeability').NE.0 ) ) THEN
          I_INVP(3,N) = 5
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Y-Dir. Hydraulic Conductivity / Intrinsic Permeability'
          INVUM = 1
          INVUKG = 0
          INVUS = -1
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'z-dir').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'hydraulic').NE.0 .OR.
     &      INDEX(ADUM(1:),'permeability').NE.0 ) ) THEN
          I_INVP(3,N) = 6
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Z-Dir. Hydraulic Conductivity / Intrinsic Permeability'
          INVUM = 1
          INVUKG = 0
          INVUS = -1
          SIUNIT = 'm/s'
        ELSEIF( INDEX(ADUM(1:),'hydraulic').NE.0 .OR.
     &      INDEX(ADUM(1:),'permeability').NE.0 ) THEN
          I_INVP(3,N) = 7
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Hydraulic Conductivity / Intrinsic Permeability'
          INVUM = 1
          INVUKG = 0
          INVUS = -1
          SIUNIT = 'm/s'
        ELSEIF( INDEX(ADUM(1:),'alpha param').NE.0 ) THEN
          I_INVP(3,N) = 8
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': van Genuchten ''alpha'' Parameter'
          INVUM = -1
          INVUKG = 0
          INVUS = 0
          SIUNIT = '1/m'
        ELSEIF( INDEX(ADUM(1:),'psi param').NE.0 ) THEN
          I_INVP(3,N) = 9
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Brooks-Corey ''psi'' Parameter'
          INVUM = 1
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'm'
        ELSEIF( INDEX(ADUM(1:),'n param').NE.0 ) THEN
          I_INVP(3,N) = 10
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': van Genuchten ''n'' Parameter'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'lambda param').NE.0 ) THEN
          I_INVP(3,N) = 11
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Brooks-Corey ''lambda'' Parameter'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSEIF( INDEX(ADUM(1:),'residual sat').NE.0 ) THEN
          I_INVP(3,N) = 12
          WRITE(IWR,'(2X,2A)') VARB(1:IVR),
     &      ': Residual Aqueous Saturation'
          INVUM = 0
          INVUKG = 0
          INVUS = 0
          SIUNIT = 'null'
        ELSE
          INDX = 2
          CHMSG = 'Unrecognized Inverse Parameter Type: '//ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read rock/soil name or scaling group name ---
!
        VARB = 'Rock/Soil or Scaling Group Name'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
!
!---    Search known rock/soil types for a matching type ---
!
        DO 10 M = 1,NROCK
           IF( ADUM.EQ.ROCK(M) ) THEN
              I_INVP(4,N) = M
              I_INVP(5,N) = 0
              WRITE(IWR,'(2X,2A)') 'Rock/Soil Name: ',ROCK(M)
              GOTO 30
           ENDIF
   10   CONTINUE
!
!---    Search known scaling groups for a matching type ---
!
        IF( ISLC(19).EQ.1 ) THEN
          DO 20 M = 1,NSCALE
             IF( ADUM.EQ.SCALNM(M) ) THEN
                I_INVP(5,N) = M
                I_INVP(4,N) = 0
                WRITE(IWR,'(2X,2A)') 'Scaling Group Name: ',SCALNM(M)
                GOTO 30
             ENDIF
   20     CONTINUE
          INDX = 4
          CHMSG = 'Unrecognized Rock/Soil Type or Scaling Group: '
     &      // ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
        INDX = 4
        CHMSG = 'Unrecognized Rock/Soil Type: ' // ADUM(1:NCH)
        CALL WRMSGS( INDX )
   30   CONTINUE
!
!---    Read a new line  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---    Read minimum parameter value  ---
!
        VARB = 'Minimum Parameter Value'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_INVP(1,N))
        VARB = 'Minimum Parameter Value Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( SIUNIT.NE.'null' ) THEN
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',R_INVP(1,N)
          INDX = 0
          IUNM = INVUM
          IUNKG = INVUKG
          IUNS = INVUS
          CALL RDUNIT(UNTS,R_INVP(1,N),INDX)
          NCH = INDEX( SIUNIT(1:),'  ')-1
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_INVP(1,N),', ',
     &      SIUNIT(1:NCH),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_INVP(1,N)
        ENDIF
!
!---    Read minimum parameter value  ---
!
        VARB = 'Maximum Parameter Value'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_INVP(2,N))
        VARB = 'Maximum Parameter Value Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( SIUNIT.NE.'null' ) THEN
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',R_INVP(2,N)
          INDX = 0
          IUNM = INVUM
          IUNKG = INVUKG
          IUNS = INVUS
          CALL RDUNIT(UNTS,R_INVP(2,N),INDX)
          NCH = INDEX( SIUNIT(1:),'  ')-1
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_INVP(2,N),', ',
     &      SIUNIT(1:NCH),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_INVP(2,N)
        ENDIF
!
!---    Read perturbation parameter value  ---
!
        VARB = 'Perturbation Parameter Value'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_INVP(2,N))
        VARB = 'Perturbation Parameter Value Units'
        CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
        IF( SIUNIT.NE.'null' ) THEN
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',
     &      UNTS(1:NCH),': ',R_INVP(2,N)
          INDX = 0
          IUNM = INVUM
          IUNKG = INVUKG
          IUNS = INVUS
          CALL RDUNIT(UNTS,R_INVP(2,N),INDX)
          NCH = INDEX( SIUNIT(1:),'  ')-1
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_INVP(2,N),', ',
     &      SIUNIT(1:NCH),')'
        ELSE
          WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_INVP(2,N)
        ENDIF
!
!---    Read inverse output format  ---
!
        VARB = 'Inverse Output Format'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_INVP(2,N))
        WRITE(IWR,'(2X,3A)') VARB(1:IVR),': ',C_INVP(2,N)(1:NCH)
!
!---    Read inverse log-transform flag  ---
!
        VARB = 'Inverse Log-Transform Flag'
        CALL RDINT(ISTART,ICOMMA,CHDUM,I_INVP(1,N))
        WRITE(FORM1(9:9),'(I1)') ICOUNT( I_INVP(1,N) )
        WRITE(IWR,FORM1) VARB(1:IVR),': ',I_INVP(1,N)
        IF( I_INVP(1,N).LT.0 .OR. I_INVP(1,N).GT.1 ) THEN
          INDX = 7
          CHMSG = 'Out of Range Inverse Log-Transform Flag: '
          IMSG = I_INVP(1,N)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read inverse estimate flag  ---
!
        VARB = 'Inverse Estimate Flag'
        CALL RDINT(ISTART,ICOMMA,CHDUM,I_INVP(2,N))
        WRITE(FORM1(9:9),'(I1)') ICOUNT( I_INVP(2,N) )
        WRITE(IWR,FORM1) VARB(1:IVR),': ',I_INVP(2,N)
        IF( I_INVP(2,N).LT.0 .OR. I_INVP(2,N).GT.1 ) THEN
          INDX = 7
          CHMSG = 'Out of Range Inverse Estimate Flag: '
          IMSG = I_INVP(2,N)
          CALL WRMSGS( INDX )
        ENDIF
  100 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDINVR group ---
!
      RETURN
      END

