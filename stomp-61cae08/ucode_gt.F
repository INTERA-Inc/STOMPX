!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDOBDA
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
!     Read observed data card.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 May 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE TRNSPT
      USE SOLTN
      USE OUTPU
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
      CHARACTER*64 ADUM,BDUM,CDUM,FDUM,SOLNM,UNTS,FLNM
      CHARACTER*512 CHDUM
      CHARACTER*10 FORM1
      CHARACTER*4 FORM2,FORM4
      CHARACTER*6 FORM3
      LOGICAL FLG_CHK
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1,FORM2,FORM3,FORM4
      DATA FORM1 / '(2X,2A,I2)' /
      DATA FORM2 / '(I2)' /
      DATA FORM3 / '(I6,$)' /
      DATA FORM4 / '(I6)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDOBDA'
!
!---  Write card information to output file  ---
!
      CARD = 'Observed-Data Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Read number of observed-data types  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
      VARB = 'Number of Observed-Data Types'
      CALL RDINT(ISTART,ICOMMA,CHDUM,NOBDT)
      IF( NOBDT.GT.LOBDT ) THEN
        INDX = 6
        CHMSG = 'Number of Observed-Data Types > Parameter LOBDT'
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Loop over the number of observed-data types  ---
!
      DO 300 NT = 1,NOBDT
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
        VARB = 'Observed-Data Type'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
!
!---  Field-observation variable  ---
!
        IF( INDEX(ADUM(1:),'field').NE.0 ) THEN
          VARB = 'Field-Observation Variable'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,BDUM )
          IVRX = IVR
          IF( INDEX( BDUM(1:),'solute' ).NE.0 ) THEN
            VARB = 'Field-Observation Variable: Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
            DO 10 NSL = 1,NSOLU
              IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 20
   10       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Field-Observation Variable ' //
     &        'Solute Name: ' // SOLNM
            CALL WRMSGS( INDX )
   20       CONTINUE
          ENDIF
          I_OBDT(1,NT) = 1
!
!---      Convert field-observation variable to an index  ---
!
          VARB = 'Field-Observation Variable'
          IVR = IVRX
          IF( INDEX(BDUM(1:),'atmospheric temp').NE.0 ) THEN
            I_OBDT(2,NT) = 206
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Temperature'
          ELSEIF( INDEX(BDUM(1:),'atmospheric rel').NE.0 ) THEN
            I_OBDT(2,NT) = 207
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Relative Humidity'
          ELSEIF( INDEX(BDUM(1:),'atmospheric solar').NE.0 ) THEN
            I_OBDT(2,NT) = 208
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Solar Radiation'
          ELSEIF( INDEX(BDUM(1:),'atmospheric wind').NE.0 ) THEN
            I_OBDT(2,NT) = 209
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Wind Speed'
          ELSEIF( INDEX(BDUM(1:),'surface temp').NE.0 ) THEN
            I_OBDT(2,NT) = 213
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Ground Surface Temperature'
          ELSEIF( INDEX(BDUM(1:),'surface vapor press').NE.0 ) THEN
            I_OBDT(2,NT) = 214
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Ground Surface Water-Vapor Pressure'
          ELSEIF( INDEX(BDUM(1:),'actual evapor').NE.0 ) THEN
            I_OBDT(2,NT) = 215
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Actual Evaporation Rate'
          ELSEIF( INDEX(BDUM(1:),'potential evapor').NE.0 ) THEN
            I_OBDT(2,NT) = 216
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Potential Evaporation Rate'
          ELSEIF( INDEX(BDUM(1:),'actual trans').NE.0 ) THEN
            I_OBDT(2,NT) = 217
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Actual Transpiration Rate'
          ELSEIF( INDEX(BDUM(1:),'potential trans').NE.0 ) THEN
            I_OBDT(2,NT) = 218
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Potential Transpiration Rate'
          ELSEIF( INDEX(BDUM(1:),'aqueous pres').NE.0 ) THEN
            I_OBDT(2,NT) = 1
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Pressure'
          ELSEIF( INDEX(BDUM(1:),'gas pres').NE.0 ) THEN
            I_OBDT(2,NT) = 2
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Pressure'
          ELSEIF( INDEX(BDUM(1:),'aqueous sat').NE.0 ) THEN
            I_OBDT(2,NT) = 11
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Saturation'
          ELSEIF( INDEX(BDUM(1:),'gas sat').NE.0 ) THEN
            I_OBDT(2,NT) = 12
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Saturation'
          ELSEIF( INDEX(BDUM(1:),'aqueous moisture cont').NE.0 ) THEN
            I_OBDT(2,NT) = 15
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Moisture Content'
          ELSEIF( INDEX(BDUM(1:),'water gas mass').NE.0 ) THEN
            I_OBDT(2,NT) = 21
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Gas Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'air gas mass').NE.0 ) THEN
            I_OBDT(2,NT) = 22
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Gas Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'water aqueous mass').NE.0 ) THEN
            I_OBDT(2,NT) = 23
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Aqueous Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'air aqueous mass').NE.0 ) THEN
            I_OBDT(2,NT) = 24
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Aqueous Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'aqueous hydraulic head').NE.0 ) THEN
            I_OBDT(2,NT) = 27
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Hydrualic Head'
          ELSEIF( INDEX(BDUM(1:),'gas hydraulic head').NE.0 ) THEN
            I_OBDT(2,NT) = 28
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Hydrualic Head'
          ELSEIF( INDEX(BDUM(1:),'aqueous relative perm').NE.0 ) THEN
            I_OBDT(2,NT) = 31
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Relative Permeability'
          ELSEIF( INDEX(BDUM(1:),'gas relative perm').NE.0 ) THEN
            I_OBDT(2,NT) = 32
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Relative Permeability'
          ELSEIF( INDEX(BDUM(1:),'aqueous density').NE.0 ) THEN
            I_OBDT(2,NT) = 34
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Density'
          ELSEIF( INDEX(BDUM(1:),'gas density').NE.0 ) THEN
            I_OBDT(2,NT) = 35
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Density'
          ELSEIF( INDEX(BDUM(1:),'x thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 44
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'y thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 45
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'z thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 46
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'x aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 51
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'y aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 52
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'z aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 53
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'x gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 54
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'y gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 55
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'z gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 56
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'x heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 60
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'y heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 61
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'z heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 62
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'matric potential').NE.0 ) THEN
            I_OBDT(2,NT) = 63
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Matric Potential'
          ELSEIF( INDEX(BDUM(1:),'water gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 73
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Gas Concentration'
          ELSEIF( INDEX(BDUM(1:),'air gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 74
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Gas Concentration'
          ELSEIF( INDEX(BDUM(1:),'water aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 76
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Aqueous Concentration'
          ELSEIF( INDEX(BDUM(1:),'air aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 77
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Aqueous Concentration'
          ELSEIF( INDEX(BDUM(1:),'solute volumetric conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 1
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Volumetric Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 2
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Aqueous Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 3
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Gas Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute aqueous mol').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 5
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Aqueous Mole Fracton'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute gas mol').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 6
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Gas Mole Fracton'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'x solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 8
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'y solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 9
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'z solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 10
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Field-Observation Variable: '
     &        // BDUM
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Read field-observation output units  ---
!
          VARB = 'Field-Observation Output Units'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_OBDT(NT) )
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),
     &      ': ',C_OBDT(NT)(1:NCH)
!
!---      Read field-observation x-dir. coordinate and units  ---
!
          VARB = 'Field-Observation X-Dir. Coordinate'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(1,NT))
          VARB = 'Field-Observation X-Dir. Coordinate Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',R_OBDT(1,NT)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,R_OBDT(1,NT),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_OBDT(1,NT),', m)'
          IF( R_OBDT(1,NT).LT.XE(1,1) .OR.
     &      R_OBDT(1,NT).GT.XE(8,NFLD) ) THEN
            INDX = 4
            CHMSG = 'Field-Observation X-Dir. Coordinate ' //
     &        'Outside of Computational Domain'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Read field-observation y-dir. coordinate and units  ---
!
          VARB = 'Field-Observation Y-Dir. Coordinate'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(2,NT))
          VARB = 'Field-Observation Y-Dir. Coordinate Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',R_OBDT(2,NT)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,R_OBDT(2,NT),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_OBDT(2,NT),', m)'
          IF( R_OBDT(2,NT).LT.YE(1,1) .OR.
     &      R_OBDT(2,NT).GT.YE(8,NFLD) ) THEN
            INDX = 4
            CHMSG = 'Field-Observation Y-Dir. Coordinate ' //
     &        'Outside of Computational Domain'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Read field-observation z-dir. coordinate and units  ---
!
          VARB = 'Field-Observation Z-Dir. Coordinate'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(3,NT))
          VARB = 'Field-Observation Z-Dir. Coordinate Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &      ': ',R_OBDT(3,NT)
          INDX = 0
          IUNM = 1
          CALL RDUNIT(UNTS,R_OBDT(3,NT),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_OBDT(3,NT),', m)'
          IF( R_OBDT(3,NT).LT.ZE(1,1) .OR.
     &      R_OBDT(3,NT).GT.ZE(8,NFLD) ) THEN
            INDX = 4
            CHMSG = 'Field-Observation Z-Dir. Coordinate ' //
     &        'Outside of Computational Domain'
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Determine field-observation bounding I,J,K indices  ---
!
          ILO = 1
          IHI = IFLD
          JLO = 1
          JHI = JFLD
          KLO = 1
          KHI = KFLD
          IF( R_OBDT(3,NT).GT.ZP(ND(1,1,KFLD)) ) THEN
            KLO = KFLD
          ELSEIF( R_OBDT(3,NT).LT.ZP(ND(1,1,1)) ) THEN
            KHI = 1
          ENDIF
          IF( R_OBDT(2,NT).GT.YP(ND(1,JFLD,1)) ) THEN
            JLO = JFLD
          ELSEIF( R_OBDT(2,NT).LT.YP(ND(1,1,1)) ) THEN
            JHI = 1
          ENDIF
          IF( R_OBDT(1,NT).GT.XP(ND(IFLD,1,1)) ) THEN
            ILO = IFLD
          ELSEIF( R_OBDT(1,NT).LT.XP(ND(1,1,1)) ) THEN
            IHI = 1
          ENDIF
   30     CONTINUE
          IF( KHI-KLO.GT.1 ) THEN
            K = (KHI+KLO)/2
            N = ND(1,1,K)
            IF( ZP(N).GT.R_OBDT(3,NT) ) THEN
              KHI = K
            ELSE
              KLO = K
            ENDIF
            GOTO 30
          ENDIF
   40     CONTINUE
          IF( JHI-JLO.GT.1 ) THEN
            J = (JHI+JLO)/2
            N = ND(1,J,1)
            IF( YP(N).GT.R_OBDT(2,NT) ) THEN
              JHI = J
            ELSE
              JLO = J
            ENDIF
            GOTO 40
          ENDIF
   50     CONTINUE
          IF( IHI-ILO.GT.1 ) THEN
            I = (IHI+ILO)/2
            N = ND(I,1,1)
            IF( XP(N).GT.R_OBDT(1,NT) ) THEN
              IHI = I
            ELSE
              ILO = I
            ENDIF
            GOTO 50
          ENDIF
          DO K = KLO,KHI
          DO J = JLO,JHI
          DO I = ILO,IHI
            N = ND(I,J,K)
            IF( IXP(N).EQ.0 ) THEN
              INDX = 7
              IMSG = N
              CHMSG = 'Field-Observation ' //
     &          'Spacial Interpolation: Inactive Node:'
              CALL WRMSGS( INDX )
            ENDIF
          ENDDO
          ENDDO
          ENDDO
          I_OBDT(4,NT) = ILO
          I_OBDT(5,NT) = IHI
          I_OBDT(6,NT) = JLO
          I_OBDT(7,NT) = JHI
          I_OBDT(8,NT) = KLO
          I_OBDT(9,NT) = KHI
!
!---  Reference-observation variable  ---
!
        ELSEIF( INDEX(ADUM(1:),'reference').NE.0 ) THEN
          I_OBDT(1,NT) = 2
          VARB = 'Reference-Observation Variable'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          IVRX = IVR
          IF( INDEX( ADUM(1:),'solute' ).NE.0 ) THEN
            VARB = 'Reference-Observation Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
            DO 70 NSL = 1,NSOLU
              IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 72
   70       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Reference-Observation Solute Name: '
     &        // SOLNM
            CALL WRMSGS( INDX )
   72       CONTINUE
          ENDIF
!
!---      Convert reference-observation variable
!         to an index  ---
!
          VARB = 'Reference-Observation Variable'
          IVR = IVRX
          IF( INDEX(BDUM(1:),'atmospheric temp').NE.0 ) THEN
            I_OBDT(2,NT) = 206
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Temperature'
          ELSEIF( INDEX(BDUM(1:),'atmospheric rel').NE.0 ) THEN
            I_OBDT(2,NT) = 207
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Relative Humidity'
          ELSEIF( INDEX(BDUM(1:),'atmospheric solar').NE.0 ) THEN
            I_OBDT(2,NT) = 208
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Solar Radiation'
          ELSEIF( INDEX(BDUM(1:),'atmospheric wind').NE.0 ) THEN
            I_OBDT(2,NT) = 209
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Atmospheric Wind Speed'
          ELSEIF( INDEX(BDUM(1:),'surface temp').NE.0 ) THEN
            I_OBDT(2,NT) = 213
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Ground Surface Temperature'
          ELSEIF( INDEX(BDUM(1:),'surface vapor press').NE.0 ) THEN
            I_OBDT(2,NT) = 214
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Ground Surface Water-Vapor Pressure'
          ELSEIF( INDEX(BDUM(1:),'actual evapor').NE.0 ) THEN
            I_OBDT(2,NT) = 215
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Actual Evaporation Rate'
          ELSEIF( INDEX(BDUM(1:),'potential evapor').NE.0 ) THEN
            I_OBDT(2,NT) = 216
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Potential Evaporation Rate'
          ELSEIF( INDEX(BDUM(1:),'actual trans').NE.0 ) THEN
            I_OBDT(2,NT) = 217
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Actual Transpiration Rate'
          ELSEIF( INDEX(BDUM(1:),'potential trans').NE.0 ) THEN
            I_OBDT(2,NT) = 218
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Potential Transpiration Rate'
          ELSEIF( INDEX(BDUM(1:),'aqueous pres').NE.0 ) THEN
            I_OBDT(2,NT) = 1
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Pressure'
          ELSEIF( INDEX(BDUM(1:),'gas pres').NE.0 ) THEN
            I_OBDT(2,NT) = 2
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Pressure'
          ELSEIF( INDEX(BDUM(1:),'aqueous sat').NE.0 ) THEN
            I_OBDT(2,NT) = 11
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Saturation'
          ELSEIF( INDEX(BDUM(1:),'gas sat').NE.0 ) THEN
            I_OBDT(2,NT) = 12
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Saturation'
          ELSEIF( INDEX(BDUM(1:),'aqueous moisture cont').NE.0 ) THEN
            I_OBDT(2,NT) = 15
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Moisture Content'
          ELSEIF( INDEX(BDUM(1:),'water gas mass').NE.0 ) THEN
            I_OBDT(2,NT) = 21
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Gas Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'air gas mass').NE.0 ) THEN
            I_OBDT(2,NT) = 22
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Gas Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'water aqueous mass').NE.0 ) THEN
            I_OBDT(2,NT) = 23
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Aqueous Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'air aqueous mass').NE.0 ) THEN
            I_OBDT(2,NT) = 24
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Aqueous Mass Fraction'
          ELSEIF( INDEX(BDUM(1:),'aqueous hydraulic head').NE.0 ) THEN
            I_OBDT(2,NT) = 27
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Hydrualic Head'
          ELSEIF( INDEX(BDUM(1:),'gas hydraulic head').NE.0 ) THEN
            I_OBDT(2,NT) = 28
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Hydrualic Head'
          ELSEIF( INDEX(BDUM(1:),'aqueous relative perm').NE.0 ) THEN
            I_OBDT(2,NT) = 31
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Relative Permeability'
          ELSEIF( INDEX(BDUM(1:),'gas relative perm').NE.0 ) THEN
            I_OBDT(2,NT) = 32
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Relative Permeability'
          ELSEIF( INDEX(BDUM(1:),'aqueous density').NE.0 ) THEN
            I_OBDT(2,NT) = 34
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Aqueous Density'
          ELSEIF( INDEX(BDUM(1:),'gas density').NE.0 ) THEN
            I_OBDT(2,NT) = 35
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Gas Density'
          ELSEIF( INDEX(BDUM(1:),'x thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 44
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'y thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 45
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'z thermal cond').NE.0 ) THEN
            I_OBDT(2,NT) = 46
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Effective Thermal Conductivity'
          ELSEIF( INDEX(BDUM(1:),'x aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 51
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'y aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 52
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'z aqueous vol').NE.0 ) THEN
            I_OBDT(2,NT) = 53
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Aqueous Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'x gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 54
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'y gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 55
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'z gas vol').NE.0 ) THEN
            I_OBDT(2,NT) = 56
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Gas Volumetric Flux'
          ELSEIF( INDEX(BDUM(1:),'x heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 60
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'y heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 61
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'z heat flux').NE.0 ) THEN
            I_OBDT(2,NT) = 62
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z-Dir. Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'matric potential').NE.0 ) THEN
            I_OBDT(2,NT) = 63
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Matric Potential'
          ELSEIF( INDEX(BDUM(1:),'water gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 73
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Gas Concentration'
          ELSEIF( INDEX(BDUM(1:),'air gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 74
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Gas Concentration'
          ELSEIF( INDEX(BDUM(1:),'water aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 76
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Water Aqueous Concentration'
          ELSEIF( INDEX(BDUM(1:),'air aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 77
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Air Aqueous Concentration'
          ELSEIF( INDEX(BDUM(1:),'integrated water mass').NE.0 ) THEN
            I_OBDT(2,NT) = 191
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Water Mass'
          ELSEIF( INDEX(BDUM(1:),'integrated air mass').NE.0 ) THEN
            I_OBDT(2,NT) = 192
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Air Mass'
          ELSEIF( INDEX(BDUM(1:),'integrated aqueous water').NE.0 ) THEN
            I_OBDT(2,NT) = 194
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Aqueous Water Mass'
          ELSEIF( INDEX(BDUM(1:),'integrated aqueous air').NE.0 ) THEN
            I_OBDT(2,NT) = 195
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Aqueous Air Mass'
          ELSEIF( INDEX(BDUM(1:),'integrated gas water').NE.0 ) THEN
            I_OBDT(2,NT) = 197
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Gas Water Mass'
          ELSEIF( INDEX(BDUM(1:),'integrated gas air').NE.0 ) THEN
            I_OBDT(2,NT) = 198
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Integrated Gas Air Mass'
          ELSEIF( INDEX(BDUM(1:),'solute volumetric conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 1
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Volumetric Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute aqueous conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 2
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Aqueous Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute gas conc').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 3
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Gas Concentration'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute aqueous mol').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 5
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Aqueous Mole Fracton'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'solute gas mol').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 6
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Gas Mole Fracton'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'x solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 8
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': X Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'y solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 9
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Y Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSEIF( INDEX(BDUM(1:),'z solute flux').NE.0 ) THEN
            I_OBDT(2,NT) = 400 + (NSL-1)*33 + 10
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Z Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Reference-Observation Variable: '
     &        // BDUM
            CALL WRMSGS( INDX )
          ENDIF
!
!---      Read and check reference-observation output units  ---
!
          IDFLT = 1
          VARB = 'Reference-Observation Output Units'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_OBDT(NT) )
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),
     &      ': ',C_OBDT(NT)(1:NCH)
          CALL RDOUUN( I_OBDT(2,NT) )
          VAR = 0.D+0
          INDX = 0
          CALL RDUNIT( C_OBDT(NT),VAR,INDX )
!
!---      Read IJK Indices  ---
!
          VARB = 'Reference-Observation IJK Indices'
          CALL RDINT(ISTART,ICOMMA,CHDUM,I_OBDT(4,NT))
          CALL RDINT(ISTART,ICOMMA,CHDUM,I_OBDT(5,NT))
          CALL RDINT(ISTART,ICOMMA,CHDUM,I_OBDT(6,NT))
          IF( I_OBDT(4,NT).LT.1 .OR. I_OBDT(4,NT).GT.IFLD .OR.
     &      I_OBDT(5,NT).LT.1 .OR. I_OBDT(5,NT).GT.JFLD. OR.
     &      I_OBDT(6,NT).LT.1 .OR. I_OBDT(6,NT).GT.KFLD) THEN
            INDX = 4
            CHMSG = 'Out-of-Range Reference-Observation Index'
            CALL WRMSGS( INDX )
          ENDIF
          WRITE(IWR,'(2X,A,$)') 'Reference-Observation Indices: '
          WRITE(FORM3(3:3),'(I1)') ICOUNT(I_OBDT(4,NT))
          WRITE(IWR,'(2X,A,$)') 'I = '
          WRITE(IWR,FORM3) I_OBDT(4,NT)
          WRITE(FORM3(3:3),'(I1)') ICOUNT(I_OBDT(5,NT))
          WRITE(IWR,'(2X,A,$)') 'J = '
          WRITE(IWR,FORM3) I_OBDT(5,NT)
          WRITE(FORM4(3:3),'(I1)') ICOUNT(I_OBDT(6,NT))
          WRITE(IWR,'(2X,A,$)') 'K = '
          WRITE(IWR,FORM4) I_OBDT(6,NT)
!
!---  Surface-rate-observation variable  ---
!
        ELSEIF( INDEX(ADUM(1:),'surface').NE.0 .AND.
     &    ( INDEX(ADUM(1:),'rate').NE.0 .OR.
     &    INDEX(ADUM(1:),'flux').NE.0 ) ) THEN
          I_OBDT(1,NT) = 3
          NSF = NSF + 1
          IF( NSF.GT.LSF ) THEN
            INDX = 5
            CHMSG = 'Number of Surface Flux Domains > Parameter LSF'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Surface-Rate-Observation Variable'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          IVRX = IVR
          IF( INDEX( BDUM(1:),'solute' ).NE.0 ) THEN
            VARB = 'Surface-Rate-Observation Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
            DO 80 NSL = 1,NSOLU
              IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 82
   80       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Surface-Rate-Observ. Solute Name: '
     &        // SOLNM
            CALL WRMSGS( INDX )
   82       CONTINUE
          ENDIF
!
!---      Convert surface-rate-observation variable
!         to an index  ---
!
          VARB = 'Surface-Rate-Observation Variable'
          IVR = IVRX
          I_OBDT(4,NT) = NSF
          IF( INDEX(BDUM(1:),'aqueous').NE.0 ) THEN
            IF( INDEX(BDUM(1:),'volum').NE.0 ) THEN
              I_OBDT(2,NT) = 2
              C_OBDT(NT) = 'm^3/s'
              IUNM = 3
              IUNS = -1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Aqueous Volumetric Flux'
            ELSE
              I_OBDT(2,NT) = 5
              C_OBDT(NT) = 'kg/s'
              IUNKG = 1
              IUNS = -1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Aqueous Mass Flux'
            ENDIF
          ELSEIF( INDEX(BDUM(1:),'gas').NE.0 ) THEN
            IF( INDEX(BDUM(1:),'volum').NE.0 ) THEN
              I_OBDT(2,NT) = 3
              C_OBDT(NT) = 'm^3/s'
              IUNM = 3
              IUNS = -1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Gas Volumetric Flux'
            ELSEIF( INDEX(BDUM(1:),'advective').NE.0 ) THEN
              IF( INDEX(BDUM(1:),'heat').NE.0) THEN
                I_OBDT(2,NT) = 20
                C_OBDT(NT) = 'w'
                IUNM = 2
                IUNKG = 1
                IUNS = -3
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Heat Flux'
              ELSEIF( INDEX(BDUM(1:),'water').NE.0) THEN
                I_OBDT(2,NT) = 21
                C_OBDT(NT) = 'kg/s'
                IUNKG = 1
                IUNS = -1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Water-Mass Flux'
              ELSEIF( INDEX(BDUM(1:),'air').NE.0) THEN
                I_OBDT(2,NT) = 22
                C_OBDT(NT) = 'kg/s'
                IUNKG = 1
                IUNS = -1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Air-Mass Flux'
              ENDIF
            ELSEIF( INDEX(BDUM(1:),'diffusive').NE.0 ) THEN
              IF( INDEX(BDUM(1:),'heat').NE.0) THEN
                I_OBDT(2,NT) = 25
                C_OBDT(NT) = 'w'
                IUNM = 2
                IUNKG = 1
                IUNS = -3
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Heat Flux'
              ELSEIF( INDEX(BDUM(1:),'water').NE.0) THEN
                I_OBDT(2,NT) = 26
                C_OBDT(NT) = 'kg/s'
                IUNKG = 1
                IUNS = -1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Water-Mass Flux'
              ELSEIF( INDEX(BDUM(1:),'air').NE.0) THEN
                I_OBDT(2,NT) = 27
                C_OBDT(NT) = 'kg/s'
                IUNKG = 1
                IUNS = -1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Air-Mass Flux'
              ENDIF
            ELSE
              I_OBDT(2,NT) = 6
              C_OBDT(NT) = 'kg/s'
              IUNKG = 1
              IUNS = -1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Gas Mass Flux'
            ENDIF
          ELSEIF( INDEX(BDUM(1:),'heat').NE.0 ) THEN
            I_OBDT(2,NT) = 2
            C_OBDT(NT) = 'w'
            IUNM = 2
            IUKG = 1
            IUNS = -3
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Heat Flux'
          ELSEIF( INDEX(BDUM(1:),'solute').NE.0 ) THEN
            I_OBDT(2,NT) = 100 + NSL
            C_OBDT(NT) = 'sol/s'
            IUNS = -1
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Flux'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Surface-Rate-Observation Variable: '
     &        // BDUM
            CALL WRMSGS( INDX )
          ENDIF
          ISFT(NSF) = I_OBDT(2,NT)
!
!---      Read and check units  ---
!
          IDFLT = 1
          VARB = 'Surface-Rate-Observation Output Units'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_OBDT(NT) )
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),
     &      ': ',C_OBDT(NT)(1:NCH)
          VAR = 0.D+0
          INDX = 0
          CALL RDUNIT( C_OBDT(NT),VAR,INDX )
!
!---      Read surface-rate-observation orientation  ---
!
          VARB = 'Surface-Rate-Observation Orientation'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          IF( INDEX(ADUM(1:),'west').NE.0) THEN
            ISFD(NSF) = -1
            WRITE(IWR,'(A)') 'X-Direction: West Surface'
          ELSEIF( INDEX(ADUM(1:),'east').NE.0) THEN
            ISFD(NSF) = 1
            WRITE(IWR,'(A)') 'X-Direction: East Surface'
          ELSEIF( INDEX(ADUM(1:),'south').NE.0) THEN
            ISFD(NSF) = -2
            WRITE(IWR,'(A)') 'Y-Direction: South Surface'
          ELSEIF( INDEX(ADUM(1:),'north').NE.0) THEN
            ISFD(NSF) = 2
            WRITE(IWR,'(A)') 'Y-Direction: North Surface'
          ELSEIF( INDEX(ADUM(1:),'bottom').NE.0) THEN
            ISFD(NSF) = -3
            WRITE(IWR,'(A)') 'Z-Direction: Bottom Surface'
          ELSEIF( INDEX(ADUM(1:),'top').NE.0) THEN
            ISFD(NSF) = 3
            WRITE(IWR,'(A)') 'Z-Direction: Top Surface'
          ENDIF
!
!---      Read surface-rate-observation domain  ---
!
          VARB = 'Surface-Rate-Observation Domain: '
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(4,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(5,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(6,NSF))
          ISFC(1,NSF) = MAX( 1,ISFC(1,NSF) )
          ISFC(1,NSF) = MIN( IFLD,ISFC(1,NSF),ISFC(2,NSF) )
          ISFC(2,NSF) = MAX( 1,ISFC(1,NSF),ISFC(2,NSF) )
          ISFC(2,NSF) = MIN( IFLD,ISFC(2,NSF) )
          ISFC(3,NSF) = MAX( 1,ISFC(3,NSF) )
          ISFC(3,NSF) = MIN( JFLD,ISFC(3,NSF),ISFC(4,NSF) )
          ISFC(4,NSF) = MAX( 1,ISFC(3,NSF),ISFC(4,NSF) )
          ISFC(4,NSF) = MIN( JFLD,ISFC(4,NSF) )
          ISFC(5,NSF) = MAX( 1,ISFC(5,NSF) )
          ISFC(5,NSF) = MIN( KFLD,ISFC(5,NSF),ISFC(6,NSF) )
          ISFC(6,NSF) = MAX( 1,ISFC(5,NSF),ISFC(6,NSF) )
          ISFC(6,NSF) = MIN( KFLD,ISFC(6,NSF) )
          WRITE(IWR,'(T3,A)') VARB(1:IVR)
          WRITE(IWR, '(T5,2(A,I6))') 'I = ',ISFC(1,NSF),' to ',
     &      ISFC(2,NSF)
          WRITE(IWR, '(T5,2(A,I6))') 'J = ',ISFC(3,NSF),' to ',
     &      ISFC(4,NSF)
          WRITE(IWR, '(T5,2(A,I6))') 'K = ',ISFC(5,NSF),' to ',
     &      ISFC(6,NSF)
!
!---  Surface-integral-observation variable  ---
!
        ELSEIF( INDEX(ADUM(1:),'surface').NE.0 .AND.
     &    INDEX(ADUM(1:),'integral').NE.0 ) THEN
          I_OBDT(1,NT) = 4
          NSF = NSF + 1
          IF( NSF.GT.LSF ) THEN
            INDX = 5
            CHMSG = 'Number of Surface Flux Domains > Parameter LSF'
            CALL WRMSGS( INDX )
          ENDIF
          VARB = 'Surface-Rate-Observation Variable'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,BDUM)
          IVRX = IVR
          IF( INDEX( BDUM(1:),'solute' ).NE.0 ) THEN
            VARB = 'Surface-Rate-Observation Solute Name'
            CALL RDCHR(ISTART,ICOMMA,NCS,CHDUM,SOLNM)
            DO 90 NSL = 1,NSOLU
              IF( SOLNM.EQ.SOLUT(NSL) ) GOTO 92
   90       CONTINUE
            INDX = 4
            CHMSG = 'Unrecognized Surface-Rate-Observ. Solute Name: '
     &        // SOLNM
            CALL WRMSGS( INDX )
   92       CONTINUE
          ENDIF
!
!---      Convert surface-integral-observation variable
!         to an index  ---
!
          I_OBDT(4,NT) = NSF
          VARB = 'Surface-Integral-Observation Variable'
          IVR = IVRX
          IF( INDEX(BDUM(1:),'aqueous').NE.0 ) THEN
            IF( INDEX(BDUM(1:),'volum').NE.0 ) THEN
              I_OBDT(2,NT) = 2
              C_OBDT(NT) = 'm^3'
              IUNM = 3
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Aqueous Volumetric Flux Integral'
            ELSE
              I_OBDT(2,NT) = 5
              C_OBDT(NT) = 'kg'
              IUNKG = 1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Aqueous Mass Flux Integral'
            ENDIF
          ELSEIF( INDEX(BDUM(1:),'gas').NE.0 ) THEN
            IF( INDEX(BDUM(1:),'volum').NE.0 ) THEN
              I_OBDT(2,NT) = 3
              C_OBDT(NT) = 'm^3'
              IUNM = 3
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Gas Volumetric Flux Integral'
            ELSEIF( INDEX(BDUM(1:),'advective').NE.0 ) THEN
              IF( INDEX(BDUM(1:),'heat').NE.0) THEN
                I_OBDT(2,NT) = 20
                C_OBDT(NT) = 'j'
                IUNM = 2
                IUNKG = 1
                IUNS = -2
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Heat Flux Integral'
              ELSEIF( INDEX(BDUM(1:),'water').NE.0) THEN
                I_OBDT(2,NT) = 21
                C_OBDT(NT) = 'kg'
                IUNKG = 1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Water-Mass Flux Integral'
              ELSEIF( INDEX(BDUM(1:),'air').NE.0) THEN
                I_OBDT(2,NT) = 22
                C_OBDT(NT) = 'kg'
                IUNKG = 1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Air-Mass Flux Integral'
              ENDIF
            ELSEIF( INDEX(BDUM(1:),'diffusive').NE.0 ) THEN
              IF( INDEX(BDUM(1:),'heat').NE.0) THEN
                I_OBDT(2,NT) = 25
                C_OBDT(NT) = 'j'
                IUNM = 2
                IUNKG = 1
                IUNS = -2
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Heat Flux Integral'
              ELSEIF( INDEX(BDUM(1:),'water').NE.0) THEN
                I_OBDT(2,NT) = 26
                C_OBDT(NT) = 'kg'
                IUNKG = 1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Water-Mass Flux Integral'
              ELSEIF( INDEX(BDUM(1:),'air').NE.0) THEN
                I_OBDT(2,NT) = 27
                C_OBDT(NT) = 'kg'
                IUNKG = 1
                WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &            ': Gas Advective Air-Mass Flux Integral'
              ENDIF
            ELSE
              I_OBDT(2,NT) = 6
              C_OBDT(NT) = 'kg'
              IUNKG = 1
              WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &          ': Gas Mass Flux Integral'
            ENDIF
          ELSEIF( INDEX(BDUM(1:),'heat').NE.0 ) THEN
            I_OBDT(2,NT) = 2
            C_OBDT(NT) = 'j'
            IUNM = 2
            IUKG = 1
            IUNS = -2
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Heat Flux Integral'
          ELSEIF( INDEX(BDUM(1:),'solute').NE.0 ) THEN
            I_OBDT(2,NT) = 100 + NSL
            C_OBDT(NT) = 'sol'
            WRITE(IWR,'(/,2X,2A)') VARB(1:IVR),
     &        ': Solute Flux Integral'
            WRITE(IWR,'(2X,2A)') 'Solute Name: ',SOLUT(NSL)
          ELSE
            INDX = 4
            CHMSG = 'Unrecognized Surface-Rate-Observation Variable: '
     &        // BDUM
            CALL WRMSGS( INDX )
          ENDIF
          ISFT(NSF) = I_OBDT(2,NT)
!
!---      Read and check units  ---
!
          IDFLT = 1
          VARB = 'Surface-Integral-Observation Output Units'
          CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,C_OBDT(NT) )
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),
     &      ': ',C_OBDT(NT)(1:NCH)
           VAR = 0.D+0
          INDX = 0
          CALL RDUNIT( C_OBDT(NT),VAR,INDX )
!
!---      Read surface-integral-observation orientation  ---
!
          VARB = 'Surface-Integral-Observation Orientation'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          WRITE(IWR,'(2X,2A,$)') VARB(1:IVR),': '
          IF( INDEX(ADUM(1:),'west').NE.0) THEN
            ISFD(NSF) = -1
            WRITE(IWR,'(A)') 'X-Direction: West Surface'
          ELSEIF( INDEX(ADUM(1:),'east').NE.0) THEN
            ISFD(NSF) = 1
            WRITE(IWR,'(A)') 'X-Direction: East Surface'
          ELSEIF( INDEX(ADUM(1:),'south').NE.0) THEN
            ISFD(NSF) = -2
            WRITE(IWR,'(A)') 'Y-Direction: South Surface'
          ELSEIF( INDEX(ADUM(1:),'north').NE.0) THEN
            ISFD(NSF) = 2
            WRITE(IWR,'(A)') 'Y-Direction: North Surface'
          ELSEIF( INDEX(ADUM(1:),'bottom').NE.0) THEN
            ISFD(NSF) = -3
            WRITE(IWR,'(A)') 'Z-Direction: Bottom Surface'
          ELSEIF( INDEX(ADUM(1:),'top').NE.0) THEN
            ISFD(NSF) = 3
            WRITE(IWR,'(A)') 'Z-Direction: Top Surface'
          ENDIF
!
!---      Read surface-rate-observation domain  ---
!
          VARB = 'Surface-Rate-Observation Domain: '
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(1,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(2,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(3,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(4,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(5,NSF))
          CALL RDINT(ISTART,ICOMMA,CHDUM,ISFC(6,NSF))
          ISFC(1,NSF) = MAX( 1,ISFC(1,NSF) )
          ISFC(1,NSF) = MIN( IFLD,ISFC(1,NSF),ISFC(2,NSF) )
          ISFC(2,NSF) = MAX( 1,ISFC(1,NSF),ISFC(2,NSF) )
          ISFC(2,NSF) = MIN( IFLD,ISFC(2,NSF) )
          ISFC(3,NSF) = MAX( 1,ISFC(3,NSF) )
          ISFC(3,NSF) = MIN( JFLD,ISFC(3,NSF),ISFC(4,NSF) )
          ISFC(4,NSF) = MAX( 1,ISFC(3,NSF),ISFC(4,NSF) )
          ISFC(4,NSF) = MIN( JFLD,ISFC(4,NSF) )
          ISFC(5,NSF) = MAX( 1,ISFC(5,NSF) )
          ISFC(5,NSF) = MIN( KFLD,ISFC(5,NSF),ISFC(6,NSF) )
          ISFC(6,NSF) = MAX( 1,ISFC(5,NSF),ISFC(6,NSF) )
          ISFC(6,NSF) = MIN( KFLD,ISFC(6,NSF) )
          WRITE(IWR,'(T3,A)') VARB(1:IVR)
          WRITE(IWR, '(T5,2(A,I6))') 'I = ',ISFC(1,NSF),' to ',
     &      ISFC(2,NSF)
          WRITE(IWR, '(T5,2(A,I6))') 'J = ',ISFC(3,NSF),' to ',
     &      ISFC(4,NSF)
          WRITE(IWR, '(T5,2(A,I6))') 'K = ',ISFC(5,NSF),' to ',
     &      ISFC(6,NSF)
!
!---    Unrecognized Observed-Data Type  ---
!
        ELSE
          INDX = 2
          CHMSG = 'Unrecognized Observed-Data Type: '//ADUM(1:NCH)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read observed data statistical index  ---
!
        VARB = 'Observed Data Statistical Index'
        CALL RDINT(ISTART,ICOMMA,CHDUM,I_OBDT(3,NT))
        WRITE(FORM1(9:9),'(I1)') ICOUNT( I_OBDT(3,NT) )
        WRITE(IWR,FORM1) VARB(1:IVR),': ',I_OBDT(3,NT)
        IF( I_OBDT(3,NT).LT.1 .OR. I_OBDT(3,NT).GT.1 ) THEN
          INDX = 7
          CHMSG = 'Out of Range Observed Data Statistical Index: '
          IMSG = I_OBDT(3,NT)
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Read observed data statistical parameters  ---
!
        VARB = 'Observed Data Statistic'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(4,NT))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_OBDT(4,NT)
        VARB = 'Observed Data Time Weighting Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(5,NT))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_OBDT(5,NT)
        VARB = 'Observed Data Space Weighting Factor'
        CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDT(6,NT))
        WRITE(IWR,'(2X,2A,1PE11.4)') VARB(1:IVR),': ',R_OBDT(6,NT)
!
!---    Read number of observed data samples
!       or an external file name  ---
!
        ISTART = 1
        CALL RDINPL( CHDUM )
        CALL LCASE( CHDUM )
!
!---    Read observed-data samples from an external file  ---
!
        IF( INDEX( CHDUM(1:),'file').NE.0 ) THEN
          VARB = 'Observed-Data External File Name'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,ADUM)
          CALL RDCHR(ISTART,ICOMMA,NCHF,CHDUM,FDUM)
          WRITE(IWR,'(2X,3A)') VARB(1:IVR),': ',FDUM(1:NCHF)
!
!---      Read external file time units  ---
!
          VARB = 'Observed-Data External File Time Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),': ',UNTS(1:NCH)
          INDX = 0
          IUNS = 1
          VART = 1.D+0
          CALL RDUNIT(UNTS,VART,INDX)
!
!---      Read external file variable units  ---
!
          VARB = 'Observed-Data External File Variable Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),': ',UNTS(1:NCH)
          INDX = 0
          CALL OBDAUNT( NT )
          VARV = 1.D+0
          CALL RDUNIT(UNTS,VARV,INDX)
!
!---      Check that external file exists  ---
!
          INQUIRE( FILE=FDUM(1:NCHF), FORM=CDUM, EXIST=FLG_CHK )
          IF( .NOT.FLG_CHK ) THEN
            INDX = 4
            CHMSG = 'Missing Observed-Data External File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ELSEIF( CDUM.EQ.'UNFORMATTED' ) THEN
            INDX = 4
            CHMSG = 'Unformatted Observed-Data External File: '
     &        // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ENDIF
          OPEN(UNIT=43,FILE=FDUM(1:NCHF),STATUS='OLD',FORM='FORMATTED')
          NS = 0
  100     READ(43,'(A)',END=110) CHDUM
          IF( CHDUM(1:1).EQ.'#' .OR. CHDUM(1:1).EQ.'!' ) GOTO 100
          BACKSPACE(43)
          NS = NS + 1
          IF( NS.GT.LOBDS ) THEN
            INDX = 6
            CHMSG = 'Number of Observed Data Samples > Parameter LOBDS'
     &        // 'External File: ' // FDUM(1:NCHF)
            CALL WRMSGS( INDX )
          ENDIF
          READ(43,*,END=110) R_OBDS(2,NS,NT),R_OBDS(1,NS,NT)
          R_OBDS(2,NS,NT) = R_OBDS(2,NS,NT)*VART
          R_OBDS(1,NS,NT) = R_OBDS(1,NS,NT)*VARV
          GOTO 100
  110     CONTINUE
          NOBDS(NT) = NS
          CLOSE(UNIT=43)
          GOTO 300
        ENDIF
!
!---    Read observed-data samples from input file  ---
!
        VARB = 'Number of Observed Data Samples'
        CALL RDINT(ISTART,ICOMMA,CHDUM,NOBDS(NT))
        IF( NOBDS(NT).GT.LOBDS ) THEN
          INDX = 6
          CHMSG = 'Number of Observed Data Samples > Parameter LOBDS'
          CALL WRMSGS( INDX )
        ENDIF
!
!---    Loop over number of observed data samples  ---
!
        DO 200 NS = 1,NOBDS(NT)
          ISTART = 1
          CALL RDINPL( CHDUM )
          CALL LCASE( CHDUM )
!
!---      Read observed data time and units  ---
!
          VARB = 'Observed-Data Time'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDS(2,NS,NT))
          VARB = 'Observed-Data Time Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4,$)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',R_OBDS(2,NS,NT)
          INDX = 0
          IUNS = 1
          CALL RDUNIT(UNTS,R_OBDS(2,NS,NT),INDX)
          WRITE(IWR,'(A,1PE11.4,A)') ' (',R_OBDS(2,NS,NT),', s)'
          TMOB = MIN( R_OBDS(2,NS,NT),TMOB )
!
!---      Read observed data value and units  ---
!
          VARB = 'Observed Data Value'
          CALL RDDPR(ISTART,ICOMMA,CHDUM,R_OBDS(1,NS,NT))
          VARB = 'Observed Data Value Units'
          CALL RDCHR(ISTART,ICOMMA,NCH,CHDUM,UNTS)
          WRITE(IWR,'(2X,4A,1PE11.4)') VARB(1:IVR),', ',UNTS(1:NCH),
     &    ': ',R_OBDS(1,NS,NT)
          INDX = 0
          CALL OBDAUNT( NT )
          CALL RDUNIT(UNTS,R_OBDS(1,NS,NT),INDX)
  200   CONTINUE
  300 CONTINUE
!
!---  Create out_uc.sto file for UCODE ---
!
      FLNM = 'out_uc'
      IX = 7
      JX = IX+ICOUNT(IOM)-1
      WRITE(FORM2(3:3),'(I1)') ICOUNT(IOM)
      WRITE(FLNM(IX:JX),FORM2) IOM
      IX = JX+1
      JX = JX+4
      FLNM(IX:JX) = '.sto'
      OPEN(UNIT=IOBDSF, FILE=FLNM(1:JX), STATUS='UNKNOWN',
     &  FORM='FORMATTED')
      CLOSE(UNIT=IOBDSF, STATUS='DELETE')
      OPEN(UNIT=IOBDSF, FILE=FLNM(1:JX), STATUS='NEW', FORM='FORMATTED')
!
!---  Create stompx.uni file for UCODE ---
!
      FLNM(IX:JX) = '.uni'
      INQUIRE( FILE=FLNM(1:JX), EXIST=FLG_CHK )
      IF( .NOT. FLG_CHK ) THEN
        FLG_UNI = .TRUE.
      ELSE
        FLG_UNI = .FALSE.
      ENDIF
      IF( FLG_UNI ) THEN
        OPEN(UNIT=IOBDUF, FILE=FLNM(1:JX), STATUS='UNKNOWN',
     &    FORM='FORMATTED')
        CLOSE(UNIT=IOBDUF, STATUS='DELETE')
        OPEN(UNIT=IOBDUF, FILE=FLNM(1:JX), STATUS='NEW',
     &    FORM='FORMATTED')
      ENDIF
!
!---  Create stompx.ext file for UCODE ---
!
      FLNM(IX:JX) = '.ext'
      INQUIRE( FILE=FLNM(1:JX), EXIST=FLG_CHK )
      IF( .NOT. FLG_CHK ) THEN
        FLG_EXT = .TRUE.
      ELSE
        FLG_EXT = .FALSE.
      ENDIF
      IF( FLG_EXT ) THEN
        OPEN(UNIT=IOBDEF, FILE=FLNM(1:JX), STATUS='UNKNOWN',
     &    FORM='FORMATTED')
        CLOSE(UNIT=IOBDEF, STATUS='DELETE')
        OPEN(UNIT=IOBDEF, FILE=FLNM(1:JX), STATUS='NEW',
     &    FORM='FORMATTED')
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDOBDA group ---
!
      RETURN
      END


!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE RDUCODE
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
!     Read UCode control card.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 11 October 2001.
!     Last Modified by MD White, PNNL, 11 October 2001.
!     Last Modified by MD White, PNNL, 10 December 2003.
!     ucode_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE SOLTN
      USE FILES
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 ADUM
      CHARACTER*512 CHDUM
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*10 FORM1
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1
      DATA FORM1 / '(I2,T24,A)' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/RDUCODE'
!
!---  Write card information to output file  ---
!
      CARD = 'UCode Control Card'
      ICD = INDEX( CARD,'  ' )-1
      WRITE(IWR,'(//,3A)') ' ~ ',CARD(1:ICD),': '
!
!---  Write header data to stomp.uni file  ---
!
      IF( FLG_UNI ) THEN
        WRITE(IOBDUF,9001) '#STOMP.UNI FILE FOR UCODE'
        WRITE(IOBDUF,9001) '#'
      ENDIF
!
!---  Read input and write output to stomp.uni file  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
!
!---  Phase  ---
!
      VARB = 'UCode Phase'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Forward Modeling'
      ELSEIF( IVAR.EQ.11 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Calculates Sum-of-Squares'
      ELSEIF( IVAR.EQ.2 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Sensitivities at ' //
     &    'Starting Parameters'
      ELSEIF( IVAR.EQ.22 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Sensitivities at ' //
     &    'Starting Parameters using Centeral Differences'
      ELSEIF( IVAR.EQ.3 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Peform Regression'
      ELSEIF( IVAR.EQ.33 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Calculate Model ' //
     &    'Linearity'
      ELSEIF( IVAR.EQ.44 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Calculate Prediction ' //
     &    'Intervals'
      ELSEIF( IVAR.EQ.45 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Calculate Differences ' //
     &    'and Prediction Intervals'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Phase'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,'#phase'
!
!---  Differencing  ---
!
      VARB = 'UCode Differencing Index'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Forward Differencing'
      ELSEIF( IVAR.EQ.2 ) THEN
        WRITE(IWR,'(/,T2,2A)') VARB(1:IVR),': Central Differencing'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Differencing Index'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#differencing (1=forward [recommended], 2=central)'
!
!---  Tolerance  ---
!
      VARB = 'UCode Tolerance'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR)
      WRITE(IWR,'(T2,2A,1PE11.4)') VARB(1:IVR),': ',VAR
      IF( FLG_UNI ) WRITE(IOBDUF,9002) VAR,
     &  '#tol (0.01 recommended)'
!
!---  Sum-of-Squared Residual Factor  ---
!
      VARB = 'UCode Sum-of-Squared Residual Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR)
      WRITE(IWR,'(T2,2A,1PE11.4)') VARB(1:IVR),': ',VAR
      IF( FLG_UNI ) WRITE(IOBDUF,9002) VAR,
     &  '#tolerance sosr (0.01 or 0.1 [recommended])'
!
!---  Quasi-Newton Updating Index  ---
!
      VARB = 'UCode Quasi-Newton Updating Index'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.0 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': No Quasi-Newton Updating'
      ELSEIF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Apply Quasi-Newton Updating'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Quasi-Newton Updating Index'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#nopt (0=no quasi-Newton updating, ' //
     &  '1=quasi-Newton updating)'
!
!---  Maximum number of iterations  ---
!
      VARB = 'UCode Maximum Number of Iterations'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      WRITE(IWR,'(T2,2A,I4)') VARB(1:IVR),': ',IVAR
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#maximum number of iterations'
!
!---  Maximum change factor  ---
!
      VARB = 'UCode Maximum Change Factor'
      CALL RDDPR(ISTART,ICOMMA,CHDUM,VAR)
      WRITE(IWR,'(T2,2A,1PE11.4)') VARB(1:IVR),': ',VAR
      IF( FLG_UNI ) WRITE(IOBDUF,9002) VAR,
     &  '#maximum fractional parameter change'
!
!---  Read new input line  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
!
!---  Path and Name of Inverse Code  ---
!
      VARB = 'UCode Path and Name of Inverse Code'
      CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
      WRITE(IWR,'(T2,3A)') VARB(1:IVR),': ',ADUM(1:NCH)
      IF( FLG_UNI ) WRITE(IOBDUF,9003) ADUM(1:NCH),
     &  '#path and name of inverse code'
!
!---  Read new input line  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
!
!---  Number of Application Models  ---
!
      VARB = 'UCode Number of Application Models'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      WRITE(IWR,'(T2,2A,I4)') VARB(1:IVR),': ',IVAR
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#number of application models'
!
!---  Application Model Execution Commands  ---
!
      DO 100 N = 1,IVAR
        ISTART = 1
        CALL RDINPL( CHDUM )
        VARB = 'UCode Application Model Execution Commands'
        CALL RDCHR( ISTART,ICOMMA,NCH,CHDUM,ADUM )
        WRITE(IWR,'(T2,3A)') VARB(1:IVR),': ',ADUM(1:NCH)
        IF( FLG_UNI ) WRITE(IOBDUF,9003) ADUM(1:NCH),
     &    '#application model execution commands'
  100 CONTINUE
!
!---  Read new input line  ---
!
      ISTART = 1
      CALL RDINPL( CHDUM )
      CALL LCASE( CHDUM )
!
!---  Scale Sensitivities  ---
!
      VARB = 'UCode Scale Sensitivities Index'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.0 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': No scaling is applied, ' //
     &    'and unscaled sensitivities are printed.'
      ELSEIF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Dimensionless scaled ' //
     &    'sensitivities are printed.'
      ELSEIF( IVAR.EQ.2 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': One-percent scaled ' //
     &    'sensitivities are printed.'
      ELSEIF( IVAR.EQ.3 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Both dimensionless and ' //
     &    'one-percent scaled sensitivities are printed.'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Scale Sensitivities Index'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#scale-sensitivities ( 0=no scaling, 1=dimensionless, ' //
     &  '2=1%, and 3=both 1 and 2)'
!
!---  Print Intermediate Index  ---
!
      VARB = 'UCode Print Intermediate Index'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.0 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': No printing for ' //
     &    'intermediate iterations.'
      ELSEIF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Printing for ' //
     &    'intermediate iterations.'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Print Intermediate Index'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#print intermediate ( 0=no printing, 1=print )'
!
!---  Graph Index  ---
!
      VARB = 'UCode Graph Index'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      IF( IVAR.EQ.0 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Do not print ' //
     &    'post-processing files.'
      ELSEIF( IVAR.EQ.1 ) THEN
        WRITE(IWR,'(T2,2A)') VARB(1:IVR),': Print post-' //
     &    'processing files.'
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized UCode Graph Index'
        IMSG = IVAR
        CALL WRMSGS( INDX )
      ENDIF
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#graph ( 0=no printing, 1=print )'
!
!---  Number of Residual Sets  ---
!
      VARB = 'Number of Residual Sets'
      CALL RDINT(ISTART,ICOMMA,CHDUM,IVAR)
      WRITE(FORM1(3:3),'(I1)') ICOUNT(IVAR)
      IF( FLG_UNI ) WRITE(IOBDUF,FORM1) IVAR,
     &  '#number-residual-sets'
!
!---  Observation header  ---
!
      IF( FLG_UNI ) THEN
        WRITE(IOBDUF,9001) '#'
        WRITE(IOBDUF,9001) '# Observations'
        WRITE(IOBDUF,9001) '# Stat-Flag (0=variance, 1=standard ' //
     &    'deviation, 2=coefficient of variation)'
        WRITE(IOBDUF,9004) '# Obs-Name','Obs-Value','Stat.',
     &    'Stat-Flag','Plot-Symbol'
        WRITE(IOBDUF,9001) '#'
      ENDIF
!
!---  Format statements ---
!
 9001 FORMAT(A)
 9002 FORMAT(F7.4,T24,A)
 9003 FORMAT(A,T24,A)
 9004 FORMAT(A,T12,A,T24,A,T36,A,T48,A)
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of RDUCODE group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WROBDA
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
!     Write computed observed data.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 23 May 2001.
!     Last Modified by MD White, PNNL, 23 May 2001.
!     Last Modified by MD White, PNNL, 10 December 2003.
!     ucode_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE TRNSPT
      USE SOLTN
      USE OUTPU
      USE GRID
      USE FLUXT
      USE FLUXP
      USE FDVT
      USE FDVP
      USE FDVG
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*64 FLNM
      CHARACTER*16 OBDSNM
      CHARACTER*4 FORM1
      REAL*8 VAR(8),XVAR(8),YVAR(8),ZVAR(8)
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1
      DATA FORM1 / '(I )' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WROBDA'
!
!---  Write ".sto" file name at begining of ".ext" file  ---
!
      IF( NOBDP.EQ.0 .AND. FLG_EXT ) THEN
        FLNM = '<out_uc'
        IX = 8
        JX = IX+ICOUNT(IOM)-1
        WRITE(FORM1(3:3),'(I1)') ICOUNT(IOM)
        WRITE(FLNM(IX:JX),FORM1) IOM
        IX = JX+1
        JX = JX+4
        FLNM(IX:JX) = '.sto'
        WRITE(IOBDEF,'(A)') FLNM(1:JX)
      ENDIF
!
!---  Find observed data samples with matching times ---
!
      DO NT = 1,NOBDT
      DO NS = 1,NOBDS(NT)
        IF( ABS(TM-R_OBDS(2,NS,NT)).LT.(1.D+1*EPSL) ) THEN
          NOBDP = NOBDP + 1
          OBDSNM = '                '
!
!---      Field-Observation Output ---
!
          IF( I_OBDT(1,NT).EQ.1 ) THEN
            XOB = R_OBDT(1,NT)
            YOB = R_OBDT(2,NT)
            ZOB = R_OBDT(3,NT)
            NC = 0
            DO KX = 1,2
            DO JX = 1,2
            DO IX = 1,2
              NC = NC + 1
              I = I_OBDT(3+IX,NT)
              J = I_OBDT(5+JX,NT)
              K = I_OBDT(7+KX,NT)
              N = ND(I,J,K)
              XVAR(NC) = XP(N)
              YVAR(NC) = YP(N)
              ZVAR(NC) = ZP(N)
            ENDDO
            ENDDO
            ENDDO
            IF( I_OBDT(2,NT).GT.0 .AND. I_OBDT(2,NT).LE.32 ) THEN
            IF( I_OBDT(2,NT).EQ.1 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PL(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'pl_'
            ELSEIF( I_OBDT(2,NT).EQ.2 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PG(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'pg_'
            ELSEIF( I_OBDT(2,NT).EQ.11 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = SL(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'sl_'
            ELSEIF( I_OBDT(2,NT).EQ.12 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = SG(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'sg_'
            ELSEIF( I_OBDT(2,NT).EQ.15 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PORD(2,N)*SL(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'mcl'
            ELSEIF( I_OBDT(2,NT).EQ.21 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XGW(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'xgw'
            ELSEIF( I_OBDT(2,NT).EQ.22 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XGA(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'xga'
            ELSEIF( I_OBDT(2,NT).EQ.23 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XLW(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'xlw'
            ELSEIF( I_OBDT(2,NT).EQ.24 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XLA(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'xla'
            ELSEIF( I_OBDT(2,NT).EQ.27 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PL(2,N)/RHORL/GRAV + ZP(N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'hhl'
            ELSEIF( I_OBDT(2,NT).EQ.28 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PG(2,N)/RHORL/GRAV + ZP(N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'hhg'
            ELSEIF( I_OBDT(2,NT).EQ.31 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(1./3.)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'krl'
            ELSEIF( I_OBDT(2,NT).EQ.32 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RKG(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'krg'
            ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.32 .AND. I_OBDT(2,NT).LE.64 ) THEN
            IF( I_OBDT(2,NT).EQ.34 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOL(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rhl'
            ELSEIF( I_OBDT(2,NT).EQ.35 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOG(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rhg'
            ELSEIF( I_OBDT(2,NT).EQ.44 ) THEN
              NC = 0
              INDX = 1
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &            THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &            PORD(2,N),PORT(2,N),INDX )
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'tkx'
            ELSEIF( I_OBDT(2,NT).EQ.45 ) THEN
              NC = 0
              INDX = 2
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &            THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &            PORD(2,N),PORT(2,N),INDX )
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'tky'
            ELSEIF( I_OBDT(2,NT).EQ.46 ) THEN
              NC = 0
              INDX = 3
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &            THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &            PORD(2,N),PORT(2,N),INDX )
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'tkz'
            ELSEIF( I_OBDT(2,NT).EQ.51 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPX = NSX(N)
                NQX = NSX(N)+1
                VAR(NC) = 0.5D+0*(UL(1,NPX)+UL(1,NQX))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'ul_'
            ELSEIF( I_OBDT(2,NT).EQ.52 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPY = NSY(N)
                NQY = NSY(N)+IFLD
                VAR(NC) = 0.5D+0*(VL(1,NPY)+VL(1,NQY))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'vl_'
            ELSEIF( I_OBDT(2,NT).EQ.53 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPZ = NSZ(N)
                NQZ = NSZ(N)+IJFLD
                VAR(NC) = 0.5D+0*(WL(1,NPZ)+WL(1,NQZ))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'wl_'
            ELSEIF( I_OBDT(2,NT).EQ.54 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPX = NSX(N)
                NQX = NSX(N)+1
                VAR(NC) = 0.5D+0*(UG(1,NPX)+UG(1,NQX))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'ug_'
            ELSEIF( I_OBDT(2,NT).EQ.55 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPY = NSY(N)
                NQY = NSY(N)+IFLD
                VAR(NC) = 0.5D+0*(VG(1,NPY)+VG(1,NQY))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'vg_'
            ELSEIF( I_OBDT(2,NT).EQ.56 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPZ = NSZ(N)
                NQZ = NSZ(N)+IJFLD
                VAR(NC) = 0.5D+0*(WG(1,NPZ)+WG(1,NQZ))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'wg_'
            ELSEIF( I_OBDT(2,NT).EQ.60 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPX = NSX(N)
                NQX = NSX(N)+1
                VAR(NC) = 0.5D+0*(UQ(1,NPX)+UQ(1,NQX))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'uq_'
            ELSEIF( I_OBDT(2,NT).EQ.61 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPY = NSY(N)
                NQY = NSY(N)+IFLD
                VAR(NC) = 0.5D+0*(VQ(1,NPY)+VQ(1,NQY))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'vq_'
            ELSEIF( I_OBDT(2,NT).EQ.62 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                NPZ = NSZ(N)
                NQZ = NSZ(N)+IJFLD
                VAR(NC) = 0.5D+0*(WQ(1,NPZ)+WQ(1,NQZ))
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'wq_'
            ELSEIF( I_OBDT(2,NT).EQ.63 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = (PL(2,N)-PG(2,N))/RHORL/GRAV
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'mph'
            ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.64 .AND. I_OBDT(2,NT).LE.96 ) THEN
            IF( I_OBDT(2,NT).EQ.73 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOG(2,N)*XGW(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rgw'
            ELSEIF( I_OBDT(2,NT).EQ.74 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOG(2,N)*XGA(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rga'
            ELSEIF( I_OBDT(2,NT).EQ.75 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOL(2,N)*XLW(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rlw'
            ELSEIF( I_OBDT(2,NT).EQ.76 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHOL(2,N)*XLA(2,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rla'
            ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.192 .AND. I_OBDT(2,NT).LE.224 ) THEN
            IF( I_OBDT(2,NT).EQ.206 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PN(1,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'ta_'
            ELSEIF( I_OBDT(2,NT).EQ.207 ) THEN
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PN(3,N)
              ENDDO
              ENDDO
              ENDDO
              OBDSNM(1:3) = 'rha'
            ELSEIF( I_OBDT(2,NT).EQ.208 ) THEN
              OBDSNM(1:3) = 'rna'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PN(5,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.209 ) THEN
              OBDSNM(1:3) = 'uza'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = PN(4,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.213 ) THEN
              OBDSNM(1:3) = 'ts'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XLO(1,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.214 ) THEN
              OBDSNM(1:3) = 'pvs'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = XLO(5,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.215 ) THEN
              OBDSNM(1:3) = 'aev'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHON(2,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.216 ) THEN
              OBDSNM(1:3) = 'pev'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHON(3,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.217 ) THEN
              OBDSNM(1:3) = 'atr'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHON(4,N)
              ENDDO
              ENDDO
              ENDDO
            ELSEIF( I_OBDT(2,NT).EQ.218 ) THEN
              OBDSNM(1:3) = 'ptr'
              NC = 0
              DO KX = 1,2
              DO JX = 1,2
              DO IX = 1,2
                NC = NC + 1
                I = I_OBDT(3+IX,NT)
                J = I_OBDT(5+JX,NT)
                K = I_OBDT(7+KX,NT)
                N = ND(I,J,K)
                VAR(NC) = RHON(5,N)
              ENDDO
              ENDDO
              ENDDO
            ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.400 ) THEN
              NSL = ((I_OBDT(2,NT)-400)/33)+1
              INDX = MOD((I_OBDT(2,N)-400),33)+400
              IF( INDX.EQ.401 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  VAR(NC) = C(N,NSL)
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'c  '
              ELSEIF( INDX.EQ.402 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  IF( (SL(2,N)*PORD(2,N))/EPSL.GT.EPSL ) THEN
                    VAR(NC) = C(N,NSL)*YL(N,NSL)/(SL(2,N)*PORD(2,N))
                  ELSE
                    VAR(NC) = 0.D+0
                  ENDIF
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'cl '
              ELSEIF( INDX.EQ.403 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  IF( (SG(2,N)*PORD(2,N))/EPSL.GT.EPSL ) THEN
                    VAR(NC) = C(N,NSL)*YG(N,NSL)/(SG(2,N)*PORD(2,N))
                  ELSE
                    VAR(NC) = 0.D+0
                  ENDIF
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'cg '
              ELSEIF( INDX.EQ.405 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  VAR(NC) = YL(N,NSL)
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'yl '
              ELSEIF( INDX.EQ.406 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  VAR(NC) = YG(N,NSL)
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'yg '
              ELSEIF( INDX.EQ.408 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  NPX = NSX(N)
                  NQX = NSX(N)+1
                  VAR(NC) = 0.5D+0*(UC(NPX,NSL)+UC(NQX,NSL))
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'uc_'
              ELSEIF( INDX.EQ.409 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  NPY = NSY(N)
                  NQY = NSY(N)+IFLD
                  VAR(NC) = 0.5D+0*(VC(NPY,NSL)+VC(NQY,NSL))
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'vc_'
              ELSEIF( INDX.EQ.410 ) THEN
                NC = 0
                DO KX = 1,2
                DO JX = 1,2
                DO IX = 1,2
                  NC = NC + 1
                  I = I_OBDT(3+IX,NT)
                  J = I_OBDT(5+JX,NT)
                  K = I_OBDT(7+KX,NT)
                  N = ND(I,J,K)
                  NPZ = NSZ(N)
                  NQZ = NSZ(N)+IJFLD
                  VAR(NC) = 0.5D+0*(WC(NPZ,NSL)+WC(NQZ,NSL))
                ENDDO
                ENDDO
                ENDDO
                OBDSNM(1:3) = 'wc_'
              ENDIF
            ENDIF
!
!---        Interpolate field-observation to observation point ---
!
            CALL TRI_LIN( VAR,XVAR,YVAR,ZVAR,XOB,YOB,ZOB,VARX )
!
!---      Reference-observation variable ---
!
          ELSEIF( I_OBDT(1,NT).EQ.2 ) THEN
            I = I_OBDT(4,NT)
            J = I_OBDT(5,NT)
            K = I_OBDT(6,NT)
            N = ND(I,J,K)
            IRNV = I_OBDT(2,NT)
            IF( I_OBDT(2,NT).GT.0 .AND. I_OBDT(2,NT).LE.32 ) THEN
              IF( I_OBDT(2,NT).EQ.1 ) THEN
                OBDSNM(1:3) = 'pl_'
                VARX = PL(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.2 ) THEN
                OBDSNM(1:3) = 'pg_'
                VARX = PG(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.11 ) THEN
                OBDSNM(1:3) = 'sl_'
                VARX = SL(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.12 ) THEN
                OBDSNM(1:3) = 'sg_'
                VARX = SG(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.15 ) THEN
                OBDSNM(1:3) = 'mcl'
                VARX = SL(2,N)*PORD(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.21 ) THEN
                OBDSNM(1:3) = 'xgw'
                VARX = XGW(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.22 ) THEN
                OBDSNM(1:3) = 'xga'
                VARX = XGA(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.23 ) THEN
                OBDSNM(1:3) = 'xlw'
                VARX = XLW(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.24 ) THEN
                OBDSNM(1:3) = 'xlw'
                VARX = XLW(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.27 ) THEN
                OBDSNM(1:3) = 'hhl'
                VARX = PL(2,N)/RHORL/GRAV + ZP(N)
              ELSEIF( I_OBDT(2,NT).EQ.28 ) THEN
                OBDSNM(1:3) = 'hhg'
                VARX = PG(2,N)/RHORL/GRAV + ZP(N)
              ELSEIF( I_OBDT(2,NT).EQ.31 ) THEN
                OBDSNM(1:3) = 'krl'
                VARX = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(1./3.)
              ELSEIF( I_OBDT(2,NT).EQ.32 ) THEN
                OBDSNM(1:3) = 'krg'
                VARX = RKG(2,N)
              ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.32 .AND. I_OBDT(2,NT).LE.64 ) THEN
              IF( I_OBDT(2,NT).EQ.44 ) THEN
                OBDSNM(1:3) = 'tkx'
                INDX = 1
                VARX = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &             THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &             PORD(2,N),PORT(2,N),INDX )
              ELSEIF( I_OBDT(2,NT).EQ.45 ) THEN
                OBDSNM(1:3) = 'tky'
                INDX = 2
                VARX = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &             THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &             PORD(2,N),PORT(2,N),INDX )
              ELSEIF( I_OBDT(2,NT).EQ.46 ) THEN
                OBDSNM(1:3) = 'tkz'
                INDX = 3
                VARX = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &             THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &             PORD(2,N),PORT(2,N),INDX )
              ELSEIF( I_OBDT(2,NT).EQ.51 ) THEN
                OBDSNM(1:3) = 'ul_'
                VARX = UL(1,NSX(N))
              ELSEIF( I_OBDT(2,NT).EQ.52 ) THEN
                OBDSNM(1:3) = 'vl_'
                VARX = VL(1,NSY(N))
              ELSEIF( I_OBDT(2,NT).EQ.53 ) THEN
                OBDSNM(1:3) = 'wl_'
                VARX = WL(1,NSZ(N))
              ELSEIF( I_OBDT(2,NT).EQ.54 ) THEN
                OBDSNM(1:3) = 'ug_'
                VARX = UG(1,NSX(N))
              ELSEIF( I_OBDT(2,NT).EQ.55 ) THEN
                OBDSNM(1:3) = 'vg_'
                VARX = VG(1,NSY(N))
              ELSEIF( I_OBDT(2,NT).EQ.56 ) THEN
                OBDSNM(1:3) = 'wg_'
                VARX = WG(1,NSZ(N))
              ELSEIF( I_OBDT(2,NT).EQ.60 ) THEN
                OBDSNM(1:3) = 'uq_'
                VARX = UQ(1,NSX(N))
              ELSEIF( I_OBDT(2,NT).EQ.61 ) THEN
                OBDSNM(1:3) = 'vq_'
                VARX = VQ(1,NSY(N))
              ELSEIF( I_OBDT(2,NT).EQ.62 ) THEN
                OBDSNM(1:3) = 'wq_'
                VARX = WQ(1,NSZ(N))
              ELSEIF( I_OBDT(2,NT).EQ.63 ) THEN
                OBDSNM(1:3) = 'mph'
                VARX = (PL(2,N)-PG(2,N))/RHORL/GRAV
              ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.64 .AND. I_OBDT(2,NT).LE.96 ) THEN
              IF( I_OBDT(2,NT).EQ.73 ) THEN
                OBDSNM(1:3) = 'rgw'
                VARX = XGW(2,N)*RHOG(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.74 ) THEN
                OBDSNM(1:3) = 'rga'
                VARX = XGA(2,N)*RHOG(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.76 ) THEN
                OBDSNM(1:3) = 'rla'
                VARX = XLA(2,N)*RHOL(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.77 ) THEN
                OBDSNM(1:3) = 'rla'
                VARX = XLA(2,N)*RHOL(2,N)
              ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.160 .AND. I_OBDT(2,NT).LE.192 ) THEN
              IF( I_OBDT(2,NT).EQ.191 ) THEN
                OBDSNM(1:3) = 'imw'
                VARX = 0.D+0
                DO 1191 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1191
                  VARX = VARX + PORD(2,L)*VOL(L)*(XLW(2,L)*
     &              SL(2,L)*RHOL(2,L) + XGW(2,L)*SG(2,L)*RHOG(2,L))
 1191           CONTINUE
              ELSEIF( I_OBDT(2,NT).EQ.192 ) THEN
                OBDSNM(1:3) = 'ima'
                VARX = 0.D+0
                DO 1192 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1192
                  VARX = VARX + PORD(2,L)*VOL(L)*(XLA(2,L)*
     &              SL(2,L)*RHOL(2,L) + XGA(2,L)*SG(2,L)*RHOG(2,L))
 1192           CONTINUE
              ELSEIF( I_OBDT(2,NT).EQ.194 ) THEN
                OBDSNM(1:3) = 'ilw'
                VARX = 0.D+0
                DO 1194 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1194
                  VARX = VARX + PORD(2,L)*VOL(L)*XLW(2,L)*
     &              SL(2,L)*RHOL(2,L)
 1194           CONTINUE
              ELSEIF( I_OBDT(2,NT).EQ.195 ) THEN
                OBDSNM(1:3) = 'ila'
                VARX = 0.D+0
                DO 1195 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1195
                  VARX = VARX + PORD(2,L)*VOL(L)*XLA(2,L)*
     &              SL(2,L)*RHOL(2,L)
 1195           CONTINUE
              ELSEIF( I_OBDT(2,NT).EQ.197 ) THEN
                OBDSNM(1:3) = 'igw'
                VARX = 0.D+0
                DO 1197 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1197
                  VARX = VARX + PORD(2,L)*VOL(L)*XGW(2,L)*
     &              SG(2,L)*RHOG(2,L)
 1197           CONTINUE
              ELSEIF( I_OBDT(2,NT).EQ.198 ) THEN
                OBDSNM(1:3) = 'iga'
                VARX = 0.D+0
                DO 1198 L = 1,NFLD
                  IF( IXP(L).EQ.0 ) GOTO 1198
                  VARX = VARX + PORD(2,L)*VOL(L)*XGA(2,L)*
     &              SG(2,L)*RHOG(2,L)
 1198           CONTINUE
              ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.192 .AND. I_OBDT(2,NT).LE.224 ) THEN
              IF( I_OBDT(2,NT).EQ.206 ) THEN
                OBDSNM(1:3) = 'ta_'
                VARX = PN(1,N)
              ELSEIF( I_OBDT(2,NT).EQ.207 ) THEN
                OBDSNM(1:3) = 'rha'
                VARX = PN(3,N)
              ELSEIF( I_OBDT(2,NT).EQ.208 ) THEN
                OBDSNM(1:3) = 'rna'
                VARX = PN(5,N)
              ELSEIF( I_OBDT(2,NT).EQ.209 ) THEN
                OBDSNM(1:3) = 'uza'
                VARX = PN(4,N)
              ELSEIF( I_OBDT(2,NT).EQ.213 ) THEN
                OBDSNM(1:3) = 'ts'
                VARX = XLO(1,N)
              ELSEIF( I_OBDT(2,NT).EQ.214 ) THEN
                OBDSNM(1:3) = 'pvs'
                VARX = XLO(5,N)
              ELSEIF( I_OBDT(2,NT).EQ.215 ) THEN
                OBDSNM(1:3) = 'aev'
                VARX = RHON(2,N)
              ELSEIF( I_OBDT(2,NT).EQ.216 ) THEN
                OBDSNM(1:3) = 'pev'
                VARX = RHON(3,N)
              ELSEIF( I_OBDT(2,NT).EQ.217 ) THEN
                OBDSNM(1:3) = 'atr'
                VARX = RHON(4,N)
              ELSEIF( I_OBDT(2,NT).EQ.218 ) THEN
                OBDSNM(1:3) = 'ptr'
                VARX = RHON(5,N)
              ENDIF
            ENDIF
            IF( I_OBDT(2,NT).GT.400 ) THEN
              NSL = ((I_OBDT(2,NT)-400)/33)+1
              INDX = MOD((I_OBDT(2,N)-400),33)+400
              IF( INDX.EQ.401 ) THEN
                VARX = C(N,NSL)
                OBDSNM(1:3) = 'c  '
              ELSEIF( INDX.EQ.402 ) THEN
                IF( (SL(2,N)*PORD(2,N))/EPSL.GT.EPSL ) THEN
                  VARX = C(N,NSL)*YL(N,NSL)/(SL(2,N)*PORD(2,N))
                ELSE
                  VARX = 0.D+0
                ENDIF
                OBDSNM(1:3) = 'cl '
              ELSEIF( INDX.EQ.403 ) THEN
                IF( (SG(2,N)*PORD(2,N))/EPSL.GT.EPSL ) THEN
                  VARX = C(N,NSL)*YG(N,NSL)/(SG(2,N)*PORD(2,N))
                ELSE
                  VARX = 0.D+0
                ENDIF
                OBDSNM(1:3) = 'cg '
              ELSEIF( INDX.EQ.405 ) THEN
                VARX = YL(N,NSL)
                OBDSNM(1:3) = 'yl '
              ELSEIF( INDX.EQ.406 ) THEN
                VARX = YG(N,NSL)
                OBDSNM(1:3) = 'yg '
              ELSEIF( INDX.EQ.408 ) THEN
                VARX = 0.5D+0*(UC(NPX,NSL)+UC(NQX,NSL))
                OBDSNM(1:3) = 'uc_'
              ELSEIF( INDX.EQ.409 ) THEN
                VARX = 0.5D+0*(VC(NPY,NSL)+VC(NQY,NSL))
                OBDSNM(1:3) = 'vc_'
              ELSEIF( INDX.EQ.410 ) THEN
                VARX = 0.5D+0*(WC(NPZ,NSL)+WC(NQZ,NSL))
                OBDSNM(1:3) = 'wc_'
              ENDIF
            ENDIF
!
!---      Surface-rate-observation variable ---
!
          ELSEIF( I_OBDT(1,NT).EQ.3 ) THEN
            IF( I_OBDT(2,NT).EQ.2 ) THEN
              OBDSNM(1:3) = 'vlr'
            ELSEIF( I_OBDT(2,NT).EQ.2 ) THEN
              OBDSNM(1:3) = 'qr'
            ELSEIF( I_OBDT(2,NT).EQ.3 ) THEN
              OBDSNM(1:3) = 'vgr'
            ELSEIF( I_OBDT(2,NT).EQ.5 ) THEN
              OBDSNM(1:3) = 'mlr'
            ELSEIF( I_OBDT(2,NT).EQ.6 ) THEN
              OBDSNM(1:3) = 'mgr'
            ELSEIF( I_OBDT(2,NT).EQ.20 ) THEN
              OBDSNM(1:3) = 'aqr'
            ELSEIF( I_OBDT(2,NT).EQ.21 ) THEN
              OBDSNM(1:3) = 'awr'
            ELSEIF( I_OBDT(2,NT).EQ.22 ) THEN
              OBDSNM(1:3) = 'aar'
            ELSEIF( I_OBDT(2,NT).EQ.25 ) THEN
              OBDSNM(1:3) = 'dqr'
            ELSEIF( I_OBDT(2,NT).EQ.26 ) THEN
              OBDSNM(1:3) = 'dwr'
            ELSEIF( I_OBDT(2,NT).EQ.27 ) THEN
              OBDSNM(1:3) = 'dar'
            ELSEIF( I_OBDT(2,NT).GT.(100) ) THEN
              OBDSNM(1:3) = 'cr_'
            ENDIF
            VARX = SF(1,I_OBDT(4,NT))
!
!---      Surface-integral-observation variable ---
!
          ELSEIF( I_OBDT(1,NT).EQ.4 ) THEN
            IF( I_OBDT(2,NT).EQ.2 ) THEN
              OBDSNM(1:3) = 'vli'
            ELSEIF( I_OBDT(2,NT).EQ.2 ) THEN
              OBDSNM(1:3) = 'qi'
            ELSEIF( I_OBDT(2,NT).EQ.3 ) THEN
              OBDSNM(1:3) = 'vgi'
            ELSEIF( I_OBDT(2,NT).EQ.5 ) THEN
              OBDSNM(1:3) = 'mli'
            ELSEIF( I_OBDT(2,NT).EQ.6 ) THEN
              OBDSNM(1:3) = 'mgi'
            ELSEIF( I_OBDT(2,NT).EQ.20 ) THEN
              OBDSNM(1:3) = 'aqi'
            ELSEIF( I_OBDT(2,NT).EQ.21 ) THEN
              OBDSNM(1:3) = 'awi'
            ELSEIF( I_OBDT(2,NT).EQ.22 ) THEN
              OBDSNM(1:3) = 'aai'
            ELSEIF( I_OBDT(2,NT).EQ.25 ) THEN
              OBDSNM(1:3) = 'dqi'
            ELSEIF( I_OBDT(2,NT).EQ.26 ) THEN
              OBDSNM(1:3) = 'dwi'
            ELSEIF( I_OBDT(2,NT).EQ.27 ) THEN
              OBDSNM(1:3) = 'dai'
            ELSEIF( I_OBDT(2,NT).GT.(100) ) THEN
              OBDSNM(1:3) = 'ci_'
            ENDIF
            VARX = SF(2,I_OBDT(4,NT))
          ENDIF
!
!---      Convert to STOMP observation data to output units ---
!
          IF( C_OBDT(NT).NE.'null' ) THEN
            CALL OBDAUNT( NT )
            INDX = 1
            CALL RDUNIT(C_OBDT(NT),VARX,INDX)
          ENDIF
!
!---      Create a unique sample name ---
!
          NCH = INDEX( OBDSNM(1:),'  ')
          IC = ICOUNT(NOBDP)
          WRITE(FORM1(3:3),'(I1)') IC
          IF( IC.EQ.1 ) THEN
            OBDSNM(NCH:NCH+4) = '00000'
            WRITE(OBDSNM(NCH+5:),FORM1) NOBDP
          ELSEIF( IC.EQ.2 ) THEN
            OBDSNM(NCH:NCH+3) = '0000'
            WRITE(OBDSNM(NCH+4:),FORM1) NOBDP
          ELSEIF( IC.EQ.3 ) THEN
            OBDSNM(NCH:NCH+2) = '000'
            WRITE(OBDSNM(NCH+3:),FORM1) NOBDP
          ELSEIF( IC.EQ.4 ) THEN
            OBDSNM(NCH:NCH+1) = '00'
            WRITE(OBDSNM(NCH+2:),FORM1) NOBDP
          ELSEIF( IC.EQ.5 ) THEN
            OBDSNM(NCH:NCH) = '0'
            WRITE(OBDSNM(NCH+1:),FORM1) NOBDP
          ELSEIF( IC.EQ.6 ) THEN
            WRITE(OBDSNM(NCH:),FORM1) NOBDP
          ENDIF
          NCH = INDEX( OBDSNM(1:),'  ')-1
!
!---      Write STOMP observation data to out_ucx.sto file  ---
!
          WRITE(IOBDSF,9001) OBDSNM(1:NCH),VARX
!
!---      Convert to field observation data to output units ---
!
          VARX = R_OBDS(1,NS,NT)
          IF( C_OBDT(NT).NE.'null' ) THEN
            CALL OBDAUNT( NT )
            INDX = 1
            CALL RDUNIT(C_OBDT(NT),VARX,INDX)
          ENDIF
!
!---      Write field observation data to out_ucx.uni file  ---
!
          IF( FLG_UNI ) WRITE(IOBDUF,9002) OBDSNM(1:NCH),VARX,
     &      R_OBDT(4,NT),I_OBDT(3,NT),NT
!
!---      Write correlation between field and STOMP observation data
!         to out_ucx.ext file  ---
!
          IF( FLG_EXT ) THEN
            IF( NOBDP.EQ.1 ) WRITE(IOBDEF,'(A)') '#'
            WRITE(IOBDEF,'(A,T5,A)') 'o',OBDSNM(1:NCH)
            WRITE(IOBDEF,'(3A)') '/',OBDSNM(1:NCH),'/'
            WRITE(IOBDEF,'(A)') 'C12_23'
            WRITE(IOBDEF,'(A)') '#'
          ENDIF
        ENDIF
      ENDDO
      ENDDO
!
!---  Format statements ---
!
 9001 FORMAT(A,T12,1PE11.4)
 9002 FORMAT(A,T12,1PE11.4,T24,1PE11.4,T36,I2,T48,I4)
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WROBDA group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TRI_LIN( VAR,XVAR,YVAR,ZVAR,XOB,YOB,ZOB,VARX )
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
!     Tri-linear interpolation on a structured grid.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 24 May 2001.
!     Last Modified by MD White, PNNL, 24 May 2001.
!     Last Modified by MD White, PNNL, 10 December 2003.
!     ucode_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
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
      REAL*8 VAR(8),XVAR(8),YVAR(8),ZVAR(8)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TRI_LIN'
      UX = (XOB-XVAR(1))/(XVAR(2)-XVAR(1)+EPSL)
      VX = (YOB-YVAR(1))/(YVAR(3)-YVAR(1)+EPSL)
      WX = (ZOB-ZVAR(1))/(ZVAR(5)-ZVAR(1)+EPSL)
      VARX = (1.D+0-UX)*(1.D+0-VX)*(1.D+0-WX)*VAR(1) +
     &  UX*(1.D+0-VX)*(1.D+0-WX)*VAR(2) +
     &  (1.D+0-UX)*VX*(1.D+0-WX)*VAR(3) +
     &  UX*VX*(1.D+0-WX)*VAR(4) +
     &  (1.D+0-UX)*(1.D+0-VX)*WX*VAR(5) +
     &  UX*(1.D+0-VX)*WX*VAR(6) +
     &  (1.D+0-UX)*VX*WX*VAR(7) +
     &  UX*VX*WX*VAR(8)
!
!---  Reset subroutine string sequence  ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TRI_LIN group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE OBDAUNT( NT )
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
!     Observed data units.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 31 May 2001.
!     Last Modified by MD White, PNNL, 31 May 2001.
!     Last Modified by MD White, PNNL, 10 December 2003.
!     ucode_gt.F https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
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
      SUB_LOG(ISUB_LOG) = '/OBDAUNT'
!
!---  Field-observation units  ---
!
      IF( I_OBDT(1,NT).EQ.1 ) THEN
        CALL RDOUUN( I_OBDT(2,NT) )
!
!---  Reference node  ---
!
      ELSEIF( I_OBDT(1,NT).EQ.2 ) THEN
        CALL RDOUUN( I_OBDT(2,NT) )
!
!---  Surface flux rate  ---
!
      ELSEIF( I_OBDT(1,NT).EQ.3 ) THEN
        CALL RDSFUN( I_OBDT(2,NT) )
!
!---  Surface flux integral  ---
!
      ELSEIF( I_OBDT(1,NT).EQ.4 ) THEN
        INDX = -I_OBDT(2,NT)
        CALL RDSFUN( INDX )
      ELSE
        INDX = 7
        CHMSG = 'Unrecognized Observed-Data Type: '
        IMSG = I_OBDT(1,NT)
        CALL WRMSGS( INDX )
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of OBDAUNT group ---
!
      RETURN
      END

