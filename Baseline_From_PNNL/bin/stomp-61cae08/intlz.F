!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE INTLZ
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
!     Initialize all variables contained in common blocks.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, December 1992.
!     Last Modified by MD White, PNNL, October 15, 1997.
!     Last Modified by MD White, PNNL, 1 August 2002.
!     Last Modified by SK Wurstner, PNNL, December 04, 2007.
#ifdef sac
!
!     SAC Modification by WE Nichols, PNNL, June 27, 2003.
#endif
!     intlz.F 1344 2020-07-28 22:36:16Z d3c002 https://stash.pnnl.gov/scm/stomp/stomp.git V3.0
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE TRNSPT
      USE TABL
      USE SPILL
      USE SOURC
      USE SOLTN
#ifdef sac
      USE SACTM
#endif
#ifdef sac
      USE SACIN
#endif
#ifdef sac
      USE SACEV
#endif
      USE REACT
      USE PORMED
      USE POINTE
      USE NCG_PT
      USE NAPL
      USE JACOB
      USE HYST
      USE GRID
      USE FLUXT
      USE FLUXS
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FLUXD
      USE FILES
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVI
      USE FDVH
      USE FDVGC
      USE FDVG
      USE FDVD
      USE FDVA
      USE DUAL_POR
      USE CONST
      USE CCP
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      LOGICAL FLG_EX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/INTLZ'
!
!---  [ BCV, BCVP, BCVT, BCVG, BCVN, BCVI, BCVS, BCVA ]
!     Boundary condition variables  ---
!
      CALL IN_BOUN
!
!---  [ CONST ] Constants  ---
!
      SMALL = 1.D-20
      BIG = 1.D+20
      ZERO = 0.D+0
      THIRD = 1.D+0/3.D+0
      ONE = 1.D+0
      EPSL = 1.D-14
      GRAV = 9.81D+0
      TABS = 273.15D+0
      TMX = 374.14D+0
      TMN = 0.01D+0
      PATM = 101325.D+0
      PMX = 2.212D+8
      PMN = 6.1125D+2
      RCU = 8.31434D+3
      RCW = 461.52D+0
      RCA = 287.0D+0
      RHORL = 998.32142721500441D+0
      RHORG = 1.199D+0
      VISRL = 1.0176489259594200D-3
      VISRG = 1.82D-5
      TSPRF = 293.15D+0
      WTMW = 18.015D+0
      WTMA = 28.97D+0
      WTMN = 28.01348D+0
      WTMS = 58.4428D+0
      TCRW = 647.096D+0
      PCRW = 22.064D+6
      THKRW = 0.6068
      THKRA = 26.1D-3
      ZCRW = 0.235D+0
      VCRW = 57.1075D+0
      VCRA = 86.2269D+0
      PAFW = 0.344D+0
      DPMW = 1.8D+0
      RHOGI = 0.D+0
      RHOLI = 0.D+0
      VISGI = 0.D+0
      VISLI = 0.D+0
      VISNI = 0.D+0
      FSFLA = 1.D+0
      TBW = 373.2D+0
      TBA = 83.35D+0
      DFGAC = 0.D+0
      DFGOC = 0.D+0
      DFGNC = 0.D+0
      DFGWC = 0.D+0
      DFLAC = 0.D+0
      DFLOC = 0.D+0
      DFLNC = 0.D+0
      DFLSC = 0.D+0
      DFNAC = 0.D+0
      DFNOC = 0.D+0
      DFNNC = 0.D+0
      DFNWC = 0.D+0
      PTPS(1) = 0.0765D+0
      PTPS(2) = 0.2664D+0
      PTPS(3) = 0.D+0
      PTPS(4) = 0.00127D+0
      PTPS(5) = 0.10898D+0
      HCAW = 6.7419D+9
      SUFW = 72.8D-3
      HDOD = 1.D+9/RHORL/9.81D+0
      GPI = 3.1415926536D+0
      TENTH = 1.D-1
      TOLN = LOG(1.D+1)
      ISMALL = -32000
      IBIG = 32000
      IPTPS(1) = 1
      IPTPS(2) = 1
      IPTPS(3) = 1
      IPTPS(4) = 1
!
!---  [ CRTCL ] Compound critical properties  ---
!
      DO 400 L = 1,LCMP
        TC(L) = 0.D+0
        PC(L) = 0.D+0
        ZRA(L) = 0.D+0
        VMC(L) = 0.D+0
        PAF(L) = 0.D+0
        WTM(L) = 0.D+0
  400 CONTINUE
      DO L = 1,LNGC
      DO K = 1,70
        GCPP(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO 404 K = 1,70
        GWPP(K) = 0.D+0
  404 CONTINUE
!
!---  [ FILES ] External file names and unit numbers ---
!
      IRD = 21
      IWR = 22
      IPL = 23
      IRS = 24
      IRS_BF = 34
      ISF(1) = 25
#ifdef sac
      IBL = 26
      IVD = 27
      IRL = 28
      IRC = 29
      IRG = 30
#endif
      M = 2
      N = LSF
      IF( LSF.GE.2 ) THEN
        DO 410 L = M,N
          ISF(L) = 49 + L
  410   CONTINUE
      ENDIF
      ISC = 6
      FNRD = 'input'
      FNWR = 'output'
      FNPL = 'plot'
      FNRS = 'restart'
      FNRS_BF = 'restart_bf'
      DO 500 L = 1,LSF
        FNSF(L) = 'surface'
  500 CONTINUE
      FNSR = 'screen'
#ifdef sac
      FNBL = 'balance'
      FNVD = 'remediate'
      FNRL = 'release'
      FNRC = 'recharge'
      FNRG = 'airborne'
#endif
      OPEN(UNIT=IWR, FILE=FNWR, STATUS='UNKNOWN', FORM='FORMATTED')
      CLOSE(UNIT=IWR,STATUS='DELETE')
      OPEN(UNIT=IWR, FILE=FNWR, STATUS='NEW', FORM='FORMATTED')
      INQUIRE( FILE=FNRD, EXIST=FLG_EX )
      IF( .NOT.FLG_EX ) THEN
        INDX = 4
        NCH = INDEX( FNRD(1:),'  ' )-1
        CHMSG = 'Nonexistent Input File: ' // FNRD(1:NCH)
        CARD = 'Simulation Title Card'
        CALL WRMSGS( INDX )
      ELSE
        OPEN(UNIT=IRD, FILE=FNRD, STATUS='OLD', FORM='FORMATTED')
      ENDIF
!
!---  [ FLUXP ] Primary flux variables  ---
!
      DO L = 1,LSX
      DO K = 1,LSFV
        UL(K,L) = 0.D+0
        UDLO(K,L) = 0.D+0
        UDLA(K,L) = 0.D+0
        UDLW(K,L) = 0.D+0
        UG(K,L) = 0.D+0
        UGA(K,L) = 0.D+0
        UDGA(K,L) = 0.D+0
        UDGO(K,L) = 0.D+0
        UDGW(K,L) = 0.D+0
        ULA(K,L) = 0.D+0
        ULW(K,L) = 0.D+0
        UDS(K,L) = 0.D+0
        UGW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSXN2
      DO K = 1,LSFV
        UGN(K,L) = 0.D+0
        ULN(K,L) = 0.D+0
        UDGN(K,L) = 0.D+0
        UDLN(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSY
      DO K = 1,LSFV
        VL(K,L) = 0.D+0
        VDLO(K,L) = 0.D+0
        VDLA(K,L) = 0.D+0
        VDLW(K,L) = 0.D+0
        VG(K,L) = 0.D+0
        VGA(K,L) = 0.D+0
        VDGA(K,L) = 0.D+0
        VDGO(K,L) = 0.D+0
        VDGW(K,L) = 0.D+0
        VLA(K,L) = 0.D+0
        VLW(K,L) = 0.D+0
        VDS(K,L) = 0.D+0
        VGW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYN2
      DO K = 1,LSFV
        VGN(K,L) = 0.D+0
        VLN(K,L) = 0.D+0
        VDGN(K,L) = 0.D+0
        VDLN(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZ
      DO K = 1,LSFV
        WL(K,L) = 0.D+0
        WDLO(K,L) = 0.D+0
        WDLA(K,L) = 0.D+0
        WDLW(K,L) = 0.D+0
        WG(K,L) = 0.D+0
        WGA(K,L) = 0.D+0
        WDGA(K,L) = 0.D+0
        WDGO(K,L) = 0.D+0
        WDGW(K,L) = 0.D+0
        WLA(K,L) = 0.D+0
        WLW(K,L) = 0.D+0
        WDS(K,L) = 0.D+0
        WGW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZN2
      DO K = 1,LSFV
        WGN(K,L) = 0.D+0
        WLN(K,L) = 0.D+0
        WDGN(K,L) = 0.D+0
        WDLN(K,L) = 0.D+0
      ENDDO
      ENDDO
!
!---  [ FLUXT ] Energy flux variables  ---
!
      DO L = 1,LSXT
      DO K = 1,LSFV
        UQ(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYT
      DO K = 1,LSFV
        VQ(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZT
      DO K = 1,LSFV
        WQ(K,L) = 0.D+0
      ENDDO
      ENDDO
!
!---  [ FLUXN ] NAPL flux variables  ---
!
      DO L = 1,LSXN
      DO K = 1,LSFV
        UN(K,L) = 0.D+0
        UNA(K,L) = 0.D+0
        UDNA(K,L) = 0.D+0
        UDNO(K,L) = 0.D+0
        UDNW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYN
      DO K = 1,LSFV
        VN(K,L) = 0.D+0
        VNA(K,L) = 0.D+0
        VDNA(K,L) = 0.D+0
        VDNO(K,L) = 0.D+0
        VDNW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZN
      DO K = 1,LSFV
        WN(K,L) = 0.D+0
        WNA(K,L) = 0.D+0
        WDNA(K,L) = 0.D+0
        WDNO(K,L) = 0.D+0
        WDNW(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSXN2
      DO K = 1,LSFV
        UDNN(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYN2
      DO K = 1,LSFV
        VDNN(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZN2
      DO K = 1,LSFV
        WDNN(K,L) = 0.D+0
      ENDDO
      ENDDO
!
!---  [ FLUXS ] Salt/Surfactant flux variables  ---
!
      DO L = 1,LSXS
      DO K = 1,LSFV
        US(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYS
      DO K = 1,LSFV
        VS(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZS
      DO K = 1,LSFV
        WS(K,L) = 0.D+0
      ENDDO
      ENDDO
!
!---  [ FLUXD ] Dissolved-oil flux variables  ---
!
      DO L = 1,LSXN
      DO K = 1,LSFV
        UGO(K,L) = 0.D+0
        ULO(K,L) = 0.D+0
        UNO(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSYN
      DO K = 1,LSFV
        VGO(K,L) = 0.D+0
        VLO(K,L) = 0.D+0
        VNO(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO L = 1,LSZN
      DO K = 1,LSFV
        WGO(K,L) = 0.D+0
        WLO(K,L) = 0.D+0
        WNO(K,L) = 0.D+0
      ENDDO
      ENDDO
!
!---  [ FLUXGC ] Gas component flux variables  ---
!
      DO L = 1,LSXGC
      DO K = 1,LSFV
      DO J = 1,LNGC
        UGC(J,K,L) = 0.D+0
        UDGC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSYGC
      DO K = 1,LSFV
      DO J = 1,LNGC
        VGC(J,K,L) = 0.D+0
        VDGC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSZGC
      DO K = 1,LSFV
      DO J = 1,LNGC
        WGC(J,K,L) = 0.D+0
        WDGC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSXLC
      DO K = 1,LSFV
      DO J = 1,LNGC
        ULC(J,K,L) = 0.D+0
        UDLC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSYLC
      DO K = 1,LSFV
      DO J = 1,LNGC
        VLC(J,K,L) = 0.D+0
        VDLC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSZLC
      DO K = 1,LSFV
      DO J = 1,LNGC
        WLC(J,K,L) = 0.D+0
        WDLC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSXNC
      DO K = 1,LSFV
      DO J = 1,LNGC
        UNC(J,K,L) = 0.D+0
        UDNC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSYNC
      DO K = 1,LSFV
      DO J = 1,LNGC
        VNC(J,K,L) = 0.D+0
        VDNC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
      DO L = 1,LSZNC
      DO K = 1,LSFV
      DO J = 1,LNGC
        WNC(J,K,L) = 0.D+0
        WDNC(J,K,L) = 0.D+0
      ENDDO
      ENDDO
      ENDDO
!
!---  [ GRID ] Grid variables  ---
!
      DO 715 L = 1,LFX
        RP(L) = 0.D+0
  715 CONTINUE
      DO 730 L = 1,LSX
        DO K = 1,3
          SRFNX(K,L) = 0.D+0
        ENDDO
        AFX(L) = 0.D+0
        DXGP(L) = 0.D+0
        GRVX(L) = 0.D+0
  730 CONTINUE
      DO 735 L = 1,LSY
        DO K = 1,3
          SRFNY(K,L) = 0.D+0
        ENDDO
        AFY(L) = 0.D+0
        DYGP(L) = 0.D+0
        GRVY(L) = 0.D+0
  735 CONTINUE
      DO 740 L = 1,LSZ
        DO K = 1,3
          SRFNZ(K,L) = 0.D+0
        ENDDO
        AFZ(L) = 0.D+0
        DZGP(L) = 0.D+0
        GRVZ(L) = 0.D+0
  740 CONTINUE
      DO L = 1,LBR
        INP(L) = 0
      ENDDO
      DO 746 L = 1,LFD
        DO 742 K = 1,8
          XE(K,L) = 0.D+0
          YE(K,L) = 0.D+0
          ZE(K,L) = 0.D+0
  742   CONTINUE
        DO K = 1,6
          INBS(K,L) = 0
        ENDDO
        DO K = 1,5
          IBR(K,L) = 0
        ENDDO
        DO K = 1,6
          DO J = 1,5
            ICM(J,K,L) = 0
            ICMS(J,K,L) = 0
          ENDDO
        ENDDO
        XP(L) = 0.D+0
        YP(L) = 0.D+0
        ZP(L) = 0.D+0
        DXGF(L) = 1.D+0
        DYGF(L) = 1.D+0
        DZGF(L) = 1.D+0
        VOL(L) = 0.D+0
        IXP(L) = 0
        IZ(L) = 0
        IZ2(L) = 0
        ID(L) = 0
        JD(L) = 0
        KD(L) = 0
        NSX(L) = 0
        NSY(L) = 0
        NSZ(L) = 0
        NSSX(L) = 0
        NSSY(L) = 0
        NSSZ(L) = 0
        IXF(L) = 0
  746 CONTINUE
      MFX_MPI = LFX_MPI+2
      IF( LFX_MPI.EQ.1 ) MFX_MPI = 1
      MFY_MPI = LFY_MPI+2
      IF( LFY_MPI.EQ.1 ) MFY_MPI = 1
      MFZ_MPI = LFZ_MPI+2
      IF( LFZ_MPI.EQ.1 ) MFZ_MPI = 1
      DO L = 1,LP_MPI
        NC_MPI(L) = 0
        NCN_MPI(L) = 0
        DO K = 1,(MFX_MPI*MFY_MPI*MFZ_MPI)
          ND_MPI(K,L) = 1
        ENDDO
      ENDDO
      GRAVX = 0.D+0
      GRAVY = 0.D+0
      GRAVZ = 9.81D+0
      THXZ = 0.D+0
      THYZ = 0.D+0
      IFLD = 1
      JFLD = 1
      KFLD = 1
      IJFLD = 1
      JKFLD = 1
      KIFLD = 1
      NBRN = 0
      NFLD = 1
      NFBN = 1
      NRFN = 0
      ICS = 1
      NXP = 0
      NROCK = 0
      NROCK2 = 0
      NDIM = 1
      N_DB = 0
      DO 748 L = 1,3
        XREF(L) = 0.D+0
        XREFU(L) = 'null'
        IXREF(L) = 0
  748 CONTINUE
      DO 762 L = 1,6
        MDIM(L) = 0
  762 CONTINUE
      DO LZ = 1,LFZ
        DO LY = 1,LFY
          DO LX = 1,LFX
            ND(LX,LY,LZ) = 0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LRC
        ROCK(L) = 'null'
        ROCK2(L) = 'null'
      ENDDO
      DO L = 1,3
        DO K = 1,3
          ROTMAT(K,L) = 0.D+0
          IF( K.EQ.L ) ROTMAT(K,L) = 1.D+0
        ENDDO
      ENDDO
!
!---  [ SPILL ] Spill variables  ---
!
      DO L = 1,LFX*LFY
        XSP(L) = 0.D+0
        YSP(L) = 0.D+0
        ZSP(L) = 0.D+0
        AFZSP(L) = 0.D+0
      ENDDO
      DO 772 L = 1,(LFX+1)*LFY
        TXSP(L) = 0.D+0
        DYSP(L) = 0.D+0
        PDXSP(L) = 0.D+0
  772 CONTINUE
      DO 773 L = 1,(LFY+1)*LFX
        TYSP(L) = 0.D+0
        DXSP(L) = 0.D+0
        PDYSP(L) = 0.D+0
  773 CONTINUE
      DO 776 K = 1,LSFV
        DO 774 L = 1,(LFY+1)*LFX
          VNSP(K,L) = 0.D+0
          VLSP(K,L) = 0.D+0
  774   CONTINUE
        DO 775 L = 1,(LFX+1)*LFY
          UNSP(K,L) = 0.D+0
          ULSP(K,L) = 0.D+0
  775   CONTINUE 
  776 CONTINUE
      DO 778 L = 1,LFX*LFY
        DO 777 K = 1,LSV
          RHONSP(K,L) = 0.D+0
          VISNSP(K,L) = 0.D+0
          HNSP(K,L) = 0.D+0
          RHOLSP(K,L) = 0.D+0
          VISLSP(K,L) = 0.D+0
          HLSP(K,L) = 0.D+0
          SRCOSP(K,L) = 0.D+0
          SRCWSP(K,L) = 0.D+0
  777   CONTINUE
        SPNORM(L) = 0.D+0
  778 CONTINUE
      DO 780 L = 1,LFX*LFY
        DO 779 K = 1,LUK
          DNRSP(K,L) = 0.D+0
  779   CONTINUE 
  780 CONTINUE
      DO 782 L = 1,LRC
        SPHMIN(L) = 0.D+0
  782 CONTINUE
!
!---  [ JACOB ] Linear system solution variables  ---
!
      DO 805 L = 1,LJE
        DO 800 K = 1,LJD
          ALU(K,L) = 0.D+0
  800   CONTINUE
  805 CONTINUE
      DO 810 L = 1,LJF
        BLU(L) = 0.D+0
  810 CONTINUE
      DO 815 L = 1,LJA
        CLU(L) = 0.D+0
  815 CONTINUE
      DO 825 L = 1,LJB
        DLU(L) = 0.D+0
  825 CONTINUE
      DO 835 L = 1,LFD
        DO 830 K = 1,LUKW
          RSDL(K,L) = 0.D+0
  830   CONTINUE
  835 CONTINUE
#ifdef lis
      DO L = 1,LJI
        ILU(L) = L
      ENDDO
#elif petsc
      DO L = 1,LJI
        ILU(L) = L - 1
      ENDDO
#else
      DO L = 1,LJI
        ILU(L) = 0
      ENDDO
#endif
      DO 855 L = 1,LJJ
        JLU(L) = 0
  855 CONTINUE
      DO 865 L = 1,LJH
        DO 860 K = 1,LJG
          KLU(K,L) = 0
  860   CONTINUE
  865 CONTINUE
      DO 866 K = 1,LJO
        MLU(K) = 0
  866 CONTINUE
      DO 867 K = 1,LJC
        NLU(K) = 0
  867 CONTINUE
      DO 869 L = 1,LJL
        DO 868 K = 1,LJK
          KLUC(K,L) = 0
  868   CONTINUE
  869 CONTINUE
      DO 870 K = 1,LJM
        MLUC(K) = 0
  870 CONTINUE
      DO 871 K = 1,LJN
        NLUC(K) = 0
  871 CONTINUE
      DO 873 L = 1,LFD
        DO 872 K = 1,LUKW
          IM(K,L) = 0
  872   CONTINUE
  873 CONTINUE
      DO 885 L = 1,LUK
        DO 880 K = 1,LUK
          IF( L.EQ.K ) THEN
            JM(L,K) = 1
          ELSEIF( K.GT.L ) THEN
            JM(L,K) = K-L+1
          ELSEIF( K.LT.L ) THEN
            JM(L,K) = K-L+LUK+1
          ENDIF
  880   CONTINUE
  885 CONTINUE
      ISVC = 0
      ISVT = 0
      ISVF = 0
      MUC = 0
      MLC = 0
      MDC = 0
      MUT = 0
      MLT = 0
      MDT = 0
      MKC = 0
      MKT = 0
      NNZFR = 0
      NNZTR = 0
!
!---  [ NAPL ] Volatile organic compound and surfactant
!     critical properties and coefficients  ---
!
      PCRO = 0.D+0
      TCRO = 0.D+0
      ZCRO = 0.D+0
      VCRO = 0.D+0
      TFPO = 0.D+0
      WTMO = 153.823D+0
      RCO = 0.D+0
      PAFO = 0.D+0
      TBO = 0.D+0
      TTPO = 0.D+0
      DPMO = 0.D+0
      PCCVO = 0.D+0
      PCCVA = 0.D+0
      WHBTO = 0.D+0
      WHBTA = 0.D+0
      RHORO = 0.D+0
      TDRO = 0.D+0
      ZRAO = 0.D+0
      DO 900 L = 1,4
        VISCO(L) = 0.D+0
        SATOC(L) = 0.D+0
        SATAC(L) = 0.D+0
        CPOC(L) = 0.D+0
        CPAC(L) = 0.D+0
        CIMTC(L) = 0.D+0
        SFCSF(L) = 0.D+0
        VISCS(L) = 0.D+0
  900 CONTINUE
      DO 901 L = 1,6
        SFCSF(L) = 0.D+0
  901 CONTINUE
      DO 902 L = 1,11
        TERDC(L) = 0.D+0
  902 CONTINUE
      CIMTC(1) = 5.5D-1
      CIMTC(2) = 2.5D-1
      CIMTC(3) = 1.5D+0
      CIMTC(4) = 1.5D+0
      VISRN = 0.D+0
      RHORN = 0.D+0
      HCOW = 0.D+0
      SUFO = 0.D+0
      RCS = 0.D+0
      TCRS = 0.D+0
      PCRS = 0.D+0
      VCRS = 0.D+0
      ZCRS = 0.D+0
      PAFS = 0.D+0
      RHORA = 0.D+0
      DPMA = 0.D+0
      TDRA = 0.D+0
      TFPA = 0.D+0
      ZRAA = 0.D+0
      SUFA = 0.D+0
      TCRA = 132.63D+0
      PCRA = 3.7858D+6
      TTPA = 0.D+0
      PTPA = 0.D+0
      DO 905 L = 1,LRCS
        PCSLS(1,L) = 0.D+0
        PCSLS(2,L) = 0.D+0
        IPCSLS(L) = 0
  905 CONTINUE
      IVISO = 0
      IVISS = 0
      IVAPO = 0
      IRHOO = 0
      IMTC = 1
      IVAPS = 0
      IRHOS = 0
      NQO = 0
      NQA = 0
      ITERDC = 0
!
!---  [ OUTPU ] Output control variables  ---
!
      CALL IN_OUTP
!
!---  [ PORMED ] Porous media properties  ---
!
      DO 1130 L = 1,LRC
        DO 1100 K = 1,6
          POR(K,L) = 0.D+0
 1100   CONTINUE
        DO 1102 K = 1,13
          GAMMA(K,L) = 1.D+0
 1102   CONTINUE
        DO 1104 K = 1,6
          TOR(K,L) = 0.D+0
 1104   CONTINUE
        CMP(1,L) = 1.D-7
        CMP(2,L) = 1.D-7
        CMP(3,L) = 0.D+0
        CMP(4,L) = 0.D+0
        CHML(L) = 0.D+0
        DO 1110 K = 1,9
          THKS(K,L) = 0.D+0
 1110   CONTINUE
        DO 1112 K = 1,5
          DFEF(K,L) = 0.D+0
 1112   CONTINUE
        RHOS(L) = 2650.D+0
        CPS(L) = 0.D+0
        DO 1116 K = 1,9
          PERM(K,L) = 0.D+0
 1116   CONTINUE
        DO 1120 K = 1,LSCHR
          SCHR(K,L) = 0.D+0
 1120   CONTINUE
        DO 1122 K = 1,LRPL
          RPLT(1,K,L) = 0.D+0
          RPLT(2,K,L) = 0.D+0
          RPLT(3,K,L) = 0.D+0
          RPLT(4,K,L) = 0.D+0
 1122   CONTINUE
        DO 1123 K = 1,LRPGC
          RPGC(K,L) = 0.D+0
 1123   CONTINUE
        DO 1124 K = 1,LRPLC
          RPLC(K,L) = 0.D+0
 1124   CONTINUE
        DO 1125 K = 1,LRPNC
          RPNC(K,L) = 0.D+0
 1125   CONTINUE
        DPLGA(L) = 0.D+0
        DPLLA(L) = 0.D+0
        DPTGA(L) = 0.D+0
        DPTLA(L) = 0.D+0
        SCALNM(L) = 'null'
        ISCALE(L) = 0
        HCMWE(L) = 0.D+0
        ITOR(L) = 0
        IRPG(L) = 0
        IRPL(L) = 0
        IRPLT(1,L) = 0
        IRPLT(2,L) = 0
        IRPLT(3,L) = 0
        IRPN(L) = 0
        ISCHR(L) = 0
        ISM(L) = 0
        ITHK(L) = 0
        IDP(L) = 0
        ISKP(L) = 0
        IPRF(L) = 0
        INHYD(L) = 0
        DO 1126 K = 1,2
          IRLTBLT(K,L,1) = 0
          IRLTBLT(K,L,2) = 0
          IRLTBLT(K,L,3) = 0
 1126   CONTINUE
        DO 1128 K = 1,4
          ISLTBL(K,L) = 0
          IRLTBL(K,L) = 0
          IRGTBL(K,L) = 0
          IRNTBL(K,L) = 0
 1128   CONTINUE
 1130 CONTINUE
      NSCALE = 0
      IGAMMA(1) = 0
      IGAMMA(2) = 0
      IGAMMA(3) = 0
      IGAMMA(4) = 0
      IGAMMA(5) = 0
      IGAMMA(6) = 0
      IGAMMA(7) = 0
      IGAMMA(8) = 0
!
!---  [ SOLTN ] Solution control variables  ---
!
      TM = 1.D+20
      TMMX = -1.D+20
      TMPR = 1.D+20
      DT = 0.D+0
      DTI = 1.D+0
      DTMX = 0.D+0
      DTMN = 1.D+20
      DTAF = 0.D+0
      DTCF = 2.D-1
      DTO = 0.D+0
      DTSO = 0.D+0
      RSDMX = 0.D+0
      RLXF = 1.D+0
      RLMSG = 0.D+0
      CPUMX = 0.D+0
      CLKMX = 0.D+0
      CPUSEC = 0.D+0
      CLKSEC = 0.D+0
      USER = 'null'
      CMPNY = 'null'
      TITLE = 'null'
      INPDAT = 'null'
      INPTIM = 'null'
      CARD = 'null'
      CHMSG = 'null'
      CH_VRSN = '1087'
      DO L = 1,LNOTES
       NOTES(L) = ' '
      ENDDO
      ICD = 1
      IVR = 1
      IVRSN = 1
      ISIC = 0
      ICNV = 3
      JSTOP = 0
      IEO = 0
      ILES = 1
      IOM = 0
      IUNM = 0
      IUNKG = 0
      IUNS = 0
      IUNK = 0
      IUNMOL = 0
      IMSG = 0
      IEQT = 0
      IEQW = 0
      IEQA = 0
      IEQDA = 0
      IEQDO = 0
      IEQO = 0
      IEQC = 0
      IEQS = 0
      IEQD = 0
      IEQALC = 0
      DO L = 1,LNGC
        IEQGC(L) = 0
      ENDDO
      DO L = 1,100
        ISLC(L) = 0
        VSLC(L) = 0.D+0
      ENDDO
      RTOL_PETSC = 1.D-5
      ATOL_PETSC = 1.D-50
      DTOL_PETSC = 1.D+5
      MAXITS_PETSC = 100000
!
!---  Upwind oil transport default  ---
!
      ISLC(8) = 3
!
!---  Supersaturation factor default  ---
!
      VSLC(1) = 1.5D+0
      ISLC(52) = 1
      CRNTMXC = 0.D+0
      DO L = 1,20
        IDMN(L) = 1
        WFMN(L) = 0.D+0
      ENDDO
      IDMN(2) = 4
      IDMN(3) = 4
      IDMN(4) = 4
      IDMN(8) = 4
      IDMN(9) = 4
      IDMN(10) = 4
      IDFLT = 0
      IDFLTD = 0
      MXSTEP = 0
      NRSD = 0
      NSTEP = 0
      NRST = 0
      NITER = 0
      NTSR = 0
      NGC = 0
      NEPD = 0
      MEPD = 0
      IEPD = 0
      DO L = 1,(LUK*(1+LWELL+LSPILL))
        RSD(L) = 0.D+0
        NSD(L) = 0
      ENDDO
      DO L = 1,LEPD
        TMPS(L) = 0.D+0
        TMPE(L) = 0.D+0
        TMPD(L) = 0.D+0
        TMPX(L) = 0.D+0
        TMPN(L) = 1.D+20
        TMPA(L) = 0.D+0
        TMPC(L) = 2.D-1
      ENDDO
      DO L = ISUB_LOG+1,32
        SUB_LOG(L) = 'null'
      ENDDO
      IPID = 0
      PID_CHAR = ' '
      PIDFNM = ' '
      IFXST = .FALSE.
!
!---  [ SOURC ] Source variables  ---
!
      DO J = 1,8+LSOLU
        DO K = 1,LSTM
          DO L = 1,LSR
            SRC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LSR
        DO K = 1,8
          SRCP(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          QNW(K,L) = 0.D+0
          QLW(K,L) = 0.D+0
          QTW(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          PLWB(K,L) = 0.D+0
          PGW(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDT
        SRCIT(L) = 0.D+0
        DO K = 1,LSV
          SRCT(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDD
        SRCID(L) = 0.D+0
        DO K = 1,LSV
          SRCD(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFD
        SRCIW(L) = 0.D+0
        DO K = 1,LSV
          SRCW(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDG
        SRCIA(L) = 0.D+0
        DO K = 1,LSV
          SRCA(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDGC
        DO K = 1,LNGC
          SRCIGC(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDGC
        DO K = 1,LSV
          DO J = 1,LNGC
            SRCGC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LFDNH
        SRCIO(L) = 0.D+0
        DO K = 1,LSV
          SRCO(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDS
        SRCIS(L) = 0.D+0
        DO K = 1,LSV
          SRCS(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSR
        DO K = 1,13
          ISRDM(K,L) = 0
        ENDDO
        ISRT(L) = 0
        ISRM(L) = 0
        NSOLSR(L) = 0
      ENDDO
      DO L = 1,LSR
        DO K = 1,LWSI
          IWSI(K,L) = 0
          GWSI(1,K,L) = 0.D+0
          GWSI(2,K,L) = 0.D+0
          GWSI(3,K,L) = 0.D+0
          SWSI(1,K,L) = 0.D+0
          SWSI(2,K,L) = 0.D+0
          SWSI(3,K,L) = 0.D+0
        ENDDO
      ENDDO
      NSR = 0
!
!---  [ TABL ] Tabular data variables  ---
!
      DO L = 1,LTBL
        TBLX(L) = 0.D+0
        TBLY(L) = 0.D+0
        TBLDDX(L) = 0.D+0
        TBLDDY(L) = 0.D+0
      ENDDO
      NTBL = 0
!
!---  [ TRNSPT ] Solute transport variables  ---
!
      DO L = 1,LSOLU
        DO K = 1,LFDC
          CNL(K,L) = 0.D+0
          CNLO(K,L) = 0.D+0
          ICTN(K,L) = 0
        ENDDO
        DO K = 1,LFDCR
          SRCIC(K,L) = 0.D+0
          YL(K,L) = 0.D+0
          YG(K,L) = 0.D+0
          YN(K,L) = 0.D+0
        ENDDO
        DO K = 1,LFDCR
          C(K,L) = 0.D+0
          CO(K,L) = 0.D+0
          ICT(K,L) = 0
        ENDDO
        DO K = 1,LSXC
          UC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSYC
          VC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSZC
          WC(K,L) = 0.D+0
        ENDDO
        HLF(L) = 1.D+20
        DO K = 1,LCHEM
          RHLF(L,K) = 1.D+20
          RHLFL(L,K) = 1.D+20
          RHLFN(L,K) = 1.D+20
        ENDDO
        SMDL(L) = 0.D+0
        SMDG(L) = 0.D+0
        SMDN(L) = 0.D+0
        DO K = 1,LRC
          SDCL(1,K,L) = 0.D+0
          SDCL(2,K,L) = 0.D+0
          SDCL(3,K,L) = 0.D+0
          PCSL(1,K,L) = 1.D-20
          PCSL(2,K,L) = 0.D+0
          PCSL(3,K,L) = 0.D+0
          PCSL(4,K,L) = 0.D+0
          PCSL(5,K,L) = 0.D+0
          PCSN(1,K,L) = 1.D-20
          PCSN(2,K,L) = 0.D+0
          PCSN(3,K,L) = 0.D+0
          PCSN(4,K,L) = 0.D+0
          PCSN(5,K,L) = 0.D+0
          IPCSL(K,L) = 0
          SMDEF(K,L) = 1.D+0
        ENDDO
        DO K = 1,5
          PCLN(K,L) = 1.D+20
          PCGL(K,L) = 1.D-20
          PCGN(K,L) = 1.D-20
        ENDDO
        DO K = 1,LBCC
          CB(K,L) = 0.D+0
          CBO(K,L) = 0.D+0
          YLB(K,L) = 0.D+0
          YGB(K,L) = 0.D+0
          YNB(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSOLU
          CHDF(K,L) = 0.D+0
          DO J = 1,LCHEM
            RCHDF(K,L,J) = 0.D+0
            RCHDFL(K,L,J) = 0.D+0
            RCHDFN(K,L,J) = 0.D+0
          ENDDO
        ENDDO
        SOLUT(L) = 'null'
        IEDL(L) = 1
        IPCL(L) = 0
        IPCN(L) = 0
        IPCGL(L) = 0
        IPCGN(L) = 0
        IPCLN(L) = 0
        IMTLN(L) = 0
        NCHEM(L) = 0
        CMTLN(1,L) = 0.D+0
        CMTLN(2,L) = 0.D+0
        CMTLN(3,L) = 0.D+0
        CMTLN(4,L) = 0.D+0
        N_CRN(L) = 0
        CCL_CRN(L) = 0.D+0
        SOLML(L) = 0.D+0
      ENDDO
      DO L = 1,LCDC+LSOLU
        IBCDS(L) = 0
      ENDDO
      DO L = 1,LCDC
        NBCDP(L) = 0
        DO K = 1,LCDP
          DO J = 1,LCDS
            IBCDP(J,K,L) = 0
          ENDDO
        ENDDO
      ENDDO
      NBCDS = 0
!
!---  [ TRNSPT ] Species transport variables  ---
!
      IF( LSPT.GT.1 ) THEN
        I = LSOLU+1
        J = LSOLU+LSPT
        DO L = I,J
          SOLUT(L) = 'null'
          DO K = 1,LFDC
            SRCIC(K,L) = 0.D+0
            YL(K,L) = 0.D+0
            YG(K,L) = 0.D+0
            YN(K,L) = 0.D+0
          ENDDO
          DO K = 1,LFDCR
            C(K,L) = 0.D+0
            CO(K,L) = 0.D+0
          ENDDO
          DO K = 1,LSXC
            UC(K,L) = 0.D+0
          ENDDO
          DO K = 1,LSYC
            VC(K,L) = 0.D+0
          ENDDO
          DO K = 1,LSZC
            WC(K,L) = 0.D+0
          ENDDO
          HLF(L) = 1.D+20
          SMDL(L) = 0.D+0
          SMDG(L) = 0.D+0
          SMDN(L) = 0.D+0
          DO K = 1,LRC
            SDCL(1,K,L) = 0.D+0
            SDCL(2,K,L) = 0.D+0
            SDCL(3,K,L) = 0.D+0
            SMDEF(K,L) = 1.D+0
          ENDDO
          DO K = 1,LBCC
            CB(K,L) = 0.D+0
            CBO(K,L) = 0.D+0
            YLB(K,L) = 0.D+0
            YGB(K,L) = 0.D+0
            YNB(K,L) = 0.D+0
          ENDDO
          IEDL(L) = 1
        ENDDO
      ENDIF
      ELC_SOL = 'null'
      ELC_DUN = 0.D+0
      ELC_VUN = 0.D+0
      DT_CRN = 0.D+0
      DTI_CRN = 1.D+0
      TM_CRN = 0.D+0
      IDF_ELC = 0
      IVF_ELC = 0
      NSL_ELC = 0
      DO L = 1,4
        ELC_DCF(L) = 0.D+0
        ELC_VCF(L) = 0.D+0
        ELC_SCF(L) = 0.D+0
      ENDDO
      SMDLS = 0.D+0
      DO L = 1,LRC
        SDCLS(1,L) = 0.D+0
        SDCLS(2,L) = 0.D+0
        SDCLS(3,L) = 0.D+0
        DISPL(L) = 0.D+0
        DISPT(L) = 0.D+0
      ENDDO
      DO K = 1,LRCN
        PCSLD(1,K) = 0.D+0
        PCSLD(2,K) = 0.D+0
        DPTD(K) = 0.D+0
        DPLD(K) = 0.D+0
        D50(K) = 0.D+0
        PMDD(K) = 1.D-4
        IPCSLD(K) = 0
      ENDDO
      DO L = 1,LRCS
        DPLGS(L) = 0.D+0
        DPTRS(L) = 0.D+0
      ENDDO
      DO L = 1,LFD
        CRNTL(L) = 0.D+0
        CRNTG(L) = 0.D+0
        CRNTN(L) = 0.D+0
      ENDDO
      CRNTMXT = 1.D+0
      NSOLU = 0
      IDISP = 0
      IEDLS = 0
      IDSPS = 0
      IDSPD = 0
      ICRNT = 0
!
!---  [ FDVP ] Primary field variables  ---
!
      DO L = 1,LFD
        DO K = 1,LSV
          T(K,L) = 20.D+0
          PG(K,L) = 0.D+0
          PN(K,L) = 0.D+0
          PSO(K,L) = 0.D+0
          SG(K,L) = 0.D+0
          SN(K,L) = 0.D+0
          PORD(K,L) = 0.D+0
          PORT(K,L) = 0.D+0
          PSW(K,L) = 0.D+0
          PVA(K,L) = 0.D+0
          PVO(K,L) = 0.D+0
          PVW(K,L) = 0.D+0
          XGA(K,L) = 0.D+0
          XGO(K,L) = 0.D+0
          XGW(K,L) = 0.D+0
          XMGA(K,L) = 0.D+0
          XMGO(K,L) = 0.D+0
          XMGW(K,L) = 0.D+0
          RHOG(K,L) = 0.D+0
          RHON(K,L) = 0.D+0
          RHOMG(K,L) = 0.D+0
          RHOMN(K,L) = 0.D+0
          SI(K,L) = 0.D+0
          PI(K,L) = 0.D+0
          POSM(K,L) = -1.D+0
          BTGL(K,L) = 1.D+0
          BTGN(K,L) = 1.D+0
          BTNL(K,L) = 1.D+0
          TMS(K,L) = 0.D+0
          PERMRF(K,L) = 1.D+0
        ENDDO
        PCMP(L) = 0.D+0
        TCMP(L) = 0.D+0
        DO K = 1,3
          POR0(K,L) = 0.D+0
        ENDDO
        DO K = 1,LUK
          DNR(K,L) = 0.D+0
        ENDDO
        SDPM(L) = 0.D+0
        SDPF(L) = 0.D+0
        DO K = 1,LSU
          TRPNL(K,L) = 0.D+0
          TRPGL(K,L) = 0.D+0
        ENDDO
        PG(1,L) = -2.D+0*EPSL
        PERMV(1,L) = 0.D+0
        PERMV(2,L) = 0.D+0
        PERMV(3,L) = 0.D+0
        SCHRV(1,L) = 0.D+0
      ENDDO
!
!---  [ FDVP ] Aqueous primary field variables  ---
!
      DO 2090 L = 1,LFDL
        DO 2082 K = 1,LSV
          PL(K,L) = 0.D+0
          SL(K,L) = 0.D+0
          XLA(K,L) = 0.D+0
          XLO(K,L) = 0.D+0
          XLW(K,L) = 0.D+0
          XMLA(K,L) = 0.D+0
          XMLO(K,L) = 0.D+0
          XMLW(K,L) = 0.D+0
          RHOL(K,L) = 0.D+0
          VISL(K,L) = 0.D+0
          TORL(K,L) = 0.D+0
          RKL(1,K,L) = 0.D+0
          RKL(2,K,L) = 0.D+0
          RKL(3,K,L) = 0.D+0
          DFLO(K,L) = 0.D+0
          DFLA(K,L) = 0.D+0
          DFLS(K,L) = 0.D+0
          RHOML(K,L) = 0.D+0
 2082   CONTINUE
        PL(1,L) = -2.D+0*EPSL
        SL(1,L) = -2.D+0*EPSL
 2090 CONTINUE
!
!---  Special initialization for STOMP-WNE  ---
!
      IF( ICODE.EQ.30 ) THEN
        DO 2112 L = 1,LFD
          PVA(2,L) = -EPSL
 2112   CONTINUE
      ENDIF
!
!---  Special initialization for STOMP-OS  ---
!
      IF( ICODE.EQ.50 .OR. ICODE.EQ.51 ) THEN
        DO 2114 L = 1,LFD
          PG(2,L) = -EPSL
          SN(2,L) = -EPSL
 2114   CONTINUE
      ENDIF
!
!---  [ FDVD ] Dissolved-oil field variables  ---
!
      DO 2116 L = 1,LFDD
        AVPVG(L) = 0.D+0
        AVPVL(L) = 0.D+0
 2116 CONTINUE
!
!---  [ FDVT ] Energy field variables  ---
!
      DO 2120 L = 1,LFDT
        DO 2115 K = 1,LSV
          THKL(K,L) = 0.D+0
          THKG(K,L) = 0.D+0
          THKN(K,L) = 0.D+0
          HL(K,L) = 0.D+0
          HG(K,L) = 0.D+0
          HN(K,L) = 0.D+0
          HGA(K,L) = 0.D+0
          HGD(K,L) = 0.D+0
          HGO(K,L) = 0.D+0
          HGW(K,L) = 0.D+0
          HLA(K,L) = 0.D+0
          HLO(K,L) = 0.D+0
          HLS(K,L) = 0.D+0
          HLW(K,L) = 0.D+0
          UEG(K,L) = 0.D+0
          UEL(K,L) = 0.D+0
          UEN(K,L) = 0.D+0
 2115   CONTINUE
 2120 CONTINUE
!
!---  [ FDVG ] Gas field variables  ---
!
      DO 2130 L = 1,LFDG
        DO 2125 K = 1,LSV
          GNIFT(K,L) = 0.D+0
          VISG(K,L) = 0.D+0
          VISDG(K,L) = 0.D+0
          TORG(K,L) = 0.D+0
          RKG(K,L) = 0.D+0
          DFGW(K,L) = 0.D+0
          DFGO(K,L) = 0.D+0
          DFGA(K,L) = 0.D+0
          TMBP_A(K,L) = 0.D+0
          TMBP_O(K,L) = 0.D+0
 2125   CONTINUE
        IF( ICODE.EQ.34 ) THEN
          ICAIR(L) = 0
          ICCO2(L) = 0
        ELSE
          ICAIR(L) = 1
          ICCO2(L) = 1
        ENDIF
 2130 CONTINUE
!
!---  [ FDVN ] NAPL field variables  ---
!
      DO 2134 L = 1,LFD
        DO 2132 K = 1,LSV
          VISN(K,L) = 0.D+0
          TORN(K,L) = 0.D+0
          RKN(K,L) = 0.D+0
          XSO(K,L) = 0.D+0
 2132   CONTINUE
 2134 CONTINUE
      DO 2140 L = 1,LFDN
        ANLT(L) = 0.D+0
        ANLF(L) = 0.D+0
        HKL(L) = 0.D+0
        HKNT(L) = 0.D+0
        HKNF(L) = 0.D+0
 2140 CONTINUE
!
!---  [ FDVI ] Ice field variables  ---
!
      DO 2150 L = 1,LFDI
        DO 2145 K = 1,LSV
          RHOI(K,L) = 0.D+0
          THKI(K,L) = 0.D+0
          HI(K,L) = 0.D+0
 2145   CONTINUE
 2150 CONTINUE
!
!---  [ FDVS ] Salt/Surfactant field variables  ---
!
      DO 2160 L = 1,LFDS
        DO 2155 K = 1,LSV
          XLS(K,L) = 0.D+0
          OEC(K,L) = 0.D+0
          YLS(K,L) = 0.D+0
          RHOSP(K,L) = 0.D+0
          HSP(K,L) = 0.D+0
 2155   CONTINUE
        ICBRN(L) = 0
 2160 CONTINUE
      ISALT = 1
!
!---  [ FDVA ] Nonaqueous liquid field variables  ---
!
      DO L = 1,LFDA
        DO K = 1,LSV
          XNA(K,L) = 0.D+0
          XNO(K,L) = 0.D+0
          XNW(K,L) = 0.D+0
          XMNA(K,L) = 0.D+0
          XMNO(K,L) = 0.D+0
          XMNW(K,L) = 0.D+0
          SUGL(K,L) = 0.D+0
          SUGN(K,L) = 0.D+0
          SUNL(K,L) = 0.D+0
          DFNA(K,L) = 0.D+0
          DFNO(K,L) = 0.D+0
          DFNW(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDN2
        DO K = 1,LSV
          XNN(K,L) = 0.D+0
          XMNN(K,L) = 0.D+0
          DFNN(K,L) = 0.D+0
          DFLN(K,L) = 0.D+0
          SRCN(K,L) = 0.D+0
          TMBP_N(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  [ FDVGC ] Gas component field variables  ---
!
      IFK = 0
      DO 2167 M = 1,LFDGC
        DO 2164 L = 1,LNGC
          FK(L,M) = 0.D+0
 2164   CONTINUE
        DO L = 1,6
          BETA(L,M) = 0.D+0
        ENDDO
        DO 2166 L = 1,LSV
          DO 2165 K = 1,LNGC
            PVC(K,L,M) = 0.D+0
            XGC(K,L,M) = 0.D+0
            XMGC(K,L,M) = 0.D+0
            XNC(K,L,M) = 0.D+0
            XMNC(K,L,M) = 0.D+0
            XLC(K,L,M) = 0.D+0
            XMLC(K,L,M) = 0.D+0
            DFLC(K,L,M) = 0.D+0
            DFGC(K,L,M) = 0.D+0
            DFNC(K,L,M) = 0.D+0
            HGC(K,L,M) = 0.D+0
            TMC(K,L,M) = 0.D+0
            ZMC(K,L,M) = 0.D+0
 2165     CONTINUE
          ZG(L,M) = 0.D+0
          ZN(L,M) = 0.D+0
 2166   CONTINUE
        IZMC(M) = 0
        IBETA(M) = 0
 2167 CONTINUE
!
!---  [ FDVH ] Hydrate field variables  ---
!
      DO L = 1,LNHC
      DO K = 1,30
        HCPP(K,L) = 0.D+0
      ENDDO
      ENDDO
      DO 2172 L = 1,LFDH
        DO 2170 K = 1,LSV
          XHO(K,L) = 0.D+0
          XHW(K,L) = 0.D+0
          XHA(K,L) = 0.D+0
          RHOH(K,L) = 0.D+0
          THKH(K,L) = 0.D+0
          HH(K,L) = 0.D+0
          SH(K,L) = 0.D+0
          PH(K,L) = 0.D+0
          UEGA(K,L) = 0.D+0
          YMGA(K,L) = 0.D+0
          YMHGA(K,L) = 0.D+0
          YMGO(K,L) = 0.D+0
          YMHGO(K,L) = 0.D+0
          PVHA(K,L) = 0.D+0
          PVHO(K,L) = 0.D+0
          TMHA(K,L) = 0.D+0
          TMHO(K,L) = 0.D+0
          ZMCA(K,L) = 0.D+0
          ZMCO(K,L) = 0.D+0
 2170   CONTINUE
        PPEL(L) = 0.D+0
        PPEU(L) = 0.D+0
        TCR(L) = 0.D+0
        TCT(L) = 0.D+0
 2172 CONTINUE
      DO 2175 L = 1,LFD
        DO 2174 K = 1,15
          IC_OPT(K,L) = 0
 2174   CONTINUE
 2175 CONTINUE
      DO 2178 L = 1,LFDN2
        DO 2176 K = 1,LSV
          XHN(K,L) = 0.D+0
          PVHN(K,L) = 0.D+0
          PVN(K,L) = 0.D+0
          XGN(K,L) = 0.D+0
          XMGN(K,L) = 0.D+0
          XLN(K,L) = 0.D+0
          XMLN(K,L) = 0.D+0
          TMHN(K,L) = 0.D+0
          ZMCN(K,L) = 0.D+0
          YMHGN(K,L) = 0.D+0
          YMGN(K,L) = 0.D+0
          DFGN(K,L) = 0.D+0
 2176   CONTINUE
 2178 CONTINUE
      CHKN(1) = 1.D+0
      CHKN(2) = 1.D-2
      CHKN(3) = 1.D-2
      CHKN(4) = 1.D-4
      CHKN(5) = 2.D+0
      CHKN(6) = 1.D+0
      CHKN(7) = 1.D+0
      CHKN(8) = 1.D+0
      CHKN(9) = 1.D+0
      CHKN(10) = 1.D+0
!
!---  [ HYST ] Hysteretic k-s-P function variables ---
!
      DO 2220 L = 1,LFD
        DO 2210 K = 1,LSV
          SGT(K,L) = 0.D+0
          SNT(K,L) = 0.D+0
          SNR(K,L) = 0.D+0
 2210   CONTINUE
        SGTL(L) = 0.D+0
        SGTN(L) = 0.D+0
        DO 2212 K = 1,LSU
          ASLMIN(K,L) = 1.D+0
          ASTMIN(K,L) = 1.D+0
          ASTMAX(K,L) = 0.D+0
          NPHAZ(K,L) = 2
          IPH(K,L) = -1
 2212   CONTINUE
        ASNT(L) = 0.D+0
        ASNR(L) = 0.D+0
        ASGT(L) = 0.D+0
        ASGTL(L) = 0.D+0
        ASGTN(L) = 0.D+0
        ASL(L) = 0.D+0
        AST(L) = 0.D+0
 2220 CONTINUE
      BGN = 0.D+0
      BGL = 0.D+0
      BHL = 0.D+0
      BIL = 0.D+0
      BNL = 0.D+0
      CA_GN = 1.D+0
      CA_GL = 1.D+0
      CA_NL = 1.D+0
      SIG_GN = 0.D+0
      SIG_GL = 0.D+0
      SIG_NL = 0.D+0
      SIG_HL = 0.D+0
      SIG_IL = 0.D+0
      INSR = 0
!
!---  [ WELLS ] Well variables ---
!
      CALL IN_WELL
!
!---  [ COUP_WELL ] Coupled-well variables ---
!
      CALL IN_COUP_WELL
!
!---  [ LEAK_WELL ] Leaky-well variables ---
!
      CALL IN_LEAK_WELL
!
!---  [ HYDT ] Ternary hydrate variables ---
!
      CALL IN_HYDT
!
!---  [ EOR ] Enhanced oil recovery variables ---
!
      CALL IN_EOR
!
!---  [ PLT_ATM ] Surface cover, plant, and atmospheric variables  ---
!
      CALL IN_PLT_ATM
!
!---  [ UCODE ] Inverse (UCode) variables ---
!
      DO 2304 M = 1,LOBDT
        DO 2302 L = 1,LOBDS
          DO 2300 K = 1,2
            R_OBDS(K,L,M) = 0.D+0
 2300     CONTINUE
 2302   CONTINUE
 2304 CONTINUE
      DO 2340 L = 1,LOBDT
        DO 2320 K = 1,6
          R_OBDT(K,L) = 0
 2320   CONTINUE
        DO 2330 K = 1,9
          I_OBDT(K,L) = 0
 2330   CONTINUE
        C_OBDT(L) = 'null'
        NOBDS(L) = 0
 2340 CONTINUE
      NOBDT = 0
      NOBDP = 0
      IOBDEF = 40
      IOBDSF = 41
      IOBDUF = 42
      TMOB = 1.D+20
      FLG_UNI = .FALSE.
      FLG_EXT = .FALSE.
!
!---  [ NCG_PT ] Noncondensible gas property tables  ---
!
      DO INCG = 1,LNNGC
        IP_TA(INCG) = 0
        I_LV(INCG) = 0
        IPTP(INCG) = 0
        IPCR(INCG) = 0
        DO IPX = 1,LP_TA
          P_TA(IPX,INCG) = 0.D+0
          IT_TA(IPX,INCG) = 0
          IV_TA(IPX,INCG) = 0
          DO ITX = 1,LT_TA
            T_TA(ITX,IPX,INCG) = 0.D+0
            RHO_TA(ITX,IPX,INCG) = 0.D+0
            H_TA(ITX,IPX,INCG) = 0.D+0
            U_TA(ITX,IPX,INCG) = 0.D+0
            FUG_TA(ITX,IPX,INCG) = 0.D+0
            S_TA(ITX,IPX,INCG) = 0.D+0
            RHO_ST(ITX,IPX,INCG) = 0.D+0
            H_ST(ITX,IPX,INCG) = 0.D+0
            U_ST(ITX,IPX,INCG) = 0.D+0
            FUG_ST(ITX,IPX,INCG) = 0.D+0
            S_ST(ITX,IPX,INCG) = 0.D+0
          ENDDO
        ENDDO
        DO IPX = 1,L_LV
          T_LV(IPX,INCG) = 0.D+0
          P_LV(IPX,INCG) = 0.D+0
          RHOL_LV(IPX,INCG) = 0.D+0
          HL_LV(IPX,INCG) = 0.D+0
          UL_LV(IPX,INCG) = 0.D+0
          RHOV_LV(IPX,INCG) = 0.D+0
          HV_LV(IPX,INCG) = 0.D+0
          UV_LV(IPX,INCG) = 0.D+0
          FUG_LV(IPX,INCG) = 0.D+0
          SL_LV(IPX,INCG) = 0.D+0
          SV_LV(IPX,INCG) = 0.D+0
        ENDDO
      ENDDO
      INCG = 1
      DO ITX = 1,LT_PH
        DO IOX = 1,LO_PH
          P_PH(ITX,IOX) = 0.D+0
          YMHO_PH(ITX,IOX) = 0.D+0
          XSCA_PH(ITX,IOX) = 0.D+0
          XLCA_PH(ITX,IOX) = 0.D+0
          XSCO_PH(ITX,IOX) = 0.D+0
          XLCO_PH(ITX,IOX) = 0.D+0
        ENDDO    
      ENDDO
      DO ITX = 1,LT_PH
        T_PH(ITX) = 0.D+0
      ENDDO
      DO IOX = 1,LO_PH
        YMGO_PH(IOX) = 0.D+0
      ENDDO
      IT_TH = 0
      IS_TH = 0
      IO_TH = 0
      I_INH = 0
      DO ITX = 1,LT_TH
        DO IOX = 1,LO_TH
          T_TH(ITX,IOX) = 0.D+0
          YMHO_TH(ITX,IOX) = 0.D+0
          XSCA_TH(ITX,IOX) = 0.D+0
          XLCA_TH(ITX,IOX) = 0.D+0
          XSCO_TH(ITX,IOX) = 0.D+0
          XLCO_TH(ITX,IOX) = 0.D+0
        ENDDO
      ENDDO
      DO ITX = 1,LT_TH
        P_TH(ITX) = 0.D+0
      ENDDO
      DO IOX = 1,LO_TH
        YMGO_TH(IOX) = 0.D+0
      ENDDO
      DO L = 1,LFDG
        DO K = 1,LNNGC
          IC_NCG(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LINH
        INHNM(L) = ''
        XLIMX(L) = 0.D+0
        DO K = 1,6
          CINH(K,L) = 0.D+0
        ENDDO
      ENDDO
#ifdef sac
!
!---  [ BALAN ] global mass/activity balance reporting variables
!
      INMBLT = 0
      TMBL = 1.D+20
      DO 2600 L = 1,LBAL
        DBLTIM(L) = 0.D+0
 2600 CONTINUE
!
!---  [ SACRM ] SAC remediation function variables
!
      INMREM = 0
      TMRM = 1.D+20
      DO 2700 L = 1,LREM
        DRMTIM(L) = 0.D+0
        IRMSOL(L) = 0
        DRMDEP(L) = 0.D+0
        DRMEFF(L) = 0.D+0
        CRMSIT(L) = ''
 2700 CONTINUE
!
!---  [ SACRL ] SAC release flux plane location variables
!
      INMREL = 0
      TMRL = 1.D+20
      DO 2800 L = 1,LREL
        DRLTIM(L) = 0.D+0
        IRLK(L) = 0
 2800 CONTINUE
#endif
!
!---  [ REACT ] Reaction Variables  ---
!
      DO 2904 L = 1,LEQC
        DO 2900 K = 1,LSEC
          EQ_C(K,L) = 0.D+0
 2900   CONTINUE
        DO 2902 K = 1,LSEC+1
          IEQ_C(K,L) = 0
 2902   CONTINUE
 2904 CONTINUE
      DO 2906 L = 1,6
        ACTV16(L) = 0.D+0
 2906 CONTINUE
      DO 2914 L = 1,LEQE
        DO 2910 K = 1,LSEE
          EQ_E(K,L) = 0.D+0
 2910   CONTINUE
        DO 2912 K = 1,LSEE+1
          IEQ_E(K,L) = 0
 2912   CONTINUE
 2914 CONTINUE
      DO 2924 L = 1,LEQK
        DO 2920 K = 1,LSEK+LREK
          EQ_K(K,L) = 0.D+0
 2920   CONTINUE
        DO 2922 K = 1,LSEK+LREK+2
          IEQ_K(K,L) = 0
 2922   CONTINUE
 2924 CONTINUE
      DO 2934 L = 1,LRCE
        DO 2930 K = 1,5
          RC_E(K,L) = 0.D+0
 2930   CONTINUE
        RCNME(L) = ' '
 2934 CONTINUE
      DO 2944 L = 1,LRCK
        DO 2940 J = 1,LCKN
          DO 2938 K = 1,LSPK+11
            RC_K(K,J,L) = 0.D+0
            IRCKN(K) = 0
 2938     CONTINUE
 2940   CONTINUE
        DO 2942 K = 1,LSPK+3
          IRC_K(K,L) = 0
 2942   CONTINUE
        RCNMK(L) = ' '
        IRCKT(L) = 0
 2944 CONTINUE
      DO 2946 L = 1,LSPR
        IEQ_S(L) = 0
        CHARG(L) = 0.D+0
        FACTV(L) = 0.D+0
        ACTVS(L) = 0.D+0
 2946 CONTINUE
      DO 2948 L = 1,LEQE+LEQC+LEQK
        ISP_S(L) = 0
 2948 CONTINUE
      DO 2952 L = 1,LSPL
        SPNML(L) = ' '
        SP_L(1,L) = 0.D+0
        SP_L(2,L) = 0.D+0
        SP_L(3,L) = 0.D+0
 2952 CONTINUE
      DO 2956 L = 1,(LEQC+LEQK)
      DO 2954 K = 1,LFDRL
        YSPL(K,L) = 0.D+0
 2954 CONTINUE
 2956 CONTINUE
      SP_MDL = 0.D+0
      DO 2960 L = 1,LSPS
        SPNMS(L) = ' '
        SP_S(1,L) = 0.D+0
        SP_S(2,L) = 0.D+0
 2960 CONTINUE
      DO 2962 L = 1,LSPE
        SPNME(L) = ' '
        ISP_E(L) = 0
        IEL_LK(L) = 0
 2962 CONTINUE
      DO L = 1,LFDR
      DO K = 1,LSPS
        RS_S(1,K,L) = 0.D+0
        RS_S(2,K,L) = 0.D+0
        RS_S(3,K,L) = 0.D+0
        ISP_OW(K,L) = 0
      ENDDO
      ENDDO
      DO 2966 L = 1,LSOLU+LSPT
        IMMB(L) = 0
 2966 CONTINUE
      DO 2972 L = 1,LSPG
        SPNMG(L) = ' '
 2972 CONTINUE
      DO 2976 L = 1,(LEQC+LEQK)
      DO 2974 K = 1,LFDRG
        YSPG(K,L) = 0.D+0
 2974 CONTINUE
 2976 CONTINUE
      SP_MDG = 0.D+0
      DO 2982 L = 1,LSPN
        SPNMN(L) = ' '
 2982 CONTINUE
      DO 2986 L = 1,(LEQC+LEQK)
      DO 2984 K = 1,LFDRN
        YSPN(K,L) = 0.D+0
 2984 CONTINUE
 2986 CONTINUE
      SP_MDN = 0.D+0
      DO 2988 L = 1,LEQC
        SPNMC(L) = 'null'
 2988 CONTINUE
      DO 2990 L = 1,LEQK
        SPNMK(L) = 'null'
 2990 CONTINUE
      DO 2992 L=1,LSPE
        IEL_LK(L) = 0
 2992 CONTINUE
      DO 2994 L = 1,LSPR+14
        ISPLK(L) = 0
 2994 CONTINUE
      ACTVC = 1.D+0
      IACTV = 0
      IACTEX = 0
      NEQC = 0
      NEQE = 0
      NEQK = 0
      NRCE = 0
      NRCK = 0
      NSPC = 0
      NSPG = 0
      NSPK = 0
      NSPL = 0
      NSPLK = 0
      NSPN = 0
      NSPS = 0
      NSPE = 0
      NESITE = 0
      NRTSI = 0
      N_RST = 10
      ECKE_ER = .FALSE.
      DO 2997 L = 1,LSPR
        DO 2995 K = 1,LFDR
          SP_C(K,L) = 0.D+0
          SP_CO(K,L) = 0.D+0
          SP_CI(K,L) = 0.D+0
          IC_SP(K,L) = 0
 2995   CONTINUE
        DO 2996 K = 1,LBCC
          SP_CBO(K,L) = 0.D+0
 2996   CONTINUE
        ISP_MN(L) = 0
 2997 CONTINUE
      DO 2999 K = 1,LFDR
        C_PH(K) = 0.D+0
        DO 2998 L = 1,LSPS
          SP_CMN(K,L) = 0.D+0
          SP_RATE(K,L) = 0.D+0
          SP_AREA(K,L) = 0.D+0
 2998   CONTINUE
 2999 CONTINUE
      DO 3000 L = 1,LMC
        CFMX(L) = 1.D+0
 3000 CONTINUE
      DT_RST = 0.D+0
      DTI_RST = 0.D+0
      TM_RST = 0.D+0
!
!---  [ GEO_MECH ] Geomechanical Variables  ---
!
      CALL IN_GEO_MECH
!
!---  [ FDVP_BH ] Borehole Field Variables  ---
!
      CALL IN_FDVP_BH
!
!---  [ FDVP_FRC ] Fracture/Fault Field Variables  ---
!
      CALL IN_FDVP_FRC
!
!---  [ FDVT_FRC ] Fracture/Fault Energy Variables  ---
!
      CALL IN_FDVT_FRC
!
!---  [ FDVG_FRC ] Fracture/Fault Gas Variables  ---
!
      CALL IN_FDVG_FRC
!
!---  [ FDVN_FRC ] Fracture/Fault NAPL Variables  ---
!
      CALL IN_FDVN_FRC
!
!---  [ FDVS_FRC ] Fracture/Fault Salt Variables  ---
!
      CALL IN_FDVS_FRC
!
!---  [ FDVGC_FRC ] Fracture/Fault Gas Component Variables  ---
!
      CALL IN_FDVGC_FRC
!
!---  [ FLUX_BH ] Borehole Flux Variables  ---
!
      CALL IN_FLUX_BH
!
!---  [ FLUX_FRC ] Fracture/Fault Flux Variables  ---
!
      CALL IN_FLUX_FRC
!
!---  [ FLUXC_FRC ] Component Fracture/Fault Flux Variables  ---
!
      CALL IN_FLUXC_FRC
!
!---  [ GEOM_BH ] Borehole Geometry Variables  ---
!
      CALL IN_GEOM_BH
!
!---  [ GEOM_FRC ] Fracture/Fault Geometry Variables  ---
!
      CALL IN_GEOM_FRC
!
!---  [ PARM_BH ] Borehole Parameter Variables  ---
!
      CALL IN_PARM_BH
!
!---  [ PARM_FRC ] Fracture/Fault Parameter Variables  ---
!
      CALL IN_PARM_FRC
!
!---  [ TRNS_BH ] Borehole Transport Variables  ---
!
      CALL IN_TRNS_BH
!
!---  [ TRNS_FRC ] Fracture/Fault Transport Variables  ---
!
      CALL IN_TRNS_FRC
!
!---  [ POINTE ] Numerical scheme pointer arrays  ---
!
!      MNOD /2,3,2,4,2,5,2,6,2,7,2,8,2,9,2/
!      MADJ /2,2,3,2,4,2,5,2,6,2,7,2,8,2,9/
!      MPOS /2,3,2,4,2,5,2,6,2,7,2,8,2,9,2/
!      MNEG /2,2,3,2,4,2,5,2,6,2,7,2,8,2,9/
!      MFLX /1,3,2,5,4,7,6,9,8,11,10,13,12,15,14/
!      MPOSB /1,2,4,6,8,10,12,14/
!      MNEGB /1,3,5,7,9,11,13,15/
!
      DO 3400 L = 1,LSFV
        IF( MOD(L,2).EQ.1 ) THEN
          MNOD(L) = 2
        ELSE
          MNOD(L) = L/2 + 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MADJ(L) = (L-1)/2 + 2
        ELSE
          MADJ(L) = 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MPOS(L) = 2
        ELSE
          MPOS(L) = L/2 + 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MNEG(L) = (L-1)/2 + 2
        ELSE
          MNEG(L) = 2
        ENDIF
        IF( MOD(L,2).EQ.1 ) THEN
          MFLX(L) = MAX( 1,(L-1) )
        ELSE
          MFLX(L) = L + 1
        ENDIF
 3400 CONTINUE
      DO 3410 L = 1,LSV
        MPOSB(L) = MAX( 1,((2*L)-2) )
        MNEGB(L) = (2*L) - 1
 3410 CONTINUE
#ifdef dual_porosity
!
!---  [ DUAL_POR ] Dual porosity variables  ---
!
      DO L = 1,LRC
        DO K = 1,17
          FRAC_P = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFD_EC
        DO K = 1,LSV
          SG_F(K,L) = 0.D+0
          SL_F(K,L) = 0.D+0
          SN_F(K,L) = 0.D+0
          SG_M(K,L) = 0.D+0
          SL_M(K,L) = 0.D+0
          SN_M(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFD_DP
        VOL_M(L) = 0.D+0
        DO K = 1,LUK
          DNR_M(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSFV
          DO J = 1,LNGC
            QGC_FM(J,K,L) = 0.D+0
            QNC_FM(J,K,L) = 0.D+0
          ENDDO
          QLA_FM(K,L) = 0.D+0
          QGW_FM(K,L) = 0.D+0
          QLW_FM(K,L) = 0.D+0
          QS_FM(K,L) = 0.D+0
          QQ_FM(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          DO J = 1,LNGC
            DFNC_M(J,K,L) = 0.D+0
            TMC_M(J,K,L) = 0.D+0
            XGC_M(J,K,L) = 0.D+0
            XMGC_M(J,K,L) = 0.D+0
            XNC_M(J,K,L) = 0.D+0
            XMNC_M(J,K,L) = 0.D+0
            ZMC_M(J,K,L) = 0.D+0
          ENDDO
          DFGW_M(K,L) = 0.D+0
          DFLA_M(K,L) = 0.D+0
          DFLS_M(K,L) = 0.D+0
          HG_M(K,L) = 0.D+0
          HGA_M(K,L) = 0.D+0
          HGW_M(K,L) = 0.D+0
          HL_M(K,L) = 0.D+0
          HLW_M(K,L) = 0.D+0
          HN_M(K,L) = 0.D+0
          HSP_M(K,L) = 0.D+0
          UEG_M(K,L) = 0.D+0
          UEL_M(K,L) = 0.D+0
          UEN_M(K,L) = 0.D+0
          THKG_M(K,L) = 0.D+0
          THKL_M(K,L) = 0.D+0
          THKN_M(K,L) = 0.D+0
          PERMRF_M(K,L) = 1.D+0
          PG_M(K,L) = 0.D+0
          PL_M(K,L) = 0.D+0
          PN_M(K,L) = 0.D+0
          POSM_M(K,L) = 0.D+0
          PSO_M(K,L) = 0.D+0
          PORD_M(K,L) = 0.D+0
          PORT_M(K,L) = 0.D+0
          PVA_M(K,L) = 0.D+0
          PVW_M(K,L) = 0.D+0
          RHOG_M(K,L) = 0.D+0
          RHOL_M(K,L) = 0.D+0
          RHOMG_M(K,L) = 0.D+0
          RHOML_M(K,L) = 0.D+0
          RHOMN_M(K,L) = 0.D+0
          RHON_M(K,L) = 0.D+0
          RHOSP_M(K,L) = 0.D+0
          RKG_M(K,L) = 0.D+0
          RKL_M(1,K,L) = 0.D+0
          RKL_M(2,K,L) = 0.D+0
          RKL_M(3,K,L) = 0.D+0
          RKN_M(K,L) = 0.D+0
          SS_M(K,L) = 0.D+0
          T_M(K,L) = 0.D+0
          TMS_M(K,L) = 0.D+0
          TORG_M(K,L) = 0.D+0
          TORL_M(K,L) = 0.D+0
          TORN_M(K,L) = 0.D+0
          VISG_M(K,L) = 0.D+0
          VISL_M(K,L) = 0.D+0
          VISN_M(K,L) = 0.D+0
          XGW_M(K,L) = 0.D+0
          XLA_M(K,L) = 0.D+0
          XLS_M(K,L) = 0.D+0
          XLW_M(K,L) = 0.D+0
          XMGW_M(K,L) = 0.D+0
          XMLA_M(K,L) = 0.D+0
          XMLS_M(K,L) = 0.D+0
          XMLW_M(K,L) = 0.D+0
          YLS_M(K,L) = 0.D+0
          ZG_M(K,L) = 0.D+0
          ZN_M(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          POR0(K,L) = 0.D+0
        ENDDO
        DO K = 1,LUK
          IM_M(K,L) = 0
        ENDDO
        DO K = 1,LSU
          NPHAZ_M(K,L) = 0
        ENDDO
      ENDDO
#endif
      ISUB_LOG = ISUB_LOG-1
!
!---  End of INTLZ group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_BOUN
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
!     Initialize boundary condition variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 June 2001.
!     Last Modified by MD White, PNNL, 19 June 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE BCVT
      USE BCVS
      USE BCVP
      USE BCVN
      USE BCVI
      USE BCVH
      USE BCVGC
      USE BCVG
      USE BCVA
      USE BCV
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_BOUN'
!
!---  [ BCV ] Global boundary condition variables  ---
!
      DO L = 1,LBCIN
        TSBC(L) = 0.D+0
        XSBC(L) = 0.D+0
        YSBC(L) = 0.D+0
        ZSBC(L) = 0.D+0
        DO K = 1,LBTM
          BCXYZG(K,1) = 0.D+0
          BCXYZG(K,2) = 0.D+0
          BCXYZG(K,3) = 0.D+0
          DO J = 1,LBCV
            BC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,LBC
        PHDL(1,L) = 0.D+0
        PHDL(2,L) = 0.D+0
        PHDN(1,L) = 0.D+0
        PHDN(2,L) = 0.D+0
        XPBC(L) = 0.D+0
        YPBC(L) = 0.D+0
        ZPBC(L) = 0.D+0
        IBCC(L) = 0
        IBCD(L) = 0
        IBCM(L) = 0
        IBCN(L) = 0
        IBCIN(L) = 0
        DO K = 1,LUK+LPH*LSOLU*LC+LR+LL+LG+LN
          IBCT(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LBC
        DO K = 1,LPH
          NBHG(K,L) = 0
        ENDDO
      ENDDO
      NBC = 0
!
!---  [ BCVP ] Primary boundary condition variables  ---
!
      DO L = 1,LBC
        DO K = 1,LSV
          TB(K,L) = 0.D+0
          PLB(K,L) = 0.D+0
          PGB(K,L) = 0.D+0
          PNB(K,L) = 0.D+0
          PSOB(K,L) = 0.D+0
          PSWB(K,L) = 0.D+0
          SLB(K,L) = 0.D+0
          SGB(K,L) = 0.D+0
          SNB(K,L) = 0.D+0
          PORDB(K,L) = 0.D+0
          PORTB(K,L) = 0.D+0
          RHOMGB(K,L) = 0.D+0
          RHOMLB(K,L) = 0.D+0
          RHOMNB(K,L) = 0.D+0
          PVAB(K,L) = 0.D+0
          PVOB(K,L) = 0.D+0
          PVWB(K,L) = 0.D+0
          XGAB(K,L) = 0.D+0
          XGOB(K,L) = 0.D+0
          XGWB(K,L) = 0.D+0
          XMGAB(K,L) = 0.D+0
          XMGOB(K,L) = 0.D+0
          XMGWB(K,L) = 0.D+0
          XLAB(K,L) = 0.D+0
          XLOB(K,L) = 0.D+0
          XLWB(K,L) = 0.D+0
          XMLAB(K,L) = 0.D+0
          XMLOB(K,L) = 0.D+0
          XMLWB(K,L) = 0.D+0
          RHOLB(K,L) = 0.D+0
          RHOGB(K,L) = 0.D+0
          RHONB(K,L) = 0.D+0
          TORLB(K,L) = 0.D+0
          VISLB(K,L) = 0.D+0
          RKLB(1,K,L) = 0.D+0
          RKLB(2,K,L) = 0.D+0
          RKLB(3,K,L) = 0.D+0
          DFLOB(K,L) = 0.D+0
          DFLAB(K,L) = 0.D+0
          DFLSB(K,L) = 0.D+0
          SIB(K,L) = 0.D+0
          PIB(K,L) = 0.D+0
          TMSB(K,L) = 0.D+0
          POSB(K,L) = 0.D+0
          BTGLB(K,L) = 0.D+0
          BTGNB(K,L) = 0.D+0
          BTNLB(K,L) = 0.D+0
        ENDDO
        SDPMB(L) = 0.D+0
        SDPFB(L) = 0.D+0
!
!---    Detection of missing initial boundary conditions  ---
!
        TB(1,L) = -1.D+9
        PLB(1,L) = -1.D+9
        PGB(1,L) = -1.D+9
        YLSB(1,L) = -1.D+9
      ENDDO
      DO L = 1,LBC
        DO K = 1,LSPBC+1
          IBCSP(K,L) = 0
        ENDDO
      ENDDO
!
!---  [ BCVT ] Energy boundary condition variables  ---
!
      DO K = 1,LBCT
        DO L = 1,LSV
          THKLB(L,K) = 0.D+0
          THKGB(L,K) = 0.D+0
          THKNB(L,K) = 0.D+0
          HLB(L,K) = 0.D+0
          HGB(L,K) = 0.D+0
          HNB(L,K) = 0.D+0
          UEGB(L,K) = 0.D+0
          HGAB(L,K) = 0.D+0
          HGOB(L,K) = 0.D+0
          HGWB(L,K) = 0.D+0
          HLAB(L,K) = 0.D+0
          HLOB(L,K) = 0.D+0
          HLSB(L,K) = 0.D+0
          HLWB(L,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  [ BCVG ] Gas boundary condition variables  ---
!
      DO K = 1,LBCG
        DO L = 1,LSV
          TORGB(L,K) = 0.D+0
          VISGB(L,K) = 0.D+0
          RKGB(L,K) = 0.D+0
          DFGWB(L,K) = 0.D+0
          DFGOB(L,K) = 0.D+0
          DFGAB(L,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  [ BCVN ] NAPL boundary condition variables  ---
!
      DO 520 K = 1,LBC
        DO 510 L = 1,LSV
          TORNB(L,K) = 0.D+0
          VISNB(L,K) = 0.D+0
          RKNB(L,K) = 0.D+0
  510   CONTINUE
  520 CONTINUE
      DO 540 K = 1,LBCN
        DO 530 L = 1,LSV
          XSOB(L,K) = 0.D+0
  530   CONTINUE
  540 CONTINUE
!
!---  [ BCVI ] Ice boundary condition variables  ---
!
      DO 700 K = 1,LBCI
        DO 610 L = 1,LSV
          RHOIB(L,K) = 0.D+0
          THKIB(L,K) = 0.D+0
          HIB(LSV,K) = 0.D+0
  610   CONTINUE
  700 CONTINUE
!
!---  [ BCVS ] Salt/Surfactant boundary condition variables  ---
!
      DO 800 K = 1,LBCS
        DO 710 L = 1,LSV
          XLSB(L,K) = 0.D+0
          XMLSB(L,K) = 0.D+0
          OECB(L,K) = 0.D+0
          YLSB(L,K) = 0.D+0
  710   CONTINUE
  800 CONTINUE
!
!---  [ BCVA ] Nonaqueous liquid boundary condition variables  ---
!
      DO K = 1,LBCA
        DO L = 1,LSV
          XNAB(L,K) = 0.D+0
          XNOB(L,K) = 0.D+0
          XNWB(L,K) = 0.D+0
          XMNAB(L,K) = 0.D+0
          XMNOB(L,K) = 0.D+0
          XMNWB(L,K) = 0.D+0
          DFNAB(L,K) = 0.D+0
          DFNOB(L,K) = 0.D+0
          DFNWB(L,K) = 0.D+0
        ENDDO
      ENDDO
      DO K = 1,LBCN2
        DO L = 1,LSV
          XNNB(L,K) = 0.D+0
          XMNNB(L,K) = 0.D+0
          DFNNB(L,K) = 0.D+0
          DFLNB(L,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  [ BCVH ] Hydrate boundary condition variables  ---
!
      DO 920 K = 1,LBCH
        DO 910 L = 1,LSV
          XHWB(L,K) = 0.D+0
          XHAB(L,K) = 0.D+0
          XHOB(L,K) = 0.D+0
          RHOHB(L,K) = 0.D+0
          THKHB(L,K) = 0.D+0
          HHB(L,K) = 0.D+0
          SHB(L,K) = 0.D+0
          PHB(L,K) = 0.D+0
          ZLAB(L,K) = 0.D+0
          PVHOB(L,K) = 0.D+0
          PVHAB(L,K) = 0.D+0
  910   CONTINUE
  920 CONTINUE
      DO 940 L = 1,LBCN2
        DO 930 K = 1,LSV
          DFGNB(K,L) = 0.D+0
          XHNB(K,L) = 0.D+0
          PVHNB(K,L) = 0.D+0
          PVNB(K,L) = 0.D+0
          XGNB(K,L) = 0.D+0
          XMGNB(K,L) = 0.D+0
          XLNB(K,L) = 0.D+0
          XMLNB(K,L) = 0.D+0
  930   CONTINUE
  940 CONTINUE
!
!---  [ BCVGC ] Gas component boundary condition variables  ---
!
      DO 1104 L = 1,LBCGC
        DO 1102 K = 1,LSV
          DO 1100 J = 1,LNGC
            PVCB(J,K,L) = 0.D+0
            XGCB(J,K,L) = 0.D+0
            XMGCB(J,K,L) = 0.D+0
            XNCB(J,K,L) = 0.D+0
            XMNCB(J,K,L) = 0.D+0
            XLCB(J,K,L) = 0.D+0
            XMLCB(J,K,L) = 0.D+0
            DFLCB(J,K,L) = 0.D+0
            DFGCB(J,K,L) = 0.D+0
            DFNCB(J,K,L) = 0.D+0
            HGCB(J,K,L) = 0.D+0
            ZMCB(J,K,L) = 0.D+0
 1100     CONTINUE
          ZGB(K,L) = 0.D+0
          ZNB(K,L) = 0.D+0
 1102   CONTINUE
 1104 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_BOUN group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_GEO_MECH
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
!     Initialize geomechanical variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 8 September 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GEO_MECH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      INTEGER NKX(8,8)
!
!----------------------Data Statements---------------------------------!
!
!      DATA NKX /14,13,11,10,5,4,2,1,
!     &          15,14,12,11,6,5,3,2,
!     &          17,16,14,13,8,7,5,4,
!     &          18,17,15,14,9,8,6,5,
!     &          23,22,20,19,14,13,11,10,
!     &          24,23,21,20,15,14,12,11,
!     &          26,25,23,22,17,16,14,13,
!     &          27,26,24,23,18,17,15,14 /
      DATA NKX /14,15,17,18,23,24,26,27,
     &          13,14,16,17,22,23,25,26,
     &          11,12,14,15,20,21,23,24,
     &          10,11,13,14,19,20,22,23,
     &          5,6,8,9,14,15,17,18,
     &          4,5,7,8,13,14,16,17,
     &          2,3,5,6,11,12,14,15,
     &          1,2,4,5,10,11,13,14 /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_GEO_MECH'
      DO L = 1,LBCIN_GM
        DO K = 1,LBTM_GM
          DO J = 1,LBCV_GM
            BC_GM(J,K,L) = 0.D+0
          ENDDO
        ENDDO
        IBCM_GM(L) = 0
        IBCR_GM(L) = 0
      ENDDO
      DO L = 1,LRC
        DO K = 1,9
          PROP_GM(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          IPROP_GM(K,L) = 0
        ENDDO
        IHCM_GM(L) = 0
      ENDDO
      DO L = 1,LFEN
        DO K = 1,2
          U_GM(K,L) = 0.D+0
          V_GM(K,L) = 0.D+0
          W_GM(K,L) = 0.D+0
        ENDDO
        DO K = 1,8
          NE_GM(K,L) = 0
        ENDDO
        IM_GM(L) = 0
      ENDDO
      DO L = 1,LFD
        DO K = 1,6
          EPS_GM(K,L) = 0.D+0
          SIG_GM(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LFDM
        IXP_GM(L) = 1
        DO K = 1,8
          ND_GM(K,L) = 0
        ENDDO
        DO K = 1,3
          P_GM(K,L) = 0.D+0
          EPSV_GM(K,L) = 0.D+0
          SIGV_GM(K,L) = 0.D+0
        ENDDO
        EPSV_CMP(L) = 0.D+0
        SIGV_CMP(L) = 0.D+0
      ENDDO
      DO K = 1,8
        DO L = 1,8
          NK_GM(K,L) = NKX(K,L)
        ENDDO
      ENDDO
      DO L = 1,LJH_GM
        DO K = 1,LJG_GM
          KLU_GM(K,L) = 0
        ENDDO
      ENDDO
      DO K = 1,LJO_GM
        MLU_GM(K) = 0
      ENDDO
      DO K = 1,LJC_GM
        NLU_GM(K) = 0
      ENDDO
      DO L = 1,LBC_GM
        IBCC_GM(L) = 0
        IBCN_GM(L) = 0
        IBCIN_GM(L) = 0
        IBCD_GM(L) = 0
        DO K = 1,3
          IBCT_GM(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LEPD
        RSDM_GM(L) = 1.D-6
      ENDDO
      IREF_GM = 0
      NNZR_GM = 0
      MK_GM = 0
      ML_GM = 0
      MU_GM = 0
      MD_GM = 0
      NFEN_GM = 0
      NBC_GM = 0
      GRAV_GM = 9.81D+0
      RSD_GM = 0.D+0
      NSD_GM = 0
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_GEO_MECH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVP_BH
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
!     Initialize borehole field variables.
!
!     DNR_BH(LUK,LBN_BH) - primary variable increments of borehole node
!     DFGW_BH(LUK,LBN_BH) - gas water diffusion coefficient of 
!       borehole node
!     DFLA_BH(LUK,LBN_BH) - aqueous air diffusion coefficient of 
!       borehole node
!     DFLS_BH(LUK,LBN_BH) - aqueous salt diffusion coefficient of 
!       borehole node
!     HGW_BH(LUK,LBN_BH) - water vapor enthalpy of borehole node
!     HGA_BH(LUK,LBN_BH) - air enthalpy of borehole node
!     HG_BH(LUK,LBN_BH) - gas enthalpy of borehole node
!     HL_BH(LUK,LBN_BH) - aqueous enthalpy of borehole node
!     HSP_BH(LUK,LBN_BH) - precipitated salt enthalpy of borehole node
!     PERMRF_BH(LUK,LBN_BH) - permeability reduction facture of 
!       borehole node
!     PL_BH(LUK,LBN_BH) - aqueous pressure of borehole node
!     PG_BH(LUK,LBN_BH) - gas pressure of borehole node
!     PORD_BH(LUK,LBN_BH) - diffusive porosity of borehole node
!     PSW_BH(LUK,LBN_BH) - saturated water pressure of borehole node
!     PVA_BH(LUK,LBN_BH) - partial pressure of air of borehole node
!     PVW_BH(LUK,LBN_BH) - partial pressure of water of borehole node
!     RHOG_BH(LUK,LBN_BH) - gas density of borehole node
!     RHOL_BH(LUK,LBN_BH) - aqueous density of borehole node
!     RHOMG_BH(LUK,LBN_BH) - gas molar density of borehole node
!     RHOML_BH(LUK,LBN_BH) - aqueous molar density of borehole node
!     RHOSP_BH(LUK,LBN_BH) - precipitated salt density of borehole node
!     RKG_BH(LUK,LBN_BH) - gas relative permeability of borehole node
!     RKL_BH(LUK,LBN_BH) - aqueous rel. permeability of borehole node
!     SG_BH(LUK,LBN_BH) - gas saturation of borehole node
!     SL_BH(LUK,LBN_BH) - aqueous saturation of borehole node
!     SS_BH(LUK,LBN_BH) - precipitated salt saturation of borehole node
!     T_BH(LUK,LBN_BH) - temperature of borehole node
!     THKG_BH(LUK,LBN_BH) - gas thermal conductivity of borehole node
!     THKL_BH(LUK,LBN_BH) - aqueous thermal conduc. of borehole node
!     T_BH(LUK,LBN_BH) - temperature of borehole node
!     TMS_BH(LUK,LBN_BH) - total salt mass of borehole node
!     UEG_BH(LUK,LBN_BH) - gas internal energy of borehole node
!     VISG_BH(LUK,LBN_BH) - gas viscosity of borehole node
!     VISL_BH(LUK,LBN_BH) - aqueous viscosity of borehole node
!     XGA_BH(LUK,LBN_BH) - gas air mass fraction of borehole node
!     XGW_BH(LUK,LBN_BH) - gas water mass fraction of borehole node
!     XLA_BH(LUK,LBN_BH) - aqueous air mass fraction of borehole node
!     XLS_BH(LUK,LBN_BH) - aqueous salt mass fraction of borehole node
!     XLW_BH(LUK,LBN_BH) - aqueous water mass fraction of borehole node
!     YLS_BH(LUK,LBN_BH) - aqueous total salt mass fraction of 
!       borehole node
!     XLW_BH(LUK,LBN_BH) - aqueous water mass fraction of borehole node
!     XMGA_BH(LUK,LBN_BH) - gas air mole fraction of borehole node
!     XMGW_BH(LUK,LBN_BH) - gas water mole fraction of borehole node
!     XMLA_BH(LUK,LBN_BH) - aqueous air mole fraction of borehole node
!     XMLS_BH(LUK,LBN_BH) - aqueous salt mole fraction of borehole node
!     XMLW_BH(LUK,LBN_BH) - aqueous water mole fraction of borehole node
!     NPHAZ_BH(LUK,LBN_BH) - phase condition of borehole node
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVP_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVP_BH'
!
!---  Loop over the number of borehole nodes ---
!
      DO L = 1,LBN_BH
        DO K = 1,LUK
          DNR_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          POR0_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          DFGW_BH(K,L) = 0.D+0
          DFLA_BH(K,L) = 0.D+0
          DFLS_BH(K,L) = 0.D+0
          HGW_BH(K,L) = 0.D+0
          HGA_BH(K,L) = 0.D+0
          HG_BH(K,L) = 0.D+0
          HL_BH(K,L) = 0.D+0
          HSP_BH(K,L) = 0.D+0
          PERMRF_BH(K,L) = 1.D+0
          PL_BH(K,L) = 0.D+0
          PG_BH(K,L) = 0.D+0
          PORD_BH(K,L) = 0.D+0
          PSW_BH(K,L) = 0.D+0
          PVA_BH(K,L) = 0.D+0
          PVW_BH(K,L) = 0.D+0
          RHOG_BH(K,L) = 0.D+0
          RHOL_BH(K,L) = 0.D+0
          RHOMG_BH(K,L) = 0.D+0
          RHOML_BH(K,L) = 0.D+0
          RHOSP_BH(K,L) = 0.D+0
          RKG_BH(K,L) = 0.D+0
          RKL_BH(K,L) = 0.D+0
          SG_BH(K,L) = 0.D+0
          SL_BH(K,L) = 0.D+0
          SS_BH(K,L) = 0.D+0
          T_BH(K,L) = 0.D+0
          THKG_BH(K,L) = 0.D+0
          THKL_BH(K,L) = 0.D+0
          TMS_BH(K,L) = 0.D+0
          UEG_BH(K,L) = 0.D+0
          VISG_BH(K,L) = 0.D+0
          VISL_BH(K,L) = 0.D+0
          XGA_BH(K,L) = 0.D+0
          XGW_BH(K,L) = 0.D+0
          XLA_BH(K,L) = 0.D+0
          XLS_BH(K,L) = 0.D+0
          YLS_BH(K,L) = 0.D+0
          XLW_BH(K,L) = 0.D+0
          XMGA_BH(K,L) = 0.D+0
          XMGW_BH(K,L) = 0.D+0
          XMLA_BH(K,L) = 0.D+0
          XMLS_BH(K,L) = 0.D+0
          XMLW_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSU
          NPHAZ_BH(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVP_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVP_FRC
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
!     Initialize fracture/fault field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVP_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVP_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LUK
          DNR_FRC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          T_FRC(K,L) = 0.D+0
          PL_FRC(K,L) = 0.D+0
          PG_FRC(K,L) = 0.D+0
          PN_FRC(K,L) = 0.D+0
          PSO_FRC(K,L) = 0.D+0
          SG_FRC(K,L) = 0.D+0
          SL_FRC(K,L) = 0.D+0
          SN_FRC(K,L) = 0.D+0
          PORD_FRC(K,L) = 0.D+0
          PORT_FRC(K,L) = 0.D+0
          PSW_FRC(K,L) = 0.D+0
          PVA_FRC(K,L) = 0.D+0
          PVO_FRC(K,L) = 0.D+0
          PVW_FRC(K,L) = 0.D+0
          XGA_FRC(K,L) = 0.D+0
          XGO_FRC(K,L) = 0.D+0
          XGW_FRC(K,L) = 0.D+0
          XMGA_FRC(K,L) = 0.D+0
          XMGO_FRC(K,L) = 0.D+0
          XMGW_FRC(K,L) = 0.D+0
          XLA_FRC(K,L) = 0.D+0
          XLO_FRC(K,L) = 0.D+0
          XLW_FRC(K,L) = 0.D+0
          XMLA_FRC(K,L) = 0.D+0
          XMLO_FRC(K,L) = 0.D+0
          XMLW_FRC(K,L) = 0.D+0
          RHOG_FRC(K,L) = 0.D+0
          RHOL_FRC(K,L) = 0.D+0
          RHON_FRC(K,L) = 0.D+0
          VISL_FRC(K,L) = 0.D+0
          TORL_FRC(K,L) = 0.D+0
          RKL_FRC(K,L) = 0.D+0
          DFLO_FRC(K,L) = 0.D+0
          DFLA_FRC(K,L) = 0.D+0
          DFLS_FRC(K,L) = 0.D+0
          SI_FRC(K,L) = 0.D+0
          PI_FRC(K,L) = 0.D+0
          RHOMG_FRC(K,L) = 0.D+0
          RHOML_FRC(K,L) = 0.D+0
          RHOMN_FRC(K,L) = 0.D+0
          TMS_FRC(K,L) = 0.D+0
          POSM_FRC(K,L) = 0.D+0
          PERMRF_FRC(K,L) = 1.D+0
          APH_FRC(K,L) = 0.D+0
          APM_FRC(K,L) = 0.D+0
          PERM_FRC(K,L) = 0.D+0
          ASLMIN_FRC(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          POR0_FRC = 0.D+0
        ENDDO
        DO K = 1,LSU
          NPHAZ_FRC(K,L) = 0
        ENDDO
        PCMP_FRC(L) = 0.D+0
        TCMP_FRC(L) = 0.D+0
        VOL_FRC(L) = 0.D+0
        ASL_FRC(L) = 0.D+0
        AST_FRC(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVP_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVT_FRC
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
!     Initialize fracture/fault energy field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVT_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVT_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LSV
          THKG_FRC(K,L) = 0.D+0
          THKL_FRC(K,L) = 0.D+0
          THKN_FRC(K,L) = 0.D+0
          HL_FRC(K,L) = 0.D+0
          HG_FRC(K,L) = 0.D+0
          HN_FRC(K,L) = 0.D+0
          HGA_FRC(K,L) = 0.D+0
          HGO_FRC(K,L) = 0.D+0
          HGW_FRC(K,L) = 0.D+0
          HLA_FRC(K,L) = 0.D+0
          HLO_FRC(K,L) = 0.D+0
          HLS_FRC(K,L) = 0.D+0
          HLW_FRC(K,L) = 0.D+0
          UEG_FRC(K,L) = 0.D+0
          UEL_FRC(K,L) = 0.D+0
          UEN_FRC(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVT_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVG_FRC
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
!     Initialize fracture/fault gas field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVG_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVG_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LSV
          GNIFT_FRC(K,L) = 0.D+0
          VISG_FRC(K,L) = 0.D+0
          TORG_FRC(K,L) = 0.D+0
          RKG_FRC(K,L) = 0.D+0
          DFGW_FRC(K,L) = 0.D+0
          DFGO_FRC(K,L) = 0.D+0
          DFGA_FRC(K,L) = 0.D+0
        ENDDO
        ICAIR_FRC(L) = 0
        ICCO2_FRC(L) = 0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVG_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVN_FRC
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
!     Initialize fracture/fault NAPL field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVN_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVN_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LSV
          VISN_FRC(K,L) = 0.D+0
          TORN_FRC(K,L) = 0.D+0
          RKN_FRC(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVN_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVS_FRC
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
!     Initialize fracture/fault salt field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVS_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVS_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LSV
          XLS_FRC(K,L) = 0.D+0
          XMLS_FRC(K,L) = 0.D+0
          YLS_FRC(K,L) = 0.D+0
          RHOSP_FRC(K,L) = 0.D+0
          HSP_FRC(K,L) = 0.D+0
          SS_FRC(K,L) = 0.D+0
        ENDDO
        ICBRN_FRC(L) = 0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVS_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FDVGC_FRC
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
!     Initialize fracture/fault gas component field variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 March 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVGC_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FDVGC_FRC'
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,6
          BETA_FRC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LNGC
          FK_FRC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          DO J = 1,LNGC
            XGC_FRC(J,K,L) = 0.D+0
            XMGC_FRC(J,K,L) = 0.D+0
            XNC_FRC(J,K,L) = 0.D+0
            XMNC_FRC(J,K,L) = 0.D+0
            XLC_FRC(J,K,L) = 0.D+0
            XMLC_FRC(J,K,L) = 0.D+0
            DFGC_FRC(J,K,L) = 0.D+0
            DFNC_FRC(J,K,L) = 0.D+0
            DFLC_FRC(J,K,L) = 0.D+0
            HGC_FRC(J,K,L) = 0.D+0
            XMVGC_FRC(J,K,L) = 0.D+0
            XMVLC_FRC(J,K,L) = 0.D+0
            TMC_FRC(J,K,L) = 0.D+0
            ZMC_FRC(J,K,L) = 0.D+0
          ENDDO
          XMVGW_FRC(K,L) = 0.D+0
          XMVLW_FRC(K,L) = 0.D+0
          ZG_FRC(K,L) = 0.D+0
          ZN_FRC(K,L) = 0.D+0
        ENDDO
        IZMC_FRC(L) = 0.D+0
        IBETA_FRC(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FDVGC_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FLUX_BH
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
!     Initialize borehole node to borehole node flux variables
!
!     UBBC(LBC_BH,LSOLU+LSPT) solute flux 1/m^2 s
!     UBBDGW(LSFV,LBC_BH) diffusive water vapor flux through gas, 
!       kmol/m^2 s
!     UBBDLA(LSFV,LBC_BH) diffusive air flux through aqu., kmol/m^2 s
!     UBBDS(LSFV,LBC_BH) diffusive salt flux through aqu., kg/m^2 s
!     UBBS(LSFV,LBC_BH) salt flux through aqu., kg/m^2 s
!     UBBG(LSFV,LBC_BH) gas flux, m/s
!     UBBGW(LSFV,LBC_BH) water vapor flux through gas, kg/m^2 s
!     UBBL(LSFV,LBC_BH) aqueous flux, m/s
!     UBBL(LSFV,LBC_BH) nonaqueous flux, m/s
!     UBBQ(LSFV,LBC_BH) heat flux, W/m^2
!     UBBQC(LSFV,LBN_BH) coaxial heat flux, W/m^2
!
!     Initialize borehole node to fracture/fault triangle flux variables
!
!     UBFC(LFC_BH,LSOLU+LSPT) solute flux 1/m^2 s
!     UBFDGW(LSFV,LFC_BH) diffusive water vapor flux through gas, 
!       kmol/m^2 s
!     UBFDLA(LSFV,LFC_BH) diffusive air flux through aqu., kmol/m^2 s
!     UBFDS(LSFV,LFC_BH) diffusive salt flux through aqu., kg/m^2 s
!     UBFS(LSFV,LFC_BH) salt flux through aqu., kg/m^2 s
!     UBFG(LSFV,LFC_BH) gas flux, m/s
!     UBFGW(LSFV,LFC_BH) water vapor flux through gas, kg/m^2 s
!     UBFL(LSFV,LFC_BH) aqueous flux, m/s
!     UBFL(LSFV,LFC_BH) nonaqueous flux, m/s
!     UBFQ(LSFV,LFC_BH) heat flux, W/m^2
!
!     Borehole node to field node flux variables
!
!     UBML(LBN_BH) aqueous flow rate, m^3/s
!     UBMG(LBN_BH) gas flow rate, m^3/s
!     UBMT(LBN_BH) heat flow rate, W
!     UBMC(LBN_BH,LSOLU+LSPT) tracer flow rate, 1/s
!     TRNSA_BH(LSFV,LBN_BH) air mass transfer rate between borehole node
!       and grid cell, kg/s
!     TRNSS_BH(LSFV,LBN_BH) salt mass transfer rate between borehole 
!       node and grid cell, kg/s
!     TRNSQ_BH(LSFV,LBN_BH) energy transfer rate between borehole node
!       and grid cell, kg/s
!     TRNSW_BH(LSFV,LBN_BH) water mass transfer rate between borehole 
!       node and grid cell, kg/s
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FLUX_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FLUX_BH'
!
!---  Loop over the number of borehole node to borehole node
!     connections ---
!
      DO L = 1,LBC_BH
        DO K = 1,LSOLU+LSPT
          UBBC(L,K) = 0.D+0
        ENDDO
        DO K = 1,LSFV
          UBBDGW(K,L) = 0.D+0
          UBBDLA(K,L) = 0.D+0
          UBBS(K,L) = 0.D+0
          UBBDS(K,L) = 0.D+0
          UBBG(K,L) = 0.D+0
          UBBGW(K,L) = 0.D+0
          UBBL(K,L) = 0.D+0
          UBBN(K,L) = 0.D+0
          UBBQ(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,(LBN_BH**LCOAX_BH)
        DO K = 1,LSFV
          UBBQC(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over the number of borehole node to fracture/fault triangle
!     connections ---
!
      DO L = 1,LFC_BH
        DO K = 1,LSOLU+LSPT
          UBFC(L,K) = 0.D+0
        ENDDO
        DO K = 1,LSFV
          UBFDGW(K,L) = 0.D+0
          UBFDLA(K,L) = 0.D+0
          UBFS(K,L) = 0.D+0
          UBFDS(K,L) = 0.D+0
          UBFG(K,L) = 0.D+0
          UBFGW(K,L) = 0.D+0
          UBFL(K,L) = 0.D+0
          UBFN(K,L) = 0.D+0
          UBFQ(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over the number of borehole nodes ---
!
      DO L = 1,LBN_BH
        UBML(L) = 0.D+0
        UBMG(L) = 0.D+0
        UBMT(L) = 0.D+0
        DO K = 1,LSFV
          TRNSA_BH(K,L) = 0.D+0
          TRNSS_BH(K,L) = 0.D+0
          TRNSQ_BH(K,L) = 0.D+0
          TRNSW_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSOLU+LSPT 
          UBMC(L,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FLUX_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FLUX_FRC
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
!     Initialize fracture/fault flux variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 12 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FLUX_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FLUX_FRC'
!
!---  Loop over the number of fracture/fault triangle connections ---
!
      DO L = 1,LFC_FRC
        DO K = 1,LSFV
          UFFDGW(K,L) = 0.D+0
          UFFDLA(K,L) = 0.D+0
          UFFDLW(K,L) = 0.D+0
          UFFDS(K,L) = 0.D+0
          UFFS(K,L) = 0.D+0
          UFFG(K,L) = 0.D+0
          UFFGW(K,L) = 0.D+0
          UFFLA(K,L) = 0.D+0
          UFFLW(K,L) = 0.D+0
          UFFL(K,L) = 0.D+0
          UFFN(K,L) = 0.D+0
          UFFQ(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LNC_FRC
        UFMG(L) = 0.D+0
        UFML(L) = 0.D+0
        UFMN(L) = 0.D+0
        UFMT(L) = 0.D+0
        DO K = 1,LSFV
          TRNSA_FRC(K,L) = 0.D+0
          TRNSS_FRC(K,L) = 0.D+0
          TRNSQ_FRC(K,L) = 0.D+0
          TRNSW_FRC(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LSOLU+LSPT
        DO K = 1,LFC_FRC
          UFFC(K,L) = 0.D+0
        ENDDO
        DO K = 1,LNC_FRC
          UFMC(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FLUX_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_FLUXC_FRC
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
!     Initialize component fracture/fault flux variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 16 April 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FLUXC_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_FLUXC_FRC'
!
!---  Loop over the number of fracture/fault triangle connections ---
!
      DO L = 1,(LFC_FRC**LGC)
        DO K = 1,LSFV
          DO J = 1,LNGC
            UFFDGC(J,K,L) = 0.D+0
            UFFDNC(J,K,L) = 0.D+0
            UFFGC(J,K,L) = 0.D+0
            UFFNC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
      DO L = 1,(LNC_FRC**LGC)
        DO K = 1,LSFV
          DO J = 1,LNGC
            TRNSC_FRC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_FLUXC_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_GEOM_BH
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
!     Initialize borehole geometry variables
!
!     VOL_BH(:) volume of borehole node
!     DBF_BH(:) distance of borehole node to fracture/fault triangle 
!       connection
!     DBB_BH(:) distance of borehole node to borehole node connection
!     DBN_BH(LBN_BH) distance of borehole node to field node
!     PLX_BH(LBN_BH) - x-direction well projection on local grid
!     PLY_BH(LBN_BH) - y-direction well projection on local grid
!     PLZ_BH(LBN_BH) - z-direction well projection on local grid
!     XE_BH(:,:) x vertices (10x radius) of borehole node, m
!     YE_BH(:,:) y vertices (10x radius) of borehole node, m
!     ZE_BH(:,:) z vertices (10x radius) of borehole node, m
!     XP_BH(:,:) x starting and ending point of borehole node
!     YP_BH(:,:) y starting and ending point of borehole node
!     ZP_BH(:,:) z starting and ending point of borehole node
!     XTP_BH(2,LI_BH) - x-transition points for borehole interval, m
!     YTP_BH(2,LI_BH) - y-transition points for borehole interval, m
!     ZTP_BH(2,LI_BH) - z-transition points for borehole interval, m
!     XP_BF(LBC_FRC) x coor. borehole node to frac. tri. connection, m
!     YP_BF(LBC_FRC) y coor. borehole node to frac. tri. connection, m
!     ZP_BF(LBC_FRC) z coor. borehole node to frac. tri. connection, m
!     NM_BH(:) borehole name
!     IFCM_BH(LBN_BH) borehole node to grid cell connection map
!     IBCM_BH(2,LBN_BH) borehole node to borehole node connection map
!     IT_BH(2,LN_BH) starting and ending borehole node type indices
!     Borehole types
!     1 - gas injection w/ zero water
!     2 - gas injection w/ water mass frac.
!     3 - gas injection w/ water rel. humidity
!     4 - aqueous injection w/ zero air and salt
!     5 - aqueous injection w/ air and salt mass frac.
!     6 - aqueous injection w/ air and salt rel. sat.
!     -1 - production
!    IS_BH(LI_BH) borehole interval type
!    Borehole interval types
!    0 - cased filled borehole interval
!    1 - open filled borehole interval
!    2 - screened pipe interval
!    3 - closed pipe interval
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GEOM_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_GEOM_BH'
!
!---  Loop over the number of borehole node to borehole node
!     connections ---
!
      DO L = 1,LBC_BH
        IBCM_BH(L) = 0
        DBB_BH(L) = 0.D+0
        DBBM_BH(L) = 0.D+0
      ENDDO
!
!---  Loop over the number of borehole node to fracture/fault triangle
!     connections ---
!
      DO L = 1,LFC_BH
        DO M = 1,3
          SRFN_BH(3,L) = 0.D+0
        ENDDO
        DBF_BH(L) = 0.D+0
        XP_BF(L) = 0.D+0
        YP_BF(L) = 0.D+0
        ZP_BF(L) = 0.D+0
      ENDDO
!
!---  Loop over the number of boreholes ---
!
      DO L = 1,LN_BH
        IT_BH(1,L) = 0
        IT_BH(2,L) = 0
        NM_BH(L) = 'null'
      ENDDO
!
!---  Loop over the number of borehole nodes ---
!
      DO L = 1,LBN_BH
        DBN_BH(L) = 0.D+0
        DIST_BH(L) = 0.D+0
        PLX_BH(L) = 0.D+0
        PLY_BH(L) = 0.D+0
        PLZ_BH(L) = 0.D+0
        DO K = 1,8
          XE_BH(K,L) = 0.D+0
          YE_BH(K,L) = 0.D+0
          ZE_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,2
          XP_BH(K,L) = 0.D+0
          YP_BH(K,L) = 0.D+0
          ZP_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,3
          IPB_BH(K,L) = 0
        ENDDO
        INV_BH(L) = 0
        IBN_BH(L) = 0
        IZ_BH(L) = 0
      ENDDO
!
!---  Loop over the number of borehole intervals ---
!
      DO L = 1,LI_BH
        DO K = 1,2
          XTP_BH(K,L) = 0.D+0
          YTP_BH(K,L) = 0.D+0
          ZTP_BH(K,L) = 0.D+0
        ENDDO
        IS_BH(L) = 0
      ENDDO
      N_BH = 0
      NBN_BH = 0
      NXP_BH = 0
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_GEOM_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_GEOM_FRC
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
!     Initialize fracture/fault geometry variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 1 August 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE GEOM_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*4 FORM1
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 /'(I1)'/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_GEOM_FRC'
      NF_FRC = 0
      NT_FRC = 0
      NXP_FRC = 0
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,3
          SIGN_FRC(K,L) = 0.D+0
          SIGS_FRC(K,L) = 0.D+0
          XE_FRC(K,L) = 0.D+0
          YE_FRC(K,L) = 0.D+0
          ZE_FRC(K,L) = 0.D+0
          SFNT_FRC(K,L) = 0.D+0
          ND_FRC(K,L) = 0
        ENDDO
        AF_FRC(L) = 0.D+0
        IZ_FRC(L) = 0
        IXP_FRC(L) = 0
      ENDDO    
!
!---  Loop over the number of fracture/fault triangle to grid cell
!     connections  ---
!
      DO K = 1,LNC_FRC
        AFN_FRC(K) = 0.D+0
        XPFN_FRC(K) = 0.D+0
        YPFN_FRC(K) = 0.D+0
        ZPFN_FRC(K) = 0.D+0
        DFN_FRC(K) = 0.D+0
        DFNM_FRC(K) = 0.D+0
        INCM_FRC(K) = 0
      ENDDO
!
!---  Loop over the number of fracture/fault triangle to 
!     fracture/fault triangle connections  ---
!
      DO K = 1,LFC_FRC
        AFF_FRC(K) = 0.D+0
        DFF_FRC(K) = 0.D+0
        DFFM_FRC(K) = 0.D+0
        ITCM_FRC(K) = 0
        DO J = 1,3
          SFNC_FRC(J,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over the number of borehole to fracture/fault triangle
!     connections  ---
!
      DO K = 1,LBC_FRC
        IBHT_FRC(K) = 0
        IBHN_FRC(K) = 0
      ENDDO
!
!---  Loop over the number of fractures/faults ---
!
      DO L = 1,LF_FRC
        DO K = 1,2
          IP_FRC(K,L) = 0
        ENDDO
        ICX = ICOUNT(LF_FRC)
        WRITE(FORM1(3:3),'(I1)') ICOUNT(L)
        NM_FRC(L) = 'Fracture/Fault # '
        WRITE(NM_FRC(L)(12:11+ICX),FORM1) L
        NTP_FRC(L) = 0
      ENDDO
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,2
          IPF_FRC(K,L) = 0
          IPN_FRC(K,L) = 0
        ENDDO
      ENDDO
!
!---  Loop over the MPI processors ---
!
      DO NP = 1,LP_MPI
        NCFGT_MPI(NP) = 0
        NCFT_MPI(NP) = 0
        NCGT_MPI(NP) = 0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_GEOM_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_PARM_BH
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
!     Initialize borehole parameter variables
!
!     HR_BH(LBN_BH) - heat of thermocatalytic reaction, J/m^3 aqu
!     PCMP_BH(LBN_BH) - reference pressure for compressibility
!     PAR_BH(1,LI_BH) - skin factor for borehole interval
!     PAR_BH(2,LI_BH) - borehole outer radius for interval, m
!     PAR_BH(3,LI_BH) - borehole inner radius for interval, m
!     PAR_BH(4,LI_BH) - borehole wall thermal conductivity, W/m ˚K
!     PAR_BH(5,LI_BH) - borehole wall volumetric heat capacity, J/m^3 ˚K
!     PAR_BH(6,LI_BH) - borehole Colebrook roughness, m
!     PAR_BH(7,LI_BH) - borehole Brooks and Corey psi, m
!     PAR_BH(8,LI_BH) - borehole Brooks and Corey lambda
!     PAR_BH(9,LI_BH) - borehole Brooks and Corey slr
!     PAR_BH(10,LI_BH) - borehole Webb oven dried head, m
!     PAR_BH(11,LI_BH) - borehole Webb saturation matching point
!     PAR_BH(12,LI_BH) - borehole Webb head matching point, m
!     VAR_BH(1,LT_BH,LN_BH) - borehole time for time point, s
!     VAR_BH(2,LT_BH,LN_BH) - mass rate for time point, kg/s
!     VAR_BH(3,LT_BH,LN_BH) - production pressure, Pa
!     VAR_BH(4,LT_BH,LN_BH) - aqueous air conc. for time point
!     VAR_BH(4,LT_BH,LN_BH) - gas water conc. for time point
!     VAR_BH(5,LT_BH,LN_BH) - temperature for time point
!     VAR_BH(6,LT_BH,LN_BH) - aqu. salt conc. for time point
!     VIC_BH(21,NF) borehole initial condition variables
!     ICC_BH(LN_BH) - cyclic borehole time index for borehole
!     IC_BH(NF) fracture/fault initial condition state index
!       1 - aqueous saturation, gas pressure
!       2 - aqueous saturation, aqueous pressure
!       3 - aqueous pressure, gas pressure
!       4 - hydrostatic
!     ID_BH(1,NBH) - starting borehole interval index
!     ID_BH(2,NBH) - ending borehole interval index
!     ID_BH(3,NBH) - starting borehole node index
!     ID_BH(4,NBH) - ending borehole node index
!     ID_BH(5,NBH) - starting frac. triangle index, pointer to 
!       IBHT_FRC and IBHN_FRC
!     ID_BH(6,NBH) - ending frac. triangle index, pointer to 
!       IBHT_FRC and IBHN_FRC
!     ITC_BH - thermo-catalytic fluid equation form
!       22 - Zero-order reaction without catalyst with salt tracking
!       23 - First-order reaction without catalyst with salt tracking
!       24 - Zero-order reaction with catalyst with salt tracking
!       25 - First-order reaction with catalyst with salt tracking
!     ISPH_BH - thermo-catalytic enthalpy tracking solute
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_BH
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_PARM_BH'
!
!---  Loop over the number of boreholes ---
!
      DO L = 1,LN_BH
        DO K = 1,LT_BH
          DO J = 1,11
            VAR_BH(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,LSOLU_BH
            VARC_BH(J,K,L) = 0.D+0
          ENDDO
        ENDDO
        DO K = 1,21
          VIC_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSOLU_BH
          ISOLU_BH(K,L) = 0
        ENDDO
        ITM_BH(L) = 0
        ICC_BH(L) = 0
        IC_BH(L) = 0
        DO K = 1,6
          ID_BH(K,L) = 0
        ENDDO
      ENDDO
!
!---  Loop over the number of unknowns ---
!
      DO L = 1,LUK
        RSD_BH(L) = 0.D+0
        NSD_BH(L) = 0.D+0
      ENDDO
!
!---  Loop over the number of borehole nodes ---
!
      DO L = 1,LBN_BH
        DO K = 1,LUK
          IM_BH(K,L) = 0
          RSDL_BH(K,L) = 0.D+0
        ENDDO
        DO K = 1,LSV
          SRCA_BH(K,L) = 0.D+0
          SRCS_BH(K,L) = 0.D+0
          SRCT_BH(K,L) = 0.D+0
          SRCW_BH(K,L) = 0.D+0
        ENDDO
        HR_BH(L) = 0.D+0
        PCMP_BH(L) = 0.D+0
      ENDDO
!
!---  Loop over the number of borehole intervals ---
!
      DO L = 1,LI_BH
        DO K = 1,16
          PAR_BH(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LREF
        NDREF_BH(L) = 0
      ENDDO
!
!---  Loop over the number of borehole source ---
!
      DO L = 1,LSR_BH
        DO K = 1,LSTM_BH
          DO J = 1,7+LSOLU+LSPT
            SRC_BH(J,K,L) = 0.D+0
          ENDDO
        ENDDO
        ISRM_BH(L) = 0
        ISRT_BH(L) = 0
        DO K = 1,2
          ISRDM_BH(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_BH
        DO K = 1,LJG_BH
          KLU_BH(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_BH
        DO K = 1,LJK_BH
          KLUC_BH(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_MCB
        DO K = 1,LJG_MCB
          KLU_MCB(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_MCB
        DO K = 1,LJK_MCB
          KLUC_MCB(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_BCM
        DO K = 1,LJG_BCM
          KLU_BCM(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_BCM
        DO K = 1,LJK_BCM
          KLUC_BCM(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_BCF
        DO K = 1,LJG_BCF
          KLU_BCF(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_BCF
        DO K = 1,LJK_BCF
          KLUC_BCF(K,L) = 0
        ENDDO
      ENDDO
!
!---  Thermo-catalytic parameters ---
!
      ITC_BH = 0
      ISPH_BH = 0
      DO L = 1,3
        PTC_BH(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_PARM_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_PARM_FRC
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
!     Initialize fracture/fault parameter variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 6 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PARM_FRC
      USE GEOM_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_PARM_FRC'
      NF_FRC = 0
!
!---  Loop over the number of fractures/faults ---
!
      DO L = 1,LF_FRC
        DO K = 1,7
          DSBB_FRC(K,L) = 0.D+0
          RKSP_FRC(K,L) = 0.D+0
        ENDDO
        DTHCP_FRC(L) = 0.D+0
        SKF_FRC(L) = 0.D+0
        DO K = 1,21
          VIC_FRC(K,L) = 0.D+0
        ENDDO
        IC_FRC(L) = 0
      ENDDO
      DO M = 1,LSR_FRC
        DO L = 1,LSTM_FRC
          DO K = 1,7+LSOLU+LSPT
            SRC_FRC(K,L,M) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO L = 1,LT_FRC
        DO K = 1,LUK
          IM_FRC(K,L) = 0
          RSDL_FRC(K,L) = 0.D+0
        ENDDO
      ENDDO
      DO L = 1,LT_FRC
        JRC_FRC(L) = 0.D+0
        JCS_FRC(L) = 0.D+0
        UCS_FRC(L) = 0.D+0
        THCP_FRC(L) = 0.D+0
        IJM_FRC(L) = 0
        DO K = 1,LSV
          DO J = 1,LNGC
            SRCGC_FRC(J,K,L) = 0.D+0
          ENDDO
          SRCA_FRC(K,L) = 0.D+0
          SRCS_FRC(K,L) = 0.D+0
          SRCT_FRC(K,L) = 0.D+0
          SRCW_FRC(K,L) = 0.D+0
        ENDDO
        SRCG_FRC(L) = 0.D+0
        SRCL_FRC(L) = 0.D+0
      ENDDO
      DO L = 1,LJH_FRC
        DO K = 1,LJG_FRC
          KLU_FRC(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_FRC
        DO K = 1,LJK_FRC
          KLUC_FRC(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_MCF
        DO K = 1,LJG_MCF
          KLU_MCF(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_MCF
        DO K = 1,LJK_MCF
          KLUC_MCF(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_FCM
        DO K = 1,LJG_FCM
          KLU_FCM(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_FCM
        DO K = 1,LJK_FCM
          KLUC_FCM(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJH_FCB
        DO K = 1,LJG_FCB
          KLU_FCB(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LJL_FCB
        DO K = 1,LJK_FCB
          KLUC_FCB(K,L) = 0
        ENDDO
      ENDDO
      NREF_FRC = 0
      DO L = 1,LREF
        NDREF_FRC(L) = 0
      ENDDO
      DO L = 1,LUK
        NSD_FRC(L) = 0
        RSD_FRC(L) = 0.D+0
      ENDDO
      DO L = 1,LSR_FRC
        ISRM_FRC(L) = 0
        ISRT_FRC(L) = 0
        DO K = 1,2
          ISRDM_FRC(K,L) = 0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_PARM_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_TRNS_BH
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
!     Initialize borehole transport variables
!
!     C_BH(LBN_BH,LSOLU+LSPT) volumetric solute or total 
!       species concentration of borehole node
!     CO_BH(LBN_BH,LSOLU+LSPT) volumetric solute or total 
!       species concentration of borehole node at old time step
!     CRNTG_BH(LBN_BH) gas Courant number for fracture/fault
!     CRNTL_BH(LBN_BH) aqueous Courant number for fracture/fault
!     CRNTN_BH(LBN_BH) nonaqueous-liq. Courant number for fracture/fault
!     ICT_BH(LBN_BH,LSOLU+LSPT) initial condition index of solute or 
!       total parameters species of borehole node
!     SP_C_BH(LBN_BH,LSPR) - fracture/fault species concentration, 
!       mol/m^3 node volume
!     SP_CO_BH(LBN_BH,LSPR) - old fracture/fault species concentration, 
!       mol/m^3 node volume
!     YG_BH(LBN_BHC,LSOLU+LSPT) gas mole fraction of solute or total 
!       species of borehole node
!     YL_BH(LBN_BHC,LSOLU+LSPT) aqueous mole fraction of solute or 
!       total species of borehole node
!     YN_BH(LBN_BHC,LSOLU+LSPT) nonaqueous liquid mole fraction of 
!       solute or total species of borehole node
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 18 April 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
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
      SUB_LOG(ISUB_LOG) = '/IN_TRNS_BH'
!
!---  Loop over the number of borehole nodes ---
!
      DO L = 1,LBN_BHC
!
!---    Loop over the number of solutes and transported species ---
!
        DO K = 1,LSOLU+LSPT
          C_BH(L,K) = 0.D+0
          CO_BH(L,K) = 0.D+0
          YG_BH(L,K) = 0.D+0
          YL_BH(L,K) = 0.D+0
          YN_BH(L,K) = 0.D+0
        ENDDO
!
!---    Loop over the number of reactive species ---
!
        DO K = 1,LSPR
          SP_C_BH(L,K) = 0.D+0
          SP_CO_BH(L,K) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over the number of boreholes ---
!
      DO L = 1,LN_BHC
        DO K = 1,LSOLU+LSPT
          ICT_BH(L,K) = 0
        ENDDO
      ENDDO
!
!---  Loop over the number of borehole nodes  ---
!
      DO L = 1,LBN_BH
        CRNTG_BH(L) = 0.D+0
        CRNTL_BH(L) = 0.D+0
        CRNTN_BH(L) = 0.D+0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_TRNS_BH group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_TRNS_FRC
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
!     Initialize fracture/fault transport variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 2 July 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE SOLTN
      USE GEOM_FRC
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_TRNS_FRC'
      NF_FRC = 0
!
!---  Loop over the number of fractures/faults ---
!
      DO K = 1,LF_FRCC
        DO L = 1,LSOLU+LSPT
          ICT_FRC(K,L) = 0
        ENDDO
        DO L = 1,LSPR
          IC_SP_FRC(K,L) = 0
        ENDDO
      ENDDO
!
!---  Loop over the number of fracture/fault triangles ---
!
      DO K = 1,LT_FRCC
        DO L = 1,LSOLU+LSPT
          C_FRC(K,L) = 0.D+0
          CO_FRC(K,L) = 0.D+0
          YG_FRC(K,L) = 0.D+0
          YL_FRC(K,L) = 0.D+0
          YN_FRC(K,L) = 0.D+0
        ENDDO
        DO L = 1,LSPR
          SP_C_FRC(K,L) = 0.D+0
          SP_CO_FRC(K,L) = 0.D+0
        ENDDO
        RHOS_FRC(K) = 0.D+0
        CRNTG_FRC(K) = 0.D+0
        CRNTL_FRC(K) = 0.D+0
        CRNTN_FRC(K) = 0.D+0
      ENDDO
      DO L = 1,LSOLU
        DO K = 1,LT_FRCC
          DO J = 1,5
            PCSL_FRC(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,3
            SDCL_FRC(J,K,L) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_TRNS_FRC group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_OUTP
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
!     Initialize output control variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 June 2001.
!     Last Modified by MD White, PNNL, 19 June 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE TRNSPT
      USE SPILL
      USE SOURC
      USE SOLTN
      USE REACT
      USE OUTPU
      USE NAPL
      USE HYST
      USE GRID
      USE FLUX_FRC
      USE FLUX_BH
      USE FLUXT
      USE FLUXS
      USE FLUXP
      USE FLUXN
      USE FLUXGC
      USE FLUXD
      USE FDVT
      USE FDVS
      USE FDVP
      USE FDVN
      USE FDVH
      USE FDVG
      USE FDVA
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL ICOUNT
!
!----------------------Type Declarations-------------------------------!
!
      CHARACTER*4 FORM1
!
!----------------------Data Statements---------------------------------!
!
      DATA FORM1 / '(I )' /
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_OUTP'
!
!---  [ OUTPU ] Output control variables  ---
!
      DO 10 L = 1,LPTM
        PRTM(L) = 0.D+0
   10 CONTINUE
      DO 30 L = 1,LSF
        SF(1,L) = 0.D+0
        SF(2,L) = 0.D+0
        UNSF(1,L) = 'null'
        UNSF(2,L) = 'null'
        ISFT(L) = 0
        ISFF(L) = 0
        ISFD(L) = 0
        ISFGC(L) = 0
        NSFDOM(L) = 0
        DO 20 K = 1,6
          ISFC(K,L) = 0
   20   CONTINUE
        DO 24 M = 1,LSFDOM
          DO 22 K = 1,4
            ISFDOM(K,M,L) = 0
   22     CONTINUE
   24   CONTINUE
        ISFGP(L) = 0
        ISFSN(L) = 0
   30 CONTINUE
      UNTM = 's'
      UNLN = 'm'
      UNAR = 'rad'
      DO 50 L = 1,LOUPV
        UNPLOT(L) = 'null'
        UNREF(L) = 'null'
        CHREF(L) = 'null'
        DO 40 K = 1,3
          CHSF(L,K) = 'null'
   40   CONTINUE
   50 CONTINUE
      DO 52 L = 1,LVREF
        IREF(L) = 0
        IREFGC(L) = 0
   52 CONTINUE
      DO 54 L = 1,LVPLOT
        IPLOT(L) = 0
        IPLOTGC(L) = 0
   54 CONTINUE
      DO 60 L = 1,LREF
        NDREF(L) = 0
   60 CONTINUE
      NSF = 0
      NSFGP = 1
      IHSF = 0
      NPRTM = 0
      NVPLOT = 0
      NREF = 0
      NVREF = 0
      ICNO = -1
      ICNS = -1
      IFQS = 0
      IFQO = 0
      ISGNS = 4
      ISGNO = 4
      ISGNP = 5
      CHREF(1) = 'PL'
      UNPLOT(1) = 'pa'
      UNREF(1) = 'pa'
      CHREF(2) = 'PG'
      UNPLOT(2) = 'pa'
      UNREF(2) = 'pa'
      CHREF(3) = 'PN'
      UNPLOT(3) = 'pa'
      UNREF(3) = 'pa'
      CHREF(4) = 'T'
      UNPLOT(4) = 'c'
      UNREF(4) = 'c'
      CHREF(5) = 'PHCN'
      CHREF(6) = 'GPL'
      UNPLOT(6) = 'pa'
      UNREF(6) = 'pa'
      CHREF(7) = 'GPG'
      UNPLOT(7) = 'pa'
      UNREF(7) = 'pa'
      CHREF(8) = 'GPN'
      UNPLOT(8) = 'pa'
      UNREF(8) = 'pa'
      CHREF(9) = 'ASL'
      CHREF(10) = 'AST'
      CHREF(11) = 'SL'
      CHREF(12) = 'SG'
      CHREF(13) = 'SN'
      CHREF(14) = 'ST'
      CHREF(15) = 'MCL'
      CHREF(16) = 'MCN'
      CHREF(17) = 'MCT'
      CHREF(18) = 'ESNT'
      CHREF(19) = 'ESGT'
      CHREF(20) = 'PORD'
      CHREF(21) = 'XGW'
      CHREF(22) = 'XGA'
      CHREF(23) = 'XGO'
      CHREF(24) = 'XLW'
      CHREF(25) = 'XLA'
      CHREF(26) = 'XLO'
      CHREF(27) = 'HHL'
      UNPLOT(27) = 'm'
      UNREF(27) = 'm'
      CHREF(28) = 'HHG'
      UNPLOT(28) = 'm'
      UNREF(28) = 'm'
      CHREF(29) = 'HHN'
      UNPLOT(29) = 'm'
      UNREF(29) = 'm'
      CHREF(30) = 'RSZN'
      CHREF(31) = 'RPL'
      CHREF(32) = 'RPG'
      CHREF(33) = 'RPN'
      CHREF(34) = 'RHOL'
      UNPLOT(34) = 'kg/m^3'
      UNREF(34) = 'kg/m^3'
      CHREF(35) = 'RHOG'
      UNPLOT(35) = 'kg/m^3'
      UNREF(35) = 'kg/m^3'
      CHREF(36) = 'RHON'
      UNPLOT(36) = 'kg/m^3'
      UNREF(36) = 'kg/m^3'
      CHREF(37) = 'TMW'
      UNPLOT(37) = 'kg'
      UNREF(37) = 'kg'
      CHREF(38) = 'TMA'
      UNPLOT(38) = 'kg'
      UNREF(38) = 'kg'
      CHREF(39) = 'TMO'
      UNPLOT(39) = 'kg'
      UNREF(39) = 'kg'
      CHREF(40) = 'SRIW'
      UNPLOT(40) = 'kg'
      UNREF(40) = 'kg'
      CHREF(41) = 'SRIA'
      UNPLOT(41) = 'kg'
      UNREF(41) = 'kg'
      CHREF(42) = 'SRIO'
      UNPLOT(42) = 'kg'
      UNREF(42) = 'kg'
      CHREF(43) = 'SRIT'
      UNPLOT(43) = 'j'
      UNREF(43) = 'j'
      CHREF(44) = 'THKX'
      UNPLOT(44) = 'w/m k'
      UNREF(44) = 'w/m k'
      CHREF(45) = 'THKY'
      UNPLOT(45) = 'w/m k'
      UNREF(45) = 'w/m k'
      CHREF(46) = 'THKZ'
      UNPLOT(46) = 'w/m k'
      UNREF(46) = 'w/m k'
      CHREF(47) = 'CS'
      UNPLOT(47) = 'kg/m^3'
      UNREF(47) = 'kg/m^3'
      CHREF(48) = 'CSL'
      UNPLOT(48) = 'kg/m^3'
      UNREF(48) = 'kg/m^3'
      CHREF(49) = 'CRNL'
      CHREF(50) = 'TMS'
      UNPLOT(50) = 'kg'
      UNREF(50) = 'kg'
      CHREF(51) = 'UL'
      UNPLOT(51) = 'm/s'
      UNREF(51) = 'm/s'
      CHREF(52) = 'VL'
      UNPLOT(52) = 'm/s'
      UNREF(52) = 'm/s'
      CHREF(53) = 'WL'
      UNPLOT(53) = 'm/s'
      UNREF(53) = 'm/s'
      CHREF(54) = 'UG'
      UNPLOT(54) = 'm/s'
      UNREF(54) = 'm/s'
      CHREF(55) = 'VG'
      UNPLOT(55) = 'm/s'
      UNREF(55) = 'm/s'
      CHREF(56) = 'WG'
      UNPLOT(56) = 'm/s'
      UNREF(56) = 'm/s'
      CHREF(57) = 'UN'
      UNPLOT(57) = 'm/s'
      UNREF(57) = 'm/s'
      CHREF(58) = 'VN'
      UNPLOT(58) = 'm/s'
      UNREF(58) = 'm/s'
      CHREF(59) = 'WN'
      UNPLOT(59) = 'm/s'
      UNREF(59) = 'm/s'
      CHREF(60) = 'UQ'
      UNPLOT(60) = 'w/m^2 s'
      UNREF(60) = 'w/m^2 s'
      CHREF(61) = 'VQ'
      UNPLOT(61) = 'w/m^2 s'
      UNREF(61) = 'w/m^2 s'
      CHREF(62) = 'WQ'
      UNPLOT(62) = 'w/m^2 s'
      UNREF(62) = 'w/m^2 s'
      CHREF(63) = 'MPH'
      UNPLOT(63) = 'm'
      UNREF(63) = 'm'
      CHREF(64) = 'US'
      UNPLOT(64) = 'kg/m^2 s'
      UNREF(64) = 'kg/m^2 s'
      CHREF(65) = 'VS'
      UNPLOT(65) = 'kg/m^2 s'
      UNREF(65) = 'kg/m^2 s'
      CHREF(66) = 'WS'
      UNPLOT(66) = 'kg/m^2 s'
      UNREF(66) = 'kg/m^2 s'
      CHREF(67) = 'USNC'
      UNPLOT(67) = 'kg/m^2 s'
      UNREF(67) = 'kg/m^2 s'
      CHREF(68) = 'VSNC'
      UNPLOT(68) = 'kg/m^2 s'
      UNREF(68) = 'kg/m^2 s'
      CHREF(69) = 'WSNC'
      UNPLOT(69) = 'kg/m^2 s'
      UNREF(69) = 'kg/m^2 s'
      CHREF(70) = 'XMGW'
      CHREF(71) = 'XMGA'
      CHREF(72) = 'XMGO'
      CHREF(73) = 'CGW'
      UNPLOT(73) = 'kg/m^3'
      UNREF(73) = 'kg/m^3'
      CHREF(74) = 'CGA'
      UNPLOT(74) = 'kg/m^3'
      UNREF(74) = 'kg/m^3'
      CHREF(75) = 'CGO'
      UNPLOT(75) = 'kg/m^3'
      UNREF(75) = 'kg/m^3'
      CHREF(76) = 'CLW'
      UNPLOT(76) = 'kg/m^3'
      UNREF(76) = 'kg/m^3'
      CHREF(77) = 'CLA'
      UNPLOT(77) = 'kg/m^3'
      UNREF(77) = 'kg/m^3'
      CHREF(78) = 'CLO'
      UNPLOT(78) = 'kg/m^3'
      UNREF(78) = 'kg/m^3'
      CHREF(79) = 'CRNG'
      CHREF(80) = 'PI'
      UNPLOT(80) = 'pa'
      UNREF(80) = 'pa'
      CHREF(81) = 'SI'
      CHREF(82) = 'RHOF'
      UNPLOT(82) = 'kg/m^3'
      UNREF(82) = 'kg/m^3'
      CHREF(83) = 'DSLF'
      CHREF(84) = 'DSLM'
      CHREF(85) = 'DSGF'
      CHREF(86) = 'DSGM'
      CHREF(87) = 'ULNC'
      UNPLOT(87) = 'm/s'
      UNREF(87) = 'm/s'
      CHREF(88) = 'VLNC'
      UNPLOT(88) = 'm/s'
      UNREF(88) = 'm/s'
      CHREF(89) = 'WLNC'
      UNPLOT(89) = 'm/s'
      UNREF(89) = 'm/s'
      CHREF(90) = 'UGNC'
      UNPLOT(90) = 'm/s'
      UNREF(90) = 'm/s'
      CHREF(91) = 'VGNC'
      UNPLOT(91) = 'm/s'
      UNREF(91) = 'm/s'
      CHREF(92) = 'WGNC'
      UNPLOT(92) = 'm/s'
      UNREF(92) = 'm/s'
      CHREF(93) = 'UNNC'
      UNPLOT(93) = 'm/s'
      UNREF(93) = 'm/s'
      CHREF(94) = 'VNNC'
      UNPLOT(94) = 'm/s'
      UNREF(94) = 'm/s'
      CHREF(95) = 'WNNC'
      UNPLOT(95) = 'm/s'
      UNREF(95) = 'm/s'
      CHREF(96) = 'UQNC'
      UNPLOT(96) = 'w/m^2 s'
      UNREF(96) = 'w/m^2 s'
      CHREF(97) = 'VQNC'
      UNPLOT(97) = 'w/m^2 s'
      UNREF(97) = 'w/m^2 s'
      CHREF(98) = 'WQNC'
      UNPLOT(98) = 'w/m^2 s'
      UNREF(98) = 'w/m^2 s'
      CHREF(99) = 'CRNN'
      CHREF(100) = 'NODE'
      CHREF(101) = 'POSM'
      UNPLOT(101) = 'pa'
      UNREF(101) = 'pa'
      CHREF(102) = 'OEC'
      CHREF(103) = 'CMPL'
      UNPLOT(103) = '1/pa'
      UNREF(103) = '1/pa'
      CHREF(104) = 'CMPG'
      UNPLOT(104) = '1/pa'
      UNREF(104) = '1/pa'
      CHREF(105) = 'SGT'
      CHREF(106) = 'SNT'
      CHREF(107) = 'SGTL'
      CHREF(108) = 'SGTN'
      CHREF(109) = 'CLO'
      UNREF(109) = 'kg/m^3'
      UNPLOT(109) = 'kg/m^3'
      CHREF(110) = 'XLS'
      CHREF(111) = 'CS'
      UNREF(111) = 'kg/m^3'
      UNPLOT(111) = 'kg/m^3'
      CHREF(112) = 'CLS'
      UNREF(112) = 'kg/m^3'
      UNPLOT(112) = 'kg/m^3'
      CHREF(113) = 'XLS'
      CHREF(114) = 'US'
      UNREF(114) = 'kg/m^2 s'
      UNPLOT(114) = 'kg/m^2 s'
      CHREF(115) = 'VS'
      UNREF(115) = 'kg/m^2 s'
      UNPLOT(115) = 'kg/m^2 s'
      CHREF(116) = 'WS'
      UNREF(116) = 'kg/m^2 s'
      UNPLOT(116) = 'kg/m^2 s'
      CHREF(117) = 'USNC'
      UNREF(117) = 'kg/m^2 s'
      UNPLOT(117) = 'kg/m^2 s'
      CHREF(118) = 'VSNC'
      UNREF(118) = 'kg/m^2 s'
      UNPLOT(118) = 'kg/m^2 s'
      CHREF(119) = 'WSNC'
      UNREF(119) = 'kg/m^2 s'
      UNPLOT(119) = 'kg/m^2 s'
      CHREF(120) = 'ULO'
      UNREF(120) = 'kg/m^2 s'
      UNPLOT(120) = 'kg/m^2 s'
      CHREF(121) = 'VLO'
      UNREF(121) = 'kg/m^2 s'
      UNPLOT(121) = 'kg/m^2 s'
      CHREF(122) = 'WLO'
      UNREF(122) = 'kg/m^2 s'
      UNPLOT(122) = 'kg/m^2 s'
      CHREF(123) = 'ULOC'
      UNREF(123) = 'kg/m^2 s'
      UNPLOT(123) = 'kg/m^2 s'
      CHREF(124) = 'VLOC'
      UNREF(124) = 'kg/m^2 s'
      UNPLOT(124) = 'kg/m^2 s'
      CHREF(125) = 'WLOC'
      UNREF(125) = 'kg/m^2 s'
      UNPLOT(125) = 'kg/m^2 s'
      CHREF(126) = 'TPNL'
      CHREF(127) = 'ESLM'
      CHREF(128) = 'PVW'
      UNPLOT(128) = 'pa'
      UNREF(128) = 'pa'
      CHREF(129) = 'PVA'
      UNPLOT(129) = 'pa'
      UNREF(129) = 'pa'
      CHREF(130) = 'PVO'
      UNPLOT(130) = 'pa'
      UNREF(130) = 'pa'
      CHREF(131) = 'BGL'
      CHREF(132) = 'ANFL'
      UNPLOT(132) = 'm^2'
      UNREF(132) = 'm^2'
      CHREF(132) = 'ANTL'
      UNPLOT(132) = 'm^2'
      UNREF(132) = 'm^2'
      CHREF(134) = 'HKL'
      UNPLOT(134) = '1/s'
      UNREF(134) = '1/s'
      CHREF(135) = 'HKNF'
      UNPLOT(135) = '1/s'
      UNREF(135) = '1/s'
      CHREF(136) = 'HKNT'
      UNPLOT(136) = '1/s'
      UNREF(136) = '1/s'
      CHREF(137) = 'RHW'
      CHREF(138) = 'PLWB'
      UNREF(138) = 'pa'
      CHREF(140) = 'SRCW'
      UNPLOT(140) = 'kg/s'
      UNREF(140) = 'kg/s'
      CHREF(141) = 'SRCA'
      UNPLOT(141) = 'kg/s'
      UNREF(141) = 'kg/s'
      CHREF(142) = 'SRCO'
      UNPLOT(142) = 'kg/s'
      UNREF(142) = 'kg/s'
      CHREF(143) = 'SRCQ'
      UNPLOT(143) = 'w'
      UNREF(143) = 'w'
      CHREF(144) = 'PLWB'
      UNREF(144) = 'm'
      CHREF(145) = 'QLW'
      UNREF(145) = 'm^3'
      CHREF(146) = 'QLWI'
      UNREF(146) = 'm^3/s'
      CHREF(147) = 'SRCS'
      UNPLOT(147) = 'kg/s'
      UNREF(147) = 'kg/s'
      CHREF(148) = 'SRIS'
      UNPLOT(148) = 'kg'
      UNREF(148) = 'kg'
      CHREF(149) = 'PATH'
      CHREF(150) = 'DAPS'
      CHREF(151) = 'BVF'
      CHREF(152) = 'XBA'
      CHREF(153) = 'MCO2'
      UNPLOT(153) = 'kg/m^3'
      UNREF(153) = 'kg/m^3'
      CHREF(154) = 'QNW'
      UNREF(154) = 'm^3'
      CHREF(155) = 'QNWI'
      UNREF(155) = 'm^3/s'
      CHREF(156) = 'QTW'
      UNREF(156) = 'm^3'
      CHREF(157) = 'QTWI'
      UNREF(157) = 'm^3/s'
      CHREF(161) = 'CNW'
      UNPLOT(161) = 'kg/m^3'
      UNREF(161) = 'kg/m^3'
      CHREF(162) = 'XMNW'
      CHREF(163) = 'XNW'
      CHREF(164) = 'CNO'
      UNPLOT(164) = 'kg/m^3'
      UNREF(164) = 'kg/m^3'
      CHREF(165) = 'XMNO'
      CHREF(166) = 'XNO'
      CHREF(167) = 'TMA'
      UNPLOT(167) = 'kg/m^3'
      UNREF(167) = 'kg/m^3'
      CHREF(168) = 'XLW'
      CHREF(169) = 'SRIA'
      UNPLOT(169) = 'kg'
      UNREF(169) = 'kg'
      CHREF(170) = 'SRCA'
      UNPLOT(170) = 'kg/s'
      UNREF(170) = 'kg/s'
      CHREF(171) = 'IMA'
      UNREF(171) = 'kg'
      CHREF(172) = 'IMLA'
      UNREF(172) = 'kg'
      CHREF(173) = 'IMNW'
      UNREF(173) = 'kg'
      CHREF(174) = 'IMNO'
      UNREF(174) = 'kg'
      CHREF(175) = 'IMNA'
      UNREF(175) = 'kg'
      CHREF(176) = 'VISL'
      UNPLOT(176) = 'pa s'
      UNREF(176) = 'pa s'
      CHREF(177) = 'WDT'
      UNPLOT(177) = 'm'
      UNREF(177) = 'm'
      CHREF(178) = 'WDL'
      UNPLOT(178) = 'm'
      UNREF(178) = 'm'
      CHREF(179) = 'WDN'
      UNPLOT(179) = 'm'
      UNREF(179) = 'm'
      CHREF(180) = 'PW'
      UNPLOT(180) = 'pa'
      UNREF(180) = 'pa'
      CHREF(181) = 'SLW'
      CHREF(182) = 'XGWW'
      CHREF(183) = 'XLAW'
      CHREF(184) = 'UL_W'
      UNPLOT(184) = 'm/s'
      UNREF(184) = 'm/s'
      CHREF(185) = 'UG_W'
      UNPLOT(185) = 'm/s'
      UNREF(185) = 'm/s'
      CHREF(186) = 'WL_W'
      UNPLOT(186) = 'm/s'
      UNREF(186) = 'm/s'
      CHREF(187) = 'WG_W'
      UNPLOT(187) = 'm/s'
      UNREF(187) = 'm/s'
      CHREF(188) = 'IPLW'
      UNREF(188) = 'kg'
      CHREF(189) = 'IPNW'
      UNREF(189) = 'kg'
      CHREF(190) = 'IMGT'
      UNREF(190) = 'kg'
      CHREF(191) = 'IMW'
      UNREF(191) = 'kg'
      CHREF(192) = 'IMA'
      UNREF(192) = 'kg'
      CHREF(193) = 'IMO'
      UNREF(193) = 'kg'
      CHREF(194) = 'IMLW'
      UNREF(194) = 'kg'
      CHREF(195) = 'IMLA'
      UNREF(195) = 'kg'
      CHREF(196) = 'IMLO'
      UNREF(196) = 'kg'
      CHREF(197) = 'IMGW'
      UNREF(197) = 'kg'
      CHREF(198) = 'IMGA'
      UNREF(198) = 'kg'
      CHREF(199) = 'IMGO'
      UNREF(199) = 'kg'
      CHREF(201) = 'RKLX'
      CHREF(202) = 'RKLY'
      CHREF(203) = 'RKLZ'
      CHREF(204) = 'XMLA'
      CHREF(205) = 'XMLS'
      CHREF(206) = 'TA'
      UNREF(206) = 'k'
      UNPLOT(206) = 'k'
      CHREF(207) = 'RH'
      CHREF(208) = 'RN'
      UNREF(208) = 'w/m^2'
      UNPLOT(208) = 'w/m^2'
      CHREF(209) = 'WS'
      UNREF(209) = 'm/s'
      UNPLOT(209) = 'm/s'
      CHREF(210) = 'SNR'
      CHREF(211) = 'SNM'
      CHREF(212) = 'SNF'
      CHREF(213) = 'T_S'
      UNREF(213) = 'k'
      CHREF(214) = 'PV_S'
      UNREF(214) = 'pa'
      CHREF(215) = 'E_SA'
      UNREF(215) = 'kg/s'
      CHREF(216) = 'PESA'
      UNREF(216) = 'kg/s'
      CHREF(217) = 'T_PA'
      UNREF(217) = 'kg/s'
      UNPLOT(217) = 'kg/s'
      CHREF(218) = 'PTPA'
      UNREF(218) = 'kg/s'
      UNPLOT(218) = 'kg/s'
      CHREF(219) = 'SXLA'
      CHREF(220) = 'XMLA'
      CHREF(221) = 'XMNA'
      CHREF(222) = 'XLA'
      CHREF(223) = 'XNA'
      CHREF(224) = 'PA'
      UNREF(224) = 'pa'
      UNPLOT(224) = 'pa'
      CHREF(225) = 'PL_S'
      UNREF(225) = 'pa'
      UNPLOT(225) = 'pa'
      CHREF(226) = 'PG_S'
      UNREF(226) = 'pa'
      UNPLOT(226) = 'pa'
      CHREF(227) = 'SL_S'
      CHREF(228) = 'QL_S'
      UNREF(228) = 'w/m^2'
      CHREF(229) = 'QH_S'
      UNREF(229) = 'w/m^2'
      CHREF(230) = 'RL_S'
      UNREF(230) = 'w/m^2'
      CHREF(231) = 'RS_S'
      UNREF(231) = 'w/m^2'
      CHREF(232) = 'QG_S'
      UNREF(232) = 'w/m^2'
      CHREF(233) = 'WB_S'
      UNREF(233) = 'kg/s'
      CHREF(234) = 'T_P'
      UNREF(234) = 'k'
      CHREF(235) = 'TP'
      UNREF(235) = 'k'
      CHREF(236) = 'TP'
      UNREF(236) = 'k'
      CHREF(237) = 'TP'
      UNREF(237) = 'k'
      CHREF(238) = 'TP'
      UNREF(238) = 'k'
      CHREF(239) = 'RFIM'
      UNREF(239) = 'kg'
      CHREF(240) = 'TSO'
      UNREF(240) = 'kg'
      UNPLOT(240) = 'kg'
      CHREF(241) = 'XSO'
      CHREF(242) = 'CSO'
      UNREF(242) = 'kg/m^3'
      UNPLOT(242) = 'kg/m^3'
      CHREF(243) = 'RABS'
      UNREF(243) = 's/m'
      UNPLOT(243) = 's/m'
      CHREF(244) = 'PV_S'
      UNREF(244) = 'm^3/s'
      UNPLOT(244) = 'm^3/s'
      CHREF(245) = 'PM_S'
      UNREF(245) = 'kg/s'
      UNPLOT(245) = 'kg/s'
      CHREF(246) = 'PM_A'
      UNREF(246) = 'kg/s'
      UNPLOT(246) = 'kg/s'
      CHREF(247) = 'UK'
      UNREF(247) = 'm^2'
      UNPLOT(247) = 'm^2'
      CHREF(248) = 'VK'
      UNREF(248) = 'm^2'
      UNPLOT(248) = 'm^2'
      CHREF(249) = 'WK'
      UNREF(249) = 'm^2'
      UNPLOT(249) = 'm^2'
      CHREF(250) = 'XHW'
      CHREF(251) = 'XHA'
      CHREF(252) = 'XHO'
      CHREF(253) = 'RHOH'
      UNREF(253) = 'kg/m^3'
      UNPLOT(253) = 'kg/m^3'
      CHREF(254) = 'SH'
      CHREF(255) = 'PH'
      UNREF(255) = 'Pa'
      UNPLOT(255) = 'Pa'
      CHREF(256) = 'IMHW'
      UNREF(256) = 'kg'
      UNPLOT(256) = 'kg'
      CHREF(257) = 'IMHA'
      UNREF(257) = 'kg'
      UNPLOT(257) = 'kg'
      CHREF(258) = 'IMHO'
      UNREF(258) = 'kg'
      UNPLOT(258) = 'kg'
      CHREF(259) = 'IMLO'
      UNREF(259) = 'kg'
      UNPLOT(259) = 'kg'
      CHREF(260) = 'IMGO'
      UNREF(260) = 'kg'
      UNPLOT(260) = 'kg'
      CHREF(261) = 'IMSW'
      UNREF(261) = 'kg'
      UNPLOT(261) = 'kg'
      CHREF(262) = 'IMSA'
      UNREF(262) = 'kg'
      UNPLOT(262) = 'kg'
      CHREF(263) = 'IMSO'
      UNREF(263) = 'kg'
      UNPLOT(263) = 'kg'
      CHREF(264) = 'SS'
      CHREF(265) = 'XMHW'
      CHREF(266) = 'XMHA'
      CHREF(267) = 'XMHO'
      CHREF(274) = 'RWRO'
      UNREF(274) = 'm^3/s'
      UNPLOT(274) = 'm^3/s'
      CHREF(275) = 'WNTP'
      UNREF(275) = 'Pa'
      UNPLOT(275) = 'Pa'
      CHREF(276) = 'WNBP'
      UNREF(276) = 'Pa'
      UNPLOT(276) = 'Pa'
      CHREF(277) = 'DMHW'
      UNREF(277) = 'kg'
      UNPLOT(277) = 'kg'
      CHREF(278) = 'DMHA'
      UNREF(278) = 'kg'
      UNPLOT(278) = 'kg'
      CHREF(279) = 'DMHO'
      UNREF(279) = 'kg'
      UNPLOT(279) = 'kg'
      CHREF(280) = 'DMLO'
      UNREF(280) = 'kg'
      UNPLOT(280) = 'kg'
      CHREF(281) = 'DMGO'
      UNREF(281) = 'kg'
      UNPLOT(281) = 'kg'
      CHREF(282) = 'DMW'
      UNREF(282) = 'kg'
      UNPLOT(282) = 'kg'
      CHREF(283) = 'DMA'
      UNREF(283) = 'kg'
      UNPLOT(283) = 'kg'
      CHREF(284) = 'DMO'
      UNREF(284) = 'kg'
      UNPLOT(284) = 'kg'
      CHREF(285) = 'SWP'
      UNREF(285) = 'Pa'
      UNPLOT(285) = 'Pa'
      CHREF(286) = 'SWT'
      UNREF(286) = 'C'
      UNPLOT(286) = 'C'
      CHREF(287) = 'HLSP'
      UNREF(287) = 'm'
      UNPLOT(287) = 'm'
      CHREF(288) = 'HNSP'
      UNREF(288) = 'm'
      UNPLOT(288) = 'm'
      CHREF(289) = 'VISG'
      UNPLOT(289) = 'pa s'
      UNREF(289) = 'pa s'
      CHREF(290) = 'VISN'
      UNPLOT(290) = 'pa s'
      UNREF(290) = 'pa s'
      CHREF(291) = 'XP'
      UNPLOT(291) = 'm'
      UNREF(291) = 'm'
      CHREF(292) = 'YP'
      UNPLOT(292) = 'm'
      UNREF(292) = 'm'
      CHREF(293) = 'ZP'
      UNPLOT(293) = 'm'
      UNREF(293) = 'm'
      CHREF(294) = 'YMGO'
      CHREF(295) = 'YMHGO'
      CHREF(296) = 'HL'
      UNPLOT(296) = 'j/kg'
      UNREF(296) = 'j/kg'
      CHREF(297) = 'HG'
      UNPLOT(297) = 'j/kg'
      UNREF(297) = 'j/kg'
      CHREF(298) = 'HN'
      UNPLOT(298) = 'j/kg'
      UNREF(298) = 'j/kg'
      CHREF(299) = 'SIMV'
      UNPLOT(299) = 'm^2/s'
      UNREF(299) = 'm^2/s'
      CHREF(300) = 'TAI'
      UNREF(300) = 'C'
      UNPLOT(300) = 'C'
      CHREF(301) = 'TAO'
      UNREF(301) = 'C'
      UNPLOT(301) = 'C'
      CHREF(302) = 'TFI'
      UNREF(302) = 'C'
      UNPLOT(302) = 'C'
      CHREF(303) = 'TFO'
      UNREF(303) = 'C'
      UNPLOT(303) = 'C'
      CHREF(304) = 'TAIW'
      UNREF(304) = 'C'
      UNPLOT(304) = 'C'
      CHREF(305) = 'TAOW'
      UNREF(305) = 'C'
      UNPLOT(305) = 'C'
      CHREF(306) = 'TFIW'
      UNREF(306) = 'C'
      UNPLOT(306) = 'C'
      CHREF(307) = 'TFOW'
      UNREF(307) = 'C'
      UNPLOT(307) = 'C'
      CHREF(308) = 'TFC'
      UNREF(308) = 'C'
      UNPLOT(308) = 'C'
      CHREF(309) = 'TSF'
      UNREF(309) = 'C'
      UNPLOT(309) = 'C'
      CHREF(310) = 'TIC'
      UNREF(310) = 'C'
      UNPLOT(310) = 'C'
      CHREF(311) = 'TIN'
      UNREF(311) = 'C'
      UNPLOT(311) = 'C'
      CHREF(312) = 'TOC'
      UNREF(312) = 'C'
      UNPLOT(312) = 'C'
      CHREF(313) = 'PAI'
      UNREF(313) = 'Pa'
      UNPLOT(313) = 'Pa'
      CHREF(314) = 'PAO'
      UNREF(314) = 'Pa'
      UNPLOT(314) = 'Pa'
      CHREF(315) = 'PFI'
      UNREF(315) = 'Pa'
      UNPLOT(315) = 'Pa'
      CHREF(316) = 'PFO'
      UNREF(316) = 'Pa'
      UNPLOT(316) = 'Pa'
      CHREF(317) = 'WAI'
      UNREF(317) = 'm^3/s'
      UNPLOT(317) = 'm^3/s'
      CHREF(318) = 'WAS'
      UNREF(318) = 'm^3/s'
      UNPLOT(318) = 'm^3/s'
      CHREF(319) = 'WAO'
      UNREF(319) = 'm^3/s'
      UNPLOT(319) = 'm^3/s'
      CHREF(320) = 'WFI'
      UNREF(320) = 'm^3/s'
      UNPLOT(320) = 'm^3/s'
      CHREF(321) = 'WFS'
      UNREF(321) = 'm^3/s'
      UNPLOT(321) = 'm^3/s'
      CHREF(322) = 'WFO'
      UNREF(322) = 'm^3/s'
      UNPLOT(322) = 'm^3/s'
      CHREF(323) = 'MAI'
      UNREF(323) = 'kg/s'
      UNPLOT(323) = 'kg/s'
      CHREF(324) = 'MAS'
      UNREF(324) = 'kg/s'
      UNPLOT(324) = 'kg/s'
      CHREF(325) = 'MAO'
      UNREF(325) = 'kg/s'
      UNPLOT(325) = 'kg/s'
      CHREF(326) = 'MFI'
      UNREF(326) = 'kg/s'
      UNPLOT(326) = 'kg/s'
      CHREF(327) = 'MFS'
      UNREF(327) = 'kg/s'
      UNPLOT(327) = 'kg/s'
      CHREF(328) = 'MFO'
      UNREF(328) = 'kg/s'
      UNPLOT(328) = 'kg/s'
      CHREF(329) = 'QGM'
      UNREF(329) = 'w'
      UNPLOT(329) = 'w'
      CHREF(330)= 'QSPH'
      UNREF(330) = 'w'
      UNPLOT(330) = 'w'
      CHREF(331) = 'QFCP'
      UNREF(331) = 'w'
      UNPLOT(331) = 'w'
      CHREF(332) = 'QNAH'
      UNREF(332) = 'w'
      UNPLOT(332) = 'w'
      CHREF(333) = 'QNFH'
      UNREF(333) = 'w'
      UNPLOT(333) = 'w'
      CHREF(334) = 'V'
      CHREF(335) = 'LCD'
      CHREF(341) = 'UEL'
      UNPLOT(341) = 'j/kg'
      UNREF(341) = 'j/kg'
      CHREF(342) = 'UEG'
      UNPLOT(342) = 'j/kg'
      UNREF(342) = 'j/kg'
      CHREF(343) = 'UEN'
      UNPLOT(343) = 'j/kg'
      UNREF(343) = 'j/kg'
      CHREF(344) = 'THKL'
      UNPLOT(344) = 'w/m K'
      UNREF(344) = 'w/m K'
      CHREF(345) = 'THKG'
      UNPLOT(345) = 'w/m K'
      UNREF(345) = 'w/m K'
      CHREF(346) = 'THKN'
      UNPLOT(346) = 'w/m K'
      UNREF(346) = 'w/m K'
      CHREF(347) = 'DFLA'
      UNPLOT(347) = 'm^2/s'
      UNREF(347) = 'm^2/s'
      CHREF(348) = 'DFGW'
      UNPLOT(348) = 'm^2/s'
      UNREF(348) = 'm^2/s'
      CHREF(349) = 'QMRA'
      UNPLOT(349) = 'kg/s'
      UNREF(349) = 'kg/s'
      CHREF(350) = 'QMIA'
      UNPLOT(350) = 'kg'
      UNREF(350) = 'kg'
      CHREF(351) = 'QMRW'
      UNPLOT(351) = 'kg/s'
      UNREF(351) = 'kg/s'
      CHREF(352) = 'QMIW'
      UNPLOT(352) = 'kg'
      UNREF(352) = 'kg'
      CHREF(361) = 'VIA'
      UNPLOT(361) = 'kg'
      UNREF(361) = 'kg'
      CHREF(362) = 'VIAPA'
      UNPLOT(362) = 'kg/m^2'
      UNREF(362) = 'kg/m^2'
      CHREF(363) = 'VIGA'
      UNPLOT(363) = 'kg'
      UNREF(363) = 'kg'
      CHREF(364) = 'VIGAPA'
      UNPLOT(364) = 'kg/m^2'
      UNREF(364) = 'kg/m^2'
      CHREF(365) = 'VILA'
      UNPLOT(365) = 'kg'
      UNREF(365) = 'kg'
      CHREF(366) = 'VILAPA'
      UNPLOT(366) = 'kg/m^2'
      UNREF(366) = 'kg/m^2'
      CHREF(367) = 'IMPS'
      UNPLOT(367) = 'kg'
      UNREF(367) = 'kg'
      CHREF(368) = 'STRS-M'
      UNPLOT(368) = 'Pa'
      UNREF(368) = 'Pa'
      CHREF(369) = 'EPS-XX'
      CHREF(370) = 'EPS-YY'
      CHREF(371) = 'EPS-ZZ'
      CHREF(372) = 'EPS-YZ'
      CHREF(373) = 'EPS-XZ'
      CHREF(374) = 'EPS-XY'
      CHREF(375) = 'DISPX'
      UNPLOT(375) = 'm'
      UNREF(375) = 'm'
      CHREF(376) = 'DISPY'
      UNPLOT(376) = 'm'
      UNREF(376) = 'm'
      CHREF(377) = 'DISPZ'
      UNPLOT(377) = 'm'
      UNREF(377) = 'm'
      CHREF(378) = 'IQ'
      UNPLOT(378) = 'J'
      UNREF(378) = 'J'
      CHREF(379) = 'YMGA'
      CHREF(380) = 'YMHGA'
      CHREF(381) = 'YMGN'
      CHREF(382) = 'YMHGN'
      CHREF(383) = 'ZMCA'
      CHREF(384) = 'ZMCO'
      CHREF(385) = 'ZMCN'
      CHREF(386) = 'PEQH'
      UNPLOT(386) = 'pa'
      UNREF(386) = 'pa'
      CHREF(387) = 'PVTF'
      UNPLOT(387) = 'pa'
      UNREF(387) = 'pa'
      CHREF(388) = 'I'
      CHREF(389) = 'J'
      CHREF(390) = 'K'
      CHREF(391) = 'AFX'
      UNPLOT(391) = 'm^2'
      UNREF(391) = 'm^2'
      CHREF(392) = 'AFY'
      UNPLOT(392) = 'm^2'
      UNREF(392) = 'm^2'
      CHREF(393) = 'AFZ'
      UNPLOT(393) = 'm^2'
      UNREF(393) = 'm^2'
      CHREF(394) = 'SFZN'
      CHREF(395) = 'PID'
      CHREF(396) = 'IFT'
      UNPLOT(396) = 'n/m'
      UNREF(396) = 'n/m'
      CHREF(397) = 'APM'
      UNPLOT(397) = 'm'
      UNREF(397) = 'm'
      CHREF(398) = 'APH'
      UNPLOT(398) = 'm'
      UNREF(398) = 'm'
      CHREF(399) = 'DABH'
      UNPLOT(399) = 'm'
      UNREF(399) = 'm'
      CHSF(1,1) = 'UQV'
      CHSF(1,2) = 'VQV'
      CHSF(1,3) = 'WQV'
      CHSF(2,1) = 'ULV'
      CHSF(2,2) = 'VLV'
      CHSF(2,3) = 'WLV'
      CHSF(3,1) = 'UGV'
      CHSF(3,2) = 'VGV'
      CHSF(3,3) = 'WGV'
      CHSF(4,1) = 'UNV'
      CHSF(4,2) = 'VNV'
      CHSF(4,3) = 'WNV'
      CHSF(5,1) = 'ULM'
      CHSF(5,2) = 'VLM'
      CHSF(5,3) = 'WLM'
      CHSF(6,1) = 'UGM'
      CHSF(6,2) = 'VGM'
      CHSF(6,3) = 'WGM'
      CHSF(7,1) = 'UNM'
      CHSF(7,2) = 'VNM'
      CHSF(7,3) = 'WNM'
      CHSF(8,1) = 'USM'
      CHSF(8,2) = 'VSM'
      CHSF(8,3) = 'WSM'
      CHSF(9,1) = 'ULO'
      CHSF(9,2) = 'VLO'
      CHSF(9,3) = 'WLO'
      CHSF(10,1) = 'UWM'
      CHSF(10,2) = 'VWM'
      CHSF(10,3) = 'WWM'
      CHSF(11,1) = 'UGOM'
      CHSF(11,2) = 'VGOM'
      CHSF(11,3) = 'WGOM'
      CHSF(12,1) = 'ULOM'
      CHSF(12,2) = 'VLOM'
      CHSF(12,3) = 'WLOM'
      CHSF(13,1) = 'UOM'
      CHSF(13,2) = 'VOM'
      CHSF(13,3) = 'WOM'
      CHSF(20,1) = 'UGAQ'
      CHSF(20,2) = 'VGAQ'
      CHSF(20,3) = 'WGAQ'
      CHSF(21,1) = 'UGAW'
      CHSF(21,2) = 'VGAW'
      CHSF(21,3) = 'WGAW'
      CHSF(22,1) = 'UGAA'
      CHSF(22,2) = 'VGAA'
      CHSF(22,3) = 'WGAA'
      CHSF(25,1) = 'UGDQ'
      CHSF(25,2) = 'VGDQ'
      CHSF(25,3) = 'WGDQ'
      CHSF(26,1) = 'UGDW'
      CHSF(26,2) = 'VGDW'
      CHSF(26,3) = 'WGDW'
      CHSF(27,1) = 'UGDA'
      CHSF(27,2) = 'VGDA'
      CHSF(27,3) = 'WGDA'
      CHSF(28,1) = 'UGAM'
      CHSF(28,2) = 'VGAM'
      CHSF(28,3) = 'WGAM'
      CHSF(29,1) = 'ULAM'
      CHSF(29,2) = 'VLAM'
      CHSF(29,3) = 'WLAM'
      CHSF(30,1) = 'UAM'
      CHSF(30,2) = 'VAM'
      CHSF(30,3) = 'WAM'
      CHSF(31,1) = 'UGAO'
      CHSF(31,2) = 'VGAO'
      CHSF(31,3) = 'WGAO'
      CHSF(32,1) = 'UGDO'
      CHSF(32,2) = 'VGDO'
      CHSF(32,3) = 'WGDO'
      CHSF(33,1) = 'UGO'
      CHSF(33,2) = 'VGO'
      CHSF(33,3) = 'WGO'
      CHSF(34,1) = 'AE'
      CHSF(34,2) = 'AE'
      CHSF(34,3) = 'AE'
      CHSF(35,1) = 'PE'
      CHSF(35,2) = 'PE'
      CHSF(35,3) = 'PE'
      CHSF(36,1) = 'AT'
      CHSF(36,2) = 'AT'
      CHSF(36,3) = 'AT'
      CHSF(37,1) = 'PT'
      CHSF(37,2) = 'PT'
      CHSF(37,3) = 'PT'
      CHSF(38,1) = 'NTR'
      CHSF(38,2) = 'NTR'
      CHSF(38,3) = 'NTR'
      CHSF(39,1) = 'NSWR'
      CHSF(39,2) = 'NSWR'
      CHSF(39,3) = 'NSWR'
      CHSF(40,1) = 'NLWR'
      CHSF(40,2) = 'NLWR'
      CHSF(40,3) = 'NLWR'
      CHSF(41,1) = 'WMB'
      CHSF(41,2) = 'WMB'
      CHSF(41,3) = 'WMB'
      CHSF(42,1) = 'RWRO'
      CHSF(42,2) = 'RWRO'
      CHSF(42,3) = 'RWRO'
      CHSF(43,1) = 'ULWM'
      CHSF(43,2) = 'VLWM'
      CHSF(43,3) = 'WLWM'
      CHSF(44,1) = 'UGWM'
      CHSF(44,2) = 'VGWM'
      CHSF(44,3) = 'WGWM'
      CHSF(45,1) = 'UWM'
      CHSF(45,2) = 'VWM'
      CHSF(45,3) = 'WWM'
      CHSF(46,1) = 'ULAC'
      CHSF(46,2) = 'VLAC'
      CHSF(46,3) = 'WLAC'
      CHSF(47,1) = 'ULDC'
      CHSF(47,2) = 'VLDC'
      CHSF(47,3) = 'WLDC'
      CHSF(48,1) = 'UGAC'
      CHSF(48,2) = 'VGAC'
      CHSF(48,3) = 'WGAC'
      CHSF(49,1) = 'UGDC'
      CHSF(49,2) = 'VGDC'
      CHSF(49,3) = 'WGDC'
      CHSF(50,1) = 'UTAC'
      CHSF(50,2) = 'VTAC'
      CHSF(50,3) = 'WTAC'
      CHSF(51,1) = 'UTDC'
      CHSF(51,2) = 'VTDC'
      CHSF(51,3) = 'WTDC'
      CHSF(52,1) = 'UNAC'
      CHSF(52,2) = 'VNAC'
      CHSF(52,3) = 'WNAC'
      CHSF(53,1) = 'UNDC'
      CHSF(53,2) = 'VNDC'
      CHSF(53,3) = 'WNDC'
      CHSF(54,1) = 'UTC'
      CHSF(54,2) = 'VTC'
      CHSF(54,3) = 'WTC'
      CHSF(55,1) = 'UNAM'
      CHSF(55,2) = 'VNAM'
      CHSF(55,3) = 'WNAM'
      CHSF(56,1) = 'UNOM'
      CHSF(56,2) = 'VNOM'
      CHSF(56,3) = 'WNOM'
      CHSF(57,1) = 'UNWM'
      CHSF(57,2) = 'VNWM'
      CHSF(57,3) = 'WNWM'
      CHSF(58,1) = 'UN2M'
      CHSF(58,2) = 'VN2M'
      CHSF(58,3) = 'WN2M'
      CHSF(59,1) = 'ULNM'
      CHSF(59,2) = 'VLNM'
      CHSF(59,3) = 'WLNM'
      CHSF(60,1) = 'UGNM'
      CHSF(60,2) = 'VGNM'
      CHSF(60,3) = 'WGNM'
      CHSF(61,1) = 'UNNM'
      CHSF(61,2) = 'VNNM'
      CHSF(61,3) = 'WNNM'
      CHSF(62,1) = 'UGC'
      CHSF(62,2) = 'VGC'
      CHSF(62,3) = 'WGC'
      CHSF(63,1) = 'ULC'
      CHSF(63,2) = 'VLC'
      CHSF(63,3) = 'WLC'
      CHSF(64,1) = 'UNC'
      CHSF(64,2) = 'VNC'
      CHSF(64,3) = 'WNC'
      CHSF(65,1) = 'UTC'
      CHSF(65,2) = 'VTC'
      CHSF(65,3) = 'WTC'
      CHSF(66,1) = 'UFML'
      CHSF(66,2) = 'UFML'
      CHSF(66,3) = 'UFML'
      CHSF(67,1) = 'UBBL'
      CHSF(67,2) = 'UBBL'
      CHSF(67,3) = 'UBBL'
      CHSF(100,1) = 'UCM'
      CHSF(100,2) = 'VCM'
      CHSF(100,3) = 'WCM'
      CHSF(101,1) = 'UCNM'
      CHSF(101,2) = 'VCNM'
      CHSF(101,3) = 'WCNM'
      CHSF(102,1) = 'UCC'
      CHSF(102,2) = 'VCC'
      CHSF(102,3) = 'WCC'
      CHSF(103,1) = 'UKC'
      CHSF(103,2) = 'VKC'
      CHSF(103,3) = 'WKC'
      DO 70 L = 1,((LSOLU**LC)+((LSPR+LSPT)**LR)+2-LC-LR)
        M = 400 + (L-1)*33
        WRITE(FORM1(3:3),'(I1)') ICOUNT(L)
        CHREF(M+1) = 'C     '
        UNPLOT(M+1) = '1/m^3'
        UNREF(M+1) = '1/m^3'
        WRITE( CHREF(M+1)(2:),FORM1) L
        CHREF(M+2) = 'CL    '
        UNPLOT(M+2) = '1/m^3'
        UNREF(M+2) = '1/m^3'
        WRITE( CHREF(M+2)(3:),FORM1) L
        CHREF(M+3) = 'CG    '
        UNPLOT(M+3) = '1/m^3'
        UNREF(M+3) = '1/m^3'
        WRITE( CHREF(M+3)(3:),FORM1) L
        CHREF(M+4) = 'CN    '
        UNPLOT(M+4) = '1/m^3'
        UNREF(M+4) = '1/m^3'
        WRITE( CHREF(M+4)(3:),FORM1) L
        CHREF(M+5) = 'YL    '
        WRITE( CHREF(M+5)(3:),FORM1) L
        CHREF(M+6) = 'YG    '
        WRITE( CHREF(M+6)(3:),FORM1) L
        CHREF(M+7) = 'YN    '
        WRITE( CHREF(M+7)(3:),FORM1) L
        CHREF(M+8) = 'UC    '
        UNPLOT(M+8) = '1/m^2 s'
        UNREF(M+8) = '1/m^2 s'
        WRITE( CHREF(M+8)(3:),FORM1) L
        CHREF(M+9) = 'VC    '
        UNPLOT(M+9) = '1/m^2 s'
        UNREF(M+9) = '1/m^2 s'
        WRITE( CHREF(M+9)(3:),FORM1) L
        CHREF(M+10) = 'WC    '
        UNPLOT(M+10) = '1/m^2 s'
        UNREF(M+10) = '1/m^2 s'
        WRITE( CHREF(M+10)(3:),FORM1) L
        CHREF(M+11) = 'SRC   '
        UNPLOT(M+11) = '1/s'
        UNREF(M+11) = '1/s'
        WRITE( CHREF(M+11)(4:),FORM1) L
        CHREF(M+12) = 'C     '
        UNPLOT(M+12) = '1/m^3'
        UNREF(M+12) = '1/m^3'
        WRITE( CHREF(M+12)(2:),FORM1) L
        CHREF(M+13) = 'CL    '
        UNPLOT(M+13) = '1/m^3'
        UNREF(M+13) = '1/m^3'
        WRITE( CHREF(M+13)(3:),FORM1) L
        CHREF(M+14) = 'CN    '
        UNPLOT(M+14) = '1/m^3'
        UNREF(M+14) = '1/m^3'
        WRITE( CHREF(M+14)(3:),FORM1) L
        CHREF(M+15) = 'CS    '
        UNPLOT(M+15) = '1/m^3'
        UNREF(M+15) = '1/m^3'
        WRITE( CHREF(M+15)(3:),FORM1) L
        CHREF(M+16) = 'UCL   '
        UNPLOT(M+16) = '1/m^2 s'
        UNREF(M+16) = '1/m^2 s'
        WRITE( CHREF(M+16)(4:),FORM1) L
        CHREF(M+17) = 'VCL   '
        UNPLOT(M+17) = '1/m^2 s'
        UNREF(M+17) = '1/m^2 s'
        WRITE( CHREF(M+17)(4:),FORM1) L
        CHREF(M+18) = 'WCL   '
        UNPLOT(M+18) = '1/m^2 s'
        UNREF(M+18) = '1/m^2 s'
        WRITE( CHREF(M+18)(4:),FORM1) L
        CHREF(M+19) = 'UCN   '
        UNPLOT(M+19) = '1/m^2 s'
        UNREF(M+19) = '1/m^2 s'
        WRITE( CHREF(M+19)(4:),FORM1) L
        CHREF(M+20) = 'VCN   '
        UNPLOT(M+20) = '1/m^2 s'
        UNREF(M+20) = '1/m^2 s'
        WRITE( CHREF(M+20)(4:),FORM1) L
        CHREF(M+21) = 'WCN   '
        UNPLOT(M+21) = '1/m^2 s'
        UNREF(M+21) = '1/m^2 s'
        WRITE( CHREF(M+21)(4:),FORM1) L
        CHREF(M+22) = 'FCB   '
        UNPLOT(M+22) = '1/s'
        UNREF(M+22) = '1/s'
        WRITE( CHREF(M+22)(4:),FORM1) L
        CHREF(M+23) = 'ICM   '
        WRITE( CHREF(M+23)(4:),FORM1) L
        CHREF(M+30) = 'CM   '
        WRITE( CHREF(M+30)(4:),FORM1) L
        UNPLOT(M+30) = '1/kg'
        UNREF(M+30) = '1/kg'
        CHREF(M+31) = 'CA    '
        UNPLOT(M+31) = 'Bq'
        UNREF(M+31) = 'Bq'
        CHREF(M+32) = 'CIA   '
        UNPLOT(M+32) = 'Bq'
        UNREF(M+32) = 'Bq'
   70 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_OUTP group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_WELL
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
!     Initialize well variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 19 June 2001.
!     Last Modified by MD White, PNNL, 19 June 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
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
      SUB_LOG(ISUB_LOG) = '/IN_WELL'
!
!---  Well state variables  ---
!
      DO 100 L = 1,LNWN
        DO 10 K = 1,LSV
          PW(K,L) = 0.D+0
          PWLW(K,L) = 0.D+0
          PWGW(K,L) = 0.D+0
          PWNW(K,L) = 0.D+0
          RHOLW(K,L) = 0.D+0
          RHOGW(K,L) = 0.D+0
          RHONW(K,L) = 0.D+0
          SLW(K,L) = 0.D+0
          SNW(K,L) = 0.D+0
          STW(K,L) = 0.D+0
          XLWW(K,L) = 0.D+0
          XLAW(K,L) = 0.D+0
          XLOW(K,L) = 0.D+0
          XGWW(K,L) = 0.D+0
          XGAW(K,L) = 0.D+0
          XMLAW(K,L) = 0.D+0
          XMLOW(K,L) = 0.D+0
          XMGWW(K,L) = 0.D+0
          VISLW(K,L) = 0.D+0
          VISGW(K,L) = 0.D+0
          VISNW(K,L) = 0.D+0
          DFLAW(K,L) = 0.D+0
          DFLOW(K,L) = 0.D+0
          DFGWW(K,L) = 0.D+0
          YMLOW(K,L) = 0.D+0
          SRCW_W(K,L) = 0.D+0
          SRCO_W(K,L) = 0.D+0
   10   CONTINUE
!
!---    Well phase condition  ---
!
        DO 20 K = 1,LSU
          NPHAZW(K,L) = 2
   20   CONTINUE
!
!---    Newton-Raphson increments  ---
!
        DO 30 K = 1,LUK
          DNRW(K,L) = 0.D+0
   30   CONTINUE
!
!---    Second-order time differencing state variables  ---
!
        RHOLW_O(L) = 0.D+0
        RHOGW_O(L) = 0.D+0
        RHONW_O(L) = 0.D+0
        SLW_O(L) = 0.D+0
        SNW_O(L) = 0.D+0
        STW_O(L) = 0.D+0
        XLWW_O(L) = 0.D+0
        XLOW_O(L) = 0.D+0
        XGWW_O(L) = 0.D+0
        XLAW_O(L) = 0.D+0
        XGAW_O(L) = 0.D+0
!
!---    Integer pointers  ---
!
        IWL(L) = 0
        IWN(L) = 0
        NSZW(L) = 0
  100 CONTINUE
!
!---  Well vertical fluxes  ---
!
      DO 200 L = 1,LSZW
        DO 110 K = 1,LSFV
          WL_W(K,L) = 0.D+0
          WG_W(K,L) = 0.D+0
          WN_W(K,L) = 0.D+0
          WDLA_W(K,L) = 0.D+0
          WDLO_W(K,L) = 0.D+0
          WDGW_W(K,L) = 0.D+0
  110   CONTINUE
  200 CONTINUE
!
!---  Well cross-screen fluxes  ---
!
      DO 300 L = 1,LNWN
        DO 210 K = 1,LSFV
          UL_W(K,L) = 0.D+0
          UG_W(K,L) = 0.D+0
          UN_W(K,L) = 0.D+0
          UDLA_W(K,L) = 0.D+0
          UDLO_W(K,L) = 0.D+0
          UDGW_W(K,L) = 0.D+0
  210   CONTINUE
  300 CONTINUE
!
!---  Well boundary conditions  ---
!
      DO 400 L = 1,LNW
        DO 320 K = 1,LNWT
          DO 310 J = 1,LNWV
            WLVR(J,K,L) = 0.D+0
  310     CONTINUE
  320   CONTINUE
!
!---    Well boundary conditions  ---
!
        WBR(L) = 0.D+0
        WBRS(L) = 0.D+0
        WWD(L) = 0.D+0
        WWDL(L) = 0.D+0
        WWDN(L) = 0.D+0
        WHP(L) = 0.D+0
        WRH(L) = 0.D+0
        WIDA(L) = 0.D+0
        WIDO(L) = 0.D+0
        IWM(L) = 0
        IWT(L) = 0
        IWCC(L) = 0
!
!---    Well domain indices  ---
!
        DO 330 K = 1,2*LNWS+4
          IWLDM(K,L) = 0
  330   CONTINUE
        DO 340 K = 1,4
          QL_W(K,L) = 0.D+0
          QN_W(K,L) = 0.D+0
          QT_W(K,L) = 0.D+0
  340   CONTINUE
  400 CONTINUE
!
!---  Well node type index  ---
!
      DO 500 L = 1,LFD
        IXW(L) = 0
  500 CONTINUE
!
!---  Well counters  ---
!
      NWLN = 0
      NWLS = 0
!
!---  Well boundary condition variables  ---
!
      DO 600 L = 1,LNW
        DO 510 K = 1,LSV
          PWB(K,L) = 0.D+0
          RHOLWB(K,L) = 0.D+0
          RHONWB(K,L) = 0.D+0
          RHOGWB(K,L) = 0.D+0
          XLWWB(K,L) = 0.D+0
          XLOWB(K,L) = 0.D+0
          XGWWB(K,L) = 0.D+0
          XLAWB(K,L) = 0.D+0
          XGAWB(K,L) = 0.D+0
          XMLAWB(K,L) = 0.D+0
          XMLOWB(K,L) = 0.D+0
          XMGWWB(K,L) = 0.D+0
          VISLWB(K,L) = 0.D+0
          VISGWB(K,L) = 0.D+0
          VISNWB(K,L) = 0.D+0
          DFLOWB(K,L) = 0.D+0
          DFLAWB(K,L) = 0.D+0
          DFGWWB(K,L) = 0.D+0
  510   CONTINUE
  600 CONTINUE
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_WELL group  ---
!
      RETURN
      END
      
!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_COUP_WELL
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
!     Initialize CO2 well variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 January 2011.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_COUP_WELL'
!
!---  Coupled well variables  ---
!
      DO 100 L = 1,LN_CW
        DNR_CW(L) = 0.D+0
        FX_CW(L) = 0.D+0
        TML_CW(L) = 0.D+0
        PL_CW(L) = 0.D+0
        RHOF_CW(L) = 0.D+0
        WNM_CW(L) = 'null'
        T_CW(1,L) = 0.D+0
        T_CW(2,L) = 0.D+0
        DO K = 1,LWT_CW
          DO J = 1,7+LNGC
            VAR_CW(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,LSOLU_CW
            VARC_CW(J,K,L) = 0.D+0
          ENDDO
          DO J = 1,LSPC_CW
            VARSP_CW(J,K,L) = 0.D+0
          ENDDO
        ENDDO
        DO 30 K = 1,LWTP_CW
          IMP_CW(K,L) = 0
          ITS_CW(K,L) = 0
   30   CONTINUE
        ICC_CW(L) = 0
        DO 40 K = 1,10
          ID_CW(K,L) = 0
   40   CONTINUE
        DO 50 K = 1,3
          FF_CW(K,L) = 0.D+0
          P_CW(K,L) = -1.D+20
   50   CONTINUE
        DO 60 K = 1,(LUK_CW+1)
          RS_CW(K,L) = 0.D+0
   60   CONTINUE
        IM_CW(L) = 0
        IT_CW(L) = 0
        JM_CW(L) = 0
        NSP_CW(L) = 0
        DO 70 K = 1,(2*(3+LNGC))
          QM_CW(K,L) = 0.D+0
   70   CONTINUE
  100 CONTINUE
      DO 200 L = 1,LWI_CW
        DO 110 K = 1,2
          XTP_CW(K,L) = 0.D+0
          YTP_CW(K,L) = 0.D+0
          ZTP_CW(K,L) = 0.D+0
  110   CONTINUE
        DO 112 K = 1,5
          PAR_CW(K,L) = 0.D+0
  112   CONTINUE
        IS_CW(L) = 0
  200 CONTINUE
      DO 300 L = 1,LWN_CW
        DO 210 K = 1,2
          XP_CW(K,L) = 0.D+0
          YP_CW(K,L) = 0.D+0
          ZP_CW(K,L) = 0.D+0
  210   CONTINUE
        DO 230 K = 1,(LUK+2)
          FXA_CW(K,L) = 0.D+0
          FXO_CW(K,L) = 0.D+0
          FXE_CW(K,L) = 0.D+0
          FXS_CW(K,L) = 0.D+0
          FXW_CW(K,L) = 0.D+0
          DO 220 J = 1,LNGC
            FXC_CW(J,K,L) = 0.D+0
  220     CONTINUE
  230   CONTINUE
        DO K = 1,4
          Q_CW(K,L) = 0.D+0
        ENDDO
        PLX_CW(L) = 0.D+0
        PLY_CW(L) = 0.D+0
        PLZ_CW(L) = 0.D+0
        IWN_CW(L) = 0
        IWP_CW(L) = 0
        INV_CW(L) = 0
  300 CONTINUE
      DO 320 K = 1,LN_CW
        DO J = 1,LSOLU_CW
          ISOLU_CW(J,K) = 0
        ENDDO
        DO J = 1,LSPC_CW
          ISPC_CW(J,K) = 0
        ENDDO
        DO J = 1,LEQC
          ISOLC_CW(J,K) = 0
        ENDDO
        DO J = 1,LEQK
          ISOLK_CW(J,K) = 0
        ENDDO
        DO 310 J = 1,LUK_CW
          KLU_CW(J,K) = 0
  310   CONTINUE
  320 CONTINUE
      DO 400 L = 1,LWF_CW
        IWF_CW(L) = 0
        PF_CW(L) = 0.D+0
  400 CONTINUE
      DO 410 L = 1,LVREF
        IREF_CW(L) = 0
  410 CONTINUE
      N_CW = 0
      N_GLW = 0
      NWF_CW = 0
      NWN_CW = 0
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_COUP_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_LEAK_WELL
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
!     Initialize leaky well variables.
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 10 October 2018.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE LEAK_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_LEAK_WELL'
!
!---  Coupled well variables  ---
!
      DO L = 1,LN_LW
        WNM_LW(L) = 'null'
        DO K = 1,8
          ID_LW(K,L) = 0
        ENDDO
        DO K = 1,3
          FF_LW(K,L) = 0.D+0
        ENDDO
        IT_LW(L) = 0
      ENDDO
      DO L = 1,LWI_LW
        DO K = 1,2
          XTP_LW(K,L) = 0.D+0
          YTP_LW(K,L) = 0.D+0
          ZTP_LW(K,L) = 0.D+0
        ENDDO
        DO K = 1,5
          PAR_LW(K,L) = 0.D+0
        ENDDO
        IS_LW(L) = 0
      ENDDO
      DO L = 1,LWN_LW
        DO K = 1,2
          XP_LW(K,L) = 0.D+0
          YP_LW(K,L) = 0.D+0
          ZP_LW(K,L) = 0.D+0
        ENDDO
        PLX_LW(L) = 0.D+0
        PLY_LW(L) = 0.D+0
        PLZ_LW(L) = 0.D+0
        INV_LW(L) = 0
        NF_LW(L) = 0
        IWP_LW(L) = 0
        ND_LW(L) = 0
      ENDDO
      DO L = 1,LWF_LW
        IWF_LW(L) = 0
      ENDDO
      DO L = 1,LUK
        DO K = 1,LWN_LW*2*LUK
          KLU_LW(K,L) = 0
        ENDDO
      ENDDO
      DO L = 1,LWN_LW*2
        KLUC_LW(L) = 0
      ENDDO
      N_LW = 0
      NWN_LW = 0
      NWF_LW = 0
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_LEAK_WELL group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_HYDT
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
!     Initialize ternary hydrate variables.
!
!     ZPC_HT - mole fraction of hydrate formers for phase equilibria
!     ZHC_HT - mole fraction of hydrate formers for hydrate equilibria
!     DZPC_HT - mole fraction spacing for phase equilibria 
!     DZHC_HT - mole fraction spacing for hydrate equilibria 
!     TCR_HT - critical temperature, K
!     PCR_HT - critical pressure, Pa
!     TCB_HT - cricondenbar temperature, K
!     PCB_HT - cricondenbar pressure, Pa
!     TCT_HT - cricondenterm temperature, K
!     PCT_HT - cricondenterm pressure, Pa
!     TLE_HT - temperature for lower phase envelope, K
!     PLE_HT - pressure for lower phase envelope, Pa
!     D2PLE_HT - second derivative of pressure for lower phase envelope
!     FKLE_HT - K-factor for lower phase envelope, Pa
!     D2FKLE_HT - second derivative of K-factor for lower phase envelope
!     TUE_HT - temperature for upper phase envelope, K
!     PUE_HT - pressure for upper phase envelope, Pa
!     D2PUE_HT - second derivative of pressure for upper phase envelope
!     FKUE_HT - K-factor for upper phase envelope, Pa
!     D2FKUE_HT - second derivative of K-factor for upper phase envelope
!     T2P_HT - temperature for two-phase region, K
!     P2P_HT - pressure for two-phase region, Pa
!     B2P_HT - beta for two-phase region
!     D2B2P_HT - second derivative of beta for two-phase region
!     FK2P_HT - K-factor for two-phase region
!     D2FK2P_HT - second derivative of K-factor for two-phase region
!     THE_HT - temperature for hydrate equilibrium table, K
!     PHE_HT - pressure for hydrate equilibrium table, Pa
!     XSCA_HT - CO2 small cage occupany for hydrate equilibrium table
!     XSCO_HT - CH4 small cage occupany for hydrate equilibrium table
!     XSCN_HT - N2 small cage occupany for hydrate equilibrium table
!     XLCA_HT - CO2 large cage occupany for hydrate equilibrium table
!     XLCO_HT - CH4 large cage occupany for hydrate equilibrium table
!     XLCN_HT - N2 large cage occupany for hydrate equilibrium table
!     THE2P_HT - second derivative of temperature for hydrate equil.
!     XSCA2T_HT - second derivative of CO2 small cage occupany for hyd.
!     XSCO2T_HT - second derivative of CH4 small cage occupany for hyd.
!     XSCN2T_HT - second derivative of N2 small cage occupany for hyd.
!     XLCA2T_HT - second derivative of CO2 large cage occupany for hyd.
!     XLCO2T_HT - second derivative of CH4 large cage occupany for hyd.
!     XLCN2T_HT - second derivative of N2 large cage occupany for hyd.
!     PHE2T_HT - second derivative of pressure for hydrate equilibrium
!     NPEP_HT - number of phase envelope table points
!     NHEP_HT - number of hydrate equilibrium table points
!     NTP_HT - number of two-phase temperature points
!     NPP_HT - number of two-phase presure points
!     NHF_HT - number of hydrate formers
!     NZP_HT - dimension of phase envelope table
!     NZH_HT - dimension of hydrate table
!     IZP_HT - bounding concentrations indices for phase equilibria
!     IZH_HT - bounding concentrations indices for hydrate equilibrium
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 July 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE HYDT
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_HYDT'
!
!---  Phase equilibria variables  ---
!
      DZPC_HT = 1.D-1
      DZHC_HT = 1.D-1
      DO 400 L = 1,LCP_HT
        DO 100 K = 1,LHF_HT
          ZPC_HT(K,L) = 0.D+0
  100   CONTINUE
        TCR_HT(L) = 0.D+0
        PCR_HT(L) = 0.D+0
        TCB_HT(L) = 0.D+0
        PCB_HT(L) = 0.D+0
        TCT_HT(L) = 0.D+0
        PCT_HT(L) = 0.D+0
        NPEP_HT(L) = 0
        DO 202 K = 1,LPE_HT
          TLE_HT(K,L) = 0.D+0
          PLE_HT(K,L) = 0.D+0
          D2PLE_HT(K,L) = 0.D+0
          TUE_HT(K,L) = 0.D+0
          PUE_HT(K,L) = 0.D+0
          D2PUE_HT(K,L) = 0.D+0
          DO 200 M = 1,LHF_HT
            FKLE_HT(K,L,M) = 0.D+0
            D2FKLE_HT(K,L,M) = 0.D+0
            FKUE_HT(K,L,M) = 0.D+0
            D2FKUE_HT(K,L,M) = 0.D+0
  200     CONTINUE
  202   CONTINUE
        NTP_HT(L) = 0
        DO 210 K = 1,LTP_HT
          T2P_HT(K,L) = 0.D+0
  210   CONTINUE
        NPP_HT(L) = 0
        DO 220 K = 1,LPP_HT
          P2P_HT(K,L) = 0.D+0
  220   CONTINUE
        DO 320 K = 1,LPP_HT
          DO 310 J = 1,LTP_HT
            B2P_HT(J,K,L) = 0.D+0
            D2B2P_HT(J,K,L) = 0.D+0
  310     CONTINUE
  320   CONTINUE
  400 CONTINUE
      DO 440 L = 1,LHF_HT
        DO 430 K = 1,LCP_HT
          DO 420 J = 1,LPP_HT
            DO 410 I = 1,LTP_HT
              FK2P_HT(I,J,K,L) = 0.D+0
              D2FK2P_HT(I,J,K,L) = 0.D+0
  410       CONTINUE
  420     CONTINUE
  430   CONTINUE
  440 CONTINUE
      DO 520 L = 1,LCN_HT
        DO 510 K = 1,LCN_HT
          IZP_HT(K,L) = 0
  510   CONTINUE
  520 CONTINUE
!
!---  Hydrate equilibrium variables  ---
!
      DO 610 L = 1,LCH_HT
        DO 600 K = 1,LHF_HT
          ZHC_HT(K,L) = 0.D+0
  600   CONTINUE
        NHEP_HT(L) = 0
  610 CONTINUE
      DO 630 L = 1,LCN_HT
        DO 620 K = 1,LCN_HT
        IZH_HT(K,L) = 0
  620   CONTINUE
  630 CONTINUE
      DO 650 L = 1,LCH_HT
        DO 640 K = 1,LHE_HT
          THE_HT(K,L) = 0.D+0
          PHE_HT(K,L) = 0.D+0
          XSCA_HT(K,L) = 0.D+0
          XSCO_HT(K,L) = 0.D+0
          XSCN_HT(K,L) = 0.D+0
          XLCA_HT(K,L) = 0.D+0
          XLCO_HT(K,L) = 0.D+0
          XLCN_HT(K,L) = 0.D+0
          THE2P_HT(K,L) = 0.D+0
          XSCA2T_HT(K,L) = 0.D+0
          XSCO2T_HT(K,L) = 0.D+0
          XSCN2T_HT(K,L) = 0.D+0
          XLCA2T_HT(K,L) = 0.D+0
          XLCO2T_HT(K,L) = 0.D+0
          XLCN2T_HT(K,L) = 0.D+0
          PHE2T_HT(K,L) = 0.D+0
  640   CONTINUE
  650 CONTINUE
      DO 660 L = 1,LHF_HT
        ZMIH_HT(L) = 0.D+0
  660 CONTINUE
      NZP_HT = 0
      NZH_HT = 0
      NHF_HT = 0
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_HYDT group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_EOR
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
!     Initialize enhanced oil recovery variables.
!
!     BIPC - binary interaction parameters of petroleum components
!     BIPF - binary interaction parameters of petroleum fractions
!     GC_MMP - gas composition for MMP
!     PFNM - petroleum fraction name
!     PFPP - petroleum fraction properties
!     XMPCF - mole fraction of petroleum component fractions
!     PR0LKR - reduced pressure for the Lee-Kesler H0 residual enthalpy
!     TR0LKR - reduced temp. for the Lee-Kesler H0 residual enthalpy
!     H0LKR - enthalpy for the Lee-Kesler H0 residual enthalpy
!     PR1LKR - reduced pressure for the Lee-Kesler H1 residual enthalpy
!     TR1LKR - reduced temp. for the Lee-Kesler H1 residual enthalpy
!     H1LKR - enthalpy for the Lee-Kesler H0 residual enthalpy
!     IPCF - indices of petroleum component fractions
!     NPF - number of petroleum fractions
!     NPCF - number of petroleum component fractions
!
!       PFPP(1,I) molecular weight, kg/kmol
!       PFPP(2,I) normal boiling point, K
!       PFPP(3,I) critical temperature, K
!       PFPP(4,I) critical pressure, Pa
!       PFPP(5,I) critical molar volume, cm^3/mole
!       PFPP(6,I) critical compressibility
!       PFPP(7,I) Pitzer acentric factor
!       PFPP(8,I)  ---
!       PFPP(9,I)  ---
!       PFPP(13,I) Peng-Robinson pure component parameter a
!       PFPP(14,I) Peng-Robinson pure component parameter b
!       PFPP(15,I) Peng-Robinson pure component parameter fw
!       PFPP(16,I)  ---
!       PFPP(17,I)  ---
!       PFPP(18,I)  ---
!       PFPP(19,I)  ---
!       PFPP(20,I)  ---
!       PFPP(21,I)  ---
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 26 July 2012.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE EOR
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_EOR'
!
!---  Loop over number of petroleum fractions  ---
!
      DO L = 1,LPF_EOR
        PFNM(L) = 'null'
!
!---    Loop over number of petroleum fraction parameters  ---
!
        DO K = 1,10
          PFPP(K,L) = 0.D+0
        ENDDO
      ENDDO
      NPF = 0
!
!---  Loop over number of petroleum components  ---
!
      DO L = 1,LNGC
        NPCF(L) = 0
        GC_MMP(L) = 0.D+0
!
!---    Loop over number of petroleum-component fractions  ---
!
        DO K = 1,LPCF
          XMPCF(K,L) = 0.D+0
          IPCF(K,L) = 0
        ENDDO
!
!---    Loop over number of petroleum components  ---
!
        DO K = 1,LNGC
          BIPC(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Loop over CO2, CH4 and number of petroleum fractions  ---
!
      DO L = 1,LPF_EOR+2
!
!---    Loop over CO2, CH4 and number of petroleum fractions   ---
!
        DO K = 1,LPF_EOR+2
          BIPF(K,L) = 0.D+0
        ENDDO
      ENDDO
!
!---  Tables for the Lee-Kesler residual enthalpy   ---
!
      DO N = 1,40
        TR0LKR(N) = 0.D+0
        TR1LKR(N) = 0.D+0
        DO M = 1,15
          H0LKR(M,N) = 0.D+0
          H1LKR(M,N) = 0.D+0
        ENDDO
      ENDDO
      DO M = 1,15
        IL0LKR(M) = 0
        IL1LKR(M) = 0
        PR0LKR(M) = 0.D+0
        PR1LKR(M) = 0.D+0
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_EOR group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE IN_PLT_ATM
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
!     Initialize surface cover variables.
!
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 24 November 2015.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE PLT_ATM
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/IN_PLT_ATM'
!
!---  [ PLT_ATM ] Surface cover, plant, and atmospheric variables  ---
!
      NPLANT = 0
      ISWBCF = 0
      NSFCA = 0
      NSFCC = 0
      NSFCN = 0
      ISCVF = 0
!
!---  Loop over number of surface cover areas  ---
!
      DO K = 1,LSFCA
        NAME_SFC(K) = 'null'
        NSFCT(K) = 0
        IALB(K) = 0
        PHMX(K) = 0.D+0
        DO J = 1,5
          ALBEDO(J,K) = 0.D+0
        ENDDO
!
!---    Loop over number of surface cover times  ---
!
        DO J = 1,LSFCT
!
!---      Loop over number of plant species times 2 plus 1  ---
!
          DO I = 1,2*LPLANT+1
            PARMS_SFC(I,J,K) = 0.D+0
          ENDDO
        ENDDO
      ENDDO
!
!---  Loop over number of surface cover nodes  ---
!
      DO J = 1,LSFCN
        ID_SFC(J) = 0
!
!---    Loop over number of ground-surface equation residual 
!       increments  ---
!
        DO I = 1,LUK+6
          RE_GS(I,J) = 0.D+0
          RW_GS(I,J) = 0.D+0
        ENDDO
!
!---    Loop over number of plant and canopy equation residual 
!       increments  ---
!
        DO I = 1,LUK*LSFCC+6
          RE_PL(I,J) = 0.D+0
          RE_CP(I,J) = 0.D+0
          RW_CP(I,J) = 0.D+0
        ENDDO
!
!---    Loop over number of surface cover equations  ---
!
        DO I = 1,5
          DNR_SFC(I,J) = 0.D+0
        ENDDO
!
!---    Loop over number primary unknown indices  ---
!
        DO I = 1,7
          ACW_PL(I,J) = 0.D+0
          DM_PL(I,J) = 0.D+0
          WM_PL(I,J) = 0.D+0
          DFGW_GS(I,J) = 0.D+0
          DFLA_GS(I,J) = 0.D+0
          HG_GS(I,J) = 0.D+0
          HGA_GS(I,J) = 0.D+0
          HGW_GS(I,J) = 0.D+0
          HL_GS(I,J) = 0.D+0
          HLW_GS(I,J) = 0.D+0
          HLW_PL(I,J) = 0.D+0
          PVA_GS(I,J) = 0.D+0
          PVW_CP(I,J) = 0.D+0
          PVW_GS(I,J) = 0.D+0
          RHOG_GS(I,J) = 0.D+0
          RHOL_GS(I,J) = 0.D+0
          RHOMG_GS(I,J) = 0.D+0
          RHOML_GS(I,J) = 0.D+0
          RKG_GS(I,J) = 0.D+0
          RKL_GS(I,J) = 0.D+0
          PL_GS(I,J) = 0.D+0
          SG_GS(I,J) = 0.D+0
          SL_GS(I,J) = 0.D+0
          T_CP(I,J) = 0.D+0
          T_GS(I,J) = 0.D+0
          T_PL(I,J) = 0.D+0
          THKG_GS(I,J) = 0.D+0
          THKL_GS(I,J) = 0.D+0
          TORG_GS(I,J) = 0.D+0
          TORL_GS(I,J) = 0.D+0
          UEG_GS(I,J) = 0.D+0
          VISG_GS(I,J) = 0.D+0
          VISL_GS(I,J) = 0.D+0
          XGA_GS(I,J) = 0.D+0
          XGW_GS(I,J) = 0.D+0
          XLA_GS(I,J) = 0.D+0
          XLW_GS(I,J) = 0.D+0
          XMGA_GS(I,J) = 0.D+0
          XMGW_GS(I,J) = 0.D+0
          XMLA_GS(I,J) = 0.D+0
          XMLW_GS(I,J) = 0.D+0
        ENDDO
!
!---    Loop over number of surface cover connections  ---
!
        DO I = 1,LSFCC
          ICM_SFC(I,J) = 0
        ENDDO
!
!---    Loop over number of surface cover equations  ---
!
        DO I = 1,5
          JM_SFC(I,J) = 0
        ENDDO
!
!---    Loop over number of surface cover equation matrix elements  ---
!
        DO I = 1,LUK_SFC
          KLU_SFC(I,J) = 0
        ENDDO
      ENDDO
!
!---  Loop over number of plant species  ---
!
      DO J = 1,LPLANT
        DO I = 1,29
          PARMS_P(I,J) = 0.D+0
        ENDDO
!
!---    Minimum maximum plant root depth of 1 mm  ---
!
        PARMS_P(1,J) = 1.D-3
!
!---    Minimum plant canopy height of 1 mm  ---
!
        PARMS_P(11,J) = 1.D-3     
        IRSM_P(J) = 0
        ITMP_P(J) = 0
        IALB_P(J) = 0
        ISRM_P(J) = 0
        PLANT(J) = 'null'
      ENDDO
      NATM_T = 1
      IATM_C = 0
      DO J = 1,7
        DO I = 1,LATM
          ATMOS(I,J) = 0.D+0
        ENDDO
      ENDDO
      ATMC(1) = 1.D+0
      ATMC(2) = 1.D+0
      ATMC(3) = 0.D+0
      ATMC(4) = 0.D+0
      ATMC(5) = 120.D+0
      ATMC(6) = 1.3D-3
      ATMC(7) = 1.3D-3
      ATMST = 0.D+0
!
!---  Transpiration  ---
!
      DO J = 1,LSFCN
        DO I = 1,6
          TP_SFC(I,J) = 0.D+0
        ENDDO
      ENDDO
!
!---  Root water uptake  ---
!
      DO J = 1,LSFCN*LSFCC
        DO I = 1,6+LUK
          RWU_SFC(I,J) = 0.D+0
        ENDDO
      ENDDO
!
!---  Stored water mass on plant leaves, kg/m^2 ground  ---
!
      DO J = 1,LSFCN
        DO I = 1,7
          SWM_PL(I,J) = 0.D+0
        ENDDO
      ENDDO
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of IN_PLT_ATM group  ---
!
      RETURN
      END

