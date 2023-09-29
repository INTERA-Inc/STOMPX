!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REFNOD
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
!     Prints convergence and variable information to the screen and
!     output file.
!
!   1 aqueous pressure (Pa), ' PL '
!   2 gas pressure (Pa), ' PG '
!   3 NAPL pressure (Pa), ' PN '
!   4 temperature (C), ' T  '
!   5 phase condition (null) 'PHCN'
!   6 aqueous gauge pressure (Pa), 'GPL '
!   7 gas gauge pressure (Pa), 'GPG '
!   8 NAPL gauge pressure (Pa), 'GPN '
!   9 apparent aqueous saturation (null), 'ASL ' or
!     hydrate area (m^2)
!   10 apparent total-liquid saturation (null), 'AST ' or
!      hydrate dissociation/formation rate (kmol/s m^2)
!      formation is positive and dissociation is negative
!   11 aqueous saturation (null), ' SL '
!   12 gas saturation (null), ' SG '
!   13 NAPL saturation (null), ' SN '
!   14 total-liquid saturation (null), ' ST '
!   15 aqueous moisture content (null), 'MCL '
!   16 NAPL moisture content (null), 'MCN '
!   17 total-liquid moisture content (null), 'MCT '
!   18 effective trapped NAPL (null), 'ESNT'
!   19 effective trapped gas (null), 'ESGT'
!   20 diffusive porosity (null), 'PORD'
!   21 gas water mass fraction (null), 'XGW '
!       or volumetric component mass fraction (null), 'XVC'
!   22 gas air mass fraction (null), 'XGA '
!       or gas component mass fraction (null), 'XGC'
!   23 gas oil mass fraction (null), 'XGO '
!       or NAPL component mass fraction (null), 'XNC'
!   24 aqueous water mass fraction (null), 'XLW '
!   25 aqueous air mass fraction (null), 'XLA '
!   26 aqueous oil mass fraction (null), 'XLO '
!   27 aqueous hydraulic head (m), 'HHL '
!   28 gas hydraulic head (m), 'HHG '
!   29 NAPL hydraulic head (m), 'HHN '
!   30 rock/soil type (null), 'RSZN'
!   31 aqueous relative permeability (null), 'RKL '
!   32 gas relative permeability (null), 'RKG '
!   33 NAPL relative permeability (null), 'RKN '
!       or liquid-CO2 relative permeability (null), 'RKN '
!   34 aqueous density (kg/m^3), 'RHOL'
!   35 gas density (kg/m^3), 'RHOG'
!   36 NAPL density (kg/m^3), 'RHON'
!   37 total water mass (kg), 'TMW '
!   38 total air mass (kg), 'TMA '
!       or total component mass (kg), 'TMC'
!   39 total oil mass (kg), 'TMO '
!   40 water mass source integral (kg), 'SRIW'
!   41 air mass source integral (kg), 'SRIA'
!       or component mass source integral (kg), 'SRIC'
!   42 oil mass source integral (kg), 'SRIO'
!   43 energy source integral (W), 'SRIT'
!   44 x thermal conductivity (W/m K), 'THKX'
!   45 y thermal conductivity (W/m K), 'THKY'
!   46 z thermal conductivity (W/m K), 'THKZ'
!   47 salt volumetric concentration (kg/m^3), ' CS '
!   48 salt aqueous concentration (kg/m^3), 'CSL '
!   49 aqueous courant (null), 'CRNL'
!   50 total salt mass (kg), 'TMS '
!   51 x aqueous volumetric flux (m/s), ' UL '
!   52 y aqueous volumetric flux (m/s), ' VL '
!   53 z aqueous volumetric flux (m/s), ' WL '
!   54 x gas volumetric flux (m/s), ' UG '
!   55 y gas volumetric flux (m/s), ' VG '
!   56 z gas volumetric flux (m/s), ' WG '
!   57 x NAPL volumetric flux (m/s), ' UN '
!   58 y NAPL volumetric flux (m/s), ' VN '
!   59 z NAPL volumetric flux (m/s), ' WN '
!   60 x heat flux (W/m^2), ' UQ '
!   61 y heat flux (W/m^2), ' VQ '
!   62 z heat flux (W/m^2), ' WQ '
!   63 matric potential head (m), 'MPH '
!       or Webb matching point head (m), 'WMPH'
!   64 x salt flux (kg/m^2 s), ' US '
!   65 y salt flux (kg/m^2 s), ' VS '
!   66 z salt flux (kg/m^2 s), ' WS '
!   67 xnc salt flux (kg/m^2 s), 'USNC'
!   68 ync salt flux (kg/m^2 s), 'VSNC'
!   69 znc salt flux (kg/m^2 s), 'WSNC'
!   70 gas water mole fraction (null), 'XMGW'
!   71 gas air mole fraction (null), 'XMGA'
!       or gas component mole fraction (null), 'XMGC'
!   72 gas oil mole fraction (null), 'XMGO'
!       or NAPL component mole fraction (null), 'XMNC'
!   73 gas water concentration (kg/m^3), 'CGW '
!       or volumetric component concentration (kg/m^3), 'CVC'
!   74 gas air concentration (kg/m^3), 'CGA '
!       or gas component concentration (kg/m^3), 'CGC'
!   75 gas oil concentration (kg/m^3), 'CGO '
!       or NAPL component concentration (kg/m^3), 'CNC'
!   76 aqueous water concentration (kg/m^3), 'CLW '
!   77 aqueous air concentration (kg/m^3), 'CLA '
!   78 aqueous oil concentration (kg/m^3), 'CLO '
!   79 gas courant (null), 'CRNG'
!   80 ice pressure (Pa), 'PI '
!   80 system pressure (Pa), 'P'
!   80 stress-xx (Pa), 'SIGXX'
!   81 ice saturation (null), 'SI '
!      or beta factor (null), 'BETA'
!   82 ice density (kg/m^3), 'RHOF'
!   83 aqueous matrix (null), 'DSLM'
!       or NAPL matrix (null), 'DSNM'
!       or Eclipse gas saturation (null), 'ESG'
!   84 aqueous fracture (null), 'DSLF'
!       or NAPL fracture (null), 'DSNF'
!       or Eclipse oil saturation (null), 'ESO'
!       or Eclipse pressure (null), 'EPR'
!   85 gas matrix (null), 'DSGM'
!   86 gas fracture (null), 'DSGF'
!   87 xnc aqueous volumetric flux (m/s), 'ULNC'
!   88 ync aqueous volumetric flux (m/s), 'VLNC'
!   89 znc aqueous volumetric flux (m/s), 'WLNC'
!   90 xnc gas volumetric flux (m/s), 'UGNC'
!   91 ync gas volumetric flux (m/s), 'VGNC'
!   92 znc gas volumetric flux (m/s), 'WGNC'
!   93 xnc NAPL volumetric flux (m/s), 'UNNC'
!   94 ync NAPL volumetric flux (m/s), 'VNNC'
!   95 znc NAPL volumetric flux (m/s), 'WNNC'
!   96 xnc heat flux (W/m^2), 'UQNC'
!   97 ync heat flux (W/m^2), 'VQNC'
!   98 znc heat flux (W/m^2), 'WQNC'
!   99 NAPL courant (null), 'CRNN'
!   100 node number (null), 'NODE'
!   101 osmotic pressure (Pa), 'POSM'
!   101 stress-yy (Pa), 'SIGXX'
!   102 osmotic efficiency factor (null), 'OEC '
!   103 aqueous compressibility (1/Pa), CMPL
!   104 gas compressibility (1/Pa), CMPG
!   105 trapped gas saturation (null), 'SGT '
!   106 trapped NAPL saturation (null), 'SNT '
!   107 aqueous trapped gas (null), 'SGTL'
!   108 NAPL trapped gas (null), 'SGTN'
!   109 dissolved-aqueous oil concentration (kg/m^3), 'CLO '
!   110 salt aqueous mass fraction (null), 'XLS '
!   111 surfactant volumetric concentration (kg/m^3), 'CS '
!   112 surfactant aqueous concentration (kg/m^3), 'CLS '
!   113 surfactant aqueous mass fraction (null), 'XLS '
!       or gas N2 mass fraction (null), 'XGN'
!   114 x surfactant flux (kg/m^2 s), ' US '
!   115 y surfactant flux (kg/m^2 s), ' VS '
!   116 z surfactant flux (kg/m^2 s), ' WS '
!   117 xnc surfactant flux (kg/m^2 s), 'USNC'
!   118 ync surfactant flux (kg/m^2 s), 'VSNC'
!   119 znc surfactant flux (kg/m^2 s), 'WSNC'
!   120 x dissolved-oil flux (kg/m^2 s), 'ULO '
!   121 y dissolved-oil flux (kg/m^2 s), 'VLO '
!   122 z dissolved-oil flux (kg/m^2 s), 'WLO '
!   123 xnc dissolved-oil flux (kg/m^2 s), 'ULOC'
!   124 ync dissolved-oil flux (kg/m^2 s), 'VLOC'
!   125 znc dissolved-oil flux (kg/m^2 s), 'WLOC'
!   126 NAPL-aqueous trapping number (null), 'TPNL'
!   127 minimum effect aqueous saturation (null), 'ESLM'
!   128 water vapor partial pressure (Pa), 'PVW '
!   129 air partial pressure (Pa), 'PVA '
!   130 oil vapor partial pressure (Pa), 'PVO '
!   130 stress-zz (Pa), 'SIGZZ'
!   131 gas-aqueous scaling factor (null), 'BGL '
!   132 free-NAPL aqueous interfacial area (m^2), 'ANFL'
!      or X CO2 displacement variance (m^2), 'SIGXX'
!   133 trapped-NAPL aqueous interfacial area (m^2), 'ANTL'
!      or Y CO2 displacement variance (m^2), 'SIGYY'
!   134 aqueous solute coefficient 'HKL '
!      or Z CO2 displacement variance (m^2), 'SIGZZ'
!   135 free-NAPL solute coefficient 'HKNF'
!   136 trapped-NAPL solute coefficient 'HKNT'
!   137 water relative humidity (null), 'RHW '
!   138 well pressure (Pa), 'PLWB'
!       or coupled-well pressure (Pa), 'PCW'
!       or injection well pressure (Pa), 'IWP'
!   139 volumetric molar density (kg/m^3), 'RHMV'
!   140 water mass source rate (kg/s), 'SRCW'
!       or total water mass source rate (kg/s), 'TSRW'
!   141 air mass source rate (kg/s), 'SRCA'
!       or component mass source rate (kg/s), 'SRCC'
!   142 oil mass source rate (kg/s), 'SRCO'
!       or coupled-well CO2 mass nodal rate (kg/s), 'QNRA'
!       or coupled-well volumetric nodal rate (m^3/s), 'QNVR'
!   143 energy source rate (W), 'SRCQ'
!   144 aqueous well depth (m), 'PLWB'
!   145 well mass flow rate (kg/s), 'QLW '
!       or coupled-well water mass nodal rate (kg/s), 'QNRW'
!   146 well flow integral (kg), 'QLWI'
!   147 salt mass source rate (kg/s), 'SRCS'
!   148 salt mass source integral (kg), 'SRIS'
!   149 scanning path (null), 'PATH'
!   150 aqueous air or co2 saturation (null), 'DAPS'
!       or Webb matching point saturation (null), 'WMPS'
!   151 bubble void fraction (null), 'BVF '
!       or aqueous N2 mole fraction (null), 'XMLN'
!   152 bubble air mass fraction (null), 'XBA '
!       or aqueous N2 mass fraction (null), 'XLN'
!   153 mineralized co2 'MCO2 '
!   154 napl well flow rate 'QNW '
!   155 napl well flow integral 'QNWI'
!   156 total well flow rate 'QTW '
!   157 total well flow integral 'QTWI'
!   160 gas N2 mole fraction (null), 'XMGN'
!   161 NAPL dissolved water concentration 'CNW '
!   162 NAPL dissolved water mole fraction (null), 'XMNW'
!       or N2 nonaqueous liquid mole fraction (null), 'XMNN'
!   163 NAPL dissolved water mass fraction (null), 'XNW '
!       or water liquid-CO2 mass fraction (null), 'XNW '
!   164 NAPL dissolved oil concentration 'CNO '
!   165 NAPL dissolved oil mole fraction (null), 'XMNO'
!       or CH4 nonaqueous liquid mole fraction (null), 'XMNO'
!   166 NAPL dissolved oil mass fraction (null), 'XNO '
!       or CH4 liquid-CO2 mass fraction (null), 'XNO '
!   167 total alcohol mass 'TMA '
!   168 aqueous dissolved water mass fraction (null), 'XLW '
!   169 alcohol mass source integral 'SRIA'
!   170 alcohol mass source rate 'SRCA'
!   171 integrated NAPL and aqueous dissolved alcohol 'IMA '
!       or integrated kerogen mass 'IMK'
!       or integrated hydrate N2 mass 'IMHN'
!   172 integrated aqueous dissolved alcohol 'IMLA'
!       or integrated char mass 'IMCH'
!   173 integrated NAPL dissolved water 'IMNW'
!       or integrated coke mass 'IMCK'
!   174 integrated NAPL dissolved oil 'IMNO'
!   175 integrated NAPL dissolved alcohol 'IMNA'
!   176 aqueous viscosity 'VISL'
!   177 monitoring well water depth or total-liquid well depth 'WDT '
!   178 aqueous well depth 'WDL '
!   179 NAPL well depth 'WDN '
!   180 monitoring well pressure ' PW '
!       or coupled-well nodal pressure 'PNCW'
!   181 monitoring well aqueous saturation 'SLW '
!   182 monitoring well water-vapor mass fraction or well NAPL saturation 'XGWW'
!   183 monitoring well dissolved-air mass fraction or dissolved-oil mass fraction 'XLAW'
!   184 monitoring well axial aqueous flux or well total-liquid pumping rate 'UL_W'
!   185 monitoring well axial gas flux or well aqueous pumping rate 'UG_W'
!   186 monitoring well vertical aqueous flux or well NAPL pumping rate 'WL_W'
!   187 monitoring well vertical gas flux or well total-liquid pumping integral 'WG_W'
!   188 integrated well aqueous pumping 'IPLW'
!   189 integrated mineral CO2 or well NAPL pumping integral 'IPNW'
!   190 integrated trapped gas air 'IMGT'
!       or integrated N2 mass
!   191 integrated water mass (kg), 'IMW '
!   192 integrated air mass (kg), 'IMA '
!       or integrated petroleum component mass (kg), IMC
!   193 integrated oil mass (kg), 'IMO '
!   194 integrated aqueous water mass (kg), 'IMLW'
!   195 integrated aqueous air mass (kg), 'IMLA'
!   196 integrated aqueous oil mass (kg), 'IMLO'
!   197 integrated gas water mass (kg), 'IMGW'
!       or integrated volumetric component mass (kg), 'IMVC'
!   198 integrated gas air mass (kg), 'IMGA'
!       or integrated gas component mass (kg), 'IMGC'
!   199 integrated gas oil mass (kg), 'IMGO'
!       or integrated NAPL component mass (kg), 'IMNC'
!   200 reserved to control plot file output
!   201 x aqueous relative permeability (null), 'RKLX'
!   202 y aqueous relative permeability (null), 'RKLY'
!   203 z aqueous relative permeability (null), 'RKLZ'
!   204 aqueous co2 mole fraction (null), 'XMLA'
!       or aqueous component mole fraction (null), 'XMLC' 
!   205 aqueous salt mole fraction (null), 'XMLS'
!   206 atmospheric temperature (C), ' TA '
!       or critical temperature (C), ' TCR '
!       or ground-loop well head temperature (C), 'TGLW'
!   207 atmospheric relative humidity (null), ' RH '
!   208 atmospheric solar radiation ' RN '
!   209 atmospheric wind speed (m/s), ' WS '
!   210 residual NAPL saturation (null), 'SNR '
!   211 mobile NAPL saturation (null), 'SNM '
!   212 free NAPL saturation (null), 'SNF '
!   213 surface temperature (C), 'T_S'
!       or cricondentherm temperature (C), ' TCT '
!   214 surface vapor pressure (Pa), 'PV_S'
!   214 stress-yz (Pa), 'SIGYZ'
!   215 actual evaporation rate 'E_SA'
!   216 potential evaporation rate 'PE_SA'
!   217 actual transpiration rate 'T_SA'
!   218 potential transpiration rate 'PT_SA'
!   219 saturated co2 aqueous mass fraction 'SXLA'
!   220 aqueous alcohol mole fraction (null), 'XMLA'
!       or aqueous N2 mole fraction (null), 'XMLN'
!   221 NAPL alcohol mode fraction (null), 'XMNA'
!       or CO2 nonaqueous liquid mole fraction (null), 'XMNA'
!   222 aqueous alcohol mass fraction (null), 'XLA '
!       or nonaqueous N2 liquid mass fraction (null), 'XNN'
!   223 NAPL alcohol mass fraction (null), 'XNA '
!       or CO2 liquid-CO2 mass fraction (null), 'XNA '
!   224 atmospheric pressure (Pa), ' PA '
!   224 stress-xz (Pa), 'SIGXZ'
!   225 surface aqueous pressure (Pa), 'PL_S'
!   226 surface gas pressure (Pa), 'PG_S'
!   226 stress-xy (Pa), 'SIGXY'
!   227 surface aqueous saturation (null), 'SL_S
!   228 surface latent heat flux (W/m^2), 'QL_S'
!   229 surface sensible heat flux (W/m^2), 'QH_S'
!   230 surface net long-wave radiation, 'RL_S'
!   231 surface net short-wave radiation, 'RS_S'
!   232 surface ground heat flux, 'QG_S'
!   233 surface water mass balance (kg/s), 'WB_S'
!   234 plant temperature (C), 'T_P' or 'TPXX'
!   235 plant temperature (C), 'TPXX'
!   236 plant temperature (C), 'TPXX'
!   237 plant temperature (C), 'TPXX'
!   238 plant temperature (C), 'TPXX'
!   239 rainfall interception mass (kg), 'RFIM'
!       or integrated aqueous N2 mass (kg), 'IMLN'
!   240 sorbed oil mass, (kg oil), 'TSO '
!       or integrated gas N2 mass (kg), 'IMGN'
!   241 sorbed oil mass fraction (kg oil/kg soil), 'XSO '
!       or hydrate N2 mass fraction (null), 'XHN'
!   242 sorbed oil volumetric concentration (kg oil/m^3), 'CSO '
!   243 bare-soil aerodynamic resistance (s/m), 'RABS'
!   244 surface volumetric precipitation rate (m^3/s), 'PV_S'
!   245 surface mass precipitation rate (kg/s), 'PM_S'
!       or fracture source gas mass rate (kg/s), 'SRCG'
!   246 atmospheric precipitation rate (kg/s), 'PM_A'
!       or injection well petroleum component mass rate (kg/s), 'IWCR'
!       or production well petroleum component mass rate (kg/s), 'PWCR'
!       or fracture source aqueous mass rate (kg/s), 'SRCL'
!   247 x-direction intrinsic permeability (m^2), ' UK '
!   247 matrix permeability (m^2), 'KM' (STOMP-OS)
!   248 y-direction intrinsic permeability (m^2), ' VK '
!   248 fracture permeability (m^2), 'FM' (STOMP-OS)
!   249 z-direction intrinsic permeability (m^2), ' WK '
!   250 hydrate water mass fraction (null), 'XHW '
!   251 hydrate CO2 mass fraction (null), 'XHA '
!   252 hydrate CH4 mass fraction (null), 'XHO '
!   253 hydrate density (kg/m^3), 'RHOH'
!       or production well fluid density (kg/m^3), 'PWFD'
!   254 hydrate saturation (null), ' SH '
!   255 hydrate pressure (Pa), ' PH '
!   256 integrated hydrate water mass (kg), 'IMHW'
!   257 integrated hydrate CO2 mass (kg), 'IMHA'
!       or injection well petroleum comp. mass integral (kg), 'IWCI'
!       or production well petroleum comp. mass integral (kg), 'PWCI'
!   258 integrated hydrate CH4 mass (kg), 'IMHO'
!   259 integrated aqueous CH4 mass (kg), 'IMLO'
!   260 integrated gas CH4 mass (kg), 'IMGO'
!   261 integrated water mass source (kg), 'IMSW'
!   262 integrated air or CO2 mass source (kg), 'IMSA'
!   263 integrated oil or CH4 mass source (kg), 'IMSO'
!   264 precipitated salt saturation (null), ' SS '
!   265 hydrate water mole fraction (null), 'XMHW'
!   266 hydrate CO2 mole fraction (null), 'XMHA'
!   267 hydrate CH4 mole fraction (null), 'XMHO'
!   268 plant stomatal resistance (s/m), 'RS_P' or 'RSPX'
!   269 plant stomatal resistance (s/m), 'RSPX'
!   270 plant stomatal resistance (s/m), 'RSPX'
!   271 plant stomatal resistance (s/m), 'RSPX'
!   272 plant stomatal resistance (s/m), 'RSPX'
!   273 fracture adjacent index (null), 'FAI '
!   274 rain-water runoff volumetric rate (m^3/s), 'RWRO'
!   275 well-node pressure (Pa), 'WNP'
!       or MCStress (Pa), 'MCSTR'
!   276 well-node temperature (C), 'WNT'
!   277 differential integrated hydrate water mass (kg), 'DMHW'
!   278 differential integrated hydrate CO2 mass (kg), 'DMHA'
!   279 differential integrated hydrate CH4 mass (kg), 'DMHO'
!   280 differential integrated aqueous CH4 mass (kg), 'DMLO'
!       or differential integrated hydrate N2 mass (kg), 'DMHN'
!   281 differential integrated gas CH4 mass (kg), 'DMGO'
!       or differential integrated mobile water mass (kg), 'DMMW'
!   282 differential integrated water or water mass (kg), 'DMW'
!       or differential integrated mobile CO2 mass (kg), 'DMMA'
!   283 differential integrated air or CO2 mass (kg), 'DMA'
!       or differential integrated mobile CH4 mass (kg), 'DMMO'
!   284 differential integrated oil or CH4 mass (kg), 'DMO'
!       or differential integrated mobile N2 mass (kg), 'DMMN'
!   285 source-well pressure Pa, 'SWP'
!   286 source-well temperature C, 'SWT'
!   287 aqueous spill head/height m, 'HLSP'
!   288 NAPL spill head/height m, 'HNSP'
!   289 gas viscosity, 'VISG'
!   290 NAPL viscosity, 'VISN'
!   291 x node position, 'XP'
!   292 y node position, 'YP'
!   293 z node position, 'ZP'
!   294 gas CH4 mole fraction of formers, 'YMGO'
!   295 hydrate CH4 mole fraction of formers, 'YMHGO'
!   296 aqueous enthalpy, 'HL'
!   297 gas enthalpy, 'HG'
!   298 nonaqueous liquid phase enthalpy, 'HN'
!   299 similitude parameter (m^2/s), 'SIMV'
!   341 aqueous internal energy, 'UL'
!   342 gas internal energy, 'UG'
!   343 nonaqueous liquid phase internal energy, 'UN'
!   344 aqueous thermal conductivity, 'THKL'
!   345 gas thermal conductivity, 'THKG'
!   346 nonaqueous liquid phase thermal conductivity, 'THKN'
!   347 CO2 aqueous diffusion coefficient, 'DFLA' 
!   348 water gas diffusion coefficient, 'DFGW' 
!   349 coupled-well CO2 mass rate (kg/s), 'QMRA'
!       or injection well total petroleum mass rate (kg/s), 'IWPR'
!       or production well total petroleum mass rate (kg/s), 'PWPR'
!   350 coupled-well CO2 mass integral (kg), 'QMIA'
!       or injection well total petroleum mass integral (kg), 'IWPI'
!       or production well total petroleum mass integral (kg), 'PWPI'
!   351 coupled-well water mass rate (kg/s), 'QMRW'
!       or coupled-well mass rate (kg/s), 'QWMR'
!       or injection well water mass rate (kg/s), 'IWWR'
!       or production well water mass rate (kg/s), 'PWWR'
!   352 coupled-well water mass integral (kg), 'QMIW'
!       or coupled-well mass integral (kg), 'QWMI'
!       or injection well water mass integral (kg), IWWI
!       or production well water mass integral (kg), PWWI
!   353 coupled-well air mass rate (kg/s), 'QMRW'
!       or injection well air mass rate (kg/s), 'IWWR'
!       or production well air mass rate (kg/s), 'PWWR'
!   354 coupled-well air mass integral (kg) 'QMIW'
!       or injection well air mass integral (kg), IWWI
!       or production well air mass integral (kg), PWWI
!   361 vertically-integrated CO2 mass (kg), VIA
!       or differential integrated N2 gas mass (kg), DMGN
!   362 vertically-integrated CO2 mass, (kg/m^2), VIAPA
!   363 vertically-integrated CO2 mass (kg), VIGA
!       or differential integrated water aqueous mass (kg), DMLW
!   364 vertically-integrated CO2 mass, kg/m^2, VIGAPA
!   365 vertically-integrated CO2 mass (kg), VILA
!   366 vertically-integrated CO2 mass, kg/m^2, VILAPA
!   367 integrated precipitated salt mass (kg), IMPS
!   368 mean effective stress, Pa, STRS-M
!   369 strain xx, EPSXX
!   370 strain yy, EPSYY
!   371 strain zz, EPSZZ
!   372 strain yz, EPSYZ
!   373 strain xz, EPSXZ
!   374 strain xy, EPSXY
!   375 x displacement, m, DISPX
!   376 y displacement, m, DISPY
!   377 z displacement, m, DISPZ
!   378 integrated energy, J, 'IQ'
!   379 gas CO2 mole fraction of formers, 'YMGA'
!   380 hydrate-gas CO2 mole fraction of former,s 'YMHGA'
!   381 gas N2 mole fraction of formers, 'YMGN'
!   382 hydrate-gas N2 mole fraction of formers, 'YMHGN'
!   383 total nonaqueous CO2 mole fraction of formers, 'ZMCA'
!       or total petroleum component mole fraction, 'ZMC'
!   384 total nonaqueous CH4 mole fraction of formers, 'ZMCO'
!   385 total nonaqueous N2 mole fraction of formers, 'ZMCN'
!   386 hydrate equilibrium pressure, 'PEQH'
!   387 total mobile formers vapor pressure, 'PVTF'
!   388 x-direction index of node, 'I'
!   389 y-direction index of node, 'J'
!   390 z-direction index of node, 'K'
!   391 xnc surface area, m^2 'AFX'
!   392 ync surface area, m^2 'AFY'
!   393 znc surface area, m^2 'AFZ'
!   394 saturation function index, 'SFZN'
!   395 processor id, 'PID'
!   396 gas-nonaqueous liquid interfacial tension, N/m 'IFT'
!   397 fracture mechanical aperture (m), 'APM'
!   398 fracture hydraulic aperture (m), 'APH'
!   399 distance along borehole (m), 'DABH'
!   401 solute volumetric concentration 'C   '
!   401 total_species volumetric concentration, mol/m^3, 'SP  '
!   402 solute aqueous concentration 'CL  '
!   402 total_species aqueous concentration 'SPL  '
!   403 solute gas concentration 'CG  '
!   403 total_species gas concentration 'SPG  '
!   404 solute NAPL concentration 'CN  '
!   404 total_species nonaqueus-liquid concentration 'SPN  '
!   405 solute aqueous mole fraction 'YL  '
!   406 solute gas mole fraction 'YG  '
!   407 solute NAPL mole fraction 'YN  '
!   408 xnc solute flux 'UC  '
!   409 ync solute flux 'VC  '
!   410 znc solute flux 'WC  '
!   411 solute source 'SRC '
!   411 total_species source 'SPSR'
!   412 exchange total_species volumetric concentration, mol/m^3, 'SPX '
!   413 exchange total_species aqueous concentration mol/m^3, 'SPLX'
!   414 exchange total_species gas concentration mol/m^3, 'SPGX'
!   415 exchange total_species nonaqueus-liquid concentration 'SPNX'
!   423 solute integrated mass 'ICM'
!   423 total_species integrated mass 'SPIM'
!   426 exchange total_species volumetric concentration, mol/m^3, 'SPS '
!   427 exchange total_species aqueous concentration mol/m^3, 'SPLS'
!   428 exchange total_species gas concentration mol/m^3, 'SPGS'
!   429 exchange total_species nonaqueus-liquid concentration 'SPNS'
!   430 solute mass concentration, 1/kg soil, 'CM'
!   431 solute activity, Bq, 'CA'
!   432 integrated solute activity, Bq, 'CIA'
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, PNL, February, 1993.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_FRC
      USE PARM_BH
      USE OUTPU
      USE LEAK_WELL
      USE GRID
      USE GEOM_FRC
      USE GEOM_BH
      USE GEO_MECH
      USE FLUXN
      USE FILES
      USE FDVP
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      EXTERNAL ICOUNT
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(:), ALLOCATABLE ::  VSKP
      INTEGER, DIMENSION(:), ALLOCATABLE ::  JSKP
      CHARACTER*12 CHREFX
      CHARACTER*7  FORM1
      CHARACTER*27 FORM2
      CHARACTER*9  FORM3
      CHARACTER*12 FORM4
      CHARACTER*28 FORM5
      CHARACTER*14 FORM6
      CHARACTER*16 FORM7
      CHARACTER*14 FORM8
      CHARACTER*5  FORM9
      CHARACTER*28 FORM10
      CHARACTER*43 FORM11
      CHARACTER*39 FORM12
      CHARACTER*49 FORM13
      CHARACTER*12 FORM14
      CHARACTER*40 FORM15
      CHARACTER*14 FORM16
      CHARACTER*4  FORM17
      CHARACTER*14 FORM18
      CHARACTER*5  FORM19
      CHARACTER*10 FORM20
      CHARACTER*6 FORM21
      CHARACTER*42 FORM22
      CHARACTER*33 FORM23
      CHARACTER*49 FORM24
      CHARACTER*22 FORM25
      CHARACTER*38 FORM26
      CHARACTER*48 FORM27
      CHARACTER*42 FORM28
      CHARACTER*20 FORM29
      CHARACTER*64 SPACES
      CHARACTER*64 TXTX,STRX
      CHARACTER*9 TEXT
      CHARACTER*4096 N_IWR,N_ISC
!
!----------------------Data Statements---------------------------------!
!
      SAVE FORM1,FORM2,FORM3,FORM4,FORM5,FORM6,FORM7,FORM8,FORM9
      SAVE FORM10,FORM11,FORM12,FORM13,FORM14,FORM15,FORM16,FORM17
      SAVE FORM18,FORM19,FORM20,FORM21,FORM22,FORM23,FORM24,FORM25
      SAVE FORM26,FORM27,FORM28,FORM29,SPACES
      DATA FORM1 /'(/,A,$)'/
      DATA FORM2 /'(3X,A,6X,A,5X,A,7X,A,2X,$)'/
      DATA FORM3 /'(3X,2A,$)'/
      DATA FORM4 /'(4X,A4,A,$)'/
      DATA FORM5 /'(A,I3,A,I3,A,I3,A,I8,A,$)'/
      DATA FORM6 /'(1X,1PE10.3,$)'/
      DATA FORM7 /'(20X,3A,5X,4A,$)'/
      DATA FORM8 /'(1X,A,A7,2A,$)'/
      DATA FORM9 /'(A,$)'/
      DATA FORM10 /'(3X,A,6X,A,5X,A,7X,A,3X,A,$)'/
      DATA FORM11 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I2,$)'/
      DATA FORM12 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,A,$)'/
      DATA FORM13 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I2,1X,I6,$)'/
      DATA FORM14 /'(4X,A4,A,$)'/
      DATA FORM15 /'(3X,A,6X,A,5X,A,7X,A,3X,A,1X,A,$)'/
      DATA FORM16 /'(1X,1PE10.3,$)'/
      DATA FORM17 /'(1X)'/
      DATA FORM18 /'(1X,A,A7,2A,$)'/
      DATA FORM19 /'(A,$)'/
      DATA FORM20 /'(10X,2A,$)'/
      DATA FORM21 /'(I1)'/
      DATA FORM22 /'(A,I3,A,I3,A,I3,A,I3,A,I3,A,I3,A,I8,A,$)'/
      DATA FORM23 /'(3X,A,6X,A,5X,A,7X,A,4X,A,2X,A,$)'/
      DATA FORM24 /'(1X,I6,1X,I9,1X,1PE12.5,1X,1PE12.5,1X,I3,2X,I3,$)'/
      DATA FORM25 /'(A,I4,A,I4,A,I8,A,$)'/
      DATA FORM26 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,A,$)'/
      DATA FORM27 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,1X,I2,2X,I2,$)'/
      DATA FORM28 /'(1X,I6,1X,A,1X,1PE12.5,1X,1PE12.5,1X,I2,$)'/
      DATA FORM29 /'(A,I4,A,I8,A,I8,A,$)'/
      DATA SPACES /'                                                    
     &      '/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REFNOD'
      ALLOCATE( VSKP(1:LVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: VSKP'
        CALL WRMSGS( INDX )
      ENDIF
      ALLOCATE( JSKP(1:LVREF),STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Allocation Error: JSKP'
        CALL WRMSGS( INDX )
      ENDIF
      DO 2 N = 1,LVREF
        JSKP(N) = 0
        VSKP(N) = 0.D+0
    2 CONTINUE
!
!---  Compute format spacings according to the number of significant
!     digits requested  ---
!
      IVAR = (ISGNS+3)/2
      WRITE(FORM4(2:2),'(I1)') MAX( 4,ISGNS+3-IVAR )
      WRITE(FORM6(8:9),'(I2)') MAX( 10,ISGNS+6 )
      WRITE(FORM6(11:11),'(I1)') MIN( 9,ISGNS-1 )
      IVAR = (ISGNS-2)/2
      WRITE(FORM8(2:2),'(I1)') MAX( 1,IVAR )
      IVAR = (ISGNO+3)/2
      WRITE(FORM14(2:2),'(I1)') MAX( 4,ISGNO+3-IVAR )
      WRITE(FORM16(8:9),'(I2)') MAX( 10,ISGNO+6 )
      WRITE(FORM16(11:11),'(I1)') MIN( 9,ISGNO-1 )
      IVAR = (ISGNO-2)/2
      WRITE(FORM18(2:2),'(I1)') MAX( 1,ISGNO-2-IVAR )
!
!---  Write main header line  ---
!
      IF( ICNO.EQ.-1 ) THEN
        WRITE(IWR,'(//,A)') ' ---  Reference Node Output Record  ---'
        ICNO = 10
      ENDIF
      IF( ICNS.EQ.-1 ) THEN
        WRITE(ISC,'(//,A)') ' ---  Reference Node Output Record  ---'
        ICNS = 10
      ENDIF
!
!---  Write subheader lines to the output file ---
!
      IF( ICNO.EQ.10 ) THEN
        ICNO = 0
        WRITE(IWR,FORM1) 'Reference Node(s)'
        DO 10 M = 1,NREF
          IF( NDREF(M).LT.0 ) THEN
            NWN = -NDREF(M)
            N = ND_LW(NWN)
          ELSE
            N = NDREF(M)
          ENDIF
 !
 !---     Block refined node, 
 !        print additional reference node indices  ---
 !
          IF( N.GT.NFLD .AND. NDREF(M).GT.0 ) THEN
            IRX = IBR(1,N)
            JRX = IBR(2,N)
            KRX = IBR(3,N)
            I = ID(INP(N-NFLD))
            J = JD(INP(N-NFLD))
            K = KD(INP(N-NFLD))
            WRITE(FORM22(5:5),'(I1)') ICOUNT(I)
            WRITE(FORM22(10:10),'(I1)') ICOUNT(J)
            WRITE(FORM22(15:15),'(I1)') ICOUNT(K)
            WRITE(FORM22(20:20),'(I1)') ICOUNT(IRX)
            WRITE(FORM22(25:25),'(I1)') ICOUNT(JRX)
            WRITE(FORM22(30:30),'(I1)') ICOUNT(KRX)
            WRITE(FORM22(35:35),'(I1)') ICOUNT(N)
            WRITE(IWR,FORM22) ' (',I,',',J,',',K,':',
     &        IRX,',',JRX,',',KRX,':',N,')'
 !
 !---     Leaky well node  ---
 !
          ELSEIF( N.GT.NFLD .AND. NDREF(M).LT.0 ) THEN
            NLWX = 0
            DO NLW = 1,N_LW
              DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
                IF( N.EQ.ND_LW(NWN) ) THEN
                  NLWX = NLW
                  NWNX = NWN
                ENDIF
              ENDDO
            ENDDO
            WRITE(FORM29(5:5),'(I1)') ICOUNT(NLWX)
            WRITE(FORM29(10:10),'(I1)') ICOUNT(NWNX)
            WRITE(FORM29(15:15),'(I1)') ICOUNT(N)
            WRITE(IWR,FORM29) ' (LW',NLWX,',WN',NWNX,':',N,')'           
          ELSE
            I = ID(N)
            J = JD(N)
            K = KD(N)
            WRITE(FORM5(5:5),'(I1)') ICOUNT(I)
            WRITE(FORM5(10:10),'(I1)') ICOUNT(J)
            WRITE(FORM5(15:15),'(I1)') ICOUNT(K)
            WRITE(FORM5(20:20),'(I1)') ICOUNT(N)
            WRITE(IWR,FORM5) ' (',I,',',J,',',K,':',N,')'
          ENDIF
   10   CONTINUE
        WRITE(IWR,FORM17)
!
!---    Fracture triangles  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 .AND. NF_FRC.GT.0 ) THEN
          WRITE(IWR,FORM9) 'Reference Fracture Triangle(s)'
          DO M = 1,NREF_FRC
            NTX = NDREF_FRC(M)
            DO NFX = 1,NF_FRC
              IF( NTX.GE.IP_FRC(1,NFX) .AND. NTX.LE.IP_FRC(2,NFX) ) EXIT
            ENDDO
            NLX = NTX - IP_FRC(1,NFX) + 1
            WRITE(FORM25(5:5),'(I1)') ICOUNT(NFX)
            WRITE(FORM25(10:10),'(I1)') ICOUNT(NLX)
            WRITE(FORM25(15:15),'(I1)') ICOUNT(NTX)
            WRITE(IWR,FORM25) ' (F',NFX,':L',NLX,':T',NTX,')'
          ENDDO
          WRITE(IWR,FORM17)
        ENDIF
!
!---    Borehole nodes  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 .AND. N_BH.GT.0 ) THEN
          WRITE(IWR,FORM9) 'Reference Borehole Node(s)'
          DO M = 1,NREF_BH
            NBN = NDREF_BH(M)
            DO NBH = 1,N_BH
              IF( NBN.GE.ID_BH(3,NBH) .AND. NBN.LE.ID_BH(4,NBH) ) EXIT
            ENDDO
            NBL = NBN - ID_BH(3,NBH) + 1
            WRITE(FORM25(5:5),'(I1)') ICOUNT(NBH)
            WRITE(FORM25(10:10),'(I1)') ICOUNT(NBL)
            WRITE(FORM25(15:15),'(I1)') ICOUNT(NBN)
            WRITE(IWR,FORM25) ' (B',NBH,':L',NBL,':N',NBN,')'
          ENDDO
          WRITE(IWR,FORM17)
        ENDIF
        IN_IWR = 0
        IF( ISLC(17).NE.0 ) THEN
          WRITE(IWR,FORM2) 'Step','Node','Time','Timestep'
          INX = 1
          WRITE(N_IWR(INX:INX+44),'(A)') SPACES(1:45)
          INX = INX+45
          TXTX = 'Itr'
          ITX = 3
          ISX = 14
          CALL CNTRTXT( TXTX,STRX,ITX,ISX )
          WRITE(IWR,'(A,$)') STRX(1:ISX)
          WRITE(N_IWR(INX:INX+ISX-1),'(A)') SPACES(1:ISX)
          INX = INX + ISX
        ELSEIF( ISLC(50).NE.0 ) THEN
          WRITE(IWR,FORM23) 'Step','Node','Time','Timestep','Itr','MLp'
          INX = 1
          WRITE(N_IWR(INX:INX+51),'(A)') SPACES(1:52)
          INX = INX+54
        ELSE
          WRITE(IWR,FORM10) 'Step','Node','Time','Timestep','Itr'
          INX = 1
          WRITE(N_IWR(INX:INX+46),'(A)') SPACES(1:47)
          INX = INX+47
        ENDIF
        DO 20 NV = 1,NVREF
          IRNV = IREF(NV)
          IRNV_CW = IREF_CW(NV)
          IRNV_GC = IREFGC(NV)
          ISX = ISGNO+7
          CHREFX = '            '
          CHREFX(1:6) = CHREF(IRNV)
!
!---      Insert well number after the letter 'W'  ---
!
          IF( IRNV_CW.GT.0 ) THEN
            ICX = ICOUNT(IRNV_CW)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            IWX = INDEX(CHREFX(1:),'W') + 1
            CHREFX(IWX+ICX:ICH+ICX) = CHREFX(IWX:ICH)
            WRITE(CHREFX(IWX:IWX+ICX-1),FORM21) IRNV_CW
          ENDIF
!
!---      Insert component number after the letter 'C'  ---
!
          IF( IRNV_GC.GT.0 ) THEN
            ICX = ICOUNT(IRNV_GC)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            IWX = INDEX(CHREFX(1:),'C') + 1
            CHREFX(IWX+ICX:ICH+ICX) = CHREFX(IWX:ICH)
            WRITE(CHREFX(IWX:IWX+ICX-1),FORM21) IRNV_GC
          ENDIF
!
!---      Insert reactive species number after reference
!         character string  ---
!
          INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
          IF( IRNV.GT.INDX ) THEN
            IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
              NSP = ((IRNV-INDX)/33)
            ELSE
              NSP = ((IRNV-INDX)/33) + 1
            ENDIF
!
!---        Skip for pH output  ---
!
            IF( MOD((IRNV-INDX),33).NE.27 ) THEN
              ICX = ICOUNT(NSP)
              WRITE(FORM21(3:3),'(I1)') ICX
              ICH = INDEX(CHREFX(1:),' ') - 1
              WRITE(CHREFX(ICH+1:ICH+ICX),FORM21) NSP
            ENDIF
!
!---      Insert solutes, conservation component species, and 
!         kinetic component species number after reference
!         character string  ---
!
          ELSEIF( IRNV.GT.400 ) THEN
            IF( MOD((IRNV-400),33).EQ.0 ) THEN
              NSL = ((IRNV-400)/33)
            ELSE
              NSL = ((IRNV-400)/33) + 1
            ENDIF
            ICX = ICOUNT(NSL)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            WRITE(CHREFX(ICH+1:ICH+ICX),FORM21) NSL
          ENDIF 
          TXTX = CHREFX
          ITX = INDEX(TXTX(1:),' ') - 1
          CALL CNTRTXT( TXTX,STRX,ITX,ISX )
          WRITE(IWR,'(A,$)') STRX(1:ISX)
!
!---      Reactive species  ---
!
          INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
          IF( IRNV.GT.INDX ) THEN
            IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
              NSP = ((IRNV-INDX)/33)
            ELSE
              NSP = ((IRNV-INDX)/33) + 1
            ENDIF
            IF( NSP.GT.NSPL+NSPS+NSPE ) THEN
              NCH = INDEX(SPNMG(NSP-NSPL-NSPS-NSPE)(1:),'  ')-1 
              TXTX = SPNMG(NSP-NSPL-NSPS-NSPE)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSEIF( NSP.GT.NSPL+NSPS ) THEN
              NCH = INDEX(SPNME(NSP-NSPL-NSPS)(1:),'  ')-1 
              TXTX = SPNME(NSP-NSPL-NSPS)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSEIF( NSP.GT.NSPL ) THEN
              NCH = INDEX(SPNMS(NSP-NSPL)(1:),'  ')-1 
              TXTX = SPNMS(NSP-NSPL)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSE
              NCH = INDEX(SPNML(NSP)(1:),'  ')-1 
              TXTX = SPNML(NSP)(1:NCH)
              ITX = MIN(NCH,ISX)
            ENDIF
            CALL CNTRTXT( TXTX,STRX,ITX,ISX )
            WRITE(N_IWR(INX:INX+ISX-1),'(A)') STRX(1:ISX)
            INX = INX+ISX
            IN_IWR = 1
!
!---      Solutes, conservation component species, and 
!         kinetic component species  ---
!
          ELSEIF( IRNV.GT.400 ) THEN
            IF( MOD((IRNV-400),33).EQ.0 ) THEN
              NSL = ((IRNV-400)/33)
            ELSE
              NSL = ((IRNV-400)/33) + 1
            ENDIF
            NCH = INDEX(SOLUT(NSL)(1:),'  ')-1 
            TXTX = SOLUT(NSL)(1:NCH)
            ITX = MIN(NCH,ISX)
            CALL CNTRTXT( TXTX,STRX,ITX,ISX )
            WRITE(N_IWR(INX:INX+ISX-1),'(A)') STRX(1:ISX)
            INX = INX+ISX
            IN_IWR = 1
          ELSE
            WRITE(N_IWR(INX:INX+ISX-1),'(A)') SPACES(1:ISX)
            INX = INX+ISX
          ENDIF
   20   CONTINUE
        WRITE(IWR,FORM17)
        IF( IN_IWR.EQ.1 ) WRITE(IWR,'(A)') N_IWR(1:INX-1)
        IF( ISLC(17).NE.0 ) THEN
          IC = 14
          WRITE(IWR,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:IC)
        ELSEIF( ISLC(50).NE.0 ) THEN
          WRITE(IWR,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:10)
        ELSE
          WRITE(IWR,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:6)
        ENDIF
        DO 30 NV = 1,NVREF
          IRNV = IREF(NV)
          ITX = 9
          IF( UNREF(IRNV).NE.'null' ) THEN
             TXTX = '[' // UNREF(IRNV)(1:7) // ']'
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(IWR,'(A,$)') STRX(1:ISX)
          ELSEIF( IREFGC(NV).NE.0 ) THEN
             TXTX = '[' // GCNM(IREFGC(NV))(1:7) // ']'
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(IWR,'(A,$)') STRX(1:ISX)
          ELSE
             TXTX = SPACES(1:9)
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(IWR,'(A,$)') STRX(1:ISX)
          ENDIF
   30   CONTINUE
        WRITE(IWR,FORM17)
      ENDIF
!
!---  Write subheader lines to the screen ---
!
      IF( ICNS.EQ.10 ) THEN
        ICNS = 0
        WRITE(ISC,FORM1) 'Reference Node(s)'
        DO 40 M = 1,NREF
          IF( NDREF(M).LT.0 ) THEN
            NWN = -NDREF(M)
            N = ND_LW(NWN)
          ELSE
            N = NDREF(M)
          ENDIF
 !
 !---     Block refined node, 
 !        print additional reference node indices  ---
 !
          IF( N.GT.NFLD .AND. NDREF(M).GT.0 ) THEN
            IRX = IBR(1,N)
            JRX = IBR(2,N)
            KRX = IBR(3,N)
            I = ID(INP(N-NFLD))
            J = JD(INP(N-NFLD))
            K = KD(INP(N-NFLD))
            WRITE(FORM22(5:5),'(I1)') ICOUNT(I)
            WRITE(FORM22(10:10),'(I1)') ICOUNT(J)
            WRITE(FORM22(15:15),'(I1)') ICOUNT(K)
            WRITE(FORM22(20:20),'(I1)') ICOUNT(IRX)
            WRITE(FORM22(25:25),'(I1)') ICOUNT(JRX)
            WRITE(FORM22(30:30),'(I1)') ICOUNT(KRX)
            WRITE(FORM22(35:35),'(I1)') ICOUNT(N)
            WRITE(ISC,FORM22) ' (',I,',',J,',',K,':',
     &        IRX,',',JRX,',',KRX,':',N,')'
 !
 !---     Leaky well node  ---
 !
          ELSEIF( N.GT.NFLD .AND. NDREF(M).LT.0 ) THEN
            NLWX = 0
            DO NLW = 1,N_LW
              DO NWN = ID_LW(3,NLW),ID_LW(4,NLW)
                IF( N.EQ.ND_LW(NWN) ) THEN
                  NLWX = NLW
                  NWNX = NWN
                ENDIF
              ENDDO
            ENDDO
            WRITE(FORM29(5:5),'(I1)') ICOUNT(NLWX)
            WRITE(FORM29(10:10),'(I1)') ICOUNT(NWNX)
            WRITE(FORM29(15:15),'(I1)') ICOUNT(N)
            WRITE(ISC,FORM29) ' (LW',NLWX,',WN',NWNX,':',N,')'           
          ELSE
            I = ID(N)
            J = JD(N)
            K = KD(N)
            WRITE(FORM5(5:5),'(I1)') ICOUNT(I)
            WRITE(FORM5(10:10),'(I1)') ICOUNT(J)
            WRITE(FORM5(15:15),'(I1)') ICOUNT(K)
            WRITE(FORM5(20:20),'(I1)') ICOUNT(N)
            WRITE(ISC,FORM5) ' (',I,',',J,',',K,':',N,')'
          ENDIF
   40   CONTINUE
        IN_ISC = 0
        WRITE(ISC,FORM17)
!
!---    Fracture triangles  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 .AND. NF_FRC.GT.0 ) THEN
          WRITE(ISC,FORM9) 'Reference Fracture Triangles(s)'
          DO M = 1,NREF_FRC
            NTX = NDREF_FRC(M)
            DO NFX = 1,NF_FRC
              IF( NTX.GE.IP_FRC(1,NFX) .AND. NTX.LE.IP_FRC(2,NFX) ) EXIT
            ENDDO
            NLX = NTX - IP_FRC(1,NFX) + 1
            WRITE(FORM25(5:5),'(I1)') ICOUNT(NFX)
            WRITE(FORM25(10:10),'(I1)') ICOUNT(NLX)
            WRITE(FORM25(15:15),'(I1)') ICOUNT(NTX)
            WRITE(ISC,FORM25) ' (F',NFX,':L',NLX,':T',NTX,')'
          ENDDO
          WRITE(ISC,FORM17)
        ENDIF
!
!---    Borehole nodes  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 .AND. N_BH.GT.0 ) THEN
          WRITE(ISC,FORM9) 'Reference Borehole Node(s)'
          DO M = 1,NREF_BH
            NBN = NDREF_BH(M)
            DO NBH = 1,N_BH
              IF( NBN.GE.ID_BH(3,NBH) .AND. NBN.LE.ID_BH(4,NBH) ) EXIT
            ENDDO
            NBL = NBN - ID_BH(3,NBH) + 1
            WRITE(FORM25(5:5),'(I1)') ICOUNT(NBH)
            WRITE(FORM25(10:10),'(I1)') ICOUNT(NBL)
            WRITE(FORM25(15:15),'(I1)') ICOUNT(NBN)
            WRITE(ISC,FORM25) ' (B',NBH,':L',NBL,':N',NBN,')'
          ENDDO
          WRITE(ISC,FORM17)
        ENDIF
        IF( ISLC(17).NE.0 ) THEN
          WRITE(ISC,FORM2) 'Step','Node','Time','Timestep'
          INX = 1
          WRITE(N_ISC(INX:INX+41),'(A)') SPACES(1:42)
          INX = INX+42
          TXTX = 'Itr'
          ITX = 3
          ISX = 14
          CALL CNTRTXT( TXTX,STRX,ITX,ISX )
          WRITE(ISC,'(A,$)') STRX(1:ISX)
          WRITE(N_ISC(INX:INX+ISX-1),'(A)') SPACES(1:ISX)
          INX = INX + ISX
        ELSEIF( ISLC(50).NE.0 ) THEN
          WRITE(ISC,FORM23) 'Step','Node','Time','Timestep','Itr','MLp'
          INX = 1
          WRITE(N_ISC(INX:INX+48),'(A)') SPACES(1:49)
          INX = INX+51
        ELSE
          WRITE(ISC,FORM10) 'Step','Node','Time','Timestep','Itr'
          INX = 1
          WRITE(N_ISC(INX:INX+43),'(A)') SPACES(1:44)
          INX = INX+44
        ENDIF
        DO 50 NV = 1,NVREF
          IRNV = IREF(NV)
          IRNV_CW = IREF_CW(NV)
          IRNV_GC = IREFGC(NV)
          ISX = ISGNS+7
          CHREFX = '            '
          CHREFX(1:6) = CHREF(IRNV)
!
!---      Insert well number after the letter 'W'  ---
!
          IF( IRNV_CW.GT.0 ) THEN
            ICX = ICOUNT(IRNV_CW)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            IWX = INDEX(CHREFX(1:),'W') + 1
            CHREFX(IWX+ICX:ICH+ICX) = CHREFX(IWX:ICH)
            WRITE(CHREFX(IWX:IWX+ICX-1),FORM21) IRNV_CW
          ENDIF
!
!---      Insert component number after the letter 'C'  ---
!
          IF( IRNV_GC.GT.0 ) THEN
            ICX = ICOUNT(IRNV_GC)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            IWX = INDEX(CHREFX(1:),'C') + 1
            CHREFX(IWX+ICX:ICH+ICX) = CHREFX(IWX:ICH)
            WRITE(CHREFX(IWX:IWX+ICX-1),FORM21) IRNV_GC
          ENDIF
!
!---      Insert reactive species number after reference
!         character string  ---
!
          INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
          IF( IRNV.GT.INDX ) THEN
            IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
              NSP = ((IRNV-INDX)/33)
            ELSE
              NSP = ((IRNV-INDX)/33) + 1
            ENDIF
!
!---        Skip for pH output  ---
!
            IF( MOD((IRNV-INDX),33).NE.27 ) THEN
              ICX = ICOUNT(NSP)
              WRITE(FORM21(3:3),'(I1)') ICX
              ICH = INDEX(CHREFX(1:),' ') - 1
              WRITE(CHREFX(ICH+1:ICH+ICX),FORM21) NSP
            ENDIF
!
!---      Insert solutes, conservation component species, and 
!         kinetic component species number after reference
!         character string  ---
!
          ELSEIF( IRNV.GT.400 ) THEN
            IF( MOD((IRNV-400),33).EQ.0 ) THEN
              NSL = ((IRNV-400)/33)
            ELSE
              NSL = ((IRNV-400)/33) + 1
            ENDIF
            ICX = ICOUNT(NSL)
            WRITE(FORM21(3:3),'(I1)') ICX
            ICH = INDEX(CHREFX(1:),' ') - 1
            WRITE(CHREFX(ICH+1:ICH+ICX),FORM21) NSL
          ENDIF 
          TXTX = CHREFX
          ITX = INDEX(TXTX(1:),' ') - 1
          CALL CNTRTXT( TXTX,STRX,ITX,ISX )
          WRITE(ISC,'(A,$)') STRX(1:ISX)
!
!---      Reactive species  ---
!
          INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
          IF( IRNV.GT.INDX ) THEN
            IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
              NSP = ((IRNV-INDX)/33)
            ELSE
              NSP = ((IRNV-INDX)/33) + 1
            ENDIF
            IF( NSP.GT.NSPL+NSPS+NSPE ) THEN
              NCH = INDEX(SPNMG(NSP-NSPL-NSPS-NSPE)(1:),'  ')-1 
              TXTX = SPNMG(NSP-NSPL-NSPS-NSPE)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSEIF( NSP.GT.NSPL+NSPS ) THEN
              NCH = INDEX(SPNME(NSP-NSPL-NSPS)(1:),'  ')-1 
              TXTX = SPNME(NSP-NSPL-NSPS)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSEIF( NSP.GT.NSPL ) THEN
              NCH = INDEX(SPNMS(NSP-NSPL)(1:),'  ')-1 
              TXTX = SPNMS(NSP-NSPL)(1:NCH)
              ITX = MIN(NCH,ISX)
            ELSE
              NCH = INDEX(SPNML(NSP)(1:),'  ')-1 
              TXTX = SPNML(NSP)(1:NCH)
              ITX = MIN(NCH,ISX)
            ENDIF
            CALL CNTRTXT( TXTX,STRX,ITX,ISX )
            WRITE(N_ISC(INX:INX+ISX-1),'(A)') STRX(1:ISX)
            INX = INX+ISX
            IN_ISC = 1
!
!---      Solutes, conservation component species, and 
!         kinetic component species  ---
!
          ELSEIF( IRNV.GT.400 ) THEN
            IF( MOD((IRNV-400),33).EQ.0 ) THEN
              NSL = ((IRNV-400)/33)
            ELSE
              NSL = ((IRNV-400)/33) + 1
            ENDIF
            NCH = INDEX(SOLUT(NSL)(1:),'  ')-1 
            TXTX = SOLUT(NSL)(1:NCH)
            ITX = MIN(NCH,ISX)
            CALL CNTRTXT( TXTX,STRX,ITX,ISX )
            WRITE(N_ISC(INX:INX+ISX-1),'(A)') STRX(1:ISX)
            INX = INX+ISX
            IN_ISC = 1
          ELSE
            WRITE(N_ISC(INX:INX+ISX-1),'(A)') SPACES(1:ISX)
            INX = INX+ISX
          ENDIF
   50   CONTINUE
        WRITE(ISC,FORM17)
        IF( IN_ISC.EQ.1 ) WRITE(ISC,'(A)') N_ISC(1:INX-1)
        IF( ISLC(17).NE.0 ) THEN
          IC = 14
          WRITE(ISC,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:IC)
        ELSEIF( ISLC(50).NE.0 ) THEN
          WRITE(ISC,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:10)
        ELSE
          WRITE(ISC,FORM7) '[',UNTM(1:6),']','[',UNTM(1:6),']',
     &      SPACES(1:6)
        ENDIF
        DO 60 NV = 1,NVREF
          IRNV = IREF(NV)
          ITX = 9
          IF( UNREF(IRNV).NE.'null' ) THEN
             TXTX = '[' // UNREF(IRNV)(1:7) // ']'
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(ISC,'(A,$)') STRX(1:ISX)
          ELSEIF( IREFGC(NV).NE.0 ) THEN
             TXTX = '[' // GCNM(IREFGC(NV))(1:7) // ']'
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(ISC,'(A,$)') STRX(1:ISX)
          ELSE
             TXTX = SPACES(1:9)
             CALL CNTRTXT( TXTX,STRX,ITX,ISX )
             WRITE(ISC,'(A,$)') STRX(1:ISX)
          ENDIF
   60   CONTINUE
        WRITE(ISC,FORM17)
      ENDIF
!
!---  Write reference number information  ---
!
      IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 .OR.
     &  MOD( (NSTEP-NRST),IFQS ).EQ.0 ) THEN
        NINAC = NFLD-NXP
        DO 900 M = 1,NREF
          IF( NDREF(M).LT.0 ) THEN
            NWN = -NDREF(M)
            N = ND_LW(NWN)
          ELSE
            N = NDREF(M)
          ENDIF
          TMY = TM
          DTY = DT
          IF( UNTM.NE.'null' ) THEN
            INDX = 4
            IUNS = 1
            CALL RDUNIT(UNTM,TMY,INDX)
            IUNS = 1
            CALL RDUNIT(UNTM,DTY,INDX)
          ENDIF
!
!---      Output file  ---
!
          IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) THEN
!
!---        Courant number control output  ---
!
            IF( ISLC(17).NE.0 ) THEN
              WRITE(IWR,FORM12) NSTEP,N,TMY,DTY,' ['
              WRITE(IWR,'(I4,A,$)') NITER,'/'
!              NSOLUX = NSOLU
!              IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!              DO 62 MM = 1,NSOLUX-1
!                WRITE(IWR,'(I3,A,$)') N_CRN(MM),'/'
!   62         CONTINUE
!              WRITE(IWR,'(I3,A,$)') N_CRN(NSOLUX),']'
              WRITE(IWR,'(I4,A,$)') N_CRN(1),']'
            ELSEIF( ISLC(50).NE.0 ) THEN
              WRITE(IWR,FORM24) NSTEP,N,TMY,DTY,K_GM(2),K_GM(1)
            ELSE
              WRITE(IWR,FORM11) NSTEP,N,TMY,DTY,NITER
            ENDIF
          ENDIF
!
!---      Standard output (screen)  ---
!
          IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) THEN
!
!---        Courant number control output  ---
!
            IF( ISLC(17).NE.0 ) THEN
              WRITE(ISC,FORM12) NSTEP,N,TMY,DTY,' ['
              WRITE(ISC,'(I4,A,$)') NITER,'/'
!              NSOLUX = NSOLU
!              IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!              DO 64 MM = 1,NSOLUX-1
!                WRITE(ISC,'(I3,A,$)') N_CRN(MM),'/'
!   64         CONTINUE
!              WRITE(ISC,'(I3,A,$)') N_CRN(NSOLUX),']'
              WRITE(ISC,'(I4,A,$)') N_CRN(1),']'
            ELSEIF( ISLC(50).NE.0 ) THEN
              WRITE(ISC,FORM24) NSTEP,N,TMY,DTY,K_GM(2),K_GM(1)
            ELSE
              WRITE(ISC,FORM11) NSTEP,N,TMY,DTY,NITER
            ENDIF
          ENDIF
          DO 800 NV = 1,NVREF
            IRNV = IREF(NV)
            IRNV_GC = IREFGC(NV)
            IRNV_CW = IREF_CW(NV)
            IF( M.EQ.1 ) JSKP(NV) = 0
            CALL REFVAR( VAR,VSKP(NV),IRNV,IRNV_GC,IRNV_CW,JSKP(NV),N )
!
!---        Skip unit conversion for fracture-only outputs  ---
!
            IF( IRNV.GT.0 ) THEN
              IF( UNREF(IRNV).NE.'null' ) THEN
                INDX = 4
                CALL RDUNIT(UNREF(IRNV),VAR,INDX)
              ENDIF
            ENDIF
            IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM16) VAR
            IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM6) VAR
  800     CONTINUE
          IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM17)
          IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM17)
  900   CONTINUE
!
!---    Fracture triangle output  ---
!
        IF( ISLC(74).EQ.1 .OR. ISLC(74).EQ.3 .AND. NF_FRC.GT.0 ) THEN
          DO M = 1,NREF_FRC
            NTX = NDREF_FRC(M)
            TMY = TM
            DTY = DT
            IF( UNTM.NE.'null' ) THEN
              INDX = 4
              IUNS = 1
              CALL RDUNIT(UNTM,TMY,INDX)
              IUNS = 1
              CALL RDUNIT(UNTM,DTY,INDX)
            ENDIF
            TEXT(1:9) = SPACES(1:9)
            I1X = 10 - ICOUNT(NTX)
            WRITE(FORM21(3:3),'(I1)') ICOUNT(NTX)
            WRITE(TEXT(I1X:9),FORM21) NTX
            I1X = I1X - 1
            TEXT(I1X:I1X) = 'T'
!
!---        Output file  ---
!
            IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) THEN
!
!---          Courant number control output  ---
!
              IF( ISLC(17).NE.0 ) THEN
                WRITE(IWR,FORM26) NSTEP,TEXT,TMY,DTY,' ['
                WRITE(IWR,'(I4,A,$)') NITER,'/'
!                NSOLUX = NSOLU
!                IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!                DO MM = 1,NSOLUX-1
!                  WRITE(IWR,'(I3,A,$)') N_CRN(MM),'/'
!                ENDDO
!                WRITE(IWR,'(I3,A,$)') N_CRN(NSOLUX),']'
                WRITE(IWR,'(I4,A,$)') N_CRN(1),']'
              ELSEIF( ISLC(50).NE.0 ) THEN
                WRITE(IWR,FORM27) NSTEP,TEXT,TMY,DTY,K_GM(2),
     &            K_GM(1)
              ELSE
                WRITE(IWR,FORM28) NSTEP,TEXT,TMY,DTY,NITER
              ENDIF
            ENDIF
!
!---        Standard output (screen)  ---
!
            IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) THEN
!
!---          Courant number control output  ---
!
              IF( ISLC(17).NE.0 ) THEN
                WRITE(ISC,FORM26) NSTEP,TEXT,TMY,DTY,' ['
                WRITE(ISC,'(I4,A,$)') NITER,'/'
!                NSOLUX = NSOLU
!                IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!                DO MM = 1,NSOLUX-1
!                  WRITE(ISC,'(I3,A,$)') N_CRN(MM),'/'
!                ENDDO
!                WRITE(ISC,'(I3,A,$)') N_CRN(NSOLUX),']'
                WRITE(ISC,'(I4,A,$)') N_CRN(1),']'
              ELSEIF( ISLC(50).NE.0 ) THEN
                WRITE(ISC,FORM27) NSTEP,TEXT,TMY,DTY,K_GM(2),
     &            K_GM(1)
              ELSE
                WRITE(ISC,FORM28) NSTEP,TEXT,TMY,DTY,NITER
              ENDIF
            ENDIF
            DO NV = 1,NVREF
              IRNV = IREF(NV)
              IRNV_GC = IREFGC(NV)
              IF( M.EQ.1 ) JSKP(NV) = 0
              CALL REFVAR_FRC( VAR,VSKP(NV),IRNV,IRNV_GC,JSKP(NV),NTX )
!
!---          Skip unit conversion for non-fracture outputs  ---
!
              IF( IRNV.GT.0 ) THEN
                IF( UNREF(IRNV).NE.'null' ) THEN
                  INDX = 4
                  CALL RDUNIT(UNREF(IRNV),VAR,INDX)
                ENDIF
              ENDIF
              IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM16) VAR
              IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM6) VAR
            ENDDO
            IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM17)
            IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM17)
          ENDDO
        ENDIF
!
!---    Borehole node output  ---
!
        IF( ISLC(74).EQ.2 .OR. ISLC(74).EQ.3 .AND. N_BH.GT.0 ) THEN
          DO M = 1,NREF_BH
            NBN = NDREF_BH(M)
            TMY = TM
            DTY = DT
            IF( UNTM.NE.'null' ) THEN
              INDX = 4
              IUNS = 1
              CALL RDUNIT(UNTM,TMY,INDX)
              IUNS = 1
              CALL RDUNIT(UNTM,DTY,INDX)
            ENDIF
            TEXT(1:9) = SPACES(1:9)
            I1X = 10 - ICOUNT(NBN)
            WRITE(FORM21(3:3),'(I1)') ICOUNT(NBN)
            WRITE(TEXT(I1X:9),FORM21) NBN
            I1X = I1X - 1
            TEXT(I1X:I1X) = 'B'
!
!---        Output file  ---
!
            IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) THEN
!
!---          Courant number control output  ---
!
              IF( ISLC(17).NE.0 ) THEN
                WRITE(IWR,FORM26) NSTEP,TEXT,TMY,DTY,' ['
                WRITE(IWR,'(I4,A,$)') NITER,'/'
!                NSOLUX = NSOLU
!                IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!                DO MM = 1,NSOLUX-1
!                  WRITE(IWR,'(I3,A,$)') N_CRN(MM),'/'
!                ENDDO
!                WRITE(IWR,'(I3,A,$)') N_CRN(NSOLUX),']'
                WRITE(IWR,'(I4,A,$)') N_CRN(1),']'
              ELSEIF( ISLC(50).NE.0 ) THEN
                WRITE(IWR,FORM27) NSTEP,TEXT,TMY,DTY,K_GM(2),
     &            K_GM(1)
              ELSE
                WRITE(IWR,FORM28) NSTEP,TEXT,TMY,DTY,NITER
              ENDIF
            ENDIF
!
!---        Standard output (screen)  ---
!
            IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) THEN
!
!---          Courant number control output  ---
!
              IF( ISLC(17).NE.0 ) THEN
                WRITE(ISC,FORM26) NSTEP,TEXT,TMY,DTY,' ['
                WRITE(ISC,'(I4,A,$)') NITER,'/'
!                NSOLUX = NSOLU
!                IF( NEQC+NEQK.GT.0 ) NSOLUX = NSOLU+1
!                DO MM = 1,NSOLUX-1
!                  WRITE(ISC,'(I3,A,$)') N_CRN(MM),'/'
!                ENDDO
!                WRITE(ISC,'(I3,A,$)') N_CRN(NSOLUX),']'
                WRITE(ISC,'(I4,A,$)') N_CRN(1),']'
              ELSEIF( ISLC(50).NE.0 ) THEN
                WRITE(ISC,FORM27) NSTEP,TEXT,TMY,DTY,K_GM(2),
     &            K_GM(1)
              ELSE
                WRITE(ISC,FORM28) NSTEP,TEXT,TMY,DTY,NITER
              ENDIF
            ENDIF
            DO NV = 1,NVREF
              IRNV = IREF(NV)
              IF( M.EQ.1 ) JSKP(NV) = 0
              CALL REFVAR_BH( VAR,VSKP(NV),IRNV,JSKP(NV),NBN )
!
!---          Skip unit conversion for non-fracture outputs  ---
!
              IF( IRNV.GT.0 ) THEN
                IF( UNREF(IRNV).NE.'null' ) THEN
                  INDX = 4
                  CALL RDUNIT(UNREF(IRNV),VAR,INDX)
                ENDIF
              ENDIF
              IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM16) VAR
              IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM6) VAR
            ENDDO
            IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) WRITE(IWR,FORM17)
            IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) WRITE(ISC,FORM17)
          ENDDO
        ENDIF
      ENDIF
      IF( MOD( (NSTEP-NRST),IFQO ).EQ.0 ) ICNO = ICNO + 1
      IF( MOD( (NSTEP-NRST),IFQS ).EQ.0 ) ICNS = ICNS + 1
      IF( ALLOCATED(VSKP) ) THEN
      DEALLOCATE( VSKP,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: VSKP'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      IF( ALLOCATED(JSKP) ) THEN
      DEALLOCATE( JSKP,STAT=ISTAT )
      IF( ISTAT.NE.0 ) THEN
        INDX = 3
        CHMSG = 'Deallocation Error: JSKP'
        CALL WRMSGS( INDX )
      ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REFNOD group
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REFVAR( VAR,VSKPX,IRNV,IRNV_GC,IRNV_CW,JSKPX,N )
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
!     Reference node variable conversion.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 29 May 2001.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE WELL_FX
      USE WELL_FD
      USE WELL_CL
      USE TRNSPT
      USE SPILL
      USE SOURC
      USE SOLTN
      USE REACT
      USE PORMED
      USE NAPL
      USE HYST
      USE GRID
      USE GEO_MECH
      USE FLUXT
      USE FLUXS
      USE FLUXP
      USE FLUXN
      USE FLUXD
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
      USE COUP_WELL
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Data Statements---------------------------------!
!
      SAVE DMHWX,DMHAX,DMHOX,DMLOX,DMGOX,DMWX,DMAX,DMOX
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REFVAR'
!
!---  Assign surface indices ---
!
      NPX = NSX(N)
      NPY = NSY(N)
      NPZ = NSZ(N)
      NQX = NSX(N)+1
      IF( INBS(4,N).GT.0 ) NQX = INBS(4,N)
      NQY = NSY(N)+IFLD
      IF( INBS(5,N).GT.0 ) NQY = INBS(5,N)
      NQZ = NSZ(N)+IJFLD
      IF( INBS(6,N).GT.0 ) NQZ = INBS(6,N)
!
!---  Convert reference node variable ---
!
      IF( IRNV.LE.32 ) THEN
!
!---  Aqueous pressure (absolute)  ---
!
      IF( IRNV.EQ.1 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PL(2,N) + PATM
!
!---  Gas pressure (absolute)  ---
!
      ELSEIF( IRNV.EQ.2 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PG(2,N) + PATM
!
!---  NAPL or Liquid-CO2 pressure (absolute)  ---
!
      ELSEIF( IRNV.EQ.3 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PN(2,N) + PATM
!
!---  Temperature  ---
!
      ELSEIF( IRNV.EQ.4 ) THEN
        IUNK = 1
        VAR = T(2,N)
!
!---  Phase condition  ---
!
      ELSEIF( IRNV.EQ.5 ) THEN
        VAR = REAL( NPHAZ(2,N) )
!
!---  Aqueous pressure (gauge)  ---
!
      ELSEIF( IRNV.EQ.6 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PL(2,N)
!
!---  Gas pressure (gauge)  ---
!
      ELSEIF( IRNV.EQ.7 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PG(2,N)
!
!---  NAPL or Liquid-CO2 pressure (gauge)  ---
!
      ELSEIF( IRNV.EQ.8 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PN(2,N)
!
!---  Apparent aqueous saturation or
!     hydrate specific area (m^2/m^3)  ---
!
      ELSEIF( IRNV.EQ.9 ) THEN
        VAR = ASL(N)
!
!---  Apparent total-liquid saturation
!     or hydrate dissociation/formation rate (kmol/s m^2)
!     formation is positive and dissociation is negative  ---
!
      ELSEIF( IRNV.EQ.10 ) THEN
        VAR = AST(N)
!
!---  Actual aqueous saturation  ---
!
      ELSEIF( IRNV.EQ.11 ) THEN
        VAR = SL(2,N)
!
!---  Actual gas saturation  ---
!
      ELSEIF( IRNV.EQ.12 ) THEN
        VAR = SG(2,N)
!
!---  Actual NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.13 ) THEN
        VAR = SN(2,N)
!
!---  Actual total-liquid saturation  ---
!
      ELSEIF( IRNV.EQ.14 ) THEN
        VAR = SL(2,N) + SN(2,N)
!
!---  Aqueous moisture content  ---
!
      ELSEIF( IRNV.EQ.15 ) THEN
        VAR = SL(2,N)*PORD(2,N)
!
!---  NAPL moisture content  ---
!
      ELSEIF( IRNV.EQ.16 ) THEN
        VAR = SN(2,N)*PORD(2,N)
!
!---  Total-liquid moisture content  ---
!
      ELSEIF( IRNV.EQ.17 ) THEN
        VAR = (SL(2,N)+SN(2,N))*PORD(2,N)
!
!---  Apparent trapped-NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.18 ) THEN
        VAR = ASNT(N)
!
!---  Apparent trapped-gas saturation  ---
!
      ELSEIF( IRNV.EQ.19 ) THEN
        VAR = ASGT(N)
!
!---  Diffusive porosity  ---
!
      ELSEIF( IRNV.EQ.20 ) THEN
        VAR = PORD(2,N)
!
!---  Gas water mass fraction
!     or volumetric component mass fraction  ---
!
      ELSEIF( IRNV.EQ.21 ) THEN
        VAR = XGW(2,N)
!
!---  Gas air mass fraction
!     or gas component mass fraction  ---
!
      ELSEIF( IRNV.EQ.22 ) THEN
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = XGC(IRNV_GC,2,N)
        ELSE
          VAR = XGA(2,N)
        ENDIF
!
!---  Gas oil mass fraction
!     or NAPL component mass fraction ---
!
      ELSEIF( IRNV.EQ.23 ) THEN
        IF( IOM.EQ.43 ) THEN
          VAR = XNC(IRNV_GC,2,N)
        ELSE
          VAR = XGO(2,N)
        ENDIF
!
!---  Aqueous water mass fraction  ---
!
      ELSEIF( IRNV.EQ.24 ) THEN
        VAR = XLW(2,N)
!
!---  Aqueous air or gas component mass fraction  ---
!
      ELSEIF( IRNV.EQ.25 ) THEN
        IF( IOM.EQ.30 .OR. IOM.EQ.40 ) THEN
          VAR = XLC(IRNV_GC,2,N)
        ELSE
          VAR = XLA(2,N)
        ENDIF
!
!---  Aqueous oil mass fraction  ---
!
      ELSEIF( IRNV.EQ.26 ) THEN
        VAR = XLO(2,N)
!
!---  Aqueous hydraulic head  ---
!
      ELSEIF( IRNV.EQ.27 ) THEN
        IUNM = 1
        VAR = PL(2,N)/RHORL/GRAV + ZP(N)
!
!---  Gas hydraulic head  ---
!
      ELSEIF( IRNV.EQ.28 ) THEN
        IUNM = 1
        VAR = PG(2,N)/RHORL/GRAV + ZP(N)
!
!---  NAPL hydraulic head  ---
!
      ELSEIF( IRNV.EQ.29 ) THEN
        IUNM = 1
        VAR = PN(2,N)/RHORL/GRAV + ZP(N)
!
!---  Rock/soil type  ---
!
      ELSEIF( IRNV.EQ.30 ) THEN
        IF( IXP(N).EQ.0 ) THEN
          VAR = 0.D+0
        ELSE
          VAR = REAL( IZ(N) )
        ENDIF
!
!---  Aqueous relative permeability ---
!
      ELSEIF( IRNV.EQ.31 ) THEN
        VAR = (RKL(1,2,N)*RKL(2,2,N)*RKL(3,2,N))**(1./3.)
!
!---  Gas relative permeability
!     or residual-NAPL saturation ---
!
      ELSEIF( IRNV.EQ.32 ) THEN
        VAR = RKG(2,N)
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.32 .AND. IRNV.LE.64 ) THEN
!
!---  NAPL relative permeability
!     or liquid-CO2 relative permeability  ---
!
      IF( IRNV.EQ.33 ) THEN
        VAR = RKN(2,N)
!
!---  Aqueous density ---
!
      ELSEIF( IRNV.EQ.34 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = RHOL(2,N)
!
!---  Gas density ---
!
      ELSEIF( IRNV.EQ.35 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = RHOG(2,N)
!
!---  NAPL density ---
!
      ELSEIF( IRNV.EQ.36 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = RHON(2,N)
!
!---  Total water mass ---
!
      ELSEIF( IRNV.EQ.37 ) THEN
        IUNKG = 1
        VAR = 0.D+0
!
!---    Aqueous, gas, and NAPL water  ---
!
        VAR = PORD(2,N)*VOL(N)*(XLW(2,N)*SL(2,N)*RHOL(2,N) +
     &    XGW(2,N)*SG(2,N)*RHOG(2,N))
!
!---    NAPL water  ---
!
        IF( INAPL.NE.0 ) THEN
          VAR = VAR + PORD(2,N)*VOL(N)*XNW(2,N)*SN(2,N)*RHON(2,N)
        ENDIF
!
!---    Hydrate water  ---
!
        IF( LFDH.EQ.LFD ) THEN
          VAR = VAR + PORD(2,N)*VOL(N)*XHW(2,N)*SH(2,N)*RHOH(2,N)
        ENDIF
!
!---    Ice water  ---
!
        IF( LFDI.EQ.LFD ) THEN
          VAR = VAR + PORD(2,N)*VOL(N)*SI(2,N)*RHOI(2,N)
        ENDIF
!
!---  Total air mass, 
!     or total CO2 mass, 
!     or total gas-component mass
!     or total component mass  ---
!
      ELSEIF( IRNV.EQ.38 ) THEN
        IUNKG = 1
        VAR = 0.D+0
        IF( IOM.EQ.30 .OR. IOM.EQ.40 ) THEN
          VAR = PORD(2,N)*VOL(N)*(XLC(IRNV_GC,2,N)*SL(2,N)*RHOL(2,N) +
     &      XGC(IRNV_GC,2,N)*SG(2,N)*RHOG(2,N))
        ELSEIF( IOM.EQ.43 ) THEN
          IF( IRNV_GC.EQ.1 ) THEN
            VAR = PORD(2,N)*VOL(N)*(XLA(2,N)*SL(2,N)*RHOL(2,N) +
     &        XGC(IRNV_GC,2,N)*SG(2,N)*RHOG(2,N) + 
     &        XNC(IRNV_GC,2,N)*SN(2,N)*RHON(2,N))
          ELSE
            VAR = PORD(2,N)*VOL(N)*(XGC(IRNV_GC,2,N)*SG(2,N)*RHOG(2,N) + 
     &        XNC(IRNV_GC,2,N)*SN(2,N)*RHON(2,N))
          ENDIF
        ELSEIF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
          VAR = PORD(2,N)*VOL(N)*(XLA(2,N)*SL(2,N)*RHOL(2,N) +
     &      XGA(2,N)*SG(2,N)*RHOG(2,N) + XHA(2,N)*SH(2,N)*RHOH(2,N) +
     &      XNA(2,N)*SN(2,N)*RHON(2,N))
        ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = PORD(2,N)*VOL(N)*(XLA(2,N)*SL(2,N)*RHOL(2,N) +
     &      XGA(2,N)*SG(2,N)*RHOG(2,N))
        ELSE
          VAR = PORD(2,N)*VOL(N)*(XLA(2,N)*SL(2,N)*RHOL(2,N) +
     &      XGA(2,N)*SG(2,N)*RHOG(2,N) + XNA(2,N)*SN(2,N)*RHON(2,N))
        ENDIF
!
!---  Total oil mass (kg),
!     or CH4 mass (kg)  ---
!
      ELSEIF( IRNV.EQ.39 ) THEN
        IUNKG = 1
        VAR = 0.D+0
!
!---    Total CH4 mass (kg)
!
        IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
          VAR = PORD(2,N)*VOL(N)*(XLO(2,N)*SL(2,N)*RHOL(2,N) +
     &      XGO(2,N)*SG(2,N)*RHOG(2,N) + XHO(2,N)*SH(2,N)*RHOH(2,N) +
     &      XNO(2,N)*SN(2,N)*RHON(2,N)) +
     &      (1.D+0-PORT(2,N))*VOL(N)*XSO(2,N)*RHOS(IZ(N))
!
!---    Total oil mass (kg)
!
        ELSE
          VAR = PORD(2,N)*VOL(N)*(XLO(2,N)*SL(2,N)*RHOL(2,N) +
     &      XGO(2,N)*SG(2,N)*RHOG(2,N) +
     &      XNO(2,N)*SN(2,N)*RHON(2,N)) +
     &      (1.D+0-PORT(2,N))*VOL(N)*XSO(2,N)*RHOS(IZ(N))
        ENDIF
!
!---  Water mass source integral ---
!
      ELSEIF( IRNV.EQ.40 ) THEN
        IUNKG = 1
        VAR = SRCIW(N)
!
!---  Air mass source integral 
!     or gas-component mass source integral
!     or component mass source integral  ---
!
      ELSEIF( IRNV.EQ.41 ) THEN
        IUNKG = 1
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = SRCIGC(IRNV_GC,N)
        ELSE
          VAR = SRCIA(N)
        ENDIF
!
!---  Oil mass source integral  ---
!
      ELSEIF( IRNV.EQ.42 ) THEN
        IUNKG = 1
        VAR = SRCIO(N)
!
!---  Energy source integral  ---
!
      ELSEIF( IRNV.EQ.43 ) THEN
        IUNM = 2
        IUNKG = 1
        IUNS = -2
        VAR = SRCIT(N)
!
!---  X-dir. equivalent thermal conductivity  ---
!
      ELSEIF( IRNV.EQ.44 ) THEN
        INDX = 1
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        IF( ISLC(5).EQ.1 ) THEN
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),SI(2,N),
     &      THKL(2,N),THKG(2,N),THKN(2,N),THKI(2,N),
     &      PORD(2,N),PORT(2,N),INDX )
        ELSE
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &      THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &      PORD(2,N),PORT(2,N),INDX )
        ENDIF
!
!---  Y-dir. equivalent thermal conductivity  ---
!
      ELSEIF( IRNV.EQ.45 ) THEN
        INDX = 2
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        IF( ISLC(5).EQ.1 ) THEN
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),SI(2,N),
     &      THKL(2,N),THKG(2,N),THKN(2,N),THKI(2,N),
     &      PORD(2,N),PORT(2,N),INDX )
        ELSE
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &      THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &      PORD(2,N),PORT(2,N),INDX )
        ENDIF
!
!---  Z-dir. equivalent thermal conductivity  ---
!
      ELSEIF( IRNV.EQ.46 ) THEN
        INDX = 3
        IUNM = 1
        IUNKG = 1
        IUNS = -3
        IUNK = -1
        IF( ISLC(5).EQ.1 ) THEN
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),SI(2,N),
     &      THKL(2,N),THKG(2,N),THKN(2,N),THKI(2,N),
     &      PORD(2,N),PORT(2,N),INDX )
        ELSE
          VAR = THKE( IZ(N),SL(2,N),SN(2,N),ZERO,
     &      THKL(2,N),THKG(2,N),THKN(2,N),ZERO,
     &      PORD(2,N),PORT(2,N),INDX )
        ENDIF
!
!---  Salt volumetric concentration  ---
!
      ELSEIF( IRNV.EQ.47 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = TMS(2,N)
!
!---  Salt aqueous concentration  ---
!
      ELSEIF( IRNV.EQ.48 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XLS(2,N)*RHOL(2,N)
!
!---  Aqueous-flux Courant number  ---
!
      ELSEIF( IRNV.EQ.49 ) THEN
        VAR = CRNTL(N)
!
!---  Total salt mass  ---
!
      ELSEIF( IRNV.EQ.50 ) THEN
        IUNKG = 1
        VAR = TMS(2,N)*VOL(N)
!
!---  X aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.51 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = UL(1,NPX)
!
!---  Y aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.52 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = VL(1,NPY)
!
!---  Z aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.53 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = WL(1,NPZ)
!
!---  X gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.54 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = UG(1,NPX)
!
!---  Y gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.55 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = VG(1,NPY)
!
!---  Z gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.56 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = WG(1,NPZ)
!
!---  X NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.57 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = UN(1,NPX)
!
!---  Y NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.58 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = VN(1,NPY)
!
!---  Z NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.59 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = WN(1,NPZ)
!
!---  X heat flux  ---
!
      ELSEIF( IRNV.EQ.60 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = UQ(1,NPX)
!
!---  Y heat flux  ---
!
      ELSEIF( IRNV.EQ.61 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = VQ(1,NPY)
!
!---  Z heat flux  ---
!
      ELSEIF( IRNV.EQ.62 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = WQ(1,NPZ)
!
!---  Matric potential head or Webb matching point head  ---
!
      ELSEIF( IRNV.EQ.63 ) THEN
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          IUNM = 1
          VAR = SCHR(9,IZ(N))
        ELSE
          IUNM = 1
          VAR = MIN(PL(2,N)-PG(2,N),0.D+0)/RHORL/GRAV
        ENDIF
!
!---  X salt flux  ---
!
      ELSEIF( IRNV.EQ.64 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = US(1,NPX)
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.64 .AND. IRNV.LE.96 ) THEN
!
!---  Y salt flux  ---
!
      IF( IRNV.EQ.65 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = VS(1,NPY)
!
!---  Z salt flux  ---
!
      ELSEIF( IRNV.EQ.66 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = WS(1,NPZ)
!
!---  X node-centered salt flux  ---
!
      ELSEIF( IRNV.EQ.67 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(US(1,NPX)+US(1,NQX))
!
!---  Y node-centered salt flux  ---
!
      ELSEIF( IRNV.EQ.68 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(VS(1,NPY)+VS(1,NQY))
!
!---  Z node-centered salt flux  ---
!
      ELSEIF( IRNV.EQ.69 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(WS(1,NPZ)+WS(1,NQZ))
!
!---  Gas water mole fraction
!     or volumetric component mole fraction  ---
!
      ELSEIF( IRNV.EQ.70 ) THEN
        VAR = XMGW(2,N)
!
!---  Gas air mole fraction
!     or gas component mole fraction  ---
!
      ELSEIF( IRNV.EQ.71 ) THEN
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = XMGC(IRNV_GC,2,N)
        ELSE
          VAR = XMGA(2,N)
        ENDIF
!
!---  Gas oil mole fraction
!     or NAPL component mole fraction  ---
!
      ELSEIF( IRNV.EQ.72 ) THEN
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = XMNC(IRNV_GC,2,N)
        ELSE
          VAR = XMGO(2,N)
        ENDIF
!
!---  Gas water concentration
!     or volumetric component concentration  ---
!
      ELSEIF( IRNV.EQ.73 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XGW(2,N)*RHOG(2,N)
!
!---  Gas air concentration
!     or gas component concentration  ---
!
      ELSEIF( IRNV.EQ.74 ) THEN
        IUNM = -3
        IUNKG = 1
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = XGC(IRNV_GC,2,N)*RHOG(2,N)
        ELSE
          VAR = XGA(2,N)*RHOG(2,N)
        ENDIF
!
!---  Gas oil concentration
!     or NAPL component concentration  ---
!
      ELSEIF( IRNV.EQ.75 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XGO(2,N)*RHOG(2,N)
!
!---  Aqueous water concentration  ---
!
      ELSEIF( IRNV.EQ.76 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XLW(2,N)*RHOL(2,N)
!
!---  Aqueous air or gas component concentration  ---
!
      ELSEIF( IRNV.EQ.77 ) THEN
        IUNM = -3
        IUNKG = 1
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = XLC(IRNV_GC,2,N)*RHOL(2,N)
        ELSE
          VAR = XLA(2,N)*RHOL(2,N)
        ENDIF
!
!---  Aqueous oil concentration or nonaqueous liquid
!     component concentration  ---
!
      ELSEIF( IRNV.EQ.78 ) THEN
        IUNM = -3
        IUNKG = 1
        IF( IOM.EQ.43 ) THEN
          VAR = XNC(IRNV_GC,2,N)*RHON(2,N)
        ELSE
          VAR = XLO(2,N)*RHOL(2,N)
        ENDIF
!
!---  Gas Courant number  ---
!
      ELSEIF( IRNV.EQ.79 ) THEN
        VAR = CRNTG(N)
!
!---  Ice pressure (absolute), Pa,
!     or system pressure (absolute), Pa
!     or stress-xx, Pa  ---
!
      ELSEIF( IRNV.EQ.80 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          ISLC50X = ABS(ISLC(50))
          IMODX = MOD(ISLC50X,10)
!
!---      Effective stress  ---
!
          VAR = SIG_GM(1,N)
!
!---      Total stress  ---
!
          IF( IRNV_GC.EQ.0 ) THEN
!
!---        Poroelasticity ---
!
            IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
              PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
              VAR = VAR + PROP_GM(3,IZ(N))*(PX-PCMP(N))
            ENDIF
!
!---        Thermoelasticity ---
!
            IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
              BLK3X = PROP_GM(1,IZ(N))/(1.0D+0-2.0D+0*PROP_GM(2,IZ(N)))
              VAR = VAR + PROP_GM(4,IZ(N))*BLK3X*(T(2,N)-TCMP(N))
            ENDIF
          ENDIF
!
!---    STOMP-EOR  ---
!
        ELSEIF( IOM.EQ.43 ) THEN
          VAR = PSO(2,N) + PATM
        ELSE
          VAR = PI(2,N) + PATM
        ENDIF
!
!---  Ice saturation  ---
!
      ELSEIF( IRNV.EQ.81 ) THEN
        IF( IOM.EQ.34 ) THEN
          VAR = PI(2,N)
        ELSE
          VAR = SI(2,N)
        ENDIF
!
!---  Ice density  ---
!
      ELSEIF( IRNV.EQ.82 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = RHOI(2,N)
!
!---  Aqueous matrix saturation
!     or NAPL matrix saturation
!     or Eclipse gas saturation  ---
!
      ELSEIF( IRNV.EQ.83 ) THEN
        IF( IOM.EQ.43 ) THEN
          VAR = SI(1,N)
        ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = SI(1,N)
        ELSE
          VAR = SDPM(N)
        ENDIF
!
!---  Aqueous fracture saturation
!     or NAPL fracture saturation
!     or Eclipse oil saturation  ---
!     or Eclipse pressure  ---
!
      ELSEIF( IRNV.EQ.84 ) THEN
        IF( IOM.EQ.43 ) THEN
          VAR = SI(2,N)
        ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = SI(2,N)
       ELSE
          VAR = SDPF(N)
        ENDIF 
!
!---  Gas matrix  ---
!
      ELSEIF( IRNV.EQ.85 ) THEN
        VAR = MAX( 1.D+0-SDPM(N),ZERO )
!
!---  Gas fracture  ---
!
      ELSEIF( IRNV.EQ.86 ) THEN
        VAR = MAX( 1.D+0-SDPF(N),ZERO )
!
!---  X node-centered aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.87 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(UL(1,NPX)+UL(1,NQX))
!
!---  Y node-centered aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.88 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(VL(1,NPY)+VL(1,NQY))
!
!---  Z node-centered aqueous volumetric flux  ---
!
      ELSEIF( IRNV.EQ.89 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(WL(1,NPZ)+WL(1,NQZ))
!
!---  X node-centered gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.90 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(UG(1,NPX)+UG(1,NQX))
!
!---  Y node-centered gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.91 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(VG(1,NPY)+VG(1,NQY))
!
!---  Z node-centered gas volumetric flux  ---
!
      ELSEIF( IRNV.EQ.92 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(WG(1,NPZ)+WG(1,NQZ))
!
!---  X node-centered NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.93 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(UN(1,NPX)+UN(1,NQX))
!
!---  Y node-centered NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.94 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(VN(1,NPY)+VN(1,NQY))
!
!---  Z node-centered NAPL volumetric flux  ---
!
      ELSEIF( IRNV.EQ.95 ) THEN
        IUNM = 1
        IUNS = -1
        VAR = 0.5D+0*(WN(1,NPZ)+WN(1,NQZ))
!
!---  X node-centered heat flux  ---
!
      ELSEIF( IRNV.EQ.96 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = 0.5D+0*(UQ(1,NPX)+UQ(1,NQX))
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.96 .AND. IRNV.LE.128 ) THEN
!
!---  Y node-centered heat flux  ---
!
      IF( IRNV.EQ.97 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = 0.5D+0*(VQ(1,NPY)+VQ(1,NQY))
!
!---  Z node-centered heat flux  ---
!
      ELSEIF( IRNV.EQ.98 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = 0.5D+0*(WQ(1,NPZ)+WQ(1,NQZ))
!
!---  NAPL Courant number  ---
!
      ELSEIF( IRNV.EQ.99 ) THEN
        VAR = CRNTN(N)
!
!---  Node number  ---
!
      ELSEIF( IRNV.EQ.100 ) THEN
        VAR = REAL(N)
!
!---  Osmotic pressure, Pa 
!     or stress-yy, Pa  ---
!
      ELSEIF( IRNV.EQ.101 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          ISLC50X = ABS(ISLC(50))
          IMODX = MOD(ISLC50X,10)
!
!---      Effective stress  ---
!
          VAR = SIG_GM(2,N)
!
!---      Total stress  ---
!
          IF( IRNV_GC.EQ.0 ) THEN
!
!---        Poroelasticity ---
!
            IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
              PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
              VAR = VAR + PROP_GM(3,IZ(N))*(PX-PCMP(N))
            ENDIF
!
!---        Thermoelasticity ---
!
            IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
              BLK3X = PROP_GM(1,IZ(N))/(1.0D+0-2.0D+0*PROP_GM(2,IZ(N)))
              VAR = VAR + PROP_GM(4,IZ(N))*BLK3X*(T(2,N)-TCMP(N))
            ENDIF
          ENDIF
        ELSE
          VAR = POSM(2,N)
        ENDIF
!
!---  Osmotic efficiency factor  ---
!
      ELSEIF( IRNV.EQ.102 ) THEN
        VAR = OEC(2,N)
!
!---  Aqueous compressibility, 1/Pa  ---
!
      ELSEIF( IRNV.EQ.103 ) THEN
        IUNM = 1
        IUNKG = -1
        IUNS = 2
        VAR = TORL(1,N)
!
!---  Gas compressibility, 1/Pa  ---
!
      ELSEIF( IRNV.EQ.104 ) THEN
        IUNM = 1
        IUNKG = -1
        IUNS = 2
        VAR = TORG(1,N)
!
!---  Trapped-gas saturation ---
!
      ELSEIF( IRNV.EQ.105 ) THEN
        VAR = SGT(2,N)
!
!---  Trapped-NAPL saturation ---
!
      ELSEIF( IRNV.EQ.106 ) THEN
        VAR = SNT(2,N)
!
!---  Aqueous trapped-gas saturation  ---
!
      ELSEIF( IRNV.EQ.107 ) THEN
        VAR = SGTL(N)
!
!---  NAPL trapped-gas saturation  ---
!
      ELSEIF( IRNV.EQ.108 ) THEN
        VAR = SGTN(N)
!
!---  Aqueous oil concentration  ---
!
      ELSEIF( IRNV.EQ.109 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XLO(2,N)*RHOL(2,N)
!
!---  Salt aqueous mass fraction  ---
!
      ELSEIF( IRNV.EQ.110 ) THEN
        VAR = XLS(2,N)
!
!---  Surfactant volumetric concentration  ---
!
      ELSEIF( IRNV.EQ.111 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = TMS(2,N)
!
!---  Surfactant aqueous concentration  ---
!
      ELSEIF( IRNV.EQ.112 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XLS(2,N)*RHOL(2,N)
!
!---  Surfactant aqueous mass fraction
!     or gas N2 mass fraction  ---
!
      ELSEIF( IRNV.EQ.113 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          VAR = XGN(2,N)
        ELSE
          VAR = XLS(2,N)
        ENDIF
!
!---  X surfactant flux  ---
!
      ELSEIF( IRNV.EQ.114 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = US(1,NPX)
!
!---  Y surfactant flux  ---
!
      ELSEIF( IRNV.EQ.115 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = VS(1,NPY)
!
!---  Z surfactant flux  ---
!
      ELSEIF( IRNV.EQ.116 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = WS(1,NPZ)
!
!---  X node-centered surfactant flux  ---
!
      ELSEIF( IRNV.EQ.117 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(US(1,NPX)+US(1,NQX))
!
!---  Y node-centered surfactant flux  ---
!
      ELSEIF( IRNV.EQ.118 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(VS(1,NPY)+VS(1,NQY))
!
!---  Z node-centered surfactant flux  ---
!
      ELSEIF( IRNV.EQ.119 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(WS(1,NPZ)+WS(1,NQZ))
!
!---  X dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.120 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = ULO(1,NPX)
!
!---  Y dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.121 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = VLO(1,NPY)
!
!---  Z dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.122 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = WLO(1,NPZ)
!
!---  X node-centered dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.123 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(ULO(1,NPX)+ULO(1,NQX))
!
!---  Y node-centered dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.124 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(VLO(1,NPY)+VLO(1,NQY))
!
!---  Z node-centered dissolved-oil flux  ---
!
      ELSEIF( IRNV.EQ.125 ) THEN
        IUNKG = 1
        IUNM = -2
        IUNS = -1
        VAR = 0.5D+0*(WLO(1,NPZ)+WLO(1,NQZ))
!
!---  Aqueous-NAPL trapping number  ---
!
      ELSEIF( IRNV.EQ.126 ) THEN
        VAR = TRPNL(2,N)
!
!---  Minimum effective aqueous saturation  ---
!
      ELSEIF( IRNV.EQ.127 ) THEN
        VAR = ASLMIN(2,N)
!
!---  Water-vapor partial pressure  ---
!
      ELSEIF( IRNV.EQ.128 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PVW(2,N)
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.128 .AND. IRNV.LE.160 ) THEN
!
!---  Air partial pressure or component partial pressure  ---
!
      IF( IRNV.EQ.129 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = PVC(IRNV_GC,2,N)
        ELSE
          VAR = PVA(2,N)
        ENDIF
!
!---  Oil-vapor partial pressure, Pa
!     or stress-zz, Pa  ---
!
      ELSEIF( IRNV.EQ.130 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          ISLC50X = ABS(ISLC(50))
          IMODX = MOD(ISLC50X,10)
!
!---      Effective stress  ---
!
          VAR = SIG_GM(3,N)
!
!---      Total stress  ---
!
          IF( IRNV_GC.EQ.0 ) THEN
!
!---        Poroelasticity ---
!
            IF( IMODX.NE.2 .AND. IMODX.NE.4 ) THEN
              PX = MAX( PG(2,N),PL(2,N),PN(2,N) ) + PATM
              VAR = VAR + PROP_GM(3,IZ(N))*(PX-PCMP(N))
            ENDIF
!
!---        Thermoelasticity ---
!
            IF( IMODX.NE.3 .AND. IMODX.NE.4 ) THEN
              BLK3X = PROP_GM(1,IZ(N))/(1.0D+0-2.0D+0*PROP_GM(2,IZ(N)))
              VAR = VAR + PROP_GM(4,IZ(N))*BLK3X*(T(2,N)-TCMP(N))
            ENDIF
          ENDIF
        ELSE
          VAR = PVO(2,N)
        ENDIF
!
!---  Gas-aqueous scaling factor  ---
!
      ELSEIF( IRNV.EQ.131 ) THEN
        VAR = BTGL(2,N)
!
!---  Free-NAPL aqueous interfacial area
!     or X CO2 displacement variance (m^2), 'SIGXX'  ---
!
      ELSEIF( IRNV.EQ.132 ) THEN
!
!---    STOMP-CO2 and STOMP-CO2E  ---
!
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          IUNM = 2
          IF( JSKPX.EQ.0 ) THEN
            VAR0X = 0.D+0
            VAR1X = 0.D+0
            DO L = 1,NFBN
              IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) CYCLE
              VARX = (SL(2,L)*XLA(2,L)*RHOL(2,L) +
     &          SG(2,L)*XGA(2,L)*RHOG(2,L))*VOL(N)*PORD(2,L)
              VAR0X = VAR0X + VARX
              VAR1X = VAR1X + VARX*XP(L)
              VAR2X = VAR2X + VARX*(XP(L)**2)
            ENDDO
            IF( VAR0X.GT.EPSL ) THEN
              VAR = VAR2X/VAR0X - (VAR1X/VAR0X)**2
            ELSE
              VAR = 0.D+0
            ENDIF
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          IUNM = -1
          VAR = ANLF(N)
        ENDIF
!
!---  Trapped-NAPL aqueous interfacial area  ---
!     or Y CO2 displacement variance (m^2), 'SIGYY'  ---
!
      ELSEIF( IRNV.EQ.133 ) THEN
!
!---    STOMP-CO2 and STOMP-CO2E  ---
!
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          IUNM = 2
          IF( JSKPX.EQ.0 ) THEN
            VAR0X = 0.D+0
            VAR1X = 0.D+0
            DO L = 1,NFBN
              IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) CYCLE
              VARX = (SL(2,L)*XLA(2,L)*RHOL(2,L) +
     &          SG(2,L)*XGA(2,L)*RHOG(2,L))*VOL(N)*PORD(2,L)
              VAR0X = VAR0X + VARX
              VAR1X = VAR1X + VARX*YP(L)*RP(ID(L))
              VAR2X = VAR2X + VARX*((YP(L)*RP(ID(L)))**2)
            ENDDO
            IF( VAR0X.GT.EPSL ) THEN
              VAR = VAR2X/VAR0X - (VAR1X/VAR0X)**2
            ELSE
              VAR = 0.D+0
            ENDIF
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          IUNM = -1
          VAR = ANLT(N)
        ENDIF
!
!---  Aqueous solute coefficient  ---
!     or Z CO2 displacement variance (m^2), 'SIGZZ'  ---
!
      ELSEIF( IRNV.EQ.134 ) THEN
!
!---    STOMP-CO2 and STOMP-CO2E  ---
!
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          IUNM = 2
          IF( JSKPX.EQ.0 ) THEN
            VAR0X = 0.D+0
            VAR1X = 0.D+0
            DO L = 1,NFBN
              IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) CYCLE
              VARX = (SL(2,L)*XLA(2,L)*RHOL(2,L) +
     &          SG(2,L)*XGA(2,L)*RHOG(2,L))*VOL(N)*PORD(2,L)
              VAR0X = VAR0X + VARX
              VAR1X = VAR1X + VARX*ZP(L)
              VAR2X = VAR2X + VARX*(ZP(L)**2)
            ENDDO
            IF( VAR0X.GT.EPSL ) THEN
              VAR = VAR2X/VAR0X - (VAR1X/VAR0X)**2
            ELSE
              VAR = 0.D+0
            ENDIF
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          IUNS = -1
          VAR = HKL(N)
        ENDIF
!
!---  Free-NAPL solute coefficient  ---
!
      ELSEIF( IRNV.EQ.135 ) THEN
        IUNS = -1
        VAR = HKNF(N)
!
!---  Trapped-NAPL solute coefficient  ---
!
      ELSEIF( IRNV.EQ.136 ) THEN
        IUNS = -1
        VAR = HKNT(N)
!
!---  Water relative humidity  ---
!
      ELSEIF( IRNV.EQ.137 ) THEN
        VAR = PVW(2,N)/PSW(2,N)
!
!---  Aqueous well pressure, or coupled-well pressure    ---
!
      ELSEIF( IRNV.EQ.138 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Coupled-well pressure    ---
!
        IF( L_CW.EQ.1 ) THEN
          VAR = P_CW(2,IRNV_CW) + PATM
!
!---     Aqueous well pressure    ---
!
        ELSE
          VAR = 0.D+0
          DO 62 NS = 1,NSR
            IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
              DO I = ISRDM(1,NS),ISRDM(2,NS)
              DO J = ISRDM(3,NS),ISRDM(4,NS)
              DO K = ISRDM(5,NS),ISRDM(6,NS)
                IF( ND(I,J,K).EQ.N ) THEN
                  VAR = PLWB(2,NS)
                  GOTO 64
                ENDIF
              ENDDO
              ENDDO
              ENDDO
            ENDIF
   62     CONTINUE
   64     CONTINUE
        ENDIF
!
!---  Water mass source rate or total water mass source rate (kg/s)  ---
!
      ELSEIF( IRNV.EQ.140 ) THEN
        IUNKG = 1
        IUNS = -1
        IF( IRNV_CW.EQ.-1 ) THEN
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
            DO L = 1,NFBN
              IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) CYCLE
              VAR = VAR + SRCW(2,L)
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          VAR = SRCW(2,N)
        ENDIF
!
!---  Air mass source rate
!     or gas component mass source rate  ---
!
      ELSEIF( IRNV.EQ.141 ) THEN
        IUNKG = 1
        IUNS = -1
        IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
          VAR = SRCGC(IRNV_GC,2,N)
        ELSE
          VAR = SRCA(2,N)
        ENDIF
!
!---  Oil mass source rate
!     or coupled-well volumetric nodal rate (m^3/s),
!     or coupled-well mass nodal rate (kg/s),
!     or EOR injection well mass rate (kg/s),
!
      ELSEIF( IRNV.EQ.142 ) THEN
!
!---    Coupled-well mass rate    ---
!
        IF( L_CW.EQ.1 ) THEN
!
!---      Coupled-well volumetric rate    ---
!
          IF( IOM.EQ.1 ) THEN
            IUNM = 3
            IUNS = -1
            VAR = 0.D+0
            DO NWN = ID_CW(3,IRNV_CW),ID_CW(4,IRNV_CW)
              IF( N.EQ.IWN_CW(NWN) ) THEN
                VAR = VAR + Q_CW(1,NWN)
              ENDIF
            ENDDO
!
!---      Coupled-well CO2 mass rate    ---
!
          ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
            IUNKG = 1
            IUNS = -1
            VAR = 0.D+0
            DO NWN = ID_CW(3,IRNV_CW),ID_CW(4,IRNV_CW)
              IF( N.EQ.IWN_CW(NWN) ) THEN
                VAR = VAR + FXA_CW(1,NWN)
              ENDIF
            ENDDO
!
!---      EOR injection well mass rate (kg/s)    ---
!
          ELSEIF( IOM.EQ.43 ) THEN
            IUNKG = 1
            IUNS = -1
            VAR = QM_CW(1,IRNV_CW)
          ENDIF
!
!---    Oil mass source rate  ---
!
        ELSE
          IUNM = -3
          IUNKG = 1
          IUNS = -1
          VAR = SRCO(2,N)
        ENDIF
!
!---  Energy source rate  ---
!
      ELSEIF( IRNV.EQ.143 ) THEN
        IUNM = 2
        IUNKG = 1
        IUNS = -3
        VAR = SRCT(2,N)
!
!---  Aqueous well depth  ---
!
      ELSEIF( IRNV.EQ.144 ) THEN
        VAR = 0.D+0
        DO 72 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = (PLWB(2,NS)-PGW(2,NS))/RHORL/GRAV
                GOTO 74
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
   72   CONTINUE
   74   CONTINUE
        IUNM = 1
!
!---  Aqueous well flow rate or coupled-well water mass nodal rate  ---
!
      ELSEIF( IRNV.EQ.145 ) THEN
!
!---    Coupled-well water mass rate    ---
!
        IF( L_CW.EQ.1 ) THEN
          IUNKG = 1
          IUNS = -1
          VAR = 0.D+0
          DO NWN = ID_CW(3,IRNV_CW),ID_CW(4,IRNV_CW)
            IF( N.EQ.IWN_CW(NWN) ) THEN
              VAR = VAR + FXW_CW(1,NWN)
            ENDIF
          ENDDO
!
!---    Aqueous well flow rate  ---
!
        ELSE
          VAR = 0.D+0
          DO NS = 1,NSR
            IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
              DO I = ISRDM(1,NS),ISRDM(2,NS)
              DO J = ISRDM(3,NS),ISRDM(4,NS)
              DO K = ISRDM(5,NS),ISRDM(6,NS)
                IF( ND(I,J,K).EQ.N ) THEN
                  VAR = QLW(3,NS)
                  EXIT
                ENDIF
              ENDDO
              ENDDO
              ENDDO
            ENDIF
          ENDDO
          IUNM = 3
          IUNS = -1
        ENDIF
!
!---  Aqueous well flow integral  ---
!
      ELSEIF( IRNV.EQ.146 ) THEN
        VAR = 0.D+0
        DO 92 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = QLW(1,NS)
                GOTO 94
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
   92   CONTINUE
   94   CONTINUE
        IUNM = 3
!
!---  Salt mass source rate  ---
!
      ELSEIF( IRNV.EQ.147 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = SRCS(2,N)
!
!---  Salt mass source integral  ---
!
      ELSEIF( IRNV.EQ.148 ) THEN
        IUNKG = 1
        VAR = SRCIS(N)
!
!---  Scanning path  ---
!
      ELSEIF( IRNV.EQ.149 ) THEN
        VAR = REAL( IPH(2,N) )
!
!---  Aqueous air or CO2 saturation
!     or Webb matching point saturation  ---
!
      ELSEIF( IRNV.EQ.150 ) THEN
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = SCHR(8,IZ(N))
        ELSE
          VAR = MIN( XMLA(2,N)*HCAW/(PVA(2,N)+SMALL),1.D+0 )
        ENDIF
!
!---  Bubble void fraction
!     or aqueous N2 mole fraction  ---
!
      ELSEIF( IRNV.EQ.151 ) THEN
        IF( IOM.EQ.39 ) THEN
          VAR = XMLN(2,N)
        ELSE
          VAR = SN(2,N)
        ENDIF
!
!---  Bubble air mass fraction
!     or aqueous N2 mass fraction  ---
!
      ELSEIF( IRNV.EQ.152 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          VAR = XLN(2,N)
        ELSE
          VAR = XGO(2,N)
        ENDIF
!
!---  Mineralized CO2   ---
!
      ELSEIF( IRNV.EQ.153 ) THEN
        IUNKG = 1
        IUNM = -3
        VAR = XLO(2,N)
!
!---  NAPL well flow rate   ---
!
      ELSEIF( IRNV.EQ.154 ) THEN
        VAR = 0.D+0
        IUNM = 3
        IUNS = -1
        DO 108 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = QNW(3,NS)
                GOTO 110
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
  108   CONTINUE
  110   CONTINUE
!
!---  NAPL well flow integral  ---
!
      ELSEIF( IRNV.EQ.155 ) THEN
        VAR = 0.D+0
        IUNM = 3
        DO 118 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = QNW(1,NS)
                GOTO 120
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
  118   CONTINUE
  120   CONTINUE
!
!---  Total well flow rate   ---
!
      ELSEIF( IRNV.EQ.156 ) THEN
        VAR = 0.D+0
        IUNM = 3
        IUNS = -1
        DO 128 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = QTW(3,NS)
                GOTO 130
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
  128   CONTINUE
  130   CONTINUE
!
!---  Total well flow integral  ---
!
      ELSEIF( IRNV.EQ.157 ) THEN
        VAR = 0.D+0
        IUNM = 3
        DO 138 NS = 1,NSR
          IF( ISRT(NS).GE.20 .AND. ISRT(NS).LE.29 ) THEN
            DO I = ISRDM(1,NS),ISRDM(2,NS)
            DO J = ISRDM(3,NS),ISRDM(4,NS)
            DO K = ISRDM(5,NS),ISRDM(6,NS)
              IF( ND(I,J,K).EQ.N ) THEN
                VAR = QTW(1,NS)
                GOTO 140
              ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF
  138   CONTINUE
  140   CONTINUE
!
!---  Gas N2 mole fraction  ---
!
      ELSEIF( IRNV.EQ.160 ) THEN
          VAR = XMGN(2,N)
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.160 .AND. IRNV.LE.200 ) THEN
!
!---  NAPL dissolved water concentration  ---
!
      IF( IRNV.EQ.161 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XNW(2,N)*RHON(2,N)
!
!---  NAPL dissolved water mole fraction
!     nonaqueous liquid mole fraction of N2 ---
!
      ELSEIF( IRNV.EQ.162 ) THEN
        IF( IOM.EQ.39 ) THEN
          VAR = XMNN(2,N)
        ELSE
          VAR = XMNW(2,N)
        ENDIF
!
!---  NAPL dissolved water mass fraction  ---
!
      ELSEIF( IRNV.EQ.163 ) THEN
        VAR = XNW(2,N)
!
!---  NAPL dissolved oil concentration  ---
!
      ELSEIF( IRNV.EQ.164 ) THEN
        IUNM = -3
        IUNKG = 1
        VAR = XNO(2,N)*RHON(2,N)
!
!---  NAPL dissolved oil mole fraction
!     or CH4 nonaqueous liquid mole fraction  ---
!
      ELSEIF( IRNV.EQ.165 ) THEN
        VAR = XMNO(2,N)
!
!---  NAPL dissolved oil mass fraction  ---
!
      ELSEIF( IRNV.EQ.166 ) THEN
        VAR = XNO(2,N)
!
!---  Total alcohol mass  ---
!
      ELSEIF( IRNV.EQ.167 ) THEN
        IUNKG = 1
        VAR = 0.D+0
        IF( IOM.EQ.7 ) VAR = PORD(2,N)*VOL(N)*
     &    (XLA(2,N)*SL(2,N)*RHOL(2,N)+XNA(2,N)*SN(2,N)*RHON(2,N))
!
!---  Aqueous dissolved water mass fraction  ---
!
      ELSEIF( IRNV.EQ.168 ) THEN
        VAR = XMLW(2,N)
!
!---  Alcohol mass source integral  ---
!
      ELSEIF( IRNV.EQ.169 ) THEN
        IUNKG = 1
        VAR = SRCIA(N)
!
!---  Alcohol mass source rate  ---
!
      ELSEIF( IRNV.EQ.170 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = SRCA(2,N)
!
!---  Total NAPL and aqueous dissolved alcohol
!     or integrated kerogen mass
!     or integrated hydrate N2 mass  ---
!
      ELSEIF( IRNV.EQ.171 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 171 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 171
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*XHN(2,L)*SH(2,L)*RHOH(2,L)
            ELSE
              VAR = VAR + PORD(2,L)*VOL(L)*(XLA(2,L)*SL(2,L)*RHOL(2,L)
     &          + XNA(2,L)*SN(2,L)*RHON(2,L))
            ENDIF
  171     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Total aqueous dissolved alcohol
!     or integrated char mass  ---
!
      ELSEIF( IRNV.EQ.172 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 172 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 172
            VAR = VAR + PORD(2,L)*VOL(L)*XLA(2,L)*SL(2,L)*RHOL(2,L)
  172     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Total NAPL dissolved water
!     or integrated coke mass  ---
!
      ELSEIF( IRNV.EQ.173 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 173 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 173
            VAR = VAR + PORD(2,L)*VOL(L)*XNW(2,L)*SN(2,L)*RHON(2,L)
  173     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Total NAPL dissolved oil  ---
!
      ELSEIF( IRNV.EQ.174 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 174 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 174
            VAR = VAR + PORD(2,L)*VOL(L)*XNO(2,L)*SN(2,L)*RHON(2,L)
  174     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Total NAPL dissolved alcohol  ---
!
      ELSEIF( IRNV.EQ.175 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 175 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 175
            VAR = VAR + PORD(2,L)*VOL(L)*XNA(2,L)*SN(2,L)*RHON(2,L)
  175     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Aqueous viscosity  ---
!
      ELSEIF( IRNV.EQ.176 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -1
        VAR = VISL(2,N)
!
!---  Monitoring well water depth
!     or total-liquid well depth  ---
!
      ELSEIF( IRNV.EQ.177 ) THEN
        IUNM = 1
        VAR = PI(2,N)
!
!---  Matric potential  ---
!
      ELSEIF( IRNV.EQ.178 ) THEN
        IUNM = 1
        IF( IOM.EQ.2 ) THEN
          VAR = (PL(2,N)-PG(2,N))/RHORL/GRAV
        ELSEIF( IOM.EQ.4 ) THEN
          VAR = PI(IEQW+2,N)
        ELSE
          VAR = (PL(2,N)-PG(2,N))/RHORL/GRAV
        ENDIF
!
!---  NAPL well depth  ---
!
      ELSEIF( IRNV.EQ.179 ) THEN
        IUNM = 1
        IF( IOM.EQ.4 ) THEN
          VAR = PI(IEQO+2,N)
        ENDIF
!
!---  Monitoring well pressure or coupled-well nodal pressure  ---
!
      ELSEIF( IRNV.EQ.180 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Coupled-well nodal pressure    ---
!
        IF( L_CW.EQ.1 ) THEN
          VAR = 0.D+0
          DO 180 NWF = ID_CW(5,IRNV_CW),ID_CW(6,IRNV_CW)
            IF( N.EQ.IWF_CW(NWF) ) VAR = PF_CW(NWF) + PATM
  180     CONTINUE
!
!---    Monitoring well pressure    ---
!
        ELSE
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSEIF( IOM.EQ.2 ) THEN
            VAR = MAX(PWLW(2,NW),PWGW(2,NW)) + PATM
          ELSEIF( IOM.EQ.4 ) THEN
            VAR = PW(2,NW) + PATM
          ENDIF
        ENDIF
!
!---  Monitoring well aqueous saturation  ---
!
      ELSEIF( IRNV.EQ.181 ) THEN
        NW = ABS(IXW(N))
        IF( NW.EQ.0 ) THEN
          VAR = 0.D+0
        ELSE
          VAR = SLW(2,NW)
        ENDIF
!
!---  Monitoring well water-vapor mass fraction or
!     well NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.182 ) THEN
        NW = ABS(IXW(N))
        IF( NW.EQ.0 ) THEN
          VAR = 0.D+0
        ELSEIF( IOM.EQ.2 ) THEN
          VAR = XGWW(2,NW)
        ELSEIF( IOM.EQ.4 ) THEN
          VAR = SNW(2,NW)
        ENDIF
!
!---  Monitoring well dissolved-air mass fraction
!     or dissolved-oil mass fraction  ---
!
      ELSEIF( IRNV.EQ.183 ) THEN
        NW = ABS(IXW(N))
        IF( NW.EQ.0 ) THEN
          VAR = 0.D+0
        ELSEIF( IOM.EQ.2 ) THEN
          VAR = XLAW(2,NW)
        ELSEIF( IOM.EQ.4 ) THEN
          VAR = XLOW(2,NW)
        ENDIF
!
!---  Monitoring well axial aqueous flux
!     or well total-liquid pumping rate  ---
!
      ELSEIF( IRNV.EQ.184 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            NWL = IWL(NW)
            VAR = QT_W(3,NWL)
          ENDIF
        ELSE
          IUNM = 1
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            VAR = UL_W(1,NW)
          ENDIF
        ENDIF
!
!---  Monitoring well axial gas flux or
!     or well aqueous pumping rate  ---
!
      ELSEIF( IRNV.EQ.185 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            NWL = IWL(NW)
            VAR = QL_W(3,NWL)
          ENDIF
        ELSE
          IUNM = 1
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            VAR = UG_W(1,NW)
          ENDIF
        ENDIF
!
!---  Monitoring well vertical aqueous flux
!     or well NAPL pumping rate  ---
!
      ELSEIF( IRNV.EQ.186 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            NWL = IWL(NW)
            VAR = QN_W(3,NWL)
          ENDIF
        ELSE
          IUNM = 1
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            VAR = WL_W(1,NW)
          ENDIF
        ENDIF
!
!---  Monitoring well vertical gas flux
!     or well total-liquid pumping integral  ---
!
      ELSEIF( IRNV.EQ.187 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            NWL = IWL(NW)
            VAR = QT_W(4,NWL)
          ENDIF
        ELSE
          IUNM = 1
          IUNS = -1
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            VAR = WG_W(1,NW)
          ENDIF
        ENDIF
!
!---  Well aqueous pumping integral  ---
!
      ELSEIF( IRNV.EQ.188 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          NW = ABS(IXW(N))
          IF( NW.EQ.0 ) THEN
            VAR = 0.D+0
          ELSE
            NWL = IWL(NW)
            VAR = QL_W(4,NWL)
          ENDIF
        ENDIF
!
!---  Integrated mineral CO2
!     or well NAPL pumping integral  ---
!
      ELSEIF( IRNV.EQ.189 ) THEN
        IF( IOM.EQ.4 ) THEN
          IUNM = 3
          NW = ABS(IXW(N))
          NWL = IWL(NW)
          VAR = QN_W(4,NWL)
        ELSE
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
            DO 189 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 189
              VAR = VAR + VOL(L)*XLO(2,L)
  189       CONTINUE
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ENDIF
!
!---  Integrated trapped gas air mass, or
!     integrated trapped gas CO2 mass, or
!     integrated N2 mass  ---
!
      ELSEIF( IRNV.EQ.190 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 190 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 190
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLN(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGN(2,L)*SG(2,L)*RHOG(2,L) +
     &           XHN(2,L)*SH(2,L)*RHOH(2,L) +
     &           XNN(2,L)*SN(2,L)*RHON(2,L))
            ELSE
              VAR = VAR + PORD(2,L)*VOL(L)*XGA(2,L)*SGT(2,L)*RHOG(2,L)
            ENDIF
  190     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated water mass  ---
!
      ELSEIF( IRNV.EQ.191 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 191 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 191
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLW(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGW(2,L)*SG(2,L)*RHOG(2,L) +
     &           XNW(2,L)*SN(2,L)*RHON(2,L) +
     &           XHW(2,L)*SH(2,L)*RHOH(2,L) +
     &           SI(2,L)*RHOI(2,L))
!
!---        STOMP-HYD through STOMP-HYD-KE  ---
!
            ELSEIF( IOM.GE.36 .AND. IOM.LE.38 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLW(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGW(2,L)*SG(2,L)*RHOG(2,L) +
     &           XHW(2,L)*SH(2,L)*RHOH(2,L) +
     &           SI(2,L)*RHOI(2,L))
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLW(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGW(2,L)*SG(2,L)*RHOG(2,L)
              IF( ISLC(5).EQ.1 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*SI(2,L)*RHOI(2,L)
              IF( IOM.EQ.7 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XNW(2,L)*SN(2,L)*RHON(2,L)
            ENDIF
  191     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated air, CO2, or gas component mass  ---
!
      ELSEIF( IRNV.EQ.192 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 192 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 192
!
!---        STOMP-HYD through STOMP-HYDT-KE  ---
!
            IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLA(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGA(2,L)*SG(2,L)*RHOG(2,L) +
     &           XHA(2,L)*SH(2,L)*RHOH(2,L) +
     &           XNA(2,L)*SN(2,L)*RHON(2,L))
            ELSEIF( IOM.EQ.30 .OR. IOM.EQ.40 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLC(IRNV_GC,2,L)*SL(2,L)*RHOL(2,L) +
     &           XGC(IRNV_GC,2,L)*SG(2,L)*RHOG(2,L))
!
!---        STOMP-EOR  ---
!
            ELSEIF( IOM.EQ.43 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XLC(IRNV_GC,2,L)*SL(2,L)*RHOL(2,L) +
     &           XGC(IRNV_GC,2,L)*SG(2,L)*RHOG(2,L) +
     &           XNC(IRNV_GC,2,L)*SN(2,L)*RHON(2,L) )
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLA(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGA(2,L)*SG(2,L)*RHOG(2,L)
              IF( ISLC(13).EQ.1 ) VAR = VAR + XGO(2,L)
            ENDIF
  192     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated oil or CH4 mass  ---
!
      ELSEIF( IRNV.EQ.193 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 193 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 193
!
!---        STOMP-HYD through STOMP-HYDT-KE  ---
!
            IF( IOM.GE.36 .AND. IOM.LE.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLO(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGO(2,L)*SG(2,L)*RHOG(2,L) +
     &           XNO(2,L)*SN(2,L)*RHON(2,L) +
     &           XHO(2,L)*SH(2,L)*RHOH(2,L))
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLO(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGO(2,L)*SG(2,L)*RHOG(2,L)
              IF( IEQO.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*SN(2,L)*RHON(2,L) +
     &          (1.D+0-PORT(2,L))*VOL(L)*XSO(2,L)*RHOS(IZ(L))
              IF( IOM.EQ.7 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XNO(2,L)*SN(2,L)*RHON(2,L)
            ENDIF
  193     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated aqueous water mass  ---
!
      ELSEIF( IRNV.EQ.194 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 194 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 194
            VAR = VAR +
     &        PORD(2,L)*VOL(L)*XLW(2,L)*SL(2,L)*RHOL(2,L)
            IF( ISLC(5).EQ.1 ) VAR = VAR +
     &        PORD(2,L)*VOL(L)*SI(2,L)*RHOI(2,L)
  194     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated aqueous air, integrated aqueous CO2, or 
!     integrated gas component mass  ---
!
      ELSEIF( IRNV.EQ.195 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 195 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 195
            IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLC(IRNV_GC,2,L)*SL(2,L)*RHOL(2,L)
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLA(2,L)*SL(2,L)*RHOL(2,L)
            ENDIF
  195     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated aqueous oil mass  ---
!
      ELSEIF( IRNV.EQ.196 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 196 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 196
            VAR = VAR + PORD(2,L)*VOL(L)*XLO(2,L)*SL(2,L)*RHOL(2,L)
  196     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated gas water mass  ---
!
      ELSEIF( IRNV.EQ.197 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 197 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 197
            VAR = VAR + PORD(2,L)*VOL(L)*XGW(2,L)*SG(2,L)*RHOG(2,L)
  197     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated gas air, integrated gas CO2 or integrated gas 
!     component mass  ---
!
      ELSEIF( IRNV.EQ.198 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 198 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 198
            IF( IOM.EQ.30 .OR. IOM.EQ.40 .OR. IOM.EQ.43 ) THEN
              VAR = VAR +
     &        PORD(2,L)*VOL(L)*XGC(IRNV_GC,2,L)*SG(2,L)*RHOG(2,L)
            ELSE
              VAR = VAR + PORD(2,L)*VOL(L)*XGA(2,L)*SG(2,L)*RHOG(2,L)
              IF( ISLC(13).EQ.1 ) VAR = VAR + XGO(2,L)
            ENDIF
  198     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated gas oil mass  ---
!
      ELSEIF( IRNV.EQ.199 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 199 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 199
            VAR = VAR + PORD(2,L)*VOL(L)*XGO(2,L)*SG(2,L)*RHOG(2,L)
  199     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.200 .AND. IRNV.LE.240 ) THEN
!
!---  X aqueous relative permeability ---
!
      IF( IRNV.EQ.201 ) THEN
        VAR = RKL(1,2,N)
!
!---  Y aqueous relative permeability ---
!
      ELSEIF( IRNV.EQ.202 ) THEN
        VAR = RKL(2,2,N)
!
!---  Z aqueous relative permeability ---
!
      ELSEIF( IRNV.EQ.203 ) THEN
        VAR = RKL(3,2,N)
!
!---  Aqueous CO2 mole fraction  ---
!     or aqueous component mole fraction  ---
!
      ELSEIF( IRNV.EQ.204 ) THEN
        IF( IOM.EQ.30 .OR. IOM.EQ.40 ) THEN
          VAR = XMLC(IRNV_GC,2,N)
        ELSE
        VAR = XMLA(2,N)
        ENDIF
!
!---  Aqueous salt mole fraction  ---
!
      ELSEIF( IRNV.EQ.205 ) THEN
        VAR = XMLS(2,N)
!
!---  Atmospheric temperature  ---
!
      ELSEIF( IRNV.EQ.206 ) THEN
        IUNK = 1
        IF( IRNV_CW.GT.0 ) THEN
          VAR = P_CW(2,IRNV_CW)
        ELSE
          VAR = PN(1,N)
        ENDIF
!
!---  Atmospheric relative humidity  ---
!
      ELSEIF( IRNV.EQ.207 ) THEN
        VAR = PN(3,N)
!
!---  Atmospheric solar radiation  ---
!
      ELSEIF( IRNV.EQ.208 ) THEN
        IUNKG = 1
        IUNS = -3
        INDX = 5
        VAR = PN(INDX,N)
!
!---  Atmospheric wind speed  ---
!
      ELSEIF( IRNV.EQ.209 ) THEN
        IUNM = 1
        IUNS = -1
        INDX = 4
        VAR = PN(INDX,N)
!
!---  Residual NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.210 ) THEN
        VAR = SNR(2,N)
!
!---  Mobile NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.211 ) THEN
        VAR = SN(2,N)-SNR(2,N)-SNT(2,N)
!
!---  Free NAPL saturation  ---
!
      ELSEIF( IRNV.EQ.212 ) THEN
        VAR = SN(2,N)-SNT(2,N)
!
!---  Surface temperature  ---
!
      ELSEIF( IRNV.EQ.213 ) THEN
        IUNK = 1
        VAR = XLO(1,N)
!
!---  Surface vapor pressure, Pa
!     or stress-yz, Pa  ---
!
      ELSEIF( IRNV.EQ.214 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          VAR = SIG_GM(4,N)
        ELSE
          INDX = 5
          VAR = XLO(INDX,N)
        ENDIF
!
!---  Actual evaporation rate  ---
!
      ELSEIF( IRNV.EQ.215 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = RHON(2,N)
!
!---  Potential evaporation rate  ---
!
      ELSEIF( IRNV.EQ.216 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = RHON(3,N)
!
!---  Actual transpiration rate  ---
!
      ELSEIF( IRNV.EQ.217 ) THEN
        IUNKG = 1
        IUNS = -1
        INDX = 4
        VAR = RHON(INDX,N)
!
!---  Potential transpiration rate  ---
!
      ELSEIF( IRNV.EQ.218 ) THEN
        IUNKG = 1
        IUNS = -1
        INDX = 5
        VAR = RHON(INDX,N)
!
!---  Saturated CO2 aqueous mass fraction  ---
!
      ELSEIF( IRNV.EQ.219 ) THEN
        VAR = PVO(2,N)
!
!---  Aqueous alcohol mole fraction
!     or aqueous N2 mole fraction  ---
!
      ELSEIF( IRNV.EQ.220 ) THEN
        IF( IOM.EQ.39 ) THEN
          VAR = XMLN(2,N)
        ELSE
          VAR = XMLA(2,N)
        ENDIF
!
!---  NAPL alcohol mole fraction
!     or CO2 nonaqueous liquid mole fraction  ---
!
      ELSEIF( IRNV.EQ.221 ) THEN
        VAR = XMNA(2,N)
!
!---  Aqueous alcohol mass fraction
!     or N2 nonaqueous liquid mass fraction  ---
!
      ELSEIF( IRNV.EQ.222 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          VAR = XNN(2,N)
        ELSE
          VAR = XLA(2,N)
        ENDIF
!
!---  NAPL alcohol mass fraction  ---
!
      ELSEIF( IRNV.EQ.223 ) THEN
        VAR = XNA(2,N)
!
!---  Atmospheric pressure, Pa
!     or stress-xz, Pa  ---
!
      ELSEIF( IRNV.EQ.224 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          VAR = SIG_GM(5,N)
        ELSE
          VAR = PN(2,N)
        ENDIF
!
!---  Surface aqueous pressure  ---
!
      ELSEIF( IRNV.EQ.225 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = XLO(2,N)
!
!---  Surface gas pressure, Pa
!     or stress-xy, Pa  ---
!
      ELSEIF( IRNV.EQ.226 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
!
!---    Geomechanics  ---
!
        IF( ISLC(50).NE.0 ) THEN
          VAR = SIG_GM(5,N)
        ELSE
          VAR = XLO(3,N)
        ENDIF
!
!---  Surface aqueous saturation  ---
!
      ELSEIF( IRNV.EQ.227 ) THEN
        INDX = 4
        VAR = XLO(INDX,N)
!
!---  Surface latent heat flux  ---
!
      ELSEIF( IRNV.EQ.228 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = SN(1,N)
!
!---  Surface sensible heat flux  ---
!
      ELSEIF( IRNV.EQ.229 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = SN(2,N)
!
!---  Surface net long-wave radiation flux  ---
!
      ELSEIF( IRNV.EQ.230 ) THEN
        IUNKG = 1
        IUNS = -3
        VAR = SN(3,N)
!
!---  Surface net short-wave radiation flux  ---
!
      ELSEIF( IRNV.EQ.231 ) THEN
        IUNKG = 1
        IUNS = -3
        INDX = 4
        VAR = SN(INDX,N)
!
!---  Surface ground heat flux  ---
!
      ELSEIF( IRNV.EQ.232 ) THEN
        IUNKG = 1
        IUNS = -3
        INDX = 5
        VAR = SN(INDX,N)
!
!---  Surface water mass flux  ---
!
      ELSEIF( IRNV.EQ.233 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = RHON(1,N)
!
!---  Plant temperature  ---
!
      ELSEIF( IRNV.EQ.234 ) THEN
        IUNK = 1
        VAR = XMLO(1,N)
!
!---  Plant temperature  ---
!
      ELSEIF( IRNV.EQ.235 ) THEN
        IUNK = 1
        VAR = XMLO(2,N)
!
!---  Plant temperature  ---
!
      ELSEIF( IRNV.EQ.236 ) THEN
        IUNK = 1
        VAR = XMLO(3,N)
!
!---  Plant temperature  ---
!
      ELSEIF( IRNV.EQ.237 ) THEN
        IUNK = 1
        INDX = 4
        VAR = XMLO(INDX,N)
!
!---  Plant temperature  ---
!
      ELSEIF( IRNV.EQ.238 ) THEN
        IUNK = 1
        INDX = 5
        VAR = XMLO(INDX,N)
!
!---  Rainfall interception mass or integrated aqueous N2, mass, kg  ---
!
      ELSEIF( IRNV.EQ.239 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
            DO 239 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 239
              VAR = VAR + PORD(2,L)*VOL(L)*XLN(2,L)*SL(2,L)*RHOL(2,L)
  239       CONTINUE
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          IUNKG = 1
          VAR = XGO(1,N)
        ENDIF
!
!---  Sorbed oil mass, kg or integrated gas N2 mass, kg  ---
!
      ELSEIF( IRNV.EQ.240 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
            DO 240 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 240
              VAR = VAR + PORD(2,L)*VOL(L)*XGN(2,L)*SG(2,L)*RHOG(2,L)
  240       CONTINUE
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ELSE
          IUNKG = 1
          VAR = VOL(N)*(1.D+0-PORT(2,N))*XSO(2,N)*RHOS(IZ(N))
        ENDIF
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.240 .AND. IRNV.LE.280 ) THEN
!
!---  Sorbed Oil Mass Fraction
!     or N2 Hydrate Mass Fraction  ---
!
      IF( IRNV.EQ.241 ) THEN
!
!---    STOMP-HYDT-KE  ---
!
        IF( IOM.EQ.39 ) THEN
          VAR = XHN(2,N)
        ELSE
          VAR = XSO(2,N)
        ENDIF
!
!---  Sorbed Oil Concentration ---
!
      ELSEIF( IRNV.EQ.242 ) THEN
        IUNKG = 1
        IUNM = -3
        VAR = (1.D+0-PORT(2,N))*XSO(2,N)*RHOS(IZ(N))
!
!---  Bare-Soil Aerodynamic Resistance ---
!
      ELSEIF( IRNV.EQ.243 ) THEN
        IUNS = 1
        IUNM = -1
        VAR = XGO(2,N)
!
!---  Surface volumetric precipitation rate  ---
!
      ELSEIF( IRNV.EQ.244 ) THEN
        IUNS = -1
        IUNM = 3
        IX = 3
        VAR = XGO(IX,N)
!
!---  Surface Mass Precipitation Rate (kg/s)  ---
!
      ELSEIF( IRNV.EQ.245 ) THEN
        IUNS = -1
        IUNKG = 1
!
!---    Fracture flow, use null value  ---
!
        IF( ISLC(74).NE.0 ) THEN
          VAR = -9999.D+22
        ELSE
          IX = 4
          VAR = XGO(IX,N)
        ENDIF
!
!---  Atmospheric precipitation rate (kg/s),
!     or injection well petroleum component mass rate (kg/s),
!     or production well petroleum component mass rate (kg/s)  ---
!
      ELSEIF( IRNV.EQ.246 ) THEN
        IUNS = -1
        IUNKG = 1
!
!---    Fracture flow, use null value  ---
!
        IF( ISLC(74).NE.0 ) THEN
          VAR = -9999.D+22
!
!---    Injection well petroleum component mass rate (kg/s),
!       or production well petroleum component mass rate (kg/s)  ---
!
        ELSEIF( L_CW.EQ.1 ) THEN
          IVAR = 7 + (IRNV_GC-1)*2
          VAR = QM_CW(IVAR,IRNV_CW)
!
!---    Atmospheric precipitation rate (kg/s)  ---
!
        ELSE
          IX = 5
          VAR = XGO(IX,N)
        ENDIF
!
!---  X-Direction Intrinsic Permeability
!     or Matrix Permeability  ---
!
      ELSEIF( IRNV.EQ.247 ) THEN
        IUNM = 2
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = PERMV(1,N)*PERMRF(2,N)
        ELSE
          VAR = PERM(1,IZ(N))*PERMRF(2,N)
        ENDIF
!
!---  Y-Direction Intrinsic Permeability
!     or Fracture Permeability  ---
!
      ELSEIF( IRNV.EQ.248 ) THEN
        IUNM = 2
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = PERMV(2,N)*PERMRF(2,N)
        ELSE
          VAR = PERM(2,IZ(N))*PERMRF(2,N)
        ENDIF
!
!---  Z-Direction Intrinsic Permeability  ---
!
      ELSEIF( IRNV.EQ.249 ) THEN
        IUNM = 2
        IF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
          VAR = PERMV(3,N)*PERMRF(2,N)
        ELSE
          VAR = PERM(3,IZ(N))*PERMRF(2,N)
        ENDIF
!
!---  Hydrate water mass fraction  ---
!
      ELSEIF( IRNV.EQ.250 ) THEN
        VAR = XHW(2,N)
!
!---  Hydrate CO2 mass fraction  ---
!
      ELSEIF( IRNV.EQ.251 ) THEN
        VAR = XHA(2,N)
!
!---  Hydrate CH4 mass fraction  ---
!
      ELSEIF( IRNV.EQ.252 ) THEN
        VAR = XHO(2,N)
!
!---  Hydrate density or production well fluid density ---
!
      ELSEIF( IRNV.EQ.253 ) THEN
        IUNM = -3
        IUNKG = 1
!
!---    Production well fluid density  ---
!
        IF( L_CW.EQ.1 ) THEN
          VAR = RHOF_CW(IRNV_CW)
!
!---    Hydrate density  ---
!
        ELSE
          VAR = RHOH(2,N)
        ENDIF
!
!---  Hydrate saturation ---
!
      ELSEIF( IRNV.EQ.254 ) THEN
        VAR = SH(2,N)
!
!---  Hydrate pressure ---
!
      ELSEIF( IRNV.EQ.255 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PH(2,N)
!
!---  Integrated hydrate water mass  ---
!
      ELSEIF( IRNV.EQ.256 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 256 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 256
            VAR = VAR + PORD(2,L)*VOL(L)*XHW(2,L)*SH(2,L)*RHOH(2,L)
  256     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated hydrate CO2 mass (kg),
!     or injection well petroleum component mass integral (kg),
!     or production well petroleum component mass integral (kg)  ---
!
      ELSEIF( IRNV.EQ.257 ) THEN
        IUNKG = 1
!
!---    Injection well petroleum component mass integral (kg)
!       or Production well petroleum component mass integral (kg)---
!
        IF( L_CW.EQ.1 ) THEN
          IVAR = 8 + (IRNV_GC-1)*2
          VAR = QM_CW(IVAR,IRNV_CW)
!
!---    Integrated hydrate CO2 mass (kg)  ---
!
        ELSE
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
            DO 257 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 257
              VAR = VAR + PORD(2,L)*VOL(L)*XHA(2,L)*SH(2,L)*RHOH(2,L)
  257       CONTINUE
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ENDIF
!
!---  Integrated hydrate CH4 mass  ---
!
      ELSEIF( IRNV.EQ.258 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 258 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 258
            VAR = VAR + PORD(2,L)*VOL(L)*XHO(2,L)*SH(2,L)*RHOH(2,L)
  258     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated aqueous CH4 mass  ---
!
      ELSEIF( IRNV.EQ.259 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 259 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 259
            VAR = VAR + PORD(2,L)*VOL(L)*XLO(2,L)*SL(2,L)*RHOL(2,L)
  259     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated gas CH4 mass  ---
!
      ELSEIF( IRNV.EQ.260 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 260 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 260
            VAR = VAR + PORD(2,L)*VOL(L)*XGO(2,L)*SG(2,L)*RHOG(2,L)
  260     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated water mass source  ---
!
      ELSEIF( IRNV.EQ.261 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 261 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 261
            VAR = VAR + SRCIW(L)
  261     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated air or CO2 mass source  ---
!
      ELSEIF( IRNV.EQ.262 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 262 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 262
            VAR = VAR + SRCIA(L)
  262     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Integrated oil or CH4 mass source  ---
!
      ELSEIF( IRNV.EQ.263 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 263 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 263
            VAR = VAR + SRCIO(L)
  263     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Precipitated salt saturation  ---
!
      ELSEIF( IRNV.EQ.264 ) THEN
        VAR = SS(2,N)
!
!---  Hydrate water mole fraction  ---
!
      ELSEIF( IRNV.EQ.265 ) THEN
        VAR = (XHW(2,N)/WTMW)/
     &  ((XHA(2,N)/WTMA)+(XHO(2,N)/WTMO)+(XHW(2,N)/WTMW))
!
!---  Hydrate CO2 mole fraction  ---
!
      ELSEIF( IRNV.EQ.266 ) THEN
        VAR = (XHA(2,N)/WTMA)/
     &  ((XHA(2,N)/WTMA)+(XHO(2,N)/WTMO)+(XHW(2,N)/WTMW))
!
!---  Hydrate CH4 mole fraction  ---
!
      ELSEIF( IRNV.EQ.267 ) THEN
        VAR = (XHO(2,N)/WTMO)/
     &  ((XHA(2,N)/WTMA)+(XHO(2,N)/WTMO)+(XHW(2,N)/WTMW))
!
!---  Rain-Water Runoff Mass Rate  ---
!
      ELSEIF( IRNV.EQ.274 ) THEN
        IUNS = -1
        IUNM = 3
        INDX = 1
        VAR = XGO(INDX,N)
!
!---  Well-Node Pressure or MCStress  ---
!
      ELSEIF( IRNV.EQ.275 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PVO(1,N)
!
!---  Well-Node Temperature  ---
!
      ELSEIF( IRNV.EQ.276 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PSO(1,N)
!
!---  Differential integrated hydrate water mass  ---
!
      ELSEIF( IRNV.EQ.277 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 277 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 277
            VAR = VAR + PORD(2,L)*VOL(L)*XHW(2,L)*SH(2,L)*RHOH(2,L)
  277     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMHWX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMHWX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated hydrate CO2 mass  ---
!
      ELSEIF( IRNV.EQ.278 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 278 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 278
 !           VAR = VAR + PORD(2,L)*VOL(L)*XHA(2,L)*SH(2,L)*RHOH(2,L)
            VAR = VAR + PORD(2,L)*VOL(L)*TMHA(2,L)
  278     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMHAX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMHAX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated hydrate CH4 mass  ---
!
      ELSEIF( IRNV.EQ.279 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 279 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 279
!            VAR = VAR + PORD(2,L)*VOL(L)*XHO(2,L)*SH(2,L)*RHOH(2,L)
            VAR = VAR + PORD(2,L)*VOL(L)*TMHO(2,L)
  279     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMHOX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMHOX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated hydrate N2 mass  ---
!
      ELSEIF( IRNV.EQ.280 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 280 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 280
!            VAR = VAR + PORD(2,L)*VOL(L)*XHN(2,L)*SH(2,L)*RHOH(2,L)
            VAR = VAR + PORD(2,L)*VOL(L)*TMHN(2,L)
  280     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMLOX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMLOX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.280 .AND. IRNV.LE.299 ) THEN
!
!---  Differential integrated mobile water mass  ---
!
      IF( IRNV.EQ.281 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 281 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 281
            VAR = VAR + PORD(2,L)*VOL(L)*(XGW(2,L)*SG(2,L)*RHOG(2,L)
     &        + XNW(2,L)*SN(2,L)*RHON(2,L) + XLW(2,L)*SL(2,L)*RHOL(2,L)
     &        + SI(2,L)*RHOI(2,L))
  281     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMGOX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMGOX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated mobile CO2 mass  ---
!
      ELSEIF( IRNV.EQ.282 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 282 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 282
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XLA(2,L)*SL(2,L)*RHOL(2,L) +
     &          XGA(2,L)*SG(2,L)*RHOG(2,L) +
     &          XNA(2,L)*SN(2,L)*RHON(2,L) )
!
!---        STOMP-HYD through STOMP-HYD-KE  ---
!
            ELSEIF( IOM.GE.36 .AND. IOM.LE.38 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLW(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGW(2,L)*SG(2,L)*RHOG(2,L) +
     &           XHW(2,L)*SH(2,L)*RHOH(2,L))
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLW(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGW(2,L)*SG(2,L)*RHOG(2,L)
              IF( ISLC(5).EQ.1 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*SI(2,L)*RHOI(2,L)
              IF( IOM.EQ.7 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XNW(2,L)*SN(2,L)*RHON(2,L)
            ENDIF
  282     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMWX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMWX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated air, CO2 mass, or mobile CH4 mass  ---
!
      ELSEIF( IRNV.EQ.283 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 283 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 283
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XLO(2,L)*SL(2,L)*RHOL(2,L) +
     &          XGO(2,L)*SG(2,L)*RHOG(2,L) +
     &          XNO(2,L)*SN(2,L)*RHON(2,L) )
!
!---        STOMP-HYD through STOMP-HYD-KE  ---
!
            ELSEIF( IOM.GE.36 .AND. IOM.LE.38 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XLA(2,L)*SL(2,L)*RHOL(2,L) +
     &          XGA(2,L)*SG(2,L)*RHOG(2,L) +
     &          XHA(2,L)*SH(2,L)*RHOH(2,L) )
!
!---        STOMP-EOR  ---
!
            ELSEIF( IOM.GE.43 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XGC(IRNV_GC,2,L)*SG(2,L)*RHOG(2,L) +
     &          XLC(IRNV_GC,2,L)*SL(2,L)*RHOL(2,L) +
     &          XNC(IRNV_GC,2,L)*SN(2,L)*RHON(2,L) )
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLA(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGA(2,L)*SG(2,L)*RHOG(2,L)
              IF( ISLC(13).EQ.1 ) VAR = VAR + XGO(2,L)
            ENDIF
  283     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMAX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMAX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Differential integrated oil mass or mobile N2 mass  ---
!
      ELSEIF( IRNV.EQ.284 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 284 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 284
!
!---        STOMP-HYDT-KE  ---
!
            IF( IOM.EQ.39 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          ( XLN(2,L)*SL(2,L)*RHOL(2,L) +
     &          XGN(2,L)*SG(2,L)*RHOG(2,L) +
     &          XNN(2,L)*SN(2,L)*RHON(2,L) )
!
!---        STOMP-HYD through STOMP-HYD-KE  ---
!
            ELSEIF( IOM.GE.36 .AND. IOM.LE.38 ) THEN
              VAR = VAR + PORD(2,L)*VOL(L)*
     &          (XLO(2,L)*SL(2,L)*RHOL(2,L) +
     &           XGO(2,L)*SG(2,L)*RHOG(2,L) +
     &           XHO(2,L)*SH(2,L)*RHOH(2,L))
            ELSE
              VAR = VAR +
     &          PORD(2,L)*VOL(L)*XLO(2,L)*SL(2,L)*RHOL(2,L)
              IF( IEQA.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XGO(2,L)*SG(2,L)*RHOG(2,L)
              IF( IEQO.GT.0 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*SN(2,L)*RHON(2,L) +
     &          (1.D+0-PORT(2,L))*VOL(L)*XSO(2,L)*RHOS(IZ(L))
              IF( IOM.EQ.7 ) VAR = VAR +
     &          PORD(2,L)*VOL(L)*XNO(2,L)*SN(2,L)*RHON(2,L)
            ENDIF
  284     CONTINUE
          IF( NSTEP-NRST.EQ.0 ) THEN
            DMOX = VAR
            VAR = 0.D+0
          ELSE
            VAR = VAR-DMOX
          ENDIF
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Source-well pressure  ---
!
      ELSEIF( IRNV.EQ.285 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PVO(2,N)
!
!---  Source-well temperature  ---
!
      ELSEIF( IRNV.EQ.286 ) THEN
        IUNK = 1
        VAR = PSO(2,N)
!
!---  Aqueous surface spill height  ---
!
      ELSEIF( IRNV.EQ.287 ) THEN
        IUNM = 1
        NSP = (JD(N)-1)*IFLD + ID(N)
        VAR = HLSP(2,NSP)
!
!---  NAPL surface spill height  ---
!
      ELSEIF( IRNV.EQ.288 ) THEN
        IUNM = 1
        NSP = (JD(N)-1)*IFLD + ID(N)
        VAR = HNSP(2,NSP)
!
!---  Gas viscosity  ---
!
      ELSEIF( IRNV.EQ.289 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -1
        VAR = VISG(2,N)
!
!---  NAPL viscosity  ---
!
      ELSEIF( IRNV.EQ.290 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -1
        VAR = VISN(2,N)
!
!---  X-Node Centroid Position  ---
!
      ELSEIF( IRNV.EQ.291 ) THEN
        IUNM = 1
        VAR = XP(N)
!
!---  Y-Node Centroid Position  ---
!
      ELSEIF( IRNV.EQ.292 ) THEN
        IUNM = 1
        VAR = YP(N)
!
!---  Z-Node Centroid Position  ---
!
      ELSEIF( IRNV.EQ.293 ) THEN
        IUNM = 1
        VAR = ZP(N)
!
!---  Gas CH4 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.294 ) THEN
        VAR = YMGO(2,N)
!
!---  Hydrate-gas CH4 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.295 ) THEN
        VAR = YMHGO(2,N)
!
!---  Aqueous enthalpy  ---
!
      ELSEIF( IRNV.EQ.296 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = HL(2,N)
!
!---  Gas enthalpy  ---
!
      ELSEIF( IRNV.EQ.297 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = HG(2,N)
!
!---  NAPL enthalpy  ---
!
      ELSEIF( IRNV.EQ.298 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = HN(2,N)
!
!---  Similarity variable  ---
!
      ELSEIF( IRNV.EQ.299 ) THEN
        IUNM = 2
        IUNS = -1
        VAR = (XP(N)**2)/(TM+SMALL)
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
      IF( IRNV.GT.340 .AND. IRNV.LE.400 ) THEN
!
!---  Aqueous internal energy  ---
!
      IF( IRNV.EQ.341 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = HL(2,N)-(PG(2,N)+PATM)/RHOL(2,N)
!
!---  Gas internal energy  ---
!
      ELSEIF( IRNV.EQ.342 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = UEG(2,N)
!
!---  NAPL internal energy  ---
!
      ELSEIF( IRNV.EQ.343 ) THEN
        IUNM = 2
        IUNS = -2
        VAR = HN(2,N)-(PG(2,N)+PATM)/RHON(2,N)
!
!---  Aqueous thermal concductivity  ---
!
      ELSEIF( IRNV.EQ.344 ) THEN
        IUNKG = 1
        IUNM = 1
        IUNS = -3
        IUNK = -1
        VAR = THKL(2,N)
!
!---  Gas thermal concductivity  ---
!
      ELSEIF( IRNV.EQ.345 ) THEN
        IUNKG = 1
        IUNM = 1
        IUNS = -3
        IUNK = -1
        VAR = THKG(2,N)
!
!---  NAPL thermal concductivity  ---
!
      ELSEIF( IRNV.EQ.346 ) THEN
        IUNKG = 1
        IUNM = 1
        IUNS = -3
        IUNK = -1
        VAR = THKN(2,N)
!
!---  CO2 aqueous diffusion coefficient, 'DFLA'  ---
!
      ELSEIF( IRNV.EQ.347 ) THEN
        IUNM = 2
        IUNS = -1
        VAR = DFLA(2,N)
!
!---  Water gas diffusion coefficient, 'DFGW'  ---
!
      ELSEIF( IRNV.EQ.348 ) THEN
        IUNM = 2
        IUNS = -1
        VAR = DFGW(2,N)
!
!---  Coupled-well CO2 mass rate (kg/s),
!     or coupled-well volumetric rate (m^3/s),
!     or EOR injection well total petroleum mass rate (kg/s),
!     or EOR production well total petroleum mass rate (kg/s)  ---
!
      ELSEIF( IRNV.EQ.349 ) THEN
        IF( IOM.EQ.1 ) THEN
          IUNM = 3
          IUNS = -1
          VAR = QM_CW(3,IRNV_CW)
        ELSE
          IUNKG = 1
          IUNS = -1
          VAR = QM_CW(1,IRNV_CW)
        ENDIF
!
!---  Coupled-well CO2 mass integral (kg),
!     or coupled-well volumetric integral (m^3),
!     or EOR injection well total petroleum mass integral (kg),
!     or EOR production well total petroleum mass integral (kg)  ---
!
      ELSEIF( IRNV.EQ.350 ) THEN
        IF( IOM.EQ.1 ) THEN
          IUNM = 3
          VAR = QM_CW(4,IRNV_CW)
        ELSE
          IUNKG = 1
          VAR = QM_CW(2,IRNV_CW)
        ENDIF
!
!---  Coupled-well water mass rate (kg/s),
!     or coupled-well mass rate (kg/s),
!     or EOR injection well water mass rate (kg/s),
!     or EOR production well water mass rate (kg/s)   ---
!
      ELSEIF( IRNV.EQ.351 ) THEN
        IUNKG = 1
        IUNS = -1
        IF( IOM.EQ.1 ) THEN
          VAR = QM_CW(1,IRNV_CW)
        ELSE
          VAR = QM_CW(3,IRNV_CW)
        ENDIF
!
!---  Coupled-well water mass integral (kg),
!     or coupled-well mass integral (kg),
!     or EOR injection well water mass integral (kg),
!     or EOR production well water mass integral (kg)  ---
!
      ELSEIF( IRNV.EQ.352 ) THEN
        IUNKG = 1
        IF( IOM.EQ.1 ) THEN
          VAR = QM_CW(2,IRNV_CW)
        ELSE
          VAR = QM_CW(4,IRNV_CW)
        ENDIF
!
!---  Coupled-well air mass rate (kg/s),
!     or EOR injection well air mass rate (kg/s),
!     or EOR production well air mass rate (kg/s)   ---
!
      ELSEIF( IRNV.EQ.353 ) THEN
        IUNKG = 1
        IUNS = -1
        VAR = QM_CW(7,IRNV_CW)
!
!---  Coupled-well air mass integral (kg),
!     or EOR injection well air mass integral (kg),
!     or EOR production well air mass integral (kg)  ---
!
      ELSEIF( IRNV.EQ.354 ) THEN
        IUNKG = 1
        VAR = QM_CW(8,IRNV_CW)
!
!---  Vertically-integrated CO2 mass  ---
!
      ELSEIF( IRNV.EQ.361 ) THEN
        IUNKG = 1
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 361 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 361
          VAR = VAR + VOL(NX)*PORD(2,NX)*
     &      (SL(2,NX)*RHOL(2,NX)*XLA(2,NX) +
     &      SG(2,NX)*RHOG(2,NX)*XGA(2,NX) + 
     &      SN(2,NX)*RHON(2,NX)*XNA(2,NX))
  361   CONTINUE
!
!---  Vertically-integrated CO2 mass per area  ---
!
      ELSEIF( IRNV.EQ.362 ) THEN
        IUNKG = 1
        IUNM = -2
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 362 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 362
          VAR = VAR + VOL(NX)*PORD(2,NX)*
     &      (SL(2,NX)*RHOL(2,NX)*XLA(2,NX) +
     &      SG(2,NX)*RHOG(2,NX)*XGA(2,NX) + 
     &      SN(2,NX)*RHON(2,NX)*XNA(2,NX))/AFZ(NSZ(NX))
  362   CONTINUE
!
!---  Vertically-integrated gas CO2 mass  ---
!
      ELSEIF( IRNV.EQ.363 ) THEN
        IUNKG = 1
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 363 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 363
          VAR = VAR + VOL(NX)*PORD(2,NX)*SG(2,NX)*RHOG(2,NX)*XGA(2,NX)
  363   CONTINUE
!
!---  Vertically-integrated gas CO2 mass per area  ---
!
      ELSEIF( IRNV.EQ.364 ) THEN
        IUNKG = 1
        IUNM = -2
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 364 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 364
          VAR = VAR + VOL(NX)*PORD(2,NX)*SG(2,NX)*RHOG(2,NX)*XGA(2,NX)
     &      /AFZ(NSZ(NX))
  364   CONTINUE
!
!---  Vertically-integrated aqueous CO2 mass  ---
!
      ELSEIF( IRNV.EQ.365 ) THEN
        IUNKG = 1
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 365 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 365
          VAR = VAR + VOL(NX)*PORD(2,NX)*SL(2,NX)*RHOL(2,NX)*XLA(2,NX)
  365   CONTINUE
!
!---  Vertically-integrated aqueous CO2 mass per area  ---
!
      ELSEIF( IRNV.EQ.366 ) THEN
        IUNKG = 1
        IUNM = -2
        VAR = 0.D+0
        IX = ID(N)
        JX = JD(N)
        DO 366 KX = 1,KFLD
          NX = ND(IX,JX,KX)
          IF( IXP(NX).EQ.0 ) GOTO 366
          VAR = VAR + VOL(NX)*PORD(2,NX)*SL(2,NX)*RHOL(2,NX)*XLA(2,NX)
     &      /AFZ(NSZ(NX))
  366   CONTINUE
!
!---  Integrated precipitated salt mass  ---
!
      ELSEIF( IRNV.EQ.367 ) THEN
        IUNKG = 1
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 367 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 367
            VAR = VAR +
     &        PORD(2,L)*VOL(L)*SS(2,L)*RHOSP(2,L)
  367     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Strain-xx  ---
!
      ELSEIF( IRNV.EQ.369 ) THEN
        VAR = EPS_GM(1,N)
!
!---  Strain-yy  ---
!
      ELSEIF( IRNV.EQ.370 ) THEN
        VAR = EPS_GM(2,N)
!
!---  Strain-zz  ---
!
      ELSEIF( IRNV.EQ.371 ) THEN
        VAR = EPS_GM(3,N)
!
!---  Reference volumetric strain  ---
!
      ELSEIF( IRNV.EQ.372 ) THEN
        VAR = EPSV_CMP(N)
!
!---  Volumetric strain  ---
!
      ELSEIF( IRNV.EQ.373 ) THEN
        VAR = EPSV_GM(2,N)
!!
!!---  Strain-yz  ---
!!
!      ELSEIF( IRNV.EQ.372 ) THEN
!        VAR = EPS_GM(4,N)
!!
!!---  Strain-xz  ---
!!
!      ELSEIF( IRNV.EQ.373 ) THEN
!        VAR = EPS_GM(5,N)
!
!---  Strain-xy  ---
!
      ELSEIF( IRNV.EQ.374 ) THEN
        VAR = EPS_GM(6,N)
!
!---  X Displacement  ---
!
      ELSEIF( IRNV.EQ.375 ) THEN
        IUNM = 1
        IF( IRNV_GC.EQ.0 ) THEN
          VAR = 0.D+0
          DO L = 1,8
            NFEN = ND_GM(L,N)
            VAR = VAR + U_GM(2,NFEN) - U_GM(1,NFEN)
          ENDDO
          VAR = 1.25D-1*VAR
        ELSE
          L = -IRNV_GC
          NFEN = ND_GM(L,N)
          VAR = U_GM(2,NFEN) - U_GM(1,NFEN)
        ENDIF
!
!---  Y Displacement  ---
!
      ELSEIF( IRNV.EQ.376 ) THEN
        IUNM = 1
        IF( IRNV_GC.EQ.0 ) THEN
          VAR = 0.D+0
          DO L = 1,8
            NFEN = ND_GM(L,N)
            VAR = VAR + V_GM(2,NFEN) - V_GM(1,NFEN)
          ENDDO
          VAR = 1.25D-1*VAR
        ELSE
          L = -IRNV_GC
          NFEN = ND_GM(L,N)
          VAR = V_GM(2,NFEN) - V_GM(1,NFEN)
        ENDIF
!
!---  Z Displacement  ---
!
      ELSEIF( IRNV.EQ.377 ) THEN
        IUNM = 1
        IF( IRNV_GC.EQ.0 ) THEN
          VAR = 0.D+0
          DO L = 1,8
            NFEN = ND_GM(L,N)
            VAR = VAR + W_GM(2,NFEN) - W_GM(1,NFEN)
          ENDDO
          VAR = 1.25D-1*VAR
        ELSE
          L = -IRNV_GC
          NFEN = ND_GM(L,N)
          VAR = W_GM(2,NFEN) - W_GM(1,NFEN)
        ENDIF
!
!---  Integrated energy, J  ---
!
      ELSEIF( IRNV.EQ.378 ) THEN
        IUNKG = 1
        IUNM = 2
        IUNS = -2
        IF( JSKPX.EQ.0 ) THEN
          VAR = 0.D+0
          DO 378 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 378
            VAR = VAR + (PORD(2,L)*(SL(2,L)*RHOL(2,L)*UEL(2,L) +
     &        SG(2,L)*RHOG(2,L)*UEG(2,L) + SN(2,L)*RHON(2,L)*UEN(2,L) +
     &        SH(2,L)*RHOH(2,L)*HH(2,L) + SI(2,L)*RHOI(2,L)*HI(2,L) +
     &        SS(2,L)*RHOSP(2,L)*HSP(2,L)) +
     &        (1.D+0-PORD(2,L))*RHOS(IZ(L))*CPS(IZ(L))*T(2,L))*VOL(L)
  378     CONTINUE
          VSKPX = VAR
          JSKPX = 1
        ELSE
          VAR = VSKPX
        ENDIF
!
!---  Gas CO2 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.379 ) THEN
        VAR = YMGA(2,N)
!
!---  Hydrate-gas CO2 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.380 ) THEN
        VAR = YMHGA(2,N)
!
!---  Gas N2 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.381 ) THEN
        VAR = YMGN(2,N)
!
!---  Hydrate-gas N2 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.382 ) THEN
        VAR = YMHGN(2,N)
!
!---  Total nonaqueous CO2 mole fraction of formers
!     or total petroleum mole fraction  ---
!
      ELSEIF( IRNV.EQ.383 ) THEN
        IF( IOM.EQ.43 ) THEN
          VAR = ZMC(IRNV_GC,2,N)
        ELSE
          VAR = ZMCA(2,N)
        ENDIF
!
!---  Total nonaqueous CH4 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.384 ) THEN
        VAR = ZMCO(2,N)
!
!---  Total nonaqueous N2 mole fraction of formers  ---
!
      ELSEIF( IRNV.EQ.385 ) THEN
        VAR = ZMCN(2,N)
!
!---  Hydrate equilibrium pressure  ---
!
      ELSEIF( IRNV.EQ.386 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PVHA(2,N)+PVHO(2,N)+PVHN(2,N)
!
!---  Total mobile formers vapor pressure  ---
!
      ELSEIF( IRNV.EQ.387 ) THEN
        IUNM = -1
        IUNKG = 1
        IUNS = -2
        VAR = PVA(2,N)+PVO(2,N)+PVN(2,N)
!
!---  Saturation function index  ---
!
      ELSEIF( IRNV.EQ.394 ) THEN
        IF( IXP(N).EQ.0 ) THEN
          VAR = 0.D+0
        ELSE
          VAR = REAL( IZ2(N) )
        ENDIF
!
!---  Gas-nonaqueous liquid interfacial tension  ---
!
      ELSEIF( IRNV.EQ.396 ) THEN
        IUNKG = 1
        IUNS = -2
        VAR = GNIFT(2,N)
!
!---  Fracture mechanical aperture (m), 'APM'  ---
!
      ELSEIF( IRNV.EQ.397 ) THEN
        IRNV = -1
        VAR = -9999.D+22
!
!---  Fracture mechanical aperture (m), 'APH'  ---
!
      ELSEIF( IRNV.EQ.398 ) THEN
        IRNV = -1
        VAR = -9999.D+22
!
!---  Null Value  ---
!
      ELSE
        IRNV = -1
        VAR = -9999.D+22
      ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Solute, conservation-component species, and
!     kinetic-component species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.400 .AND. IRNV.LE.INDX ) THEN
        IF( MOD((IRNV-400),33).EQ.0 ) THEN
          NSL = ((IRNV-400)/33)
          IRNVX = 33
        ELSE
          NSL = ((IRNV-400)/33) + 1
          IRNVX = MOD((IRNV-400),33)
        ENDIF
        IF( NSL.GT.NSOLU ) IUNMOL = 1
        IF( IRNVX.EQ.1 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 600 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  600       CONTINUE
          ELSE
            VAR = C(N,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.2 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 602 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  602       CONTINUE
          ELSE
            VAR = C(N,NSL)
          ENDIF
          IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = VAR*YL(N,NSL)/(SL(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.3 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 604 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  604       CONTINUE
          ELSE
            VAR = C(N,NSL)
          ENDIF
          IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = VAR*YG(N,NSL)/(SG(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.4 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 606 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  606       CONTINUE
          ELSE
            VAR = C(N,NSL)
          ENDIF
          IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = VAR*YN(N,NSL)/(SN(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.5 ) THEN
          VAR = YL(N,NSL)
        ELSEIF( IRNVX.EQ.6 ) THEN
          IF( IOM.EQ.1 ) THEN
            VAR = 0.D+0
            DO 610 L = 1,NFBN
              IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 610
              VAR = VAR + C(L,NSL)*YL(L,NSL)*VOL(L)
  610       CONTINUE
          ELSE
            VAR = YG(N,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.7 ) THEN
          VAR = YN(N,NSL)
        ELSEIF( IRNVX.EQ.8 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.5D+0*(UC(NPX,NSL)+UC(NQX,NSL))
        ELSEIF( IRNVX.EQ.9 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.5D+0*(VC(NPY,NSL)+VC(NQY,NSL))
        ELSEIF( IRNVX.EQ.10 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.5D+0*(WC(NPZ,NSL)+WC(NQZ,NSL))
        ELSEIF( IRNVX.EQ.11 ) THEN
          VAR = SRCIC(N,NSL)
        ELSEIF( IRNVX.EQ.12 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 612 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  612       CONTINUE
          ELSE
            VAR = YL(N,NSL) + YN(N,NSL) + YG(N,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.13 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 614 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  614       CONTINUE
            IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YL(N,NSL)/(SL(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = C(N,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.14 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 616 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  616       CONTINUE
            IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YG(N,NSL)/(SG(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = CNL(N,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.15 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 618 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  618       CONTINUE
            IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YN(N,NSL)/(SN(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = YG(N,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.16 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = UC(NPX,NSL)
        ELSEIF( IRNVX.EQ.17 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = VC(NPY,NSL)
        ELSEIF( IRNVX.EQ.18 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = WC(NPZ,NSL)
        ELSEIF( IRNVX.EQ.22 ) THEN
          IUNS = -1
          VAR = YN(NPZ,NSL)
        ELSEIF( IRNVX.EQ.23 ) THEN
          VAR = 0.D+0
          DO 620 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 620
            VAR = VAR + C(L,NSL)*VOL(L)
            IF( ISLC(13).EQ.1 ) VAR = VAR + CNL(L,NSL)*VOL(L)*SN(2,L)
  620     CONTINUE
        ELSEIF( IRNVX.EQ.26 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 622 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  622       CONTINUE
          ENDIF     
        ELSEIF( IRNVX.EQ.27 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 624 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  624       CONTINUE
            IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YL(N,NSL)/(SL(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ENDIF     
        ELSEIF( IRNVX.EQ.28 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 626 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  626       CONTINUE
            IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YG(N,NSL)/(SG(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ENDIF     
        ELSEIF( IRNVX.EQ.29 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 628 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  628       CONTINUE
            IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
              VAR = VAR*YN(N,NSL)/(SN(2,N)*PORD(2,N))
            ELSE
              VAR = 0.D+0
            ENDIF
          ENDIF     
!
!---    Mass concentration ---
!
        ELSEIF( IRNVX.EQ.30 ) THEN
          IUNKG = -1
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 630 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C(N,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C(N,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  630       CONTINUE
          ELSE
            VAR = C(N,NSL)
          ENDIF
          DIVX = (1.D+0-PORD(2,N))*RHOS(IZ(N))
          VAR = VAR/(DIVX+SMALL)
!
!---    Solute activity ---
!
        ELSEIF( IRNVX.EQ.31 ) THEN
          VAR = C(N,NSL)*VOL(N)*(LOG(2.D+0)/HLF(NSL))
!
!---    Integrated solute activity ---
!
        ELSEIF( IRNVX.EQ.32 ) THEN
          VAR = 0.D+0
          DO L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) CYCLE
            VAR = VAR + C(N,NSL)*VOL(N)*(LOG(2.D+0)/HLF(NSL))
          ENDDO
        ENDIF
        IF( NSL.GT.NSOLU ) VAR = VAR*1.D-3
      ENDIF
!
!---  Reactive species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.INDX .AND. 
     &  IRNV.LE.(INDX+NSPR*33) ) THEN
        IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
          NSP = ((IRNV-INDX)/33)
          IRNVX = 33
        ELSE
          NSP = ((IRNV-INDX)/33) + 1
          IRNVX = MOD((IRNV-INDX),33)
        ENDIF
        IF( IRNVX.EQ.1 ) THEN
          IUNMOL = 1
          IUNM = -3
          VAR = SP_C(N,NSP)*1.D-3
        ELSEIF( IRNVX.EQ.2 ) THEN
          IUNMOL = 1
          IUNM = -3
          IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = SP_C(N,NSP)/(SL(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
          VAR = VAR*1.D-3
        ELSEIF( IRNVX.EQ.3 ) THEN
          IUNMOL = 1
          IUNM = -3
          IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = SP_C(N,NSP)/(SG(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
          VAR = VAR*1.D-3
        ELSEIF( IRNVX.EQ.4 ) THEN
          IUNMOL = 1
          IUNM = -3
          IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
            VAR = SP_C(N,NSP)/(SN(2,N)*PORD(2,N))
          ELSE
            VAR = 0.D+0
          ENDIF
          VAR = VAR*1.D-3
        ELSEIF( IRNVX.EQ.5 ) THEN
          VAR = YL(N,NSP)
        ELSEIF( IRNVX.EQ.6 ) THEN
          IF( IOM.EQ.1 ) THEN
            VAR = 0.D+0
            DO 710 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 710
              VAR = VAR + C(L,NSP)*YL(L,NSP)*VOL(L)
  710       CONTINUE
          ELSE
            VAR = YG(N,NSP)
          ENDIF
        ELSEIF( IRNVX.EQ.7 ) THEN
          VAR = YN(N,NSP)
        ELSEIF( IRNVX.EQ.8 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(UC(NPX,NSP)+UC(NQX,NSP))
        ELSEIF( IRNVX.EQ.9 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(VC(NPY,NSP)+VC(NQY,NSP))
        ELSEIF( IRNVX.EQ.10 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(WC(NPZ,NSP)+WC(NQZ,NSP))
        ELSEIF( IRNVX.EQ.11 ) THEN
          VAR = SRCIC(N,NSP)
        ELSEIF( IRNVX.EQ.12 ) THEN
          IUNM = -3
          VAR = YL(N,NSP) + YN(N,NSP) + YG(N,NSP)
        ELSEIF( IRNVX.EQ.13 ) THEN
          IUNM = -3
          VAR = SP_C(N,NSP)
        ELSEIF( IRNVX.EQ.14 ) THEN
          IUNM = -3
          VAR = CNL(N,NSP)
        ELSEIF( IRNVX.EQ.15 ) THEN
          IUNM = -3
          VAR = YG(N,NSP)
        ELSEIF( IRNVX.EQ.16 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(UC(NPX,NSP)+UC(NQX,NSP))
        ELSEIF( IRNVX.EQ.17 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(VC(NPY,NSP)+VC(NQY,NSP))
        ELSEIF( IRNVX.EQ.18 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 5.D-1*(WC(NPZ,NSP)+WC(NQZ,NSP))
        ELSEIF( IRNVX.EQ.22 ) THEN
          IUNS = -1
          VAR = YN(NPZ,NSP)
        ELSEIF( IRNVX.EQ.23 ) THEN
          IUNMOL = 1
          VAR = 0.D+0
          DO 720 L = 1,NFBN
            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 720
            VAR = VAR + SP_C(L,NSP)*VOL(L)
  720     CONTINUE
          VAR = VAR*1.D-3
        ELSEIF( IRNVX.EQ.24 ) THEN
          IUNM = 2
          NSP_M = NSP-NSPL
          VAR = SP_AREA(N,NSP_M)
        ELSEIF( IRNVX.EQ.25 ) THEN
          IUNMOL = 1
          IUNS = -1
          NSP_M = NSP-NSPL
          VAR = SP_RATE(N,NSP_M)
          VAR = VAR*1.D-3
        ELSEIF( IRNVX.EQ.26 ) THEN
          NSP_M = NSP-NSPL
          VAR = RS_S(3,NSP_M,N)
        ELSEIF( IRNVX.EQ.27 ) THEN
          VAR = C_PH(N)
        ENDIF
      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REFVAR group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REFVAR_BH( VAR,VSKPX,IRNV,JSKPX,NBN )
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
!     Reference node variable conversion for borehole nodes.
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 14 May 2019.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_BH
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PORMED
      USE PARM_BH
      USE GEOM_BH
      USE FLUX_BH
      USE FDVP_BH
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
!      REAL*8, DIMENSION(3) :: P1X,P2X,VNX
      REAL*8, DIMENSION(3) :: VXNX,VYNX,VZNX
!
!----------------------Data Statements---------------------------------!
!
      DATA VXNX /1.D+0,0.D+0,0.D+0/
      DATA VYNX /0.D+0,1.D+0,0.D+0/
      DATA VZNX /0.D+0,0.D+0,1.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REFVAR_BH'
!
!---  Convert reference node variables (indices 1 through 32) ---
!
      IF( IRNV.LE.32 ) THEN
!
!---    Aqueous pressure (absolute)  ---
!
        IF( IRNV.EQ.1 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PL_BH(2,NBN) + PATM
!
!---    Gas pressure (absolute)  ---
!
        ELSEIF( IRNV.EQ.2 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PG_BH(2,NBN) + PATM
!
!---    Temperature  ---
!
        ELSEIF( IRNV.EQ.4 ) THEN
          IUNK = 1
          VAR = T_BH(2,NBN)
!
!---    Phase condition  ---
!
        ELSEIF( IRNV.EQ.5 ) THEN
          VAR = REAL( NPHAZ_BH(2,NBN) )
!
!---    Aqueous pressure (gauge)  ---
!
        ELSEIF( IRNV.EQ.6 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PL_BH(2,NBN)
!
!---    Gas pressure (gauge)  ---
!
        ELSEIF( IRNV.EQ.7 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PG_BH(2,NBN)
!
!---    Actual aqueous saturation  ---
!
        ELSEIF( IRNV.EQ.11 ) THEN
          VAR = SL_BH(2,NBN)
!
!---    Actual gas saturation  ---
!
        ELSEIF( IRNV.EQ.12 ) THEN
          VAR = SG_BH(2,NBN)
!
!---    Actual total-liquid saturation  ---
!
        ELSEIF( IRNV.EQ.14 ) THEN
          VAR = SL_BH(2,NBN)
!
!---    Aqueous moisture content  ---
!
        ELSEIF( IRNV.EQ.15 ) THEN
          VAR = SL_BH(2,NBN)
!
!---    Total-liquid moisture content  ---
!
        ELSEIF( IRNV.EQ.17 ) THEN
          VAR = SL_BH(2,NBN)
!
!---    Diffusive porosity  ---
!
        ELSEIF( IRNV.EQ.20 ) THEN
          VAR = PORD_BH(2,NBN)
!
!---    Gas water mass fraction
!       or volumetric component mass fraction  ---
!
        ELSEIF( IRNV.EQ.21 ) THEN
          VAR = XGW_BH(2,NBN)
!
!---    Gas air mass fraction  ---
!
        ELSEIF( IRNV.EQ.22 ) THEN
          VAR = XGA_BH(2,NBN)
!
!---    Aqueous water mass fraction  ---
!
        ELSEIF( IRNV.EQ.24 ) THEN
          VAR = XLW_BH(2,NBN)
!
!---    Aqueous air mass fraction  ---
!
        ELSEIF( IRNV.EQ.25 ) THEN
            VAR = XLA_BH(2,NBN)
!
!---    Aqueous hydraulic head  ---
!
        ELSEIF( IRNV.EQ.27 ) THEN
          IUNM = 1
          ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
          VAR = PL_BH(2,NBN)/RHORL/GRAV + ZPX
!
!---    Gas hydraulic head  ---
!
        ELSEIF( IRNV.EQ.28 ) THEN
          IUNM = 1
          ZPX = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
          VAR = PG_BH(2,NBN)/RHORL/GRAV + ZPX
!
!---    Aqueous relative permeability ---
!
        ELSEIF( IRNV.EQ.31 ) THEN
          VAR = RKL_BH(2,NBN)
!
!---    Gas relative permeability
!       or residual-NAPL saturation ---
!
        ELSEIF( IRNV.EQ.32 ) THEN
          VAR = RKG_BH(2,NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 33 through 64) ---
!
      IF( IRNV.GT.32 .AND. IRNV.LE.64 ) THEN
!
!---    Aqueous density ---
!
        IF( IRNV.EQ.34 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = RHOL_BH(2,NBN)
!
!---    Gas density ---
!
        ELSEIF( IRNV.EQ.35 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = RHOG_BH(2,NBN)
!
!---    Total water mass ---
!
        ELSEIF( IRNV.EQ.37 ) THEN
          IUNKG = 1
          VAR = VOL_BH(NBN)*
     &      (XLW_BH(2,NBN)*SL_BH(2,NBN)*RHOL_BH(2,NBN) +
     &      XGW_BH(2,NBN)*SG_BH(2,NBN)*RHOG_BH(2,NBN))
!
!---    Total air mass  ---
!
        ELSEIF( IRNV.EQ.38 ) THEN
          IUNKG = 1
          VAR = VOL_BH(NBN)*
     &      (XLA_BH(2,NBN)*SL_BH(2,NBN)*RHOL_BH(2,NBN) +
     &      XGA_BH(2,NBN)*SG_BH(2,NBN)*RHOG_BH(2,NBN))
!
!---    Salt volumetric concentration  ---
!
        ELSEIF( IRNV.EQ.47 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = TMS_BH(2,NBN)
!
!---    Salt aqueous concentration  ---
!
        ELSEIF( IRNV.EQ.48 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XLS_BH(2,NBN)*RHOL_BH(2,NBN)
!
!---    Total salt mass  ---
!
        ELSEIF( IRNV.EQ.50 ) THEN
          IUNKG = 1
          VAR = TMS_BH(2,NBN)*VOL_BH(NBN)
!
!---    Matric potential head or Webb matching point head  ---
!
        ELSEIF( IRNV.EQ.63 ) THEN
          IUNM = 1
          VAR = MIN(PL_BH(2,NBN)-PG_BH(2,NBN),0.D+0)/RHORL/GRAV
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 65 through 96) ---
!
      IF( IRNV.GT.64 .AND. IRNV.LE.96 ) THEN
!
!---    X node-centered salt flux  ---
!
        IF( IRNV.EQ.67 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Y node-centered salt flux  ---
!
        ELSEIF( IRNV.EQ.68 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered salt flux  ---
!
        ELSEIF( IRNV.EQ.69 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Gas water mole fraction
!       or volumetric component mole fraction  ---
!
        ELSEIF( IRNV.EQ.70 ) THEN
          VAR = XMGW_BH(2,NBN)
!
!---    Gas air mole fraction
!       or gas component mole fraction  ---
!
        ELSEIF( IRNV.EQ.71 ) THEN
          VAR = XMGA_BH(2,NBN)
!
!---    Gas water concentration
!       or volumetric component concentration  ---
!
        ELSEIF( IRNV.EQ.73 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XGW_BH(2,NBN)*RHOG_BH(2,NBN)
!
!---    Gas air concentration
!       or gas component concentration  ---
!
        ELSEIF( IRNV.EQ.74 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XGA_BH(2,NBN)*RHOG_BH(2,NBN)
!
!---    Aqueous water concentration  ---
!
        ELSEIF( IRNV.EQ.76 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XLW_BH(2,NBN)*RHOL_BH(2,NBN)
!
!---    Aqueous air or gas component concentration  ---
!
        ELSEIF( IRNV.EQ.77 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XLA_BH(2,NBN)*RHOL_BH(2,NBN)
!
!---    X node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.87 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
!
!---      Loop over borehole node to node connections ---
!
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBL(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Y node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.88 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
!
!---      Loop over borehole node to node connections ---
!
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBL(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.89 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBL(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    X node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.90 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBG(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Y node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.91 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBG(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.92 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBG(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    X node-centered heat flux  ---
!
        ELSEIF( IRNV.EQ.96 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBQ(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 97 through 128) ---
!
      IF( IRNV.GT.96 .AND. IRNV.LE.128 ) THEN
!
!---    Y node-centered heat flux  ---
!
        IF( IRNV.EQ.97 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBQ(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered heat flux  ---
!
        ELSEIF( IRNV.EQ.98 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBQ(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Node number  ---
!
        ELSEIF( IRNV.EQ.100 ) THEN
          VAR = REAL(NBN)
!
!---    Salt aqueous mass fraction  ---
!
        ELSEIF( IRNV.EQ.110 ) THEN
          VAR = XLS_BH(2,NBN)
!
!---    X node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.117 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Y node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.118 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.119 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBS(1,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Water-vapor partial pressure  ---
!
        ELSEIF( IRNV.EQ.128 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PVW_BH(2,NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 129 through 160) ---
!
      IF( IRNV.GT.128 .AND. IRNV.LE.160 ) THEN
!
!---    Air partial pressure or component partial pressure  ---
!
        IF( IRNV.EQ.129 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PVA_BH(2,NBN)
!
!---    Water relative humidity  ---
!
        ELSEIF( IRNV.EQ.137 ) THEN
          VAR = PVW_BH(2,NBN)/PSW_BH(2,NBN)
!
!---    Aqueous air saturation  ---
!
        ELSEIF( IRNV.EQ.150 ) THEN
          VAR = MIN( XMLA_BH(2,NBN)*HCAW/(PVA_BH(2,NBN)+SMALL),1.D+0 )
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 161 through 200) ---
!
      IF( IRNV.GT.160 .AND. IRNV.LE.200 ) THEN
!
!---    Aqueous dissolved water mass fraction  ---
!
        IF( IRNV.EQ.168 ) THEN
          VAR = XMLW_BH(2,NBN)
!
!---    Aqueous viscosity  ---
!
        ELSEIF( IRNV.EQ.176 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -1
          VAR = VISL_BH(2,NBN)
!
!---    Matric potential  ---
!
        ELSEIF( IRNV.EQ.178 ) THEN
          IUNM = 1
          VAR = (PL_BH(2,NBN)-PG_BH(2,NBN))/RHORL/GRAV
!
!---    Integrated water mass  ---
!
        ELSEIF( IRNV.EQ.191 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + 
     &            (XLW_BH(2,MBN)*SL_BH(2,MBN)*RHOL_BH(2,MBN)
     &            + XGW_BH(2,MBN)*SG_BH(2,MBN)*RHOG_BH(2,MBN))*
     &            VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated air, CO2, or gas component mass  ---
!
        ELSEIF( IRNV.EQ.192 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + 
     &            (XLA_BH(2,MBN)*SL_BH(2,MBN)*RHOL_BH(2,MBN)
     &            + XGA_BH(2,MBN)*SG_BH(2,MBN)*RHOG_BH(2,MBN))*
     &            VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated aqueous water mass  ---
!
        ELSEIF( IRNV.EQ.194 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + XLW_BH(2,MBN)*SL_BH(2,MBN)*RHOL_BH(2,MBN)
     &            *VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated aqueous air mass  ---
!
        ELSEIF( IRNV.EQ.195 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + XLA_BH(2,MBN)*SL_BH(2,MBN)*RHOL_BH(2,MBN)
     &            *VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated gas water mass  ---
!
        ELSEIF( IRNV.EQ.197 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + XGW_BH(2,MBN)*SG_BH(2,MBN)*RHOG_BH(2,MBN)
     &            *VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated gas air mass  ---
!
        ELSEIF( IRNV.EQ.198 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + XGA_BH(2,MBN)*SG_BH(2,MBN)*RHOG_BH(2,MBN)
     &            *VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 201 through 240) ---
!
      IF( IRNV.GT.200 .AND. IRNV.LE.240 ) THEN
!
!---    Aqueous CO2 mole fraction  ---
!
        IF( IRNV.EQ.204 ) THEN
          VAR = XMLA_BH(2,NBN)
!
!---    Aqueous salt mole fraction  ---
!
        ELSEIF( IRNV.EQ.205 ) THEN
          VAR = XMLS_BH(2,NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 241 through 280) ---
!
      IF( IRNV.GT.240 .AND. IRNV.LE.280 ) THEN
!
!---    X-Direction Intrinsic Permeability  ---
!
        IF( IRNV.EQ.247 ) THEN
          IUNM = 2
          VAR = PERM(1,IZ_BH(NBN))
!
!---    Y-Direction Intrinsic Permeability  ---
!
        ELSEIF( IRNV.EQ.248 ) THEN
          IUNM = 2
          VAR = PERM(1,IZ_BH(NBN))
!
!---    Z-Direction Intrinsic Permeability  ---
!
        ELSEIF( IRNV.EQ.249 ) THEN
          IUNM = 2
          VAR = PERM(1,IZ_BH(NBN))
!
!---    Precipitated salt saturation  ---
!
        ELSEIF( IRNV.EQ.264 ) THEN
          VAR = SS_BH(2,NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 281 through 299) ---
!
      IF( IRNV.GT.280 .AND. IRNV.LE.299 ) THEN
!
!---    Gas viscosity  ---
!
        IF( IRNV.EQ.289 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -1
          VAR = VISG_BH(2,NBN)
!
!---    X-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.291 ) THEN
          IUNM = 1
          VAR = 5.D-1*(XP_BH(1,NBN)+XP_BH(2,NBN))
!
!---    Y-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.292 ) THEN
          IUNM = 1
          VAR = 5.D-1*(YP_BH(1,NBN)+YP_BH(2,NBN))
!
!---    Z-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.293 ) THEN
          IUNM = 1
          VAR = 5.D-1*(ZP_BH(1,NBN)+ZP_BH(2,NBN))
!
!---    Aqueous enthalpy  ---
!
        ELSEIF( IRNV.EQ.296 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HL_BH(2,NBN)
!
!---    Gas enthalpy  ---
!
        ELSEIF( IRNV.EQ.297 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HG_BH(2,NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 341 through 400) ---
!
      IF( IRNV.GT.340 .AND. IRNV.LE.400 ) THEN
!
!---    Aqueous internal energy  ---
!
        IF( IRNV.EQ.341 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HL_BH(2,NBN)-(PG_BH(2,NBN)+PATM)/RHOL_BH(2,NBN)
!
!---    Gas internal energy  ---
!
        ELSEIF( IRNV.EQ.342 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = UEG_BH(2,NBN)
!
!---    Aqueous thermal concductivity  ---
!
        ELSEIF( IRNV.EQ.344 ) THEN
          IUNKG = 1
          IUNM = 1
          IUNS = -3
          IUNK = -1
          VAR = THKL_BH(2,NBN)
!
!---    Gas thermal concductivity  ---
!
        ELSEIF( IRNV.EQ.345 ) THEN
          IUNKG = 1
          IUNM = 1
          IUNS = -3
          IUNK = -1
          VAR = THKG_BH(2,NBN)
!
!---    CO2 aqueous diffusion coefficient, 'DFLA'  ---
!
        ELSEIF( IRNV.EQ.347 ) THEN
          IUNM = 2
          IUNS = -1
          VAR = DFLA_BH(2,NBN)
!
!---    Water gas diffusion coefficient, 'DFGW'  ---
!
        ELSEIF( IRNV.EQ.348 ) THEN
          IUNM = 2
          IUNS = -1
          VAR = DFGW_BH(2,NBN)
!
!---    Integrated precipitated salt mass  ---
!
        ELSEIF( IRNV.EQ.367 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + SS_BH(2,MBN)*RHOSP_BH(2,MBN)*
     &            VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated energy, J  ---
!
        ELSEIF( IRNV.EQ.378 ) THEN
          IUNKG = 1
          IUNM = 2
          IUNS = -2
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + 
     &          (SL_BH(2,MBN)*RHOL_BH(2,MBN)*HL_BH(2,MBN)
     &          + SG_BH(2,MBN)*RHOG_BH(2,MBN)*UEG_BH(2,MBN))*
     &          VOL_BH(MBN)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Fracture mechanical aperture (m), 'APM'  ---
!
        ELSEIF( IRNV.EQ.397 ) THEN
          IRNV = -1
          VAR = -9999.D+22
!
!---    Fracture hydraulic aperture (m), 'APH'  ---
!
        ELSEIF( IRNV.EQ.398 ) THEN
          IRNV = -1
          VAR = -9999.D+22
!
!---    Distance along borehole (m), 'DABH' ---
!
        ELSEIF( IRNV.EQ.399 ) THEN
          IUNM = 1
          VAR = DIST_BH(NBN)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Solute, conservation-component species, and
!     kinetic-component species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.400 .AND. IRNV.LE.INDX ) THEN
        IF( MOD((IRNV-400),33).EQ.0 ) THEN
          NSL = ((IRNV-400)/33)
          IRNVX = 33
        ELSE
          NSL = ((IRNV-400)/33) + 1
          IRNVX = MOD((IRNV-400),33)
        ENDIF
        IF( NSL.GT.NSOLU ) IUNMOL = 1
        IF( IRNVX.EQ.1 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 600 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_BH(NBN,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_BH(NBN,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  600       CONTINUE
          ELSE
            VAR = C_BH(NBN,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.2 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 602 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_BH(NBN,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_BH(NBN,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  602       CONTINUE
          ELSE
            VAR = C_BH(NBN,NSL)
          ENDIF
          IF( SL_BH(2,NBN).GT.SMALL ) THEN
            VAR = VAR*YL_BH(NBN,NSL)/(SL_BH(2,NBN)*PORD_BH(2,NBN))
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.3 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 604 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_BH(NBN,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_BH(NBN,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  604       CONTINUE
          ELSE
            VAR = C_BH(NBN,NSL)
          ENDIF
          IF( SG_BH(2,NBN).GT.SMALL ) THEN
            VAR = VAR*YG_BH(NBN,NSL)/(SG_BH(2,NBN)*PORD_BH(2,NBN))
          ELSE
            VAR = 0.D+0
          ENDIF
!        ELSEIF( IRNVX.EQ.4 ) THEN
!          IUNM = -3
!          IF( NSL.GT.NSOLU )THEN
!            NEQ = NSL-NSOLU
!            VAR = 0.D+0
!            IUNMOL = 1
!            DO 606 M = 1,IEQ_C(1,NEQ)
!              NSP = IEQ_C(M+1,NEQ)
!              IF( NSP.LE.NSPL ) THEN
!                IF( ABS(SP_C_BH(NBN,NSP)).LT.1.D-30 ) THEN
!                  SP_CX = 0.D+0
!                ELSE
!                  SP_CX = SP_C_BH(NBN,NSP)
!                ENDIF
!                VAR = VAR + EQ_C(M,NEQ)*SP_CX
!              ENDIF
!  606       CONTINUE
!          ELSE
!            VAR = C_BH(NBN,NSL)
!          ENDIF
!          IF( SN_BH(2,NBN).GT.SMALL ) THEN
!            VAR = VAR*YN_BH(NBN,NSL)/(SN_BH(2,NBN)*PORD_BH(2,NBN))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
        ELSEIF( IRNVX.EQ.5 ) THEN
          VAR = YL_BH(NBN,NSL)
        ELSEIF( IRNVX.EQ.6 ) THEN
          IF( IOM.EQ.1 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + C_BH(MBN,NSL)*YL_BH(MBN,NSL)*VOL_BH(MBN)
              ENDDO
            ENDDO
          ELSE
            VAR = YG_BH(NBN,NSL)
          ENDIF
!        ELSEIF( IRNVX.EQ.7 ) THEN
!          VAR = YN_BH(NBN,NSL)
!
!---    X node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.8 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBC(NSL,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Y node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.9 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBC(NSL,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Z node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.10 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          MC = 0
          DO ICX = 1,2
            NCX = IPB_BH(ICX,NBN)
            IF( NCX.EQ.0 ) CYCLE
            MC = MC + 1
            NBN2 = IBCM_BH(NCX)
            VAR = VAR + UBBC(NSL,NCX)*REAL(SIGN(1,(NBN2-NBN)))
          ENDDO
          IF( MC.GT.0 ) VAR = VAR/REAL(MC)
!
!---    Solute integrated mass  ---
!
        ELSEIF( IRNVX.EQ.23 ) THEN
          VAR = 0.D+0
!
!---      Loop over boreholes  ---
!
          DO MBH = 1,N_BH
!
!---        Loop over borehole nodes  ---
!
            DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
              VAR = VAR + C_BH(MBN,NSL)*VOL_BH(MBN)
            ENDDO
          ENDDO
!
!---    Solute activity ---
!
        ELSEIF( IRNVX.EQ.31 ) THEN
          VAR = C_BH(NBN,NSL)*VOL_BH(NBN)*(LOG(2.D+0)/HLF(NSL))
!
!---    Integrated solute activity ---
!
        ELSEIF( IRNVX.EQ.32 ) THEN
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over boreholes  ---
!
            DO MBH = 1,N_BH
!
!---          Loop over borehole nodes  ---
!
              DO MBN = ID_BH(3,MBH),ID_BH(4,MBH)
                VAR = VAR + C_BH(MBN,NSL)*VOL_BH(MBN)*
     &            (LOG(2.D+0)/HLF(NSL))
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ENDIF
        IF( NSL.GT.NSOLU ) VAR = VAR*1.D-3
      ENDIF
!
!---  Reactive species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.INDX .AND. 
     &  IRNV.LE.(INDX+NSPR*33) ) THEN
        IRNV = -1
        VAR = -9999.D+22
      ENDIF
!!
!!---  Reactive species reference-node output ---
!!
!      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
!      IF( IRNV.GT.INDX .AND. 
!     &  IRNV.LE.(INDX+NSPR*33) ) THEN
!        IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
!          NSP = ((IRNV-INDX)/33)
!          IRNVX = 33
!        ELSE
!          NSP = ((IRNV-INDX)/33) + 1
!          IRNVX = MOD((IRNV-INDX),33)
!        ENDIF
!        IF( IRNVX.EQ.1 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          VAR = SP_C(N,NSP)*1.D-3
!        ELSEIF( IRNVX.EQ.2 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SL(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.3 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SG(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.4 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SN(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.5 ) THEN
!          VAR = YL(N,NS)
!        ELSEIF( IRNVX.EQ.6 ) THEN
!          IF( IOM.EQ.1 ) THEN
!            VAR = 0.D+0
!            DO 710 L = 1,NFBN
!            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 710
!              VAR = VAR + C(L,NS)*YL(L,NS)*VOL(L)
!  710       CONTINUE
!          ELSE
!            VAR = YG(N,NS)
!          ENDIF
!        ELSEIF( IRNVX.EQ.7 ) THEN
!          VAR = YN(N,NS)
!        ELSEIF( IRNVX.EQ.8 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = UC(NPX,NS)
!        ELSEIF( IRNVX.EQ.9 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = VC(NPY,NS)
!        ELSEIF( IRNVX.EQ.10 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = WC(NPZ,NS)
!        ELSEIF( IRNVX.EQ.11 ) THEN
!          VAR = SRCIC(N,NS)
!        ELSEIF( IRNVX.EQ.12 ) THEN
!          IUNM = -3
!          VAR = YL(N,NS) + YN(N,NS) + YG(N,NS)
!        ELSEIF( IRNVX.EQ.13 ) THEN
!          IUNM = -3
!          VAR = SP_C(N,NSP)
!        ELSEIF( IRNVX.EQ.14 ) THEN
!          IUNM = -3
!          VAR = CNL(N,NS)
!        ELSEIF( IRNVX.EQ.15 ) THEN
!          IUNM = -3
!          VAR = YG(N,NS)
!        ELSEIF( IRNVX.EQ.16 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = UC(NPX,NS)
!        ELSEIF( IRNVX.EQ.17 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = VC(NPY,NS)
!        ELSEIF( IRNVX.EQ.18 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = WC(NPZ,NS)
!        ELSEIF( IRNVX.EQ.22 ) THEN
!          IUNS = -1
!          VAR = YN(NPZ,NS)
!        ELSEIF( IRNVX.EQ.23 ) THEN
!          IUNMOL = 1
!          VAR = 0.D+0
!          DO 720 L = 1,NFBN
!            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 720
!            VAR = VAR + SP_C(L,NSP)*VOL(L)
!  720     CONTINUE
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.24 ) THEN
!          IUNM = 2
!          NSP_M = NSP-NSPL
!          VAR = SP_AREA(N,NSP_M)
!        ELSEIF( IRNVX.EQ.25 ) THEN
!          IUNMOL = 1
!          IUNS = -1
!          NSP_M = NSP-NSPL
!          VAR = SP_RATE(N,NSP_M)
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.26 ) THEN
!          NSP_M = NSP-NSPL
!          VAR = RS_S(3,NSP_M,N)
!        ELSEIF( IRNVX.EQ.27 ) THEN
!          VAR = C_PH(N)
!        ENDIF
!      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REFVAR_BH group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE REFVAR_FRC( VAR,VSKPX,IRNV,IRNV_GC,JSKPX,NTX )
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
!     Reference node variable conversion for fracture triangles
!
!----------------------Authors-----------------------------------------!
!
!     Written by Mark White, PNNL, 21 September 2017.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNS_FRC
      USE TRNSPT
      USE SOLTN
      USE REACT
      USE PARM_FRC
      USE GEOM_FRC
      USE FLUX_FRC
      USE FDVT_FRC
      USE FDVS_FRC
      USE FDVP_FRC
      USE FDVP
      USE FDVN_FRC
      USE FDVG_FRC
      USE FDVGC_FRC
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8, DIMENSION(3) :: P1X,P2X,VNX
      REAL*8, DIMENSION(3) :: VXNX,VYNX,VZNX
!
!----------------------Data Statements---------------------------------!
!
      SAVE DMWX,DMAX
      DATA VXNX /1.D+0,0.D+0,0.D+0/
      DATA VYNX /0.D+0,1.D+0,0.D+0/
      DATA VZNX /0.D+0,0.D+0,1.D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/REFVAR_FRC'
!
!---  Convert reference node variables (indices 1 through 32) ---
!
      IF( IRNV.LE.32 ) THEN
!
!---    Aqueous pressure (absolute)  ---
!
        IF( IRNV.EQ.1 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PL_FRC(2,NTX) + PATM
!
!---    Gas pressure (absolute)  ---
!
        ELSEIF( IRNV.EQ.2 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PG_FRC(2,NTX) + PATM
!
!---    NAPL pressure (absolute)  ---
!
        ELSEIF( IRNV.EQ.3 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PN_FRC(2,NTX) + PATM
!
!---    Temperature  ---
!
        ELSEIF( IRNV.EQ.4 ) THEN
          IUNK = 1
          VAR = T_FRC(2,NTX)
!
!---    Phase condition  ---
!
        ELSEIF( IRNV.EQ.5 ) THEN
          VAR = REAL( NPHAZ_FRC(2,NTX) )
!
!---    Aqueous pressure (gauge)  ---
!
        ELSEIF( IRNV.EQ.6 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PL_FRC(2,NTX)
!
!---    Gas pressure (gauge)  ---
!
        ELSEIF( IRNV.EQ.7 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PG_FRC(2,NTX)
!
!---    NAPL or Liquid-CO2 pressure (gauge)  ---
!
        ELSEIF( IRNV.EQ.8 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PN_FRC(2,NTX)
!
!---    Actual aqueous saturation  ---
!
        ELSEIF( IRNV.EQ.11 ) THEN
          VAR = SL_FRC(2,NTX)
!
!---    Actual gas saturation  ---
!
        ELSEIF( IRNV.EQ.12 ) THEN
          VAR = SG_FRC(2,NTX)
!
!---    Actual NAPL saturation  ---
!
        ELSEIF( IRNV.EQ.13 ) THEN
          VAR = SN_FRC(2,NTX)
!
!---    Actual total-liquid saturation  ---
!
        ELSEIF( IRNV.EQ.14 ) THEN
          VAR = SL_FRC(2,NTX) + SN_FRC(2,NTX)
!
!---    Aqueous moisture content  ---
!
        ELSEIF( IRNV.EQ.15 ) THEN
          VAR = SL_FRC(2,NTX)*PORD_FRC(2,NTX)
!
!---    NAPL moisture content  ---
!
        ELSEIF( IRNV.EQ.16 ) THEN
          VAR = SN_FRC(2,NTX)*PORD_FRC(2,NTX)
!
!---    Total-liquid moisture content  ---
!
        ELSEIF( IRNV.EQ.17 ) THEN
          VAR = (SL_FRC(2,NTX)+SN_FRC(2,NTX))*PORD_FRC(2,NTX)
!
!---    Diffusive porosity  ---
!
        ELSEIF( IRNV.EQ.20 ) THEN
          VAR = PORD_FRC(2,NTX)
!
!---    Gas water mass fraction
!       or volumetric component mass fraction  ---
!
        ELSEIF( IRNV.EQ.21 ) THEN
          VAR = XGW_FRC(2,NTX)
!
!---    Gas air mass fraction
!       or gas component mass fraction  ---
!
        ELSEIF( IRNV.EQ.22 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = XGC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = XGA_FRC(2,NTX)
          ENDIF
!
!---    Gas oil mass fraction
!       or NAPL component mass fraction ---
!
        ELSEIF( IRNV.EQ.23 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = XNC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = XGO_FRC(2,NTX)
          ENDIF
!
!---    Aqueous water mass fraction  ---
!
        ELSEIF( IRNV.EQ.24 ) THEN
          VAR = XLW_FRC(2,NTX)
!
!---    Aqueous air or gas component mass fraction  ---
!
        ELSEIF( IRNV.EQ.25 ) THEN
          IF( IOM.EQ.30 .OR. IOM.EQ.40 ) THEN
            VAR = XLC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = XLA_FRC(2,NTX)
          ENDIF
!
!---    Aqueous hydraulic head  ---
!
        ELSEIF( IRNV.EQ.27 ) THEN
          IUNM = 1
          VAR = PL_FRC(2,NTX)/RHORL/GRAV + ZP_FRC(NTX)
!
!---    Gas hydraulic head  ---
!
        ELSEIF( IRNV.EQ.28 ) THEN
          IUNM = 1
          VAR = PG_FRC(2,NTX)/RHORL/GRAV + ZP_FRC(NTX)
!
!---    NAPL hydraulic head  ---
!
        ELSEIF( IRNV.EQ.29 ) THEN
          IUNM = 1
          VAR = PN_FRC(2,NTX)/RHORL/GRAV + ZP_FRC(NTX)
!
!---    Rock/soil type  ---
!
        ELSEIF( IRNV.EQ.30 ) THEN
          VAR = REAL( IZ_FRC(NTX) )
!
!---    Aqueous relative permeability ---
!
        ELSEIF( IRNV.EQ.31 ) THEN
          VAR = RKL_FRC(2,NTX)
!
!---    Gas relative permeability
!       or residual-NAPL saturation ---
!
        ELSEIF( IRNV.EQ.32 ) THEN
          VAR = RKG_FRC(2,NTX)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 33 through 64) ---
!
      IF( IRNV.GT.32 .AND. IRNV.LE.64 ) THEN
!
!---    NAPL relative permeability
!       or liquid-CO2 relative permeability  ---
!
        IF( IRNV.EQ.33 ) THEN
          VAR = RKN_FRC(2,NTX)
!
!---    Aqueous density ---
!
        ELSEIF( IRNV.EQ.34 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = RHOL_FRC(2,NTX)
!
!---    Gas density ---
!
        ELSEIF( IRNV.EQ.35 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = RHOG_FRC(2,NTX)
!
!---    NAPL density ---
!
        ELSEIF( IRNV.EQ.36 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = RHON_FRC(2,NTX)
!
!---    Total water mass ---
!
        ELSEIF( IRNV.EQ.37 ) THEN
          IUNKG = 1
          VAR = AF_FRC(NTX)*APM_FRC(2,NTX)*
     &      (XLW_FRC(2,NTX)*SL_FRC(2,NTX)*RHOL_FRC(2,NTX) +
     &      XGW_FRC(2,NTX)*SG_FRC(2,NTX)*RHOG_FRC(2,NTX))
!
!---    Total air mass, 
!       or total CO2 mass, 
!       or total gas-component mass
!       or total component mass  ---
!
        ELSEIF( IRNV.EQ.38 ) THEN
          IUNKG = 1
          VAR = 0.D+0
          IF( IOM.EQ.43 ) THEN
            IF( IRNV_GC.EQ.1 ) THEN
              VAR = PORD_FRC(2,NTX)*VOL_FRC(NTX)*
     &          (XLA_FRC(2,NTX)*SL_FRC(2,NTX)*RHOL_FRC(2,NTX) +
     &          XGC_FRC(IRNV_GC,2,NTX)*SG_FRC(2,NTX)*RHOG_FRC(2,NTX) + 
     &          XNC_FRC(IRNV_GC,2,NTX)*SN_FRC(2,NTX)*RHON_FRC(2,NTX))
            ELSE
              VAR = PORD_FRC(2,NTX)*VOL_FRC(NTX)*
     &          (XGC_FRC(IRNV_GC,2,NTX)*SG_FRC(2,NTX)*RHOG_FRC(2,NTX) + 
     &          XNC_FRC(IRNV_GC,2,NTX)*SN_FRC(2,NTX)*RHON_FRC(2,NTX))
            ENDIF
          ELSE
            VAR = PORD_FRC(2,NTX)*VOL_FRC(NTX)*
     &        (XLA_FRC(2,NTX)*SL_FRC(2,NTX)*RHOL_FRC(2,NTX) +
     &        XGA_FRC(2,NTX)*SG_FRC(2,NTX)*RHOG_FRC(2,NTX) + 
     &        XNA_FRC(2,NTX)*SN_FRC(2,NTX)*RHON_FRC(2,NTX))
          ENDIF
!
!---    Salt volumetric concentration  ---
!
        ELSEIF( IRNV.EQ.47 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = TMS_FRC(2,NTX)
!
!---    Salt aqueous concentration  ---
!
        ELSEIF( IRNV.EQ.48 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XLS_FRC(2,NTX)*RHOL_FRC(2,NTX)
!
!---    Aqueous-flux Courant number  ---
!
        ELSEIF( IRNV.EQ.49 ) THEN
          VAR = CRNTL_FRC(NTX)
!
!---    Total salt mass  ---
!
        ELSEIF( IRNV.EQ.50 ) THEN
          IUNKG = 1
          VAR = TMS_FRC(2,NTX)
!
!---    X aqueous volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.51 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y aqueous volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.52 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z aqueous volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.53 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    X gas volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.54 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y gas volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.55 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z gas volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.56 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    X NAPL volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.60 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFN(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y NAPL volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.61 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFN(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z NAPL volumetric flux (node centered)  ---
!
        ELSEIF( IRNV.EQ.62 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFN(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    Matric potential head or Webb matching point head  ---
!
        ELSEIF( IRNV.EQ.63 ) THEN
          IUNM = 1
          VAR = MIN(PL_FRC(2,NTX)-PG_FRC(2,NTX),0.D+0)/RHORL/GRAV
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 65 through 96) ---
!
      IF( IRNV.GT.64 .AND. IRNV.LE.96 ) THEN
!
!---    X node-centered salt flux  ---
!
        IF( IRNV.EQ.64 .OR. IRNV.EQ.67 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y node-centered salt flux  ---
!
        ELSEIF( IRNV.EQ.65 .OR. IRNV.EQ.68 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered salt flux  ---
!
        ELSEIF( IRNV.EQ.66 .OR. IRNV.EQ.69 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    Gas water mole fraction
!       or volumetric component mole fraction  ---
!
        ELSEIF( IRNV.EQ.70 ) THEN
          VAR = XMGW_FRC(2,NTX)
!
!---    Gas air mole fraction
!       or gas component mole fraction  ---
!
        ELSEIF( IRNV.EQ.71 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = XMGC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = XMGA_FRC(2,NTX)
          ENDIF
!
!---    Gas oil mole fraction
!       or NAPL component mole fraction  ---
!
        ELSEIF( IRNV.EQ.72 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = XMNC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = XMGO_FRC(2,NTX)
          ENDIF
!
!---    Gas water concentration
!       or volumetric component concentration  ---
!
        ELSEIF( IRNV.EQ.73 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XGW_FRC(2,NTX)*RHOG_FRC(2,NTX)
!
!---    Gas air concentration
!       or gas component concentration  ---
!
        ELSEIF( IRNV.EQ.74 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XGA_FRC(2,NTX)*RHOG_FRC(2,NTX)
!
!---    Aqueous water concentration  ---
!
        ELSEIF( IRNV.EQ.76 ) THEN
          IUNM = -3
          IUNKG = 1
          VAR = XLW_FRC(2,NTX)*RHOL_FRC(2,NTX)
!
!---    Aqueous air or gas component concentration  ---
!
        ELSEIF( IRNV.EQ.77 ) THEN
          IUNM = -3
          IUNKG = 1
          IF( IOM.EQ.43 ) THEN
            VAR = XLC_FRC(IRNV_GC,2,NTX)*RHOL_FRC(2,NTX)
          ELSE
            VAR = XLA_FRC(2,NTX)*RHOL_FRC(2,NTX)
          ENDIF
!
!---    Aqueous oil concentration or nonaqueous liquid
!       component concentration  ---
!
        ELSEIF( IRNV.EQ.78 ) THEN
          IUNM = -3
          IUNKG = 1
          IF( IOM.EQ.43 ) THEN
            VAR = XNC_FRC(IRNV_GC,2,NTX)*RHON_FRC(2,NTX)
          ELSE
            VAR = XLO_FRC(2,NTX)*RHOL_FRC(2,NTX)
          ENDIF
!
!---    Gas Courant number  ---
!
        ELSEIF( IRNV.EQ.79 ) THEN
          VAR = CRNTG_FRC(NTX)
!
!---    System pressure (absolute), Pa  ---
!
        ELSEIF( IRNV.EQ.80 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
!
!---      STOMP-EOR  ---
!
          IF( IOM.EQ.43 ) THEN
            VAR = PSO_FRC(2,NTX) + PATM
          ENDIF
!
!---    Aqueous matrix saturation
!       or NAPL matrix saturation
!       or Eclipse gas saturation  ---
!
        ELSEIF( IRNV.EQ.83 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = SI_FRC(1,NTX)
          ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
            VAR = SI_FRC(1,NTX)
          ELSE
            IRNV = -1
            VAR = -9999.D+22
          ENDIF
!
!---    Aqueous fracture saturation
!       or NAPL fracture saturation
!       or Eclipse oil saturation  ---
!       or Eclipse pressure  ---
!
        ELSEIF( IRNV.EQ.84 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = SI_FRC(2,NTX)
          ELSEIF( IOM.EQ.32 .OR. IOM.EQ.33 ) THEN
            VAR = SI_FRC(2,NTX)
         ELSE
           IRNV = -1
           VAR = -9999.D+22
         ENDIF 
!
!---    X node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.87 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.88 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered aqueous volumetric flux  ---
!
        ELSEIF( IRNV.EQ.89 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFL(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    X node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.90 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.91 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered gas volumetric flux  ---
!
        ELSEIF( IRNV.EQ.92 ) THEN
          IUNM = 1
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFG(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    X node-centered heat flux  ---
!
        ELSEIF( IRNV.EQ.96 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFQ(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 97 through 128) ---
!
      IF( IRNV.GT.96 .AND. IRNV.LE.128 ) THEN
!
!---    Y node-centered heat flux  ---
!
        IF( IRNV.EQ.97 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFQ(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered heat flux  ---
!
        ELSEIF( IRNV.EQ.98 ) THEN
          IUNKG = 1
          IUNS = -3
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFQ(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    Node number  ---
!
        ELSEIF( IRNV.EQ.100 ) THEN
          VAR = REAL(NTX)
!
!---    Salt aqueous mass fraction  ---
!
        ELSEIF( IRNV.EQ.110 ) THEN
          VAR = XLS_FRC(2,NTX)
!
!---    X node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.117 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.118 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered surfactant flux  ---
!
        ELSEIF( IRNV.EQ.119 ) THEN
          IUNKG = 1
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFS(1,NCX)*VDOTP(VNX,VZNX)
          ENDDO
!
!---    Water-vapor partial pressure  ---
!
        ELSEIF( IRNV.EQ.128 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PVW_FRC(2,NTX)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 129 through 160) ---
!
      IF( IRNV.GT.128 .AND. IRNV.LE.160 ) THEN
!
!---    Air partial pressure or component partial pressure  ---
!
        IF( IRNV.EQ.129 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -2
          VAR = PVA_FRC(2,NTX)
!
!---    Water relative humidity  ---
!
        ELSEIF( IRNV.EQ.137 ) THEN
          VAR = PVW_FRC(2,NTX)/PSW_FRC(2,NTX)
!
!---    Air mass source rate
!       or gas component mass source rate  ---
!
        ELSEIF( IRNV.EQ.141 ) THEN
          IUNKG = 1
          IUNS = -1
          IF( IOM.EQ.43 ) THEN
            VAR = SRCGC_FRC(IRNV_GC,2,NTX)
          ELSE
            VAR = SRCA_FRC(2,NTX)
          ENDIF
!
!---    Salt mass source rate  ---
!
        ELSEIF( IRNV.EQ.147 ) THEN
          IUNKG = 1
          IUNS = -1
          VAR = SRCS_FRC(2,NTX)
!
!---    Aqueous air saturation  ---
!
        ELSEIF( IRNV.EQ.150 ) THEN
          VAR = MIN( XMLA_FRC(2,NTX)*HCAW/(PVA_FRC(2,NTX)+SMALL),1.D+0 )
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 161 through 200) ---
!
      IF( IRNV.GT.160 .AND. IRNV.LE.200 ) THEN
!
!---    Aqueous dissolved water mass fraction  ---
!
        IF( IRNV.EQ.168 ) THEN
          VAR = XMLW_FRC(2,NTX)
!
!---    Aqueous viscosity  ---
!
        ELSEIF( IRNV.EQ.176 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -1
          VAR = VISL_FRC(2,NTX)
!
!---    Matric potential  ---
!
        ELSEIF( IRNV.EQ.178 ) THEN
          IUNM = 1
          VAR = (PL_FRC(2,NTX)-PG_FRC(2,NTX))/RHORL/GRAV
!
!---    Integrated water mass  ---
!
        ELSEIF( IRNV.EQ.191 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            (XLW_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX)
     &            + XGW_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX))
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated air, CO2, or gas component mass  ---
!
        ELSEIF( IRNV.EQ.192 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
!
!---            STOMP-EOR  ---
!
                IF( IOM.EQ.43 ) THEN
                  VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &          ( XLC_FRC(IRNV_GC,2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX) +
     &            XGC_FRC(IRNV_GC,2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX) +
     &            XNC_FRC(IRNV_GC,2,MTX)*SN_FRC(2,MTX)*RHON_FRC(2,MTX) )
                ELSE
                  VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &          ( XLA_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX) + 
     &            XGA_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX) )
                ENDIF
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated aqueous water mass  ---
!
        ELSEIF( IRNV.EQ.194 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            XLW_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated aqueous air mass  ---
!
        ELSEIF( IRNV.EQ.195 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            XLA_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated gas water mass  ---
!
        ELSEIF( IRNV.EQ.197 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            XGW_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated gas air mass  ---
!
        ELSEIF( IRNV.EQ.198 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            XGA_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 201 through 240) ---
!
      IF( IRNV.GT.200 .AND. IRNV.LE.240 ) THEN
!
!---    Aqueous CO2 mole fraction  ---
!
        IF( IRNV.EQ.204 ) THEN
          VAR = XMLA_FRC(2,NTX)
!
!---    Aqueous salt mole fraction  ---
!
        ELSEIF( IRNV.EQ.205 ) THEN
          VAR = XMLS_FRC(2,NTX)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 241 through 280) ---
!
      IF( IRNV.GT.240 .AND. IRNV.LE.280 ) THEN
!
!---    Fracture gas mass source rate, kg/s  ---
!
        IF( IRNV.EQ.245 ) THEN
          IUNKG = 1
          IUNS = -1
          VAR = SRCG_FRC(NTX)
!
!---    Fracture aqueous mass source rate, kg/s  ---
!
        ELSEIF( IRNV.EQ.246 ) THEN
          IUNKG = 1
          IUNS = -1
          VAR = SRCL_FRC(NTX)
!
!---    X-Direction Intrinsic Permeability  ---
!
        ELSEIF( IRNV.EQ.247 ) THEN
          IUNM = 2
          VAR = PERM_FRC(2,NTX)
!
!---    Y-Direction Intrinsic Permeability  ---
!
        ELSEIF( IRNV.EQ.248 ) THEN
          IUNM = 2
          VAR = PERM_FRC(2,NTX)
!
!---    Z-Direction Intrinsic Permeability  ---
!
        ELSEIF( IRNV.EQ.249 ) THEN
          IUNM = 2
          VAR = PERM_FRC(2,NTX)
!
!---    Precipitated salt saturation  ---
!
        ELSEIF( IRNV.EQ.264 ) THEN
          VAR = SS_FRC(2,NTX)
!
!---    Differential integrated mobile water mass  ---
!
        ELSEIF( IRNV.EQ.282 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR +
     &            PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            (XLW_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX) +
     &             XGW_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX))
              ENDDO
            ENDDO
            IF( NSTEP-NRST.EQ.0 ) THEN
              DMWX = VAR
              VAR = 0.D+0
            ELSE
              VAR = VAR-DMWX
            ENDIF
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Differential integrated air, CO2 mass, or mobile CH4 mass  ---
!
        ELSEIF( IRNV.EQ.283 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
!
!---            STOMP-EOR  ---
!
                IF( IOM.EQ.43 ) THEN
                  VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &          ( XGC_FRC(IRNV_GC,2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX) +
     &          XLC_FRC(IRNV_GC,2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX) +
     &          XNC_FRC(IRNV_GC,2,MTX)*SN_FRC(2,MTX)*RHON_FRC(2,MTX) )
                ELSE
                  VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &          ( XGA_FRC(2,MTX)*SG_FRC(2,MTX)*RHOG_FRC(2,MTX) +
     &          XLA_FRC(2,MTX)*SL_FRC(2,MTX)*RHOL_FRC(2,MTX) )
                ENDIF
              ENDDO
            ENDDO
            IF( NSTEP-NRST.EQ.0 ) THEN
              DMAX = VAR
              VAR = 0.D+0
            ELSE
              VAR = VAR-DMAX
            ENDIF
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 281 through 299) ---
!
      IF( IRNV.GT.280 .AND. IRNV.LE.299 ) THEN
!
!---    Gas viscosity  ---
!
        IF( IRNV.EQ.289 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -1
          VAR = VISG_FRC(2,NTX)
!
!---    NAPL viscosity  ---
!
        ELSEIF( IRNV.EQ.290 ) THEN
          IUNM = -1
          IUNKG = 1
          IUNS = -1
          VAR = VISN_FRC(2,NTX)
!
!---    X-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.291 ) THEN
          IUNM = 1
          VAR = XP_FRC(NTX)
!
!---    Y-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.292 ) THEN
          IUNM = 1
          VAR = YP_FRC(NTX)
!
!---    Z-Node Centroid Position  ---
!
        ELSEIF( IRNV.EQ.293 ) THEN
          IUNM = 1
          VAR = ZP_FRC(NTX)
!
!---    Aqueous enthalpy  ---
!
        ELSEIF( IRNV.EQ.296 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HL_FRC(2,NTX)
!
!---    Gas enthalpy  ---
!
        ELSEIF( IRNV.EQ.297 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HG_FRC(2,NTX)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Convert reference node variables (indices 341 through 400) ---
!
      IF( IRNV.GT.340 .AND. IRNV.LE.400 ) THEN
!
!---    Aqueous internal energy  ---
!
        IF( IRNV.EQ.341 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = HL_FRC(2,NTX)-(PG_FRC(2,NTX)+PATM)/RHOL_FRC(2,NTX)
!
!---    Gas internal energy  ---
!
        ELSEIF( IRNV.EQ.342 ) THEN
          IUNM = 2
          IUNS = -2
          VAR = UEG_FRC(2,NTX)
!
!---    Aqueous thermal concductivity  ---
!
        ELSEIF( IRNV.EQ.344 ) THEN
          IUNKG = 1
          IUNM = 1
          IUNS = -3
          IUNK = -1
          VAR = THKL_FRC(2,NTX)
!
!---    Gas thermal concductivity  ---
!
        ELSEIF( IRNV.EQ.345 ) THEN
          IUNKG = 1
          IUNM = 1
          IUNS = -3
          IUNK = -1
          VAR = THKG_FRC(2,NTX)
!
!---    CO2 aqueous diffusion coefficient, 'DFLA'  ---
!
        ELSEIF( IRNV.EQ.347 ) THEN
          IUNM = 2
          IUNS = -1
          VAR = DFLA_FRC(2,NTX)
!
!---    Water gas diffusion coefficient, 'DFGW'  ---
!
        ELSEIF( IRNV.EQ.348 ) THEN
          IUNM = 2
          IUNS = -1
          VAR = DFGW_FRC(2,NTX)
!
!---    Integrated precipitated salt mass  ---
!
        ELSEIF( IRNV.EQ.367 ) THEN
          IUNKG = 1
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &            SS_FRC(2,MTX)*RHOSP_FRC(2,MTX)
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Integrated energy, J  ---
!
        ELSEIF( IRNV.EQ.378 ) THEN
          IUNKG = 1
          IUNM = 2
          IUNS = -2
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + PORD_FRC(2,MTX)*VOL_FRC(MTX)*
     &           (SL_FRC(2,MTX)*RHOL_FRC(2,MTX)*UEL_FRC(2,MTX) +
     &            SG_FRC(2,MTX)*RHOG_FRC(2,MTX)*UEG_FRC(2,MTX) + 
     $            SN_FRC(2,MTX)*RHON_FRC(2,MTX)*UEN_FRC(2,MTX) +
     &            SS_FRC(2,MTX)*RHOSP_FRC(2,MTX)*HSP_FRC(2,MTX))
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
!
!---    Total nonaqueous CO2 mole fraction of formers
!       or total petroleum mole fraction  ---
!
        ELSEIF( IRNV.EQ.383 ) THEN
          IF( IOM.EQ.43 ) THEN
            VAR = ZMC_FRC(IRNV_GC,2,NTX)
          ENDIF
!
!---    Fracture mechanical aperture (m), 'APM'  ---
!
        ELSEIF( IRNV.EQ.397 ) THEN
          IUNM = 1
          VAR = APM_FRC(2,NTX)
!
!---    Fracture hydraulic aperture (m), 'APH'  ---
!
        ELSEIF( IRNV.EQ.398 ) THEN
          IUNM = 1
          VAR = APH_FRC(2,NTX)
!
!---    Null Value  ---
!
        ELSE
          IRNV = -1
          VAR = -9999.D+22
        ENDIF
        ISUB_LOG = ISUB_LOG-1
        RETURN
      ENDIF
!
!---  Solute, conservation-component species, and
!     kinetic-component species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.400 .AND. IRNV.LE.INDX ) THEN
        IF( MOD((IRNV-400),33).EQ.0 ) THEN
          NSL = ((IRNV-400)/33)
          IRNVX = 33
        ELSE
          NSL = ((IRNV-400)/33) + 1
          IRNVX = MOD((IRNV-400),33)
        ENDIF
        IF( NSL.GT.NSOLU ) IUNMOL = 1
        IF( IRNVX.EQ.1 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 600 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  600       CONTINUE
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.2 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 602 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  602       CONTINUE
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF
          IF( SL_FRC(2,NTX).GT.SMALL ) THEN
            VAR = VAR*YL_FRC(NTX,NSL)/SL_FRC(2,NTX)
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.3 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 604 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  604       CONTINUE
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF
          IF( SG_FRC(2,NTX).GT.SMALL ) THEN
            VAR = VAR*YG_FRC(NTX,NSL)/SG_FRC(2,NTX)
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.4 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO 606 M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.LE.NSPL ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
  606       CONTINUE
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF
          IF( SN_FRC(2,NTX).GT.SMALL ) THEN
            VAR = VAR*YN_FRC(NTX,NSL)/SN_FRC(2,NTX)
          ELSE
            VAR = 0.D+0
          ENDIF
        ELSEIF( IRNVX.EQ.5 ) THEN
          VAR = YL_FRC(NTX,NSL)
        ELSEIF( IRNVX.EQ.6 ) THEN
          IF( IOM.EQ.1 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO MFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,MFX),IP_FRC(2,MFX)
                VAR = VAR + VOL_FRC(MTX)*
     &            C_FRC(MTX,NSL)*YL_FRC(MTX,NSL)
              ENDDO
            ENDDO
          ELSE
            VAR = YG_FRC(NTX,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.7 ) THEN
          VAR = YN_FRC(NTX,NSL)
!
!---    X node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.8 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFC(NSL,NCX)*VDOTP(VNX,VXNX)
          ENDDO
!
!---    Y node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.9 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFC(NSL,NCX)*VDOTP(VNX,VYNX)
          ENDDO
!
!---    Z node-centered solute flux  ---
!
        ELSEIF( IRNVX.EQ.10 ) THEN
          IUNM = -2
          IUNS = -1
          VAR = 0.D+0
          P1X(1) = XP_FRC(NTX)
          P1X(2) = YP_FRC(NTX)
          P1X(3) = ZP_FRC(NTX)
          DO NCX = IPF_FRC(1,NTX),IPF_FRC(2,NTX)
            NT2X = ITCM_FRC(NCX)
            P2X(1) = XP_FRC(NT2X)
            P2X(2) = YP_FRC(NT2X)
            P2X(3) = ZP_FRC(NT2X)
            CALL VSUB( P2X,P1X,VNX )
            VLX = SQRT(VNX(1)**2 + VNX(2)**2 + VNX(3)**2)
            DO M = 1,3
              VNX(M) = VNX(M)/VLX
            ENDDO
            VAR = VAR + UFFC(NSL,NCX)*VDOTP(VNX,VZNX)
          ENDDO
        ELSEIF( IRNVX.EQ.12 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
          ELSE
            VAR = YL_FRC(NTX,NSL) + YN_FRC(NTX,NSL) + YG_FRC(NTX,NSL)
          ENDIF
        ELSEIF( IRNVX.EQ.13 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SL_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YL_FRC(NTX,NSL)/(SL_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.14 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SG_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YG_FRC(NTX,NSL)/(SG_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            IRNV = -1
            VAR = -9999.D+22
          ENDIF     
        ELSEIF( IRNVX.EQ.15 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL+NSPS .AND. NSP.LE.NSPL+NSPS+NSPE ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SN_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YN_FRC(NTX,NSL)/(SN_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = YG_FRC(NTX,NSL)
          ENDIF     
!
!---    Solute integrated mass  ---
!
        ELSEIF( IRNVX.EQ.23 ) THEN
          VAR = 0.D+0
!
!---      Loop over fractures  ---
!
          DO NFX = 1,NF_FRC
!
!---        Loop over fracture triangles  ---
!
            DO MTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
              VAR = VAR + VOL_FRC(MTX)*C_FRC(MTX,NSL)
            ENDDO
          ENDDO
        ELSEIF( IRNVX.EQ.26 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
          ENDIF     
        ELSEIF( IRNVX.EQ.27 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SL_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YL_FRC(NTX,NSL)/(SL_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.28 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SG_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YG_FRC(NTX,NSL)/(SG_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF     
        ELSEIF( IRNVX.EQ.29 ) THEN
          IUNM = -3
          IF( NSL.GT.NSOLU )THEN
            NEQ = NSL-NSOLU
            VAR = 0.D+0
            IUNMOL = 1
            DO M = 1,IEQ_C(1,NEQ)
              NSP = IEQ_C(M+1,NEQ)
              IF( NSP.GT.NSPL .AND. NSP.LE.NSPL+NSPS ) THEN
                IF( ABS(SP_C_FRC(NTX,NSP)).LT.1.D-30 ) THEN
                  SP_CX = 0.D+0
                ELSE
                  SP_CX = SP_C_FRC(NTX,NSP)
                ENDIF
                VAR = VAR + EQ_C(M,NEQ)*SP_CX
              ENDIF
            ENDDO
            IF( SN_FRC(2,NTX)*PORD_FRC(2,NTX).GT.SMALL ) THEN
              VAR = VAR*YN_FRC(NTX,NSL)/(SN_FRC(2,NTX)*PORD_FRC(2,NTX))
            ELSE
              VAR = 0.D+0
            ENDIF
          ELSE
            VAR = C_FRC(NTX,NSL)
          ENDIF     
!
!---    Solute activity ---
!
        ELSEIF( IRNVX.EQ.31 ) THEN
          VAR = C_FRC(NTX,NSL)*VOL_FRC(NTX)*(LOG(2.D+0)/HLF(NSL))
!
!---    Integrated solute activity ---
!
        ELSEIF( IRNVX.EQ.32 ) THEN
          IF( JSKPX.EQ.0 ) THEN
            VAR = 0.D+0
!
!---        Loop over fractures  ---
!
            DO NFX = 1,NF_FRC
!
!---          Loop over fracture triangles  ---
!
              DO MTX = IP_FRC(1,NFX),IP_FRC(2,NFX)
                VAR = VAR + C_FRC(MTX,NSL)*VOL_FRC(MTX)*
     &            (LOG(2.D+0)/HLF(NSL))
              ENDDO
            ENDDO
            VSKPX = VAR
            JSKPX = 1
          ELSE
            VAR = VSKPX
          ENDIF
        ENDIF
        IF( NSL.GT.NSOLU ) VAR = VAR*1.D-3
      ENDIF
!
!---  Reactive species reference-node output ---
!
      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
      IF( IRNV.GT.INDX .AND. 
     &  IRNV.LE.(INDX+NSPR*33) ) THEN
        IRNV = -1
        VAR = -9999.D+22
      ENDIF
!!
!!---  Reactive species reference-node output ---
!!
!      INDX = (400+(NSOLU*33)+((NEQC+NEQK)*33))
!      IF( IRNV.GT.INDX .AND. 
!     &  IRNV.LE.(INDX+NSPR*33) ) THEN
!        IF( MOD((IRNV-INDX),33).EQ.0 ) THEN
!          NSP = ((IRNV-INDX)/33)
!          IRNVX = 33
!        ELSE
!          NSP = ((IRNV-INDX)/33) + 1
!          IRNVX = MOD((IRNV-INDX),33)
!        ENDIF
!        IF( IRNVX.EQ.1 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          VAR = SP_C(N,NSP)*1.D-3
!        ELSEIF( IRNVX.EQ.2 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SL(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SL(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.3 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SG(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SG(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.4 ) THEN
!          IUNMOL = 1
!          IUNM = -3
!          IF( SN(2,N)*PORD(2,N).GT.SMALL ) THEN
!            VAR = SP_C(N,NSP)/(SN(2,N)*PORD(2,N))
!          ELSE
!            VAR = 0.D+0
!          ENDIF
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.5 ) THEN
!          VAR = YL(N,NS)
!        ELSEIF( IRNVX.EQ.6 ) THEN
!          IF( IOM.EQ.1 ) THEN
!            VAR = 0.D+0
!            DO 710 L = 1,NFBN
!            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 710
!              VAR = VAR + C(L,NS)*YL(L,NS)*VOL(L)
!  710       CONTINUE
!          ELSE
!            VAR = YG(N,NS)
!          ENDIF
!        ELSEIF( IRNVX.EQ.7 ) THEN
!          VAR = YN(N,NS)
!        ELSEIF( IRNVX.EQ.8 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = UC(NPX,NS)
!        ELSEIF( IRNVX.EQ.9 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = VC(NPY,NS)
!        ELSEIF( IRNVX.EQ.10 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = WC(NPZ,NS)
!        ELSEIF( IRNVX.EQ.11 ) THEN
!          VAR = SRCIC(N,NS)
!        ELSEIF( IRNVX.EQ.12 ) THEN
!          IUNM = -3
!          VAR = YL(N,NS) + YN(N,NS) + YG(N,NS)
!        ELSEIF( IRNVX.EQ.13 ) THEN
!          IUNM = -3
!          VAR = SP_C(N,NSP)
!        ELSEIF( IRNVX.EQ.14 ) THEN
!          IUNM = -3
!          VAR = CNL(N,NS)
!        ELSEIF( IRNVX.EQ.15 ) THEN
!          IUNM = -3
!          VAR = YG(N,NS)
!        ELSEIF( IRNVX.EQ.16 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = UC(NPX,NS)
!        ELSEIF( IRNVX.EQ.17 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = VC(NPY,NS)
!        ELSEIF( IRNVX.EQ.18 ) THEN
!          IUNM = -2
!          IUNS = -1
!          VAR = WC(NPZ,NS)
!        ELSEIF( IRNVX.EQ.22 ) THEN
!          IUNS = -1
!          VAR = YN(NPZ,NS)
!        ELSEIF( IRNVX.EQ.23 ) THEN
!          IUNMOL = 1
!          VAR = 0.D+0
!          DO 720 L = 1,NFBN
!            IF( IXP(L).EQ.0 .OR. IBR(4,L).NE.L ) GOTO 720
!            VAR = VAR + SP_C(L,NSP)*VOL(L)
!  720     CONTINUE
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.24 ) THEN
!          IUNM = 2
!          NSP_M = NSP-NSPL
!          VAR = SP_AREA(N,NSP_M)
!        ELSEIF( IRNVX.EQ.25 ) THEN
!          IUNMOL = 1
!          IUNS = -1
!          NSP_M = NSP-NSPL
!          VAR = SP_RATE(N,NSP_M)
!          VAR = VAR*1.D-3
!        ELSEIF( IRNVX.EQ.26 ) THEN
!          NSP_M = NSP-NSPL
!          VAR = RS_S(3,NSP_M,N)
!        ELSEIF( IRNVX.EQ.27 ) THEN
!          VAR = C_PH(N)
!        ENDIF
!      ENDIF
!
!---  Reset subroutine character string ---
!
      ISUB_LOG = ISUB_LOG-1
!
!---  End of REFVAR_FRC group ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE WROBDA_WALLUA
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
!     Write pressure differential for well #3 at observation times
!     to the observation data file (wallula_sim.dat) for PEST output
!     to determine parameters for the CO2 mineralization at Wallula
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, PNNL, 11 August 2020.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE UCODE
      USE SOLTN
      USE FILES
      USE COUP_WELL
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL(KIND=DP), DIMENSION(164) ::  OBSTM
!
!----------------------Data Statements---------------------------------!
!
      DATA OBSTM /0.00000000D+0, 0.00504108D+0, 0.00740319D+0, 
     &  0.00974609D+0, 0.01454713D+0, 0.01474877D+0, 0.02465810D+0, 
     &  0.02708742D+0, 0.03211890D+0, 0.03715038D+0, 0.04219146D+0, 
     &  0.04726135D+0, 0.06273028D+0, 0.08348034D+0, 0.09905488D+0, 
     &  0.11201767D+0, 0.12224387D+0, 0.13267170D+0, 0.16895790D+0, 
     &  0.18720183D+0, 0.20540734D+0, 0.20958424D+0, 0.21167748D+0,
     &  0.21668975D+0, 0.22953732D+0, 0.25273590D+0, 0.27333233D+0,
     &  0.31203824D+0, 0.34560704D+0, 0.38178763D+0, 0.42059916D+0,
     &  0.46456700D+0, 0.50335940D+0, 0.55000615D+0, 0.59924555D+0,
     &  0.64069766D+0, 0.67433370D+0, 0.71057190D+0, 0.74681010D+0,
     &  0.77009505D+0, 0.79080670D+0, 0.80749510D+0, 0.81797093D+0,
     &  0.83100090D+0, 0.84658510D+0, 0.85437240D+0, 0.85868370D+0,
     &  0.86632690D+0, 0.86882347D+0, 0.88432120D+0, 0.90240186D+0,
     &  0.92306554D+0, 0.94371957D+0, 0.96698534D+0, 0.98507565D+0,
     &  1.01610950D+0, 1.04972640D+0, 1.10146220D+0, 1.13767160D+0,
     &  1.15317890D+0, 1.17127880D+0, 1.19716600D+0, 1.22564580D+0,
     &  1.24115310D+0, 1.24885390D+0, 1.27474100D+0, 1.30579410D+0,
     &  1.32907920D+0, 1.33700080D+0, 1.34739030D+0, 1.36283040D+0,
     &  1.36792910D+0, 1.39901100D+0, 1.43008330D+0, 1.45854380D+0,
     &  1.48959680D+0, 1.52842760D+0, 1.56206360D+0, 1.60090400D+0,
     &  1.63713250D+0, 1.68115810D+0, 1.72259090D+0, 1.75624610D+0,
     &  1.79508650D+0, 1.82874170D+0, 1.85203640D+0, 1.85732710D+0,
     &  1.88576840D+0, 1.88581630D+0, 1.90124690D+0, 1.92710530D+0,
     &  1.96077010D+0, 2.00219350D+0, 2.03844120D+0, 2.06949420D+0,
     &  2.09538150D+0, 2.11607380D+0, 2.15232160D+0, 2.17819930D+0,
     &  2.19632800D+0, 2.22740000D+0, 2.25585100D+0, 2.28432130D+0,
     &  2.31023720D+0, 2.33869770D+0, 2.37236260D+0, 2.40861030D+0,
     &  2.44485800D+0, 2.48368880D+0, 2.50448680D+0, 2.51235100D+0,
     &  2.53564550D+0, 2.55375500D+0, 2.56933930D+0, 2.57173010D+0,
     &  2.59253790D+0, 2.61844420D+0, 2.65210900D+0, 2.66760680D+0,
     &  2.67527900D+0, 2.69602900D+0, 2.71931390D+0, 2.73482100D+0,
     &  2.76070830D+0, 2.76084280D+0, 2.79452680D+0, 2.82560870D+0,
     &  2.83855220D+0, 2.84622430D+0, 2.85930230D+0, 2.87742140D+0,
     &  2.87990830D+0, 2.89542530D+0, 2.89811400D+0, 2.90576670D+0,
     &  2.92389540D+0, 2.93176910D+0, 2.95510220D+0, 2.97818540D+0,
     &  2.97838700D+0, 3.03255250D+0, 3.03533700D+0, 3.04053160D+0,
     &  3.04817490D+0, 3.06888650D+0, 3.09994940D+0, 3.13103100D+0,
     &  3.16990020D+0, 3.20096300D+0, 3.23461820D+0, 3.24765780D+0,
     &  3.26314600D+0, 3.27085640D+0, 3.29933620D+0, 3.32523270D+0,
     &  3.34336160D+0, 3.35894560D+0, 3.37705520D+0, 3.40034990D+0,
     &  3.42365400D+0, 3.44955090D+0, 3.47026250D+0, 3.48585600D+0,
     &  3.49881890D+0/
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/WROBDA_WALLUA'
!
!---  Store old and current time step values of the pressure
!     differential at well #3
!
      R_OBDT(5,1) = R_OBDT(6,1)
      R_OBDT(6,1) = P_CW(2,3)-R_OBDT(4,1)
!
!---  Interpolate well pressure differential (MPa) at observation 
!     times (day)  ---
!
      TMX = (TM-R_OBDT(1,1))/24.D+0/3600.D+0
      TMOX = (TM-DT-R_OBDT(1,1))/24.D+0/3600.D+0
      DO N = 1,164
        IF( TMOX.LE.OBSTM(N) .AND. TMX.GE.OBSTM(N) ) THEN
          VARX = (OBSTM(N)-TMOX)*(R_OBDT(6,1)-R_OBDT(5,1))/(TMX-TMOX)
     &      + R_OBDT(5,1)
          WRITE(46,'(1PE14.7)') 1.D-6*VARX 
        ENDIF
      ENDDO
      ISUB_LOG = ISUB_LOG-1
!
!---  End of WROBDA_WALLUA group  ---
!
      RETURN
      END

