#-------------------------Disclaimer------------------------------------
#
#     This material was prepared as an account of work sponsored by
#     an agency of the United States Government. Neither the
#     United States Government nor the United States Department of
#     Energy, nor Battelle, nor any of their employees, makes any
#     warranty, express or implied, or assumes any legal liability or
#     responsibility for the accuracy, completeness, or usefulness
#     of any information, apparatus, product, software or process
#     disclosed, or represents that its use would not infringe
#     privately owned rights.
#
#----------------------Acknowledgement----------------------------------
#
#     This software and its documentation were produced with Government
#     support under Contract Number DE-AC06-76RLO-1830 awarded by the
#     United Department of Energy. The Government retains a paid-up
#     non-exclusive, irrevocable worldwide license to reproduce,
#     prepare derivative works, perform publicly and display publicly
#     by or for the Government, including the right to distribute to
#     other Government contractors.
#
#---------------------Copyright Notices---------------------------------
#
#            Copyright © Battelle Memorial Institute, 1996.
#                    All Rights Reserved.
#     makefile https://stash.pnnl.gov/scm/stomp/stomp.git V3.0

#
#---  Command-line Directives  ---
#

decay = radioactive
storage = nada
wellbore = nada
refinement = nada
special = nada
compiler_opt = nada
chem = nada
barrier = nada
fracture = nada
mech = nada
guest_exchange = equilibrium
library = dynamic

#
#---  Environment Variables  ---
#
#  The user must set the following environment variables, otherwise Mark's default paths will be used
#
# $STOMP_COMPILER [gnu,intel,absoft]
# $STOMP_SRC_DIR [path to stomp source code]
# $STOMP_BLAS_DIR [path to blas libraries - see compiler sections below for examples]
# $STOMP_LIS_DIR or $STOMP_SP_DIR depending on solver [path to solver libraries]
# if using the petsc solver $PETSC_DIR and $PETSC_ARCH must be set
#
STOMP_COMPILER ?= gnu
compiler ?= $(STOMP_COMPILER)
STOMP_SRC_DIR ?= /Users/d3c002/stomp/src_brev
#STOMP_SRC_DIR ?= /Users/d3e144/Documents/Projects/sToMp/src_brev
DIR = $(STOMP_SRC_DIR)
STOMP_BLAS_DIR ?= -L /usr/lib -llapack -lblas
BLASLIB = $(STOMP_BLAS_DIR)
STOMP_SP_DIR ?= /Users/d3c002/splib
PETSC_DIR ?= /Users/d3c002/petsc
PETSC_ARCH ?= arch-darwin-serial-opt
#
#---  Definitions  ---
#

ifeq ($(solver),splib)
  solver = sp
endif
ifeq ($(solver),PETSc)
  solver = petsc
endif
ifeq ($(solver),Lis)
  solver = lis
endif
ifeq ($(solver),banded)
  solver = bd
endif
ifndef solver
  $(error Error: Unrecognized solver: Please use [solver=sp, solver=lis, solver=bd or solver=petsc])
endif
ifneq ($(solver),bd)
  ifneq ($(solver),sp)
    ifneq ($(solver),lis)
      ifneq ($(solver),petsc)
        $(error Error: Unrecognized solver: Please use [solver=sp, solver=lis, solver=bd or solver=petsc])
      endif
    endif
  endif
endif
ifneq ($(decay),radioactive)
  ifneq ($(decay),reactive)
    $(error Error: Unrecognized decay: Please use [decay=radioactive or decay=reactive])
  endif
endif
ifeq ($(geochem),ecke)
  chem = ecke
endif
ifeq ($(geochemistry),ecke)
  chem = ecke
endif
ifeq ($(chemistry),ecke)
  chem = ecke
endif
ifeq ($(chem),ecke)
  chem = ecke
endif
ifeq ($(geochem),eckechem)
  chem = ecke
endif
ifeq ($(geochemistry),eckechem)
  chem = ecke
endif
ifeq ($(chemistry),eckechem)
  chem = ecke
endif
ifeq ($(chem),eckechem)
  chem = ecke
endif
ifeq ($(chem),coupled_oil_shale)
  chem = coupled-oil-shale
endif
ifeq ($(chem),lagged_oil_shale)
  chem = lagged-oil-shale
endif
ifeq ($(guest_exchange),equilibrium)
  guest_exchange = eq_hyd
endif
ifeq ($(guest_exchange),equil)
  guest_exchange = eq_hyd
endif
ifeq ($(guest_exchange),eq)
  guest_exchange = eq_hyd
endif
ifeq ($(guest_exchange),equil_hyd)
  guest_exchange = eq_hyd
endif
ifeq ($(guest_exchange),k1)
  guest_exchange = k1_hyd
endif
ifeq ($(guest_exchange),kinetic 1)
  guest_exchange = k1_hyd
endif
ifeq ($(guest_exchange),k2)
  guest_exchange = k2_hyd
endif
ifeq ($(guest_exchange),kinetic 2)
  guest_exchange = k2_hyd
endif

#
#---  Set Petsc flags ---
#

ifeq ($(solver),petsc)

  ifeq ($(compiler),absoft)
    ifeq ($(library),static)
      GNU_LINKING = /usr/local/lib/libgfortran.a /usr/local/lib/libquadmath.a /usr/local/lib/gcc/x86_64-apple-darwin15.6.0/7.1.0/libgcc.a
    else
      GNU_LINKING = -L /usr/local/lib -lgfortran
    endif
  endif

  ifeq ($(compiler),intel)
    ifeq ($(library),static)
      GNU_LINKING = /usr/local/lib/libgfortran.a /usr/local/lib/libquadmath.a /usr/local/lib/gcc/x86_64-apple-darwin15.6.0/7.1.0/libgcc.a
    else
      GNU_LINKING = -L /usr/local/lib -lgfortran
    endif
  endif

  ifeq ($(compiler),intel_duo)
    ifeq ($(library),static)
      GNU_LINKING = /usr/local/lib/libgfortran.a /usr/local/lib/libquadmath.a /usr/local/lib/gcc/x86_64-apple-darwin15.6.0/7.1.0/libgcc.a
    else
      GNU_LINKING = -L /usr/local/lib -lgfortran
    endif
  endif
  include $(PETSC_DIR)/lib/petsc/conf/variables
  include $(PETSC_DIR)/lib/petsc/conf/rules
  PETSC_CPP_FLAGS = $(FCPPFLAGS) -I$(PETSC_DIR) -I$(PETSC_DIR)/$(PETSC_ARCH)/include
endif

#
#---  Set splib directory  ---
#
ifeq ($(solver),sp)
  ifeq ($(compiler),gnu)
     STOMP_SP_DIR = /Users/d3c002/splib_gnu
     SPDIR = $(STOMP_SP_DIR)
  else ifeq ($(compiler),intel)
     STOMP_SP_DIR = /Users/d3c002/splib_intel
     SPDIR = $(STOMP_SP_DIR)
  else ifeq ($(compiler),absoft)
     STOMP_SP_DIR = /Users/d3c002/splib_absoft
     SPDIR = $(STOMP_SP_DIR)
  else 
     SPDIR = $(STOMP_SP_DIR)
  endif    
endif
#
#
#---  Set lis directory  ---
#
ifeq ($(solver),lis)
  ifeq ($(compiler),gnu)
    STOMP_LIS_DIR ?= /Users/d3c002/lis-2.0.24/lis_gnu
    LIS_DIR = $(STOMP_LIS_DIR)
  else
    STOMP_LIS_DIR ?= /Users/d3c002/lis-2.0.24/lis_gnu
    LIS_DIR = $(STOMP_LIS_DIR)
    ifeq ($(library),static)
      GNU_LINKING = /usr/local/lib/libgfortran.a /usr/local/lib/libquadmath.a /usr/local/lib/gcc/x86_64-apple-darwin15.6.0/7.1.0/libgcc.a
    else
      GNU_LINKING = -L /usr/local/lib -lgfortran
    endif
  endif
  LIS_CPP_FLAGS = $(FCPPFLAGS) -I$(LIS_DIR) 
endif

#
#---  C Preprocessor  ---
#

CPP_FLAGS = -w -P -traditional -D$(solver) -D$(library)
ifneq ($(decay),nada)
  CPP_FLAGS += -D$(decay)
endif
ifneq ($(storage),nada)
  CPP_FLAGS += -D$(storage)
endif
ifneq ($(refinement),nada)
  CPP_FLAGS += -D$(refinement)
endif
ifneq ($(wellbore),nada)
  CPP_FLAGS += -D$(wellbore)
endif
ifneq ($(special),nada)
  CPP_FLAGS += -D$(special)
endif
ifneq ($(compiler),nada)
  CPP_FLAGS += -D$(compiler)
endif
ifneq ($(compiler_opt),nada)
  CPP_FLAGS += -D$(compiler_opt)
endif
ifneq ($(chem),nada)
  CPP_FLAGS += -D$(chem)
endif
ifneq ($(barrier),nada)
  CPP_FLAGS += -D$(barrier)
endif
ifneq ($(fracture),nada)
  CPP_FLAGS += -D$(fracture)
endif
ifneq ($(mech),nada)
  CPP_FLAGS += -D$(mech)
endif
ifneq ($(guest_exchange),nada)
  CPP_FLAGS += -D$(guest_exchange)
endif
ifeq ($(solver),petsc)
  CPP_FLAGS += $(PETSC_CPP_FLAGS)
endif
ifeq ($(solver),lis)
  CPP_FLAGS += $(LIS_CPP_FLAGS)
endif
CPP_FORM = $< $@

#
#---  Source Paths  ---
#

SRCH_PATH =  -I . -I $(DIR)
ifeq ($(solver),petsc)
  SRCH_PATH += -I $(PETSC_DIR)/$(PETSC_ARCH)/include
endif

ifneq ($(solver),petsc)
#
#---  Fortran Compilers  ---
#

  FC_INTEL = ifort -m64
  FC_ABSOFT = f90
  FC_LAHEY = lf95
  FC_IBM = xlf -qextname -qinitauto
  FC_GNU = gfortran
  FC_F95 = f95
  FC_UTAH = /usr/bin/gfortran-4.2
  
#
#--- Intel Fortran Compiler Options  ---
#
  
  ifeq ($(compiler),intel)
    FC = $(FC_INTEL)
    DEBUG_FLAG = -g -C -check -traceback -fpe0
    OPTIMIZE_FLAG = -O2
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = -L/share/apps/intel/xe/composer_xe_2011_sp1.9.293/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lX11 -lpthread
#    STOMP_BLAS_DIR = /opt/intel/mkl/lib/32/libmkl_lapack.a /opt/intel/mkl/lib/32/libmkl_p3.a
  endif
  
#
#--- Absoft Fortran Compiler Options  ---
#
  
  ifeq ($(compiler),absoft)
    FC = $(FC_ABSOFT)
    STOMP_BLAS_DIR = /Applications/Absoft20.0/lib64/librefblas.a
    BLASLIB = $(STOMP_BLAS_DIR)
    DEBUG_FLAG = -g -m64
    STANDARD_FLAG = -m64
    OPTIMIZE_FLAG = -O2 -m64
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = /Users/d3c002/splib/absoft/libblas.a
  endif
  
#
#--- Lahey Fortran Compiler Options  ---
#
  
  ifeq ($(compiler),lahey)
    FC = $(FC_LAHEY)
    DEBUG_FLAG = -g -C
    STANDARD_FLAG = -f
    OPTIMIZE_FLAG = -O
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = /usr/local/lf9561/lib/libblasmt.a
  endif
  
#
#--- IBM Fortran Compiler Options  ---
#
  
  ifeq ($(compiler),ibm)
    FC = $(FC_IBM)
    DEBUG_FLAG = -g -C -qflttrap -qsigtrap
    OPTIMIZE_FLAG = -O -qmaxmem=-1
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = /System/Library/Frameworks/vecLib.framework/vecLib $(SP_DIR)/blas/lsame.o
    CPP_FORM = $< | sed '/^\#pragma/d' > $@
  endif
  
# 
#--- Intel_Duo Fortran Compiler Options  ---
# 
    
  ifeq ($(compiler),intel_duo)
    FC = $(FC_INTEL)
    DEBUG_FLAG = -g -C -check -traceback -fpe0
    OPTIMIZE_FLAG = -O2
    CPP_FORM = $< | sed '/^\#pragma/d' > $@
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = -L /usr/lib -llapack -lblas
  endif
  
# 
#--- University of Utah GNU Fortran Compiler Options  ---
# 
    
  ifeq ($(compiler),utah_gnu)
    FC = $(FC_UTAH)
    DEBUG_FLAG = -g -C
    OPTIMIZE_FLAG = -O
    SPDIR ?= /home/yem/stomp/splib
    SPDIR = $(STOMP_SP_DIR)
    DIR = /home/yem/stomp/src
    CPP_FORM = $< | sed '/^\#pragma/d' > $@
    BLASLIB ?= -L /usr/lib -llapack -lblas
    BLASLIB = $(STOMP_BLAS_DIR)
  endif
  
# 
#--- GNU Fortran Compiler Options  ---
# 
    
  ifeq ($(compiler),gnu)
    FC = $(FC_GNU) -w
#  DEBUG_FLAG = -g -gdwarf-3 -gstrict-dwarf -m64 -fomit-frame-pointer
#  OPTIMIZE_FLAG = -O2 -m64 -fomit-frame-pointer
#  DEBUG_FLAG = -g -gdwarf-3 -gstrict-dwarf -m64
    DEBUG_FLAG = -g -m64 -fcheck= -fbounds-check
    OPTIMIZE_FLAG = -O2 -m64 -fcheck= -fbounds-check
#    CPP_FORM = $< | sed '/^\#pragma/d' > $@
    CPP_FORM = $< $@
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = -L /usr/lib -llapack -lblas
  endif
  
# 
#--- F95 (GNU Linux) Fortran Compiler Options  ---
# 
    
  ifeq ($(compiler),f95)
    FC = $(FC_F95)
    DEBUG_FLAG = -g -C
    OPTIMIZE_FLAG = -O
    CPP_FORM = $< | sed '/^\#pragma/d' > $@
#  in system environment variables set $STOMP_BLAS_DIR to 
#    STOMP_BLAS_DIR = -L /usr/lib -llapack -lblas
  endif
  
endif
#
#---  Fortran Compiler Options  ---
#

FFLAGS = $(STANDARD_FLAG)
SUFFIX_FLAG = .e
SP_DIR = $(SPDIR)
ifeq ($(mode),debug)
  FFLAGS = $(DEBUG_FLAG)
  SUFFIX_FLAG = .d
  SP_DIR = $(SPDIR)_debug
endif
ifeq ($(mode),standard)
  FFLAGS = $(STANDARD_FLAG)
  SUFFIX_FLAG = .e
  SP_DIR = $(SPDIR)
endif
ifeq ($(mode),optimize)
  FFLAGS = $(OPTIMIZE_FLAG)
  SUFFIX_FLAG = .x
  SP_DIR = $(SPDIR)
endif
ifeq ($(mode),optimized)
  FFLAGS = $(OPTIMIZE_FLAG)
  SUFFIX_FLAG = .x
  SP_DIR = $(SPDIR)
endif
ifeq ($(solver),petsc)
  FCFLAGS = $(FC_FLAGS) 
  LDFLAGS = $(FC_LINKER_FLAGS)
endif

#
#---  Solver Libraries  ---
#

ifeq ($(solver),sp)
  SP_LIB = $(SP_DIR)/splib.a
  BLAS_LIB = $(BLASLIB)
endif

ifeq ($(solver),lis)
  LIS_LIB = $(LIS_DIR)/lib/liblis.a
endif

ifeq ($(solver),petsc)
  PETSC_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
endif

#
#---  Object File Groups  ---
#

ifeq ($(special),sac)
  sac_objs = balanc.o rdbala.o rdsrel.o rdsrem.o sacout.o sacrem.o saczip.o
endif

ifeq ($(solver),petsc)
  petsc_objs = petsc.o
  petsc_srcs = petsc.f
endif

ifeq ($(solver),lis)
  lis_objs = lis.o
  lis_srcs = lis.f
endif

allo_objs = allo.o $(petsc_objs) $(lis_objs) step.o
global_objs = alloc.o advbb.o band.o banner.o crntnb.o dbgbfa.o dbgbsl.o difmn.o \
            fnitp.o frac.o frac_geom.o fsplnx.o fsplny.o fntblx.o fntbly.o \
            intlz.o jcbbm.o jcbbmb.o jcbp.o jcbz.o lcase.o \
            geomech.o perm_i.o psplib.o rdinpl.o \
            rdchr.o rddpr.o rdijk.o rdgrid.o rdhydr.o rdinac.o rdini.o \
            rdint.o rdmech.o rdouun.o rdrock.o rdrst.o rdsimu.o \
            rdsolu.o rdunit.o refnod.o shdp.o splinx.o spliny.o srchm.o thke.o \
            tmstep.o tortu.o vplwr.o wrmsg.o wrplot.o wrrst.o xyz.o \
            $(sac_objs) $(splib_objs) $(blas_objs)

ifeq ($(chem),ecke)
  eckechem_objs = eckechem.o
  eckechem_srcs = eckechem.f
endif

mode_w_objs = advb.o drcvl.o drcvlb.o eos_w.o \
             jcbwl.o jcbwlb.o mpipp_w.o porsty.o rdlrp.o rombint.o \
             sfin.o sfxl.o sfxlb.o sfxz.o shdpl.o sjcbl.o updtc.o \
             ucode_w.o well_w.o

mode_ws_objs = stomp_w.o $(mode_w_objs) $(eckechem_objs)

mode_gt_objs =  advb.o eos_gt.o evptrns_gt.o flsh_gt.o flux_gt.o \
              frac_gt.o frac_trans_gt.o gasvis.o jcb_gt.o \
              masfrc.o molfrc.o mpipp_gt.o porsty.o rdgrp.o rdlrp.o rdther.o \
              sfin.o sfxg.o sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o \
              shdpl.o sjcb_gt.o ucode_gt.o updtc.o vocgsd.o \
              watdfg.o watgsd.o watgsh.o watgsv.o waticd.o watich.o \
              watick.o watlqd.o watlqh.o watlqk.o watlqv.o well_gt.o

mode_gts_objs = stomp_gt.o $(mode_gt_objs) $(eckechem_objs)

mode_wo_objs = adsorp.o advb.o airgsd.o airgsv.o airdfl.o \
             dffla.o dfflab.o dfflo.o dfflob.o ddflo.o ddflob.o \
             drcvl.o drcvlb.o drcvn.o drcvnb.o \
             fnhgbl.o fnhgbn.o gasvis.o isotherm.o \
             jcbnob.o jcblwb.o jcblob.o jcbwln.o jcboln.o \
             masfrc.o molfrc.o porsty.o rdnrp.o rdvoc.o rdlrp.o\
             rt_napl.o sfin.o sfxl.o sfxlb.o sfxn.o sfxnb.o sfxz.o \
             shdpl.o shdpn.o sjcbl.o sjcbn.o spill.o updtc.o \
             vocdfl.o vocgsd.o vocgsv.o vocsp.o voclqd.o voclqv.o \
             watgsd.o watgsv.o watsp.o watlqd.o watlqv.o well_wo.o

mode_wos_objs = stomp_wo.o $(mode_wo_objs)

mode_woa_objs = adsorp.o advb.o airgsd.o airgsv.o airdfl.o \
             bndfao.o bndfaw.o bndfow.o dff_woa.o \
             drcvg.o drcvgb.o drcvl.o drcvlb.o drcvn.o drcvnb.o \
             fnhgbg.o fnhgbl.o fnhgbn.o gasvis.o jcb_woa.o \
             masfrc.o molfrc.o porsty.o \
             rdvoc.o rdgrp.o rdnrp.o rdlrp.o rt_napl.o \
             sfin.o sfxg.o sfxgb.o sfxl.o sfxlb.o sfxn.o sfxnb.o \
             sfxz.o shdpg.o shdpl.o shdpn.o sjcbg.o sjcbl.o sjcbn.o \
             updtc.o vocdfg.o vocdfl.o vocgsd.o vocgsv.o vocsp.o \
             voclqd.o voclqv.o \
             watdfg.o watgsd.o watgsv.o watsp.o watlqd.o watlqv.o

mode_woas_objs = stomp_woa.o $(mode_woa_objs)

mode_co2_objs = advb.o dff_co2e.o drcv_co2e.o eos_co2e.o flsh_co2e.o \
              jcb_co2e.o mpipp_co2.o porsty.o rdgrp.o rdlrp.o sfin.o sfxb_co2.o \
              sfxg.o sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o updtc.o well_co2.o

mode_co2s_objs = stomp_co2.o $(mode_co2_objs) $(eckechem_objs)

mode_co2e_objs = advb.o dff_co2e.o drcv_co2e.o eos_co2e.o flsh_co2e.o \
              jcb_co2e.o porsty.o rdgrp.o rdlrp.o rdther.o sfin.o \
              sfxb_co2.o sfxg.o sfxl.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_co2e.o updtc.o \
              well_co2e.o

mode_co2es_objs = stomp_co2e.o $(mode_co2e_objs) $(eckechem_objs)

mode_seq_objs = advb.o dff_seq.o drcv_seq.o eos_seq.o flsh_seq.o \
              jcb_seq.o porsty.o rdgrp.o rdlrp.o rdnrp.o rdther.o \
              sfin.o sfx_seq.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_seq.o updtc.o \
              well_seq.o

mode_seqs_objs = stomp_seq.o $(mode_seq_objs) $(eckechem_objs)

mode_comp_objs = advb.o dff_comp.o drcv_comp.o eos_comp.o hgbg_comp.o \
              hgbl_comp.o jcb_comp.o porsty.o rdgrp.o rdlrp.o rdther.o \
              sfin.o sfxb_co2.o sfxg.o sfxl.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_comp.o updtc.o 

mode_comps_objs = stomp_comp.o $(mode_comp_objs) $(eckechem_objs)

mode_hyd_objs = advb.o dff_hyd.o drcv_hyd.o eos_hyd.o flsh_hyd.o \
              hgb_hyd.o jcb_hyd.o porsty.o \
              rdgrp.o rdlrp.o rdnrp.o rdther.o sfin.o sfxg.o \
              sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_hyd.o updtc.o

mode_hyds_objs = stomp_hyd.o $(mode_hyd_objs)

mode_hyd_ke_objs = advb.o dff_hyd.o drcv_hyd.o eos_hyd_ke.o flsh_hyd_ke.o \
              hgb_hyd.o jcb_hyd_ke.o porsty.o \
              rdgrp.o rdlrp.o rdnrp.o rdther.o sfin.o sfxg.o \
              sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_hyd.o updtc.o

mode_hyd_kes_objs = stomp_hyd_ke.o $(mode_hyd_ke_objs)

mode_hydt_ke_objs = advb.o dffg_hydt_ke.o dffl_hydt_ke.o dffn_hydt_ke.o \
              drcv_hydt_ke.o eos_hydt_ke.o flsh_hydt_ke.o \
              jcba_hydt_ke.o jcbn_hydt_ke.o jcbo_hydt_ke.o \
              jcbs_hydt_ke.o jcbt_hydt_ke.o jcbw_hydt_ke.o mpipp_hydt_ke.o \
              porsty.o rdgrp.o rdlrp.o rdnrp.o rdther.o sfin.o sfxg.o \
              sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_hydt_ke.o updtc.o

mode_hydt_kes_objs = stomp_hydt_ke.o $(mode_hydt_ke_objs)

mode_hyd_o2_ch4_n2_objs = advb.o dffg_hydt_ke.o dffl_hydt_ke.o dffn_hydt_ke.o \
              drcv_hydt_ke.o eos_hyd_o2_ch4_n2.o flsh_hyd_o2_ch4_n2.o \
              jcba_hydt_ke.o jcbn_hydt_ke.o jcbo_hydt_ke.o \
              jcbs_hydt_ke.o jcbt_hydt_ke.o jcbw_hydt_ke.o mpipp_hyd_o2_ch4_n2.o\
              porsty.o rdgrp.o rdlrp.o rdnrp.o rdther.o sfin.o sfxg.o \
              sfxgb.o sfxl.o sfxlb.o sfxz.o shdpg.o shdpl.o \
              sjcbg.o sjcbl.o thd_hydt_ke.o updtc.o

mode_hyd_o2_ch4_n2s_objs = stomp_hyd_o2_ch4_n2.o $(mode_hyd_o2_ch4_n2_objs)

mode_eor_objs = advb.o dffg_eor.o dffl_eor.o dffn_eor.o \
              drcv_eor.o dual_por_eor.o eos_eor.o flsh_eor.o \
              frac_eor.o frac_trans_eor.o frac_rdrock.o \
              jcbc_eor.o jcbs_eor.o jcbt_eor.o jcbw_eor.o mpipp_eor.o \
              porsty.o rdfrac.o rdther.o sfin.o sfx_eor.o sfxz.o \
              shdpg.o shdpl.o shdpn.o \
              sjcb_eor.o thd_eor.o updtc.o well_eor.o

mode_eors_objs = stomp_eor.o $(mode_eor_objs) $(eckechem_objs)

#
#---  Source File Groups  ---
#

ifeq ($(special),sac)
  sac_srcs = balanc.f rdbala.f rdsrel.f rdsrem.f sacout.f sacrem.f saczip.f
endif

allo_srcs = allo.f
global_srcs = $(petsc_srcs) $(lis_srcs) alloc.f step.f \
            advbb.f band.f banner.f \
            crntnb.f dbgbfa.f dbgbsl.f difmn.f fnitp.f \
            fsplnx.f fsplny.f fntblx.f fntbly.f \
            frac.f frac_geom.f intlz.f jcbbm.f jcbbmb.f jcbp.f \
            jcbz.f lcase.f geomech.f perm_i.f psplib.f \
            rdinpl.f rdijk.f \
            rdchr.f rddpr.f rdgrid.f rdhydr.f \
            rdinac.f rdini.f rdint.f \
            rdmech.f rdouun.f rdrock.f rdrst.f \
            rdsimu.f rdsolu.f rdunit.f refnod.f \
            shdp.f splinx.f spliny.f srchm.f thke.f \
            tmstep.f tortu.f vplwr.f wrmsg.f \
            wrplot.f wrrst.f xyz.f \
            $(sac_srcs) $(splib_srcs) $(blas_srcs)

mode_w_srcs = advb.f drcvl.f drcvlb.f \
             eos_w.f jcbwl.f jcbwlb.f mpipp_w.f \
             porsty.f rdlrp.f rombint.f sfin.f \
             sfxl.f sfxlb.f sfxz.f shdpl.f \
             sjcbl.f updtc.f ucode_w.f well_w.f $(eckechem_srcs)

mode_gt_srcs = advb.f eos_gt.f evptrns_gt.f flsh_gt.f \
              flux_gt.f frac_gt.f frac_trans_gt.f gasvis.f jcb_gt.f \
              masfrc.f molfrc.f mpipp_gt.f porsty.f rdgrp.f rdlrp.f rdther.f \
              sfin.f sfxg.f sfxgb.f sfxl.f sfxlb.f sfxz.f \
              shdpg.f shdpl.f sjcb_gt.f \
              updtc.f ucode_gt.f vocgsd.f \
              watdfg.f watgsd.f watgsh.f watgsv.f \
              waticd.f watich.f watick.f watlqd.f \
              watlqh.f watlqk.f watlqv.f well_gt.f \
              $(eckechem_srcs)

mode_wo_srcs = adsorp.f advb.f airgsd.f airgsv.f \
             airdfl.f ddflo.f ddflob.f \
             dfflo.f dfflob.f dffla.f dfflab.f \
             drcvl.f drcvlb.f drcvn.f isotherm.f \
             drcvnb.f fnhgbl.f fnhgbn.f \
             gasvis.f jcbnob.f jcblwb.f \
             jcblob.f jcbwln.f jcboln.f \
             masfrc.f molfrc.f porsty.f rdnrp.f \
             rdvoc.f rdlrp.f rt_napl.f sfin.f sfxl.f \
             sfxlb.f sfxn.f sfxnb.f sfxz.f \
             shdpl.f shdpn.f sjcbl.f sjcbn.f spill.f \
             updtc.f vocdfl.f vocgsd.f vocgsv.f \
             vocsp.f voclqd.f voclqv.f watgsd.f \
             watgsv.f watsp.f watlqd.f watlqv.f \
             well_wo.f

mode_woa_srcs = adsorp.f advb.f airgsd.f airgsv.f airdfl.f \
             bndfao.f bndfaw.f bndfow.f dff_woa.f \
             drcvg.f drcvgb.f drcvl.f drcvlb.f \
             drcvn.f drcvnb.f fnhgbg.f fnhgbl.f \
             fnhgbn.f gasvis.f jcb_woa.f \
             masfrc.f molfrc.f porsty.f rdvoc.f rdgrp.f \
             rdnrp.f rdlrp.f rt_napl.f sfin.f sfxg.f \
             sfxgb.f sfxl.f sfxlb.f sfxn.f \
             sfxnb.f sfxz.f shdpg.f \
             shdpl.f shdpn.f sjcbg.f sjcbl.f \
             sjcbn.f  updtc.f vocdfg.f vocdfl.f \
             vocgsd.f vocgsv.f vocsp.f voclqd.f \
             voclqv.f watdfg.f watgsd.f \
             watgsv.f watsp.f watlqd.f watlqv.f

mode_co2_srcs = advb.f dff_co2e.f drcv_co2e.f eos_co2e.f flsh_co2e.f \
              jcb_co2e.f mpipp_co2.f porsty.f rdgrp.f rdlrp.f sfin.f sfxb_co2.f \
              sfxg.f sfxgb.f sfxl.f sfxlb.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f updtc.f well_co2.f \
              $(eckechem_srcs)

mode_co2e_srcs = advb.f dff_co2e.f drcv_co2e.f eos_co2e.f flsh_co2e.f \
              jcb_co2e.f porsty.f rdgrp.f rdlrp.f rdther.f sfin.f \
              sfxb_co2.f sfxg.f sfxl.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_co2e.f updtc.f \
              well_co2e.f $(eckechem_srcs)

mode_seq_srcs = advb.f dff_seq.f drcv_seq.f eos_seq.f flsh_seq.f \
              jcb_seq.f porsty.f rdgrp.f rdlrp.f rdnrp.f rdther.f sfin.f \
              sfx_seq.f sfxg.f sfxl.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_seq.f updtc.f \
              well_seq.f $(eckechem_srcs)

mode_comp_srcs = advb.f dff_comp.f drcv_comp.f eos_comp.f hgbg_comp.f \
              hgbl_comp.f jcb_comp.f jcbgcb_comp.f jcbgwb_comp.f \
              jcblcb_comp.f jcblwb_comp.f jcbsb_comp.f jcbtb_comp.f \
              porsty.f rdgrp.f rdlrp.f rdther.f sfin.f sfxb_co2.f \
              sfxg.f sfxl.f sfxz.f shdpg.f shdpl.f sjcbg.f sjcbl.f \
              thd_comp.f updtc.f $(eckechem_srcs)

mode_hyd_srcs = advb.f dff_hyd.f drcv_hyd.f eos_hyd.f flsh_hyd.f \
              hgb_hyd.f jcb_hyd.f porsty.f \
              rdgrp.f rdlrp.f rdnrp.f rdther.f sfin.f sfxg.f \
              sfxgb.f sfxl.f sfxlb.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_hyd.f updtc.f

mode_hyd_ke_srcs = advb.f dff_hyd.f drcv_hyd.f eos_hyd_ke.f flsh_hyd_ke.f \
              hgb_hyd.f jcb_hyd_ke.f porsty.f \
              rdgrp.f rdlrp.f rdnrp.f rdther.f sfin.f sfxg.f \
              sfxgb.f sfxl.f sfxlb.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_hyd.f updtc.f

mode_hydt_ke_srcs = advb.f dffg_hydt_ke.f dffl_hydt_ke.f dffn_hydt_ke.f \
              drcv_hydt_ke.f eos_hydt_ke.f flsh_hydt_ke.f \
              jcba_hydt_ke.f jcbn_hydt_ke.f jcbo_hydt_ke.f \
              jcbs_hydt_ke.f jcbt_hydt_ke.f jcbw_hydt_ke.f \
              mpipp_hydt_ke.f porsty.f \
              rdgrp.f rdlrp.f rdnrp.f rdther.f sfin.f sfxg.f \
              sfxgb.f sfxl.f sfxlb.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_hydt_ke.f updtc.f

mode_hyd_o2_ch4_n2_srcs = advb.f dffg_hydt_ke.f dffl_hydt_ke.f dffn_hydt_ke.f \
              drcv_hydt_ke.f eos_hyd_o2_ch4_n2.f flsh_hyd_o2_ch4_n2.f \
              jcba_hydt_ke.f jcbn_hydt_ke.f jcbo_hydt_ke.f \
              jcbs_hydt_ke.f jcbt_hydt_ke.f jcbw_hydt_ke.f \
              mpipp_hyd_o2_ch4_n2.f porsty.f \
              rdgrp.f rdlrp.f rdnrp.f rdther.f sfin.f sfxg.f \
              sfxgb.f sfxl.f sfxlb.f sfxz.f shdpg.f shdpl.f \
              sjcbg.f sjcbl.f thd_hydt_ke.f updtc.f

mode_eor_srcs = advb.f dffg_eor.f dffl_eor.f dffn_eor.f \
              drcv_eor.f dual_por_eor.f eos_eor.f flsh_eor.f \
              frac_eor.f frac_trans_eor.f frac_rdrock.f \
              jcbc_eor.f jcbs_eor.f jcbt_eor.f jcbw_eor.f mpipp_eor.f \
              porsty.f rdfrac.f rdther.f sfin.f sfx_eor.f \
              sfxz.f shdpg.f shdpl.f shdpn.f sjcb_eor.f thd_eor.f \
              updtc.f well_eor.f $(eckechem_srcs)

#
#---  Targets  ---
#

prb1:
	echo "Copying STOMP Short Course Problem #1."
	cp $(DIR)/../short_course/prb1/input .
	echo "Creating parameters file for banded solver."
	$(DIR)/../bin/step_bd.e
	echo "Making static-memory, banded-solver version of STOMP1."
	make stomp1 mode=optimize solver=bd
	rm *.o parameters
	echo "Executing static-memory, banded-solver version of STOMP1 for Short Course Problem #1."
	stomp1_bd.x > screen
	echo "Testing static-memory, banded-solver version of STOMP1 for Short Course Problem #1."
	rm stomp1_bd.x screen
	echo "Creating parameters file for SPLIB solver."
	$(DIR)/../bin/step_bd.e
	echo "Making static-memory, SPLIB-solver version of STOMP1."
	make stomp1 mode=optimize solver=sp
	rm *.o parameters
	echo "Executing static-memory, SPLIB-solver version of STOMP1 for Short Course Problem #1."
	stomp1_sp.x > screen
	echo "Testing static-memory, SPLIB-solver version of STOMP1 for Short Course Problem #1."
	
print:
	echo "FFLAGS = $(FFLAGS)"

ifeq ($(mode),source)
ifeq ($(solver),sp)
step: allo.f step.f
	cat allo.f step.f > step_sp.f
	echo "Source step_sp has been updated."
endif
ifeq ($(solver),lis)
step: allo.f step.f
	cat allo.f step.f > step_lis.f
	echo "Source step_lis has been updated."
endif
ifeq ($(solver),bd)
step: allo.f step.f
	cat allo.f step.f > step_bd.f
	echo "Source step_bd has been updated."
endif
else
ifeq ($(solver),sp)
step: allo.o step.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o step_sp$(SUFFIX_FLAG) \
      allo.o step.o
	echo "Executable step_sp has been updated."
endif
ifeq ($(solver),lis)
step: allo.o step.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o step_lis$(SUFFIX_FLAG) \
      allo.o step.o
	echo "Executable step_lis has been updated."
endif
ifeq ($(solver),bd)
step: allo.o step.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o step_bd$(SUFFIX_FLAG) \
      allo.o step.o
	echo "Executable step_bd has been updated."
endif
endif
#endif
#endif

ifeq ($(solver),sp)
ifeq ($(mode),source)
splib: $(splib_srcs)
	zip splib.zip $(splib_srcs)
	echo "Source splib has been updated."
endif
endif

ifeq ($(solver),sp)
ifeq ($(mode),source)
blas: $(blas_srcs)
	zip blas.zip $(blas_srcs)
	echo "Source blas has been updated."
endif
endif

ifeq ($(mode),source)
stomp_w stomp-w: $(global_srcs) $(allo_srcs) $(mode_w_srcs) stomp_w.f
	cp stomp_w.f stomp-w_$(solver).f
	cat $(mode_w_srcs) $(global_srcs) >> stomp-w_$(solver).f
	echo "Source stomp-w_$(solver) has been updated."
else
stomp_w stomp-w:  $(allo_objs) $(global_objs) $(mode_ws_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-w_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_ws_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-w_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_gt stomp-gt: $(global_srcs) $(allo_srcs) $(mode_gt_srcs) stomp_gt.f
	cp stomp_gt.f stomp-gt_$(solver).f
	cat $(mode_gt_srcs) $(global_srcs) >> stomp-gt_$(solver).f
	echo "Source stomp-gt_$(solver) has been updated."
else
stomp_gt stomp-gt:  $(allo_objs) $(global_objs) $(mode_gts_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-gt_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_gts_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-gt_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_wo stomp-wo: $(global_srcs) $(allo_srcs) $(mode_wo_srcs) stomp_wo.f
	cp stomp_wo.f stomp-wo_$(solver).f
	cat $(mode_wo_srcs) $(global_srcs) >> stomp-wo_$(solver).f
	echo "Source stomp-wo_$(solver) has been updated."
else
stomp_wo stomp-wo:  $(allo_objs) $(global_objs) $(mode_wos_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-wo_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_wos_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-wo_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_woa stomp-woa stomp-voc: $(global_srcs) $(allo_srcs) $(mode_woa_srcs) stomp_woa.f
	cp stomp_woa.f stomp-woa_$(solver).f
	cat $(mode_woa_srcs) $(global_srcs) >> stomp-woa_$(solver).f
	echo "Source stomp-woa_$(solver) has been updated."
else
stomp_woa stomp-woa stomp-voc:  $(allo_objs) $(global_objs) $(mode_woas_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-woa_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_woas_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-woa_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_hyd stomp-hyd: $(global_srcs) $(allo_srcs) $(mode_hyd_srcs) stomp_hyd.f
	cp stomp_hyd.f stomp-hyd_$(solver).f
	cat $(mode_hyd_srcs) $(global_srcs) >> stomp-hyd_$(solver).f
	echo "Source stomp-hyd_$(solver) has been updated."
else
stomp_hyd stomp-hyd:  $(allo_objs) $(global_objs) $(mode_hyds_objs) 
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-hyd_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_hyds_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-hyd_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_co2 stomp-co2: $(global_srcs) $(allo_srcs) $(mode_co2_srcs) stomp_co2.f
	cp stomp_co2.f stomp-co2_$(solver).f
	cat $(mode_co2_srcs) $(global_srcs) >> stomp-co2_$(solver).f
	echo "Source stomp-co2_$(solver) has been updated."
else
stomp_co2 stomp-co2:  $(allo_objs) $(global_objs) $(mode_co2s_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-co2_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_co2s_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-co2_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_co2e stomp-co2e: $(global_srcs) $(allo_srcs) $(mode_co2e_srcs) stomp_co2e.f
	cp stomp_co2e.f stomp-co2e_$(solver).f
	cat $(mode_co2e_srcs) $(global_srcs) >> stomp-co2e_$(solver).f
	echo "Source stomp-co2e_$(solver) has been updated."
else
stomp_co2e stomp-co2e:  $(allo_objs) $(global_objs) $(mode_co2es_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-co2e_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_co2es_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-co2e_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_seq stomp-seq: $(global_srcs) $(allo_srcs) $(mode_seq_srcs) stomp_seq.f
	cp stomp_seq.f stomp-seq_$(solver).f
	cat $(mode_seq_srcs) $(global_srcs) >> stomp-seq_$(solver).f
	echo "Source stomp-seq_$(solver) has been updated."
else
stomp_seq stomp-seq:  $(allo_objs) $(global_objs) $(mode_seqs_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-seq_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_seqs_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-seq_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_hyd_ke stomp-hyd-ke: $(global_srcs) $(allo_srcs) $(mode_hyd_ke_srcs) stomp_hyd_ke.f
	cp stomp_hyd_ke.f stomp-hyd-ke_$(solver).f
	cat $(mode_hyd_ke_srcs) $(global_srcs) >> stomp-hyd-ke_$(solver).f
	echo "Source stomp-hyd-ke_$(solver) has been updated."
else
stomp_hyd_ke stomp-hyd-ke:  $(allo_objs) $(global_objs) $(mode_hyd_kes_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-hyd-ke_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_hyd_kes_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING) 
	echo "Executable stomp-hyd-ke_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_hydt_ke stomp-hydt-ke: $(global_srcs) $(allo_srcs) $(mode_hydt_ke_srcs) stomp_hydt_ke.f
	cp stomp_hydt_ke.f stomp-hydt-ke_$(solver).f
	cat $(mode_hydt_ke_srcs) $(global_srcs) >> stomp-hydt-ke_$(solver).f
	echo "Source stomp-hydt-ke_$(solver) has been updated."
else
stomp_hydt_ke stomp-hydt-ke:  $(allo_objs) $(global_objs) $(mode_hydt_kes_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-hydt-ke_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_hydt_kes_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-hydt-ke_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_hyd_o2_ch4_n2 stomp-hyd-o2-ch4-n2: $(global_srcs) $(allo_srcs) $(mode_hyd_o2_ch4_n2_srcs) stomp_hyd_o2_ch4_n2.f
	cp stomp_hyd_o2_ch4_n2.f stomp-hyd-o2-ch4-n2_$(solver).f
	cat $(mode_hyd_o2_ch4_n2_srcs) $(global_srcs) >> stomp-hyd-o2-ch4-n2_$(solver).f
	echo "Source stomp-hyd-o2-ch4-n2_$(solver) has been updated."
else
stomp_hyd_o2_ch4_n2 stomp-hyd-o2-ch4-n2:  $(allo_objs) $(global_objs) $(mode_hyd_o2_ch4_n2s_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-hyd-o2-ch4-n2_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_hyd_o2_ch4_n2s_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-hyd-o2-ch4-n2_$(solver) has been updated."
endif

ifeq ($(mode),source)
stomp_eor stomp-eor: $(global_srcs) $(allo_srcs) $(mode_eor_srcs) stomp_eor.f
	cp stomp_eor.f stomp-eor_$(solver).f
	cat $(mode_eor_srcs) $(global_srcs) >> stomp-eor_$(solver).f
	echo "Source stomp-eor_$(solver) has been updated."
else
stomp_eor stomp-eor:  $(allo_objs) $(global_objs) $(mode_eors_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stomp-eor_$(solver)$(SUFFIX_FLAG) $(allo_objs) \
      $(global_objs) $(mode_eors_objs) $(BLAS_LIB) $(LIS_LIB) $(PETSC_LIB) $(SP_LIB) $(GNU_LINKING)
	echo "Executable stomp-eor_$(solver) has been updated."
endif


ifeq ($(mode),debug)
.SECONDARY: $(sac_srcs) $(global_srcs) $(allo_srcs) $(mode_w_srcs) \
  $(mode_gt_srcs) $(mode_wo_srcs) $(mode_woa_srcs) \
  $(mode_co2_srcs) \
  $(mode_co2e_srcs) \
  $(mode_seq_srcs) $(mode_hyd_srcs) \
  $(mode_hyd_ke_srcs) $(mode_hydt_ke_srcs) $(mode_hyd_o2_ch4_n2_srcs) \
  $(mode_comp_srcs) $(mode_eor_srcs) $(eckechem_srcs) \
  stomp_w.f stomp_gt.f stomp_wo.f stomp_woa.f \
  stomp_co2.f stomp_co2e.f stomp_seq.f stomp_comp.f \
  stomp_hyd.f stomp_hyd_ke.f stomp_hydt_ke.f stomp_hyd_o2_ch4_n2.f \
  stomp_eor.f step.f eos_comp.f *.f petsc.f lis.f
endif


#
#---  Dependencies  ---
#

step.o:			 allo.o
$(global_objs):    allo.o
$(mode_ws_objs):    allo.o
$(mode_gts_objs):    allo.o
$(mode_eors_objs):    allo.o
$(mode_wos_objs):    allo.o
$(mode_woas_objs):    allo.o
$(mode_co2s_objs):   allo.o
$(mode_co2es_objs):   allo.o
$(mode_seq_objs):   allo.o
$(mode_comps_objs):   allo.o
$(mode_hyd_kes_objs):  allo.o
$(mode_hydt_kes_objs):  allo.o
$(mode_hyd_o2_ch4_n2s_objs):  allo.o
$(mode_hyds_objs):   allo.o

#
#---  Pattern Rules  ---
#

%.o:  %.F

%.f:  $(DIR)/%.F
	/usr/bin/cpp $(CPP_FLAGS) $(CPP_FORM)

%.f:  $(DIR)/splib/%.F
	/usr/bin/cpp $(CPP_FLAGS) $(CPP_FORM)

%.o:  %.f
	$(FC) $(FFLAGS) $(FCFLAGS) $(SRCH_PATH) -c $<

#
#---  End of makefile  ---
#

