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
solver = nada

#
#---  Environment Variables  ---
#
#  The user must set the following environment variables, otherwise Mark's default paths will be used
#
# $STOMP_COMPILER [gnu]
# $STOMP_SRC_DIR [path to stomp source code]
# $STOMP_LIS_DIR [path to solver libraries]
# if using the petsc solver $PETSC_DIR and $PETSC_ARCH must be set
#
STOMP_COMPILER ?= mpi
compiler ?= $(STOMP_COMPILER)
STOMP_SRC_DIR ?= /Users/d3c002/stomp/src_x
DIR = $(STOMP_SRC_DIR)
PETSC_DIR ?= /Users/d3c002/petsc
PETSC_ARCH ?= arch-darwin-c-opt

#
#---  Definitions  ---
#

ifeq ($(solver),PETSc)
  solver = petsc
endif
ifeq ($(solver),Lis)
  solver = lis
endif
ifndef solver
  $(error Error: Unrecognized solver: Please use [solver=lis or solver=petsc])
endif
ifneq ($(solver),lis)
  ifneq ($(solver),petsc)
    ifneq ($(solver),nada)
      $(error Error: Unrecognized solver: Please use [solver=lis or solver=petsc])
    endif
  endif
endif
ifneq ($(decay),radioactive)
  ifneq ($(decay),reactive)
    $(error Error: Unrecognized decay: Please use [decay=radioactive or decay=reactive])
  endif
endif
ifeq ($(chem),NAPL)
  chem = napl
endif
ifeq ($(chem),napl)
  chem = napl
endif
#
#---  Set Petsc flags ---
#

ifeq ($(solver),petsc)
  ifeq ($(mode),debug )
    PETSC_ARCH = arch-darwin-c-debug
  endif
  include $(PETSC_DIR)/lib/petsc/conf/variables
  include $(PETSC_DIR)/lib/petsc/conf/rules
  ifeq ($(compiler),mpi)
    PETSC_LINKING = -L /usr/local/lib -lgfortran
  endif
  PETSC_CPP_FLAGS = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
endif

#
#---  Set lis directory  ---
#
ifeq ($(solver),lis)
  ifeq ($(compiler),mpi)
    STOMP_LIS_DIR ?= /Users/d3c002/lis-2.0.24/lis_mpi
    LIS_DIR = $(STOMP_LIS_DIR)
    LIS_LINKING = -L /usr/local/lib -lgfortran
  endif
  LIS_CPP_FLAGS = $(FCPPFLAGS) -I$(LIS_DIR)/include 
endif

#
#---  C Preprocessor  ---
#

CPP_FLAGS = -w -P -traditional -D$(solver) 
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
ifeq ($(compiler_opt),omp)
  CPP_FLAGS += $(OMP_CPP_FLAGS)
endif
CPP_FORM = $< $@

#
#---  Source Paths  ---
#

SRCH_PATH =  -I . -I $(DIR)
ifeq ($(solver),petsc)
  SRCH_PATH += -I $(PETSC_DIR)/$(PETSC_ARCH)/include
endif

#
#---  Fortran Compilers  ---
#
  FC_MPI = mpif90
# 
#--- OpenMPI Fortran Compiler Options  ---
# 
  FC = $(FC_MPI)
  DEBUG_FLAG = -g
  OPTIMIZE_FLAG = -O2
  CPP_FORM = $< $@

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

#
#---  Solver Libraries  ---
#

ifeq ($(solver),lis)
  LIS_LIB = $(LIS_DIR)/lib/liblis.a
endif

ifeq ($(solver),petsc)
  PETSC_LIB = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
endif

#
#---  Object File Groups  ---
#

ifeq ($(solver),petsc)
  petsc_objs = petsc_x.o
  petsc_srcs = petsc_x.f
endif

ifeq ($(solver),lis)
  lis_objs = lis_x.o
  lis_srcs = lis_x.f
endif

allo_hydt_ke_objs = allo_hydt_ke_x.o $(petsc_objs) $(lis_objs)

allo_hyd_o2_ch4_n2_objs = allo_hyd_o2_ch4_n2_x.o $(petsc_objs) $(lis_objs)

allo_co2_objs = allo_co2_x.o $(petsc_objs) $(lis_objs)

allo_eor_objs = allo_eor_x.o $(petsc_objs) $(lis_objs)

allo_gt_objs = allo_gt_x.o $(petsc_objs) $(lis_objs)

allo_w_objs = allo_w_x.o $(petsc_objs) $(lis_objs)

mode_hydt_ke_objs = stomp_hydt_ke_x.o alloc_hydt_ke_x.o \
	difmn_x.o eos_hydt_ke_x.o flsh_hydt_ke_x.o flux_hydt_ke_x.o \
	fluxb_hydt_ke_x.o geomech_x.o intlz_hydt_ke_x.o jcb_hydt_ke_x.o jcbb_hydt_ke_x.o \
	output_hydt_ke_x.o restart_hydt_ke_x.o tmstep_hydt_ke_x.o wrmsg_x.o

mode_hyd_o2_ch4_n2_objs = stomp_hyd_o2_ch4_n2_x.o alloc_hyd_o2_ch4_n2_x.o \
	difmn_x.o eos_hyd_o2_ch4_n2_x.o flsh_hyd_o2_ch4_n2_x.o flux_hyd_o2_ch4_n2_x.o \
	fluxb_hyd_o2_ch4_n2_x.o geomech_x.o intlz_hyd_o2_ch4_n2_x.o jcb_hyd_o2_ch4_n2_x.o jcbb_hyd_o2_ch4_n2_x.o \
	output_hyd_o2_ch4_n2_x.o restart_hyd_o2_ch4_n2_x.o tmstep_hyd_o2_ch4_n2_x.o wrmsg_x.o

mode_co2_objs = stomp_co2_x.o alloc_co2_x.o difmn_x.o eckechem_x.o \
	eos_co2e_x.o flsh_co2e_x.o flux_co2_x.o fluxb_co2_x.o geomech_x.o intlz_co2_x.o \
	jcb_co2_x.o jcbb_co2_x.o output_co2_x.o restart_co2_x.o tmstep_co2_x.o \
	tport_co2_x.o well_co2_x.o wrmsg_x.o

mode_eor_objs = stomp_eor_x.o alloc_eor_x.o difmn_x.o eckechem_x.o \
	eos_eor_x.o flsh_eor_x.o flux_eor_x.o fluxb_eor_x.o geomech_x.o intlz_eor_x.o \
	jcb_eor_x.o jcbb_eor_x.o output_eor_x.o restart_eor_x.o tmstep_eor_x.o \
	tport_eor_x.o well_eor_x.o wrmsg_x.o

mode_gt_objs = stomp_gt_x.o alloc_gt_x.o difmn_x.o eckechem_x.o \
	eos_gt_x.o flsh_gt_x.o flux_gt_x.o fluxb_gt_x.o frac_gt_x.o \
	frac_tport_gt_x.o geomech_x.o intlz_gt_x.o \
	jcb_gt_x.o jcbb_gt_x.o output_gt_x.o restart_gt_x.o tmstep_gt_x.o \
	tport_gt_x.o wrmsg_x.o

mode_w_objs = stomp_w_x.o alloc_w_x.o difmn_x.o eckechem_x.o \
	eos_w_x.o flux_w_x.o fluxb_w_x.o geomech_x.o intlz_w_x.o \
	jcb_w_x.o jcbb_w_x.o output_w_x.o restart_w_x.o tmstep_w_x.o \
	tport_w_x.o well_w_x.o wrmsg_x.o

#
#---  Source File Groups  ---
#

allo_hydt_ke_srcs = allo_hydt_ke_x.f $(petsc_srcs) $(lis_srcs) 

allo_hyd_o2_ch4_n2_srcs = allo_hyd_o2_ch4_n2_x.f $(petsc_srcs) $(lis_srcs) 

allo_co2_srcs = allo_co2_x.f $(petsc_srcs) $(lis_srcs) 

allo_eor_srcs = allo_eor_x.f $(petsc_srcs) $(lis_srcs) 

allo_gt_srcs = allo_gt_x.f $(petsc_srcs) $(lis_srcs) 

allo_w_srcs = allo_w_x.f $(petsc_srcs) $(lis_srcs) 

mode_hydt_ke_srcs = stomp_hydt_ke_x.f alloc_hydt_ke_x.f \
	difmn_x.f eos_hydt_ke_x.f flsh_hydt_ke_x.f flux_hydt_ke_x.f \
	fluxb_hydt_ke_x.f geomech_x.f intlz_hydt_ke_x.f jcb_hydt_ke_x.f jcbb_hydt_ke_x.f \
	output_hydt_ke_x.f restart_hydt_ke_x.f tmstep_hydt_ke_x.f wrmsg_x.f

mode_hyd_o2_ch4_n2_srcs = stomp_hyd_o2_ch4_n2_x.f alloc_hyd_o2_ch4_n2_x.f \
	difmn_x.f eos_hyd_o2_ch4_n2_x.f flsh_hyd_o2_ch4_n2_x.f flux_hyd_o2_ch4_n2_x.f \
	fluxb_hyd_o2_ch4_n2_x.f geomech_x.f intlz_hyd_o2_ch4_n2_x.f jcb_hyd_o2_ch4_n2_x.f jcbb_hyd_o2_ch4_n2_x.f \
	output_hyd_o2_ch4_n2_x.f restart_hyd_o2_ch4_n2_x.f tmstep_hyd_o2_ch4_n2_x.f wrmsg_x.f

mode_co2_srcs = stomp_co2_x.f alloc_co2_x.f difmn_x.f eckechem_x.f\
	eos_co2e_x.f flsh_co2e_x.f flux_co2_x.f fluxb_co2_x.f geomech_x.f intlz_co2_x.f \
	jcb_co2_x.f jcbb_co2_x.f output_co2_x.f restart_co2_x.f tmstep_co2_x.f \
	tport_co2_x.f well_co2_x.f wrmsg_x.f

mode_eor_srcs = stomp_eor_x.f alloc_eor_x.f difmn_x.f eckechem_x.f\
	eos_eor_x.f flsh_eor_x.f flux_eor_x.f fluxb_eor_x.f geomech_x.f intlz_eor_x.f \
	jcb_eor_x.f jcbb_eor_x.f output_eor_x.f restart_eor_x.f tmstep_eor_x.f \
	tport_eor_x.f well_eor_x.f wrmsg_x.f

mode_gt_srcs = stomp_gt_x.f alloc_gt_x.f difmn_x.f eckechem_x.f\
	eos_gt_x.f flsh_gt_x.f flux_gt_x.f fluxb_gt_x.f frac_gt_x.f \
	frac_tport_gt_x.f geomech_x.f intlz_gt_x.f \
	jcb_gt_x.f jcbb_gt_x.f output_gt_x.f restart_gt_x.f tmstep_gt_x.f \
	tport_gt_x.f well_gt_x.f wrmsg_x.f

mode_w_srcs = stomp_w_x.f alloc_w_x.f difmn_x.f eckechem_x.f\
	eos_w_x.f flux_w_x.f fluxb_w_x.f geomech_x.f intlz_w_x.f \
	jcb_w_x.f jcbb_w_x.f output_w_x.f restart_w_x.f tmstep_w_x.f \
	tport_w_x.f well_w_x.f wrmsg_x.f

#
#---  Targets  ---
#

ifeq ($(mode),source)
stompx_hydt_ke stompx-hydt-ke: $(allo_hydt_ke_srcs) $(mode_hydt_ke_srcs) stomp_hydt_ke_x.f
	cp stomp_hydt_ke_x.f stompx-hydt-ke_$(solver).f
	cat $(mode_hydt_ke_srcs) >> stompx-hydt-ke_$(solver).f
	echo "Source stompx-hydt-ke-mpi_$(solver) has been updated."
stompx_hyd_o2_ch4_n2 stompx-hyd-o2-ch4-n2: $(allo_hyd_o2_ch4_n2_srcs) $(mode_hyd_o2_ch4_n2_srcs) stomp_hyd_o2_ch4_n2_x.f
	cp stomp_hyd_o2_ch4_n2_x.f stompx-hyd-o2-ch4_n2_$(solver).f
	cat $(mode_hyd_o2_ch4_n2_srcs) >> stompx-hyd_o2_ch4_n2_$(solver).f
	echo "Source stompx-hyd-o2-ch4-n2-mpi_$(solver) has been updated."
stompx_co2 stompx-co2: $(allo_co2_srcs) $(mode_co2_srcs) stomp_co2_x.f
	cp stomp_co2_x.f stompx-co2_$(solver).f
	cat $(mode_co2_srcs) >> stompx-co2_$(solver).f
	echo "Source stompx-co2_$(solver) has been updated."
stompx_eor stompx-eor: $(allo_eor_srcs) $(mode_eor_srcs) stomp_eor_x.f
	cp stomp_eor_x.f stompx-eor_$(solver).f
	cat $(mode_eor_srcs) >> stompx-eor_$(solver).f
	echo "Source stompx-eor_$(solver) has been updated."
stompx_gt stompx-gt: $(allo_gt_srcs) $(mode_gt_srcs) stomp_gt_x.f
	cp stomp_gt_x.f stompx-gt_$(solver).f
	cat $(mode_gt_srcs) >> stompx-gt_$(solver).f
	echo "Source stompx-gt_$(solver) has been updated."
stompx_w stompx-w: $(allo_w_srcs) $(mode_w_srcs) stomp_w_x.f
	cp stomp_w_x.f stompx-w_$(solver).f
	cat $(mode_w_srcs) >> stompx-w_$(solver).f
	echo "Source stompx-w_$(solver) has been updated."
else
stompx_hydt_ke stompx-hydt-ke:  $(allo_hydt_ke_objs) $(mode_hydt_ke_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-hydt-ke_$(solver)$(SUFFIX_FLAG) $(allo_hydt_ke_objs) \
      $(mode_hydt_ke_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-hydt-ke_$(solver) has been updated."
stompx_hyd_o2_ch4_n2 stompx-hyd-o2-ch4-n2:  $(allo_hyd_o2_ch4_n2_objs) $(mode_hyd_o2_ch4_n2_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-hyd-o2-ch4-n2_$(solver)$(SUFFIX_FLAG) $(allo_hyd_o2_ch4_n2_objs) \
      $(mode_hyd_o2_ch4_n2_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-hyd-o2-ch4-n2_$(solver) has been updated."
stompx_co2 stompx-co2:  $(allo_co2_objs) $(mode_co2_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-co2_$(solver)$(SUFFIX_FLAG) $(allo_co2_objs) \
      $(mode_co2_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-co2_$(solver) has been updated."
stompx_eor stompx-eor:  $(allo_eor_objs) $(mode_eor_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-eor_$(solver)$(SUFFIX_FLAG) $(allo_eor_objs) \
      $(mode_eor_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-eor_$(solver) has been updated."
stompx_gt stompx-gt:  $(allo_gt_objs) $(mode_gt_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-gt_$(solver)$(SUFFIX_FLAG) $(allo_gt_objs) \
      $(mode_gt_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-gt_$(solver) has been updated."
stompx_w stompx-w:  $(allo_w_objs) $(mode_w_objs)
	$(FC) $(FFLAGS) $(LDFLAGS) -o stompx-w_$(solver)$(SUFFIX_FLAG) $(allo_w_objs) \
      $(mode_w_objs) $(LIS_LIB) $(PETSC_LIB) $(LIS_LINKING)
	echo "Executable stompx-w_$(solver) has been updated."
endif

ifeq ($(mode),debug)
.SECONDARY: $(allo_hydt_ke_srcs) $(mode_hydt_ke_srcs) \
  $(allo_hyd_o2_ch4_n2_srcs) $(mode_hyd_o2_ch4_n2_srcs) \
  $(allo_co2_srcs) $(mode_co2_srcs) \
  $(allo_eor_srcs) $(mode_eor_srcs) \
  $(allo_gt_srcs) $(mode_gt_srcs) \
  $(allo_w_srcs) $(mode_w_srcs) \
  stomp_hydt_ke_x.f stomp_hyd_o2_ch4_n2_x.f stomp_co2_x.f  \
  stomp_eor_x.f stomp_gt_x.f stomp_w_x.f \
  *.f petsc.f lis.f
endif

#
#---  Pattern Rules  ---
#

%.o:  %.F

%.f:  $(DIR)/%.FF
	/usr/bin/cpp $(CPP_FLAGS) $(CPP_FORM)

%.o:  %.f
	$(FC) $(FFLAGS) $(FCFLAGS) $(SRCH_PATH) -c $<

#
#---  End of makefile  ---
#

