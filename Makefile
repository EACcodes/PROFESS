#--------------------------------------------------------------------------------
# PROFESS 3.0 MAKEFILE 
# THIS MAKEFILE CONTAINS FOUR PARTS:
# (1) PLATFORM SPECIFIC VARIABLES.
# (2) PLATFORM INDEPENDENT VARIABLES.
# (3) OBJECT FILES.
# (4) IMPLICIT (PATTERN) RULES DEFINITIONS.
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# (1) PLATFORM SPECIFIC VARIABLES.
#--------------------------------------------------------------------------------
SHELL = /bin/sh
MACHINE = $(shell uname -m)
OS = $(shell uname)

# Linux - Opterons and Xeons
ifeq ($(MACHINE),x86_64)

#--------------------------------------------------------------------------------
# SHORT NOTE FOR FFTW PACKAGE INSTALLATION:
# Installation document for FFTW 3 (3.3.4 version is tested) can be found here: 
# http://www.fftw.org/doc/
# For PROFESS 3 version, one header file named 'fftw3.h' file is needed for
# serial version of PROFESS while 'fftw3-mpi.f03' is needed for parallel version.
# These files can be found in the 'FFTW_DIR/include/' directory after the 
# installation of FFTW 3, where 'FFTW_DIR' is your chosen direcotory of FFTW 3.
# The detail procedures to install FFTW 3 are: 
# 1) download FFTW source file from http://www.fftw.org/download.html 
# 2) enter the FFTW package and type: ./configure --prefix=/FFTW_DIR/
# (for mpi version, try ./configure --prefix=/FFTW_DIR/ --enable-mpi)
# 3) type: make instal
# 4) type: make
# (Make sure you use the same compiler for PROFESS and FFTW3)
#--------------------------------------------------------------------------------

  #FFTW_DIR=/home/your_fftw3_directory/

#--------------------------------------------------------------------------------
# SHORT NOTE FOR LIBXC PACKAGE INSTALLATION:
# PROFESS 3 depends on Libxc library (2.0.1 version is tested), 
# Important note: the interfaces change since the new version Libxc 2.1.0,
# so please use version before that.
# which can be found here:
# http://www.tddft.org/programs/octopus/wiki/index.php/Libxc
# The detail procedures to install libxc are:
# 1) download Libxc library from: 
# http://www.tddft.org/programs/octopus/wiki/index.php/Libxc:download
# 2) enter the FFTW package and type: ./configure --prefix=/XC_DIR/
# 3) type: make install
# 4) type: make
# (Make sure you use the same compiler for PROFESS and Libxc)
#----------------------------------------------------------------------------------
 
#  XC_DIR=/home/your_libxc_directory/

#--------------------------------------------------------------------------------
# SHORT NOTE FOR LAPACK PACKAGE INSTALLATION:
# If you do not plan to use Intel Math Kernel Library (Intel MKL), 
# you need to install LAPACK package by yourself.
# PROFESS 3 depends on LAPACK package (3.5.0 version is tested), 
# which can be found here:
# http://www.netlib.org/lapack/ 
# The detail procedures to install LAPACK are:
# 1) download LAPACK library from: 
# http://www.netlib.org/lapack/#_lapack_version_3_5_0_2
# 2) enter this package and copy "make.inc.example" to "make.inc" 
# 3) set the detail compiler information in "make.inc" 
# (the default compiler in make.inc.examples is gnu compiler 
# include gfortran and gcc). 
# 4) enter BLAS/SRC/ subdirectory and type: make
# this procedure will help you compile the BLAS package first,
# which is inside LAPACK package.
# 5) go back to the main directory where you have make.inc,
# and type make
# 6) then you will find both: liblapack.a and librefblas.a in 
# the main directory, then it is finished.
# P.S. For compiling LAPACK package using intel compiler,
# every procedure is the same except you need to update the compiler information 
# "FORTRAN  = ifort" and "LOADER   = ifort" in make.inc,
# and also change "TIMER    = INT_CPU_TIME" to "TIMER    = EXT_ETIME" in make.inc
# (Make sure you use the same compiler for PROFESS and LAPACK)
#----------------------------------------------------------------------------------

  #LAPACK_DIR=/home/your_lapack_directory/

#----------------------------------------------------------------------------------
# For gnu compiler (tested under version 4.4.7) 
#----------------------------------------------------------------------------------

  #---------------------
  # compiler option 
  #---------------------
  FC = gfortran

  #---------------------
  # FFLAGS flag
  #---------------------
  FFLAGS=-O3 -I${XC_DIR}include -I${FFTW_DIR}/include -cpp -J${WORKDIR}

  #---------------------
  # LDLIBS flag
  #---------------------
  # for those who want to use compiled LAPACK, Libxc and FFTW3 packages.
  LDLIBS = -L${LAPACK_DIR} -llapack -lrefblas -L${XC_DIR}/lib/ -lxcf90 -lxc -L${FFTW_DIR}/lib -lfftw3 -lm 


#----------------------------------------------------------------------------------
# For intel compiler (tested under version 13.0.1)
#----------------------------------------------------------------------------------
  #---------------------
  # compiler option 
  #---------------------
  # serial version
  #FC = ifort
  #CC = icc
  # parallel version
  #FC_par = mpif90
  #CC_par = mpicc

  #---------------------
  # FFLAGS flag
  #---------------------
  # for serial version (you can choose one of them for differernt purpose)
  # for debugging
  #FFLAGS = -g -warn all -cpp -r8 -I${XC_DIR}include -I${FFTW_DIR}/include -module ${WORKDIR}
  #FFLAGS = -g -traceback -check all -warn all -cpp -r8 -I${XC_DIR}include -I${FFTW_DIR}/include -module ${WORKDIR}
  #FFLAGS = -g -traceback  -warn all -cpp -r8 -I${XC_DIR}include -I${FFTW_DIR}/include -module ${WORKDIR}
  # for high performance
  #FFLAGS = -warn all -O3 -ipo  -cpp -r8 -I${XC_DIR}include -I${FFTW_DIR}/include -module ${WORKDIR}

  # for parallel version (you can choose one of them for differernt purpose)
  # for debugging
  #FFLAGS_PAR = -g -warn -module ${pWORKDIR} -I${XC_DIR}include -I${FFTW_DIR}/include -r8 -i_dynamic -cpp $(PARAFLAG) 
  #FFLAGS_PAR = -g -traceback -check all -warn -module ${pWORKDIR} -I${XC_DIR}/include -I${FFTW_DIR}/include -r8 -i_dynamic -cpp $(PARAFLAG) 
  # for high performance
  #FFLAGS_PAR = -warn all -O3 -xHost -ipo -module ${pWORKDIR} -I${XC_DIR}/include -I${FFTW_DIR}/include -r8 -i_dynamic -cpp $(PARAFLAG) 

  #---------------------
  # LDLIBS flag
  #---------------------
  # serial version
  # for those who want to use Intel Math Kernel Library (Intel MKL) plus compiled Libxc package 
  #LDLIBS = -shared-intel -mkl=sequential ${XC_DIR}/lib/libxc.a 
  # for those who want to use compiled LAPACK, Libxc and FFTW3 packages.
  #LDLIBS = -L${LAPACK_DIR} -llapack -lrefblas -L${XC_DIR}/lib/ -lxcf90 -lxc -L${FFTW_DIR}/lib -lfftw3 -lm
  
  # parallel version
  #LDLIBS_PAR = -L${FFTW_DIR}/lib -lfftw3_mpi -lfftw3 -lm 


endif


#----------------------------------------------------------------------------------
# (2) PLATFORM INDEPENDENT VARIABLES.
#----------------------------------------------------------------------------------

GOAL = PROFESS
PARA_GOAL = pPROFESS
CINEB_GOAL = CINEB
CONVERT_GOAL = RhoConvert

# containing object files for serial version
WORKDIR = obj
# containing object files for parallel version
pWORKDIR = pobj
# source files
SOURCEDIR = Source

VPATH = $(SOURCEDIR)

PARAFLAG = -D__USE_PARALLEL


#----------------------------------------------------------------------------------
# (3) OBJECT FILES.
#----------------------------------------------------------------------------------
# Definitions
# The relative path of the object files
FFT_OBJS= Constants.o Fourier_v1.o Fourier.o
FFT_POBJS= Constants.o Fourier_para_v1.o Fourier.o

COMMON_OBJS = MathNMS.o MathLBFGS.o MathRKsuite.o MathDCstep.o MathDCsrch.o \
       OutputFiles.o Timer.o MPI_Functions.o \
	   MathSplines.o MathFunctions.o MathLNsrch.o \
	   CellInfo.o LocalPseudoPot.o SetupFFT.o \
	   System.o Output.o KEDF_KernelODE.o PlaneWave.o CBSpline.o \
	   XC_LDA.o XC_PBE.o Hartree.o \
	   IonElectronSpline.o IonElectron.o ReadAtomDen.o Ewald.o \
	   KEDF_TF.o KEDF_VW.o KEDF_GGA.o KEDF_WTkernel.o KEDF_WGCkernel.o \
	   KEDF_WGC.o KEDF_CAT.o KEDF_WGCD.o KEDF_WT.o KEDF_Q.o \
	   KEDF_HC10.o KEDF_DenDec.o KEDF_EvW.o \
       KEDF_MGPkernel.o KEDF_MGP.o \
	   SetupKEDF.o Report.o Neighbor.o \
	   CalForce.o RefreshIons.o RefreshCell.o CalStress.o CalPotential.o\
       RhoOptimizers.o RhoDirCG.o RhoDirBFGS.o RhoDirNew.o RhoLineSearch.o \
	   RhoOptN.o RhoOptSTN.o RhoOptSCG.o RhoOptLOG.o \
	   IonOptimizers.o IonOptQui.o IonOptCON.o IonOptCG2.o IonOptBFGS.o CellOptimizers.o \
       MolecularDynamics.o MolecularDynamicsNVT.o MolecularDynamicsNPT.o MolecularDynamicsNVE.o\
	   Optimizer.o ReadIonFile.o ReadInputFile.o Initializer.o PROFESS.o 

OBJS = $(FFT_OBJS) $(COMMON_OBJS)
POBJS =$(FFT_POBJS) $(COMMON_OBJS) 

CINEBS = Constants.o Timer.o Fourier_v1.o Fourier.o MathSplines.o MathFunctions.o\
         CINEB.o

CONVS = RhoConvert.o


#----------------------------------------------------------------------------------
# (4) IMPLICIT (PATTERN) RULES DEFINITIONS.
#----------------------------------------------------------------------------------
# for serial version
$(WORKDIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(WORKDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@
    
$(WORKDIR)/%.o : %.F90
	$(FC) -c $(FFLAGS) $< -o $@

$(WORKDIR)/%.o : %.c
	$(CC) -c  $(CFLAGS) $< -o $@

# for parallel version
$(pWORKDIR)/%.o : %.f
	$(FC_par) -c $(FFLAGS_PAR) $< -o $@

$(pWORKDIR)/%.o : %.f90
	$(FC_par) -c $(FFLAGS_PAR) $< -o $@

$(pWORKDIR)/%.o : %.c
	$(CC_par) -c $< -o $@

# Absolute path of the object files (with working dir attached)
SER_ABS_OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(OBJS))
PAR_ABS_OBJS = $(patsubst %.o, $(pWORKDIR)/%.o, $(POBJS))
CINEB_ABS_OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(CINEBS))
CONV_ABS_OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(CONVS))

# Rules
default :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi 
	@ chmod a+rwX $(WORKDIR)       # this is for pathscale compiler
	@ $(MAKE) $(GOAL)                               

cineb :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi 
	@ chmod a+rwX $(WORKDIR)       # this is for pathscale compiler
	@ $(MAKE) $(CINEB_GOAL) 

parallel :
	@ if [ ! -d $(pWORKDIR) ]; then mkdir $(pWORKDIR); fi
	@ chmod a+rwX $(pWORKDIR)       # this is for pathscale compiler
	@ $(MAKE) $(PARA_GOAL)

convert :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi 
	@ chmod a+rwX $(WORKDIR)       # this is for pathscale compiler
	@ $(MAKE) $(CONVERT_GOAL)

$(GOAL) : $(SER_ABS_OBJS)
	$(FC) $(FFLAGS) -o $(GOAL) $(SER_ABS_OBJS) $(LDLIBS)

$(CINEB_GOAL) : $(CINEB_ABS_OBJS) 
	$(FC) $(FFLAGS) -o CINEB $(CINEB_ABS_OBJS) $(LDLIBS)

$(PARA_GOAL) : $(PAR_ABS_OBJS)
	$(FC_par) $(FFLAGS_PAR) -o $(PARA_GOAL) $(PAR_ABS_OBJS) $(LDLIBS_PAR) $(LDLIBS)

$(CONVERT_GOAL) : $(CONV_ABS_OBJS)
	$(FC) $(FFLAGS) -o $(CONVERT_GOAL) $(CONV_ABS_OBJS)

.PHONY : clean
clean : 
	@ if [ -d $(WORKDIR) ]; then rm -rf $(WORKDIR); fi
	@ if [ -d $(pWORKDIR) ]; then rm -rf $(pWORKDIR); fi
	@ if [ -f $(GOAL) ] ; then rm $(GOAL); fi
	@ if [ -f $(PARA_GOAL) ]; then rm $(PARA_GOAL); fi
	@ if [ -f $(CINEB_GOAL) ] ; then rm $(CINEB_GOAL); fi
	@ if [ -f $(CONVERT_GOAL) ] ; then rm $(CONVERT_GOAL) ; fi
