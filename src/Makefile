# fortran compiler
FC = mpif90
#FC = ftn
#FC = mpiifort
#
BIG    := #-mcmodel=medium
DBG    := #-O0 -g -fbacktrace -ffpe-trap=invalid -Wall -Wextra -pedantic -fcheck=all #-finit-real=snan
OMP    := #-fopenmp
OPT    := -O3 --fast-math -march=native 
#OPT    := -O3 -ipo -xHost
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_ITER_PI -DLOW_MACH 
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_ITER_MO -DLOW_MACH 
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_ITER_LA -DLOW_MACH 
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_DIRECT -DTIMING #-DLOW_MACH 
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_DIRECT -DLOW_MACH -DEVAP #-DENERGY 
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_DIRECT -DLOW_MACH -DENERGY -DVAR_MUG -DVAR_KAG -DEVAP -DVAP_MASS -DVAR_D_LG
OTH    :=  -fdefault-real-8 -fdefault-double-8 -DFFT_DIRECT -DUSE_VOF -DENERGY -DEVAP -DVAP_MASS -DTWOD #-DVAR_D_LG -DLOW_MACH #-DTURB_FORCING
#OTH    := -r8 -fpconstant -O3 -132 -DTWOD -DFFT_DIRECT -DLOW_MACH -DENERGY #-DEVAP -DVAP_MASS -DVAR_D_LG
#OTH    := -r8 -fpconstant -O3 -132 -DTWOD -DFFT_DIRECT -DENERGY #-DEVAP -DVAP_MASS -DVAR_D_LG
#OTH    := -fdefault-real-8 -fdefault-double-8 -DTWOD -DFFT_DIRECT -DEVAP #-DVAP_MASS #-DENERGY 
#OTH    := -DFFT_DIRECT -DTIMING #-DLOW_MACH 
PROF   := #-O0 -g -pg
FFLAGS := -cpp $(OPT) $(BIG) $(DBG) $(PROF) $(OMP) $(OTH)
LFLAGS := -fopenmp
#
TARGET = cans
#
OBJ = bound.o chkdiv.o chkdt.o cmpt_divth.o common_mpi.o correc.o debug.o fft.o fftw.o fillps.o funcs.o gradls.o initvof.o initflow.o initgrid.o initmpi.o initsolver.o inte_vel.o load.o main.o mom.o moms.o output.o param.o rk.o rks.o sanity.o solver.o tagging.o thermo.o vof.o vof_to_ls.o decomp_2d.o io.o
#
#LIBS = -L /scratch/nicolos/nicolos_phd_files/software/fftw-3.3.7/lib -lfftw3 -L /scratch/nicolos/nicolos_phd_files/software/hypre-2.15.1/src/lib -lHYPRE 
#LIBS = -L /software/sse/manual/HYPRE/2.11.2/nsc1-intel-2018a-eb/lib -lHYPRE -L /software/sse/easybuild/prefix/software/FFTW/3.3.6-intel-2018a-nsc1/lib -lfftw3
#LIBS = -L /pdc/vol/fftw/3.3.7/gcc/7.2.0/openmpi/3.0/double/lib/ -lfftw3 -L /pdc/vol/hypre/2.11.2/intel/15.0.2/lib -lHYPRE
#LIBS =  -L /cluster/work/users/salarzs/fftw-3.3.9/lib -lfftw3 -L /cluster/work/users/salarzs/hypre/src/lib -lHYPRE  
LIBS =  -L /home/salarz/scratch/fftw-3.3.9/lib -lfftw3
#LIBS = -L /cluster/work/users/nicolos/software/fftw-3.3.9/lib -lfftw3 -L /cluster/work/users/nicolos/software/hypre/src/lib -lHYPRE 
#LIBS = -L /home/nicolo/Documents/numerical_libraries/fftw-3.3.8/lib -lfftw3 -L /home/nicolo/Documents/numerical_libraries/hypre/src/lib -lHYPRE
# nb: the user should change the line LIBS accordingly with the path for FFTW3 and HYPRE libraries of his/her machine
#
all: $(TARGET)
#
$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(LFLAGS) $(OBJ) $(LIBS) -o $(TARGET)
#
.PHONY: clean
clean:
	rm -rf *.o *.mod *dSYM $(TARGET)
#
bound.o: bound.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  bound.f90
chkdiv.o: chkdiv.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  chkdiv.f90
chkdt.o: chkdt.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  chkdt.f90
cmpt_divth.o: cmpt_divth.f90 bound.o param.o thermo.o
	$(FC) $(FFLAGS) -c  cmpt_divth.f90
common_mpi.o: common_mpi.f90 param.o
	$(FC) $(FFLAGS) -c  common_mpi.f90
correc.o: correc.f90 param.o
	$(FC) $(FFLAGS) -c  correc.f90
debug.o: debug.f90 common_mpi.o
	$(FC) $(FFLAGS) -c  debug.f90
fillps.o: fillps.f90 param.o
	$(FC) $(FFLAGS) -c  fillps.f90
fft.o: fft.f90 common_mpi.o fftw.o
	$(FC) $(FFLAGS) -c  fft.f90
fftw.o: fftw.f90
	$(FC) $(FFLAGS) -c  fftw.f90
funcs.o: funcs.f90 param.o
	$(FC) $(FFLAGS) -c  funcs.f90
gradls.o: gradls.f90 
	$(FC) $(FFLAGS) -c  gradls.f90
initflow.o: initflow.f90 common_mpi.o decomp_2d.o output.o rk.o rks.o
	$(FC) $(FFLAGS) -c  initflow.f90
initgrid.o: initgrid.f90 param.o
	$(FC) $(FFLAGS) -c  initgrid.f90
initmpi.o: initmpi.f90 decomp_2d.o common_mpi.o
	$(FC) $(FFLAGS) -c  initmpi.f90
initsolver.o: initsolver.f90 common_mpi.o fft.o param.o
	$(FC) $(FFLAGS) -c  initsolver.f90 
initvof.o: initvof.f90 decomp_2d.o param.o common_mpi.o
	$(FC) $(FFLAGS) -c  initvof.f90
inte_vel.o: inte_vel.f90 param.o 
	$(FC) $(FFLAGS) -c  inte_vel.f90
load.o: load.f90 common_mpi.o decomp_2d.o io.o param.o
	$(FC) $(FFLAGS) -c  load.f90 
main.o: main.f90 bound.o chkdiv.o chkdt.o cmpt_divth.o debug.o fft.o fillps.o initflow.o initflow.o initgrid.o initmpi.o initsolver.o inte_vel.o load.o output.o param.o rk.o rks.o sanity.o solver.o solver_vc.o tagging.o thermo.o vof.o vof_to_ls.o out1d.h90 out2d.h90 out3d.h90 vol_based/time_averaging_liquid.h90 mass_based/time_averaging_gas.h90 vol_based/time_averaging_gas.h90 mass_based/time_averaging_heat_transfer.h90 vol_based/time_averaging_heat_transfer.h90 mass_based/time_averaging_mass_transfer.h90 vol_based/time_averaging_mass_transfer.h90
	$(FC) $(FFLAGS) -c  main.f90
mom.o: mom.f90 param.o funcs.o common_mpi.o 
	$(FC) $(FFLAGS) -c  mom.f90
moms.o: moms.f90 gradls.o common_mpi.o thermo.o
	$(FC) $(FFLAGS) -c  moms.f90
output.o: output.f90 common_mpi.o io.o
	$(FC) $(FFLAGS) -c  output.f90
param.o: param.f90 bc.h90
	$(FC) $(FFLAGS) -c  param.f90
rk.o: rk.f90 debug.o mom.o param.o 
	$(FC) $(FFLAGS) -c  rk.f90
rks.o: rks.f90 debug.o moms.o param.o solver_vc.o
	$(FC) $(FFLAGS) -c  rks.f90
sanity.o: sanity.f90 bound.o chkdiv.o common_mpi.o correc.o debug.o decomp_2d.o fft.o fillps.o initflow.o initmpi.o initsolver.o param.o solver.o
	$(FC) $(FFLAGS) -c  sanity.f90
solver_vc.o: solver_vc.f90 common_mpi.o param.o decomp_2d.o thermo.o
	$(FC) $(FFLAGS) -c  solver_vc.f90
solver.o: solver.f90 fft.o param.o
	$(FC) $(FFLAGS) -c  solver.f90
tagging.o: tagging.f90 bound.o common_mpi.o param.o thermo.o
	$(FC) $(FFLAGS) -c  tagging.f90
thermo.o: common_mpi.o thermo.f90 param.o
	$(FC) $(FFLAGS) -c  thermo.f90
vof.o: vof.f90 bound.o initvof.o param.o
	$(FC) $(FFLAGS) -c  vof.f90
vof_to_ls.o: vof_to_ls.f90 bound.o
	$(FC) $(FFLAGS) -c  vof_to_ls.f90
#
FFLAGS_2DECOMP := -DDOUBLE_PREC -DOVERWRITE -DEVEN
2decompdir=2decomp/
decomp_2d.o: $(2decompdir)decomp_2d.f90
	$(FC) $(FFLAGS) $(FFLAGS_2DECOMP) -c $(2decompdir)decomp_2d.f90
io.o: $(2decompdir)io.f90 decomp_2d.o
	$(FC) $(FFLAGS) $(FFLAGS_2DECOMP) -c $(2decompdir)io.f90
