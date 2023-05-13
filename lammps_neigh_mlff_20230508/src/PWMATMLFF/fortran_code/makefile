# ugly makefile, but do not change this makefile unless you have full test it.
.PHONY:all
#all: cmain.x mlff_predict.x main_MD.x cmain.a
#all: cmain.x mlff_predict.x main_MD.x pair_qcad.a
all: f2c_calc_energy_force.a

MPIFC = mpiifort
MPICXX = mpiicpc

is_mpiifort := $(shell which mpiifort 1>&2 2> /dev/null; echo $$?)
ifneq ($(is_mpiifort),0)
	MPIFC=mpif90
endif

FFLAG = -O3 -I/opt/iotk/src/ -mkl
#FFLAG = -g -traceback -C -I/opt/iotk/src/ -mkl -debug all
CFLAG = -O3 --std=c++11 -mkl -I../
FLINKFLAG = -cxxlib -nofor_main
#SRCS = $(wildcard *.f90)
#OBJ = $(patsubst %.f90,%.o,$(SRCS))
# priority-sensitive
OBJ =   dp_ff_mod.o mod_matom.o mod_data.o \
        gen_deepMD_feature.o calc_deepMD_forceM.o\
        ML_FF_EF.o \
        f2c_calc_energy_force.o  

# priority-sensitive
MODS = calc_deepmd_f.mod calc_deepMD.mod mod_data.mod mod_parameter.mod 

%.o:  %.f90
	$(MPIFC) $(FFLAG) -c $<
%.o:  %.f
	$(MPIFC) $(FFLAG) -c $<
%.o: %.cpp
	$(MPICXX) $(CFLAG) -c $<

calc_deepMD.mod: calc_deepMD_forceM.f90
	$(MPIFC) $(FFLAG) -c $^
calc_deepMD_f.mod: gen_deepMD_feature.f90
	$(MPIFC) $(FFLAG) -c $^
mod_data.mod: mod_data.f90
	$(MPIFC) $(FFLAG) -c $^
mod_parameter.mod: mod_matom.f90
	$(MPIFC) $(FFLAG) -c $^

mods: mod_parameter.mod calc_deepMD.mod \
      mod_data.mod mod_md.mod calc_deepmd_f.mod

#mlff_predict.x: mods $(OBJ) mlff_predict.o
#	make -C iotk
#	\cp iotk/src/libiotk.a libiotk.a
#	$(MPIFC) $(FFLAG) libiotk.a -o $@ $(OBJ) mlff_predict.o

#cmain.x: mods $(OBJ) cmain.o f2c_calc_energy_force.a
#	\cp iotk/src/libiotk.a libiotk.a
#	$(MPIFC) $(FFLAG) $(FLINKFLAG) libiotk.a -o $@ $(OBJ) cmain.o
#cmain.x: cmain.o f2c_calc_energy_force.a
#	$(MPIFC) $(FFLAG) $(FLINKFLAG) f2c_calc_energy_force.a -o $@ $(OBJ) cmain.o

libiotk.a:
	make -C iotk/
	\cp iotk/src/libiotk.a libiotk.a

# the $(OBJ) contains f2c_calc_energy_force.o
# this f2c_calc_energy_force.a is used in pair_qcad.cpp
f2c_calc_energy_force.a: $(mods) $(OBJ) libiotk.a
	ar ruv $@ $(OBJ) libiotk.a
	\cp f2c_calc_energy_force.a ../
	
# this the fortran MD/MD100(inference) main program
#main_MD.x: $(mods) $(OBJ) libiotk.a main_MD.o
#	$(MPIFC) $(FFLAG) $(FLINKFLAG) libiotk.a -o $@ $(OBJ) main_MD.o

clean:
	#make clean -C iotk
	rm -f *.a *.o *.mod *.x 
