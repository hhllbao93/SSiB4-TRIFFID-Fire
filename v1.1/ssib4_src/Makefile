EXECUTABLE = ../a.out

OBJ_SSIB= \
    module_comsdef.o          \
    ssib4_module_comsconst.o        \
    ssib4_module_veg9table.o        \
    ssib4_module_comsveg.o          \
    ssib4_module_comssib4sfc.o      \
    ssib4_module_trifparms.o        \
    ssib4_module_fireparms.o	    \
    ssib4_module_obsco2.o           \
    ssib4_module_vegcnvt.o          \
    module_time.o             \
    module_string.o           \
    module_netcdf.o           \
    module_radc2.o            \
    module_records.o          \
    ssib4_module_sitriffid.o        \
    ssib4_module_ssibco2.o          \
    ssib4_module_ssibsub.o          \
    ssib4_module_updatepara.o       \
    ssib4_module_weighting.o        \
    ssib4_module_ssibsfcbio.o       \
    ssib4_module_pressibtrif.o      \
    ssib4_module_fireseason.o       \
    ssib4_module_fireimp.o          \
    ssib4_module_sfctrif.o          \
    ssib4_module_sfcssib2.o         \
    ssib4_module_ssib4sfc.o         \
    module_ssibio.o           \
    module_read1deg_clim.o    \
    module_writenc.o          \
    module_ssib4_main_netcdf.o  

OBJS = ${OBJ_SSIB}
CP = cp
RM = rm
CD = cd
CPP = /usr/bin/cpp
FC = ifort
LIB = -L/glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1/lib -lnetcdff -lnetcdf
FINCS=-I/glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1/include
CPP_FLAGS = -C -nostdinc
FORT_FLAGS =  
LOAD_FLAGS =  

INCLUDE_MYDIR = -I../.. -I../include

${EXECUTABLE} : ${OBJS} 
	${F77} ${LOAD_FLAGS} -o $@ ${OBJS} $(LIB) 

.SUFFIXES : .o .f .F90


.SUFFIXES : .o .f .F90

.F90.o :
	${CPP}  $(CPP_FLAGS) ${INCLUDE_MYDIR} $(FINCS) $*.F90 >$*.f90
#	sed -i '4,45d' $*.f90 
	${FC} ${FORT_FLAGS} $(FINCS) -c $*.f90

.f.o :
	${FC} ${FORT_FLAGS} $(FINCS) -c $*.f $(LIB)
	 ${FC} -qsuffix=f=f90 ${FORT_FLAGS} -c $*.f90  $(LIB)

clean :
	 - ${RM} -f  *.o *.f90  *.mod  ${EXECUTABLE}


