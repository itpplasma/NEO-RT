#FC       = f95f
#FC       = f95-lah
FC = gfortran
#DEBUGFLAG = -g -ggdb -C -p -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check
DEBUGFLAG = -g -ggdb -C -p -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fbounds-check
#LDFLAGS = -L/usr/local/opt/lapack/lib
#CPPFLAGS = -I/usr/local/opt/lapack/include
#FLIBS = -L/usr/local/Cellar/gcc/9.2.0_1/lib/gcc/9
#OPTS= -M OBJS --chk a,e,s,u,x --trace --trap -g
#OPTS= -M OBJS -O 
#OPTS= -std=gnu -M OBJS -O3
#OPTS= -M OBJS -Waliasing -Wampersand  -Wline-truncation  -Wnonstd-intrinsics  -Wsurprising -Wno-tabs  -Wunderflow #-Wall# -Wunused-parameter -Wconversion -Wimplicit-interface -Wcharacter-truncation
#OPTS=-Wunused
OPTS = -J OBJS $(DEBUGFLAG) -O
	

OBJS =  OBJS/mag_surf_module.o \
	OBJS/search_grid_mod.o \
	OBJS/file_reader_module.o \
	OBJS/intsrc.o \
	OBJS/plag_coeff.o \
	OBJS/integrator_stw.o \
	OBJS/spl_three_to_five_mod.o \
	OBJS/binsrc.o \
	OBJS/magdata_in_symfluxcoord.o \
	OBJS/boozer_data_in_symfluxcoord.o \
	OBJS/window_filter.o \
	OBJS/extender.o \
	OBJS/write_eqfile.o \
	OBJS/boozer_to_efit.o 
	
boozer_transform.x: $(OBJS) Boozer_make.mk
	$(FC) $(OPTS) -o boozer_transform.x $(OBJS) -llapack
OBJS/mag_surf_module.o: SRC/mag_surf_module.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/mag_surf_module.f90
	mv mag_surf_module.o OBJS/
OBJS/search_grid_mod.o: SRC/search_grid_mod.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/search_grid_mod.f90
	mv search_grid_mod.o OBJS/
OBJS/file_reader_module.o: SRC/file_reader_module.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/file_reader_module.f90
	mv file_reader_module.o OBJS/
OBJS/intsrc.o: SRC/intsrc.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/intsrc.f90
	mv intsrc.o OBJS/
OBJS/plag_coeff.o: SRC/plag_coeff.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/plag_coeff.f90
	mv plag_coeff.o OBJS/
OBJS/integrator_stw.o: SRC/integrator_stw.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/integrator_stw.f90
	mv integrator_stw.o OBJS/
OBJS/spl_three_to_five_mod.o: SRC/spl_three_to_five_mod.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/spl_three_to_five_mod.f90
	mv spl_three_to_five_mod.o OBJS/
OBJS/binsrc.o: SRC/binsrc.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/binsrc.f90
	mv binsrc.o OBJS/
OBJS/magdata_in_symfluxcoord.o: SRC/magdata_in_symfluxcoord.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/magdata_in_symfluxcoord.f90
	mv magdata_in_symfluxcoord.o OBJS/
OBJS/boozer_data_in_symfluxcoord.o: SRC/boozer_data_in_symfluxcoord.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/boozer_data_in_symfluxcoord.f90
	mv boozer_data_in_symfluxcoord.o OBJS/
OBJS/window_filter.o: SRC/window_filter.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/window_filter.f90
	mv window_filter.o OBJS/
OBJS/extender.o: SRC/extender.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/extender.f90
	mv extender.o OBJS/
OBJS/write_eqfile.o: SRC/write_eqfile.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/write_eqfile.f90
	mv write_eqfile.o OBJS/
OBJS/boozer_to_efit.o: SRC/boozer_to_efit.f90 Boozer_make.mk
	$(FC) $(OPTS) -c SRC/boozer_to_efit.f90
	mv boozer_to_efit.o OBJS/
#OBJS/.o: SRC/.f90 Boozer_make.mk
#	$(FC) $(OPTS) -c SRC/.f90
#	mv .o OBJS/
