MPIF90 = mpif90
FCFLAGS = -fdefault-double-8 -fdefault-integer-8 -finit-integer=0 -finit-real=zero -c -fallow-argument-mismatch
FCFLAGS_F77 = -ffixed-form -ffixed-line-length-none -fallow-argument-mismatch
FCFLAGS_F90 = -fdefault-real-8 -ffree-form -ffree-line-length-none -fallow-argument-mismatch
LDFLAGS = -fdefault-real-8 -fdefault-integer-8 -finit-integer=0 -finit-real=zero -lmpi -lmpi -lfftw3 -lm

# Program name
PROG    = hit36

# Modules
MODULES = m_openmpi.o\
	m_io.o\
	m_parameters.o\
	m_work.o\
	m_fields.o\
	m_timing.o\
	x_fftw.o\
	m_filter_xfftw.o\
	m_particles.o\
	m_stats.o\
	m_force.o\
	m_rand_knuth.o\
	m_les.o\
	RANDu.o

# Objects
OBJ     = main.o\
	begin_new.o\
	begin_restart.o\
	dealias_all.o\
	get_file_ext.o\
	init_velocity.o\
	init_scalars.o\
	io_write_4.o\
	my_dt.o\
	my_exit.o\
	pressure.o\
	restart_io.o\
	rhs_velocity.o\
	rhs_scalars.o\
	velocity_rescale.o\
	write_tmp4.o


# -------------------------------------------------------
# link

$(PROG):  $(MODULES) $(OBJ)
	$(MPIF90) $(MODULES) $(OBJ) -o $(PROG) $(LDFLAGS)
# -------------------------------------------------------
# compile

$(OBJ): $(MODULES)

%.o: %.f
	$(MPIF90) $(FCFLAGS) $(FCFLAGS_F77)  $<

%.o: %.f90
	$(MPIF90) $(FCFLAGS) $<

clean:
	rm *.o *.mod $(PROG)
