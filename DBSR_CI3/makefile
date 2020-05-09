include ${HOME}/DEFS_03

N  = dbsr_ci3

# Libraries to be searched:
 
L  =  $(DBS) $(ZCONFJJ) $(ZCOM) $(LAPACK)

S  =	mod_dbsr_ci.f90	   mod_c_data.f90                      \
	add_integral.f90   dbsr_ci.f90        diag.f90         \
	gen_matrix.f90     int_corr.f90       inf_dbsr_ci.f90  \
	read_conf_jj.f90   read_data.f90                       \
	shift.f90          rk_evaluate.f90 

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod *~
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<

# -fbacktrace -fbounds-check 

