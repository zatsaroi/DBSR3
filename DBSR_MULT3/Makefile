include ${HOME}/DEFS_03
 
N  = dbsr_mult3
 
# Libraries to be searched:
 
L  =  $(ZCONFJJ) $(ZCOM)

S  = mod_dbsr_mult.f90          mod_nljm.f90               mod_term_exp.f90           \
	mod_coef.f90               check_mult_bnk.f90         add_res.f90                \
	conf_loop.f90              dbsr_mult.f90              det_exp.f90                \
	det_mult.f90               read_arg.f90               term_loop.f90              \
	Idet_fact.f90              me_jj.f90

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
	$(FC) $(comp) $<
 

