include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_breit3_mpi

# Libraries to be searched:
 
L  =  $(ZCONFJJ) $(ZCOM) 

S  =	mod_param.f90            mod_boef.f90             mod_coef.f90             \
	mod_nljm.f90             mod_term_exp.f90                                  \
	add_res.f90              bj_mpi.f90               br_conf_jj.f90           \
	conf_calc.f90            conf_loop.f90            det_exp.f90              \
	det_jq.f90               det_me.f90               get_res.f90              \
	open_jj.f90              read_conf_jj.f90         term_loop.f90            \
 
O = $(S:.f90=.o)

$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC_mpi) $(comp) $<

 
