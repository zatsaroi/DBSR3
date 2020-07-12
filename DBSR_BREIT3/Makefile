include ${HOME}/DEFS_03
 
N  = dbsr_breit3 
 
# Libraries to be searched: 
 
L  = $(ZCONFJJ) $(ZCOM) 

S  = mod_dbsr_breit.f90   mod_nljm.f90         mod_term_exp.f90    \
     mod_coef.f90         mod_boef_list.f90    read_arg.f90        \
     dbsr_breit.f90       add_res.f90          conf_loop.f90       \
     det_exp.f90          det_me.f90           read_conf_jj.f90    \
     term_loop.f90        inf_dbsr_breit.f90   me_jj.f90

O = $(S:.f90=.o) 
 
$(N): $(O) 
	$(FC) -o $(N) $(O) $(link) $(L) 
	mv $(N) $(BINDIR) 
	rm *.o *.mod  
 
clean:
	rm -f *.o *.mod  
 
.SUFFIXES: 
.SUFFIXES: .f90 .mod .o 
 
.f90.o: 
	$(FC) $(comp)  $< 
 

