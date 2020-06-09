include ${HOME}/DEFS_03
 
#  make  PROGRAM with name:
 
N  = jj2ls
 
# Libraries to be searched:

L  =  $(ZCONFJJ) $(ZCONFLS) $(ZCOM) 
 
S  =	mod_jj2ls.f90            lq_jj_ls.f90             check_jj_shells.f90      \
	check_ovl.f90            cm2j.f90                 get_LS_conf.f90          \
	get_LS_terms.f90         get_tr_coef.f90          g_factor.f90             \
	jj2ls.f90                make_couplings.f90       transform_jfile.f90      \
	recup_3j.f90             recup_jj_shells.f90      sub_nlLS.f90             

O = $(S:.f90=.o) recup_a.o  
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
#	rm *.o *.mod 

recup_a.o: recup_a.f
	$(FC) $(comp) recup_a.f
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<
 
.f.o:
	$(FC) $(comp)  $<
 
