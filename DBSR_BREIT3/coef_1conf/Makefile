include ${HOME}/DEFS_03
 
N  = test_coef_1conf
 
# Libraries to be searched:
 
# L  =  $(ZCONFJJ) $(ZCOM)

S  =	COEF_MOD_boef.f90    COEF_MOD_det_jq.f90  COEF_1conf.f90     \
	COEF_decode_cj.f90   COEF_jdef_ncfg.f90   COEF_Jterm.f90     \
	COEF_me_jj.f90       COM_3j_symbol.f90    COM_check_file.f90 \
	COM_Cjkj.f90         COM_EL_nljk.f90      COM_kappa.f90      \
	COM_read_arg.f90     test_coef_1conf.f90

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) # $(L)
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
 

