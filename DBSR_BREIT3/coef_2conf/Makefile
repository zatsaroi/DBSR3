include ${HOME}/DEFS_03
 
N  = test_coef_2conf_jj
 
# Libraries to be searched:
 
L  =  $(ZCONFJJ) $(ZCOM)

S  =	mod_boef_list.f90 test_coef_2conf_jj.f90  coef_2conf_jj.f90 

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
 
# -fbacktrace -fbounds-check   (gfortran)

