include ${HOME}/DEFS_03
 
N  = coef_jj
 
# Libraries to be searched:
 
L  =  $(ZCONFJJ) $(ZCOM)

S  =	 coef_jj.f90   

O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o  
 
clean:
	rm -f *.o *.mod *~
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) -fcheck=all -Wall -fbacktrace -fcheck=bounds $<
 

