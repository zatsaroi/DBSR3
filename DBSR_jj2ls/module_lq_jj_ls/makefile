
include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = ls_jj_nl
 
# Libraries to be searched:
 
L  =  $(ZCONFJJ) $(ZCONFLS) $(SHELLS) $(ZCOM) 

S  =	read_ln.f90  ls_jj_nl.f90  
 
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
	$(FC) $(comp)  $<

# -fcheck=bounds -fbacktrace 
