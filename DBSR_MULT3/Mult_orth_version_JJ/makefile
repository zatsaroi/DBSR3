FC = gfortran 
link = -O2 -I.    -fbacktrace    -fcheck=bounds
comp = -c -O2 -I.  -fbacktrace    -fcheck=bounds
 
#  make  PROGRAM with name:
 
N  = test_mult_2conf_jj
 
S  = mod_det_jq.f90 mult_2conf_jj.f90  test_mult_2conf_jj.f90
	
O = $(S:.f90=.o)   
 
$(N): $(O)
	$(FC) -o $(N) $(link) $(O)  
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<
 
# -fbacktrace  -fcheck=bounds
 
