include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_hd3
 
# Libraries to be searched:

L  =  $(DBS) $(ZCONFJJ) $(ZCOM) $(LAPACK)
 
S  =	mod_dbsr_hd.f90          add_exp.f90              b_out.f90                \
	dbsr_hd.f90              diag_mat.f90             h_out.f90                \
	check_dbsr_mat.f90       rsol_out.f90             read_arg.f90             \
	sub1_hd.f90              w_out.f90                channel_label.f90 
 
O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O)  $(link) $(L)
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
