include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_hd3_mpi
 
# Libraries to be searched:

LIB = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential \
      -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl

L  =  $(DBS) $(ZCONFJJ) $(ZCOM)

S  =	mod_blacs.f90            mod_hd.f90               add_exp.f90              \
	b_out.f90                dbsr_hdb.f90             diag_mat.f90             \
	h_out.f90                check_dbsr_mat.f90       read_overlaps.f90        \
	read_matrix.f90          rsol_out.f90             read_arg.f90             \
	sub1_hd.f90              w_out.f90                channel_label.f90
 
O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) -mkl $(link) $(L) $(LIB)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) -I${MKLROOT}/include $<
 
#-CB -traceback 