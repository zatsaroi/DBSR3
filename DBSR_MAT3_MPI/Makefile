include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_mat3a

# Libraries to be searched:
 
L  =  $(DBS) $(ZCONFJJ) $(ZCOM) $(LAPACK)
 
S  = mod_dbsr_mat.f90        dbsr_mat.f90            L_data.f90      \
     O_data.f90              R_data.f90              S_data.f90      \
     add_integral.f90        br_fake.f90             check_mat.f90   \
     dbs_orth.f90            diag.f90                gen_matrix.f90  \
     idef_type.f90           read_arg.f90            read_data.f90   \
     record_matrix.f90       state_res.f90           sub1.f90        \
     target_subs.f90         transform_matrix.f90    updates.f90     \
     zero_cond.f90           

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
 

# $(FC) $(comp) -fcheck=bounds -fbacktrace  $<

#-fno-automatic -fcheck=bounds -fbacktrace -fcheck=mem
#-traceback -CB