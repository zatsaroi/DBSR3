include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_prep3
 
# Libraries to be searched:

L  =  $(DBS) $(ZCONFJJ) $(ZCOM) 

S  =	mod_dbsr_prep.f90       dbsr_prep.f90                                      \
     add_sub_orbital.f90     assign_index.f90      check_bsr_name.f90           \
     dbsr_prep_inf.f90       read_arg.f90          sub_check_orb.f90            \
     sub_phys_orb.f90        sub_phys_pert.f90     write_targ_example.f90

O = $(S:.f90=.o)
 
$(N): $(O) 
	$(FC) -o $(N)  $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) $<

# -fbounds-check  -fbacktrace
#-fbounds-check  -fbacktrace -fcheck=all -Wall 