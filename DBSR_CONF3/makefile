include ${HOME}/DEFS_03
 
#  make PROGRAM with name:
 
N  = dbsr_conf3 
M  = dbsr_conf3_mpi 

L  = $(ZCONFJJ) $(ZCOM) 
 
S  = mod_dbsr_conf.f90     gen_conf.f90       make_coupling.f90                     \
     check_perturber.f90   def_phys_targ.f90  read_arg.f90      sub_jj.f90          \
     record_orth.f90       def_orth_cond.f90  def_del.f90       check_cfg.f90       \
     dbsr_conf_inf.f90     
 
SN = dbsr_conf.f90                 
SM = dbsr_conf_mpi.f90   br_arg.f90  br_conf_jj.f90  br_phys_orb_jj.f90  br_target_jj.f90
	
O  = $(S:.f90=.o)  
ON = $(SN:.f90=.o)  
OM = $(SM:.f90=.o)  

$(N): $(O) $(ON)
	$(FC) -o $(N) $(ON) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 

$(M): $(O) $(OM)
	$(FC_mpi) -o $(M) $(OM) $(O) $(link) $(L) $(BSR_MPI)
	mv $(M) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<

# -fbounds-check  -fbacktrace
# -traceback -CB 