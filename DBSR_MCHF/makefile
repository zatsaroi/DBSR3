include ${HOME}/DEFS_03

#  make  PROGRAM with name:
 
N  = dbsr_mchf
 
# Libraries to be searched:
 
L  =   $(DBS) $(ZCONFJJ) $(ZCOM) $(LAPACK)
 
S  = MOD_dbsr_mchf.f90         MOD_df_orbitals.f90       MOD_average_energy.f90     \
	MOD_hd_core.f90           MOD_c_data.f90                                       \
	apply_orthogonality.f90   boundary_conditions.f90   check_tails.f90     \
	core_energy.f90           dbsr_mchf.f90             def_blocks.f90      \
	def_orbitals.f90          diag.f90                  get_case.f90        \
	get_estimates.f90         get_spline_param.f90                          \
	read_ang_coef.f90         read_conf_jj.f90          read_pqbs.f90       \
	rotate.f90                solve_eiv.f90             mchf_matrix.f90     \
	solve_nr.f90              solve_direct.f90          scf_mchf.f90        \
	scf_mchf_all.f90          summry.f90                rk_df.f90

	
O = $(S:.f90=.o)
 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
#	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp) $<
 
#	$(FC) $(comp) -fcheck=bounds -fcheck=all -Wall $<
