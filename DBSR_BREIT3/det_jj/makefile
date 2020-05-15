#         make  det_exp_jj

include ${HOME}/DEFS_03

det_exp_jj: det_exp_jj.f90
	$(FC) -o det_exp_jj $(link) det_exp_jj.f90 $(ZCONFJJ) $(ZCOM)
	mv det_exp_jj $(BINDIR)

