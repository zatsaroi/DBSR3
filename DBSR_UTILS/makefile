#       utilities DBSR

include ${HOME}/DEFS_03

all: coef_jj           dbound_bsw          dbound_tab        \
	dbsr_merge        mult_jj             bsw_rw            \
     bsw123            rw_bsw                                \
     rw_dat            rw123

coef_jj: coef_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

dbound_bsw: dbound_bsw.f90    
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

dbound_tab: dbound_tab.f90
	$(FC) -o $@ $(link) $< $(ZCOM)
	mv $@ $(BINDIR)

dbsr_merge: dbsr_merge.f90        
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

dbsw_tab: dbsw_tab.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

mult_jj: mult_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)


photo_jj: photo_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

bsw_rw: bsw_rw.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

bsrw123: bsrw123.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_bsw: rw_bsw.f90
	$(FC) -o $@ $(link) $< $(DBS)  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw123: rw123.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_bb: zf_bb.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_cb: zf_cb.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_cc: zf_cc.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_dat: rw_dat.f90    
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)
