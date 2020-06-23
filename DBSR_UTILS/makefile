#       utilities DBSR

include ${HOME}/DEFS_03

all: bsw_rw       cm_j       coef_jj     dbound_bsw    dbound_tab   \
	dbsr_merge   dbsw123                genjconf      genjterm     \
	jcfile       jfile      jj_list     mult_jj       order_jj     \
	rw123        rw_bsw     rw_dat      zf_bb         zf_cc

# dbsw_tab

bsw_rw: bsw_rw.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

cm_j: cm_j.f90
	$(FC) -o $@ $(link)  $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

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

dbsw123: dbsw123.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

#dbsw_tab: dbsw_tab.f90
#	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
#	mv $@ $(BINDIR)

genjconf: genjconf.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

genjterm: genjterm.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jcfile: jcfile.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jfile: jfile.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jj_list: jj_list.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

mult_jj: mult_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

order_jj: order_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

rw123: rw123.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_bsw: rw_bsw.f90
	$(FC) -o $@ $(link) $< $(DBS)  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_dat: rw_dat.f90    
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_bb: zf_bb.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_cc: zf_cc.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

# -fbacktrace -fcheck=bounds 
