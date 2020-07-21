!======================================================================
      Subroutine read_shift
!======================================================================
!     fine-turning of interaction matrix by shifting the configuration
!     energies. List of shifting configurations are given in input file
!     (name.inp_ci) as:
!     
!     shift:
!     met  S  unit  jot1  jot2
!     
!     At the moment, one option: shifting all atomic configuration
!     states (ACS) with the same "pure" configuration.
!----------------------------------------------------------------------
      Use dbsr_ci
      Use conf_jj
        
      Implicit none
      Integer :: ic,j, jot1,jot2
      Real(8) :: S
      Character :: unit*2, met*12
      Integer, external :: Icheck_file, Iadd_symc

      if(inp.eq.0) Return
      if(allocated(shift)) Deallocate(shift)
      Allocate(shift(ncfg)); shift = 0.d0; nshift = 1

      rewind(inp)
    1 read(inp,'(a)',end=2) AS
      if(AS(1:6).ne.'shift:') go to 1
      read(AS(7:),*) met,S,unit,jot1,jot2

      if(unit.eq.'eV') S = S / au_eV
      if(unit.eq.'cm') S = S / au_cm

      Select case(met)

      Case('conf')
       read(inp,'(a)') CONFIG
       Call Decode_cj
       no1 = no; iq1=iq; kn1=kn
       Do j = 1,njbl; jot = JJc(j)
        if(jot1.ge.0.and.jot1.gt.jot) Cycle
        if(jot2.ge.0.and.jot2.lt.jot) Cycle
        iconf1 = Iadd_symc(jot,no1,iq1,kn1)
        Do ic=JTc1(j),JTc2(j)
         Call Get_cfg_jj(ic)
         if(iconf.ne.iconf1) Cycle
         shift(ic) = shift(ic) + S
        End do
       End do
      
      Case default
       Stop 'Stop in read_shift: unknoun case '
      End Select

      go to 1
    2 Continue

      End Subroutine read_shift


