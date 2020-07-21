!======================================================================
      Subroutine Gen_zf
!======================================================================
! ... provide oscillator strengths for given dipole matrix
!----------------------------------------------------------------------
      Use dbsr_dmat
      Use DBS_nuclear, only: atomic_number, atomic_weight

      Implicit none
      Integer :: i,j, j1,j2, ic1,ic2, jc1,jc2, kp, &
                 nhm,kch,kpert,ns1,nc,parity
      Real(8) :: SL,SV, FL,FV, GFL,GFV, WL,WV, GWL,GWV, &
                 alfaL, alfaV, alfL, alfV, E1,E2, de,   &
                 au_cm, au_eV, S, g1,g2, ANGS, ANGSA , dmatL,dmatV
      Character( 1) :: par1, par2
      Character( 5) :: AA
      Character(64) :: Label1, Label2
      Integer, external :: Ifind_position, ITRA
      Real(8), external :: ANGS_AIR

! ... check the initial state:

      nstate1=1
      if(ctype1.eq.'b') then

       i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS1
       Call Check_file(AF); Open(nub1,file=AF,form='UNFORMATTED')
       rewind(nub1)

       read(nub1) nhm,kch,kpert,ns1,jot1,parity,nstate1

       if(ns1 .ne.ns )       Stop 'dbsr_dmat: ns1 <> ns '
       if(nch1.ne.kch)       Stop 'dbsr_dmat: nch1 --> ?'
       if(npert1.ne.kpert)   Stop 'dbsr_dmat: npert1 --> ?'
       if(kdm1.ne.nhm)       Stop 'dbsr_dmat: nhm1 --> ?'
       if(parity1.ne.parity) Stop 'dbsr_dmat: parity1 --> ?'

      elseif(ctype1.eq.'j') then

       AF = name1(1:iname1)//'j'
       Call Check_file(AF);  Open(nub1,file=AF)
       Call Read_ipar(nub1,'ncfg',nc     ) 
       if(nc.ne.kdm1) Stop 'dbsr_dmat: nc in j-file <> ncfg1'
       Call Read_ipar(nub1,'nsol',nstate1) 

      elseif(ctype1.eq.'c') then;   
      
       nstate1=1; Call Jdef_JPE(nuc1,jot1,parity1,E1) 

      else;   Stop 'Gen_zf: ctype1 --> ?'

      end if
      par1 = '+'; if(parity1.le.0) par1 = '-' 

! ... check the final state:

      nstate2=1
      if(ctype2.eq.'b') then

       i=INDEX(BF_b,'.'); AF=BF_b(1:i)//ALS2
       Call Check_file(AF); Open(nub2,file=AF,form='UNFORMATTED')
       rewind(nub2)

       read(nub2) nhm,kch,kpert,ns1,jot2,parity,nstate2

       if(ns1 .ne.ns )       Stop 'dbsr_dmat: ns2 <> ns '
       if(nch2.ne.kch)       Stop 'dbsr_dmat: nch2 --> ?'
       if(npert2.ne.kpert)   Stop 'dbsr_dmat: npert2 --> ?'
       if(kdm2.ne.nhm)       Stop 'dbsr_dmat: nhm2 --> ?'
       if(parity2.ne.parity) Stop 'dbsr_dmat: parity2 --> ?'

      elseif(ctype2.eq.'j') then

       AF = name2(1:iname2)//'j'
       Call Check_file(AF);  Open(nub2,file=AF)
       Call Read_ipar(nub2,'ncfg',nc     ) 
       if(nc.ne.kdm2) Stop 'dbsr_dmat: nc in j-file <> ncfg2'
       Call Read_ipar(nub2,'nsol',nstate2) 

      elseif(ctype2.eq.'c') then; 
      
       nstate2=1; Call Jdef_JPE(nuc2,jot2,parity2,E2)  

      else;  Stop 'Gen_zf: ctype2 --> ?'

      end if
      par2 = '+'; if(parity2.le.0) par2 = '-' 

      write(pri,'(/a,i5 )') 'nstate1 = ',nstate1
      write(pri,'( a,i5/)') 'nstate2 = ',nstate2
      if(mstate1.le.0.or.mstate1.gt.nstate1) mstate1=nstate1
      if(mstate2.le.0.or.mstate2.gt.nstate2) mstate2=nstate2

!----------------------------------------------------------------------
! ... generation of f-values:

      Open(nuf,file=AF_zf,position='APPEND')

      ! Rydberg constants:
      Call Conv_au (atomic_number,atomic_weight,au_cm,au_eV,pri)

      Deallocate(C1); Allocate(C1(kdm1)); C1=1.d0
      Deallocate(C2); Allocate(C2(kdm2)); C2=1.d0

      label1=name1(1:iname1-1); ic1=1; ic2=kdm1
      label2=name2(1:iname2-1); jc1=1; jc2=kdm2

! ... loop over initial set:

      if(ctype1.eq.'b') then
       rewind(nub1);  read(nub1) i 
      elseif(ctype1.eq.'j') then
       i=Ifind_position(nub1,'Solutions');  read(nub1,*)  
      end if

      Do j1=1,mstate1
       if(ctype1.eq.'j') then
        read(nub1,*) j,label1
        read(nub1,*) E1,jot1,ic1,ic2
        read(nub1,*) C1(ic1:ic2)
       elseif(ctype1.eq.'b') then
        read(nub1) j,label1
        read(nub1) E1
        read(nub1) C1
       end if

       if(istate1.gt.0.and.istate1.ne.j1) Cycle
       if(istate1.gt.0.and.istate1.lt.j1) Exit


       alfaL=0.0; alfaV=0.0

!... loop over final set:

      if(ctype2.eq.'b') then
       rewind(nub2); read(nub2) i
      elseif(ctype2.eq.'j') then 
       i=Ifind_position(nub2,'Solutions'); read(nub2,*)  
      end if

      Do j2=1,mstate2
       if(ctype2.eq.'j') then
        read(nub2,*) j,label2
        read(nub2,*) E2,jot2,jc1,jc2
        read(nub2,*) C2(jc1:jc2)
       elseif(ctype2.eq.'b') then
        read(nub2) j,label2
        read(nub2) E2
        read(nub2) C2
       end if

       if(istate2.gt.0.and.istate2.ne.j2) Cycle
       if(istate2.gt.0.and.istate2.lt.j2) Exit

! ... angular orthogonality: 

      if(ITRA(jot1,kpol+kpol,jot2).eq.0) Cycle

! ... transition energy:

      de = abs(E2-E1); if(de.lt.2.d-8) Cycle

! ... line strengths:

      SL=0.d0; SV=0.d0
      Do j=jc1,jc2
       S = SUM(DL(ic1:ic2,j)*C1(ic1:ic2));  SL = SL + C2(j)*S
       if(ktype.eq.'M') Cycle
       S = SUM(DV(ic1:ic2,j)*C1(ic1:ic2));  SV = SV + C2(j)*S
      End do

      if(SL.eq.0.0) Cycle

      if(kpol.eq.0) then
       write(nur,'(/i4,f14.8,2x,a)') jot1-1,E1,Label1
       write(nur,'( i4,f14.8,2x,a)') jot2-1,E2,Label2
       write(nur,'(a,2E12.5,f10.2)')  'S = ',SL,SL*SL,SL*SL*100
       Cycle
      end if

      dmatL = SL
      dmatV = SV/de*c_au

      SL=dmatL**2 
      SV=dmatV**2

      kp = kpol+kpol+1
      S = 1.d0;  Do i = 1,kp,2;  S = S * i;  End do
      S = 2*kp*(kpol+1)/(S*S)/kpol * (de/c_au)**kp
       
      GWL = S*SL; GFL = c_au**3/2*GWL/(de*de)
      GWV = S*SV; GFV = c_au**3/2*GWV/(de*de)

      g1 = jot1+1; g2 = jot2+1
 
      alfL = 2*dmatL*dmatL/(E2-E1)/(kp*g1); alfaL = alfaL + alfL   
      alfV = 2*dmatV*dmatV/(E2-E1)/(kp*g1); alfaV = alfaV + alfV

      if(E2.gt.E1) then
       FL = GFL/g1; WL = GWL/g2/time_au
       FV = GFV/g1; WV = GWV/g2/time_au
       write(nuf,'(/a,i2,f15.8,2x,a)') par1,jot1,E1,trim(Label1) 
       write(nuf,'( a,i2,f15.8,2x,a)') par2,jot2,E2,trim(Label2) 
      else
       FL = GFL/g2; WL = GWL/g1/time_au
       FV = GFV/g2; WV = GWV/g1/time_au
       write(nuf,'(/a,i2,f15.8,2x,a)') par2,jot2,E2,trim(Label2) 
       write(nuf,'( a,i2,f15.8,2x,a)') par1,jot1,E1,trim(Label1) 
      end if

      de = abs(E1-E2) * au_cm
      ANGS = 1.d+8 / de;  ANGSA = ANGS_AIR(ANGS)

      write(nuf,'(f11.2,a5,2(3x,f10.2,a10))') &
                DE,' CM-1',ANGS,' ANGS(VAC)',ANGSA,' ANGS(AIR)'

      AA='FIK= '
      if(GF.ne.'f') then; FL=GFL; FV=GFV; AA='GF = '; end if

      if(ialpha.eq.0) then
       write(nuf,'(1x,a1,i1,2x,a4,1PD12.5,3x,a5,D12.5,3x,a6,D12.5)') &
                   ktype,kpol,'S = ',SL,AA,FL,'AKI = ',WL
       if(SV.ne.0.d0) & 
       write(nuf,'(9x,1PD12.5,8x,D12.5,9x,D12.5)') SV,FV,WV
      else
       write(nuf,'(1x,a1,i1,2x,a4,1PD12.5,3x,a5,D12.5,3x,a6,D12.5,a8,d16.8,a8,0P2f10.5)') &
                   ktype,kpol,'S = ',SL,AA,FL,'AKI = ',WL, &
                   '   RME =',dmatL, '  alpha=',alfL,alfaL
       if(SV.ne.0.d0) & 
       write(nuf,'(9x,1PD12.5,8x,D12.5,9x,D12.5,8x,d16.8,8x,0P2f10.5)') &
                   SV,FV,WV, dmatV, alfV,alfaV 
      end if      

      End do; End do   !  over solutions

      End Subroutine Gen_zf
