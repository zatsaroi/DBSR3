!======================================================================
!     PROGRAM  test_coef_2conf_jj                      
!
!               C O P Y R I G H T -- 2015
!
!     Written by:   Oleg Zatsarinny
!                   email: oleg_zoi@yahoo.com
!======================================================================
!
!     it is a debug program to check subroutine   "coef_2conf_jj" 
!
!----------------------------------------------------------------------
!
!     INPUT ARGUMENTS:  name  or   cfg=...  tab=...  
!    
!----------------------------------------------------------------------
!
!     INPUT FILE:     AF_cfg    (default  - rcsl.inp) 
!     OUTPUT FILES:   AF_tab    (default  - coef.tab)
!    
!----------------------------------------------------------------------     
      Use conf_jj
      Use rk4_data

      Implicit none 

      Integer :: nuc=1; Character(40) :: AF_cfg = 'rcsl.inp'
      Integer :: out=2; Character(40) :: AF_tab = 'coef.tab'
      Character(80) :: name = ' ', line, int

      Character(5) :: EL1, EL2, EL3, EL4, ELj
      Real(8) :: C, eps_C = 1.d-7, S(8), SMU
      Integer :: i,j, k,v, k1,k2,k3,k4, ic,jc, i1,i2,j1,j2, m
      
! ... result coefficients:

      Integer, parameter :: mcoef = 1000, ibi = 2**16
      Integer :: ncoef
      Real(8) :: coefs(mcoef)
      Integer :: icoefs(5,mcoef)

!----------------------------------------------------------------------
      Call Read_name(name)
      if(len_trim(name).ne.0) then
       AF_cfg=trim(name)//'.c'
       AF_tab=trim(name)//'.tab'
      else
       Call Read_aarg('cfg',AF_cfg)
       Call Read_aarg('tab',AF_tab)
      end if

      Call Check_file(AF_cfg)
      Open(nuc,file=AF_cfg)
      Call Read_conf_jj(nuc,0,'detect','check')

      Open(out,file=AF_tab)

      Do ic = 1,ncfg
       Call Get_cfg_jj(ic)
       Call Save_cfg_jj(1)

      Do jc = ic,ncfg
       Call Get_cfg_jj(jc)
       Call Save_cfg_jj(2)

       Call coef_2conf_jj(no1,nn1,ln1,jn1,iq1,Jshell1,Vshell1,Jintra1,  &
                          no2,nn2,ln2,jn2,iq2,Jshell2,Vshell2,Jintra2,  &
                          mcoef,ncoef,icoefs,coefs)

      if(ncoef.eq.0) Cycle

      write(out,'(///70(''=''))')
      write(out,'(12x,''< state'',i2,'' || (H-E) || state'',i2,''>'')') &
        ic,jc
      write(out,'(70(''=''))') 

       write(out,'(72("-"))')
       Call print_conf_jj (out,ic,0.d0)
       if(ic.ne.jc) Call print_conf_jj (out,jc,0.d0)
       write(out,'(72("-"))')

! ... one-electron integrals

      m = 0
      Do i=1,ncoef
       k =icoefs(1,i) 
       if(k.ne.-1) Cycle
       C = coefs(i); if(abs(C).lt.eps_c) Cycle

       if(m.eq.0) then
        write(out,'(/12x,'' One-electron integrals:''/)')
        m=1
       end if

       i1=icoefs(2,i); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0)
       j1=icoefs(4,i); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0)

       Call Num(C,i1,i2,999999,1.d-9);  i1=iabs(i1)
       write(line,'(a,i6,a,i6,a)') '[',i1,':',i2,']'
       j1=sqrt(1.*i1)+0.1; j2=sqrt(1.*i2)+0.1
       if(i1.eq.j1*j1.and.i2.eq.j2*j2) &
       write(line,'(a,i6,a,i6,a)') '(',j1,':',j2,')'

       write(out,'(f12.6,2x,a,4x,a,a,a,a,a,a)') &
                   C,trim(line),'L','(',EL1,';',EL3,')'        

      End do

      write(out,'(/12x,'' two-electron integrals:''/)')

      Do i=1,ncoef
       C = coefs(i); if(abs(C).lt.eps_c) Cycle
       k =icoefs(1,i) 
       if(k.lt.0) Cycle

       i1=icoefs(2,i); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0)
       i2=icoefs(3,i); EL2 = ELj(nn1(i2),ln1(i2),jn1(i2),0)
       j1=icoefs(4,i); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0)
       j2=icoefs(5,i); EL4 = ELj(nn2(j2),ln2(j2),jn2(j2),0)

       if(mod(k+ln1(i1)+ln2(j1),2).ne.0) Cycle
       if(mod(k+ln1(i2)+ln2(j2),2).ne.0) Cycle


       if(ic.eq.jc.and.k.ge.0) then
        if(j1.le.j2) write(int,'(a,i2,a,a,a,a,a,f10.5)') &
                               'F',k,'(',EL1,',',EL2,')'        
        if(j1.gt.j2) write(int,'(a,i2,a,a,a,a,a,f10.5)') &
                               'G',k,'(',EL1,',',EL2,')'        
       elseif(k.ge.0) then
        write(int,'(a,i2,a,a,a,a,a,a,a,a,a,f10.5)') &
                  'R',k,'(',EL1,',',EL2,';',EL3,',',EL4,')'        
       end if


       Call Num(C,i1,i2,999999,1.d-9);  i1=iabs(i1)
       write(line,'(a,i6,a,i6,a)') '[',i1,':',i2,']'
       j1=sqrt(1.*i1)+0.1; j2=sqrt(1.*i2)+0.1
       if(i1.eq.j1*j1.and.i2.eq.j2*j2) &
       write(line,'(a,i6,a,i6,a)') '(',j1,':',j2,')'

       write(out,'(f12.6,2x,a,4x,a)') C,trim(line),trim(int)        

      End do

!--------------------------------------------------------------------------
      write(out,'(/12x,'' Breit adding:''/)')

      Do i=1,ncoef
       C = coefs(i); if(abs(C).lt.eps_c) Cycle
       k =icoefs(1,i) 
       if(k.lt.0) Cycle

       i1=icoefs(2,i); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0)
       i2=icoefs(3,i); EL2 = ELj(nn1(i2),ln1(i2),jn1(i2),0)
       j1=icoefs(4,i); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0)
       j2=icoefs(5,i); EL4 = ELj(nn2(j2),ln2(j2),jn2(j2),0)

       if(mod(k+ln1(i1)+ln2(j1),2).eq.0.and. &
          mod(k+ln1(i2)+ln2(j2),2).eq.0) Cycle

       if(ic.eq.jc.and.k.ge.0) then
        if(j1.le.j2) write(int,'(a,i2,a,a,a,a,a,f10.5)') &
                               'F',k,'(',EL1,',',EL2,')'        
        if(j1.gt.j2) write(int,'(a,i2,a,a,a,a,a,f10.5)') &
                               'G',k,'(',EL1,',',EL2,')'        
       elseif(k.ge.0) then
        write(int,'(a,i2,a,a,a,a,a,a,a,a,a,f10.5)') &
                  'R',k,'(',EL1,',',EL2,';',EL3,',',EL4,')'        
       end if


       Call Num(C,i1,i2,999999,1.d-9);  i1=iabs(i1)
       write(line,'(a,i6,a,i6,a)') '[',i1,':',i2,']'
       j1=sqrt(1.*i1)+0.1; j2=sqrt(1.*i2)+0.1
       if(i1.eq.j1*j1.and.i2.eq.j2*j2) &
       write(line,'(a,i6,a,i6,a)') '(',j1,':',j2,')'

       write(out,'(f12.6,2x,a,4x,a)') C,trim(line),trim(int)        

      End do

!---------------------------------------------------------------------------
! ... generate  Sk integrals:

      Do i=1,ncoef
       C = coefs(i); if(abs(C).lt.eps_c) Cycle
       k =icoefs(1,i) 
       if(k.lt.0) Cycle

       i1=icoefs(2,i); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0); k1=kn1(i1)
       i2=icoefs(3,i); EL2 = ELj(nn1(i2),ln1(i2),jn1(i2),0); k2=kn1(i2)
       j1=icoefs(4,i); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0); k3=kn2(j1)
       j2=icoefs(5,i); EL4 = ELj(nn2(j2),ln2(j2),jn2(j2),0); k4=kn2(j2)

!       if(mod(k+ln1(i1)+ln2(j1),2).ne.0) Cycle
!       if(mod(k+ln1(i2)+ln2(j2),2).ne.0) Cycle

       Do v = k-1,k+1
         if(v.lt.0) Cycle
         if(SMU(k1,k2,k3,k4,k,v,S).eq.0.d0) Cycle
         S = S * C
         Call Add_rk4_data(1,v,i1+ibi*i2,j1+ibi*j2,S(1))
         Call Add_rk4_data(1,v,i2+ibi*i1,j2+ibi*j1,S(2))
         Call Add_rk4_data(1,v,j1+ibi*j2,i1+ibi*i2,S(3))
         Call Add_rk4_data(1,v,j2+ibi*j1,i2+ibi*i1,S(4))
         Call Add_rk4_data(1,v,i1+ibi*j2,j1+ibi*i2,S(5))
         Call Add_rk4_data(1,v,j2+ibi*i1,i2+ibi*j1,S(6))
         Call Add_rk4_data(1,v,j1+ibi*i2,i1+ibi*j2,S(7))
         Call Add_rk4_data(1,v,i2+ibi*j1,j2+ibi*i1,S(8))
       End do

      End do

!--------------------------------------------------------------------------
      write(out,'(/12x,'' Sk integrals (ppqq):''/)')

      Do i=1,nrk
       
       C = crk(i); if(abs(C).lt.eps_c) Cycle
                   if(kr1(i).ne.1) Cycle
       k = kr2(i) 
       i1 = mod(kr3(i),ibi); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0)
       i2 = kr3(i)/ibi;      EL2 = ELj(nn1(i2),ln1(i2),jn1(i2),0)
       j1 = mod(kr4(i),ibi); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0)
       j2 = kr4(i)/ibi;      EL4 = ELj(nn2(j2),ln2(j2),jn2(j2),0)

       write(int,'(a,i2,a,a,a,a,a,a,a,a,a,f10.5)') &
                  'S1',k,'(',EL1,',',EL2,';',EL3,',',EL4,')'        

       Call Num(C,i1,i2,999999,1.d-9);  i1=iabs(i1)
       write(line,'(a,i6,a,i6,a)') '[',i1,':',i2,']'
       j1=sqrt(1.*i1)+0.1; j2=sqrt(1.*i2)+0.1
       if(i1.eq.j1*j1.and.i2.eq.j2*j2) &
       write(line,'(a,i6,a,i6,a)') '(',j1,':',j2,')'

       write(out,'(f12.6,2x,a,4x,a)') C,trim(line),trim(int)        
      End do

      write(out,'(/12x,'' Sk integrals (pqqp):''/)')

      Do i=1,nrk
       
       C = crk(i); if(abs(C).lt.eps_c) Cycle
                   if(kr1(i).ne.2) Cycle
       k = kr2(i) 
       i1 = mod(kr3(i),ibi); EL1 = ELj(nn1(i1),ln1(i1),jn1(i1),0)
       i2 = kr3(i)/ibi;      EL2 = ELj(nn1(i2),ln1(i2),jn1(i2),0)
       j1 = mod(kr4(i),ibi); EL3 = ELj(nn2(j1),ln2(j1),jn2(j1),0)
       j2 = kr4(i)/ibi;      EL4 = ELj(nn2(j2),ln2(j2),jn2(j2),0)
       write(int,'(a,i2,a,a,a,a,a,a,a,a,a,f10.5)') &
                  'S2',k,'(',EL1,',',EL2,';',EL3,',',EL4,')'        

       Call Num(C,i1,i2,999999,1.d-9);  i1=iabs(i1)
       write(line,'(a,i6,a,i6,a)') '[',i1,':',i2,']'
       j1=sqrt(1.*i1)+0.1; j2=sqrt(1.*i2)+0.1
       if(i1.eq.j1*j1.and.i2.eq.j2*j2) &
       write(line,'(a,i6,a,i6,a)') '(',j1,':',j2,')'

       write(out,'(f12.6,2x,a,4x,a)') C,trim(line),trim(int)        
      End do


      End do; End do   ! over ic

      End  ! PROGRAM  test_coef_2conf_jj

