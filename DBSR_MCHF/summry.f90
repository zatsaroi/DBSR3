!====================================================================== 
      Subroutine Summry
!====================================================================== 
!     The results of a calculation are summarized in "SUMMRY" file 
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer :: i,j
      Real(8) :: r1,r2,rm1,ar1,ar2,am1,am2,am3 
      Real(8), external :: quadr

      write(log,'(/80(''-'')/)')
      write(log,'(a,a)') 'ATOM:  ',ATOM 

      write(log,'(/a)') 'Convergence (latest difference):'
      write(log,'(a,T40,1Pd10.2)') 'Orbital diff. =', orb_diff
      write(log,'(a,T40,1Pd10.2)') 'SCF diff. =', scf_diff

! ... Compute Moments

      write(log,'(/2x,a,7x,a,8x,a,7x,a,7x,a,7x,a,7x,a,5x,a/)') &
           'nl','E(nl)','1/R','R','R**2','dmp','ns','max_R'

      Do i = 1, nbf 
       r1  = quadr(p(1,1,i),p(1,1,i), 1) 
       r2  = quadr(p(1,1,i),p(1,1,i), 2)
       rm1 = quadr(p(1,1,i),p(1,1,i),-1) 
       Call DCME(nbs(i),kbs(i),z,ar1,ar2,am1,am2,am3)
       am1 = z - ar1/r1   ! ???
       write(log,'(a5,f15.8,3f9.3,1PE12.2,i6,0Pf10.2)') &
        ebs(i), e(i), rm1, r1, r2, dpm(i), mbs(i), t(mbs(i)+ks)  
      End do 
      
      write(log,'(/a,T36,f25.16)') 'Average energy:', Etotal

      write(log,'(/a/)') 'Optimized states:'

      Call SortR(nlevels,elevel,ip_level)
      Do j = 1,nlevels; i = ip_level(j)
       write(log,'(a34,T36,f25.16)') labeln(i),elevel(i)
      End do

      End Subroutine Summry 


!======================================================================
      Subroutine Output_jfile
!======================================================================
!     output mixing coefficients in j-file (DBSR) format:
!----------------------------------------------------------------------
      Use dbsr_mchf
     
      Implicit none
      Integer :: i,j,ib,ip,nc,k
      Real(8) :: CM
      Character(mlab) :: LAB
    
      AF_j = trim(name)//'.j'
      Call Read_aarg('j',AF_j)
      open(nuj,file=AF_j)
    
      write(nuj,'(a,i8)') 'ncfg =',ncfg
      write(nuj,*)
      write(nuj,'(a,i8)') 'nsol =',nlevels
      write(nuj,*)
      write(nuj,'(a)') 'Solutions:'
      Do i=1,nlevels; ib=block(i); nc=Jncfg(ib); ip=ip_level(i)
        CM = 0.d0; k = 1
        Do j=ip+1,ip+nc
         if(abs(coefs(j)).lt.CM) Cycle
         k = j-ip; CM = abs(coefs(j))
        End do
        Call Get_cfg_jj(JTc1(ib)+k-1)
        Call Label_jj (mlab,Lab,0); Labeln(i) = Lab
        write(nuj,'(i8,2x,a)') i,trim(LAB)
        write(nuj,'(f16.8,3i8)') elevel(i), JJc(ib), JTc1(ib),JTc2(ib)
        write(nuj,'(6f12.8)') (coefs(j),j=ip+1,ip+nc)
      End do
    
      End Subroutine Output_jfile


