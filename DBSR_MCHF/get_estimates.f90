!======================================================================
      Subroutine get_estimates
!======================================================================
! ... Get initial estimates:
! ... (1) reading bsw-file
! ... (2) screened hydrogenic
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
      
      Implicit none
      Real(8) :: s,snl(nbf)
      Real(8), external :: quadr, Ecore_df
      Integer :: i, j

      write(log,'(/80(''-''))')
      write(log,'(/a/)') 'Initial estimations:'

! ... Read AF_inp for fixed orbitals or initial estimates:

      AF_inp = trim(name)//'.bsw'
      Call Read_apar(inp,'inp',AF_inp)
      Call Read_aarg('inp',AF_inp)

      if(Icheck_file(AF_inp).ne.0) then
       open(nuw,file=AF_inp,form='UNFORMATTED')
       Call Read_pqbs(nuw)
      end if

      snl = 0.d0
      Do i = 1, nbf
       if(mbs(i).ne.0) Cycle  ! We already have an initial estimate
       Call Screen_nl(i,snl(i))
      End do 
      
      Do i=1,nbf; if(mbs(i).ne.0) Cycle
                  if(snl(i).eq.0.d0) Cycle
       Do j=i+1,nbf; if(mbs(j).ne.0) Cycle
                     if(snl(j).eq.0.d0) Cycle
        if(nbs(i).ne.nbs(j)) Cycle
        if(lbs(i).ne.lbs(j)) Cycle
        S = (snl(i)+snl(j))/2
        snl(i) = s; snl(j) = s
       End do
      End do

      Do i = 1, nbf;   if(mbs(i).ne.0) Cycle  

       Call bdcwf_pq(nbs(i),kbs(i),z-snl(i),p(1,1,i),p(1,2,i))

! ... set orthogonality constraints:

       Do j = 1,i-1
        if (kbs(i) /= kbs(j)) Cycle
        S=quadr(p(1,1,i),p(1,1,j),0)
        if(abs(S).lt.orb_tol) Cycle
        p(:,:,i) = p(:,:,i) - S * p(:,:,j)
       End do
       S=quadr(p(1,1,i),p(1,1,i),0);  p(:,:,i) = p(:,:,i)/ sqrt(S)

       write(log,'(a,a,f5.2,a,i2)') ebs(i),&
       ' - hydrogenic orbital with screening ',snl(i)

      End do

      Call Check_tails(0)

      Call Gen_hd_core (ncore,mbreit,0)    

      Ecore = Ecore_df (ncore) 

      End Subroutine get_estimates


!======================================================================
      Subroutine Screen_nl(io,S)
!======================================================================
! ... find screenning parameter for orbital i
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer, intent(in) :: io
      Real(8), intent(out) :: S
      Integer :: i, ic,jc, ip
      Real(8) :: C

! ... core orbitals:

      if(io.le.ncore) then
       S = 0.d0
       Do i = 1, io-1
        S = S + qsum(i)
       End do
       S = S + qsum(io)/2      
       Return
      end if

! ... find leading configuration for this orbital:

      S = 0.d0
      Do ic = 1,ncfg; C = WC(ic)*WC(ic); if(C.le.S) Cycle
       Call Get_cfg_jj(ic)     
       ip = ip_state(ic)
       Do i=1,no; ip=ip+1
        if(IP_orb(ip).ne.io) Cycle
        ip = -1; Exit  
       End do
       if(ip.gt.0) Cycle;  S = C; jc = ic
      End do

! ... find screening papameter from given configuration:

      S = 0.d0
      if(ncore.gt.0) S = SUM(qsum(1:ncore))
      Call Get_cfg_jj (jc)     
      ip = ip_state(jc)
      Do i=1,no; ip=ip+1
       if(IP_orb(ip).ne.io) then
        S = S + iq(i)
       else
        S = S + iq(i)/2 
        Exit
       end if
      End do

      End Subroutine Screen_nl


