!======================================================================
      Subroutine Read_pqbs(nu)
!======================================================================
!     read B-spline w.f. from bsw-file (unit nu)
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer, intent(in) :: nu
      Integer :: i,j,k,ikp,ikq,ka,l,n,m,itype,nsw,ksw,mw,kp,kq
      Character(5) :: elw
      Integer, external :: Ifind_orb_df 
      Real(8), allocatable :: tw(:),pw(:),qw(:)
      Real(8) :: ee

! ... read the written B-spline grid and check if it matches 

      rewind(nu)
      read(nu) itype,nsw,ksw
      BACKSPACE(nu)
      allocate(tw(nsw+ksw),pw(nsw),qw(nsw))
      read(nu) itype,nsw,ksw,tw,kp,kq
      ikp = ksw - kp
      ikq = ksw - kq
      k=1
      if(ksw.ne.ks) k=0
      if(nsw.ne.ns) k=0
      if(ksp.ne.kp) k=0
      if(ksq.ne.kq) k=0
      Do i=1,ns+ks; if(k.eq.0) Exit
       if(abs(t(i)-tw(i)).lt.1.d-12) Cycle; k=0; Exit
      End do    

! ... read radial functions and converte them if necessary 

    1 read(nu,end=2) elw,mw,ee
      pw=0.d0; read(nu) pw(1:mw)
      qw=0.d0; read(nu) qw(1:mw)
      Call EL_NLJK(elw,n,ka,l,j,i)
      m = Ifind_orb_df(n,ka,i)
      if(m.eq.0) go to 1    ! skip that orbital 

      if(k.eq.1) then
       mbs(m)=mw
       p(:,1,m) = pw 
       p(:,2,m) = qw 
       e(m)=ee
       write(log,'(a,a,a)') ebs(m),' - read from file ', trim(AF_inp)
      else
       Call Convert_pq(nsw-ikp,kp,tw(1+ikp),pw,nsp,ksp,p(1,1,m),pbsp,fpbs,1,0)
       Call Convert_pq(nsw-ikq,kq,tw(1+ikq),qw,nsq,ksq,p(1,2,m),qbsp,fqbs,1,0)
       mbs(m) = nsp-1
       e(m)=ee
       write(log,'(a,a,a,a)') ebs(m),' - read from file ', trim(AF_inp),' and converted'
      end if

      go to 1
    2 Close(nu)

      if(allocated(tw)) deallocate(tw,pw,qw)

      End subroutine Read_pqbs


!======================================================================
      Subroutine Write_pqbs
!======================================================================
!     read B-spline w.f. from bsw-file (unit nu)
!----------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Integer :: i

      AF_out = trim(name)//'.bsw'
      Call Read_apar(inp,'out',AF_out)
      Call Read_aarg('out',AF_out)

      open(nuw,file=AF_out,form='UNFORMATTED')

      rewind(nuw)
      write(nuw) grid_type,ns,ks,t(1:ns+ks),ksp,ksq
      Do i=1,nbf
       write(nuw) ebs(i),mbs(i),e(i)
       write(nuw) p(1:mbs(i),1,i)
       write(nuw) p(1:mbs(i),2,i)
      End do

      End Subroutine Write_pqbs

