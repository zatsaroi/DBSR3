!======================================================================
      Subroutine read_int_corr
!======================================================================
!     read semi-empirical corrections for integrals
!----------------------------------------------------------------------
      Use dbsr_ci
      Implicit none
      Integer :: i,n,kappa,l,j,k
      Character :: AS*80
      Integer, external :: Ifind_position, Ifind_jjorb

      ncorr = 0
      if(inp.le.0) Return
      if(Ifind_position(inp,'integrals corrections:').eq.0) Return

    1 read(inp,'(a)',end=2) AS
      if(AS(1:1).eq.'*') go to 2
      ncorr=ncorr+1
      go to 1
    2 if(ncorr.eq.0) Return

      Allocate(icorr(6,ncorr), scorr(ncorr))
      i = Ifind_position(inp,'integrals corrections:')
      i = 0
      rewind(inp)
    3 read(inp,'(a)',end=4) AS
      if(AS(1:1).eq.'*') go to 4
      if(AS(1:1).eq.'R') then
       i = i + 1
       icorr(1,i) = 2      
       read(AS(2:3),*) icorr(2,i)
       Call EL_NLJK(AS(6:10),n,kappa,l,j,k)
       icorr(3,i) = Ifind_jjorb(n,kappa,k,0)
       Call EL_NLJK(AS(13:17),n,kappa,l,j,k)
       icorr(4,i) = Ifind_jjorb(n,kappa,k,0)
       Call EL_NLJK(AS(20:24),n,kappa,l,j,k)
       icorr(5,i) = Ifind_jjorb(n,kappa,k,0)
       Call EL_NLJK(AS(27:31),n,kappa,l,j,k)
       icorr(6,i) = Ifind_jjorb(n,kappa,k,0)
       read(AS(34:),*) scorr(i)
      end if
      go to 3
    4 Continue

      End Subroutine read_int_corr


!======================================================================
      Real(8) Function S_corr(met,k,i1,i2,i3,i4)
!======================================================================
!     apply the semi-empirical corrections for integrals
!----------------------------------------------------------------------
      Use dbsr_ci
      Implicit none
      Integer :: i,met,k,i1,i2,i3,i4

      S_corr = 1.d0
      if(ncorr.eq.0) Return

      Do i=1,ncorr
       if(icorr(1,i).ne.met) Cycle
       if(icorr(2,i).ne.k) Cycle
       if(icorr(3,i).ne.i1) Cycle
       if(icorr(4,i).ne.i2) Cycle
       if(icorr(5,i).ne.i3) Cycle
       if(icorr(6,i).ne.i4) Cycle
       S_corr = scorr(i)
       Exit
      End do

      End Function S_corr
