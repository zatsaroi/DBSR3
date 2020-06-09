!=====================================================================
      Subroutine Scf_mchf_all
!=====================================================================
!     Solve the DF equations in turns 
!---------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Real(8) :: hfm(ms,ms,nbf), hx(ms,ms,nbf), rhs(ms,nbf), et ! hx - ???
      Real(8) :: t1,t2
      Real(8), external :: quadr
      Integer :: ip, i,j, it

      Call Write_run_par(log)

      etotal = -1.d0;  Call DIAG (etotal);  et = etotal
      write(log,'(/a,F16.8)') 'Etotal = ',etotal

      if(.not.allocated(pp)) Allocate(pp(ns,2,nbf))
      hfm = 0.d0; hx = 0.d0; rhs = 0.d0; pp = 0.d0; dpm=0.d0

      Do it = 1,max_it

       write(log,'(/80(''-''))')
       write(log,'(a,i5)') 'Iteration',it
       write(log,'(80(''-''))')

       write(log,'(/3x,a,7x,a,10x,a,8x,a/)') 'nl', 'E(nl)','dpm', 'ns'

       ! ... set up of matrixes:

       Do ip = 1,nbf; i=iord(ip); if(i.eq.0) Cycle
        Call mchf_matrix (i,hfm(1,1,i),rhs(1,i),hx(1,1,i))
       End do

       ! .. diagonalize the hf matrix

       Call CPU_time(t1)

       Do ip = 1,nbf; i=iord(ip); if(i.eq.0) Cycle

        if(it <= 1 .or. newton == 0  .or. dpm(i) > 0.1 ) then
         Call solve_eiv (i,hfm(1,1,i),pp(1,1,i),rhs(1,i)) 
        else
         Call solve_nr  (i,hfm(1,1,i),pp(1,1,i),rhs(1,i),hx(1,1,i))
        end if

        dpm(i)=maxval(abs(p(1:ns,1,i)-pp(1:ns,1,i)))/maxval(abs(p(1:ns,1,i)))
        write(log,'(1x,a5,E14.5,E13.3,i7)') ebs(i), e(i), dpm(i), mbs(i)

        if(debug.gt.0) then
         Do j = 1,nbf 
          if(i.eq.j) Cycle 
          if(kbs(i).ne.kbs(j)) Cycle
          write(log,'(7x,a,a,a,a,a,E12.3)') &
          '<',ebs(i),'|',ebs(j),'> ',quadr(pp(1,1,i),p(1,1,j),0)
         End do
        end if

        if(abs(dpm(i)).lt.orb_tol) iord(i) = 0 

       End do ! over functions 

       Call CPU_time(t2)
       time_solve = time_solve + (t2-t1)

       Do ip = 1,nbf; i=iord(ip); if(i.eq.0) Cycle
        p(:,:,i) = pp(:,:,i);  Call Check_tails(i)
       End do

       Call DIAG(etotal)

       ! .. convergence test 

       orb_diff =  maxval(abs(dpm(1:nbf)))
       scf_diff =  abs(et-etotal)/abs(etotal)
       et = etotal

       write(log,*)
       write(log,'(a,F16.8)') 'Etotal = ',etotal
       write(log,*)
       write(log,'(A,T40,1P2D10.2)') 'Maximum orbital diff (relative) ', orb_diff
       write(log,'(A,T40,1P2D10.2)') 'Energy difference (relative) ', scf_diff
       

       write(scr,'(a,i5,f22.15,1P2d12.2)') 'it,etotal,scf_diff,orb_diff:', &
                                            it,etotal,scf_diff,orb_diff

       if ( orb_diff < orb_tol .and. scf_diff  < scf_tol ) Exit

       icore = 0;  if(ncore.gt.0) icore = sum(iord(:))
    
      End do ! over itterations

      End Subroutine Scf_mchf_all


