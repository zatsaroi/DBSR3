!=====================================================================
      Subroutine Scf_mchf
!=====================================================================
!     Solve the DF equations in turns 
!---------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals

      Implicit none
      Real(8) :: hfm(ms,ms), v(ms), et, rhs(ms), hx(ms,ms)
      Real(8) :: t1,t2
      Real(8), external :: quadr
      Integer :: ip, i,j, it

      Call Write_run_par(log)

      Call diag(etotal);    et = etotal

      dpm = 0.d0

      Do it = 1,max_it

       write(log,'(/80(''-''))')
       write(log,'(a,i5)') 'Iteration',it
       write(log,'(80(''-''))')

       write(log,'(/3x,a,7x,a,10x,a,11x,a,11x,a/)') 'nl','e(nl)','qsum', 'dpm', 'ns'

       Do ip = 1,nbf; i=iord(ip); if(i.eq.0) Cycle

        Do j=i+1,nbf
         if(rotate.eq.0) Cycle; if(it <= 1) Cycle;  Call Rotate_ij(i,j)
        End do
   
        Call mchf_matrix (i,hfm,rhs,hx)

        ! .. diagonalize the hf matrix

        Call CPU_time(t1)
        if(it <= 1 .or. newton == 0  .or. dpm(i) > 0.1 ) then
         Call solve_eiv (i,hfm,v,rhs) 
        else
         Call solve_nr  (i,hfm,v,rhs,hx)
        end if
        Call CPU_time(t2)
        time_solve = time_solve + (t2-t1)

        dpm(i)=maxval(abs(p(1:ns,1,i)-v(1:ns)))/maxval(abs(p(:,1,i)))
        p(1:ns,1,i) = v(1:ns)
        p(1:ns,2,i) = v(ns+1:ms)       

        ! .. remove tail zero

        Call Check_tails(i)

        write(log,'(1x,a5,3e15.5,i7)') &
                   ebs(i),  e(i), qsum(i),dpm(i), mbs(i)

       End do ! over functions 


       Do i = 1,nbf 
        Do j = 1,nbf 
         if(i.eq.j) Cycle 
         if(kbs(i).ne.kbs(j)) Cycle
         if(debug.gt.0) write(log,'(7x,a,a,a,a,a,E12.3)') &
         '<',ebs(i),'|',ebs(j),'> ',quadr(p(1,1,i),p(1,1,j),0)
        End do
       End do

       Call diag(etotal)

       ! .. test convergence

       orb_diff =  maxval(abs(dpm(1:nbf)))
       scf_diff =  abs(et-etotal)/abs(etotal)

       write(log,*)
       write(log,'(a,F22.8)') 'etotal = ',etotal
       write(log,*)
       write(log,'(A,T40,1P2d10.2)') 'Maximum orbital diff (relative): ', orb_diff,orb_tol
       write(log,'(A,T40,1P2d10.2)') 'Energy difference (relative): ', scf_diff,scf_tol

       write(scr,'(a,i5,f22.15,1P4d12.2)') 'it,etotal,scf_diff,orb_diff:', &
                                            it,etotal,scf_diff,orb_diff
        
       et = etotal

       if ( orb_diff < orb_tol .and. scf_diff  < scf_tol) Exit
    
       icore = 0;  if(ncore.gt.0) icore = sum(iord(:))

      End do ! over itterations

      End subroutine Scf_mchf



