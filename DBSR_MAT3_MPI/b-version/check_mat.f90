!======================================================================
      Subroutine Check_mat(met)
!======================================================================
!     check overlap matrix for existence of big elements ( > S_ovl)
!----------------------------------------------------------------------
      Use dbsr_mat
      Use phys_orb_jj

      Implicit none
      Real(8) :: S, SM, C, CM, v(ms)
      Integer :: met, i,j,  j1,j2, jj, ich,jch, it,jt, js, ip, &
                 ksol, jm, k
      Integer, external :: IBORT

      met = 0
      if(SUM(overlaps).eq.0.d0) Return

      Close(nuc);  Open(nuc,file=AF_cfg,position='APPEND'); write(nuc,*)

! ... channel-channel blocks:

      Do ich=2,nch
       Do jch=1,ich-1
!        if(icc(ich,jch).eq.0) Cycle
        if(overlaps(ich,jch).lt.s_ovl) Cycle 
        met = met + 1
        ! ... find bigest overlap and put orth.condition
        Call Def_SM(ich,jch)       
        if(SM.eq.0.d0) Call Def_SM(jch,ich)
        if(SM.eq.0.d0) then
         write(pri,'(2a10,f10.5)') elc(ich),elc(jch),overlaps(ich,jch)
         Stop 'I cannot find orth.condition'
        end if
       End do  ! over jch
       if(met.ne.0) Exit     
      End do

!----------------------------------------------------------------------
      if(npert.eq.0.or.met.ne.0) Return

! ... channel-perturber blocks:

      Do ich=1,nch

       Do ip=1,npert
!        if(icb(ich,ip).eq.0) Cycle
        if(overlaps(nch+ip,ich).lt.s_ovl) Cycle 
        met = met + 1 

        ! ... find bigest overlap and put orth.condition

        SM = 0.d0; jm = 0
        Do jj=1,npert_sub; j=np_phy(jj); js=np_sub(jj)
         if(kch(ich).ne.kbs(js)) Cycle
         if(IBORT(ipch(ich),js).eq.0) Cycle
         Call Get_qv(js,v,ns)
          ksol=ipsol(ich); CM = 0.d0 
          Do i = 1,ksol
           C = abs( SUM(v(:)*diag(:,i,ich)) )        
           if(C.lt.CM) Cycle; CM = C
          End do
          if(CM.lt.SM) Cycle
          SM = CM; jm = js
         End do                                                               
 
         if(SM.eq.0.d0) Stop 'I cannot find orth.condition'

         it = iptar(ich)
         write(nuc,'(a1,a5,a1,a5,a3,5x,a12,a6,5x,a12,i6,3f10.3)') &
            '<',elc(ich),'|',ebs(jm),'>=0', AFT(it),elc(ich), 'perturber   ',ip,SM, &
            overlaps(nch+ip,ich), s_ovl 
         Call Nulify_BORT(ipch(ich),jm)  
         overlaps(:,ich) = 0.d0

       End do ! over perturbers

      End do ! over ich

! ... perturber-perturber blocks:

      if(npert.le.1) Return

      Do i=2,npert; Do j=1,i-1                   !; k=ibb(i,j); if(k.eq.0) Cycle
       S = abs(overlaps(i+nch,j+nch))
       if(S.lt.S_ovl) Cycle
       write(nuc,'(f10.5,2i5,a)') S, i,j , ' - suspicious perturber overlap '
      End do; End do

 CONTAINS

!======================================================================
      Subroutine Def_SM(ich,jch)
!======================================================================
      Implicit none
      Integer :: ich,jch

      jt = iptar(jch)
      j1=1; if(jt.gt.1) j1=ip_tar(jt-1)+1; j2=ip_tar(jt)
      SM = 0.d0; jm = 0
      Do jj=j1,j2; j=ip_phy(jj); js=ip_sub(jj)
       if(kch(ich).ne.kbs(j)) Cycle
       if(IBORT(ipch(ich),js).eq.0) Cycle
       Call Get_qv(js,v,ns)
       ksol=ipsol(ich); CM = 0.d0 
       Do i = 1,ksol
        C = abs( SUM( v(1:ms)*diag(1:ms,i,ich) ) )        
        if(C.lt.CM) Cycle; CM = C
       End do
       if(CM.lt.SM) Cycle
       SM = CM; jm = js
      End do                                                               

      it = iptar(ich)
      write(nuc,'(a1,a5,a1,a5,a3,5x,a12,a6,5x,a12,a6,3f10.3)') &
         '<',elc(ich),'|',ebs(jm),'>=0', AFT(it),elc(ich), AFT(jt),elc(jch), SM, &
         overlaps(ich,jch), s_ovl 
      Call Nulify_bort(ipch(ich),jm) 
      overlaps(:,jch) = 0.d0

      End Subroutine Def_SM


      End Subroutine Check_mat


