!====================================================================
      Subroutine Boundary_conditions
!====================================================================
! ... apply zero conditions at r=0 and on boundary r=a  which are 
! ... recorded in the array "iprm"
!--------------------------------------------------------------------
      Use dbsr_mchf
      Use df_orbitals
      
      Implicit none
      Integer :: i,j

      iprm=1
      Do i=1,nbf

       j=lbs(i)+1
       if(j.gt.ksp-1) j=1 
       if(nuclear.eq.'point') j=1
       if(ipzero.gt.0) j=ipzero 
       iprm(1:j,i)=0

       if(kbs(i).lt.0) j=lbs(i)+2
       if(kbs(i).ge.1) j=lbs(i)
       if(j.gt.ksq-1) j=1
       if(nuclear.eq.'point') j=1
       if(iqzero.gt.0) j=iqzero
       iprm(ns+1:ns+j,i)=0

       j=nsp-jpzero+1; iprm(j:ns,i)=0  
       j=nsq-jqzero+1; iprm(j+ns:ms,i)=0

      End do

      End Subroutine Boundary_conditions
       
