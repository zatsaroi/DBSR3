!======================================================================
      Subroutine Def_blocks
!======================================================================
! ... define the J-blocks and optimized levels
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Integer :: i,ib

      Call Def_Jblocks

      write(log,*)
      write(log,'(a,i3)') 'number of J-blocks: ',njbl
      write(log,*)
      Do i=1,njbl
       write(log,'(a,i2,a,i2,a,i5)') 'block ',i,'  2J= ',JJc(i),'   ncfg= ',Jncfg(i) 
      End do

      Call Read_ipar(inp,'eol',eol)
      Call Read_iarg('eol',eol)

      if(eol.eq.9) then
       Call Read_weights
      else
       Call Read_levels
      end if 

      write(log,*)
      if(eol.eq.1) &
      write(log,'(a)') 'optimization mode: eol = 1  -  equally weighted levels'
      if(eol.eq.5) &
      write(log,'(a)') 'optimization mode: eol = 5  -  statistically weighted levels'
      if(eol.eq.9) &
      write(log,'(a)') 'optimization mode: eol = 9  -  level"s weights are defined by user' 
      write(log,*)
      write(log,'(a,i3)') 'number of levels to optimize:',nlevels
      write(log,*)
      Do i=1,nlevels
       write(log,'(a,i2,a,i2,a,f10.5)') &
        'block ',block(i),'  level= ',level(i),'   weight= ',weight(i)
      End do

      i = maxval(Jncfg)
      memory_mat = i * i * 8.d0 /(1024*1024)
      write(log,'(/a,T40,F10.2,a)') 'Memory required for Hamiltonian matrix: ', &
                                     memory_mat,' Mb' 
 
! ... allocate expansion coefficients arrays:      
      
      Allocate(ip_level(nlevels)); ip_level = 0
      Do i=2,nlevels
       ip_level(i) = ip_level(i-1) + Jncfg(block(i-1))      ! bad assignment
      End do

      coefs_size = ip_level(nlevels)+Jncfg(block(nlevels))
      Allocate(coefs(coefs_size)) 
      Allocate(elevel(nlevels))

      Allocate(nlevelb (njbl));  nlevelb  = 0
      Do i=1,nlevels
       ib = block(i);  nlevelb(ib) = nlevelb(ib)+1
      End do

      Allocate(Labeln(nlevels))
       
      End Subroutine Def_blocks


!======================================================================
      Subroutine Read_weights
!======================================================================
!     read optimized levels with given weights (from the inp-file)
!     in format:  weight  =  block, level, weight
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Character(80) :: line
      Integer :: i,j
      Real(8) :: S 

      nlevels = 0

      Call Read_weights_arg; if(nlevels.ne.0) Return

      rewind(inp)
    1 read(inp,'(a)',end=2) line
      if(line(1:5).ne.'weight') go to 1
      nlevels=nlevels+1
      go to 1
    2 Continue
      if(nlevels.eq.0) &
       Stop 'Stop in Read_weights: no weights given for eol=9'

      Allocate(block(nlevels),level(nlevels),weight(nlevels))
      j = 0
      rewind(inp)
    3 read(inp,'(a)',end=4) line
      if(line(1:5).ne.'weight') go to 3
      j=j+1
      i=INDEX(line,'=')+1
      read(line(i:),*) block(j),level(j),weight(j)
      go to 3
    4 Continue

      S = SUM(weight); weight = weight / S

      End Subroutine Read_weights 


!======================================================================
      Subroutine Read_weights_arg
!======================================================================
!     read optimized levels with given weights (from the inp-file)
!     in format:  weight=block,level,weight;block,level,weight; ...
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Character(280) :: line 
      Integer :: j,i1,i2
      Real(8) :: S 

      line = ' '; Call Read_aarg('weights',line)

      if(len_trim(line).eq.0) Return

      i1=1 
      Do
       i2=INDEX(line(i1:),';')+i1
       if(i2.le.i1) i2=len_trim(line)
       if(i2.le.i1) Exit
       nlevels=nlevels + 1
       i1 = i2
      End do

      Allocate(block(nlevels),level(nlevels),weight(nlevels))

      j=0; i1=1 
      Do
       i2=INDEX(line(i1:),';')+i1-2
       if(i2.le.i1) i2=len_trim(line)
       if(i2.le.i1) Exit
       j = j + 1
       read(line(i1:),*) block(j),level(j),weight(j)
       i1 = i2+2
      End do

      S = SUM(weight); weight = weight / S

      End Subroutine Read_weights_arg 


!======================================================================
      Subroutine Read_levels
!======================================================================
!     read optimized levels (from the inp-file)  in the format
!     levels =  block, list of levels
!     Corresponding levels will have non-zero weights 
!     = 1.d0      for   eol =  1     (equally weighted)
!     = 2J+1      for   eol =  5     (statistically weighted)
!
!     If no such information was found, by defaults,  the first
!     ACS in each block will be used for optimization
!----------------------------------------------------------------------
      Use dbsr_mchf

      Implicit none
      Character(280) :: line = ' '
      Integer :: i,j, i1,i2, ib,ip, na, iarr(ncfg+1)
      Real(8) :: S 

      WC = 0.d0

      Call Read_aarg('levels',line)

      if(len_trim(line).ne.0) then

       i1=1 
       Do
        i2=INDEX(line(i1:),';')+i1-2
        if(i2.le.i1) i2=len_trim(line)
        if(i2.le.i1) Exit

        Call Read_iarr_string(line(i1:i2),na,iarr)
        if(na.le.1) Exit
   
        ib=iarr(1); ip = JTc1(ib)-1 
        Do i=2,na; j = ip + iarr(i); nlevels = nlevels + 1
         WC(j) = JJc(ib)+1; if(eol.eq.1) WC(j)=1.d0 
        End do
        i1 = i2+2
       End do

      else

         rewind(inp)
       1 read(inp,'(a)',end=2) line
         if(line(1:6).ne.'levels') go to 1
         i=INDEX(line,'=')+1
         Call Read_iarr_string(line(i:),na,iarr)
    
         ib=iarr(1); ip = JTc1(ib)-1 
         Do i=2,na; j = ip + iarr(i); nlevels = nlevels + 1
          WC(j) = JJc(ib)+1; if(eol.eq.1) WC(j)=1.d0 
         End do
         go to 1
       2 Continue

      end if

      if(nlevels.eq.0) then 
       Call Read_ipar(inp,'nlevels',nlevels)
       Call Read_ipar(inp,'nlevels',nlevels)
      end if 

      if(nlevels.eq.0) then
       nlevels = ncfg
       Do ib=1,njbl
        Do i=JTc1(ib),JTc2(ib)
         WC(i) = JJc(ib)+1; if(eol.eq.1) WC(i)=1.d0 
        End do
       End do
      end if     

      if(nlevels.eq.-1) then
       nlevels = njbl
       Do ib=1,njbl; j=JTc1(ib)
        WC(j) = JJc(ib)+1; if(eol.eq.1) WC(j)=1.d0 
       End do
      end if     

      Allocate(block(nlevels),level(nlevels),weight(nlevels))
      j = 0
      Do ib=1,njbl
       Do i=JTc1(ib),JTc2(ib)
        if(WC(i).eq.0.d0) Cycle      
        j = j + 1
        block(j) = ib
        level(j) = i - JTc1(ib) + 1
        weight(j) = WC(i)
       End do 
      End do

      S = SUM(weight); weight = weight / S

      End Subroutine Read_levels 


