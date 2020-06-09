!======================================================================
      Subroutine cm2j(num,nuj,nblock,JBtot,JBsol)
!======================================================================
!     Convert name.mix or name.cm files from rscf2 or rci2
!     to j-file of DBSR format
!
!     arguments:  file units for cm- and j-files
!
!     list of configurations are placed in module conf_jj
!----------------------------------------------------------------------
      Use conf_jj
     
      Implicit real(8) (A-H,O-Z)

      Character(6) :: G92MIX

      Integer :: num,nuj,nblock,JBtot(*),JBsol(*)
      Real(8), allocatable :: EVAL(:), EVEC(:,:)
      Integer, allocatable :: ICCMIN(:)

      write(nuj,'(a,i8)') 'ncfg =',ncfg
      nsol=0

! ... read solution:

      read(num) G92MIX
      if(G92MIX.ne.'G92MIX') Stop 'cm2j: input file is not GRASP2K mix-file'

      read(num) NELEC, NCFTOT, NW, NVECTOT, NVECSIZ, NBLOCK
      if(NCFTOT.ne.ncfg) Stop 'ncfg in c-file <> NCFTOT in cm-file'

      write(nuj,*)
      write(nuj,'(a,i8)') 'nsol =',NVECTOT

      write(nuj,*)
      write(nuj,'(a)') 'Solutions:'

      ic2=0
      Do ib = 1,NBLOCK
       read(num) NB, NCFBLK, NEVBLK, IATJP, IASPA 
       ic1=ic2+1; ic2=ic2+NCFBLK

       JBtot(ib) = IATJP-1
       JBsol(ib) = NEVBLK

       if(allocated(ICCMIN)) Deallocate(ICCMIN) 
       Allocate(ICCMIN(NEVBLK))
       read(num) (ICCMIN(I),I=1,NEVBLK)

       if(allocated(EVAL)) Deallocate(EVAL) 
       Allocate(EVAL(NEVBLK))
       read(num) EAV, (EVAL(i), I = 1, NEVBLK)

       if(allocated(EVEC)) Deallocate(EVEC) 
       Allocate(EVEC(NCFBLK,NEVBLK))
       read(num) ((EVEC(i,j), i=1,NCFBLK), J=1,NEVBLK)

       Do k=1,NEVBLK
        nsol = nsol+1
        write(nuj,'(i8,2x,2F16.8)') nsol
        write(nuj,'(f16.8,3i8)') EVAL(k)+EAV, IATJP-1, ic1,ic2       
        write(nuj,'(6f12.8)') (EVEC(i,k),i=1,NCFBLK)
       End do

      End do

      End Subroutine cm2j
