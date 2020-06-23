!====================================================================
!     Utility    jj_list 
!====================================================================
!
!     creates list of states in set of j-files  
!     given in jj_list.inp or one j-file as argument
!
!     name.c + name.j -> name.list
!
!     additional parameters:
!
!     mconf  -   max.number of leading configurations to show [1]
!     eps_c  -   tolerance for configuration weights [0.2]           
!     unit   -   au, Ry, eV, cm  [au]                      
!     shift  -   overall energy shift [0.d0]'                  
!     msol   -   max.number of solutions [0 -> all]'         
!
!     Example:
!
!       jj_list  name.j  mconf=2  eps_c=0.2  unit=eV  shift=1.d0      
!
!     All arguments except first - optional
!
!--------------------------------------------------------------------
      Use conf_jj

      Implicit real(8) (A-H,O-Z)
      Integer :: mconf = 0
      Integer :: iev   = 0
      Integer :: msol  = 0
      Real(8) :: eps_c = 0.0d0      
      Real(8) :: shift = 0.d0
      Real(8) :: Z     = 100.d0

      Character(1) :: ctype='j'
      Character(2) :: unit='au'

      Real(8), allocatable :: EE(:), C(:)
      Real(8), allocatable :: CC(:,:)

      Integer, allocatable :: IP(:), IPE(:), JJ(:), kco(:), isol(:)
      Integer, allocatable :: IPC(:,:)

      Integer, parameter :: ma=64 

      Character(ma) :: Label1
      Character(ma) :: AF

      Character(ma), Allocatable :: LAB(:,:)
      Character(ma), Allocatable :: AFF(:)

      Integer :: inp= 5; Character(ma) :: AFi = 'jj_list.inp'
      Integer :: nuc=11; Character(ma) :: AFc = ' '
      Integer :: nuj=12; Character(ma) :: AFj = ' '
      Integer :: nur=13; Character(ma) :: AFr = 'jj_list'
      Integer :: nua=21; ! scratch file

!----------------------------------------------------------------------
! ... input parameters:

      nfile = 0
      if(Icheck_file(AFi).eq.1) then
       Open(inp,file=AFi)
       Call Read_rpar(inp,'eps_c',eps_c)
       Call Read_rpar(inp,'shift',shift)
       Call Read_ipar(inp,'mconf',mconf)
       Call Read_ipar(inp,'msol' ,msol )
       Call Read_apar(inp,'unit' ,unit )
       Call Read_rpar(inp,'Z'    ,Z    )
       Call Read_rpar(inp,'awt'  ,awt  )
       Call Read_apar(inp,'out'  ,AFr  )
       rewind(inp)
      end if

      Call Read_rarg('eps_c',eps_c)
      Call Read_rarg('shift',shift)
      Call Read_iarg('mconf',mconf)
      Call Read_iarg('msol' ,msol)
      Call Read_aarg('unit' ,unit)
      Call Read_rarg('Z'    ,Z    )
      Call Read_rarg('awt'  ,awt  )

      Call Read_aarg('inp' ,AFj)
      if(LEN_TRIM(AFj).ne.0) then
       nfile=1
       i=LEN_TRIM(AFj)
       AFr=AFj(1:i-1)//'list'
      end if

! ... find conversion factor:
       
      AWT = 0.d0  
      Call Conv_au (Z,AWT,au_cm,au_eV,0)

      Select case(unit)
       case('au'); conv = 1.d0;  iev=0
       case('Ry'); conv = 2.d0;  iev=1
       case('eV'); conv = au_eV; iev=2
       case('cm'); conv = au_cm; iev=3
       case default; conv = 1.d0; iev = 0
      End select

      Open(nua,status='SCRATCH',form='UNFORMATTED') 

!----------------------------------------------------------------------
! ... cycle over j-files:

      iaf = 0

      nsol=0
    1 Continue 


      AF = AFj
      if(nfile.ne.1) then
       read(inp,'(a)',end=2) AF
       if(AF(1:1).eq.'*') go to 2
       i=INDEX(AF,'='); if(i.ne.0) go to 1
      end if

      j = LEN_TRIM(AF); if(j.gt.iaf) iaf=j

      AFj=AF; open(nuj,file=AFj,status='old')
      i=LEN_TRIM(AFj);  AFc=AFj(1:i-1)//'c'
      open(nuc,file=AFc,status='old')

! ... find the labels for config.s:

      ncfg=0; Call Read_ipar(nuj,'ncfg',ncfg) 
      if(ncfg.eq.0) Stop ' ncfg = 0'

      Call R_label_jj(nuc,0)      

      ilab = 1
      Do i=1,ncfg
       j=LEN_TRIM(LABEL(i))
       if(j.gt.ilab) ilab=j
      End do
      if(ilab.ge.64) Stop 'ilab > 63'

      if(allocated(C )) Deallocate(C ); Allocate(C (ncfg))
      if(allocated(ip)) Deallocate(ip); Allocate(ip(ncfg))

!----------------------------------------------------------------------
! ... read the solutions:

      mco = 0

      ns=0; Call Read_ipar(nuj,'nsol',ns);  if(ns.eq.0) Stop ' ns = 0'

      if(mconf.eq.0.or.mconf.lt.ncfg) mconf=ncfg
      if(mco.lt.mconf) mco = mconf

      i=Ifind_position(nuj,'Solutions');  read(nuj,*) 

      Do is = 1,ns

       read(nuj,*) j
       read(nuj,*) es,jot,i1,i2

       C = 0.d0

       read(nuj,*) C(i1:i2)
       Call Sorta(ncfg,C,IP)

! ...  define leading weights:

       nconf=1                                 
       Do i=2,mconf
        if(i.gt.ncfg) Exit
        if(abs(C(IP(i))).lt.eps_c) Exit; nconf=i
       End do

! ... save information:
 
       write(nua) es,jot,nconf,j,AF
       Do i=1,nconf
        j=IP(i); write(nua) C(j),LABEL(j)
       End do
       nsol=nsol+1         

      End do   

      Close(nuj)

      if(nfile.ne.1) go to 1
    
    2 Continue     

      write(*,*) ' nsol = ',nsol,mco

!----------------------------------------------------------------------
! ... restore all solutions:

      ns = nsol
      Allocate(EE(ns),IPE(ns),kco(ns),JJ(ns),CC(ns,mco),IPC(ns,mco))
      Allocate(LAB(ns,mco),isol(ns),AFF(ns))

      cc = 0.d0
      rewind(nua)
      Do i=1,ns
       read(nua) EE(i),JJ(i),nco,isol(i),AFF(i)
       kco(i)=nco
       Do j=1,nco
        read(nua) CC(i,j),LAB(i,j)
       End do
      End do
      cc=cc*cc 

! ... check the labels:

      go to 5
      
      i=2
    3 j=1
    4 Continue      
      if(LAB(i,1).ne.LAB(j,1)) then
        j=j+1
        if(j.lt.i) go to 4 
        i=i+1
        if(i.le.ns) go to 3
      elseif(cc(j,1).ge.cc(i,1)) then
        ij=0
        Do k=2,kco(i)
         if(cc(i,k).le.0.d0) Cycle
         S=-cc(i,1); cc(i,1)=cc(i,k); cc(i,k)=S
         LABEL1=LAB(i,1); LAB(i,1)=LAB(i,k); lab(i,k)=LABEL1
         if(LAB(i,1).eq.LAB(j,1)) Cycle
         ij=1
         Exit
        End do
        if(ij.eq.0) Stop 'problems with LABELS' 
        go to 3
      else
        ij=0
        Do k=2,kco(j)
         if(cc(j,k).le.0.d0) Cycle
         S=-cc(j,1); cc(j,1)=cc(j,k); cc(j,k)=S
         LABEL1=LAB(j,1); LAB(j,1)=LAB(j,k); lab(j,k)=LABEL1
         if(LAB(i,1).eq.LAB(j,1)) Cycle
         ij=1
         Exit
        End do
        if(ij.eq.0) Stop 'problems with LABELS' 
        i=j; if(i.eq.1) i=2
        go to 3
      end if 

    5 Continue

! ... sort of results:

      Call SortR(ns,EE,IPE)

      if(shift.eq.1.d0) shift = EE(IPE(1)) 

      Do i=1,ns; EE(i)=(EE(i) - shift)*conv; End do

! ... print the results:

      Open(nur,file=AFr)

      ms=ns; if(ms.gt.msol.and.msol.ne.0) ms=msol

      Do ii=1,ms
       i=IPE(ii);  nco=kco(i)
       BS = LAB(i,1); AS(1:ilab) = BS(1:ilab)
!       write(AS(1:ilab),'(a)') trim(LAB(i,1)); 
       k=ilab+1

       Select case(iev)
        case(0,1);  write(AS(k:),'(2x,f16.8)') EE(i); k=k+18
        case(2);    write(AS(k:),'(2x,f16.3)') EE(i); k=k+18
        case(3);    write(AS(k:),'(2x,f16.1)') EE(i); k=k+18
       End Select

       write(AS(k:),'(f8.1)') CC(i,1)*100; k=k+8

       write(AS(k:),'(i5,i6,2x,a)') ii,isol(i),trim(AFF(i))

       write(nur,'(a)') trim(AS)

!       Do j=2,nco
!        write(nur,'(f8.3,2x,a64)') CC(i,j),LAB(i,j)
!       End do

       AS = ' '
  
      End do

      End ! utility jj_list


