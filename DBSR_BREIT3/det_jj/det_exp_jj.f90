!======================================================================
!     PROGRAM  D E T _ E X P _ J J
!======================================================================
!     Expanding the j^q-subshell w.f. into determinants
!     on the basis of the expansion for j^(q-1)-subshell
!     by using parentage coefficients
!----------------------------------------------------------------------

      Implicit real(8) (A-H,O-P,R-Z)
      Implicit integer (q)

      Real(8) :: zero=0.d0, one = 1.d0

      Integer, Allocatable, Dimension(:) :: MD, MD1, ID
      Integer, Allocatable, Dimension(:,:) :: Idet,Idet1,ICD1,ICD2
      Real(8), Allocatable, Dimension(:,:) :: CD, CD1

      nua=1; open(nua,file='det_exp_jj.out')
      nub=2; open(nub,file='det_exp_jj.tab')
      nus=3; open(nus,form='UNFORMATTED',status='SCRATCH') 

!----------------------------------------------------------------------
! ... cycle over subshells:

      Do j = 1,11,2
      Do q = 2,j+1

       if(j.gt.7.and.q.gt.2) Exit

! ... prepare the parents determinants:

       q1 = q-1
       kd1 = Ndets_jq(j,q1)

       nt1 = Jterm(j,q1,-1,JT1,JV1,JW1,JQ1)
       Allocate(Idet1(q1,kd1), MD1(kd1), CD1(kd1,nt1))

       if(q1.eq.1) then
        if(nt1.ne.1) Stop 'nt1 <> 1 for q=1'
        if(kd1.ne.j+1) Stop 'kd1 <> j+1 for q=1'
        Do i=1,kd1
         Idet1(1,i) = i
         MD1(i) = mj_value(i)
         CD1(i,1) = 1.d0
        End do
       else
        rewind(nus) 
        read(nus) kd,nt
        if(kd.ne.kd1) Stop ' kd <> kd1'
        if(nt.ne.nt1) Stop ' nt <> nt1'
        read(nus) Idet1,MD1,CD1
       end if

! ... daughter determinants:

       kd = Ndets_jq(j,q)
       nt = Jterm(j,q,-1,JT,JV,JW,JQ)
       Allocate(Idet(q,kd), MD(kd), CD(kd,nt), ID(q), &
                ICD1(kd,nt),ICD2(kd,nt))
       CD = 0.d0

! ... find possible determinants and their coefficients

      k=0;  i=1;  ID(1)=0; mq = j+1

   11 ID(i)=ID(i)+1
      if(ID(i).GT.mq-q+1) go to 15          ! end of exhaustion
   12 i=i+1                              
      ID(i)=ID(i-1)+1
   13 if(i.lt.q) go to 12

      MJ = 0                                ! total azimutal numbers
      Do iq = 1,q
       MJ = MJ + mj_value(ID(iq))
      End do

      k = k + 1; if(k.gt.kd) Stop ' k > kd'
      MD(k) = MJ; Idet(1:q,k) = ID(1:q)

      Do it=1,nt                            ! cycle on terms
       ii  = Jterm(j,q,it,JT,JV,JW,JQ)
       if(IABS(MJ).GT.JT) Cycle
       CD(k,it) = CDET() 
      End do

! ... set next determinant

   14 ID(i)=ID(i)+1             
      if(ID(i).le.mq-q+i) go to 13  
      i=i-1
      if(i.EQ.1) go to 11
      go to 14

   15 CONTINUE

! ... check of normalization:

      eps = 1.d-5
      Do id1=1,kd; Do id2=id1,kd
       if(MD(id1).ne.MD(id2)) Cycle
       C = SUM(CD(id1,:)*CD(id2,:))
       if(id1.ne.id2.and.abs(C).lt.eps) Cycle
       if(id1.eq.id2.and.abs(1.d0-C).lt.eps) Cycle
       Call NUM(C,k1,k2,1000,1.d-6)
       write(nua,'(a,i1,a,i1,3x,a,2i4,3x,a,2i4,3x,a,f10.6,2i8)') &
                'subshell [',j,'/2]^',q, &
                'id1,id2=',id1,id2,'MJ1,MJ2=',MD(id1),MD(id2), &
                'C =',C,k1,k2
       Stop
      End do; End do

! ... record new results 

      rewind(nus) 
      write(nus) kd,nt
      write(nus) Idet,MD,CD

      write(nua,'(/a,i1,a,i1,5x,a,i3,4x,a,i3/)') &
                'subshell [',j,'/2]^',q, 'kd =',kd, 'nt =',nt

      Do ik = 1,kd
      Do it = 1,nt
       Call NUM(CD(ik,it),ICD1(ik,it),ICD2(ik,it),1000,1.d-8)
      End do; End do

      Do i = 1,kd
       write(nua,'(i3,a,i4,3x,10i2)') & 
         i,'.  ',MD(i),Idet(:,i)
       write(nua,'(20x,10f12.8)') CD(i,:)
      End do

      Call Print_table

      Deallocate(ID, Idet,Idet1, MD,MD1, CD,CD1, ICD1,ICD2)

      End do   ! ovet q
      End do   ! ovet j

CONTAINS

!====================================================================
      Real(8) Function CDET()
!====================================================================
!     expansion coefficients for given detrminant and term
!--------------------------------------------------------------------

      Implicit real(8) (A-H,O-Z)

      sq = q; sq = sqrt(sq)               ! normalization factor

      CDET = 0.d0

! ... cycle on parent det.'s

      Do id1=1,kd1                       

! ... mmj - 2*M of the additional electron:

       mmj=MJ-MD1(id1);    if(IABS(mmj).GT.j) Cycle

! ... its index in the common list: 

       ip = Index_mj(mmj)  
       
! ... its index in the given determinant: 

       ia = Ipointer (q,ID,ip);  if(ia.eq.0) Cycle

! ... check if the det.id1. is the parent for det.ID

       m=0
       Do iq=1,q             
        if(iq.EQ.ia) Cycle
        m=m+1
        if(ID(iq).ne.Idet1(m,id1)) then; m=-1; Exit; end if
       End do
       if(m.eq.-1) Cycle

! ... sign factor of permutation the additional electron
!     from the end to I

       kz=(-1)**(q-ia)

!--------------------------------------------------------------------

      Do it1=1,nt1                   ! cycle on parent terms

       C = CD1(id1,it1)*kz; if(C.eq.0.d0) Cycle
                                  
       ii = Jterm(j,q-1,it1,JT1,JV1,JW1,JQ1)

       C = C * cfp_jj(j,q,JT1,JV1,JT,JV); if(C.eq.0.d0) Cycle
                                    
       C = C * CLEBSH2(JT1,MD1(id1),j,mmj,JT,MJ)
                                    

       CDET = CDET + C

      End do         ! over parents terms
      End do         ! over parent determinants

       CDET = CDET/sq


      END FUNCTION CDET


!====================================================================
  Subroutine Print_table
!====================================================================
! output results in table format
!--------------------------------------------------------------------

  Character( 8) :: akd, ant
  Character( 8) :: aMD = 'MD_j#_q#'
  Character( 8) :: aID = 'ID_j#_q#'
  Character( 8) :: aJD = 'JD_j#_q#'
  Character(10) :: aIdet = 'Idet_j#_q#'
  Character(16) :: sjq = 'subshell [j/2]^q'
  Character(80) :: AS

  write(sjq(11:11),'(i1)') j
  write(sjq(16:16),'(i1)') q
  write(akd,'(a4,i1,a2,i1)') 'kd_j',j,'_q',q
  write(ant,'(a4,i1,a2,i1)') 'nt_j',j,'_q',q
  write(aMD,'(a4,i1,a2,i1)') 'MD_j',j,'_q',q
  write(aID,'(a4,i1,a2,i1)') 'ID_j',j,'_q',q
  write(aJD,'(a4,i1,a2,i1)') 'JD_j',j,'_q',q
  write(aIdet,'(a6,i1,a2,i1)') 'Idet_j',j,'_q',q

  write(nub,'(/a6,a/)') '! ... ',sjq
  write(nub,'(6x,3(a),i3)') 'Integer(4), parameter :: ',akd,' =',kd
  write(nub,'(6x,3(a),i3)') 'Integer(4), parameter :: ',ant,' =',nt
  write(nub,*)
  write(nub,'(6x,3(a),i1,3(a))') 'Integer(4) :: ',aIdet,'(',q,',',akd,')'
  write(nub,'(6x,6(a))') 'Integer(4) :: ',aMD,'(',akd,')'
  write(nub,'(6x,7(a))') 'Integer(4) :: ',aID,'(',ant,',',akd,')'
  write(nub,'(6x,7(a))') 'Integer(4) :: ',aJD,'(',ant,',',akd,')'

  write(nub,'(/6x,a,a,a)') 'Data ',aIdet,'/ &'
  Do ikd = 1,kd
   write(AS,'(6x,10(i5,a1))') (Idet(i,ikd),',',i=1,q)  
   ip=6+q*6
   AS(ip:ip+2)=', &'; if(ikd.eq.kd) AS(ip:ip+2)='  /'
   write(nub,'(a)') AS(1:ip+2)
  End do

  write(nub,'(/6x,a,a,a)') 'Data ',aMD,'/ &'

  i1 = 1; i2 = 10; if(i2.gt.kd) i2=kd
  Do 
   write(AS,'(6x,10(i5,a1))') (MD(i),',',i=i1,i2)  
   ip=6+(i2-i1+1)*6
   AS(ip:ip+2)=', &'; if(i2.eq.kd) AS(ip:ip+2)='  /'
   write(nub,'(a)') AS(1:ip+2)
   i1=i1+10; if(i1.gt.kd) Exit; i2=i2+10; if(i2.gt.kd) i2=kd
  End do

  write(nub,'(/6x,a,a,a)') 'Data ',aID,'/ &'

  Do ikd = 1,kd
   write(AS,'(6x,10(i5,a1))') (ICD1(ikd,i),',',i=1,nt)  
   ip=6+nt*6
   AS(ip:ip+2)=', &'; if(ikd.eq.kd) AS(ip:ip+2)='  /'
   write(nub,'(a)') AS(1:ip+2)
  End do

  write(nub,'(/6x,a,a,a)') 'Data ',aJD,'/ &'

  Do ikd = 1,kd
   write(AS,'(6x,10(i5,a1))') (ICD2(ikd,i),',',i=1,nt)  
   ip=6+nt*6
   AS(ip:ip+2)=', &'; if(ikd.eq.kd) AS(ip:ip+2)='  /'
   write(nub,'(a)') AS(1:ip+2)
  End do

  END Subroutine Print_table


  END   !  PROGRAM  D E T _ E X P _ J J




