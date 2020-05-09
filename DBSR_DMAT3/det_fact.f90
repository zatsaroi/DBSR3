!======================================================================
      Subroutine Det_fact(idf,IPN1,IPN2)
!======================================================================
!     Analize the determinant overlap factor. 
!     First each overlap determinant is expended into one-electron
!     overlaps <nl|n'l'>, with overlaps included only known radial func.
!     to be directly calculated. Results are places in module new_dets. 
!     Then overlap determinant are combined to obtain the final 
!     overlap factors as:  C, C <kl|nl> or C <kl|k'l'>. 
!     Results are places in module new_defs.
!----------------------------------------------------------------------
      Use det_list;     Use new_dets 
      Use def_list;     Use new_defs 
   
      Implicit none
      Integer, intent(in) :: idf, IPN1(*),IPN2(*)
      Real(8) :: DF(jmdef)
      Integer :: NP(jmdef),NP1(jmdef),NP2(jmdef),IDPNT(jmdef)
      Integer :: i,j,ii,jj,io,jo,kd,ip,id,iext,kn,k

      nndef = 0       ! number of new det.factors in module new_defs
!----------------------------------------------------------------------
! ... trivial case:

      if(idf.eq.0) then; Call Iadd_ndefs(0,0,1.d0); Return; end if
      
! ... pick up the determinant factor:

      kd=KPF(idf); ip=IPF(idf);  NP(1:kd)=NPF(ip+1:ip+kd)

!----------------------------------------------------------------------
! ... find overlaps for specific orbitals with extracted continuum:  

      nndet = 0       ! number of new determinants in module new_dets
      Do j = 1,kd
       id=NP(j)/ibf; iext=mod(NP(j),ibf); kn=KPD(id); ip=IPD(id)
       Do i=1,kn
        k=NPD(i+ip); NP1(i)=IPN1(k/ibd); NP2(i)=IPN2(mod(k,ibd))
       End do
       Call Check_det(kn,NP1,NP2,iext)
       IDPNT(j)=nndet
       if(j.gt.1) then; if(IDPNT(j).eq.IDPNT(j-1)) nndet=0; end if
       if(nndet.eq.0) Return   ! exit if overlap = 0 
      End do 

!----------------------------------------------------------------------
! ... find new det.factors corrected with bound overlaps:

      i=1                !  index of overlap determinant 
      np(1)=1            !  index of chosen one-electron overlap 
    1 Continue

      if(i.eq.1) then
       k=np(i); DF(i)=ADET(k)
      else
       k=IDPNT(i-1)+np(i); DF(i)=DF(i-1)*ADET(k)
      end if

      if(i.lt.kd) then; i=i+1; np(i)=1; go to 1; end if

! ... find the final one-elctron overlaps with continuum functions:

      io=0; jo=0
      Do j=1,kd
       k=np(j); if(j.gt.1) k=IDPNT(j-1)+np(j) 
       ii = IZOD(k);  if(ii.eq.0) Cycle
       jj = JZOD(k)
       if(jj.eq.0) then
        if(io.eq.0) then
         io=ii
        elseif(io.gt.0.and.jo.eq.0) then
         jo=ii
        else
         Stop 'Det_fact: problems with io,jo'
        end if
       else
        if(io.eq.0) then
         io=ii;jo=jj
        else
         Stop 'Det_fact: problems with io,jo'
        end if
       end if   
      End do

      Call Iadd_ndefs(io,jo,DF(kd))

    2 np(i)=np(i)+1; k=np(i); if(i.gt.1) k=IDPNT(i-1)+np(i) 
      if(k.le.IDPNT(i)) go to 1
      i=i-1; if(i.gt.0) go to 2

      End Subroutine Det_fact
