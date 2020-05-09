!=======================================================================
      Subroutine EL_NLJK(EL,n,kappa,l,j,k)
!=======================================================================
!
!     decodes the specroscopic notation for electron orbital (n,l,j,k)
!
!     It is allowed following notations: 1s , 2s 3, 2p-30, 20p-3,
!     1s h, 20s h, kp , kp 1, kp-11, ns , ns 3, ns 33, ... .
!
!     Call:  LA, INDEX
!
!----------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4), parameter :: kset = 61
      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

      Character(5), Intent(in) :: EL
      Integer(4), Intent(out) :: n,l,j,k,kappa    

      Integer(4) :: jj, k1,k2, n1,n2  

      Integer(4), EXTERNAL :: LA, kappa_lj

!----------------------------------------------------------------------
      jj=0
      Do j=5,3,-1
       if(EL(j:j).eq.'-') then; jj=j; Exit; end if
      End do
      if(jj.eq.0) then
       Do j=5,3,-1
        if(EL(j:j).eq.' '.and.EL(j-1:j-1).ne.' ') then
         jj=j; Exit
        end if
       End do
      end if

      if(jj.eq.5) then

       k = 0
       l = LA(EL(4:4))
       n1 = INDEX(ASET,EL(3:3))
       n2 = INDEX(ASET,EL(2:2))
       n = n2*kset+n1

      elseif(jj.eq.4) then

       k = INDEX(ASET,EL(5:5))
       l = LA(EL(3:3))
       n1 = INDEX(ASET,EL(2:2))
       n2 = INDEX(ASET,EL(1:1))
       n = n2*kset+n1

      elseif(jj.eq.3) then

       k1 = INDEX(ASET,EL(5:5))
       k2 = INDEX(ASET,EL(4:4))
       k = k2*kset+k1
       l = LA(EL(2:2))
       n = INDEX(ASET,EL(1:1))

      else

       write(*,*) 'EL_NLJK: can not decode ',EL
       Stop ' '

      end if

      j = l+l+1; if(EL(jj:jj).eq.'-') j = l+l-1
      kappa = kappa_lj(l,j)

      End Subroutine EL_NLJK


!=======================================================================
      Character(5) function ELi(n,kappa,k)     
!=======================================================================
!
!     provides the specroscopic notation for electron orbital (n,l,j,k)
!
!     k must be < 61*61; if k<=61 - one character from ASET
!
!-----------------------------------------------------------------------

      IMPLICIT NONE

      Integer(4) :: n,l,j,k, ll,jj, i,k1,k2,n1,n2, kappa

      Character(5) :: EL
      Character(1), EXTERNAL :: AL
      Integer(4), External :: l_kappa, j_kappa

      Integer(4), parameter :: kset = 61
      Character(61) :: ASET = &
       '123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

       l = l_kappa(kappa) 
       j = j_kappa(kappa)

      if(n.le.0.or.l.lt.0.or.j.le.0.or.k.lt.0) then
        write(*,'(a,4i5)') &
            ' ELj: parmeters are out of limits: n,l,j,k=',n,l,j,k
        Stop 'Stop in Elj'
      end if

      EL='    '; i=5

      if(k.lt.0) Stop 'ELi: set index < 0'
      if(k.gt.0) then
       if(k.le.kset) then
        EL(i:i)=ASET(k:k); i=i-1
       else
        k1=k/kset; k2=mod(k,kset); 
        if(k2.eq.0) then; k1=k1-1; k2=kset; end if
        if(k1.gt.kset) Stop 'ELi: set index too big'
        EL(i:i)=ASET(k2:k2); i=i-1
        EL(i:i)=ASET(k1:k1); i=i-1
       end if
      end if

      if(j.eq.l+l-1) EL(i:i) = '-'; i=i-1
      
      EL(i:i)=AL(l,1);  i=i-1

      if(n.lt.0) Stop 'ELi: n < 0'
      if(n.gt.0) then
       if(n.le.kset) then
        EL(i:i)=ASET(n:n)          
       else
        n1=n/kset; n2=mod(n,kset); 
        if(n2.eq.0) then; n1=n1-1; n2=kset; end if
        if(n1.gt.kset) Stop 'ELi: n is too big'
        EL(i:i)=ASET(n2:n2); i=i-1
        EL(i:i)=ASET(n1:n1); i=i-1
       end if
      end if

      ELi = EL

      End Function ELi




!====================================================================
      CHARACTER FUNCTION AL(L,k)
!====================================================================
!        
!     provides spectroscopic symbols for L values
!
!     Limits: L <= 153
!
!--------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(4), intent(in) :: L,K
      INTEGER(4) :: I
      CHARACTER(21), SAVE :: AS, AB

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = L+1; IF(k.eq.5.or.k.eq.6) i=(L-1)/2+1
      if(i.ge.1.and.i.le.21) then
       if(k.eq.1.or.k.eq.5) AL=AS(I:I)
       if(k.eq.2.or.k.eq.6) AL=AB(I:I)
      elseif(i.ge.22.and.i.le.153) then
       AL=CHAR(i+101)  ! from {[} ...
      else
       write(*,*) 'L,k=',L,k
       i = i/0
       Stop ' AL: L out of range'
      end if

      END FUNCTION AL


!====================================================================
      Integer(4) Function LA(a)
!====================================================================
!     gives the value of L from spetroscopic symbol "a"
!--------------------------------------------------------------------

      IMPLICIT NONE  
      Character, Intent(in) :: a
      CHARACTER(21), SAVE :: AS, AB
      Integer(4) :: i

      DATA AS/'spdfghiklmnoqrtuvwxyz'/
      DATA AB/'SPDFGHIKLMNOQRTUVWXYZ'/

      i = INDEX (AS,a)
      if(i.eq.0) i = INDEX (AB,a)
      if(i.eq.0) i = ICHAR(a)-101
      LA = i-1

      End Function LA
