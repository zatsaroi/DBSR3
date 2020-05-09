!======================================================================
      Subroutine Open_jj(nu,fail)
!======================================================================
!     open (closed) files in the bj program
!----------------------------------------------------------------------

      USE param_jj

      Implicit none
      Integer, Intent(in) :: nu
      Integer :: i,iarg,fail
      Character(3) :: ALSP
      Logical :: EX
      Integer, External :: Icheck_file, IARGC

       write(ALSP,'(i3.3)') klsp

       Select case(nu)

       Case(1)               ! c-file

        AF = AF_c
        if(klsp.gt.0) then
         i=Index(BF_c,'.'); AF=BF_c(1:i)//ALSP; BF_c=AF
        else
         iarg = IARGC()
         if(iarg.gt.0) then
          Call GETARG(1,BF); i=INDEX(BF,'=')
          if(i.eq.0) then 
           i=INDEX(BF,'.')-1; if(i.lt.0) i=LEN_TRIM(BF)
           AF=BF(1:i)//'.c'; AF_b = BF(1:i)//'.bnk'
          end if
         end if
        end if  

        if(Icheck_file(AF).eq.0) then; fail=1; Return; end if         
        Open(nu,file=AF,status='OLD')

       Case(2)               ! bnk-file

        AF = AF_b
        if(klsp.gt.0) then
         i=Index(BF_b,'.'); AF=BF_b(1:i)//ALSP; BF_b=AF
        end if  
        Inquire (FILE=AF, exist=EX)

        new=.TRUE.
        if(EX) then
         new = .FALSE.
         Open(nub,file=AF,form='UNFORMATTED',STATUS='OLD')
        end if   

        Case(3)

         AF = AF_r
         if(klsp.gt.0) then
          i=Index(BF_r,'.'); AF=BF_r(1:i)//ALSP; BF_r=AF
         end if  
         Open(nur,file=AF,form='UNFORMATTED')

        Case(7)

         AF = AF_pri;  Open(nu,file=AF)

        Case(11)

         Open(nui,file=AF_i,form='UNFORMATTED')

        Case(12,13)

         Open(nu,form='UNFORMATTED',status='SCRATCH')

        Case default
        
         Stop ' open_adcf: unit value <nu> is out of list '

       End select

       End Subroutine Open_jj



