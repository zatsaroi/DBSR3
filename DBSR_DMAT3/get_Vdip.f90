!======================================================================
      Subroutine Get_Vdip(i,j,k,side,v,w)
!======================================================================
!     convolutes dipole radial matrix with one orbital:
!----------------------------------------------------------------------
      Use dbsr_dmat
      Use DBS_orbitals_pq

      Implicit none
      Character(1) :: side   ! 'left' or 'right' 
      Integer :: i,j,k
      Real(8) :: C, v(ms), w(ms)

      if(ktype.eq.'M') then

       if(side.eq.'l') then
        v(1:ns) = dbs(1:ns,2,i); v(ns+1:ms)=dbs(1:ns,1,i)
       else
        v(1:ns) = dbs(1:ns,2,j); v(ns+1:ms)=dbs(1:ns,1,j)
       end if
       C  = kbs(i)+kbs(j); C = C / (k+1); v = C * v; w = 0.d0
        
      else

       if(side.eq.'l') then
        v(1:ns) = dbs(1:ns,1,i); v(ns+1:ms)=dbs(1:ns,2,i)
        w(1:ns) = vbs(1:ns,2,i); w(ns+1:ms)=vbs(1:ns,1,i)
        C = kbs(i)-kbs(j)+k;  w(   1:ns) = C * w(   1:ns) 
        C = kbs(i)-kbs(j)-k;  w(ns+1:ms) = C * w(ns+1:ms) 

       else
        v(1:ns) = dbs(1:ns,1,j); v(ns+1:ms)=dbs(1:ns,2,j)
        w(1:ns) = vbs(1:ns,2,j); w(ns+1:ms)=vbs(1:ns,1,j)
        C = kbs(i)-kbs(j)-k;  w(   1:ns) = C * w(   1:ns) 
        C = kbs(i)-kbs(j)+k;  w(ns+1:ms) = C * w(ns+1:ms) 

       end if                                              

      end if

      End Subroutine Get_Vdip