subroutine Output
use BEM3D 
implicit none
!                                                                               
! Dam3D: Output Results                                                        
! 
 write(*,'(/a)')' Dam-3D:  Printing The Results'
!
 write(prn,'(a,f14.5/)')'w=',w
!
! writing the results
!
 if(ndyn<=3)call SOutput
!
 if(ndyn>=2)call HOutput
!
 close (inp)
 close (prn)
!
 write(*,'(a/)')' Dam-3D:  End of the Analysis'

end subroutine Output