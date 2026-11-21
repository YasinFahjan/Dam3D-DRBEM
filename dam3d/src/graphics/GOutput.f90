subroutine GOutput
use BEM3D 
use graphics
implicit none
!                                                                               
! Dam3D: Output Results                                                        
!
 write(*,'(a)')' Dam-3D: Printing The results'
!
 if(ndyn<=2)call GSOutput
!
 if(ndyn>=2)call GHOutput
!
end subroutine GOutput