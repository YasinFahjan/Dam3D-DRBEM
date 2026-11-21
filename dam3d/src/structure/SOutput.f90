subroutine SOutput
use Fem3D
implicit none
!
! Modal: Input ground motion excitation direction, magnitude, frequency
!
 print*,' Dam-3D: Structural Output'
!
 write(prn,'(//a/)')'Displacement and Accelerations' 
 do i=1,numnp
 write(prn,'(i4,5x,6g14.5)')i,node(i)%disp,node(i)%accel
 enddo
!
end subroutine SOutput
