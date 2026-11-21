subroutine GSOutput
use Fem3D
use graphics
use modal
implicit none
!
! Modal: Input ground motion excitation direction, magnitude, frequency
!
 print*,' Dam-3D: Structural Output'
!
 ii=GHY(gn1,gn2/2+1)
 if(guw==1)then
 write(prn7,'(4G12.4)')w/wo,node(ii)%disp(1:3)
 endif

 if(gaw==1)then
 write(prn8,'(4G12.4)')w/wo,node(ii)%accel(1:3)
 endif

 if(gyw==1)then
 write(prn9,'(<nf+1>G12.4)')w/wo,YY(1:nf)
 endif
!
end subroutine GSOutput
