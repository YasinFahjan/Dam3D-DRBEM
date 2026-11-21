subroutine Snodal
use fem3d
use modal
implicit none
!
! Modal: Nodal Displacement, acceleration, Pressure
!
 print*,'DRM-3D: Structural Nodal Analysis '
!
!
 Y=abs(L)
 YY=abs(-w**2*L)
!
 allocate(u(neq),uu(neq))
 u=0.0; uu=0.
 u(1:neq)=matmul(eigen_vect(:,1:nf),L(1:nf))
 uu(1:neq)=matmul(eigen_vect(:,1:nf),(-w**2*L(1:nf)))
 do i=1,numnp
 do j=1,3
 select case(node(i)%id(j))
 case(0)
 node(i)%disp(j)=0.0
 node(i)%accel(j)=0.0
 case	default
 node(i)%disp(j)=abs(u(node(i)%id(j)))
 node(i)%accel(j)=abs(uu(node(i)%id(j)))
 end select
 enddo
 enddo
!
 deallocate(u,uu)
!
end subroutine Snodal
