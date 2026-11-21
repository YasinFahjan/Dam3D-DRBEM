subroutine HOutput
use Bem3D
implicit none
!
! BEM3D: Hydrodynamic Pressure results of  boundary & internal points 
!
 print*,' DRM-3D: Hydrodynamic Output'
!
 write(prn,'(/a/)')     '[ HYDRODYNAMIC PRESSURE RESULTS ] '
 write(prn,'(//a/)')    'Dam3D: Boundary Element Results '
 write(prn,'(a,10x,<nod>(a,i1,a,11x),a,8x,<nod>(a,i1,a,11x)/)')&
 &'element',('p',i,' [real,complex]',i=1,nod),'code',&
 &('q',i,' [real,complex]',i=1,nod)
 write(prn,'(i5,3x,<2*nod>f14.3,i12,<2*nod>f14.3)')(i,element(i)%p(1:nod),&
 &element(i)%code,element(i)%q(1:nod),i=1,nelm)
 if(nint.gt.0)then
 write(prn,'(//a/)')    'Dam3D: interior nodes Results '
 write(prn,'(4x,8(a,12x)/)')'node','x','y','z','p [real,complex]'&
 &,'qx [real,complex]','qy [real,complex]','qz [real,complex]'  
 write(prn,'(i5,4x,11f14.3)')(i,node(i)%coord,node(i)%p,node(i)%q&
 &,i=ntot+1,ntot+nint)
 endif
!                                                                               
! BEM3D: hydrodynamic Pressure Output for the graph files 
!
 open(22,file='test1.dat')
 write(22,'(4x,6(a,12x)/)')'x','y','z','p[real]','p[complex]','p[abs]'
 node%p=node%p*gravity/rho_w/hwater
 do i=1,ntot+nint
! if(i>ntot.or.node(i)%code/=0.0)then
 if(node(i)%coord(1)==.0.and.node(i)%coord(2)==0.0)then
 write(22,'(6f14.6)')node(i)%coord,node(i)%p,abs(node(i)%p)
 endif
 enddo
 close(22)

 open(23,file='test2.dat')
 write(23,'(4x,6(a,12x)/)')'x','y','z','p[real]','p[complex]','p[abs]'
 do i=1,ntot+nint
 if(node(i)%coord(3)==0.0.and.node(i)%coord(2)==0.0)then
 write(23,'(6f14.6)')node(i)%coord,node(i)%p,abs(node(i)%p)
 endif
 enddo
 close(23)
!
end subroutine HOutput
 
