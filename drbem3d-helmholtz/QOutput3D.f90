subroutine Output
! Quadratic DRM3D  
 use arrays 
!
! Output of the boundary & internal points results
!
 write(prn,'(//a/)')'> boundary results'
 write(prn,'(a,8x,<nod>(a,i1,9x),a,8x,<nod>(a,i1,12x)/)')&
     &'element',('   p',i,i=1,nod),'code', ('q',i,i=1,nod)
 write(prn,'(i5,3x,<nod>f14.6,i10,<nod>f14.6)')(i,element(i)%&
     &p(1:nod),element(i)%code,element(i)%q(1:nod),i=1,nelm)
 if(nint.gt.0)then
 write(prn,'(//a/)')'> results at interior nodes'
 write(prn,'(4x,8(a2,12x)/)')'node','x','y','z','p','qx','qy','qz'  
 write(prn,'(i5,4x,7f14.6,5x,i2)')(i,node(i),i=nnod+1,nnod+nint)
 endif
 close (inp)
 close (prn)
 open(22,file='test1.dat')
 open(23,file='test2.dat')
 write(22,'(4x,8(a,12x)/)')'no','x','y','z','p','qx','qy','qz'
 write(23,'(4x,8(a,12x)/)')'no','x','y','z','p','qx','qy','qz'
 do i=nnod+1,nnod+nint
 write(23,'(i8,4f14.6\)')i,node(i)%x,node(i)%p
 if(w/=0)write(23,'(64x,3f14.6)')node(i)%q/w
 enddo
 do i=1,nnod
 if(node(i)%x(3)==0.0.and.node(i)%x(2)==.5)then
 write(22,'(i8,7f14.6)')i,node(i)%x,node(i)%p,node(i)%q
 endif
 enddo
 close(22)
 close(23)

end subroutine Output
 
