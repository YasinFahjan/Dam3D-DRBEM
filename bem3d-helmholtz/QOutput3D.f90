subroutine Output
! Quadratic BEM3D  
 use arrays 
!
! Output of the boundary & internal points results
!
 write(prn,'(//a/)')'> boundary results'
 write(prn,'(a,f14.5/)')'w=',w
 write(prn,'(a,10x,<nod>(a,i1,a,11x),a,8x,<nod>(a,i1,a,11x)/)')&
 &'element',('p',i,' [real,complex]',i=1,nod),'code',&
 &('q',i,' [real,complex]',i=1,nod)
 write(prn,'(i5,3x,<2*nod>f14.6,i12,<2*nod>f14.6)')(i,element(i)%p(1:nod),&
 &element(i)%code,element(i)%q(1:nod),i=1,nelm)
 if(nint.gt.0)then
 write(prn,'(//a/)')'> results at interior nodes'
 write(prn,'(4x,8(a,12x)/)')'node','x','y','z','p [real,complex]'&
 &,'qx [real,complex]','qy [real,complex]','qz [real,complex]'  
 write(prn,'(i5,4x,11f14.6,5x,i2)')(i,node(i),i=nnod+1,nnod+nint)
 endif
 write(prn,'(i5,4x,4f14.6)')(i,node(i)%x,abs(node(i)%p),i=nnod+1,nnod+nint)

end subroutine Output
 
