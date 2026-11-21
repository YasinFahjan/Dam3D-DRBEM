subroutine Output
! Linear BEM3D  
 use arrays 
!
! Output of the boundary & internal points results
!
  write(*,'(/a)')" Output"
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
!
! Output for graph files
!
end subroutine Output   
			