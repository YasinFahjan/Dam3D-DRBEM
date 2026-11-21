subroutine Assemble
use arrays
! Quadratic DRM3D
 allocate(a_mat(neq,neq),xy_vect(neq))
 a_mat=0.0; xy_vect=0.0
!
! DRM: building the system [a]{xy}={xy}
!
 do i=1,neq
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod
 select case (node(n(ij))%code)
 case (1)
 if(element(j)%code==1)a_mat(i,n(ij))=a_mat(i,n(ij))-g_mat(i,j,ij)
 if(element(j)%code==2)xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
 case (2)
 xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
 end select
 enddo
 enddo
 enddo
!
 h_mat=h_mat+w**2*s_mat
 do i=1,neq
 do j=1,neq
 select case (node(j)%code)
 case (1)
 xy_vect(i)=xy_vect(i)-h_mat(i,j)*node(j)%p
 case (0,2)
 a_mat(i,j)=a_mat(i,j)+h_mat(i,j)
 end select
 enddo
 enddo
!
end subroutine Assemble

subroutine Re_Assemble
use arrays
! Quadratic DRM3D  
 where(node(1:neq)%code==0)node(1:neq)%p=xy_vect
 where(node(1:neq)%code==2)node(1:neq)%p=xy_vect
 do i=1,nelm
 n(1:nod)=element(i)%node(1:nod)
 do ij=1,nod
 select case (element(i)%code)
 case(1); element(i)%q(ij)=xy_vect(n(ij))
 case(2); element(i)%p(ij)=node(n(ij))%p
 end select
 enddo
 enddo
!
end subroutine Re_Assemble

subroutine node_code
use arrays
!
! Assigning code to the nodes on the element
!
 select case(code)
 case(1)
 element(i)%p(1:nod)=value(1:nod)
 node(n(1:nod))%p=value(1:nod)
 node(n(1:nod))%code=code
 case(2)
 element(i)%q(1:nod)=value(1:nod)
 do ij=1,nod; if(node(n(ij))%code/=1)node(n(ij))%code=2; enddo
 end select
!
end subroutine node_code
