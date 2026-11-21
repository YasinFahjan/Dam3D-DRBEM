 subroutine rhsmat
! Quadratic DRM3D  
 use arrays
!
 allocate(f_inv(neq,neq),s_mat(neq,neq))
 allocate(u_hat(neq,neq),q_hat(nelm,neq,nod))
 finv=0.0; S_mat=0.; u_hat=0.0; q_hat=0.0 
!
! DRM: Building [alfa], [alfa]-1
!
 do i=1,neq
 do j=1,neq
 r=sqrt(sum((node(i)%x-node(j)%x)**2))
 f_inv(i,j)=c1+c2*r+c3*r**2+c4*r**3
 enddo; enddo
 if(neq>0)call linrg (neq, f_inv, neq, f_inv, neq)
!
! DRM: [S]=[H]*[uhat]-[G]*[qhat]
!
 call CornerFunctions
!
 do i=1,neq
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord_x(ij,:)=node(n(ij))%x(:); enddo
 do ij=1,nod
 r=sqrt(sum((node(i)%x-node(n(ij))%x)**2))
 call ShapeFunctions(eta_c(1,ij),eta_c(2,ij))
 call norm_jac
 drdx=(node(n(ij))%x-node(i)%x); drdn=sum(drdx*norm) 
 q_hat(j,i,ij)=(c1*1./2.+c2*r/3.+c3*(r**2)/4.+c4*(r**3)/5.)*drdn
 enddo; enddo
 enddo
 do ij=1,nod
 s_mat=s_mat-matmul(g_mat(:,:,ij),q_hat(:,:,ij))
 enddo
!
 do i=1,neq
 do j=1,neq
 r=sqrt(sum((node(i)%x-node(j)%x)**2))
 u_hat(i,j)=c1*(r**2)/4.+c2*(r**3)/9.+c3*(r**4)/16.+c4*(r**5)/25.
 enddo
 enddo
 s_mat=s_mat+matmul(h_mat,u_hat)
!
! DRM: Building [S]*[alfa]-1
!
 s_mat=matmul(s_mat,f_inv)
!
 deallocate(u_hat,q_hat)
!
end subroutine rhsmat
