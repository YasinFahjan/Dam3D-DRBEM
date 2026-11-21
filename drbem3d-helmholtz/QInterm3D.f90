subroutine interm
! Quadratic DRM3D  
 use arrays
!
! calculate values of fluxes at intenal points once
! potential solution is known
!
 allocate(fx(nnod+1:neq,neq,n_dim))
 fx=0.0
 do i=nnod+1,neq
 do j=1,neq
 r=sqrt(sum((node(i)%x-node(j)%x)**2))
 drdx=(node(j)%x-node(i)%x) 
 if(r/=0.0)fx(i,j,:)=-(c2/r+2.*c3+3.*c4*r)*drdx
 enddo
 enddo
 do ij=1,n_dim
 fx(:,:,ij)=matmul(fx(:,:,ij),f_inv)
 node(nnod+1:neq)%q(ij)=matmul(fx(nnod+1:neq,:,ij),node(:)%p)
 enddo
 deallocate(fx)
end subroutine interm
