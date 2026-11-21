subroutine Hfluid
use Bem3D
implicit none
!                                                                               
! BEM3D: Building  [h] and [g] and [S] Matrices                                                               
!                                                                               
 print*,'DRM-3D: Form of [G],[H],[S]'
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 allocate(h_mat(nq,nq),g_mat(nq,nelm,nod))
 allocate(funct(nod),derv_eta(3,nod))
 allocate(coord(nod,3),derv_x(3,nod))
 allocate(N_funct(3,nod),N_Der(3,nod))  
!
 print*,' Building [G] and [H] matrices '

 h_mat=0.; g_mat=0.0
 do i=1,nq
  write(*,'(i5\)')i
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord(ij,:)=node(n(ij))%coord(:); enddo
!
! BEM3D: Evaluating h,g Terms using Guass Quadrature 
!		
 g=0.; h=0.
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(nod,eta(ii),eta(jj))
 xx=0.; do ij=1,nod; xx=xx+funct(ij)*node(n(ij))%coord; enddo
 r=sqrt(sum((xx-node(lmw(i))%coord)**2))
 drdx=(xx-node(lmw(i))%coord)/r; drdn=sum(drdx*norm) 
 g(1:nod)=g(1:nod)+funct(1:nod)*ww(ii)*ww(jj)*green()*det_jac
 h(1:nod)=h(1:nod)+funct(1:nod)*ww(ii)*ww(jj)*dgdr()*drdn*det_jac
 enddo; enddo;
!
! BEM3D: Arranging the Matrices [g], [h]
!  
 n(1:nod)=node(n(1:nod))%idw 
 g_mat(i,j,1:nod)=g(1:nod)
 h_mat(i,n(1:nod))=h_mat(i,n(1:nod))+h(1:nod)
 enddo
 h_mat(i,i)=-(sum(h_mat(i,1:nq))-h_mat(i,i)) 
 enddo
!                                                                               
! DRM3D: Building  [S] Matrix                                                             
!                                                                               
 cterm=3
 c=0; c(3)=1.0
 allocate(s_mat(nq,nq))
 s_mat=0.0
!
! DRM3D: [S]=[H]*[uhat]-[G]*[qhat]
!
 write(*,'(/a)')' DRM-3D: Building [G]*[qhat] '
 allocate(q_hat(nelm,nq,nod))
 q_hat=0.0
 do i=1,nq
 write(*,'(i5\)')i
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord(ij,:)=node(n(ij))%coord(:); enddo
 do ij=1,nod
 r=sqrt(sum((node(lmw(i))%coord-node(n(ij))%coord)**2))
 norm(1:3)=element(j)%norm(ij,1:3)
 drdx=(node(n(ij))%coord-node(lmw(i))%coord); drdn=sum(drdx*norm) 
! do k=0,cterm
! q_hat(j,i,ij)=q_hat(j,i,ij)+drdn*(c(k)*(r**(k)/(k+3.0)))
! enddo
 if(r/=0.0) q_hat(j,i,ij)=drdn*(r**2*(-1.0/25.+log(r)/5.0))
 enddo; enddo
 enddo
 do ij=1,nod
 s_mat=s_mat-matmul(g_mat(:,:,ij),q_hat(:,:,ij))
 enddo
 deallocate(q_hat)
!
 write(*,'(/a)')' DRM-3D: Building [H]*[uhat] '
 allocate(u_hat(nq,nq))
 u_hat=0.0
 do i=1,nq
 write(*,'(i5\)')i
 do j=1,nq
 r=sqrt(sum((node(lmw(i))%coord-node(lmw(j))%coord)**2))
! do k=0,cterm
!  u_hat(i,j)=u_hat(i,j)+c(k)*(r**(k+2)/((k+3.0)*(k+2.0)))
! enddo
 if(r/=0.0)u_hat(i,j)=r**4*(-9.0/400.+log(r)/20.)
 enddo
 enddo
 s_mat=s_mat+matmul(h_mat,u_hat)
 deallocate(u_hat)
!
! DRM3D: Building [alfa], [alfa]-1
!
 write(*,'(/a)')' DRM-3D: Building [F]-1 '
 allocate(f_inv(nq,nq))
 f_inv=0.0
 do i=1,nq
 write(*,'(i5\)')i
 do j=1,nq
 r=sqrt(sum((node(lmw(i))%coord-node(lmw(j))%coord)**2))
! do k=0,cterm
!  f_inv(i,j)=f_inv(i,j)+c(k)*r**(k)
! enddo
 if(r/=0.0)f_inv(i,j)=r**2*log(r)
 enddo; enddo
 if(nq>0)call dlinrg (nq, f_inv, nq, f_inv, nq)
!
! DRM3D: Building [S]*[alfa]-1
!
 write(*,'(/a)')' Building [S] matrix '
 s_mat=matmul(s_mat,f_inv)
 deallocate(f_inv)
!
!
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)
!
!
 if(run_type==1)then
 open(1003,file='fluid.bin',form='binary',status='unknown',action='write')
 write(1003)h_mat
 write(1003)g_mat
 write(1003)s_mat
 close(1003)
 endif
!
 contains
!
! BEM3D: Definition of the green functions and derivatives
!
 real function green()
 green=1./(4.*pi*r)         ; end function green
 real function dgdr()
 dgdr=-1./(4.*pi*r**2)      ; end function dgdr
 real function d2gdr()
 d2gdr=1./(2.*pi*r**3)	     ; end function d2gdr
!
end subroutine Hfluid
