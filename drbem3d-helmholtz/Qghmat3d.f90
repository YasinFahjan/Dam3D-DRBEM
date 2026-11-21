subroutine ghmat
! Quadratic DRM3D  
 use arrays
!
! classical Gauss-Legendre quadrature formula using imsl
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
! Assembling [h] and [g] elements 
!
 allocate(h_mat(neq,neq),g_mat(neq,nelm,nod))
 allocate(fi(nod),fi_eta(n_dim,nod),eta_c(n_dim,nod))
 allocate(coord_x(nod,n_dim),jac_mat(n_dim,n_dim),norm(n_dim))
!
 h_mat=0.; g_mat=0.0
 do i=1,neq
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord_x(ij,:)=node(n(ij))%x(:); enddo
!
! Evaluating h,g terms using Guass quadrature 
!		
 g=0.; h=0.
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(eta(ii),eta(jj))
 call norm_jac
 xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
 r=sqrt(sum((xx-node(i)%x)**2))
 drdx=(xx-node(i)%x)/r; drdn=sum(drdx*norm) 
 g(1:nod)=g(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*green(r)*det_jac
 h(1:nod)=h(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*dgdr(r)*drdn*det_jac
 enddo; enddo;
!
! arranging the matrices [g], [h]
!  
 g_mat(i,j,1:nod)=g(1:nod)
 h_mat(i,n(1:nod))=h_mat(i,n(1:nod))+h(1:nod)
 enddo
 h_mat(i,i)=-(sum(h_mat(i,1:neq))-h_mat(i,i)) 
 enddo
! 
end subroutine ghmat
