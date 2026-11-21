subroutine interm
! Quadratic BEM3D  
 use arrays
!
! solution at intenal points once boundary solution solution is known
!
 do i=nnod+1,nnod+nint
 node(i)%p=0.; node(i)%q=0.
 zeta=node(i)%x
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord_x(ij,:)=node(n(ij))%x(:); enddo
!    
! Evaluating the g, h, dg/d(zeta), dh/d(zeta)  
!   		
 g=0.; h=0.; dgdz=0.; dhdz=0.;
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(eta(ii),eta(jj))
 jac_mat=matmul(fi_eta,coord_x)
 call cross_product(jac_mat(1,:),jac_mat(2,:),norm)
 det_jac=sqrt(sum(norm*norm)); norm=norm/det_jac
 xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
 r=sqrt(sum((zeta-xx)**2))
 drdx=(xx-zeta)/r; drdn=sum(drdx*norm); drdnx=-norm+drdn*drdx
 g(1:nod)=g(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*green(r)*det_jac
 h(1:nod)=h(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*dgdr(r)*drdn*det_jac
 do ij=1,nod
 dgdz(ij,1:n_dim)=dgdz(ij,1:n_dim)+fi(ij)*ww(ii)*ww(jj)*&
  &dgdr(r)*(-drdx(1:n_dim))*det_jac
 dhdz(ij,1:n_dim)=dhdz(ij,1:n_dim)+fi(ij)*ww(ii)*ww(jj)*&
  &(d2gdr(r)*drdn*(-drdx(1:n_dim))+dgdr(r)/r*drdnx(1:n_dim))*det_jac
 enddo
 enddo; enddo;
!
! Evaluating p, q  
!		
 do ij=1,nod
 node(i)%p=node(i)%p-element(j)%p(ij)*h(ij)+element(j)%q(ij)*g(ij)
 node(i)%q(1:n_dim)=node(i)%q(1:n_dim)-element(j)%p(ij)*dhdz(ij,1:n_dim)&
  &+element(j)%q(ij)*dgdz(ij,1:n_dim)
 enddo
 enddo
 enddo

end subroutine interm
