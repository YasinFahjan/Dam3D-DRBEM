subroutine ghmat
! Quadratic BEM3D  
 use arrays
 w=1.0
!
! classical Gauss-Legendre quadrature formula using imsl
!
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
! Assembling [h] and [g] elements 
!
 a_mat=0.; xy_vect=0.; f_inv=0.0; s_mat=0.0
 do i=1,neq
 zeta=node(i)%x
 h_diag=0.
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod; coord_x(ij,:)=node(n(ij))%x(:); enddo
!
! Evaluating h,g terms using Guass quadrature 
!		
 g=0.; h=0.; dgdz=0.; dhdz=0.;
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(eta(ii),eta(jj))
! Jacobian Matrix & the Normal to the Element Surface
 jac_mat=matmul(fi_eta,coord_x)
 call cross_product(jac_mat(1,:),jac_mat(2,:),norm)
 det_jac=sqrt(sum(norm*norm)); norm=norm/det_jac
!
 xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
 r=sqrt(sum((xx-zeta)**2))
 drdx=(xx-zeta)/r; drdn=sum(drdx*norm) 
 g(1:nod)=g(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*green(r)*det_jac
 h(1:nod)=h(1:nod)+fi(1:nod)*ww(ii)*ww(jj)*dgdr(r)*drdn*det_jac
 enddo; enddo;
! Diagonal elements of h matrix
 h_diag=h_diag-sum(h(1:nod)) 
!
! arranging the system 	[a]{xy}={xy}
!  
 do ij=1,nod
 select case (node(n(ij))%code)
 case (1)
 xy_vect(i)=xy_vect(i)-h(ij)*node(n(ij))%p
 if(element(j)%code==1)a_mat(i,n(ij))=a_mat(i,n(ij))-g(ij)
 if(element(j)%code==2)xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
 case (2)
 xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
 a_mat(i,n(ij))=a_mat(i,n(ij))+h(ij)
 end select
 enddo
 enddo;
!
! Modification for Diagonal Term
!
! if(node(i)%code==0)h_diag=1.0
 select case (node(i)%code)
 case (0); a_mat(i,i)=a_mat(i,i)+1.0
 case (1); xy_vect(i)=xy_vect(i)-h_diag*node(i)%p
 case (2); a_mat(i,i)=a_mat(i,i)+h_diag
 end select
!
 enddo
!
end subroutine ghmat
