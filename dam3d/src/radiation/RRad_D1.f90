subroutine Rrad
use Rad3D
use Bem3D
implicit none
!
! DRM2D\Rad: Building  [h] and [g] and [S] Radiation Matrices
!
 print*,'DRM-2D: Form of  Radiation [G],[H],[S]'
!
 allocate(funct(ndr),derv_eta(3,ndr))
 allocate(coord(ndr,3),derv_x(3,ndr))
 allocate(N_funct(3,ndr),N_Der(3,ndr))  
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 allocate(Hrad(neqr,neqr),Grad(neqr,nelr,ndr))
!
 Hrad=0.0; Grad=0.0
!
 do i=1,neqr
 do j=1,nelr
 n(1:ndr)=elementr(j)%node(1:ndr) 
 do ij=1,ndr; coord(ij,:)=node(n(ij))%coord(:); enddo
!
! DRM2D\Rad: Evaluating h,g Terms using Guass Quadrature 
!		
 g=0.; h=0.
 do ii=1,ng
 call ShapeFunctions(ndr,eta(ii),.0)
 xx=0.; do ij=1,ndr; xx=xx+funct(ij)*node(n(ij))%coord; enddo
 r=sqrt(sum((xx-node(lmr(i))%coord)**2))
 drdx=(xx-node(lmr(i))%coord)/r; drdn=sum(drdx*norm) 
 g(1:ndr)=g(1:ndr)+funct(1:ndr)*ww(ii)*green()*det_jac
 h(1:ndr)=h(1:ndr)+funct(1:ndr)*ww(ii)*dgdr()*drdn*det_jac
 enddo;
!
! DRM2D\Rad: Arranging the Matrices [g], [h]
!  
 n(1:ndr)=node(n(1:ndr))%idr 
 grad(i,j,1:ndr)=g(1:ndr)
 hrad(i,n(1:ndr))=hrad(i,n(1:ndr))+h(1:ndr)
 enddo
 hrad(i,i)=-(sum(hrad(i,1:neqr))-hrad(i,i)) 
 enddo
!                                                                               
! DRM2D\Rad: Building  [S]                                                               
!                                                                               
 allocate(f_inv(neqr,neqr),srad(neqr,neqr))
 allocate(u_hat(neqr,neqr),q_hat(nelr,neqr,ndr))
!
! DRM2D\Rad: Building [alfa], [alfa]-1
!
 f_inv=0.0; Srad=0.; u_hat=0.0; q_hat=0.0 
 do i=1,neqr
 do j=1,neqr
 r=sqrt(sum((node(lmr(i))%coord-node(lmr(j))%coord)**2))
 f_inv(i,j)=c1+c2*r+c3*r**2+c4*r**3
 enddo; enddo
 if(neqr>0)call linrg (neqr, f_inv, neqr, f_inv, neqr)
!
! DRM2D\Rad: [S]=[H]*[uhat]-[G]*[qhat]
!
 do i=1,neqr
 do j=1,nelr
 n(1:ndr)=elementr(j)%node(1:ndr) 
 do ij=1,ndr; coord(ij,:)=node(n(ij))%coord(:); enddo
 do ij=1,ndr
 r=sqrt(sum((node(lmr(i))%coord-node(n(ij))%coord)**2))
 norm(1:3)=elementr(j)%norm(ij,1:3)
 drdx=(node(n(ij))%coord-node(lmr(i))%coord); drdn=sum(drdx*norm) 
 q_hat(j,i,ij)=(c1*1./2.+c2*r/3.+c3*(r**2)/4.+c4*(r**3)/5.)*drdn
 enddo; enddo
 enddo
 do ij=1,ndr
 srad=srad-matmul(grad(:,:,ij),q_hat(:,:,ij))
 enddo
!
 do i=1,neqr
 do j=1,neqr
 r=sqrt(sum((node(lmr(i))%coord-node(lmr(j))%coord)**2))
 u_hat(i,j)=c1*(r**2)/4.+c2*(r**3)/9.+c3*(r**4)/16.+c4*(r**5)/25.
 enddo
 enddo
 srad=srad+matmul(hrad,u_hat)
!
! DRM2D\Rad: Building [S]*[alfa]-1
!
 srad=matmul(srad,f_inv)
!
!                                                                               
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)
!
 deallocate(f_inv)
 deallocate(u_hat,q_hat)
!
 contains
!
! DRM2D\RAD : Definition of the green functions and derivatives
!
 real function green()
  green=-1./(2.*pi)*alog(1./r)   ; end function green
 real function dgdr()
  dgdr=1./(2.*pi*r)            ; end function dgdr
 real function d2gdr()
  d2gdr=-1./(2.*pi*r**2)	        ; end function d2gdr
!
end subroutine Rrad