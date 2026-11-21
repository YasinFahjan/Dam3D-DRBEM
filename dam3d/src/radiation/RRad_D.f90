subroutine Rrad
use Rad3D
use Bem3D
implicit none
!
! DRM2D\Rad: Building [h],[g] and [S] Radiation Matrices
!
  write(*,'(/a)')'DRM-2D: Form of  Radiation [G],[H],[S]'
!
 allocate(funct(ndr),derv_eta(3,ndr))
 allocate(coord(ndr,3),derv_x(3,ndr))
 allocate(N_funct(3,ndr),N_Der(3,ndr))  
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 allocate(Hrad(neqr,neqr),G_rad(neqr,nelr,ndr),Srad(neqr,neqr))
 allocate(Grad(neqr,neqr),G11(nk1,nk1),drad(neqr))
!
 Hrad=0.0; Grad=0.0; Srad=0.
 G_rad=0.0; G11=0.0; Drad=0.
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
! DRM2D\Rad: Arranging the Matrices [grad],[hrad],[G_rad],{D_rad}
!  
 do ij=1,ndr
 jj=node(ij)%idr
! if(jj>nk1.and.node(n(ij))%codr==3)Drad(jj)=Drad(jj)+elementr(j)%q(ij)
 if(jj>nk1)Drad(jj)=Drad(jj)+elementr(j)%q(ij)
 enddo
 n(1:ndr)=node(n(1:ndr))%idr 
 g_rad(i,j,1:ndr)=g(1:ndr)
 hrad(i,n(1:ndr))=hrad(i,n(1:ndr))+h(1:ndr)
 grad(i,n(1:ndr))=grad(i,n(1:ndr))+g(1:ndr)
 enddo
 hrad(i,i)=-(sum(hrad(i,1:neqr))-hrad(i,i)) 
 enddo
!
! DRM2D\Rad: Building G11 from Grad
!
 g11(1:nk1,1:nk1)=grad(1:nk1,1:nk1)
 if(nk1>0)call linrg(nk1,g11,nk1,g11,nk1)
!
! DRM2D\Rad: [S]=-[G]*[qhat]
!                                                                               
 c=0
 c(3)=1.0
 allocate(q_hat(nelr,neqr,ndr))
 q_hat=0.0 
 do i=1,neqr
 do j=1,nelr
 n(1:ndr)=elementr(j)%node(1:ndr) 
 do ij=1,ndr; coord(ij,:)=node(n(ij))%coord(:); enddo
 do ij=1,ndr
 r=sqrt(sum((node(lmr(i))%coord-node(n(ij))%coord)**2))
 norm(1:3)=elementr(j)%norm(ij,1:3)
 drdx=(node(n(ij))%coord-node(lmr(i))%coord); drdn=sum(drdx*norm) 
! q_hat(j,i,ij)=(c(0)*1./2.+c(1)*r/3.+c(2)*(r**2)/4.+c(3)*(r**3)/5.)*drdn
 if(r/=0.0) q_hat(j,i,ij)=drdn*(r**2*(-1.0/16.+log(r)/4.0))
 enddo; enddo
 enddo
 do ij=1,ndr
 srad=srad-matmul(g_rad(:,:,ij),q_hat(:,:,ij))
 enddo
 deallocate(q_hat)
!
! DRM2D\Rad: [S]=[H]*[uhat]-[G]*[qhat]
!
 allocate(u_hat(neqr,neqr))
 u_hat=0.0 
 do i=1,neqr
 do j=1,neqr
 r=sqrt(sum((node(lmr(i))%coord-node(lmr(j))%coord)**2))
! u_hat(i,j)=c(0)*(r**2)/4.+c(1)*(r**3)/9.+c(2)*(r**4)/16.+c(3)*(r**5)/25.
 if(r/=0.0)u_hat(i,j)=r**4*(-1.0/32.+log(r)/16.)
 enddo
 enddo
 srad=srad+matmul(hrad,u_hat)
 deallocate(u_hat)
!
! DRM2D\Rad: Building [alfa], [alfa]-1
!
 allocate(f_inv(neqr,neqr))
 f_inv=0.0
 do i=1,neqr
 do j=1,neqr
 r=sqrt(sum((node(lmr(i))%coord-node(lmr(j))%coord)**2))
! f_inv(i,j)=c(0)+c(1)*r+c(2)*r**2+c(3)*r**3
 if(r/=0.0)f_inv(i,j)=r**2*log(r)
 enddo; enddo
 if(neqr>0)call dlinrg (neqr, f_inv, neqr, f_inv, neqr)
!
! DRM2D\Rad: Building [S]*[alfa]-1
!
 srad=matmul(srad,f_inv)
 deallocate(f_inv)
!                                                                               
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)

!     write(prn,'(/a/)')' Hrad :::: '
!     write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!     write(prn,'(i8,<neqr>f20.7)')(i,(Hrad(i,j),j=1,neqr),i=1,neqr)

!     write(prn,'(/a/)')' Grad :::: '
!     write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!     write(prn,'(i8,<neqr>f20.7)')(i,(Grad(i,j),j=1,neqr),i=1,neqr)

!     write(prn,'(/a/)')' srad :::: '
!     write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!     write(prn,'(i8,<neqr>f20.7)')(i,(srad(i,j),j=1,neqr),i=1,neqr)

!     write(prn,'(/a/)')' g11 :::: '
!     write(prn,'(4x,<nk1>i20/)')(i,i=1,nk1)
!     write(prn,'(i8,<nk1>f20.7)')(i,(g11(i,j),j=1,nk1),i=1,nk1)

!
! DRM2D\Rad: Building [H]-lamda[S]={d}
!
 Hrad(1:nqr,1:nqr)= Hrad(nk1+1:neqr,nk1+1:neqr)-matmul(&
 &  matmul(Grad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  Hrad(1:nk1,nk1+1:neqr))
!
 Srad(1:nqr,1:nqr)= Srad(nk1+1:neqr,nk1+1:neqr)-matmul(&
 &  matmul(Grad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  Srad(1:nk1,nk1+1:neqr))
!
 drad(1:nqr)=matmul((Grad(nk1+1:neqr,nk1+1:neqr)-matmul(&
 &  matmul(Grad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  Grad(1:nk1,nk1+1:neqr))),drad(nk1+1:neqr))
!
 Grad(1:nqr,1:nqr)= Grad(nk1+1:neqr,nk1+1:neqr)-matmul(&
 &  matmul(Grad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  Grad(1:nk1,nk1+1:neqr))
!
!
 if(run_type==1)then
 open(1004,file='drad.bin',form='binary',status='unknown',action='write')
 write(1004)Hrad
 write(1004)Grad
 write(1004)Srad
 close(1004)
 endif
!
 contains
!
! DRM2D\RAD : Definition of the green functions and derivatives
!
 real function green()
  green=1./(2.*pi)*alog(1./r)   ; end function green
 real function dgdr()
  dgdr=-1./(2.*pi*r)            ; end function dgdr
 real function d2gdr()
  d2gdr=1./(2.*pi*r**2)	        ; end function d2gdr
!
end subroutine Rrad