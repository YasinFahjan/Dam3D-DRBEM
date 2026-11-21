subroutine Rrad
use Rad3D
use Bem3D
implicit none
!
! FEM2D\Rad: [B],{D} 
!
  write(*,'(/a)')'FEM-2D: Form of  Radiation [G],[H],[B]'
!
 allocate(funct(ndr),derv_eta(3,ndr))
 allocate(coord(ndr,3),derv_x(3,ndr))
 allocate(N_funct(3,ndr),N_Der(3,ndr))  
 allocate(eta_c(3,ndr))
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 allocate(Brad(ndr,ndr),Drad(ndr))
 allocate(B_rad(nqr,nqr),D_rad(nqr))
!
 B_rad=0.0; D_rad=0.0
!
 do j=1,nelr
 Brad=0.0; Drad=0.0; coord=0.0
 n(1:ndr)=elementr(j)%node(1:ndr) 
 do ij=1,ndr; coord(ij,:)=node(n(ij))%coord(:); enddo
!
! FEM2D\Rad: Building [B],{D} 
!
 do ii=1,ndr
 do jj=1,ndr
 do i=1,ng			     
 call ShapeFunctions(ndr,eta(i),.0)
 Brad(ii,1:ndr)=Brad(ii,1:ndr)+ww(i)*det_jac*funct(ii)*funct(1:ndr)
 Drad(1:ndr)=Drad(1:ndr)+ww(i)*det_jac*funct(1:ndr)*elementr(i)%q(ii)
 enddo
 enddo; enddo
!
! FEM2D\Rad: Global [B],{D}
!
 n(1:ndr)=node(n(1:ndr))%idr-nk1 
 do ii=1,ndr
 do jj=1,ndr
 if(n(ii)>0.and.n(jj)>0)then
 B_rad(n(ii),n(jj))=B_rad(n(ii),n(jj))+Brad(ii,jj)
 endif
 enddo
 if(n(ii)>0.and.n(jj)>0)D_rad(n(ii))=D_rad(n(ii))+Drad(ii)
 enddo
 enddo
!
 deallocate(Brad,Drad)
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)
 deallocate(eta_c)
!
! FEM2D\Rad: Building [G],[H] 
!
 deallocate(lmr); allocate(lmr(nod))
 allocate(funct(nod),derv_eta(3,nod))
 allocate(coord(nod,3),derv_x(3,nod))
 allocate(N_funct(3,nod),N_Der(3,nod))  
 allocate(eta_c(3,nod))
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
 call CornerFunctions(nod)
!
 allocate(Hrad(nod,nod),Grad(nod,nod))
 allocate(H_rad(nqr,nqr),G_rad(nqr,nqr))
!
 H_rad=0.0; G_rad=0.00
!
 do j=1,nelm
 if(element(j)%code==4)then
 Hrad=0.0; Grad=0.0;  coord=0.0; lmr=0
 n(1:nod)=element(j)%node(1:nod) 
 do ij=1,nod
 coord(ij,1:3)=node(n(ij))%coord(1:3)
 if(node(n(ij))%idr>nk1)lmr(ij)=node(n(ij))%idr-nk1
 enddo
 do ii=1,ng
 do jj=1,ng
 jac_inv=0.0; N_funct=0.0; N_Der=0.0; jac_mat=0.0;
 call ShapeFunctions(nod,eta(ii),eta(jj))
 N_funct(1,1:nod)=funct(1:nod)
 jac_mat(1:3,1:3)=matmul(derv_eta(1:3,1:nod),coord(1:nod,1:3))
 jac_inv(2,1)= jac_mat(2,3)
 jac_inv(2,2)=-jac_mat(1,3)
 jac_inv(3,1)=-jac_mat(2,2)
 jac_inv(3,2)= jac_mat(1,2)
 det_jac=(jac_mat(1,2)*jac_mat(2,3)-jac_mat(1,3)*jac_mat(2,2))
 jac_inv(1:3,1:3)=jac_inv(1:3,1:3)/det_jac
 N_Der(1:3,1:nod)=matmul(jac_inv(1:3,1:3),derv_eta(1:3,1:nod)) 
 Grad(1:nod,1:nod)=Grad(1:nod,1:nod)+ww(ii)*ww(jj)*det_jac&
 &*matmul(transpose(N_funct),N_funct)
 Hrad(1:nod,1:nod)=Hrad(1:nod,1:nod)+ww(ii)*ww(jj)*det_jac&
 &*matmul(transpose(N_Der(1:3,1:nod)),N_Der(1:3,1:nod))
 enddo
 enddo
!
! FEM2D\Rad: Global [G],[H]
!
 do ii=1,nod
 do jj=1,nod
 if(lmr(ii)/=0.and.lmr(jj)/=0)then
 H_rad(lmr(ii),lmr(jj))=H_rad(lmr(ii),lmr(jj))+Hrad(ii,jj)
 G_rad(lmr(ii),lmr(jj))=G_rad(lmr(ii),lmr(jj))+Grad(ii,jj)
 endif
 enddo
 enddo
!
 endif
 enddo
!
 deallocate(Hrad,Grad)
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)
 deallocate(eta_c)
!
 if(run_type==1)then
 open(1004,file='frad.bin',form='binary',status='unknown',action='write')
 write(1004)H_rad
 write(1004)G_rad
 write(1004)b_rad
 close(1004)
 endif
!
end subroutine Rrad
