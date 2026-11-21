subroutine HNodal
use Bem3D
use Modal
implicit none
!
! DRM-3D: Nodal Hydrýodynamic Response 
!
 print*,' DRM-3D: Hydrodynamic Nodal Analysis '
!
 xy_vect=0
!
 do j=1,nelm
 if(element(j)%code==2)then
 fi_r=0.0 
 n(1:nod)=element(j)%node(1:nod)
!
! DRM-3D: Normal Component of the mode shapes
!
 do ij=1,nod
 norm(1:3)=element(j)%norm(ij,1:3)
 id1=node(n(ij))%id(1); id2=node(n(ij))%id(2); id3=node(n(ij))%id(3)
 if(id1/=0)fi_r(ij,1:nf)=fi_r(ij,1:nf)+eigen_vect(id1,1:nf)*norm(1)
 if(id2/=0)fi_r(ij,1:nf)=fi_r(ij,1:nf)+eigen_vect(id2,1:nf)*norm(2)
 if(id3/=0)fi_r(ij,1:nf)=fi_r(ij,1:nf)+eigen_vect(id3,1:nf)*norm(3)
 enddo
!
! DRM-3D: Modification of xy-vect to include the modal acceleration
!
 xy_vect(1:nq)=xy_vect(1:nq)+rho_w/gravity*(w)**2*&
 & matmul(matmul(G_mat(1:nq,j,1:nod),fi_r(1:nod,1:nf)),L(1:nf))
 endif
 enddo
!
! DRM-3D: Solving for Ptot
!
 xy_vect=matmul(a_mat,xy_vect)
!
! DRM3D: Re-assembling {xy} --> {p},{q}
!
 where(node(lmw(1:nq))%code==0)node(lmw(1:nq))%p=node(lmw(1:nq))%p+xy_vect
 where(node(lmw(1:nq))%code==2)node(lmw(1:nq))%p=node(lmw(1:nq))%p+xy_vect
 where(node(lmw(1:nq))%code==3)node(lmw(1:nq))%p=node(lmw(1:nq))%p+xy_vect
 where(node(lmw(1:nq))%code==4)node(lmw(1:nq))%p=node(lmw(1:nq))%p+xy_vect

 do i=1,nelm
 n(1:nod)=element(i)%node(1:nod) 
 n(1:nod)=node(n(1:nod))%idw 
 do ij=1,nod
 select case (element(i)%code)
 case(1); element(i)%q(ij)=xy_vect(n(ij))
 case(2); element(i)%p(ij)=node(lmw(n(ij)))%p
 case(3); element(i)%p(ij)=node(lmw(n(ij)))%p
!          element(i)%q(ij)=element(i)%q(ij)-(0,1)*(w/sound)*gama*node(lmw(n(ij)))%p
 case(4); element(i)%p(ij)=node(lmw(n(ij)))%p 
 end select
 enddo
 enddo
!
! DRM-3D: Hydrodynamic Force on the Upstream face of the dam
!
!
 allocate(funct(nod),derv_eta(3,nod))
 allocate(coord(nod,3),derv_x(3,nod))
 allocate(N_funct(3,nod),N_Der(3,nod))  
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 Force=0.0
 Force_static=0.0
 do j=1,nelm
 if(element(j)%code==2)then
 n(1:nod)=element(j)%node(1:nod)
 do ij=1,nod; coord(ij,:)=node(n(ij))%coord(:); enddo
 Po(1:nod)=gravity*node(n(1:nod))%p
 Po_static(1:nod)=rho_w*node(n(1:nod))%coord(3)
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(nod,eta(ii),eta(jj))
 Force=Force+ww(ii)*ww(jj)*det_jac*sum(funct(1:nod)*Po(1:nod))*norm
 Force_static=Force_static+ww(ii)*ww(jj)*det_jac*sum(funct(1:nod)*&
              &Po_static(1:nod))*norm
 enddo
 enddo
!
 endif
 enddo
!
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta,ww)
!
 deallocate(a_mat,xy_vect)
 deallocate(Bmod,fi,fi_r)
 deallocate(Po,Po_static,Pj,Qo,Qj)  
!
end subroutine HNodal
