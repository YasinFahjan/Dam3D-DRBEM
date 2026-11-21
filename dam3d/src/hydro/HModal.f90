Subroutine HModal
use Bem3D
use modal
implicit none
!
! DRM3D: Assembling the system [H]{p}-[G]{q}=-w^2*[S]{p}
!
 print*,'DRM-3D: Hydrodynamic Modal Analysis '
 allocate(a_mat(nq,nq),xy_vect(nq))
!
 a_mat=0.0; xy_vect=0.0
!
! DRM3D: Assembling [G]{q} --> [a]{xy}={xy}
!
 do i=1,nq
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod) 
 n(1:nod)=node(n(1:nod))%idw 
 do ij=1,nod
 select case (node(lmw(n(ij)))%code)
 case (1)
 select case(element(j)%code)
! case(1); a_mat(i,n(ij))=a_mat(i,n(ij))-g_mat(i,j,ij)
 case(1,4); a_mat(i,n(ij))=a_mat(i,n(ij))-g_mat(i,j,ij)
 case(2); xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
 case(3); xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
!          xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*(0.,1.)*(w/sound)*gama*node(lmw(n(ij)))%p
 end select
 case (2:4)
 select case(element(j)%code)
 case(2); xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
 case(3); xy_vect(i)=xy_vect(i)+g_mat(i,j,ij)*element(j)%q(ij)
          a_mat(i,n(ij))=a_mat(i,n(ij))-(0.,1.)*(w/sound)*gama*g_mat(i,j,ij)
 endselect
 end select
 enddo
 enddo
 enddo
!
! DRM3D: Assembling [H+w^2*S]{p} --> [a]{xy}={xy}
!
 do i=1,nq
 do j=1,nq
 select case (node(lmw(j))%code)
 case (1)
 xy_vect(i)=xy_vect(i)-(h_mat(i,j)+(w/sound)**2*s_mat(i,j))*node(lmw(j))%p
 case (0,2:4)
 a_mat(i,j)=a_mat(i,j)+h_mat(i,j)+(w/sound)**2*s_mat(i,j)
 end select
 enddo
 enddo
!
! DRM-3D: Formulation for  Po ( rigid dam system )
!
 node%p=0.0
!
 if(ndyn==2.or.ndyn==4)then
!
! DRM3D: Radiation condition Modification to [a]{xy}={xy}
!
 if(nqr>0.and.ixyz==1)call RModal
!
! DRM-3D: Solving for Po 
!
 call dlslcg (nq,a_mat,nq,xy_vect,1,xy_vect)
 where(node(lmw(1:nq))%code==0)node(lmw(1:nq))%p=xy_vect
 where(node(lmw(1:nq))%code==2)node(lmw(1:nq))%p=xy_vect
 where(node(lmw(1:nq))%code==3)node(lmw(1:nq))%p=xy_vect
 where(node(lmw(1:nq))%code==4)node(lmw(1:nq))%p=xy_vect
 endif
!
! DRM-3D: Formulation for Pj ( deformational dam motions)
!
 if(ndyn==2.or.ndyn==3)then
!
! DRM3D: Radiation condition Modification to [a]{xy}={xy}
!
 if(nqr>0)call RModal
!
! DRM-3D: Inverting the matrix [a] using IMSL library
!
 call dlincg (nq,a_mat,nq,a_mat,nq)
 endif
!
! DRM-3D: Modification of the Hydrodynamic Force to [S] and {L}
!
 allocate(Bmod(nod,nod),fi(nod,nf),fi_r(nod,nf))
 allocate(Po(nod),Po_static(nod),Pj(nod,nf),Qo(nod),Qj(nod,nf))
!
 allocate(funct(nod),derv_eta(3,nod))
 allocate(coord(nod,3),derv_x(3,nod))
 allocate(N_funct(3,nod),N_Der(3,nod))  
!
 allocate(eta(ng),ww(ng))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),eta,ww)
!
 Po_=0.0; Pj_=0.0 
!
 do j=1,nelm
 if(element(j)%code==2)then
!
 fi_r=0.0 
 Po=0.0; Pj=0.0 
 n(1:nod)=element(j)%node(1:nod)
 do ij=1,nod; coord(ij,:)=node(n(ij))%coord(:); enddo
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
! DRM-3D: Solving for Po, Pj
!
 Po(1:nod)=node(n(1:nod))%p
 Po_(n(1:nod))=Po(1:nod)
! 
 n(1:nod)=node(n(1:nod))%idw 
 Pj(1:nod,1:nf)=rho_w/gravity*matmul(matmul(a_mat(lmw(n(1:nod)),1:nq)&
                &,G_mat(1:nq,j,1:nod)),fi_r(1:nod,1:nf))
 n(1:nod)=element(j)%node(1:nod) 
 Pj_(n(1:nod),1:nf)=Pj(1:nod,1:nf)
!
! DRM-3D: Forming  x,y,z components for Qo and Qj
!
 do i_xyz=1,3
!
 fi=0.0
 do ij=1,nod
 idd=node(n(ij))%id(i_xyz)
 if(idd/=0)fi(ij,1:nf)=fi(ij,1:nf)+eigen_vect(idd,1:nf)
 enddo
!
! DRM-3D: form B-mod matrix ( pressure --> Force ) 
!
 Bmod=0.0; Qo=0.0; Qj=0.0
 do ii=1,ng
 do jj=1,ng
 call ShapeFunctions(nod,eta(ii),eta(jj))
 Bmod(1:nod,1:nod)=Bmod(1:nod,1:nod)+ww(ii)*ww(jj)*det_jac&
 &*matmul(transpose(N_funct),N_funct)*norm(i_xyz)
 enddo
 enddo
!
! DRM-3D: Solving for Qo, Qj
!
 Qo(1:nod)=matmul(Bmod(1:nod,1:nod),Po(1:nod))
 Qj(1:nod,1:nf)=matmul(Bmod(1:nod,1:nod),Pj(1:nod,1:nf))
!
! DRM-3D: Modification of [S] and {L}
!
 S(1:nf,1:nf)=S(1:nf,1:nf)-(w)**2*&
              & matmul(transpose(fi(1:nod,1:nf)),Qj(1:nod,1:nf))
 L(1:nf)=L(1:nf)-matmul(transpose(fi(1:nod,1:nf)),Qo(1:nod))
!
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
end Subroutine HModal

