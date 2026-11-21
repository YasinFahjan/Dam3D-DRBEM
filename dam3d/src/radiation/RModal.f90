subroutine RModal
use Rad3D
use Bem3D
implicit none
!
! DRM2D\Rad: Eigen Values and Vectors for Radiation Model
! 
 allocate(H_rad(nqr,nqr),S_rad(nqr,nqr),D_rad(neqr))
 allocate(G_rad(neqr,neqr))
 allocate(G11(nk1,nk1))
 H_rad=0.0; G_rad=0.0; S_rad=0.0; D_rad=0.0; G11=0.0
!
! DRM2D\Rad: Building G_rad from Grad
!
 do j=1,nelr
 n(1:ndr)=elementr(j)%node(1:ndr)
 do ij=1,ndr
 if(node(n(ij))%idr/=0)then
 jj=node(n(ij))%idr
 G_rad(:,jj)=G_rad(:,jj)+Grad(:,j,ij)
 if(node(n(ij))%codr==3)D_rad(jj)=D_rad(jj)+elementr(j)%q(ij)
 endif
 enddo; enddo
!
 G11(1:nk1,1:nk1)=G_rad(1:nk1,1:nk1)
 if(nk1>0)call lincg(nk1,g11,nk1,g11,nk1)
!
! DRM2D\Rad: Building [H]-lamda[S]={d}
!
 H_rad(1:nqr,1:nqr)=Hrad(nk1+1:neqr,nk1+1:neqr)-(0.,1.)*w*gama&
 &  *G_rad(nk1+1:neqr,nk1+1:neqr)
 H_rad(1:nqr,1:nqr)= H_rad(1:nqr,1:nqr)-matmul(&
 &  matmul(G_rad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  (Hrad(1:nk1,nk1+1:neqr)&
 &  -(0.,1.)*w*gama*G_rad(1:nk1,nk1+1:neqr)))
!
 S_rad(1:nqr,1:nqr)=Srad(nk1+1:neqr,nk1+1:neqr)&
 &  -matmul(matmul(G_rad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  Srad(1:nk1,nk1+1:neqr))
!
 d_rad(nk1+1:neqr)=matmul((G_rad(nk1+1:neqr,nk1+1:neqr)&
 &  -matmul(matmul(G_rad(nk1+1:neqr,1:nk1),G11(1:nk1,1:nk1)),&
 &  G_rad(1:nk1,nk1+1:neqr))),d_rad(nk1+1:neqr))
!                                                                     
! DRM2D\Rad: Constructing DRM2D EVP 
! 
 allocate(eval(nqr),evect(nqr,nqr))
 eval=0.0; evect=0.0
! 
 call lincg(nqr,s_rad,nqr,s_rad,nqr)
 h_rad=matmul(s_rad,h_rad)
 call evccg(nqr,h_rad,nqr,eval,evect,nqr)
!
 eval=csqrt(eval)
! evect=transpose(evect)
!
 print*,eval; print*,' '; write(prn,*)eval; write(prn,*)
!
! DRM2D\Rad: Building RAD Matrix
!
 print*,'DRM-3D: Modal Radiation conditions '
 allocate(Rad(nqr,nqr),G4(nq,nqr))
 allocate(Knm(nqr,nqr))
 Rad=0.0; G4=0.0; Knm=0.0
!
 do i=1,nqr
 Knm(i,i)=csqrt(eval(i)**2-(w)**2)
 enddo
 call lincg (nqr,evect,nqr,Rad,nqr)

      write(prn,'(/a/)')' Knm :::: '
      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
      write(prn,'(i8,<2*nqr>f20.7)')(i,(Knm(i,j),j=1,nqr),i=1,nqr)
      write(prn,'(/a/)')' evect :::: '
      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
      write(prn,'(i8,<2*nqr>f20.7)')(i,(evect(i,j),j=1,nqr),i=1,nqr)

 Rad=matmul(matmul(evect,Knm),Rad)

      write(prn,'(/a/)')' Rad :::: '
      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
      write(prn,'(i8,<2*nqr>f20.7)')(i,(Rad(i,j),j=1,nqr),i=1,nqr)


!
! DRM2D\Rad: Building G4 from G_mat
!
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod)
 do ij=1,nod
 if(node(n(ij))%idr/=0)then
 jj=node(n(ij))%idr
 G4(:,jj)=G4(:,jj)+G_mat(:,j,ij)
 endif
 enddo; enddo
!
! DRM2D\Rad: Modification to DRM System [a]{xy}={xy}
!
 G4=matmul(G4,Rad)
 do i=1,nq
 if(node(i)%idr/=0)then
 ii=node(i)%idw
 jj=node(i)%idr
 a_mat(:,ii)=a_mat(:,ii)-G4(:,jj)
 endif
 enddo
! 
 !deallocate(eval,evect)
 !deallocate(Rad,G4,Knm)
 !deallocate(H_rad,S_rad,D_rad)
 !deallocate(G_rad,G11)
!
end subroutine RModal
