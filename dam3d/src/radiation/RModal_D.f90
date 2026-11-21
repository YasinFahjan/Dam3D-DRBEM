subroutine RModal
use Rad3D
use Bem3D
implicit none
!                                                                     
! DRM2D\Rad: Constructing DRM2D EVP 
! 
 write(*,'(a)')' DRM2D\Rad: Constructing DRM2D EVP '
!
 allocate(Rad(nqr,nqr))
 allocate(eval(nqr),evect(nqr,nqr))
 rad=0.0; eval=0.0; evect=0.0
!
! DRM2D\Rad: Inverse [S]
!
 rad(1:nqr,1:nqr)=-srad(1:nqr,1:nqr)
 call dlincg(nqr, rad, nqr, rad, nqr)
 rad(1:nqr,1:nqr)=matmul(rad(1:nqr,1:nqr),(hrad(1:nqr,1:nqr)&
 & +(0.,1.0)*(w/sound)*gama*grad(1:nqr,1:nqr)))
 call devccg(nqr, rad, nqr, eval, evect, nqr)
 eval=cdsqrt(eval)
! write(prn,*)eval; write(prn,*)
!
! DRM2D\Rad: Building RAD Matrix
!
 print*,'DRM-3D: Modal Radiation conditions '
 allocate(G4(nq,nqr),Knm(nqr,nqr))
 rad=0.0;G4=0.0; Knm=0.0
 do i=1,nqr
 Knm(i,i)=cdsqrt(eval(i)**2-(w/sound)**2)
 enddo
 call dlincg (nqr,evect,nqr,Rad,nqr)

 !    write(prn,'(/a/)')' Knm :::: '
 !    write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
 !    write(prn,'(i8,<2*nqr>f20.7)')(i,(Knm(i,j),j=1,nqr),i=1,nqr)
 !    write(prn,'(/a/)')' evect :::: '
 !    write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
 !    write(prn,'(i8,<2*nqr>f20.7)')(i,(evect(i,j),j=1,nqr),i=1,nqr)

 Rad=matmul(matmul(evect,Knm),Rad)

 !     write(prn,'(/a/)')' Rad :::: '
 !     write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
 !     write(prn,'(i8,<2*nqr>f20.7)')(i,(Rad(i,j),j=1,nqr),i=1,nqr)
!
! DRM2D\Rad: Building G4 from G_mat
!
 do j=1,nelm
 n(1:nod)=element(j)%node(1:nod)
 do ij=1,nod
 if(node(n(ij))%idr>nk1)then
 jj=node(n(ij))%idr-nk1
 G4(:,jj)=G4(:,jj)+G_mat(:,j,ij)
 endif
 enddo; enddo
!
! DRM2D\Rad: Modification to DRM System [a]{xy}={xy}
!
 G4=matmul(G4,Rad)
 do i=1,nq
 if(node(i)%idr>nk1)then
 ii=node(i)%idw
 jj=node(i)%idr-nk1
 a_mat(:,ii)=a_mat(:,ii)-G4(:,jj)
 endif
 enddo
! 
 deallocate(eval,evect)
 deallocate(Rad,G4,Knm)

!  write(prn,'(/a/)')' Hrad :::: '
!  write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!  write(prn,'(i8,<neqr>f20.7)')(i,(Hrad(i,j),j=1,neqr),i=1,neqr)

!  write(prn,'(/a/)')' G_rad :::: '
!  write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!  write(prn,'(i8,<neqr>f20.7)')(i,(G_rad(i,j),j=1,neqr),i=1,neqr)

!  write(prn,'(/a/)')' srad :::: '
!  write(prn,'(4x,<neqr>i20/)')(i,i=1,neqr)
!  write(prn,'(i8,<neqr>f20.7)')(i,(srad(i,j),j=1,neqr),i=1,neqr)

!  write(prn,'(/a/)')' g11 :::: '
!  write(prn,'(4x,<nk1>i20/)')(i,i=1,nk1)
!  write(prn,'(i8,<nk1>f20.7)')(i,(g11(i,j),j=1,nk1),i=1,nk1)

!  write(prn,'(/a/)')' g11 :::: '
!  write(prn,'(4x,<nk1>i20/)')(i,i=1,nk1)
!  write(prn,'(i8,<nk1>f20.7)')(i,(g11(i,j),j=1,nk1),i=1,nk1)

end subroutine RModal
