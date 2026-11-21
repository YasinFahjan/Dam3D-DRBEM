subroutine RModal
use Rad3D
use Bem3D
implicit none
!
! Analytical\Rad: Building RAD Matrix
!
 print*,'DRM-3D: Modal Radiation conditions '
 allocate(Rad(nqr,nqr),G4(nq,nqr))
 allocate(Knm(nqr,nqr))
 Rad=0.0; G4=0.0; Knm=0.0
!      write(prn,'(/a/)')' evect :::: '
!      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
!      write(prn,'(i8,<2*nqr>f20.7)')(i,(evect(i,j),j=1,nqr),i=1,nqr)

!
 do i=1,nqr
 Knm(i,i)=csqrt(eval(i)**2-(w/sound)**2)
 enddo
 call dlincg(nqr,evect,nqr,Rad,nqr)

!      write(prn,'(/a/)')' Knm :::: '
!      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
!      write(prn,'(i8,<2*nqr>f20.7)')(i,(Knm(i,j),j=1,nqr),i=1,nqr)
Rad=matmul(matmul(evect,Knm),Rad)

!      write(prn,'(/a/)')' Rad :::: '
!      write(prn,'(4x,<2*nqr>i20/)')(i,i,i=1,nqr)
!      write(prn,'(i8,<2*nqr>f20.7)')(i,(Rad(i,j),j=1,nqr),i=1,nqr)
!
! Analytical\Rad: Building G4 from G_mat
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
! Analytical\Rad: Modification to DRM System [a]{xy}={xy}
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
 deallocate(Rad,G4,Knm)
!
end subroutine RModal
