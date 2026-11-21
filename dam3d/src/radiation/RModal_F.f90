subroutine RModal
use Rad3D
use Bem3D
implicit none
!                                                                     
! FEM2D\Rad: Constructing FEM2D EVP
! 
!
 if(sum(ag*(/1.,0.,0./))/=0.0)then
! 
 print*,'DRM-3D: Modal Radiation conditions '
 allocate(Rad(nqr,nqr))
 allocate(eval(nqr),evect(nqr,nqr))
 rad=0.0; eval=0.0; evect=0.0
!
 rad(1:nqr,1:nqr)=g_rad(1:nqr,1:nqr)
 call dlincg(nqr, rad, nqr, rad, nqr)
 rad=matmul(rad,(h_rad+(0.,1.0)*(w/sound)*gama*b_rad))
 call devccg(nqr, rad, nqr, eval, evect, nqr)
!
 eval=cdsqrt(eval)
! print*,eval; print*,' '
 write(prn,*)eval; write(prn,*)
!
! FEM2D\Rad: Building RAD Matrix
!
 allocate(G4(nq,nqr))
 allocate(Knm(nqr,nqr))
 Rad=0.0; G4=0.0; Knm=0.0
!
 do i=1,nqr
 Knm(i,i)=cdsqrt(eval(i)**2-(w/sound)**2)
 enddo
 call dlincg (nqr,evect,nqr,Rad,nqr)

 Rad=matmul(matmul(evect,Knm),Rad)
!
! FEM2D\Rad: Building G4 from G_mat
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
! FEM2D\Rad: Modification to DRM System [a]{xy}={xy}
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
!
 endif
!
end subroutine RModal
