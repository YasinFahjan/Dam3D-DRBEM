subroutine Rrad
use Rad3D
use Bem3D
implicit none
!
! Analytical\Rad: Closed Form Analytical Radiation Condition 
!                 for Rectangular Cross Section
!
 print*,'DRM-3D: Analytical Radiation Condition'
!
! Analytical\Rad: Constructing Eigen Values and Vectors
! 
 allocate(eval(nqr),evect(nqr,nqr))
 allocate(evalz(nqr),evaly(nqr))
 eval=0.0;evalz=0.0;evaly=0.0; evect=0.0
!
 i=0
 do ii=1,nqr
 do jj=1,nk1
 if(nod==8.and.mod(ii,2)==0.and.mod(jj,2)==0)cycle
 i=i+1
 if(i<=nqr)then
 evalz(i)=pi*(2*ii-1)/(2*hh)
 evaly(i)=pi*(jj-1)/bb
 eval(i)=csqrt(evalz(i)**2+evaly(i)**2)
 endif
 enddo; enddo

 do i=1,nqr
 do j=1,nqr
 evect(i,j)=cos(evalz(j)*node(lmr(i+nk1))%coord(3))&
          &*cos(evaly(j)*node(lmr(i+nk1))%coord(2))
 enddo; enddo
!
 write(prn,*)eval; write(prn,*)
!
end subroutine Rrad