subroutine SEVP
use FEM3D
use modal
implicit none
!
! FEM3D\EVP: Smallest Eigenvalues & Vectors for the system
! [M]{v}=1/(e_v)**2*[K]*{v}  Using IMSL library
 print*,'Dam-3D: EVP - Free Vibration '
!
 allocate(stiff_mass(neq,neq))
!
 call linrg(neq,stiff_mat,neq,stiff_mass,neq)
 stiff_mass=matmul(mass_mat,stiff_mass)
!call gvcsp (neq,mass_mat,neq,stiff_mat,neq,eigen_value,eigen_vect,neq)
 call evesf(neq,nf_f,stiff_mass,neq,.false.,eigen_value,eigen_vect,neq)
 eigen_value(1:nf)=1./sqrt(eigen_value(nf_i:nf_f:nf_s))
 eigen_vect(:,1:nf)=eigen_vect(:,nf_i:nf_f:nf_s)
! eigen_value(1:nf)=sqrt(abs(eigen_value(1:nf)))
!
! FEM3D\EVP: Print Eigen values 
!
 write(prn,'(//10x,a//a//a/a/)')&
 &' free vibration of the dam structure',&
 &' e i g e n v a l u e s   a n d   f r e q u e n c i e s',&
 &'   mode     eigenvalue  circular freq      frequency      period',&
 &' number   (rad/sec)**2      (rad/sec)   (cycles/sec)       (sec)'

 do i=1,nf
! betaa=ww*h*sqrt(rho(1)/(ee(1)*grav))
 write(prn,'(i6,5x,4(f10.4,5x))')i,(eigen_value(i))**2,eigen_value(i),&
 & eigen_value(i)/(2*pi),(2*pi)/eigen_value(i)
 enddo
!
 deallocate(stiff_mass)
!
 if(run_type==1)then
 open(1001,file='eigen.bin',form='binary',status='unknown',action='write')
 write(1001)eigen_value
 write(1001)eigen_vect
 close(1001)
 endif
!
end subroutine SEVP
