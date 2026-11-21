subroutine SModal
use FEM3D
use modal
implicit none
!                                                                               
! Modal: Influence vecors rx,ry,rx                                                              
! 
 r_g=0.0
 do i =1,numnp                                                       
 do j =1,3                                                            
 if(node(i)%id(j)>0)r_g(node(i)%id(j),j)=1.0
 enddo; enddo                                                           
!
! FEM3D\Modal: Mass, Damp, Stiffness Matrices
!
 print*,'FEM-3D: Structural Modal Analysis '
!
 m_mat=0.; c_mat=0.; k_mat=0.

 do i=1,nf
 m_mat(i)=dot_product (eigen_vect(:,i),matmul(mass_mat,eigen_vect(:,i)))
 c_mat(i)=2*hyst_damp*eigen_value(i)*m_mat(i)
 k_mat(i)=dot_product (eigen_vect(:,i),matmul(stiff_mat,eigen_vect(:,i)))
 enddo
 L_vect=ag(ixyz)*matmul(transpose(eigen_vect(:,1:nf)),matmul(mass_mat,r_g(:,ixyz)))

! do j=1,3
! L=L-ag(ixyz)*matmul(transpose(eigen_vect(:,1:nf)),matmul(mass_mat,r_g(:,ixyz)))
! enddo
!
 deallocate(mass_mat,stiff_mat)
!
 if(run_type==1)then
 open(1002,file='solid.bin',form='binary',status='unknown',action='write')
 write(1002)m_mat
 write(1002)c_mat
 write(1002)k_mat
 write(1002)L_vect
 write(1002)r_g
 close(1002)
 endif
!
end subroutine SModal
