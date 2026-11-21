subroutine DHistory
use dam3d
use modal
use graphics
implicit none
!
! Dam3D: Form the FEM [M],[K],[C] for the Dam Structure
!
 if(ndyn<=3)call Ssolid
!
! Dam3D: Dam Free Vibration EVP
!
 allocate(eigen_value(nf),eigen_vect(neq,nf))
 allocate(Po_(ntot),Pj_(ntot,nf))
 eigen_value=0.0; eigen_vect=0.0
 if(ndyn<=3)call SEVP
!
! Dam3D: Dam and Reservoir System Response
!
 if(ndyn>0)then
!
 allocate(m_mat(nf),c_mat(nf),k_mat(nf),L_vect(nf))
 allocate(r_g(neq,3))
 allocate(S(nf,nf),L(nf)) 
 allocate(Y(nf),YY(nf))
!
 if(ndyn<=3)call SModal
!
! Dam3D: Form the DRM [H],[G],[S]
!
 if(ndyn>=2)call Hfluid
!
! Dam3D: Form the FEM2D and EVP for the radiation Problem 
!
 if(ndyn>=2.and.nqr/=0)call Rrad
!
! Dam3D: Modal Analysis Including Dam-Water Interaction
!
 write(*,'(/a)')'Dam-3D: Modal Analysis'
!
 wo=pi*sound/(2.*hwater)  ! first frequency
 w=0
! write(*,'(/a,g12.4)')'  w/wr =',w/eigen_value(1) 
 do num=1,20
 w=num*0.1*wo
!
 write(*,'(//a,f14.5/)')' Dam-3D: Analysis for w=',w
 write(prn,'(//a,f14.5)')'w=',w
!w=pi/2
! w=0.0
! do w=wi*eigen_value(1),wf*eigen_value(1),wr*eigen_value(1)
! do w=wi*eigen_value(1),wf*eigen_value(1),wr*eigen_value(1)

!
! Modal: Form [S] & {L} for the required frequency range
!
 S=0.0; L=0.0; Y=0.0; YY=0.0
!
! Modal: Solve EVP for the Dam Structure and Construct [S] & {L}
!
 if(ndyn<=3)then
 do i=1,nf
 S(i,i)=-w**2*m_mat(i)+(0,1.)*w*c_mat(i)+k_mat(i)
 enddo
 L=L_vect
 endif
!
! Modal: Modification of Hydrodynamic Terms to [S] & {L}
!
   write(prn,'(/a/)')' L :::: Structure'
   write(prn,'(i8,2f20.7)')(i,L(i),i=1,nf)
 if(ndyn>=2)call HModal
!
! Modal: Solving [S]{Y}={L} System Using IMSL Library
!
   write(prn,'(/a/)')' L :::: Structure + Hydrodynamic'
   write(prn,'(i8,2f20.7)')(i,L(i),i=1,nf)
 if(ndyn>0.and.ndyn<=3)call dlslcg (nf,S,nf,L,1,L)
   write(prn,'(/a/)')' L :::: Result'
   write(prn,'(i8,2f20.7)')(i,L(i),i=1,nf)
!
! Nodal: Modal Analysis of the Nodal Structural Response
!
 if(ndyn<=3)call SNodal
!
! Nodal: Modal Analysis of the Nodal Hydrodynamic Response
!
 if(ndyn>=2)call HNodal
!
 if(Gmode==1.and.ndyn>0)call GOutput
 enddo

 deallocate(m_mat,c_mat,k_mat,L_vect)
 deallocate(r_g)
 deallocate(S,L,Y,YY)
 endif
!
 deallocate(Po_,Pj_)  
 deallocate(eigen_value,eigen_vect)
!
end subroutine DHistory
