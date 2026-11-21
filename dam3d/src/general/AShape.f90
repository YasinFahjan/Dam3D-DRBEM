module shape_functions
implicit none
!
! Dam3D: Module Shape functions used in FEM and BEM3D/ models
!
 real,allocatable,      save :: funct(:)        ! DAM3D/shape function
 real,allocatable,      save :: derv_eta(:,:)   ! DAM3D/(df/deta)
 real,allocatable,      save :: derv_x(:,:)     ! DAM3D/(df/dx)
 real,allocatable,      save :: eta_c(:,:)      ! DAM3D/corner nodes
 real,allocatable,      save :: coord(:,:)      ! DAM3D/nodes coord.
 real,allocatable,      save :: N_funct(:,:)	! DAM3D/shape function 
 real,allocatable,      save :: N_Der(:,:)     	! DAM3D/Der.shape funct.
 real,                  save :: jac_mat(3,3)    ! DAM3D/Jacob. matrix
 real,                  save :: jac_inv(3,3)    ! DAM3D/Jacob. inverse
 real,                  save :: norm(3)         ! BEM3D/normal to element   
 real,		      save :: det_jac 		! DAM3D/Jacob. det.
!
contains

subroutine CornerFunctions(nod)
 integer nod
!
! BEM3D: Quadrilateral 8 node quadratic Corner functions 
!
 eta_c=0.0
 if(nod==2)then
 eta_c(1,1)=-1.; eta_c(1,2)=+1.
 elseif(nod==3)then
 eta_c(1,1)=-1.; eta_c(1,2)=0.0; eta_c(1,3)=+1.
 elseif(nod==4)then
 eta_c(1,1)=-1.; eta_c(1,2)=+1.; eta_c(1,3)=+1.; eta_c(1,4)=-1.
 eta_c(2,1)=-1.; eta_c(2,2)=-1.; eta_c(2,3)=+1.; eta_c(2,4)=+1.
 elseif(nod==8)then
 eta_c(1,1)=-1.; eta_c(1,2)=+1.; eta_c(1,3)=+1.; eta_c(1,4)=-1.
 eta_c(1,5)=0.0; eta_c(1,6)=+1.; eta_c(1,7)=0.0; eta_c(1,8)=-1.
 eta_c(2,1)=-1.; eta_c(2,2)=-1.; eta_c(2,3)=+1.; eta_c(2,4)=+1.
 eta_c(2,5)=-1.; eta_c(2,6)=0.0; eta_c(2,7)=+1.; eta_c(2,8)=0.0
 endif
end subroutine CornerFunctions

subroutine ShapeFunctions(nod,eta1,eta2)
 real eta1,eta2
 integer nod
!
! BEM3D: Quadrilateral 8 node quadratic shape functions 
!
 funct=0.0; derv_eta=0.0
 jac_inv=0.0;  N_Der=0.0; jac_mat=0.0; N_funct=0.0
 if(nod==2)then
 funct(1)=.5*(1-eta1)
 funct(2)=.5*(1+eta1)
 derv_eta(1,1)=-.5
 derv_eta(1,2)= .5
 elseif(nod==3)then
 funct(1)=-.5*eta1*(1-eta1)
 funct(2)=(1+eta1)*(1-eta1)
 funct(3)=.5*eta1*(1+eta1)
 derv_eta(1,1)=-.5*(1.-2.*eta1)
 derv_eta(1,2)=-2*eta1
 derv_eta(1,3)=.5*(1.+2.*eta1)
 elseif(nod==4)then
 funct(1)=.25*(1-eta1)*(1-eta2)
 funct(2)=.25*(1+eta1)*(1-eta2)
 funct(3)=.25*(1+eta1)*(1+eta2)
 funct(4)=.25*(1-eta1)*(1+eta2)
 derv_eta(1,1)=-.25*(1-eta2)
 derv_eta(1,2)= .25*(1-eta2)
 derv_eta(1,3)= .25*(1+eta2)
 derv_eta(1,4)=-.25*(1+eta2)
 derv_eta(2,1)=-.25*(1-eta1)
 derv_eta(2,2)=-.25*(1+eta1)
 derv_eta(2,3)= .25*(1+eta1)
 derv_eta(2,4)= .25*(1-eta1)
 elseif(nod==8)then
 funct(1)= .25*(1.-eta1)*(1.-eta2)*(-eta1-eta2-1.) 
 funct(2)= .25*(1.+eta1)*(1.-eta2)*( eta1-eta2-1.) 
 funct(3)= .25*(1.+eta1)*(1.+eta2)*( eta1+eta2-1.) 
 funct(4)= .25*(1.-eta1)*(1.+eta2)*(-eta1+eta2-1.) 
 funct(5)= .50*(1.-eta1)*(1.+eta1)*(1.-eta2)             
 funct(6)= .50*(1.+eta1)*(1.+eta2)*(1.-eta2)             
 funct(7)= .50*(1.+eta1)*(1.-eta1)*(1.+eta2)             
 funct(8)= .50*(1.-eta1)*(1.-eta2)*(1.+eta2)             
!
 derv_eta(1,1)= .25*(1.-eta2)*(2.*eta1+eta2)
 derv_eta(1,2)= .25*(1.-eta2)*(2.*eta1-eta2)
 derv_eta(1,3)= .25*(1.+eta2)*(2.*eta1+eta2)
 derv_eta(1,4)= .25*(1.+eta2)*(2.*eta1-eta2)
 derv_eta(1,5)= -eta1*(1.-eta2)
 derv_eta(1,6)= .5*(1.-eta2**2)
 derv_eta(1,7)= -eta1*(1.+eta2)
 derv_eta(1,8)=-.5*(1.-eta2**2)
!
 derv_eta(2,1)= .25*(1.-eta1)*(2.*eta2+eta1)
 derv_eta(2,2)= .25*(1.+eta1)*(2.*eta2-eta1)
 derv_eta(2,3)= .25*(1.+eta1)*(2.*eta2+eta1)
 derv_eta(2,4)= .25*(1.-eta1)*(2.*eta2-eta1)
 derv_eta(2,5)=-.5*(1.-eta1**2)
 derv_eta(2,6)=-eta2*(1.+eta1)
 derv_eta(2,7)= .5*(1.-eta1**2)
 derv_eta(2,8)=-eta2*(1.-eta1)
 endif
!
! BEM3D: Jacobian Matrix &  Normal to the Element Surface
!
 N_funct(1,1:nod)=funct(1:nod)
 jac_mat=matmul(derv_eta,coord)
 if(nod==2.or.nod==3)jac_mat(2,:)=(/1.,0.,0./)
 call cross_product(jac_mat(1,:),jac_mat(2,:),norm)
 det_jac=sqrt(sum(norm*norm)); norm=norm/det_jac
!
end subroutine ShapeFunctions

subroutine shell_funct(eta1,eta2,eta3)
!                                                                     
! FEM3D: shape functions for isoparametric 3D thick shell element                                                    
!                                                                     
 integer i,j,k
 real eta1,eta2,eta3
 real iperm(3)	      ! permutation to  [J]-1
 funct=0.0; derv_eta=0.0; derv_x=0.0
 jac_inv=0.0;  N_Der=0.0; jac_mat=0.0; N_funct=0.0
 funct( 1)= .125*(1.+eta1)*(1.+eta2)*(1.-eta3)*( eta1+eta2-1.)                    
 funct( 2)= .125*(1.-eta1)*(1.+eta2)*(1.-eta3)*(-eta1+eta2-1.)                           
 funct( 3)= .125*(1.-eta1)*(1.-eta2)*(1.-eta3)*(-eta1-eta2-1.)                           
 funct( 4)= .125*(1.+eta1)*(1.-eta2)*(1.-eta3)*( eta1-eta2-1.)                    
 funct( 5)= .125*(1.+eta1)*(1.+eta2)*(1.+eta3)*( eta1+eta2-1.)                    
 funct( 6)= .125*(1.-eta1)*(1.+eta2)*(1.+eta3)*(-eta1+eta2-1.)                           
 funct( 7)= .125*(1.-eta1)*(1.-eta2)*(1.+eta3)*(-eta1-eta2-1.)                           
 funct( 8)= .125*(1.+eta1)*(1.-eta2)*(1.+eta3)*( eta1-eta2-1.)                    
 funct( 9)= .25*(1.-eta1**2)*(1.+eta2)*(1.-eta3)                      
 funct(10)= .25*(1.-eta1)*(1.-eta2**2)*(1.-eta3)                             
 funct(11)= .25*(1.-eta1**2)*(1.-eta2)*(1.-eta3)                      
 funct(12)= .25*(1.+eta1)*(1.-eta2**2)*(1.-eta3)                      
 funct(13)= .25*(1.-eta1**2)*(1.+eta2)*(1.+eta3)                      
 funct(14)= .25*(1.-eta1)*(1.-eta2**2)*(1.+eta3)                             
 funct(15)= .25*(1.-eta1**2)*(1.-eta2)*(1.+eta3)                      
 funct(16)= .25*(1.+eta1)*(1.-eta2**2)*(1.+eta3)
! funct(17)= eta1*(1-eta1**2)       
! funct(18)= eta2*(1-eta2**2)                      
! funct(19)= (1-eta3**2)                      
! funct(20)= eta1*eta2*(1-eta3**2)                      
! funct(21)= eta2*eta1*(1-eta2**2)                      
!                                                                     
! derivatives of shape functions  w.r.t eta1,eta2,eta3                                    
!
 derv_eta(1, 1)= .125*(1.+eta2)*(1.-eta3)*( 2.*eta1+eta2)                           
 derv_eta(1, 2)=-.125*(1.+eta2)*(1.-eta3)*(-2.*eta1+eta2)                          
 derv_eta(1, 3)=-.125*(1.-eta2)*(1.-eta3)*(-2.*eta1-eta2)                          
 derv_eta(1, 4)= .125*(1.-eta2)*(1.-eta3)*( 2.*eta1-eta2)                          
 derv_eta(1, 5)= .125*(1.+eta2)*(1.+eta3)*( 2.*eta1+eta2)                           
 derv_eta(1, 6)=-.125*(1.+eta2)*(1.+eta3)*(-2.*eta1+eta2)                          
 derv_eta(1, 7)=-.125*(1.-eta2)*(1.+eta3)*(-2.*eta1-eta2)                          
 derv_eta(1, 8)= .125*(1.-eta2)*(1.+eta3)*( 2.*eta1-eta2)                          
 derv_eta(1, 9)=-.5*eta1*(1.+eta2)*(1.-eta3)                                    
 derv_eta(1,10)=-.25*(1.-eta2**2)*(1.-eta3)                           
 derv_eta(1,11)=-.5*eta1*(1.-eta2)*(1.-eta3)                                   
 derv_eta(1,12)= .25*(1.-eta2**2)*(1.-eta3)                           
 derv_eta(1,13)=-.5*eta1*(1.+eta2)*(1.+eta3)                                   
 derv_eta(1,14)=-.25*(1.-eta2**2)*(1.+eta3)                           
 derv_eta(1,15)=-.5*eta1*(1.-eta2)*(1.+eta3)                                   
 derv_eta(1,16)= .25*(1.-eta2**2)*(1.+eta3)                           
 derv_eta(1,17)= 1.-3.*eta1**2                                               
 derv_eta(1,18)= 0.                                                          
 derv_eta(1,19)= 0.                                                          
 derv_eta(1,20)= eta2*(1.-3.*eta1**2)                                        
 derv_eta(1,21)= eta2*(1.-eta2**2)                                            
!                                                                     
 derv_eta(2, 1)= .125*(1.+eta1)*(1.-eta3)*( 2.*eta2+eta1)                        
 derv_eta(2, 2)= .125*(1.-eta1)*(1.-eta3)*( 2.*eta2-eta1)                               
 derv_eta(2, 3)=-.125*(1.-eta1)*(1.-eta3)*(-2.*eta2-eta1)                               
 derv_eta(2, 4)=-.125*(1.+eta1)*(1.-eta3)*(-2.*eta2+eta1)                         
 derv_eta(2, 5)= .125*(1.+eta1)*(1.+eta3)*( 2.*eta2+eta1)                        
 derv_eta(2, 6)= .125*(1.-eta1)*(1.+eta3)*( 2.*eta2-eta1)                               
 derv_eta(2, 7)=-.125*(1.-eta1)*(1.+eta3)*(-2.*eta2-eta1)                               
 derv_eta(2, 8)=-.125*(1.+eta1)*(1.+eta3)*(-2.*eta2+eta1)                         
 derv_eta(2, 9)= .25*(1.-eta1**2)*(1.-eta3)                            
 derv_eta(2,10)=-.5*eta2*(1.-eta1)*(1.-eta3)                                   
 derv_eta(2,11)=-.25*(1.-eta1**2)*(1.-eta3)                           
 derv_eta(2,12)=-.5*eta2*(1.+eta1)*(1.-eta3)                            
 derv_eta(2,13)= .25*(1.-eta1**2)*(1.+eta3)                           
 derv_eta(2,14)=-.5*eta2*(1.-eta1)*(1.+eta3)                                   
 derv_eta(2,15)=-.25*(1.-eta1**2)*(1.+eta3)                           
 derv_eta(2,16)=-.5*eta2*(1.+eta1)*(1.+eta3)                            
 derv_eta(2,17)= 0.                                                          
 derv_eta(2,18)= 1.-3.*eta2**2                                               
 derv_eta(2,19)= 0.                                                          
 derv_eta(2,20)= eta1*(1.-eta1**2)                                            
 derv_eta(2,21)= eta1*(1.-3.*eta2**2)                                        
!                                                                     
 derv_eta(3, 1)=-.125*(1.+eta1)*(1.+eta2)*( eta1+eta2-1.)                          
 derv_eta(3, 2)=-.125*(1.-eta1)*(1.+eta2)*(-eta1+eta2-1.)                                 
 derv_eta(3, 3)=-.125*(1.-eta1)*(1.-eta2)*(-eta1-eta2-1.)                                 
 derv_eta(3, 4)=-.125*(1.+eta1)*(1.-eta2)*( eta1-eta2-1.)                          
 derv_eta(3, 5)= .125*(1.+eta1)*(1.+eta2)*( eta1+eta2-1.)                          
 derv_eta(3, 6)= .125*(1.-eta1)*(1.+eta2)*(-eta1+eta2-1.)                                 
 derv_eta(3, 7)= .125*(1.-eta1)*(1.-eta2)*(-eta1-eta2-1.)                                 
 derv_eta(3, 8)= .125*(1.+eta1)*(1.-eta2)*( eta1-eta2-1.)                          
 derv_eta(3, 9)=-.25*(1.-eta1**2)*(1.+eta2)                            
 derv_eta(3,10)=-.25*(1.-eta1)*(1.-eta2**2)                                  
 derv_eta(3,11)=-.25*(1.-eta1**2)*(1.-eta2)                           
 derv_eta(3,12)=-.25*(1.+eta1)*(1.-eta2**2)                           
 derv_eta(3,13)= .25*(1.-eta1**2)*(1.+eta2)                           
 derv_eta(3,14)= .25*(1.-eta1)*(1.-eta2**2)                                  
 derv_eta(3,15)= .25*(1.-eta1**2)*(1.-eta2)                           
 derv_eta(3,16)= .25*(1.+eta1)*(1.-eta2**2)                           
 derv_eta(3,17)= 0.                                                          
 derv_eta(3,18)= 0.                                                          
 derv_eta(3,19)=-2.*eta3                                                      
 derv_eta(3,20)= 0.                                                         
 derv_eta(3,21)= 0.                                                         
!                                                                     
! FEM3D: jacobian matrix [J]                                                    
!                                                                     
 jac_mat(1:3,1:3)=matmul(derv_eta(1:3,1:16),coord(1:16,1:3))
!                                                                     
! FEM3D: jacobian inverse [J]-1                                                    
!                                                                     
 iperm=(/2,3,1/)
 do i=1,3                                                             
 j=iperm(i)                                                          
 k=iperm(j)                                                          
 jac_inv(i,i)=jac_mat(j,j)*jac_mat(k,k)-jac_mat(k,j)*jac_mat(j,k)                                  
 jac_inv(i,j)=jac_mat(k,j)*jac_mat(i,k)-jac_mat(i,j)*jac_mat(k,k)                                  
 jac_inv(j,i)=jac_mat(j,k)*jac_mat(k,i)-jac_mat(j,i)*jac_mat(k,k)                                  
 enddo
 det_jac=sum(jac_mat(1,1:3)*jac_inv(1:3,1))                       
 jac_inv(1:3,1:3)=jac_inv(1:3,1:3)/det_jac
!                                                                     
! FEM3D: derivatives of shape functions  w.r.t x,y,z                                    
!                                                                     
 derv_x(1:3,1:21)=matmul(jac_inv(1:3,1:3),derv_eta(1:3,1:21)) 
!                                                                     
end subroutine shell_funct

subroutine solid_funct(eta1,eta2,eta3)
!                                                                     
! FEM3D: shape functions for isoparametric 3D 8-node solid element                                                    
!                                                                     
 integer i,j,k
 real eta1,eta2,eta3
 real iperm(3)	      ! permutation to  [J]-1
 funct=0.0; derv_eta=0.0; derv_x=0.0
 jac_inv=0.0;  N_Der=0.0; jac_mat=0.0; N_funct=0.0
 funct( 1)= .125*(1.+eta1)*(1.+eta2)*(1.-eta3)                    
 funct( 2)= .125*(1.-eta1)*(1.+eta2)*(1.-eta3)                    
 funct( 3)= .125*(1.-eta1)*(1.-eta2)*(1.-eta3)                           
 funct( 4)= .125*(1.+eta1)*(1.-eta2)*(1.-eta3)                           
 funct( 5)= .125*(1.+eta1)*(1.+eta2)*(1.+eta3)                    
 funct( 6)= .125*(1.-eta1)*(1.+eta2)*(1.+eta3)                    
 funct( 7)= .125*(1.-eta1)*(1.-eta2)*(1.+eta3)                           
 funct( 8)= .125*(1.+eta1)*(1.-eta2)*(1.+eta3)                           
! funct( 9)=  1.-eta1**2
! funct(10)=  1.-eta2**2       
! funct(11)=  1.-eta3**2                      
!                                                                     
! derivatives of shape functions  w.r.t eta1,eta2,eta3                                    
!
 derv_eta(1, 1)= .125*(1.+eta2)*(1.-eta3)
 derv_eta(1, 2)=-.125*(1.+eta2)*(1.-eta3)
 derv_eta(1, 3)=-.125*(1.-eta2)*(1.-eta3)
 derv_eta(1, 4)= .125*(1.-eta2)*(1.-eta3)
 derv_eta(1, 5)= .125*(1.+eta2)*(1.+eta3)
 derv_eta(1, 6)=-.125*(1.+eta2)*(1.+eta3)
 derv_eta(1, 7)=-.125*(1.-eta2)*(1.+eta3)
 derv_eta(1, 8)= .125*(1.-eta2)*(1.+eta3)
 derv_eta(1, 9)=-2*eta1                                    
 derv_eta(1,10)= .0                            
 derv_eta(1,11)= .0                                   
!                                                                     
 derv_eta(2, 1)= .125*(1.+eta1)*(1.-eta3)
 derv_eta(2, 2)= .125*(1.-eta1)*(1.-eta3) 
 derv_eta(2, 3)=-.125*(1.-eta1)*(1.-eta3)       
 derv_eta(2, 4)=-.125*(1.+eta1)*(1.-eta3)       
 derv_eta(2, 5)= .125*(1.+eta1)*(1.+eta3)
 derv_eta(2, 6)= .125*(1.-eta1)*(1.+eta3) 
 derv_eta(2, 7)=-.125*(1.-eta1)*(1.+eta3)       
 derv_eta(2, 8)=-.125*(1.+eta1)*(1.+eta3)       
 derv_eta(2, 9)= .0                            
 derv_eta(2,10)=-2*eta2                                   
 derv_eta(2,11)= .0                           
!                                                                     
 derv_eta(3, 1)=-.125*(1.+eta1)*(1.+eta2)
 derv_eta(3, 2)=-.125*(1.-eta1)*(1.+eta2)
 derv_eta(3, 3)=-.125*(1.-eta1)*(1.-eta2)
 derv_eta(3, 4)=-.125*(1.+eta1)*(1.-eta2)
 derv_eta(3, 5)= .125*(1.+eta1)*(1.+eta2)
 derv_eta(3, 6)= .125*(1.-eta1)*(1.+eta2)
 derv_eta(3, 7)= .125*(1.-eta1)*(1.-eta2)
 derv_eta(3, 8)= .125*(1.+eta1)*(1.-eta2)
 derv_eta(3, 9)= .0                            
 derv_eta(3,10)= .0                                  
 derv_eta(3,11)=-2*eta3                           
!                                                                     
! FEM3D: jacobian matrix [J]                                                    
!                                                                     
 jac_mat(1:3,1:3)=matmul(derv_eta(1:3,1:8),coord(1:8,1:3))
!                                                                     
! FEM3D: jacobian inverse [J]-1                                                    
!                                                                     
 iperm=(/2,3,1/)
 do i=1,3                                                             
 j=iperm(i)                                                          
 k=iperm(j)                                                          
 jac_inv(i,i)=jac_mat(j,j)*jac_mat(k,k)-jac_mat(k,j)*jac_mat(j,k)                                  
 jac_inv(i,j)=jac_mat(k,j)*jac_mat(i,k)-jac_mat(i,j)*jac_mat(k,k)                                  
 jac_inv(j,i)=jac_mat(j,k)*jac_mat(k,i)-jac_mat(j,i)*jac_mat(k,k)                                  
 enddo
 det_jac=sum(jac_mat(1,1:3)*jac_inv(1:3,1))                       
 jac_inv(1:3,1:3)=jac_inv(1:3,1:3)/det_jac
!                                                                     
! FEM3D: derivatives of shape functions  w.r.t x,y,z                                    
!                                                                     
 derv_x(1:3,1:11)=matmul(jac_inv(1:3,1:3),derv_eta(1:3,1:11)) 
!                                                                     
end subroutine solid_funct
							    
subroutine cross_product(a,b,c)
! Cross Product of Two vectors
 real :: a(3),b(3),c(3)
 c(1)=a(2)*b(3)-a(3)*b(2)
 c(2)=a(3)*b(1)-a(1)*b(3)
 c(3)=a(1)*b(2)-a(2)*b(1)
end subroutine cross_product

end module shape_functions