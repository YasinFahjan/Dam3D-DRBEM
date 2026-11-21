module shape_functions
! Quadratic DRM3D  
! Definition of the variables in the module 
  integer,             save :: nod          ! # of nodes per element
  real,allocatable,    save :: fi(:)        ! shape function
  real,allocatable,    save :: fi_eta(:,:)  ! deriv. fi w.r.t eta
  real,allocatable,    save :: eta_c(:,:)   ! source point coord.
  real,allocatable,    save :: jac_mat(:,:) ! Jacob. matrix
  real,allocatable,    save :: coord_x(:,:) ! shape function
  real,allocatable,    save :: norm(:)      ! normal to element   
  real,		     save :: det_jac 
  contains

 subroutine CornerFunctions
! Quadrilateral 8 node quadratic Corner functions 
 eta_c=0.0
 if(nod==4)then
 eta_c(1,1)=-1.; eta_c(1,2)=+1.; eta_c(1,3)=+1.; eta_c(1,4)=-1.
 eta_c(2,1)=-1.; eta_c(2,2)=-1.; eta_c(2,3)=+1.; eta_c(2,4)=+1.
 elseif(nod==8)then
 eta_c(1,1)=-1.; eta_c(1,2)=+1.; eta_c(1,3)=+1.; eta_c(1,4)=-1.
 eta_c(1,5)=0.0; eta_c(1,6)=+1.; eta_c(1,7)=0.0; eta_c(1,8)=-1.
 eta_c(2,1)=-1.; eta_c(2,2)=-1.; eta_c(2,3)=+1.; eta_c(2,4)=+1.
 eta_c(2,5)=-1.; eta_c(2,6)=0.0; eta_c(2,7)=+1.; eta_c(2,8)=0.0
 endif
 end subroutine CornerFunctions

 subroutine ShapeFunctions(eta1,eta2)
! Quadrilateral 8 node quadratic shape functions 
  fi=0.0; fi_eta=0.0
  if(nod==4)then
  fi(1)=.25*(1-eta1)*(1-eta2)
  fi(2)=.25*(1+eta1)*(1-eta2)
  fi(3)=.25*(1+eta1)*(1+eta2)
  fi(4)=.25*(1-eta1)*(1+eta2)
  fi_eta(1,1)=-.25*(1-eta2)
  fi_eta(1,2)= .25*(1-eta2)
  fi_eta(1,3)= .25*(1+eta2)
  fi_eta(1,4)=-.25*(1+eta2)
  fi_eta(2,1)=-.25*(1-eta1)
  fi_eta(2,2)=-.25*(1+eta1)
  fi_eta(2,3)= .25*(1+eta1)
  fi_eta(2,4)= .25*(1-eta1)
 elseif(nod==8)then
  fi(1)= .25*(1.-eta1)*(1.-eta2)*(-eta1-eta2-1.) 
  fi(2)= .25*(1.+eta1)*(1.-eta2)*( eta1-eta2-1.) 
  fi(3)= .25*(1.+eta1)*(1.+eta2)*( eta1+eta2-1.) 
  fi(4)= .25*(1.-eta1)*(1.+eta2)*(-eta1+eta2-1.) 
  fi(5)= .50*(1.-eta1)*(1.+eta1)*(1.-eta2)             
  fi(6)= .50*(1.+eta1)*(1.-eta2)*(1.+eta2)             
  fi(7)= .50*(1.-eta1)*(1.+eta1)*(1.+eta2)             
  fi(8)= .50*(1.-eta1)*(1.-eta2)*(1.+eta2)             
!
  fi_eta(1,1)= .25*(1.-eta2)*(2.*eta1+eta2)
  fi_eta(1,2)= .25*(1.-eta2)*(2.*eta1-eta2)
  fi_eta(1,3)= .25*(1.+eta2)*(2.*eta1+eta2)
  fi_eta(1,4)= .25*(1.+eta2)*(2.*eta1-eta2)
  fi_eta(1,5)= -eta1*(1.-eta2)
  fi_eta(1,6)= .5*(1.-eta2**2)
  fi_eta(1,7)= -eta1*(1.+eta2)
  fi_eta(1,8)=-.5*(1.-eta2**2)
!
  fi_eta(2,1)= .25*(1.-eta1)*(2.*eta2+eta1)
  fi_eta(2,2)= .25*(1.+eta1)*(2.*eta2-eta1)
  fi_eta(2,3)= .25*(1.+eta1)*(2.*eta2+eta1)
  fi_eta(2,4)= .25*(1.-eta1)*(2.*eta2-eta1)
  fi_eta(2,5)=-.5*(1.-eta1**2)
  fi_eta(2,6)=-eta2*(1.+eta1)
  fi_eta(2,7)= .5*(1.-eta1**2)
  fi_eta(2,8)=-eta2*(1.-eta1)
!
 endif
 end subroutine ShapeFunctions
!
 subroutine norm_jac
! Jacobian Matrix & the Normal to the Element Surface
 jac_mat=matmul(fi_eta,coord_x)
 call cross_product(jac_mat(1,:),jac_mat(2,:),norm)
 det_jac=sqrt(sum(norm*norm)); norm=norm/det_jac
!
 end subroutine norm_jac

! Cross Product of Two vectors
 subroutine cross_product(a,b,c)
 real :: a(3),b(3),c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
 end subroutine cross_product

end module shape_functions
