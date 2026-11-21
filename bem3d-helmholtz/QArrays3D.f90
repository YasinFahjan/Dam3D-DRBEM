module arrays
! Quadratic BEM3D  
 use msimsl
 use portlib
  integer,        parameter :: max_nod=8       ! max # of nodes per element
  integer,        parameter :: n_dim=3	     ! Problem Dimensionality
  integer,        parameter :: ng=16            ! # of g-quadrature points 
  integer,        parameter :: inp=11, prn=12  ! Input & Output Files #     

! Definition the derived type structure node & element 
 type node_record
  real                      :: x(n_dim)	     ! Geom. Coord. of the node
  complex                      :: p      	     ! potential at the node
  complex                      :: q(n_dim)	     ! Internal points flux 
  integer                   :: code 	     ! indicates wether p is known
 end type node_record
  type(node_record),allocatable,   save:: node(:) 
! 
 type element_record
  integer                   :: node(max_nod)   ! element nodes 
  integer                   :: code		     ! Element code(1 Pot.,2 Flux)
  complex                      :: q(max_nod)	     ! normal fluxes at the element nodes 
  complex                      :: p(max_nod)	     ! Potential at the middle point 
 end type element_record
  type(element_record),allocatable,  save:: element(:) 

! Definition of the variables in the program 
  real,                save :: pi
  real,                save :: w
  integer,             save :: nnod,nelm,nint  ! # nodes, elements,int.pts
  integer,             save :: neq             ! # equations
  integer,             save :: nod             ! # of nodes per element
  complex, allocatable,   save :: a_mat(:,:)	     ! system matrix
  complex, allocatable,   save :: xy_vect(:)      ! system right vector
  real, allocatable,   save :: eta(:)          ! Numerical Integration points
  real, allocatable,   save :: ww(:)           ! Numerical Integration weights
  complex                      :: g(max_nod)      ! g matrix element
  complex                      :: h(max_nod)      ! h matrix element
  complex                      :: uhat(max_nod)         ! g matrix element
  complex                      :: qhat(max_nod)         ! h matrix element
  complex                      :: dgdz(max_nod,n_dim)   ! 1st deriv. of g w.r.t zeta
  complex                      :: dhdz(max_nod,n_dim)   ! 1st deriv. of h w.r.t zeta
  real,		     save :: fi(max_nod)           ! shape function
  real,		     save :: fi_eta(n_dim,max_nod) ! deriv. fi w.r.t eta
  real 		          :: zeta(n_dim)           ! source point coord.
  real 		          :: xx(n_dim)	           ! Target point coord.
  real		          :: jac_mat(n_dim,n_dim)  ! Jacob. matrix
  real		          :: coord_x(max_nod,n_dim) ! shape function
  real 		          :: norm(n_dim)     ! the normal at element   
  real 		          :: drdx(n_dim)     ! 1st deriv. of r w.r.t x
  real 		          :: drdnx(n_dim)    ! drdn*(drdx)
  character                 :: file*12,title*72  
  integer                   :: i,ii,j,jj,ij,code,n(max_nod)              
  real                      :: value(max_nod), det_jac 
  contains
				     
! Definition of the green functions and derivatives for 3D-Helmholtz
complex function green(r)
  green=1./(4.*pi*r)*cexp((0.,-1.)*w*r)
!   green=1./(4.*pi*r)*cos(w*r)
  end function green
complex function dgdr(r)
 dgdr=-1./(4.*pi*r)*((0.,1.)*w+1./r)*cexp((0.,-1.)*w*r)
!  dgdr=-1./(4.*pi*r)*(w*sin(w*r)+cos(w*r)/r)
  end function dgdr
complex function d2gdr(r)
  d2gdr=1./(4.*pi*r)*(w**2+(0.,1.)*w/r-2./r**2)*cexp((0.,-1.)*w*r)
!  d2gdr=sin(w*r)
  end function d2gdr




! Allocating data & working  arrays  
 subroutine arrays_allocation
  allocate(node(nnod+nint),element(nelm))
  allocate(a_mat(neq,neq),xy_vect(neq))
  allocate(eta(ng),ww(ng))
 end subroutine arrays_allocation

! Quadrilateral 4, 8 node quadratic shape functions 

 subroutine ShapeFunctions(eta1,eta2)
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

! Cross Product of Two vectors
 subroutine cross_product(a,b,c)
 real :: a(3),b(3),c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
 end subroutine cross_product

end module arrays
