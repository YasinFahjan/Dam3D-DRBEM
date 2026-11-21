module arrays
! Linear BEM3D  
 use msimsl
 use portlib
  integer,        parameter :: max_nod=4       ! max # of nodes per element
  integer,        parameter :: n_dim=3	     ! Problem Dimensionality
  integer,        parameter :: ng=4           ! # of g-quadrature points 
  integer,        parameter :: inp=11, prn=12  ! Input & Output Files #     

! Definition the derived type structure node & element 
 type node_record
  real                      :: x(n_dim)	     ! Geom. Coord. of the node
  real                      :: p      	     ! potential at the node
  real                      :: q(n_dim)	     ! Internal points flux 
  integer                   :: code 	     ! indicates wether p is known
 end type node_record
  type(node_record),allocatable,   save:: node(:) 
! 
 type element_record
  integer                   :: node(max_nod)   ! element nodes 
  integer                   :: code		     ! Element code(1 Pot.,2 Flux)
  real                      :: q(max_nod)	     ! normal fluxes at the element nodes 
  real                      :: p(max_nod)	     ! Potential at the middle point 
 end type element_record
  type(element_record),allocatable,  save:: element(:) 

! Definition of the variables in the program 
  real,                save :: pi
  integer,             save :: nnod,nelm,nint   ! # nodes, elements,int.pts
  integer,             save :: nod              ! # of nodes per element
  integer,             save :: nod_eff          ! # of nodes per element
  real, allocatable,   save :: a_mat(:,:)	      ! system matrix
  real, allocatable,   save :: xy_vect(:)       ! system right vector
  real, allocatable,   save :: eta1(:),eta2(:)  ! Numerical Integration points
  real, allocatable,   save :: w(:),etafix(:)   ! Numerical Integration weights
  real                      :: g(max_nod)       ! g matrix element
  real                      :: h(max_nod)       ! h matrix element
  real                      :: dgdz(max_nod,n_dim) ! 1st deriv. of g w.r.t zeta
  real                      :: dhdz(max_nod,n_dim) ! 1st deriv. of h w.r.t zeta
  real,		     save :: fi(max_nod)	! shape function
  real 		          :: zeta(n_dim)	! source point coord.
  real 		          :: xx(n_dim)	      ! Target point coord.
  real 		          :: ds1(n_dim)       ! element surface vector 1
  real 		          :: ds2(n_dim)       ! element surface vector 2
  real 		          :: norm(n_dim)      ! the normal at element   
  real 		          :: drdx(n_dim)      ! 1st deriv. of r w.r.t x
  real 		          :: drdnx(n_dim)     ! drdn*(drdx)
  character                 :: file*12,title*20  
  integer                   :: i,ii,j,jj,ij,code,n(max_nod)              
  real                      :: value(max_nod)              
  contains
				     
! Definition of the green functions and derivatives
 function green(r)
  green=1./(4.*pi*r)         ; end function green
 function dgdr(r)
  dgdr=-1./(4.*pi*r**2)      ; end function dgdr
 function d2gdr(r)
  d2gdr=1./(2.*pi*r**3)	     ; end function d2gdr

! Allocating data & working  arrays  
 subroutine arrays_allocation
  allocate(node(nnod+nint),element(nelm))
  allocate(a_mat(nnod,nnod),xy_vect(nnod))
  allocate(eta1(ng),eta2(ng),w(ng),etafix(2))
 end subroutine arrays_allocation

! shape functions 
 subroutine ShapeFunctions_Quad4(eta1,eta2)
  fi(1)=.25*(1-eta1)*(1-eta2)
  fi(2)=.25*(1+eta1)*(1-eta2)
  fi(3)=.25*(1+eta1)*(1+eta2)
  fi(4)=.25*(1-eta1)*(1+eta2)
 end subroutine ShapeFunctions_Quad4

 subroutine ShapeFunctions_Triang3(eta1,eta2)
  fi(1)=eta1
  fi(2)=eta2
  fi(3)=1-eta1-eta2
 end subroutine ShapeFunctions_Triang3

end module arrays
