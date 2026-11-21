module arrays
! Quadratic DRM3D  
 use msimsl
 use portlib
 use shape_functions
  integer,        parameter :: max_nod=8       ! max # of nodes per element
  integer,        parameter :: n_dim=3	     ! Problem Dimensionality
  integer,        parameter :: ng=8            ! # of g-quadrature points 
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
  real,                save :: w
  integer,             save :: nnod,nelm,nint  ! # nodes, elements,int.pts
  integer,             save :: neq             ! # euations
  real, allocatable,   save :: eta(:)          ! Numerical Integration points
  real, allocatable,   save :: ww(:)           ! Numerical Integration weights
  real, allocatable,   save :: a_mat(:,:)	     ! system matrix
  real, allocatable,   save :: xy_vect(:)      ! system right vector
  real, allocatable,   save :: g_mat(:,:,:)    ! system matrix
  real, allocatable,   save :: h_mat(:,:)	     ! system matrix
  real, allocatable,   save :: s_mat(:,:)	     ! system matrix
  real, allocatable,   save :: f_inv(:,:)	     ! system matrix
  real, allocatable         :: u_hat(:,:)      ! g matrix element
  real, allocatable         :: q_hat(:,:,:)    ! h matrix element
  real, allocatable         :: fx(:,:,:)       ! h matrix element
  real                      :: g(max_nod)      ! g matrix element
  real                      :: h(max_nod)      ! h matrix element
  real                      :: dgdz(max_nod,n_dim) ! 1st deriv. of g w.r.t zeta
  real                      :: dhdz(max_nod,n_dim) ! 1st deriv. of h w.r.t zeta
  real 		          :: xx(n_dim)	     ! Target point coord.
  real 		          :: drdx(n_dim)     ! 1st deriv. of r w.r.t x
  real 		          :: drdnx(n_dim)    ! drdn*(drdx)
  real, 		     save :: c1,c2,c3,c4
  character                 :: file*12,title*72  
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
!
end module arrays
