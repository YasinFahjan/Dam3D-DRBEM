module Bem3D
use Dam3D
use shape_functions
implicit none
!
! BEM3D: Variables used in BEM3D model
!
 integer,          parameter :: max_n=8           ! BEM3D/# nodes per element
 integer,          parameter :: ng=4              ! BEM3D/quadrature points 
!
! BEM3D: Definition the derived type: element 
!
type element_record
 integer                     :: node(max_n)       ! BEM3D/element nodes 
 integer                     :: code		  ! BEM3D/code(1 Pot.,2 Flux)
 complex                     :: q(max_n)	        ! BEM3D/nodes:fluxes 
 complex                     :: p(max_n)	        ! BEM3D/nodes:potential 
 real                        :: norm(max_n,3)	  ! BEM3D/nodes:normal 
end type element_record
 type(element_record),allocatable,  save:: element(:) 
!
! BEM3D: Definition of the BEM variables  
!
 real, allocatable,     save :: eta(:)            ! BEM3D/quadrature points
 real, allocatable,     save :: ww(:)             ! BEM3D/quadrature weights

 integer,allocatable,   save :: lmw(:)            ! BEM3D/{LM}
 complex(8),allocatable,save :: a_mat(:,:)	  ! BEM3D/system matrix
 complex(8),allocatable,save :: xy_vect(:)        ! BEM3D/system right vector
 real, allocatable,     save :: g_mat(:,:,:)      ! BEM3D/G matrix
 real, allocatable,     save :: h_mat(:,:)	  ! BEM3D/H matrix
 real, allocatable,     save :: s_mat(:,:)	  ! DRM3D/S matrix
 real(8),allocatable,   save :: f_inv(:,:)        ! DRM3D/F inverse
 real, allocatable           :: u_hat(:,:)        ! DRM3D/Uhat matrix
 real, allocatable           :: q_hat(:,:,:)      ! DRM3D/Qhat matrix
 real, allocatable           :: fx(:,:,:)         ! DRM3D/(df/dx) matrix 
 real                        :: g(max_n)          ! BEM3D/g matrix element
 real                        :: h(max_n)          ! BEM3D/h matrix element
 real                        :: dgdz(max_n,3)     ! BEM3D/(dg/d.zeta)
 real                        :: dhdz(max_n,3)     ! BEM3D/(dh/d.zeta)
 real 		           :: xx(3)	        ! BEM3D/Target point coord.
 real 		           :: drdx(3)           ! BEM3D/(dr/dx)
 real 		           :: drdn              ! BEM3D/drdn*(drdx)
 real              	     :: c(0:20)           ! DRM3D/Approx.funct.Coef.
 integer              	     :: cterm             ! DRM3D/Approx.funct.Coef.
!
 integer                     :: nod              
 integer                     :: np(4),np1,np2,np3,np4              
 integer                     :: nf1,nf2,nf3,nf4             
 integer                     :: code,n(max_n)              
 real                        :: value(max_n), r
!
end module Bem3D
