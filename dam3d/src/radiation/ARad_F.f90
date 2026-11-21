module Rad3D
implicit none
!
! FEM2D: Variables used in Radiation FEM2D Model
!
 integer,          parameter :: max_r=3           ! FEM2D/# nodes per element
!
! FEM\Rad: Definition the derived type: elementr 
!
type elementr_record
 integer                     :: node(max_r)       ! DRM2D/element nodes 
 integer                     :: code		  ! DRM2D/code(1 Pot.,2 Flux)
 complex                     :: q(max_r)	        ! DRM2D/nodes:fluxes 
 complex                     :: p(max_r)	        ! DRM2D/nodes:potential 
 real                        :: norm(max_r,3)	  ! DRM2D/nodes:normal 
end type elementr_record
 type(elementr_record),allocatable,  save:: elementr(:) 

!
! FEM\Rad: Definitons of Variables used in Radiation FEM2D Model
!
 integer,allocatable,   save :: lmr(:)		  ! Radiation {LM}
 complex(8),allocatable,save :: eval(:)		  ! Eigen Values 
 complex(8),allocatable,save :: evect(:,:)	  ! Eigen Vectors
 complex,allocatable,   save :: Knm(:,:)		  ! Knm matrix
 complex(8),allocatable,save :: Rad(:,:)		  ! Radiation Matrix
 complex,allocatable,   save :: G4(:,:)		  ! G matrix partition
!
 real, allocatable,     save :: grad(:,:)         ! RAD3D/G matrix
 real, allocatable,     save :: hrad(:,:)	        ! RAD3D/H matrix
 real, allocatable,     save :: brad(:,:)	        ! RAD3D/S matrix
 real, allocatable,     save :: drad(:)	        ! RAD3D/S matrix
!
 real, allocatable,     save :: g_rad(:,:)	  ! RAD3D/EVP G matrix
 real, allocatable,     save :: h_rad(:,:)	  ! RAD3D/EVP H matrix
 real, allocatable,     save :: b_rad(:,:)        ! RAD3D/EVP S matrix
 real, allocatable,     save :: d_rad(:)          ! RAD3D/EVP d vector
!
 integer                     :: ndr, neqr             
 integer                     :: nr1,nr2,nk1,nk2             
!
 complex ,allocatable :: HB_mat(:,:),GH_mat(:,:)
!
end module Rad3D 

