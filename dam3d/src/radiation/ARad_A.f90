module Rad3D
implicit none
!
! BEM3D: Variables used in Radiation DRM2D Model
!
 integer,          parameter :: max_r=3           ! DRM2D/# nodes per element
!
! DRM\Rad: Definition the derived type: elementr 
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
! DRM\Rad: Definitons of Variables used in Radiation DRM2D Model
!
!
! Analytical\Rad:  Variables used in BEM3D Model
!
! integer,          parameter :: nn=9              ! BEM3D/# nodes per element
! integer,          parameter :: mm=9              ! BEM3D/# nodes per element
 real, 	       parameter :: bb=100. !318.19 	        ! DRM3D/Approx.funct.Coef.
 real, 	       parameter :: hh=100. 	        ! DRM3D/Approx.funct.Coef.
 complex,allocatable,   save :: evalz(:)		  ! Eigen Values 
 complex,allocatable,   save :: evaly(:)		  ! Eigen Values 
!
!
 integer,allocatable,   save :: lmr(:)		  ! Radiation {LM}
 complex,allocatable,   save :: eval(:)		  ! Eigen Values 
 complex(8),allocatable,save :: evect(:,:)	  ! Eigen Vectors
 complex,allocatable,   save :: Knm(:,:)		  ! Knm matrix
 complex(8),allocatable,save :: Rad(:,:)		  ! Radiation Matrix
 complex,allocatable,   save :: G4(:,:)		  ! G matrix partition
!
 real, allocatable,     save :: grad(:,:,:)       ! RAD3D/G matrix
 real, allocatable,     save :: hrad(:,:)	        ! RAD3D/H matrix
 real, allocatable,     save :: srad(:,:)	        ! RAD3D/S matrix
!
 complex,allocatable,   save :: g_rad(:,:)	  ! RAD3D/EVP G matrix
 complex,allocatable,   save :: G11(:,:)		  ! G matrix partition
 complex,allocatable,   save :: h_rad(:,:)	  ! RAD3D/EVP H matrix
 complex,allocatable,   save :: s_rad(:,:)        ! RAD3D/EVP S matrix
 complex,allocatable,   save :: d_rad(:)          ! RAD3D/EVP d vector
!
 integer                     :: ndr, neqr             
 integer                     :: nr1,nr2,nk1,nk2             
!
end module Rad3D 

