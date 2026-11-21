module Dam3D
use msimsl
use portlib
implicit none
!
! Dam3D: General Variables used in FEM and BEM3D/ models
!
 integer,          parameter :: inp=11          ! DAM3D/Input File #     
 integer,          parameter :: prn=12          ! DAM3D/Output File #     
!
! Dam3D: Definition of the derived type: node 
!
type node_record
 integer                     :: id(3)	      ! FEM3D/FEM ID 
 integer                     :: idw  	      ! BEM3D/BEM ID
 integer                     :: idr  	      ! BEM3D/RAD ID 
 real                        :: coord(3)	      ! DAM3D/node Geom. Coord.
 real                        :: disp(3)         ! FEM3D/node displacements
 real                        :: accel(3)        ! FEM3D/node accelerations
 complex                     :: p      	      ! BEM3D/node potential
 complex                     :: q(3)	      ! BEM3D/node fluxes 
 integer                     :: code 	      ! BEM3D/node p or q is known
 integer                     :: codr 	      ! RAD2D/node p or q is known
end type node_record
type(node_record),allocatable,   save:: node(:) 
!
! Dam3D: Definition of the General Variables
!
 integer,               save :: numel            ! FEM3D/#dam elements                 
 integer,               save :: nelm             ! BEM3D/#reservoir elements
 integer,               save :: nelr             ! DRM2D/#Radiation elements                  
 integer,               save :: neq              ! FEM3D/#structure DOF    
 integer,               save :: nq               ! BEM3D/#Reservoir NEQ  
 integer,               save :: nqr              ! RAD2D/#Radiation NEQ   
 integer,               save :: ntot             ! DAM3D/#total nodes   
 integer,               save :: numnp            ! FEM3D/#Dam nodes   
 integer,               save :: nnod             ! DAM3D/#rteservoir nodes   
 integer,               save :: nint             ! BEM3D/#DRM internal nodes  
 real,                  save :: pi
 real,                  save :: gravity          ! DAM3D/gravity 
 real,                  save :: w	             ! DAM3D/frequency
 real,                  save :: wo	             ! DAM3D/frequency
 real,                  save :: wi,wf,wr         ! DAM3D/intial,final,step relative frequency 
 real,                  save :: rho_w            ! BEM3D/water density  
 real,                  save :: alfa_r           ! BEM3D/reflection Coeffecient
 real,                  save :: hwater           ! BEM3D/reflection Coeffecient
 real,                  save :: gama             ! BEM3D/foundation Damping
 real,                  save :: sound            ! BEM3D/sound speed in the water
 real,                  save :: hyst_damp	       ! FEM3D/hystertic damping
 integer,               save :: ixyz             ! DAM3D/ground motion direction
 real,                  save :: ag(3)            ! DAM3D/ground acceleration
 character(72),         save :: heading	       ! DAM3D/heading
 character(72),         save :: Void	       ! DAM3D/heading
 character(72),         save :: file  	       ! DAM3D/file name 
 integer,               save :: nf               ! DAM3D/# req. mode shapes              
 integer,               save :: nf_i             ! DAM3D/# req. mode shapes              
 integer,               save :: nf_f             ! DAM3D/# req. mode shapes              
 integer,               save :: nf_s             ! DAM3D/# req. mode shapes              
 integer,               save :: ndyn             ! DAM3D/0:EVP, 1:D, 2:D+H, 3:H             
 integer,               save :: run_type         ! DAM3D/0:EVP, 1:D, 2:D+H, 3:H             
!
 integer                     :: i,ii,j,jj,ij,k,kk,i_xyz             
 integer		           :: nel,mat
 integer		           :: counter
!
end module Dam3D 