module Fem3D
use Dam3D
use shape_functions
implicit none
!
! FEM3D: Variables used in FEM3D model
!
 integer,          parameter :: prn1=22           ! FEM3D/Output File #     
 integer,          parameter :: max_n=16          ! BEM3D/# nodes per element
!
! FEM3D: Definition the derived type: element 
! 
type elements_record
 integer                     :: node(16)          ! FEM3D/element nodes
 integer                     :: mat	              ! FEM3D/material type
 integer                     :: ng     	        ! FEM3D/# integ.points
 end type elements_record
type(elements_record),allocatable,  save:: elements(:) 
!
! FEM3D: Definition of the variables in the program 
!
 integer,               save :: nmat              ! FEM3D/# dam materials                 
 integer,allocatable,   save :: lm(:)       	  ! FEM {LM}
 real,allocatable            :: xk1(:),xk2(:)     ! G-quad. abss.
 real,allocatable            :: xk3(:)            ! G-quad. abss.
 real,allocatable            :: wgt1(:),wgt2(:)   ! G-quad. weights
 real,allocatable            :: wgt3(:)           ! G-quad. weights
 real,allocatable,      save :: ee(:)		  ! Dam Material: E
 real,allocatable,      save :: enu(:)	        ! Dam Material: nu
 real,allocatable,      save :: rho_d(:)	        ! Dam Material: unit weight
 real,allocatable,      save :: mass(:,:) 	  ! element mass matrix
 real,allocatable,      save :: stiff(:,:)        ! element stiffness matrix
 real,allocatable,      save :: stiff_alfa(:,:)   ! element rotational stiff
 real,allocatable,      save :: mass_mat(:,:)     ! Global mass matrix 
 real,allocatable,      save :: stiff_mat(:,:)    ! Global stifness matrix 
 integer,allocatable,   save :: nn(:)             ! Element: nodes numbers
 integer,allocatable,   save :: id(:,:)	        ! Element: nodes ID
 real,allocatable,      save :: G(:,:)	        ! Element:Strain Matrix [G]
 real,allocatable,      save :: D(:,:)	        ! Element:stress-strain matrix [D]
!
 integer                     :: ng
 integer                     :: nods,rots,ndr
 integer                     :: lx,ly,lz
 real                        :: wg,ee1,ee2,ee3
 real                        :: Stiff_coeff,vol
!
end module Fem3D
