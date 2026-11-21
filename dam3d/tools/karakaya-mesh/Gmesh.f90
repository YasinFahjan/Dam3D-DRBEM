module mesh
!
! Definition of the variables used in the Mesh generation
!
 integer,         parameter:: inp=11     	 
 integer,         parameter:: prn=12       
 integer,         parameter:: prn2=13      
 integer,         parameter:: prn3=14      
 integer,         parameter:: prn4=15      
 integer,         parameter:: prn5=16      
!
! Definition the derived type structure node & element 
!
 type node_type
 real            :: radius
 real            :: angle
 real            :: height
 real            :: coord(3)
 integer         :: nnum
 integer         :: code
 integer         :: id(3)
 end type node_type
!
 type nodes_type
 real            :: radius
 real            :: angle
 real            :: height
 real            :: coord(3)
 integer         :: nnum
 integer         :: id(3)
 character(1)    :: indicator
 end type nodes_type
!
 type CornerPoints_type
 real            :: coord(3)
 real            :: radian(3)
 end type  CornerPoints_type   
!
 type elements_type
 integer         :: node(16)
 integer         :: mat
 integer         :: ng
 real            :: coord(3)
 end type elements_type
!
 type element_type
 integer         :: node(8)
 integer         :: code
 real            :: coord(3)
 real            :: value
 end type element_type
!
 type elementr_type
 integer         :: node(3)
 integer         :: code
 end type elementr_type

!
! Definition of the variables in the program
! 
 type(node_type),allocatable,    save:: node(:,:) 
 type(nodes_type),allocatable,   save:: nodes(:,:,:) 
 type(element_type),allocatable, save:: element(:) 
 type(elements_type),allocatable,save:: elements(:) 
 type(elementr_type),allocatable,save:: elementr(:) 
 type(CornerPoints_type),        save:: corner(8) 
 type(CornerPoints_type),        save:: corners(32) 
 character(12),                  save:: file*12
 character(80),                  save:: title*80
 character(1),                   save:: indicator*1
 integer,                        save:: ntot,numnp,nnod,nint
 integer,                        save:: numel,mat
 integer,                        save:: nelm
 integer,                        save:: nelr
 integer,                        save:: nx,ny,nz
 integer,                        save:: inx,iny,inz
 integer,                        save:: nn1,nn2
 integer,                        save:: n1(6),n2(6)
 integer,                        save:: nnx(6),nnz(6)
 integer,                        save:: ni(6),nj(6)
 integer,                        save:: code(6)
 real,                           save:: pi
 real,                           save:: radius
 real,                           save:: angle 
 real,                           save:: height
 real,                           save:: ri,rf
 real,                           save:: value(6)
 real,                           save:: T(3,3)
 real,                           save:: x10(3),x20(3),x30(3),x40(3)
 real,                           save:: x11(3),x21(3),x31(3),x41(3)
 real,                           save:: x1(3),x2(3),x3(3),x4(3)
 real,                           save:: xx(3),xx1(3),xx2(3),x(3)
 real,                           save:: xi_int(3),xf_int(3),dx_int(3)
 integer,                        save:: nnum
 integer,                        save:: llmm(16)
 character(1),                   save:: rff(3)

 character(8),                   save:: kind
 real,                           save:: ag(3),ee,enu,rho_d,gravity,rho_w
 real,                           save:: sound,alfa_r,hwater
 integer,                        save:: nod,ndr,nrot
 integer,                        save:: sap,porter,lusas,sap_res,Gmode
 integer,                        save:: nf,ndyn,ixyz,nmat,ng,radiation,geometry
 integer,                        save:: start_i(16),finish_i(16)
 integer,                        save:: start_j(16),finish_j(16)
 integer,                        save:: inc_i(16),inc_j(16)
 integer,                        save:: nl,nc,nr,lb,lc,lu,ldam,ndam,mdam,mb
 integer,                        save:: pts(8,4,2),pts1,pts2,pts3,pts4
 integer,                        save:: test

contains
!
! Transformation matrix 
!
 subroutine transformation_matrix(theta)
 T(1,1)=cos(theta) ; T(1,2)=-sin(theta) ; T(1,3)=0. 
 T(2,1)=sin(theta) ; T(2,2)= cos(theta) ; T(2,3)=0. 
 T(3,1)=0.         ; T(3,2)=0.          ; T(3,3)=1. 
 end subroutine transformation_matrix

end module mesh
