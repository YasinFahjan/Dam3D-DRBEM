module arrays
 integer,         parameter:: n_dim=3, nod=4      
 integer,         parameter:: inp=11, ipr=12      
! Definition the derived type structure node & element 
 type node_record
 real            :: x(n_dim)
 integer         :: nnum
 end type node_record
 type CornerPoints_record
 real            :: x(n_dim)
 end type  CornerPoints_record   
 type element_record
 integer         :: node(nod),kode
 real            :: value
 real            :: x(n_dim)
 end type element_record
! Definition of the variables in the program 
 type(node_record),allocatable,   save:: node(:,:) 
 type(element_record),allocatable,save:: element(:) 
 type(CornerPoints_record),       save:: corner(8) 
 character,            save:: file*12,title*80,indicator*1
 integer,              save:: nx,ny,nz,n1(6),n2(6),nnx(6),nnz(6),ni(6),nj(6)
 integer,              save:: nn1,nn2,nnod,nelm,kode(6)
 real,                 save:: value(6)
 real,                 save:: T(3,3)
 real,                 save:: x1(3),x2(3),x3(3),x4(3),xx1(3),xx2(3),x(3)
 contains
				     
! Allocating data & working  arrays  
 subroutine array_allocations
 allocate(node(nn1,nn2),element(nelm))
 end subroutine array_allocations
! Transformation matrix 
 subroutine transformation_matrix(theta)
  T(1,1)=cos(theta) ; T(1,2)=-sin(theta) ; T(1,3)=0. 
  T(2,1)=sin(theta) ; T(2,2)= cos(theta) ; T(2,3)=0. 
  T(3,1)=0.         ; T(3,2)=0.          ; T(3,3)=1. 
 end subroutine transformation_matrix

! Face-Shape Functions 
 subroutine Face_shape(face,x)
 integer      :: face
 real         :: x(3) 
 select case(face)
 case(1)

 case(2)

 case(3)
 
 case(4)
 
 case(5)
 
 case(6)
 
 case default
 end select
 x=x
 end subroutine Face_shape

end module arrays
	
