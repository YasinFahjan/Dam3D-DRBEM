program Gmesh
!
! Mesh generation Program for Arch Dam-Reservoir 
! Prepared by: by Yasin Fahjan   Date:  04.03.1998
! Ms Fortran 90
! 
 use mesh
!	
! data input
!
 call input
!
! Geometrical Transformation matrix
!
 pi=4.*atan(1.)
 theta=0.*pi/180
 call transformation_matrix(theta)
!
! Data Generation for Dam structure and Reservoir
!
 call dam
 call fluid
 if(code(3)==4)call rad
!
! writing the results
!
 call output
!
end program Gmesh
