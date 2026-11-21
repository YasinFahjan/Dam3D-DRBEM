! Mesh generation Program for 3D- BEM 
! Laplace Equation  with Drichlet & Neuman Conditions
! Prepared by: by Yasin Fahjan   Date:  18.08.1997
! Ms Fortran 90
! 
program main
 use arrays
!	
! data input
!
 call input
!
! data Generation 
!
 call generation
!
! writing the results
!
 call output
!
 close (inp)
 close (ipr)
end program main
