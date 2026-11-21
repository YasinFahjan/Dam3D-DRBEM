!
! 3 Dimensional FEM - BEM Fluid - Structure Interaction
!
! FEM 3D Thick-Shell and 8 solid Elements For 3D Dams Including
! Dam-Water Interaction Using Modal Analysis 
!
! DRBEM Linear 4 node and Quadratic 8 node elements to solve
! the hydrodynamic Pressure on the upstream face of the dam
!
!  The radiation condition for the reservoir
!  -DRM2D  2 linear node and Quadaratic 3 node elements
!  -FEM2D  4 linear node and Quadaratic 8 node elements 
!  -Analytical expression 2 linear node 
!  combined with continuum solution
!
!  Response spectra and time history analysis
!  utilizing FFT
!
! Using Fortran 90 and IMSL Library
! Developed By: Yasin Fahjan Feb 01, 1998
! Under Supervision : Assoc. Prof. Osman Borekci
!                     Prof. Mustafa Erdik
!
program Arch3D
use Dam3D
implicit none
!	
! Dam3D: data input for the Dam Structure	and reservoir
!
 call input
!
! Dam3D: Modal Analysis  
!
 select case(run_type)
 case(2); call DModalR
 case default; call DModal
 end select
!
! Dam3D: Time History Analysis 
!
 ! call DHistory
 !	
! Dam3D: data Output for the Dam Structure and reservoir
!
 if(ndyn>0)call Output
!
!
end program Arch3D
