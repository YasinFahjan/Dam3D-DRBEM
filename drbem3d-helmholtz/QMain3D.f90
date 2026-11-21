! 3 Dimensional Dual Reciprocity Boundary Element Program
! to Solve Helmholtze Equation
! \ Constructing all the matrices then Assembling
! Using Linear 4 node and Quadratic 8 node elements
! Prepared by: by Yasin Fahjan
! Date:  15.01.1998
! Ms Fortran 90 and IMSL library
! 
program main
! Quadratic DRM3D  
 use arrays
!	
! data input
!
 print*,'DRM-3D: input data'
 call input
!
! calculation of [g],[h] and [s] 
!
 print*,'DRM-3D: calculation of [g] and [h]'
 call ghmat
 print*,'DRM-3D: calculation of [S]'
 c1=1.0; c2=1.0; c3=0.0; c4=0.0
 w=pi
 call rhsmat
!
! arranging the system 	[a]{xy}={xy}
!
 print*,'DRM-3D: arranging the system [a]{xy}={xy}'
 call Assemble
!       write(prn,'(/a/)')' A-mat :::: '
!       write(prn,'(4x,<neq>i20/)')(i,i=1,neq)
!       write(prn,'(i8,<neq>f20.7)')(i,(A_mat(i,j),j=1,neq),i=1,neq)
!       write(prn,'(/a/)')' XY-vect :::: '
!       write(prn,'(f20.7)')(xy_vect(i),i=1,neq)
!
! solving [a]{xy}={xy} using IMSL library
!
 print*,'DRM-3D: solving [a]{xy}={xy} using IMSL library	'
 call lslrg(neq,a_mat,neq,xy_vect,1,xy_vect)
!       write(prn,'(/a/)')' XY-vect :::: '
!       write(prn,'(f20.7)')(xy_vect(i),i=1,neq)
!
! distributing {xy} to {p} or {q}
!
 print*,'DRM-3D: distributing {xy} to {p} or {q}'
 call Re_Assemble
!
! calculating internal points derivatives	{qx},{qy},{qz}
!
 print*,'DRM-3D: calculating internal points derivatives {qx},{qy},{qz}'
 call interm
!
! writing the results
!
 print*,'DRM-3D: Printing the results'
 call output
!
end program main

