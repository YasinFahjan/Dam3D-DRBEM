! 3 Dimensional Boundary Element Program to Solve 
! Laplace Equation  with Drichlet & Neuman Conditions
! Prepared by: by Yasin Fahjan
! Date:  11.08.1997
! Ms Fortran 90
! 
program main
! Linear BEM3D  
 use arrays
!	
! data input
!
 call input
!
! calculation of [g] and [h] 
!
 if(nod==3)call ghmat_triangular
 if(nod==4)call ghmat_quadrilateral
!
! solving [a]{xy}={xy} using IMSL library
!
 write(*,'(/a)')" Solving [a]{xy}={xy}"
 call lslrg(nnod,a_mat,nnod,xy_vect,1,xy_vect)
!
! distributing {xy} to {p} or {q}
!
 where(node(1:nnod)%code==2)node(1:nnod)%p=xy_vect
 do i=1,nelm
  n(1:nod)=element(i)%node(1:nod)
  do ij=1,nod
   select case (element(i)%code)
    case(1); element(i)%q(ij)=xy_vect(n(ij))
    case(2); element(i)%p(ij)=node(n(ij))%p
   end select
  enddo
 enddo
!
! calculating internal points derivatives	{qx},{qy},{qz}
!
 if(nod==3)call interm_triangular
 if(nod==4)call interm_quadrilateral
!
! writing the results
!
 call output
 close (inp)
 close (prn)
end program main


