! 3 Dimensional Dual Reciprocity Boundary Element Program
! to Solve Helmholtze Equation 
! Using Linear 4 node and Quadratic 8 node elements
! Prepared by: by Yasin Fahjan
! Date:  15.12.1997
! Ms Fortran 90 and IMSL library
! 
program main
! Quadratic BEM3D  
 use arrays
!	
! data input
!
 call input
!
! calculation of [g] and [h] 
!
 print*,'calculation of [g] and [h]'
 call ghmat
!
! solving [a]{xy}={xy} using IMSL library
!
 print*,'solving [a]{xy}={xy} using IMSL library	'
 call lslcg(neq,a_mat,neq,xy_vect,1,xy_vect)
!
! distributing {xy} to {p} or {q}
!
       write(prn,'(//a/)')' xy-vect'
       write(prn,'(i8,f20.7)')(i,abs(xy_vect(i)),i=1,neq)

 print*,'distributing {xy} to {p} or {q}'
 where(node(1:neq)%code==0)node(1:neq)%p=xy_vect
 where(node(1:neq)%code==2)node(1:neq)%p=xy_vect
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
 print*,'calculating internal points derivatives	{qx},{qy},{qz}'
 call interm
!
! writing the results
!
 call output
 close (inp)
 close (prn)
end program main

