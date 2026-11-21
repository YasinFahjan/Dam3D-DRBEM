subroutine ghmat_quadrilateral
! Linear BEM3D  
 use arrays
 use numeric
!
! classical Gauss-Legendre quadrature formula using imsl
!
call gqrul (ng,1,0.,0.,0,etafix,eta1,w) 
!
! Assembling [h] and [g] elements 
!
  a_mat=0.; xy_vect=0.
 write(*,'(/a)')" Assembling [h] and [g] elements"
 do i=1,nnod
  zeta=node(i)%x
  h_diag=0.
  write(*,'(i5,1x\)')i
  do j=1,nelm
   n(1:nod)=element(j)%node(1:nod) 
   ds1=1/2.*(-node(n(1))%x-node(n(2))%x+node(n(3))%x+node(n(4))%x)
   ds2=1/2.*( node(n(1))%x-node(n(2))%x-node(n(3))%x+node(n(4))%x)
   call cross_product(ds1,ds2,norm)
   ds1_ds2=sqrt(sum(norm*norm)); norm=norm/ds1_ds2
!
! Evaluating the h,g matrices terms using Guass quadrature 
!		
   g=0.; h=0.; dgdz=0.; dhdz=0.;
   do ii=1,ng
   do jj=1,ng
    call ShapeFunctions_Quad4(eta1(ii),eta1(jj))
    xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
    r=sqrt(sum((xx-zeta)**2))
    drdx=(xx-zeta)/r; drdn=sum(drdx*norm) 
    g(1:nod)=g(1:nod)+fi(1:nod)*w(ii)*w(jj)*green(r)*ds1_ds2
    h(1:nod)=h(1:nod)+fi(1:nod)*w(ii)*w(jj)*dgdr(r)*drdn*ds1_ds2
   enddo; enddo;
! Diagonal elements of h matrix
   h_diag=h_diag-sum(h(1:nod)) 
!
! arranging the system 	[a]{xy}={xy}
!  
   do ij=1,nod
    select case (node(n(ij))%code)
     case (1)
      xy_vect(i)=xy_vect(i)-h(ij)*node(n(ij))%p
    	if(element(j)%code==1)then
       a_mat(i,n(ij))=a_mat(i,n(ij))-g(ij)
	else
	 xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
	endif
     case (2)
      xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
      a_mat(i,n(ij))=a_mat(i,n(ij))+h(ij)
    end select
   enddo
  enddo;
!
! Modification for Diagonal Term
!
  select case (node(i)%code)
   case (1); xy_vect(i)=xy_vect(i)-h_diag*node(i)%p
   case (2); a_mat(i,i)=a_mat(i,i)+h_diag
  end select

 enddo

end subroutine ghmat_quadrilateral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ghmat_triangular
! Linear BEM3D  
 use arrays
 use numeric
!
! Gauss-quadrature Weights for Triangular Domains
!
 call GaussQuadrature_Triangle(eta1,eta2,w,ng)
!
! Assembling [h] and [g] elements 
!
  a_mat=0.; xy_vect=0.
 write(*,'(/a)')" Assembling [h] and [g] elements"
 do i=1,nnod
  zeta=node(i)%x
  h_diag=0.
  write(*,'(i5,1x\)')i
  do j=1,nelm
   n(1:nod)=element(j)%node(1:nod) 
   ds1=node(n(2))%x-node(n(1))%x
   ds2=node(n(3))%x-node(n(2))%x
   call cross_product(ds1,ds2,norm)
   ds1_ds2=sqrt(sum(norm*norm)); norm=norm/ds1_ds2
!
! Evaluating the h,g matrices terms using Guass quadrature 
!
   g=0.; h=0.; dgdz=0.; dhdz=0.;
   do ii=1,ng
    call ShapeFunctions_Triang3(eta1(ii),eta2(ii))
    xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
    r=sqrt(sum((xx-zeta)**2))
    drdx=(xx-zeta)/r; drdn=sum(drdx*norm) 
    g(1:nod)=g(1:nod)+fi(1:nod)*w(ii)*green(r)*ds1_ds2
    h(1:nod)=h(1:nod)+fi(1:nod)*w(ii)*dgdr(r)*drdn*ds1_ds2
   enddo
! Diagonal elements of h matrix
   h_diag=h_diag-sum(h(1:nod)) 
!
! arranging the system 	[a]{xy}={xy}
!  
   do ij=1,nod
    select case (node(n(ij))%code)
     case (1)
      xy_vect(i)=xy_vect(i)-h(ij)*node(n(ij))%p
    	if(element(j)%code==1)then
       a_mat(i,n(ij))=a_mat(i,n(ij))-g(ij)
	else
	 xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
	endif
     case (2)
      xy_vect(i)=xy_vect(i)+g(ij)*element(j)%q(ij)
      a_mat(i,n(ij))=a_mat(i,n(ij))+h(ij)
    end select
   enddo

  enddo
!
! Modification for Diagonal Term
!
  select case (node(i)%code)
   case (1); xy_vect(i)=xy_vect(i)-h_diag*node(i)%p
   case (2); a_mat(i,i)=a_mat(i,i)+h_diag
  end select

 enddo

end subroutine ghmat_triangular
