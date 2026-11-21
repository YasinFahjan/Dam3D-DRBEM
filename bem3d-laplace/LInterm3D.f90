subroutine interm_quadrilateral
! Linear BEM3D  
 use arrays
 use numeric
!
! solution at internal points once boundary solution solution is known
!
 write(*,'(/a)')" Solution at internal points"
 do i=nnod+1,nnod+nint
  node(i)%p=0.; node(i)%q=0.
  zeta=node(i)%x
  do j=1,nelm
   n(1:nod)=element(j)%node(1:nod) 
   ds1=1/2.*(-node(n(1))%x-node(n(2))%x+node(n(3))%x+node(n(4))%x)
   ds2=1/2.*( node(n(1))%x-node(n(2))%x-node(n(3))%x+node(n(4))%x)
   call cross_product(ds1,ds2,norm)
   ds1_ds2=sqrt(sum(norm*norm)); norm=norm/ds1_ds2
!    
! Evaluating the g,h,dg/d(zeta),dh/d(zeta)  
!   		
   g=0.; h=0.; dgdz=0.; dhdz=0.;
   do ii=1,ng
   do jj=1,ng
    call ShapeFunctions_Quad4(eta1(ii),eta1(jj))
    xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
    r=sqrt(sum((zeta-xx)**2))
    drdx=(xx-zeta)/r; drdn=sum(drdx*norm); drdnx=-norm+drdn*drdx
    g(1:nod)=g(1:nod)+fi(1:nod)*w(ii)*w(jj)*green(r)*ds1_ds2
    h(1:nod)=h(1:nod)+fi(1:nod)*w(ii)*w(jj)*dgdr(r)*drdn*ds1_ds2
    do ij=1,nod
     dgdz(ij,1:n_dim)=dgdz(ij,1:n_dim)+fi(ij)*w(ii)*w(jj)*&
        &dgdr(r)*(-drdx(1:n_dim))*ds1_ds2
     dhdz(ij,1:n_dim)=dhdz(ij,1:n_dim)+fi(ij)*w(ii)*w(jj)*&
        &(d2gdr(r)*drdn*(-drdx(1:n_dim))+dgdr(r)/r*drdnx(1:n_dim))*ds1_ds2
    enddo
   enddo; enddo;
!
! Evaluating p, q  
!		
   do ij=1,nod
    node(i)%p=node(i)%p-element(j)%p(ij)*h(ij)+element(j)%q(ij)*g(ij)
    node(i)%q(1:n_dim)=node(i)%q(1:n_dim)-element(j)%p(ij)*dhdz(ij,1:n_dim)&
       &+element(j)%q(ij)*dgdz(ij,1:n_dim)
   enddo
  enddo;
  node(i)%p=node(i)%p/4.0
  node(i)%q(:)=node(i)%q(:)/4.0
 enddo

end subroutine interm_quadrilateral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interm_triangular
! Linear BEM3D  
 use arrays
 use numeric
!
! solution at intenal points once boundary solution solution is known
!
 do i=nnod+1,nnod+nint
  node(i)%p=0.; node(i)%q=0.
  zeta=node(i)%x
  do j=1,nelm
   n=element(j)%node 
   ds1=node(n(2))%x-node(n(1))%x
   ds2=node(n(3))%x-node(n(2))%x
   call cross_product(ds1,ds2,norm)
   ds1_ds2=sqrt(sum(norm*norm)); norm=norm/ds1_ds2
!
! Evaluating the g,h,dg/d(zeta),dh/d(zeta)  
!		
   g=0.; h=0.; dgdz=0.; dhdz=0.;
   do ii=1,ng
    call ShapeFunctions_Triang3(eta1(ii),eta2(ii))
    xx=0.; do ij=1,nod; xx=xx+fi(ij)*node(n(ij))%x; enddo
    r=sqrt(sum((zeta-xx)**2))
    drdx=(xx-zeta)/r; drdn=sum(drdx*norm); drdnx=-norm+drdn*drdx
    g(1:nod)=g(1:nod)+fi(1:nod)*w(ii)*green(r)*ds1_ds2
    h(1:nod)=h(1:nod)+fi(1:nod)*w(ii)*dgdr(r)*drdn*ds1_ds2
    do ij=1,nod
     dgdz(ij,1:n_dim)=dgdz(ij,1:n_dim)+fi(ij)*w(ii)*&
       &dgdr(r)*(-drdx(1:n_dim))*ds1_ds2
     dhdz(ij,1:n_dim)=dhdz(ij,1:n_dim)+fi(ij)*w(ii)*&
       &(d2gdr(r)*drdn*(-drdx(1:n_dim))+dgdr(r)/r*drdnx(1:n_dim))*ds1_ds2
    enddo
   enddo
   do ij=1,nod
    node(i)%p=node(i)%p-element(j)%p(ij)*h(ij)+element(j)%q(ij)*g(ij)
    node(i)%q(1:n_dim)=node(i)%q(1:n_dim)-element(j)%p(ij)*dhdz(ij,1:n_dim)&
       &+element(j)%q(ij)*dgdz(ij,1:n_dim)
   enddo
  enddo;  
  node(i)%p=node(i)%p/2.0
  node(i)%q(:)=node(i)%q(:)/2.0
 enddo

end subroutine interm_triangular
