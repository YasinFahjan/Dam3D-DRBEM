subroutine generation
 use arrays
! Geometrical Transformation matrix
 pi=4.*atan(1.)
 theta=0.*pi/180
 call transformation_matrix(theta)
!
! face- node generation	 1-4 ,5-6
!
 nnum=0
 n1(1)=1        ; n2(1)=nz+1   ; nnx(1)=nx; nnz(1)=nz; nj(1)=nx; ni(1)=nz-1  
 n1(2)=nz+1     ; n2(2)=nz+1   ; nnx(2)=nx; nnz(2)=ny; nj(2)=nx; ni(2)=ny-1 
 n1(3)=ny+nz+1  ; n2(3)=nz+1   ; nnx(3)=nx; nnz(3)=nz; nj(3)=nx; ni(3)=nz-1 
 n1(4)=ny+2*nz+1; n2(4)=nz+1   ; nnx(4)=nx; nnz(4)=ny; nj(4)=nx; ni(4)=ny-1 
 n1(5)=nz+1     ; n2(5)=1      ; nnx(5)=nz; nnz(5)=ny; nj(5)=nz-1; ni(5)=ny-1  
 n1(6)=nz+1     ; n2(6)=nz+nx+1; nnx(6)=nz; nnz(6)=ny; nj(6)=nz-1; ni(6)=ny-1 
 
 do iface=1,6
  inc=(iface-1)*2
  select case(iface)
   case(4); inc1=inc+1; inc2=inc+2;inc3=1; inc4=2; i1=0; j1=0
   case(5); inc1=1; inc2=3; inc3=7; inc4=5; i1=1; j1=1
   case(6); inc1=4; inc2=2; inc3=6; inc4=8; i1=1; j1=1
   case default;  inc1=inc+1; inc2=inc+2; inc3=inc+3; inc4=inc+4; i1=0; j1=0
  end select
 do i=i1,ni(iface)
 do j=j1,nj(iface)
  i_index=n1(iface)+i; j_index=n2(iface)+j
  x1=corner(inc1)%x ;x2=corner(inc2)%x; x3=corner(inc3)%x; x4=corner(inc4)%x
  xx1=x1+(x2-x1)/nnx(iface)*j; xx2=x3+(x4-x3)/nnx(iface)*j
  node(i_index,j_index)%x=xx1+(xx2-xx1)/nnz(iface)*i
  x= node(i_index,j_index)%x
 call face_shape(iface,x)
!  if(iface==1.and.(x(1)>2.0.and.x(1)<4.0).and.x(3)/=0.)then
!  node(i_index,j_index)%nnum=0
!  else
  nnum=nnum+1
  node(i_index,j_index)%nnum=nnum
  node(i_index,j_index)%x=matmul(T,node(i_index,j_index)%x)
!  endif
 enddo;enddo;enddo
  nnod=nnum
!
! transfering the information to the associated nodes  
!
  node(nn1,:)=node(1,:)
  node(nz+1,1:nz)=node(1:nz,nz+1)
  node(nz+1,nn2:nn2-nz:-1)=node(1:nz,nz+nx+1)
 
  node(nz+ny+1,1:nz)=node(2*nz+ny+1:nz+ny+1:-1,nz+1)
  node(nz+ny+1,nn2:nn2-nz:-1)=node(2*nz+ny+1:nz+ny+1:-1,nz+nx+1)
 
  node(nz+1:nz+ny,1)=node(nn1:nn1-ny:-1,nz+1)
  node(nz+1:nz+ny,nn2)=node(nn1:nn1-ny:-1,nz+nx+1)
!
! Generation of 4- nodes quadilateral element information 
!
   n1(1)=1         ; nnx(1)=nz      ; n2(1)=nz+1     ; nnz(1)=nz+nx
   n1(2)=nz+1      ; nnx(2)=nz+ny   ; n2(2)=nz+1     ; nnz(2)=nz+nx
   n1(3)=nz+ny+1   ; nnx(3)=2*nz+ny ; n2(3)=nz+1     ; nnz(3)=nz+nx
   n1(4)=2*nz+ny+1 ; nnx(4)=nn1-1   ; n2(4)=nz+1     ; nnz(4)=nz+nx
   n1(5)=nz+1      ; nnx(5)=nz+ny   ; n2(5)=1        ; nnz(5)=nz 
   n1(6)=nz+1      ; nnx(6)=nz+ny   ; n2(6)=nz+nx+1  ; nnz(6)=2*nz+nx 
   nnum=0
  do ii=1,6
  do i=n1(ii),nnx(ii)
  do j=n2(ii),nnz(ii)
  select case(nod)
  case(3)
   if(node(i,j)%nnum*node(i,j+1)%nnum*node(i+1,j)%nnum*node(i+1,j+1)%nnum==0)then
   else
   nnum=nnum+1
   element(nnum)%node(1)=node(i,j)%nnum
   element(nnum)%node(2)=node(i,j+1)%nnum
   element(nnum)%node(3)=node(i+1,j)%nnum
   element(nnum)%kode=kode(ii); element(nnum)%value =value(ii)
   element(nnum)%x=(node(i,j)%x+node(i,j+1)%x+node(i+1,j)%x)/3.
   nnum=nnum+1
   element(nnum)%node(1)=node(i,j+1)%nnum
   element(nnum)%node(2)=node(i+1,j+1)%nnum
   element(nnum)%node(3)=node(i+1,j)%nnum
   element(nnum)%kode=kode(ii); element(nnum)%value =value(ii)
   element(nnum)%x=(node(i,j+1)%x+node(i+1,j+1)%x+node(i+1,j)%x)/3.
   endif
  case(4)
   if(node(i,j)%nnum*node(i,j+1)%nnum*node(i+1,j)%nnum*node(i+1,j+1)%nnum==0)then
   else
   nnum=nnum+1
   element(nnum)%node(1)=node(i,j)%nnum
   element(nnum)%node(4)=node(i,j+1)%nnum
   element(nnum)%node(3)=node(i+1,j+1)%nnum
   element(nnum)%node(2)=node(i+1,j)%nnum
   element(nnum)%kode=kode(ii); element(nnum)%value =value(ii)
   element(nnum)%x=(node(i,j)%x+node(i,j+1)%x+node(i+1,j+1)%x+node(i+1,j)%x)/4.
   endif
  case default
  end select 
  enddo; enddo; enddo
  nelm=nnum
end subroutine generation