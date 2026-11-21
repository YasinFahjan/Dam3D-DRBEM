subroutine fluid
 use mesh
!
! Allocating data & working  arrays  
!  
 nn1=2*(ny+nz)+1
 nn2=nx+2*nz+1
 nelm=2*nx*(nz+ny)+2*nz*ny
 nnum=numnp
 nnod=(nz+1)*(nx+1)
 allocate(node(nn1+2,nn2+2),element(nelm))
!
! Transfer the Dam element information to the Resrvoir
!
 node(1:nz+1,nz+1:nz+nx+1)%nnum=nodes(1:nz+1,1:nx+1,mdam+1)%nnum
 node(1:nz+1,nz+1:nz+nx+1)%coord(1)=nodes(1:nz+1,1:nx+1,mdam+1)%coord(1)
 node(1:nz+1,nz+1:nz+nx+1)%coord(2)=nodes(1:nz+1,1:nx+1,mdam+1)%coord(2)
 node(1:nz+1,nz+1:nz+nx+1)%coord(3)=nodes(1:nz+1,1:nx+1,mdam+1)%coord(3)
!
! face- node generation	 1-4 ,5-6
!
 ri=corner(1)%radian(1)
 rf=corner(8)%radian(1)
!
 n1(1)=1        ; n2(1)=nz+1   ; nnx(1)=nx; nnz(1)=nz; nj(1)=nx; ni(1)=nz-1  
 n1(2)=nz+1     ; n2(2)=nz+1   ; nnx(2)=nx; nnz(2)=ny; nj(2)=nx; ni(2)=ny-1 
 n1(3)=ny+nz+1  ; n2(3)=nz+1   ; nnx(3)=nx; nnz(3)=nz; nj(3)=nx; ni(3)=nz-1 
 n1(4)=ny+2*nz+1; n2(4)=nz+1   ; nnx(4)=nx; nnz(4)=ny; nj(4)=nx; ni(4)=ny-1 
 n1(5)=nz+1     ; n2(5)=1      ; nnx(5)=nz; nnz(5)=ny; nj(5)=nz-1; ni(5)=ny-1  
 n1(6)=nz+1     ; n2(6)=nz+nx+1; nnx(6)=nz; nnz(6)=ny; nj(6)=nz-1; ni(6)=ny-1 
 do iface=2,6
 incc=(iface-1)*2
 select case(iface)
 case(2);  inc1=4; inc2=3; inc3=8; inc4=7; i1=1; j1=0
 case(3);  inc1=8; inc2=7; inc3=5; inc4=6; i1=0; j1=0
 case(4); inc1=5; inc2=6;inc3=1; inc4=2; i1=0; j1=0
 case(5); inc1=1; inc2=4; inc3=5; inc4=8; i1=1; j1=1
 case(6); inc1=3; inc2=2; inc3=7; inc4=6; i1=1; j1=1
 case default;  inc1=incc+1; inc2=incc+2; inc3=incc+3; inc4=incc+4; i1=0; j1=0
 end select
 do i=i1,ni(iface)
 do j=j1,nj(iface)
  i_index=n1(iface)+i; j_index=n2(iface)+j
  x1=corner(inc1)%radian ;x2=corner(inc2)%radian
  x3=corner(inc3)%radian; x4=corner(inc4)%radian
  xx1=x1+(x3-x1)/nnz(iface)*i; xx2=x2+(x4-x2)/nnz(iface)*i
  if(geometry==1)call contours
  xx=xx1+(xx2-xx1)/nnx(iface)*j
  radius=xx(1)
  angle =xx(2)
  height=xx(3)
  if(geometry==1)call contours
!  if(code(3)==4.and.i_index>=nz+ny+1.and.i_index<=2*nz+ny+1)then
  if(code(3)==4.and.radiation==1)then
  radius=radius/cos(angle*(radius-ri)/(rf-ri))
  endif
  node(i_index,j_index)%coord(1)=radius*cos(angle)
  node(i_index,j_index)%coord(2)=radius*sin(angle)
  node(i_index,j_index)%coord(3)=height
  x= node(i_index,j_index)%coord
 select case(nod)
  node(i_index,j_index)%id(1:3)=0
  case(4)
  nnum=nnum+1
  nnod=nnod+1
  node(i_index,j_index)%nnum=nnum
  node(i_index,j_index)%coord=matmul(T,node(i_index,j_index)%coord)
  case(8)
  if(mod(i+1,2)==0.and.mod(j+1,2)==0)then
  node(i_index,j_index)%nnum=0
  else
  nnum=nnum+1
  nnod=nnod+1
  node(i_index,j_index)%nnum=nnum
  node(i_index,j_index)%coord=matmul(T,node(i_index,j_index)%coord)
  endif
  end select

 enddo;enddo;enddo
  nnod=nnod
  ntot=nnum
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
   nelm=0
  do ii=1,6
  do i=n1(ii),nnx(ii)
  do j=n2(ii),nnz(ii)
  select case(nod)
  case(4)
   if(node(i,j)%nnum*node(i,j+1)%nnum*node(i+1,j)%nnum*node(i+1,j+1)%nnum==0)then
   else
   nelm=nelm+1
   element(nelm)%node(1)=node(i,j)%nnum
   element(nelm)%node(2)=node(i,j+1)%nnum
   element(nelm)%node(3)=node(i+1,j+1)%nnum
   element(nelm)%node(4)=node(i+1,j)%nnum
   element(nelm)%code=code(ii); element(nelm)%value =value(ii)
   element(nelm)%coord=(node(i,j)%coord+node(i,j+1)%coord+node(i+1,j+1)%coord+node(i+1,j)%coord)/4.
   endif
  case(8)
   if(node(i,j)%nnum*node(i,j+1)%nnum*node(i,j+2)%nnum*node(i+1,j)%nnum&
   &*node(i+1,j+2)%nnum*node(i+2,j)%nnum*node(i+2,j+1)%nnum*node(i+2,j+2)%nnum==0)then
   else
   nelm=nelm+1
   element(nelm)%node(1)=node(i,j)%nnum
   element(nelm)%node(2)=node(i,j+2)%nnum
   element(nelm)%node(3)=node(i+2,j+2)%nnum
   element(nelm)%node(4)=node(i+2,j)%nnum
   element(nelm)%node(5)=node(i,j+1)%nnum
   element(nelm)%node(6)=node(i+1,j+2)%nnum
   element(nelm)%node(7)=node(i+2,j+1)%nnum
   element(nelm)%node(8)=node(i+1,j)%nnum
   element(nelm)%code=code(ii); element(nnum)%value =value(ii)
   element(nelm)%coord=(node(i,j)%coord+node(i,j+1)%coord+node(i,j+2)%coord+node(i+1,j)%coord&
   &+node(i+1,j+2)%coord+node(i+2,j)%coord+node(i+2,j+1)%coord+node(i+2,j+2)%coord)/8.0
   endif
  case default
  end select 
  enddo; enddo; enddo
!
! internal points generations   
!
 inx=0; iny=0; inz=0
 if(dx_int(1)/=0)inx=int((xf_int(1)-xi_int(1))/dx_int(1))
 if(dx_int(2)/=0)iny=int((xf_int(2)-xi_int(2))/dx_int(2))
 if(dx_int(3)/=0)inz=int((xf_int(3)-xi_int(3))/dx_int(3))
 if(inx+iny+inz/=0)inx=inx+1; iny=iny+1; inz=inz+1;
 nint=inx*iny*inz
!
  open(55,file='test.dat',status='unknown')
  do i=1,nn1
  do j=1,nn2
  write(55,'(3i8)')i,j,node(i,j)%nnum
  enddo;enddo
  close(55)
! 
 contains

function theta_r(yr,xr)
 real xr,yr		     
 theta_r=19.794 -0.225017*yr +0.000994457*yr**2 -2.12983e-6*yr**3 &
 &+2.20405e-9*yr**4-8.80752e-13*yr**5 -0.0156063*xr +0.00039348*xr*yr &
 &-2.35875e-6*xr*yr**2 +5.29353e-9*xr*yr**3 -3.96588e-12*xr*yr**4 &
 &-0.000268273*xr**2 +3.72533e-6*xr**2*yr -1.147e-8*xr**2*yr**2 &
 &+9.74481e-12*xr**2*yr**3 -2.93403e-6*xr**3 +9.33748e-9*xr**3*yr &
 &-3.14911e-12*xr**3*yr**2 +5.99698e-009*xr**4 -2.3342e-11*xr**4*yr &
 &+9.52922e-12*xr**5
  end function theta_r

function theta_l(yl,xl)
 real xl,yl
 theta_l=-17.1355 +0.200978*yl -0.000923077*yl**2 +2.04823e-6*yl**3 &
 &-2.19921e-9*yl**4 +9.15841e-13*yl**5 +0.0101603*xl -0.000266593*xl*yl&
 &+1.40078e-6*xl*yl**2 -2.75017e-9*xl*yl**3 +1.84224e-12*xl*yl**4 &
 &+0.000392754*xl**2 -2.17808e-6*xl**2*yl +4.03644e-9*xl**2*yl**2 &
 &-2.56469e-12*xl**2*yl**3 -2.13971e-6*xl**3 +5.85002e-9*xl**3*yl &
 &-3.8135e-12*xl**3*yl**2 +7.61562e-9*xl**4 -9.10346e-12*xl**4*yl &
 &-1.09202e-11*xl**5
 
 end function theta_l
 
 subroutine contours
  select case(iface)
  case(2)
  xx1(2)=theta_l(xx1(1),xx1(3))
  xx2(2)=theta_r(xx2(1),xx2(3))
  case(3)
  xx1(2)=theta_l(xx1(1),xx1(3))
  xx2(2)=theta_r(xx2(1),xx2(3))
  case(4)
  xx1(2)=theta_l(xx1(1),xx1(3))
  xx2(2)=theta_r(xx2(1),xx2(3))
  case(5)
  angle=theta_l(radius,height)
  case(6)
  angle=theta_r(radius,height)
  case default
  end select

 end subroutine contours


end subroutine fluid