subroutine rad
 use mesh
!
! Allocating data & working  arrays  
!  
 nelr=2*(nx+nz)
 allocate(elementr(nelr))
!
! Generation of 2-nodes quadilateral elementr information 
!
   n1(3)=nz+ny+1   ; nnx(3)=2*nz+ny ; n2(3)=nz+1     ; nnz(3)=nz+nx
   nelr=0
  do i=n1(3),nnx(3)
   j=n2(3)
  select case(ndr)
  case(2)
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i,j)%nnum
   elementr(nelr)%node(2)=node(i+1,j)%nnum
   elementr(nelr)%code=code(5)
  case(3)
   if(mod(i,2)==0)then
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i-1,j)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%node(3)=node(i+1,j)%nnum
   elementr(nelr)%code=code(5)
   endif
  end select 
  enddo 
! 
  do j=n2(3),nnz(3)
     i=nnx(3)+1
  select case(ndr)
  case(2)
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i,j)%nnum
   elementr(nelr)%node(2)=node(i,j+1)%nnum
   elementr(nelr)%code=code(5)
  case(3)
   if(mod(j,2)==0)then
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i,j-1)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%node(3)=node(i,j+1)%nnum
   elementr(nelr)%code=code(5)
   endif
  end select 
  enddo 
!
  do i=nnx(3),n1(3),-1
   j=nnz(3)+1
  select case(ndr)
  case(2)
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i+1,j)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%code=code(4)
  case(3)
   if(mod(i,2)==0)then
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i+1,j)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%node(3)=node(i-1,j)%nnum
   elementr(nelr)%code=code(4)
   endif
  end select 
  enddo 
! 
  do j=nnz(3),n2(3),-1
     i=n1(3)
  select case(ndr)
  case(2)
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i,j+1)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%code=code(2)
  case(3)
   if(mod(j,2)==0)then
   nelr=nelr+1
   elementr(nelr)%node(1)=node(i,j+1)%nnum
   elementr(nelr)%node(2)=node(i,j)%nnum
   elementr(nelr)%node(3)=node(i,j-1)%nnum
   elementr(nelr)%code=code(2)
   endif
  end select 
  enddo 
 !
end subroutine rad