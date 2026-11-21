subroutine input
! Linear BEM3D  
 use arrays
 use numeric
  pi=4.*atan(1.0)
! 
! read names and open  files for input and output    
! 
 write(*,' (a) ')' name of input file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(inp,file=file,status='old')
 write(*,' (a) ') ' name of output file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(prn,file=file,status='unknown')
! 
! read job title 
! 
 read(inp,'(a)') title
 write(prn,'(a/)') title
! 
! print date and time
! 
 write(prn,'(6hdate :,a\)')date()
 write(prn,'(10x,6hTime :,a)')clock()
! 
! reading  # nodes, elements, internal points, # nodes/element
! 
 read(inp,*)nnod,nelm,nint,nod
! 
! allocating data and working arrays
! 
 call arrays_allocation
!
! initialize vectors for nodes & elements
! 
  element%code=0;node%code=2; node%p=0.
! 
! read coordinates of boundary nodes
! 
 do i=1,nnod
  read(inp,*)node(i)%x
 enddo
! 
! read boundary element information 
!  
 do i=1,nelm
  read(inp,*)(n(j),j=1,nod),code,(value(j),j=1,nod)
   element(i)%node(1:nod)=n(1:nod) 
   element(i)%code=code
  select case(code)
   case(1)
    element(i)%p(1:nod)=value(1:nod)
    node(n(1:nod))%p=value(1:nod)
    node(n(1:nod))%code=code
   case(2)
    element(i)%q(1:nod)=value(1:nod)
  end select
 enddo
! 
! read coordinates of  internal nodes 
! 
 do i=nnod+1,nnod+nint
  read(inp,*) node(i)%x
 enddo
! 
! print input data 
! 
 write(prn,'(/a,i5)')'number of boundary nodes    =',nnod
 write(prn,'(a,i5)')'number of boundary elements =',nelm
 write(prn,'(a,i5)')'number of internal points   =',nint
 write(prn,'(a,i5)')'number of nodes per element =',nod
 write(prn,'(//a/)')'> data for boundary nodes'
 write(prn,'(2x,5(a,13x)/)')'node','x','y','z','p' 
 write(prn,'(i5,3x,4f14.6)')(i,node(i)%x,node(i)%p,i=1,nnod)
 write(prn,'(//a/)')'> data for boundary elements'
 write(prn,'(a,6x,<nod>(a,i1,3x),a,8x,<nod>(a,i1,12x)/)')&
     &'element',('node',i,i=1,nod),'code', ('q',i,i=1,nod)
 write(prn,'(i5,3x,<nod+1>i8,<nod>f14.6)')(i,element(i)%&
     &node(1:nod),element(i)%code,element(i)%q(1:nod),i=1,nelm)
 if(nint.gt.0)then
 write(prn,'(/a/)')'> coordinates of internal nodes' 
 write(prn,'(2x,4(a,13x)/)')'node','x','y','z'  
 write(prn,'(i5,3x,3f14.6)')(i,node(i)%x,i=nnod+1,nnod+nint)
 endif
end subroutine input


