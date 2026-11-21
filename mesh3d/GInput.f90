subroutine input
 use arrays
! 
! read names and open  files for input and output    
! 
 write(*,' (a) ')' name of input file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(inp,file=file,status='old')
 write(*,' (a) ') ' name of output file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(ipr,file=file,status='unknown')
 open(13,file='boxn.dat',status='unknown')
 open(14,file='boxe.dat',status='unknown')
! 
! read job title 
! 
 read(inp,'(a)') title
 write(ipr,'(a)') title
! 
! nx  : number of nodes in x- direction
! ny  : number of nodes in y- direction
! nz  : number of nodes in z- direction
! 
 read(inp,*)nx,ny,nz
  nn1=2*(ny+nz)+1
  nn2=nx+2*nz+1
  nelm=2*nx*(nz+ny)+2*nz*ny
  if(nod.eq.3)nelm=2*nelm
! 
! allocating data and working arrays
! 
 call array_allocations
! 
! read coordinates of fluid domain boundary nodes
! 
 do i=1,8
 read(inp,*)corner(i)%x
 enddo
! 
! read boundary element information 
!  kode  : boundary element type (1: p known, 2: q known) 
!  value : 
 do i=1,6
 read(inp,*)kode(i),value(i)
 enddo
end subroutine input
