subroutine Input
use Dam3D
implicit none
!                                                                               
! Dam3D: Input and control data for structure-Fluid system
! 
 pi=4.*atan(1.0)
! 
! Dam3D: input and output files names     
! 
 write(*,' (a) ')' name of input file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(inp,file=file,status='old')
 write(*,' (a) ') ' name of output file (max. 12 chart.)'
 read(*,'  (a) ')file
 open(prn,file=file,status='unknown')
! 
 print*,' Dam-3D: input data'
!
! Dam3D: read job title 
! 
 read(inp,'(a)')heading
 write(prn,'(a/)')heading
! 
! Dam3D: print date and time
! 
 write(prn,'(6hdate :,a\)')date()
 write(prn,'(10x,6hTime :,a)')clock()
!
! Dam3D: reading  General data				   
! 									   
 read(inp,*)gravity						   
 read(inp,*)rho_w
 read(inp,*)sound
 read(inp,*)hwater
 read(inp,*)alfa_r
 read(inp,*)nf_i, nf_f, nf_s
 read(inp,*)ndyn
 read(inp,*)run_type
 !
 nf=(nf_f-nf_i)/nf_s+1
 gama=(1.-alfa_r)/(1.+alfa_r) 
 !
 if(ndyn==0)heading='  Free Vibration'
 if(ndyn==1)heading='  Dam Only - no water'
 if(ndyn==2)heading='  Dam & Reservoir: Full interaction'
 if(ndyn==3)heading='  Reservoir: Dam Deformational motion  '
 if(ndyn==4)heading='  Reservoir: Rigid Dam Motion'
!
! Dam3D: reading # total, structure, fluid boundary and internal points
! 
 read(inp,*)ntot
 read(inp,*)numnp
 read(inp,*)nnod
 read(inp,*)nint
!
! Dam3D: reading input acceleration
!
 read(inp,*)ag
!
 if(sum(ag*(/1.,0.,0./))/=0.)ixyz=1
 if(sum(ag*(/0.,1.,0./))/=0.)ixyz=2
 if(sum(ag*(/0.,0.,1./))/=0.)ixyz=3
! 
! Dam3D: read coordinates of nodal points
! 
 allocate(node(ntot+nint))
 read(inp,*)void
 do i=1,ntot
 read(inp,*)node(i)%coord(1:3),node(i)%id(1:3)
 enddo
!
! Dam3D: print input data - nodes informartion
! 
 write(prn,'(//72(1h=)/)')
 write(prn,'(a//)')     '[ FlUID-STRUCTURE INTERACTION ] '
 write(prn,'(a,a/)')   'Dam3D:',heading
 write(prn,'(a,i5)')   'Dam3D: number of total nodes......=',ntot
 write(prn,'(a,i5)')   'Dam3D: number of Structure nodes..=',numnp
 write(prn,'(a,i5)')   'Dam3D: Boundary Fluid nodes.......=',nnod
 write(prn,'(a,i5/)')  'Dam3D: number of modal shapes.....=',nf
 write(prn,'(a,f10.5)')'Dam3D: accel. of the Gravity......=',gravity
 write(prn,'(a,f10.5)')'Dam3D: density of the water.......=',rho_w
 write(prn,'(a,f10.5)')'Dam3D: sound speed in the watery..=',sound
 write(prn,'(//a/)')    'Dam3D: data for nodal coordinates '
 write(prn,'(2x,4(a,13x)/)')'node','x','y','z'
 write(prn,'(i5,3x,3f14.6)')(i,node(i)%coord,i=1,ntot)
!
! Dam3D : Elements in FEM & BEM
!
 call SInput
 call HInput
 if(nqr>0)call RInput
 call GInput
!
end subroutine input