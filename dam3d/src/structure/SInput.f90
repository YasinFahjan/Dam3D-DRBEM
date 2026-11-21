subroutine SInput
use Fem3D
implicit none
!                                                                               
! FEM3D: number of DOF (Global Stiff & Mass Matrices Dimensions)                                                              
!                                                                               
 neq = 0;                                                                  
 do i = 1,numnp                                                       
 do j = 1,3                                                            
 select case (node(i)%id(j))
 case(0); neq = neq + 1; node(i)%id(j) = neq                                                             
 case default; node(i)%id(j) = 0
 end select                                                               
 enddo; enddo                                                           
!                                                                               
! FEM3D: input and control data of the structure elements                                                      
! 
 read(inp,*) void
!                                                                     
! FEM3D: read elements: nodes, material, number of integration points                                                
!                                                                     
 read(inp,*)numel
 read(inp,*)nods,rots
 read(inp,*)nmat
!
 allocate(elements(numel))
 allocate(ee(nmat),enu(nmat),rho_d(nmat))
!
 read(inp,*) (ee(i),enu(i),rho_d(i),i=1,nmat)
!
 do nel=1,numel 
 if(ndyn<=3)write(*,*)' Reading element # ',nel                                          
 read(inp,*)elements(nel)%ng,elements(nel)%mat                   
 read(inp,*)elements(nel)%node(1:nods)
 enddo                                                                    
!
! FEM3D: print input data - nodes informartion
! 
 write(prn,'(//72(1h=))')
 write(prn,'(/a//)')     '[ STRUCTURE MODEL ] '
 write(prn,'(a,i5)')   'FEM3D: number of nodes............=',numnp
 write(prn,'(a,i5)')   'FEM3D: number of modal shapes.....=',nf
 write(prn,'(a,i5)')   'FEM3D: number of DOF..............=',neq
 write(prn,'(a,i5)')   'FEM3D: number of elements ....... =',numel
 write(prn,'(a,i5)')   'FEM3D: number of materials .......=',nmat
 write(prn,'(//a/)')   'FEM3D: data for nodal points '
 write(prn,'(2x,4(a,13x),3(a,9x)/)')'node','x','y','z','idx','idy','idz'
 write(prn,'(i5,3x,3f14.6,1x,3i12)')(i,node(i)%coord,node(i)%id,i=1,numnp)
!
! FEM3D: print input data - element informartion
!
 write(prn,'(//a/)')   'FEM3D: data for material property '
 write(prn,'(2x,a,6x,a,12x,a,10x,a,/ )')'material no.','E','nu','rho'      
 write(prn,'(1x,i5,4x,3g15.4)')(mat,ee(mat),enu(mat),rho_d(mat),mat=1,nmat)

 write(prn,'(//a/)')   'FEM3D: data for elements '
 write(prn,'(a,6x,<nods>(a,i2,2x),4x,2(a,2x)/)')&
 &'element',('node',i,i=1,nods),'material type','Integration Order'
 write(prn,'(i5,3x,<nods>i8,2i16)')(nel,elements(nel)%node(1:nods),&
 &elements(nel)%mat,elements(nel)%ng,nel=1,numel)
 write(prn,'(//90(1h=)/)')
!
end subroutine SInput

