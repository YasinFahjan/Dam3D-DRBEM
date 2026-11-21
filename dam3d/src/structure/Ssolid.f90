subroutine Ssolid
use fem3D
implicit none
!                                                                     
! FEM3D : Formation of Stifness and mass matrices
!
 print*,'FEM-3D: Form the FEM [M],[K],[C]'
!
 ndr=nods+rots   ! total displacement and rotations
! 
 allocate(mass(3*nods,3*nods),stiff(3*ndr,3*ndr))
 if(rots>0)allocate(stiff_alfa(3*rots,3*rots))
 allocate(lm(3*nods),nn(nods),id(nods,3))
 allocate(N_funct(3,3*nods),G(6,3*ndr),D(6,6))
!
 allocate(funct(nods+5),derv_eta(3,nods+5))
 allocate(coord(nods,3),derv_x(3,nods+5))
 allocate(N_Der(3,nods))  
!
 allocate(mass_mat(neq,neq),stiff_mat(neq,neq))

 stiff_mat=0.0; mass_mat=0.0

 vol=0.0
 do nel=1,numel 
!                                                                     
! FEM3D\element: nodes, material, LM vector                                                
!                                                                     
 stiff=0.0; mass=0.0; N_funct=0.0; G=0.0; D=0.0 
 mat=elements(nel)%mat
 ng=elements(nel)%ng
 nn(1:nods)=elements(nel)%node(1:nods)
 do ij=1,nods
 coord(ij,1:3)=node(nn(ij))%coord(1:3)
 lm(3*ij-2:3*ij)=node(nn(ij))%id(1:3)
 enddo
!                                                                     
! FEM3D\Element: stress-strain law [D]                              
!                                                                     
 Stiff_coeff=ee(mat)/((1.-2.*enu(mat))*(1.+enu(mat)))
 D(1:3,1:3)=Stiff_coeff*(enu(mat))
 do i=1,6
 if(i<=3)D(i,i)=Stiff_coeff*(1.-enu(mat))
 if(i> 3)D(i,i)=Stiff_coeff*(.5-enu(mat))
 enddo
!
! FEM3D: classical Gauss-Legendre quadrature formula using imsl
!
 allocate(xk1(ng),xk2(ng),xk3(ng-1))
 allocate(wgt1(ng),wgt2(ng),wgt3(ng-1))
 call gqrul (ng,1,0.,0.,0,(/0.,0./),xk1,wgt1) 
 call gqrul (ng,1,0.,0.,0,(/0.,0./),xk2,wgt2) 
 call gqrul ((ng-1),1,0.,0.,0,(/0.,0./),xk3,wgt3) 
!
! FEM3D\Element: loop over ng**3 integration points
!                                                                     
 do lx = 1,ng                                                       
 do ly = 1,ng                                                       
 do lz = 1,(ng-1)                                                     
 ee1=xk1(lx); ee2=xk2(ly); ee3=xk3(lz)                                  
 wg=wgt1(lx)*wgt2(ly)*wgt3(lz)                                
!
! FEM3D\Element: shape functions for isoparametric 3D thick shell element
!
 select case(nods)
 case( 8); call solid_funct(ee1,ee2,ee3)
 case(16); call shell_funct(ee1,ee2,ee3)
 end select                                                    
!
! FEM3D\Element: shape function matrix [N] and Strain Matrix [G]
!
 do i=1,3 
 N_funct(i,i:3*nods:3)=funct(1:nods)
 G(i,i:3*ndr:3)=derv_x(i,1:ndr)
 enddo
 G(4,1:3*ndr:3)=derv_x(2,1:ndr); G(4,2:3*ndr:3)=derv_x(1,1:ndr)
 G(5,2:3*ndr:3)=derv_x(3,1:ndr); G(5,3:3*ndr:3)=derv_x(2,1:ndr)
 G(6,1:3*ndr:3)=derv_x(3,1:ndr); G(6,3:3*ndr:3)=derv_x(1,1:ndr)
!                                                                     
! FEM3D\Element: consistent mass matrix                                              
!                                                                     
 mass(1:3*nods,1:3*nods)=mass(1:3*nods,1:3*nods)+ &
   &(wg*det_jac)*rho_d(mat)/gravity*matmul(transpose(N_funct),N_funct)
 vol=vol+wg*det_jac                                                        
!                                                                     
! FEM3D\Element: General stiffness matrix (displacement+ rotations)                                
!                                                                     
 stiff=stiff+(wg*det_jac)*matmul(transpose(G),matmul(D,G)) 
 enddo; enddo; enddo
!                                                                     
! FEM3D\Element Stiffness : static condensation of the rotations                                                
!                                                                     
 if(rots>0)then
 stiff_alfa(1:3*rots,1:3*rots)=stiff(3*nods+1:3*ndr,3*nods+1:3*ndr)
 call linds (3*rots,stiff_alfa,3*rots,stiff_alfa,3*rots)
 stiff(1:3*nods,1:3*nods)= stiff(1:3*nods,1:3*nods)- &
      & matmul(stiff(1:3*nods,3*nods+1:3*ndr), &
	& matmul(stiff_alfa,stiff(3*nods+1:3*ndr,1:3*nods)))
 endif
!
! FEM3D\Global: cosistent mass and stiffness matrices
!
 do i=1,3*nods
 do j=1,3*nods
 if(lm(i)/=0.and.lm(j)/=0)then
 mass_mat(lm(i),lm(j))=mass_mat(lm(i),lm(j))+mass(i,j)
 stiff_mat(lm(i),lm(j))=stiff_mat(lm(i),lm(j))+stiff(i,j)
 endif
 enddo
 enddo
!
 deallocate(xk1,xk2,xk3,wgt1,wgt2,wgt3)
!                   
! FEM3D\print for check
!                                                                     
 write(*,*)' building element stifness matrix # ',nel                                          
! write(prn,'(i6,x,<nods>i5,2x,i9,i12)')nel,elements(nel)%node(1:nods),ng,mat 

 enddo                                                                    

 write(prn,'(//5x,e15.7)') vol
 write(prn,'(5x,e15.7)') vol*rho_d(1)

 deallocate(funct,derv_eta,coord,derv_x)
 deallocate(mass,stiff)
 deallocate(N_funct,G,D,lm,nn,id)
 deallocate(N_Der)  
 if(rots>0)deallocate(stiff_alfa)
!
!open(20,file='m.dat')
!open(30,file='k.dat')
!write(20,'(<neq>g20.7)')((mass_mat(i,j),j=1,neq),i=1,neq)
!write(30,'(<neq>g20.7)')((stiff_mat(i,j),j=1,neq),i=1,neq)
!close(20)
!close(30)

end subroutine Ssolid	                                                  
