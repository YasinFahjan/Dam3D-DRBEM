subroutine HInput
use Bem3D
implicit none
!                                                                               
! BEM3D: Input and Control Data of the Fluid Elements                                                      
! 
 read(inp,*) void
!                                                                     
! BEM3D: Read Elements: nodes, code, absorption                                                 
!                                                                               
 read(inp,*)nelm
 read(inp,*)nod
!
 allocate(element(nelm))
!
 allocate(funct(nod),derv_eta(3,nod))
 allocate(coord(nod,3),derv_x(3,nod))
 allocate(N_funct(3,nod),N_Der(3,nod))  
 allocate(eta_c(3,nod))
!  
 element%code=0;node%code=0; node%p=0; 
!  
! code=0    internal point    
! code=1    free surface b.c.    
! code=2    dam interaction b.c. 
! code=3    bottom absorbing b.c.
! code=4    scattering wave b.c. 
!
 call CornerFunctions(nod)
!
 do i=1,nelm
 read(inp,*)(n(j),j=1,nod),code,(value(j),j=1,nod)
 element(i)%node(1:nod)=n(1:nod) 
 element(i)%code=code
 element(i)%p=0.0; element(i)%q=0.0
 call Hcode
 enddo
!                                                                               
! BEM3D: Number of Equation ( G & H & S Matrices Dimensions)                                                              
!                                                                               
 nnod=0
 do i=1,ntot
 if(node(i)%code/=0)nnod=nnod+1
 enddo
 nq=nnod+nint
! 
! BEM3D: Read Coordinates of Internal Nodes 
! 
 read(inp,*)void
 do i=ntot+1,ntot+nint							 
 read(inp,*)node(i)%coord
 enddo
!
! BEM3D: Defining ID for the Reservoir Nodes 
!
 allocate(lmw(nq))
 counter=0; node%idw=0; np=0; lmw=0
 do code=1,4
 do i=1,ntot
 if(node(i)%code==code)then
 np(code)=np(code)+1
 counter=counter+1
 lmw(counter)=i
 node(i)%idw=counter
 endif
 enddo; enddo
! intenal points
 do i=ntot+1,ntot+nint
 counter=counter+1
 lmw(counter)=i
 node(i)%idw=counter
 enddo
 np1=np(1);	np2=np(2)  ; np3=np(3)  ; np4=np(4)  
 nf1=np1  ; nf2=nf1+np2; nf3=nf2+np3; nf4=nf3+np4
 nqr=np4 
!
write(prn,'(a,i10)/)')'np1',np1  
write(prn,'(a,i10)/)')'np2',np2  
write(prn,'(a,i10)/)')'np3',np3  
write(prn,'(a,i10)/)')'np4',np4
write(prn,'(a,i10)/)')'nf1',nf1  
write(prn,'(a,i10)/)')'nf2',nf2  
write(prn,'(a,i10)/)')'nf3',nf3  
write(prn,'(a,i10)/)')'nf4',nf4
  
write(prn,'(2x,5(a,13x)/)')'node','x','y','z','code'  
do counter=1,nf4+nint
i=lmw(counter)
write(prn,'(i5,3f8.3,3x,4i6)')i,node(i)%coord,node(i)%code,node(i)%idw,node(i)%idr,node(i)%codr
enddo
deallocate(funct,derv_eta)
deallocate(coord,derv_x)
deallocate(N_funct,N_Der)  
deallocate(eta_c)
! 
! BEM3D: Print Input Data  - Nodes Informartion
! 
 write(prn,'(a//)')     '[ HYDRODYNAMIC MODEL ] '
 write(prn,'(a,i5)')   'BEM3D: number of Boundary nodes...=',nnod
 write(prn,'(a,i5)')   'BEM3D: number of Internal nodes...=',nint
 write(prn,'(a,i5)')   'BEM3D: number of elements ....... =',nelm
 write(prn,'(a,i5)')   'BEM3D: number of nodes\elements.. =',nod
 write(prn,'(a,f10.5/)')'BEM3D: absorbing coeffecient......=',alfa_r
 write(prn,'(a,f10.5)')'DRM3D: Approx.Funct.Coeff.,  c0...=',c(1)
 write(prn,'(a,f10.5)')'DRM3D: Approx.Funct.Coeff.,  c1...=',c(2)
 write(prn,'(a,f10.5)')'DRM3D: Approx.Funct.Coeff.,  c2...=',c(3)
 write(prn,'(a,f10.5)')'DRM3D: Approx.Funct.Coeff.,  c3...=',c(4)
!
 write(prn,'(//a/)')    'BEM3D: data for boundary nodes'
 write(prn,'(2x,4(a,13x),a/)')'node','x','y','z','code'
 do i=1,ntot
 if(node(i)%code/=0)then
 write(prn,'(i5,3x,3f14.6,3x,i10)')i,node(i)%coord,node(i)%code
 endif; enddo
 if(nint.gt.0)then
 write(prn,'(//a/)')    'BEM3D: coordinates of internal nodes' 
 write(prn,'(2x,4(a,13x)/)')'node','x','y','z'  
 write(prn,'(i5,3x,3f14.6)')(i,node(i)%coord,i=ntot+1,ntot+nint)
 endif
!
! BEM3D: Print Input Data - Element Informartion
!
 write(prn,'(//a/)')    'BEM3D: element informartion' 
 write(prn,'(a,6x,<nod>(a,i1,3x),a,8x,<2*nod>(a,12x)/)')&
 &'element',('node',i,i=1,nod),'code',('p [real,complex]',i=1,nod),&
 &('q [real,complex]',i=1,nod)
 write(prn,'(i5,3x,<nod+1>i8,<7*nod>f14.6)')(i,element(i)%node(1:nod),&
 &element(i)%code,element(i)%p(1:nod),element(i)%q(1:nod),element(i)%norm(1:nod,1:3),i=1,nelm)
 write(prn,'(//90(1h=)/)')
!
 contains
! 
subroutine Hcode
implicit none
!
! DAM3D: Defining the Normal Acceleration on the Element  
!
 n(1:nod)=element(i)%node(1:nod) 
 coord=0.0; norm=0.0
 do ij=1,nod; coord(ij,:)=node(n(ij))%coord(:); enddo
 do ij=1,nod
 call ShapeFunctions(nod,eta_c(1,ij),eta_c(2,ij))
 element(i)%norm(ij,1:3)=norm(1:3)
 if(code==2.or.code==3)value(ij)=-rho_w/gravity*sum(ag*norm)
! if(code==2.or.code==3)value(ij)=-rho_w/gravity*sum(ag*norm)
 enddo
!
! DAM3D: Assigning Code to the Nodes on the Element  
!
 select case(code)
 case(1)
 element(i)%p(1:nod)=value(1:nod)
 node(n(1:nod))%p=value(1:nod)
 node(n(1:nod))%code=code
 case(2)
 element(i)%q(1:nod)=value(1:nod)
 do ij=1,nod; if(node(n(ij))%code/=1.and.node(n(ij))%code/=4)&
 &node(n(ij))%code=code; enddo
 case(3)
 element(i)%q(1:nod)=value(1:nod)
 do ij=1,nod; if(node(n(ij))%code==0)node(n(ij))%code=code; enddo
 case(4)
 element(i)%q(1:nod)=value(1:nod)
 do ij=1,nod; if(node(n(ij))%code/=1)node(n(ij))%code=code; enddo
 end select
end subroutine Hcode
!
end subroutine HInput
