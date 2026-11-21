subroutine RInput
use Bem3D
use Rad3D
implicit none
!                                                                               
! DRM2D\RAD: Input and Control Data of the Fluid Elements                                                      
! 
 read(inp,*) void
!                                                                     
! DRM2D\RAD: Read Elements: nodes, code, absorption                                                 
!                                                                               
 read(inp,*)nelr
 read(inp,*)ndr
!
 allocate(elementr(nelr))
!
 allocate(funct(ndr),derv_eta(3,ndr))
 allocate(coord(ndr,3),derv_x(3,ndr))
 allocate(N_funct(3,ndr),N_Der(3,ndr))  
 allocate(eta_c(3,ndr))
!  
 elementr%code=0; node%codr=0  
!
! code=0    internal point    
! code=1    free surface b.c.    
! code=3    bottom absorbing b.c.
!
 call CornerFunctions(ndr)
!
 do i=1,nelr
 read(inp,*)(n(j),j=1,ndr),code
 elementr(i)%node(1:ndr)=n(1:ndr) 
 elementr(i)%code=code
 elementr(i)%p=0.0; elementr(i)%q=0.0
 call Rcode
 enddo
!
! DRM2D\RAD: Defining ID for the Radiation Model 
!
 allocate(lmr(nq))
 lmr=0; node%idr=0; neqr=0 ; np=0 
 do code=1,3,2
 do i=1,ntot
 if(node(i)%codr==code)then
 np(code)=np(code)+1
 neqr=neqr+1
 lmr(neqr)=i
 node(i)%idr=neqr
 endif
 enddo; enddo
! intenal points
 do i=1,ntot
 if(node(i)%code==4.and.node(i)%codr==0)then
 neqr=neqr+1
 lmr(neqr)=i
 node(i)%idr=neqr
 endif
 enddo
!
 nr1=np(1);	nr2=np(3)   
 nk1=nr1  ; nk2=nk1+nr2
 nqr=neqr-nk1
!
 deallocate(funct,derv_eta)
 deallocate(coord,derv_x)
 deallocate(N_funct,N_Der)  
 deallocate(eta_c)
! 
! DRM2D\RAD: Print Input Data  - Nodes Informartion
! 
 write(prn,'(/a//)')     '[ RADIATION MODEL ] '
 write(prn,'(a,i5)')   'DRM2D\RAD: number of Boundary nodes...=',nk2
 write(prn,'(a,i5)')   'DRM2D\RAD: number of Internal nodes...=',neqr-nk2
 write(prn,'(a,i5)')   'DRM2D\RAD: number of elements ....... =',nelr
 write(prn,'(a,i5/)')  'DRM2D\RAD: number of nodes\elements.. =',ndr
!
 write(prn,'(//a/)')    'DRM2D\RAD: data for boundary and internal nodes'
 write(prn,'(2x,4(a,13x),a/)')'node','x','y','z','code'
 do i=1,neqr
 write(prn,'(i5,3x,3f14.6,3x,i10)')lmr(i),node(lmr(i))%coord,node(lmr(i))%codr
 enddo
!
! DRM2D\RAD: Print Input Data - Element Informartion
!
 write(prn,'(//a/)')    'DRM2D\RAD: element informartion' 
 write(prn,'(a,6x,<ndr>(a,i1,3x),a,8x,<2*ndr>(a,12x)/)')&
 &'element',('node',i,i=1,ndr),'code',('q [real,complex]',i=1,ndr),&
 &(' normal ',i=1,ndr)
 write(prn,'(i5,3x,<ndr+1>i8,<5*ndr>f14.6)')(i,elementr(i)%node(1:ndr),&
 &elementr(i)%code,elementr(i)%q(1:ndr),(elementr(i)%norm(j,1:3),j=1,ndr),i=1,nelr)
 write(prn,'(//90(1h=)/)')
!
 contains
! 
subroutine Rcode
implicit none
!
! DRM2D\RAD: Defining the Normal to the Element  
!
 n=0; n(1:ndr)=elementr(i)%node(1:ndr) 
 coord=0.0; norm=0.0
 do ij=1,ndr; coord(ij,:)=node(n(ij))%coord(:); enddo
 do ij=1,ndr
 call ShapeFunctions(ndr,eta_c(1,ij),eta_c(2,ij))
 elementr(i)%norm(ij,1:3)=norm(1:3)
 if(code==2.or.code==3)value(ij)=-rho_w/gravity*sum(ag*norm)
 enddo
!
! DRM2D\RAD: Assigning Code to the Nodes on the Element  
!
 select case(code)
 case(1)
 node(n(1:ndr))%codr=code
 case(3)
 elementr(i)%q(1:ndr)=value(1:ndr)
 do ij=1,ndr; if(node(n(ij))%codr==0)node(n(ij))%codr=code; enddo
 end select
end subroutine Rcode
!
end subroutine RInput
