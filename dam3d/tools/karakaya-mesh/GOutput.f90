subroutine Output
 use mesh 
!
! This subroutine to prepare the Input File for DAM3D
!
 write(prn,'(g12.4,5x,a)')gravity,': gravity'
 write(prn,'(g12.4,5x,a)')rho_w,': specific weight of the water'
 write(prn,'(g12.4,5x,a)')sound,': sound speed in water'
 write(prn,'(g12.4,5x,a)')hwater,':water height in the reservoir'
 write(prn,'(g12.4,5x,a)')alfa_r,': reflection coeffecient'
 write(prn,'(3i4,5x,a)')1,nf,1,': # req. mode shapes range'
 write(prn,'(i7,10x,a)')ndyn,': type of problem'
 write(prn,'(i7,10x,a)')0,': 0:Run 1:Run&Save 2:Load&Run'
!
 write(prn,'(/i7,10x,a)')ntot,': # total nodes'
 write(prn,'(i7,10x,a)')numnp,': # dam nodes'
 write(prn,'(i7,10x,a)')nnod,': # reservoir nodes'
 write(prn,'(i7,10x,a/)')nint,': # Reservoir Internal nodes'

 write(prn,'(3g12.4,5x,a/)')ag,': ground acceleration'
!
! print nodes coordinates for dam structure
! 
 write(prn,'(/a/)')' Dam-Reservoir System nodes coordinates:'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 write(prn,'(3f12.5,5x,3i8)')nodes(i,j,k)%coord,nodes(i,j,k)%id(1:3)
 write(prn4,'(i7,1x,3f12.5)')nnum,nodes(i,j,k)%coord
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo
!
 write(prn,*)
!
! print nodes coordinates for reservoir
! 
 do nnum=numnp+1,ntot
 indicator=' '
 do i=1,nn1
 do j=1,nn2
 if(indicator==' '.and.node(i,j)%nnum==nnum)then
 write(prn,'(3f12.5,5x,3i8)')node(i,j)%coord,node(i,j)%id(1:3)
 indicator='*'
 endif
 enddo ;enddo ;enddo

 do nnum=1,ntot
 indicator=' '
 do i=1,nn1
 do j=1,nn2
 if(indicator==' '.and.node(i,j)%nnum==nnum)then
 write(prn2,'(i7,1x,3f12.5)')nnum,node(i,j)%coord
 indicator='*'
 endif
 enddo ;enddo ;enddo
!
! print dam structure elements
!
 write(prn,'(/a)')' Dam Structure Elements'
 write(prn,'(i8,10x,a)')numel,': # Dam Elements'
 write(prn,'(2i8,10x,a)')2*nod,nrot,': # nodes per Element'
 write(prn,'(i8,10x,a/)')nmat, ': # Dam materials'
 write(prn,'(3(g14.5,5x))')(ee,enu,rho_d,mat=1,nmat)
 write(prn,*)
 do i=1,numel
 write(prn,'(2i8)')elements(i)%ng,elements(i)%mat
 write(prn,'(8i8:8i8)')(elements(i)%node(j),j=1,2*nod)
 write(prn5,'(i3,5x,16i5)')i,(elements(i)%node(j),j=1,2*nod)
 enddo
!
! print reservoir elements
!
 write(prn,'(/a)')' Resrvoir Elements'
 write(prn,'(i8,10x,a)')nelm,': # Reservoir Elements Elements'
 write(prn,'(i8,10x,a/)')nod,  ': # nodes per Element'
 do i=1,nelm
 write(prn,'(4i8\:4i8\)')(element(i)%node(j),j=1,nod)
 write(prn,'(i8,4f10.5:4f10.5)')element(i)%code,(element(i)%value,j=1,nod)
! write(prn3,'(i3,5x,4i5:4i5)')i,(element(i)%node(j),j=1,nod)
 write(prn3,'(i7,1x,3f10.5)')i,element(i)%coord
 enddo
!
! internal points generations   
!
 write(prn,'(/a/)')'Internal points'
 do k=1,inz
 do j=1,iny
 do i=1,inx
 radius=xi_int(1)+(i-1)*dx_int(1)
 angle =xi_int(2)+(j-1)*dx_int(2)
 height=xi_int(3)+(k-1)*dx_int(3)
 if(code(3)==4.and.radiation==1)then
 radius=radius/cos(angle*(radius-ri)/(rf-ri))
 endif
 x(1)=radius*cos(angle)
 x(2)=radius*sin(angle)
 x(3)=height
 write(prn,'(3f12.3)')x
 enddo; enddo; enddo
!
! print radiation elements
!
 if(code(3)==4)then
 write(prn,'(/a)')' Radiaiton Elements'
 write(prn,'(i8,10x,a)')nelr,': # Reservoir Elements Elements'
 write(prn,'(i8,10x,a/)')ndr  ,': # nodes per Element'
 do i=1,nelr
 write(prn,'(2i8\:i8\)')(elementr(i)%node(j),j=1,ndr)
 write(prn,'(i8)')elementr(i)%code
 enddo
 endif
!
!
! print Graphics nodes
!
 if(Gmode==1)then
!
 write(prn,'(/a)')'Graphs: Upstream Face Nodes '
 write(prn,'(i4,2x,i4,10x,a)')nz+1,nx+1,': # nodes across face'
 do i=1,nz+1
 write(prn,'(<nx+1>(i4,3x))')(node(i,j)%nnum,j=nz+1,nz+nx+1)
 enddo
!
 write(prn,'(/a)')'Graphs: Reservoir Bottom Nodes '
 write(prn,'(i4,2x,i4,10x,a)')ny+1,nx+1,': # nodes across face'
 do i=2*(nz+ny)+1,2*nz+ny+1,-1
 write(prn,'(<nx+1>(i4,3x))')(node(i,j)%nnum,j=nz+1,nz+nx+1)
 enddo
!
 endif
!
 close (inp)
 close (ipr)
!
! print Simple Dam input
!
 if(porter==1.and.nod==8)then
 open(25,file='porter',status='unknown')
!
 write(25,'(a)') title
 write(25,'(i7,10x,a)')nf,': # req. mode shapes '
 write(25,'(i7,10x,a)')ndyn,': type of problem'
 write(25,'(i7,10x,a)')numnp,': # dam nodes'
 write(25,'(2(f8.3,2x))')corners(8)%radian(1)&
 &,corners(8)%radian(1)/corners(8)%radian(3)
!
! print nodes coordinates for dam structure
! 
 write(25,'(/a/)')' Dam System nodes coordinates:'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 write(25,'(5f12.5,5x,3i8)')nodes(i,j,k)%angle,nodes(i,j,k)%height,&
 &nodes(i,j,k)%coord,nodes(i,j,k)%id(1:3)
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo
!
! print dam structure elements
!
 write(25,'(/a)')' Dam Structure Elements'
 write(25,'(i8,10x,a)')numel,': # Dam Elements'
 write(25,'(i8,10x,a/)')nmat, ': # Dam materials'
 write(25,'(3(g14.5,5x))')(ee,enu,rho_d,mat=1,nmat)
 write(25,*)
 do i=1,numel
 write(25,'(2i8)')elements(i)%ng,elements(i)%mat
 write(25,'(8i8:8i8)')(elements(i)%node(j),j=1,2*nod)
 enddo
!
! print dam structure elements
!
 write(25,'(/i2)')ixyz	                                                       
 write(25,'(f5.2)')0.05  ! hyst_damp
 write(25,'(3(f5.2,2x))')0.1,0.2,0.1
 write(25,'(3(f5.2,2x))')(/-1.0,-1.0,-1.0/)
 write(25,'(f8.3,2x,4(i2,2x),f5.2)')corners(8)%radian(3),ldam+1,ndam+1,3,3,0.85	                                                       
!
! print nodes for graphs
!
 write(25,'(/i3)')ldam
 write(25,'(<ldam>(i4,2x))')nodes(1:ldam,ndam/2+1,mdam+1)%nnum
 write(25,'(/i3)')ndam+1
 write(25,'(<ndam+1>(i4,2x))')nodes(ldam,1:ndam+1,mdam+1)%nnum
 
 write(25,'(/i3)')numnp/2
 do i=1,ldam
 nstp=1
 if(mod(i,2)==0)nstp=2
 write(25,'(<ndam>(i4,2x))')nodes(i,1:ndam:nstp,mdam+1)%nnum
 enddo
!
 write(25,'(/i3)')ldam
 write(25,'(<ldam>(i4,2x))')nodes(1:ldam,ndam/2+1,mdam+1)%nnum

 close(25)
 endif
!
! print SAP input
!
 if(sap==1.and.nod==4)then
 open(15,file='Sap',status='unknown')
!
 write(15,'(a)') title
 write(15,'(a)') 'SYSTEM'
 write(15,'(a,i2)') 'V=',nf
!
 write(15,'(a)') 'JOINTS'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 write(15,'(i4,3x,"X=",f9.3,3x,"Y=",f9.3,3x,"Z=",f9.3)')&
      &nnum,nodes(i,j,k)%coord
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo
!
!******** reservoir nodes
 if(sap_res==1)then
 do nnum=numnp+1,ntot
 indicator=' '
 do i=1,nn1
 do j=1,nn2
 if(indicator==' '.and.node(i,j)%nnum==nnum)then
 write(15,'(i4,3x,"X=",f9.3,3x,"Y=",f9.3,3x,"Z=",f9.3)')&
      &nnum,node(i,j)%coord
 indicator='*'
 endif
 enddo ;enddo ;enddo
 endif
!
!*******
 write(15,'(/a)') 'RESTRAINTS'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 write(15,'(3(i4,2x),2x,a,1x,5(i1,","),i1)')&
 &nnum,nnum,1,'R=',nodes(i,j,k)%id(1:3),1,1,1
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo

!*********** reservoir nodes
!
 if(sap_res==1)then
 do nnum=numnp+1,ntot
 indicator=' '
 do i=1,nn1
 do j=1,nn2
 if(indicator==' '.and.node(i,j)%nnum==nnum)then
 write(15,'(3(i4,2x),2x,a,1x,5(i1,","),i1)')&
 &nnum,nnum,1,'R=',node(i,j)%id(1:3),1,1,1
 indicator='*'
 endif
 enddo ;enddo ;enddo
 endif
!
!***********

 write(15,'(/a)') 'SOLID'
 write(15,'("NM= ",i2,4x,"maxn=",i2)')nmat,1
 write(15,'(i2,1x,"Numt=",i2,2x,"W=",f9.3,2x,"M=",f9.3)')&
 &(mat,1,rho_d*gravity,rho_d,mat=1,nmat)
 write(15,'("T=",f5.2,2x,"E=",G12.5,2x,"U=",f8.4)')&
 &(0.,ee,enu,mat=1,nmat)
 llmm(1:8)=(/2,1,3,4,6,5,7,8/)
 do i=1,numel
 write(15,'(i4,4x,"JQ=",7(i4,","),i4,4x,"M=",i2,3x,"I=",i2)')&
 &i,(elements(i)%node(llmm(j)),j=1,2*nod),elements(i)%mat,0
 enddo

!*********** reservoir nodes
!
 if(sap_res==1)then
 write(15,'(/a)') 'SHELL'
 write(15,'("NM= ",i2)')nmat
 write(15,'(i2,2x,"E=",G12.5,2x,"U=",f8.4,1x,"W=",f9.3,2x,"M=",f9.3)')&
 &(mat,ee,enu,rho_d*gravity,rho_d,mat=1,nmat)
 llmm(1:4)=(/2,1,3,4/)
 do i=1,nelm
 write(15,'(i4,4x,"JQ=",3(i4,","),i4,4x,"M=",i2,3x,"TH=",2(f5.2,2x))')&
 &i,(element(i)%node(llmm(j)),j=1,nod),1,.1,.1
 enddo
 endif
!
 endif
!
!***********

 close(15)
!
! print LUSAS input
!
 if(Lusas==1)then
 open(35,file='Lusas.dat',status='unknown')
!
 write(35,'(a,a)')'problem title ',title
! write(35,'(a)')'units n m kg'
! write(35,'(a)')' option  '
!
! Elements properties
! 
 if(nod==4)then
 write(35,'(/a)')'hx8 element topology'
 llmm(1:8)=(/1,2,3,4,5,6,7,8/)
 do i=1,numel
 write(35,'(i4,4x,8(i4,2x))')i,(elements(i)%node(llmm(j)),j=1,2*nod)
 enddo

 elseif(nod==8)then
 write(35,'(/a)')'hx16 element topology'
 llmm(1:16)=(/1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16/)
 do i=1,numel
 write(35,'(i4,4x,16(i4,2x))')i,(elements(i)%node(llmm(j)),j=1,2*nod)
 enddo

 endif
! 
! nodes properties
!
 write(35,'(/a)')'node coordinates'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 write(35,'(i4,3x,3(f9.3,2x))')nnum,nodes(i,j,k)%coord
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo
!
! material properties
!
 write(35,'(/a)')'material properties'
 write(35,'(a)')'c   E(N/mm2)  pois  massdens(kg)'
 write(35,'(i2,2x,G12.5,2x,f5.2,2x,f8.4)')(mat,ee,enu,rho_d,mat=1,nmat)
 write(35,'(a)')'material assign'
 write(35,'(4(i4,2x))')1,numel,1,1
!
! restraints
!
 write(35,'(/a)')'support nodes'
 write(35,'(a)')'c dof U V W Rx Ry Rz'
 do nnum=1,numnp
 indicator=' '
 do i=1,ldam+1
 do j=1,ndam+1
 do k=1,mdam+1
 if(indicator==' '.and.nodes(i,j,k)%nnum==nnum)then
 rff(1:3)='r'
 if(nodes(i,j,k)%id(1)==0)rff(1)='f'
 if(nodes(i,j,k)%id(2)==0)rff(2)='f'
 if(nodes(i,j,k)%id(3)==0)rff(3)='f'
 write(35,'(3(i4,2x),2x,3(a1,2x))')nnum,nnum,0,rff(1:3)
 indicator='*'
 endif
 enddo ; enddo ; enddo; enddo
!
! Eigen value
!
 write(35,'(/a)')'eigenvalue control'
 write(35,'(a,i3)')'constants ',nf
 write(35,'(a)')'convegence 1.e-5  30'
!
 write(35,'(/a)')'end'
!
 endif
!
 close(35)


end subroutine Output   
