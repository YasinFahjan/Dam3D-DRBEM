subroutine input
 use mesh
! 
! read names and open  files for input and output    
! 
 write(*,' (/28x,a/26x,a) ')'-KARAKAYA Dam-',' Mesh Generation '
 open(prn,file='kara.in',status='unknown')
 open(prn2,file='node.dat',status='unknown')
 open(prn3,file='element.dat',status='unknown')
 open(prn4,file='nodes.dat',status='unknown')
 open(prn5,file='elements.dat',status='unknown')
!
 title=' Karakaya	Dam Modeling FEM-BEM scheme'
 write(prn,'(a/)') title
! 
! Coordinates of Karakaya dam structure nodes
! 
 corners( 1)%radian(1:3)=(/175.,-0.55476,0./)
 corners( 2)%radian(1:3)=(/175.,-0.28,0./)
 corners( 3)%radian(1:3)=(/175., 0.28,0./)
 corners( 4)%radian(1:3)=(/175., 0.541666,0./)
 corners( 5)%radian(1:3)=(/175.,-0.80714,85./)
 corners( 6)%radian(1:3)=(/175.,-0.40,85./)
 corners( 7)%radian(1:3)=(/175., 0.40,85./)
 corners( 8)%radian(1:3)=(/175., 0.790476,85./)
 corners( 9)%radian(1:3)=(/215.,-1.047619,155./)
 corners(10)%radian(1:3)=(/215.,-0.40,155./)
 corners(11)%radian(1:3)=(/215., 0.40,155./)
 corners(12)%radian(1:3)=(/215., .988095,155./)
 corners(13)%radian(1:3)=(/215.,-1.1,175./)
 corners(14)%radian(1:3)=(/215.,-0.40,175./)
 corners(15)%radian(1:3)=(/215., 0.40,175./)
 corners(16)%radian(1:3)=(/215., 1.042857,175./)
!
 corners(17)%radian(1:3)=(/225.,-0.55476,0./)   
 corners(18)%radian(1:3)=(/225.,-0.28,0./)      
 corners(19)%radian(1:3)=(/225., 0.28,0./)      
 corners(20)%radian(1:3)=(/225., 0.541666,0./)  
 corners(21)%radian(1:3)=(/225.,-0.80714,85./)  
 corners(22)%radian(1:3)=(/225.,-0.40,85./)     
 corners(23)%radian(1:3)=(/225., 0.40,85./)     
 corners(24)%radian(1:3)=(/225., 0.790476,85./) 
 corners(25)%radian(1:3)=(/225.,-1.047619,155./)
 corners(26)%radian(1:3)=(/225.,-0.40,155./)    
 corners(27)%radian(1:3)=(/225., 0.40,155./)    
 corners(28)%radian(1:3)=(/225., .988095,155./) 
 corners(29)%radian(1:3)=(/225.,-1.1,175./)     
 corners(30)%radian(1:3)=(/225.,-0.40,175./)    
 corners(31)%radian(1:3)=(/225., 0.40,175./)    
 corners(32)%radian(1:3)=(/225., 1.042857,175./)
! 				    
! Coordinates of Karakaya fluid domain boundary nodes
! 
 corner(1)%radian(1:3)=(/225.,-0.52,0./)
 corner(2)%radian(1:3)=(/225., 0.52,0./)
 corner(3)%radian(1:3)=(/225., 1.00,155./)
 corner(4)%radian(1:3)=(/225.,-0.96,155./)
 corner(5)%radian(1:3)=(/450.,-0.52,0./)
 corner(6)%radian(1:3)=(/450., 0.52,0./)
 corner(7)%radian(1:3)=(/450., 1.00,155./)
 corner(8)%radian(1:3)=(/450.,-0.96,155./)
! corner(5)%radian(1:3)=(/425.,-0.345,0./)
! corner(6)%radian(1:3)=(/400., 0.0,0./)
! corner(7)%radian(1:3)=(/460., .514,155./)
! corner(8)%radian(1:3)=(/475.,-0.57,155./)
! 
! Resevoir boundary conditions 
!
 code(1)=2; value(1)=0.0	! Dam-Reservoir Interaction
 code(2)=1; value(2)=0.0	! Free surface
 code(3)=4; value(3)=0.0	! Radiation
 code(4)=3; value(4)=0.0	! Absoption Bottom and sides
 code(5)=3; value(5)=0.0	! Absoption Bottom and sides
 code(6)=3; value(6)=0.0	! Absoption Bottom and sides
!
! Reservoir Internal points limits 
!
 xi_int(1:3)=(/250.,-0.30, 20./)   ! r,theta,z  starting corner
 xf_int(1:3)=(/425., 0.35,125./)   ! r,theta,z  end corner
 dx_int(1:3)=(/ 550., 2.15, 125./)   !	r,theta,z  step
!
! number of elements in the model
!
 nl=4		     ! # dam elmement at left part in radial direction 
 nc=6		     ! # dam elmement at central part in radial direction
 nr=4		     ! # dam elmement at right part in radial direction
 lb=4		     ! # dam elmement at bottom part in vertical direction
 lc=4		     ! # dam elmement at central part in vertical direction
 lu=1		     ! # dam elmement at upper part in vertical direction
 mb=1		     ! # dam elmement at thickness of the dam
 ny=10 	     ! # reservoir elements in upstream direction

! nl=2		     ! # dam elmement at left part in radial direction 
! nc=2		     ! # dam elmement at central part in radial direction
! nr=2		     ! # dam elmement at right part in radial direction
! lb=2		     ! # dam elmement at bottom part in vertical direction
! lc=2		     ! # dam elmement at central part in vertical direction
! lu=1		     ! # dam elmement at upper part in vertical direction
! mb=1		     ! # dam elmement at thickness of the dam
! ny=4 	     ! # reservoir elements in upstream direction

 nx=nl+nc+nr     ! # reservoir elements Horizontal Direction across Dam
 lr=0		     ! # reservoir water level under the spillway
 nz=lb+lc-lr     ! # reservoir elements Vertical Direction across Dam
!
! Input Parameter for the model
!
 ag(1:3)=(/1.,0.,0./)   ! ground motion acceleration
 ee=3.5E9*9.806         ! modulus of elsticitiy
 enu=.2       	    !	poisson ratio
 rho_d=2.45e3         ! concrete density
 gravity=9.806        ! gravity
 rho_w  =1.0e3*9.806  ! water density
 sound  =1438.656     ! sound speed in water
 hwater =155.0  
 alfa_r =1.0          ! bottom reflection coeffecient
 nf     =50       	    ! number of vibration modes
 ndyn   =4       	    !	problem type
 ixyz   =1       	    !	ground water direction
 nmat=1               ! Material types
 ng=3       	    ! order of integration
 radiation=1          ! radiation condition          !    
 geometry=0           ! real reservoir geometry:1, symmetry reservir:0    
 kind='solid8r'       ! elment type: solid8- solid8r-3dshell-3dshellr
!
! Element Types
!
 select case(kind)
 case('solid8')
 nod=4; ndr=2; nrot=0     	 
 case('solid8r')
 nod=4; ndr=2; nrot=3     	 
 case('3dshell')
 nod=8; ndr=3; nrot=0     	 
 case('3dshellr')
 nod=8; ndr=3; nrot=5     	 
 end select
!
! Output Files
!
 sap=1        ! 1,0   Input for Sap
 porter=0     ! 1,0   Input for Simple Dam  
 Lusas=1      ! 1,0   Input for Lusas 
 sap_res=0	  ! 1,0   check the geometry of reservoir using SAP
 Gmode  =1    ! 1,0   node data for graphs    

!
! Modification for quadratic elements
!
 select case(nod)
 case(4)
 case(8)
 nl=2*nl; nc=2*nc; nr=2*nr
 lb=2*lb; lc=2*lc; lu=2*lu
 nx=2*nx; ny=2*ny; nz=2*nz
 mb=mb
 case default
 end select 
!
 ndam=nl+nc+nr
 ldam=lb+lc+lu
 mdam=mb
end subroutine input

!(/175.,-0.52,0./)  
!(/175.,-0.28,0./)  
!(/175., 0.28,0./)  
!(/175., 0.52,0./)  
!(/175.,-0.74,85./) 
!(/175.,-0.40,85./) 
!(/175., 0.40,85./) 
!(/175., 0.78,85./) 
!(/215.,-0.96,155./)
!(/215.,-0.40,155./)
!(/215., 0.40,155./)
!(/215., 1.00,155./)
!(/215.,-1.00,175./)
!(/215.,-0.40,175./)
!(/215., 0.40,175./)
!(/215., 1.04,175./)