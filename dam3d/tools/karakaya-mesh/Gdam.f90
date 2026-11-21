subroutine dam
 use mesh
!
! Allocating data & working arrays  
!  
 numel=ndam*ldam*mdam
 allocate(nodes(ldam+1,ndam+1,mdam+1),elements(numel))
!
! Generation of KARAKAYA Dam structure nodes
!
 inc=0
!
 start_i(1)=1       ;finish_i(1)=nl+1     ;start_j(1)=1       ;finish_j(1)=lb+1 
 start_i(2)=nl+1    ;finish_i(2)=nl+nc+1  ;start_j(2)=1       ;finish_j(2)=lb+1 
 start_i(3)=nl+nc+1 ;finish_i(3)=ndam+1   ;start_j(3)=1       ;finish_j(3)=lb+1 
 start_i(4)=1       ;finish_i(4)=nl+1     ;start_j(4)=lb+1    ;finish_j(4)=lb+lc+1 
 start_i(5)=nl+1    ;finish_i(5)=nl+nc+1  ;start_j(5)=lb+1    ;finish_j(5)=lb+lc+1  
 start_i(6)=nl+nc+1 ;finish_i(6)=ndam+1   ;start_j(6)=lb+1    ;finish_j(6)=lb+lc+1  
 start_i(7)=1       ;finish_i(7)=nl+1     ;start_j(7)=lb+lc+1 ;finish_j(7)=ldam+1 
 start_i(8)=nl+nc+1 ;finish_i(8)=ndam+1   ;start_j(8)=lb+lc+1 ;finish_j(8)=ldam+1
!
 inc_i(1)=nl        ;inc_j(1)=lb          ;inc_i(2)=nc        ;inc_j(2)=lb     
 inc_i(3)=nr        ;inc_j(3)=lb	      ;inc_i(4)=nl        ;inc_j(4)=lc     
 inc_i(5)=nc        ;inc_j(5)=lc          ;inc_i(6)=nr        ;inc_j(6)=lc	      
 inc_i(7)=nl        ;inc_j(7)=lu	      ;inc_i(8)=nr        ;inc_j(8)=lu	      
!
 pts(1,1,1)=1     ;pts(1,2,1)=2     ;pts(1,3,1)=6     ;pts(1,4,1)=5     
 pts(1,1,2)=17    ;pts(1,2,2)=18    ;pts(1,3,2)=22    ;pts(1,4,2)=21    
 pts(2,1,1)=2     ;pts(2,2,1)=3     ;pts(2,3,1)=7     ;pts(2,4,1)=6     
 pts(2,1,2)=18    ;pts(2,2,2)=19    ;pts(2,3,2)=23    ;pts(2,4,2)=22    
 pts(3,1,1)=3     ;pts(3,2,1)=4     ;pts(3,3,1)=8     ;pts(3,4,1)=7     
 pts(3,1,2)=19    ;pts(3,2,2)=20    ;pts(3,3,2)=24    ;pts(3,4,2)=23    
 pts(4,1,1)=5     ;pts(4,2,1)=6     ;pts(4,3,1)=10    ;pts(4,4,1)=9    
 pts(4,1,2)=21    ;pts(4,2,2)=22    ;pts(4,3,2)=26    ;pts(4,4,2)=25    
 pts(5,1,1)=6     ;pts(5,2,1)=7     ;pts(5,3,1)=11    ;pts(5,4,1)=10    
 pts(5,1,2)=22    ;pts(5,2,2)=23    ;pts(5,3,2)=27    ;pts(5,4,2)=26    
 pts(6,1,1)=7     ;pts(6,2,1)=8     ;pts(6,3,1)=12    ;pts(6,4,1)=11    
 pts(6,1,2)=23    ;pts(6,2,2)=24    ;pts(6,3,2)=28    ;pts(6,4,2)=27    
 pts(7,1,1)=9     ;pts(7,2,1)=10    ;pts(7,3,1)=14    ;pts(7,4,1)=13    
 pts(7,1,2)=25    ;pts(7,2,2)=26    ;pts(7,3,2)=30    ;pts(7,4,2)=29    
 pts(8,1,1)=11    ;pts(8,2,1)=12    ;pts(8,3,1)=16    ;pts(8,4,1)=15    
 pts(8,1,2)=27    ;pts(8,2,2)=28    ;pts(8,3,2)=32    ;pts(8,4,2)=31    
!
 do ij=1,8
 do k=1,mb+1
 do j=start_j(ij),finish_j(ij)
 do i=start_i(ij),finish_i(ij)
 pts1=pts(ij,1,1) ; pts2=pts(ij,2,1)
 pts3=pts(ij,3,1) ; pts4=pts(ij,4,1)
 pts5=pts(ij,1,2) ; pts6=pts(ij,2,2)
 pts7=pts(ij,3,2) ; pts8=pts(ij,4,2)
 ii=i-start_i(ij)
 jj=j-start_j(ij)
!
 x10=corners(pts1)%radian ;x20=corners(pts2)%radian
 x30=corners(pts3)%radian ;x40=corners(pts4)%radian
 x11=corners(pts5)%radian ;x21=corners(pts6)%radian
 x31=corners(pts7)%radian ;x41=corners(pts8)%radian
!
 x1=x10+(x11-x10)/mdam*(k-1)
 x2=x20+(x21-x20)/mdam*(k-1)
 x3=x30+(x31-x30)/mdam*(k-1)
 x4=x40+(x41-x40)/mdam*(k-1)
!
 xx1=x1+(x2-x1)/inc_i(ij)*ii; xx2=x4+(x3-x4)/inc_i(ij)*ii
 xx=xx1+(xx2-xx1)/inc_j(ij)*jj
 radius=xx(1)
 angle =xx(2)
 height=xx(3)
 nodes(j,i,k)%angle=angle
 nodes(j,i,k)%height=height
 nodes(j,i,k)%coord(1)=radius*cos(angle)
 nodes(j,i,k)%coord(2)=radius*sin(angle)
 nodes(j,i,k)%coord(3)=height
 nodes(j,i,k)%id=0
 nodes(j,i,k)%coord=matmul(T,nodes(j,i,k)%coord)
 nodes(j,i,k)%indicator='*'
 if(nod==8.and.mod(i,2)==0.and.mod(j,2)==0)nodes(j,i,k)%indicator=' '
 enddo;enddo;enddo;
 enddo
!
! numbering of the nodes
!
 nnum=0
 do j=1,ldam+1
 do k=1,mdam+1
 do i=1,ndam+1
 if(nodes(j,i,k)%indicator=='*')then
 nnum=nnum+1
 nodes(j,i,k)%nnum=nnum
 endif
 enddo;enddo;enddo
 numnp=nnum
!
! Define Restraints
!
 do i=1,ndam+1
 do j=1,ldam +1
 do k=1,mdam +1
 nodes(1,i,k)%id(1:3)=(/1,1,1/)
 nodes(j,1,k)%id(1:3)=(/1,1,1/)
 nodes(j,ndam+1,k)%id(1:3)=(/1,1,1/)
 enddo; enddo; enddo
!
! Generation Dam structure solid and thick shell elements
!
 numel=0
 select case(nod)
 case(4)
 do i=1,ldam
 do j=1,ndam
 do k=1,mdam
 test=nodes(i,j,k)%nnum*nodes(i,j+1,k)%nnum*nodes(i+1,j+1,k)%nnum*nodes(i+1,j,k)%nnum&
     &*nodes(i,j,k+1)%nnum*nodes(i,j+1,k+1)%nnum*nodes(i+1,j+1,k+1)%nnum*nodes(i+1,j,k+1)%nnum
 if(test/=0)then
 numel=numel+1
 elements(numel)%node( 1)=nodes(i  ,j  ,k)%nnum
 elements(numel)%node( 2)=nodes(i  ,j+1,k)%nnum
 elements(numel)%node( 3)=nodes(i+1,j+1,k)%nnum
 elements(numel)%node( 4)=nodes(i+1,j  ,k)%nnum
 elements(numel)%node( 5)=nodes(i  ,j  ,k+1)%nnum
 elements(numel)%node( 6)=nodes(i  ,j+1,k+1)%nnum
 elements(numel)%node( 7)=nodes(i+1,j+1,k+1)%nnum
 elements(numel)%node( 8)=nodes(i+1,j  ,k+1)%nnum
 elements(numel)%mat=nmat
 elements(numel)%ng =ng
 endif
 enddo; enddo; enddo

 case(8)
 do i=1,ldam,2
 do j=1,ndam,2
 do k=1,mdam
 test=nodes(i,j,k)%nnum*nodes(i,j+2,k)%nnum*nodes(i+2,j+2,k)%nnum*nodes(i+2,j,k)%nnum&
  &*nodes(i,j+1,k)%nnum*nodes(i+1,j+2,k)%nnum*nodes(i+2,j+1,k)%nnum*nodes(i+1,j,k)%nnum&
  &*nodes(i,j,k+1)%nnum*nodes(i,j+2,k+1)%nnum*nodes(i+2,j+2,k+1)%nnum*nodes(i+2,j,k+1)%nnum&
  &*nodes(i,j+1,k+1)%nnum*nodes(i+1,j+2,k+1)%nnum*nodes(i+2,j+1,k+1)%nnum*nodes(i+1,j,k+1)%nnum
 if(test/=0)then
 numel=numel+1
 elements(numel)%node( 1)=nodes(i  ,j  ,k)%nnum
 elements(numel)%node( 2)=nodes(i  ,j+2,k)%nnum
 elements(numel)%node( 3)=nodes(i+2,j+2,k)%nnum
 elements(numel)%node( 4)=nodes(i+2,j  ,k)%nnum
 elements(numel)%node( 5)=nodes(i  ,j  ,k+1)%nnum
 elements(numel)%node( 6)=nodes(i  ,j+2,k+1)%nnum
 elements(numel)%node( 7)=nodes(i+2,j+2,k+1)%nnum
 elements(numel)%node( 8)=nodes(i+2,j  ,k+1)%nnum
 elements(numel)%node( 9)=nodes(i  ,j+1,k)%nnum
 elements(numel)%node(10)=nodes(i+1,j+2,k)%nnum
 elements(numel)%node(11)=nodes(i+2,j+1,k)%nnum
 elements(numel)%node(12)=nodes(i+1,j  ,k)%nnum
 elements(numel)%node(13)=nodes(i  ,j+1,k+1)%nnum
 elements(numel)%node(14)=nodes(i+1,j+2,k+1)%nnum
 elements(numel)%node(15)=nodes(i+2,j+1,k+1)%nnum
 elements(numel)%node(16)=nodes(i+1,j  ,k+1)%nnum
 elements(numel)%mat=nmat
 elements(numel)%ng =ng
 endif
 enddo; enddo; enddo
!
 end select 
!
! print test data for elements
!
 open(55,file='tests.dat',status='unknown')
 do j=1,ldam+1
 do i=1,ndam+1
 write(55,'(3i8)')j,i,nodes(j,i,mdam+1)%nnum
 enddo;enddo
 close(55)


end subroutine dam