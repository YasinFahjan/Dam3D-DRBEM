subroutine GHOutput
use Bem3D
use graphics
use modal
implicit none
!                                                                               
! BEM3D: hydrodynamic Pressure Output for the graph files 
!
 node%p=node%p*gravity/rho_w/hwater
 if(ixyz==1)str1='U'
 if(ixyz==2)str1='C'
 if(ixyz==3)str1='V'
 if(ndyn==2)str4='I'
 if(ndyn/=2)str4='N'
!
 write(str2,'(i2.2)')int(alfa_r*10)
 write(str3,'(i2.2)')num
! 
 if(garch==1)then
!
 if(gpy==1.and.num/=0)then
 file=trim(path_py)//str1//'poy'//str3//'.dat'
 open(prn3,file=file)
 do i=1,size(GHY,1)
 do j=1,size(GHY,2)
 jj=gn(j)
 if(jj/=0)then
 ij=GHY(i,jj)
 write(prn3,'(2G12.4\)')node(ij)%coord(3)/hwater,abs(node(ij)%p)
 endif
 enddo
 write(prn3,'()')
 enddo
 close(prn3)
 endif
!
 if(gpt==1.and.num/=0)then
 file=trim(path_px)//str1//'pot'//str3//'.dat'
 open(prn4,file=file)
 do j=1,size(GHY,2)
 ij=GHY(1,j)
 gtheta=pi/4-(j-1)*pi/(2.0*(size(GHY,2)-1))
 xtheta=(gradius+abs(node(ij)%p))*cos(gtheta+pi/4.)-gradius
 ytheta=(gradius+abs(node(ij)%p))*sin(gtheta+pi/4.)-gradius
 write(prn4,'(2g12.4)')xtheta,ytheta
 enddo
 close(prn4)
 endif
!
 if(gfw==1)then
! Forces on the upstream face of the dam
 force=force/cdabs(cdsqrt(sum(force_static*force_static)))
 write(prn5,'(4G12.4)')w/wo,cdsqrt(sum(force*force)),cdabs(cdsqrt(sum(force*force)))
 endif
!
 if(gpc==1.and.num/=0)then
 file=trim(path_pc)//str1//'poc'//str3//'.dat'
 open(prn6,file=file)
 do i=1,size(GHY,1)
 do j=1,size(GHY,2)
 ij=GHY(i,j)
 write(prn6,'(3G12.4)')node(ij)%coord(2),&
    &node(ij)%coord(3),abs(node(ij)%p)
 enddo
 enddo
 close(prn6)
 endif
!
 else
!
 if(gpy==1.and.num/=0)then
 file=trim(path_py)//str1//'py'//str2//'w'//str3//'.dat'
 open(prn3,file=file)
 write(prn3,'(4x,4(a,10x))')'y','p[real]','p[complex]','p[abs]'
 do i=1,size(GHY,1)
 do j=1,size(GHY,2)
 jj=gn(j)
 if(jj/=0)then
 ij=GHY(i,jj)
 write(prn3,'(4G12.4\)')node(ij)%coord(3)/hwater,node(ij)%p,abs(node(ij)%p)
 endif
 enddo
 write(prn3,'()')
 enddo
 close(prn3)
 endif
!

 if(gpx==1.and.num/=0)then
 file=trim(path_px)//str1//'px'//str2//'w'//str3//'.dat'
 open(prn4,file=file)
 write(prn4,'(4x,4(a,10x)/)')'x','p[real]','p[complex]','p[abs]'
 do i=1,size(GHX,1)
 do j=1,size(GHX,2)
 jj=gn(j)
 if(jj/=0)then
 ij=GHX(i,j)
 write(prn4,'(4G12.4\)')node(ij)%coord(1)/hwater,node(ij)%p,abs(node(ij)%p)
 endif
 enddo
 write(prn4,'()')
 enddo
 close(prn4)
 endif
!
 if(gfw==1)then
! Forces on the upstream face of the dam
 force=force/cdabs(cdsqrt(sum(force_static*force_static)))
 write(prn5,'(4G12.4)')w/wo,cdsqrt(sum(force*force)),cdabs(cdsqrt(sum(force*force)))
 endif
!
! pressure contoure on the upstream face of the dam
!
 if(gpc==1.and.num/=0)then
 file=trim(path_pc)//str1//'pc'//str2//'w'//str3//'.dat'
 open(prn6,file=file)
 do i=1,size(GHY,1)
 do j=1,size(GHY,2)
 ij=GHY(i,j)
 write(prn6,'(3G12.4)')node(ij)%coord(2),&
    &node(ij)%coord(3),abs(node(ij)%p)
 enddo
 enddo
 close(prn6)
 endif
!
!
 endif
!
end subroutine GHOutput
 
