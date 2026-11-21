subroutine GInput
use Dam3D
use graphics
implicit none
! 
! Dam3D: input and output Graphics File   
!
 if(Gmode==1)then
 path_py='c:\yasin\phd\graphs\kara\testy\'
 path_px='c:\yasin\phd\graphs\kara\testx\'
 path_fw='c:\yasin\phd\graphs\kara\testf\'
 path_pc='c:\yasin\phd\graphs\kara\testc\'
 path_uw='c:\yasin\phd\graphs\kara\tests\'
! 
! Dam3D: Read Nodes for P Vs Y   
!
 read(inp,*,iostat=iostat)Void
 if(iostat==0)then
 read(inp,*)gn1,gn2
 allocate(GHY(gn1,gn2))
 GHY=0
 read(inp,*)((GHY(i,j),j=1,gn2),i=1,gn1)
 endif
! 
! Dam3D: Read Nodes for P Vs X   
!
 read(inp,*,iostat=iostat)Void
 if(iostat==0)then
 read(inp,*)gn1,gn2
 allocate(GHX(gn1,gn2))
 GHX=0
 read(inp,*)((GHX(i,j),j=1,gn2),i=1,gn1)
 endif
!
! Dam3D: Specify the nodes series to be printed
!
 allocate(gn(gn2))
 gn=0
 gn(1)=gn2/2+1
 gn(2)=gn2/4+1
 gn(3)=1
 if(ixyz==2)then
 gn(1)=1
 gn(2)=gn2/4+1
 gn(3)=gn2/2+1
 endif
! 
! Dam3D: Read Nodes for Mode Shapes  
!
 read(inp,*,iostat=iostat)Void
 if(iostat==0)then
 read(inp,*)gn1,gn2
 allocate(GSS(gn1,gn2))
 GSS=0
 read(inp,*)((GSS(i,j),j=1,gn2),i=1,gn1)
 endif
!
! Dam3D: Output file for F vs	w
!
 if(ixyz==1)str1='U'
 if(ixyz==2)str1='C'
 if(ixyz==3)str1='V'
 if(ndyn==2)str4='I'
 if(ndyn/=2)str4='N'
!
 write(str2,'(i2.2)')int(alfa_r*10)
! 
 if(ndyn>=2)then
 if(gfw==1)then
 file=trim(path_fw)//str1//'fw'//str2//str4//'.dat'
 open(prn5,file=file)
 write(prn5,'(4x,4(a,7x)/)')'w','F[real]','F[complex]','F[abs]'
 endif
 endif

 if(ndyn<=2)then
 if(guw==1)then
 file=trim(path_uw)//str1//'uw'//str2//str4//'.dat'
 open(prn7,file=file)
 write(prn7,'(4x,4(a,10x)/)')'w','u(x)','u(y)','u(z)'
 endif
!
 if(gaw==1)then
 file=trim(path_uw)//str1//'aw'//str2//str4//'.dat'
 open(prn8,file=file)
 write(prn8,'(4x,4(a,10x)/)')'w','a(x)','a(y)','a(z)'
 endif
!
 if(gyw==1)then
 file=trim(path_uw)//str1//'yw'//str2//str4//'.dat'
 open(prn9,file=file)
 write(prn9,'(5x,a,13x,<nf>(a,i2.2,9x)/)')'w',('y',i,i=1,nf)
 endif
 endif
!
!
 endif
!
end subroutine GInput
