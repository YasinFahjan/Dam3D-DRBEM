module Graphics
implicit none
!
! Dam3D: General Variables used in FEM and BEM3D/ models
!
 integer,          parameter :: prn3=33   ! DAM3D/Output File #     
 integer,          parameter :: prn4=44   ! DAM3D/Output File #     
 integer,          parameter :: prn5=55   ! DAM3D/Output File #     
 integer,          parameter :: prn6=66   ! DAM3D/Output File #     
 integer,          parameter :: prn7=77   ! DAM3D/Output File #     
 integer,          parameter :: prn8=88   ! DAM3D/Output File #     
 integer,          parameter :: prn9=99   ! DAM3D/Output File #     
 integer,          parameter :: prn10=100 ! DAM3D/Output File #     
 integer,          parameter :: prn11=111 ! DAM3D/Output File #     
 integer,          parameter :: prn12=112 ! DAM3D/Output File #     
 integer,          parameter :: prn13=113 ! DAM3D/Output File #     
!
 integer,          parameter :: Gmode=1   ! 0,1 Graphics Output Files      
 integer,          parameter :: garch=1   ! 0,1 arch dam represetation      
 integer,          parameter :: gpx=0     ! print P Vs x      
 integer,          parameter :: gpy=1     ! print P Vs y     
 integer,          parameter :: gpt=1     ! print P Vs theta     
 integer,          parameter :: gpm=0     ! print P Vs y (modal shape)     
 integer,          parameter :: gpc=1     ! print P across the upstream     
 integer,          parameter :: gfw=1     ! print F Vs w     
 integer,          parameter :: guw=0     ! print u Vs w     
 integer,          parameter :: gaw=0     ! print a Vs w     
 integer,          parameter :: gyw=0     ! print Y Vs w     
 real,             parameter :: gradius=3.0 
 real                        :: gtheta 
 real                        :: xtheta 
 real                        :: ytheta 
!
 integer, allocatable,  save :: GHX(:,:)              
 integer, allocatable,  save :: GHY(:,:)              
 integer, allocatable,  save :: GSS(:,:)              
 integer, allocatable,  save :: gn(:)              
!
 character(72),         save :: path_py  	       ! DAM3D/file name 
 character(72),         save :: path_px  	       ! DAM3D/file name 
 character(72),         save :: path_fw  	       ! DAM3D/file name 
 character(72),         save :: path_pc  	       ! DAM3D/file name 
 character(72),         save :: path_uw  	       ! DAM3D/file name 
!
 integer                     :: num,iostat              
 integer                     :: gn1,gn2          ! number of nodes    
 character(1),          save :: str1
 character(2),          save :: str2
 character(2),          save :: str3
 character(1),          save :: str4
!
end module Graphics 