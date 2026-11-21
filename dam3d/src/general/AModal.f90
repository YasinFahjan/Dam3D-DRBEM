module Modal
implicit none
!
! DAM3D: Variables used in Modal Analysis
!
 real,allocatable,      save :: Y(:)	        ! Generalized Coordinates 
 real,allocatable,      save :: YY(:)	        ! Generalized Coordinates 2nd Deriv.
 real,allocatable,      save :: r_g(:,:)          ! influence vecors rx,ry,rx 
 real,allocatable,      save :: stiff_mass(:,:)    ! Eigen Values {lamda}
 real,allocatable,      save :: eigen_value(:)    ! Eigen Values {lamda}
 real,allocatable,      save :: eigen_vect(:,:)   ! Eigen vectors [v]
 complex(8),allocatable,save :: S(:,:)            ! Modal S matrix
 complex(8),allocatable,save :: L(:)              ! Modal L Vector
 complex(8),allocatable,save :: L_vect(:)         ! Modal L Vector
 integer                     :: id1,id2,id3,idd
!
! FEM3D: Variables used in Modal Analysis	of Structural response
!
 real,allocatable,      save :: m_mat(:)   	  ! Modal Mass Matrix
 real,allocatable,      save :: k_mat(:)	        ! Modal Stiffness	Matrix
 real,allocatable,      save :: c_mat(:)   	  ! Modal Damping Matrix
 complex,allocatable,   save :: u(:)              ! Nodal Displacement
!
! BEM3D: Variables used in Modal Analysis of Hydrodynamic Response
!
 real, allocatable           :: Bmod(:,:)         ! BEM3D/Cosine direction matrix
 complex, allocatable        :: fi(:,:)           ! BEM3D/Vibration mode
 complex, allocatable        :: fi_r(:,:)         ! BEM3D/Vibration mode
 complex, allocatable        :: Po(:)             ! BEM3D/H.Pressure (rigid dam)
 complex, allocatable        :: Po_static(:)      ! BEM3D/H.Pressure (rigid dam)
 complex, allocatable        :: Po_(:)            ! BEM3D/H.Pressure (rigid dam)
 complex, allocatable        :: Qo(:)      	  ! BEM3D/H.Force    (rigid dam)
 complex, allocatable        :: Pj(:,:)   	  ! BEM3D/H.Pressure (dam V.Modes)
 complex, allocatable        :: Pj_(:,:)   	  ! BEM3D/H.Pressure (dam V.Modes)
 complex, allocatable        :: Qj(:,:)   	  ! BEM3D/H.Force    (dam V.Modes)
 complex(8)                  :: Force(3)   	  ! BEM3D/H.Force    (dam V.Modes)
 complex(8)                  :: Force_static(3)	  ! BEM3D/H.Force    (dam V.Modes)
 complex,allocatable,   save :: uu(:)             ! Nodal Acceleration
!
end module Modal
