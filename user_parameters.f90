module user_parameters
	use mesh
! =======================================================================================
! Geometry variables
! Computational domain extends from START_X to START_X+WIDTH in x direction
! and from START_Y to START_Y+WIDTH in y direction
	integer,parameter              :: Nex=20            ! Number of element in x
	integer,parameter              :: Ney=20            ! Number of element in y
	real,parameter                 :: LENGTH=1.         ! Length of the computational domain
	real,parameter                 :: WIDTH=1.          ! Width  of the computational domain
	real                           :: START_X=0.        ! minimum x coordinate of the domain
	real                           :: START_Y=0.        ! minimum y coordinate of the domain
	real                           :: growth_rate_x=1.  ! refinement if needed - not used in this test case (replaced with 2nd order function in geometry.f90)
	real                           :: growth_rate_y=1.  ! refinement if needed - not used in this test case 
! =======================================================================================
! Time-dependent parameters
	real,parameter                 :: DT=0.001          ! This global variable is the global time
	integer,parameter              :: NT=60000          ! Number of time steps 
	integer                        :: STEP_OUT=10000    ! create every STEP_OUT steps
	! Solver paramters
	real,parameter                 :: Tolerance_u=1E-8  ! solver ttolerance for u field
	real,parameter                 :: Tolerance_v=Tolerance_u
	real,parameter                 :: Tolerance_p=1E-6  ! pressure
	real,parameter                 :: dudt_epsilon=1E-5 ! steady-state indicator
	real,parameter                 :: w_u=0.7           ! Relaxation factor in SOR for u field
	real,parameter                 :: w_v=w_u
	real,parameter                 :: w_p=w_u
	real,parameter                 :: w_p2=w_p
	integer,parameter              :: sweep_u=1         ! Number of sweeps for SIMPLE-LIKE algorithm (typical 1)
	integer,parameter              :: sweep_p=1

! =============================== Problem specification ===============================
	real,parameter                 :: Re=1000.0 ! Reynold number
! ================================ Flags ==============================================
	logical,parameter              :: steady_state=.false.         ! not used here		      
	logical,parameter              :: is_bx=.false.                ! if we have momentum source
	logical,parameter              :: is_by=.false.                ! if we have momentum source
	logical,parameter              :: is_bx_TD=.false.             ! if bx time-dependent
	logical,parameter              :: is_by_TD=.false.             ! if by time-dependent
	logical,parameter              :: is_analytical=.false.        ! if analytical solution is provided
	logical,parameter              :: is_analytical_TD=.false.     ! if analytical solution is time-dependent
	logical,parameter              :: is_dpdx_tolerance=.true.     ! recommended  ->  true
	logical,parameter              :: is_check_steady_state=.true. ! check for steady-state
	integer,parameter              :: nim_formulation=2            ! which formulation of NIM to use (2 is recommended -> more stable)
! ==========================================================================================
	! Restart file specifications
	logical,parameter               :: is_write_restart=.true.     ! if to write restart file or not
	logical,parameter               :: is_restart=.false.          ! start from a restart file
	integer,parameter               :: write_restart_step=10000    ! write restart every write_restart_step steps
	integer,parameter               :: read_restart_step=80000     ! Which step to restart from

end module user_parameters
