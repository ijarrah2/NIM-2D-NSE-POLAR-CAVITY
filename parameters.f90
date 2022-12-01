module parameters
  use user_parameters
  use mesh
! ++++++++++++++++ GAUSS INTEGRATION VARIABLES ++++++++++++++++++++++++++++++++++++++
!  integer,parameter              :: Gauss_int_p=3
!  real, dimension(1:Gauss_int_p) :: Gauss_x_i=(/0.7745966692414834,0.,-0.7745966692414834/)
!  real, dimension(1:Gauss_int_p) :: Gauss_w_i=(/0.5555555555555556,0.8888888888888888, 0.5555555555555556/)

	integer,parameter              :: Gauss_int_p=5
	integer,parameter              :: N_Gauss=Gauss_int_p
	real, dimension(1:Gauss_int_p) :: Gauss_x_i=(/-1./3.*sqrt(5.-2.*sqrt(10./7.)), 1./3.*sqrt(5.-2.*sqrt(10./7.)), 0. , &
												  -1./3.*sqrt(5.+2.*sqrt(10./7.)), 1./3.*sqrt(5.+2.*sqrt(10./7.))/)
	real, dimension(1:Gauss_int_p) :: Gauss_w_i=(/(322.+13.*sqrt(70.))/900.,(322.+13.*sqrt(70.))/900., 128./225. , &
												  (322.-13.*sqrt(70.))/900.,(322.-13.*sqrt(70.))/900. /)
! =======================================================================================
	integer,parameter:: left=1
	integer,parameter:: right=2
	integer,parameter:: bottom=3
	integer,parameter:: top=4

! =======================================================================================

	integer,parameter:: velx1=1
	integer,parameter:: velx2=2
	integer,parameter:: PPE=3


! =======================================================================================
	integer,parameter:: no_slip=0
	integer,parameter:: constant_bc=1
	integer,parameter:: dirichlet=2
	integer,parameter:: zero_gradient=3
	integer,parameter:: neumann=4


! =======================================================================================

	integer,parameter:: analytical_velx1=1
	integer,parameter:: analytical_velx2=2
	integer,parameter:: analytical_p=3
	integer,parameter:: IC_velx1=5
	integer,parameter:: IC_velx2=6
	integer,parameter:: IC_PPE=8
	integer,parameter:: bx_fun=10
	integer,parameter:: by_fun=11

	integer,parameter:: bc_dvdt   = 10
	integer,parameter:: bc_dvdtt  = 11
	integer,parameter:: bc_dvdttt = 12

! =======================================================================================
! Geometry variables
	integer,parameter                             :: NMESH_i=Nex
	integer,parameter                             :: NMESH_j=Ney
	integer,parameter                             :: NMESH=NMESH_i*NMESH_j
	integer,parameter                             :: NEDGE=(NMESH_i+1)*NMESH_j+NMESH_i*(NMESH_j+1)
	integer,parameter                             :: NNODE=(NMESH_i+1)*(NMESH_j+1)
	type(cell),dimension(1:NMESH_i,1:NMESH_j)     :: MESHES_ARRAY
	type(edge),dimension(1:NEDGE)                 :: EDGES_ARRAY
	type(node),dimension(1:NMESH_i+1,1:NMESH_j+1) :: NODES_ARRAY
	integer,parameter                             :: tot_edges_bcs=2*NMESH_i+2*NMESH_j
! =======================================================================================
	real,parameter                             :: v_limit=1E-2   ! Limit where v is considered zero
	real,parameter                             :: f_limit=0.0    ! limit where F is considered zero
! =======================================================================================

! Time-dependent parameters
	real                           :: TI=DT
	integer                        :: STEP=0

	! RMS variables
	real                           :: RMS_time=0.
	integer                        :: RMS_time_ctr=0
	real                           :: RMS_edges_u,RMS_edges_v
    real                           :: RMS_meshes_u,RMS_meshes_v
	integer,parameter              :: max_itr=200  ! max # of iterations

	real,parameter                 :: rho=1.0 ! not used
	real,parameter                 :: nu=1/Re

! =============================== TIMING PARAMETERS ====================================
	real*8              :: timing_start
	real*8              :: timing_initialize_1
	real*8              :: timing_initialize_2
	real*8              :: timing_integration
	real*8              :: timing_integration_step, timing_integration_step_start
	real*8              :: timing_all
	real*8              :: initialization_1_time           ! Reading and preparing meshes
	real*8              :: initialization_2_time           ! Initialize the problem: bc, ic, .....
	real*8              :: integration_time                ! NSE integration step timing
	real*8              :: simulation_time                 ! Overall Simulation time

! =============================== Global Arrays Definitions =============================  
	! Arrays related to meshes
	! Unknowns variables
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_uxy_old, meshes_uxy, meshes_uxyt, & 	! Volume averaged variables
														meshes_vxy_old, meshes_vxy, meshes_vxyt, &
														meshes_u, meshes_v, meshes_p			! Previous time step variables -> converted to curvilinear
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_uxy_old_con,meshes_vxy_old_con
	! Sources in NSE -> RC, RD, RB
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_bx, meshes_by, & 
														meshes_R_C_u, meshes_R_B_u, meshes_R_D_u, &
														meshes_R_C_v, meshes_R_B_v, meshes_R_D_v, & 
														meshes_R_B_p, meshes_R_D_p
										   
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_uT, meshes_vT, meshes_pT 		! Analytical values of vel and p

	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_dila, meshes_dila_old						! residual from the continuity equation

	! Discrete equations related variables
	real,dimension(1:11,1:NMESH_i,1:NMESH_j)         :: AA_vx, AA_vy	! Arrays of the discrete coefficients - v
	real,dimension(1:11,1:NMESH_i,1:NMESH_j)         :: AA_px, AA_py ! Arrays of the discrete coefficients - p
	real,dimension(1:7,1:NMESH_i,1:NMESH_j)          :: a_t						! Array of the time-dependent velocity

	! Curvilinear related variables:
	real,dimension(1:2,1:2,1:NMESH_i,1:NMESH_j)      :: meshes_xi_x, meshes_gij  
	real,dimension(1:2,1:2,1:4,1:NMESH_i,1:NMESH_j)  :: meshes_gij_edges
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_G1,meshes_G2
	real,dimension(1:NMESH_i,1:NMESH_j)              :: meshes_area
	real,dimension(1:NMESH_i,1:NMESH_j)              :: edges_dpdx,edges_dpdx_old,edges_dpdy,edges_dpdy_old


	! Geometry related variables
	integer,dimension(1:4,1:NMESH_i,1:NMESH_j)       :: meshes_neighbors
	integer,dimension(1:4,1:NMESH_i,1:NMESH_j)       :: meshes_edges
	real   ,dimension(1:4,1:NMESH_i,1:NMESH_j)       :: meshes_nodes_x,meshes_nodes_y	! Store position of each node in the mesh
! ==========================================================================================
	! Arrays related to edges

	! Unknowns
	real,dimension(1:NEDGE)             :: edges_u, edges_v, edges_p				! Area averaged variables
	real,dimension(1:NEDGE)             :: edges_uT, edges_vT, edges_pT			! Analytical solution over edges
	real,dimension(1:NEDGE)             :: edges_u_old, edges_v_old, edges_p_old
	real,dimension(1:2,1:NEDGE)         :: edges_normal				! Normal vectors
	real,dimension(1:2,1:NEDGE)         :: edges_tangential
	real,dimension(1:2,1:NEDGE)         :: edges_mid_coordinates	! mid-point coordinates
	real,dimension(1:NEDGE)             :: edges_bx, edges_by
	real,dimension(1:NEDGE)             :: edges_length  

	integer,dimension(1:2,1:NEDGE)      :: edges_shared_meshes


	integer,dimension(1:tot_edges_bcs)     :: BCs_edges, BCs_edges_id   	! All edges that have boundary conditions, id:BC id
	integer,dimension(1:tot_edges_bcs)     :: BCs_edges_position
	integer,dimension(1:2,1:tot_edges_bcs) :: BCs_edges_mesh_num
	integer,dimension(1:4,1:tot_edges_bcs) :: BCs_edges_points
	real,   dimension(1:6,1:tot_edges_bcs) :: BCs_edges_ap
	real,   dimension(1:10,1:tot_edges_bcs):: BCs_edges_a
	real,   dimension(1:tot_edges_bcs)     :: BCs_edges_v_t, BCs_edges_v_tt, BCs_edges_v_ttt	! Derivatives with respect to the tangential on the boundaries

! ==========================================================================================
	real                                 :: Dilatation
	real                                 :: dudt_max				! to check the steady-state solution

	! Method 2 formulation variables
    real,dimension(1:NMESH_i,1:NMESH_j)              :: M_A1,M_A2,M_A3,M_A4,M_C1,M_C2,M_B1,M_B2,M_B3,M_B4,M_D1,M_D2
    real,dimension(1:NMESH_i,1:NMESH_j)              :: M_F11,M_F12,M_F13,M_F14,M_F15
    real,dimension(1:NMESH_i,1:NMESH_j)              :: M_F21,M_F22,M_F23,M_F24,M_F25
    real,dimension(1:NMESH_i,1:NMESH_j)              :: M_F31,M_F32,M_F33,M_F34,M_F35,M_F36,M_F37
	
end module parameters
