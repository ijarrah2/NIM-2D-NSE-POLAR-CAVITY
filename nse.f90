module nse
  use general
  use nim
  use user_parameters
  use parameters
  use mesh
  use geometry
  use bc
  use postprocess
  use restart
  use omp_lib
  implicit none
  contains
! =========================================================
! =========================================================
! ==================== Main Algorithm =====================
! =============== Time Integration Steps ==================
! =========================================================
! =========================================================

  subroutine solve_nse
    implicit none
    integer l
! *********************************************************
! ************    INITIALIZE THE DOMAIN    ****************
! *********************************************************
! Setting the current NIM time
	TI=DT
	STEP=1
    call initialize_nse
! Apply all Dirichlet velocity boundary conditions - if not TD
    call apply_all_dirichlet_bcs

	timing_initialize_2=omp_get_wtime()

! *********************************************************
! ************    START THE TIME INTEGRAL    **************
! *********************************************************

    do l=1,NT
! Update the time-dependent variables : bx,by,bz + analytical solution.
! call evaluate_time_dependent_quantities
! Calculate the velocity NIM constants based on the previous time step.
	if(nim_formulation==2) then
		call calculate_velocity_nim_constants_formulation2
	else
		call calculate_velocity_nim_constants
	endif 
	! Main integtation function
		call integrate_nse
		if(is_write_restart .and. mod(STEP,write_restart_step) == 0) call write_restart
		if(mod(STEP,STEP_OUT) == 0)   call write_data

		if(is_check_steady_state) then
!			call check_steady_state
			if(dudt_max<dudt_epsilon) then
				print*, "+++++++++++++++++++++++++++++++++ Steady state is reached ++++++++++++++++++++++++++++++"
				print*, "Simulation is stopped with dudt = ", dudt_max
				if(is_write_restart) call write_restart
				call write_data
				exit
			endif
		endif


! Fix time step and copy the arrays
	   call update_time_step

    enddo

! *********************************************************
! ************    END THE TIME INTEGRAL    ****************
! *********************************************************

  if(is_analytical) call find_rms
  timing_integration=omp_get_wtime()

  end subroutine solve_nse


! =========================================================
  subroutine initialize_nse
  
! *********************************************************
! ************    INITIALIZE THE DOMAIN    ****************
! *********************************************************
! evaluate all constants related to mesh: normals, areas, volume, metrics ... etc
! apply initial condition, physical sources (if steady-state)
    call evaluate_all_mesh_constants

	call evaluate_pressure_bc_array

! Evaluate all PPE constans
	call calculate_pressure_nim_constants
	call evaluate_velocity_derivatives

	if (is_restart) then
		call read_restart
		call write_data
		call update_time_step
	endif
! Write the initial condition
    if (.not. is_restart) call write_data
  end subroutine initialize_nse
! =========================================================

subroutine integrate_nse
  implicit none
  integer                              :: i,j,ff,cnt
  real                                 :: err_u_max,err_v_max,err_p_max
! cnt is iteration counter
  cnt=0
  err_u_max=1.
  err_v_max=1.
  err_p_max=1.


 do while ((err_u_max>TOLERANCE_u .or. err_v_max>TOLERANCE_v .or. err_p_max>TOLERANCE_p) .and. cnt<max_itr)
   cnt=cnt+1
! =================================== VELOCITY INTEGRATION ================================
   do ff=1,sweep_u
		call update_velocity_sources
		call integrate_velocity_field
   enddo

! =================================== PRESSURE INTEGRATION =================================

   do ff=1,sweep_p
		call update_pressure_sources
		call apply_pressure_bc
		call integrate_pressure_field
   enddo

   call find_error(err_u_max,err_v_max,err_p_max)

  enddo
	if(is_check_steady_state) call check_steady_state
	call evaluate_D
	if(.not. is_check_steady_state) then
		if(.not. is_dpdx_tolerance) write (*,'(A, I5, A, I5, A, 1p,E11.4, A, 1p,E11.4)') 'step #', STEP, '  | # iters. ', cnt, '  | err V=', real(MAX(err_u_max,err_v_max),4), '  | err p=', err_p_max 
		if(is_dpdx_tolerance) write (*,'(A, I5, A, I5, A, 1p,E11.4, A, 1p,E11.4)') 'step #', STEP, '  | # iters. ', cnt, '  | err V=', real(MAX(err_u_max,err_v_max),4), '  | err dpdx=', err_p_max 
	else
		if(.not. is_dpdx_tolerance) write (*,'(A, I5, A, I5, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4)') 'step #', STEP, '  | # iters. ', cnt, '  | err V=', real(MAX(err_u_max,err_v_max),4), '  | err p=', err_p_max & 
			 , '  | dudt=', dudt_max
		if(is_dpdx_tolerance) write (*,'(A, I5, A, I5, A, 1p,E11.4, A, 1p,E11.4,A, 1p,E11.4)') 'step #', STEP, '  | # iters. ', cnt, '  | err V=', real(MAX(err_u_max,err_v_max),4), '  | err dpdx=', err_p_max & 
			, '  | dudt=', dudt_max
	endif
end subroutine integrate_nse



subroutine integrate_velocity_field
	integer i,j
	if(nim_formulation==2) then
		!$OMP     PARALLEL DO  PRIVATE (i,j)
		do j=1,NMESH_j
			do i=1,NMESH_i
				call integrate_velocity_field_formulation2_i(i,j)
			enddo
		enddo
		!$OMP     END PARALLEL DO
	else
		!$OMP     PARALLEL DO  PRIVATE (i,j)
		do j=1,NMESH_j
			do i=1,NMESH_i
				call integrate_velocity_field_i(i,j)
			enddo
		enddo
		!$OMP     END PARALLEL DO
	endif
end subroutine integrate_velocity_field



subroutine update_velocity_sources
	integer i,j
     !$OMP     PARALLEL DO  PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
			call calculate_R_b(i,j)
			call calculate_R_C_new(i,j)
			call calculate_R_D(i,j)
		enddo
	enddo
     !$OMP     END PARALLEL DO
	 
end subroutine update_velocity_sources


subroutine integrate_pressure_field
	integer i,j
     !$OMP     PARALLEL DO  PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
        call integrate_pressure_field_i(i,j)
		enddo
	enddo
     !$OMP     END PARALLEL DO
	 
end subroutine integrate_pressure_field


subroutine update_pressure_sources
	integer i,j
     !$OMP     PARALLEL DO  PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
		call calculate_R_b_p(i,j)
		call calculate_R_D_p(i,j)
		enddo
	enddo
    !$OMP     END PARALLEL DO
	 
end subroutine update_pressure_sources



subroutine find_error(err_u_max,err_v_max,err_p_max)
	real, intent(inout)         :: err_u_max,err_v_max,err_p_max
	integer                     :: i,j
	real                        :: error_max_dpdx, error_max_dpdy,error_max_vxy,error_max_uxy
	err_u_max=0.
	err_v_max=0.
	err_p_max=0.
	error_max_dpdx=0.
	error_max_dpdy=0.
	error_max_uxy=0.
	error_max_vxy=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(max: err_u_max,err_v_max,err_p_max)
	do i=1,NEDGE
		err_u_max=MAX(err_u_max,abs(edges_u(i)-edges_u_old(i)))
		err_v_max=MAX(err_v_max,abs(edges_v(i)-edges_v_old(i)))
		err_p_max=MAX(err_p_max,abs(edges_p(i)-edges_p_old(i)))
		
		edges_u_old(i)=edges_u(i)
		edges_v_old(i)=edges_v(i)
		edges_p_old(i)=edges_p(i)
	enddo
	!$OMP     END PARALLEL DO
	if(is_dpdx_tolerance) then
		!$OMP     PARALLEL DO  PRIVATE (i,j) REDUCTION(max: error_max_dpdx,error_max_dpdy)
		do j=1,NMESH_j
			do i=1,NMESH_i
				edges_dpdx(i,j)=(edges_p(meshes_edges(2,i,j))-edges_p(meshes_edges(1,i,j)))/2.!*meshes_xi_x(1,1,i,j)
				edges_dpdy(i,j)=(edges_p(meshes_edges(4,i,j))-edges_p(meshes_edges(3,i,j)))/2.!*meshes_xi_x(2,2,i,j)
				error_max_dpdx=MAX(error_max_dpdx,abs(edges_dpdx(i,j)-edges_dpdx_old(i,j)))
				error_max_dpdy=MAX(error_max_dpdy,abs(edges_dpdy(i,j)-edges_dpdy_old(i,j)))
				edges_dpdx_old(i,j)=edges_dpdx(i,j)
				edges_dpdy_old(i,j)=edges_dpdy(i,j)
			enddo
		enddo
		!$OMP     END PARALLEL DO
		err_p_max=MAX(error_max_dpdx,error_max_dpdy)
	endif
		!$OMP     PARALLEL DO  PRIVATE (i,j) REDUCTION(max: error_max_uxy,error_max_vxy)
		do j=1,NMESH_j
			do i=1,NMESH_i
				error_max_uxy=MAX(error_max_uxy,abs(meshes_uxy(i,j)-meshes_uxy_old_con(i,j)))
				error_max_vxy=MAX(error_max_vxy,abs(meshes_vxy(i,j)-meshes_vxy_old_con(i,j)))
				meshes_uxy_old_con(i,j)=meshes_uxy(i,j)
				meshes_vxy_old_con(i,j)=meshes_vxy(i,j)
			enddo
		enddo
		!$OMP     END PARALLEL DO
		err_u_max=MAX(err_u_max,error_max_uxy)
		err_v_max=MAX(err_v_max,error_max_vxy)
end subroutine find_error


subroutine check_steady_state
	integer                     :: i,j
	real                        :: max_dudt, max_dvdt
	max_dudt=0.
	max_dvdt=0.
	!$OMP     PARALLEL DO  PRIVATE (i,j) REDUCTION(max: max_dudt,max_dvdt)
	do j=1,NMESH_j
		do i=1,NMESH_i
			max_dudt=MAX(max_dudt,abs(meshes_uxy(i,j)-meshes_uxy_old(i,j))/(2*DT))
			max_dvdt=MAX(max_dvdt,abs(meshes_vxy(i,j)-meshes_vxy_old(i,j))/(2*DT))
		enddo
	enddo
	!$OMP     END PARALLEL DO
	dudt_max=MAX(max_dudt,max_dvdt)
end subroutine check_steady_state

subroutine integrate_velocity_field_i(i,j)
	integer, intent(in)                  :: i,j
	real                                 :: B1,B2,B3,B_ij
	integer                              :: fnn
	integer, dimension(1:4)              :: n1_edges, n2_edges
!	call apply_zero_gradient_bc
! =========================== INTEGRATE U ============================================
	n1_edges(:)=meshes_edges(:,i,j)
	B_ij=a_t(5,i,j)*(meshes_R_C_u(i,j)+meshes_R_D_u(i,j)+meshes_R_B_u(i,j))
	meshes_uxy(i,j)=(1-w_u)* meshes_uxy(i,j)+w_u*( &
	             -a_t(1,i,j)* edges_u(n1_edges(1)) &
	             -a_t(2,i,j)* edges_u(n1_edges(2)) &
				 -a_t(3,i,j)* edges_u(n1_edges(3)) &
				 -a_t(4,i,j)* edges_u(n1_edges(4)) &
				 -a_t(6,i,j)* meshes_uxy_old(i,j) &
				 -B_ij)/(a_t(7,i,j))
! =========================== INTEGRATE V ============================================

	B_ij=a_t(5,i,j)*(meshes_R_C_v(i,j)+meshes_R_D_v(i,j)+meshes_R_B_v(i,j))
	meshes_vxy(i,j)=(1-w_v)* meshes_vxy(i,j)+w_v*( &
	             -a_t(1,i,j)* edges_v(n1_edges(1)) &
	             -a_t(2,i,j)* edges_v(n1_edges(2)) &
				 -a_t(3,i,j)* edges_v(n1_edges(3)) &
				 -a_t(4,i,j)* edges_v(n1_edges(4)) &
				 -a_t(6,i,j)* meshes_vxy_old(i,j) &
				 -B_ij)/(a_t(7,i,j))
	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================

	if (i<NMESH_i) THEN
	   n2_edges(:)=meshes_edges(:,i+1,j)
! =========================== INTEGRATE U ============================================
	   B1=(meshes_R_C_u(i,j)+meshes_R_D_u(i,j)+meshes_R_B_u(i,j))*AA_vx(9,i,j)*2*DT
	   B2=(meshes_R_C_u(i+1,j)+meshes_R_D_u(i+1,j)+meshes_R_B_u(i+1,j))*AA_vx(11,i,j)*2*DT
	   B3= (meshes_gij_edges(1,2,1,i+1,j)-meshes_gij_edges(1,2,2,i,j))*(edges_u(n1_edges(4))-edges_u(n1_edges(3)) + &
										                            edges_u(n2_edges(4))-edges_u(n2_edges(3)))/4.0

	   edges_u(n1_edges(2))=(1-w_u)*edges_u(n1_edges(2)) + w_u/(AA_vx(2,i,j)) * ( &
											-AA_vx(1 ,i,j)*edges_u(n1_edges(1)) &
											-AA_vx(3 ,i,j)*edges_u(n2_edges(2)) &
											-AA_vx(4 ,i,j)*edges_u(n1_edges(3)) &
											-AA_vx(5 ,i,j)*edges_u(n1_edges(4)) &
											-AA_vx(6 ,i,j)*edges_u(n2_edges(3)) &
											-AA_vx(7 ,i,j)*edges_u(n2_edges(4)) &
											-AA_vx(8,i,j)*(meshes_uxy(i,j)-meshes_uxy_old(i,j)) &
											-AA_vx(10,i,j)*(meshes_uxy(i+1,j)-meshes_uxy_old(i+1,j)) &
											+B1+B2+B3)

! =========================== INTEGRATE v ============================================
	   B1=(meshes_R_C_v(i,j)+meshes_R_D_v(i,j)+meshes_R_B_v(i,j))*AA_vx(9,i,j)*2*DT
	   B2=(meshes_R_C_v(i+1,j)+meshes_R_D_v(i+1,j)+meshes_R_B_v(i+1,j))*AA_vx(11,i,j)*2*DT
	   B3= (meshes_gij_edges(1,2,1,i+1,j)-meshes_gij_edges(1,2,2,i,j))*(edges_v(n1_edges(4))-edges_v(n1_edges(3)) + &
										                            edges_v(n2_edges(4))-edges_v(n2_edges(3)))/4.0

	   edges_v(n1_edges(2))=(1-w_v)*edges_v(n1_edges(2)) + w_v/(AA_vx(2,i,j)) * ( &
											-AA_vx(1 ,i,j)*edges_v(n1_edges(1)) &
											-AA_vx(3 ,i,j)*edges_v(n2_edges(2)) &
											-AA_vx(4 ,i,j)*edges_v(n1_edges(3)) &
											-AA_vx(5 ,i,j)*edges_v(n1_edges(4)) &
											-AA_vx(6 ,i,j)*edges_v(n2_edges(3)) &
											-AA_vx(7 ,i,j)*edges_v(n2_edges(4)) &
											-AA_vx(8,i,j)*(meshes_vxy(i,j)-meshes_vxy_old(i,j)) &
											-AA_vx(10,i,j)*(meshes_vxy(i+1,j)-meshes_vxy_old(i+1,j)) &
											+B1+B2+B3)

	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (j<NMESH_j) THEN
	   n2_edges(:)=meshes_edges(:,i,j+1)
! =========================== INTEGRATE U ============================================
	   B1=(meshes_R_C_u(i,j)+meshes_R_D_u(i,j)+meshes_R_B_u(i,j))*AA_vy(9,i,j)*2*DT
	   B2=(meshes_R_C_u(i,j+1)+meshes_R_D_u(i,j+1)+meshes_R_B_u(i,j+1))*AA_vy(11,i,j)*2*DT
	   B3= (meshes_gij_edges(2,1,3,i,j+1)-meshes_gij_edges(2,1,4,i,j))*(edges_u(n1_edges(2))-edges_u(n1_edges(1)) + &
										                            edges_u(n2_edges(2))-edges_u(n2_edges(1)))/4.0

	   edges_u(n1_edges(4))=(1-w_u)*edges_u(n1_edges(4)) + w_u/(AA_vy(6,i,j)) * ( &
											-AA_vy(1 ,i,j)*edges_u(n1_edges(1)) &
											-AA_vy(2 ,i,j)*edges_u(n1_edges(2)) &
											-AA_vy(3 ,i,j)*edges_u(n2_edges(1)) &
											-AA_vy(4 ,i,j)*edges_u(n2_edges(2)) &
											-AA_vy(5 ,i,j)*edges_u(n1_edges(3)) &
											-AA_vy(7 ,i,j)*edges_u(n2_edges(4)) &
											-AA_vy(8,i,j)*(meshes_uxy(i,j)-meshes_uxy_old(i,j)) &
											-AA_vy(10,i,j)*(meshes_uxy(i,j+1)-meshes_uxy_old(i,j+1)) &
											+B1+B2+B3)
											
! =========================== INTEGRATE V ============================================
	   B1=(meshes_R_C_v(i,j)+meshes_R_D_v(i,j)+meshes_R_B_v(i,j))*AA_vy(9,i,j)*2*DT
	   B2=(meshes_R_C_v(i,j+1)+meshes_R_D_v(i,j+1)+meshes_R_B_v(i,j+1))*AA_vy(11,i,j)*2*DT
	   B3= (meshes_gij_edges(2,1,3,i,j+1)-meshes_gij_edges(2,1,4,i,j))*(edges_v(n1_edges(2))-edges_v(n1_edges(1)) + &
										                            edges_v(n2_edges(2))-edges_v(n2_edges(1)))/4.0


	   edges_v(n1_edges(4))=(1-w_v)*edges_v(n1_edges(4)) + w_v/(AA_vy(6,i,j)) * ( &
											-AA_vy(1 ,i,j)*edges_v(n1_edges(1)) &
											-AA_vy(2 ,i,j)*edges_v(n1_edges(2)) &
											-AA_vy(3 ,i,j)*edges_v(n2_edges(1)) &
											-AA_vy(4 ,i,j)*edges_v(n2_edges(2)) &
											-AA_vy(5 ,i,j)*edges_v(n1_edges(3)) &
											-AA_vy(7 ,i,j)*edges_v(n2_edges(4)) &
											-AA_vy(8,i,j)*(meshes_vxy(i,j)-meshes_vxy_old(i,j)) &
											-AA_vy(10,i,j)*(meshes_vxy(i,j+1)-meshes_vxy_old(i,j+1)) &
											+B1+B2+B3)

	endif

end subroutine integrate_velocity_field_i


subroutine integrate_velocity_field_formulation2_i(i,j)
	integer, intent(in)                  :: i,j
	real                                 :: B1,B2,B3,B_ij
	integer                              :: fnn
	integer, dimension(1:4)              :: n1_edges, n2_edges
	n1_edges(:)=meshes_edges(:,i,j)
	B_ij=(meshes_R_C_u(i,j)+meshes_R_D_u(i,j)+meshes_R_B_u(i,j))
	
! =========================== INTEGRATE U ============================================

	meshes_uxy(i,j)=(1-w_u)*meshes_uxy(i,j)+(w_u/M_F31(i,j))*( &
					   +M_F32(i,j)*meshes_uxy_old(i,j) &
					   +M_F33(i,j)*edges_u(n1_edges(1)) &
					   +M_F34(i,j)*edges_u(n1_edges(2)) &
					   +M_F35(i,j)*edges_u(n1_edges(3)) &
					   +M_F36(i,j)*edges_u(n1_edges(4)) &
					   +B_ij)
!	print*, M_F31(i,j)
! =========================== INTEGRATE V ============================================
								 
	B_ij=(meshes_R_C_v(i,j)+meshes_R_D_v(i,j)+meshes_R_B_v(i,j))
	meshes_vxy(i,j)=(1-w_v)*meshes_vxy(i,j)+(w_v/M_F31(i,j))*( &
					   +M_F32(i,j)*meshes_vxy_old(i,j) &
					   +M_F33(i,j)*edges_v(n1_edges(1)) &
					   +M_F34(i,j)*edges_v(n1_edges(2)) &
					   +M_F35(i,j)*edges_v(n1_edges(3)) &
					   +M_F36(i,j)*edges_v(n1_edges(4)) &
					   +B_ij)
	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================

	if (i<NMESH_i) THEN
	     n2_edges(:)=meshes_edges(:,i+1,j)
! =========================== INTEGRATE U ============================================
        B3= (meshes_gij_edges(1,2,1,i+1,j)-meshes_gij_edges(1,2,2,i,j))*(edges_u(n1_edges(4))-edges_u(n1_edges(3)) + &
                                                                         edges_u(n2_edges(4))-edges_u(n2_edges(3)))/4.0

		edges_u(n1_edges(2))=(1-w_u)*edges_u(n1_edges(2))+(w_u/M_F11(i,j)) * ( &
								 +M_F12(i,j) *  edges_u(n1_edges(1)) &
								 +M_F13(i,j) *  edges_u(n2_edges(2)) &
								 +M_F14(i,j) * (meshes_uxy(i,j)+meshes_uxy_old(i,j))     &
								 +M_F15(i,j) * (meshes_uxy(i+1,j)+meshes_uxy_old(i+1,j)) &
								 -B3)

		B3= (meshes_gij_edges(1,2,1,i+1,j)-meshes_gij_edges(1,2,2,i,j))*(edges_v(n1_edges(4))-edges_v(n1_edges(3)) + &
                                                                         edges_v(n2_edges(4))-edges_v(n2_edges(3)))/4.0

		edges_v(n1_edges(2))=(1-w_v)*edges_v(n1_edges(2))+(w_v/M_F11(i,j)) * ( &
								 +M_F12(i,j) *  edges_v(n1_edges(1)) &
								 +M_F13(i,j) *  edges_v(n2_edges(2)) &
								 +M_F14(i,j) * (meshes_vxy(i,j)+meshes_vxy_old(i,j))     &
								 +M_F15(i,j) * (meshes_vxy(i+1,j)+meshes_vxy_old(i+1,j)) &
								 -B3)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (j<NMESH_j) THEN
		n2_edges(:)=meshes_edges(:,i,j+1)

        B3= (meshes_gij_edges(2,1,3,i,j+1)-meshes_gij_edges(2,1,4,i,j))*(edges_u(n1_edges(2))-edges_u(n1_edges(1)) + &
                                                                         edges_u(n2_edges(2))-edges_u(n2_edges(1)))/4.0

		edges_u(n1_edges(4))=(1-w_u)*edges_u(n1_edges(4))+(w_u/M_F21(i,j)) * ( &
								 +M_F22(i,j) * edges_u(n1_edges(3)) &
								 +M_F23(i,j) * edges_u(n2_edges(4)) &
								 +M_F24(i,j) * (meshes_uxy(i,j)+meshes_uxy_old(i,j)) &
								 +M_F25(i,j) * (meshes_uxy(i,j+1)+meshes_uxy_old(i,j+1)) &
								 -B3)

        B3= (meshes_gij_edges(2,1,3,i,j+1)-meshes_gij_edges(2,1,4,i,j))*(edges_v(n1_edges(2))-edges_v(n1_edges(1)) + &
                                                                         edges_v(n2_edges(2))-edges_v(n2_edges(1)))/4.0

		edges_v(n1_edges(4))=(1-w_v)*edges_v(n1_edges(4))+(w_v/M_F21(i,j)) * ( &
								 +M_F22(i,j) * edges_v(n1_edges(3)) &
								 +M_F23(i,j) * edges_v(n2_edges(4)) &
								 +M_F24(i,j) * (meshes_vxy(i,j)+meshes_vxy_old(i,j)) &
								 +M_F25(i,j) * (meshes_vxy(i,j+1)+meshes_vxy_old(i,j+1)) &
								 -B3)

	endif

end subroutine integrate_velocity_field_formulation2_i


subroutine integrate_pressure_field_i(i,j)
	integer, intent(in)                  :: i,j
	real                                 :: B1,B2,B3,B_ij
	integer                              :: fnn,n1,n2
	integer, dimension(1:4)              :: n1_edges, n2_edges

	n1_edges(:)=meshes_edges(:,i,j)

	!                               ======================
	! ============================== ||F||A||C||E|| ||2|| =============================
	!                               ======================
	if (i<NMESH_i) THEN
	   n2_edges(:)=meshes_edges(:,i+1,j)

	   B1=(meshes_R_D_p(i,j)+meshes_R_B_p(i,j))*AA_px(9,i,j)
	   B2=(meshes_R_D_p(i+1,j)+meshes_R_B_p(i+1,j))*AA_px(11,i,j)
	   B3= (meshes_gij_edges(1,2,1,i+1,j)-meshes_gij_edges(1,2,2,i,j))*(edges_p(n1_edges(4))-edges_p(n1_edges(3)) + &
										                            edges_p(n2_edges(4))-edges_p(n2_edges(3)))/4.0

	   edges_p(n1_edges(2))=(1-w_p)*edges_p(n1_edges(2)) + w_p/(AA_px(2,i,j)) * ( &
											-AA_px(1 ,i,j)*edges_p(n1_edges(1)) &
											-AA_px(3 ,i,j)*edges_p(n2_edges(2)) &
											-AA_px(4 ,i,j)*edges_p(n1_edges(3)) &
											-AA_px(5 ,i,j)*edges_p(n1_edges(4)) &
											-AA_px(6 ,i,j)*edges_p(n2_edges(3)) &
											-AA_px(7 ,i,j)*edges_p(n2_edges(4)) &
											+B1+B2+B3)
	endif
	!                               ======================
	! ============================== ||F||A||C||E|| ||4|| =============================
	!                               ======================
	if (j<NMESH_j) THEN
	   n2_edges(:)=meshes_edges(:,i,j+1)

	   B1=(meshes_R_D_p(i,j)+meshes_R_B_p(i,j))*AA_py(9,i,j)
	   B2=(meshes_R_D_p(i,j+1)+meshes_R_B_p(i,j+1))*AA_py(11,i,j)
	   B3= (meshes_gij_edges(2,1,3,i,j+1)-meshes_gij_edges(2,1,4,i,j))*(edges_p(n1_edges(2))-edges_p(n1_edges(1)) + &
										                            edges_p(n2_edges(2))-edges_p(n2_edges(1)))/4.0

	   edges_p(n1_edges(4))=(1-w_p)*edges_p(n1_edges(4)) + w_p/(AA_py(6,i,j)) * ( &
											-AA_py(1 ,i,j)*edges_p(n1_edges(1)) &
											-AA_py(2 ,i,j)*edges_p(n1_edges(2)) &
											-AA_py(3 ,i,j)*edges_p(n2_edges(1)) &
											-AA_py(4 ,i,j)*edges_p(n2_edges(2)) &
											-AA_py(5 ,i,j)*edges_p(n1_edges(3)) &
											-AA_py(7 ,i,j)*edges_p(n2_edges(4)) &
											+B1+B2+B3)
	endif
end subroutine integrate_pressure_field_i


subroutine update_time_step
	integer i,j
! Update time and current time-step
		STEP=STEP+1
		TI=TI+2*DT
! UPDATE TIME STEP
!$OMP     PARALLEL DO PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
			meshes_uxyt(i,j)=(meshes_uxy_old(i,j)+meshes_uxy(i,j))/2.
			meshes_uxy_old(i,j)=meshes_uxy(i,j)
			meshes_vxyt(i,j)=(meshes_vxy_old(i,j)+meshes_vxy(i,j))/2.
			meshes_vxy_old(i,j)=meshes_vxy(i,j)
			meshes_p(i,j)=(edges_p(meshes_edges(1,i,j))+edges_p(meshes_edges(2,i,j))+edges_p(meshes_edges(3,i,j))+edges_p(meshes_edges(4,i,j)))/4.0
			meshes_dila_old(i,j)=meshes_dila(i,j)
		enddo
	enddo
!$OMP     END PARALLEL DO

end subroutine update_time_step


subroutine calculate_R_c(Mi,Mj)
	integer, intent(in)         :: Mi,Mj
	real                        :: du_dxi1,du_dxi2,dv_dxi1,dv_dxi2
	real                        :: duu_dx,duv_dy,duv_dx,dvv_dy
    integer, dimension(1:4)     :: n_edges
	integer                     :: k
    real, dimension(1:4)        :: u_edges,v_edges,uu_edges,uv_edges,vv_edges

	
	n_edges(:)=meshes_edges(:,Mi,Mj)
	do k=1,4
		u_edges(k) =edges_u(n_edges(k))
		v_edges(k) =edges_v(n_edges(k))
		uu_edges(k)=edges_u(n_edges(k))*edges_u(n_edges(k))
		uv_edges(k)=edges_u(n_edges(k))*edges_v(n_edges(k))
		vv_edges(k)=edges_v(n_edges(k))*edges_v(n_edges(k))
	enddo 
! ==============================================================================================
    du_dxi1=(edges_u(n_edges(2))-edges_u(n_edges(1)))/2.
    du_dxi2=(edges_u(n_edges(4))-edges_u(n_edges(3)))/2.

    dv_dxi1=(edges_v(n_edges(2))-edges_v(n_edges(1)))/2.
    dv_dxi2=(edges_v(n_edges(4))-edges_v(n_edges(3)))/2.

! ==============================================================================================

	duu_dx=dphi_dx2(uu_edges,Mi,Mj)
	duv_dx=dphi_dx2(uv_edges,Mi,Mj)
	duv_dy=dphi_dy2(uv_edges,Mi,Mj)
	dvv_dy=dphi_dy2(vv_edges,Mi,Mj)
	
! ==============================================================================================


	meshes_R_C_u(Mi,Mj)=-(meshes_u(Mi,Mj)*du_dxi1 +meshes_v(Mi,Mj)*du_dxi2)+duu_dx+duv_dy

	meshes_R_C_v(Mi,Mj)=-(meshes_u(Mi,Mj)*dv_dxi1 +meshes_v(Mi,Mj)*dv_dxi2)+duv_dx+dvv_dy




end subroutine calculate_R_c


subroutine calculate_R_c_new(Mi,Mj)
        integer, intent(in)         :: Mi,Mj
        real                        :: du_dxi1,du_dxi2,dv_dxi1,dv_dxi2
        real                        :: uc,vc
    integer, dimension(1:4)     :: n_edges
        integer                     :: k
    real, dimension(1:4)        :: u_edges,v_edges,uu_edges,uv_edges,vv_edges


        n_edges(:)=meshes_edges(:,Mi,Mj)
        do k=1,4
                u_edges(k) =edges_u(n_edges(k))
                v_edges(k) =edges_v(n_edges(k))
        enddo
! ==============================================================================================
    du_dxi1=(edges_u(n_edges(2))-edges_u(n_edges(1)))/2.
    du_dxi2=(edges_u(n_edges(4))-edges_u(n_edges(3)))/2.

    dv_dxi1=(edges_v(n_edges(2))-edges_v(n_edges(1)))/2.
    dv_dxi2=(edges_v(n_edges(4))-edges_v(n_edges(3)))/2.

! ==============================================================================================

        uc=meshes_uxyt(Mi,Mj)*meshes_xi_x(1,1,Mi,Mj)+meshes_vxyt(Mi,Mj)*meshes_xi_x(1,2,Mi,Mj)
        vc=meshes_uxyt(Mi,Mj)*meshes_xi_x(2,1,Mi,Mj)+meshes_vxyt(Mi,Mj)*meshes_xi_x(2,2,Mi,Mj)

! ==============================================================================================


        meshes_R_C_u(Mi,Mj)=-((meshes_u(Mi,Mj)-uc)*du_dxi1 +(meshes_v(Mi,Mj)-vc)*du_dxi2)

        meshes_R_C_v(Mi,Mj)=-((meshes_u(Mi,Mj)-uc)*dv_dxi1 +(meshes_v(Mi,Mj)-vc)*dv_dxi2)



end subroutine calculate_R_c_new



subroutine calculate_R_D(Mi,Mj)
	integer, intent(in)         :: Mi,Mj
	real                        :: ux1,vx1,ux2,vx2,uxy,vxy,uy1,uy2,vy1,vy2,uyx,vyx,ux,uy,vx,vy,&
								   sum1,sum2,sum3,sum4,cnt1,cnt2,cnt3,cnt4,sum12,sum22,sum32,sum42
    cnt1=0
    sum1=0.
    cnt2=0
    sum2=0.
    cnt3=0
    sum3=0.
    cnt4=0
    sum4=0.	
	sum12=0.
	sum22=0.
	sum32=0.
	sum42=0.
	! Find ux1,2, vx1,2
	if(Mi>1) then
		sum3 =sum3+(edges_u(meshes_edges(3,Mi,Mj))-edges_u(meshes_edges(3,Mi-1,Mj)))/2.
		sum32=sum32+(edges_v(meshes_edges(3,Mi,Mj))-edges_v(meshes_edges(3,Mi-1,Mj)))/2.
		cnt3=cnt3+1
		sum4 =sum4+(edges_u(meshes_edges(4,Mi,Mj))-edges_u(meshes_edges(4,Mi-1,Mj)))/2.
		sum42=sum42+(edges_v(meshes_edges(4,Mi,Mj))-edges_v(meshes_edges(4,Mi-1,Mj)))/2.
		cnt4=cnt4+1
	endif
	if(Mi<NMESH_i) then
		sum3=sum3+(edges_u(meshes_edges(3,Mi+1,Mj))-edges_u(meshes_edges(3,Mi,Mj)))/2.
		sum32=sum32+(edges_v(meshes_edges(3,Mi+1,Mj))-edges_v(meshes_edges(3,Mi,Mj)))/2.
		cnt3=cnt3+1
		sum4=sum4+(edges_u(meshes_edges(4,Mi+1,Mj))-edges_u(meshes_edges(4,Mi,Mj)))/2.
		sum42=sum42+(edges_v(meshes_edges(4,Mi+1,Mj))-edges_v(meshes_edges(4,Mi,Mj)))/2.
		cnt4=cnt4+1
	endif
	if(Mj>1) then
		sum1=sum1+(edges_u(meshes_edges(1,Mi,Mj))-edges_u(meshes_edges(1,Mi,Mj-1)))/2.
		sum12=sum12+(edges_v(meshes_edges(1,Mi,Mj))-edges_v(meshes_edges(1,Mi,Mj-1)))/2.
		cnt1=cnt1+1
		sum2=sum2+(edges_u(meshes_edges(2,Mi,Mj))-edges_u(meshes_edges(2,Mi,Mj-1)))/2.
		sum22=sum22+(edges_v(meshes_edges(2,Mi,Mj))-edges_v(meshes_edges(2,Mi,Mj-1)))/2.
		cnt2=cnt2+1
	endif
	if(Mj<NMESH_j) then
		sum1=sum1+(edges_u(meshes_edges(1,Mi,Mj+1))-edges_u(meshes_edges(1,Mi,Mj)))/2.
		sum12=sum12+(edges_v(meshes_edges(1,Mi,Mj+1))-edges_v(meshes_edges(1,Mi,Mj)))/2.
		cnt1=cnt1+1
		sum2=sum2+(edges_u(meshes_edges(2,Mi,Mj+1))-edges_u(meshes_edges(2,Mi,Mj)))/2.
		sum22=sum22+(edges_v(meshes_edges(2,Mi,Mj+1))-edges_v(meshes_edges(2,Mi,Mj)))/2.
		cnt2=cnt2+1
	endif
	sum1=sum1/cnt1
	sum2=sum2/cnt2
	sum3=sum3/cnt3
	sum4=sum4/cnt4

	sum12=sum12/cnt1
	sum22=sum22/cnt2
	sum32=sum32/cnt3
	sum42=sum42/cnt4
	
	uxy=(sum2-sum1)/2.
	uyx=(sum4-sum3)/2.
	vxy=(sum22-sum12)/2.
	vyx=(sum42-sum32)/2.
	
	ux=(edges_u(meshes_edges(2,Mi,Mj))-edges_u(meshes_edges(1,Mi,Mj)))/2.
	uy=(edges_u(meshes_edges(4,Mi,Mj))-edges_u(meshes_edges(3,Mi,Mj)))/2.

	vx=(edges_v(meshes_edges(2,Mi,Mj))-edges_v(meshes_edges(1,Mi,Mj)))/2.
	vy=(edges_v(meshes_edges(4,Mi,Mj))-edges_v(meshes_edges(3,Mi,Mj)))/2.


	meshes_R_D_u(Mi,Mj)=1/Re*(meshes_G1(Mi,Mj)*ux+meshes_G2(Mi,Mj)*uy-meshes_gij(2,1,Mi,Mj)*uxy-meshes_gij(1,2,Mi,Mj)*uyx)

	meshes_R_D_v(Mi,Mj)=1/Re*(meshes_G1(Mi,Mj)*vx+meshes_G2(Mi,Mj)*vy-meshes_gij(2,1,Mi,Mj)*vxy-meshes_gij(1,2,Mi,Mj)*vyx)

end subroutine calculate_R_D


subroutine calculate_R_b(Mi,Mj)
	integer, intent(in)         :: Mi,Mj
	real                        :: dp_dx, dp_dy
    integer, dimension(1:4)     :: n_edges
    real, dimension(1:4)        :: p_edges
	integer                     :: k
	
	n_edges(:)=meshes_edges(:,Mi,Mj)
	do k=1,4
		p_edges(k)=edges_p(n_edges(k))
	enddo 
! ==============================================================================================
	dp_dx=dphi_dx2(p_edges,Mi,Mj)
	dp_dy=dphi_dy2(p_edges,Mi,Mj)
! ==============================================================================================
	meshes_R_B_u(Mi,Mj)=meshes_bx(Mi,Mj)+1/rho*dp_dx
	meshes_R_B_v(Mi,Mj)=meshes_by(Mi,Mj)+1/rho*dp_dy
end subroutine calculate_R_b




subroutine calculate_R_b_p(Mi,Mj)
	integer, intent(in)         :: Mi,Mj
	real                        :: du_dx, du_dy, dv_dx, dv_dy
	real                        :: dbx_dx, dby_dy, Dila
    integer, dimension(1:4)     :: n_edges
	integer                     :: k
    real, dimension(1:4)        :: u_edges,v_edges,bx_edges,by_edges

	n_edges(:)=meshes_edges(:,Mi,Mj)
	do k=1,4
		u_edges(k)=edges_u(n_edges(k))
		v_edges(k)=edges_v(n_edges(k))
		bx_edges(k)=edges_bx(n_edges(k))
		by_edges(k)=edges_by(n_edges(k))
	enddo 
	du_dx=dphi_dx2(u_edges,Mi,Mj)
	du_dy=dphi_dy2(u_edges,Mi,Mj)
	dv_dx=dphi_dx2(v_edges,Mi,Mj)
	dv_dy=dphi_dy2(v_edges,Mi,Mj)
	dbx_dx=dphi_dx2(bx_edges,Mi,Mj)
	dby_dy=dphi_dy2(by_edges,Mi,Mj)
	Dila=du_dx+dv_dy

	meshes_dila(Mi,Mj)=Dila

	meshes_R_b_p(Mi,Mj)=-(du_dx**2. + dv_dy**2.)-2.0*(du_dy*dv_dx)-dbx_dx-dby_dy+Dila/(2*DT)

end subroutine calculate_R_b_p


subroutine calculate_R_D_p(Mi,Mj)
	integer, intent(in)         :: Mi,Mj
	real                        :: px1,px2,pxy,py1,py2,pyx,px,py,sum1,sum2,sum3,sum4,cnt1,cnt2,cnt3,cnt4
    cnt1=0
    sum1=0.
    cnt2=0
    sum2=0.
    cnt3=0
    sum3=0.
    cnt4=0
    sum4=0.	
	! Find ux1,2, vx1,2
	if(Mi>1) then
		sum3=sum3+(edges_p(meshes_edges(3,Mi,Mj))-edges_p(meshes_edges(3,Mi-1,Mj)))/2.
		cnt3=cnt3+1
		sum4=sum4+(edges_p(meshes_edges(4,Mi,Mj))-edges_p(meshes_edges(4,Mi-1,Mj)))/2.
		cnt4=cnt4+1
	endif
	if(Mi<NMESH_i) then
		sum3=sum3+(edges_p(meshes_edges(3,Mi+1,Mj))-edges_p(meshes_edges(3,Mi,Mj)))/2.
		cnt3=cnt3+1
		sum4=sum4+(edges_p(meshes_edges(4,Mi+1,Mj))-edges_p(meshes_edges(4,Mi,Mj)))/2.
		cnt4=cnt4+1
	endif
	if(Mj>1) then
		sum1=sum1+(edges_p(meshes_edges(1,Mi,Mj))-edges_p(meshes_edges(1,Mi,Mj-1)))/2.
		cnt1=cnt1+1
		sum2=sum2+(edges_p(meshes_edges(2,Mi,Mj))-edges_p(meshes_edges(2,Mi,Mj-1)))/2.
		cnt2=cnt2+1
	endif
	if(Mj<NMESH_j) then
		sum1=sum1+(edges_p(meshes_edges(1,Mi,Mj+1))-edges_p(meshes_edges(1,Mi,Mj)))/2.
		cnt1=cnt1+1
		sum2=sum2+(edges_p(meshes_edges(2,Mi,Mj+1))-edges_p(meshes_edges(2,Mi,Mj)))/2.
		cnt2=cnt2+1
	endif
	sum1=sum1/cnt1
	sum2=sum2/cnt2
	sum3=sum3/cnt3
	sum4=sum4/cnt4

	
	pxy=(sum2-sum1)/2.
	pyx=(sum4-sum3)/2.

	px=(edges_p(meshes_edges(2,Mi,Mj))-edges_p(meshes_edges(1,Mi,Mj)))/2.
	py=(edges_p(meshes_edges(4,Mi,Mj))-edges_p(meshes_edges(3,Mi,Mj)))/2.

	meshes_R_D_p(Mi,Mj)=1/rho*(meshes_G1(Mi,Mj)*px+meshes_G2(Mi,Mj)*py-meshes_gij(2,1,Mi,Mj)*pxy-meshes_gij(1,2,Mi,Mj)*pyx)


end subroutine calculate_R_D_p

subroutine evaluate_D
	integer                     :: i,j,k
	real                        :: Dila
    integer, dimension(1:4)     :: n_edges
    real, dimension(1:4)        :: u_edges,v_edges
	real                        :: sum1
	
	sum1=0.0
	!$OMP PARALLEL DO PRIVATE (i,j,k,n_edges,u_edges,v_edges,Dila) REDUCTION(+: sum1)
	do j=1,NMESH_j
		do i=1,NMESH_i
			n_edges(:)=meshes_edges(:,i,j)
			do k=1,4
			u_edges(k)=edges_u(n_edges(k))
			v_edges(k)=edges_v(n_edges(k))
			enddo
			Dila=dphi_dx2(u_edges,i,j)+dphi_dy2(v_edges,i,j)
			meshes_dila(i,j)=Dila
			
			sum1=sum1+Dila**2

		enddo
	enddo
    !$OMP     END PARALLEL DO
	Dilatation=sqrt(sum1/(NMESH_i*NMESH_j))
end subroutine evaluate_D


subroutine find_rms
	integer        :: i,j
	real           :: sum1, sum2
	
	sum1=0.
	sum2=0.
	!$OMP     PARALLEL DO  PRIVATE (i,j) REDUCTION(+: sum1,sum2)
	do j=1,NMESH_j
		do i=1,NMESH_i
			sum1=sum1+abs(meshes_uT(i,j)-meshes_uxy(i,j))**2.
			sum2=sum2+abs(meshes_vT(i,j)-meshes_vxy(i,j))**2.
		enddo
	enddo
    !$OMP     END PARALLEL DO

	RMS_meshes_u=sqrt(sum1/(NMESH_i*NMESH_j))
	RMS_meshes_v=sqrt(sum2/(NMESH_i*NMESH_j))

	sum1=0.
	sum2=0.
	!$OMP     PARALLEL DO  PRIVATE (i) REDUCTION(+: sum1,sum2)
	do i=1,NEDGE
		sum1=sum1+abs(edges_uT(i)-edges_u(i))**2.
		sum2=sum2+abs(edges_vT(i)-edges_v(i))**2.
	enddo
    !$OMP     END PARALLEL DO
	RMS_edges_u=sqrt(sum1/NEDGE)
	RMS_edges_v=sqrt(sum2/NEDGE)
	print*, "================================================================================="
	write (*,'(A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4, A, 1p,E11.4)') & 
	'RMS uxy=', RMS_meshes_u , ' |  RMS vxy=', RMS_meshes_v, ' |  RMS edges u=', RMS_edges_u, ' |  RMS edges v=', RMS_edges_v
	print*, "================================================================================="

end subroutine find_rms

! This function is used to find dphi_dx in cartesian system from the derivatives in curvilinear coordinates
function dphi_dx(dphi_dxi1,dphi_dxi2,Mi,Mj) result(T)
	real, intent(in)    :: dphi_dxi1,dphi_dxi2
	integer, intent(in) :: Mi,Mj 							! mesh umber
	real                :: T
    T= meshes_xi_x(1,1,Mi,Mj)*dphi_dxi1+meshes_xi_x(2,1,Mi,Mj)*dphi_dxi2
end function dphi_dx

function dphi_dy(dphi_dxi1,dphi_dxi2,Mi,Mj) result(T)
	real, intent(in)    :: dphi_dxi1,dphi_dxi2
	integer, intent(in) :: Mi,Mj 							! mesh umber
	real                :: T
    T= meshes_xi_x(1,2,Mi,Mj)*dphi_dxi1+meshes_xi_x(2,2,Mi,Mj)*dphi_dxi2
end function dphi_dy
! To find first derivative with respect to x in cartesian coordinates using the corner points of a mesh
function dphi_dx2(phi,Mi,Mj) result(T)
	real,dimension(1:4),intent(in):: phi
	integer, intent(in)           :: Mi,Mj 							! mesh umber
	real                          :: T
    T= 1./meshes_area(Mi,Mj)*(phi(2)*(meshes_nodes_y(3,Mi,Mj)-meshes_nodes_y(2,Mi,Mj)) &
	                         -phi(1)*(meshes_nodes_y(4,Mi,Mj)-meshes_nodes_y(1,Mi,Mj)) &
	                         +phi(4)*(meshes_nodes_y(4,Mi,Mj)-meshes_nodes_y(3,Mi,Mj)) &
	                         -phi(3)*(meshes_nodes_y(1,Mi,Mj)-meshes_nodes_y(2,Mi,Mj)) )
end function dphi_dx2

function dphi_dy2(phi,Mi,Mj) result(T)
	real,dimension(1:4),intent(in):: phi
	integer, intent(in)           :: Mi,Mj 							! mesh umber
	real                          :: T
    T= 1./meshes_area(Mi,Mj)*(phi(2)*(meshes_nodes_x(2,Mi,Mj)-meshes_nodes_x(3,Mi,Mj)) &
	                         -phi(1)*(meshes_nodes_x(1,Mi,Mj)-meshes_nodes_x(4,Mi,Mj)) &
	                         +phi(4)*(meshes_nodes_x(3,Mi,Mj)-meshes_nodes_x(4,Mi,Mj)) &
	                         -phi(3)*(meshes_nodes_x(2,Mi,Mj)-meshes_nodes_x(1,Mi,Mj)) )
end function dphi_dy2
end module nse

