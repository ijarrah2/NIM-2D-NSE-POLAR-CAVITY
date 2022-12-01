module geometry
use mesh
use general
use curvilinear
use user_parameters
use parameters
implicit none
contains

subroutine create_geometry
  implicit none
	call build_geometry
	call fill_nodes_array
	call fill_edges_numbers
	call deform_grid
	!  call check_edge_average
end subroutine create_geometry

subroutine evaluate_all_mesh_constants
  implicit none
  integer :: i,j

!$OMP     PARALLEL DO PRIVATE (j,i)
	do j=1,NMESH_j
		do i=1,NMESH_i
			call calculate_mapping_coefficients_2D(MESHES_ARRAY(i,j))
			call find_area_2d(i,j)
			call find_normal(i,j)
			call calculate_defusion_tensor_over_edges(i,j)
			call convection_diffusion_constants_mesh_i(i,j)
		enddo
	enddo
!$OMP     END PARALLEL DO

end subroutine evaluate_all_mesh_constants



subroutine evaluate_time_dependent_quantities
  implicit none
  integer :: i,j

  if(.not. steady_state .and. (is_bx_TD .or. is_by_TD .or. is_analytical_TD)) then
  !$OMP     PARALLEL DO PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
! Not implemented yet

!		   call time_dependent_face_averages(i,j)
!		   call time_dependent_mesh_averages(i,j)
		enddo
	enddo
  !$OMP     END PARALLEL DO
  endif

end subroutine evaluate_time_dependent_quantities



subroutine deform_grid
  implicit none

!  call h2_perturb
!  call skew_meshes(0.,0.5)
!  call cube_to_cylinder(0.5,1.5)  
!  call h2_perturb_cylinder(0.5,1.5)
   call refine_walls_quadratic(2.,2.)
   call cube_to_cylinder2(1.0,2.0,-0.5,0.5)  
   
!  call random_perturb
end subroutine deform_grid


subroutine h2_perturb
  implicit none
  integer     i,j,k
  real  x,y,PI

  PI=4.D0*DATAN(1.D0)
  do i=1,NMESH_i+1
    do j=1,NMESH_j+1
       x=NODES_ARRAY(i,j)%x
       y=NODES_ARRAY(i,j)%y
       NODES_ARRAY(i,j)%x=x+0.1*sin(2*PI*x)*sin(2*PI*y)
       NODES_ARRAY(i,j)%y=y-0.1*sin(2*PI*x)*sin(2*PI*y)
    enddo
  enddo
  call fill_mesh_nodes_from_nodes_array
end subroutine h2_perturb

subroutine cube_to_cylinder2(R_i,R_o,theta1,theta2)
! The cube should be 1x1x1 centered at 0,0,0
  real, intent(in)                           :: R_i,R_o,theta1,theta2
  real                                       :: RR,theta,PI,x,y
  integer i,j,k
  do i=1,NMESH_i
    do j=1,NMESH_j
       do k=1,4
         PI=4.D0*DATAN(1.D0)
         x=MESHES_ARRAY(i,j)%nodes(k)%x
         y=MESHES_ARRAY(i,j)%nodes(k)%y
         call geometry_dilation(x,y,(R_o-R_i),(R_o-R_i))
         call move_geometry(x,y,(R_i),0.0)
         RR= x
         theta=(theta2-theta1)*(y)/(R_o-R_i)+theta1!*2
         MESHES_ARRAY(i,j)%nodes(k)%y=RR*sin(theta)
         MESHES_ARRAY(i,j)%nodes(k)%x=RR*cos(theta)
       enddo
    enddo
  enddo

end subroutine cube_to_cylinder2

subroutine refine_walls_quadratic(b1,b2)
! The cube should be 1x1x1 centered at 0,0,0
  real, intent(in)                           :: b1,b2
  real                                       :: x,y
  integer i,j,k
  do i=1,NMESH_i
    do j=1,NMESH_j
       do k=1,4
         x=MESHES_ARRAY(i,j)%nodes(k)%x
         y=MESHES_ARRAY(i,j)%nodes(k)%y
		 if(x<=0.5) MESHES_ARRAY(i,j)%nodes(k)%x=((2*x)**b1)/2.
		 if(x>0.5 .and. x<1.)  MESHES_ARRAY(i,j)%nodes(k)%x=1.-((-2.*(x-1))**b1)/2.
		 if(y<=0.5) MESHES_ARRAY(i,j)%nodes(k)%y=((2*y)**b2)/2.
		 if(y>0.5 .and. y<1.)  MESHES_ARRAY(i,j)%nodes(k)%y=1.-((-2.*(y-1))**b2)/2.
       enddo
    enddo
  enddo

end subroutine refine_walls_quadratic


subroutine random_perturb
  implicit none
  integer     i,j,k
  real  x,y,PI,num,dx,dy,a,b
  integer,parameter :: seed = 86456
  
  call srand(seed)
  a=0.01
  b=0.01
  do i=2,NMESH_i
    do j=2,NMESH_j
          CALL RANDOM_NUMBER(num)
          dx=MESHES_ARRAY(i,j)%a*2
          dy=MESHES_ARRAY(i,j)%b*2
          x=NODES_ARRAY(i,j)%x
          y=NODES_ARRAY(i,j)%y
          NODES_ARRAY(i,j)%x=x+(num-0.5)*a*dx
          CALL RANDOM_NUMBER(num)
          NODES_ARRAY(i,j)%y=y+(num-0.5)*b*dy
    enddo
  enddo

  call fill_mesh_nodes_from_nodes_array



end subroutine random_perturb



subroutine skew_meshes(tan_x,tan_y)
  real, intent(in)                           :: tan_x,tan_y
  real, dimension(1:3)                       :: temp_old
  real  x,y,PI
  integer i,j,k
  do i=1,NMESH_i
    do j=1,NMESH_j
       do k=1,4
          x=MESHES_ARRAY(i,j)%nodes(k)%x
          y=MESHES_ARRAY(i,j)%nodes(k)%y
          MESHES_ARRAY(i,j)%nodes(k)%x=y*tan_x+x
          MESHES_ARRAY(i,j)%nodes(k)%y=x*tan_y+y
       enddo
    enddo
  enddo

end subroutine skew_meshes


subroutine cube_to_cylinder(R_i,R_o)
! The cube should be 1x1x1 centered at 0,0,0
  real, intent(in)                           :: R_i,R_o
  real                                       :: RR,theta,PI,x,y
  integer i,j,k
  do i=1,NMESH_i
    do j=1,NMESH_j
       do k=1,4
         PI=4.D0*DATAN(1.D0)
         x=MESHES_ARRAY(i,j)%nodes(k)%x
         y=MESHES_ARRAY(i,j)%nodes(k)%y
         call geometry_dilation(x,y,(R_o-R_i),(R_o-R_i))
         call move_geometry(x,y,(R_o-R_i)/2.,(R_o-R_i)/2.)
         RR= x
         theta=PI/2.*y/(R_o-R_i)*2
         MESHES_ARRAY(i,j)%nodes(k)%y=RR*sin(theta)
         MESHES_ARRAY(i,j)%nodes(k)%x=-RR*cos(theta)
       enddo
    enddo
  enddo

end subroutine cube_to_cylinder



subroutine h2_perturb_cylinder(R_i,R_o)
  implicit none
  real, intent(in) :: R_i,R_o
  integer     i,j,k
  real  x,y,PI,r,alpha,beta

  alpha=0.06
  beta=-0.05
  PI=4.D0*DATAN(1.D0)
  do i=1,NMESH_i
    do j=1,NMESH_j
       do k=1,4
          x=MESHES_ARRAY(i,j)%nodes(k)%x
          y=MESHES_ARRAY(i,j)%nodes(k)%y
          r=sqrt(x**2+y**2)
          MESHES_ARRAY(i,j)%nodes(k)%x=x+alpha*sin(2*PI*(r-R_i))*sin(PI*x)
          MESHES_ARRAY(i,j)%nodes(k)%y=y+beta *sin(2*PI*(r-R_i))*sin(PI*x)

       enddo
    enddo
  enddo

end subroutine h2_perturb_cylinder



subroutine geometry_dilation(x,y,xx,xy)
  real, intent(inout)         :: x,y
  real, intent(in)            :: xx,xy
  x=x*xx
  y=y*xy
end subroutine geometry_dilation

subroutine move_geometry(x,y,dx,dy)
  real, intent(inout)         :: x,y
  real, intent(in)            :: dx,dy
  x=x+dx
  y=y+dy
end subroutine move_geometry




subroutine calculate_mapping_coefficients_2D(facei)
  implicit none
  type(cell), intent(inout)          :: facei
  facei%C(1)=( facei%nodes(1)%x&
              +facei%nodes(2)%x&
              +facei%nodes(3)%x&
              +facei%nodes(4)%x)/4.d0

  facei%C(2)=(-facei%nodes(1)%x&
              +facei%nodes(2)%x&
              +facei%nodes(3)%x&
              -facei%nodes(4)%x)/4.d0

  facei%C(3)=(-facei%nodes(1)%x&
              -facei%nodes(2)%x&
              +facei%nodes(3)%x&
              +facei%nodes(4)%x)/4.d0

  facei%C(4)=(+facei%nodes(1)%x&
              -facei%nodes(2)%x&
              +facei%nodes(3)%x&
              -facei%nodes(4)%x)/4.d0

  facei%D(1)=( facei%nodes(1)%y&
              +facei%nodes(2)%y&
              +facei%nodes(3)%y&
              +facei%nodes(4)%y)/4.d0

  facei%D(2)=(-facei%nodes(1)%y&
              +facei%nodes(2)%y&
              +facei%nodes(3)%y&
              -facei%nodes(4)%y)/4.d0

  facei%D(3)=(-facei%nodes(1)%y&
              -facei%nodes(2)%y&
              +facei%nodes(3)%y&
              +facei%nodes(4)%y)/4.d0

  facei%D(4)=( facei%nodes(1)%y&
              -facei%nodes(2)%y&
              +facei%nodes(3)%y&
              -facei%nodes(4)%y)/4.d0

end subroutine calculate_mapping_coefficients_2D


subroutine find_area_2d(Mi,Mj)
  implicit none
  integer, intent(in)                   :: Mi,Mj
  real,dimension(1:N_Gauss,1:N_Gauss)	:: A1
  real									:: x,y
  integer								:: i,j
  type(curvi2d)							:: Curvidata

  do i=1,N_Gauss
     do j=1,N_Gauss
        call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),Gauss_x_i(j),x,y)
        call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),curvidata,Gauss_x_i(i),Gauss_x_i(j))
        A1(i,j)= curvidata%detJ
     enddo
  enddo
  meshes_area(Mi,Mj)=Gauss_integration2d(A1)

end subroutine find_area_2d

subroutine find_normal(i,j)
! Outward normals
  implicit none
  integer,intent(in)                   :: i,j
  real                                 :: x,y, x_i, y_i,nx,ny,tx,ty
  type(curvi2d)                        :: Curvidata

  ! left x_i=-1.0
  x_i=-1.0
  y_i= 0.0
  call map_xy_2d(MESHES_ARRAY(i,j),x_i,y_i,x,y)
  call evaluate_curvilinear(MESHES_ARRAY(i,j),curvidata,x_i,y_i)
  nx=curvidata%xi_x(1,1)/sqrt(curvidata%gij(1,1))
  ny=curvidata%xi_x(1,2)/sqrt(curvidata%gij(1,1))
  tx=curvidata%x_xi(1,2)/sqrt(curvidata%g_ij(2,2))
  ty=curvidata%x_xi(2,2)/sqrt(curvidata%g_ij(2,2))

  edges_normal(:,meshes_edges(1,i,j))=(/nx,ny/)
  edges_tangential(:,meshes_edges(1,i,j))=(/tx,ty/)
  edges_mid_coordinates(1:2,meshes_edges(1,i,j))=(/x,y/)
  ! Right
  x_i=1.0
  y_i=0.0
  call map_xy_2d(MESHES_ARRAY(i,j),x_i,y_i,x,y)
  call evaluate_curvilinear(MESHES_ARRAY(i,j),curvidata,x_i,y_i)
  nx=curvidata%xi_x(1,1)/sqrt(curvidata%gij(1,1))
  ny=curvidata%xi_x(1,2)/sqrt(curvidata%gij(1,1))
  tx=curvidata%x_xi(1,2)/sqrt(curvidata%g_ij(2,2))
  ty=curvidata%x_xi(2,2)/sqrt(curvidata%g_ij(2,2))

  edges_normal(:,meshes_edges(2,i,j))=(/nx,ny/)
  edges_tangential(:,meshes_edges(2,i,j))=(/tx,ty/)
  edges_mid_coordinates(1:2,meshes_edges(2,i,j))=(/x,y/)

! Bottom
  x_i= 0.0
  y_i=-1.0
  call map_xy_2d(MESHES_ARRAY(i,j),x_i,y_i,x,y)
  call evaluate_curvilinear(MESHES_ARRAY(i,j),curvidata,x_i,y_i)


  nx=curvidata%xi_x(2,1)/sqrt(curvidata%gij(2,2))
  ny=curvidata%xi_x(2,2)/sqrt(curvidata%gij(2,2))
  tx=curvidata%x_xi(1,1)/sqrt(curvidata%g_ij(1,1))
  ty=curvidata%x_xi(2,1)/sqrt(curvidata%g_ij(1,1))

  edges_normal(:,meshes_edges(3,i,j))=(/nx,ny/)
  edges_tangential(:,meshes_edges(3,i,j))=(/tx,ty/)
  edges_mid_coordinates(1:2,meshes_edges(3,i,j))=(/x,y/)
  
  ! Top
  x_i= 0.0
  y_i= 1.0
  call map_xy_2d(MESHES_ARRAY(i,j),x_i,y_i,x,y)
  call evaluate_curvilinear(MESHES_ARRAY(i,j),curvidata,x_i,y_i)
  nx=curvidata%xi_x(2,1)/sqrt(curvidata%gij(2,2))
  ny=curvidata%xi_x(2,2)/sqrt(curvidata%gij(2,2))
  tx=curvidata%x_xi(1,1)/sqrt(curvidata%g_ij(1,1))
  ty=curvidata%x_xi(2,1)/sqrt(curvidata%g_ij(1,1))

  edges_normal(:,meshes_edges(4,i,j))=(/nx,ny/)
  edges_tangential(:,meshes_edges(4,i,j))=(/tx,ty/)
  edges_mid_coordinates(1:2,meshes_edges(4,i,j))=(/x,y/)
  
end subroutine find_normal

  subroutine calculate_defusion_tensor_over_edges(i,j)
! -----------------------------------
    integer, intent(in)             :: i,j
    real, dimension(1:Gauss_int_p)  :: ones_plus,ones_minus
    integer                         :: k

! -----------------------------------
    do k=1,Gauss_int_p
       ones_plus(k)=1.
       ones_minus(k)=-1.
    enddo
    call calculate_diffusion_tensor_edges_i(i,j,ones_minus,Gauss_x_i,1)
    call calculate_diffusion_tensor_edges_i(i,j,ones_plus ,Gauss_x_i,2)
    call calculate_diffusion_tensor_edges_i(i,j,Gauss_x_i,ones_minus,3)
    call calculate_diffusion_tensor_edges_i(i,j,Gauss_x_i,ones_plus ,4)
 
 end subroutine calculate_defusion_tensor_over_edges

  subroutine calculate_diffusion_tensor_edges_i(Mi,Mj,x_array,y_array,n)
    implicit none
! -----------------------------------
    integer, intent(in)                         :: Mi,Mj
    integer, intent(in)                         :: n
    real, intent(in), dimension(1:Gauss_int_p)  :: x_array,y_array
    real, dimension (1:Gauss_int_p,1:Gauss_int_p) :: K11,K12,K21,K22,dS,u_i, v_i, p_i,bx_i,by_i
    integer   :: i,j,k
    real :: length,x,y
    type(curvi2d)    :: C
! Calculate diffusion tensor, area of faces, and analytical over faces

! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
			call map_xy_2d(MESHES_ARRAY(Mi,Mj),x_array(i),y_array(j),x,y)
			call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,x_array(i),y_array(j))
			if(n==1 .or. n==2) then
			   dS(i,j)=1.*sqrt(C%g_ij(2,2))
			elseif(n==3 .or. n==4) then
			   dS(i,j)=1.*sqrt(C%g_ij(1,1))
			endif
			K11(i,j)=C%gij(1,1)/sqrt(C%gij(1,1))*dS(i,j)
			K12(i,j)=C%gij(1,2)/sqrt(C%gij(1,1))*dS(i,j)
			K21(i,j)=C%gij(2,1)/sqrt(C%gij(2,2))*dS(i,j)
			K22(i,j)=C%gij(2,2)/sqrt(C%gij(2,2))*dS(i,j)
			if(.not. is_analytical_TD .and. is_analytical) then
				u_i(i,j)=functions_id (x,y,0.0,velx1)*dS(i,j)
				v_i(i,j)=functions_id (x,y,0.0,velx2)*dS(i,j)
				p_i(i,j)=functions_id (x,y,0.0,PPE)  *dS(i,j)
			endif
			if(is_bx .and. .not. is_bx_TD) bx_i(i,j)=functions_id (x,y,0.0,bx_fun)*dS(i,j)
			if(is_by .and. .not. is_by_TD) by_i(i,j)=functions_id (x,y,0.0,by_fun)*dS(i,j)
       enddo
    enddo
	length=Gauss_integration1d_t(dS)
	edges_length(meshes_edges(n,Mi,Mj))=length
	if(.not. is_analytical_TD .and. is_analytical) then
		edges_uT(meshes_edges(n,Mi,Mj))= Gauss_integration1d_t(u_i)/length
		edges_vT(meshes_edges(n,Mi,Mj))= Gauss_integration1d_t(v_i)/length
		edges_pT(meshes_edges(n,Mi,Mj))= Gauss_integration1d_t(p_i)/length
	endif
	if(is_bx .and. .not. is_bx_TD) edges_bx(meshes_edges(n,Mi,Mj))= Gauss_integration1d_t(bx_i)/length
	if(is_by .and. .not. is_by_TD) edges_by(meshes_edges(n,Mi,Mj))= Gauss_integration1d_t(by_i)/length
	
    meshes_gij_edges(1,1,n,Mi,Mj)= Gauss_integration1d_t(K11)/length
    meshes_gij_edges(1,2,n,Mi,Mj)= Gauss_integration1d_t(K12)/length
    meshes_gij_edges(2,1,n,Mi,Mj)= Gauss_integration1d_t(K21)/length
    meshes_gij_edges(2,2,n,Mi,Mj)= Gauss_integration1d_t(K22)/length
  end subroutine calculate_diffusion_tensor_edges_i


  subroutine convection_diffusion_constants_mesh_i (Mi,Mj)
! -----------------------------------
    integer, intent (in)      :: Mi,Mj
    real, dimension (1:Gauss_int_p,1:Gauss_int_p) :: K11_i,K22_i,K12_i,K21_i,G_11,G_22,bx_i,by_i, &
												     J11_i,J22_i,J12_i,J21_i, u_ic, v_ic, p_ic,u_T, v_T, p_T
    real                                          :: G11,G22,bx,by,x,y
    integer                                       :: i,j,k
    type(curvi2d)                                 :: curvidata

! -----------------------------------
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
			call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),Gauss_x_i(j),x,y)
			call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),curvidata,Gauss_x_i(i),Gauss_x_i(j))
			G11=(curvidata%gamma_tensor(1,1,2)*curvidata%gij(1,2)+curvidata%gamma_tensor(1,2,1)*curvidata%gij(2,1))
			G22=(curvidata%gamma_tensor(2,1,2)*curvidata%gij(1,2)+curvidata%gamma_tensor(2,2,1)*curvidata%gij(2,1))

! ============================= Filling arrays ==================================
			G_11(i,j) =G11    *curvidata%detJ
			G_22(i,j) =G22    *curvidata%detJ
			K11_i(i,j)=curvidata%gij(1,1)*curvidata%detJ
			K12_i(i,j)=curvidata%gij(1,2)*curvidata%detJ
			K21_i(i,j)=curvidata%gij(2,1)*curvidata%detJ
			K22_i(i,j)=curvidata%gij(2,2)*curvidata%detJ


			J11_i(i,j)=curvidata%jacinv(1,1)*curvidata%detJ
			J12_i(i,j)=curvidata%jacinv(1,2)*curvidata%detJ
			J21_i(i,j)=curvidata%jacinv(2,1)*curvidata%detJ
			J22_i(i,j)=curvidata%jacinv(2,2)*curvidata%detJ

			
			! Initial conditions
			u_ic(i,j)=functions_id(x,y,0.0,IC_velx1)*curvidata%detJ
			v_ic(i,j)=functions_id(x,y,0.0,IC_velx2)*curvidata%detJ
			p_ic(i,j)=functions_id(x,y,0.0,IC_PPE)  *curvidata%detJ
			if(.not. is_analytical_TD .and. is_analytical) then
				u_T(i,j)=functions_id (x,y,0.0,analytical_velx1)*curvidata%detJ
				v_T(i,j)=functions_id (x,y,0.0,analytical_velx2)*curvidata%detJ
				p_T(i,j)=functions_id (x,y,0.0,analytical_p)*curvidata%detJ
			endif			
			
			if(is_bx .and. .not. is_bx_TD) bx_i(i,j)=functions_id (x,y,0.0,bx_fun)*curvidata%detJ
			if(is_by .and. .not. is_bx_TD) by_i(i,j)=functions_id (x,y,0.0,by_fun)*curvidata%detJ

       enddo
    enddo
! ================================================================
! initial conditions
! u
    meshes_uxy_old(Mi,Mj)=Gauss_integration2d(u_ic)/meshes_area(Mi,Mj)
	meshes_uxy(Mi,Mj)=meshes_uxy_old(Mi,Mj)
	meshes_uxyt(Mi,Mj)=meshes_uxy_old(Mi,Mj)
! v	
    meshes_vxy_old(Mi,Mj)=Gauss_integration2d(v_ic)/meshes_area(Mi,Mj)
	meshes_vxy(Mi,Mj)=meshes_vxy_old(Mi,Mj)
	meshes_vxyt(Mi,Mj)=meshes_vxy_old(Mi,Mj)
! p -> does not time-dependent, but in case I need to initialize it with a function.	
	meshes_p(Mi,Mj)=Gauss_integration2d(p_ic)/meshes_area(Mi,Mj)

	if(.not. is_analytical_TD .and. is_analytical) then
		meshes_uT(Mi,Mj)=Gauss_integration2d(u_T)/meshes_area(Mi,Mj)
		meshes_vT(Mi,Mj)=Gauss_integration2d(v_T)/meshes_area(Mi,Mj)
		meshes_pT(Mi,Mj)=Gauss_integration2d(p_T)/meshes_area(Mi,Mj)
	endif		
! ================================================================
    if(is_bx .and. .not. is_bx_TD) meshes_bx(Mi,Mj)=Gauss_integration2d(bx_i)/meshes_area(Mi,Mj)
    if(is_by .and. .not. is_by_TD) meshes_by(Mi,Mj)=Gauss_integration2d(by_i)/meshes_area(Mi,Mj)
! ================================================================
    meshes_G1(Mi,Mj)   =Gauss_integration2d(G_11)/meshes_area(Mi,Mj)
    meshes_G2(Mi,Mj)   =Gauss_integration2d(G_22)/meshes_area(Mi,Mj)
! ================================================================
    meshes_gij(1,1,Mi,Mj)=Gauss_integration2d(K11_i)/meshes_area(Mi,Mj)
    meshes_gij(1,2,Mi,Mj)=Gauss_integration2d(K12_i)/meshes_area(Mi,Mj)
    meshes_gij(2,1,Mi,Mj)=Gauss_integration2d(K21_i)/meshes_area(Mi,Mj)
    meshes_gij(2,2,Mi,Mj)=Gauss_integration2d(K22_i)/meshes_area(Mi,Mj)

! ===============================================================
	meshes_xi_x(1,1,Mi,Mj)=Gauss_integration2d(J11_i)/meshes_area(Mi,Mj)
    meshes_xi_x(2,1,Mi,Mj)=Gauss_integration2d(J12_i)/meshes_area(Mi,Mj)
    meshes_xi_x(1,2,Mi,Mj)=Gauss_integration2d(J21_i)/meshes_area(Mi,Mj)
    meshes_xi_x(2,2,Mi,Mj)=Gauss_integration2d(J22_i)/meshes_area(Mi,Mj)

  end subroutine convection_diffusion_constants_mesh_i



subroutine build_geometry
  implicit none
! -----------------------------------------------------
	integer                    :: i,j
	real           :: aa,bb
	real,dimension(1:NMESH_i) :: growth_array_x
	real,dimension(1:NMESH_j) :: growth_array_y


	if(mod(NMESH_i,2)==0)  aa=width/2./((1.-growth_rate_x**(NMESH_i/2))/(1.-growth_rate_x))
	if(mod(NMESH_i,2)/=0)  aa=width/(((1.-growth_rate_x**(NMESH_i/2+1))/(1.-growth_rate_x))+((1.-growth_rate_x**(NMESH_i/2))/(1.-growth_rate_x)))

	if(mod(NMESH_j,2)==0)  bb=length/2./((1.-growth_rate_y**(NMESH_j/2))/(1.-growth_rate_y))
	if(mod(NMESH_j,2)/=0)  bb=length/(((1.-growth_rate_y**(NMESH_j/2+1))/(1.-growth_rate_y))+((1.-growth_rate_y**(NMESH_j/2))/(1.-growth_rate_y)))
	if(growth_rate_x==1. .and. mod(NMESH_i,2)==0) aa=width/((NMESH_i))
	if(growth_rate_y==1. .and. mod(NMESH_j,2)==0) bb=length/((NMESH_j))
	if(growth_rate_x==1. .and. mod(NMESH_i,2)/=0) aa=width/((NMESH_i))
	if(growth_rate_y==1. .and. mod(NMESH_j,2)/=0) bb=length/((NMESH_j))

	 do i=1,NMESH_i/2
		growth_array_x(i)=i-1
	 enddo
	 do i=NMESH_i/2+1,NMESH_i
		growth_array_x(i)=NMESH_i-i!+1
	 enddo

	 do j=1,NMESH_j/2
		growth_array_y(j)=j-1
	 enddo
	 do j=NMESH_j/2+1,NMESH_j
		growth_array_y(j)=NMESH_j-j!+1
	 enddo



	do j=1,NMESH_j
		 do i=1,NMESH_i
			MESHES_ARRAY(i,j)%n=(j-1)*NMESH_i+i
			MESHES_ARRAY(i,j)%edges(1)=(j-1)*(NMESH_i+1)+i
			MESHES_ARRAY(i,j)%edges(2)=(j-1)*(NMESH_i+1)+i+1
			MESHES_ARRAY(i,j)%edges(3)=(NMESH_j)*(NMESH_i+1)+i+(j-1)*(NMESH_i)
			MESHES_ARRAY(i,j)%edges(4)=(NMESH_j)*(NMESH_i+1)+i+(j)*(NMESH_i)
			if(i==1) then
				MESHES_ARRAY(i,j)%nodes(1)%x=START_X
				MESHES_ARRAY(i,j)%nodes(2)%x=aa+START_X
				MESHES_ARRAY(i,j)%nodes(3)%x=aa+START_X
				MESHES_ARRAY(i,j)%nodes(4)%x=START_X
			else
				MESHES_ARRAY(i,j)%nodes(1)%x=MESHES_ARRAY(i-1,j)%nodes(3)%x
				MESHES_ARRAY(i,j)%nodes(2)%x=growth_rate_x**(growth_array_x(i))*aa+MESHES_ARRAY(i-1,j)%nodes(3)%x
				MESHES_ARRAY(i,j)%nodes(3)%x=growth_rate_x**(growth_array_x(i))*aa+MESHES_ARRAY(i-1,j)%nodes(3)%x
				MESHES_ARRAY(i,j)%nodes(4)%x=MESHES_ARRAY(i-1,j)%nodes(3)%x
			endif
			if(j==1) then
				MESHES_ARRAY(i,j)%nodes(1)%y=START_Y
				MESHES_ARRAY(i,j)%nodes(2)%y=START_Y
				MESHES_ARRAY(i,j)%nodes(3)%y=bb+START_Y
				MESHES_ARRAY(i,j)%nodes(4)%y=bb+START_Y
			else
				MESHES_ARRAY(i,j)%nodes(1)%y=MESHES_ARRAY(i,j-1)%nodes(3)%y
				MESHES_ARRAY(i,j)%nodes(2)%y=MESHES_ARRAY(i,j-1)%nodes(3)%y
				MESHES_ARRAY(i,j)%nodes(3)%y=growth_rate_y**(growth_array_y(j))*bb+MESHES_ARRAY(i,j-1)%nodes(3)%y
				MESHES_ARRAY(i,j)%nodes(4)%y=growth_rate_y**(growth_array_y(j))*bb+MESHES_ARRAY(i,j-1)%nodes(3)%y
			endif

			MESHES_ARRAY(i,j)%a=(MESHES_ARRAY(i,j)%nodes(2)%x-MESHES_ARRAY(i,j)%nodes(1)%x)/2.
			MESHES_ARRAY(i,j)%b=(MESHES_ARRAY(i,j)%nodes(3)%y-MESHES_ARRAY(i,j)%nodes(2)%y)/2.
			MESHES_ARRAY(i,j)%i=i
			MESHES_ARRAY(i,j)%j=j
		 enddo
	enddo

end subroutine build_geometry



subroutine fill_edges_numbers
  implicit none
  integer i
	do i=1,NEDGE
		EDGES_ARRAY(i)%n=i
	enddo


end subroutine fill_edges_numbers

subroutine fill_nodes_array
	implicit none
	integer i,j
	do j=1, NMESH_j
		do i=1,NMESH_i
			NODES_ARRAY(i,j)%x=MESHES_ARRAY(i,j)%nodes(1)%x
			NODES_ARRAY(i,j)%y=MESHES_ARRAY(i,j)%nodes(1)%y
		enddo
		NODES_ARRAY(NMESH_i+1,j)%x=MESHES_ARRAY(NMESH_i,j)%nodes(2)%x
		NODES_ARRAY(NMESH_i+1,j)%y=MESHES_ARRAY(NMESH_i,j)%nodes(2)%y
	enddo
	do i=1,NMESH_i
		NODES_ARRAY(i,NMESH_j+1)%x=MESHES_ARRAY(i,NMESH_j)%nodes(4)%x
		NODES_ARRAY(i,NMESH_j+1)%y=MESHES_ARRAY(i,NMESH_j)%nodes(4)%y
	enddo
	NODES_ARRAY(NMESH_i+1,NMESH_j+1)%x=MESHES_ARRAY(NMESH_i,NMESH_j)%nodes(3)%x
	NODES_ARRAY(NMESH_i+1,NMESH_j+1)%y=MESHES_ARRAY(NMESH_i,NMESH_j)%nodes(3)%y

end subroutine fill_nodes_array


subroutine fill_mesh_nodes_from_nodes_array
  implicit none
  integer i,j
  do j=1, NMESH_j
     do i=1,NMESH_i
        MESHES_ARRAY(i,j)%nodes(1)=NODES_ARRAY(i,j)
        MESHES_ARRAY(i,j)%nodes(2)=NODES_ARRAY(i+1,j)
        MESHES_ARRAY(i,j)%nodes(3)=NODES_ARRAY(i+1,j+1)
        MESHES_ARRAY(i,j)%nodes(4)=NODES_ARRAY(i,j+1)
     enddo
  enddo

end subroutine fill_mesh_nodes_from_nodes_array


subroutine evaluate_pressure_bc_array
  implicit none
  integer             :: i,j
  integer             :: M1i,M1j,M2i,M2j,M3i,M3j,M4i,M4j,fn
  real,dimension(1:4) :: n,t         !!  distance of each points from normal (n), binormal (s) and tangential (t)
  real,dimension(1:2) :: n_vec, t_vec
  real                :: D
  real,dimension(0:9) :: a

	!$OMP     PARALLEL DO PRIVATE (i,j,M1i,M1j,M2i,M2j,M3i,M3j,M4i,M4j,fn,n,t,n_vec, t_vec,D,a)
	do i=1,tot_edges_bcs
		n_vec=edges_normal      (:,BCs_edges(i))
		t_vec=edges_tangential  (:,BCs_edges(i))
		M1i=BCs_edges_mesh_num(1,i)
		M1j=BCs_edges_mesh_num(2,i)
		BCs_edges_points(1,i)=BCs_edges(i)
		! edges on the left
		if(BCs_edges_position(i)==1) then
			fn=2
			M2i=M1i+1
			M2j=M1j
			M3i=M1i
			M3j=M1j+1
			if(M1j==NMESH_J) M3j=M1j-1
		elseif(BCs_edges_position(i)==2) then
			fn=1
			M2i=M1i-1
			M2j=M1j
			M3i=M1i
			M3j=M1j+1
			if(M1j==NMESH_J) M3j=M1j-1
		elseif(BCs_edges_position(i)==3) then
			fn=4
			M2j=M1j+1
			M2i=M1i
			M3j=M1j
			M3i=M1i+1
			if(M1i==NMESH_i) M3i=M1i-1
		elseif(BCs_edges_position(i)==4) then
			fn=3
			M2j=M1j-1
			M2i=M1i
			M3j=M1j
			M3i=M1i+1
			if(M1i==NMESH_i) M3i=M1i-1
		endif
		BCs_edges_points(2,i)=meshes_edges(fn,M1i,M1j)
		BCs_edges_points(3,i)=meshes_edges(fn,M2i,M2j)
		BCs_edges_points(4,i)=meshes_edges(fn,M3i,M3j) ! j+1 or j-1
		! find the distances
		do j=1,3
			n(j)=dot((edges_mid_coordinates(:,BCs_edges_points(j+1,i))-edges_mid_coordinates(:,BCs_edges_points(1,i))),n_vec)
			t(j)=dot((edges_mid_coordinates(:,BCs_edges_points(j+1,i))-edges_mid_coordinates(:,BCs_edges_points(1,i))),t_vec)
		enddo
		D=n(3)*(t(1)-t(2))+n(1)*(t(2)-t(3))+n(2)*(t(3)-t(1))
		a(1)=2*(n(3)*t(2)-n(2)*t(3))/(n(1)**2*D)
		a(2)=-2*(n(3)*t(1)-n(1)*t(3))/(n(2)**2*D)
		a(3)=2*(n(2)*t(1)-n(1)*t(2))/(n(3)**2*D)
		a(0)=-a(1)-a(2)-a(3)
		a(4)=-a(1)*n(1)-a(2)*n(2)-a(3)*n(3)
		a(5)=-a(1)*t(1)-a(2)*t(2)-a(3)*t(3)
		a(6)=-a(1)*n(1)*t(1)-a(2)*n(2)*t(2)-a(3)*n(3)*t(3)
		a(7)=-a(1)*t(1)**2./2.-a(2)*t(2)**2./2.-a(3)*t(3)**2./2.
		a(8)=-a(1)*n(1)*t(1)**2./2.-a(2)*n(2)*t(2)**2./2.-a(3)*n(3)*t(3)**2./2.
		a(9)=-a(1)*t(1)**3./6.-a(2)*t(2)**3./6.-a(3)*t(3)**3./6.
		BCs_edges_a(1:10,i)=a(0:9)
	enddo
	!$OMP     END PARALLEL DO

end subroutine evaluate_pressure_bc_array



subroutine fill_global_arrays
    integer i,j,k,l
	do j=1,NMESH_j
		do i=1,NMESH_i
			do k=1,4
				meshes_edges(k,i,j)=MESHES_ARRAY(i,j)%edges(k)
				meshes_nodes_x(k,i,j)=MESHES_ARRAY(i,j)%nodes(k)%x
				meshes_nodes_y(k,i,j)=MESHES_ARRAY(i,j)%nodes(k)%y
			enddo
		enddo
	enddo
! ================================	
! allocate BCs_edges
	l=tot_edges_bcs
	print*, "Total # of edges that have BCs = ", tot_edges_bcs
! Fill the arrays
l=0
! left
	do j=1,NMESH_j
		l=l+1
		BCs_edges(l)=EDGES_ARRAY(MESHES_ARRAY(1,j)%edges(1))%n
		BCs_edges_id(l)=1
		BCs_edges_position(l)=1
		BCs_edges_mesh_num(:,l)=(/1,j/)
	enddo
!! right
	do j=1,NMESH_j
		l=l+1
		BCs_edges(l)=EDGES_ARRAY(MESHES_ARRAY(NMESH_i,j)%edges(2))%n
		BCs_edges_id(l)=2
		BCs_edges_position(l)=2
		BCs_edges_mesh_num(:,l)=(/NMESH_i,j/)
	enddo
! bottom
	do i=1,NMESH_i
		l=l+1
		BCs_edges(l)=EDGES_ARRAY(MESHES_ARRAY(i,1)%edges(3))%n
		BCs_edges_id(l)=3
		BCs_edges_position(l)=3
		BCs_edges_mesh_num(:,l)=(/i,1/)
	enddo
! top
	do i=1,NMESH_i
		l=l+1
		BCs_edges(l)=EDGES_ARRAY(MESHES_ARRAY(i,NMESH_j)%edges(4))%n
		BCs_edges_id(l)=4
		BCs_edges_position(l)=4
		BCs_edges_mesh_num(:,l)=(/i,NMESH_j/)
	enddo

! ================================
end subroutine fill_global_arrays


end module geometry

