! NSE class -> solve NSE using NIM
module bc
  use general
  use nim
  use user_parameters
  use parameters
  use mesh
  use geometry
  implicit none

  contains


  subroutine apply_all_dirichlet_bcs
  ! ---------------------------------
    call apply_dirichlet_velocity_bcs
!	call apply_dirichlet_pressure_bcs

  end subroutine apply_all_dirichlet_bcs


  subroutine apply_dirichlet_velocity_bcs
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_edges_bcs
! Fill BCs for fields
	   edges_u(BCs_edges(i))=apply_dirichlet_bc_i(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),velx1,BCs_edges_id(i))
	   edges_v(BCs_edges(i))=apply_dirichlet_bc_i(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),velx2,BCs_edges_id(i))
    enddo
	!$OMP     END PARALLEL DO
  end subroutine apply_dirichlet_velocity_bcs


  subroutine apply_dirichlet_pressure_bcs
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_edges_bcs
! Fill BCs for pressure
	   edges_p(BCs_edges(i))=apply_dirichlet_bc_i(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),PPE,BCs_edges_id(i))
    enddo
	!$OMP     END PARALLEL DO
  end subroutine apply_dirichlet_pressure_bcs


  function apply_dirichlet_bc_i(Mi,Mj,n,field_id,bc_id) result(TT)
  implicit none
! -----------------------------------
    integer, intent (in)                          :: Mi,Mj
    integer, intent(in)                           :: n,field_id,bc_id
    real, dimension (1:Gauss_int_p)               :: TT_i, dS
	real      :: TT, area,x,y
    integer   :: i
    type(curvi2d)    :: C
! -----------------------------------
    do i=1,Gauss_int_p
		if(n==1 .or. n==2) then
			if (n==1) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),-1.0,Gauss_x_i(i),x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,-1.0,Gauss_x_i(i))
			elseif (n==2) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),1.0,Gauss_x_i(i),x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,1.0,Gauss_x_i(i))
		   endif
		   dS(i)=sqrt(C%g_ij(2,2))
		elseif(n==3 .or. n==4) then
			if (n==3) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),-1.0,x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,Gauss_x_i(i),-1.0)
			elseif (n==4) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),1.0,x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,Gauss_x_i(i),1.0)
		   endif
		   dS(i)=sqrt(C%g_ij(1,1))
		endif
		TT_i(i)=BCs_id(x,y,field_id,bc_id)*dS(i)
    enddo
	TT=Gauss_integration1d(TT_i)/edges_length(meshes_edges(n,Mi,Mj))
 end function apply_dirichlet_bc_i



  function find_bcs_velocity_derivatives(Mi,Mj,n,bc_id,dvdt_id) result(TT)
  implicit none
! -----------------------------------
    integer, intent (in)                          :: Mi,Mj
    integer, intent(in)                           :: n,dvdt_id,bc_id
    real, dimension (1:Gauss_int_p)               :: TT_i, dS
	real      :: TT, area,x,y
    integer   :: i
    type(curvi2d)    :: C
! -----------------------------------
    do i=1,Gauss_int_p
		if(n==1 .or. n==2) then
			if (n==1) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),-1.0,Gauss_x_i(i),x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,-1.0,Gauss_x_i(i))
			elseif (n==2) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),1.0,Gauss_x_i(i),x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,1.0,Gauss_x_i(i))
		   endif
		   dS(i)=sqrt(C%g_ij(2,2))
		elseif(n==3 .or. n==4) then
			if (n==3) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),-1.0,x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,Gauss_x_i(i),-1.0)
			elseif (n==4) then
				call map_xy_2d(MESHES_ARRAY(Mi,Mj),Gauss_x_i(i),1.0,x,y)
				call evaluate_curvilinear(MESHES_ARRAY(Mi,Mj),C,Gauss_x_i(i),1.0)
		   endif
		   dS(i)=sqrt(C%g_ij(1,1))
		endif
		TT_i(i)=BCs_derivatives(x,y,bc_id,dvdt_id)*dS(i)
    enddo
	TT=Gauss_integration1d(TT_i)/edges_length(meshes_edges(n,Mi,Mj))
 end function find_bcs_velocity_derivatives


  subroutine evaluate_velocity_derivatives
    implicit none
    integer              :: i
	!$OMP     PARALLEL DO PRIVATE (i)
    do i=1,tot_edges_bcs
	! Evaluate the velocity derivatives on the boundaries
		BCs_edges_v_t  (i)=find_bcs_velocity_derivatives(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),BCs_edges_id(i),bc_dvdt)
		BCs_edges_v_tt (i)=find_bcs_velocity_derivatives(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),BCs_edges_id(i),bc_dvdtt)
		BCs_edges_v_ttt(i)=find_bcs_velocity_derivatives(BCs_edges_mesh_num(1,i),BCs_edges_mesh_num(2,i),BCs_edges_position(i),BCs_edges_id(i),bc_dvdttt)
    enddo
	!$OMP     END PARALLEL DO
  end subroutine evaluate_velocity_derivatives


subroutine apply_pressure_bc
   integer                 :: i,j,M1i,M1j,fn
   real,dimension(1:2)     :: v,b,n
   real,dimension(1:4)     :: v_n		! normal velocities for the 5 points
   real                    :: b_n,d2udn2, B1,B2		! B1 is the pseudo source terms. B2,B3 are the second derivetives
   integer,dimension(1:4)  :: n1_edges					! edges of the mesh
	!$OMP     PARALLEL DO PRIVATE (i,j,n,v,b,b_n,v_n,d2udn2,M1i,M1j,fn,B1,B2,n1_edges)
    do i=1,tot_edges_bcs
		! Get the mesh where the face exists
		M1i=BCs_edges_mesh_num(1,i)
		M1j=BCs_edges_mesh_num(2,i)
		! Get the edges of the mesh
		n1_edges=meshes_edges(:,M1i,M1j)
		! Get the normal to the face
		n=edges_normal(:,BCs_edges(i))
		! Get the b-normal
		b=(/edges_bx(BCs_edges(i)),edges_by(BCs_edges(i))/)
		b_n=dot(b,n)
		! Get the normal velocities
		do j=1,4
			v=(/edges_u(BCs_edges_points(j,i)),edges_v(BCs_edges_points(j,i))/)
			v_n(j)=dot(v,n)
		enddo
		 d2udn2=-BCs_edges_a(1,i)*v_n(1)-BCs_edges_a(2,i)*v_n(2)-BCs_edges_a(3,i)*v_n(3)-BCs_edges_a(4,i)*v_n(4) & 
		        +(BCs_edges_a(5,i))*BCs_edges_v_t(i)+(BCs_edges_a(7,i))*BCs_edges_v_tt(i) &
			    +(BCs_edges_a(9,i))*BCs_edges_v_ttt(i)
	! ---------> The integration part <----------
		if(BCs_edges_position(i)==1) then
			fn=1
			if(M1j==1) then
				B2=(-3.*edges_p(meshes_edges(fn,M1i,M1j))+4.*edges_p(meshes_edges(fn,M1i,M1j+1))-edges_p(meshes_edges(fn,M1i,M1j+2)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			elseif(M1j==NMESH_j) then
				B2=( 3.*edges_p(meshes_edges(fn,M1i,M1j))-4.*edges_p(meshes_edges(fn,M1i,M1j-1))+edges_p(meshes_edges(fn,M1i,M1j-2)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			else
				B2=(edges_p(meshes_edges(fn,M1i,M1j+1))-edges_p(meshes_edges(fn,M1i,M1j-1)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			endif
			B1=BCs_edges_ap(6,i)*(meshes_R_D_p(M1i,M1j)+meshes_R_B_p(M1i,M1j))
			edges_p(n1_edges(1))=(1-w_p2)*edges_p(n1_edges(1)) + w_p2/(BCs_edges_ap(1,i)) * ( &
								 -BCs_edges_ap(2,i)*edges_p(n1_edges(2)) &
								 -BCs_edges_ap(3,i)*edges_p(n1_edges(3)) &
								 -BCs_edges_ap(4,i)*edges_p(n1_edges(4)) &
								 +B1-B2+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_edges_position(i)==2) then
			fn=2
			if(M1j==1) then
				B2=(-3.*edges_p(meshes_edges(fn,M1i,M1j))+4.*edges_p(meshes_edges(fn,M1i,M1j+1))-edges_p(meshes_edges(fn,M1i,M1j+2)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			elseif(M1j==NMESH_j) then
				B2=( 3.*edges_p(meshes_edges(fn,M1i,M1j))-4.*edges_p(meshes_edges(fn,M1i,M1j-1))+edges_p(meshes_edges(fn,M1i,M1j-2)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			else
				B2=(edges_p(meshes_edges(fn,M1i,M1j+1))-edges_p(meshes_edges(fn,M1i,M1j-1)))/4.0*meshes_gij_edges(1,2,fn,M1i,M1j)
			endif
			B1=BCs_edges_ap(6,i)*(meshes_R_D_p(M1i,M1j)+meshes_R_B_p(M1i,M1j))
			edges_p(n1_edges(2))=(1-w_p2)*edges_p(n1_edges(2)) + w_p2/(BCs_edges_ap(2,i)) * ( &
								 -BCs_edges_ap(1,i)*edges_p(n1_edges(1)) &
								 -BCs_edges_ap(3,i)*edges_p(n1_edges(3)) &
								 -BCs_edges_ap(4,i)*edges_p(n1_edges(4)) &
								 +B1-B2+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_edges_position(i)==3) then
			fn=3
			if(M1i==1) then
				B2=(-3.*edges_p(meshes_edges(fn,M1i,M1j))+4.*edges_p(meshes_edges(fn,M1i+1,M1j))-edges_p(meshes_edges(fn,M1i+2,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			elseif(M1i==NMESH_i) then
				B2=( 3.*edges_p(meshes_edges(fn,M1i,M1j))-4.*edges_p(meshes_edges(fn,M1i-1,M1j))+edges_p(meshes_edges(fn,M1i-2,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			else
				B2=(edges_p(meshes_edges(fn,M1i+1,M1j))-edges_p(meshes_edges(fn,M1i-1,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			endif
			B1=BCs_edges_ap(6,i)*(meshes_R_D_p(M1i,M1j)+meshes_R_B_p(M1i,M1j))
			edges_p(n1_edges(3))=(1-w_p2)*edges_p(n1_edges(3)) + w_p2/(BCs_edges_ap(3,i)) * ( &
								 -BCs_edges_ap(1,i)*edges_p(n1_edges(1)) &
								 -BCs_edges_ap(2,i)*edges_p(n1_edges(2)) &
								 -BCs_edges_ap(4,i)*edges_p(n1_edges(4)) &
								 +B1-B2+(nu*rho*d2udn2-b_n*rho))
		elseif(BCs_edges_position(i)==4) then
			fn=4
			if(M1i==1) then
				B2=(-3.*edges_p(meshes_edges(fn,M1i,M1j))+4.*edges_p(meshes_edges(fn,M1i+1,M1j))-edges_p(meshes_edges(fn,M1i+2,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			elseif(M1i==NMESH_i) then
				B2=( 3.*edges_p(meshes_edges(fn,M1i,M1j))-4.*edges_p(meshes_edges(fn,M1i-1,M1j))+edges_p(meshes_edges(fn,M1i-2,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			else
				B2=(edges_p(meshes_edges(fn,M1i+1,M1j))-edges_p(meshes_edges(fn,M1i-1,M1j)))/4.0*meshes_gij_edges(2,1,fn,M1i,M1j)
			endif
			B1=BCs_edges_ap(6,i)*(meshes_R_D_p(M1i,M1j)+meshes_R_B_p(M1i,M1j))
			edges_p(n1_edges(4))=(1-w_p2)*edges_p(n1_edges(4)) + w_p2/(BCs_edges_ap(4,i)) * ( &
								 -BCs_edges_ap(1,i)*edges_p(n1_edges(1)) &
								 -BCs_edges_ap(2,i)*edges_p(n1_edges(2)) &
								 -BCs_edges_ap(3,i)*edges_p(n1_edges(3)) &
								 +B1-B2+(nu*rho*d2udn2-b_n*rho))	
		endif
    enddo
	!$OMP     END PARALLEL DO
end subroutine apply_pressure_bc

! Apply zero gradient velocity on the boundaries
subroutine apply_zero_gradient_bc
	implicit none
	integer  :: i,M1i,M1j
	!$OMP     PARALLEL DO PRIVATE (i,M1i,M1j)
    do i=1,tot_edges_bcs
		if(BCs_edges_position(i)==4) then
			M1i=BCs_edges_mesh_num(1,i)
			M1j=BCs_edges_mesh_num(2,i)
			edges_u(meshes_edges(4,M1i,M1j))=edges_u(meshes_edges(3,M1i,M1j))
			edges_v(meshes_edges(4,M1i,M1j))=edges_v(meshes_edges(3,M1i,M1j))
		
		endif
	enddo
end subroutine apply_zero_gradient_bc


end module bc

