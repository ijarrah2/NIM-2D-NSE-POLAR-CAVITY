module nim
use user_parameters
use general
use parameters
use curvilinear
use mesh
use geometry
implicit none
contains

! ++++++++++++++++++++++++++++ Velocity +++++++++++++++++++++++++++++
! *****************  Calculating average constants ******************
  subroutine calculate_velocity_nim_constants
    implicit none
    integer :: i,i2,j,j2
	real    :: B_x,B_y,a,b,Re_x,Re_y,K11,K22
	real,dimension(1:NMESH_i,1:NMESH_j)    :: Det
	real,dimension(1:2,1:NMESH_i,1:NMESH_j):: CC,EE,FF,GG,HH,II,JJ
    a=1.
    b=1.
	!$OMP     PARALLEL DO PRIVATE (i,j,Re_x,Re_y,B_x,B_y,K11,K22)
	do j=1,NMESH_j
		do i=1,NMESH_i
			meshes_u(i,j)=meshes_uxyt(i,j)*meshes_xi_x(1,1,i,j)+meshes_vxyt(i,j)*meshes_xi_x(1,2,i,j)
			meshes_v(i,j)=meshes_uxyt(i,j)*meshes_xi_x(2,1,i,j)+meshes_vxyt(i,j)*meshes_xi_x(2,2,i,j)
			K11=meshes_gij(1,1,i,j)*nu
			K22=meshes_gij(2,2,i,j)*nu
			
			Re_x=2*meshes_u(i,j)/K11
			Re_y=2*meshes_v(i,j)/K22

			B_x=B_nim(meshes_u(i,j),K11, Re_x)
			B_y=B_nim(meshes_v(i,j),K22, Re_y)

			CC(1,i,j)=C_nim(meshes_u(i,j),K11,Re_x,B_x,a)
			CC(2,i,j)=C_nim(meshes_v(i,j),K22,Re_y,B_y,b)

			Det(i,j)=1./(CC(1,i,j)+CC(2,i,j))

			EE(1,i,j)=E_nim(meshes_u(i,j),K11,Re_x,B_x,Det(i,j),a)
			EE(2,i,j)=E_nim(meshes_v(i,j),K22,Re_y,B_y,Det(i,j),b)

			FF(1,i,j)=F_nim(meshes_u(i,j),K11,Re_x,B_x,Det(i,j),a)
			FF(2,i,j)=F_nim(meshes_v(i,j),K22,Re_y,B_y,Det(i,j),b)

			GG(1,i,j)=G_nim(meshes_u(i,j),K11,Re_x,a)
			GG(2,i,j)=G_nim(meshes_v(i,j),K22,Re_y,b)

			HH(1,i,j)=H_nim(meshes_u(i,j),K11,Re_x,a)
			HH(2,i,j)=H_nim(meshes_v(i,j),K22,Re_y,b)

			II(1,i,j)=I_nim(meshes_u(i,j),K11,Re_x,a)
			II(2,i,j)=I_nim(meshes_v(i,j),K22,Re_y,b)

			JJ(1,i,j)=J_nim(meshes_u(i,j),K11,Re_x,a)
			JJ(2,i,j)=J_nim(meshes_v(i,j),K22,Re_y,b)
		enddo
    enddo
  !$OMP     END PARALLEL DO 

!$OMP     PARALLEL DO PRIVATE (i,j,i2,j2)
	do j=1,NMESH_j
		do i=1,NMESH_i
	! ==================== |F|A|C|E| |2| =====================
		if (i<NMESH_i) THEN
		   i2=i+1
		   AA_vx(1,i,j)=-(GG(1,i ,j)+EE(1,i,j)*HH(1,i,j))*meshes_gij_edges(1,1,2,i,j)
		   AA_vx(2,i,j)= (GG(1,i ,j)-FF(1,i,j)*HH(1,i,j))*meshes_gij_edges(1,1,2,i,j) &
						+(II(1,i2,j)+EE(1,i2,j)*JJ(1,i2,j))*meshes_gij_edges(1,1,1,i2,j)
		   AA_vx(3,i,j)=-(II(1,i2,j)-FF(1,i2,j)*JJ(1,i2,j))*meshes_gij_edges(1,1,1,i2,j)
		   AA_vx(4,i,j)=  EE(2,i ,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_vx(5,i,j)=  FF(2,i ,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_vx(6,i,j)= -EE(2,i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		   AA_vx(7,i,j)= -FF(2,i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		   AA_vx(8,i,j)= CC(2,i ,j)*Det(i,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)/(2*DT)
		   AA_vx(9,i,j)=-CC(2,i ,j)*Det(i,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)/(2*DT)
		   AA_vx(10,i,j)=-CC(2,i2,j)*Det(i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)/(2*DT)
		   AA_vx(11,i,j)= CC(2,i2,j)*Det(i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)/(2*DT)
		endif

	! ==================== |F|A|C|E| |4| =====================
		if (j<NMESH_j) THEN
		   j2=j+1
		   AA_vy(1,i,j)=  EE(1,i, j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_vy(2,i,j)=  FF(1,i,j )*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_vy(3,i,j)= -EE(1,i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		   AA_vy(4,i,j)= -FF(1,i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		   AA_vy(5,i,j)=-(GG(2,i,j )+EE(2,i,j)*HH(2,i,j))*meshes_gij_edges(2,2,4,i,j)
		   AA_vy(6,i,j)= (GG(2,i,j )-FF(2,i,j)*HH(2,i,j))*meshes_gij_edges(2,2,4,i,j) &
						+(II(2,i,j2)+EE(2,i,j2)*JJ(2,i,j2))*meshes_gij_edges(2,2,3,i,j2)
		   AA_vy(7,i,j)=-(II(2,i,j2)-FF(2,i,j2)*JJ(2,i,j2))*meshes_gij_edges(2,2,3,i,j2)
		   AA_vy(8,i,j)= CC(1,i,j )*Det(i,j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)/(2*DT)
		   AA_vy(9,i,j)=-CC(1,i,j )*Det(i,j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)/(2*DT)
		   AA_vy(10,i,j)=-CC(1,i,j2)*Det(i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)/(2*DT)
		   AA_vy(11,i,j)= CC(1,i,j2)*Det(i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)/(2*DT)
		endif
		   a_t(1,i,j)=2*EE(1,i,j)*CC(2,i,j)
		   a_t(2,i,j)=2*FF(1,i,j)*CC(2,i,j)
		   a_t(3,i,j)=2*EE(2,i,j)*CC(1,i,j)
		   a_t(4,i,j)=2*FF(2,i,j)*CC(1,i,j)
		   a_t(5,i,j)=2*CC(1,i,j)*CC(2,i,j)*Det(i,j)
		   a_t(6,i,j)=-a_t(5,i,j)/(2*DT)-1
		   a_t(7,i,j)= a_t(5,i,j)/(2*DT)-1
		enddo
    enddo
  !$OMP     END PARALLEL DO 

  end subroutine calculate_velocity_nim_constants



  subroutine calculate_velocity_nim_constants_formulation2
! -----------------------------------
    integer                         :: i,j
	real							::Re_x,Re_y,u,v,a,b,k11,k22
	!$OMP     PARALLEL DO PRIVATE (i,j,Re_x,Re_y,u,v,a,b,k11,k22)
	do j=1,NMESH_j
		do i=1,NMESH_i
		
			a=1.
			b=1.
			meshes_u(i,j)=meshes_uxyt(i,j)*meshes_xi_x(1,1,i,j)+meshes_vxyt(i,j)*meshes_xi_x(1,2,i,j)
			meshes_v(i,j)=meshes_uxyt(i,j)*meshes_xi_x(2,1,i,j)+meshes_vxyt(i,j)*meshes_xi_x(2,2,i,j)
			u=meshes_u(i,j)
			v=meshes_v(i,j)
			k11=meshes_gij(1,1,i,j)*nu
			k22=meshes_gij(2,2,i,j)*nu

			Re_x=2*a*u/K11
			Re_y=2*b*v/K22

			if(abs(Re_x)>v_limit) then
			  M_A1(i,j)=(u/(k11*(1./exp(Re_x)-1.)))*meshes_gij_edges(1,1,2,i,j)
			  M_A2(i,j)=(2*a/(k11*(1./exp(Re_x)-1.))+1./u)*meshes_gij_edges(1,1,2,i,j)
			  M_A3(i,j)=(u/(k11*(1.-exp(Re_x))))*meshes_gij_edges(1,1,1,i,j)
			  M_A4(i,j)=(-2*a/(k11*(1.-exp(Re_x)))-1./u)*meshes_gij_edges(1,1,1,i,j)

			  if(Re_x>0) M_C1(i,j)=(-k11+(a*u*(1./exp(Re_x)+1.)/(-1./exp(Re_x)+1.)))/(u**2)
			  if(Re_x<0) M_C1(i,j)=(-k11+(a*u*(1.+exp(Re_x))/(-1.+exp(Re_x))))/(u**2)

			  if(Re_x>0) M_C2(i,j)=1./(-1./exp(Re_x)+1.)-1./(Re_x)
			  if(Re_x<0) M_C2(i,j)=exp(Re_x)/(-1.+exp(Re_x))-1./(Re_x)
			else
			  M_A1(i,j)= (-1/(2.*a) - u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
			  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_edges(1,1,2,i,j)

			  M_A2(i,j)= (-(a/k11) - (a**2*u)/(3.*k11**2) + (a**4*u**3)/(45.*k11**4) - (2*a**6*u**5)/(945.*k11**6) + &
			  (a**8*u**7)/(4725.*k11**8) - (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_edges(1,1,2,i,j)

			  M_A4(i,j)=(-(a/k11) + (a**2*u)/(3.*k11**2) - (a**4*u**3)/(45.*k11**4) + (2*a**6*u**5)/(945.*k11**6) - &
			  (a**8*u**7)/(4725.*k11**8) + (2*a**10*u**9)/(93555.*k11**10))*meshes_gij_edges(1,1,1,i,j)

			  M_A3(i,j)=(-1/(2.*a) + u/(2.*k11) - (a*u**2)/(6.*k11**2) + (a**3*u**4)/(90.*k11**4) - &
			  (a**5*u**6)/(945.*k11**6) + (a**7*u**8)/(9450.*k11**8) - (a**9*u**10)/(93555.*k11**10))*meshes_gij_edges(1,1,1,i,j)

			  M_C1(i,j)= a**2/(3.*k11) - (a**4*u**2)/(45.*k11**3) + (2*a**6*u**4)/(945.*k11**5) - &
			  (a**8*u**6)/(4725.*k11**7) + (2*a**10*u**8)/(93555.*k11**9) - &
			  (1382*a**12*u**10)/(6.38512875e8*k11**11)

			  M_C2(i,j)= 0.5 + (a*u)/(6.*k11) - (a**3*u**3)/(90.*k11**3) + (a**5*u**5)/(945.*k11**5) - &
			  (a**7*u**7)/(9450.*k11**7) + (a**9*u**9)/(93555.*k11**9)

			endif

			if(abs(Re_y)>v_limit) then
			  M_B1(i,j)=(v/(k22*(1./exp(Re_y)-1.)))*meshes_gij_edges(2,2,4,i,j)
			  M_B2(i,j)=(2*b/(k22*(1./exp(Re_y)-1.))+1./v)*meshes_gij_edges(2,2,4,i,j)
			  M_B3(i,j)=(v/(k22*(1-exp(Re_y))))*meshes_gij_edges(2,2,3,i,j)
			  M_B4(i,j)=(-2*b/(k22*(1-exp(Re_y)))-1./v)*meshes_gij_edges(2,2,3,i,j)

			  if(Re_y>0) M_D1(i,j)=(-k22+(b*v*(1./exp(Re_y)+1.)/(-1./exp(Re_y)+1.)))/(v**2)
			  if(Re_y<0) M_D1(i,j)=(-k22+(b*v*(1+exp(Re_y))/(-1.+exp(Re_y))))/(v**2)

			  if(Re_y>0) M_D2(i,j)=1./(-1./exp(Re_y)+1.)-1./(Re_y)
			  if(Re_y<0) M_D2(i,j)=exp(Re_y)/(-1+exp(Re_y))-1./(Re_y)

			else
			  M_B1(i,j)= (-1/(2.*b) - v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
			  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_edges(2,2,4,i,j)

			  M_B2(i,j)= (-(b/k22) - (b**2*v)/(3.*k22**2) + (b**4*v**3)/(45.*k22**4) - (2*b**6*v**5)/(945.*k22**6) + &
			  (b**8*v**7)/(4725.*k22**8) - (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_edges(2,2,4,i,j)

			  M_B4(i,j)=(-(b/k22) + (b**2*v)/(3.*k22**2) - (b**4*v**3)/(45.*k22**4) + (2*b**6*v**5)/(945.*k22**6) - &
			  (b**8*v**7)/(4725.*k22**8) + (2*b**10*v**9)/(93555.*k22**10))*meshes_gij_edges(2,2,3,i,j)

			  M_B3(i,j)=(-1/(2.*b) + v/(2.*k22) - (b*v**2)/(6.*k22**2) + (b**3*v**4)/(90.*k22**4) - &
			  (b**5*v**6)/(945.*k22**6) + (b**7*v**8)/(9450.*k22**8) - (b**9*v**10)/(93555.*k22**10))*meshes_gij_edges(2,2,3,i,j)

			  M_D1(i,j)= b**2/(3.*k22) - (b**4*v**2)/(45.*k22**3) + (2*b**6*v**4)/(945.*k22**5) - &
			  (b**8*v**6)/(4725.*k22**7) + (2*b**10*v**8)/(93555.*k22**9) - &
			  (1382*b**12*v**10)/(6.38512875e8*k22**11)

			  M_D2(i,j)= 0.5 + (b*v)/(6.*k22) - (b**3*v**3)/(90.*k22**3) + (b**5*v**5)/(945.*k22**5) - &
			  (b**7*v**7)/(9450.*k22**7) + (b**9*v**9)/(93555.*k22**9)

			endif
		enddo
    enddo
  !$OMP     END PARALLEL DO 

  !$OMP     PARALLEL DO PRIVATE (i,j)
	do j=1,NMESH_j
		do i=1,NMESH_i
			M_F32(i,j)=0.5*(1./M_C1(i,j)+1./M_D1(i,j)-1./DT)
			M_F33(i,j)=-M_C2(i,j)/M_C1(i,j)
			M_F34(i,j)=(M_C2(i,j)-1.)/M_C1(i,j)
			M_F35(i,j)=-M_D2(i,j)/M_D1(i,j)
			M_F36(i,j)=(M_D2(i,j)-1.)/M_D1(i,j)
			M_F31(i,j)=M_F32(i,j)+M_F33(i,j)+M_F34(i,j)+M_F35(i,j)+M_F36(i,j)

		   if (i<NMESH_i) then
				 M_F12(i,j)=M_A1(i,j)-(M_A2(i,j)*M_C2(i,j)/M_C1(i,j))
				 M_F13(i,j)=M_A3(i+1,j)-M_A4(i+1,j)/M_C1(i+1,j)+M_A4(i+1,j)*M_C2(i+1,j)/M_C1(i+1,j)
				 M_F14(i,j)=M_A2(i,j)/(2.*M_C1(i,j))
				 M_F15(i,j)=M_A4(i+1,j)/(2.*M_C1(i+1,j))
				 M_F11(i,j)=M_F12(i,j)+M_F13(i,j)+2.*M_F14(i,j)+2.*M_F15(i,j)
			endif
			if (j<NMESH_j) then
				 M_F22(i,j)=M_B1(i,j)-(M_B2(i,j)*M_D2(i,j)/M_D1(i,j))
				 M_F23(i,j)=M_B3(i,j+1)-M_B4(i,j+1)/M_D1(i,j+1)+M_B4(i,j+1)*M_D2(i,j+1)/M_D1(i,j+1)
				 M_F24(i,j)=M_B2(i,j)/(2.*M_D1(i,j))
				 M_F25(i,j)=M_B4(i,j+1)/(2.*M_D1(i,j+1))
				 M_F21(i,j)=M_F22(i,j)+M_F23(i,j)+2.*M_F24(i,j)+2.*M_F25(i,j)
			endif
		enddo
    enddo
  !$OMP     END PARALLEL DO 
  end subroutine calculate_velocity_nim_constants_formulation2


! ++++++++++++++++++++++++++++ Velocity +++++++++++++++++++++++++++++
! *****************  Calculating average constants ******************

  subroutine calculate_pressure_nim_constants
    implicit none
    integer :: i,i2,j,j2,Mi,Mj
	real    :: B_x,B_y,a,b,Re_x,Re_y,K11,K22
	real,dimension(1:NMESH_i,1:NMESH_j)    :: Det
	real,dimension(1:2,1:NMESH_i,1:NMESH_j):: CC,EE,FF,GG,HH,II,JJ
    a=1.
    b=1.
	!$OMP     PARALLEL DO PRIVATE (i,j,Re_x,Re_y,B_x,B_y,K11,K22)
	do j=1,NMESH_j
		do i=1,NMESH_i

			K11=meshes_gij(1,1,i,j)/rho
			K22=meshes_gij(2,2,i,j)/rho
			Re_x=0.
			Re_y=0.

			B_x=B_nim(0.,K11, Re_x)
			B_y=B_nim(0.,K22, Re_y)

			CC(1,i,j)=C_nim(0.,K11,Re_x,B_x,a)
			CC(2,i,j)=C_nim(0.,K22,Re_y,B_y,b)

			Det(i,j)=1./(CC(1,i,j)+CC(2,i,j))

			EE(1,i,j)=E_nim(0.,K11,Re_x,B_x,Det(i,j),a)
			EE(2,i,j)=E_nim(0.,K22,Re_y,B_y,Det(i,j),b)

			FF(1,i,j)=F_nim(0.,K11,Re_x,B_x,Det(i,j),a)
			FF(2,i,j)=F_nim(0.,K22,Re_y,B_y,Det(i,j),b)

			GG(1,i,j)=G_nim(0.,K11,Re_x,a)
			GG(2,i,j)=G_nim(0.,K22,Re_y,b)

			HH(1,i,j)=H_nim(0.,K11,Re_x,a)
			HH(2,i,j)=H_nim(0.,K22,Re_y,b)

			II(1,i,j)=I_nim(0.,K11,Re_x,a)
			II(2,i,j)=I_nim(0.,K22,Re_y,b)

			JJ(1,i,j)=J_nim(0.,K11,Re_x,a)
			JJ(2,i,j)=J_nim(0.,K22,Re_y,b)
		enddo
    enddo
  !$OMP     END PARALLEL DO 

  !$OMP     PARALLEL DO PRIVATE (i,j,i2,j2)
	do j=1,NMESH_j
		do i=1,NMESH_i
	! ==================== |F|A|C|E| |2| =====================
		if (i<NMESH_i) THEN
		   i2=i+1
		   AA_px(1,i,j)=-(GG(1,i ,j)+EE(1,i,j)*HH(1,i,j))*meshes_gij_edges(1,1,2,i,j)
		   AA_px(2,i,j)= (GG(1,i ,j)-FF(1,i,j)*HH(1,i,j))*meshes_gij_edges(1,1,2,i,j) &
						+(II(1,i2,j)+EE(1,i2,j)*JJ(1,i2,j))*meshes_gij_edges(1,1,1,i2,j)
		   AA_px(3,i,j)=-(II(1,i2,j)-FF(1,i2,j)*JJ(1,i2,j))*meshes_gij_edges(1,1,1,i2,j)
		   AA_px(4,i,j)=  EE(2,i ,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_px(5,i,j)=  FF(2,i ,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_px(6,i,j)= -EE(2,i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		   AA_px(7,i,j)= -FF(2,i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		   AA_px(8,i,j)= CC(2,i,j)*Det(i,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_px(9,i,j)=-CC(2,i,j)*Det(i,j)*HH(1,i,j)*meshes_gij_edges(1,1,2,i,j)
		   AA_px(10,i,j)=-CC(2,i2,j)*Det(i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		   AA_px(11,i,j)= CC(2,i2,j)*Det(i2,j)*JJ(1,i2,j)*meshes_gij_edges(1,1,1,i2,j)
		endif

	! ==================== |F|A|C|E| |4| =====================
		if (j<NMESH_j) THEN
		   j2=j+1
		   AA_py(1,i,j)=  EE(1,i, j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_py(2,i,j)=  FF(1,i,j )*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_py(3,i,j)= -EE(1,i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		   AA_py(4,i,j)= -FF(1,i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		   AA_py(5,i,j)=-(GG(2,i,j )+EE(2,i,j)*HH(2,i,j))*meshes_gij_edges(2,2,4,i,j)
		   AA_py(6,i,j)= (GG(2,i,j )-FF(2,i,j)*HH(2,i,j))*meshes_gij_edges(2,2,4,i,j) &
						+(II(2,i,j2)+EE(2,i,j2)*JJ(2,i,j2))*meshes_gij_edges(2,2,3,i,j2)
		   AA_py(7,i,j)=-(II(2,i,j2)-FF(2,i,j2)*JJ(2,i,j2))*meshes_gij_edges(2,2,3,i,j2)
		   AA_py(8,i,j)= CC(1,i,j)*Det(i,j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_py(9,i,j)=-CC(1,i,j)*Det(i,j)*HH(2,i,j)*meshes_gij_edges(2,2,4,i,j)
		   AA_py(10,i,j)=-CC(1,i,j2)*Det(i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		   AA_py(11,i,j)= CC(1,i,j2)*Det(i,j2)*JJ(2,i,j2)*meshes_gij_edges(2,2,3,i,j2)
		endif
		enddo
    enddo
  !$OMP     END PARALLEL DO 

!$OMP     PARALLEL DO PRIVATE (i,Mi,Mj)
    do i=1,tot_edges_bcs
! ==================== |F|A|C|E| |2| =====================
		Mi=BCs_edges_mesh_num(1,i)
		Mj=BCs_edges_mesh_num(2,i)
		if (BCs_edges_position(i)==1) THEN
			BCs_edges_ap(1,i)=-(II(1,Mi,Mj)+EE(1,Mi,Mj)*JJ(1,Mi,Mj))*meshes_gij_edges(1,1,1,Mi,Mj)
			BCs_edges_ap(2,i)= (II(1,Mi,Mj)-FF(1,Mi,Mj)*JJ(1,Mi,Mj))*meshes_gij_edges(1,1,1,Mi,Mj)
			BCs_edges_ap(3,i)= EE(2,Mi,Mj)*JJ(1,Mi,Mj)*meshes_gij_edges(1,1,1,Mi,Mj)
			BCs_edges_ap(4,i)= FF(2,Mi,Mj)*JJ(1,Mi,Mj)*meshes_gij_edges(1,1,1,Mi,Mj)
			BCs_edges_ap(5,i)= CC(2,Mi,Mj)*Det(Mi,Mj)*JJ(1,Mi,Mj)*meshes_gij_edges(1,1,1,Mi,Mj)
			BCs_edges_ap(6,i)=-CC(2,Mi,Mj)*Det(Mi,Mj)*JJ(1,Mi,Mj)*meshes_gij_edges(1,1,1,Mi,Mj)
		elseif(BCs_edges_position(i)==2) THEN
			BCs_edges_ap(1,i)=-(GG(1,Mi,Mj)+EE(1,Mi,Mj)*HH(1,Mi,Mj))*meshes_gij_edges(1,1,2,Mi,Mj)
			BCs_edges_ap(2,i)= (GG(1,Mi,Mj)-FF(1,Mi,Mj)*HH(1,Mi,Mj))*meshes_gij_edges(1,1,2,Mi,Mj)
			BCs_edges_ap(3,i)= EE(2,Mi,Mj)*HH(1,Mi,Mj)*meshes_gij_edges(1,1,2,Mi,Mj)
			BCs_edges_ap(4,i)= FF(2,Mi,Mj)*HH(1,Mi,Mj)*meshes_gij_edges(1,1,2,Mi,Mj)
			BCs_edges_ap(5,i)= CC(2,Mi,Mj)*Det(Mi,Mj)*HH(1,Mi,Mj)*meshes_gij_edges(1,1,2,Mi,Mj)
			BCs_edges_ap(6,i)=-CC(2,Mi,Mj)*Det(Mi,Mj)*HH(1,Mi,Mj)*meshes_gij_edges(1,1,2,Mi,Mj)
		elseif(BCs_edges_position(i)==3) THEN
			BCs_edges_ap(1,i)= EE(1,Mi,Mj)*JJ(2,Mi,Mj)*meshes_gij_edges(2,2,3,Mi,Mj)
			BCs_edges_ap(2,i)= FF(1,Mi,Mj)*JJ(2,Mi,Mj)*meshes_gij_edges(2,2,3,Mi,Mj)
			BCs_edges_ap(3,i)=-(II(2,Mi,Mj)+EE(2,Mi,Mj)*JJ(2,Mi,Mj))*meshes_gij_edges(2,2,3,Mi,Mj)
			BCs_edges_ap(4,i)= (II(2,Mi,Mj)-FF(2,Mi,Mj)*JJ(2,Mi,Mj))*meshes_gij_edges(2,2,3,Mi,Mj)
			BCs_edges_ap(5,i)= CC(1,Mi,Mj)*Det(Mi,Mj)*JJ(2,Mi,Mj)*meshes_gij_edges(2,2,3,Mi,Mj)
			BCs_edges_ap(6,i)=-CC(1,Mi,Mj)*Det(Mi,Mj)*JJ(2,Mi,Mj)*meshes_gij_edges(2,2,3,Mi,Mj)
		elseif(BCs_edges_position(i)==4) THEN
			BCs_edges_ap(1,i)= EE(1,Mi,Mj)*HH(2,Mi,Mj)*meshes_gij_edges(2,2,4,Mi,Mj)
			BCs_edges_ap(2,i)= FF(1,Mi,Mj)*HH(2,Mi,Mj)*meshes_gij_edges(2,2,4,Mi,Mj)
			BCs_edges_ap(3,i)=-(GG(2,Mi,Mj)+EE(2,Mi,Mj)*HH(2,Mi,Mj))*meshes_gij_edges(2,2,4,Mi,Mj)
			BCs_edges_ap(4,i)= (GG(2,Mi,Mj)-FF(2,Mi,Mj)*HH(2,Mi,Mj))*meshes_gij_edges(2,2,4,Mi,Mj)
			BCs_edges_ap(5,i)= CC(1,Mi,Mj)*Det(Mi,Mj)*HH(2,Mi,Mj)*meshes_gij_edges(2,2,4,Mi,Mj)
			BCs_edges_ap(6,i)=-CC(1,Mi,Mj)*Det(Mi,Mj)*HH(2,Mi,Mj)*meshes_gij_edges(2,2,4,Mi,Mj)
		endif
    enddo
  !$OMP     END PARALLEL DO 
  end subroutine calculate_pressure_nim_constants


! ========================================= NIM constants =========================================
  function condition(u,F,Re) result(cond)
     real, intent(in)        :: Re,u,F
     integer                    :: cond
    if (abs(Re)<v_limit) then     ! 0.01
       cond=1
    else
       cond=3
    endif
  end function condition

  function B_nim(u,F,Re)  result(P)
     real, intent(in)        :: Re,u,F
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=0         ! For this case, it is undefined.
       case (2)
            P=0
       case default
            P=1./u*(1./Re - 1./(exp(Re)-1))
    end select
  end function B_nim

  function C_nim(u,F,Re,B,a)  result(P)
     real, intent(in)        :: u,F,Re,B,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=-a**2/(3.*F)+(a**4*u**2)/(45.*F**3)-(2*a**6*u**4)/(945.*F**5) + &
              (a**8*u**6)/(4725.*F**7) - (2*a**10*u**8)/(93555.*F**9) + &
              (1382*a**12*u**10)/(6.38512875e8*F**11) !! obtained using mathematica
       case (2)
            P=-1./u
       case default
            P=(-1./u + 2.*B)*a
    end select
  end function C_nim

  function E_nim(u,F,Re,B,D,a)  result(P)
     real, intent(in)        :: u,F,Re,B,D,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=D/2. + (a*D*u)/(6.*F) - (a**3*D*u**3)/(90.*F**3) + (a**5*D*u**5)/(945.*F**5) - &
             (a**7*D*u**7)/(9450.*F**7) + (a**9*D*u**9)/(93555.*F**9)
       case (2)
            P=D
       case default
            P=D*(1-B*u)
    end select
  end function E_nim

  function F_nim(u,FF,Re,B,D,a)  result(P)
     real, intent(in)        :: u,FF,Re,B,D,a
     real                    :: P
    select case (condition(u,FF,Re))
       case (1)
            P=D/2. - (a*D*u)/(6.*FF) + (a**3*D*u**3)/(90.*FF**3) - (a**5*D*u**5)/(945.*FF**5) + &
             (a**7*D*u**7)/(9450.*FF**7) - (a**9*D*u**9)/(93555.*FF**9)
       case (2)
            P=0
       case default
            P=D*(B*u)
    end select
  end function F_nim


  function G_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=1/(2.*a) + u/(2.*F) + (a*u**2)/(6.*F**2) - (a**3*u**4)/(90.*F**4) + &
              (a**5*u**6)/(945.*F**6) - (a**7*u**8)/(9450.*F**8) + (a**9*u**10)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=exp(Re)*u/(F*(exp(Re)-1.))
            p=u/F/(1.-1./exp(Re))
    end select
  end function G_nim

  function H_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=a/F+(a**2*u)/(3.*F**2)-(a**4*u**3)/(45.*F**4)+(2*a**6*u**5)/(945.*F**6) - &
             (a**8*u**7)/(4725.*F**8) + (2*a**10*u**9)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=(2*a*exp(Re))/(F*(exp(Re)-1))-1/u
            p=2.*a/F/(1-1/exp(Re))-1./u
    end select
  end function H_nim

  function I_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=1/(2.*a) - u/(2.*F) + (a*u**2)/(6.*F**2) - (a**3*u**4)/(90.*F**4) + &
              (a**5*u**6)/(945.*F**6) - (a**7*u**8)/(9450.*F**8) + (a**9*u**10)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=u/(F*(exp(Re)-1))
    end select
  end function I_nim

  function J_nim(u,F,Re,a)  result(P)
     real, intent(in)        :: u,F,Re,a
     real                    :: P
    select case (condition(u,F,Re))
       case (1)
            P=-(a/F) + (a**2*u)/(3.*F**2) - (a**4*u**3)/(45.*F**4) + &
               (2*a**6*u**5)/(945.*F**6) - (a**8*u**7)/(4725.*F**8) + &
               (2*a**10*u**9)/(93555.*F**10)
       case (2)
            P=0
       case default
            P=2.*a/(F*(exp(Re)-1.))-1./u
    end select
  end function J_nim

end module nim

