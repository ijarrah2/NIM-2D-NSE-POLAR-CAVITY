module general
use mesh
use user_parameters
use parameters
implicit none
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
  function analytical_xy (x,y,t,a) result (TT)
    real, intent (in)            :: x,y,t
    integer, intent(in)                      :: a
    real                         :: TT
    if (a==velx1) then
       TT=0.0
    elseif(a==velx2) then
       TT=0.0
    elseif(a==PPE) then
       TT=0.0
    endif
  end function analytical_xy

  function Initial_Condition (x,y,id) result (TT)
    real, intent (in)            :: x,y
    integer, intent(in)                      :: id
    real                         :: TT
    if (id==velx1) then
      TT=0.0
    elseif(id==velx2) then
      TT=0.0
    elseif(id==PPE) then
      TT=0.0
    endif
  end function Initial_Condition

  function Diffusion (x,y,t) result (Diff)
    real, intent (in)            :: x,y,t
    real                         :: Diff
    Diff=1./Re
  end function Diffusion

  function density (x,y,t) result (D)
    real, intent (in)            :: x,y,t
    real                         :: D
    D=1.
  end function density

  function gforce_x (x,y,t) result (S)
    real, intent (in)            :: x,y,t
    real                         :: S
    S=0.
  end function gforce_x

  function gforce_y (x,y,t) result (S)
    real, intent (in)            :: x,y,t
    real                         :: S
    S=0.

  end function gforce_y


  function boundary_condition_1 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT, theta, r
	r=1.
	theta=atan(y/x)
    if (id==velx1) then
       TT=-sin(theta)
    elseif(id==velx2) then
       TT=cos(theta)
    elseif(id==PPE) then
       TT=0.
    endif
  end function boundary_condition_1
  

  function boundary_condition_2 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    endif
  end function boundary_condition_2
  
  function boundary_condition_3 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    endif
  end function boundary_condition_3

  function boundary_condition_4 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==velx1) then
       TT=0.
    elseif(id==velx2) then
       TT=0.
    elseif(id==PPE) then
       TT=0.
    endif
  end function boundary_condition_4

  function functions_id (x,y,t,id) result (fun)
    real, intent (in)            :: x,y
    real, optional               :: t
	integer, intent(in)          :: id
    real                         :: fun

    SELECT CASE (id)
      CASE (analytical_velx1)
           fun=analytical_xy (x,y,t,velx1)
      CASE (analytical_velx2)
           fun=analytical_xy (x,y,t,velx2)
      CASE (analytical_p)
           fun=analytical_xy (x,y,t,PPE)
      CASE (IC_velx1)
           fun=Initial_Condition(x,y,velx1)
      CASE (IC_velx2)
           fun=Initial_Condition(x,y,velx2)
      CASE (IC_PPE)
           fun=Initial_Condition(x,y,PPE)
      CASE (bx_fun)
           fun=gforce_x(x,y,t)
      CASE (by_fun)
           fun=gforce_y(x,y,t)
      CASE DEFAULT
           print*, "ERROR : Unknown function call (functions_id) with id : ", id
    END SELECT

  end function functions_id




  function BCs_derivatives_1 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdt) then
       TT=0.
    elseif(id==bc_dvdtt) then
       TT=0.
    elseif(id==bc_dvdttt) then
       TT=0.
    endif
  end function BCs_derivatives_1

  function BCs_derivatives_2 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdt) then
       TT=0.
    elseif(id==bc_dvdtt) then
       TT=0.
    elseif(id==bc_dvdttt) then
       TT=0.
    endif
  end function BCs_derivatives_2

  function BCs_derivatives_3 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdt) then
       TT=0.
    elseif(id==bc_dvdtt) then
       TT=0.
    elseif(id==bc_dvdttt) then
       TT=0.
    endif
  end function BCs_derivatives_3


  function BCs_derivatives_4 (x,y,id) result (TT)
    real, intent (in)            :: x,y
	integer, intent(in)          :: id
    real                         :: TT
    if (id==bc_dvdt) then
       TT=0.
    elseif(id==bc_dvdtt) then
       TT=0.
    elseif(id==bc_dvdttt) then
       TT=0.
    endif
  end function BCs_derivatives_4

  function BCs_id (x,y,field_id,bc_id) result (fun)
    real, intent (in)            :: x,y
	integer, intent(in)          :: field_id, bc_id
    real                         :: fun

    SELECT CASE (bc_id)
      CASE (1)
           fun=boundary_condition_1(x,y,field_id)
      CASE (2)
           fun=boundary_condition_2(x,y,field_id)
      CASE (3)
           fun=boundary_condition_3(x,y,field_id)
      CASE (4)
           fun=boundary_condition_4(x,y,field_id)
      CASE DEFAULT
           print*, "ERROR : Unknown function call (BCs_id) with id : ", bc_id
    END SELECT

  end function BCs_id


  function BCs_derivatives (x,y,bc_id,dvdt_id) result (fun)
    real, intent (in)            :: x,y
	integer, intent(in)          :: dvdt_id, bc_id
    real                         :: fun

    SELECT CASE (bc_id)
      CASE (1)
           fun=BCs_derivatives_1(x,y,dvdt_id)
      CASE (2)
           fun=BCs_derivatives_2(x,y,dvdt_id)
      CASE (3)
           fun=BCs_derivatives_3(x,y,dvdt_id)
      CASE (4)
           fun=BCs_derivatives_4(x,y,dvdt_id)
      CASE DEFAULT
           print*, "ERROR : Unknown function call (BCs_derivatives) with id : ", bc_id
    END SELECT

  end function BCs_derivatives

! ======================================================================================

  function Gauss_integration2d (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,j
    Int=0.
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
             Int=Int+Gauss_w_i(i)*Gauss_w_i(j)*A(i,j)
       enddo
    enddo
  end function Gauss_integration2d


  function Gauss_integration1d (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i
    Int=0.
    do i=1,Gauss_int_p
           Int=Int+Gauss_w_i(i)*A(i)
    enddo
  end function Gauss_integration1d


  function Gauss_integration2d_t (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,j,k
    Int=0.
    do i=1,Gauss_int_p
       do j=1,Gauss_int_p
          do k=1,Gauss_int_p
             Int=Int+Gauss_w_i(i)*Gauss_w_i(j)*Gauss_w_i(k)*A(i,j,k)
          enddo
       enddo
    enddo
    Int=Int/2.
  end function Gauss_integration2d_t

 function Gauss_integration1d_t (A) result (Int)
    real, intent (in), dimension(1:Gauss_int_p,1:Gauss_int_p) :: A
    real                                      :: Int
    integer :: i,k
    Int=0.
    do i=1,Gauss_int_p
       do k=1,Gauss_int_p
           Int=Int+Gauss_w_i(i)*Gauss_w_i(k)*A(i,k)
       enddo
    enddo
    Int=Int/2.
  end function Gauss_integration1d_t

  subroutine map_ti(time,t_i)
! This function is to map the time variable to the -1 -> 1 domain
    real, intent(inout)              :: t_i
    real, intent(in)                 :: time
    real                             :: a, b
    a = TI-DT
    b = TI+DT
    t_i= (b-a)/2. *time +(b+a)/2.
  end subroutine map_ti

 function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real, intent(in) :: A(2,2)   !! Matrix
    real             :: B(2,2)   !! Inverse matrix
    real             :: detinv
    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

 function matdet2(A) result(B)
    real, intent(in) :: A(2,2)   !! Matrix
    real             :: B

    B = (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function



 function mattrans2(A) result(B)
    real, intent(in) :: A(2,2)   !! Matrix
    real             :: B(2,2)   !! Transpose matrix
    B(1,1) = A(1,1)
    B(2,1) = A(1,2)
    B(1,2) = A(2,1)
    B(2,2) = A(2,2)
  end function


subroutine map_xy_2d(facei,x,y,x2,y2)
! Inputs are face, local x and y
  implicit none
  type(cell)                          :: facei
  real,intent(in)      :: x,y
  real,intent(inout)   :: x2,y2

  x2=facei%C(1)+facei%C(2)*x+facei%C(3)*y+facei%C(4)*x*y
  y2=facei%D(1)+facei%D(2)*x+facei%D(3)*y+facei%D(4)*x*y

end subroutine map_xy_2d


function dot(a,b) result (x)
	implicit none
	real,dimension(1:2), intent (in)  :: a,b
	real                              :: x
	x=a(1)*b(1) + a(2)*b(2)
end function dot


end module general


