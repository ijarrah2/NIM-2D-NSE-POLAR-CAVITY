!
!   4____4_____3
!   |          |
!   |          | 
! 1 |          |2  
!   |          |  
!   |____3_____|          
!   1          2
module mesh
  implicit none
type curvi2d
  real, dimension(1:2)        :: g_1,g_2,g1,g2
  real                        :: g_11,g_12,g_22,g_21,g11,g12,g22,g21,detG,detJ
  real, dimension(1:2,1:2)    :: jacmat,jacinv,g_ij,gij
  real, dimension(1:2,1:2,1:2):: gamma_tensor
  real, dimension(1:2,1:2)    :: xi_x, x_xi
end type curvi2d
! ==================== NODE ==========================
  type node
    integer                   :: n=0
    real                      :: x,y
    real                      :: u,v
    real                      :: T
  end type node
! ==================== EDGE ==========================
  type edge
    integer                   :: n=0
    integer,dimension(1:2)    :: shared_faces
    type(node),dimension(1:2) :: nodes
    logical,dimension(1:3)    :: BC=(/.false.,.false.,.false./)
    integer,dimension(1:3)    :: BC_type=(/1,1,1/)  ! 1: dirichlet, 2: neumann, 3:robin
    real                      :: length, mid_x,mid_y
  end type edge
! ==================== Cell =========================
  type cell
    integer                   :: n=0   
    type(node),dimension(1:4) :: nodes
    integer,dimension(1:4)    :: edges
    real, dimension(1:4)      :: C,D         ! Mapping coefficients
    real                      :: area 
    integer,dimension(1:4)    :: neighbors
    type(node)                :: center      !  x,y of center point
    real                      :: a,b
    integer                   :: i,j         ! arrangment in grid
  end type cell
! ====================================================
  contains
	
end module mesh
