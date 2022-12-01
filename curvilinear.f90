module curvilinear
use mesh
use general
implicit none

contains

subroutine evaluate_curvilinear(M,C,x,y)
  implicit none
  real, intent(in)        :: x,y
  type(cell), intent(in)              :: M
  type(curvi2d),intent(inout)         :: C

  ! Find g1, g2 vectors
  C%jacmat= jacobian_matrix(M,x,y)
  C%g_1=(/C%jacmat(1,1),C%jacmat(1,2)/)
  C%g_2=(/C%jacmat(2,1),C%jacmat(2,2)/)
  ! Find g_1, g_2 vectors
  C%jacinv= inverse_jocobian_matrix(M,C%jacmat)
  C%g1=(/C%jacinv(1,1),C%jacinv(2,1)/)
  C%g2=(/C%jacinv(1,2),C%jacinv(2,2)/)
  ! Find metric matrices
  C%g_ij= g_ij_matrix(C%g_1,C%g_2)
  C%gij = gij_matrix (C%g1 ,C%g2 )
  ! Determinants  
  C%gamma_tensor = find_gamma_tensor(M,C%jacinv)
  C%detG = matdet2(C%g_ij)
  C%detJ = matdet2(C%jacmat)

  ! x_xi, xi_x arrays
  C%x_xi(1,1) = C%jacmat(1,1)
  C%x_xi(1,2) = C%jacmat(2,1)
  C%x_xi(2,1) = C%jacmat(1,2)
  C%x_xi(2,2) = C%jacmat(2,2)
  
  C%xi_x(1,1) = C%jacinv(1,1)
  C%xi_x(1,2) = C%jacinv(2,1)
  C%xi_x(2,1) = C%jacinv(1,2)
  C%xi_x(2,2) = C%jacinv(2,2)



end subroutine evaluate_curvilinear


function jacobian_matrix(M,x,y)  result(J)
  implicit none
  real, intent(in)        :: x,y
  type(cell), intent(in)              :: M
  real, dimension(1:2,1:2):: J

  J(1,1)=M%C(2)+M%C(4)*y
  J(1,2)=M%D(2)+M%D(4)*y

  J(2,1)=M%C(3)+M%C(4)*x
  J(2,2)=M%D(3)+M%D(4)*x

end function jacobian_matrix


function g_ij_matrix(g1,g2)  result(A)
  implicit none
  real, dimension(1:2),intent(in):: g1,g2
  real, dimension(1:2,1:2)       :: A

  A(1,1)=DOT_PRODUCT(g1,g1)
  A(1,2)=DOT_PRODUCT(g1,g2)
  A(2,1)=DOT_PRODUCT(g2,g1)
  A(2,2)=DOT_PRODUCT(g2,g2)

end function g_ij_matrix


function inverse_jocobian_matrix(M,J)  result(J2)
  implicit none
  real, dimension(1:2,1:2),intent(in):: J
  type(cell), intent(in)              :: M
  real, dimension(1:2,1:2):: J2

  J2=matinv2(J)
end function inverse_jocobian_matrix

function gij_matrix(g1,g2)  result(A)
  implicit none
  real, dimension(1:2),intent(in):: g1,g2
  real, dimension(1:2,1:2)       :: A

  A(1,1)=DOT_PRODUCT(g1,g1)
  A(1,2)=DOT_PRODUCT(g1,g2)
  A(2,1)=DOT_PRODUCT(g2,g1)
  A(2,2)=DOT_PRODUCT(g2,g2)

end function gij_matrix


function find_gamma_tensor(M,J)  result(A)
  implicit none
  real, dimension(1:2,1:2),intent(in):: J
  type(cell), intent(in)                         :: M
  real, dimension(1:2,1:2,1:2)       :: A

  A(1,1,1) = 0.
  A(1,1,2) = M%C(4)*J(1,1)+M%D(4)*J(2,1)
  A(1,2,1) = A(1,1,2)
  A(1,2,2) = 0.
  A(2,1,2) = M%C(4)*J(1,2)+M%D(4)*J(2,2)
  A(2,2,1) = A(2,1,2)
  A(2,2,2) = 0.
end function find_gamma_tensor



end module curvilinear
