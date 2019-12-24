!!****m* TOOLS/m_qtlmap_solve_linear
!!  NAME
!!    m_qtlmap_solve_linear
!!  DESCRIPTION
!!    set of functions to solve a general system of linear equations.
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!
!!  COPYRIGHT
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/
module m_qtlmap_solve_linear
    use m_qtlmap_constant
    use m_qtlmap_log
    use m_qtlmap_base
    implicit none

!solves the double precision system
!     A * X = B
    public  :: solve_lin
    private :: solve_lin_dgefs


    contains

!!****f* m_qtlmap_solve_linear/solve_lin
!! NAME
!!   solve_lin
!! DESCRIPTION
!!   general interface to solve a system of linear equations
!!   A * X = B
!! INPUTS
!!   A    : matrix
!!  LDA   : First dimension of A
!!  N     : the dimension of B
!!  B     : vector
!! OUTPUTS
!!  X     : the solution
!! ERROR  : error message
!!
!! NOTE
!!  This module is not finalized. only the slatec method is defined....
!! SOURCE
      subroutine solve_lin(A, LDA, N, B,X,ERROR)
            integer                          , intent(in) ::   LDA,N
            real(kind=dp)  , dimension(LDA,N), intent(in) ::   A
            real(kind=dp)  , dimension(N)    , intent(in) ::   B
            real(kind=dp)  ,dimension(N)    , intent(out) ::   X
            integer                         , intent(out) ::   ERROR

            call solve_lin_dgefs(A, LDA, N, B,X,ERROR)

      end subroutine solve_lin
!!***

!!****f* m_qtlmap_solve_linear/solve_lin_dgefs
!! NAME
!!   solve_lin_dgefs
!! DESCRIPTION
!!   use the slatec function DGEFS (LU method) to solve a system of linear equations
!!   A * X = B
!!
!!   *** Error Messages Printed ***
!!
!!    IND=-1  terminal   N is greater than LDA.
!!    IND=-2  terminal   N is less than 1.
!!    IND=-3  terminal   ITASK is less than 1.
!!    IND=-4  terminal   The matrix A is computationally singular.
!!                         A solution has not been computed.
!!    IND=-10 warning    The solution has no apparent significance.
!!                         The solution may be inaccurate or the matrix
!!                         A may be poorly scaled.
!! INPUTS
!!   A    : matrix
!!  LDA   : First dimension of A
!!  N     : the dimension of B
!!  B     : vector
!! OUTPUTS
!!  X     : the solution
!! ERROR  : error message
!!
!! SOURCE
      subroutine solve_lin_dgefs(A, LDA, N, B,X,ERROR)
            integer                          , intent(in) ::   LDA,N
            real(kind=dp)  , dimension(LDA,N), intent(in) ::   A
            real(kind=dp)  , dimension(N)    , intent(in) ::   B
            real(kind=dp)  ,dimension(N)    , intent(out) ::   X
            integer                         , intent(out) ::   ERROR

            real (kind=dp) , dimension(N)           :: WORK
            integer, dimension(N)                   :: IWORK

            external :: DGEFS
!
           ERROR=0
           X=B
           ! Solve a general system of linear equations.
           call DGEFS (A, LDA, N, X, 1 , ERROR, WORK, IWORK)
      end subroutine solve_lin_dgefs
!!***


end module m_qtlmap_solve_linear
