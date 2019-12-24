!!****m* TOOLS/m_qtlmap_isymmax2sat
!!  NAME
!!    m_qtlmap_isymmax2sat
!!  DESCRIPTION
!!
!!  NOTES
!!
!!  BUGS
!!
!!  HISTORY
!!
!!  SEE ALSO
!!
!!  COPYRIGHT
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/
! Interface du module d andres
module m_qtlmap_isymmax2sat
  use m_qtlmap_constant
  use m_qtlmap_base
  use m_qtlmap_log
  implicit none

!!****f* m_qtlmap_isymmax2sat/solvesymmax2sat
!! NAME
!!    solvesymmax2sat
!! DESCRIPTION
!!     Interface with C++ development of symmax2sat
!!     Calcul de h'W h
!!     resolution pour 1 famille
!!
!!     interpretation de h
!! GENOTYPE :    [A11 A12  A21 A22  A23 A23  ..... ]
!! ORIENTATION H:[  1         -1       1     ..... ]
!! SOLUTION   => [A11 A12  A22 A21  A23 A23  ..... ]
!!
!! SOURCE
  interface
            function solvesymmax2sat (n, m,posx,posy,cost,sol) result (r)
               integer ,intent(in):: n ! number of variable/markers
               integer ,intent(in):: m ! number of constraints
               integer,dimension(m),intent(in) :: posx
               integer,dimension(m),intent(in) :: posy
               double precision,dimension(m),intent(in) :: cost
               integer,dimension(m),intent(inout) :: sol
               integer :: r
            end function solveSymMax2SAT
  end interface
!!***

  public :: get_h_from_w
  public :: test_module_isymmax2sat1
  public :: test_module_isymmax2sat2

  contains
!!****f* m_qtlmap_isymmax2sat/get_h_from_w
!! NAME
!!    get_h_from_w
!! DESCRIPTION
!!
!! SOURCE
    function get_h_from_w(sW,W,H,m) result(ok)
       integer                     ,intent(in)              :: sW
       real(kind=dp),dimension(sW,sW),intent(in)            :: W
       integer                     ,intent(out)             :: m
       integer     ,dimension(sW),intent(out)        :: H
       integer     ,dimension(:)    ,allocatable    :: posx,posy ! On exit dimension value is m
       real(kind=dp),dimension(:)   ,allocatable    :: cost

       integer :: n,ret,io
       logical :: ok
       
       allocate (posx(sW*sW),posy(sW*sW),cost(sW*sW),stat=io)
       if ( io /= 0 ) call stop_application("QTLMap need more memory to run solveSymMax2SAT : not enough memory...")

       cost=0.d0
       posx=0
       posy=0
       n = size(W,1)
       ret=1
       call create_sparse_W (sW,W ,m,posx,posy,cost)
       ret = solveSymMax2SAT(n,m,posx,posy,cost,H)
       ok = ( ret /= 0 )
       deallocate (posx,posy,cost)

    end function get_h_from_w
!!***


!!****f* m_qtlmap_isymmax2sat/create_sparse_W
!! NAME
!!    create_sparse_W
!! DESCRIPTION
!!
!! SOURCE
    ! this programs serves to take full-stored W and translate it into
    ! row, column, value
    ! lf<F4>stored (i<=j) format for non-null elements
    subroutine create_sparse_W (sW,W ,m,posx,posy,cost)
    integer , intent(in)                       :: sW
    real(kind=dp),dimension(sW,sW),intent(in)  :: W
    integer                     ,intent(out) :: m
    integer     ,dimension(sW*sW)   ,intent(out) :: posx,posy ! On exit dimension value is m
    real(kind=dp),dimension(sW*sW)  ,intent(out)  :: cost ! On exit dimension value is m

    integer:: i,j,n
    real(kind=dp):: val

    if ( size(W,1) /= size(W,2)) then
      call log_mess("W,dim 1:"//str(size(W,1)),ERROR_DEF)
      call log_mess("W,dim 2:"//str(size(W,2)),ERROR_DEF)
      call stop_application("Devel error : create_sparse_W W have to be a squared matrice ")
    end if

    if ( size(W,1) <=0 ) then
      call stop_application("Devel error : create_sparse_W  W have a no dimension (size=0)")
    end if

    n = size(W,1)
    m = 0
    ! loop over values
    do i=1,n
     do j=i,n
      if(W(i,j)/=0) then
   !    if (i<=j) then
!        print *,m,W(i,j),'SIZE,W : ',size(W,1),size(cost),size(posx)
        m = m+1
        posx(m)=i
        posy(m)=j
        cost(m)=W(i,j)
    !   endif
      endif
     enddo
    enddo

   call log_mess("create_sparse_W find m:"//str(m)//' non-null values in W',INFO_DEF)


   !print *,'m= ',m,'non-null values'
   !print *,'in',count((count_row/=0)),'rows from n= ',n

   end subroutine create_sparse_W
!!***

!!****f* m_qtlmap_isymmax2sat/test_module_isymmax2sat1
!! NAME
!!    test_module_isymmax2sat1
!! DESCRIPTION
!!
!! SOURCE
   subroutine test_module_isymmax2sat1

      integer   :: n,m
      integer,dimension(3) :: posx,posy,sol
      double precision,dimension(3) :: cost
      integer :: step
      posx(1)=1
      posx(2)=1
      posx(3)=2
      posy(1)=2
      posy(2)=3
      posy(3)=3
      cost(1)=1.d0
      cost(2)=-3.d0
      cost(3)=-2.d0
      n = 3
      m = 3
      step = solveSymMax2SAT(n,m,posx,posy,cost,sol)
      print *,'res:',step
      print *,sol(1),sol(2),sol(3)
      return

   end subroutine test_module_isymmax2sat1

   subroutine test_module_isymmax2sat2

      integer   :: i,m
      real(kind=dp),dimension(20,20)  :: W
      integer      ,dimension(20)     :: H
      logical                         :: ok

      data (W(1,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(2,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,-0.46352d0,-0.39319d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.11661d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(3,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(4,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(5,i), i = 1, 20) / &
      0.00000d0,-0.46352d0,0.00000d0,0.00000d0,0.00000d0,0.73634d0,2.25498d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.14534d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(6,i), i = 1, 20) / &
      0.00000d0,-0.39319d0,0.00000d0,0.00000d0,0.73634d0,0.00000d0,0.73634d0,0.00000d0,0.00000d0,0.00000d0,0.33943d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(7,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,2.25498d0,0.73634d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,1.96596d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.20458d0,0.18216d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(8,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(9,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(10,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(11,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.33943d0,1.96596d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,1.01829d0,-0.29629d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(12,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(13,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(14,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(15,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(16,i), i = 1, 20) / &
      0.00000d0,0.11661d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.20458d0,0.00000d0,0.00000d0,0.00000d0,1.01829d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.73634d0,0.00000d0,0.00000d0,-0.39319d0/
      data (W(17,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.14534d0,0.00000d0,0.18216d0,0.00000d0,0.00000d0,0.00000d0,-0.29629d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.73634d0,0.00000d0,0.00000d0,0.00000d0,-1.85409d0/
      data (W(18,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(19,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0/
      data (W(20,i), i = 1, 20) / &
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,0.00000d0,&
      0.00000d0,0.00000d0,0.00000d0,0.00000d0,-0.39319d0,-1.85409d0,0.00000d0,0.00000d0,0.00000d0/

      print *,'N=20'
      print *,'M doit etre egale a 16'
      ok = get_h_from_w(20,W,H,m)

!      print *,'res:',step
!      print *,sol(1),sol(2),sol(3)
      return

   end subroutine test_module_isymmax2sat2
!!***

end module m_qtlmap_isymmax2sat

