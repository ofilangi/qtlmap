!!****m* TOOLS/m_qtlmap_optimization
!!  NAME
!!    m_qtlmap_optimization
!!  DESCRIPTION
!!    set of functions to find a minimum of a objective function.
!!    this module is an interface to offers :
!!     * a generic interface for all qtlmap analysis
!!     * a specific implementation for full and half sib family structure and likelihood function based on.
!!
!! - We keep a compatibility with the oldest optimization NAG method e04jyf : http://www.nag.co.uk/numeric/FL/manual/xhtml/E04/e04jyf.xml
!! The NAG optimization is only available with a NAG compiler (xlf,f95)
!!
!! source of implementation :
!!  * L-BFGS : http://users.eecs.northwestern.edu/~nocedal/lbfgs.html
!!  * NLOPT  : http://ab-initio.mit.edu/wiki/index.php/NLopt
!!
!! The module offers the following optimization :
!!
!! 1  : * NAG * E04JYF
!! 2  : * L-BFGS-B * Broyden–Fletcher–Goldfarb–Shanno update to approximate the Hessian matrix (L-BFGS stands for (limited memory BFGS))
!! 12 : * NLopt * DIRECT (global, no-derivative)
!! 13 : * NLopt * DIRECT-L (global, no-derivative)
!! 14 : * NLopt * Randomized DIRECT-L (global, no-derivative)
!! 15 : * NLopt * Unscaled DIRECT (global, no-derivative)
!! 16 : * NLopt * Unscaled DIRECT-L (global, no-derivative)
!! 17 : * NLopt * Unscaled Randomized DIRECT-L (global, no-derivative)
!! 18 : * NLopt * Original DIRECT version (global, no-derivative)
!! 19 : * NLopt * Original DIRECT-L version (global, no-derivative)
!! 23 : * NLopt * Limited-memory BFGS (L-BFGS) (local, derivative-based)
!! 24 : * NLopt * Principal-axis, praxis (local, no-derivative)
!! 25 : * NLopt * Limited-memory variable-metric, rank 1 (local, derivative-based)
!! 26 : * NLopt * Limited-memory variable-metric, rank 2 (local, derivative-based)
!! 27 : * NLopt * Truncated Newton (local, derivative-based)
!! 28 : * NLopt * Truncated Newton with restarting (local, derivative-based)
!! 29 : * NLopt * Preconditioned truncated Newton (local, derivative-based)
!! 30 : * NLopt * Preconditioned truncated Newton with restarting (local, derivative-based)
!! 31 : * NLopt * Controlled random search (CRS2) with local mutation (global, no-derivative)
!! 32 : * NLopt * Multi-level single-linkage (MLSL), random (global, no-derivative)
!! 33 : * NLopt * Multi-level single-linkage (MLSL), random (global, derivative)
!! 34 : * NLopt * Multi-level single-linkage (MLSL), quasi-random (global, no-derivative)
!! 35 : * NLopt * Multi-level single-linkage (MLSL), quasi-random (global, derivative)
!! 36 : * NLopt * Method of Moving Asymptotes (MMA) (local, derivative)
!! 37 : * NLopt * COBYLA (Constrained Optimization BY Linear Approximations) (local, no-derivative)
!! 38 : * NLopt * NEWUOA unconstrained optimization via quadratic models (local, no-derivative)
!! 39 : * NLopt * Bound-constrained optimization via NEWUOA-based quadratic models (local, no-derivative)
!! 40 : * NLopt * Nelder-Mead simplex algorithm (local, no-derivative)
!! 41 : * NLopt * Sbplx variant of Nelder-Mead (re-implementation of Rowan's Subplex) (local, no-derivative)
!! 42 : * NLopt * Augmented Lagrangian method (local, no-derivative)
!! 43 : * NLopt * Augmented Lagrangian method (local, derivative)
!! 44 : * NLopt * Augmented Lagrangian method for equality constraints (local, no-derivative)
!! 45 : * NLopt * Augmented Lagrangian method for equality constraints (local, derivative)
!! 46 : * NLopt * BOBYQA bound-constrained optimization via quadratic models (local, no-derivative)
!! 47 : * NLopt * ISRES evolutionary constrained optimization (global, no-derivative)
!!
!!
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
module m_qtlmap_optimization
      use m_qtlmap_constant
      use m_qtlmap_math
      use m_qtlmap_log
      use m_qtlmap_base
      use m_qtlmap_types

      implicit none
      private
      save

!!****d* m_qtlmap_optimization/OPTI_NAG
!! NAME
!!   OPTI_NAG
!! DESCRIPTION
!!   constant id e04jyf nag optimization
!! NOTES
!!   http://www.nag.co.uk/numeric/FL/manual/xhtml/E04/e04jyf.xml
!! SOURCE
      integer , parameter , private            :: OPTI_NAG                  = 1
!!***

!!****v* m_qtlmap_optimization/NB_OPTI_NLOPT
!! NAME
!!   NB_OPTI_NLOPT
!! DESCRIPTION
!!   Number of nplot optimization available
!! SOURCE
      integer             , private           :: NB_OPTI_NLOPT             = 0
!!***

!!****v* m_qtlmap_optimization/opti_user_fam
!! NAME
!!   opti_user_fam
!! DESCRIPTION
!!   default optimization for a specific objective function (likelihood with a family organization)
!! SOURCE
#ifdef HAVE_NLOPT
      integer ,  private                       :: opti_user_fam             = OPTI_LAST+16 ! EQUI PNET in NLOPT
#else
      integer ,  private                       :: opti_user_fam             = OPTI_LBFGS
#endif
!!***

!!****d* m_qtlmap_optimization/LIBNAG_AVAILABLE
!! NAME
!!   LIBNAG_AVAILABLE
!! DESCRIPTION
!!   flag to know if the nag e04jyf is available
!! SOURCE
#ifdef HAVE_LIBNAG
      logical , parameter,private             :: LIBNAG_AVAILABLE = .true.
#else
      logical , parameter,private             :: LIBNAG_AVAILABLE = .false.
#endif
!!***

!!****v* m_qtlmap_optimization/opti_user
!! NAME
!!   opti_user
!! DESCRIPTION
!!   default optimization for a classical optimization of a objective function
!! SOURCE
#ifdef HAVE_LIBNAG
      integer ,   public                      :: opti_user = OPTI_NAG
#else
      integer , public                        :: opti_user = OPTI_LAST+16! OPTI_LBFGS
#endif
!!***

!!****d* m_qtlmap_optimization/OPT_CLASSIC
!! NAME
!!   OPT_CLASSIC
!! DESCRIPTION
!!   Internal ID to execute a classical optimization
!! SOURCE
      integer , parameter , private            :: OPT_CLASSIC          = 1
!!***

!!****d* m_qtlmap_optimization/OPT_FAM
!! NAME
!!   OPT_FAM
!! DESCRIPTION
!!   Internal ID to execute a likelihood optimization based on full sib family
!! SOURCE
      integer , parameter , private            :: OPT_FAM              = 2
!!***

!!****d* m_qtlmap_optimization/OPT_FAM_SIRE
!! NAME
!!   OPT_FAM_SIRE
!! DESCRIPTION
!!   Internal ID to execute a likelihood optimization based on half sib family
!! SOURCE
      integer , parameter , private            :: OPT_FAM_SIRE         = 3
!!***

!!****d* m_qtlmap_optimization/OPT_FAM_MULTI
!! NAME
!!   OPT_FAM_MULTI
!! DESCRIPTION
!!   Internal ID to execute a likelihood optimization based on full sib family and multi trait
!! SOURCE
      integer , parameter , private            :: OPT_FAM_MULTI        = 4
!!***

!!****v* m_qtlmap_optimization/determ
!! NAME
!!   determ
!! DIMENSIONS
!!   np
!! DESCRIPTION
!!   pointer : a determinant vector
!! SOURCE
      real(kind=dp),dimension(:), private, pointer      ::  determ
!!***

!!****v* m_qtlmap_optimization/inv_vci
!! NAME
!!   inv_vci
!! DIMENSIONS
!!   np,ncar,ncar
!! DESCRIPTION
!!   pointer : the inverse of the covariance phenotypics trait matrix
      real(kind=dp),dimension(:,:,:), private,pointer   ::  inv_vci
!!***

!!****v* m_qtlmap_optimization/filter_vci
!! NAME
!!   filter_vci
!! DIMENSIONS
!!   np,nbnivest
!! DESCRIPTION
!!   pointer : a boolean vector to active or not the computation of inv_vci
      logical,dimension(:,:), private, pointer          :: filter_vci
!!***

      type(QTLMAP_DATASET)           ,pointer, private  :: P_DATA


      !$omp threadprivate (P_DATA,inv_vci,filter_vci,determ)
#ifdef BENCHMARK_VIEW
!!****iv* m_qtlmap_optimization/count_call_funct
!! NAME
!!   count_call_funct
!! DESCRIPTION
!!   count the number of function call
!! NOTES
!!   used in the benchmark mode
!!***
      integer                   ,private   :: count_call_funct          = 0
!!****iv* m_qtlmap_optimization/count_call_gradient
!! NAME
!!   count_call_gradient
!! DESCRIPTION
!!   count the number of gradient call
!! NOTES
!!   used in the benchmark mode
!!***
      integer                   ,private   :: count_call_gradient       = 0
!!****iv* m_qtlmap_optimization/count_call_vci
!! NAME
!!   count_call_vci
!! DESCRIPTION
!!   count the number of function (inverse of the covariance matrix) call
!! NOTES
!!   used in the benchmark mode
!!***
      integer                   ,private   :: count_call_vci            = 0
!!****iv* m_qtlmap_optimization/count_avoid_vci
!! NAME
!!   count_avoid_vci
!! DESCRIPTION
!!   count the number of function (inverse of a covariance trait matrix) call avoid
!! NOTES
!!   used in the benchmark mode
!!***
      integer                   ,private   :: count_avoid_vci           = 0
!!****iv* m_qtlmap_optimization/count_avoid_funct
!! NAME
!!   count_avoid_funct
!! DESCRIPTION
!!   count the number of function call avoid
!! NOTES
!!   used in the benchmark mode
!!***
      integer                   ,private   :: count_avoid_funct         = 0
!!****iv* m_qtlmap_optimization/total_time
!! NAME
!!   total_time
!! DESCRIPTION
!!   additional time of computation
!! NOTES
!!   used in the benchmark mode
!!***
      real(kind=dp)             ,private   :: total_time                = 0
#endif


      public :: set_optimization
      public :: minimizing_funct
      public :: minimizing_funct_family
      public :: minimizing_funct_family_sire
      public :: minimizing_funct_family_multi

      private :: name_optim

      private :: minimizing_e04jyf
      private :: minimizing_lbfgs

      contains

!!****f* m_qtlmap_optimization/set_optimization
!!  NAME
!!    set_optimization
!!  DESCRIPTION
!!    set the default optimization
!!  INPUTS
!!    value : a specific optimization id
!!  NOTES
!!   optimization available : see m_qtlmap_optimization
!!  SOURCE
           subroutine set_optimization(value)
               integer, intent(in)   :: value
               integer :: i
               external :: get_nb_opti_nlopt_qtlmap

#ifdef HAVE_NLOPT
               call get_nb_opti_nlopt_qtlmap(NB_OPTI_NLOPT)
#else
               NB_OPTI_NLOPT = 0
#endif
               if ( value == OPTI_NAG .and. .not. LIBNAG_AVAILABLE ) then
                  call log_mess('**** Optimization NAG (e04jyf) is not available, qtlmap use :'//str(opti_user))
                  return
               end if


               if ( (value >= OPTI_NAG)  .and. (value <= OPTI_LAST+NB_OPTI_NLOPT) )  then
                 opti_user = value
               end if

               if ( value /= OPTI_NAG  .and. (value <= OPTI_LAST+NB_OPTI_NLOPT)) then
                  opti_user_fam = value
               end if


           end subroutine
!!***

!!****f* m_qtlmap_optimization/name_optim
!!    name_optim
!!  DESCRIPTION
!!    get the optimization name
!!  INPUTS
!!    optim : the optimization id
!!  OUTPUTS
!!    name  : the optimization name
!!  SOURCE
           function name_optim(optim) result(name)
                integer , intent(in)                           :: optim
                character(len=300)                             :: name
              SELECT CASE (optim)
                CASE (OPTI_NAG)
                   name ="E04JYF"
                CASE (OPTI_LBFGS)
                   name = "LBFGS"
                CASE (OPTI_LAST+1)
                   name = "NLOPT_GN_DIRECT"
                CASE (OPTI_LAST+2)
                   name = "NLOPT_GN_DIRECT_L"
                CASE (OPTI_LAST+3)
                   name = "NLOPT_GN_DIRECT_L_RAND"
                CASE (OPTI_LAST+4)
                   name = "NLOPT_GN_DIRECT_NOSCAL"
                CASE (OPTI_LAST+5)
                   name = "NLOPT_GN_DIRECT_L_NOSCAL"
                CASE (OPTI_LAST+6)
                   name = "NLOPT_GN_DIRECT_L_RAND_NOSCAL"
                CASE (OPTI_LAST+7)
                   name = "NLOPT_GN_ORIG_DIRECT"
                CASE (OPTI_LAST+8)
                   name = "NLOPT_GN_ORIG_DIRECT_L"
                CASE (OPTI_LAST+9)
                   name = "NLOPT_GD_STOGO"
                CASE (OPTI_LAST+10)
                   name = "NLOPT_GD_STOGO_RAND"
                CASE (OPTI_LAST+11)
                   name = "NLOPT_LD_LBFGS_NOCEDAL"
                CASE (OPTI_LAST+12)
                   name = "NLOPT_LD_LBFGS"
                CASE (OPTI_LAST+13)
                   name = "NLOPT_LN_PRAXIS"
                CASE (OPTI_LAST+14)
                   name = "NLOPT_LD_VAR1"
                CASE (OPTI_LAST+15)
                   name = "NLOPT_LD_VAR2"
                CASE (OPTI_LAST+16)
                   name = "NLOPT_LD_TNEWTON"
                CASE (OPTI_LAST+17)
                   name = "NLOPT_LD_TNEWTON_RESTART"
                CASE (OPTI_LAST+18)
                   name = "NLOPT_LD_TNEWTON_PRECOND"
                CASE (OPTI_LAST+19)
                   name = "NLOPT_LD_TNEWTON_PRECOND_RESTART"
                CASE (OPTI_LAST+20)
                   name = "NLOPT_GN_CRS2_LM"
                CASE (OPTI_LAST+21)
                   name = "NLOPT_GN_MLSL"
                CASE (OPTI_LAST+22)
                   name = "NLOPT_GD_MLSL"
                CASE (OPTI_LAST+23)
                   name = "NLOPT_GN_MLSL_LDS"
                CASE (OPTI_LAST+24)
                   name = "NLOPT_GD_MLSL_LDS"
                CASE (OPTI_LAST+25)
                   name = "NLOPT_LD_MMA"
                CASE (OPTI_LAST+26)
                   name = "NLOPT_LN_COBYLA"
                CASE (OPTI_LAST+27)
                   name = "NLOPT_LN_NEWUOA"
                CASE (OPTI_LAST+28)
                   name = "NLOPT_LN_NEWUOA_BOUND"
                CASE (OPTI_LAST+29)
                   name = "NLOPT_LN_NELDERMEAD"
                CASE (OPTI_LAST+30)
                   name = "NLOPT_LN_SBPLX"
                CASE (OPTI_LAST+31)
                   name = "NLOPT_LN_AUGLAG"
                CASE (OPTI_LAST+32)
                   name = "NLOPT_LD_AUGLAG"
                CASE (OPTI_LAST+33)
                   name = "NLOPT_LN_AUGLAG_EQ"
                CASE (OPTI_LAST+34)
                   name = "NLOPT_LD_AUGLAG_EQ"
                CASE (OPTI_LAST+35)
                   name = "NLOPT_LN_BOBYQA"
                CASE (OPTI_LAST+36)
                   name = "NLOPT_GN_ISRES"
                CASE DEFAULT
                   call stop_application("Unknown name optimize method ["//trim(str(optim))//"]")
               END SELECT
           end function
!!***

!!****f* m_qtlmap_optimization/minimizing_funct
!!  NAME
!!    minimizing_funct
!!  DESCRIPTION
!!    apply a classical optimization. The variable opti_user determines which method is used.
!!  INPUTS
!!     N          : the number n of independent variables. N >=1
!!     IBOUND     : indicates whether the facility for dealing with bounds of special forms is to be used.
!!     FUNCT1     : SUBROUTINE, supplied by the user. External Procedure
!!                  interface : FUNCT_PART(ip,jm,n,X,FM,IUSER,USER)
!!     BL         : lower bounds
!!     BU         : upper bound
!!     IUSER      : user vector of integer
!!     USER       : user vector of real
!!   INPUTS/OUTPUTS
!!     X          : start point and the point where rthe minimum was founded
!!   OUTPUTS
!!     F          : final point stored in X
!!    IFAIL       : an error CODE
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
           subroutine minimizing_funct(DATA,& ! The dataset Genealogy
                                  N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT1, & ! SUBROUTINE, supplied by the user. External Procedure
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  IFAIL)
            type(QTLMAP_DATASET),target :: DATA
            integer,intent(in)          ::  N, IBOUND
            integer,dimension(:)        ::  IUSER
            real(kind=dp),dimension(N)  ::  BL, BU, X
            real(kind=dp)               ::  F
            real(kind=dp),dimension(:)  ::  USER
            integer,intent(out)         ::  IFAIL
            external                    :: FUNCT1

            real(kind=dp) ,dimension(DATA%GENEA%NP)   :: FP
            real(kind=dp) ,dimension(DATA%GENEA%NM)   :: FM
            logical       ,dimension(DATA%GENEA%NP,N) :: FILTER_FP
            logical       ,dimension(DATA%GENEA%NP,DATA%GENEA%NM,N) :: FILTER_FM



#ifdef BENCHMARK_VIEW
            real(kind=dp)               :: t1,t2,temp
            temp=0
            t1=second();
            count_call_funct          = 0
            count_call_gradient       = 0
            count_avoid_funct         = 0
#endif

           P_DATA => DATA
           ! correction BUG 02/11/2010 : si n=0 aucun parametre a estimer on sort
           if (N == 0 ) then
             call log_mess("minimizing_funct is call with a objective function and none parameter",WARNING_DEF)
             return
           end if

           select case (opti_user)
              case (OPTI_NAG)
                call minimizing_e04jyf(DATA,N, IBOUND, FUNCT1, BL, BU, X, F,IUSER, USER,FP,FILTER_FP,FM,FILTER_FM,IFAIL)
              case (OPTI_LBFGS)
                call minimizing_lbfgs(DATA,N, IBOUND, FUNCT1, BL, BU, X, F,IUSER, USER,FP,FILTER_FP,FM,FILTER_FM,&
                IFAIL,OPT_CLASSIC,likelihood_empty)
              case default
                if ( opti_user > OPTI_LAST .and. opti_user <= (OPTI_LAST+NB_OPTI_NLOPT) ) then
                    call minimizing_nlopt(DATA,opti_user-OPTI_LAST,N, IBOUND, FUNCT1, BL, BU, X, F,IUSER, &
                      USER,FP,FILTER_FP,FM,FILTER_FM,IFAIL,OPT_CLASSIC,likelihood_empty)
                else
                  call stop_application('bad value of opti_user['//str(opti_user)//'].')
                end if
           end select

          IF (F .NE. F) THEN
              F=INIFINY_REAL_VALUE
              IFAIL=11
          END IF

          call log_mess(trim(name_optim(opti_user))//' MAX VALUE:'//str(f),VERBOSE_DEF)

#ifdef BENCHMARK_VIEW
         t2=second()-t1
         total_time=total_time+t2
         if ((count_call_funct+count_avoid_funct)>0) temp=(count_avoid_funct/(count_call_funct+count_avoid_funct))*100
         print '(a," ",f10.8,1x,f7.3,1x,a,i7,1x,a,i7,1x,a,i7,f4.1,"%",f7.3)',trim(name_optim(opti_user)),t2,f,'call funct:',&
         count_call_funct/nm,'call gradient:',count_call_gradient,'avoid call:',count_avoid_funct/nm,temp,total_time
#endif

      end subroutine minimizing_funct
!!***

!!****f* m_qtlmap_optimization/minimizing_funct_family
!!  NAME
!!    minimizing_funct_family
!!  DESCRIPTION
!!    apply a likelihood objective function optimization using the full sib family structuration and a flag by family
!!    that determines if a parameter influence the likelihood of the family.
!!    The variable opti_user_fam determines which method is used.
!!
!!
!!
!!  INPUTS
!!     N            : the number n of independent variables. N >=1
!!     IBOUND       : indicates whether the facility for dealing with bounds of special forms is to be used.
!!     FUNCT_PART   : SUBROUTINE of the sub objective function for the family ip,jm
!!                    interface : FUNCT_PART(ip,jm,n,X,FM,IUSER,USER)
!!     filter_incid : boolean array that indicates if a parameter influence the likelihood on a specific family.
!!                    dimension : np,nm,N
!!     BL           : lower bounds
!!     BU           : upper bound
!!     IUSER        : user vector of integer
!!     USER         : user vector of real
!!   INPUTS/OUTPUTS
!!     X            : start point and the point where rthe minimum was founded
!!   OUTPUTS
!!     F            : final point stored in X
!!    fm_user       : value of the sub objective function full sib family . F = SUM FM(JM), 1<=JM<=NM
!!    fp_user       : value of the sub objective function half sib family . F = SUM FP(IP), 1<=IP<=NP
!!    IFAIL         : an error CODE
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
      subroutine minimizing_funct_family(DATA,&
                                  N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT_PART, & ! SUBROUTINE, supplied by the user Likelihood for the famili IP,JM. External Procedure
                                  filter_incid,&
                                  fm_user,     &
                                  fp_user,     &
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  IFAIL)
            type(QTLMAP_DATASET),target :: DATA
            integer,intent(in)           ::  N, IBOUND
            integer,dimension(:)         ::  IUSER
            real(kind=dp),dimension(N)   ::  BL, BU, X
            real(kind=dp) ,dimension(DATA%GENEA%NM) ::  fm_user
            real(kind=dp) ,dimension(DATA%GENEA%NP) ::  fp_user
            real(kind=dp)               ::  F
            real(kind=dp),dimension(:)  ::  USER
            integer,intent(out)         ::  IFAIL
            logical, dimension(DATA%GENEA%NP,DATA%GENEA%NM,N) , intent(in) :: filter_incid

            external                    :: FUNCT_PART

            logical, dimension(DATA%GENEA%NP,N)    :: filter_fp
            integer                     :: ip

#ifdef BENCHMARK_VIEW
            real(kind=dp)               :: t1,t2,temp
            temp=0
            t1=second();
            count_call_funct          = 0
            count_call_gradient       = 0
            count_avoid_funct         = 0
#endif
            P_DATA => DATA

            ! correction BUG 02/11/2010 : si n=0 aucun parametre a estimer on sort
            if (N == 0 ) then
             call log_mess("minimizing_funct_family is call with a objective function and none parameter",WARNING_DEF)
             return
            end if

            filter_fp = .false.

            select case (opti_user_fam)
              case (OPTI_NAG)
                 call stop_application('can not use this optimization ['//str(opti_user_fam)//'].')
              case (OPTI_LBFGS)
                call minimizing_lbfgs(DATA,N, IBOUND, get_value, BL, BU, X,&
                 F,IUSER, USER,FP_USER,FILTER_FP,FM_USER,filter_incid,&
                IFAIL,OPT_FAM,FUNCT_PART)
              case default
                if ( opti_user_fam > OPTI_LAST .and. opti_user_fam <= (OPTI_LAST+NB_OPTI_NLOPT) ) then
                    call minimizing_nlopt(DATA,opti_user_fam-OPTI_LAST,N, IBOUND, get_value, BL, BU, X, F,IUSER, &
                    USER,FP_USER,FILTER_FP,FM_USER,filter_incid,IFAIL,OPT_FAM,FUNCT_PART)
                else
                  call stop_application('bad value of opti_user['//str(opti_user_fam)//'].')
                end if
           end select

           IF (F .NE. F) THEN
              F=INIFINY_REAL_VALUE
              IFAIL=11
           END IF

           do ip=1,DATA%GENEA%np
             fp_user(ip) = sum(fm_user(DATA%GENEA%nmp(ip)+1:DATA%GENEA%nmp(ip+1)))
           end do

           call log_mess(trim(name_optim(opti_user_fam))//' MAX VALUE:'//str(f),VERBOSE_DEF)

#ifdef BENCHMARK_VIEW
         t2=second()-t1
         total_time=total_time+t2
         !print *,(count_call_funct+count_avoid_funct),count_call_funct,count_avoid_funct
         if ((count_call_funct+count_avoid_funct)>0) then
           temp=(real(count_avoid_funct)/real(count_call_funct+count_avoid_funct))*100
         end if
         !print *,temp
       print '(a," ",f10.8,1x,f7.3,1x,a,i7,1x,a,i7,1x,a,i7,1x,f4.1,"%",f7.3)',trim(name_optim(opti_user_fam)),t2,f,'call funct:',&
         count_call_funct/nm,'call gradient:',count_call_gradient,'avoid call:',count_avoid_funct/nm,temp,total_time
#endif

      end subroutine minimizing_funct_family
!!***

!!****f* m_qtlmap_optimization/minimizing_funct_family_sire
!!  NAME
!!    minimizing_funct_family_sire
!!  DESCRIPTION
!!    apply a likelihood objective function optimization using the half sib family structuration and a flag by sire family
!!    that determines if a parameter influence the likelihood of the sire family.
!!    The variable opti_user_fam determines which method is used.
!!
!!
!!
!!  INPUTS
!!     N            : the number n of independent variables. N >=1
!!     IBOUND       : indicates whether the facility for dealing with bounds of special forms is to be used.
!!     FUNCT_PART   : SUBROUTINE of the sub objective function for the family ip,jm
!!                    interface : FUNCT_PART(ip,n,X,FM,IUSER,USER)
!!     filter_incid : boolean array that indicates if a parameter influence the likelihood on a specific sire family.
!!                    dimension : np,N
!!     BL           : lower bounds
!!     BU           : upper bound
!!     IUSER        : user vector of integer
!!     USER         : user vector of real
!!   INPUTS/OUTPUTS
!!     X            : start point and the point where rthe minimum was founded
!!   OUTPUTS
!!     F            : final point stored in X
!!    fp_user       : value of the sub objective function half sib family . F = SUM FP(IP), 1<=IP<=NP
!!    IFAIL         : an error CODE
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
      subroutine minimizing_funct_family_sire(DATA,&
                                  N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT_PART, & ! SUBROUTINE, supplied by the user Likelihood for the famili IP,JM. External Procedure
                                  filter_incid,&
                                  fp_user,     &
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  IFAIL)
            type(QTLMAP_DATASET),target :: DATA
            integer,intent(in)          ::  N, IBOUND
            integer,dimension(:)        ::  IUSER
            real(kind=dp),dimension(N)  ::  BL, BU, X
            real(kind=dp) ,dimension(DATA%GENEA%np)::  fp_user
            real(kind=dp)               ::  F
            real(kind=dp),dimension(:)  ::  USER
            integer,intent(out)         ::  IFAIL
            logical, dimension(DATA%GENEA%NP,N) , intent(in) :: filter_incid

            external                    :: FUNCT_PART

            integer                     :: ip
            real(kind=dp) ,dimension(DATA%GENEA%NM) ::  FM
            logical, dimension(DATA%GENEA%NP,DATA%GENEA%NM,N)  :: filter_fm

#ifdef BENCHMARK_VIEW
            real(kind=dp)               :: t1,t2,temp
            temp=0
            t1=second();
            count_call_funct          = 0
            count_call_gradient       = 0
            count_avoid_funct         = 0
#endif
            P_DATA => DATA
            ! correction BUG 02/11/2010 : si n=0 aucun parametre a estimer on sort
            if (N == 0 ) then
             call log_mess("minimizing_funct_family_sire is call with a objective function and none parameter",WARNING_DEF)
             return
            end if

            select case (opti_user_fam)
              case (OPTI_NAG)
                 call stop_application('can not use this optimization ['//str(opti_user_fam)//'].')
              case (OPTI_LBFGS)
                call minimizing_lbfgs(DATA,N, IBOUND, get_value, BL, BU, X, F,IUSER, &
                 USER, fp_user, filter_incid, FM, filter_fm, &
                IFAIL,OPT_FAM_SIRE,FUNCT_PART)
              case default
                if ( opti_user_fam > OPTI_LAST .and. opti_user_fam <= (OPTI_LAST+NB_OPTI_NLOPT) ) then
                    call minimizing_nlopt(DATA,opti_user_fam-OPTI_LAST,N, IBOUND, get_value, BL, BU, X, F,IUSER, USER,&
                    fp_user, filter_incid, FM,filter_fm,IFAIL,OPT_FAM_SIRE,FUNCT_PART)
                else
                  call stop_application('bad value of opti_user['//str(opti_user_fam)//'].')
                end if
           end select


          IF (F .NE. F) THEN
              F=INIFINY_REAL_VALUE
              IFAIL=11
          END IF

          call log_mess(trim(name_optim(opti_user_fam))//' MAX VALUE:'//str(f),VERBOSE_DEF)

#ifdef BENCHMARK_VIEW
         t2=second()-t1
         total_time=total_time+t2
         !print *,(count_call_funct+count_avoid_funct),count_call_funct,count_avoid_funct
         if ((count_call_funct+count_avoid_funct)>0) then
           temp=(real(count_avoid_funct)/real(count_call_funct+count_avoid_funct))*100
         end if
         !print *,temp
      print '(a," ",f10.8,1x,f7.3,1x,a,i10,1x,a,i10,1x,a,i9,1x,f4.1,"%",f7.3)',trim(name_optim(opti_user_fam)),t2,f,'call funct:',&
         count_call_funct/nm,'call gradient:',count_call_gradient,'avoid call:',count_avoid_funct/nm,temp,total_time
#endif

      end subroutine minimizing_funct_family_sire
!!***

!!****f* m_qtlmap_optimization/minimizing_funct_family_multi
!!  NAME
!!    minimizing_funct_family_multi
!!  DESCRIPTION
!!    apply a likelihood objective function optimization using the half sib family structuration and a flag by sire family
!!    that determines if a parameter influence the likelihood of the sire family.
!!    The variable opti_user_fam determines which method is used.
!!
!!  INPUTS
!!     N            : the number n of independent variables. N >=1
!!     IBOUND       : indicates whether the facility for dealing with bounds of special forms is to be used.
!!     FUNCT_PART   : SUBROUTINE of the sub objective function for the family ip,jm
!!                    interface : FUNCT_PART(ip,jm,n,X,vci,determ,FM,IUSER,USER)
!!     filter_incid : boolean array that indicates if a parameter influence the likelihood on a specific sire family.
!!                    dimension : np,nm,N
!!     filter_vcinv :
!!                    dimension : ncar,ncar,np
!!     BL           : lower bounds
!!     BU           : upper bound
!!     IUSER        : user vector of integer
!!     USER         : user vector of real
!!   INPUTS/OUTPUTS
!!     X            : start point and the point where rthe minimum was founded
!!    vci_temp      : buffer array.dim : ncar,ncar,np
!!    determ_temp   : buffer array.dim : np
!!   OUTPUTS
!!     F            : final point stored in X
!!    fm_user       : value of the sub objective function full sib family . F = SUM FM(JM), 1<=JM<=NM
!!    fp_user       : value of the sub objective function half sib family . F = SUM FP(IP), 1<=IP<=NP
!!    IFAIL         : an error CODE
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
  subroutine minimizing_funct_family_multi(DATA, &
                                  DPM        , &
                                  N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT_PART, & ! SUBROUTINE, supplied by the user Likelihood for the famili IP,JM. External Procedure
                                  filter_incid,&
                                  filter_vcinv,&
                                  vci_temp,&
                                  determ_temp,&
                                  fm_user,     &
                                  fp_user,     &
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  IFAIL)
            type(QTLMAP_DATASET),target :: DATA
            type(DATAMODEL_BASE),target :: DPM
            integer,intent(in)          ::  N, IBOUND
            integer,dimension(:)        ::  IUSER
            real(kind=dp),dimension(N)  ::  BL, BU, X
            real(kind=dp) ,dimension(DATA%GENEA%NM) ::  fm_user
            real(kind=dp) ,dimension(DATA%GENEA%NP) ::  fp_user
            real(kind=dp) ,dimension(DPM%ncar,DPM%ncar,DATA%GENEA%np),target ::  vci_temp
            real(kind=dp) ,dimension(DATA%GENEA%np),target ::  determ_temp
            real(kind=dp)               ::  F
            real(kind=dp),dimension(:)  ::  USER
            integer,intent(out)         ::  IFAIL
            logical, dimension(DATA%GENEA%NP,DATA%GENEA%NM,N)      , intent(in) :: filter_incid
            logical, dimension(:,:) ,target , intent(in) ::filter_vcinv

            external                    :: FUNCT_PART
            logical, dimension(DATA%GENEA%NP,N)   :: filter_fp

            integer                     :: ip

#ifdef BENCHMARK_VIEW
            real(kind=dp)               :: t1,t2,temp
            temp=0
            t1=second();
            count_call_funct          = 0
            count_call_gradient       = 0
            count_avoid_funct         = 0
            count_call_vci            = 0
            count_avoid_vci           = 0
#endif
            P_DATA => DATA
             ! correction BUG 02/11/2010 : si n=0 aucun parametre a estimer on sort
            if (N == 0 ) then
             call log_mess("minimizing_funct_family_multi is call with a objective function and none parameter",WARNING_DEF)
             return
            end if

            filter_vci=>filter_vcinv
            determ=>determ_temp
            determ=0.d0
            inv_vci=>vci_temp
            inv_vci=0.d0

            select case (opti_user_fam)
              case (OPTI_NAG)
                 call stop_application('can not use this optimization ['//str(opti_user_fam)//'].')
              case (OPTI_LBFGS)
                call minimizing_lbfgs(DATA,N, IBOUND, get_value, BL, BU, X, F,&
                 IUSER, USER, fp_user, FILTER_FP,fm_user,filter_incid,&
                 IFAIL,OPT_FAM_MULTI,FUNCT_PART)
              case default
                 if ( opti_user_fam > OPTI_LAST .and. opti_user_fam <= (OPTI_LAST+NB_OPTI_NLOPT) ) then
                    call minimizing_nlopt(DATA,opti_user_fam-OPTI_LAST,N, IBOUND, get_value, BL, BU, X, F,IUSER, &
                    USER,fp_user, FILTER_FP,fm_user,filter_incid,IFAIL,OPT_FAM_MULTI,FUNCT_PART)
                else
                 call stop_application('bad value of opti_user['//str(opti_user_fam)//'].')
                end if
           end select


          IF (F .NE. F) THEN
              F=INIFINY_REAL_VALUE
              IFAIL=11
          END IF

           do ip=1,DATA%GENEA%np
             fp_user(ip) = sum(fm_user(DATA%GENEA%nmp(ip)+1:DATA%GENEA%nmp(ip+1)))
           end do

          call log_mess(trim(name_optim(opti_user_fam))//' MAX VALUE:'//str(f),VERBOSE_DEF)

#ifdef BENCHMARK_VIEW
         t2=second()-t1
         total_time=total_time+t2
         !print *,(count_call_funct+count_avoid_funct),count_call_funct,count_avoid_funct
         if ((count_call_funct+count_avoid_funct)>0) then
           temp=(real(count_avoid_funct)/real(count_call_funct+count_avoid_funct))*100
         end if
         !print *,temp
      print '(a," ",f10.8,1x,f7.3,1x,a,i10,1x,a,i10,1x,a,i9,1x,f4.1,"%",f7.3,i10,i10)',trim(name_optim(opti_user_fam)),t2,f,&
      'call funct:',count_call_funct/nm,'call gradient:',count_call_gradient,'avoid call:',count_avoid_funct/nm,temp,total_time,&
         count_call_vci,count_avoid_vci
#endif

      end subroutine minimizing_funct_family_multi
!!***

!!****f* m_qtlmap_optimization/minimizing_e04jyf
!!  NAME
!!    minimizing_e04jyf
!!  DESCRIPTION
!!
!!
!!  INPUTS
!!     N            : the number n of independent variables. N >=1
!!     IBOUND       : indicates whether the facility for dealing with bounds of special forms is to be used.
!!     FUNCT1       : SUBROUTINE objective function
!!                    interface : FUNCT_PART(n,X,F,IUSER,USER)
!!     BL           : lower bounds
!!     BU           : upper bound
!!     IUSER        : user vector of integer
!!     USER         : user vector of real
!!   INPUTS/OUTPUTS
!!     X            : start point and the point where rthe minimum was founded
!!   OUTPUTS
!!     F            : final point stored in X
!!    fm_user       : value of the sub objective function full sib family . F = SUM FM(JM), 1<=JM<=NM
!!    fp_user       : value of the sub objective function half sib family . F = SUM FP(IP), 1<=IP<=NP
!!    IFAIL         : an error CODE
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
        subroutine minimizing_e04jyf(DATA,N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT1, & ! SUBROUTINE, supplied by the user. External Procedure
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  FP,         &  ! Likelihood by Family of Sire
                                  FILTER_FP,  &  ! Filter on The Family Sire
                                  FM,         &  ! Likelihood by Family of SIre-Dam
                                  FILTER_FM,  &  ! Filter on family Sire-Dam
                                  IFAIL   )

            type(QTLMAP_DATASET),target :: DATA
            integer,intent(in)          ::  N, IBOUND
            integer,dimension(:)        ::  IUSER
            real(kind=dp),dimension(N)  ::  BL, BU, X
            real(kind=dp)               ::  F
            real(kind=dp),dimension(:)  ::  USER
            real(kind=dp),dimension(DATA%GENEA%NP)  ,intent(OUT) :: FP            ! Likelihood by Family of Sire
            logical       ,dimension(DATA%GENEA%NP,N),intent(IN) :: FILTER_FP     ! Filter on The Family Sire
            real(kind=dp),dimension(DATA%GENEA%NM)  ,intent(OUT) :: FM            ! Likelihood by Family of SIre-Dam
            logical    ,dimension(DATA%GENEA%NP,DATA%GENEA%NM,N) ,intent(IN):: FILTER_FM  ! Filter on family Sire-Dam
            integer,intent(out)         ::  IFAIL
            external                    :: FUNCT1

#ifdef HAVE_LIBNAG
           integer,dimension(:),allocatable           ::  IW
           integer                                    ::  LIW, LW
           real(kind=dp),dimension(:),allocatable     ::  W

           LIW = N+2
           allocate (IW(LIW))
           LW  =12*N+N*N+1

           if ( LW < 13 ) LW = 13

           allocate (W(LW))
           IFAIL=1

           call e04jyf(N, IBOUND, FUNCT1, BL, BU, X, F, IW, LIW, W, LW,IUSER, USER, IFAIL)
!IFAIL = 1
!
!
!    On entry,	N<1,
!    or	IBOUND<0,
!    or	IBOUND>3,
!    or	IBOUND=0 and  BLj>BUj for some  j,
!    or	IBOUND=3 and BL1>BU1,
!    or	LIW<N+2,
!    or
!
!       LW
!       <
!       max⁡13,12×N+N×N-1/2
!      .
!
!IFAIL = 2
!
!
!    There have been 400×n function evaluations, yet the algorithm does not seem to be converging.  The calculations can be restarted from the final point held in X.  The error may also indicate that Fx has no minimum.
!
!IFAIL = 3
!
!
!    The conditions for a minimum have not all been met but a lower point could not be found and the algorithm has failed.
!
!IFAIL = 4
!
!
!    An overflow has occurred during the computation.  This is an unlikely failure, but if it occurs you should restart at the latest point given in X.
!
!IFAIL = 5
!
!IFAIL = 6
!
!IFAIL = 7
!
!IFAIL = 8
!
!
!    There is some doubt about whether the point x found by  E04JYF  is a minimum.  The degree of confidence in the result decreases as IFAIL increases.  Thus, when IFAIL =5 it is probable that the final x gives a good estimate of the position of a minimum, but when IFAIL =8 it is very unlikely that the routine has found a minimum.
!
!IFAIL = 9
!
!
!    In the search for a minimum, the modulus of one of the variables has become very large ∼106.  This indicates that there is a mistake in FUNCT1, that your problem has no finite solution, or that the problem needs rescaling (see  Section 8).
!
!IFAIL = 10
!
!
!    The computed set of forward-difference intervals (stored in W9×N+1,W9×N+2,…,
!    W10×N) is such that Xi+W9×N+i≤Xi for some i.
!
!    This is an unlikely failure, but if it occurs you should attempt to select another starting point.
!
           if ( IFAIL /= 0 ) then
              call log_mess('E04JYF IFAIL:'//str(IFAIL),WARNING_DEF)
           end if

           deallocate(IW)
           deallocate (W)

#else
           call stop_application("No NAG environment")
#endif

      end subroutine minimizing_e04jyf
!!***

!!****f* m_qtlmap_optimization/minimizing_lbfgs
!!  NAME
!!    minimizing_lbfgs
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
      subroutine minimizing_lbfgs(DATA,N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT1, & ! SUBROUTINE, supplied by the user. External Procedure
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array IUSER must be at least 1.
                                  FP,         &  ! Likelihood by Family of Sire
                                  FILTER_FP,  &  ! Filter on The Family Sire
                                  FM,         &  ! Likelihood by Family of SIre-Dam
                                  FILTER_FM,  &  ! Filter on family Sire-Dam
                                  IFAIL,  &
                                  OPTIM_FAMILY,&
                                  FUNCT_PART)

            type(QTLMAP_DATASET),target :: DATA
            integer                     :: N, IBOUND, IFAIL
            integer,dimension(:)        :: IUSER
            real(kind=dp),dimension(N)  :: BL, BU, X
            real(kind=dp)               :: F
            real(kind=dp),dimension(:)  :: USER
            real(kind=dp),dimension(DATA%GENEA%NP)  ,intent(OUT) :: FP            ! Likelihood by Family of Sire
            logical       ,dimension(DATA%GENEA%NP,N),intent(IN) :: FILTER_FP     ! Filter on The Family Sire
            real(kind=dp),dimension(DATA%GENEA%NM)  ,intent(OUT) :: FM            ! Likelihood by Family of SIre-Dam
            logical    ,dimension(DATA%GENEA%NP,DATA%GENEA%NM,N) ,intent(IN):: FILTER_FM  ! Filter on family Sire-Dam
            integer    ,intent(in)      :: OPTIM_FAMILY
            external                    :: FUNCT1,FUNCT_PART

!        nmax is the dimension of the largest problem to be solved.
!        mmax is the maximum number of limited memory corrections.

!     Declare the variables needed by the code.
!       A description of all these variables is given at the end of
!       the driver.

!     m is an INTEGER variable that must be set by the user to the
!       number of corrections used in the limited memory matrix.
!       It is not altered by the routine.  Values of m < 3  are
!       not recommended, and large values of m can result in excessive
!       computing time. The range  3 <= m <= 20 is recommended.

      integer ,parameter      :: m = 15
      character(len=60)       :: task, csave
      logical , dimension(4)  :: lsave
      integer                 :: iprint,nbd(n), iwa(3*n), isave(44)
      double precision        :: factr, pgtol
      double precision        :: l(n), u(n), g(n), dsave(29)
      double precision        :: wa(2*m*n+4*n+12*m*m+12*m)
      integer                 :: i,c,s1,s2
      double precision        :: t1,t2

      dsave=0.d0
      IFAIL = 0
      F=0.d0
!     We wish to have output at every iteration.
      s1=size(iuser)
      s2=size(user)
      iprint = -1

!     We specify the tolerances in the stopping criteria.
! Typical values for factr: 1.d+12 for
!         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
!         high accuracy.
      factr=1.0d+10 ! par defaut 1.0d+7
      pgtol=DATA%PARAMS%optim_tolg!5.0d-2 ! par defaut 1.0d-5
!     We specify the dimension n of the sample problem and the number
!        m of limited memory corrections stored.  (n and m should not
!        exceed the limits nmax and mmax respectively.)


!     We now provide nbd which defines the bounds on the variables:
!                    l   specifies the lower bounds,
!                    u   specifies the upper bounds.

      do i=1,n
         nbd(i)=2
         l(i)=BL(i)
         u(i)=BU(i)
      end do


!     We start the iteration by initializing task.
!
      task = 'START'
!        ------- the beginning of the loop ----------
      c=0
 111  continue
     
      if ( c > DATA%PARAMS%optim_maxeval  ) then
        call log_mess("** WARNING : LBFGS : Maximum iteration was exceeded without convergence **",WARNING_DEF)
        goto 112
      end if
!     This is the call to the L-BFGS-B code.
    !  write (*,fmt='(a1)',advance='no') '-'
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)

     !  write (*,fmt='(a5)',advance='no') '*'
      if (task(1:2) .eq. 'FG') then
!        the minimization routine has returned to request the
!        function f and gradient g values at the current x.
!        Compute function value f for the sample problem.
         if ( OPTIM_FAMILY == OPT_FAM) then
          call get_value_optim_fam(FUNCT_PART,n,x,f,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_FAM_MULTI) then
          call get_value_optim_fam_multi(FUNCT_PART,n,x,f,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_FAM_SIRE) then
          call get_value_optim_fam_sire(FUNCT_PART,n,x,f,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_CLASSIC) then
          call FUNCT1(n,x,f,iuser,user)
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(np+1)+count_call_funct
#endif
         end if

         c = c+1
!        Compute gradient g for the sample problem.
         if ( OPTIM_FAMILY == OPT_FAM ) then
           call get_gradient_optim_fam(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,DATA%GENEA%NP,&
           FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_FAM_MULTI ) then
           call get_gradient_optim_fam_multi(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,DATA%GENEA%NP,&
           FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_FAM_SIRE) then
          call get_gradient_optim_fam_sire(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,DATA%GENEA%NP,&
          FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM)
         else if ( OPTIM_FAMILY == OPT_CLASSIC) then
           call get_gradient(f,n,x,FUNCT1,s1,iuser,s2,user,u,l,g,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,&
           FM,FILTER_FM)
         end if
!        go back to the minimization routine.
         goto 111
!         print *,'ENCORE'
      endif

!
!      print *,'test'
      if (task(1:5) .eq. 'NEW_X')  goto 111
   !   call log_mess('LBFGD : task:'//trim(task),DEBUG_DEF)
!        the minimization routine has returned with a new iterate,
!         and we have opted to continue the iteration.

112   continue

!           ---------- the end of the loop -------------
!     If task is neither FG nor NEW_X we terminate execution.
      return
      end subroutine minimizing_lbfgs
!!***

!!****f* m_qtlmap_optimization/minimizing_nlopt
!!  NAME
!!    minimizing_nlopt
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!    Note: the dimension of the array IUSER and USER must be at least 1.
!!  SOURCE
       subroutine minimizing_nlopt(DATA,NLOPT_OPTI,&   !optimisateur nlopt
                                  N,      & ! the number n of independent variables. N >=1
                                  IBOUND, & ! indicates whether the facility for dealing with bounds of special forms is to be used.
                                  FUNCT1, & ! SUBROUTINE, supplied by the user. External Procedure
                                  BL,     & ! lower bounds
                                  BU,     & ! upper bound
                                  X,      & ! vector W
                                  F,      & ! final point stored in X
                                  IUSER,  & ! Note: the dimension of the array IUSER must be at least 1.
                                  USER,   & ! Note: the dimension of the array USER must be at least 1.
                                  FP,         &  ! Likelihood by Family of Sire
                                  FILTER_FP,  &  ! Filter on The Family Sire
                                  FM,         &  ! Likelihood by Family of SIre-Dam
                                  FILTER_FM,  &  ! Filter on family Sire-Dam
                                  IFAIL,  &
                                  OPTIM_FAMILY,&
                                  FUNCT_PART)

            type(QTLMAP_DATASET),target :: DATA
            integer,intent(in)          :: NLOPT_OPTI,N, IBOUND
            integer,dimension(:)        :: IUSER
            real(kind=dp),dimension(N)  :: BL, BU, X
            real(kind=dp)               :: F
            real(kind=dp),dimension(:)  :: USER
            real(kind=dp),dimension(DATA%GENEA%NP)  ,intent(OUT) :: FP            ! Likelihood by Family of Sire
            logical       ,dimension(DATA%GENEA%NP,N),intent(IN) :: FILTER_FP     ! Filter on The Family Sire
            real(kind=dp),dimension(DATA%GENEA%NM)  ,intent(OUT) :: FM            ! Likelihood by Family of SIre-Dam
            logical    ,dimension(DATA%GENEA%NP,DATA%GENEA%NM,N) ,intent(IN):: FILTER_FM  ! Filter on family Sire-Dam
            integer,intent(out)         :: IFAIL
            integer    ,intent(in)      :: OPTIM_FAMILY
            external                    :: FUNCT1,FUNCT_PART

#ifdef HAVE_NLOPT
            external :: interface_minimizing_nlopt_partial
#endif
            integer :: s1,s2

            s1=size(iuser)
            s2=size(user)

           if ( OPTIM_FAMILY == OPT_FAM ) then
#ifdef HAVE_NLOPT
            call interface_minimizing_nlopt_partial(NLOPT_OPTI,get_value_optim_fam,get_gradient_optim_fam,FUNCT_PART,n,x, &
                BL,BU,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM,F,&
                DATA%PARAMS%optim_tolx,DATA%PARAMS%optim_tolf,DATA%PARAMS%optim_maxeval,&
                 DATA%PARAMS%optim_maxtime,IFAIL)
           else if ( OPTIM_FAMILY == OPT_FAM_MULTI ) then
            call interface_minimizing_nlopt_partial(NLOPT_OPTI,get_value_optim_fam_multi,get_gradient_optim_fam_multi,&
                            FUNCT_PART,n,x,BL,BU,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM,F,&
                            DATA%PARAMS%optim_tolx,DATA%PARAMS%optim_tolf,DATA%PARAMS%optim_maxeval,DATA%PARAMS%optim_maxtime,IFAIL)
           else if ( OPTIM_FAMILY == OPT_FAM_SIRE ) then
              call interface_minimizing_nlopt_partial(NLOPT_OPTI,get_value_optim_fam_sire,get_gradient_optim_fam_sire,&
                            FUNCT_PART,n,x,BL,BU,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,&
                            FILTER_FM,F,DATA%PARAMS%optim_tolx,DATA%PARAMS%optim_tolf,&
                            DATA%PARAMS%optim_maxeval,DATA%PARAMS%optim_maxtime,IFAIL)

           else if ( OPTIM_FAMILY == OPT_CLASSIC ) then
              call interface_minimizing_nlopt_partial(NLOPT_OPTI,get_value,get_gradient,&
                             FUNCT1,n,x,BL,BU,s1,iuser,s2,user,DATA%GENEA%NP,FP,FILTER_FP,DATA%GENEA%NM,FM,FILTER_FM,F,&
                               DATA%PARAMS%optim_tolx,DATA%PARAMS%optim_tolf,&
                               DATA%PARAMS%optim_maxeval,DATA%PARAMS%optim_maxtime,IFAIL)
#else
           call stop_application("** Can not perform a nlopt optimization ** ")
#endif
           end if

        if ( IFAIL /= 0 ) then
         select case (IFAIL)
           case (-4)
             call log_mess("NLOPT_FAILURE",WARNING_DEF);!stop;
           case (-3)
             call log_mess("NLOPT_INVALID_ARGS",WARNING_DEF);!stop;
           case (-2)
             call log_mess("NLOPT_OUT_OF_MEMORY",WARNING_DEF);!stop;
           case (-1)
             call log_mess("NLOPT_ROUNDOFF_LIMITED",WARNING_DEF);!stop;
           case (0)
              call log_mess("NLOPT SUCCESS ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//") : "//trim(str(F)),INFO_DEF)
           case (1)
             call log_mess("NLOPT_MINF_MAX_REACHED ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//")  : "//trim(str(F)),VERBOSE_DEF)
           case (2)
             call log_mess("NLOPT_FTOL_REACHED ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//")  : "//trim(str(F)),VERBOSE_DEF)
           case (3)
             call log_mess("NLOPT_XTOL_REACHED ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//") : "//trim(str(F)),VERBOSE_DEF)
           case (4)
             call log_mess("NLOPT_MAXEVAL_REACHED ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//") : "//trim(str(F)),WARNING_DEF)
           case (5)
             call log_mess("NLOPT_MAXTIME_REACHED ("//trim(name_optim(NLOPT_OPTI+OPTI_LAST))//") : "//trim(str(F)),WARNING_DEF)
           case default
             call log_mess("THE METHOD FAILED :"//trim(str(ifail)),WARNING_DEF)
          end select
       end if

      end subroutine minimizing_nlopt
!!***

!***********************************************************************************************
!
!             FUNCTION FOR OPTIMISATION CLASSICAL
!
!


!!****f* m_qtlmap_optimization/get_gradient
!!  NAME
!!    get_gradient
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_gradient(f,n,x,FUNCT1,s1,iuser,s2,user,u,l,g,np,fp,filter_fp,nm,fm,filter_fm)
         real(kind=dp) ,intent(in)               :: f ! value in X
         integer ,intent(in)                     :: n,s1,s2,np,nm
         real(kind=dp) , dimension(n),intent(in) :: x,u,l
         real(kind=dp) ,intent(in),dimension(s2)  :: user
         integer   , intent(in),dimension(s1)     :: iuser
         real(kind=dp) , dimension(n),intent(out):: g ! gradient result
         real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
         logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm

        external                    :: FUNCT1

         !Local variable
         integer                         :: i
         double precision , dimension(n) :: xprime
         double precision                :: fplusH,fmoinH,h2
         logical                         :: bordGaucheOk,bordDroitOk
         character(len=1)                :: o
#ifdef BENCHMARK_VIEW
          count_call_gradient=1+count_call_gradient
#endif
         xprime = x
!         do i=1,n
!            print *,'x(',i,')=',x(i)
!         end do
!         stop
         ! Precision gradiant df/dx = f(x+h)-f(x-h)/2*h
         !calcul dérivé symetrique
         h2 = P_DATA%PARAMS%optim_H_PRECISION * 2
         g = 0.d0
         do i=1,n
            bordGaucheOk = l(i)+P_DATA%PARAMS%optim_H_PRECISION < x(i)
            bordDroitOk  = u(i)-P_DATA%PARAMS%optim_H_PRECISION > x(i)

           if ( bordGaucheOk .and. bordDroitOk ) then
	           xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION
	        !   print *,n
	         !  print *,xprime
	         !  print *,fmoinH
	           call FUNCT1(n,xprime,fmoinH,iuser,user)
	           xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
	           call FUNCT1(n,xprime,fplusH,iuser,user)
	           g(i)=(fplusH-fmoinH)/h2
	           xprime(i) = x(i)
#ifdef BENCHMARK_VIEW
          count_call_funct=2*nmp(np+1)+count_call_funct
#endif

           else if ( bordGaucheOk ) then
            !   print *,'-- bord droit --'
               xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION
               call FUNCT1(n,xprime,fmoinH,iuser,user)
               g(i)=(f-fmoinH)/P_DATA%PARAMS%optim_H_PRECISION
               xprime(i) = x(i)
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(np+1)+count_call_funct
#endif
           else if ( bordDroitOk ) then
               xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
	           call FUNCT1(n,xprime,fplusH,iuser,user)
	           g(i)=(fplusH-f)/P_DATA%PARAMS%optim_H_PRECISION
	           xprime(i) = x(i)
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(np+1)+count_call_funct
#endif
           end if
         end do

      end subroutine get_gradient
!!***


!!****f* m_qtlmap_optimization/get_value
!!  NAME
!!    get_value
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
     subroutine get_value(FUNCT1,n,x,f,s1,iuser,s2,user,np,fp,filter_fp,nm,fm,filter_fm)

        integer         , intent(in)                   :: n,s1,s2,np,nm
        real (kind=dp)      ,dimension(n), intent(in)  :: x
        real (kind=dp)  , intent(inout)                :: f
        integer ,       dimension(s1), intent(in)      :: iuser
        real (kind=dp)      ,dimension(s2), intent(in) :: user
        real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
        logical  ,dimension(np,n)   , intent(in)  :: filter_fp
        logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm
        external :: FUNCT1
        integer :: jm,ip
        real (kind=dp) :: ft

!        print *," ** get_value **"
!        stop
        call FUNCT1(n,x,f,iuser,user)

#ifdef BENCHMARK_VIEW
          count_call_funct=nm+count_call_funct
#endif
      end subroutine get_value
!!***

!***********************************************************************************************
!
!             FUNCTION FOR OPTIMISATION IP-JM
!
!

!!****f* m_qtlmap_optimization/get_value_optim_fam
!!  NAME
!!    get_value_optim_fam
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_value_optim_fam(FUNCT_PART,n,x,f,s1,iuser,s2,user,np,fp,filter_fp,nm,fm,filter_fm)

        integer         , intent(in)                   :: n,s1,s2,np,nm
        real (kind=dp)      ,dimension(n), intent(in)  :: x
        real (kind=dp)  , intent(inout)                :: f
        integer ,       dimension(s1), intent(inout)   :: iuser
        real (kind=dp)  ,dimension(s2), intent(inout)  :: user
        real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
        logical  ,dimension(np,n)   , intent(in)  :: filter_fp
        logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm

        external :: FUNCT_PART
        integer :: jm,ip

        f=0
        do ip=1,np
           do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
               call FUNCT_PART(ip,jm,n,x,fm(jm),iuser,user)
               f=fm(jm)+f
           end do
        end do


#ifdef BENCHMARK_VIEW
          count_call_funct=nm+count_call_funct
#endif

      end subroutine get_value_optim_fam
!!***

!!****f* m_qtlmap_optimization/get_gradient_optim_fam
!!  NAME
!!    get_gradient_optim_fam
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_gradient_optim_fam(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,np,fp,filter_fp,nm,fm,filter_fm)
         real(kind=dp) ,intent(in)                   :: f ! value in X
         integer ,intent(in)                         :: n,s1,s2,np,nm
         real(kind=dp) , dimension(n),intent(in)     :: x,u,l
         real(kind=dp) ,intent(inout),dimension(s1)  :: user
         integer   , intent(inout),dimension(s2)     :: iuser
         real(kind=dp) , dimension(n),intent(out)    :: g ! gradient result
         real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
         logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm

         external                    :: FUNCT_PART

         !Local variable
         integer                         :: i,ip,jm
         double precision , dimension(n) :: xprime
         double precision                :: fplusH,fmoinH,h2,ftemp
         logical                         :: bordGaucheOk,bordDroitOk
        ! print *,"get_gradient_optim_fam:",size(iuser),size(user)
#ifdef BENCHMARK_VIEW
          count_call_gradient=1+count_call_gradient
#endif
         xprime = x
         h2 = P_DATA%PARAMS%optim_H_PRECISION * 2
         g = 0.d0

         do i=1,n
            bordGaucheOk = l(i)+P_DATA%PARAMS%optim_H_PRECISION < x(i)
            bordDroitOk  = u(i)-P_DATA%PARAMS%optim_H_PRECISION > x(i)

           if ( bordGaucheOk .and. bordDroitOk ) then
	           xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION

               fmoinH=0
	           do ip=1,np
	            do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	   !          print *,"**>",ip,jm
	             if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la faimme ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,ftemp,iuser,user)
                  fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	             else
                  fmoinH=fmoinH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	             end if
	            end do
	           end do

	           xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
	           fplusH=0
               do ip=1,np
	            do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	             if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la faimme ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,ftemp,iuser,user)
                  fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
                  count_call_funct=1+count_call_funct
#endif
	             else
                  fplusH=fplusH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	             end if
	            end do
	           end do

	           g(i)=(fplusH-fmoinH)/h2
	           xprime(i) = x(i)

           else if ( bordGaucheOk ) then
            !   print *,'-- bord droit --'
               xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION
               fmoinH=0
	           do ip=1,np
	            do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	             if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,ftemp,iuser,user)
                  fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	             else
                  fmoinH=fmoinH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	             end if
	            end do
	           end do
               g(i)=(f-fmoinH)/P_DATA%PARAMS%optim_H_PRECISION
               xprime(i) = x(i)
           else if ( bordDroitOk ) then
               xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
               fplusH=0
               do ip=1,np
	            do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	             if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,ftemp,iuser,user)
                  fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	             else
                  fplusH=fplusH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	             end if
	            end do
	           end do

	           g(i)=(fplusH-f)/P_DATA%PARAMS%optim_H_PRECISION
	           xprime(i) = x(i)

           end if
         end do

      end subroutine get_gradient_optim_fam
!!***

!!****f* m_qtlmap_optimization/get_value_optim_fam_sire
!!  NAME
!!    get_value_optim_fam_sire
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_value_optim_fam_sire(FUNCT_PART,n,x,f,s1,iuser,s2,user,np,fp,filter_fp,nm,fm,filter_fm)

        integer         , intent(in)                   :: n,s1,s2,np,nm
        real (kind=dp)      ,dimension(n), intent(in)  :: x
        real (kind=dp)  , intent(inout)                :: f
        integer ,       dimension(s1), intent(inout)   :: iuser
        real (kind=dp)  ,dimension(s2), intent(inout)  :: user
        real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
        logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm
        external :: FUNCT_PART
        integer :: ip,i
        logical,dimension(N) :: keep_parameter

        f=0
        do ip=1,np
            call FUNCT_PART(ip,n,x,fp(ip),iuser,user)
            f=fp(ip)+f
        end do

#ifdef BENCHMARK_VIEW
          count_call_funct=nm+count_call_funct
#endif


      end subroutine get_value_optim_fam_sire
!!***

!!****f* m_qtlmap_optimization/get_gradient_optim_fam_sire
!!  NAME
!!    get_gradient_optim_fam_sire
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_gradient_optim_fam_sire(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,np,fp,filter_fp,nm,fm,filter_fm)

         real(kind=dp) ,intent(in)                   :: f ! value in X
         integer ,intent(in)                         :: n,s1,s2,np,nm
         real(kind=dp) , dimension(n),intent(in)     :: x,u,l
         real(kind=dp) ,intent(inout),dimension(s1)  :: user
         integer   , intent(inout),dimension(s2)     :: iuser
         real(kind=dp) , dimension(n),intent(out)    :: g ! gradient result
         real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
         logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm
         external                    :: FUNCT_PART

         !Local variable
         integer                         :: i,ip
         double precision , dimension(n) :: xprime
         double precision                :: fplusH,fmoinH,h2,ftemp
         logical                         :: bordGaucheOk,bordDroitOk
        ! print *,"get_gradient_optim_fam:",size(iuser),size(user)
#ifdef BENCHMARK_VIEW
          count_call_gradient=1+count_call_gradient
#endif
         xprime = x
         h2 = P_DATA%PARAMS%optim_H_PRECISION * 2
         g = 0.d0


         do i=1,n
            bordGaucheOk = l(i)+P_DATA%PARAMS%optim_H_PRECISION < x(i)
            bordDroitOk  = u(i)-P_DATA%PARAMS%optim_H_PRECISION > x(i)

           if ( bordGaucheOk .and. bordDroitOk ) then
	           xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION

               fmoinH=0

	           do ip=1,np
	             if ( filter_fp(ip,i) ) then ! cet element affecte la vraissemblance de la faimme ip
                  call FUNCT_PART(ip,n,xprime,ftemp,iuser,user)
                  fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(ip+1)-nmp(ip)+count_call_funct
#endif
	             else
                  fmoinH=fmoinH+fp(ip)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=nmp(ip+1)-nmp(ip)+count_avoid_funct
#endif
	             end if
	           end do


	           xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
	           fplusH=0
               do ip=1,np
	             if ( filter_fp(ip,i) ) then ! cet element affecte la vraissemblance de la faimme ip,jm
                  call FUNCT_PART(ip,n,xprime,ftemp,iuser,user)
                  fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
                  count_call_funct=nmp(ip+1)-nmp(ip)+count_call_funct
#endif
	             else
                  fplusH=fplusH+fp(ip)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=nmp(ip+1)-nmp(ip)+count_avoid_funct
#endif
	             end if
	           end do


	           g(i)=(fplusH-fmoinH)/h2
	           xprime(i) = x(i)

           else if ( bordGaucheOk ) then
            !   print *,'-- bord droit --'
               xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION
               fmoinH=0
	           do ip=1,np
	             if ( filter_fp(ip,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,n,xprime,ftemp,iuser,user)
                  fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(ip+1)-nmp(ip)+count_call_funct
#endif
	             else
                  fmoinH=fmoinH+fp(ip)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=nmp(ip+1)-nmp(ip)+count_avoid_funct
#endif
	             end if
	           end do
               g(i)=(f-fmoinH)/P_DATA%PARAMS%optim_H_PRECISION
               xprime(i) = x(i)
           else if ( bordDroitOk ) then
               xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
               fplusH=0
               do ip=1,np
	             if ( filter_fp(ip,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,n,xprime,ftemp,iuser,user)
                  fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=nmp(ip+1)-nmp(ip)+count_call_funct
#endif
	             else
                  fplusH=fplusH+fp(ip)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=nmp(ip+1)-nmp(ip)+count_avoid_funct
#endif
	             end if
	            end do

	           g(i)=(fplusH-f)/P_DATA%PARAMS%optim_H_PRECISION
	           xprime(i) = x(i)

           end if
         end do



      end subroutine get_gradient_optim_fam_sire
!!***

!!****f* m_qtlmap_optimization/get_inv_residual_covariance_matrix
!!  NAME
!!    get_inv_residual_covariance_matrix
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!   FUNCTION FOR OPTIMISATION IP-JM-MULTITRAITS
!!  SOURCE
    subroutine get_inv_residual_covariance_matrix(ip,n,x,vci,det,ok)
      integer         , intent(in)                   :: ip,n
      real (kind=dp)  ,dimension(n)   , intent(in)   :: x
      real(kind=dp),dimension(P_DATA%PHENOMODEL%ncar,P_DATA%PHENOMODEL%ncar) ,intent(out):: vci
      real(kind=dp)                      ,intent(out):: det
      logical                           ,intent(out) :: ok

      real(kind=dp),dimension(P_DATA%PHENOMODEL%ncar+1,P_DATA%PHENOMODEL%ncar)    :: VCIInverse

      integer :: ic,jc,k,irank
      real(kind=dp) :: rh

      ok=.true.

      !VCi : residual covariance matrix from sire I
      !         | SIGi1^2  SIGi1*SIGi2*RHO(1,2) ....  SIGi1*SIGip*RHO(1,p)
      !         | ..       SIGi2^2
      !    VCi  | ..                ...
      !         | ..                       ...
      !         | ..
      !         | SIGi1*SIGip*RHO(1,p)        SIGip^2
      k=0
      do ic=1,P_DATA%PHENOMODEL%ncar
        do jc=ic,P_DATA%PHENOMODEL%ncar
            VCIInverse(ic,jc)=x((ic-1)*P_DATA%GENEA%np+ip)*x((jc-1)*P_DATA%GENEA%np+ip)
            if ( ic /= jc ) then
              k=k+1
              !VCIInverse(ic,jc)=VCIInverse(ic,jc)*x(ncar*np+k)!rhoc(ic,jc)
              ! ** RESCALE OFI ***
              rh = x(P_DATA%PHENOMODEL%ncar*P_DATA%GENEA%np+k)
             ! print *,rh
              rh = ((dexp(rh)-1.d0)/(dexp(rh)+1.d0))
              VCIInverse(ic,jc)=VCIInverse(ic,jc)*rh
              VCIInverse(jc,ic)=VCIInverse(ic,jc)
            end if
        end do
      end do

      !Calcul Determinant de la matrice de covariance residuelle du pere ip
      !--------------------------------------------------------------------
      irank=0
      det=0
      !Calcul Inverse de la matrice de covariance residuelle du pere ip
      !-----------------------------------------------------------------
      CALL MATH_QTLMAP_INVDETMATSYM(P_DATA%PHENOMODEL%ncar,VCIInverse,P_DATA%PHENOMODEL%ncar+1,det,irank)
    !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
       if (irank /= 0) then
        ok=.false.
        return
      end if

      do ic=1,P_DATA%PHENOMODEL%ncar
        do jc=1,ic
            VCI(ic,jc)=VCIInverse(ic+1,jc)
            VCI(jc,ic)=VCI(ic,jc)
        end do
      end do



   end subroutine get_inv_residual_covariance_matrix
!!***

!!****f* m_qtlmap_optimization/get_value_optim_fam_multi
!!  NAME
!!    get_value_optim_fam_multi
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_value_optim_fam_multi(FUNCT_PART,n,x,f,s1,iuser,s2,user,np,fp,filter_fp,nm,fm,filter_fm)

        integer         , intent(in)                   :: n,s1,s2,np,nm
        real (kind=dp)      ,dimension(n), intent(in)  :: x
        real (kind=dp)  , intent(inout)                :: f
        integer ,       dimension(s1), intent(inout)   :: iuser
        real (kind=dp)  ,dimension(s2), intent(inout)  :: user
        real (kind=dp)  ,dimension(np)   , intent(inout)  :: fp
        real (kind=dp)  ,dimension(nm)   , intent(inout)  :: fm
        logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm

        external :: FUNCT_PART
        integer :: jm,ip
        logical :: valide
      !  real (kind=dp)  ,dimension(ncar,ncar) :: vct

        f=0
        do ip=1,np
           call get_inv_residual_covariance_matrix(ip,n,x,inv_vci(:,:,ip),determ(ip),valide)
#ifdef BENCHMARK_VIEW
          count_call_vci=1+count_call_vci
#endif
           if (.not. valide) then
              determ(ip)=0.d0
              f=INIFINY_REAL_VALUE
            !  fp=INIFINY_REAL_VALUE
              fm=INIFINY_REAL_VALUE
              return
           end if

           do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
               call FUNCT_PART(ip,jm,n,x,inv_vci(:,:,ip),determ(ip),fm(jm),iuser,user)
               f=fm(jm)+f
           end do
        end do

#ifdef BENCHMARK_VIEW
          count_call_funct=nm+count_call_funct
#endif

      end subroutine get_value_optim_fam_multi
!!***

!!****f* m_qtlmap_optimization/get_gradient_optim_fam_multi
!!  NAME
!!    get_gradient_optim_fam_multi
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
      subroutine get_gradient_optim_fam_multi(f,n,x,FUNCT_PART,s1,iuser,s2,user,u,l,g,np,fp,filter_fp,nm,fm,filter_fm)
         real(kind=dp) ,intent(in)                   :: f ! value in X
         integer ,intent(in)                         :: n,s1,s2,np,nm
         real(kind=dp) , dimension(n),intent(in)     :: x,u,l
         real(kind=dp) ,intent(inout),dimension(s1)  :: user
         integer   , intent(inout),dimension(s2)     :: iuser
         real(kind=dp) , dimension(n),intent(out)    :: g ! gradient result
         real (kind=dp)  ,dimension(np)   , intent(in)  :: fp
         real (kind=dp)  ,dimension(nm)   , intent(in)  :: fm
         logical  ,dimension(np,n)   , intent(in)  :: filter_fp
         logical  ,dimension(np,nm,n)   , intent(in)  :: filter_fm
         logical :: valide
         external                    :: FUNCT_PART

         !Local variable
         integer                         :: i,ip,jm,j
         real(kind=dp) , dimension(n) :: xprime
         real(kind=dp)                :: fplusH,fmoinH,h2,ftemp,det
         real(kind=dp) , dimension(:,:),pointer :: vc
         real(kind=dp) , dimension(P_DATA%PHENOMODEL%ncar,P_DATA%PHENOMODEL%ncar),target :: vct
         logical                      :: bordGaucheOk,bordDroitOk
        ! print *,"get_gradient_optim_fam:",size(iuser),size(user)
#ifdef BENCHMARK_VIEW
          count_call_gradient=1+count_call_gradient
#endif
         xprime = x
         h2 = P_DATA%PARAMS%optim_H_PRECISION * 2
         g = 0.d0

         do i=1,n
            bordGaucheOk = l(i)+P_DATA%PARAMS%optim_H_PRECISION < x(i)
            bordDroitOk  = u(i)-P_DATA%PARAMS%optim_H_PRECISION > x(i)

            if ( bordGaucheOk .and. bordDroitOk ) then
	           xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION

               fmoinH=0
	           do ip=1,np
	            valide=.true.
	            if (filter_vci(ip,i)) then
                  call get_inv_residual_covariance_matrix(ip,n,xprime,vct,det,valide)
#ifdef BENCHMARK_VIEW
          count_call_vci=1+count_call_vci
#endif
                  vc=>vct
                  if (.not. valide) then
                   fmoinH=INIFINY_REAL_VALUE
                  end if
	            else
#ifdef BENCHMARK_VIEW
          count_avoid_vci=1+count_avoid_vci
#endif
	              vc => inv_vci(:,:,ip)
	              det = determ(ip)
	              if ( det == 0.d0 ) then
                     fmoinH=INIFINY_REAL_VALUE
                     valide=.false.
                  end if
	            end if

                if (valide) then
	             do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	   !          print *,"**>",ip,jm
	              if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la faimme ip,jm
                   call FUNCT_PART(ip,jm,n,xprime,vc,det,ftemp,iuser,user)
                   fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	              else
                   fmoinH=fmoinH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	              end if
	             end do
	            end if
	           end do

	           xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
	           fplusH=0
               do ip=1,np
                valide=.true.
                if (filter_vci(ip,i)) then
                  call get_inv_residual_covariance_matrix(ip,n,xprime,vct,det,valide)
#ifdef BENCHMARK_VIEW
          count_call_vci=1+count_call_vci
#endif
                  vc=>vct
                  if (.not. valide) then
                   fmoinH=INIFINY_REAL_VALUE
                  end if
	            else
#ifdef BENCHMARK_VIEW
          count_avoid_vci=1+count_avoid_vci
#endif
	              vc => inv_vci(:,:,ip)
	              det = determ(ip)
	              if ( det == 0 ) then
                     fmoinH=INIFINY_REAL_VALUE
                     valide=.false.
                  end if
	            end if
	            if (valide) then
	             do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	              if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la faimme ip,jm
                   call FUNCT_PART(ip,jm,n,xprime,vc,det,ftemp,iuser,user)
                   fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
                   count_call_funct=1+count_call_funct
#endif
	              else
                   fplusH=fplusH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	              end if
	             end do
	            end if
	           end do

	           g(i)=(fplusH-fmoinH)/h2
	           xprime(i) = x(i)
	          !  print *,g(i)
             !  stop

           else if ( bordGaucheOk ) then
            !   print *,'-- bord droit --'
               xprime(i) = x(i)-P_DATA%PARAMS%optim_H_PRECISION
               fmoinH=0
	           do ip=1,np
	            valide=.true.
	            if (filter_vci(ip,i)) then
                  call get_inv_residual_covariance_matrix(ip,n,xprime,vct,det,valide)
#ifdef BENCHMARK_VIEW
          count_call_vci=1+count_call_vci
#endif
                  vc=>vct
                  if (.not. valide) then
                   fmoinH=INIFINY_REAL_VALUE
                  end if
	            else
#ifdef BENCHMARK_VIEW
          count_avoid_vci=1+count_avoid_vci
#endif
	              vc => inv_vci(:,:,ip)
	              det = determ(ip)
	            end if
	            if (valide) then
	             do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	              if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,vc,det,ftemp,iuser,user)
                  fmoinH=fmoinH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	             else
                  fmoinH=fmoinH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	              end if
	             end do
	            end if
	           end do
               g(i)=(f-fmoinH)/P_DATA%PARAMS%optim_H_PRECISION
               xprime(i) = x(i)
           else if ( bordDroitOk ) then
               xprime(i) = x(i)+P_DATA%PARAMS%optim_H_PRECISION
               fplusH=0
               do ip=1,np
                valide=.true.
                if (filter_vci(ip,i)) then
                  call get_inv_residual_covariance_matrix(ip,n,xprime,vct,det,valide)
#ifdef BENCHMARK_VIEW
          count_call_vci=1+count_call_vci
#endif
                  vc=>vct
                  if (.not. valide) then
                   fmoinH=INIFINY_REAL_VALUE
                  end if
	            else
#ifdef BENCHMARK_VIEW
          count_avoid_vci=1+count_avoid_vci
#endif
	              vc => inv_vci(:,:,ip)
	              det = determ(ip)
	            end if
	            if (valide) then
	             do jm=P_DATA%GENEA%nmp(ip)+1,P_DATA%GENEA%nmp(ip+1)
	              if ( filter_fm(ip,jm,i) ) then ! cet element affecte la vraissemblance de la famille ip,jm
                  call FUNCT_PART(ip,jm,n,xprime,vc,det,ftemp,iuser,user)
                  fplusH=fplusH+ftemp
#ifdef BENCHMARK_VIEW
          count_call_funct=1+count_call_funct
#endif
	              else
                  fplusH=fplusH+fm(jm)
#ifdef BENCHMARK_VIEW
          count_avoid_funct=1+count_avoid_funct
#endif
	              end if
	             end do
	            end if
	           end do

	           g(i)=(fplusH-f)/P_DATA%PARAMS%optim_H_PRECISION
	           xprime(i) = x(i)

           end if
         end do


      end subroutine get_gradient_optim_fam_multi
!!***

!!****f* m_qtlmap_optimization/likelihood_empty
!!  NAME
!!    likelihood_empty
!!  DESCRIPTION
!!
!!
!!
!!  NOTES
!!
!!  SOURCE
    subroutine likelihood_empty(ip,jm,n,x,f,iuser,user)
!      use OMP_LIB

      integer         , intent(in)                  :: ip,jm,n
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user


    end subroutine likelihood_empty
!!***



end module m_qtlmap_optimization
