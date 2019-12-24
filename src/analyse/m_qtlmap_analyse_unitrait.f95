!!****m* ANALYSE/m_qtlmap_analyse_unitrait
!!  NAME
!!    m_qtlmap_analyse_unitrait
!!  DESCRIPTION
!!
!!  NOTES
!!   optinit,opti_0qtl,opti_1qtl
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
module m_qtlmap_analyse_unitrait
     use m_qtlmap_types
     use m_qtlmap_log
     use m_qtlmap_optimization
     use m_qtlmap_analyse_gen
     implicit none
     save

     type(GENEALOGY_BASE) , pointer :: p_dg
     type(PDD_BUILD)      , pointer :: p_spt
     type(PHENOTYPE_BASE) , pointer :: p_dpa


!!****v* m_qtlmap_analyse_unitrait/std
!! NAME
!!  std
!! DESCRIPTION
!!  the standart deviation (residual) found under H1
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: std
!!****v* m_qtlmap_analyse_unitrait/xmoyp
!! NAME
!!  xmoyp
!! DESCRIPTION
!!  The polygenic sire effect found under H1
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: xmoyp
!!****v* m_qtlmap_analyse_unitrait/ap
!! NAME
!!  ap
!! DESCRIPTION
!!  The qtl sire effect found under H1
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: ap
     !$omp threadprivate (std,xmoyp,ap)

!!****v* m_qtlmap_analyse_unitrait/am
!! NAME
!!  am
!! DESCRIPTION
!!  The qtl dam effect found under H1
!! DIMENSIONS
!!  nm
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: am
!!****v* m_qtlmap_analyse_unitrait/xmoym
!! NAME
!!  xmoym
!! DESCRIPTION
!!  The polygenic dam effect found under H1
!! DIMENSIONS
!!  nm
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: xmoym
!!****v* m_qtlmap_analyse_unitrait/ap2
!! NAME
!!  ap2
!! DESCRIPTION
!!  The qtl sire effects found under H2, 1-->first qtl, 2-->2nd sqtl
!! DIMENSIONS
!!  np,2
!!***
     real (kind=dp)       ,dimension(:,:),allocatable,public :: ap2
     !$omp threadprivate (am,xmoym,ap2)

!!****v* m_qtlmap_analyse_unitrait/xmoyp2
!! NAME
!!  xmoyp2
!! DESCRIPTION
!!  TThe polygenic sire effect found under H2
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: xmoyp2
!!****v* m_qtlmap_analyse_unitrait/std2
!! NAME
!!  std2
!! DESCRIPTION
!!  the standart deviation (residual) found under H2
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: std2
!!****v* m_qtlmap_analyse_unitrait/am2
!! NAME
!!  am2
!! DESCRIPTION
!!  The qtl dam effects found under H2, 1-->first qtl, 2-->2nd sqtl
!! DIMENSIONS
!!  nm,2
!!***
     real (kind=dp)       ,dimension(:,:),allocatable,public :: am2
!!****v* m_qtlmap_analyse_unitrait/xmoym2
!! NAME
!!  xmoym2
!! DESCRIPTION
!!  The polygenic dam effect found under H2
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,public   :: xmoym2
     !$omp threadprivate (xmoyp2,std2,am2,xmoym2)


!!****v* m_qtlmap_analyse_unitrait/sig1
!! NAME
!!  sig1
!! DESCRIPTION
!!  the standart deviation (residual) found under H0
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: sig1
!!****v* m_qtlmap_analyse_unitrait/xmu1p
!! NAME
!!  xmu1p
!! DESCRIPTION
!!  The polygenic sire effect found under H0
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1p
!!****v* m_qtlmap_analyse_unitrait/xmu1m
!! NAME
!!  xmu1m
!! DESCRIPTION
!!  The polygenic dam effect found under H0
!! DIMENSIONS
!!  nm
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1m
     !$omp threadprivate (sig1,xmu1p,xmu1m)


!!****v* m_qtlmap_analyse_unitrait/fp0
!! NAME
!!  fp0
!! DESCRIPTION
!!  The likelihood by half sib family under H0
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),pointer    ,private   :: fp0
!!****v* m_qtlmap_analyse_unitrait/fm0
!! NAME
!!  fm0
!! DESCRIPTION
!!  The likelihood by full sib family under H0
!! DIMENSIONS
!!  nm
!!***
     real (kind=dp)       ,dimension(:),pointer    ,private   :: fm0
!!****v* m_qtlmap_analyse_unitrait/f0
!! NAME
!!  f0
!! DESCRIPTION
!!  The likelihood under H0
!!***
     real (kind=dp)                                ,private   :: f0
     !$omp threadprivate (fp0,fm0,f0)


!!****v* m_qtlmap_analyse_unitrait/fp1
!! NAME
!!  fp1
!! DESCRIPTION
!!  buffer to get the likelihood founded by half sib family under H1 at the current position
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),pointer   ,private   :: fp1
!!****v* m_qtlmap_analyse_unitrait/fm1
!! NAME
!!  fm1
!! DESCRIPTION
!!  buffer to get the likelihood founded by full sib family under H1 at the current position
!! DIMENSIONS
!!  nm
!!***
     real (kind=dp)       ,dimension(:),pointer   ,private   :: fm1
     !$omp threadprivate (fp1,fm1)


!!****v* m_qtlmap_analyse_unitrait/sompp
!! NAME
!!  sompp
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: sompp
!!****v* m_qtlmap_analyse_unitrait/sompm
!! NAME
!!  sompm
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: sompm
!!****v* m_qtlmap_analyse_unitrait/sompmy
!! NAME
!!  sompmy
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: sompmy
!!****v* m_qtlmap_analyse_unitrait/carpp
!! NAME
!!  carpp
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: carpp
!!****v* m_qtlmap_analyse_unitrait/carpm
!! NAME
!!  carpm
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: carpm
!!****v* m_qtlmap_analyse_unitrait/somppy
!! NAME
!!  somppy
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: somppy
!!****v* m_qtlmap_analyse_unitrait/somppm
!! NAME
!!  somppm
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: somppm
     !$omp threadprivate (sompp,sompm,sompmy,carpp,carpm,somppy,somppm)

!!****v* m_qtlmap_analyse_unitrait/somppdf
!! NAME
!!  somppdf
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: somppdf
!!****v* m_qtlmap_analyse_unitrait/carppdf
!! NAME
!!  carppdf
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf
!!****v* m_qtlmap_analyse_unitrait/somppydf
!! NAME
!!  somppydf
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1qtl
!! DIMENSIONS
!!
!!***
     real (kind=dp)       ,dimension(:),allocatable,private   :: somppydf

     !$omp threadprivate (somppdf,carppdf,somppydf)

!!****v* m_qtlmap_analyse_unitrait/fp_1
!! NAME
!!  fp_1
!! DESCRIPTION
!!  The maximum likelihood founded by sire family under H1
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),pointer,private       :: fp_1
!!****v* m_qtlmap_analyse_unitrait/fm_1
!! NAME
!!  fm_1
!! DESCRIPTION
!!  The maximum likelihood founded by full sib family under H1
!! DIMENSIONS
!!  np
!!***
     real (kind=dp)       ,dimension(:),pointer,private       :: fm_1
!!****v* m_qtlmap_analyse_unitrait/f_1
!! NAME
!!  f_1
!! DESCRIPTION
!!  The maximum likelihood founded under H1
!!***
     real (kind=dp)                                ,private   :: f_1
     !$omp threadprivate (fp_1,fm_1,f_1)

    ! 2QTL....
!!****v* m_qtlmap_analyse_unitrait/sompp1
!! NAME
!!  sompp1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp1
!!****v* m_qtlmap_analyse_unitrait/sompp2
!! NAME
!!  sompp2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp2
!!****v* m_qtlmap_analyse_unitrait/sompm1
!! NAME
!!  sompm1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompm1
!!****v* m_qtlmap_analyse_unitrait/sompm2
!! NAME
!!  sompm2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompm2
!!****v* m_qtlmap_analyse_unitrait/carpp1
!! NAME
!!  carpp1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carpp1
!!****v* m_qtlmap_analyse_unitrait/carpp2
!! NAME
!!  carpp2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carpp2
!!****v* m_qtlmap_analyse_unitrait/carpm1
!! NAME
!!  carpm1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carpm1
!!****v* m_qtlmap_analyse_unitrait/carpm2
!! NAME
!!  carpm2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carpm2
!!****v* m_qtlmap_analyse_unitrait/sompp1y
!! NAME
!!  sompp1y
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp1y
!!****v* m_qtlmap_analyse_unitrait/sompp2y
!! NAME
!!  sompp2y
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp2y
!!****v* m_qtlmap_analyse_unitrait/sompm1y
!! NAME
!!  sompm1y
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompm1y
!!****v* m_qtlmap_analyse_unitrait/sompm2y
!! NAME
!!  sompm2y
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompm2y
!!****v* m_qtlmap_analyse_unitrait/sompp1p2
!! NAME
!!  sompp1p2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:,:),allocatable,private :: sompp1p2
!!****v* m_qtlmap_analyse_unitrait/sompp1m1
!! NAME
!!  sompp1m1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp1m1
!!****v* m_qtlmap_analyse_unitrait/sompp1m2
!! NAME
!!  sompp1m2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:,:),allocatable,private :: sompp1m2
!!****v* m_qtlmap_analyse_unitrait/sompp2m1
!! NAME
!!  sompp2m1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:,:),allocatable,private :: sompp2m1
!!****v* m_qtlmap_analyse_unitrait/sompp2m2
!! NAME
!!  sompp2m2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: sompp2m2
!!****v* m_qtlmap_analyse_unitrait/sompm1m2
!! NAME
!!  sompm1m2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:,:),allocatable,private :: sompm1m2

      !$omp threadprivate (sompp1,sompp2,sompm1,sompm2,carpp1,carpp2,carpm1,carpm2,sompp1y,sompp2y,sompm1y,sompm2y)
      !$omp threadprivate (sompp1p2,sompp1m1,sompp1m2,sompp2m1,sompp2m2,sompm1m2)

!!****v* m_qtlmap_analyse_unitrait/somppdf1
!! NAME
!!  somppdf1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: somppdf1
!!****v* m_qtlmap_analyse_unitrait/somppdf2
!! NAME
!!  somppdf2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: somppdf2

!!****v* m_qtlmap_analyse_unitrait/carppdf11
!! NAME
!!  carppdf11
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf11
!!****v* m_qtlmap_analyse_unitrait/carppdf12
!! NAME
!!  carppdf12
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf12
!!****v* m_qtlmap_analyse_unitrait/carppdf22
!! NAME
!!  carppdf22
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf22
!!****v* m_qtlmap_analyse_unitrait/somppydf1
!! NAME
!!  somppydf1
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: somppydf1
!!****v* m_qtlmap_analyse_unitrait/somppydf2
!! NAME
!!  somppydf2
!! DESCRIPTION
!!  initialization : loop before the minimization of the likelihood : opti_1car_2qtl
!! DIMENSIONS
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: somppydf2

      !$omp threadprivate (somppdf1,somppdf2,carppdf11,carppdf12,carppdf22,somppydf1,somppydf2)

!!****v* m_qtlmap_analyse_unitrait/fp_2
!! NAME
!!  fp_2
!! DESCRIPTION
!!
!! DIMENSIONS
!!  np
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: fp_2

!!****v* m_qtlmap_analyse_unitrait/fm_2
!! NAME
!!  fm_2
!! DESCRIPTION
!!
!! DIMENSIONS
!!  nm
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: fm_2
      !$omp threadprivate (fp_2,fm_2)

!!****v* m_qtlmap_analyse_unitrait/fp01
!! NAME
!!  fp01
!! DESCRIPTION
!!
!! DIMENSIONS
!!  np
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: fp01
!!****v* m_qtlmap_analyse_unitrait/fm01
!! NAME
!!  fm01
!! DESCRIPTION
!!
!! DIMENSIONS
!!  nm
!!***
      real (kind=dp)       ,dimension(:),allocatable,private   :: fm01
      !$omp threadprivate (fp01,fm01)

!!****v* m_qtlmap_analyse_unitrait/current_ic
!! NAME
!!  current_ic
!! DESCRIPTION
!!  the current trait
!!***
      integer                                      , private   :: current_ic
!!****v* m_qtlmap_analyse_unitrait/current_chr
!! NAME
!!  current_chr
!! DESCRIPTION
!!  the current chromosome while the likelihood calculus
!!***
      integer                                      , private   :: current_chr
!!****v* m_qtlmap_analyse_unitrait/current_chr2
!! NAME
!!  current_chr2
!! DESCRIPTION
!!  the 2nd chromosome while the likelihood calculus (2QTL)
!!***
      integer                                      , private   :: current_chr2
      !$omp threadprivate (current_ic,current_chr,current_chr2)


     !---------------------------------------------------------------
     ! INTERFACE
     !---------------------------------------------------------------
     public :: init_analyse_unitrait_1QTL
     public :: init_analyse_unitrait_2QTL
     private :: init_sub_1qtl
     private :: end_sub_1qtl
     private :: init_sub_2qtl
     private :: end_sub_2qtl


     public :: opti_0qtl
     public :: opti_1qtl
     public :: opti_1car_2qtl
     public :: set_solution_hypothesis0
     public :: set_solution_hypothesis1
     public :: set_solution_hypothesis2
     public :: end_analyse_unitrait_1QTL
     public :: end_analyse_unitrait_2QTL


     contains
!!****f* m_qtlmap_analyse_unitrait/init_analyse_unitrait_1QTL
!! NAME
!!    init_analyse_unitrait_1QTL
!! DESCRIPTION
!!
!! SOURCE
       subroutine init_analyse_unitrait_1QTL(dataset)
           type(QTLMAP_DATASET)       ,intent(in) :: dataset
           integer           :: stat
           type(GENEALOGY_BASE) , pointer :: dg

           dg => dataset%genea

           allocate (sig1(dg%np),STAT=stat)
           call check_allocate(stat,'sig1 [m_qtlmap_analyse_unitrait]')
           sig1=0.D0
           allocate (xmu1p(dg%np),STAT=stat)
           call check_allocate(stat,'xmu1p [m_qtlmap_analyse_unitrait]')
           xmu1p=0.d0
           allocate (xmu1m(dg%nm),STAT=stat)
           call check_allocate(stat,'xmu1m [m_qtlmap_analyse_unitrait]')
           xmu1m=0.d0
           allocate (fp0(dg%np),STAT=stat)
           call check_allocate(stat,'fp0 [m_qtlmap_analyse_unitrait]')

           allocate (fm0(dg%nm),STAT=stat)
           call check_allocate(stat,'fm0 [m_qtlmap_analyse_unitrait]')

           allocate (ap(dg%np),STAT=stat)
           call check_allocate(stat,'ap [m_qtlmap_analyse_unitrait]')
           ap=0.d0
           allocate (xmoyp(dg%np),STAT=stat)
           call check_allocate(stat,'xmoyp [m_qtlmap_analyse_unitrait]')
           xmoyp=0.d0
           allocate (std(dg%np),STAT=stat)
           call check_allocate(stat,'std [m_qtlmap_analyse_unitrait]')
           std=0.d0
           allocate (am(dg%nm),STAT=stat)
           call check_allocate(stat,'am [m_qtlmap_analyse_unitrait]')
           am=0.d0
           allocate (xmoym(dg%nm),STAT=stat)
           call check_allocate(stat,'xmoym [m_qtlmap_analyse_unitrait]')
           xmoym=0.d0
           allocate (fp_1(dg%np),STAT=stat)
           call check_allocate(stat,'fp_1 [m_qtlmap_analyse_unitrait]')
           allocate (fm_1(dg%nm),STAT=stat)
           call check_allocate(stat,'fm_1 [m_qtlmap_analyse_unitrait]')

       end subroutine init_analyse_unitrait_1QTL
!!***

!!****f* m_qtlmap_analyse_unitrait/init_sub_1qtl
!! NAME
!!    init_sub_1qtl
!! DESCRIPTION
!!
!! SOURCE
        subroutine init_sub_1qtl(dataset,spt)
           type(QTLMAP_DATASET)       ,intent(in) :: dataset
           type(PDD_BUILD)            ,intent(in) :: spt

           integer           :: stat,maxng
           type(GENEALOGY_BASE) , pointer :: dg

           dg => dataset%genea

           maxng = spt%get_maxnbgenotypedam(dataset)

           allocate (sompp(maxng),STAT=stat)
           call check_allocate(stat,'sompp [m_qtlmap_analyse_unitrait]')
           allocate (sompm(maxng),STAT=stat)
           call check_allocate(stat,'sompm [m_qtlmap_analyse_unitrait]')
           allocate (sompmy(maxng),STAT=stat)
           call check_allocate(stat,'sompmy [m_qtlmap_analyse_unitrait]')
           allocate (carpp(maxng),STAT=stat)
           call check_allocate(stat,'carpp [m_qtlmap_analyse_unitrait]')
           allocate (carpm(maxng),STAT=stat)
           call check_allocate(stat,'carpm [m_qtlmap_analyse_unitrait]')
           allocate (somppy(maxng),STAT=stat)
           call check_allocate(stat,'somppy [m_qtlmap_analyse_unitrait]')
           allocate (somppm(maxng),STAT=stat)
           call check_allocate(stat,'somppm [m_qtlmap_analyse_unitrait]')
           allocate (somppdf(dg%np),STAT=stat)
           call check_allocate(stat,'somppdf [m_qtlmap_analyse_unitrait]')
           allocate (carppdf(dg%np),STAT=stat)
           call check_allocate(stat,'carppdf [m_qtlmap_analyse_unitrait]')
           allocate (somppydf(dg%np),STAT=stat)
           call check_allocate(stat,'somppydf [m_qtlmap_analyse_unitrait]')
           allocate (fp1(dg%np),STAT=stat)
           call check_allocate(stat,'fp1 [m_qtlmap_analyse_unitrait]')
           allocate (fm1(dg%nm),STAT=stat)
           call check_allocate(stat,'fm1 [m_qtlmap_analyse_unitrait]')

        end subroutine init_sub_1qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/end_sub_1qtl
!! NAME
!!    end_sub_1qtl
!! DESCRIPTION
!!
!! SOURCE
        subroutine end_sub_1qtl
           deallocate (sompp)
           deallocate (sompm)
           deallocate (sompmy)
           deallocate (carpp)
           deallocate (carpm)
           deallocate (somppy)
           deallocate (somppm)
           deallocate (somppdf)
           deallocate (carppdf)
           deallocate (somppydf)
           deallocate (fp1)
           deallocate (fm1)
        end subroutine end_sub_1qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/init_analyse_unitrait_2QTL
!! NAME
!!    init_analyse_unitrait_2QTL
!! DESCRIPTION
!!
!! SOURCE
        subroutine init_analyse_unitrait_2QTL(dataset)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         integer             :: stat,maxng
         type(GENEALOGY_BASE) , pointer :: dg

         dg => dataset%genea

         allocate (ap2(dg%np,2),STAT=stat)
         call check_allocate(stat,'ap2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (xmoyp2(dg%np),STAT=stat)
         call check_allocate(stat,'xmoyp2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (std2(dg%np),STAT=stat)
         call check_allocate(stat,'std2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (fp_2(dg%np),STAT=stat)
         call check_allocate(stat,'fp_2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (am2(dg%nm,2),STAT=stat)
         call check_allocate(stat,'am2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (xmoym2(dg%nm),STAT=stat)
         call check_allocate(stat,'xmoym2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (fm_2(dg%nm),STAT=stat)
         call check_allocate(stat,'fm_2 [m_qtlmap_analyse_unitrait_2qtl]')

     end subroutine init_analyse_unitrait_2QTL
!!***

!!****f* m_qtlmap_analyse_unitrait/init_sub_2qtl
!! NAME
!!    init_sub_2qtl
!! DESCRIPTION
!!
!! SOURCE
      subroutine init_sub_2qtl(dataset,spt)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         type(PDD_BUILD)            ,intent(in) :: spt

         integer           :: stat,maxng
         type(GENEALOGY_BASE) , pointer :: dg

         dg => dataset%genea

         maxng = spt%get_maxnbgenotypedam(dataset)

         allocate (sompp1(maxng),STAT=stat)
         call check_allocate(stat,'sompp1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp2(maxng),STAT=stat)
         call check_allocate(stat,'sompp2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompm1(maxng),STAT=stat)
         call check_allocate(stat,'sompm1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompm2(maxng),STAT=stat)
         call check_allocate(stat,'sompm2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carpp1(maxng),STAT=stat)
         call check_allocate(stat,'carpp1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carpp2(maxng),STAT=stat)
         call check_allocate(stat,'carpp2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carpm1(maxng),STAT=stat)
         call check_allocate(stat,'carpm1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carpm2(maxng),STAT=stat)
         call check_allocate(stat,'carpm2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp1y(maxng),STAT=stat)
         call check_allocate(stat,'sompp1y [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp2y(maxng),STAT=stat)
         call check_allocate(stat,'sompp2y [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompm1y(maxng),STAT=stat)
         call check_allocate(stat,'sompm1y [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompm2y(maxng),STAT=stat)
         call check_allocate(stat,'sompm2y [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp1p2(maxng,maxng),STAT=stat)
         call check_allocate(stat,'sompp1p2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp1m1(maxng),STAT=stat)
         call check_allocate(stat,'sompp1m1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp1m2(maxng,maxng),STAT=stat)
         call check_allocate(stat,'sompp1m2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp2m1(maxng,maxng),STAT=stat)
         call check_allocate(stat,'sompp2m1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompp2m2(maxng),STAT=stat)
         call check_allocate(stat,'sompp2m2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (sompm1m2(maxng,maxng),STAT=stat)
         call check_allocate(stat,'sompm1m2 [m_qtlmap_analyse_unitrait_2qtl]')
           allocate (somppdf1(dg%np),STAT=stat)
         call check_allocate(stat,'somppdf1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (somppdf2(dg%np),STAT=stat)
         call check_allocate(stat,'somppdf2 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carppdf11(dg%np),STAT=stat)
         call check_allocate(stat,'carppdf11 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carppdf12(dg%np),STAT=stat)
         call check_allocate(stat,'carppdf12 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (carppdf22(dg%np),STAT=stat)
         call check_allocate(stat,'carppdf22 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (somppydf1(dg%np),STAT=stat)
         call check_allocate(stat,'somppydf1 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (somppydf2(dg%np),STAT=stat)
         call check_allocate(stat,'somppydf2 [m_qtlmap_analyse_unitrait_2qtl]')

         allocate (fp01(dg%np),STAT=stat)
         call check_allocate(stat,'fp01 [m_qtlmap_analyse_unitrait_2qtl]')
         allocate (fm01(dg%nm),STAT=stat)
         call check_allocate(stat,'fm01 [m_qtlmap_analyse_unitrait_2qtl]')

      end subroutine init_sub_2qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/end_sub_2qtl
!! NAME
!!    end_sub_2qtl
!! DESCRIPTION
!!
!! SOURCE
      subroutine end_sub_2qtl
          deallocate (sompp1)
         deallocate (sompp2)
         deallocate (sompm1)
         deallocate (sompm2)
         deallocate (carpp1)
         deallocate (carpp2)
         deallocate (carpm1)
         deallocate (carpm2)
         deallocate (sompp1y)
         deallocate (sompp2y)
         deallocate (sompm1y)
         deallocate (sompm2y)
         deallocate (sompp1p2)
         deallocate (sompp1m1)
         deallocate (sompp1m2)
         deallocate (sompp2m1)
         deallocate (sompp2m2)
         deallocate (sompm1m2)
         deallocate (somppdf1)
         deallocate (somppdf2)
         deallocate (carppdf11)
         deallocate (carppdf12)
         deallocate (carppdf22)
         deallocate (somppydf1)
         deallocate (somppydf2)
         deallocate (fp01)
         deallocate (fm01)
     end subroutine end_sub_2qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/end_analyse_unitrait_1QTL
!! NAME
!!    end_analyse_unitrait_1QTL
!! DESCRIPTION
!!
!! SOURCE
       subroutine end_analyse_unitrait_1QTL
           deallocate (sig1)
           deallocate (xmu1p)
           deallocate (xmu1m)
           deallocate (fp0)
           deallocate (fm0)
           deallocate (ap)
           deallocate (xmoyp)
           deallocate (std)
           deallocate (xmoym)
           deallocate (am)
           deallocate (fp_1)
           deallocate (fm_1)
       end subroutine end_analyse_unitrait_1QTL
!!***

!!****f* m_qtlmap_analyse_unitrait/end_analyse_unitrait_2qtl
!! NAME
!!    end_analyse_unitrait_2qtl
!! DESCRIPTION
!!
!! SOURCE

        subroutine end_analyse_unitrait_2qtl
         deallocate (ap2)
         deallocate (xmoyp2)
         deallocate (std2)
         deallocate (fp_2)
         deallocate (am2)
         deallocate (xmoym2)
         deallocate (fm_2)

     end subroutine end_analyse_unitrait_2qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/opti_0qtl
!! NAME
!!    opti_0qtl
!! DESCRIPTION
!!    Calcul de la vraisemblance 0 QTL, 1 caractere
!! SOURCE
      subroutine opti_0qtl(dataset,ic)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
!
! Divers

      integer   , intent(in)                 :: ic
      integer iuser(1)
      double precision  ,dimension(:),allocatable :: borni,borns,par
      real (kind=dp) :: user(1)
      integer :: npar,ibound,ip,jm,ifail,nest,nestim
      logical  , dimension(:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa

      dpa => dataset%phenoAnimal
      p_dpa => dataset%phenoAnimal
      dg => dataset%genea
      p_dg => dataset%genea


!
      current_ic = ic
!******************************************************************************
! Parametres de maximisation
      npar=(2*dg%np)+dpa%nmumest(ic)

      allocate (borni(npar))
      allocate (borns(npar))
      allocate (par(npar))
      allocate (filter_inc(dg%np,npar))

      ibound=0
      filter_inc=.false.
      nestim=0
      do ip=1,dg%np
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
        filter_inc(ip,ip)=.true.
        borni(ip+dg%np)=XMU_MIN
        borns(ip+dg%np)=XMU_MAX
        filter_inc(ip,ip+dg%np)=.true.

        nest=0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if(dpa%estime(current_ic,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
             filter_inc(ip,2*dg%np+nestim)=.true.
            end if
           end if
        end do
      end do

      do jm=1,dpa%nmumest(ic)
        borni(2*dg%np+jm)=XMU_MIN
        borns(2*dg%np+jm)=XMU_MAX
      end do

!
! Point de depart
      do ip=1,dg%np
        par(ip)=sig0(ip)
        par(ip+dg%np)=xmu0p(ip)
      end do
      do jm=1,dpa%nmumest(ic)
        par(2*dg%np+jm)=xmu0m(jm)
      end do
!
! Optimisation de la vraisemblance
      ifail=1
     ! call minimizing_funct(dg,npar,ibound,funct_0qtl,borni,borns,par,f0,iuser,user,ifail)
      call minimizing_funct_family_sire(dataset,npar,ibound,funct_0qtl_family,filter_inc,&
         fp0,borni,borns,par,f0,iuser,user,ifail)

      do ip = 1,dg%np
        sig1(ip)=par(ip)
        xmu1p(ip)=par(ip+dg%np)
      end do
      do jm=1,dpa%nmumest(ic)
        xmu1m(jm)=par(2*dg%np+jm)
      end do
!
      deallocate (borni)
      deallocate (borns)
      deallocate (par)
      deallocate (filter_inc)

      end subroutine opti_0qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_0qtl
!! NAME
!!    funct_0qtl
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous H0
!! SOURCE
      subroutine funct_0qtl(n,x,f,iuser,user)
      use m_qtlmap_analyse_gen, only : somy,eff,cary,effdf,somydf,carydf,estmum,somcd,somcddf

      implicit none
      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      integer :: ip,nm1,nm2,jm
      integer :: nestim,nest

      real (kind=dp) :: sig,var,vpf,xmup,xmup2
      real (kind=dp) :: vdf,somxmu,xmum,xmum2,xmupm

!************nestim=sum(estmum(:ip-1))+nest******************************************************************
      f=0.d0
      nestim=0
      do ip=1,p_dg%np
        sig=x(ip)
        var=sig*sig
        xmup=x(ip+p_dg%np)
        xmup2=xmup*xmup
        vdf=carydf(ip)+somcddf(ip)*xmup2-2.d0*xmup*somydf(ip)
        vdf=0.5d0*vdf/var
        fp0(ip)=vdf+(dble(effdf(ip))*dlog(sig))
        f=f+fp0(ip)
  !      print *,f

        nm1=p_dg%nmp(ip)+1
        nm2=p_dg%nmp(ip+1)
        somxmu=0.d0
        nest=0
        do jm=nm1,nm2
          if(p_dpa%estime(current_ic,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
             ! print *,jm,nest,nestim,sum(estmum(:ip-1))+nest,estmum(ip)
              xmum=x(2*p_dg%np+nestim)
              somxmu=somxmu+xmum
 !             print *,somxmu
            else
              xmum=-somxmu
            end if
            xmum2=xmum*xmum
            xmupm=2.d0*(xmup+xmum)
            vpf=cary(jm)+somcd(jm)*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
            vpf=0.5d0*vpf/var
            fm0(jm)=vpf+(dble(eff(jm))*dlog(sig))
            f=f+fm0(jm)
            fp0(ip)=fp0(ip)+fm0(jm)
!            print *,ip,jm,fm0(jm),f
          !  stop
          else
            fm0(jm)=0.d0
          end if
        end do
      end do

      return
      end subroutine funct_0qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_0qtl_family
!! NAME
!!    funct_0qtl_family
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous H0
!! SOURCE
  subroutine funct_0qtl_family(ip,n,x,f,iuser,user)
      implicit none
      integer         , intent(in)                  :: ip,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      integer :: jm,nestim,nest

      real (kind=dp) :: sig,var,vpf,xmup,xmup2
      real (kind=dp) :: vdf,somxmu,xmum,xmum2,xmupm

      !  nest=count(estime(current_ic,nmp(ip)+1:jm))
      f=0.d0
      sig=x(ip)
      var=sig*sig
      xmup=x(ip+p_dg%np)
      xmup2=xmup*xmup
      vdf=carydf(ip)+somcddf(ip)*xmup2-2.d0*xmup*somydf(ip)
      vdf=0.5d0*vdf/var
      f=vdf+(dble(effdf(ip))*dlog(sig))
      somxmu=0.d0
      nest=0
      do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
        if(p_dpa%estime(current_ic,jm)) then
          nest=nest+1
          if(nest.le.estmum(ip)) then
              if (ip>1) then
               nestim=sum(estmum(:ip-1))+nest
              else
               nestim=nest
              end if
!              nestim=nestim+1
             ! print *,jm,nest,nestim,sum(estmum(:ip-1))+nest,estmum(ip)
              xmum=x(2*p_dg%np+nestim)
              somxmu=somxmu+xmum
             ! print *,somxmu
            else
              xmum=-somxmu
            end if
            xmum2=xmum*xmum
            xmupm=2.d0*(xmup+xmum)
            vpf=cary(jm)+somcd(jm)*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
            vpf=0.5d0*vpf/var
            fm0(jm)=vpf+(dble(eff(jm))*dlog(sig))
            f=f+fm0(jm)
          else
            fm0(jm)=0.d0
          end if
      end do

      end subroutine funct_0qtl_family
!!***

!!****f* m_qtlmap_analyse_unitrait/opti_1qtl
!! NAME
!!    opti_1qtl
!! DESCRIPTION
!!     Computing statistical test across the chromosome
!! HISTORY
!! 01/04/2009 : -- imported this function into this module
!!            : -- export all printing outside from this module
!! SOURCE
      subroutine opti_1qtl(dataset,spt,ic,lrtsol)
      type(QTLMAP_DATASET)       ,intent(in)              :: dataset
      type(PDD_BUILD),target     ,intent(in)              :: spt
      integer , intent(in)                                :: ic
      type(TYPE_LRT_SOLUTION)      , intent(inout)        :: lrtsol

!
! Divers
      integer :: iuser(1)
      real (kind=dp),dimension(:),allocatable :: borni,borns,par
      double precision    :: user(1)

      integer :: chr,ibound,n,ip,jm,ngeno1,ngeno2,ix,ilong,nm1,nm2,ig
      integer :: npar,nd1,nd2,kd,kkd,ifail,ii
      integer :: nestim,nest,ifem,indam,nprime,ntotal,indexchr(dataset%map%nchr),rang
      real (kind=dp) :: pp, pm, somxmu,f1,fsave0
      logical  , dimension(:,:),pointer        :: filter_inc
      real (kind=dp) , dimension(:,:,:),pointer  :: maxpar
      real (kind=dp) , dimension(:,:)  ,pointer  :: fmax
      real (kind=dp) , dimension(:,:,:),pointer  :: fpsave1
      real (kind=dp) , dimension(:,:,:),pointer  :: fmsave1

      !Trick pour recuperer la valeur de ces tableaux
      !La parallelisation sur les caracteres empeche la recopie de ces tableaux (deja instancier NTHREADS fois)
      !on recupere le pointeur de ces tableau avant la parallelisation sur le groupe de liaison
      !puis on initialise tous les thread avec ces valeurs sauvegarde
      real (kind=dp)       ,dimension(:),pointer :: save_carydf
      real (kind=dp)       ,dimension(:),pointer :: save_somcddf
      real (kind=dp)       ,dimension(:),pointer :: save_somydf
      integer              ,dimension(:),pointer :: save_effdf
      integer      ,dimension(:),pointer         :: save_estmum
      real (kind=dp)       ,dimension(:),pointer :: save_somcd
      real (kind=dp)       ,dimension(:),pointer :: save_somy
      real (kind=dp)       ,dimension(:),pointer :: save_cary
      integer      ,dimension(:),pointer         :: save_eff
      real (kind=dp)       ,dimension(:),pointer :: save_fp0
      real (kind=dp)       ,dimension(:),pointer :: save_fm0
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      !type(DATAMODEL_BASE) , pointer :: dpm

      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp,xlrm
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

      !dpm => dataset%phenoModel

      dg => dataset%genea
      dpa => dataset%phenoAnimal
      p_dg => dataset%genea
      p_dpa => dataset%phenoAnimal
      p_spt => spt

      call lrtsol%new(dataset,1)

      allocate (lrt1(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (xlrp(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))

!******************************************************************************
! Calcul de la vraisemblance sous H1
! Parametres de maximisation

      ibound=0
      npar=(3*dg%np)+dpa%nmumest(ic)+dpa%namest(ic)

      allocate (filter_inc(dg%np,npar))
      allocate (borni(npar),borns(npar),par(npar))

      save_carydf => carydf
      save_somcddf => somcddf
      save_somydf => somydf
      save_effdf => effdf
      save_estmum => estmum
      save_somcd => somcd
      save_somy => somy
      save_cary => cary
      save_eff => eff
      fsave0 = f0
      save_fp0 => fp0
      save_fm0 => fm0

      filter_inc=.false.
      nestim=0
      do ip=1,dg%np
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
        filter_inc(ip,ip)=.true.
        borni(ip+dg%np)=XMU_MIN
        borns(ip+dg%np)=XMU_MAX
        filter_inc(ip,ip+dg%np)=.true.
        borni(2*dg%np+ip)=AP_MIN
        borns(2*dg%np+ip)=AP_MAX
        filter_inc(ip,2*dg%np+ip)=.true.

        nest=0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if(dpa%estime(current_ic,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
              filter_inc(ip,3*dg%np+nestim)=.true.
            end if
            ifem=dg%repfem(jm)
            indam=dpa%iam(ic,ifem)
            filter_inc(ip,3*dg%np+dpa%nmumest(ic)+indam)=.true.
           end if
        end do
      end do
      !filter_inc(:,3*np+nmumest(ic):)=.true.

      do jm=1,dpa%nmumest(ic)
        borni(3*dg%np+jm)=XMU_MIN
        borns(3*dg%np+jm)=XMU_MAX
      end do
      do jm=1,dpa%namest(ic)
        borni(3*dg%np+dpa%nmumest(ic)+jm)=AM_MIN
        borns(3*dg%np+dpa%nmumest(ic)+jm)=AM_MAX
      end do
!
! Point de depart
      do ip=1,dg%np
        par(ip)=sig1(ip)
        par(ip+dg%np)=xmu1p(ip)
        par(2*dg%np+ip)=0.d0
      end do


      do jm=1,dpa%nmumest(ic)
        par(3*dg%np+jm)=xmu1m(jm)
      end do
      do jm=1,dpa%namest(ic)
        par(3*dg%np+dpa%nmumest(ic)+jm)=0.d0
      end do
!
! Marche le long du chromosome

      lrtsol%lrtmax(0)=-1.d75
      lrtsol%chrmax(0)=1

      ntotal=0
      do chr=1,dataset%map%nchr
        ntotal=ntotal+dataset%map%get_npo(chr)
        indexchr(chr)=ntotal
      end do

      allocate (maxpar(dataset%map%nchr,dataset%map%get_maxnpo(),npar))
      allocate (fmax(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (fpsave1(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
      allocate (fmsave1(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))

      !$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(par)  &
      !$OMP PRIVATE(ip,jm,chr,n,nd1,nd2,pp,pm,ngeno1,ngeno2,kd,kkd,ig,ifail,nestim,ii,f1)
      current_ic = ic
      call init_sub_1qtl(dataset,spt)
      carydf=>save_carydf
      somcddf => save_somcddf
      somydf => save_somydf
      effdf => save_effdf
      estmum => save_estmum
      somcd => save_somcd
      somy => save_somy
      cary => save_cary
      eff => save_eff
      f0 = fsave0
      fp0 => save_fp0
      fm0 => save_fm0

      !$OMP DO
      do nprime=1,ntotal
       n=nprime
       do chr=dataset%map%nchr,1,-1
        if ( indexchr(chr) >= nprime) then
          current_chr = chr
        else
          n=n-dataset%map%get_npo(chr)
        end if
       end do
       chr=current_chr
!       print *,chr,n
       do ip=1,dg%np
          nm1=dg%nmp(ip)+1
          nm2=dg%nmp(ip+1)
          somppdf(ip)=0.d0
          carppdf(ip)=0.d0
          somppydf(ip)=0.d0
          do jm=nm1,nm2
            ngeno1=spt%ngenom(chr,jm)+1
            ngeno2=spt%ngenom(chr,jm+1)
            do ig=ngeno1,ngeno2
              sompp(ig)=0.d0
              carpp(ig)=0.d0
              somppy(ig)=0.d0
              sompm(ig)=0.d0
              carpm(ig)=0.d0
              sompmy(ig)=0.d0
              somppm(ig)=0.d0
              nd1=spt%ngend(chr,ig)+1
              nd2=spt%ngend(chr,ig+1)
              do kd=nd1,nd2
                kkd=spt%ndesc(chr,kd)
                if(dpa%presentc(ic,kkd)) then
                  pp=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  pm=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  if(dpa%estime(ic,jm)) then
                    sompp(ig)=sompp(ig)+pp*dpa%cd(ic,kkd)
                    carpp(ig)=carpp(ig)+pp*pp*dpa%cd(ic,kkd)
                    somppy(ig)=somppy(ig)+pp*dpa%y(ic,kkd)*dpa%cd(ic,kkd)
                    sompm(ig)=sompm(ig)+pm*dpa%cd(ic,kkd)
                    carpm(ig)=carpm(ig)+pm*pm*dpa%cd(ic,kkd)
                    sompmy(ig)=sompmy(ig)+pm*dpa%y(ic,kkd)*dpa%cd(ic,kkd)
                    somppm(ig)=somppm(ig)+pp*pm*dpa%cd(ic,kkd)
                  else
                    somppdf(ip)=somppdf(ip)+pp*dpa%cd(ic,kkd)
                    carppdf(ip)=carppdf(ip)+pp*pp*dpa%cd(ic,kkd)
                    somppydf(ip)=somppydf(ip)+pp*dpa%y(ic,kkd)*dpa%cd(ic,kkd)
                  end if
                end if
              end do
            end do
          end do
        end do

!
! Optimisation de la vraisemblance a la position dx
        ifail=1
        !call minimizing_funct(npar,ibound,funct_1qtl,borni,borns,par,f1,iuser,user,ifail)

        call minimizing_funct_family_sire(dataset,npar,ibound,funct_1qtl_family,filter_inc,&
          fp1,borni,borns,par,f1,iuser,user,ifail)

        fpsave1(chr,n,:) = fp1
        fmsave1(chr,n,:) = fm1
        fmax(chr,n) = f1
        if ( f1 < INIFINY_REAL_VALUE ) then
         lrt1(chr,n)=-2.d0*(f1-f0)
        else
         lrt1(chr,n)=0
        end if

        do ii=1,dg%np
          xlrp(chr,n,ii)=-2.d0*(fp1(ii)-fp0(ii))
          lrtsol%pater_eff(chr,ii,n)=par(2*dg%np+ii)
        end do

        do ii=1,dg%nm
          xlrm(chr,n,ii)=-2.d0*(fm1(ii)-fm0(ii))
        end do

        do ii=1,dpa%namest(ic)
          lrtsol%mater_eff(chr,ii,n)=par(3*dg%np+dpa%nmumest(ic)+ii)
        end do

        maxpar(chr,n,:)=par
      end do
      !$OMP END DO
      call end_sub_1qtl
      !$OMP END PARALLEL
      do chr=1,dataset%map%nchr
        do n=1,dataset%map%get_npo(chr)
          if (lrt1(chr,n)> lrtsol%lrtmax(0)) then
             par=maxpar(chr,n,:)
             lrtsol%lrtmax(0)=lrt1(chr,n)
             lrtsol%nxmax(0)=n
             lrtsol%chrmax(0)=chr
             f_1=fmax(chr,n)
             nestim=0
             do ip=1,dg%np
              ap(ip)=par(2*dg%np+ip)
              xmoyp(ip)=par(dg%np+ip)
              std(ip)=par(ip)
              fp_1(ip)=fpsave1(chr,n,ip)
              nest=0
              somxmu=0.d0
              do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
              fm_1(jm)=fmsave1(chr,n,jm)
              if (dpa%estime(ic,jm)) then
                nest=nest+1
                ifem=dg%repfem(jm)
                indam=dpa%iam(ic,ifem)
                am(jm)=par(3*dg%np+dpa%nmumest(ic)+indam)
                if (nest.le.estmum(ip)) then
                  nestim=nestim+1
                  xmoym(jm)=par(3*dg%np+nestim)
                  somxmu=somxmu+xmoym(jm)
                else
                  xmoym(jm)=-somxmu
                end if
              else
                am(jm)=0.d0
                xmoym(jm)=0.d0
              end if
            end do
          end do
          end if
        end do
      end do

      ! 04/2013 - New structure for LRT
      call lrtsol%LRT%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),lrt1)
      do ii=1,dg%np
       call lrtsol%LRT_SIRES(ii)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrp(:,:,ii))
      end do
      do ii=1,dg%nm
       call lrtsol%LRT_DAMS(ii)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrm(:,:,ii))
      end do

      deallocate (lrt1)
      deallocate (xlrp)
      deallocate (xlrm)


      deallocate (fpsave1)
      deallocate (fmsave1)
      deallocate (fmax)
      deallocate (maxpar)
      deallocate (filter_inc)
      deallocate (borni,borns,par)

      end subroutine opti_1qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_1qtl
!! NAME
!!    funct_1qtl
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous Hypothese 1QTL 1carac
!! HISTORY
!!
!! SOURCE
      subroutine funct_1qtl(n,x,f,iuser,user)

      integer       , intent(in)                 :: n    ! number of variables
      double precision, dimension(n), intent(in) :: x ! Position
      double precision, intent(inout)            :: f    ! Result value at the position
      double precision,dimension(1)              :: user
      integer         ,dimension(1)              :: iuser


!
! Divers
      integer :: ip,nestim,nm1,nm2,nest,jm,ifem,indam,i
      integer :: ngeno1,ngeno2,ig
      real (kind=dp)  :: sig,var,xmup,xmup2,vdf,somxmu,xmum,xmum2
      real (kind=dp)  :: xmupm,vmere,z,vpf,x_ap,x_ap2,x_am,x_am2


!******************************************************************************
      f=0.d0
      nestim=0
      do ip=1,p_dg%np
         sig=x(ip)
         var=sig*sig
         xmup=x(ip+p_dg%np)
         xmup2=xmup*xmup
         x_ap=x(2*p_dg%np+ip)
         x_ap2=x_ap*x_ap
         vdf=carydf(ip)+somcddf(ip)*xmup2-2.d0*xmup*somydf(ip)
         vdf=vdf+x_ap2*carppdf(ip)+2.d0*xmup*x_ap*somppdf(ip)-2.d0*x_ap*somppydf(ip)
         vdf=0.5d0*vdf/var
         fp1(ip)=vdf+dble(effdf(ip))*dlog(sig)
         f=f+fp1(ip)
         nm1=p_dg%nmp(ip)+1
         nm2=p_dg%nmp(ip+1)
         somxmu=0.d0
         nest=0
         do jm=nm1,nm2
            if(p_dpa%estime(current_ic,jm)) then
               nest=nest+1
               if(nest.le.estmum(ip)) then
                  nestim=nestim+1
                  xmum=x(3*p_dg%np+nestim)
                  somxmu=somxmu+xmum
               else
                  xmum=-somxmu
               end if
               xmum2=xmum*xmum
               xmupm=2.d0*(xmup+xmum)
               ifem=p_dg%repfem(jm)
               indam=p_dpa%iam(current_ic,ifem)
               x_am=x(3*p_dg%np+p_dpa%nmumest(current_ic)+indam)
               x_am2=x_am*x_am
               z=cary(jm)+somcd(jm)*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
               vmere=0.d0
               ngeno1=p_spt%ngenom(current_chr,jm)+1
               ngeno2=p_spt%ngenom(current_chr,jm+1)

               do ig=ngeno1,ngeno2
                  vpf=z+x_ap2*carpp(ig)+x_am2*carpm(ig)            &
                      +2.d0*x_ap*x_am*somppm(ig)                   &
                      +xmupm*(x_ap*sompp(ig)+x_am*sompm(ig))       &
                      -2.d0*x_ap*somppy(ig)-2.d0*x_am*sompmy(ig)
                  vpf=-0.5d0*vpf/var
                  vmere=vmere+p_spt%probg(current_chr,ig)*dexp(vpf)

               end do

               if (vmere == 0) then
                  fm1(jm)=INIFINY_REAL_VALUE
               else
                 fm1(jm)=-dlog(vmere)+(dble(eff(jm))*dlog(sig))
               end if

    !         end if
               f=f+fm1(jm)
               fp1(ip)=fp1(ip)+fm1(jm)
            else
               fm1(jm)=0.d0
            end if
         end do
      end do
      return
      end subroutine funct_1qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_1qtl_family
!! NAME
!!    funct_1qtl_family
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous Hypothese 1QTL 1carac
!! HISTORY
!!
!! SOURCE
   subroutine funct_1qtl_family(ip,n,x,f,iuser,user)

      integer       , intent(in)                 :: ip,n    ! number of variables
      double precision, dimension(n), intent(in) :: x ! Position
      double precision, intent(inout)            :: f    ! Result value at the position
      double precision,dimension(1)              :: user
      integer         ,dimension(1)              :: iuser


!
! Divers
      integer :: nestim,nest,jm,ifem,indam,i
      integer :: ngeno1,ngeno2,ig
      real (kind=dp)  :: sig,var,xmup,xmup2,vdf,somxmu,xmum,xmum2
      real (kind=dp)  :: xmupm,vmere,z,vpf,x_ap,x_ap2,x_am,x_am2


!******************************************************************************
      f=0.d0
      nestim=0
      sig=x(ip)
      var=sig*sig
      xmup=x(ip+p_dg%np)
      xmup2=xmup*xmup
      x_ap=x(2*p_dg%np+ip)
      x_ap2=x_ap*x_ap
      vdf=carydf(ip)+somcddf(ip)*xmup2-2.d0*xmup*somydf(ip)
      vdf=vdf+x_ap2*carppdf(ip)+2.d0*xmup*x_ap*somppdf(ip)-2.d0*x_ap*somppydf(ip)
      vdf=0.5d0*vdf/var
      f=vdf+dble(effdf(ip))*dlog(sig)
      somxmu=0.d0
      nest=0
      do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
         if(p_dpa%estime(current_ic,jm)) then
               nest=nest+1
               if(nest.le.estmum(ip)) then
                  if (ip>1) then
                     nestim=sum(estmum(:ip-1))+nest
                  else
                     nestim=nest
                  end if
                  xmum=x(3*p_dg%np+nestim)
                  somxmu=somxmu+xmum
               else
                  xmum=-somxmu
               end if
               xmum2=xmum*xmum
               xmupm=2.d0*(xmup+xmum)
               ifem=p_dg%repfem(jm)
               indam=p_dpa%iam(current_ic,ifem)
               x_am=x(3*p_dg%np+p_dpa%nmumest(current_ic)+indam)
               x_am2=x_am*x_am
               z=cary(jm)+somcd(jm)*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
               vmere=0.d0
               ngeno1=p_spt%ngenom(current_chr,jm)+1
               ngeno2=p_spt%ngenom(current_chr,jm+1)

               do ig=ngeno1,ngeno2
                  vpf=z+x_ap2*carpp(ig)+x_am2*carpm(ig)            &
                      +2.d0*x_ap*x_am*somppm(ig)                   &
                      +xmupm*(x_ap*sompp(ig)+x_am*sompm(ig))       &
                      -2.d0*x_ap*somppy(ig)-2.d0*x_am*sompmy(ig)
                  vpf=-0.5d0*vpf/var
                  vmere=vmere+p_spt%probg(current_chr,ig)*dexp(vpf)
               end do

               if (vmere == 0) then
                  fm1(jm)=INIFINY_REAL_VALUE
               else
                 fm1(jm)=-dlog(vmere)+(dble(eff(jm))*dlog(sig))
               end if
               f=f+fm1(jm)
            else
               fm1(jm)=0.d0
            end if
      end do

      end subroutine funct_1qtl_family
!!***

!!****f* m_qtlmap_analyse_unitrait/opti_1car_2qtl
!! NAME
!!    opti_1car_2qtl
!! DESCRIPTION
!!     Computing statistical test across the chromosome
!! HISTORY
!!  MODIF OFI 01/04/2009 : -- imported this function into this module
!!                       : -- export all printing outside from this module
!! SOURCE
      subroutine opti_1car_2qtl(dataset,spt,ic,lrtsol)
      type(QTLMAP_DATASET)       ,intent(in)         :: dataset
      type(PDD_BUILD)  ,target   ,intent(in)         :: spt

      integer                             , intent(in)      :: ic
      type(TYPE_LRT_SOLUTION)            , intent(inout)    :: lrtsol

!
! Divers
      real (kind=dp) :: par(4*dataset%genea%np+3*dataset%genea%nm)
      real (kind=dp) :: borni(4*dataset%genea%np+3*dataset%genea%nm)
      real (kind=dp) :: borns(4*dataset%genea%np+3*dataset%genea%nm)
      integer        :: iuser(1),chr
      real (kind=dp) :: user(1)
      real (kind=dp)  :: pddp(dataset%genea%nd*4,2),pddm(dataset%genea%nd*4,2)

      integer :: ip,ibound,ngeno1,ngeno2,ig,nd1,nd2,kd,kkd,ilong,ntotal,nprime,indexchr(dataset%map%nchr)
      integer :: npar,jm,nestim,nm1,nm2,nest,ifem,indam,n,ifail,kd2,nstart
      integer :: n1,ix1,ii,iq,npar2,chr2,ngeno21,ngeno22,ilong2,ig1,ig2,nd11,nd12,nd21,nd22
      real (kind=dp) :: xmin,xmax,somxmu,f01_t,save_f_1,save_f0
      logical  , dimension(:,:),pointer        :: filter_inc
      real (kind=dp) , dimension(:,:,:,:,:),pointer  :: maxpar
      real (kind=dp) , dimension(:,:,:,:,:),pointer  :: fpsave1
      real (kind=dp) , dimension(:,:,:,:,:),pointer  :: fmsave1

      !Trick pour recuperer la valeur de ces tableaux
      !La parallelisation sur les caracteres empeche la recopie de ces tableaux (deja instancier NTHREADS fois)
      !on recupere le pointeur de ces tableau avant la parallelisation sur le groupe de liaison
      !puis on initialise tous les thread avec ces valeurs sauvegarde
      real (kind=dp)       ,dimension(:),pointer :: save_carydf
      real (kind=dp)       ,dimension(:),pointer :: save_somcddf
      real (kind=dp)       ,dimension(:),pointer :: save_somydf
      integer              ,dimension(:),pointer :: save_effdf
      integer      ,dimension(:),pointer         :: save_estmum
      real (kind=dp)       ,dimension(:),pointer :: save_somcd
      real (kind=dp)       ,dimension(:),pointer :: save_somy
      real (kind=dp)       ,dimension(:),pointer :: save_cary
      integer      ,dimension(:),pointer         :: save_eff
      real (kind=dp)       ,dimension(:),pointer :: save_fp_1,save_fp0
      real (kind=dp)       ,dimension(:),pointer :: save_fm_1,save_fm0
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa

      real(kind=dp)       ,dimension(:,:,:,:,:)   ,pointer    :: xlrp2,xlrm2,xlrp2_0,xlrm2_0
      real(kind=dp)       ,dimension(:,:,:,:)   ,pointer      :: lrt0_2,lrt1_2

      dg => dataset%genea
      dpa => dataset%phenoAnimal
      p_dg => dataset%genea
      p_dpa => dataset%phenoAnimal
      p_spt => spt

!
!****************************************************************************
!
      call lrtsol%new(dataset,2)
      allocate (lrt0_2(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo()))
      allocate (lrt1_2(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo()))
      allocate (xlrp2(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm2(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%nm))
      allocate (xlrp2_0(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm2_0(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%nm))

      npar2=4*dg%np+3*dg%nm

      current_ic = ic

      borni=0.d0
      borns=0.d0
      par=0.d0

      save_carydf => carydf
      save_somcddf => somcddf
      save_somydf => somydf
      save_effdf => effdf
      save_estmum => estmum
      save_somcd => somcd
      save_somy => somy
      save_cary => cary
      save_eff => eff
      save_fp_1 => fp_1
      save_fm_1 => fm_1
      save_fp0 => fp0
      save_fm0 => fm0
      save_f_1 = f_1
      save_f0  = f0
!
! Calcul de la vraisemblance sous H1
! Parametres de maximisation
      npar=(4*dg%np)+dpa%nmumest(ic)+(2*dpa%namest(ic))
      ibound=0
      allocate (filter_inc(dg%np,npar))

      filter_inc=.false.
      nestim=0
      do ip=1,dg%np
        xmin=-5.d0*std(ip)
        xmax=-xmin
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
        filter_inc(ip,ip)=.true.
        borni(ip+dg%np)=XMU_MIN
        borns(ip+dg%np)=XMU_MAX
        borni(2*dg%np+ip)=xmin
        borns(2*dg%np+ip)=xmax
        filter_inc(ip,2*dg%np+ip)=.true.
        borni(3*dg%np+ip)=xmin
        borns(3*dg%np+ip)=xmax
        filter_inc(ip,3*dg%np+ip)=.true.

        nest=0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if(dpa%estime(current_ic,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
              filter_inc(ip,4*dg%np+nestim)=.true.
            end if
            ifem=dg%repfem(jm)
            indam=dpa%iam(ic,ifem)
            filter_inc(ip,4*dg%np+dpa%nmumest(ic)+indam)=.true.
            filter_inc(ip,4*dg%np+dpa%nmumest(ic)+dpa%namest(ic)+indam)=.true.
           end if
        end do

      end do
      do jm=1,dpa%nmumest(ic)
         borni(4*dg%np+jm)=XMU_MIN
         borns(4*dg%np+jm)=XMU_MAX
      end do
      do jm=1,dpa%namest(ic)
          borni(4*dg%np+dpa%nmumest(ic)+jm)=xmin
          borns(4*dg%np+dpa%nmumest(ic)+jm)=xmax
          borni(4*dg%np+dpa%nmumest(ic)+dpa%namest(ic)+jm)=xmin
          borns(4*dg%np+dpa%nmumest(ic)+dpa%namest(ic)+jm)=xmax
      end do
!
! Point de depart
      nestim=0
      do ip=1,dg%np
        par(ip)=std(ip)
        par(ip+dg%np)=xmoyp(ip)
        par(2*dg%np+ip)=ap(ip)/2.d0
        par(3*dg%np+ip)=ap(ip)/2.d0
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        nest=0
        somxmu=0.d0
        do jm=nm1,nm2
         if (dpa%estime(ic,jm)) then
          nest=nest+1
          ifem=dg%repfem(jm)
          indam=dpa%iam(ic,ifem)
          if (nest.le.estmum(ip))then
            nestim=nestim+1
            par(4*dg%np+nestim)=xmoym(jm)
            somxmu=somxmu+xmoym(jm)
          else
            xmoym(jm)=-somxmu
          end if
          par(4*dg%np+dpa%nmumest(ic)+indam)= am(jm)/2.d0
          par(4*dg%np+dpa%nmumest(ic)+dpa%namest(ic)+indam)= am(jm)/2.d0
         end if
        end do
      end do
!
! Marche le long du chromosome
      ifail=0
      lrtsol%lrtmax=-1.d75

      ntotal=0
      do chr=1,dataset%map%nchr
        ntotal=ntotal+dataset%map%get_npo(chr)
        indexchr(chr)=ntotal
      end do

      allocate (maxpar(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),npar))
      allocate (fpsave1(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%np))
      allocate (fmsave1(dataset%map%nchr,dataset%map%nchr,dataset%map%get_maxnpo(),dataset%map%get_maxnpo(),dg%nm))

      !$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(par)  &
      !$OMP PRIVATE(ip,jm,chr,chr2,ilong2,n1,n,nstart,ig1,ig2,pddp,pddm,ngeno1,ngeno2,kd,kkd,kd2,ig,ifail,nestim,ii,f01_t)
      current_ic = ic
      call init_sub_2qtl(dataset,spt)
      carydf=>save_carydf
      somcddf => save_somcddf
      somydf => save_somydf
      effdf => save_effdf
      estmum => save_estmum
      somcd => save_somcd
      somy => save_somy
      cary => save_cary
      eff => save_eff
      fp_1 => save_fp_1
      fm_1 => save_fm_1
      fp0 => save_fp0
      fm0 => save_fm0
      f_1 = save_f_1
      f0 = save_f0

      !$OMP DO
      do nprime=1,ntotal
!       print *,'start:',nprime
       n=nprime
       do chr=dataset%map%nchr,1,-1
        if ( indexchr(chr) >= nprime) then
          current_chr = chr
        else
          n=n-dataset%map%get_npo(chr)
        end if
       end do
       chr=current_chr
       do chr2=chr,dataset%map%nchr
         if (chr2 == chr ) then
           nstart = n +1
         else
           nstart = 1
         end if

         current_chr2 = chr2
         do n1=nstart,dataset%map%get_npo(chr2)
         ! print *,chr,chr2,n,n1
          do ip=1,dg%np
           somppdf1(ip)=0.d0
           somppdf2(ip)=0.d0
           carppdf11(ip)=0.d0
           carppdf12(ip)=0.d0
           carppdf22(ip)=0.d0
           somppydf1(ip)=0.d0
           somppydf2(ip)=0.d0
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
             sompp1(ig)=0.d0
             sompm1(ig)=0.d0
             carpp1(ig)=0.d0
             carpm1(ig)=0.d0
             sompp1y(ig)=0.d0
             sompm1y(ig)=0.d0
             sompp1m1(ig)=0.d0
             do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
              kkd=spt%ndesc(chr,kd)
                if(dpa%presentc(ic,kkd)) then
                  pddp(kd,1)=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  pddm(kd,1)=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
!
               if (dpa%estime(ic,jm)) then
                 sompp1(ig)=sompp1(ig)+pddp(kd,1)
                 carpp1(ig)=carpp1(ig)+pddp(kd,1)*pddp(kd,1)
                 sompp1y(ig)=sompp1y(ig)+pddp(kd,1)*dpa%y(ic,kkd)
                 sompm1(ig)=sompm1(ig)+pddm(kd,1)
                 carpm1(ig)=carpm1(ig)+pddm(kd,1)*pddm(kd,1)
                 sompm1y(ig)=sompm1y(ig)+pddm(kd,1)*dpa%y(ic,kkd)
                 sompp1m1(ig)=sompp1m1(ig)+pddp(kd,1)*pddm(kd,1)
               else
                 somppdf1(ip)=somppdf1(ip)+pddp(kd,1)
                 carppdf11(ip)=carppdf11(ip)+pddp(kd,1)*pddp(kd,1)
                 somppydf1(ip)=somppydf1(ip)+pddp(kd,1)*dpa%y(ic,kkd)
               end if
              end if
             end do !kd
            end do !ig
            ! second chromosome
            do ig=spt%ngenom(chr2,jm)+1,spt%ngenom(chr2,jm+1)
               sompp2(ig)=0.d0
               sompm2(ig)=0.d0
               carpp2(ig)=0.d0
               carpm2(ig)=0.d0
               sompp2y(ig)=0.d0
               sompm2y(ig)=0.d0
               sompp2m2(ig)=0.d0

               do kd=spt%ngend(chr2,ig)+1,spt%ngend(chr2,ig+1)
                kkd=spt%ndesc(chr2,kd)
                if(dpa%presentc(ic,kkd)) then
                  pddp(kd,2)=-spt%pdd(chr2,kd,1,n1)-spt%pdd(chr2,kd,2,n1)+spt%pdd(chr2,kd,3,n1)+spt%pdd(chr2,kd,4,n1)
                  pddm(kd,2)=-spt%pdd(chr2,kd,1,n1)+spt%pdd(chr2,kd,2,n1)-spt%pdd(chr2,kd,3,n1)+spt%pdd(chr2,kd,4,n1)
                  if (dpa%estime(ic,jm)) then
                   sompp2(ig)=sompp2(ig)+pddp(kd,2)
                   carpp2(ig)=carpp2(ig)+pddp(kd,2)*pddp(kd,2)
                   sompp2y(ig)=sompp2y(ig)+pddp(kd,2)*dpa%y(ic,kkd)
                   sompm2(ig)=sompm2(ig)+pddm(kd,2)
                   carpm2(ig)=carpm2(ig)+pddm(kd,2)*pddm(kd,2)
                   sompm2y(ig)=sompm2y(ig)+pddm(kd,2)*dpa%y(ic,kkd)
                   sompp2m2(ig)=sompp2m2(ig)+pddp(kd,2)*pddm(kd,2)
                  else
                   somppdf2(ip)=somppdf2(ip)+pddp(kd,2)
                   carppdf22(ip)=carppdf22(ip)+pddp(kd,2)*pddp(kd,2)
                   somppydf2(ip)=somppydf2(ip)+pddp(kd,2)*dpa%y(ic,kkd)
                  end if
              end if
             end do ! kd
            end do ! ig

            !!*********************************************************************
            do ig1=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
             do ig2=spt%ngenom(chr2,jm)+1,spt%ngenom(chr2,jm+1)
              sompp1p2(ig1,ig2)=0.d0
              sompp1m2(ig1,ig2)=0.d0
              sompp2m1(ig1,ig2)=0.d0
              sompm1m2(ig1,ig2)=0.d0

              if ( (spt%ngend(chr,ig1+1) - spt%ngend(chr,ig1)+1) &
               /= (spt%ngend(chr2,ig2+1) - spt%ngend(chr2,ig2)+1) ) then
                call log_mess("Number of desc for "//trim(str(ig1))//" chrom ("//trim(dataset%map%chromo(chr))&
                //" :"//trim(str(spt%ngend(chr,ig1+1) - spt%ngend(chr,ig1)+1)),ERROR_DEF)
                call log_mess("Number of desc for "//trim(str(ig2))//" chrom ("//trim(dataset%map%chromo(chr2))&
                //" :"//trim(str(spt%ngend(chr2,ig2+1) - spt%ngend(chr2,ig2)+1)),ERROR_DEF)
                call stop_application("** Error dev **")
              end if

              kd2=spt%ngend(chr2,ig2)+1-1
              do kd=spt%ngend(chr,ig1)+1,spt%ngend(chr,ig1+1)
                kd2 = kd2 + 1
                kkd=spt%ndesc(chr,kd)
                if(dpa%presentc(ic,kkd)) then
                  pddp(kd,1)=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  pddp(kd,2)=-spt%pdd(chr2,kd2,1,n1)-spt%pdd(chr2,kd2,2,n1)+spt%pdd(chr2,kd2,3,n1)+spt%pdd(chr2,kd2,4,n1)
                  pddm(kd,1)=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  pddm(kd,2)=-spt%pdd(chr2,kd2,1,n1)+spt%pdd(chr2,kd2,2,n1)-spt%pdd(chr2,kd2,3,n1)+spt%pdd(chr2,kd2,4,n1)
!
               if (dpa%estime(ic,jm)) then
                 sompp1p2(ig1,ig2)=sompp1p2(ig1,ig2)+pddp(kd,1)*pddp(kd,2)
                 sompp1m2(ig1,ig2)=sompp1m2(ig1,ig2)+pddp(kd,1)*pddm(kd,2)
                 sompp2m1(ig1,ig2)=sompp2m1(ig1,ig2)+pddp(kd,2)*pddm(kd,1)
                 sompm1m2(ig1,ig2)=sompm1m2(ig1,ig2)+pddm(kd,1)*pddm(kd,2)
               else
                 carppdf12(ip)=carppdf12(ip)+pddp(kd,1)*pddp(kd,2)
               end if
              end if
             end do ! kd
            end do ! ig2
           end do !ig1
!
           end do
          end do
!
! Optimisation de la vraisemblance a la position dx
          ifail=1
          !call e04jyf(npar,ibound,funct_1car_2qtl,borni,borns,par,f01_t,iw,liw,w,lw,iuser,user,ifail)
    !      call minimizing_funct(npar,ibound,funct_1car_2qtl,borni,borns,par,f01_t,iuser,user,ifail)

       call minimizing_funct_family_sire(dataset,npar,ibound,funct_1car_2qtl_family,filter_inc,&
          fp01,borni,borns,par,f01_t,iuser,user,ifail)

!          if (ifail.ne.0) call log_mess('Error code e04jyf H1 : '//char(str(ifail))//&
!          ', opti_1car_2qtl, LRT('//char(str(ix))//','//char(str(ix1)//')'),WARNING_DEF)
       fpsave1(chr,chr2,n,n1,:) = fp01
       fmsave1(chr,chr2,n,n1,:) = fm01
       maxpar(chr,chr2,n,n1,:) = par(:npar)

!
! Test 1/2 QTL
          lrt1_2(chr,chr2,n,n1)=-2.d0*(f01_t-f_1)

          do ii=1,dg%np
            xlrp2(chr,chr2,n,n1,ii)=-2.d0*(fp01(ii)-fp_1(ii))
            xlrp2_0(chr,chr2,n,n1,ii)=-2.d0*(fp01(ii)-fp0(ii))
            lrtsol%pater_eff2(chr,chr2,ii,n,n1,1)=par(2*dg%np+ip)
            lrtsol%pater_eff2(chr,chr2,ii,n,n1,2)=par(3*dg%np+ip)
          end do
          do ii=1,dg%nm
            xlrm2(chr,chr2,n,n1,ii)=-2.d0*(fm01(ii)-fm_1(ii))
            xlrm2_0(chr,chr2,n,n1,ii)=-2.d0*(fm01(ii)-fm0(ii))
          end do

          do ii=1,dpa%namest(ic)
           lrtsol%mater_eff2(chr,chr2,ii,n,n1,1)=par(4*dg%np+dpa%nmumest(ic)+ii)
           lrtsol%mater_eff2(chr,chr2,ii,n,n1,2)=par(4*dg%np+dpa%nmumest(ic)+dpa%namest(ic)+ii)
          end do

! Test 0/2 QTL
          lrt0_2(chr,chr,n,n1)=-2.d0*(f01_t-f0)

         end do ! n1
       end do !chr2
      end do !nprime
      !$OMP END DO
      call end_sub_2qtl
      !$OMP END PARALLEL
      do chr=1,dataset%map%nchr
        do n=1,dataset%map%get_npo(chr)
          do chr2=chr,dataset%map%nchr
            if (chr2 == chr ) then
              nstart = n+1
            else
              nstart=1
            end if
            do n1=nstart,dataset%map%get_npo(chr2)
                 if(lrtsol%lrtmax(1) < lrt1_2(chr,chr2,n,n1)) then
                      lrtsol%chrmax(0)=chr
                      lrtsol%chrmax(1)=chr2
                      lrtsol%nxmax(0)=n
                      lrtsol%nxmax(1)=n1
                      lrtsol%lrtmax(1)=lrt1_2(chr,chr2,n,n1)
                      lrtsol%lrtmax(0)=lrt0_2(chr,chr2,n,n1)
                      nestim=0
                      par(:npar) = maxpar(chr,chr2,n,n1,:npar)
                      do ip=1,dg%np
                         std2(ip)=par(ip)
                         xmoyp2(ip)=par(dg%np+ip)
                         do iq=1,2
                          ap2(ip,iq)=par(2*dg%np+(iq-1)*dg%np+ip)
                         end do
                         nest=0
                         somxmu=0.d0
                         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                          if (dpa%estime(ic,jm)) then
                           nest=nest+1
                           ifem=dg%repfem(jm)
                           indam=dpa%iam(ic,ifem)
                           do iq=1,2
                             am2(jm,iq)=par(4*dg%np+dpa%nmumest(ic)+(iq-1)*dpa%namest(ic)+indam)
                           end do
                           if (nest.le.estmum(ip)) then
                            nestim=nestim+1
                            xmoym2(jm)=par(4*dg%np+nestim)
                            somxmu=somxmu+xmoym2(jm)
                           else
                            xmoym2(jm)=-somxmu
                           end if
                          else
                           do iq=1,2
                            am2(jm,iq)=0.d0
                           end do
                           xmoym2(jm)=0.d0
                          end if
                       end do ! np
                    end do
               end if
            end do
          end do
        end do
      end do

      ! 04/2013 - New structure for LRT
      call lrtsol%LRT%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),lrt1_2,2)
      call lrtsol%LRT%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),lrt0_2,1)
      do ii=1,dg%np
       call lrtsol%LRT_SIRES(ii)%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrp2(:,:,:,:,ii),2)
       call lrtsol%LRT_SIRES(ii)%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrp2_0(:,:,:,:,ii),1)
      end do
      do ii=1,dg%nm
       call lrtsol%LRT_DAMS(ii)%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrm2(:,:,:,:,ii),2)
       call lrtsol%LRT_DAMS(ii)%add2(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrm2_0(:,:,:,:,ii),1)
      end do

      deallocate (lrt0_2,lrt1_2,xlrp2,xlrm2,xlrp2_0,xlrm2_0)

      deallocate (filter_inc)
      deallocate (maxpar)
      deallocate (fpsave1)
      deallocate (fmsave1)

      end subroutine opti_1car_2qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_1car_2qtl
!! NAME
!!    funct_1car_2qtl
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous H1
!! HISTORY
!!
!! SOURCE
      subroutine funct_1car_2qtl(n,x,f,iuser,user)
      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      integer :: ip,nm1,nm2,jm,ngeno1,ngeno2,ig
      integer :: nestim,nest,ifem,indam,ngeno21,ngeno22,ig2

      real (kind=dp) :: sig,var,vmere,vpf,xmup,xmup2,a1p,a2p,a1p2
      real (kind=dp) :: a2p2,vdf,somxmu,xmum,xmum2,xmupm,a1m,a2m
      real (kind=dp) :: a1m2,a2m2,z
!****************************************************************************
      f=0.d0
      fm01=0.d0
      fp01=0.d0

      nestim=0
      do ip=1,p_dg%np
        sig=x(ip)
        var=sig*sig
        xmup=x(ip+p_dg%np)
        xmup2=xmup*xmup
        a1p=x(2*p_dg%np+ip)
        a2p=x(3*p_dg%np+ip)
        a1p2=a1p*a1p
        a2p2=a2p*a2p
        vdf=carydf(ip)+dble(effdf(ip))*xmup2-2.d0*xmup*somydf(ip)
        vdf=vdf+a1p2*carppdf11(ip)+a2p2*carppdf22(ip)        &
            -2.d0*a1p*somppydf1(ip)-2.d0*a2p*somppydf2(ip)   &
            +2.d0*xmup*(a1p*somppdf1(ip)+a2p*somppdf2(ip))   &
            +2.d0*a1p*a2p*carppdf12(ip)
        vdf=0.5d0*vdf/var
        fp01(ip)=vdf+dble(effdf(ip))*dlog(sig)
        f=f+fp01(ip)

        nm1=p_dg%nmp(ip)+1
        nm2=p_dg%nmp(ip+1)
        somxmu=0.d0
        nest=0
        do jm=nm1,nm2
            if(p_dpa%estime(current_ic,jm)) then
               nest=nest+1
               if(nest.le.estmum(ip)) then
                  nestim=nestim+1
                  xmum=x(4*p_dg%np+nestim)
                  somxmu=somxmu+xmum
               else
                  xmum=-somxmu
               end if
               xmum2=xmum*xmum
               xmupm=2.d0*(xmup+xmum)
               ifem=p_dg%repfem(jm)
               indam=p_dpa%iam(current_ic,ifem)
               a1m=x(4*p_dg%np+p_dpa%nmumest(current_ic)+indam)
               a2m=x(4*p_dg%np+p_dpa%nmumest(current_ic)+p_dpa%namest(current_ic)+indam)
               a1m2=a1m*a1m
               a2m2=a2m*a2m
               z=cary(jm)+dble(eff(jm))*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
               vmere=0.d0
               ngeno1=p_spt%ngenom(current_chr,jm)+1
               ngeno2=p_spt%ngenom(current_chr,jm+1)
               ngeno21=p_spt%ngenom(current_chr2,jm)+1
               ngeno22=p_spt%ngenom(current_chr2,jm+1)

            if( (ngeno2-ngeno1 == 0) .and. (ngeno22-ngeno21 == 0) )then
               ig=ngeno1
               ig2=ngeno21
               vpf=z +a1p2*carpp1(ig)+a2p2*carpp2(ig)+a1m2*carpm1(ig    &
                      ) +a2m2*carpm2(ig) +2.d0* (a1p*a2p*sompp1p2(ig,ig2)+  &
                      a1p *a1m*sompp1m1(ig)+a1p*a2m*sompp1m2(ig,ig2) +a2p   &
                      *a1m *sompp2m1(ig,ig2)+a2p*a2m*sompp2m2(ig)+a1m*a2m   &
                      *sompm1m2(ig,ig2)) +xmupm* (a1p*sompp1(ig)+a2p        &
                      *sompp2(ig)+a1m*sompm1(ig)+a2m*sompm2(ig)) -2.d0  &
                      *(a1p*sompp1y(ig)+ a2p*sompp2y(ig)+ a1m           &
                      *sompm1y(ig)+a2m*sompm2y(ig))
               !print *,"vpf:",vpf
               vpf=0.5d0*vpf/var

               fm01(jm)=vpf+dble(eff(jm))*dlog(sig)
               !print *,"fm01:",fm01(jm)
           else
               do ig=ngeno1,ngeno2
                do ig2=ngeno21,ngeno22
                  vpf=z +a1p2*carpp1(ig)+a2p2*carpp2(ig)+a1m2*carpm1(ig   &
                      ) +a2m2*carpm2(ig) +2.d0* (a1p*a2p*sompp1p2(ig,ig2)+   &
                      a1p *a1m*sompp1m1(ig)+a1p*a2m*sompp1m2(ig,ig2) +a2p    &
                      *a1m *sompp2m1(ig,ig2)+a2p*a2m*sompp2m2(ig)+a1m*a2m    &
                      *sompm1m2(ig,ig2)) +xmupm* (a1p*sompp1(ig)+a2p         &
                      *sompp2(ig)+a1m*sompm1(ig)+a2m*sompm2(ig)) -2.d0   &
                      *(a1p*sompp1y(ig)+ a2p*sompp2y(ig)+ a1m            &
                      *sompm1y(ig)+a2m*sompm2y(ig))
                  vpf=-0.5d0*vpf/var
                  ! note ofi : a voir avec pascale...il y a 2 proba de genotype sur 2 chromosomes....
                  if ( current_chr == current_chr2 ) then
                    vmere=vmere+p_spt%probg(current_chr,ig)*dexp(vpf)
                  else
                    vmere=vmere+p_spt%probg(current_chr,ig)*p_spt%probg(current_chr2,ig2)*dexp(vpf)
                  end if
               end do
              end do
               fm01(jm)=-dlog(vmere)+(dble(eff(jm))*dlog(sig))
           end if


               f=f+fm01(jm)
               fp01(ip)=fp01(ip)+fm01(jm)
            else
               fm01(jm)=0.d0
            end if
       end do
      end do

      end subroutine funct_1car_2qtl
!!***

!!****f* m_qtlmap_analyse_unitrait/funct_1car_2qtl_family
!! NAME
!!    funct_1car_2qtl_family
!! DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous H1
!! HISTORY
!!
!! SOURCE
   subroutine funct_1car_2qtl_family(ip,n,x,f,iuser,user)
      integer         , intent(in)                  :: ip,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      integer :: jm,ngeno1,ngeno2,ig
      integer :: nestim,nest,ifem,indam,ngeno21,ngeno22,ig2

      real (kind=dp) :: sig,var,vmere,vpf,xmup,xmup2,a1p,a2p,a1p2
      real (kind=dp) :: a2p2,vdf,somxmu,xmum,xmum2,xmupm,a1m,a2m
      real (kind=dp) :: a1m2,a2m2,z
!****************************************************************************
      f=0.d0
      fm01=0.d0
      sig=x(ip)
      var=sig*sig
      xmup=x(ip+p_dg%np)
      xmup2=xmup*xmup
      a1p=x(2*p_dg%np+ip)
      a2p=x(3*p_dg%np+ip)
      a1p2=a1p*a1p
      a2p2=a2p*a2p
      vdf=carydf(ip)+dble(effdf(ip))*xmup2-2.d0*xmup*somydf(ip)
      vdf=vdf+a1p2*carppdf11(ip)+a2p2*carppdf22(ip)        &
            -2.d0*a1p*somppydf1(ip)-2.d0*a2p*somppydf2(ip)   &
            +2.d0*xmup*(a1p*somppdf1(ip)+a2p*somppdf2(ip))   &
            +2.d0*a1p*a2p*carppdf12(ip)
      vdf=0.5d0*vdf/var
      f=vdf+dble(effdf(ip))*dlog(sig)
      somxmu=0.d0
      nest=0
      do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
            if(p_dpa%estime(current_ic,jm)) then
               nest=nest+1
               if(nest.le.estmum(ip)) then
                  if (ip>1) then
                     nestim=sum(estmum(:ip-1))+nest
                  else
                     nestim=nest
                  end if
                  xmum=x(4*p_dg%np+nestim)
                  somxmu=somxmu+xmum
               else
                  xmum=-somxmu
               end if
               xmum2=xmum*xmum
               xmupm=2.d0*(xmup+xmum)
               ifem=p_dg%repfem(jm)
               indam=p_dpa%iam(current_ic,ifem)
               a1m=x(4*p_dg%np+p_dpa%nmumest(current_ic)+indam)
               a2m=x(4*p_dg%np+p_dpa%nmumest(current_ic)+p_dpa%namest(current_ic)+indam)
               a1m2=a1m*a1m
               a2m2=a2m*a2m
               z=cary(jm)+dble(eff(jm))*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
               vmere=0.d0
               ngeno1=p_spt%ngenom(current_chr,jm)+1
               ngeno2=p_spt%ngenom(current_chr,jm+1)
               ngeno21=p_spt%ngenom(current_chr2,jm)+1
               ngeno22=p_spt%ngenom(current_chr2,jm+1)

            if( (ngeno2-ngeno1 == 0) .and. (ngeno22-ngeno21 == 0) )then
               ig=ngeno1
               ig2=ngeno21
               vpf=z +a1p2*carpp1(ig)+a2p2*carpp2(ig)+a1m2*carpm1(ig    &
                      ) +a2m2*carpm2(ig) +2.d0* (a1p*a2p*sompp1p2(ig,ig2)+  &
                      a1p *a1m*sompp1m1(ig)+a1p*a2m*sompp1m2(ig,ig2) +a2p   &
                      *a1m *sompp2m1(ig,ig2)+a2p*a2m*sompp2m2(ig)+a1m*a2m   &
                      *sompm1m2(ig,ig2)) +xmupm* (a1p*sompp1(ig)+a2p        &
                      *sompp2(ig)+a1m*sompm1(ig)+a2m*sompm2(ig)) -2.d0  &
                      *(a1p*sompp1y(ig)+ a2p*sompp2y(ig)+ a1m           &
                      *sompm1y(ig)+a2m*sompm2y(ig))
               !print *,"vpf:",vpf
               vpf=0.5d0*vpf/var

               fm01(jm)=vpf+dble(eff(jm))*dlog(sig)
               !print *,"fm01:",fm01(jm)
           else
               do ig=ngeno1,ngeno2
                do ig2=ngeno21,ngeno22
                  vpf=z +a1p2*carpp1(ig)+a2p2*carpp2(ig)+a1m2*carpm1(ig   &
                      ) +a2m2*carpm2(ig) +2.d0* (a1p*a2p*sompp1p2(ig,ig2)+   &
                      a1p *a1m*sompp1m1(ig)+a1p*a2m*sompp1m2(ig,ig2) +a2p    &
                      *a1m *sompp2m1(ig,ig2)+a2p*a2m*sompp2m2(ig)+a1m*a2m    &
                      *sompm1m2(ig,ig2)) +xmupm* (a1p*sompp1(ig)+a2p         &
                      *sompp2(ig)+a1m*sompm1(ig)+a2m*sompm2(ig)) -2.d0   &
                      *(a1p*sompp1y(ig)+ a2p*sompp2y(ig)+ a1m            &
                      *sompm1y(ig)+a2m*sompm2y(ig))
                  vpf=-0.5d0*vpf/var
                  ! note ofi : a voir avec pascale...il y a 2 proba de genotype sur 2 chromosomes....
                  if ( current_chr == current_chr2 ) then
                    vmere=vmere+p_spt%probg(current_chr,ig)*dexp(vpf)
                  else
                    vmere=vmere+p_spt%probg(current_chr,ig)*p_spt%probg(current_chr2,ig2)*dexp(vpf)
                  end if
               end do
              end do
               fm01(jm)=-dlog(vmere)+(dble(eff(jm))*dlog(sig))
           end if
               f=f+fm01(jm)
            else
               fm01(jm)=0.d0
            end if
       end do

     end subroutine funct_1car_2qtl_family
!!***

!!****f* m_qtlmap_analyse_unitrait/set_solution_hypothesis0
!! NAME
!!    set_solution_hypothesis0
!! DESCRIPTION
!!
!! HISTORY
!!
!! SOURCE
       subroutine set_solution_hypothesis0(dataset,ic,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       integer                            ,intent(in)       :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=0
       !  Mean family
       nteff = 1
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 2

       maxNbPar = max(dg%np,count(dpa%estime(ic,:)))
       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.
       do ip=1,dg%np
            incsol%sig(1,ip) = sig1(ip)*dpm%sigt(ic)
       end do

       ieff=1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmu1p(ip)*dpm%sigt(ic) + dpm%xmut(ic)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0

          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = xmu1m(ifem)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end do

         end if

       end subroutine set_solution_hypothesis0
!!***

!!****f* m_qtlmap_analyse_unitrait/set_solution_hypothesis1
!! NAME
!!    set_solution_hypothesis1
!! DESCRIPTION
!!
!! HISTORY
!!
!! SOURCE
      subroutine set_solution_hypothesis1(dataset,spt,ic,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(PDD_BUILD)            ,intent(in) :: spt
       integer                        , intent(in)          :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=1
       !  Mean family , Qtl effect
       nteff = 2
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 4

       maxNbPar = max(dg%np,count(dpa%estime(ic,:)))
       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,1))

       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = std(ip)*dpm%sigt(ic)
       end do

       ieff=1
         incsol%qtl_groupeName(1,1)=ieff
         incsol%groupeName(ieff) = 'Sire Qtl effect'
         incsol%nbParameterGroup(ieff)=dg%np

         do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = ap(ip)*dpm%sigt(ic)
         end do

         if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam Qtl effect'
           incsol%nbParameterGroup(ieff)=dpa%namest(ic)
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%paramaterValue(ieff,ifem) = am(jm)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end if


       ieff = ieff +1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmoyp(ip)*dpm%sigt(ic) + dpm%xmut(ic)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
         do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm) ) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = xmoym(jm)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end do
         end if

       end subroutine set_solution_hypothesis1
!!***

!!****f* m_qtlmap_analyse_unitrait/set_solution_hypothesis2
!! NAME
!!    set_solution_hypothesis2
!! DESCRIPTION
!!
!! HISTORY
!!
!! SOURCE
      subroutine set_solution_hypothesis2(dataset,spt,ic,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(PDD_BUILD)            ,intent(in) :: spt
       integer                        , intent(in)          :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff,iq
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=2
       !  Mean family , Qtl effect 1, QTL effect2
       nteff = 3
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 6

       maxNbPar = max(dg%np,count(dpa%estime(ic,:)))
       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,2))

       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%nbParameterGroup=0.d0
       incsol%eqtl_print=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = std2(ip)*dpm%sigt(ic)
       end do

       ieff=0
       do iq=1,2
         ieff = ieff +1
         incsol%qtl_groupeName(1,iq)=ieff
         incsol%groupeName(ieff) = 'Sire Qtl effect '//trim(str(iq))
         incsol%nbParameterGroup(ieff)=dg%np
         do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = ap2(ip,iq)*dpm%sigt(ic)
         end do


         if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam Qtl effect '//trim(str(iq))
           incsol%nbParameterGroup(ieff)=dpa%namest(ic)
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%paramaterValue(ieff,ifem) = am2(jm,iq)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end if
       end do

       ieff = ieff +1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np


       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmoyp2(ip)*dpm%sigt(ic) + dpm%xmut(ic)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = xmoym2(jm)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
          end do

         end if

       end subroutine set_solution_hypothesis2
!!***

     end module m_qtlmap_analyse_unitrait
