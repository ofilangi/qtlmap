!!****m* INPUT/m_qtlmap_parameter_file
!!  NAME
!!    m_qtlmap_parameter_file -- This modules manages the parameter user file (p_analyse).
!!  SYNOPSIS

!!  DESCRIPTION
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
module m_qtlmap_type_parameter
      use m_qtlmap_constant
      use m_qtlmap_base
      use m_qtlmap_log
      use m_qtlmap_type_cli

      implicit none
      save
!!//   map file
      character(len=LEN_L)    ,parameter    , public :: K_MAP                ='in_map'
!!//   genealogy file
      character(len=LEN_L)    ,parameter    , public :: K_GENEA              ='in_genealogy'
!!//   breed file of fonder
      character(len=LEN_L)    ,parameter    , public :: K_RACE               ='in_pop'
!!//   phenotypics user file
      character(len=LEN_L)    ,parameter    , public :: K_TRAITS             ='in_traits'
!!//   genotypics user file
      character(len=LEN_L)    ,parameter    , public :: K_GENOTYPE           ='in_genotype'
!!//   model phenotypic user file
      character(len=LEN_L)    ,parameter    , public :: K_MODEL              ='in_model'
!!//   parameter simulation user file
      character(len=LEN_L)    ,parameter    , public :: K_PARAMSIM           ='in_paramsimul'
!!//   output result file
      character(len=LEN_L)    ,parameter    , public :: K_OUTSIM             ='out_maxlrt'
!!//   output probabilities of transmission file
      character(len=LEN_L)    ,parameter    , public :: K_PDED               ='out_pded'
!!//   output probabilities of join transmission file
      character(len=LEN_L)    ,parameter    , public :: K_PDECPLE            ='out_pdedjoin'
!!//   output paternal effect file
      character(len=LEN_L)    ,parameter    , public :: K_PATEFF             ='out_pateff'
!!//   output maternal effect file
      character(len=LEN_L)    ,parameter    , public :: K_MATEFF             ='out_mateff'
!!//   LRT curve output file (General LRT and by half-sib family)
      character(len=LEN_L)    ,parameter    , public :: K_LRTSIRE            ='out_lrtsires'
!!//    LRT curve output file
      character(len=LEN_L)    ,parameter    , public :: K_LRTDAM             ='out_lrtdams'
!!//   Allele frequency information
      character(len=LEN_L)    ,parameter    , public :: K_FREQALL            ='out_freqall'
!!//   output coefficent discriminante analysis result file
      character(len=LEN_L)    ,parameter    , public :: K_COEFFDA            ='out_coeffda'
!!//   phases probabilities of reproductor file
      character(len=LEN_L)    ,parameter    , public :: K_PHASES             ='out_phases'
!!//   The most likely phase haplotypes for offspring
      character(len=LEN_L)    ,parameter    , public :: K_PHASES_OFFSPRING    ='out_phases_offspring'
!!    haplotypes files
      character(len=LEN_L)    ,parameter    , public :: K_HAPLOTYPES         ='out_haplotypes'
!!    informativity file
      character(len=LEN_L)    ,parameter    , public :: K_INFORMATIVITY       ='out_informativity'
!!   main output result file
      character(len=LEN_L)    ,parameter    , public :: K_OUTPUT             ='out_output'
!!   summary result file
      character(len=LEN_L)    ,parameter    , public :: K_SUMM               ='out_summary'
!!   chromosom list to explore
      character(len=LEN_L)    ,parameter    , public :: K_CHROM              ='opt_chromosome'
!!   Minimal number of progeny by dam to considere full sib family in the analysis
      character(len=LEN_L)    ,parameter    , public :: K_NDMIN              ='opt_ndmin'
!!   step definition for LA analysis (in Morgan)
      character(len=LEN_L)    ,parameter    , public :: K_STEP               ='opt_step'
!!   Minimal dam phase probability
      character(len=LEN_L)    ,parameter    , public :: K_MINDAMPHASEPROB    ='opt_mindamphaseproba'
!!   Minimal sire phase probability
      character(len=LEN_L)    ,parameter    , public :: K_MINSIREPHASEPROB   ='opt_minsirephaseproba'
!!   Unknown genotype value
      character(len=LEN_L)    ,parameter    , public :: K_UNKNOWN_GENO       ='opt_unknown_char'
!!   coefficent used in the cholesky decomposition to establish the parameters estimabilities
      character(len=LEN_L)    ,parameter    , public :: K_CHOLESKY           ='opt_eps_cholesky'
!!   Threshold to test confusion betwwen level inside a contingence matrix
      character(len=LEN_L)    ,parameter    , public :: K_THRES_CONFUSION    ='opt_eps_confusion'
!!   Threshold to check the equilibrium of marker transmission within each family
      character(len=LEN_L)    ,parameter    , public :: K_PSEUILHWE          ='opt_eps_hwe'
!!   Threshold for convergence in the linear mode heteroscedastic
      character(len=LEN_L)    ,parameter    , public :: K_LINEAR_CONV        ='opt_eps_linear_heteroscedastic'
!!   number maximum of iteration (heteroscedastic linear analysis)
      character(len=LEN_L)    ,parameter    , public :: K_MAX_LINEAR_ITERATION = 'opt_max_iteration_linear_heteroscedastic'
!!   Minimal paternal phase probability
      character(len=LEN_L)    ,parameter    , public :: K_PROB_SEUIL_RECOMB  ='opt_eps_recomb'
!!   Length of the haplotypes to be followed in LD analyses
      character(len=LEN_L)    ,parameter    , public :: K_LONGHAP            ='opt_longhap'
!!   Number maximum of haplotype allowed
      character(len=LEN_L)    ,parameter    , public :: K_NB_HAPLO_PRIOR     ='opt_nb_haplo_prior'
!!   miminum probability of a haplotype in the LD analysis
      character(len=LEN_L)    ,parameter    , public :: K_PROB_HAPLO_MIN     ='opt_prob_haplo_min'
!!
      character(len=LEN_L)    ,parameter    , public :: K_LONG_MIN_IBS       ='opt_long_min_ibs'
!!  number mmaximum of objective function evaluation
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_MAXEVAL      ='opt_optim_maxeval'
!!   maximum elapsed time in a optimisation of a objective function
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_MAXTIME      ='opt_optim_maxtime'
!!   tolerance precison of the paramter X (m_qtlmap_optimization)
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_TOLX         ='opt_optim_tolx'
!!   tolerance precison of the objective function (m_qtlmap_optimization)
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_TOLF         ='opt_optim_tolf'
!!   tolerance precison of gradient (m_qtlmap_optimization)
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_TOLG         ='opt_optim_tolg'
!!   precision of numerical computation of the derivate F'(x) = (F(x-h)+F(x+h))/2*h (m_qtlmap_optimization)
      character(len=LEN_L)    ,parameter    , public :: K_OPTIM_H_PREC       ='opt_optim_h_precision'
!!   First marker to print the phased haplotypes
      character(len=LEN_L)    ,parameter    , public :: K_PHASES_OFFSPRING_MARK_START       ='opt_phases_offspring_marker_start'
!!   Last marker to print the phased haplotypes
      character(len=LEN_L)    ,parameter    , public :: K_PHASES_OFFSPRING_MARK_END       ='opt_phases_offspring_marker_end'
!!   id of gpu device to use with cuda version
      character(len=LEN_L)    ,parameter    , public :: K_GPU_DEVICE_ID       ='opt_gpu_device_id'
!!  DESCRIPTION
!!   For experimental methods, this module defines a generic key (opt_key_dev1,opt_key_dev2,...)
!!   The interpretation of theses keys is in charge of the experimental analysis
      character(len=LEN_L)    ,parameter    , public :: K_KEY_DEV_GEN       = 'opt_key_dev'
      !!   optional keys
      integer    , private,parameter              :: NUMBER_OPT_KEYS=23

      ! NUMBER OF AUTOMATIC KEYS DEFINED.......
!!   number of essential keys to defined by the user in the paramter file
      integer    , private,parameter              :: NUMBER_AUTO_KEYS=24
      ! ALL KEYS
      ! --------
!!   Number of qtlmap keys
      integer    , private,parameter              :: NUMBER_ALL_KEYS=NUMBER_AUTO_KEYS+NUMBER_OPT_KEYS

      !!   Index keys with a default value
      character(len=LEN_L)  , dimension(NUMBER_AUTO_KEYS) , private :: index_key     = (/ &
                                                                         K_OUTPUT,  &
                                                                         K_SUMM,    &
                                                                         K_CHROM,   &
                                                                         K_NDMIN,   &
                                                                         K_STEP,    &
                                                                         K_UNKNOWN_GENO,    &
                                                                         K_MINDAMPHASEPROB,    &
                                                                         K_MINSIREPHASEPROB,  &
                                                                         K_CHOLESKY,          &
                                                                         K_THRES_CONFUSION,   &
                                                                         K_PSEUILHWE,         &
                                                                         K_LINEAR_CONV,       &
                                                                         K_MAX_LINEAR_ITERATION,&
                                                                         K_PROB_SEUIL_RECOMB, &
                                                                         K_NB_HAPLO_PRIOR,&
                                                                         K_PROB_HAPLO_MIN,&
                                                                         K_LONG_MIN_IBS,&
                                                                         K_LONGHAP,&
                                                                         K_OPTIM_MAXEVAL,&
                                                                         K_OPTIM_MAXTIME,&
                                                                         K_OPTIM_TOLX,&
                                                                         K_OPTIM_TOLF,&
                                                                         K_OPTIM_TOLG,&
                                                                         K_OPTIM_H_PREC/)

!!   Default values for key predefined by qtlmap.
      character(len=LEN_L)  , dimension(NUMBER_AUTO_KEYS) , private :: default_values = (/ &
                                                                         "       ",  &
                                                                         "       ",  &
                                                                         "       ",  &
                                                                         "10000  ",  &
                                                                         "0.05   ",  &
                                                                         "0      ",  &
                                                                         "0.10   ",  &
                                                                         "0.90   ",  &
                                                                         "0.01   ",  &
                                                                         "0.70   ",  &
                                                                         "0.01   ",  &
                                                                         "0.5    ",  &
                                                                         "30     ",  &
                                                                         "0.05   ",  &
                                                                         "200    ",  &
                                                                         "0.01000",  &
                                                                         "4      ",  &
                                                                         "4      ",  &
                                                                         "1000000",  &
                                                                         "1000000",  &
                                                                         "0.00005",  &
                                                                         "0.00005",  &
                                                                         "0.00005",  &
                                                                         "0.00005"/)
!!   The entirely list of possibility key used by qtlmap
      character(len=LEN_L)  , dimension(NUMBER_ALL_KEYS) , private :: all_key = (/ &
                                                                         K_OUTPUT,  &
                                                                         K_LRTSIRE, &
                                                                         K_LRTDAM,  &
                                                                         K_PDED,    &
                                                                         K_PDECPLE, &
                                                                         K_PATEFF,  &
                                                                         K_MATEFF,  &
                                                                         K_FREQALL, &
                                                                         K_COEFFDA, &
                                                                         K_SUMM,    &
                                                                         K_CHROM,   &
                                                                         K_NDMIN,   &
                                                                         K_STEP,    &
                                                                         K_UNKNOWN_GENO,    &
                                                                         K_MINDAMPHASEPROB,    &
                                                                         K_MINSIREPHASEPROB,  &
                                                                         K_CHOLESKY, &
                                                                         K_THRES_CONFUSION,   &
                                                                         K_PSEUILHWE,         &
                                                                         K_LINEAR_CONV,&
                                                                         K_MAP,     &
                                                                         K_GENEA,   &
                                                                         K_RACE,   &
                                                                         K_TRAITS,  &
                                                                         K_GENOTYPE, &
                                                                         K_MODEL,    &
                                                                         K_PARAMSIM, &
                                                                         K_OUTSIM,  &
                                                                         K_MAX_LINEAR_ITERATION, &
                                                                         K_PROB_SEUIL_RECOMB,    &
                                                                         K_PHASES, &
                                                                         K_PHASES_OFFSPRING, &
                                                                         K_HAPLOTYPES,&
                                                                         K_NB_HAPLO_PRIOR,&
                                                                         K_PROB_HAPLO_MIN,&
                                                                         K_LONG_MIN_IBS,&
                                                                         K_LONGHAP,&
                                                                         K_OPTIM_MAXEVAL,&
                                                                         K_OPTIM_MAXTIME,&
                                                                         K_OPTIM_TOLX,&
                                                                         K_OPTIM_TOLF,&
                                                                         K_OPTIM_TOLG,&
                                                                         K_OPTIM_H_PREC,&
                                                                         K_PHASES_OFFSPRING_MARK_START,&
                                                                         K_PHASES_OFFSPRING_MARK_END,&
                                                                         K_GPU_DEVICE_ID,&
                                                                         K_INFORMATIVITY/)


  type ,public :: PARAMETER_BASE

      !// number maximum of evaluation (module optimization)
      integer                                         , public      :: optim_maxeval        = 100000
      !// maxtime for each optimization  (module optimization)
      real(kind=dp)                                   , public      :: optim_maxtime        = 100000
      !//  tolerance for the variable x  (module optimization)
      real(kind=dp)                                   , public      :: optim_tolx           = 5.d-5
      !//  tolerance for the function evaluation result f  (module optimization)
      real(kind=dp)                                   , public      :: optim_tolf           = 5.d-5
      !// tolerance for the gradient  (module optimization)
      real(kind=dp)                                   , public      :: optim_tolg           = 5.d-5
      !// Precision gradiant df/dx = f(x+h)-f(x-h)/2*h  (module optimization)
      real(kind=dp)                                   , public      :: optim_H_PRECISION    = 5.d-5
      !// Estimabilité / nb de descendant phenotype : (1er dim indicie caractere/ 2 eme dim : indice mere)
      integer                              ,            public      :: NDMIN  = -1
      !// option de traitement des familles en demi germains ou melange plein/ demi germains
      integer                              ,            public      :: OPT_SIB = -1
      !// Threshold to check the equilibrium of marker transmission within each family
      real (kind=dp)                       ,            public      :: PSEUILHWE
      !// threshold of the test to identify an abnormal recombination rate
      real (kind=dp)                       ,            public      :: PROB_SEUIL_RECOMB
      !// The analysis is interrupted if for a sire, none of its phases reach this threshold
      real (kind=dp)                       ,            public      :: PHPSEUIL
      !// Threshold above which the probable maternal phases will be considered in the analysis
      real (kind=dp)                       ,            public      :: PRSEUIL
      !// Cholesky threshold to build X'X matrix
      real (kind=dp)                       ,            public      :: SEUIL_CHO
      !// Minimal probability to consider a predicted phase using flanking phased markers
      real   (kind=dp)                     ,            public      :: PROB_PHASE_DESC=95.d-2
      !// Threshold to test confusion betwwen level inside a contingence matrix
      real(kind=dp)                        ,            public      :: THRES_CONFUSION
      !// Threshold for convergence in the linear mode heteroscedastic
      real (kind=dp)                                  , public      :: EPS_LINEAR_HETEROSCEDASTIC
      !// Maximum iteration in the linear mode heteroscedastic to avoid infinity loop
      integer                                         , public      :: MAX_LINEAR_ITERATION
      !// get the id of the gpu device to used
      integer                                         , public      :: gpu_device_id = 0
      !// Minimal gamete probability (LDA Haplotype option)
      real(kind=dp)                                   , public      :: prob_seuil_gam = 0.25d0
      !// Number maximum of haplotype allowed
      integer                                         , public      :: NB_HAPLO_PRIOR     = 0
      !// miminum probability of a haplotype in the LD analysis
      real (kind=dp)                                  , public      :: PROB_HAPLO_MIN     = 0.0d0
      !
      integer                                         , public      :: LONG_MIN_IBS        = 0
      !//  Length of the haplotypes to be followed in LD analyses
      integer                                         , public      :: LONGHAP             = 0
      !// array of values initialized from generic keys
      character(len=LEN_L) , dimension(MAXNB_KEY_DEV_GEN)   ,public :: tabDevKey


      logical                     , public        :: simul    ! the context is a simulation or not
      character(len=LENGTH_MAX_FILE), public      :: analyse_file    ! analyse file



!!//   The values of the p_analyses file (referenced by the array keys)
      character(len=LENGTH_MAX_FILE)   , dimension(:), private ,pointer :: values
!!//   The keys defined of the p_analyses file
      character(len=LEN_L)             , dimension(:), private ,pointer :: keys


     contains

      procedure ,public  :: set     => set_parameter_base    !// init the object
      procedure ,private :: read_analyse_file                !// init keys and values arrays
      procedure ,private :: check_unknown_keys
      procedure ,private :: initialize_dev_generic_keys

      procedure ,public  :: get_string_val
      procedure ,public  :: get_int_val
      procedure ,public  :: get_real_val
      procedure ,public  :: get_file_val
      procedure ,public  :: key_exist

      procedure ,public  :: get_summary_panalyse
      procedure ,public  :: help_panalyse

      procedure ,public  :: link    => link_parameter_base
      procedure ,public  :: copy    => copy_parameter_base
      procedure ,public  :: release => release_parameter_base

   end type PARAMETER_BASE


      contains

      subroutine check_unknown_keys(this)
        class(PARAMETER_BASE)       ,intent(in) :: this
        integer :: i,j
        logical :: find
         !check if key are present but not interpreted by qtlmap
        do i=1,size(this%keys)
          if ( trim(this%keys(i))=='' ) cycle
          find = .false.
          do j=1,size(all_key)
            if (this%keys(i) == all_key(j)) then
               find = .true.
               exit
            end if
            ! manage generic keys
            if ( this%keys(i)(1:len(K_KEY_DEV_GEN)) ==  K_KEY_DEV_GEN ) then
               find = .true.
               exit
            end if
          end do

          if ( .not. find ) then
            call log_mess(" UNKNOWN KEY ["//trim(this%keys(i))//"] ",WARNING_DEF)
          end if
        end do

      end subroutine check_unknown_keys


   subroutine set_parameter_base(this,p_analyse,listKeys,listValues,ndim,nkeys)
         class(PARAMETER_BASE)       ,intent(inout) :: this
         character(len=LENGTH_MAX_FILE),intent(in)  :: p_analyse
         integer                ,intent(in)         :: ndim
         character(len=LEN_DEF) ,intent(in), dimension(ndim) :: listKeys,listValues
         integer                ,intent(in)   :: nkeys

         character(len=LEN_W)   :: buf

         integer :: i
         logical :: ok

         !//surcharge des parametres
         call this%read_analyse_file(p_analyse,listKeys,listValues,30,nkeys)

             ! check present keys
        do i=1,size(index_key)
          if (.not. this%key_exist(index_key(i))) then
            call stop_application("Please defined the key:"//trim(index_key(i)))
          end if
        end do

        call this%initialize_dev_generic_keys()

        call this%check_unknown_keys()


        if ( this%key_exist(K_GPU_DEVICE_ID) ) then
            call this%get_string_val(K_GPU_DEVICE_ID,buf)
            this%gpu_device_id  = get_int(buf,ok)
            if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_GPU_DEVICE_ID)//"] expecting int value.")
        end if


        call this%get_string_val(K_MINDAMPHASEPROB,buf)
        this%PRSEUIL  = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_MINDAMPHASEPROB)//"] expecting real value.")

        call this%get_string_val(K_NDMIN,buf)
        this%ndmin    =  get_int(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_NDMIN)//"] expecting int value.")

        call this%get_string_val(K_MINSIREPHASEPROB,buf)
        this%PHPSEUIL = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_MINSIREPHASEPROB)//"] expecting real value.")

        call this%get_string_val(K_CHOLESKY,buf)
        this%SEUIL_CHO = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_CHOLESKY)//"] expecting real value.")

        call this%get_string_val(K_THRES_CONFUSION,buf)
        this%THRES_CONFUSION = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_THRES_CONFUSION)//"] expecting real value.")

        call this%get_string_val(K_LINEAR_CONV,buf)
        this%EPS_LINEAR_HETEROSCEDASTIC = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_LINEAR_CONV)//"] expecting real value.")

        call this%get_string_val(K_MAX_LINEAR_ITERATION,buf)
        this%MAX_LINEAR_ITERATION    =  get_int(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_MAX_LINEAR_ITERATION)//"] expecting int value.")

        call this%get_string_val(K_PSEUILHWE,buf)
        this%PSEUILHWE = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_PSEUILHWE)//"] expecting real value.")

        call this%get_string_val(K_PROB_SEUIL_RECOMB,buf)
        this%PROB_SEUIL_RECOMB = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_PROB_SEUIL_RECOMB)//"] expecting real value.")

        call this%get_string_val(K_NB_HAPLO_PRIOR,buf)
        this%NB_HAPLO_PRIOR = get_int(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_NB_HAPLO_PRIOR)//"] expecting integer value.")

        call this%get_string_val(K_PROB_HAPLO_MIN,buf)
        this%PROB_HAPLO_MIN = get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_PROB_HAPLO_MIN)//"] expecting real value.")

        if ( this%PROB_HAPLO_MIN <1d-8 )  this%PROB_HAPLO_MIN = 1.d-4

        call this%get_string_val(K_LONG_MIN_IBS,buf)
        this%LONG_MIN_IBS = get_int(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_LONG_MIN_IBS)//"] expecting integer value.")

        call this%get_string_val(K_LONGHAP,buf)
        this%LONGHAP = get_int(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_LONGHAP)//"] expecting integer value.")

     end subroutine set_parameter_base

   subroutine copy_parameter_base(mapParamsBase,copyMapParamsBase)
        class(PARAMETER_BASE), intent(in)     :: mapParamsBase
        type(PARAMETER_BASE), intent(inout)  :: copyMapParamsBase

        copyMapParamsBase%NDMIN        = mapParamsBase%NDMIN
        copyMapParamsBase%OPT_SIB           = mapParamsBase%OPT_SIB
        copyMapParamsBase%PSEUILHWE         = mapParamsBase%PSEUILHWE
        copyMapParamsBase%PROB_SEUIL_RECOMB = mapParamsBase%PROB_SEUIL_RECOMB
        copyMapParamsBase%PHPSEUIL          = mapParamsBase%PHPSEUIL
        copyMapParamsBase%PRSEUIL           = mapParamsBase%PRSEUIL
        copyMapParamsBase%SEUIL_CHO         = mapParamsBase%SEUIL_CHO
        copyMapParamsBase%PROB_PHASE_DESC   = mapParamsBase%PROB_PHASE_DESC
        copyMapParamsBase%THRES_CONFUSION   = mapParamsBase%THRES_CONFUSION
        copyMapParamsBase%EPS_LINEAR_HETEROSCEDASTIC    = mapParamsBase%EPS_LINEAR_HETEROSCEDASTIC
        copyMapParamsBase%MAX_LINEAR_ITERATION          = mapParamsBase%MAX_LINEAR_ITERATION
        copyMapParamsBase%tabDevKey                     = mapParamsBase%tabDevKey
        copyMapParamsBase%prob_seuil_gam = mapParamsBase%prob_seuil_gam
        copyMapParamsBase%NB_HAPLO_PRIOR = mapParamsBase%NB_HAPLO_PRIOR
        copyMapParamsBase%PROB_HAPLO_MIN = mapParamsBase%PROB_HAPLO_MIN
        copyMapParamsBase%LONG_MIN_IBS   = mapParamsBase%LONG_MIN_IBS
        copyMapParamsBase%LONGHAP        = mapParamsBase%LONGHAP
        copyMapParamsBase%optim_maxeval        = mapParamsBase%optim_maxeval
        copyMapParamsBase%optim_maxtime        = mapParamsBase%optim_maxtime
        copyMapParamsBase%optim_tolx           = mapParamsBase%optim_tolx
        copyMapParamsBase%optim_tolf           = mapParamsBase%optim_tolf
        copyMapParamsBase%optim_tolg           = mapParamsBase%optim_tolg
        copyMapParamsBase%optim_H_PRECISION    = mapParamsBase%optim_H_PRECISION

        allocate(copyMapParamsBase%values(size(mapParamsBase%values)))
        copyMapParamsBase%values = mapParamsBase%values

        allocate(copyMapParamsBase%keys(size(mapParamsBase%keys)))
        copyMapParamsBase%keys = mapParamsBase%keys

    end subroutine copy_parameter_base

    subroutine link_parameter_base(mapParamsBase,copyMapParamsBase)
        class(PARAMETER_BASE), intent(in)     :: mapParamsBase
        type(PARAMETER_BASE), intent(inout)  :: copyMapParamsBase

        copyMapParamsBase%NDMIN        = mapParamsBase%NDMIN
        copyMapParamsBase%OPT_SIB           = mapParamsBase%OPT_SIB
        copyMapParamsBase%PSEUILHWE         = mapParamsBase%PSEUILHWE
        copyMapParamsBase%PROB_SEUIL_RECOMB = mapParamsBase%PROB_SEUIL_RECOMB
        copyMapParamsBase%PHPSEUIL          = mapParamsBase%PHPSEUIL
        copyMapParamsBase%PRSEUIL           = mapParamsBase%PRSEUIL
        copyMapParamsBase%SEUIL_CHO         = mapParamsBase%SEUIL_CHO
        copyMapParamsBase%PROB_PHASE_DESC   = mapParamsBase%PROB_PHASE_DESC
        copyMapParamsBase%THRES_CONFUSION   = mapParamsBase%THRES_CONFUSION
        copyMapParamsBase%EPS_LINEAR_HETEROSCEDASTIC    = mapParamsBase%EPS_LINEAR_HETEROSCEDASTIC
        copyMapParamsBase%MAX_LINEAR_ITERATION          = mapParamsBase%MAX_LINEAR_ITERATION
        copyMapParamsBase%tabDevKey                     = mapParamsBase%tabDevKey
        copyMapParamsBase%prob_seuil_gam = mapParamsBase%prob_seuil_gam
        copyMapParamsBase%NB_HAPLO_PRIOR = mapParamsBase%NB_HAPLO_PRIOR
        copyMapParamsBase%PROB_HAPLO_MIN = mapParamsBase%PROB_HAPLO_MIN
        copyMapParamsBase%LONG_MIN_IBS   = mapParamsBase%LONG_MIN_IBS
        copyMapParamsBase%LONGHAP        = mapParamsBase%LONGHAP
        copyMapParamsBase%optim_maxeval        = mapParamsBase%optim_maxeval
        copyMapParamsBase%optim_maxtime        = mapParamsBase%optim_maxtime
        copyMapParamsBase%optim_tolx           = mapParamsBase%optim_tolx
        copyMapParamsBase%optim_tolf           = mapParamsBase%optim_tolf
        copyMapParamsBase%optim_tolg           = mapParamsBase%optim_tolg
        copyMapParamsBase%optim_H_PRECISION    = mapParamsBase%optim_H_PRECISION

        copyMapParamsBase%values => mapParamsBase%values
        copyMapParamsBase%keys => mapParamsBase%keys

    end subroutine link_parameter_base


    subroutine release_parameter_base(this)
      class(PARAMETER_BASE)            ,intent(inout)    :: this

      deallocate (this%values)
      deallocate (this%keys)

    end subroutine release_parameter_base

!!****f* m_qtlmap_parameter_file/read_analyse_file
!!  NAME
!!    read_analyse_file
!!  DESCRIPTION
!!    read the parameter analyse user file. fill the internal buffer keys and values.
!!    keys and values are used like a stack (only the last value for the one specific key is considered)
!!  INPUTS
!!   param      : name of the parameter file
!!  NOTES
!!  SOURCE
      subroutine read_analyse_file(this,param,overloadedkeys,overloadedvalues,sizelist,overloadednkeys)
        class(PARAMETER_BASE)                       , intent(inout) :: this
        character(len=LENGTH_MAX_FILE),intent(in) :: param
        integer , intent(in) :: sizelist,overloadednkeys
        character(len=LEN_DEF) ,dimension(sizelist),intent(in):: overloadedkeys,overloadedvalues


        integer  :: ios_id = 20
        integer  :: n ,ios,cln,i,j,chr , nfromfile
        character(len=LEN_W)   :: buf
        character(len=LEN_LINE)   :: buf2
        logical :: ok,find
        real(kind=dp) :: pas_temp

         call file_exist(param)
         open(unit=ios_id,file=trim(param),action="read",form="formatted")

         n = 0
         ios = 0
         ! count the number of line to record
         do while ( ios == 0)
            read(ios_id,*,iostat=ios) buf2
            i = index (buf2,'#')
            if ( i == 1 ) cycle
            if ( ios == 0 .and. ( trim(buf2)/='') ) then
              n = n + 1
            end if
         end do

         nfromfile = n
         if ( overloadednkeys > 0 ) n = n + overloadednkeys
         !allocate structure
         allocate (this%values(n))
         allocate(this%keys(n))
         this%keys=''

         ! init structures
         rewind(ios_id)
         cln = 0! it s the current line
         n = 0
         ios=-2
         do while ( ios /= -1 )
            call GET (ios_id,buf2,maxlen=LEN_LINE,iostat=ios)
            !remove comments
            i = index (buf2,'#')

            if ( i /= 0) then
              if ( i == 1 ) buf2=''
              if ( i > 1) buf2=trim(buf2(:i-1))
            end if
            !probleme case
            if (trim(buf2)=='' .or. ios /= 0 ) cycle

            cln = cln+1
            i = index (buf2,'=')! where is =
            if ( i == 0) then
                call stop_on_error(ios,param,cln,'bad expression of "key=value"');
            end if

              ! search if the key if exist...
              do j=1,n
                 if ( this%keys(j) == trim(buf2(1:i-1)) ) exit
              end do

              if ( j <= n ) then
                this%keys(j) = trim(buf2(1:i-1))
                this%values(j) = trim(buf2(i+1:))
              else
                n = n+1
                this%keys(n) = trim(buf2(1:i-1))
                this%values(n) = trim(buf2(i+1:))
              end if

              call log_mess("KEY ["//trim(this%keys(n))//"] = ["//trim(this%values(n))//"]",DEBUG_DEF)

         end do

        close(ios_id)

        do i=1,overloadednkeys
         do j=1,n
              if ( this%keys(j) == overloadedkeys(i) ) exit
         end do
         if ( j <= n ) then
             this%keys(j) = overloadedkeys(i)
             this%values(j) = overloadedvalues(i)
         else
             n = n+1
             this%keys(n) = overloadedkeys(i)
             this%values(n) = overloadedvalues(i)
         end if

         call log_mess("KEY ["//trim(this%keys(n))//"] = ["//trim(this%values(n))//"]",DEBUG_DEF)

        end do

      end subroutine read_analyse_file
!!***


!!****f* m_qtlmap_parameter_file/initialize_qtlmap_context
!!  NAME
!!    initialize_qtlmap_context
!!  DESCRIPTION
!!    Initialize all public variables option come from the module m_qtlmap_data and file name to read
!!
!!  INPUTS
!!   sim         : if the context is in the simulation case
!!   is_permute  : if dynamic option permutation mode is activate
!!
!!  OUTPUTS
!!   parsim      : the name of the simulation file
!!  NOTES
!!  SOURCE
!      subroutine initialize_qtlmap_context(this,dataset,sim,is_permute)
!        type(QTLMAP_DATASET)       ,intent(inout)    :: dataset
!        integer                       ,intent(in)    :: sim            ! 0 -> analyse, 1 -> simul, 2-> analyse+simul (cuda)
!        logical                      ,intent(in)     :: is_permute
!
!        integer  :: ios_id = 20
!        integer  :: n ,ios,cln,i,j,chr
!        character(len=LEN_W)   :: buf
!        character(len=LEN_S),dimension(100)   :: chromo_t
!        character(len=LEN_LINE)   :: buf2
!        logical :: ok,find
!        real(kind=dp) :: pas_temp
!
!         ! check present keys
!        do i=1,size(index_key)
!          if (.not. key_exist(index_key(i))) then
!            call stop_application("Please defined the key:"//trim(index_key(i)))
!          end if
!        end do
!
!        call initialize_dev_generic_keys(dataset)
!
!        !check if key are present but not interpreted by qtlmap
!        do i=1,size(keys)
!          if ( trim(keys(i))=='' ) cycle
!          find = .false.
!          do j=1,size(all_key)
!            if (keys(i) == all_key(j)) then
!               find = .true.
!               exit
!            end if
!            ! manage generic keys
!            if ( keys(i)(1:len(K_KEY_DEV_GEN)) ==  K_KEY_DEV_GEN ) then
!               find = .true.
!               exit
!            end if
!          end do
!
!          if ( .not. find ) then
!            call log_mess(" ---------------------------------- ",WARNING_DEF)
!            call log_mess(" UNKNOWN KEY ["//trim(keys(i))//"] ",WARNING_DEF)
!            call log_mess(" ---------------------------------- ",WARNING_DEF)
!          end if
!        end do

!        call get_string_val(K_OUTPUT,dataset%files%resul)
!        call get_string_val(K_SUMM,dataset%files%summary)
!
!        if ( sim == 0 .or. sim == 2 ) then
!          call get_string_val(K_LRTSIRE,dataset%files%resp)
!          call get_string_val(K_LRTDAM,dataset%files%resm)
!          if ( key_exist(K_PDED) ) then
!            call get_string_val(K_PDED,dataset%files%pdedf)
!          end if
!          if ( key_exist(K_PDECPLE) ) then
!            call get_string_val(K_PDECPLE,dataset%files%pdecouple)
!          end if
!          if ( key_exist(K_PATEFF) ) then
!            call get_string_val(K_PATEFF,dataset%files%pateff)
!          end if
!          if ( key_exist(K_FREQALL) ) then
!            call get_string_val(K_FREQALL,dataset%files%out_freqall)
!          end if
!          if ( key_exist(K_MATEFF) ) then
!            call get_string_val(K_MATEFF,dataset%files%mateff)
!          end if
!          if ( key_exist(K_COEFFDA) ) then
!            call get_string_val(K_COEFFDA,dataset%files%coeffda)
!          end if
!          if ( key_exist(K_GRID2QTL) ) then
!            call get_string_val(K_GRID2QTL,dataset%files%grid2qtl)
!          end if
!          if ( key_exist(K_PHASES) ) then
!            call get_string_val(K_PHASES,dataset%files%out_phases)
!          end if
!          if ( key_exist(K_HAPLOTYPES) ) then
!            call get_string_val(K_HAPLOTYPES,dataset%files%out_haplotypes)
!          end if
!          if ( key_exist(K_INFORMATIVITY) ) then
!            call get_string_val(K_INFORMATIVITY,dataset%files%file_informativity)
!          end if
!        end if

        !! Optimization initialisation
        !! -----------------------

!        call get_string_val(K_OPTIM_MAXEVAL,buf)
!        optim_maxeval = get_int(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_MAXEVAL)//"] expecting integer value.")
!
!        call get_string_val(K_OPTIM_MAXTIME,buf)
!        optim_maxtime    =  get_real(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_MAXTIME)//"] expecting int value.")
!
!        call get_string_val(K_OPTIM_TOLX,buf)
!        optim_tolx    =  get_real(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_TOLX)//"] expecting int value.")
!
!        call get_string_val(K_OPTIM_TOLF,buf)
!        optim_tolf    =  get_real(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_TOLF)//"] expecting int value.")
!
!        call get_string_val(K_OPTIM_TOLG,buf)
!        optim_tolg    =  get_real(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_TOLG)//"] expecting int value.")
!
!        call get_string_val(K_OPTIM_H_PREC,buf)
!        optim_H_PRECISION    =  get_real(buf,ok)
!        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_OPTIM_H_PREC)//"] expecting int value.")


        !! Chromome initialisation
        !! -----------------------


        ! INPUT FILES
        !------------
!        if ( key_exist(K_MAP) ) then
!             call get_string_val(K_MAP,dataset%files%in_carte)
!             call file_exist(dataset%files%in_carte)
!             dataset%files%mapFileDefined = .true.
!         else
!              dataset%files%mapFileDefined = .false.
!         end if
!        if ( key_exist(K_GENEA) ) then
!            call get_string_val(K_GENEA,dataset%files%in_genea)
!            call file_exist(dataset%files%in_genea)
!            dataset%files%geneaFileDefined = .true.
!        else
!            dataset%files%geneaFileDefined = .false.
!        end if
!
!        if ( key_exist(K_RACE) ) then
!            dataset%files%raceFileDefined=.true.
!            call get_string_val(K_RACE,dataset%files%in_race)
!            call file_exist(dataset%files%in_race)
!        end if
!        dataset%files%in_perfs = ''
!        if ( key_exist(K_TRAITS) ) then
!            call get_string_val(K_TRAITS,dataset%files%in_perfs(1))
!            call file_exist(dataset%files%in_perfs(1))
!            dataset%files%traitsFileDefined = .true.
!            !The model have to be declared...
!            if ( key_exist(K_MODEL) ) then
!              call get_string_val(K_MODEL,dataset%files%in_param_ef)
!              call file_exist(dataset%files%in_param_ef)
!            else
!              call stop_application("key ["//trim(K_TRAITS)//"] is defined. The key ["//trim(K_MODEL)//"] have to be set.")
!            end if
!        else
!            dataset%files%traitsFileDefined = .false.
!        end if
!
!        if ( key_exist(K_GENOTYPE) ) then
!            call get_string_val(K_GENOTYPE,dataset%files%in_typage)
!            call file_exist(dataset%files%in_typage)
!            dataset%files%genotypeFileDefined = .true.
!        else
!            dataset%files%genotypeFileDefined = .false.
!        end if

!       if ( sim == 1 .or. sim == 2  ) then
!        if ( key_exist(K_PARAMSIM) ) then
!            call get_string_val(K_PARAMSIM,dataset%files%in_parsim)
!            call file_exist(dataset%files%in_parsim)
!        else if ( .not. is_permute) then
!            !call stop_application("key ["//trim(K_PARAMSIM)//"] is not defined.")
!            dataset%files%in_parsim = "" ! by default we consider all traits defined in model
!        end if
!
!      ! output lrt max file is now optional
!
!        if ( key_exist(K_OUTSIM) ) then
!            call get_string_val(K_OUTSIM,dataset%files%simula)
!      !  else
!      !      call stop_application("key ["//trim(K_PARAMSIM)//"] is defined. The key ["//trim(K_OUTSIM)//"] have to be set.")
!        end if


!
!       else
!        !    t_imp = .false.
!       end if

        ! check base....

!        if ( sim == 0 .or. sim == 2 ) then
!           if ( .not. dataset%files%mapFileDefined ) &
!             call stop_application("Analysis case : You have to defined the key ["//trim(K_MAP)//"]")
!           if ( .not. dataset%files%geneaFileDefined ) &
!             call stop_application("Analysis case : You have to defined the key ["//trim(K_GENEA)//"]")
!           if ( .not. dataset%files%traitsFileDefined ) &
!             call stop_application("Analysis case : You have to defined the key ["//trim(K_TRAITS)//"]")
!           if ( .not. dataset%files%genotypeFileDefined ) &
!             call stop_application("Analysis case : You have to defined the key ["//trim(K_GENOTYPE)//"]")
!        end if
!      end subroutine initialize_qtlmap_context
!!***


!!****f* m_qtlmap_parameter_file/initialize_dev_generic_keys
!!  NAME
!!    initialize_dev_generic_keys
!!  DESCRIPTION
!!    initilialized the array m_qtlmap_data/tabDevKey from gebneric key K_KEY_DEV_GEN
!!  NOTES
!!  SOURCE
      subroutine initialize_dev_generic_keys(this)
        class(PARAMETER_BASE) , intent(inout)        :: this
        integer :: i
        character(len=LEN_L) :: buf

        this%tabDevKey=""
        do i=1,MAXNB_KEY_DEV_GEN
           buf=trim(K_KEY_DEV_GEN)//trim(str(i))
           if ( this%key_exist( buf ) ) then
               call this%get_string_val(buf,this%tabDevKey(i))
           end if
        end do

      end subroutine initialize_dev_generic_keys

!!    Get value according to the key give in parameter
      subroutine get_string_val(this,keysearch,value)
         class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)     :: keysearch
         character(len=*) ,intent(out)    :: value

         integer                :: i

         value=''

         if (associated(this%keys)) then
          do i=1,size(this%keys)
           if ( trim(this%keys(i)) == trim(keysearch)) then
             value = this%values(i)
             return
           end if
          end do
         end if

         do i=1,size(index_key)
           if ( trim(index_key(i)) == trim(keysearch)) then
             if ( trim(default_values(i)) /= '') then
                 value = default_values(i)
                 return
             end if
           end if
         end do

      end subroutine get_string_val

!!
      function get_int_val(this,keysearch) result(value)
        class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)     :: keysearch
         integer    :: value
         character(len=LEN_W)   :: buf
         logical :: ok

         call this%get_string_val(keysearch,buf)
         value  = get_int(buf,ok)
      end function get_int_val

!!
      function get_real_val(this,keysearch) result(value)
        class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)     :: keysearch
         real(kind=dp)    :: value
         character(len=LEN_W)   :: buf
         logical :: ok

         call this%get_string_val(keysearch,buf)
         value  = get_real(buf,ok)
      end function get_real_val

!!
      function get_file_val(this,keysearch) result(value)
        class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)     :: keysearch
         character(len=LENGTH_MAX_FILE)    :: value

         call this%get_string_val(keysearch,value)

      end function get_file_val

!!****f* m_qtlmap_parameter_file/key_exist
!!  NAME
!!    key_exist
!!  DESCRIPTION
!!    True if the key exist in the parameter user file or if the key have a default value
!!  NOTES
!!  SOURCE
      function key_exist(this,keysearch) result(vexist)
         class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)   :: keysearch

         logical                :: vexist
         integer                :: i

         if (associated(this%keys)) then
          do i=1,size(this%keys)
           if ( trim(this%keys(i)) == trim(keysearch) ) then
               if ( trim(this%values(i)) /= '' ) then
                 vexist=.true.
                 return
               end if

           end if
          end do
         end if

         !if a default value exist...
          do i=1,size(index_key)
           if ( trim(index_key(i)) == trim(keysearch)) then
               if ( trim(default_values(i)) /= '') then
                  vexist=.true.
               else
                  vexist=.false.
               end if
               return
           end if
         end do

         vexist=.false.
         return

      end function key_exist
!!***

!!****f* m_qtlmap_parameter_file/key_exist_file
!!  NAME
!!    key_exist_file
!!  DESCRIPTION
!!    True if the key exist in the parameter user file
!!  NOTES
!!  SOURCE
      function key_exist_file(this,keysearch) result(vexist)
         class(PARAMETER_BASE) , intent(in) :: this
         character(len=*) ,intent(in)   :: keysearch

         logical                :: vexist
         integer                :: i

         if (associated(this%keys)) then
          do i=1,size(this%keys)
           if ( trim(this%keys(i)) == trim(keysearch) ) then
               if ( trim(this%values(i)) /= '' ) then
                 vexist=.true.
                 return
               end if

           end if
          end do
         end if

         vexist=.false.
         return

      end function key_exist_file
!!***

!!****f* m_qtlmap_parameter_file/help_panalyse
!!  NAME
!!    help_panalyse
!!  DESCRIPTION
!!    print a help (all keys availables).
!!  NOTES
!!  SOURCE
      subroutine help_panalyse(this)
              class(PARAMETER_BASE) , intent(in) :: this
              character(len=LEN_L) :: headkey ,headvalues,form,buf
              character(len=LEN_S) :: opt

              integer :: i

              headkey = "KEY"
              headvalues="                     DEFINITION"
              opt = "DEFAULT"
              form='(a60,a50,4x,a10)'

              print form,headkey,headvalues,opt
              print form,"------------------------------","--------------------------------","--------------"
              print form,K_OUTPUT,"Output main file"
              print form,K_LRTSIRE,"Output file paternal effects"
              print form,K_LRTDAM,"Output file maternal effects"

              print form,K_PDED,"Grand parental segment transmission marginal probabilities"
              print form,K_PDECPLE,"Grand parental segment transmission joint probabilities"
              print form,K_PATEFF,"Sire QTL effect estimations"
              print form,K_MATEFF,"Dam QTL effect estimations"

              print form,K_SUMM,"Output summary file"
              print form,K_CHROM,"Chromosomes used in the analysis : chr1,chr2,.."
              print form,K_NDMIN,'Minimal number of progeny by dam'

              call this%get_string_val(K_STEP,buf)
              print form,K_STEP,"Chromosomic segment exploration steps in Morgan",buf

              call this%get_string_val(K_MINDAMPHASEPROB,buf)
              print form,K_MINDAMPHASEPROB, 'Minimal dam phase probability',buf

              call this%get_string_val(K_MINSIREPHASEPROB,buf)
              print form,K_MINSIREPHASEPROB,'Minimal sire phase probability',buf

              call this%get_string_val(K_UNKNOWN_GENO,buf)
              print form,K_UNKNOWN_GENO,'Unknown genotype value',buf

              call this%get_string_val(K_CHOLESKY,buf)
              print form,K_CHOLESKY,"coeff cholesky decomposition",buf
              print form,K_MAP,"Input map file"
              print form,K_GENEA,"Input genealogy file"
              print form,K_TRAITS,"Input traits file"
              print form,K_GENOTYPE,"Input genotype file"
              print form,K_MODEL,"Input model description of traits"
              print form,K_PARAMSIM,"Input simulation parameters"
              print form,K_OUTSIM,"Output for each simulation (Position and maxLRT for each traits)"
              print form,K_PHASES,"Output phases file"
              print form,K_HAPLOTYPES,"Output haplotype file"

              call this%get_string_val(K_THRES_CONFUSION,buf)
              print form,K_THRES_CONFUSION,"Threshold to test confusion betwwen level inside a contingence matrix",buf
              call this%get_string_val(K_PSEUILHWE,buf)
              print form,K_PSEUILHWE,"Threshold to check the equilibrium of marker transmission within each family",buf
              call this%get_string_val(K_LINEAR_CONV,buf)
              print form,K_LINEAR_CONV,"Threshold for convergence in the linear mode heteroscedastic",buf
              call this%get_string_val(K_MAX_LINEAR_ITERATION,buf)
              print form,K_MAX_LINEAR_ITERATION,"Maximum iteration in the linear mode heteroscedastic to avoid infinity loop",buf

              call this%get_string_val(K_PROB_SEUIL_RECOMB,buf)
              print form,K_PROB_SEUIL_RECOMB,"",buf

           end subroutine help_panalyse

!!    get the both vector of key and values defined by the user
          subroutine get_summary_panalyse(this,my_index_key,values,n)
             class(PARAMETER_BASE) , intent(in) :: this
             character(len=LEN_L) ,dimension(NUMBER_ALL_KEYS),intent(out) :: my_index_key
             character(len=LEN_L) ,dimension(NUMBER_ALL_KEYS),intent(out) :: values
             integer                                         ,intent(out) :: n

             character(len=LEN_L) :: buf
             integer :: i

             n=0
             do i=1,size(all_key)
                 if ( this%key_exist(all_key(i)) ) then
                   n=n+1
                   call this%get_string_val(all_key(i),buf)
                   my_index_key(n)=trim(all_key(i))
                   values(n)=trim(buf)
                 end if
             end do

          end subroutine get_summary_panalyse
!!***

end module m_qtlmap_type_parameter
