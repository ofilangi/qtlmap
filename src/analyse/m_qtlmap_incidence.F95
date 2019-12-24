!!****m* ANALYSE/m_qtlmap_incidence
!!  NAME
!!    m_qtlmap_incidence
!!
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
module m_qtlmap_incidence
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_solve_linear
    use m_qtlmap_optimization
    use m_qtlmap_output_handler
    use m_qtlmap_math
    use m_qtlmap_analyse_gen , only : set_ntnivmax, set_phase_pdd
    use m_qtlmap_haplotype_ldla

    implicit none

    integer   , parameter   ,public    :: MASK_INIT_MOD_INC             = INT(b'0000000000000000')
    !Exclusifs
    integer   , parameter   ,public    :: KIND_HOMOSCEDASTIC            = INT(b'0001000000000000')
    integer   , parameter   ,public    :: KIND_HETEROSCEDASTIC          = INT(b'0010000000000000')

    !Exclusifs
    integer   , parameter   ,public    :: MODEL_UNITRAIT_POLYGENIC      = INT(b'0000000100000000')
    integer   , parameter   ,public    :: MODEL_MULTITRAIT_POLYGENIC    = INT(b'0000001000000000')
    integer   , parameter   ,public    :: MODEL_MULTITRAIT_POLYGENIC_LU = INT(b'0000010000000000')
    integer   , parameter   ,public    :: MODEL_UNITRAIT_RELATIONSHIP   = INT(b'0000100000000000')

    !Exclusifs
    integer   , parameter   ,public    :: RESOLUTION_LINEAR             = INT(b'0000000000010000')
    integer   , parameter   ,public    :: RESOLUTION_OPTIMISATION       = INT(b'0000000000100000')


!!****t* m_qtlmap_incidence/type_effect_contingence
!!  NAME
!!    type_effect_contingence
!!  DESCRIPTION
!!    Determines the effect to include inside contingence matrix (by defaut qtl only for a LA analysis) with the an order of insertion
!!  SOURCE
    type type_effect_contingence
      integer :: general   = 1
      integer :: polygenic = 2
      integer :: qtl       = 3
      integer :: haplo     = 0
      integer :: nuis      = 4
      integer :: traits    = 0
      integer :: inter     = 0

      integer :: nbeff     = 7 ! have to be inscrease if other kind of effects

    end type type_effect_contingence
!!***

!!****d* m_qtlmap_constant/CONTINGENCE_QTL
!!  NAME
!!   CONTINGENCE_QTL
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_QTL        = type_effect_contingence(1,3,2,0,4,0,0)
!!****d* m_qtlmap_constant/CONTINGENCE_HAPLO
!!  NAME
!!   CONTINGENCE_HAPLO
!!  DESCRIPTION
!!   Description of a contigence matrix for LD analysis
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_HAPLO      = type_effect_contingence(1,3,0,2,4,0,0) !(1,3,0,2,4,0,0) !(1,0,0,2,3,0,0)
!!****d* m_qtlmap_constant/CONTINGENCE_HAPLO_QTL
!!  NAME
!!   CONTINGENCE_HAPLO_QTL
!!  DESCRIPTION
!!   Description of a contigence matrix for LDLA analysis
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_HAPLO_QTL  = type_effect_contingence(1,4,3,2,5,0,0)
!!****d* m_qtlmap_constant/CONTINGENCE_QTL_TRAITS
!!  NAME
!!   CONTINGENCE_QTL_TRAITS
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis with covariates traits effects
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_QTL_TRAITS = type_effect_contingence(1,3,2,0,4,5,0)
!!****d* m_qtlmap_constant/CONTINGENCE_QTL_INTER
!!  NAME
!!   CONTINGENCE_QTL_INTER
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis with an interaction
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_QTL_INTER  = type_effect_contingence(1,5,2,0,4,0,3)

!!****d* m_qtlmap_constant/CONTINGENCE_QTL_WITHOUT_POLY
!!  NAME
!!   CONTINGENCE_QTL_WITHOUT_POLY
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis without polygenic effects (model animal)
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_QTL_WITHOUT_POLY  = type_effect_contingence(1,0,2,0,3,0,0)

!!****d* m_qtlmap_constant/CONTINGENCE_QTL_WITHOUT_POLY
!!  NAME
!!   CONTINGENCE_QTL_WITHOUT_POLY
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis without polygenic effects (model animal)
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_HAPLO_WITHOUT_POLY  = type_effect_contingence(1,0,0,2,3,0,0)

!!****d* m_qtlmap_constant/CONTINGENCE_RACE
!!  NAME
!!   CONTINGENCE_QTL
!!  DESCRIPTION
!!   Description of a contigence matrix for LA analysis
!!  NOTES
!!
!!***
   type(type_effect_contingence) ,parameter ::  CONTINGENCE_RACE        = type_effect_contingence(1,3,2,0,4,0,0)

!!****t* m_qtlmap_incidence/family_type
!!  NAME
!!    family_type
!!
!!  DESCRIPTION
!!   Basic type with the proba of transmission of a family
!!  SOURCE
    type family_type
      integer :: firstKd = -1 ! correspond a l index du premier descendant (equivalent ligne mat incidence)
      integer :: lastKd  = -2 ! correspond a l index du dernier descedant (equivalent ligne mat incidence)
      real(kind=dp) , dimension(:,:,:),pointer :: ppt => null()!proba de trsanmition recu par le pere
      real(kind=dp) , dimension(:,:,:),pointer :: pmt => null()! proba de transmision recu par la mere si estim
     ! integer   ,dimension(:) ,pointer         :: ntniv_dam !indice dans la matrice d incidence correspondant a l effet pere et le qtl
    end type family_type
!!***

!!****t* m_qtlmap_incidence/SireFamily_type
!!  NAME
!!    SireFamily_type
!!
!!  DESCRIPTION
!!   Complex type with the proba of transmission of a sire family
!!  SOURCE
    type SireFamily_type
      type(family_type)                        :: half_sib ! famille complete de demi-frere (ligne ref incidence)
      type(family_type) , dimension(:),pointer :: full_sib => null()! famille de plein frere (ligne ref incidence)
      real(kind=dp)     , dimension(:),pointer :: sig0 => null()!variance du caractere dans la famille de pere
      real(kind=dp)     , dimension(:),pointer :: xmu0 => null()!moyenne du caractere dans la famille de pere
      real(kind=dp) , dimension(:,:),pointer   :: ppt => null()!proba de transmition recu par le pere
      !integer       ,dimension(:) ,pointer     :: ntniv_sire !indice dans la matrice d incidence correspondant a l effet pere et le qtl
    end type SireFamily_type
!!***


!!****t* m_qtlmap_incidence/DATASET_TYPE
!!  NAME
!!    DATASET_TYPE
!!
!!  DESCRIPTION
!!   description of the user dataset
!!  SOURCE
    type DATASET_TYPE
        ! dec 2012 => link to use QTLMap datatype
        type(QTLMAP_DATASET)                 ,pointer            :: data             => NULL()
        type(PDD_BUILD)                      ,pointer            :: spt              => NULL()
        !Index trait
        integer                                                  :: ncar             = 0
        !Struct of Families (missing data are removes) with pdd for full and half sib families
        type(SireFamily_type) , dimension(:) , pointer           :: lSires           => NULL()
        !Number of animals referenced in the contingence matrix (missing data removed)
        integer                                                  :: nkd              = 0
        !Maximum number recorded in a family full sib
        integer                                                  :: nkd_max_by_fam   = 0
        !Performance without missing data : size : nkd
        real(kind=dp)      , dimension(:,:),   pointer           :: Y                => NULL()
        integer            , dimension(:,:),   pointer           :: YDISCRETORD      => NULL()
        !Censored data without missing data : size nkd
        real(kind=dp)      , dimension(:,:),   pointer           :: CD               => NULL()
        !Covariance of traits (H0) dim : ncar,ncar
        real(kind=dp)      , dimension(:,:),   pointer           :: cov              => NULL()
        !Residual correlation (H0) dim : ncar,ncar
        real(kind=dp)      , dimension(:,:),   pointer           :: rhoi             => NULL()
        !Standart deviation (H0), dim:ncar
        real(kind=dp)      , dimension(:)  ,   pointer           :: sig              => NULL()
        !Mean (H0), dim:ncar
        real(kind=dp)      , dimension(:)  ,   pointer           :: xmu              => NULL()

    end type DATASET_TYPE
!!***

!!****t* m_qtlmap_incidence/INCIDENCE_TYPE
!!  NAME
!!    INCIDENCE_TYPE
!!
!!  DESCRIPTION
!!   description of an effect enumered in the INCIDENCE_TYPE
!!  SOURCE
    type INCIDENCE_TYPE
        !Index trait
        integer                                                :: ic       = 0
        !The associated dataset
        type(DATASET_TYPE)                      ,pointer       :: dataset  => NULL()
        !number of qtl to estim
        integer                                                :: nqtl     = 0
        ! Number of level
        integer                                                :: ntniv    = 0
        ! Number of effect
        integer                                                :: nteff    = 0
        ! tab of effect to print in the eqtl format
        logical            ,dimension(20)                      :: eqtl_print    = .true.
        integer                                                :: ntnivmax=0
        !Maximum of effect defined in the model : General mean, polygenic sire,polygenic dam, qtl sire,....
        integer                                                :: nteffmax=0
        ! Number of level which can be estimate
        integer                                                :: nbnivest = 0
        ! Number of qtl*interaction level by sire or dam
        integer                                                :: ntlev    = 0
        ! Number of level for fixed effect
        integer                                                :: nbniv    = 0
        !description of each effect
        type(DESC_EFFECT_TYPE),dimension(:),pointer            :: desc     => NULL()
        !estimability of each level ntniv
        logical        , dimension(:) , pointer                :: vecsol   => NULL()
        !precision of each level estim
        real (kind=dp)       ,dimension(:),pointer             :: precis   => NULL()
        !dim:nqtl : referenced qtl sire ntniv in the contingence matrix
        integer              ,dimension(:),pointer             :: ntniv_qtlsires => NULL()
        !dim:nqtl : referenced qtl dams ntniv in the contingence matrix
        integer              ,dimension(:),pointer             :: ntniv_qtldams  => NULL()
        real (kind=dp)       ,dimension(:) ,pointer            :: pp             => NULL() !buf pour les proba de transmissions
        real (kind=dp)       ,dimension(:) ,pointer            :: pm             => NULL() !buf pour les proba de transmissions
        real (kind=dp)       ,dimension(:) ,pointer            :: pp2             => NULL() !buf pour les 2eme proba de transmissions dans interaction
        real (kind=dp)       ,dimension(:) ,pointer            :: pm2             => NULL() !buf pour les 2eme proba de transmissions dans interaction
        !dim:ntnivmax correspondance si estimable de la colonne dans la matrice xinc reduite
        integer              ,dimension(:),pointer             :: corr_niv_nivb  => NULL()
        !optim struc
        real (kind=dp) , dimension(:) ,pointer                 :: borni            => null()
        real (kind=dp) , dimension(:) ,pointer                 :: borns            => null()! les bornes de l optimisation
        real (kind=dp) , dimension(:) ,pointer                 :: par              => null()
        real(kind=dp)  , dimension(:) ,pointer                 :: fperemax             => null()
        real(kind=dp)  , dimension(:) ,pointer                 :: fmeremax             => null()

        logical                                                :: isBiAllellic     = .false.

        !optimisation qui fait office de corniv dans le modlin originale...
        ! on stocke pour chque kd, les indices des elements des elements non nul....
        integer        , dimension(:,:) ,pointer               :: nonull           => null()
        integer        , dimension(:)   ,pointer               :: nbnonull         => null()! nombre d element non nul
        character(len=LEN_W) ,dimension(:,:,:),pointer         :: incidenceSireHaplo => null()
        character(len=LEN_W) ,dimension(:,:,:),pointer         :: incidenceSireHaploRace => null()

    end type INCIDENCE_TYPE
!!***

    type WORK_CUDA
        type(QTLMAP_DATASET)           ,pointer  :: data                 => null()
        type(PDD_BUILD)                ,pointer  :: spt                  => null()
        integer       , dimension(:,:) ,pointer  :: lchr => null() ,lnpos  => null()
        type(INCIDENCE_GEN_STRUCT)     ,pointer  :: workstruct           => null()
        real(kind=dp) , dimension(:,:) ,pointer  :: xinc                 => null()
        type(INCIDENCE_TYPE)           ,pointer  :: incidenceDesc        => null()
        integer                                  :: nLevelVar,nLevelFix
        integer , dimension(:)         ,pointer  :: corrLevel            => null()
    end type WORK_CUDA

!!****t* m_qtlmap_incidence/INCIDENCE_GEN_STRUCT
!!  NAME
!!    INCIDENCE_GEN_STRUCT
!!
!!  DESCRIPTION
!!    structure pour developper des methodes generiques
!!      ce type est un agrega des structures (non commune au mode lineaire et optim) necessaires pour la solution de la vraissemblance a un point n
!!  SOURCE
    type INCIDENCE_GEN_STRUCT

       !!!!!  initializing by the caller !!!!!
       !**************************************
       integer                                      :: type_model       = MASK_INIT_MOD_INC ! linear homo,hetero,optim,...
       type(type_effect_contingence)                :: effects
       logical                                      :: performConfusion = .false.
       logical                                      :: performTestNuis  = .false.
       logical                                      :: linear           = .false.
       logical                                      :: hdam             = .false.
       logical                                      :: biq              = .false. ! biallelic effect - no by default
 !      logical                                      :: nopoly            = .false.
       integer                                      :: nqtl             = -1
       integer                                      :: hypothesis       = -1
       !id nteff inside xinc
       integer          ,dimension(:), pointer      :: listnteff        => null()
       ! index of the interaction effect
       integer                                      :: nteffinter       = -1
       !id nteff inside xinc
       integer          ,dimension(:), pointer      :: listnteffhap     => null()
       !linear information
       !------------------
       !0-->maxiqtl,np
       real (kind=dp) ,dimension(:,:,:) ,pointer      :: sigsquare        => null()
       !maxqtl,ncar,ncar
       real (kind=dp)   ,dimension(:,:,:) ,pointer    :: rhoi  => null()

       !F0,F1,etc... dim:nqtl
       real (kind=dp)   ,dimension(:,:) ,pointer      :: fnqtlsires         => null()
       real (kind=dp)   ,dimension(:,:) ,pointer      :: fnqtldams          => null()
       !ADDITIONAL OUTPUTS INFORMATIONS
       !-------------------------------
       ! alert confusion (fill if performConfusion is true )
       type(CORR_ALERT_TYPE) ,dimension(:,:),pointer    :: alertQtl         => null()
       ! size of alertQtl
       integer               ,dimension(:),pointer      :: nalert           => null()
       ! Maximum correlation finding between qtl and other effect
       real (kind=dp)        ,dimension(:),pointer      :: corrmax

       ! list of nuisance effect and test (fill if performTestNuis is true)
       type(TEST_NUISANCES_TYPE) ,dimension(:),pointer  :: listtestnuis => null()
       ! size of listtestnuis
       integer                                          :: ntest        = 0
       ! index to each trait to init as a covariate
       integer                                          :: startTraitAsCov=0,endTraitAsCov=-1
       ! Quantity : M = (I - (I + lambda*A^-1) = Animal Modele =
       real (kind=dp)   ,dimension(:,:) ,pointer        :: M             => null()

    end type INCIDENCE_GEN_STRUCT
!!***

!!****t* m_qtlmap_incidence/POSITION_LRT_INCIDENCE
!!  NAME
!!    POSITION_LRT_INCIDENCE
!!
!!  DESCRIPTION
!!    structure pour developper des methodes generiques
!!      ce type est un agrega des structures (non commune au mode lineaire et optim) necessaires pour la solution de la vraissemblance a un point n
!!  SOURCE
    type POSITION_LRT_INCIDENCE
      !Information bind to the current position
      !----------------------------------------
      !nqtl
      integer        ,dimension(:) ,pointer        :: listChr      => null()
      integer        ,dimension(:) ,pointer        :: listN        => null()
      !output info fill by the model
      real(kind=dp) ,dimension(:,:)  ,pointer      :: lrtSires     => null()
      real(kind=dp) ,dimension(:,:)  ,pointer      :: lrtDams      => null()

    end type POSITION_LRT_INCIDENCE
!!***

!!****d* m_qtlmap_incidence/MAX_KD_PRINT
!!  NAME
!!   MAX_KD_PRINT
!!  DESCRIPTION
!!   Constant : maximum of progeny to print for the debug mode
!!  NOTES
!!   debug_write_incidence
!!***
    !Constante : maximum of progeny to print
    integer , parameter  , private                            :: MAX_KD_PRINT=20

    save

    public :: init_analyse_linear
    public :: end_analyse_linear

    public :: init_incidence
    public :: end_incidence

    public :: init_dataset_1car
    public :: end_dataset

    public :: init_position
    public :: end_position

    public :: fill_nonull_elements
    public :: matmul_incidence

    !routines to build contingence matrix
    public :: add_general_mean
    public :: add_polygenic_effect
    public :: add_qtleffect
    public :: add_qtlinteraction
    public :: add_effcov

    !changing parameter at the position
    public :: change_qtleffect
    public :: change_interaction_effect

    public :: copy_incidence_desc
    public :: release_copy_incidence_desc

    public :: set_parameter
    public :: is_parameter


    contains

!!****f* m_qtlmap_incidence/set_parameter
!!  NAME
!!    set_parameter
!!  DESCRIPTION
!!    Set a characteristic for the current analysis
!!  NOTES
!!  SOURCE
    subroutine set_parameter(newparameter,allparameter)
       integer , intent(in)    :: newparameter
       integer , intent(inout) :: allparameter
       allparameter = ior(newparameter,allparameter) ! ou logic bit à bit pour les caracteristiques de l'analyse courante (homoscedastic,hetero,lineaire,optim,model animal,....)
    end subroutine set_parameter
!!***
!!****f* m_qtlmap_incidence/is_parameter
!!  NAME
!!    is_parameter
!!  DESCRIPTION
!!    give the characteristique of the current analysis
!!  NOTES
!!  SOURCE
    function is_parameter(testparam,allparameter) result(res)
      integer , intent(in)    :: testparam,allparameter
      logical :: res
      res = iand(testparam,allparameter) /= MASK_INIT_MOD_INC
    end function
!!***

!!****f* m_qtlmap_incidence/init_analyse_linear
!!  NAME
!!    init_analyse_linear
!!  DESCRIPTION
!!    initialize the linear mode
!!  NOTES
!!  SOURCE
    subroutine init_analyse_linear(dataset,spt)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         type(PDD_BUILD)            ,intent(in) :: spt

         integer :: chr,ip,jm,ngeno1,ngeno2,igprob,ig,ic
         integer :: test
         type(GENEALOGY_BASE) , pointer :: dg
         type(PHENOTYPE_BASE) , pointer :: dpa
         type(MAP_BASE) , pointer    :: map

         call log_mess("module incidence : init_analyse_linear",DEBUG_DEF)

         map => dataset%map
         dg => dataset%genea
         dpa => dataset%phenoAnimal

     end subroutine init_analyse_linear
!!***


!!****f* m_qtlmap_incidence/end_analyse_linear
!!  NAME
!!    end_analyse_linear
!!  DESCRIPTION
!!    close the linear mode
!!  NOTES
!!  SOURCE
     subroutine end_analyse_linear

     end subroutine end_analyse_linear
!!***


!!****f* m_qtlmap_incidence/init_incidence
!!  NAME
!!    init_incidence
!!  DESCRIPTION
!!    initilise the incidence structure (INCIDENCE_TYPE)
!!  INPUTS
!!    ic             : index of the trait
!!    nqtl           : number of hypothesis to test
!!    type_incidence : type of contingence matrix construction (INCIDENCE_TYPE_CLASSIC, INCIDENCE_TYPE_LD, INCIDENCE_TYPE_LDLA, INCIDENCE_TRAIT_COV)
!!    dataset        : use dataset (type DATASET_TYPE)
!!  OUTPUTS
!!    incidenceDesc  : the description of the contingence matrix
!!  NOTES
!!    see copy_dataset,init_dataset_1car
!!    if dataset are not created init_dataset_1car is called otherwise copy_dataset
!!  SOURCE
   subroutine init_incidence(dataset,spt,ic,nqtl,incidenceDesc,workstruct,datasetIncidence)
      type(QTLMAP_DATASET)                , intent(inout)     :: dataset !dec 2012 - dataset generalisé à l'ensemble des sources qtlmap
      type(PDD_BUILD)                     ,intent(inout)      :: spt
      integer                             , intent(in)        :: ic
      type(INCIDENCE_GEN_STRUCT)            ,intent(in)       :: workstruct
      integer                             , intent(in)        :: nqtl
      type(INCIDENCE_TYPE)                , intent(inout)     :: incidenceDesc
      type(DATASET_TYPE)                , intent(in),optional :: datasetIncidence
      integer :: npar
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      call log_mess("module incidence : init_incidence",DEBUG_DEF)

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      incidenceDesc%isBiAllellic=workstruct%biq

      incidenceDesc%nqtl=nqtl
      call set_ntnivmax(dataset,ic,nqtl,incidenceDesc%ntnivmax,incidenceDesc%nteffmax,incidenceDesc%ntlev,incidenceDesc%nbniv)

      if ( workstruct%effects%inter /= 0 ) then
         !! pour l instant que 2 qtl...
         if ( incidenceDesc%isBiAllellic ) then
          incidenceDesc%ntnivmax=incidenceDesc%ntnivmax+4
          if ( count(dpa%estime(ic,:)) > 0 ) incidenceDesc%ntnivmax=incidenceDesc%ntnivmax+16
          incidenceDesc%nteffmax=incidenceDesc%nteffmax+1
         else
          incidenceDesc%ntnivmax=incidenceDesc%ntnivmax+dg%np!*4+count(estime(ic,:))*16
          incidenceDesc%nteffmax=incidenceDesc%nteffmax+1
         end if

         allocate ( incidenceDesc%pp2(dg%nd))
         allocate ( incidenceDesc%pm2(dg%nd))

      end if

      if ( workstruct%effects%traits /= 0 ) then
         incidenceDesc%ntnivmax=incidenceDesc%ntnivmax+dpm%ncar
         incidenceDesc%nteffmax=incidenceDesc%nteffmax+1
      end if

      incidenceDesc%ic       = ic
      incidenceDesc%ntniv    = 0
      incidenceDesc%nbnivest = 0
      incidenceDesc%nteff    = 0

      call log_mess("********************* INIT *****************************",DEBUG_DEF)
      call log_mess("NTNIVMAX:"//str(incidenceDesc%ntnivmax),DEBUG_DEF)
      call log_mess("NTEFFMAX:"//str(incidenceDesc%nteffmax),DEBUG_DEF)
      call log_mess("NTLEVQTL / Reproductor:"//str(incidenceDesc%ntlev),DEBUG_DEF)
      call log_mess("NBLEV FIXED EFFECT:"//str(incidenceDesc%nbniv),DEBUG_DEF)
      call log_mess("*********************************************************",DEBUG_DEF)

      allocate(incidenceDesc%desc(incidenceDesc%nteffmax))
      allocate(incidenceDesc%vecsol(incidenceDesc%ntnivmax))
      incidenceDesc%vecsol=.false.
    !  allocate(incidenceDesc%corniv(ntnivmax))
      allocate(incidenceDesc%precis(incidenceDesc%ntnivmax))
      incidenceDesc%precis=0.d0
      if ( nqtl > 0) then
       allocate (incidenceDesc%ntniv_qtlsires(nqtl))
       allocate (incidenceDesc%ntniv_qtldams(nqtl))
       allocate ( incidenceDesc%pp(dg%nd))
       allocate ( incidenceDesc%pm(dg%nd))
       incidenceDesc%ntniv_qtlsires=0
       incidenceDesc%ntniv_qtldams=0
      end if
      allocate (incidenceDesc%corr_niv_nivb(incidenceDesc%ntnivmax))

      allocate (incidenceDesc%dataset)
      if ( present(datasetIncidence) ) then
        call copy_dataset(datasetIncidence,incidenceDesc%dataset)
      else
        call init_dataset_1car(dataset,spt,ic,nqtl,incidenceDesc%dataset)
      end if

      npar = dg%np+incidenceDesc%ntnivmax

      if ( associated (dpm%nmod) ) then
       if ( dpm%nmod(ic) > 1 ) npar = npar + dpm%nmod(ic)-1
      end if

      allocate (incidenceDesc%borni(npar))
      allocate (incidenceDesc%borns(npar))
      allocate (incidenceDesc%par(npar))

      allocate (incidenceDesc%fperemax(dataset%genea%np))
      allocate (incidenceDesc%fmeremax(dataset%genea%nm))

      incidenceDesc%borni(1:dg%np)=SIG_MIN
      incidenceDesc%borns(1:dg%np)=SIG_MAX
      incidenceDesc%borni(dg%np+1:)=DEFAULT_PARAM_MIN
      incidenceDesc%borns(dg%np+1:)=DEFAULT_PARAM_MAX

      if ( associated (dpm%nmod) ) then
       if (dpm%nmod(ic)>1) then
          incidenceDesc%borni(dg%np+incidenceDesc%ntnivmax+1:)=SIG_MIN
          incidenceDesc%borns(dg%np+incidenceDesc%ntnivmax+1:)=DEFAULT_PARAM_MAX
       end if
      end if

      incidenceDesc%par=0.d0

      allocate (incidenceDesc%nonull(incidenceDesc%dataset%nkd,incidenceDesc%ntnivmax))
      allocate (incidenceDesc%nbnonull(incidenceDesc%dataset%nkd))

      allocate (incidenceDesc%incidenceSireHaplo(dg%np,nqtl,2))
      incidenceDesc%incidenceSireHaplo=''
      allocate (incidenceDesc%incidenceSireHaploRace(dg%np,nqtl,2))
      incidenceDesc%incidenceSireHaploRace=''

    end subroutine init_incidence
!!***

!!****f* m_qtlmap_incidence/end_incidence
!!  NAME
!!    end_incidence
!!  DESCRIPTION
!!    release the incidence structure (INCIDENCE_TYPE)
!!  INPUTS
!!    incidenceDesc  : the description of the contingence matrix
!!  NOTES
!!  SOURCE
    subroutine end_incidence(incidenceDesc)
      type(INCIDENCE_TYPE)                     , intent(inout) :: incidenceDesc
      integer :: i,ip,jm,ifem

      call log_mess("module incidence : end_incidence",DEBUG_DEF)
      !release structure recursive of description of each effect
      do i=1,incidenceDesc%nteff
        if ( incidenceDesc%desc(i)%haveSubDesc) then
          deallocate (incidenceDesc%desc(i)%listSubDesc)
        end if
      end do
      !release all alocated array
      deallocate(incidenceDesc%desc)
      deallocate(incidenceDesc%vecsol)
      deallocate(incidenceDesc%precis)
      deallocate (incidenceDesc%corr_niv_nivb)

      if ( incidenceDesc%nqtl > 0) then
       deallocate (incidenceDesc%ntniv_qtlsires)
       deallocate (incidenceDesc%ntniv_qtldams)
       if (associated(incidenceDesc%pp)) deallocate (incidenceDesc%pp,incidenceDesc%pm)
      end if

      if ( associated(incidenceDesc%pp2) ) then
        deallocate (incidenceDesc%pp2,incidenceDesc%pm2)
      end if

     if (associated(incidenceDesc%dataset)) then
      call end_dataset(incidenceDesc%dataset)
      deallocate (incidenceDesc%dataset)
     end if

     if (associated(incidenceDesc%par) ) then
      deallocate (incidenceDesc%borni)
      deallocate (incidenceDesc%borns)
      deallocate (incidenceDesc%par)
     end if

     if (associated(incidenceDesc%fperemax) ) then
       deallocate (incidenceDesc%fperemax)
     end if
     if (associated(incidenceDesc%fmeremax) ) then
       deallocate(incidenceDesc%fmeremax)
     end if

     if ( associated (incidenceDesc%nonull) ) deallocate(incidenceDesc%nonull)
     if ( associated (incidenceDesc%nbnonull) ) deallocate(incidenceDesc%nbnonull)

     if ( associated (incidenceDesc%incidenceSireHaplo) ) deallocate (incidenceDesc%incidenceSireHaplo)
     if ( associated (incidenceDesc%incidenceSireHaploRace) )deallocate (incidenceDesc%incidenceSireHaploRace)

    end subroutine end_incidence
!!***

!!****f* m_qtlmap_incidence/init_workstruct
!!  NAME
!!    init_workstruct
!!  DESCRIPTION
!!    initialize the workstruct structure (INCIDENCE_GEN_STRUCT)
!!    The following information is given for each position on the chromosome :
!!    * variance family
!!    * residual correlation matrix (multi trait analysis case)
!!    * value of the likelihood
!!  INPUTS
!!   type_incidence  : type of contingence matrix construction (INCIDENCE_TYPE_CLASSIC, INCIDENCE_TYPE_LD, INCIDENCE_TYPE_LDLA, INCIDENCE_TRAIT_COV)
!!   type_resolution : type of resolution (LINEAR_HOMOSCEDASTIC, LINEAR_HETEROSCEDASTIC,LINEAR_MODANIM_HOMOSCEDASTIC, OPTIM_HETEROSCEDASTIC )
!!   maxHypothesisToTest    : number of hypothesis to test
!!   pc              : need to compute the confusion step
!!   pt              : need to compute the precision step
!!   multi           : boolean for multi trait analysis
!!  OUTPUTS
!!    workstruct  : buffer structure to save information about result among the chromosome
!!  NOTES
!!  SOURCE
   subroutine init_workstruct(dataset,ic,workstruct,effects,type_resolution,nqtl,maxHypothesisToTest,pc,pt,multi)
          type(QTLMAP_DATASET)                    ,intent(in)       :: dataset
          integer                                 ,intent(in)       :: ic
          type(INCIDENCE_GEN_STRUCT)              ,intent(inout)    :: workstruct
          type(type_effect_contingence)           ,intent(in)       :: effects
          integer                                 ,intent(in)       :: type_resolution,nqtl,maxHypothesisToTest
          logical                                 ,intent(in)       :: pc,pt,multi

          type(GENEALOGY_BASE) , pointer :: dg
          type(DATAMODEL_BASE) , pointer :: dpm

          call log_mess("module incidence : init_workstruct",DEBUG_DEF)

          dpm => dataset%phenoModel
          dg => dataset%genea

          workstruct%type_model=type_resolution
          workstruct%effects=effects

          workstruct%linear=.false.
          if (is_parameter( RESOLUTION_LINEAR,workstruct%type_model) ) then
            workstruct%linear=.true.
          end if

          workstruct%performConfusion=pc
          workstruct%performTestNuis=pt
          workstruct%nqtl=nqtl

          ! Traits as covariate effect. by default all traits define in the model
           workstruct%startTraitAsCov=1
           workstruct%endTraitAsCov=dpm%ncar

          if (multi) then
           allocate (workstruct%sigsquare(0:maxHypothesisToTest,dg%np,dpm%ncar))
           allocate (workstruct%rhoi(0:maxHypothesisToTest,dpm%ncar,dpm%ncar))
           workstruct%rhoi=0.d0
          else
           allocate (workstruct%sigsquare(0:maxHypothesisToTest,dg%np,1))
          end if
          allocate (workstruct%fnqtlsires(0:maxHypothesisToTest,dataset%genea%np))
          allocate (workstruct%fnqtldams(0:maxHypothesisToTest,dataset%genea%nm))

          workstruct%sigsquare=0.d0
          workstruct%fnqtlsires=0.d0
          workstruct%fnqtldams=0.d0

          allocate (workstruct%listnteff(nqtl))
          if (workstruct%effects%haplo > 0) then
             allocate (workstruct%listnteffhap(nqtl))
          end if

          workstruct%nteffinter = -1

          allocate (workstruct%nalert(maxHypothesisToTest))
          workstruct%nalert=0
          allocate (workstruct%alertQtl(maxHypothesisToTest,500))
          allocate (workstruct%corrmax(maxHypothesisToTest))


          if ( is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) then
            call build_M_modanim(dataset,ic,workstruct)
          end if


      end subroutine init_workstruct
!!***

!!****f* m_qtlmap_incidence/release_ws
!!  NAME
!!    release_ws
!!  DESCRIPTION
!!    release the workstruct structure (INCIDENCE_GEN_STRUCT)
!!  INPUTS
!!    workstruct  : buffer structure to save information about result among the chromosome
!!  NOTES
!!  SOURCE
     subroutine release_ws(workstruct)
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)  :: workstruct

        call log_mess("module incidence : release_ws",DEBUG_DEF)

         deallocate(workstruct%sigsquare)
         deallocate (workstruct%listnteff)
         deallocate (workstruct%fnqtlsires)
         deallocate (workstruct%fnqtldams)
         if (associated(workstruct%rhoi)) deallocate(workstruct%rhoi)
         if (associated(workstruct%listtestnuis)) deallocate(workstruct%listtestnuis)
         if (associated(workstruct%listnteffhap)) deallocate(workstruct%listnteffhap)
         deallocate (workstruct%nalert)
         deallocate (workstruct%alertQtl)
         deallocate (workstruct%corrmax)
         if (associated(workstruct%M)) deallocate(workstruct%M)

     end subroutine release_ws
!!***


  !constition des famille de plein et demi frere via des index
  !pas de notion d estimable...on comptabilise seulement les effectifs et
  ! on repere les lignes de la matrice d incidence

!!    * build information family without unkwnown trait (trait and data censured/ discrete trait)
!!      (algorithm dot not use presentc and estime arrays)
!!    build an subtype array : lSire of type SireFamily_type
!!      * variance and mean for each family (sig0, xmu0)
!!      * allocation of probabilities of transmission arrays (by sire family and full family)
!!  INPUTS
!!    ic      : index of the trait
!!    nqtl    : number of qtl (hypothesis)
!!  OUTPUTS
!!    dataset : the user datatset (DATASET_TYPE)
  subroutine init_dataset_1car(data,spt,ic,nqtl,dataset)
    type(QTLMAP_DATASET) ,target             ,intent(inout) :: data    ! modif inout car on point sur la structure ...bug revele avec ifort
    type(PDD_BUILD)      ,target             ,intent(inout) :: spt     ! meme chose  OFI 27/09/2013
    type(DATASET_TYPE)                      , intent(inout) :: dataset
    integer                               , intent(in)      :: ic,nqtl
    logical ,dimension(data%genea%nd) :: pres
    integer :: ip,jm,kd,ifem,kdd,s1,s2,chr,iq,kkd,effp
    real(kind=dp) :: somyp
    logical :: isContinue
    type(GENEALOGY_BASE) , pointer :: dg
    type(PHENOTYPE_BASE) , pointer :: dpa
    type(DATAMODEL_BASE) , pointer :: dpm

    call log_mess("module incidence : init_dataset_1car",DEBUG_DEF)

    dpm => data%phenoModel
    dg => data%genea
    dpa => data%phenoAnimal

    dataset%data => data
    dataset%spt  => spt

    allocate(dataset%lSires(dg%np))
    kdd=0

    do ip=1,dg%np
      allocate(dataset%lSires(ip)%full_sib(dg%nmp(ip+1)-dg%nmp(ip)))
      ifem=0
      do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
        ifem=ifem+1
        kd=dg%ndm(jm)+1
        isContinue = dpa%presentc(ic,kd)
        do while ( kd<=dg%ndm(jm+1) .and. .not. isContinue )
          kd=kd+1
          if ( kd <= dg%nd ) then
           isContinue = dpa%presentc(ic,kd)
          end if
        end do

        if (kd>dg%ndm(jm+1)) cycle
        kdd=kdd+1
        ! on a trouve le premier kd valide du pere ip dans la matrice d incidence
        dataset%lSires(ip)%full_sib(ifem)%firstKd=kdd
       ! print *,trim(pere(ip)),"-",trim(mere(jm)),ndm(jm)+1,ndm(jm+1),count(presentc(ic,kd+1:ndm(jm+1)))
       ! print *,kd+1,ndm(jm+1)

        kdd=kdd + count(dpa%presentc(ic,kd+1:dg%ndm(jm+1)))
        dataset%lSires(ip)%full_sib(ifem)%lastKd=kdd
        call log_mess("Full Sib ["//trim(dg%pere(ip))//"-"//trim(dg%mere(jm))//"] :"//&
        trim(str(dataset%lSires(ip)%full_sib(ifem)%firstKd))//"->"&
        //trim(str(dataset%lSires(ip)%full_sib(ifem)%lastKd)),VERBOSE_DEF)

        if ( nqtl >= 1) then
         s1=0
         do chr=1,data%map%nchr
          s1=max(s1,spt%ngenom(chr,jm+1)-spt%ngenom(chr,jm))
         end do

         s2=dataset%lSires(ip)%full_sib(ifem)%lastKd-dataset%lSires(ip)%full_sib(ifem)%firstKd+1

          allocate(dataset%lSires(ip)%full_sib(ifem)%ppt(nqtl,s1,s2))
          allocate(dataset%lSires(ip)%full_sib(ifem)%pmt(nqtl,s1,s2))
        end if
      end do
      ifem=0
      do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
        ifem=ifem+1
        if (dataset%lSires(ip)%full_sib(ifem)%firstKd<0) cycle
        dataset%lSires(ip)%half_sib%firstKd= dataset%lSires(ip)%full_sib(ifem)%firstKd
        exit
      end do
      ifem=dg%nmp(ip+1)-dg%nmp(ip)+1
      do jm=dg%nmp(ip+1),dg%nmp(ip)+1,-1
        ifem=ifem-1
        if (dataset%lSires(ip)%full_sib(ifem)%lastKd<0) cycle
        dataset%lSires(ip)%half_sib%lastKd=dataset%lSires(ip)%full_sib(ifem)%lastKd
        exit
      end do

      call log_mess("Half Sib ["//trim(dg%pere(ip))//"] :"//&
        trim(str(dataset%lSires(ip)%half_sib%firstKd))//"->"//trim(str(dataset%lSires(ip)%half_sib%lastKd)),VERBOSE_DEF)
      !s1=lSires(ip)%half_sib%lastKd-lSires(ip)%half_sib%firstKd+1
      if ( nqtl >= 1) then
        allocate(dataset%lSires(ip)%ppt(nqtl,dataset%lSires(ip)%half_sib%firstKd:dataset%lSires(ip)%half_sib%lastKd))
      end if
    end do ! ip

    isContinue = dpm%natureY(ic) /= 'i'
    !pour ne pas a tester la presence des caracteres on garde ce qui est seulement present
    if ( isContinue ) then
        allocate (dataset%Y(1,count(dpa%presentc(ic,:))))
        allocate (dataset%CD(1,count(dpa%presentc(ic,:))))
      else
        allocate (dataset%YDISCRETORD(1,count(dpa%presentc(ic,:))))
    end if

    dataset%nkd=0
    dataset%nkd_max_by_fam=0

    do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
         dataset%nkd_max_by_fam=max(dataset%nkd_max_by_fam,count(dpa%presentc(ic,dg%ndm(jm)+1:dg%ndm(jm+1))))
         do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
          if (dpa%presentc(ic,kd)) then
            dataset%nkd=dataset%nkd+1
            if ( isContinue ) then
              dataset%Y(1,dataset%nkd)=dpa%y(ic,kd)
              dataset%CD(1,dataset%nkd)=dpa%cd(ic,kd)
            else
              dataset%YDISCRETORD(1,dataset%nkd)=dpa%ydiscretord(ic,kd)
            end if
          end if
         end do !kd
        end do ! jm
    end do ! ip

    do ip=1,dg%np
        allocate(dataset%lSires(ip)%sig0(1))
        allocate(dataset%lSires(ip)%xmu0(1))
        effp=0
        somyp=0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if (dpa%presentc(ic,kd)) then
               effp=effp+1
               somyp=somyp+dpa%y(ic,kd)*dpa%cd(ic,kd)
            end if
          end do
        end do

        if (effp == 0.d0) then
           call stop_application('Father ['//trim(dg%pere(ip))//&
            '] has got no child with trait value ['//trim(dpm%carac(ic))//"]")
        end if
        !moyenne du caractere dans la famille du pere ip
        dataset%lSires(ip)%xmu0(1)=somyp/dble(effp)

        !ecart type du caractere dans la famille du pere ip
        dataset%lSires(ip)%sig0(1)=0.d0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if (dpa%presentc(ic,kd)) then
               dataset%lSires(ip)%sig0(1)=dataset%lSires(ip)%sig0(1)+(dpa%y(ic,kd)-&
                dataset%lSires(ip)%xmu0(1))*(dpa%y(ic,kd)-dataset%lSires(ip)%xmu0(1))
            end if
          end do !kd
        end do !jm
        dataset%lSires(ip)%sig0(1)=dsqrt(dataset%lSires(ip)%sig0(1)/dble(effp-1))
      end do !ip
  end subroutine init_dataset_1car
!!***

!********************************************************************************************************************************************************************
! Animal Model - February 2012

 subroutine init_parent_recod(dataset,id,v,idname,recod,idRecod,sizeRecod)
      type(QTLMAP_DATASET)       ,intent(in) :: dataset
      integer , intent(in) :: id
      integer , intent(in) :: v !1,2 ou 3 1=>prog,2=>size,3=>dam
      integer , intent(in) :: idRecod,sizeRecod
      character(len=LEN_DEF) , dimension(:),intent(in) :: idname
      integer , dimension(:,:) , intent(out)           :: recod

      integer :: idp,idm
      call log_mess("incidence module: init_parent_recod",DEBUG_DEF)

      if ( v == 1 ) then
            if ( dataset%genea%OldGenealogy(id,2) /= '0') then
             idp=1
             do while (idname(idp) /= dataset%genea%OldGenealogy(id,2) )
              idp=idp+1
              if ( idp == sizeRecod ) then
               call stop_application("build_genealogy_recod: does not find sire:"//trim(dataset%genea%OldGenealogy(id,2)))
              end if
             end do
             recod(idRecod,2)=idp
            end if

            if ( dataset%genea%OldGenealogy(id,3) /= '0' ) then
             idm=1
             do while(idname(idm) /= dataset%genea%OldGenealogy(id,3) )
              idm=idm+1
              if (idm == sizeRecod ) then
                call stop_application("build_genealogy_recod: does not find dam:"//trim(dataset%genea%OldGenealogy(id,3)))
              end if
             end do
             recod(idRecod,3)=idm
            end if
      end if

 end subroutine init_parent_recod

!!****f* m_qtlmap_incidence/build_genealogy_recod
!! NAME
!!   build_genealogy_recod
!! DESCRIPTION

subroutine build_genealogy_recod(dataset,idname,recod,sizeRecod)
   type(QTLMAP_DATASET)                  ,intent(in) :: dataset
   character(len=LEN_DEF) , dimension(:),intent(out) :: idname
   integer , dimension(:,:) , intent(out)            :: recod
   integer , intent(out)                             :: sizeRecod

   integer :: i,v,idp,idm,index,id
   logical :: new

   call log_mess("incidence module:build_genealogy_recod",DEBUG_DEF)
  ! print *," build_genealogy_recod ************* "
     ! recodage des individus
     recod = 0
     sizeRecod=0

     do index=2,4
      v=mod(index,4)
      if ( v == 0 ) v=1
      new=.true.
      do i=1,dataset%genea%OldGenealogySize
       if ( dataset%genea%OldGenealogy(i,v) /= '0') then
        if ( i > 1 ) then
         if ( dataset%genea%OldGenealogy(i-1,v) == dataset%genea%OldGenealogy(i,v) ) then
           recod(i,v) = recod(i-1,v)
           new=.false.
         else
           new = .true.
         end if
        end if
        !on regarde si on ne le trouve pas dans la genealogy qu'on a deja enregistre
        if ( new ) then
         do id=1,sizeRecod
           if ( trim(idname(id)) == dataset%genea%OldGenealogy(i,v) ) then
             new = .false.
             recod(i,v) = recod(id,1)
             call init_parent_recod(dataset,i,v,idname,recod,id,sizeRecod)
             exit
           end if
         end do
        end if

        if ( new ) then
           sizeRecod=sizeRecod+1
           idname(sizeRecod)=dataset%genea%OldGenealogy(i,v)
           recod(sizeRecod,1)=sizeRecod
           call init_parent_recod(dataset,i,v,idname,recod,sizeRecod,sizeRecod)
        end if
       end if
      end do
     end do

     call log_matrix("****** Old    Genealogy   ******** ",dataset%genea%OldGenealogy,1,&
      dataset%genea%OldGenealogySize,1,3,DEBUG_DEF)
     call log_matrix("****** Genealogy recoded  ******** ",recod,1,sizeRecod,1,3,DEBUG_DEF)

    ! print *," build_genealogy_recod ************* FIN "
end subroutine build_genealogy_recod

!!***

!!****f* m_qtlmap_incidence/build_A_modanim
!! NAME
!!   build_A_modanim
!! DESCRIPTION
!!   build the relationship and inbreeding matrix A (only progeny are used)
!!
!!   with the Henderson recursive method
!!   A = L.L' Aij = SUM k=1 ( Lik.Ljk )
!!
!!   i>j : Lij = (Aij - Sum Lik.Ljk) / Ljj
!!
!!  i==j : Lii = SQRT(Aii - SUM Lik^2)
!!
  subroutine build_A_modanim(dataset,ic,A,nkd)
     type(QTLMAP_DATASET)                 ,intent(in) :: dataset
     integer   ,intent(in)                            :: ic
     integer   ,intent(in)                            :: nkd
     real(kind=dp)  , dimension(nkd,nkd), intent(out) :: A


     ! the inbreeding matrix for all animal
     real(kind=dp), dimension(:,:),allocatable :: bigA
     integer , dimension(:,:) , allocatable :: recod
     integer :: i,j,k,p,q,ip,jm,kd,idSize,id,kkd,nBigA
     character(len=LEN_DEF) , dimension(:),allocatable :: idname
     integer , dimension(:),allocatable :: keepA
     logical :: new
     real(kind=dp) :: som
!     real(kind=dp)               :: t1,t2
     character , parameter       :: SIDE='R',UPLO='L',TRANSA='T',DIAG='N'
     real(kind=dp)               :: ALPHA=1.d0
     external :: DTRMM

     call log_mess("module incidence : build_A_modanim",DEBUG_DEF)

     allocate (recod(dataset%genea%OldGenealogySize*3,3))
     allocate (idname(dataset%genea%OldGenealogySize*3))

     call build_genealogy_recod(dataset,idname,recod,idSize)
     ! les peres et meme doivent etre dans la genealogy old
     allocate (bigA(idSize,idSize))

     ! fill the bigA with the OldGenealogy
     bigA=0.d0
     do id=1,idSize
        !print *,id,'/',idSize
        !fondateur
        if ( recod(id,2) == 0 .and. recod(id,3) == 0 ) then
          bigA(id,id)=1.d0
          cycle ! suivant....
        end if
       
        ! un des parent inconnue
        if ( recod (id,2) == 0 .or. recod(id,3) == 0 ) then
          if ( recod(id,2) /= 0 ) then
            p = recod(id,2)
          else
            p = recod(id,3)
          end if

          do j=1,p
            bigA(id,j)=0.5*bigA(p,j) 
          end do
      !    print *,id,j,bigA(id,j)
          som=0.d0
          do k=1,p
            som=bigA(id,k)*bigA(id,k)+som
          end do
          bigA(id,id)=sqrt(1.d0-som)
         cycle ! suivant...
        end if

       ! cas ou les fondateurs sont connues
       ! p < q
       if ( recod(id,2) < recod(id,3) ) then
        p=recod(id,2)
        q=recod(id,3)
       else
        p=recod(id,3)
        q=recod(id,2)
       end if
 !      print *,"-- p ",p,'  q:',q
       do j=1,p !ce n'est pas decrit comme ca normalement c est jusqu a p-1 
  !      print *,bigA(p,j),bigA(j,p)
        bigA(id,j) = 0.5*(bigA(p,j) + bigA(q,j))
   !     print *,id,j,bigA(id,j),' p>j>=1'
       end do
       do j=p+1,q
         bigA(id,j)=0.5*bigA(q,j)
    !     print *,id,j,bigA(id,j),' q>=j>p+1'
       end do
       som=0.d0
       do k=1,p
        som=som+bigA(p,k)*bigA(q,k)
       end do
       som=som*0.5
       do k=1,q
        som=som - bigA(id,k)*bigA(id,k)
       end do
       bigA(id,id)=sqrt(1.d0+som)
     !  print *,id,bigA(id,id)
     end do
!    t1=time()
    call log_mess("computation of A=L.L size="//trim(str(idSize)),INFO_DEF)
    !matrice triangulaire


!*  DTRMM  performs one of the matrix-matrix operations
!*
!*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!*
!*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!*
!*     op( A ) = A   or   op( A ) = A**T.
!************** DEBUG :  bigA = matmul(bigA,transpose(bigA))
    call DTRMM(SIDE,UPLO,TRANSA,DIAG,idSize,idSize,ALPHA,bigA,idSize,bigA,idSize)
!    t1=time()-t1
!    call log_mess("time of computation A=L.L :"//trim(str(t1)),DEBUG_DEF)

  !  call log_matrix_sparse(" The matrix relationship inbredding matrix A  :",bigA,1,idSize,1,idSize,INFO_DEF)

    kkd=0
    A=0.d0
    deallocate (recod)
    allocate (keepA(dataset%genea%nd))
    
     do ip=1,dataset%genea%np
       do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
         do kd=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
          if ( dataset%phenoAnimal%presentc(ic,kd) ) then
            !on cherche l'id le plus petit
            kkd=kkd+1
            id=1
            do while ( trim(idname(id)) /= trim(dataset%genea%animal(kd)) )
             !print *,"id:",id," name :",trim(idname(id))
             id=id+1
             if (id > idSize ) then
               call stop_application("can not find animal:"//trim(dataset%genea%animal(kd))//" (build_A_modelanim)")
             end if
            end do
            !print *,"ok"
            keepA(kkd)=id
            call log_mess("Keep line animal:"//trim(dataset%genea%animal(kd))//" id:"//trim(str(id)),DEBUG_DEF)
        !    print *,'keep index nkd=',nkd,"=> ligne de A:",id
          end if
         end do
       end do
     end do

     deallocate (idname)

     if ( kkd /= nkd ) then
       call stop_application("Internal error in incidence module (build_A_modanim). kkd:"//trim(str(kkd))//" nkd:"//trim(str(nkd)))
     end if

     do i=1,nkd
       do j=1,nkd
        A(i,j) = bigA(keepA(i),keepA(j))
       end do  
     end do 

     deallocate (keepA)
     deallocate (bigA)

  end subroutine build_A_modanim
!!***

!!****f* m_qtlmap_incidence/build_M_modanim
!! NAME
!!   build_M_modanim
!! DESCRIPTION
!!
!! Animal Model : Y = Xb + a + e
!!
!! VarA = Additive genetic variance / variance genetique additive
!!
!! a Gaussienne (0,A*VarA)
!!
!!
!!
!! M = (I - (I + lambda*A^-1)
!! Solution bHat = (X'.M.X)^-1 . X'.M.Y
!!
!! lambda = VarE / VarA
!!
!! H2 = VarA / ( VarA + VarE )
!!
!! Lambda is assumed fixed => lambda = (1 - H2) / H2
!!  lambda =
  subroutine build_M_modanim(data,ic,workstruct)
   type(QTLMAP_DATASET)          ,intent(in) :: data
   integer ,intent(in)                       :: ic
   type(INCIDENCE_GEN_STRUCT), intent(inout) :: workstruct
   !local
   integer       :: nkd,irank,i,ok
   real(kind=dp) :: ups
  ! type(GENEALOGY_BASE) , pointer :: dg
   type(PHENOTYPE_BASE) , pointer :: dpa
   type(DATAMODEL_BASE) , pointer :: dpm

   call log_mess("module incidence : build_M_modanim",DEBUG_DEF)

   dpm => data%phenoModel
  ! dg => dataset%genea
   dpa => data%phenoAnimal

   call log_mess("initializing animal model",INFO_DEF)
   nkd = count(dpa%presentc(ic,:))
   allocate (workstruct%M(nkd,nkd),stat=ok)
   if ( ok /= 0 ) then
    call stop_application("Not enough memory to allocate the relationship and inbreeding matrix A nd:"&
      //trim(str(nkd)))
   end if

   ! build A : the relationship and inbreeding matrix
   call log_mess("build the relationship and inbreeding matrix",INFO_DEF)
   call build_A_modanim(data,ic,workstruct%M,nkd)

   call log_matrix_sparse(" The matrix relationship inbredding matrix (only progeny) A  :",workstruct%M,1,nkd,1,nkd,DEBUG_DEF)

   call log_mess("inverse the relationship and inbreeding matrix",INFO_DEF)
   ! A^-1
   ups=1.d-15
   irank=0

   call ginv1(workstruct%M,nkd,size(workstruct%M,1),ups,irank)
   if ( irank /= nkd ) then
     call stop_application("Inversion of the relationship and inbreeding matrix failed rank:"//trim(str(irank)))
   end if

   call log_mess("H2 for the trait("//trim(dpm%carac(ic))//") = "//trim((str(dpm%h2(ic)))),INFO_DEF)
   ! Quantity : lambda.(A^-1)
   workstruct%M = workstruct%M * ((1.d0-dpm%h2(ic))/ dpm%h2(ic) )

   ! (I + lambda.A^-1)
   do i=1,nkd
    workstruct%M(i,i) = 1.d0+workstruct%M(i,i)
   end do
   ! (I + lambda.A^-1)^-1
   call ginv1(workstruct%M,nkd,size(workstruct%M,1),ups,irank)
   if ( irank /= nkd ) then
     call stop_application("Inversion of the relationship and inbreeding matrix failed rank:"//trim(str(irank)))
   end if

  ! - (I + lambda.A^-1)^-1
  workstruct%M = - workstruct%M
  ! I - (I + lambda.A^-1)^-1
  do i=1,nkd
    workstruct%M(i,i) = 1.d0-workstruct%M(i,i)
  end do

  end subroutine build_M_modanim
!!***


!!****f* m_qtlmap_incidence/end_dataset
!!  NAME
!!    end_dataset
!!  DESCRIPTION
!!    - release the variable dataset of type DATASET_TYPE
!!  OUTPUTS
!!    dataset : the user datatset (DATASET_TYPE)
!!
!!  NOTES
!!  SOURCE
   subroutine end_dataset(dataset)
      type(DATASET_TYPE)                      , intent(inout) :: dataset
      integer :: ip,jm,ifem

     call log_mess("module incidence : end_dataset",DEBUG_DEF)

     do ip=1,dataset%data%genea%np
      if ( associated(dataset%lSires(ip)%ppt) ) deallocate (dataset%lSires(ip)%ppt)
      if ( associated(dataset%lSires(ip)%sig0) ) deallocate (dataset%lSires(ip)%sig0)
      if ( associated(dataset%lSires(ip)%xmu0) ) deallocate (dataset%lSires(ip)%xmu0)
      ifem=0
      do jm=dataset%data%genea%nmp(ip)+1,dataset%data%genea%nmp(ip+1)
       ifem=ifem+1
       if (dataset%lSires(ip)%full_sib(ifem)%firstkd<0) cycle
       if ( associated(dataset%lSires(ip)%full_sib(ifem)%ppt) ) &
          deallocate(dataset%lSires(ip)%full_sib(ifem)%ppt)
       if ( associated(dataset%lSires(ip)%full_sib(ifem)%pmt) ) &
          deallocate(dataset%lSires(ip)%full_sib(ifem)%pmt)
      end do
      deallocate(dataset%lSires(ip)%full_sib)
    end do
    deallocate (dataset%lSires)

    if (associated (dataset%Y)) then
      deallocate(dataset%Y)
      deallocate(dataset%CD)
    else
      deallocate(dataset%YDISCRETORD)
    end if

    if (associated (dataset%rhoi) ) deallocate(dataset%rhoi)
    if (associated (dataset%cov) ) deallocate(dataset%cov)
    if (associated (dataset%sig) ) deallocate(dataset%sig)
    if (associated (dataset%xmu) ) deallocate(dataset%xmu)


   end subroutine end_dataset
!!***

!!****f* m_qtlmap_incidence/copy_dataset
!!  NAME
!!    copy_dataset
!!  DESCRIPTION
!!    Allocate new arrays and copy a variable of type DATASET_TYPE.
!!  INPUTS
!!    dataset : the original variable
!!  OUTPUTS
!!    out_dataset : the copy of dataset
!!
!!  NOTES
!!  SOURCE
    subroutine copy_dataset(dataset,out_dataset)
      type(DATASET_TYPE)                      , intent(in) :: dataset
      type(DATASET_TYPE)                      , intent(out) :: out_dataset
      integer :: ip,jm,ifem,s1,s2,s3,nfem

     call log_mess("module incidence : copy_dataset",DEBUG_DEF)

     out_dataset%data => dataset%data
     out_dataset%spt => dataset%spt

     out_dataset%nkd = dataset%nkd
     out_dataset%nkd_max_by_fam = dataset%nkd_max_by_fam

     allocate (out_dataset%lSires(dataset%data%genea%np))

     do ip=1,dataset%data%genea%np
      out_dataset%lSires(ip)%half_sib%firstKd = dataset%lSires(ip)%half_sib%firstKd
      out_dataset%lSires(ip)%half_sib%lastKd = dataset%lSires(ip)%half_sib%lastKd

      if ( associated(dataset%lSires(ip)%ppt) ) then
        s1=size(dataset%lSires(ip)%ppt,1)
        s2=size(dataset%lSires(ip)%ppt,2)
        allocate (out_dataset%lSires(ip)%ppt(s1,s2))
        out_dataset%lSires(ip)%ppt = dataset%lSires(ip)%ppt
      end if

      if ( associated(dataset%lSires(ip)%full_sib) ) then

        nfem = size(dataset%lSires(ip)%full_sib)
        allocate ( out_dataset%lSires(ip)%full_sib(nfem) )
        do ifem=1,nfem
           out_dataset%lSires(ip)%full_sib(ifem)%firstKd = dataset%lSires(ip)%full_sib(ifem)%firstKd
           out_dataset%lSires(ip)%full_sib(ifem)%lastKd = dataset%lSires(ip)%full_sib(ifem)%lastKd
           if ( associated(dataset%lSires(ip)%full_sib(ifem)%ppt) ) then
             s1=size(dataset%lSires(ip)%full_sib(ifem)%ppt,1)
             s2=size(dataset%lSires(ip)%full_sib(ifem)%ppt,2)
             s3=size(dataset%lSires(ip)%full_sib(ifem)%ppt,3)
             allocate ( out_dataset%lSires(ip)%full_sib(ifem)%ppt(s1,s2,s3) )
             out_dataset%lSires(ip)%full_sib(ifem)%ppt = dataset%lSires(ip)%full_sib(ifem)%ppt
           end if
            if ( associated(dataset%lSires(ip)%full_sib(ifem)%pmt) ) then
             s1=size(dataset%lSires(ip)%full_sib(ifem)%pmt,1)
             s2=size(dataset%lSires(ip)%full_sib(ifem)%pmt,2)
             s3=size(dataset%lSires(ip)%full_sib(ifem)%pmt,3)
             allocate ( out_dataset%lSires(ip)%full_sib(ifem)%pmt(s1,s2,s3) )
             out_dataset%lSires(ip)%full_sib(ifem)%pmt = dataset%lSires(ip)%full_sib(ifem)%pmt
           end if
        end do

      end if

      if ( associated(dataset%lSires(ip)%sig0) ) then
        allocate (out_dataset%lSires(ip)%sig0(size(dataset%lSires(ip)%sig0)))
        out_dataset%lSires(ip)%sig0 = dataset%lSires(ip)%sig0
      end if

       if ( associated(dataset%lSires(ip)%xmu0) ) then
        allocate (out_dataset%lSires(ip)%xmu0(size(dataset%lSires(ip)%xmu0)))
        out_dataset%lSires(ip)%xmu0 = dataset%lSires(ip)%xmu0
      end if

     end do


    if (associated (dataset%Y)) then
      s1=size(dataset%Y,1)
      s2=size(dataset%Y,2)
      allocate(out_dataset%Y(s1,s2))
      allocate(out_dataset%CD(s1,s2))
      out_dataset%Y = dataset%Y
      out_dataset%CD = dataset%CD
    else
      s1=size(dataset%YDISCRETORD,1)
      s2=size(dataset%YDISCRETORD,2)
      allocate(out_dataset%YDISCRETORD(s1,s2))
      out_dataset%YDISCRETORD = dataset%YDISCRETORD
    end if

    if (associated (dataset%rhoi) ) then
      s1=size(dataset%rhoi,1)
      s2=size(dataset%rhoi,2)
      allocate(out_dataset%rhoi(s1,s2))
      out_dataset%rhoi = dataset%rhoi
    end if

    if (associated (dataset%cov) ) then
      s1=size(dataset%cov,1)
      s2=size(dataset%cov,2)
      allocate(out_dataset%cov(s1,s2))
      out_dataset%cov = dataset%cov
    end if

    if (associated (dataset%sig) ) then
      s1=size(dataset%sig,1)
      allocate(out_dataset%sig(s1))
      out_dataset%sig = dataset%sig
    end if

    if (associated (dataset%xmu) ) then
      s1=size(dataset%xmu,1)
      allocate(out_dataset%xmu(s1))
      out_dataset%xmu = dataset%xmu
    end if

   end subroutine copy_dataset
!!***

!!****f* m_qtlmap_incidence/copy_incidence_desc
!!  NAME
!!    copy_incidence_desc
!!  DESCRIPTION
!!    Allocate partial arrays and copy a variable of type DATASET_TYPE.
!!    Static informations are referenced from the original variable in_inc.
!!
!!  INPUTS
!!    in_inc : the original variable
!!  OUTPUTS
!!    out_inc : the copy of incidence description
!!
!!  NOTES
!!  SOURCE
   subroutine copy_incidence_desc (in_inc,out_inc)
       use m_qtlmap_type_incidence_solution
       type(INCIDENCE_TYPE)  , intent(in) :: in_inc
       type(INCIDENCE_TYPE)  , intent(out) :: out_inc
       integer :: ip,ifem,jm,s1,s2,n,i
       type(GENEALOGY_BASE) , pointer :: dg

       call log_mess("module incidence : copy_incidence_desc",DEBUG_DEF)

       dg => in_inc%dataset%data%genea

       out_inc%ic  = in_inc%ic
       out_inc%nqtl     = in_inc%nqtl
       out_inc%ntniv    = in_inc%ntniv
       out_inc%nteff    = in_inc%nteff
       out_inc%eqtl_print = in_inc%eqtl_print
       out_inc%nbnivest = in_inc%nbnivest
       out_inc%ntlev    = in_inc%ntlev
       out_inc%nbniv    = in_inc%nbniv

       if ( associated (in_inc%desc) ) then
         allocate (out_inc%desc(size(in_inc%desc)))
         do i=1,size(in_inc%desc)
          call in_inc%desc(i)%copy(out_inc%desc(i))
         end do
       else
         out_inc%desc => null()
       end if

       allocate(out_inc%vecsol(size(in_inc%vecsol)))
       allocate(out_inc%precis(size(in_inc%precis)))

       if ( out_inc%nqtl > 0) then
        allocate (out_inc%ntniv_qtlsires(out_inc%nqtl))
        allocate (out_inc%ntniv_qtldams(out_inc%nqtl))
        allocate (out_inc%pp(dg%nd))
        allocate (out_inc%pm(dg%nd))

        out_inc%ntniv_qtlsires=0
        out_inc%ntniv_qtldams=0
       end if

        if ( associated(in_inc%pp2) ) then
          allocate (out_inc%pp2(dg%nd))
          allocate (out_inc%pm2(dg%nd))
        end if

       allocate (out_inc%corr_niv_nivb(size(in_inc%corr_niv_nivb)))

       out_inc%vecsol   = in_inc%vecsol
       out_inc%precis   = in_inc%precis
       out_inc%ntniv_qtlsires = in_inc%ntniv_qtlsires
       out_inc%ntniv_qtldams  = in_inc%ntniv_qtldams
       out_inc%corr_niv_nivb  = in_inc%corr_niv_nivb

       allocate (out_inc%dataset)
       !call init_dataset_1car(out_inc%ic,out_inc%nqtl,out_inc%dataset)
       allocate(out_inc%dataset%lSires(in_inc%dataset%data%genea%np))
       do ip=1,dg%np
         allocate(out_inc%dataset%lSires(ip)%full_sib(dg%nmp(ip+1)-dg%nmp(ip)))
         ifem=0
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            ifem=ifem+1
            out_inc%dataset%lSires(ip)%full_sib(ifem)%firstKd=in_inc%dataset%lSires(ip)%full_sib(ifem)%firstKd
            out_inc%dataset%lSires(ip)%full_sib(ifem)%lastKd=in_inc%dataset%lSires(ip)%full_sib(ifem)%lastKd

            if ( associated(in_inc%dataset%lSires(ip)%full_sib(ifem)%ppt) ) then
             n=size(in_inc%dataset%lSires(ip)%full_sib(ifem)%ppt,1)
             s1=size(in_inc%dataset%lSires(ip)%full_sib(ifem)%ppt,2)
             s2=size(in_inc%dataset%lSires(ip)%full_sib(ifem)%ppt,3)
             allocate(out_inc%dataset%lSires(ip)%full_sib(ifem)%ppt(n,s1,s2))
             allocate(out_inc%dataset%lSires(ip)%full_sib(ifem)%pmt(n,s1,s2))
             out_inc%dataset%lSires(ip)%full_sib(ifem)%ppt=in_inc%dataset%lSires(ip)%full_sib(ifem)%ppt
             out_inc%dataset%lSires(ip)%full_sib(ifem)%pmt=in_inc%dataset%lSires(ip)%full_sib(ifem)%pmt
            end if

         end do

         out_inc%dataset%lSires(ip)%half_sib%firstKd=in_inc%dataset%lSires(ip)%half_sib%firstKd
         out_inc%dataset%lSires(ip)%half_sib%lastKd=in_inc%dataset%lSires(ip)%half_sib%lastKd
         if ( associated(in_inc%dataset%lSires(ip)%ppt) ) then
          n=size(in_inc%dataset%lSires(ip)%ppt,1)
          s1=in_inc%dataset%lSires(ip)%half_sib%firstKd
          s2=in_inc%dataset%lSires(ip)%half_sib%lastKd
          allocate(out_inc%dataset%lSires(ip)%ppt(n,s1:s2))
          out_inc%dataset%lSires(ip)%ppt=in_inc%dataset%lSires(ip)%ppt
         end if

         out_inc%dataset%lSires(ip)%sig0=>in_inc%dataset%lSires(ip)%sig0
         out_inc%dataset%lSires(ip)%xmu0=>in_inc%dataset%lSires(ip)%xmu0

       end do

       out_inc%dataset%nkd = in_inc%dataset%nkd
       out_inc%dataset%nkd_max_by_fam = in_inc%dataset%nkd_max_by_fam
       out_inc%dataset%data  => in_inc%dataset%data ! attention on ne recopie pas le jeu de donnees qui est considerer comme statique et non modifiable
       out_inc%dataset%spt  => in_inc%dataset%spt

       out_inc%ntnivmax =in_inc%ntnivmax
       out_inc%nteffmax =in_inc%nteffmax

       out_inc%dataset%YDISCRETORD => in_inc%dataset%YDISCRETORD
       out_inc%dataset%Y           => in_inc%dataset%Y
       out_inc%dataset%CD          => in_inc%dataset%CD
       out_inc%dataset%cov         => in_inc%dataset%cov
       out_inc%dataset%rhoi        => in_inc%dataset%rhoi

       ! --- Fin copy dataset

       if (associated(in_inc%nonull) ) then
        allocate (out_inc%nonull(size(in_inc%nonull,1),size(in_inc%nonull,2)))
        allocate (out_inc%nbnonull(size(in_inc%nbnonull)))
        out_inc%nonull = in_inc%nonull
        out_inc%nbnonull = in_inc%nbnonull
       end if

       !les bornes sont sense ne pas change...
       out_inc%borni=>in_inc%borni
       out_inc%borns=>in_inc%borns

       allocate(out_inc%fperemax(dg%np))
       allocate(out_inc%fmeremax(dg%nm))

       allocate (out_inc%par(size(in_inc%par)))
       out_inc%par = in_inc%par

       if ( associated(in_inc%incidenceSireHaplo) ) then
        allocate(out_inc%incidenceSireHaplo(dg%np,in_inc%nqtl,2))
        out_inc%incidenceSireHaplo = in_inc%incidenceSireHaplo
       end if

       if ( associated(in_inc%incidenceSireHaploRace) ) then
        allocate(out_inc%incidenceSireHaploRace(dg%np,in_inc%nqtl,2))
        out_inc%incidenceSireHaploRace = in_inc%incidenceSireHaploRace
       end if

       out_inc%isBiAllellic = in_inc%isBiAllellic

   end subroutine copy_incidence_desc
!!***

!!****f* m_qtlmap_incidence/release_copy_incidence_desc
!!  NAME
!!    release_copy_incidence_desc
!!  DESCRIPTION
!!    release only arrays created during a copy (see copy_incidence_desc)
!!
!!  OUTPUTS
!!    copy : the copy of the incidence description
!!
!!  NOTES
!!  SOURCE
   subroutine release_copy_incidence_desc (copy)
     type(INCIDENCE_TYPE)  , intent(inout) :: copy
     integer :: i,ip,jm

     call log_mess("module incidence : release_copy_incidence_desc",DEBUG_DEF)

     if ( associated (copy%desc) ) then
         do i=1,size(copy%desc)
          call copy%desc(i)%release()
         end do
         deallocate (copy%desc)
     end if

     deallocate (copy%vecsol)
     deallocate (copy%precis)

     if ( copy%nqtl > 0 ) then
      deallocate (copy%ntniv_qtlsires)
      deallocate (copy%ntniv_qtldams)
      deallocate (copy%pp,copy%pm)
      if ( associated(copy%pp2)) then
        deallocate (copy%pp2,copy%pm2)
      end if
     end if

     deallocate (copy%corr_niv_nivb)
!     call end_dataset(copy%dataset)
     if ( associated(copy%dataset)) then

      do ip=1,copy%dataset%data%genea%np
        do jm=1,size(copy%dataset%lSires(ip)%full_sib)
         if ( associated(copy%dataset%lSires(ip)%full_sib(jm)%ppt) ) then
            deallocate ( copy%dataset%lSires(ip)%full_sib(jm)%ppt )
           deallocate ( copy%dataset%lSires(ip)%full_sib(jm)%pmt )
         end if
       end do
       deallocate ( copy%dataset%lSires(ip)%full_sib)
       if ( associated(copy%dataset%lSires(ip)%ppt) ) deallocate(copy%dataset%lSires(ip)%ppt )
      end do
      deallocate ( copy%dataset%lSires )
      deallocate ( copy%dataset )

      deallocate (copy%par)
      if (associated(copy%nonull) ) deallocate (copy%nonull)
      if (associated(copy%nbnonull) ) deallocate (copy%nbnonull)

      if (associated(copy%incidenceSireHaploRace) ) deallocate(copy%incidenceSireHaploRace)
      if (associated(copy%incidenceSireHaplo) ) deallocate(copy%incidenceSireHaplo)
     end if

     deallocate(copy%fperemax)
     deallocate(copy%fmeremax)

   end subroutine release_copy_incidence_desc
!!***

!!****f* m_qtlmap_incidence/init_position
!!  NAME
!!    init_position
!!  DESCRIPTION
!!    Initialize a position for a multi qtl hypothesis test (POSITION_LRT_INCIDENCE), e.g :
!!     * allocation of a LRT array value (comparaison with Hi-1, Hi-2,...)
!!     * allocation of the chromosome list (genome position of all qtl)
!!     * allocation of a discrete position list (index of absi array for each qtl position)
!!  INPUTS
!!    nqtl        : the number of qtl to test
!!  OUTPUTS
!!    position    : the variable poisiton
!!
!!  NOTES
!!  SOURCE
   subroutine init_position (dataset,hypothesis,nqtl,position)
     type(QTLMAP_DATASET)          ,intent(in)     :: dataset
     integer                       , intent(in)    :: hypothesis,nqtl
     type(POSITION_LRT_INCIDENCE)  , intent(inout) :: position

     integer :: i,j,n,chr

     call log_mess("module incidence : init_position",DEBUG_DEF)

     allocate (position%lrtSires(hypothesis,dataset%genea%np))
     allocate (position%lrtDams(hypothesis,dataset%genea%nm))
     allocate (position%listChr(nqtl))
     allocate (position%listN(nqtl))

     position%lrtSires=0.d0;
     position%lrtDams=0.d0;

     do i=1,nqtl
         n=0
         j=i
         do chr=1,dataset%map%nchr
           n=n+dataset%map%get_ilong(chr)
           if (i<=n) exit
           j=1
         end do

         if (chr>dataset%map%nchr) call stop_application("Not enough sampling point to detect QTLs")

         position%listChr(i)=1 ! attention si pas assez d echantillonage on doit passer au chromosome suivant....
         position%listN(i)=i

         do while (dataset%map%absi(1,position%listN(i))<dataset%map%posi(1,ceiling(dataset%params%longhap/2.d0)))
          position%listN(i)=position%listN(i)+1
         end do
      end do

   end subroutine init_position
!!***

!!****f* m_qtlmap_incidence/end_position
!!  NAME
!!    end_position
!!  DESCRIPTION
!!    release position of type POSITION_LRT_INCIDENCE
!!  OUTPUTS
!!    position    : the variable position
!!
!!  NOTES
!!  SOURCE
   subroutine end_position (position)
     type(POSITION_LRT_INCIDENCE)  , intent(inout) :: position

     call log_mess("module incidence : end_position",DEBUG_DEF)

     deallocate (position%lrtSires)
     deallocate (position%lrtDams)
     deallocate (position%listChr)
     deallocate (position%listN)

   end subroutine end_position
!!***

!!****f* m_qtlmap_incidence/add_general_mean
!!  NAME
!!    add_general_mean
!!  DESCRIPTION
!!    Add a general mean effect inside a contingence matrix with a description
!!
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_general_mean(xinc,incidenceDesc)
       type(INCIDENCE_TYPE)                                  , intent(inout) :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(inout) :: xinc
       integer :: ntniv

       call log_mess("incidence : ** Add General mean ** ",DEBUG_DEF)
       incidenceDesc%ntniv = incidenceDesc%ntniv+1

       ntniv = incidenceDesc%ntniv
       incidenceDesc%nteff = incidenceDesc%nteff+1
       incidenceDesc%desc(incidenceDesc%nteff)%name="General Mean"
       incidenceDesc%desc(incidenceDesc%nteff)%start=incidenceDesc%nteff
       incidenceDesc%desc(incidenceDesc%nteff)%end=incidenceDesc%nteff
       incidenceDesc%desc(incidenceDesc%nteff)%isVar = .false.
       xinc(:,ntniv) = 1
       !incidenceDesc%nivdir(:,incidenceDesc%nteff)   = 1
    end subroutine add_general_mean
!!***

!!****f* m_qtlmap_incidence/add_polygenic_effect
!!  NAME
!!    add_polygenic_effect
!!  DESCRIPTION
!!    Add a polygenic effect inside a contingence matrix with a description
!!    * An effect is adding for each Sire family
!!    * An effect is adding for each full sib family (depends ndmin variable)
!!
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_polygenic_effect(xinc,incidenceDesc)
       type(INCIDENCE_TYPE)                    , intent(inout) :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(inout) :: xinc
       integer :: ip,jm,kd,nt,j,fem,nteff,nkd
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       call log_mess("incidence : ** Add Polygenic effect to estimate **",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
       dpa => incidenceDesc%dataset%data%phenoAnimal

       nt=incidenceDesc%ntniv
       nteff = incidenceDesc%nteff

       incidenceDesc%desc(incidenceDesc%nteff+1)%name="Sire polygenic effects"
       incidenceDesc%desc(incidenceDesc%nteff+1)%start=nt+1
       incidenceDesc%desc(incidenceDesc%nteff+1)%end=nt+dg%np
       incidenceDesc%desc(incidenceDesc%nteff+1)%haveSubDesc=.true.
       incidenceDesc%desc(incidenceDesc%nteff+1)%isVar = .false.
       allocate (incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(dg%np))
       do ip=1,dg%np
         incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(ip)%name="Sire "//trim(dg%pere(ip))
         incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(ip)%start=nt+ip
         incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(ip)%end=nt+ip
       end do
       incidenceDesc%nteff = incidenceDesc%nteff+1
       incidenceDesc%ntniv = incidenceDesc%ntniv+dg%np

       if ( count(dpa%estime(incidenceDesc%ic,:))>0) then
         incidenceDesc%desc(incidenceDesc%nteff+1)%name="Dam polygenic effects"
         incidenceDesc%desc(incidenceDesc%nteff+1)%start=incidenceDesc%desc(incidenceDesc%nteff)%end+1
         incidenceDesc%desc(incidenceDesc%nteff+1)%end=nt+dg%np+count(dpa%estime(incidenceDesc%ic,:))!nbfemLoc(ic)
         incidenceDesc%desc(incidenceDesc%nteff+1)%haveSubDesc=.true.
         incidenceDesc%desc(incidenceDesc%nteff+1)%isVar = .false.
         allocate (incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(count(dpa%estime(incidenceDesc%ic,:))))
         fem=0
         do ip=1,dg%np
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if (dpa%estime(incidenceDesc%ic,jm) ) then
             fem=fem+1
             incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(fem)%name="Dam "//trim(dg%mere(jm))//&
              " [Sire "//trim(dg%pere(ip))//"]"
             incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(fem)%start=nt+dg%np+fem
             incidenceDesc%desc(incidenceDesc%nteff+1)%listSubDesc(fem)%end=nt+dg%np+fem
           end if
         end do
        end do
         incidenceDesc%nteff = incidenceDesc%nteff+1
         incidenceDesc%ntniv = incidenceDesc%ntniv+count(dpa%estime(incidenceDesc%ic,:))!nbfemLoc(ic)
       else

       end if
       j=0

       !initialisation
       xinc(:,nt+1:nt+dg%np+count(dpa%estime(incidenceDesc%ic,:dg%nm))) = 0

       nkd=0
       do ip=1,dg%np
         xinc(:,nt+ip) = 0
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
             if ( dpa%presentc(incidenceDesc%ic,kd) ) then
              nkd =nkd +1
              xinc(nkd,nt+ip) = 1
         !    incidenceDesc%nivdir(kd,nteff+1)= nt+ip
              if (dpa%estime(incidenceDesc%ic,jm)) then
               fem=count(dpa%estime(incidenceDesc%ic,:jm)) ! on compte le nombre estimable dans la famille de pere ip
               xinc(nkd,nt+dg%np+fem) = 1
          !    incidenceDesc%nivdir(kd,nteff+2)= nt+np+fem
             end if
             end if
           end do
         end do
       end do
    end subroutine add_polygenic_effect
!!***

!!****f* m_qtlmap_incidence/add_qtleffect
!!  NAME
!!    add_qtleffect
!!  DESCRIPTION
!!    Add a qtl effect inside a contingence matrix with a description
!!    * An effect is adding for each Sire family
!!    * An effect is adding for each dam with enough progenies (ndmin)
!!
!!  INPUTS
!!    nteff       : index of the effect
!!    numqtl      : index of the qtl effect
!!    chr         : chromosome localisation
!!    n           : the position
!!    mint        : index of the interaction to remove (otherwise mint <= 0)
!!    linear      : use the likely probability genotype dam or not
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_qtleffect(numqtl,chr,n,xinc,incidenceDesc,mint,linear)
       integer                                 ,intent(in)     :: numqtl,chr,n
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc
       integer                                 ,intent(in)     :: mint ! index of interaction to remove from incidence
       logical                                 ,intent(in)              :: linear

       call log_mess("incidence : ** Add QTL effect to estimate **",DEBUG_DEF)

       if ( incidenceDesc%isBiAllellic ) then
         call add_qtleffect_biallelic(numqtl,chr,n,xinc,incidenceDesc,mint,linear)
       else
         call add_qtleffect_sire(numqtl,chr,n,xinc,incidenceDesc,mint,linear)
       end if

    end subroutine add_qtleffect
!!***

!!****f* m_qtlmap_incidence/add_qtleffect_sire
!!  NAME
!!    add_qtleffect_sire
!!  DESCRIPTION
!!    Add a qtl effect inside a contingence matrix with a description
!!    * An effect is adding for each Sire family
!!    * An effect is adding for each dam with enough progenies (ndmin)
!!
!!  INPUTS
!!    nteff       : index of the effect
!!    numqtl      : index of the qtl effect
!!    chr         : chromosome localisation
!!    n           : the position
!!    mint        : index of the interaction to remove (otherwise mint <= 0)
!!    linear      : use the likely probability genotype dam or not
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_qtleffect_sire(numqtl,chr,n,xinc,incidenceDesc,mint,linear)
       integer                                 ,intent(in)     :: numqtl,chr,n
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc
       integer                                 ,intent(in)     :: mint ! index of interaction to remove from incidence
       logical                                 ,intent(in)              :: linear

       integer :: ntniv,nteffstart,ip,fem,jm,l,s,u,nbint,jef,ceff,scalei
       logical ,dimension(incidenceDesc%dataset%data%genea%nm) :: femInit
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       call log_mess("incidence : ** Add Qtl Sire effect to estimate **",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
       dpa => incidenceDesc%dataset%data%phenoAnimal
       dpm => incidenceDesc%dataset%data%phenoModel

       nteffstart = incidenceDesc%nteff
       ceff= incidenceDesc%nteff + 1
       if ( .not. incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL) ) then
        incidenceDesc%nteff = incidenceDesc%nteff+1
       end if

       if ( count(dpa%estime(incidenceDesc%ic,:)) > 0 ) incidenceDesc%nteff = incidenceDesc%nteff+1

       nbint=dpm%listModelTrait(incidenceDesc%ic)%nbint(numqtl) !modele(incidenceDesc%ic,3)
       incidenceDesc%ntlev=1

       do jef=1,nbint
          if ( mint == jef) cycle ! this interaction is not adding
          incidenceDesc%ntlev=incidenceDesc%ntlev*&
           dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(numqtl,jef))
       end do

       scalei=1

       ! pour hetero chromosome femelle, on peut ne pas estimer les QTL chez les peres
       if ( .not. incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL) ) then
        incidenceDesc%desc(nteffstart+scalei)%isVar = .true.
        incidenceDesc%desc(nteffstart+scalei)%name="Sire QTL effects ["//trim(str(numqtl))//"]"
        incidenceDesc%desc(nteffstart+scalei)%start=incidenceDesc%ntniv+1
        incidenceDesc%ntniv_qtlsires(numqtl)=incidenceDesc%desc(nteffstart+scalei)%start
        incidenceDesc%desc(nteffstart+scalei)%end=incidenceDesc%ntniv+incidenceDesc%ntlev*dg%np

        incidenceDesc%desc(nteffstart+scalei)%haveSubDesc=.true.
        allocate (incidenceDesc%desc(nteffstart+scalei)%listSubDesc(dg%np*incidenceDesc%ntlev))

        u=0
        do ip=1,dg%np
         l=0
         do s=incidenceDesc%ntniv+(ip-1)*incidenceDesc%ntlev+1,incidenceDesc%ntniv+ip*incidenceDesc%ntlev
          l=l+1
          u=u+1

          incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%name="Sire "//trim(dg%pere(ip))//" level "//trim(str(l))

          incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%start=s
          incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%end=s
         end do
        end do
        incidenceDesc%ntniv = incidenceDesc%desc(nteffstart+scalei)%end
        scalei=scalei+1
       end if


       if ( count(dpa%estime(incidenceDesc%ic,:)) > 0 ) then
         incidenceDesc%desc(nteffstart+scalei)%isVar = .true.
         incidenceDesc%desc(nteffstart+scalei)%name="Dam QTL effects ["//trim(str(numqtl))//"]"
         incidenceDesc%desc(nteffstart+scalei)%start=incidenceDesc%desc(nteffstart+scalei-1)%end+1
         incidenceDesc%ntniv_qtldams(numqtl)=incidenceDesc%desc(nteffstart+scalei)%start
         femInit=.false.
         l=0
         do jm=1,dg%nm
           if (dpa%estime(incidenceDesc%ic,jm) ) then
             fem=dpa%iam(incidenceDesc%ic,dg%repfem(jm))
             if (femInit(fem)) then
                cycle
             end if
             femInit(fem)=.true.
             l=l+1
           end if
         end do

         allocate (incidenceDesc%desc(nteffstart+scalei)%listSubDesc(l*incidenceDesc%ntlev))
         incidenceDesc%desc(nteffstart+scalei)%end=incidenceDesc%desc(nteffstart+scalei)%start+incidenceDesc%ntlev*l - 1
         incidenceDesc%desc(nteffstart+scalei)%haveSubDesc=.true.
         !fem=0
         femInit=.false.
         u=0
         do jm=1,dg%nm
           if (dpa%estime(incidenceDesc%ic,jm) ) then
             fem=dpa%iam(incidenceDesc%ic,dg%repfem(jm))
             if (femInit(fem)) then
                cycle
             end if
             femInit(fem)=.true.
             l=0
             do s=incidenceDesc%desc(nteffstart+scalei)%start+(fem-1)*incidenceDesc%ntlev,&
             incidenceDesc%desc(nteffstart+scalei)%start-1+fem*incidenceDesc%ntlev
              l=l+1
              u=u+1
              !fem=fem+1
              incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%name="Dam "//trim(dg%mere(jm))//" Level "//trim(str(l))
              incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%start=s
              incidenceDesc%desc(nteffstart+scalei)%listSubDesc(u)%end=s
             end do
           end if
         end do
         incidenceDesc%ntniv = incidenceDesc%desc(nteffstart+scalei)%end
       end if

       call change_qtleffect(ceff,numqtl,chr,n,xinc,incidenceDesc,mint,linear)

    end subroutine add_qtleffect_sire
!!***

!!****f* m_qtlmap_incidence/add_qtleffect_biallelic
!!  NAME
!!    add_qtleffect_biallelic
!!  DESCRIPTION
!!    Add a qtl bi-allelic effect inside a contingence matrix with a description
!!
!!  INPUTS
!!    nteff       : index of the effect
!!    numqtl      : index of the qtl effect
!!    chr         : chromosome localisation
!!    n           : the position
!!    mint        : index of the interaction to remove (otherwise mint <= 0)
!!    linear      : use the likely probability genotype dam or not
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_qtleffect_biallelic(numqtl,chr,n,xinc,incidenceDesc,mint,linear)
       integer                                 ,intent(in)     :: numqtl,chr,n
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc
       integer                                 ,intent(in)     :: mint ! index of interaction to remove from incidence
       logical                                 ,intent(in)              :: linear

       integer :: ntniv,nteffstart,ip,fem,jm,l,s,u,nbint,jef,ceff
       logical ,dimension(incidenceDesc%dataset%data%genea%nm) :: femInit
      ! type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       call log_mess("incidence : ** Add Qtl Bi-allelic effect to estimate **",DEBUG_DEF)
     !  dg => incidenceDesc%dataset%data%genea
       dpa => incidenceDesc%dataset%data%phenoAnimal
       dpm => incidenceDesc%dataset%data%phenoModel

       nteffstart = incidenceDesc%nteff
       incidenceDesc%nteff = incidenceDesc%nteff+1
       ceff=incidenceDesc%nteff
       nbint=dpm%listModelTrait(incidenceDesc%ic)%nbint(numqtl) !modele(incidenceDesc%ic,3)
       incidenceDesc%ntlev=1

       do jef=1,nbint
          if ( mint == jef) cycle ! this interaction is not adding
          incidenceDesc%ntlev=incidenceDesc%ntlev*&
           dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(numqtl,jef))
       end do

       incidenceDesc%desc(nteffstart+1)%isVar = .true.
       incidenceDesc%desc(nteffstart+1)%name="QTL effects ["//trim(str(numqtl))//"]"
       incidenceDesc%desc(nteffstart+1)%start=incidenceDesc%ntniv+1
       incidenceDesc%ntniv_qtlsires(numqtl)=incidenceDesc%desc(nteffstart+1)%start

       if ( count(dpa%estime(incidenceDesc%ic,:)) > 0 ) then
        incidenceDesc%ntniv_qtldams(numqtl)=incidenceDesc%desc(nteffstart+1)%start
       end if

       incidenceDesc%desc(nteffstart+1)%end=incidenceDesc%ntniv+incidenceDesc%ntlev

       incidenceDesc%desc(nteffstart+1)%haveSubDesc=.true.
       allocate (incidenceDesc%desc(nteffstart+1)%listSubDesc(incidenceDesc%ntlev))

       l=0
       do s=incidenceDesc%desc(nteffstart+1)%start,incidenceDesc%desc(nteffstart+1)%end
        l=l+1
        incidenceDesc%desc(nteffstart+1)%listSubDesc(l)%name="level "//trim(str(l))
        incidenceDesc%desc(nteffstart+1)%listSubDesc(l)%start=s
        incidenceDesc%desc(nteffstart+1)%listSubDesc(l)%end=s
       end do

       incidenceDesc%ntniv = incidenceDesc%desc(nteffstart+1)%end

       call change_qtleffect(ceff,numqtl,chr,n,xinc,incidenceDesc,mint,linear)

    end subroutine add_qtleffect_biallelic
!!***

!!****f* m_qtlmap_incidence/change_name_qtleffect
!!  NAME
!!    change_name_qtleffect
!!  DESCRIPTION
!!    change the name of the qtl effect
!!  NOTES
!!  SOURCE
    subroutine change_name_qtleffect(nteff,numqtl,chr,n,incidenceDesc)
       integer                                 ,intent(in)     :: nteff,numqtl,chr,n
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc

       integer :: ip,u,s,l,jhr
       character(LEN=100)  :: nhr(2)
       type(HAPLOTYPE_POSITION_BUILD)  :: shp
       type(GENEALOGY_BASE) , pointer  :: dg

       call log_mess("module incidence : change_name_qtleffect",DEBUG_DEF)

       if ( incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL) ) then
          return
       end if

       if ( incidenceDesc%isBiAllellic ) then
          return
       end if

       dg => incidenceDesc%dataset%data%genea
       call shp%set(incidenceDesc%dataset%data,incidenceDesc%dataset%spt)
       call shp%set_haplo_for_ldla(chr,incidenceDesc%dataset%data%map%absi(chr,n),n,.true.,.false.)
       u=0
       do ip=1,dg%np
           l=0
           do s=incidenceDesc%ntniv+(ip-1)*incidenceDesc%ntlev+1,incidenceDesc%ntniv+ip*incidenceDesc%ntlev
               l=l+1
               u=u+1
               nhr= ' '
               do jhr=1,2
                   if(shp%num_haplo_pere(ip,jhr,1) /= 0) nhr(jhr) =shp%name_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
               end do

               incidenceDesc%desc(nteff)%listSubDesc(u)%name="Sire "//trim(dg%pere(ip))//" l"//trim(str(l))&
                   //'['//trim(nhr(1))//'/'//trim(nhr(2))//']'

           end do
       end do
       call shp%free()

   end subroutine change_name_qtleffect
!!***


!!****f* m_qtlmap_incidence/add_haplotype_effect
!!  NAME
!!    add_haplotype_effect
!!  DESCRIPTION
!!    Add a haplotype effect inside a contingence matrix with a description
!!    An effect is adding by haplotype finded in the population at the specific position
!!    The number of haplotype is growing depending hdam optional variable.
!!
!!  INPUTS
!!    nteff       : index of the effect
!!    numhap      : index of the haplotype effect
!!    chr         : chromosome localisation
!!    n           : the position
!!    hdam        :
!!    start       : using in a inialize context
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_haplotype_effect(nteff,numhap,chr,n,xinc,incidenceDesc,shp,hdam,start)
       integer                                 ,intent(in)     :: nteff,numhap,chr,n
       type(HAPLOTYPE_POSITION_BUILD)         ,intent(inout)  :: shp
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc
       logical                                 ,intent(in)     :: hdam,start

       !local
       integer :: ip,jm,kd,ntniv,nkd,j_gam,i_haplo,j,i_gam,geno,nteffstart,jhr
       real (kind=dp) ,dimension(:,:),allocatable        :: pp_ldla,pm_ldla,pb_haplo

       real (kind=dp)  :: su
       integer :: nbHaploIncidence
       logical :: disableSire

       type(GENEALOGY_BASE) , pointer :: dg

       call log_mess("incidence : ** Add haplotype effect to estimate **",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea

       allocate (pb_haplo(dg%nd,incidenceDesc%dataset%data%params%NB_HAPLO_PRIOR))
       allocate (pp_ldla(dg%nd,2))
       allocate (pm_ldla(dg%nd,2))

       disableSire = incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL)

       ntniv = incidenceDesc%desc(nteff-1)%end
       nkd=0

       call shp%set_haplo_for_ldla(chr,incidenceDesc%dataset%data%map%absi(chr,n),n,.true.,hdam)
       !
       ! Attention : le nombre de niveau doit etre fixe. normalement il ne l est pas car le nombre d haplotype change a la position (nb_haplo_reduit)
       ! Pour le fonctionnement de incidence il est necessaire de fixer un maximum (2**longhap) * nb_races ( le *2 pour les races => 1 meme haplotype est evaluer dans chaque race...)
       nbHaploIncidence=incidenceDesc%dataset%data%geneaRace%nb_races*(2**incidenceDesc%dataset%data%params%longhap)

       if (start) then
         nteffstart = incidenceDesc%nteff
         incidenceDesc%nteff = incidenceDesc%nteff+1
 !        incidenceDesc%ntniv = ntniv + nb_haplo_reduit!ntniv + 2**longhap
         incidenceDesc%ntniv = ntniv + nbHaploIncidence
       else
         nteffstart = nteff - 1
       end if

       incidenceDesc%desc(nteffstart+1)%isVar = .true.
       incidenceDesc%desc(nteffstart+1)%name="Haplotypes effects["//trim(str(numhap))//"]"
       incidenceDesc%desc(nteffstart+1)%start=ntniv+1
!       incidenceDesc%desc(nteffstart+1)%end=ntniv+nb_haplo_reduit ! 2**longhap ! nb_haplo_reduit
       incidenceDesc%desc(nteffstart+1)%end=ntniv+nbHaploIncidence
       incidenceDesc%desc(nteffstart+1)%haveSubDesc=.true.
       if ( associated(incidenceDesc%desc(nteffstart+1)%listSubDesc) ) &
         deallocate (incidenceDesc%desc(nteffstart+1)%listSubDesc)

       allocate (incidenceDesc%desc(nteffstart+1)%listSubDesc(shp%nb_haplo_reduit))

       do j=1,shp%nb_haplo_reduit
        ! if (j<=nb_haplo_reduit) then
          incidenceDesc%desc(nteffstart+1)%listSubDesc(j)%name=trim(shp%name_haplo_reduit(j))//&
          "(race="//trim(shp%race_haplo_reduit(j)) //&
           ", freq="//trim(str(shp%pb_haplo_reduit(j)))//")"
        ! else
         ! incidenceDesc%desc(nteffstart+1)%listSubDesc(j)%name="none("//trim(str(j))//")"
        ! end if
         incidenceDesc%desc(nteffstart+1)%listSubDesc(j)%start=incidenceDesc%desc(nteffstart+1)%start+j-1
         incidenceDesc%desc(nteffstart+1)%listSubDesc(j)%end=incidenceDesc%desc(nteffstart+1)%start+j-1
       end do

       xinc(:,ntniv+1:ntniv+nbHaploIncidence) = 0.d0
  !     xinc(:,ntniv+1:ntniv+nb_haplo_reduit) = 0.d0
       do ip=1,dg%np
         !ajout traitement affichage des haplotypes et origines pour chaque pere
         do jhr=1,2
          incidenceDesc%incidenceSireHaplo(ip,numhap,jhr) = shp%name_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
          incidenceDesc%incidenceSireHaploRace(ip,numhap,jhr) = shp%race_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
         end do

         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          !if ( linear ) then
            call getprobld_lin(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr,n,incidenceDesc%ic,ip,jm,&
             dg%ndm(jm)+1,dg%ndm(jm+1),pp_ldla,pm_ldla)
!           else
!            call getglobalprob(chr,n,incidenceDesc%ic,ip,jm,ndm(jm)+1,ndm(jm+1),pp_ldla,pm_ldla)
           !end if
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if ( incidenceDesc%dataset%data%phenoAnimal%presentc(incidenceDesc%ic,kd) ) then
            nkd=nkd+1

            do i_haplo=1,shp%nb_haplo_reduit
             pb_haplo(kd,i_haplo)=0.d0
           end do

! haplotype du pere
           !if(hsire)then
             do j=1,2
               su = dble(sum(shp%pb_haplo_reduit(shp%num_haplo_pere(ip,j,1:shp%nb_gam_pere(ip,j)))))
               if ( su == 0 ) cycle
               do i_gam=1,shp%nb_gam_pere(ip,j)
                 pb_haplo(kd,shp%num_haplo_pere(ip,j,i_gam))= pb_haplo(kd,shp%num_haplo_pere(ip,j,i_gam)) &
                                                       +pp_ldla(kd,j) * ( shp%pb_haplo_reduit(shp%num_haplo_pere(ip,j,i_gam)) / su )
                end do !i_gam
              end do !j



           !end if ! opt_model
! haplotype de la mere si phase connue / on force si analyse hetero-chromosome
           if(hdam .or. disableSire)then
            if(incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,jm)) then
               do geno=incidenceDesc%dataset%spt%ngenom(chr,jm)+1,incidenceDesc%dataset%spt%ngenom(chr,jm+1)
                do j=1,2
                 su = dble(sum(shp%pb_haplo_reduit(shp%num_haplo_mere(geno,j,1:shp%nb_gam_mere(geno,j)))))
                 if ( su == 0 ) cycle
                 do i_gam=1,shp%nb_gam_mere(geno,j)
                   pb_haplo(kd,shp%num_haplo_mere(geno,j,i_gam))  = pb_haplo(kd,shp%num_haplo_mere(geno,j,i_gam)) &
                                                 +pm_ldla(kd,j) * ( shp%pb_haplo_reduit(shp%num_haplo_mere(geno,j,i_gam)) / su )
                 end do !i_gam
                end do !j
               end do
! haplotype de la mere si phase inconnue
           else
!             do j_gam=1,nb_gam(kd)
!               su = dble(sum(pb_haplo_reduit(num_haplo_desc(kd,j_gam,1:nb_gam_desc(kd,j_gam)))))
!               if ( su == 0 ) cycle
!               do i_gam=1,nb_gam_desc(kd,j_gam)
!                   pb_haplo(kd,num_haplo_desc(kd,j_gam,i_gam))  = pb_haplo(kd,num_haplo_desc(kd,j_gam,i_gam)) &
!                                                 +prob_gam(kd,j_gam) * (pb_haplo_reduit(num_haplo_desc(kd,j_gam,i_gam)) /su)
 !              end do !i_gam
 !            end do !j_gam

 !           write(111,*)'n,kd,nb_gam_desc(kd)',n,kd,nb_gam_desc(kd)
 !           do i_gam=1,nb_gam_desc(kd)
 !              pb_haplo(kd,num_haplo_desc(kd,i_gam))=pb_haplo(kd,num_haplo_desc(kd,i_gam)) + &
 !                       pb_haplo_desc(kd,i_gam)
              do i_gam=1,shp%nb_haplo_reduit
               pb_haplo(kd,i_gam)=pb_haplo(kd,i_gam) + shp%pb_haplo_desc(kd,i_gam)
            end do  !i_gam


           end if !estime
           end if !opt_model (hdam)
!
! ligne d'incidence
!
           do i_haplo=1,shp%nb_haplo_reduit
             xinc(nkd,ntniv+i_haplo)=pb_haplo(kd,i_haplo)
            end do
          end if !presentc
         end do
        end do
       end do

       !ntniv=ntniv+nb_haplo_reduit
      ! nteff=nteff+nb_haplo_reduit

       deallocate (pp_ldla,pm_ldla,pb_haplo)
      ! deallocate (pm_ldla)

    end subroutine add_haplotype_effect
!!***

!!****f* m_qtlmap_incidence/change_qtleffect
!!  NAME
!!    change_qtleffect
!!  DESCRIPTION
!!    Change the contingence matrix (probabilities of transmission information) during the chromosome scan .
!!
!!  INPUTS
!!    nteff       : index of the effect
!!    iq          : index of the qtl effect
!!    chr         : chromosome localisation
!!    n           : the position
!!    mint        : index of the interaction to remove (otherwise mint <= 0)
!!    linear      : use the likely probability genotype dam or not
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine change_qtleffect(nteff,iq,chr,n,xinc,incidenceDesc,mint,linear)
       integer                                 ,intent(in)     :: nteff,iq,chr,n
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc
       integer                                 ,intent(in)     :: mint ! index of interaction to remove from incidence
       logical                                 ,intent(in)     :: linear

       !local
       integer :: nbint,jef,nivint,ntlev,ip,jm,kd,ntt,km,ntniv,nkd
       logical :: disab_qtl_sire
      ! real (kind=dp) ,dimension(np+nm)         :: phases
       type(GENEALOGY_BASE) , pointer :: dg
       type(DATAMODEL_BASE) , pointer :: dpm
       type(PHENOTYPE_BASE) , pointer :: dpa

       call log_mess("module incidence : change_qtleffect",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
       dpm => incidenceDesc%dataset%data%phenoModel
       dpa => incidenceDesc%dataset%data%phenoAnimal

       disab_qtl_sire = incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL)

       nbint=dpm%listModelTrait(incidenceDesc%ic)%nbint(iq)
       ntniv = incidenceDesc%desc(nteff)%start -1
       nkd=0

      if ( .not. disab_qtl_sire ) then
       do ip=1,dg%np
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if ( linear ) then
            call getprob_lin(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr,n,incidenceDesc%ic,ip,jm,&
             dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp,incidenceDesc%pm)
           else
            call getglobalprob(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr,n,incidenceDesc%ic,ip,jm,&
             dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp,incidenceDesc%pm)
           end if
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if ( incidenceDesc%dataset%data%phenoAnimal%presentc(incidenceDesc%ic,kd) ) then
            nkd=nkd+1

           !  si des effets fixes sont en interaction avec l'allele paternel au qtl
           !  il faut estimer les effets qtl pour chacun des niveaux concernes
           !
           !  on commence donc par creer un effet composites regroupant l'ensemble de
           !  ces effets
           !
           nivint=1
           ntlev=1
           !nt=ntniv
           do jef=1,nbint
             if ( mint == jef) cycle ! this interaction is not adding
             nivint=nivint+&
              ntlev*(dpa%nivx(kd,dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(iq,jef))-1)
             ntlev=ntlev*dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(iq,jef))
           end do
           ! mode biall or sire
           if ( incidenceDesc%isBiAllellic ) then
             ntt=ntniv+nivint
             xinc(nkd,ntt)=incidenceDesc%pp(kd)
           else
             ntt=ntniv+nivint+ntlev*(ip-1)
             xinc(nkd,ntt)=incidenceDesc%pp(kd)
           end if

           if(incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,jm)) then
              km=incidenceDesc%dataset%data%phenoAnimal%iam(incidenceDesc%ic,dg%repfem(jm))
              if ( incidenceDesc%isBiAllellic ) then
               ntt= ntniv+nivint
               xinc(nkd,ntt)=incidenceDesc%pm(kd)
              else
               ntt= ntniv+ntlev*dg%np+nivint+ntlev*(km-1)
               xinc(nkd,ntt)=incidenceDesc%pm(kd)
              end if
           end if
          end if
         end do
        end do
       end do
       else
        do ip=1,dg%np
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if ( linear ) then
            call getprob_lin(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr,n,incidenceDesc%ic,ip,jm,&
             dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp,incidenceDesc%pm)
           else
            call getglobalprob(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr,n,incidenceDesc%ic,ip,jm,&
             dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp,incidenceDesc%pm)
           end if
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if ( incidenceDesc%dataset%data%phenoAnimal%presentc(incidenceDesc%ic,kd) ) then
            nkd=nkd+1
            nivint=1
            ntlev=1
            !nt=ntniv
            do jef=1,nbint
             if ( mint == jef) cycle ! this interaction is not adding
             nivint=nivint+&
              ntlev*(dpa%nivx(kd,dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(iq,jef))-1)
             ntlev=ntlev*dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffectWithInteraction(iq,jef))
            end do

            if(incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,jm)) then
               km=incidenceDesc%dataset%data%phenoAnimal%iam(incidenceDesc%ic,dg%repfem(jm))
             !  km = 1 !!!! MEME niveau !!!!! !  print *,"A ENLEVER....KM=1 change_qtleffect"
               ntt= ntniv+nivint+ntlev*(km-1)
               xinc(nkd,ntt)=incidenceDesc%pm(kd)

           end if
          end if
         end do
        end do
       end do
         end if ! DISABLE_SIRE_QTL

       ! a modifier
       ! --> on devrait integerer l intialisation des pdds dans la boucle precedente
       ! attention les interaction qtl ne sont pas encore pris en compte
       if (.not. linear) then
         call pdd_at_position(incidenceDesc%ic,incidenceDesc%dataset,iq,chr,n)
       end if
    end subroutine change_qtleffect
!!***

!!****f* m_qtlmap_incidence/pdd_at_position
!!  NAME
!!    pdd_at_position
!!  DESCRIPTION
!!    Compute the probalities to receive the qtl from the sire
!!
!!  INPUTS
!!    ic          : index of the trait
!!    iq          : index of the qtl
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!    dataset     : The user dataset
!!
!!  NOTES
!!  SOURCE
 subroutine pdd_at_position(ic,dataset,iq,chr,n)
          type (DATASET_TYPE)  , intent(inout)    :: dataset
          integer   ,intent(in)   :: ic,iq,chr,n

          integer :: ip,jm,kd,ig,kkd,igg,ifem,kkkd,kds
          real (kind=dp) :: ppt,pmt

          call log_mess("module incidence : pdd_at_position",DEBUG_DEF)

          do ip=1,dataset%data%genea%np
           dataset%lSires(ip)%ppt(iq,:)=0.d0
           ifem=0
           kkkd=0
           do jm=dataset%data%genea%nmp(ip)+1,dataset%data%genea%nmp(ip+1)
             ifem=ifem+1
             if (dataset%lSires(ip)%full_sib(ifem)%firstkd<0) cycle
             dataset%lSires(ip)%full_sib(ifem)%ppt(iq,:,:)=0.d0
             dataset%lSires(ip)%full_sib(ifem)%pmt(iq,:,:)=0.d0
             igg=0
             kds=kkkd
            ! print *,"ip,jm,firstkd:",ip,jm,lSires(ip)%half_sib%firstKd
             do ig=dataset%spt%ngenom(chr,jm)+1,dataset%spt%ngenom(chr,jm+1)
              igg=igg+1
              kkd=0
              kkkd=kds
              do kd=dataset%spt%ngend(chr,ig)+1,dataset%spt%ngend(chr,ig+1)
                if(dataset%data%phenoAnimal%presentc(ic,dataset%spt%ndesc(chr,kd))) then
                 kkd=kkd+1
	             kkkd=kkkd+1
	             if( dataset%data%phenoAnimal%estime(ic,jm) )then
                  dataset%lSires(ip)%full_sib(ifem)%ppt(iq,igg,kkd)=&
                  -dataset%spt%pdd(chr,kd,1,n)-dataset%spt%pdd(chr,kd,2,n)+dataset%spt%pdd(chr,kd,3,n)+dataset%spt%pdd(chr,kd,4,n)
                  dataset%lSires(ip)%full_sib(ifem)%pmt(iq,igg,kkd)=&
                  -dataset%spt%pdd(chr,kd,1,n)+dataset%spt%pdd(chr,kd,2,n)-dataset%spt%pdd(chr,kd,3,n)+dataset%spt%pdd(chr,kd,4,n)
                 else
                  dataset%lSires(ip)%ppt(iq,dataset%lSires(ip)%half_sib%firstKd+kkkd-1)=&
                  -dataset%spt%pdd(chr,kd,1,n)+dataset%spt%pdd(chr,kd,3,n)
                 end if
               end if
              end do !! kd
             end do !! ig
           end do ! jm
          end do ! ip
      end subroutine pdd_at_position
!!***
!!****f* m_qtlmap_incidence/add_qtlinteraction_sire
!!  NAME
!!    add_qtlinteraction_sire
!!  DESCRIPTION
!!    <DEV> modele interaction
!!  NOTES
!!  SOURCE
    subroutine add_qtlinteraction(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       integer                                 ,intent(in)     :: nteff,chr1,chr2,n1,n2
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       logical                                 ,intent(in)     :: linear
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc

       call log_mess("module incidence : add_qtlinteraction",DEBUG_DEF)

       if ( incidenceDesc%isBiAllellic ) then
        call add_qtlinteraction_biallelic(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       else
        call add_qtlinteraction_sire(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       end if

    end subroutine add_qtlinteraction

!!    <DEV> modele interaction
    subroutine add_qtlinteraction_sire(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       integer                                 ,intent(in)     :: nteff,chr1,chr2,n1,n2
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       logical                                 ,intent(in)     :: linear
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc

       integer :: ntniv,nteffstart,ip,fem,jm,start,nkd
       logical ,dimension(incidenceDesc%dataset%data%genea%nfem) :: femInit
       type(GENEALOGY_BASE) , pointer :: dg

       call log_mess("incidence : ** Add Qtl Interaction effect to estimate **",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea

       nteffstart = incidenceDesc%nteff
       incidenceDesc%nteff = incidenceDesc%nteff+1

       incidenceDesc%desc(nteffstart+1)%name="QTL Interaction effects"
       incidenceDesc%desc(nteffstart+1)%start=incidenceDesc%ntniv+1
       incidenceDesc%desc(nteffstart+1)%end=incidenceDesc%ntniv

       incidenceDesc%desc(nteffstart+1)%haveSubDesc=.true.

       !allocate ( incidenceDesc%desc(nteffstart+1)%listSubDesc(np*4+count(estime(incidenceDesc%ic,:))*16) )
       ! adding only one level by sire for estimate the interaction
       allocate ( incidenceDesc%desc(nteffstart+1)%listSubDesc(dg%np) )

       start=incidenceDesc%desc(nteffstart+1)%start

       do ip=1,dg%np
         ! on estime forcement les effets des haplotype pere QQ Qq qQ qq
          incidenceDesc%desc(nteffstart+1)%listSubDesc(ip)%name="Interaction Sire "//trim(dg%pere(ip))
          incidenceDesc%desc(nteffstart+1)%listSubDesc(ip)%start=start
          incidenceDesc%desc(nteffstart+1)%listSubDesc(ip)%end=start
          start=start+1
       end do
       incidenceDesc%desc(nteffstart+1)%end=incidenceDesc%desc(nteffstart+1)%end+dg%np

       incidenceDesc%ntniv = incidenceDesc%ntniv + dg%np

       incidenceDesc%pp=0
       incidenceDesc%pp2=0

       call change_interaction_effect(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)

    end subroutine add_qtlinteraction_sire

!!    <DEV> modele interaction
    subroutine add_qtlinteraction_biallelic(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       integer                                 ,intent(in)     :: nteff,chr1,chr2,n1,n2
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       logical                                 ,intent(in)     :: linear
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc

       integer :: ntniv,nteffstart,ip,fem,jm,nkd,nef
       type(GENEALOGY_BASE) , pointer :: dg

       dg => incidenceDesc%dataset%data%genea

       call log_mess("incidence : ** Add Qtl Interaction (biallelic) effect to estimate **",DEBUG_DEF)

       nteffstart = incidenceDesc%nteff
       incidenceDesc%nteff = incidenceDesc%nteff+1

       incidenceDesc%desc(nteffstart+1)%name="QTL Interaction effects"
       incidenceDesc%desc(nteffstart+1)%start=incidenceDesc%ntniv+1
      ! incidenceDesc%desc(nteffstart+1)%end=incidenceDesc%ntniv

       incidenceDesc%desc(nteffstart+1)%haveSubDesc=.true.

       if ( count(incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,:)) >0) then
        allocate ( incidenceDesc%desc(nteffstart+1)%listSubDesc(16) )
       else
        allocate ( incidenceDesc%desc(nteffstart+1)%listSubDesc(4) )
       end if

       ntniv=incidenceDesc%desc(nteffstart+1)%start

       ! on estime forcement les effets des haplotype pere QQ Qq qQ qq
       incidenceDesc%desc(nteffstart+1)%listSubDesc(1)%name="i-QQ"
       incidenceDesc%desc(nteffstart+1)%listSubDesc(1)%start=ntniv
       incidenceDesc%desc(nteffstart+1)%listSubDesc(1)%end=ntniv
       ntniv=ntniv+1
       incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%name="i-Qq"
       incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%start=ntniv
       incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%end=ntniv
       ntniv=ntniv+1
       incidenceDesc%desc(nteffstart+1)%listSubDesc(3)%name="i-qQ"
       incidenceDesc%desc(nteffstart+1)%listSubDesc(3)%start=ntniv
       incidenceDesc%desc(nteffstart+1)%listSubDesc(3)%end=ntniv
       ntniv=ntniv+1
       incidenceDesc%desc(nteffstart+1)%listSubDesc(4)%name="i-qq"
       incidenceDesc%desc(nteffstart+1)%listSubDesc(4)%start=ntniv
       incidenceDesc%desc(nteffstart+1)%listSubDesc(4)%end=ntniv

        ! on estime les effets combines pere et mere QQWW QQWw QQww QQwW.....
       if ( count(incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,:)) >0) then
!   A FAIRE.....
!          incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%name="Sire - Dam interactions "
!          incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%start=start
!          incidenceDesc%desc(nteffstart+1)%listSubDesc(2)%end=start+16-1
!          start=start+16
       end if


       incidenceDesc%desc(nteffstart+1)%end=ntniv

       incidenceDesc%ntniv = ntniv

       call change_interaction_effect(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)

    end subroutine add_qtlinteraction_biallelic

!!    <DEV> modele interaction
    subroutine change_interaction_effect(nteff,chr1,chr2,n1,n2,xinc,incidenceDesc,linear)
       integer                                 ,intent(in)     :: nteff,chr1,chr2,n1,n2
       type(INCIDENCE_TYPE)                    ,intent(inout)  :: incidenceDesc
       logical                                 ,intent(in)     :: linear
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)  :: xinc

       !local
       integer :: ip,jm,kd,ntniv,ig1,ig2,kkd,kd2,fem,nkd,scaleXinc
       type(GENEALOGY_BASE) , pointer :: dg

       call log_mess("module incidence : change_interaction_effect",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea

       ntniv = incidenceDesc%desc(nteff)%start -1
       nkd=0

       do ip=1,dg%np
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if ( linear ) then
            call getprobReveive1_lin(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr1,&
              n1,incidenceDesc%ic,ip,jm,dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp,incidenceDesc%pm)
            call getprobReveive1_lin(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,chr2,&
              n2,incidenceDesc%ic,ip,jm,dg%ndm(jm)+1,dg%ndm(jm+1),incidenceDesc%pp2,incidenceDesc%pm2)
           end if
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
           if ( incidenceDesc%dataset%data%phenoAnimal%presentc(incidenceDesc%ic,kd) ) then
            nkd=nkd+1
            if ( incidenceDesc%isBiAllellic ) then
             ! tout les cas d interaction en famille de pere
             xinc(nkd,ntniv+1) = incidenceDesc%pp(kd)*incidenceDesc%pp2(kd)
             xinc(nkd,ntniv+2) = incidenceDesc%pp(kd)*(1.d0-incidenceDesc%pp2(kd))
             xinc(nkd,ntniv+3) = (1d0-incidenceDesc%pp(kd))*incidenceDesc%pp2(kd)
             xinc(nkd,ntniv+4) = (1.d0-incidenceDesc%pp(kd))*(1d0-incidenceDesc%pp2(kd))
!              xinc(nkd,ntniv+1) = (1-incidenceDesc%pp(kd))
            else
             ! une seul interaction pour chaque pere.....
               xinc(nkd,ntniv+ip) = (incidenceDesc%pp(kd)*(1.d0-incidenceDesc%pp2(kd))) + &
              ((1.d0-incidenceDesc%pp(kd))*(1.d0-incidenceDesc%pp2(kd)))
            end if
           end if
         end do ! kd
        end do ! jm
       end do ! ip


    end subroutine change_interaction_effect
!!***


!!****f* m_qtlmap_incidence/add_effcov
!!  NAME
!!    add_effcov
!!  DESCRIPTION
!!    Add the nuisances effects inside a contingence matrix with a description
!!    * An effect is adding for each covariate and fixed effect
!!  INPUTS
!!   meff         : fixed effect to not added in the contingence matrix (if meff <= 0, all fixed effects are adding)
!!   mcov         : covariate effect to not added in the contingence matrix (if mcov <= 0, all covariates effects are adding)
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_effcov(xinc,incidenceDesc,meff,mcov)
       type(INCIDENCE_TYPE)                    , intent(inout) :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(inout) :: xinc
       integer                                 , intent(in)    :: meff   ! indice on effet fixed list to remove (<=0 otherwise)
       integer                                 , intent(in)    :: mcov   ! indice on covariate list to remove (<=0 otherwise)
       integer :: nkd,ntniv,ip,jm,kd,ntnifix
       integer :: nbef,nbco,nteff,jef,jcov,ieff
       integer :: nbniv,nivd,i,j,nbef2,nbco2
       type(GENEALOGY_BASE) , pointer :: dg
       type(DATAMODEL_BASE) , pointer :: dpm
       type(PHENOTYPE_BASE) , pointer :: dpa

       call log_mess("incidence : ** Add eff and cov ** ",DEBUG_DEF)
       dg => incidenceDesc%dataset%data%genea
       dpm => incidenceDesc%dataset%data%phenoModel
       dpa => incidenceDesc%dataset%data%phenoAnimal

       nbef=dpm%listModelTrait(incidenceDesc%ic)%nbfe
       nbco=dpm%listModelTrait(incidenceDesc%ic)%nbco

       nbef2=nbef
       if ( meff /= 0 ) nbef2=nbef-1
       nbco2=nbco
       if ( mcov /=0 ) nbco2=nbco-1


       if (nbef2+nbco2 <=0)return

       nbniv=0
       ieff=incidenceDesc%nteff

       if ( nbef2 > 0 ) then
        ieff=ieff+1
        incidenceDesc%eqtl_print(ieff)=.false.
        do i=1,nbef
         if ( i == meff ) cycle
         nbniv=dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffect(i))+nbniv
        end do
        incidenceDesc%desc(ieff)%isVar = .false.
        incidenceDesc%desc(ieff)%name="Fixed effects"
        incidenceDesc%desc(ieff)%start=incidenceDesc%ntniv+1
        incidenceDesc%desc(ieff)%end=incidenceDesc%ntniv+nbniv
        incidenceDesc%desc(ieff)%haveSubDesc=.true.

        !init of xinc
        xinc(:,incidenceDesc%desc(ieff)%start:incidenceDesc%desc(ieff)%end)=0.d0

        allocate ( incidenceDesc%desc(ieff)%listSubDesc(nbniv) )

        j=0
        do jef=1,nbef
         if ( jef == meff ) cycle
         do i=1,dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffect(jef))
         j=j+1
         incidenceDesc%desc(ieff)%listSubDesc(j)%name=&
          trim(dpm%namefix(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffect(jef)))&
         //' level '//trim(str(i))
         incidenceDesc%desc(ieff)%listSubDesc(j)%start=incidenceDesc%ntniv+j
         incidenceDesc%desc(ieff)%listSubDesc(j)%end=incidenceDesc%ntniv+j
         end do
        end do
        incidenceDesc%nteff = incidenceDesc%nteff+1
       end if

       if ( nbco2 > 0 ) then
        ieff=ieff+1
        incidenceDesc%eqtl_print(ieff)=.false.
        incidenceDesc%desc(ieff)%isVar = .false.
        incidenceDesc%desc(ieff)%name="Covariates"
        incidenceDesc%desc(ieff)%start=incidenceDesc%ntniv+nbniv+1
        incidenceDesc%desc(ieff)%end=incidenceDesc%ntniv+nbniv+nbco2
        incidenceDesc%desc(ieff)%haveSubDesc=.true.

        !init of xinc
        xinc(:,incidenceDesc%desc(ieff)%start:incidenceDesc%desc(ieff)%end)=0.d0

        allocate ( incidenceDesc%desc(ieff)%listSubDesc(nbco2) )

        j=0
        do jcov=1,nbco
         if ( jcov == mcov ) cycle
         j=j+1
         incidenceDesc%desc(ieff)%listSubDesc(j)%name=&
          trim(dpm%namecov(dpm%listModelTrait(incidenceDesc%ic)%indexCovariate(jcov)))
         incidenceDesc%desc(ieff)%listSubDesc(j)%start=incidenceDesc%ntniv+nbniv+j
         incidenceDesc%desc(ieff)%listSubDesc(j)%end=incidenceDesc%ntniv+nbniv+j
        end do
        incidenceDesc%nteff = incidenceDesc%nteff+1
       end if


       nkd = 0
       nbniv=0
       do ip=1,dg%np
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if ( .not. dpa%presentc(incidenceDesc%ic,kd) ) cycle
            ntniv=incidenceDesc%ntniv
            nkd=nkd+1
            !
            !  autres effets fixes
            !
            nbniv=0
            do jef=1,nbef
             if ( jef == meff ) cycle ! this effet fixed is not adding
            ! incidenceDesc%nivdir(kd,nteff+jef)=ntniv+nbniv+nivx(kd,modele(ic,3+jef))
             nivd=ntniv+nbniv+dpa%nivx(kd,dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffect(jef))
             xinc(nkd,nivd)=1
             nbniv=nbniv+dpm%nlev(dpm%listModelTrait(incidenceDesc%ic)%indexFixedEffect(jef))
            end do

            ntnifix=ntniv+nbniv
            !covariate
            j=0
            do jcov=1,nbco
             if ( jcov == mcov ) cycle ! this covariate is not adding
             j=j+1
             ! incidenceDesc%covdir(kd,jcov)=covar(kd,modele(ic,nbt+jcov))
             xinc(nkd,ntnifix+j)=dpa%covar(kd,dpm%listModelTrait(incidenceDesc%ic)%indexCovariate(jcov))
            end do
            !ntniv=ntnifix+nbco
            !nbt=nbt+nbco
       end do
      end do
     end do

     incidenceDesc%ntniv=incidenceDesc%ntniv+nbniv+nbco2
     call log_mess("END incidence : ** Add eff and cov ** ",DEBUG_DEF)
    end subroutine add_effcov
!!***


!!****f* m_qtlmap_incidence/add_other_continue_traits_as_cov
!!  NAME
!!    add_other_continue_traits_as_cov
!!  DESCRIPTION
!!    Add traits (that are not analysing in the current session) as covariate effects inside a contingence matrix with a description
!!    * An effect is adding for each covariate and fixed effect
!!  INPUTS
!!   mcov         : covariate effect to not added in the contingence matrix (if mcov <= 0, all covariates effects are adding)
!!  OUTPUTS
!!    xinc        : The contingence matrix
!!  incidenceDesc : The description of contingence matrix
!!
!!  NOTES
!!  SOURCE
    subroutine add_other_continue_traits_as_cov(xinc,incidenceDesc,mcov,startTraitAsCov,endTraitAsCov)
       type(INCIDENCE_TYPE)                    , intent(inout) :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(inout) :: xinc
       integer                                 , intent(in)    :: startTraitAsCov,endTraitAsCov,mcov   ! indice on covariate list to remove (<=0 otherwise)

       integer :: nbntniv,i,j,jj,nkd,ieff,ip,jm,kd,jcov
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       call log_mess("module incidence : add_other_continue_traits_as_cov",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
       dpa => incidenceDesc%dataset%data%phenoAnimal
       dpm => incidenceDesc%dataset%data%phenoModel

       ieff=incidenceDesc%nteff
       ieff=ieff+1

       nbntniv=0

       do i=startTraitAsCov,endTraitAsCov
         if ( i == incidenceDesc%ic )cycle
         if ( dpm%natureY(i) /= 'r') cycle
         nbntniv=nbntniv+1
       end do
       if ( nbntniv <= 0 ) return
       incidenceDesc%eqtl_print(ieff)=.false.
       incidenceDesc%desc(ieff)%isVar = .false.
       incidenceDesc%desc(ieff)%name="Traits as covariates"
       incidenceDesc%desc(ieff)%start=incidenceDesc%ntniv+1
       incidenceDesc%desc(ieff)%end=incidenceDesc%ntniv+nbntniv
       incidenceDesc%desc(ieff)%haveSubDesc=.true.

       !init of xinc
       xinc(:,incidenceDesc%desc(ieff)%start:incidenceDesc%desc(ieff)%end)=0.d0

       allocate ( incidenceDesc%desc(ieff)%listSubDesc(nbntniv) )

       jj=0;j=0
       do jcov=startTraitAsCov,endTraitAsCov
         if ( jcov == incidenceDesc%ic )cycle ! the current trait...
         if ( dpm%natureY(jcov) /= 'r') cycle ! can not be treat...
         jj=jj+1
         if ( jj == mcov ) cycle
         j=j+1
         incidenceDesc%desc(ieff)%listSubDesc(j)%name=trim(dpm%carac(jcov))
         incidenceDesc%desc(ieff)%listSubDesc(j)%start=incidenceDesc%ntniv+j
         incidenceDesc%desc(ieff)%listSubDesc(j)%end=incidenceDesc%ntniv+j

         nkd=0
         do ip=1,dg%np
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if ( .not. dpa%presentc(incidenceDesc%ic,kd) ) cycle
             nkd=nkd+1
             if ( dpa%presentc(jcov,kd) ) then
               xinc(nkd,incidenceDesc%desc(ieff)%listSubDesc(j)%start)= dpa%y(jcov,kd)
             else
               xinc(nkd,incidenceDesc%desc(ieff)%listSubDesc(j)%start)= 0.d0
             end if
           end do !kd
          end do !jm
         end do !ip
       end do ! jcov
       incidenceDesc%nteff = incidenceDesc%nteff+1
       incidenceDesc%ntniv=incidenceDesc%ntniv+nbntniv

    end subroutine add_other_continue_traits_as_cov

!remove the last effect describe by the INCIDENCE_TYPE
    subroutine remove_last_effect(incidenceDesc,xinc)
       type(INCIDENCE_TYPE)                    , intent(inout) :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(inout) :: xinc

       xinc(:,incidenceDesc%desc(incidenceDesc%nteff)%start:incidenceDesc%desc(incidenceDesc%nteff)%end)=0.d0
       incidenceDesc%ntniv=incidenceDesc%desc(incidenceDesc%nteff)%start-1
       ! deallocation of sub effect
       if (associated(incidenceDesc%desc(incidenceDesc%nteff)%listSubDesc)) then
         deallocate ( incidenceDesc%desc(incidenceDesc%nteff)%listSubDesc )
       end if

       incidenceDesc%nteff=incidenceDesc%nteff-1
    end subroutine remove_last_effect
!!***

!!****f* m_qtlmap_incidence/getprob_lin
!!  NAME
!!    getprob_lin
!!  DESCRIPTION
!!    Get the probabilities to receive qtl in a linear context depending the probabilities of transmission
!!  INPUTS
!!    ic          : index of the trait
!!    ip          : index of sire
!!    jm          : index of dam
!!    nd1         : first index of progenies of family (ip,jm)
!!    nd2         : last index of progenies of family (ip,jm)
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!    probp       : The probabilities to receive the qtl from the sire
!!    probm       : The probabilities to receive the qtl from the dam
!!
!!  NOTES
!!  SOURCE
      subroutine getprob_lin(dataset,spt,chr,n,ic,ip,jm,nd1,nd2,probp,probm)
          type(QTLMAP_DATASET)       ,intent(in) :: dataset
          type(PDD_BUILD)            ,intent(in) :: spt
          integer                 ,intent(in)                 :: chr,n,ic,ip,jm,nd1,nd2
          real (kind=dp)   ,dimension(dataset%genea%nd)   ,intent(out) :: probp,probm

          integer :: kd,ig,kkd
          real (kind=dp) :: p

          call log_mess("module incidence : getprob_lin",DEBUG_DEF)

          probp=0.d0
          probm=0.d0

          !do kkd=nd1,nd2
           do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
            do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
              kkd = spt%ndesc(chr,kd)
              if ( kkd < nd1 .or. kkd > nd2 ) cycle
               if(dataset%phenoAnimal%presentc(ic,spt%ndesc(chr,kkd))) then
	            if( dataset%phenoAnimal%estime(ic,jm) )then
                 probp(kkd)=probp(kkd) + (-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n))
                 probm(kkd)=probm(kkd) + (-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n))
                 else
                 probp(kkd)=probp(kkd) + (-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n))
                end if
               end if
            end do
           end do

      end subroutine getprob_lin
!!***

!!****f* m_qtlmap_incidence/getprobReveive1_lin
!!  NAME
!!    getprobReveive1_lin
!!  DESCRIPTION
!!    Get the probabilities to receive the chromosome one in a linear context
!!  INPUTS
!!    ic          : index of the trait
!!    ip          : index of sire
!!    jm          : index of dam
!!    nd1         : first index of progenies of family (ip,jm)
!!    nd2         : last index of progenies of family (ip,jm)
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!    probp       : The probabilities to receive the qtl from the sire
!!    probm       : The probabilities to receive the qtl from the dam
!!
!!  NOTES
!!  SOURCE
      subroutine getprobReveive1_lin(dataset,spt,chr,n,ic,ip,jm,nd1,nd2,probp,probm)
          type(QTLMAP_DATASET)       ,intent(in) :: dataset
          type(PDD_BUILD)            ,intent(in) :: spt
          integer                 ,intent(in)                 :: chr,n,ic,ip,jm,nd1,nd2
          real (kind=dp)   ,dimension(dataset%genea%nd)   ,intent(out) :: probp,probm

          integer :: kd,ig,kkd
          real (kind=dp) :: p

          call log_mess("module incidence : getprobReveive1_lin",DEBUG_DEF)

          probp=0.d0
          probm=0.d0

       !   do kkd=nd1,nd2
           do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
            do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
              kkd = spt%ndesc(chr,kd)
              if ( kkd < nd1 .or. kkd > nd2 ) cycle
              if(dataset%phenoAnimal%presentc(ic,spt%ndesc(chr,kkd))) then
                if( dataset%phenoAnimal%estime(ic,jm) )then
                 probp(kkd)=probp(kkd) + spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)
                 probm(kkd)=probm(kkd) + spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n)
                else
                 probp(kkd)=probp(kkd) + spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)
                end if
              end if
           end do
          end do

      end subroutine getprobReveive1_lin
!!***

!!****f* m_qtlmap_incidence/getprobld_lin
!!  NAME
!!    getprobld_lin
!!  DESCRIPTION
!!    Get the probabilities to receive an haplotype
!!  INPUTS
!!    ic          : index of the trait
!!    ip          : index of sire
!!    jm          : index of dam
!!    nd1         : first index of progenies of family (ip,jm)
!!    nd2         : last index of progenies of family (ip,jm)
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!    pp_ldla     :
!!    pm_ldla     :
!!
!!  NOTES
!!  SOURCE
 !attention pp_ldla,pm_ldla de taille nd car lineaire (1 seul genotype pour toute les meres)
      subroutine getprobld_lin(dataset,spt,chr,n,ic,ip,jm,nd1,nd2,pp_ldla,pm_ldla)
          type(QTLMAP_DATASET)       ,intent(in) :: dataset
          type(PDD_BUILD)            ,intent(in) :: spt
          integer                 ,intent(in)                 :: chr,n,ic,ip,jm,nd1,nd2
          real (kind=dp)   ,dimension(dataset%genea%nd,2)   ,intent(out) :: pp_ldla,pm_ldla

          integer :: kd,ig,kkd
          real (kind=dp) :: p

          call log_mess("module incidence : getprobld_lin",DEBUG_DEF)

          pp_ldla=0.d0
          pm_ldla=0.d0
        !  ig=likelyprob(chr,ip,jm)
        !  do kkd=nd1,nd2
            do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
              do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
                kkd = spt%ndesc(chr,kd)
                if ( kkd < nd1 .or. kkd > nd2 ) cycle
                 if( dataset%phenoAnimal%presentc(ic,kkd) ) then
	              if( dataset%phenoAnimal%estime(ic,jm) )then
                   pp_ldla(kkd,1)=pp_ldla(kkd,1) + pp_ldla(kkd,1)+spt%probg(chr,ig)*(spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n))
                   pp_ldla(kkd,2)=pp_ldla(kkd,2) + pp_ldla(kkd,2)+spt%probg(chr,ig)*(spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n))
                   pm_ldla(kkd,1)=pm_ldla(kkd,1) + spt%probg(chr,ig)*(spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n))
                   pm_ldla(kkd,2)=pm_ldla(kkd,2) + spt%probg(chr,ig)*(spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,4,n))

                 else
                   pp_ldla(kkd,1)=pp_ldla(kkd,1) + spt%pdd(chr,kd,1,n)
                   pp_ldla(kkd,2)=pp_ldla(kkd,2) + spt%pdd(chr,kd,3,n)

                 end if
               end if
           end do
          end do

      end subroutine getprobld_lin
!!***


!!****f* m_qtlmap_incidence/getglobalprob
!!  NAME
!!    getglobalprob
!!  DESCRIPTION
!!    Get the probabilities to receive a qtleffect in a general context
!!  INPUTS
!!    ic          : index of the trait
!!    ip          : index of sire
!!    jm          : index of dam
!!    nd1         : first index of progenies of family (ip,jm)
!!    nd2         : last index of progenies of family (ip,jm)
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!    probp       :
!!    probm       :
!!
!!  NOTES
!!  SOURCE
     subroutine getglobalprob(dataset,spt,chr,n,ic,ip,jm,nd1,nd2,probp,probm)
          type(QTLMAP_DATASET)       ,intent(in) :: dataset
          type(PDD_BUILD)            ,intent(in) :: spt
          integer                 ,intent(in)                 :: chr,n,ic,ip,jm,nd1,nd2
          real (kind=dp)    ,dimension(dataset%genea%nd)   ,intent(out) :: probp,probm

          integer :: kd,ig,kkd
          real (kind=dp) :: ppt,pmt

          call log_mess("module incidence : getglobalprob",DEBUG_DEF)

          probp=0.d0
          probm=0.d0
          do kkd=nd1,nd2
            do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
             do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
              if ( kkd == spt%ndesc(chr,kd)) then
               if(dataset%phenoAnimal%presentc(ic,spt%ndesc(chr,kkd))) then
	             if( dataset%phenoAnimal%estime(ic,jm) )then
                  ppt=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  probp(kkd)=probp(kkd)+spt%probg(chr,ig)*ppt
                  pmt=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                  probm(kkd)=probm(kkd)+spt%probg(chr,ig)*pmt
                else
                  ppt=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n)
                  probp(kkd)=ppt
                end if
                exit
              end if
             end if
           end do !! kd
          end do !! ig
        end do !! kkd

      end subroutine getglobalprob
!!***


!!****f* m_qtlmap_incidence/get_precision
!!  NAME
!!    get_precision
!!  DESCRIPTION
!!    Compute the precision for each level of the incidence matrix.
!!  INPUTS
!!    xx          : X'.X ,X contingence matrix
!!    xxx         : buffer array
!!    jm          : index of dam
!!    nd1         : first index of progenies of family (ip,jm)
!!    nd2         : last index of progenies of family (ip,jm)
!!    chr         : chromosome localisation
!!    n           : the position
!!  OUTPUTS
!!  incidenceDesc%precis : Precision of each level of the incidence matrix.
!!
!!  NOTES
!!    ginv1
!!  SOURCE
  subroutine get_precision(xx,xxx,incidenceDesc)
      type(INCIDENCE_TYPE)                                ,intent(inout) :: incidenceDesc
      real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)       , intent(in)   :: xx
      real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)       , intent(out)  :: xxx

      integer                                             :: iniv,kniv,jniv,lniv,irank
      real (kind=dp)                                      :: ups

      call log_mess("module incidence : get_precision",DEBUG_DEF)
!
!  Reduction de la matrice d'incidence
!
       xxx(:incidenceDesc%ntniv,:incidenceDesc%ntniv)=0.d0

       kniv=0
       do iniv =1,incidenceDesc%ntniv
         if(incidenceDesc%vecsol(iniv)) then
           kniv=kniv+1
           lniv=0
           do jniv=1,incidenceDesc%ntniv
             if(incidenceDesc%vecsol(jniv)) then
               lniv=lniv+1
               xxx(kniv,lniv)=xx(iniv,jniv)
             end if
           end do
          end if
       end do
!
!  inversion de la matrice d'incidence reduite pour calcul des precisions
!

      ups=1.d-15
      call ginv1(xxx,lniv,size(xxx,1),ups,irank)

      jniv=0
      do iniv=1,incidenceDesc%ntniv
        incidenceDesc%precis(iniv)=9999.d0
        if(incidenceDesc%vecsol(iniv)) then
          jniv=jniv+1
         incidenceDesc%precis(iniv)=xxx(jniv,jniv)
        end if
      end do

      end subroutine get_precision
!!***

!!****f* m_qtlmap_incidence/get_correlations
!!  NAME
!!    get_correlations
!!  DESCRIPTION
!!    Compute the correlation between level of the incidence matrix.
!!  INPUTS
!!  incidenceDesc : the description of the incidence matrix.
!!    xxx         : buffer array provide by get_precision function
!!  OUTPUTS
!!    xcor        : the correlation
!!
!!  NOTES
!!    xxx come from the get_precision function.
!!  SOURCE
     subroutine get_correlations(incidenceDesc,xxx,xcor)
       type(INCIDENCE_TYPE)                          ,intent(in)  :: incidenceDesc
       real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(in)  :: xxx
       real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out) :: xcor
       integer ::ief,jef,i,j

       call log_mess("module incidence : get_correlations",DEBUG_DEF)

       xcor=0.d0
       ief=0
       do i=1,incidenceDesc%ntniv
         if(incidenceDesc%vecsol(i))then
            ief=ief+1
            jef=0
            do j=1,incidenceDesc%ntniv
              if(incidenceDesc%vecsol(j))then
                jef=jef+1
                if (xxx(ief,ief) == 0 .or. xxx(jef,jef) == 0) then
                  xcor(i,j)=0.d0
                else
                  xcor(i,j)=xxx(ief,jef)/dsqrt(xxx(ief,ief)*xxx(jef,jef))
                end if
               end if
            end do
         end if
       end do
      end subroutine get_correlations
!!***

!!****f* m_qtlmap_incidence/confusion_between_qtl_and_effect
!!  NAME
!!    confusion_between_qtl_and_effect
!!  DESCRIPTION
!!   Detect a confusion betwwen a qtl effect and the other effect described by the model.
!!  INPUTS
!!  incidenceDesc : the description of the incidence matrix.
!!    xcor        : the correlation matrix
!!  index_effect  : the index qtl effect to test with the other effect
!!  OUTPUTS
!!    sizealert     : number of alert detected
!!    alertCorrQtl  : the correlation
!!    corrmax       : the maximum correlation finded succesively
!!
!!  NOTES
!!    xxx come from the get_precision function.
!!  HISTORY
!!    14/09/2010 : correction bug dimension de femInit
!!  SOURCE
      subroutine confusion_between_qtl_and_effect(incidenceDesc,xcor,index_effect,&
      sizealert,alertCorrQtl,nalert,corrmax)
        type(INCIDENCE_TYPE)                          ,intent(in)  :: incidenceDesc
        real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(in)  :: xcor
        integer                                       ,intent(in)  :: index_effect
        integer                                       ,intent(in)  :: sizealert
        type(CORR_ALERT_TYPE) ,dimension(sizealert)   ,intent(out) :: alertCorrQtl
        integer                                       ,intent(out) :: nalert
        real (kind=dp)                              ,intent(inout) :: corrmax

        integer :: numqtl,ip,jm,iqtl,jniv,ntlev,id,s,e,i,fem,nnp
        logical ,dimension(:),allocatable :: femInit
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa

        call log_mess("module incidence : confusion_between_qtl_and_effect",DEBUG_DEF)

        dg => incidenceDesc%dataset%data%genea
        dpa => incidenceDesc%dataset%data%phenoAnimal

        nalert=0
        allocate (femInit(maxval(dpa%iam)))
        !
        !  on caclule les correlations entre les colonnes de X'X des effets qtl pere et de
        !  tous les autres effets
        !
        do numqtl=1,incidenceDesc%nqtl
         if ( incidenceDesc%isBiAllellic ) then
           nnp = 1
         else
           nnp = dg%np
         end if

         if ( .not. incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL) ) then
         do ip=1,nnp
           ntlev=0
           !print *,"ip:",ip
           if ( incidenceDesc%isBiAllellic ) then
            s=incidenceDesc%ntniv_qtlsires(numqtl)
           else
            s=incidenceDesc%ntniv_qtlsires(numqtl)+(ip-1)*incidenceDesc%ntlev
           end if

           e=s+incidenceDesc%ntlev-1

           do iqtl=s,e
             ntlev=ntlev+1

             if (incidenceDesc%vecsol(iqtl))then
              !
              !  test d'une confusion avec l'effect concerné
              !

              do jniv=incidenceDesc%desc(index_effect)%start,incidenceDesc%desc(index_effect)%end
               if (abs(xcor(iqtl,jniv))>=corrmax) corrmax=xcor(iqtl,jniv)
               !  si la correlation est superieure a corMax, on imprime une alerte

               if(dabs(xcor(iqtl,jniv)).ge.incidenceDesc%dataset%data%params%THRES_CONFUSION)then
                nalert=nalert+1
                alertCorrQtl(nalert)%qtl=numqtl
                alertCorrQtl(nalert)%ip=ip
                alertCorrQtl(nalert)%corr=xcor(iqtl,jniv)
                alertCorrQtl(nalert)%ntlev=ntlev
                id=jniv-incidenceDesc%desc(index_effect)%start+1
                alertCorrQtl(nalert)%name_effect=incidenceDesc%desc(index_effect)%name
                alertCorrQtl(nalert)%name_level=incidenceDesc%desc(index_effect)%listSubDesc(id)%name
               end if
              end do ! jniv
            end if
           end do ! iqtl
!
!  fin pour le niveau du qtl du pere
!
        end do ! ip

        end if

     if ( .not. incidenceDesc%isBiAllellic .and. count (dpa%estime(incidenceDesc%ic,:))>0 ) then
        femInit=.false.
         do jm=1,dg%nm
          if ( .not. incidenceDesc%isBiAllellic ) then
           if (dpa%estime(incidenceDesc%ic,jm) ) then
             fem=dpa%iam(incidenceDesc%ic,dg%repfem(jm))
             if (femInit(fem)) then
                cycle
             end if
             femInit(fem)=.true.
          else
             fem=1
          end if

           ntlev=0
           s=incidenceDesc%ntniv_qtldams(numqtl)+(fem-1)*incidenceDesc%ntlev
           e=s+incidenceDesc%ntlev-1
           do iqtl=s,e
        !     print *,"iqtl:",iqtl
             ntlev=ntlev+1
             if (incidenceDesc%vecsol(iqtl))then
              !
              !  test d'une confusion avec l'effect concerné
              !
              do jniv=incidenceDesc%desc(index_effect)%start,incidenceDesc%desc(index_effect)%end
          !     print *,'jniv:',jniv
           !    print *,abs(xcor(iqtl,jniv))
               if (abs(xcor(iqtl,jniv))>=corrmax) corrmax=xcor(iqtl,jniv)
               !  si la correlation est superieure a corMax, on imprime une alerte
               !
               if(dabs(xcor(iqtl,jniv)).ge.incidenceDesc%dataset%data%params%THRES_CONFUSION)then
                nalert=nalert+1
                alertCorrQtl(nalert)%qtl=numqtl
                alertCorrQtl(nalert)%jm=jm
                alertCorrQtl(nalert)%corr=xcor(iqtl,jniv)
                alertCorrQtl(nalert)%ntlev=ntlev
                id=jniv-incidenceDesc%desc(index_effect)%start+1
                alertCorrQtl(nalert)%name_effect=incidenceDesc%desc(index_effect)%name
                alertCorrQtl(nalert)%name_level=incidenceDesc%desc(index_effect)%listSubDesc(id)%name
               end if
              end do ! jniv
             end if
            end do ! iqtl

!
!  fin pour le niveau du qtl du pere
!
         end if
        end do ! jm

        end if

      end do ! numqtl
      deallocate (femInit)

      end subroutine confusion_between_qtl_and_effect
!!***

!!****f* m_qtlmap_incidence/confusion_type1
!!  NAME
!!    confusion_type1
!!  DESCRIPTION
!!
!!  INPUTS
!!  incidenceDesc : the description of the incidence matrix.
!!    xxx         :
!!  OUTPUTS
!!    workstruct  :
!!
!!  NOTES
!!    xxx come from the get_precision function.
!!  SOURCE
      subroutine confusion_type1(hypothesis,incidenceDesc,xxx,workstruct)
        integer                                       ,intent(in)  :: hypothesis
        type(INCIDENCE_TYPE)                          ,intent(in)  :: incidenceDesc
        real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(in)  :: xxx
        type(INCIDENCE_GEN_STRUCT)            ,intent(inout)       :: workstruct

        real (kind=dp)  ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)                     :: xcor
        type(CORR_ALERT_TYPE) ,dimension(:),allocatable  :: alerts
        type(CORR_ALERT_TYPE) ,dimension(:),allocatable  :: as
        integer :: ialert,sizealert,eff,i,j,starteff
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa

        call log_mess("module incidence : confusion_type1",DEBUG_DEF)

        dg => incidenceDesc%dataset%data%genea
        dpa => incidenceDesc%dataset%data%phenoAnimal

        call get_correlations(incidenceDesc,xxx,xcor)
        workstruct%corrmax(hypothesis)=0.d0

        sizealert=( dg%np+count(dpa%estime(incidenceDesc%ic,:)) )*incidenceDesc%ntlev*incidenceDesc%nbnivest

        allocate (alerts(sizealert))
        allocate (as(sizealert*incidenceDesc%nqtl))

        workstruct%nalert(hypothesis)=0
        starteff = workstruct%listnteff(workstruct%nqtl)+1
        if (count(dpa%estime(incidenceDesc%ic,:))>0) starteff = starteff+1
        do eff=starteff,incidenceDesc%nteff
           call log_mess("computing confusion beetween qtls and:"//incidenceDesc%desc(eff)%name,VERBOSE_DEF)
           call confusion_between_qtl_and_effect(incidenceDesc,xcor,&
           eff,sizealert,alerts,ialert,workstruct%corrmax(hypothesis))
           if (ialert>0) then
      !       print *,"ialert:",alerts(1)%corr,alerts(2)%corr,alerts(3)%corr
             as(workstruct%nalert(hypothesis)+1:workstruct%nalert(hypothesis)+ialert)=alerts(:ialert)
             workstruct%nalert(hypothesis)=workstruct%nalert(hypothesis)+ialert
           end if
        end do

        if ( workstruct%nalert(hypothesis) <= size(workstruct%alertQtl,2) ) then
         workstruct%alertQtl(hypothesis,:workstruct%nalert(hypothesis))=as(:workstruct%nalert(hypothesis))
        else
          workstruct%alertQtl(hypothesis,:size(workstruct%alertQtl,2))=as(:size(workstruct%alertQtl,2))
          workstruct%nalert(hypothesis)=size(workstruct%alertQtl,2)
        end if

        deallocate (alerts)
        deallocate (as)

      end subroutine confusion_type1
!!***

!!****f* m_qtlmap_incidence/get_incidence_reduit
!!  NAME
!!    get_incidence_reduit
!!  DESCRIPTION
!!
!!  INPUTS
!!  incidenceDesc : the description of the incidence matrix.
!!    xinc        : the contingence matrix
!!  OUTPUTS
!!    xincred     : the incontingence matrixithout level which are not estimable (the reduction contingence matrix)
!!
!!  NOTES
!!  SOURCE
    subroutine get_incidence_reduit(incidenceDesc,xinc,xincred)
      type(INCIDENCE_TYPE)                               ,intent(in)          :: incidenceDesc
      real (kind=dp)   ,dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax),intent(in)  :: xinc
      real (kind=dp)   ,dimension(incidenceDesc%dataset%nkd,incidenceDesc%nbnivest)       ,intent(inout)  :: xincred

      integer :: i,j

      call log_mess("module incidence : get_incidence_reduit",DEBUG_DEF)

      j=0
      do i=1,incidenceDesc%ntniv
         if ( incidenceDesc%vecsol(i) ) then
           j=j+1
           xincred(:,j)=xinc(:incidenceDesc%dataset%nkd,i)
         end if
      end do

    end subroutine get_incidence_reduit
!!***

!!****f* m_qtlmap_incidence/fill_nonull_elements
!!  NAME
!!    fill_nonull_elements
!!  DESCRIPTION
!!    Compacte a sparse matrice.
!!  INPUTS
!!    xincred              : the reduction contingence matrix
!!  OUTPUTS
!!    incidenceDesc%nonull : the compacted array with elements not null.
!!  incidenceDesc%nbnonull : index in the xincred of no null element by line
!!
!!  NOTES
!!  SOURCE
    subroutine fill_nonull_elements(incidenceDesc,s1,s2,xincred)
      integer                     ,intent(in)         :: s1,s2
      type(INCIDENCE_TYPE)       ,intent(inout)       :: incidenceDesc
      real (kind=dp)   ,dimension(s1,s2) ,intent(in)  :: xincred
      integer :: i,j,k

      call log_mess("module incidence : fill_nonull_elements",DEBUG_DEF)

      do i=1,incidenceDesc%dataset%nkd
       k=0
       do j=1,incidenceDesc%nbnivest
        if ( xincred(i,j) /= 0.d0 ) then
         k=k+1
         incidenceDesc%nonull(i,k)=j
        end if
       end do
       incidenceDesc%nbnonull(i)=k
      end do

    end subroutine fill_nonull_elements
!!***

!!****f* m_qtlmap_incidence/matmul_incidence
!!  NAME
!!    matmul_incidence
!!  DESCRIPTION
!!    multiplication matrix routine using the compacted sparse matrix structure (incidenceDesc%nonull,incidenceDesc%nbnonull)
!!  INPUTS
!!    kd1                  : first line index
!!    kd2                  : last line index
!!    incidenceDesc%nonull : the compacted array with elements not null.
!!  incidenceDesc%nbnonull : index in the xincred of no null element by line
!!   sizeInXinc1           : size of the first dimension of xinc
!!   sizeInXinc2           : size of the second dimension of xinc
!!   in_X                  : the vector to multiplied
!!   SizeOutVec            : size of the result vector
!!   startRes              : first index in the result vector to set result
!!
!!  OUTPUTS
!!   out_vec               : the result vector
!!
!!  NOTES
!!  SOURCE
    subroutine matmul_incidence(kd1,kd2,incidenceDesc,sizeInXinc1,sizeInXinc2,in_xinc,in_X,SizeOutVec,startRes,out_vec)
      integer                                      ,intent(in)          :: kd1,kd2,SizeOutVec,startRes,sizeInXinc1,sizeInXinc2
      type(INCIDENCE_TYPE)                         ,intent(in)          :: incidenceDesc
      real (kind=dp)   ,dimension(sizeInXinc1,sizeInXinc2) ,intent(in)  :: in_xinc
      real (kind=dp)   ,dimension(incidenceDesc%nbnivest)  ,intent(in)  :: in_X
      real (kind=dp)   ,dimension(SizeOutVec)              ,intent(out) :: out_vec

      integer :: i,j,k

      call log_mess("module incidence : matmul_incidence",DEBUG_DEF)

      out_vec(:) = 0.d0
      k=startRes
      do i=kd1,kd2
        do j=1,incidenceDesc%nbnonull(i)
           out_vec(k) = out_vec(k) + in_X(incidenceDesc%nonull(i,j))*in_xinc(i,incidenceDesc%nonull(i,j))
        end do
        k=k+1
      end do

    end subroutine matmul_incidence
!!***

!!****f* m_qtlmap_incidence/model_XT_X
!!  NAME
!!    model_XT_X
!!  DESCRIPTION
!! X contingence matrix
!! Build X' X : Model linear homoscedastic without cd
!!              or
!!              Non linear model (modlin)
!!
!! Dimension of X : ND line NTNIV Column
!!
!!                Mu sire1 ...sireN ..MatEff1 Mateff
!!     Animal1    1
!!     ....       1
!!     AnimalN    1
!!  INPUTS
!!    xinc                 : the contingence matrix
!!  OUTPUTS
!!    XX                   : X' X
!!
!!  NOTES
!!    Context homoscedastic
!!  SOURCE
      subroutine model_XT_X(xinc,incidenceDesc,XX)
         type(INCIDENCE_TYPE)                         ,intent(in)   :: incidenceDesc
         real (kind=dp)       , dimension(:,:), intent(in)  :: xinc
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax), intent(out) :: XX

         real (kind=dp)       , dimension(incidenceDesc%dataset%nkd,incidenceDesc%ntniv) :: temp
         integer :: i

         call log_mess("module incidence : model_XT_X",DEBUG_DEF)

         temp = xinc(:incidenceDesc%dataset%nkd,:incidenceDesc%ntniv)
         !
         !  construction de X'X
         !
         XX(:incidenceDesc%ntniv,:incidenceDesc%ntniv)=matmul(transpose(temp),temp)

      end subroutine model_XT_X
!!***

!!****f* m_qtlmap_incidence/set_VInv
!!  NAME
!!    set_VInv
!!  DESCRIPTION
!!  Build V and VInv diagonal matrix with :
!!    V(k,k)    = CD(i,j,k) / var(i)^2
!!    VInv(k,k) = var(i)^2 / CD(i,j,k)
!!                                        i sire,j dam ,k progeny
!!
!!  INPUTS
!!    incidenceDesc       : description of the incidence matrix
!!  OUTPUTS
!!    V                   : CD(i,j,k) / var(i)^2 , i sire,j dam ,k progeny.
!!    VInv                : var(i)^2 / CD(i,j,k) , i sire,j dam ,k progeny.
!!  NOTES
!!    Context heteroscedastic
!!  SOURCE
      subroutine set_VInv( incidenceDesc,sig, VInv )
         type(INCIDENCE_TYPE)          , intent(in)      :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np), intent(in)      :: sig
         real (kind=dp) , dimension(incidenceDesc%dataset%nkd,incidenceDesc%dataset%nkd), intent(out)  :: VInv
         integer :: l,ip,jm,kd
         real (kind=dp) :: sigsquare

         call log_mess("module incidence : set_VInv",DEBUG_DEF)

         VInv=0.d0

         do ip=1,incidenceDesc%dataset%data%genea%np
          sigsquare = sig(ip)*sig(ip)
          do kd=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd,incidenceDesc%dataset%lSires(ip)%half_sib%lastKd
           VInv(kd,kd)=incidenceDesc%dataset%cd(1,kd)/sigsquare
          end do
         end do

      end subroutine set_VInv
!!***

!!****f* m_qtlmap_incidence/model_XT_V_X
!!  NAME
!!    model_XT_V_X
!!  DESCRIPTION
!!  Build X' V-1 X :
!!
!!   X : the contingence matrix
!!   V-1 : diagonal matrix / VInv(k,k) = var(i)^2 / CD(i,j,k) , i sire,j dam ,k progeny
!!  INPUTS
!!    incidenceDesc       : the description of the contingence matrix
!!    xinc                : the contingence matrix
!!    VInv                : diagonal matrix / VInv(k,k) = var(i)^2 / CD(i,j,k) , i sire,j dam ,k progeny
!!  OUTPUTS
!!    XVX                 : X' V-1 X
!!
!!  NOTES
!!    Context heteroscedastic
!!  SOURCE
      subroutine model_XT_V_X(xinc,incidenceDesc,VInv,XVX)
         type(INCIDENCE_TYPE)                         ,intent(in)  :: incidenceDesc
         real (kind=dp)       , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         real (kind=dp)       , dimension(incidenceDesc%dataset%nkd,incidenceDesc%dataset%nkd) , intent(in):: VInv
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax), intent(out)    :: XVX


         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%dataset%nkd)    :: XVXTemp
         real (kind=dp) , dimension(incidenceDesc%dataset%nkd,incidenceDesc%ntniv)    :: temp

       !  real(kind=dp) :: t1
         integer :: i,j,k
         !
         !  construction de X' V-1 X
         !
         temp = xinc(:incidenceDesc%dataset%nkd,:incidenceDesc%ntniv)

        ! t1=time()
         XVXTemp=matmul(transpose(temp),VInv)
         XVX(:incidenceDesc%ntniv,:incidenceDesc%ntniv)=matmul(XVXTemp,temp)
        ! t1=t1-time()
        ! call log_mess("model_XT_V_X sec:"//trim(str(t1)))

      end subroutine model_XT_V_X
!!***


!!****f* m_qtlmap_incidence/estim_cholesky
!!  NAME
!!    estim_cholesky
!!  DESCRIPTION
!!      Estimability with a decomposition of cholesky method
!!
!!      X'.X = L . L*
!!      where L is a lower triangular matrix with strictly positive diagonal entries, and L* denotes the conjugate transpose of L
!!
!!      output : incidenceDesc%vecsol  : array of boolean : true if estim , false otherwise
!!              incidenceDesc%nbnivest : number of effect to estimate according the decompostion
!!  INPUTS
!!    XX                         : X'.X
!!    sizeTriangComp             : size of the vector result
!!  OUTPUTS
!!    incidenceDesc%vecsol       : a logical array with the estimability of each level incidence matrix
!!    incidenceDesc%nbnivest     : the number of estimable level
!!    triangComp                 : X' V-1 X
!!
!!  NOTES
!!    triangComp is using in the resolution for the linear context
!!
!!  SOURCE
       subroutine estim_cholesky(XX,incidenceDesc,sizeTriangComp,triangComp)
         type(INCIDENCE_TYPE)                                                     , intent(inout) :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax), intent(in)    :: XX
         integer                                                                   ,intent(in)    :: sizeTriangComp
         real (kind=dp) , dimension(sizeTriangComp,sizeTriangComp),intent(out) :: triangComp
         integer                                 :: i,j,k,ief,jef
         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%ntniv) :: triang

         call log_mess("module incidence : estim_cholesky",DEBUG_DEF)

         !  application de la decomposition de choleski pour determiner les contraintes
         !
         !  Dans le vecteur vecsol les effets a estimer sont indiques a TRUE et ceux
         !  qui ne sont pas estimables a FALSE
         !
         !  La matrice triang est la matrice triangulaire suparieure M telle que
         !  M'M = X'X
         !
         if ( sizeTriangComp < incidenceDesc%ntniv ) then
           call stop_application("Devel error in m_incidence:estim_cholesky")
         end if

         incidenceDesc%vecsol=.true.
         triang=0.d0
         triangComp=0.d0
         do j=1,incidenceDesc%ntniv
           triang(j,j)=xx(j,j)
           do k=1,j-1
            triang(j,j)=triang(j,j)- triang(j,k)*triang(j,k)
           end do
           if(triang(j,j).gt. incidenceDesc%dataset%data%params%SEUIL_CHO) then
             incidenceDesc%nbnivest=incidenceDesc%nbnivest+1
             triang(j,j)=dsqrt(triang(j,j))
             do i=j+1,incidenceDesc%ntniv
               triang(i,j)=xx(i,j)
               do k=1,j-1
                triang(i,j)=triang(i,j)-triang(i,k)*triang(j,k)
               end do
               triang(i,j)=triang(i,j)/triang(j,j)
               triang(j,i)=triang(i,j)
             end do
           else
            incidenceDesc%vecsol(j)=.false.
           end if
         end do

         i=0
         incidenceDesc%nbnivest=0
         do ief=1,incidenceDesc%ntniv
          if(incidenceDesc%vecsol(ief))then
           incidenceDesc%nbnivest=incidenceDesc%nbnivest+1
           i=i+1
           j=0
           do jef=1,incidenceDesc%ntniv
             if (incidenceDesc%vecsol(jef)) then
               j=j+1
               triangComp(i,j) = triang(ief,jef)
             end if
           end do
       !    corniv(ief)=nbnivest
          end if
         end do

       end subroutine estim_cholesky
!!***


!!****f* m_qtlmap_incidence/set_MR
!!  NAME
!!    set_MR
!!  DESCRIPTION
!!     reduction Matrix :
!!       Mr = Xr' Vestim -1 Xr
!!
!!  INPUTS
!!    XX                         : X'.X
!!    incidenceDesc              : the description of contingence matrix
!!    sizeTriangComp             : size of the vector result
!!  OUTPUTS
!!    Mr                         : the number of estimable level
!!
!!  NOTES
!!    This function are not using....Mr was using for the resolution of Bestim. We are
!!    using the LU resolution (with L  lower triangular matrix from cholesky decomposition)
!!
!!  SOURCE
        subroutine set_MR(XX,incidenceDesc,Mr)
         type(INCIDENCE_TYPE)                         , intent(in)    :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax), intent(in)    :: XX
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) , intent(out)  :: Mr

         integer :: ki,kj,i,j,ix

         call log_mess("module incidence : set_MR",DEBUG_DEF)
         !
         !  compactage de triang
         !
         Mr=0.d0
         ki=0

         do i=1,incidenceDesc%ntniv
           if(incidenceDesc%vecsol(i)) then
             ki=ki+1
             kj=0
             do j=1,incidenceDesc%ntniv
              if(incidenceDesc%vecsol(j))then
               kj=kj+1
               Mr(ki,kj)=XX(i,j)
              end if
             end do
           end if
         end do

       end subroutine set_MR
!!***

!!****f* m_qtlmap_incidence/set_RHS
!!  NAME
!!    set_RHS
!!  DESCRIPTION
!!   Build reduction of matrix : RHS = Xr'  Xr
!!  INPUTS
!!   incidenceDesc              : the description of contingence matrix
!!   xincreduit                 : the reduction of the contingence matrix
!!  OUTPUTS
!!   RHS                        : Xr'  Xr
!!  NOTES
!!
!!  SOURCE
       subroutine set_RHS(xincreduit,incidenceDesc,RHS)
         type(INCIDENCE_TYPE)                         , intent(in) :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xincreduit
         real (kind=dp) , dimension(incidenceDesc%ntnivmax) , intent(out)        :: RHS

         real (kind=dp) , dimension(incidenceDesc%dataset%nkd,incidenceDesc%nbnivest) :: temp

         call log_mess("module incidence : set_RHS",DEBUG_DEF)

         RHS=0.d0
         temp = xincreduit(:incidenceDesc%dataset%nkd,:incidenceDesc%nbnivest)
         RHS(:incidenceDesc%nbnivest)=matmul(transpose(temp),incidenceDesc%dataset%Y(1,:))

       end subroutine set_RHS
!!***

!!****f* m_qtlmap_incidence/set_RHS_V
!!  NAME
!!    set_RHS
!!  DESCRIPTION
!!   Let the diagonal matrix V -1 = var(i)^2 / CD(i,j,k) , i sire,j dam ,k progeny.
!!
!!   Build reduction of matrix : RHS = Xr' V -1  Xr
!!
!!  INPUTS
!!   incidenceDesc              : the description of contingence matrix
!!   xincreduit                 : the reduction of the contingence matrix
!!   VInv                       : var(i)^2 / CD(i,j,k) , i sire,j dam ,k progeny.
!!  OUTPUTS
!!   RHS                        : Xr' V -1  Xr
!!  NOTES
!!
!!  SOURCE
       subroutine set_RHS_V(xincreduit,VInv,incidenceDesc,RHS)
         type(INCIDENCE_TYPE)                         , intent(in) :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xincreduit
         real (kind=dp) , dimension(incidenceDesc%dataset%nkd,incidenceDesc%dataset%nkd) , intent(in):: VInv
         real (kind=dp) , dimension(incidenceDesc%ntnivmax) , intent(out)        :: RHS

         real (kind=dp) , dimension(incidenceDesc%dataset%nkd,incidenceDesc%nbnivest) :: temp
         real (kind=dp) , dimension(incidenceDesc%dataset%nkd)::II
         real (kind=dp)       , dimension(incidenceDesc%nbnivest,incidenceDesc%dataset%nkd) :: RHSTemp
         integer :: i,kd

         call log_mess("module incidence : set_RHS_V",DEBUG_DEF)

         RHS=0.d0
         temp = xincreduit(:incidenceDesc%dataset%nkd,:incidenceDesc%nbnivest)

         RHSTemp=matmul(transpose(temp),VInv)
         RHS(:incidenceDesc%nbnivest)=matmul(RHSTemp,incidenceDesc%dataset%Y(1,:))

       end subroutine set_RHS_V
!!***

!!****f* m_qtlmap_incidence/set_corrxinc
!!  NAME
!!    set_corrxinc
!!  DESCRIPTION
!!  remove parameter which are not estimable from the original contingence matrix
!!
!!  Sample
!!
!!
!!         ( 1 1 0 )         ( T )          ( 1 1 )
!!   X'X = ( 1 1 0 )  Vecsol ( T )  ==>     ( 1 1 )
!!         ( 1 0 1 )         ( F )          ( 1 0 )
!!
!!  INPUTS
!!   incidenceDesc              : the description of contingence matrix
!!   xinc                       : the contingence matrix
!!  OUTPUTS
!!   XXOUT                      : xinc with column parameter not estimable removed
!!
!!  NOTES
!!
!!  SOURCE
       subroutine set_corrxinc(Xinc,incidenceDesc,XXOUT)
         type(INCIDENCE_TYPE)                         , intent(inout) :: incidenceDesc
         real (kind=dp) , dimension(:,:), intent(in)          :: Xinc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(out) :: XXOUT
         integer :: i,inb

         call log_mess("module incidence : set_corrxinc",DEBUG_DEF)

         inb=0
         incidenceDesc%corr_niv_nivb=-1
         do i=1,incidenceDesc%ntniv
           if(incidenceDesc%vecsol(i)) then
            inb=inb+1
            incidenceDesc%corr_niv_nivb(i)=inb
            XXOUT(:,inb)=Xinc(:,i)
           end if
         end do

       end subroutine set_corrxinc
!!***

!!****f* m_qtlmap_incidence/debug_write_incidence
!!  NAME
!!    debug_write_incidence
!!  DESCRIPTION
!!    Print the contingence matrix and additionnal information.
!!    The user have the possibility to print a part of the contingence matrix
!!
!!  INPUTS
!!   incidenceDesc              : the description of contingence matrix
!!   xinc                       : the contingence matrix
!!  NOTES
!!
!!***
     subroutine debug_write_incidence(xinc,incidenceDesc)
       type(INCIDENCE_TYPE) , intent(in)                    :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) , intent(in) :: xinc

       integer  :: i,maxPrint,kd,j,start,end,ip,jm,jjm
       character(len=500)                        :: mf1,myfmt,myfmt2
       character(len=10) :: c,c2
       logical ::is_ok
       type(GENEALOGY_BASE) , pointer :: dg

       call log_mess("module incidence : debug_write_incidence",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
    !   print *,trim(additional_info)
       print *,""
       print *," **** SUMMARY **** "
       print *,' NAME             ', ' | START INDEX', ' | END INDEX'

       print *,"nombre de niveau         :",incidenceDesc%ntniv
       print *,"nombre de niveau estime  :",incidenceDesc%nbnivest
       print *,"interaction effet fix-qtl:",incidenceDesc%ntlev
       if ( associated (incidenceDesc%ntniv_qtlsires) ) then
         print *,"indice des qtls pere     :",incidenceDesc%ntniv_qtlsires
       end if
       if ( associated(incidenceDesc%ntniv_qtldams) ) then
          print *,"indice des qtls mere     :",incidenceDesc%ntniv_qtldams
       end if
       print *,"famille de pere avec indice de descendant"
       do ip=1,dg%np
          jjm=0
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            jjm=jjm+1
            print *,"ip:",ip," jm:",jm," kd1:",incidenceDesc%dataset%lSires(ip)%full_sib(jjm)%firstKd,&
            " kd2:",incidenceDesc%dataset%lSires(ip)%full_sib(jjm)%lastKd,'( ndm:',dg%ndm(jm)+1,dg%ndm(jm+1),')'
          end do
       end do

       !for each effect in the incidence
       do i=1,incidenceDesc%nteff
          print *,"--",i,'---'
          write (*,FMT="("//trim(str(incidenceDesc%nteff))//"(a22,i4,i4))") &
          trim(incidenceDesc%desc(i)%name), incidenceDesc%desc(i)%start,incidenceDesc%desc(i)%end
          if ( incidenceDesc%desc(i)%haveSubDesc) then
             do j=1,size(incidenceDesc%desc(i)%listSubDesc)
               write (*,*) "   **sub**   ",incidenceDesc%desc(i)%listSubDesc(j)%name,&
                                           incidenceDesc%desc(i)%listSubDesc(j)%start,&
                                           incidenceDesc%desc(i)%listSubDesc(j)%end, &
                                           incidenceDesc%vecsol(incidenceDesc%desc(i)%listSubDesc(j)%start)
             end do
          end if
       end do

       !header
       myfmt="(i5,':',1x"
       do i=1,incidenceDesc%nteff
          myfmt=trim(myfmt)//",'|',f5.1"
          if ( incidenceDesc%desc(i)%end-incidenceDesc%desc(i)%start > 0) then
            myfmt=trim(myfmt)//','//trim(str(incidenceDesc%desc(i)%end-incidenceDesc%desc(i)%start))//'(f5.1)'
          else
            if ( incidenceDesc%desc(i)%end-incidenceDesc%desc(i)%start < 0 ) then
              print *,"**ERREUR** END PLUS PETIT QUE START"
              print *,incidenceDesc%desc(i)%name,incidenceDesc%desc(i)%start,incidenceDesc%desc(i)%end
              stop
            end if
          end if
       end do
       myfmt=trim(myfmt)//")"

       do
         write (*,*) "size xinc  (nkd):",incidenceDesc%dataset%nkd
         write (*,*) "[0-9999] [1-9999] : indice start et indice end de la matrice d incidence "
         write (*,*) "c: sortir "
         read (*,fmt="(a12)") c
         if (trim(c) == "c") exit
         start=1
         end=size(xinc,1)
         c2=trim(next_word(c,is_ok))
         if (is_ok) start=get_int(c2)
         c2=trim(next_word(c,is_ok))
         if (is_ok) end=get_int(c2)
         start=max(1,start)
         end=min(end,size(xinc,1))
         write (*,*) "INDICE",('|'//trim(incidenceDesc%desc(i)%name),i=1,incidenceDesc%nteff)
         write (*,*) "------",(('-',j=1,len(trim(incidenceDesc%desc(i)%name))+1),i=1,incidenceDesc%nteff)

         do kd=start,end
           write (*,FMT=myfmt) kd,(xinc(kd,i),i=1,incidenceDesc%ntniv)
         end do
         write (*,*) "..."
       end do

     end subroutine debug_write_incidence

!!****f* m_qtlmap_incidence/set_solution
!!  NAME
!!    set_solution
!!  DESCRIPTION
!!    Fill incsol the solution from the estimimation give in parameters
!!  INPUTS
!!    hypothesis      : the hypothesis
!!    workstruct      :
!!    sigsquareEstime : residual variance
!!    Bestim          : solution
!!    incidenceDesc   : description of the incidence
!!    lennteffQtl     : length of listnteffQtl
!!    listnteffQtl    : index of each qtl effect inside the contingence matrix
!!  OUTPUTS
!!    incsol          : description of the solution
!!  NOTES
!!
!!  SOURCE
    subroutine set_solution(hypothesis,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,lennteffQtl,listnteffQtl)
       integer                        , intent(in)          :: hypothesis
       type(INCIDENCE_GEN_STRUCT)            ,intent(in)    :: workstruct
       type(INCIDENCE_TYPE) , intent(in)                    :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np), intent(in):: sigsquareEstime
       real (kind=dp) , dimension(incidenceDesc%ntnivmax) , intent(in)    :: Bestim      ! estimation of parameter on the maximum finded
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol
       integer                               , intent(in)   ::lennteffQtl
       integer     , dimension(lennteffQtl)  , intent(in)   ::listnteffQtl
       integer :: i,j,k,ip,maxNbPar,g,nb,INDEX,ntniv,nt,iq
       type(GENEALOGY_BASE) , pointer :: dg
       type(DATAMODEL_BASE) , pointer :: dpm

       call log_mess("module incidence : set_solution",DEBUG_DEF)

       dg => incidenceDesc%dataset%data%genea
       dpm => incidenceDesc%dataset%data%phenoModel

       incsol%hypothesis=hypothesis
       call log_mess( "" , VERBOSE_DEF)
       call log_mess( "-------------------------------", VERBOSE_DEF)
       call log_mess( "Estimation of parameters under "//trim(str(incsol%hypothesis)), VERBOSE_DEF)
       call log_mess( "-------------------------------", VERBOSE_DEF)
       call log_mess( "rank of incidence:"//trim(str(incidenceDesc%ntniv)), VERBOSE_DEF)
       call log_mess( "rank max to resolve :"//trim(str(incidenceDesc%ntnivmax)), VERBOSE_DEF)

       allocate (incsol%sig(1,dg%np))
       allocate (incsol%groupeName(incidenceDesc%nteff))
       allocate (incsol%eqtl_print(incidenceDesc%nteff))
       allocate (incsol%nbParameterGroup(incidenceDesc%nteff))
       if (hypothesis > 0) then
         allocate (incsol%qtl_groupeName(1,hypothesis))

         if ( lennteffQtl < hypothesis ) then
           call stop_application("Devel error : bad init of array ** listnteff **")
         end if

         do i=1,hypothesis
           incsol%qtl_groupeName(1,i)=listnteffQtl(i)
         end do

       end if

       maxNbPar=1
       do i=1,incidenceDesc%nteff
         if ( incidenceDesc%desc(i)%haveSubDesc) then
           nb=0
           do j=1,size(incidenceDesc%desc(i)%listSubDesc)
              nb = nb + incidenceDesc%desc(i)%listSubDesc(j)%end - incidenceDesc%desc(i)%listSubDesc(j)%start + 1
           end do
           maxNbPar = max(maxNbPar,nb)
         end if
       end do

       allocate (incsol%parameterName(incidenceDesc%nteff,maxNbPar))
       allocate (incsol%paramaterValue(incidenceDesc%nteff,maxNbPar))
       allocate (incsol%parameterVecsol(incidenceDesc%nteff,maxNbPar))
       allocate (incsol%parameterPrecis(incidenceDesc%nteff,maxNbPar))

       if (hypothesis > 0) then
        if ( trim(incidenceDesc%incidenceSireHaplo(1,1,1)) /= '' ) then
         allocate (incsol%haplotypes(dg%np,incsol%hypothesis,2))
         allocate (incsol%races_haplotypes(dg%np,incsol%hypothesis,2))
         do ip=1,dg%np
          do iq=1,hypothesis
           do j=1,2
            incsol%haplotypes(ip,iq,j) =  incidenceDesc%incidenceSireHaplo(ip,iq,j)
            incsol%races_haplotypes(ip,iq,j) =  incidenceDesc%incidenceSireHaploRace(ip,iq,j)
           end do
          end do
         end do
        end if
       end if

       do ip=1,dg%np
        incsol%sig(1,ip)=sqrt(sigsquareEstime(ip))*dpm%sigt(incidenceDesc%ic)
       end do

       if ( is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model) ) then
         ! pour l'instant homoscedastic
         ! VARA = (H2 / 1 - H2) * VarE
         allocate (incsol%siga(1))
         incsol%siga(1) = dpm%h2(incidenceDesc%ic) / (1.d0 - dpm%h2(incidenceDesc%ic))
         incsol%siga(1) = incsol%siga(1) * sigsquareEstime(1)

         !est ce bon....
         incsol%siga(1) = sqrt(incsol%siga(1))*dpm%sigt(incidenceDesc%ic)
       end if

       nt=0
       do i=1,incidenceDesc%nteff
          incsol%eqtl_print(i) = incidenceDesc%eqtl_print(i)
          incsol%groupeName(i) = incidenceDesc%desc(i)%name
  !        print *,incidenceDesc%desc(i)%name,incidenceDesc%desc(i)%haveSubDesc,size(incidenceDesc%desc(i)%listSubDesc)
          if ( incidenceDesc%desc(i)%haveSubDesc) then
            incsol%nbParameterGroup(i)=0
            do j=1,size(incidenceDesc%desc(i)%listSubDesc)
              INDEX=0
              do k=incidenceDesc%desc(i)%listSubDesc(j)%start,incidenceDesc%desc(i)%listSubDesc(j)%end
                 INDEX=INDEX+1
                 incsol%nbParameterGroup(i) = incsol%nbParameterGroup(i) +1
                 g=incsol%nbParameterGroup(i)
                 incsol%parameterName(i,g)=incidenceDesc%desc(i)%listSubDesc(j)%name//" "//adjustl(STR(INDEX))
                 ntniv=incidenceDesc%desc(i)%listSubDesc(j)%start+INDEX -1

!                 print *,incsol%parameterName(i,g),i,j,k,ntniv,&
!                 incidenceDesc%desc(i)%listSubDesc(j)%start,incidenceDesc%desc(i)%listSubDesc(j)%end

                 incsol%parameterVecsol(i,g)=incidenceDesc%vecsol(ntniv)
                 if ( incsol%parameterVecsol(i,g) ) then
                    nt = nt + 1
                    incsol%paramaterValue(i,g)= Bestim(nt)*dpm%sigt(incidenceDesc%ic)
                    incsol%parameterPrecis(i,g) = incidenceDesc%precis(ntniv)
                 end if
              end do
            end do
          else
            incsol%nbParameterGroup(i)=1
            incsol%parameterName(i,1)=incidenceDesc%desc(i)%name
            incsol%parameterVecsol(i,1)=incidenceDesc%vecsol(incidenceDesc%desc(i)%start)
            if ( incsol%parameterVecsol(i,1) ) then
              nt = nt + 1
              incsol%paramaterValue(i,1)  = Bestim(nt)*dpm%sigt(incidenceDesc%ic)
              incsol%parameterPrecis(i,1) = incidenceDesc%precis(incidenceDesc%desc(i)%start)
            end if
          end if
       end do

       !specific treatment for general mean
       if ( workstruct%effects%general > 0 ) then
        if ( incsol%parameterVecsol(workstruct%effects%general,1) ) then
         incsol%paramaterValue(workstruct%effects%general,1)  = &
         incsol%paramaterValue(workstruct%effects%general,1)+dpm%xmut(incidenceDesc%ic)
        end if
       end if
    end subroutine set_solution
!!***


!!****f* m_qtlmap_incidence/build_incidence_matrix
!!  NAME
!!    build_incidence_matrix
!!  DESCRIPTION
!!    Build the contingence matrix
!!  INPUTS
!!    hypothesis      : the hypothesis
!!    workstruct      :
!!    curPos          : the initial position
!!    meff            : index of the fixed effect to no add
!!    mcov            : index of the covariate to no add
!!    mint            : index of the fixed effect-qtl interaction to not add
!!    numqtl          : number of qtl to test confusion with the mint interaction
!!  OUTPUTS
!!    incidenceDesc   : description of the incidence
!!    xinc            : the contingence matrix to build
!!  NOTES
!!
!!  SOURCE
    subroutine build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,meff,mcov,mint,numqtl)
       integer                                 ,intent(in)       :: meff,mcov,mint,numqtl
       type(POSITION_LRT_INCIDENCE)            ,intent(in)       :: curPos
       type(INCIDENCE_GEN_STRUCT)              ,intent(inout)    :: workstruct
       type(INCIDENCE_TYPE)                    ,intent(inout)    :: incidenceDesc
       real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) ,intent(inout)    :: Xinc

       !local
       integer :: i,ic2,mint2,idEffect
       type(DATAMODEL_BASE) , pointer :: dpm
       type(HAPLOTYPE_POSITION_BUILD) :: shp

       call log_mess("module incidence : build_incidence_matrix",DEBUG_DEF)

       dpm => incidenceDesc%dataset%data%phenoModel
    !   dg => incidenceDesc%dataset%data%genea

       do idEffect=1,workstruct%effects%nbeff
         call log_mess("incidence : ** Add effect "//trim(str(idEffect)),DEBUG_DEF)
         ! GENERAL MEAN
         ! ------------
         if ( idEffect == workstruct%effects%general .and. dpm%natureY(incidenceDesc%ic)/='i') then
           call add_general_mean(xinc,incidenceDesc)

         ! QTLs
         ! ------------
         else if ( idEffect == workstruct%effects%qtl) then

           do i=1,workstruct%nqtl
            mint2=0
            !add qtl effect/interaction at position n to estim
            if ( numqtl == i ) mint2=mint
            workstruct%listnteff(i) = incidenceDesc%nteff+1
            call add_qtleffect(i,curPos%listChr(i),curPos%listN(i),&
            xinc,incidenceDesc,mint2,workstruct%linear)
           end do

         ! HAPLOTYPES
         ! ------------
         else if ( idEffect == workstruct%effects%haplo) then
           if ( workstruct%nqtl > 0 ) then
            call shp%set(incidenceDesc%dataset%data,incidenceDesc%dataset%spt)
            !haplotype effect
            do i=1,workstruct%nqtl
             workstruct%listnteffhap(i)=incidenceDesc%nteff+1
             call add_haplotype_effect(workstruct%listnteffhap(i),i,curPos%listChr(i),curPos%listN(i),&
               xinc,incidenceDesc,shp,workstruct%hdam,.true.)
            end do
            call shp%free()
           end if
         ! INTERACTION
         ! ------------
         else if ( idEffect == workstruct%effects%inter) then
            print *,curPos%listChr(1),size(curPos%listChr)
            print *,"..."
            !add interaction
            workstruct%nteffinter=incidenceDesc%nteff+1
            call add_qtlinteraction(workstruct%nteffinter,curPos%listChr(1),curPos%listChr(2),&
            curPos%listN(1),curPos%listN(2),xinc,incidenceDesc,workstruct%linear)

        ! POLYGENIC
        ! ------------
        else if ( idEffect == workstruct%effects%polygenic) then
          call add_polygenic_effect(xinc,incidenceDesc)

        ! NUISANCES
        ! ------------
        else if ( idEffect == workstruct%effects%nuis) then
           !add fixed effect and covariate
          call add_effcov(xinc,incidenceDesc,meff,mcov)

        ! TRAITS AS COVARIATE
        ! ------------
        else if ( idEffect == workstruct%effects%traits) then
           call add_other_continue_traits_as_cov(xinc,incidenceDesc,mcov,&
                workstruct%startTraitAsCov,workstruct%endTraitAsCov)
        end if
       end do

    end subroutine build_incidence_matrix
!!***


!!****f* m_qtlmap_incidence/opti_0qtl
!!  NAME
!!    opti_0qtl
!!  DESCRIPTION
!!    computation of the LRT under H0 1 trait
!!  INPUTS
!!    ic              : index of the trait
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!  OUTPUTS
!!    incsol          : the solution of the estimation
!!    workstruct      : residual variance and maximum Likelihood under H0
!!  NOTES
!!
!!  SOURCE
      subroutine opti_0qtl(    dataset, &
                               ic,       &  ! the trait
                               incsol   ,       &  !
                               workstruct, &
                              FUNCT_MODEL)  ! function with a model specific
      type(QTLMAP_DATASET)                    ,intent(inout)    :: dataset
      integer                                 ,intent(in)       :: ic
      type(INCIDENCE_GEN_STRUCT)              ,intent(inout)    :: workstruct
      type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)      :: incsol
      external :: FUNCT_MODEL   ! function with a model specific

      integer :: nkd,ip,i,listnteff(1)

      real (kind=dp) ,dimension(:,:),pointer      :: xinc,xxx
      real (kind=dp) , dimension(:) ,pointer      :: Bestim
      real (kind=dp) , dimension(dataset%genea%np):: sigsquareEstime
      type(INCIDENCE_TYPE)                        :: incidenceDesc
      type(POSITION_LRT_INCIDENCE)                :: curPos
      type(PDD_BUILD)               ,pointer      :: spt => null()

      call log_mess("module incidence : opti_0qtl",DEBUG_DEF)

      !initialisation of contingence matrix
      call init_incidence(dataset,spt,ic,0,incidenceDesc,workstruct)
      workstruct%nqtl=0
      allocate (Bestim(incidenceDesc%ntnivmax))
      allocate (xinc(dataset%genea%nd,incidenceDesc%ntnivmax))
      xinc=0.d0
      Bestim=0.d0
      call build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,0,0,0,0)
      !call debug_write_incidence(xinc,incidenceDesc)
      !call model analysis
      call FUNCT_MODEL(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,.true.)
      call set_solution(0,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,1,listnteff)
      workstruct%sigsquare(0,:,1)=sigsquareEstime
      workstruct%fnqtlsires(0,:)=incidenceDesc%fperemax
      workstruct%fnqtldams(0,:)=incidenceDesc%fmeremax

      call end_incidence(incidenceDesc)
      deallocate (xinc)
      deallocate (Bestim)

      end subroutine opti_0qtl
!!***

#ifdef MANAGE_CUDA


     subroutine get_xinc_fixed_effect(w_cuda,xinc2)
        type(WORK_CUDA)                 ,intent(inout)           :: w_cuda
#ifdef CUDA_SP
        real                     ,dimension(w_cuda%data%genea%nd,w_cuda%nLevelFix) :: xinc2
#else
        real (kind=dp)           ,dimension(w_cuda%data%genea%nd,w_cuda%nLevelFix) :: xinc2
#endif
        integer :: i

        call log_mess("module incidence : get_xinc_fixed_effect",DEBUG_DEF)

        do i=1,w_cuda%nLevelFix
         xinc2(:,i) = w_cuda%xinc(:,w_cuda%corrLevel(i)+1)
        end do

     end subroutine get_xinc_fixed_effect

     subroutine get_partial_xinc0(w_cuda,pos1,pos2,xinc2)
        type(WORK_CUDA)                 ,intent(inout)                        :: w_cuda
        integer                                              , intent(in)     :: pos1,pos2
#ifdef CUDA_SP
        real                     ,dimension(w_cuda%data%genea%nd,w_cuda%incidenceDesc%ntnivmax) :: xinc2
#else
        real (kind=dp)           ,dimension(w_cuda%data%genea%nd,w_cuda%incidenceDesc%ntnivmax) :: xinc2
#endif
       call log_mess("module incidence : get_partial_xinc0",DEBUG_DEF)
       ! call debug_write_incidence(w_cuda%xinc,w_cuda%incidenceDesc)
       ! xinc2 = w_cuda%xinc

     end subroutine get_partial_xinc0

!!****f* m_qtlmap_incidence/opti_0qtl_cuda
!!  NAME
!!    opti_0qtl
!!  DESCRIPTION
!!    compute the LRT under H0 1 trait
!!  INPUTS
!!    ic              : index of the trait
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!  OUTPUTS
!!    incsol          : the solution of the estimation
!!    workstruct      : residual variance and maximum Likelihood under H0
!!  NOTES
!!
!!  SOURCE
      subroutine opti_0qtl_cuda(    data,     &
                                    nsim,     &
                                    ic,       &  ! the trait
                               incsol   ,     &  !
                               workstruct,    &  !
                               ySIMUL,        &  !
                               sizeF,         &  !
                               oAllsig,       &
                              FUNCT_MODEL)  ! function with a model specific
      type(QTLMAP_DATASET)         ,target    ,intent(inout)    :: data
      integer                                 ,intent(in)       :: nsim,ic
      type(INCIDENCE_GEN_STRUCT)   ,target    ,intent(inout)    :: workstruct
      type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)      :: incsol
      external :: FUNCT_MODEL   ! function with a model specific
      real (kind=dp), dimension (nsim+1,data%genea%nd)   , intent(in) :: ySIMUL
      integer ,dimension(data%genea%np)                , intent(out) :: sizeF
      real (kind=dp), dimension (nsim+1,data%genea%np)   , intent(out) :: oAllsig

      integer :: nkd,ip,i,listnteff(1),kd1,kd2,nLevelFix,nLevelVar,j

      real (kind=dp) ,dimension(:,:),pointer      :: xinc,xxx
      real (kind=dp) , dimension(:) ,pointer      :: Bestim
      real (kind=dp) , dimension(data%genea%np)   :: sigsquareEstime
      type(INCIDENCE_TYPE)    ,target             :: incidenceDesc
      type(POSITION_LRT_INCIDENCE)               :: curPos
      integer , dimension(:),allocatable ,target         :: corrLevel
      integer                                     :: nqtl,u,nsim2,nGLTotal,ikd,jm,kd
      integer        ,dimension(:),pointer        :: allmaxpos
      real (kind=dp)                              :: seuil_cholesky
      type(PDD_BUILD)             ,pointer        :: spt => null()

#ifdef CUDA_SP
      real                       , dimension(:,:,:)  ,pointer      :: allLrtSires
      real                       , dimension(:,:)    ,pointer      :: allmaxlrt
      real                       , dimension(:,:)   ,pointer       :: isigsq
      real                       , dimension(:,:,:) ,pointer       :: allbestim

      real                       , dimension (:,:)  ,pointer       :: ySIMULR,Mr
      real                       , dimension (:)    ,pointer       :: CDR
      real                       , dimension (:,:)  ,pointer       :: oAllsigR
#else
      real (kind=dp)           ,dimension(:,:,:,:),pointer         :: allLrtSires
      real (kind=dp)           ,dimension(:,:)   ,pointer          :: allmaxlrt
      real (kind=dp)           , dimension(:,:)   ,pointer         :: isigsq
      real (kind=dp)         , dimension(:,:,:)   ,pointer         :: allbestim
#endif


      type(WORK_CUDA)                             :: w_cuda
      integer :: mode

      call log_mess("module incidence : opti_0qtl_cuda",DEBUG_DEF)
      !initialisation of contingence matrix
      call init_incidence(data,spt,ic,0,incidenceDesc,workstruct)
      workstruct%nqtl=0
      allocate (Bestim(incidenceDesc%ntnivmax))
      allocate (xinc(data%genea%nd,incidenceDesc%ntnivmax))

      xinc=0.d0
      Bestim=0.d0
      call build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,0,0,0,0)

      w_cuda%xinc=>xinc
      w_cuda%incidenceDesc=>incidenceDesc
      w_cuda%workstruct=>workstruct
      w_cuda%data => data
      w_cuda%spt => null()

      !call debug_write_incidence(xinc,incidenceDesc)
      !call model analysis
      !call FUNCT_MODEL(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,.true.)


      do ip=1,data%genea%np
           kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
           kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
           if ( kd2>kd1) then
           sizeF(ip) = kd2-kd1+1
           else
           sizeF(ip) = 0
           end if
      end do

      nLevelFix=incidenceDesc%ntniv
      nLevelVar=0

      w_cuda%nLevelFix=nLevelFix
      w_cuda%nLevelVar=nLevelVar

      nqtl=0
      nGLTotal=1
      nsim2 = nsim+1

      allocate (corrLevel(incidenceDesc%ntniv))
!      allocate (allvecsol(nGLTotal,incidenceDesc%ntniv))
      allocate (allbestim(nGLTotal,nsim2,incidenceDesc%ntniv))
      allbestim=0
      w_cuda%corrLevel=>corrLevel

      ! ** Normalement pas utilise sous H0 **
      ! *************************************
      allocate (isigsq(nsim2,data%genea%np))
      ! INITIALISER SEULEMENT POUR LE CAS HETEROSCEDASTIC
      !---------------------------------------------------
!      do ip=1,np
!        isigsq(:,ip)=incidenceDesc%dataset%lSires(ip)%sig0(1)
!      end do

      ! FIN HETERO
      !---------------------------------------------------

      allocate (allLrtSires(nGLTotal,1,1,data%genea%np))
      allocate (allmaxlrt(1,1))
      allocate (allmaxpos(1))
      allLrtSires = 0.d0


      do u=1,incidenceDesc%ntniv
        corrLevel(u)=u-1
      end do

      mode = 0

      if (is_parameter(RESOLUTION_LINEAR,workstruct%type_model) .and.&
          is_parameter(KIND_HETEROSCEDASTIC,workstruct%type_model).and. &
          is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) mode=1

      if (is_parameter(RESOLUTION_LINEAR,workstruct%type_model) .and.&
          is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and. &
          is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) mode=2

#ifdef CUDA_SP
      allocate (ySIMULR(nsim2,incidenceDesc%dataset%nkd))
      allocate (CDR(incidenceDesc%dataset%nkd))


      CDR(:)=incidenceDesc%dataset%CD(1,:)
      ySIMULR=ySIMUL(:,:incidenceDesc%dataset%nkd)
      if ( is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) then
       allocate (Mr(size(workstruct%M,1),size(workstruct%M,2)))
       Mr = workstruct%M
      else
       allocate (Mr(1,1))
      end if

      allocate (oAllsigR(nsim2,data%genea%np))
      seuil_cholesky=0.3
      oAllsigR=0.d0

      call cuda_model_resolv_genome(data%params%gpu_device_id,     &
                                      mode,                        &
                                      nqtl,                        &
                                      isigsq,                      &
                                      nLevelFix,                   &
                                      nLevelVar,                   &
                                      corrLevel,                   &
                                      w_cuda,                      &
                                      get_xinc_fixed_effect,       &
                                      get_partial_xinc0,           &
                                      ySIMULR,                      &
                                      CDR,                         &
                                      Mr,                          &
                                      nsim2,                       &
                                      data%genea%nd,               &
                                      incidenceDesc%dataset%nkd,   &
                                      data%genea%np,               &
                                      sizeF,                       &
                                      incidenceDesc%ntnivmax,      &
                                      nGLTotal,                    &
                                      seuil_cholesky,              &
                                      allbestim,                   &
                                      oAllsigR,                    &
                                      allLrtSires,                 &
                                      allmaxlrt,                   &
                                      allmaxpos )
      deallocate (ySIMULR)
      deallocate (CDR)
      oAllsig = oAllsigR
      deallocate (oAllsigR)
      deallocate (Mr)

#else
      seuil_cholesky=data%params%SEUIL_CHO

      call cuda_model_resolv_genome(data%params%gpu_device_id,     &
                                      mode,                        &
                                      nqtl,                        &
                                      isigsq,                      &
                                      nLevelFix,                   &
                                      nLevelVar,                   &
                                      corrLevel,                   &
                                      w_cuda,                      &
                                      get_xinc_fixed_effect,       &
                                      get_partial_xinc0,           &
                                      ySIMUL,                      &
                                      incidenceDesc%dataset%CD,    &
                                      workstruct%M,                &
                                      nsim2,                       &
                                      data%genea%nd,               &
                                      incidenceDesc%dataset%nkd,   &
                                      data%genea%np,               &
                                      sizeF,                       &
                                      incidenceDesc%ntnivmax,      &
                                      nGLTotal,                    &
                                      seuil_cholesky,              &
                                      allbestim,                   &
                                      oAllsig,                     &
                                      allLrtSires,                 &
                                      allmaxlrt,                   &
                                      allmaxpos )

#endif

      sigsquareEstime = oAllsig(1,:)

      incidenceDesc%vecsol = .false.
      j=1
      do i=1,incidenceDesc%ntniv
         if ( allbestim(1,1,i) /= 0.d0 ) then
           Bestim(j) = allbestim(1,1,i)
           incidenceDesc%vecsol(i) = .true.
           j=j+1
         end if
      end do

    !  call FUNCT_MODEL(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,.true.)


      call set_solution(0,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,1,listnteff)
      workstruct%sigsquare(0,:,1)=sigsquareEstime
      workstruct%fnqtlsires(0,:)=incidenceDesc%fperemax
      workstruct%fnqtldams(0,:)=incidenceDesc%fmeremax
      !stop

      call end_incidence(incidenceDesc)

      deallocate (xinc)
      deallocate (Bestim)
      deallocate (corrLevel)
      deallocate (isigsq)
!      deallocate (allvecsol)
      deallocate (allbestim)
      deallocate (allLrtSires)
      deallocate (allmaxlrt)
      deallocate (allmaxpos)

      end subroutine opti_0qtl_cuda
!!***

!!****f* m_qtlmap_incidence/opti_nqtl_cuda
!!  NAME
!!    opti_nqtl
!!  DESCRIPTION
!!    compute the LRT under H(nqtl) 1 trait    nqtl>=1
!!  INPUTS
!!    nqtl            : the hypothesis
!!    ic              : index of the trait
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl-1)
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!  OUTPUTS
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl)
!!    incsol          : the solution of the estimation at the maximum finded
!!    lrtsol          : LRT : curves and maximum under H(nqtl)
!!  NOTES
!!    use the iterative scan gen_opti_nqtl function
!!  SOURCE
     subroutine opti_nqtl_cuda(      data,        &
                                     spt,         &
                                     nsim,        &
                                     nqtl,        &  ! Number of qtl
                                        ic,       &  ! the trait
                                workStruct,       &  ! variance / likelihood estimated under NQTL-1 hypothesis
                                 incsol   ,       &  ! incidence solution (estimation of each effect)
                                 lrtsol   ,         &  ! maximum lrt finding at a position
                                 ySIMUL   ,         &
                                 sizeNqtl,          & ! for dim of listLrtSolSim : get the number of QTL asked by the user
                                 listLrtSolSim,     &
                                 sizeF,             &
                                 AllsigQtlMoinUn,           &
                               FUNCT_MODEL)
      type(QTLMAP_DATASET)                    ,intent(inout)       :: data
      type(PDD_BUILD)                         ,intent(inout)       :: spt
      integer                                 ,intent(in)          :: nsim,nqtl,ic,sizeNqtl
      type(INCIDENCE_GEN_STRUCT)              ,intent(inout)       :: workstruct
      type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)         :: incsol
      type(TYPE_LRT_SOLUTION)                 ,intent(out)         :: lrtsol
      external                                                     :: FUNCT_MODEL
      type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(data%phenoModel%ncar,sizeNqtl,nsim)    :: listLrtSolSim
      real (kind=dp), dimension (nsim+1,data%genea%nd)             , intent(in) :: ySIMUL
      integer ,dimension(data%genea%np)                   , intent(in)   :: sizeF
      real (kind=dp), dimension (nsim+1,data%genea%np)   , intent(inout) :: AllsigQtlMoinUn


      integer :: i

      real (kind=dp) ,dimension(:,:),pointer      :: xinc ,xxx    ! incidence matrix and temp for testing nuisance effect
      real (kind=dp) , dimension(:) ,pointer      :: Bestim       ! solution of the estimation
      real (kind=dp) , dimension(data%genea%np)   :: sigsquareEstime
      type(INCIDENCE_TYPE)                        :: incidenceDesc ! description incidence matrix
      type(POSITION_LRT_INCIDENCE)                :: curPos        ! the position with additional info for the model
      integer :: n,j,chr
      type(HAPLOTYPE_POSITION_BUILD)   :: shp

      call log_mess("module incidence : opti_nqtl_cuda",DEBUG_DEF)

      if ( nqtl <= 0 ) then
        call stop_application("Error dev: can not call opti_nqtl_linear with nqtl:"//trim(str(nqtl)))
      end if

      !Position Allocation
      call init_position (data,nqtl,nqtl,curPos)

      !initialisation of incidence matrix
      call init_incidence(data,spt,ic,nqtl,incidenceDesc,workstruct)

      allocate (xinc(data%genea%nd,incidenceDesc%ntnivmax))
      xinc=0.d0
      allocate (Bestim(incidenceDesc%ntnivmax))
      allocate (xxx(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax))

      Bestim=0.d0

      call build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,0,0,0,0)
     ! call debug_write_incidence(xinc,incidenceDesc)
      !initalizing lrtsolution
      call lrtsol%new(data,nqtl)

      lrtsol%lrtmax=-9999.9
      lrtsol%chrmax=1
      lrtsol%nxmax=1

      call gen_opti_nqtl_cuda(data,spt,nsim,nqtl,curPos,workstruct,xinc,incidenceDesc,sigsquareEstime,&
                             Bestim,lrtsol,ySIMUL,sizeF,AllsigQtlMoinUn,sizeNqtl,listLrtSolSim)

      if (  workstruct%effects%haplo>0 ) then
       call shp%set(data,spt)

       do i=1,nqtl
         call add_haplotype_effect(workstruct%listnteffhap(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),&
              xinc,incidenceDesc,shp,workstruct%hdam,.false.)
       end do

        call shp%free()
      end if

      if (  workstruct%effects%qtl>0 ) then
      ! Compute precision at the maximum LRT in position posx
       do i=1,nqtl
        call change_qtleffect(workstruct%listnteff(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),&
          xinc,incidenceDesc,0,workstruct%linear)
        call change_name_qtleffect(workstruct%listnteff(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),incidenceDesc)
       end do

       if ( workstruct%nteffinter /= -1 ) then
          call change_interaction_effect(workstruct%nteffinter,lrtsol%chrmax(0),lrtsol%chrmax(1),&
            lrtsol%nxmax(0),lrtsol%nxmax(1),xinc,incidenceDesc,workstruct%linear)
       end if
      end if

      do i=1,nqtl
       curPos%listChr(i)=lrtsol%chrmax(i-1)
       curPos%listN(i)=lrtsol%nxmax(i-1)
      end do

      !call debug_write_incidence(xinc,incidenceDesc)
      call FUNCT_MODEL(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.true.,workstruct%performConfusion,xxx,.false.)

      if ( workstruct%performConfusion) then
         call confusion_type1(nqtl,incidenceDesc,xxx,workstruct)
       end if

      if ( workstruct%performTestNuis ) then
        call test_nuisances(data,spt,ic,incidenceDesc,lrtsol,curPos,workstruct,sigsquareEstime,FUNCT_MODEL)
      end if
      call set_solution(nqtl,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,nqtl,workstruct%listnteff)

      workstruct%fnqtlsires(nqtl,:)=incidenceDesc%fperemax
      workstruct%fnqtldams(nqtl,:)=incidenceDesc%fmeremax
      workstruct%sigsquare(nqtl,:,1)=sigsquareEstime

      call end_incidence(incidenceDesc)
      deallocate (xinc)
      deallocate (xxx)
      deallocate (Bestim)
      call end_position(curPos)

      end subroutine opti_nqtl_cuda
!!***

      subroutine get_partial_xinc(w_cuda,pos1,pos2,xinc2)
        type(WORK_CUDA)                 ,intent(inout)                    :: w_cuda
        integer                                              , intent(in) :: pos1,pos2
#ifdef CUDA_SP
        real             ,dimension(w_cuda%data%genea%nd,(w_cuda%nLevelVar)*(pos2-pos1+1)):: xinc2
#else
        real (kind=dp)   ,dimension(w_cuda%data%genea%nd,(w_cuda%nLevelVar)*(pos2-pos1+1))  :: xinc2
#endif
        integer                                       :: d,nlin,checkIndice,iqtl,i,scalei
        type(INCIDENCE_TYPE)             ,pointer     :: incidenceDesc
        type(INCIDENCE_GEN_STRUCT)       ,pointer     :: workstruct
        real (kind=dp)  ,dimension(:,:)  ,pointer     :: xinc

        logical :: isLD,disable_sire_qtl
        type(HAPLOTYPE_POSITION_BUILD)   :: shp

        call log_mess("module incidence : get_partial_xinc",DEBUG_DEF)

        call log_mess("get_partial_xinc : build all incidence matrix.....POS1="//trim(str(pos1))&
                         //" POS2="//trim(str(pos2)),INFO_DEF)

       incidenceDesc=>w_cuda%incidenceDesc
       xinc => w_cuda%xinc
       workstruct => w_cuda%workstruct

       isLD = workstruct%effects%haplo>0
       disable_sire_qtl = incidenceDesc%dataset%data%cli%key_exist(incidenceDesc%dataset%data%cli%DISABLE_SIRE_QTL)

       if ( workstruct%effects%haplo > 0 ) call shp%set(w_cuda%data,w_cuda%spt)
       d=1
       do nlin=pos1,pos2
        checkIndice = 0
        if ( isLD ) then
         do iqtl=1,workstruct%nqtl
            call add_haplotype_effect(workstruct%listnteffhap(iqtl),iqtl,w_cuda%lchr(nlin,iqtl),w_cuda%lnpos(nlin,iqtl),&
                  xinc,incidenceDesc,shp,workstruct%hdam,.false.)
            checkIndice = checkIndice + incidenceDesc%desc(workstruct%listnteffhap(iqtl))%end &
                        - incidenceDesc%desc(workstruct%listnteffhap(iqtl))%start + 1
            do i=incidenceDesc%desc(workstruct%listnteffhap(iqtl))%start,incidenceDesc%desc(workstruct%listnteffhap(iqtl))%end
               xinc2(:,d)=xinc(:,i)
               d=d+1
            end do
          end do
         ! print *,d," hap CheckIndice:",checkIndice
         end if

         if ( workstruct%effects%qtl > 0 ) then
          do iqtl=1,workstruct%nqtl
            scalei=workstruct%listnteff(iqtl)
            call change_qtleffect(workstruct%listnteff(iqtl),iqtl,&
              w_cuda%lchr(nlin,iqtl),w_cuda%lnpos(nlin,iqtl),xinc,incidenceDesc,0,.true.)
             ! call debug_write_incidence(xinc,incidenceDesc)
             if (.not. disable_sire_qtl) then
              checkIndice = checkIndice + incidenceDesc%desc(scalei)%end - incidenceDesc%desc(scalei)%start + 1
              do i=incidenceDesc%desc(scalei)%start,incidenceDesc%desc(scalei)%end
               xinc2(:,d)=xinc(:,i)
               d=d+1
              end do
              scalei = scalei +1
            end if

            if ( count(w_cuda%data%phenoAnimal%estime(incidenceDesc%ic,:))>0) then
             checkIndice = checkIndice + incidenceDesc%desc(scalei)%end - incidenceDesc%desc(scalei)%start + 1
             do i=incidenceDesc%desc(scalei)%start,incidenceDesc%desc(scalei)%end
              xinc2(:,d)=xinc(:,i)
              d=d+1
             end do
            end if
           ! print *,d," qtl CheckIndice:",checkIndice
          end do
          if ( workstruct%nteffinter /= -1 ) then
             call change_interaction_effect(workstruct%nteffinter,w_cuda%lchr(nlin,1),w_cuda%lchr(nlin,2),w_cuda%lnpos(nlin,1),&
               w_cuda%lnpos(nlin,2),xinc,incidenceDesc,workstruct%linear)
          end if

         end if

        ! Le nombre de niveaux variable à la position doit etre constant pour chaque position pour etre utiliser sur un GPU
        ! => probleme avec l'analyse LD => on fixe le nombre de niveau pour les haplotype sachant que l'etape cholesky va enlever les niveaux qui servent a rien
        if ( checkIndice /= w_cuda%nLevelVar ) then
          call stop_application("Error Dev : indice at the position "//trim(str(nlin))//" has a bad number of variates ["&
             //trim(str(checkIndice))//"] against ["//trim(str(w_cuda%nLevelVar))//"]")
        end if
       end do

       if ( workstruct%effects%haplo > 0 ) call shp%free()

      end subroutine get_partial_xinc

!!****f* m_qtlmap_incidence/gen_opti_nqtl_cuda
!!  NAME
!!    gen_opti_nqtl
!!  DESCRIPTION
!!    scan the genome under the hypothesis nqtl
!!  INPUTS
!!    nqtl            : the hypothesis
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl-1)
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!    xinc            : the contingence matrix
!!  OUTPUTS
!!    curPosMax          : maximum position finded
!!    workstruct         : residual variance and maximum Likelihood under H(nqtl)
!!    incidenceDescMax   : the description of the contingence matrix at the maximum position
!!    sigsquareEstimeMax : the residual variance
!!    BestimMax          : the solution
!!    lrtsol             : LRT : curves and maximum under H(nqtl)
!!  NOTES
!!
!!  SOURCE
     subroutine gen_opti_nqtl_cuda(      data,&
                                         spt,&
                                         nsim,      &
                                         nqtl,      &
                                    curPosMax,    &
                                    workstruct,&
                                    xinc,      &
                                    incidenceDescMax, &
                                    sigsquareEstimeMax, &
                                    BestimMax,     &
                                    lrtsol,&
                                    ySIMUL,&
                                    sizeF, &
                                    AllsigQtlMoinUn,&
                                    sizeNqtl,&
                                    listLrtSolSim)
          type(QTLMAP_DATASET)          ,target         ,intent(in) :: data
          type(PDD_BUILD)               ,target         ,intent(in) :: spt
          integer                                 , intent(in)    :: nsim,nqtl,sizeNqtl  ! under hypothesis nqtl
          type(POSITION_LRT_INCIDENCE)             ,intent(inout) :: curPosMax
          type(INCIDENCE_GEN_STRUCT)   ,target   ,intent(inout)  :: workstruct
          type(INCIDENCE_TYPE) , intent(inout)    ,target                :: incidenceDescMax  ! description of the contingence matrix
          real (kind=dp) , dimension(data%genea%nd,incidenceDescMax%ntnivmax) ,target, intent(inout) :: xinc  ! contingence matrix
          real (kind=dp) , dimension(data%genea%np) , intent(inout)          :: sigsquareEstimeMax      ! estimation of parameter on the maximum finded
          real (kind=dp) , dimension(incidenceDescMax%ntnivmax) , intent(inout)    :: BestimMax      ! estimation of parameter on the maximum finded
          type(TYPE_LRT_SOLUTION)               ,intent(inout)    :: lrtsol
          type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(data%phenoModel%ncar,sizeNqtl,nsim)    :: listLrtSolSim
          real (kind=dp), dimension (nsim+1,data%genea%nd)                       , intent(in) :: ySIMUL
          integer ,dimension(data%genea%np)                         , intent(in)  :: sizeF
          real (kind=dp), dimension (nsim+1,data%genea%np)   , intent(inout) :: AllsigQtlMoinUn


          real (kind=dp) , dimension(incidenceDescMax%ntnivmax)        :: Bestim
          real (kind=dp) , dimension(data%genea%np)              :: sigsquareEstime

          integer :: iip,ifem,ii,chr,ip,iqtl,i,nlinReal,nsim2,ieff,stat_ok,j,iq,iq2,jjm
          logical hypotheseOne, hypotheseTwo,ok
          real (kind=dp) ,dimension(:,:),pointer      :: xxx

          integer :: nGL,nGLTotal,nGLFact,InGL,nlinMax,nlin,nmax(nqtl),chmax(nqtl),nlinValide
          type(WORK_CUDA)                              :: w_cuda
          real(kind=dp) , dimension(:,:) ,allocatable  :: lrt_save,bestim_save,sigsq_save
          type(POSITION_LRT_INCIDENCE)                 :: curPos
          integer , dimension(:),allocatable ,target   :: corrLevel

      integer                                         :: isim,c,d,u,nteffS,nteffE,nkd,jm,kd,ikd
      integer                                         :: n,kd1,kd2,checkIndice,mode
      integer                  ,dimension(:)        ,pointer      :: allmaxpos
      real (kind=dp)                                 :: seuil_cholesky
      type(MAP_BASE)        ,pointer                 :: map
      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

#ifdef CUDA_SP
      real                     ,dimension(:,:,:)    ,pointer      :: allLrtSires
      real                     ,dimension(:,:)      ,pointer      :: allmaxlrt
      real                     ,dimension(:,:,:)    ,pointer      :: allbestim,allosig

      real                     ,dimension (:,:)     ,pointer      :: AllsigQtlMoinUnR
      real                     ,dimension (:,:)     ,pointer      :: ySIMULR,Mr
      real                     ,dimension (:)       ,pointer      :: CDR
#else
      real (kind=dp)           ,dimension(:,:,:,:)  ,pointer      :: allLrtSires
      real (kind=dp)           ,dimension(:,:)      ,pointer      :: allmaxlrt
      real (kind=dp)           ,dimension(:,:,:)    ,pointer      :: allbestim,allosig
#endif

       call log_mess("module incidence : gen_opti_nqtl_cuda",DEBUG_DEF)

       map => data%map

       hypotheseOne = ( lrtsol%hypothesis == 1 )

       nGL = 0

       !Nombre de point sur le groupe de liaison
       !on cherche le nombre de point dependant de l hyp nqtl du numbre de point sur le groupe de liason

          do chr=1,map%nchr
              nGL=nGL + map%get_npo(chr)
          end do
          
          nGLTotal=1
          do iqtl=1,nqtl
            nGLTotal=nGLTotal*nGL
          end do

          curPosMax%listN(nqtl)=0
          curPosMax%listN(1:nqtl-1)=1
          curPosMax%listChr=1

          ! structure pour transporter l info du code C cuda a Fortran
          allocate (w_cuda%lchr(nGLTotal,nqtl))
          allocate (w_cuda%lnpos(nGLTotal,nqtl))
          w_cuda%workstruct=>workstruct
          w_cuda%xinc=>xinc
          w_cuda%incidenceDesc=>incidenceDescMax
          w_cuda%data => data
          w_cuda%spt => spt


          nlinValide=0
          do nlin=1,nGLTotal
             i=nqtl
             ok=.false.
             do while ( (.not. ok) .and. (i>=1))
               curPosMax%listN(i) = curPosMax%listN(i)+1
               if ( curPosMax%listN(i) > map%get_npo(curPosMax%listChr(i)) ) then
                  if ( map%nchr > curPosMax%listChr(i) ) then
                    curPosMax%listN(i)=1
                    curPosMax%listChr(i)=curPosMax%listChr(i)+1
                    ok=.true.
                  else
                    curPosMax%listN(i)=1
                    curPosMax%listChr(i)=1
                    i=i-1
                  end if
               else
                  ok = .true.
               end if
             end do
             ok=.true.

             do iqtl=2,nqtl
               ok = ok .and. ( curPosMax%listN(iqtl-1)<=curPosMax%listN(iqtl) )
             end do

             if ( workstruct%effects%haplo > 0 ) then
                 ok = ok .and. (&
                        .not. (map%absi(curPosMax%listChr(1),curPosMax%listN(1))<=&
                                 map%posi(curPosMax%listChr(1),ceiling(data%params%longhap/2.d0))) ) &
                       .and. &
                        .not. (map%absi(curPosMax%listChr(nqtl),curPosMax%listN(nqtl))>&
                                 map%posi(curPosMax%listChr(nqtl),&
                                  map%nmk(curPosMax%listChr(nqtl))-ceiling(data%params%longhap/2.d0)+1))
             end if

             if ( ok ) then
              nlinValide=nlinValide+1
              w_cuda%lchr(nlinValide,:)=curPosMax%listChr
              w_cuda%lnpos(nlinValide,:)=curPosMax%listN
             end if
          end do

          nGLTotal=nlinValide


       nsim2 = nsim + 1

 !      allocate (allvecsol(nGLTotal,incidenceDesc%ntniv))
       allocate (allbestim(nGLTotal,nsim2,incidenceDescMax%ntniv),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 1/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif
       allocate (allosig(nGLTotal,nsim2,data%genea%np),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 2/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif
       allocate (allLrtSires(nGLTotal,nsim2,nqtl,data%genea%np),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 3/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif

       allocate (allmaxlrt(nsim2,nqtl),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 4/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif

       allocate (allmaxpos(nsim2),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 5/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif
       allocate (corrLevel(incidenceDescMax%ntniv),stat=stat_ok)
       if ( stat_ok /= 0 ) then
          call stop_application("** 6/6 -> Not enough host memory ** number of position to test:"//trim(str(nGLTotal)))
       endif

       w_cuda%corrLevel=>corrLevel

       u=1
       w_cuda%nLevelFix = 0
       !les effets fixe a la position sont en premier
       do ieff=1,incidenceDescMax%nteff
         if ( incidenceDescMax%desc(ieff)%isVar ) cycle
        ! print *,"fix:",ieff,incidenceDescMax%desc(ieff)%start,incidenceDescMax%desc(ieff)%end
         do i=incidenceDescMax%desc(ieff)%start,incidenceDescMax%desc(ieff)%end
             corrLevel(u)=i-1
             u=u+1
      !       print *,i,"->",corrLevel(i)
             w_cuda%nLevelFix = w_cuda%nLevelFix + 1
         end do
       end do

       w_cuda%nLevelVar = 0
        !puis les effets variables
       do ieff=1,incidenceDescMax%nteff
         if ( .not. incidenceDescMax%desc(ieff)%isVar ) cycle
        !  print *,"VAR:",ieff,incidenceDescMax%desc(ieff)%start,incidenceDescMax%desc(ieff)%end
         do i=incidenceDescMax%desc(ieff)%start,incidenceDescMax%desc(ieff)%end
             corrLevel(u)=i-1
             u=u+1
             w_cuda%nLevelVar = w_cuda%nLevelVar + 1
         !    print *,i,"->",corrLevel(i)
         end do
       end do

      call log_mess("call CUDA analysis.....",INFO_DEF)

      mode = 0

      if (is_parameter(RESOLUTION_LINEAR,workstruct%type_model) .and.&
          is_parameter(KIND_HETEROSCEDASTIC,workstruct%type_model).and. &
          is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) mode=1

      if (is_parameter(RESOLUTION_LINEAR,workstruct%type_model) .and.&
          is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and. &
          is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) mode=2



#ifdef CUDA_SP

     allocate (AllsigQtlMoinUnR(nsim2,data%genea%np))
     AllsigQtlMoinUnR = AllsigQtlMoinUn
     allocate (ySIMULR(nsim2,incidenceDescMax%dataset%nkd))
     allocate (CDR(incidenceDescMax%dataset%nkd))

     if ( is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) then
      allocate (Mr(size(workstruct%M,1),size(workstruct%M,2)))
      Mr=workstruct%M
     else
      allocate (Mr(1,1))
     end if

     CDR=incidenceDescMax%dataset%CD(1,:)
     ySIMULR=ySIMUL(:,:incidenceDescMax%dataset%nkd)
     seuil_cholesky=0.3


     call log_mess("** simple precision ** ",INFO_DEF)
     call cuda_model_resolv_genome(data%params%gpu_device_id,      &
                                      mode,                        &
                                      nqtl,                        &
                                      AllsigQtlMoinUnR,            &
                                      w_cuda%nLevelFix,            &
                                      w_cuda%nLevelVar,            &
                                      corrLevel,                   &
                                      w_cuda,                      &
                                      get_xinc_fixed_effect,       &
                                      get_partial_xinc,            &
                                      ySIMULR,                      &
                                      CDR,                          &
                                      Mr,                           &
                                      nsim2,                        &
                                      data%genea%nd,                &
                                      incidenceDescMax%dataset%nkd,   &
                                      data%genea%np,                  &
                                      sizeF,                       &
                                      incidenceDescMax%ntnivmax,      &
                                      nGLTotal,                    &
                                      seuil_cholesky,              &
                                      allbestim,                   &
                                      allosig,                     &
                                      allLrtSires,                 &
                                      allmaxlrt,                   &
                                      allmaxpos )

      AllsigQtlMoinUn = AllsigQtlMoinUnR
      deallocate (AllsigQtlMoinUnR)
      deallocate (ySIMULR)
      deallocate (CDR)
      deallocate (Mr)


#else
     seuil_cholesky=data%params%SEUIL_CHO
     call log_mess("** double precision ** ",INFO_DEF)
     call cuda_model_resolv_genome(data%params%gpu_device_id,      &
                                      mode,                        &
                                      nqtl,                        &
                                      AllsigQtlMoinUn,             &
                                      w_cuda%nLevelFix,            &
                                      w_cuda%nLevelVar,            &
                                      corrLevel,                   &
                                      w_cuda,                      &
                                      get_xinc_fixed_effect,       &
                                      get_partial_xinc,            &
                                      ySIMUL,                      &
                                      incidenceDescMax%dataset%CD,    &
                                      workstruct%M,                    &
                                      nsim2,                        &
                                      data%genea%nd,                &
                                      incidenceDescMax%dataset%nkd, &
                                      data%genea%np,                &
                                      sizeF,                       &
                                      incidenceDescMax%ntnivmax,      &
                                      nGLTotal,                    &
                                      seuil_cholesky,              &
                                      allbestim,                   &
                                      allosig,                     &
                                      allLrtSires,                 &
                                      allmaxlrt,                   &
                                      allmaxpos )
#endif

!       print *,"***************LRT***********"
!       do ip=1,data%genea%np
!         print *,"ip=",ip
!         print *,allLrtSires(:,1,1,ip)
!       end do

       isim=1

       do nlin=1,nGLTotal
         do iq=1,workstruct%nqtl
          call lrtsol%LRT%add(data,workstruct%nqtl,w_cuda%lchr(nlin,:),w_cuda%lnpos(nlin,:),&
                              sum(allLrtSires(nlin,isim,iq,:)),iq)
          do ip=1,data%genea%np
           call lrtsol%LRT_SIRES(ip)%add(data,workstruct%nqtl,w_cuda%lchr(nlin,:),w_cuda%lnpos(nlin,:),&
                                      allLrtSires(nlin,isim,iq,ip),iq)
          end do
          do jm=1,data%genea%nm
           call lrtsol%LRT_DAMS(jm)%add(data,workstruct%nqtl,w_cuda%lchr(nlin,:),w_cuda%lnpos(nlin,:),&
                                  0.d0,iq)
          end do
         end do
       end do


!       if ( workstruct%nqtl == 1 ) then
!         do nlin=1,nGLTotal
!          lrtsol%lrt1(w_cuda%lchr(nlin,1),w_cuda%lnpos(nlin,1))=sum(allLrtSires(nlin,isim,1,:))
!          lrtsol%xlrp(w_cuda%lchr(nlin,1),:,w_cuda%lnpos(nlin,1))=allLrtSires(nlin,isim,1,:)
!         end do
!       end if
!       if ( workstruct%nqtl == 2 ) then
!         lrtsol%lrt0_2 = 0.d0 ! pas implemente....
!         do nlin=1,nGLTotal
!          lrtsol%lrt1_2(w_cuda%lchr(nlin,1),w_cuda%lchr(nlin,2),&
!           w_cuda%lnpos(nlin,1),w_cuda%lnpos(nlin,2))=sum(allLrtSires(nlin,isim,2,:))
!          lrtsol%xlrp2(w_cuda%lchr(nlin,1),w_cuda%lchr(nlin,2),:,&
!           w_cuda%lnpos(nlin,1),w_cuda%lnpos(nlin,2))=allLrtSires(nlin,isim,2,:)
!         end do
!       end if

       if ( allmaxpos(isim) > 0 ) then
        sigsquareEstimeMax=allosig(allmaxpos(isim),isim,:)
        incidenceDescMax%vecsol=.false.

        j=1
        do i=1,incidenceDescMax%ntniv
          if ( allbestim(allmaxpos(isim),isim,i) /= 0.d0 ) then
               BestimMax(j) = allbestim(allmaxpos(isim),isim,i)
               incidenceDescMax%vecsol(i)=.true.
               j=j+1
          end if
        end do

        AllsigQtlMoinUn(isim,:) = allosig(allmaxpos(isim),isim,:)

        do iqtl=1,nqtl
         lrtsol%nxmax(iqtl-1)=w_cuda%lnpos(allmaxpos(isim),iqtl)
         lrtsol%chrmax(iqtl-1)=w_cuda%lchr(allmaxpos(isim),iqtl)
         lrtsol%lrtmax(iqtl-1)=allmaxlrt(isim,iqtl)
         curPosMax%listN=w_cuda%lnpos(allmaxpos(isim),iqtl)
         curPosMax%listChr=w_cuda%lchr(allmaxpos(isim),iqtl)
        end do
!        print *,lrtsol%lrtmax
       end if

       do isim=1,nsim

          AllsigQtlMoinUn(isim+1,:) = allosig(allmaxpos(isim+1),isim+1,:)
          allocate (listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%lrtmax(0:nqtl-1))
          allocate (listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%nxmax(0:nqtl-1))
          allocate (listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%chrmax(0:nqtl-1))

          listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%nxmax=0
          listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%chrmax=0
          listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%lrtmax=0

          if ( allmaxpos(isim+1) <= 0 .or. allmaxpos(isim+1) > nGLTotal) cycle

          do iqtl=1,nqtl
            listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%nxmax(iqtl-1)=w_cuda%lnpos(allmaxpos(isim+1),iqtl)
            listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%chrmax(iqtl-1)=w_cuda%lchr(allmaxpos(isim+1),iqtl)
            listLrtSolSim(incidenceDescMax%ic,nqtl,isim)%lrtmax(iqtl-1)=allmaxlrt(isim+1,iqtl)
         end do
       end do

       !TODO : Modifier les structure pateff et mateff pour ajouter les effet qtls
       if ( workstruct%nqtl == 1 ) then
        lrtsol%pater_eff=0.d0
        lrtsol%mater_eff=0.d0
        !if ( workstruct%effects%qtl == 0 ) cycle ! pas d'effet qtl...

        do iq=1,workstruct%nqtl
         !recuperation des effets paternels
         iip = incidenceDescMax%ntniv_qtlsires(iq)
         if ( iip == 0 ) cycle
         do ip=1,data%genea%np
          do nlin=1,nGLTotal
           do iq2=1,iq
            if ( incidenceDescMax%vecsol(iip) ) then
             lrtsol%pater_eff(w_cuda%lchr(nlin,iq2),ip,w_cuda%lnpos(nlin,iq2))=allbestim(nlin,1,iip)
            end if
           end do ! iq
          end do ! nlin
          iip = iip + 1
         end do !ip
         !recuperation des effets maternels
         jjm = incidenceDescMax%ntniv_qtldams(iq)
         if ( jjm == 0 ) cycle
         ifem=0
         do jm=1,data%genea%nm
           if ( data%phenoAnimal%estime(incidenceDescMax%ic,jm)) then
              ifem = ifem + 1
              do nlin=1,nGLTotal
                do iq2=1,iq
                 if ( incidenceDescMax%vecsol(jjm) ) then
                   lrtsol%mater_eff(w_cuda%lchr(nlin,iq2),ifem,w_cuda%lnpos(nlin,iq2))=allbestim(nlin,1,jjm)
                 end if
                end do ! iq
              end do ! nlin
             jjm = jjm+1
            end if ! estime
          end do ! jm
        end do
       end if


       deallocate (corrLevel)
!       deallocate (xinc2)
!       deallocate (allvecsol)
       deallocate (allbestim)
       deallocate (allosig)
       deallocate (w_cuda%lchr)
       deallocate (w_cuda%lnpos)
       deallocate (allLrtSires)
       deallocate (allmaxlrt)
       deallocate (allmaxpos)

     end subroutine gen_opti_nqtl_cuda
!!***
#endif

!!****f* m_qtlmap_incidence/opti_nqtl
!!  NAME
!!    opti_nqtl
!!  DESCRIPTION
!!    compute the LRT under H(nqtl) 1 trait    nqtl>=1
!!  INPUTS
!!    nqtl            : the hypothesis
!!    ic              : index of the trait
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl-1)
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!  OUTPUTS
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl)
!!    incsol          : the solution of the estimation at the maximum finded
!!    lrtsol          : LRT : curves and maximum under H(nqtl)
!!  NOTES
!!    use the iterative scan gen_opti_nqtl function
!!  SOURCE
     subroutine opti_nqtl(      dataset,   &
                                spt,       &
                                nqtl,      &  ! under hypothesis NQTL
                                hyp,       &
                                        ic,       &  ! the trait
                                workStruct,       &  ! variance / likelihood estimated under NQTL-1 hypothesis
                                 incsol   ,       &  ! incidence solution (estimation of each effect)
                                 lrtsol   ,       &  ! maximum lrt finding at a position
                               FUNCT_MODEL         ) ! function with a model specific

      type(QTLMAP_DATASET)                    ,intent(inout)       :: dataset
      type(PDD_BUILD)                         ,intent(inout)       :: spt
      integer                                 ,intent(in)          :: nqtl,ic,hyp
      type(INCIDENCE_GEN_STRUCT)              ,intent(inout)       :: workstruct
      type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)         :: incsol
      type(TYPE_LRT_SOLUTION)                 ,intent(out)         :: lrtsol
      external                                                     :: FUNCT_MODEL


      integer :: i

      real (kind=dp) ,dimension(:,:),pointer      :: xinc ,xxx    ! incidence matrix and temp for testing nuisance effect
      real (kind=dp) , dimension(:) ,pointer      :: Bestim       ! solution of the estimation
      real (kind=dp) , dimension(dataset%genea%np):: sigsquareEstime
      type(INCIDENCE_TYPE)                        :: incidenceDesc ! description incidence matrix
      type(POSITION_LRT_INCIDENCE)                :: curPos        ! the position with additional info for the model
      integer :: n,j,chr
      type(HAPLOTYPE_POSITION_BUILD)   :: shp

      call log_mess("module incidence : opti_nqtl",DEBUG_DEF)

      if ( nqtl <= 0 ) then
        call stop_application("Error dev: can not call opti_nqtl_linear with nqtl:"//trim(str(nqtl)))
      end if

      !Position Allocation
      call init_position (dataset,hyp,nqtl,curPos)

      !initialisation of incidence matrix
      call init_incidence(dataset,spt,ic,nqtl,incidenceDesc,workstruct)

      allocate (xinc(dataset%genea%nd,incidenceDesc%ntnivmax))
      xinc=0.d0
      allocate (Bestim(incidenceDesc%ntnivmax))
      allocate (xxx(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax))

      Bestim=0.d0

      call build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,0,0,0,0)

      call FUNCT_MODEL(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.false.,.false.,xxx,.false.)

      !initalizing lrtsolution
      call lrtsol%new(incidenceDesc%dataset%data,hyp,nqtl)

      do i=0,nqtl-1
        lrtsol%lrtmax(i)=sum(curPos%lrtSires(i+1,:))
      end do

      lrtsol%chrmax(0:nqtl-1)=curPos%listChr
      lrtsol%nxmax(0:nqtl-1)=curPos%listN

      call gen_opti_nqtl(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,nqtl,hyp,curPos,workstruct,xinc,&
       incidenceDesc,sigsquareEstime,Bestim,lrtsol,FUNCT_MODEL)


      if ( workstruct%effects%haplo>0 ) then
       if ( nqtl > 0 ) then
        call shp%set(dataset,spt)
        do i=1,nqtl
         call add_haplotype_effect(workstruct%listnteffhap(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),&
              xinc,incidenceDesc,shp,workstruct%hdam,.false.)
        end do
        call shp%free()
       end if
      end if

      if (  workstruct%effects%qtl>0 ) then
      ! Compute precision at the maximum LRT in position posx
       do i=1,nqtl
        call change_qtleffect(workstruct%listnteff(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),&
          xinc,incidenceDesc,0,workstruct%linear)
        call change_name_qtleffect(workstruct%listnteff(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),incidenceDesc)
       end do

       if ( workstruct%nteffinter /= -1 ) then
            call change_interaction_effect(workstruct%nteffinter,lrtsol%chrmax(0),lrtsol%chrmax(1),lrtsol%nxmax(0),&
               lrtsol%nxmax(1),xinc,incidenceDesc,workstruct%linear)
       end if
      end if

      do i=1,nqtl
       curPos%listChr(i)=lrtsol%chrmax(i-1)
       curPos%listN(i)=lrtsol%nxmax(i-1)
      end do

      !*********************************************
      ! DEBUG
      !*********************************************
      !call debug_write_incidence(xinc,incidenceDesc)

      call FUNCT_MODEL(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.true.,workstruct%performConfusion,xxx,.false.)

      if ( workstruct%performConfusion) then
         call confusion_type1(hyp,incidenceDesc,xxx,workstruct)
       end if

      if ( workstruct%performTestNuis ) then
        call test_nuisances(dataset,spt,ic,incidenceDesc,lrtsol,curPos,workstruct,sigsquareEstime,FUNCT_MODEL)
      end if

      call set_solution(nqtl,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,nqtl,workstruct%listnteff)

      workstruct%fnqtlsires(hyp,:)=incidenceDesc%fperemax
      workstruct%fnqtldams(hyp,:)=incidenceDesc%fmeremax
      workstruct%sigsquare(hyp,:,1)=sigsquareEstime

      call end_incidence(incidenceDesc)
      deallocate (xinc)
      deallocate (xxx)
      deallocate (Bestim)
      call end_position(curPos)

      end subroutine opti_nqtl
!!***

!!****f* m_qtlmap_incidence/gen_opti_nqtl
!!  NAME
!!    gen_opti_nqtl
!!  DESCRIPTION
!!    scan the genome under the hypothesis nqtl
!!  INPUTS
!!    nqtl            : the hypothesis
!!    workstruct      : residual variance and maximum Likelihood under H(nqtl-1)
!!    FUNCT_MODEL     : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!    xinc            : the contingence matrix
!!  OUTPUTS
!!    curPosMax          : maximum position finded
!!    workstruct         : residual variance and maximum Likelihood under H(nqtl)
!!    incidenceDescMax   : the description of the contingence matrix at the maximum position
!!    sigsquareEstimeMax : the residual variance
!!    BestimMax          : the solution
!!    lrtsol             : LRT : curves and maximum under H(nqtl)
!!  NOTES
!!
!!  SOURCE
     subroutine gen_opti_nqtl(      data,         &
                                    spt,          &
                                    nqtl,         &
                                    hyp,          &
                                    curPosMax,    &
                                    workstruct,&
                                    xincOrig,      &
                                    incidenceDescMax, &
                                    sigsquareEstimeMax, &
                                    BestimMax,     &
                                    lrtsol,        &
                                    FUNCT_MODEL)
          type(QTLMAP_DATASET)                    ,intent(in)     :: data
          type(PDD_BUILD)                         ,intent(in)     :: spt
          integer                                 , intent(in)    :: nqtl,hyp  ! under hypothesis nqtl
          type(POSITION_LRT_INCIDENCE)             ,intent(inout) :: curPosMax
          type(INCIDENCE_GEN_STRUCT)              ,intent(inout)  :: workstruct
          type(INCIDENCE_TYPE) , intent(inout)                    :: incidenceDescMax  ! description of the contingence matrix
          real (kind=dp) , dimension(data%genea%nd,incidenceDescMax%ntnivmax), target, intent(inout) :: xincOrig  ! contingence matrix
          real (kind=dp) , dimension(data%genea%np) , intent(inout)  :: sigsquareEstimeMax      ! estimation of parameter on the maximum finded
          real (kind=dp) , dimension(incidenceDescMax%ntnivmax) , intent(inout)    :: BestimMax      ! estimation of parameter on the maximum finded
          type(TYPE_LRT_SOLUTION)               ,intent(inout)    :: lrtsol
          external                                                :: FUNCT_MODEL


          real (kind=dp) , dimension(incidenceDescMax%ntnivmax)        :: Bestim
          real (kind=dp) , dimension(:,:) , pointer :: xincSub
          real (kind=dp) , dimension(data%genea%np)              :: sigsquareEstime

          integer :: iip,ifem,ii,chr,ip,iqtl,i,nlinReal,status
          logical hypotheseOne, ok
          real (kind=dp) ,dimension(:,:),pointer      :: xxx

          integer :: nGL,nGLTotal,nGLFact,InGL,nlinMax,nlin,nmax(nqtl),chmax(nqtl),nlinValide,ipos
          integer       , dimension(:,:) ,pointer :: lchr,lnpos
          real(kind=dp) , dimension(:,:) ,allocatable :: lrt_save,bestim_save,sigsq_save
          type(POSITION_LRT_INCIDENCE)             :: curPos
          type(INCIDENCE_TYPE)                     :: incidenceDesc
          logical                                  :: loopFalse,isLD
         ! logical , dimension(:) ,allocatable      :: bufvecsol
          type(HAPLOTYPE_POSITION_BUILD)           :: shp

          call log_mess("module incidence : gen_opti_nqtl",DEBUG_DEF)
          ! G.Salle
          !type(TYPE_INCIDENCE_SOLUTION)            :: tempincsol

         !   print *,"gen_opti_nqtl*******************"
          hypotheseOne = ( lrtsol%hypothesis == 1 )

          !nGL = 0

          call get_iterator_index(data,nqtl,lchr,lnpos,nGLTotal)

          allocate (bestim_save(nGLTotal,incidenceDescMax%ntnivmax))
          allocate (sigsq_save(nGLTotal,data%genea%np))
          allocate (lrt_save(nGLTotal,hyp))

          lrt_save=-INIFINY_REAL_VALUE

          isLD = workstruct%effects%haplo>0

          !call debug_write_incidence(xincSub,incidenceDescMax)

          !$OMP PARALLEL DEFAULT(SHARED)  &
!!!          !$OMP PRIVATE(xincSub,curPos,i,iqtl,ok,sigsquareEstime,Bestim,xxx,ip,iip,ifem,incidenceDesc,loopFalse,shp,bufvecsol)
          !$OMP PRIVATE(xincSub,curPos,i,iqtl,ok,sigsquareEstime,Bestim,xxx,ip,iip,ifem,incidenceDesc,loopFalse,shp)

          allocate (xincSub(data%genea%nd,incidenceDescMax%ntnivmax),stat=status)
!          allocate (bufvecsol(incidenceDescMax%ntnivmax))
          xincSub=xincOrig
          if ( status /= 0 ) call stop_application("Can not allocate incidence matrix ** not enough memory ** ")
          call init_position (data,hyp,nqtl,curPos)
          allocate (xxx(incidenceDescMax%ntnivmax,incidenceDescMax%ntnivmax),stat=status)
          if ( status /= 0 ) call stop_application("Can not allocate xxx array ** not enough memory ** ")
          !initialisation of contingence matrix
          call copy_incidence_desc (incidenceDescMax,incidenceDesc)

          curPos%lrtSires(nqtl,:)=0
          curPos%lrtDams(nqtl,:)=0
          curPos%listN(nqtl)=0
          curPos%listN(1:nqtl-1)=1
          curPos%listChr=1

          if (workstruct%effects%haplo > 0) call shp%set(data,spt)

          !$OMP DO
          do nlin=1,nGLTotal
             call log_mess( trim(str((float(nlin)/float(nGLTotal))*100.d0))//"%", DEBUG_DEF )
             curPos%listN=lnpos(nlin,:)
             curPos%listChr=lchr(nlin,:)

             !traitement specifique pour le LD
             if ( workstruct%effects%haplo > 0 ) then
               i=1
               loopFalse = .false.
               do while (i<=nqtl .and. .not. loopFalse )
               if (data%map%absi(curPos%listChr(i),curPos%listN(i))<=&
                data%map%posi(curPos%listChr(i),ceiling(data%params%longhap/2.d0)) ) &
                 loopFalse = .true.
               if (data%map%absi(curPos%listChr(i),curPos%listN(i))>=&
                data%map%posi(curPos%listChr(i),data%map%nmk(curPos%listChr(i))-ceiling(data%params%longhap/2.d0)+1)) &
                 loopFalse = .true.
                i = i + 1
               end do
               if ( loopFalse ) cycle
               if ( isLD ) then
                do i=1,nqtl
                 call add_haplotype_effect(workstruct%listnteffhap(i),i,curPos%listChr(i),curPos%listN(i),&
                  xincSub,incidenceDesc,shp,workstruct%hdam,.false.)
                end do
               end if
             end if

              if ( workstruct%effects%qtl > 0 ) then
               do iqtl=1,nqtl
                call change_qtleffect(workstruct%listnteff(iqtl),iqtl,curPos%listChr(iqtl),curPos%listN(iqtl),&
                                     xincSub,incidenceDesc,0,workstruct%linear)
               end do

                if ( workstruct%nteffinter /= -1 ) then
                 call change_interaction_effect(workstruct%nteffinter,curPos%listChr(1),curPos%listChr(2),curPos%listN(1),&
               curPos%listN(2),xincSub,incidenceDesc,workstruct%linear)
                end if
              end if

              !if (nqtl>=1) call debug_write_incidence(xincSub,incidenceDesc)
!              if ( hypotheseTwo ) then
!                if ( curPos%listN(1)==4 .and. curPos%listN(2)==19 ) then
!                   call debug_write_incidence(xincSub,incidenceDesc)
!                end if
!              end if
              ! Need to phase pdd ?
!              bufvecsol=incidenceDesc%vecsol
!              do iqtl=1,nqtl
!               if ( incidenceDesc%isBiAllellic ) then
              !  print *,"ntniv qtl=",incidenceDesc%desc(workstruct%listnteff(iqtl))%start
!                 bufvecsol(incidenceDesc%desc(workstruct%listnteff(iqtl))%start)=.false.
!              else
!                 bufvecsol(incidenceDesc%desc(workstruct%listnteff(iqtl))%start:&
!                  incidenceDesc%desc(workstruct%listnteff(iqtl))%start+data%genea%np)=.false.
!               end if
!              end do

              ! NOTE OFI...
              ! Normalement on devrait le faire qu'une seule fois
              ! Sous 2QTL less pdd ont deja été inversé sous HYP 1 QTL,...a voir pour deplacer ce traitement
              if ( incidenceDesc%isBiAllellic ) then
              do ipos=1,nqtl
                if ( curPos%listN(ipos) > 1 ) then
                 call set_phase_pdd(incidenceDesc%dataset%data,incidenceDesc%dataset%spt,&
                            curPos%listChr(ipos),curPos%listN(ipos),incidenceDesc%ic,&
                            incidenceDesc%dataset%data%genea%nd)!,&
!                            incidenceDesc%ntnivmax,incidenceDesc%ntniv,&
!                            xincSub,bufvecsol,bestim_save(nlin-1,:))
                end if
              end do
              end if

              call FUNCT_MODEL(xincSub,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.false.,.false.,xxx,.false.)

!              if (workstruct%effects%inter /= 0 ) then
!              if (workstruct%nqtl > 1 ) then
!                call debug_write_incidence(xincSub,incidenceDesc)
!              end if
              ! Demande G.Salle 20/09/2011
              !   call set_solution(nqtl,workstruct,sigsquareEstime,Bestim,incidenceDesc,tempincsol,nqtl,workstruct%listnteff)
              !   call debug_solution(nqtl,curPos%listChr,curPos%listN,tempincsol)
              !   call release(tempincsol)

              !save value
              bestim_save(nlin,:incidenceDesc%nbnivest)=Bestim(:incidenceDesc%nbnivest)
              sigsq_save(nlin,:data%genea%np)=sigsquareEstime(:data%genea%np)
              do i=1,hyp
                lrt_save(nlin,i)=sum(curPos%lrtSires(i,:))
                call lrtsol%LRT%add(data,nqtl,curPos%listChr,curPos%listN,lrt_save(nlin,i),i)
                do ip=1,data%genea%np
                   call lrtsol%LRT_SIRES(ip)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtSires(i,ip),i)
                end do

                do ii=1,data%genea%nm
                   call lrtsol%LRT_DAMS(ii)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtDams(i,ii),i)
                end do

              end do
          !    print *,"LRT=",lrt_save(nlin,:)!,OMP_GET_THREAD_NUM()

              if ( hypotheseOne ) then
                iip=1
                do ip=1,data%genea%np
                  if ( incidenceDesc%vecsol(1+ip)) then
                    iip=iip+1
                    ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                    lrtsol%pater_eff(curPos%listChr(1),ip,curPos%listN(1))=Bestim(iip)
                  end if
                end do
                ifem=0
                do ii=1,data%genea%nm
                  if ( data%phenoAnimal%estime(incidenceDesc%ic,ii)) then
                    ifem=ifem+1
                    if ( incidenceDesc%vecsol(data%genea%np+1+ifem) ) then
                      iip=iip+1
                      ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                      lrtsol%mater_eff(curPos%listChr(1),ifem,curPos%listN(1))=Bestim(iip)
                    end if
                  end if
                end do
              end if
        !    end if !ok
         end do ! nlin
         !$OMP END DO
         if (workstruct%effects%haplo > 0) call shp%free()
         call end_position (curPos)
         deallocate (xxx)
         call release_copy_incidence_desc(incidenceDesc)
         deallocate (xincSub)
!         deallocate (bufvecsol)
         !$OMP END PARALLEL

         nlinMax=1
         lrtsol%lrtmax=lrt_save(nlinMax,:)
         nmax=1
         chmax=1

         call init_position (data,hyp,nqtl,curPos)

         curPos%listN(nqtl)=0
         curPos%listN(1:nqtl-1)=1
         curPos%listChr=1

         do nlin=1,nGLTotal
           curPos%listN=lnpos(nlin,:)
           curPos%listChr=lchr(nlin,:)

           if(lrtsol%lrtmax(hyp-1) < lrt_save(nlin,hyp)) then
             nmax=curPos%listN
             chmax=curPos%listChr
             nlinMax=nlin
             lrtsol%lrtmax=lrt_save(nlin,:)
         !    print *,nmax
           end if
         end do

         do iqtl=1,nqtl
            lrtsol%nxmax(iqtl-1)=nmax(iqtl)
            lrtsol%chrmax(iqtl-1)=chmax(iqtl)
            BestimMax = bestim_save(nlinMax,:)
            sigsquareEstimeMax=sigsq_save(nlinMax,:)
            curPosMax%listN=curPos%listN
            curPosMax%listChr=curPos%listChr
         end do

         call end_position (curPos)

         deallocate (bestim_save)
         deallocate (sigsq_save)
         deallocate (lrt_save)
         deallocate (lchr)
         deallocate (lnpos)

     end subroutine gen_opti_nqtl
!!***

!!****f* m_qtlmap_incidence/test_nuisances
!!  NAME
!!    test_nuisances
!!  DESCRIPTION
!!    Test des differents effets de nuisance du modele par une LRT compare a une chi2
!!  INPUTS
!!    ic                     : index of the trait
!!    incidenceDescOrigine   : description of the contingence matrix at the maximum
!!    curPosMax              : maximum position
!!    sigsquareEstime        : the residual variance at the maximum
!!    FUNCT_MODEL            : the function for the resolution (model_lin_h0,model_lin_hn,model_optim_h0,model_optim_hn)
!!  OUTPUTS
!!    workstruct             : residual variance and maximum Likelihood under H(nqtl)
!!  NOTES
!!
!!  SOURCE
      subroutine test_nuisances(data,spt,ic,incidenceDescOrigine,lrtsol,curPosMax,workstruct,sigsquareEstime,FUNCT_MODEL)
        type(QTLMAP_DATASET)                    ,intent(inout)       :: data
        type(PDD_BUILD)                         ,intent(inout)       :: spt
        integer                                 ,intent(in)          :: ic!,nbnivestmax
        type(INCIDENCE_TYPE)                    ,intent(in)          :: incidenceDescOrigine
        type(TYPE_LRT_SOLUTION)                 ,intent(in)          :: lrtsol
        type(POSITION_LRT_INCIDENCE)            ,intent(in)          :: curPosMax
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)       :: workstruct
        real (kind=dp)      , dimension(data%genea%np) ,intent(in)   :: sigsquareEstime
        external                                                     :: FUNCT_MODEL

        real (kind=dp) ,dimension(:,:),pointer      :: xinc,xxx
        real (kind=dp) , dimension(:) ,pointer      :: Bestim,SIGest
        type(INCIDENCE_TYPE)                        :: incidenceDesc
        real (kind=dp) :: prob,sigEstSave(data%genea%np),fperemaxSave(data%genea%np),fmeremaxSave(data%genea%nm)
        integer :: itest,i,j,nbreduit,nbef,nbco,nkd,ip,ifail,nqtl,nbnivestmax,nbeffnuis,numqtl
        type(DATAMODEL_BASE) , pointer :: dpm

        call log_mess("module incidence : test_nuisances",DEBUG_DEF)

        dpm => data%phenoModel

        nbnivestmax = incidenceDescOrigine%nbnivest
        nbef=dpm%listModelTrait(ic)%nbfe
        nbco=dpm%listModelTrait(ic)%nbco

      !  nbint=listModelTrait(incidenceDesc%ic)%nbint(numqtl)
        nqtl=workstruct%nqtl
        nbeffnuis = nbef+nbco
        do numqtl=1,nqtl
          nbeffnuis = nbeffnuis + dpm%listModelTrait(ic)%nbint(numqtl)
        end do

        if (nbeffnuis==0) then
          workstruct%ntest = 0
          return
        end if

        sigEstSave=workstruct%sigsquare(workstruct%hypothesis-1,:,1)
        fperemaxSave=workstruct%fnqtlsires(workstruct%hypothesis-1,:)
        fmeremaxSave=workstruct%fnqtldams(workstruct%hypothesis-1,:)

        workstruct%sigsquare(workstruct%hypothesis-1,:,1)=sigsquareEstime
        !workstruct%fnqtl(nqtl)=incidenceDesc%fmax

        workstruct%ntest=nbeffnuis
        allocate (workstruct%listtestnuis(workstruct%ntest))
        itest=0

        !initialisation of contingence matrix
        call init_incidence(data,spt,ic,nqtl,incidenceDesc,workstruct,incidenceDescOrigine%dataset)

        allocate (xinc(data%genea%nd,incidenceDescOrigine%ntnivmax))
        xinc=0.d0
        allocate (Bestim(incidenceDescOrigine%ntnivmax))
        allocate (SIGest(data%genea%np))
        allocate (xxx(incidenceDescOrigine%ntnivmax,incidenceDescOrigine%ntnivmax))

        !on laisse seulement les variances comme point de depart
        incidenceDesc%par(data%genea%np+1:)=0.d0

        Bestim=0.d0
       ! workStruct%

        do i=1,nbef+nbco
          if ( i <= nbef ) then
            call build_incidence_matrix(workstruct,curPosMax,xinc,incidenceDesc,i,0,0,0)
          end if
          if ( i > nbef .and. i <= nbef+nbco ) then
            call build_incidence_matrix(workstruct,curPosMax,xinc,incidenceDesc,0,i-nbef,0,0)
          end if

          !call debug_write_incidence(xinc,incidenceDesc)
          !call model analysis

          call FUNCT_MODEL(xinc,incidenceDesc,curPosMax,workstruct,SIGest,Bestim,.false.,.false.,xxx,.true.)

          ifail=0
          if ( i > nbef .and. i <= nbco ) then
            nbreduit=1
          else
            nbreduit = nbnivestmax -  incidenceDesc%nbnivest
          end if

          itest=itest+1

          if ( i <= nbef ) then
           workStruct%listtestnuis(itest)%name=&
            trim(data%phenoModel%namefix(dpm%listModelTrait(ic)%indexFixedEffect(i)))
          else
           workStruct%listtestnuis(itest)%name=&
            trim(data%phenoModel%namefix(dpm%listModelTrait(ic)%indexCovariate(i-nbef)))
          end if

          if ( nbreduit == 0 ) then
            call log_mess("The effect ["//trim(workStruct%listtestnuis(itest)%name)//"]"//&
                          " might be confused with another effect !",WARNING_DEF)
            prob=0.d0
          else
           prob=MATH_QTLMAP_G01ECF('U',sum(curPosMax%lrtSires(workstruct%hypothesis,:)),dble(nbreduit),ifail)
          end if

          workStruct%listtestnuis(itest)%directeffect=.true.
          workStruct%listtestnuis(itest)%df=nbreduit
          workStruct%listtestnuis(itest)%lrt=sum(curPosMax%lrtSires(workstruct%hypothesis,:))
          workStruct%listtestnuis(itest)%pvalue=prob

          do while(incidenceDesc%nteff>0)
            call remove_last_effect(incidenceDesc,xinc)
          end do
        end do

        do numqtl=1,nqtl
         do i=1,dpm%listModelTrait(incidenceDesc%ic)%nbint(numqtl)
            call build_incidence_matrix(workstruct,curPosMax,xinc,incidenceDesc,0,0,i,numqtl)
          !call debug_write_incidence(xinc,incidenceDesc)
          !call model analysis

          call FUNCT_MODEL(xinc,incidenceDesc,curPosMax,workstruct,SIGest,Bestim,.false.,.false.,xxx,.true.)

          ifail=0
          nbreduit = nbnivestmax -  incidenceDesc%nbnivest

          prob=MATH_QTLMAP_G01ECF('U',sum(curPosMax%lrtSires(workstruct%hypothesis,:)),dble(nbreduit),ifail)

          itest=itest+1
          workStruct%listtestnuis(itest)%directeffect=.false.
          workStruct%listtestnuis(itest)%name=trim(data%phenoModel%namefix(&
           dpm%listModelTrait(ic)%indexFixedEffectWithInteraction(numqtl,i)))//&
            " QTL["//trim(str(numqtl))//"]"
          workStruct%listtestnuis(itest)%df=nbreduit
          workStruct%listtestnuis(itest)%lrt=sum(curPosMax%lrtSires(workstruct%hypothesis,:))
          workStruct%listtestnuis(itest)%pvalue=prob

          do while(incidenceDesc%nteff>0)
            call remove_last_effect(incidenceDesc,xinc)
          end do
         end do
        end do

        workstruct%sigsquare(workstruct%hypothesis-1,:,1)=sigEstSave
        workstruct%fnqtlsires(workstruct%hypothesis-1,:)=fperemaxSave
        workstruct%fnqtldams(workstruct%hypothesis-1,:)=fmeremaxSave

        call end_incidence(incidenceDesc)

        deallocate (SIGest)
        deallocate (xinc)
        deallocate (Bestim)
        deallocate (xxx)

       end subroutine test_nuisances
!!***

end module m_qtlmap_incidence
