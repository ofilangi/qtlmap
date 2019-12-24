!!****m* DATA/m_qtlmap_constant
!!  NAME
!!    m_qtlmap_constant -- Definition of constantes
!!  SYNOPSIS
!!    All definition useb by QTLMap are defined here....
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

module m_qtlmap_constant
   implicit none



   !precision a 15 chiffre pour les reel
   !pour une constante  ajoute _dp ex:real, public, parameter :: pi = 3.1415926535897932384_dp
   !pour un real (kind=dp) selected_real_kind(P=15)

   !integer, public, parameter             :: PRECISION_REAL     = 8
!!****d* m_qtlmap_constant/DP
!!  NAME
!!   DP
!!  DESCRIPTION
!!   Precision (kind) of real defined : 15 digits and exponent range greater at least 307.
!!  NOTES
!!
!!***
   integer, public, parameter             :: DP                 = SELECTED_REAL_KIND( 15, 307 )
    !***************************************************************************
   ! CONSTANTE PART
   !***************************************************************************

   integer, public, parameter             :: MAX_ANIMAL                 = 15000

   integer, public, parameter             :: MIN_NDMIN_PERMUTATION      = 10

   !-----------------------------------------------------------------------------
   ! dim : nb marker,nmes,2 : description of each allele by animal and by marker
   ! Modif olivier, optimisation -128<pheno<127
   ! Optimisation : codage des allele sur 256 entier
   ! 1) ==> amelioration des tps de calcul pour la comparaison d 'allele
   ! 2) ==> gain memoire
   !


!!****d* m_qtlmap_constant/KIND_PHENO
!!  NAME
!!   KIND_PHENO
!!  DESCRIPTION
!!   Precision (kind) of integer id allele defined
!!  NOTES
!!
!!***
   integer                    ,parameter     ,public :: KIND_PHENO    =    2
!!****d* m_qtlmap_constant/VAL_MIN_INDEX_PHENO
!!  NAME
!!   VAL_MIN_INDEX_PHENO
!!  DESCRIPTION
!!   Minimum value of an index allele
!!  NOTES
!!
!!***
   integer(kind=KIND_PHENO)   ,parameter     ,public :: VAL_MIN_INDEX_PHENO=1-(2**bit_size(VAL_MIN_INDEX_PHENO))/2
!!****d* m_qtlmap_constant/VAL_MAX_INDEX_PHENO
!!  NAME
!!   VAL_MAX_INDEX_PHENO
!!  DESCRIPTION
!!   Maximum value of an index allele
!!  NOTES
!!
!!***
   integer(kind=KIND_PHENO)   ,parameter     ,public :: VAL_MAX_INDEX_PHENO=(2**bit_size(VAL_MIN_INDEX_PHENO))/2 - 1
!!****d* m_qtlmap_constant/MAX_FILES_INPUT
!!  NAME
!!   MAX_FILES_INPUT
!!  DESCRIPTION
!!   Maximum files to read for input data
!!  NOTES
!!
!!***
   integer, public, parameter             :: MAX_FILES_INPUT    = 6
!!****d* m_qtlmap_constant/LENGTH_MAX_FILE
!!  NAME
!!   LENGTH_MAX_FILE
!!  DESCRIPTION
!!   Maximum length of a name file
!!  NOTES
!!
!!***
   integer, public, parameter             :: LENGTH_MAX_FILE    = 1024
!!****d* m_qtlmap_constant/LEN_S
!!  NAME
!!   LEN_S
!!  DESCRIPTION
!!   Length size of a short word
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_S              = 20
!!****d* m_qtlmap_constant/LEN_W
!!  NAME
!!   LEN_W
!!  DESCRIPTION
!!   Length size of a word
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_W              = 40
!!****d* m_qtlmap_constant/LEN_L
!!  NAME
!!   LEN_L
!!  DESCRIPTION
!!   Length size of a long word
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_L              = 80
!!****d* m_qtlmap_constant/LEN_LINE
!!  NAME
!!   LEN_LINE
!!  DESCRIPTION
!!   Maximum size length of a line readed
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_LINE           = 100000
!!****d* m_qtlmap_constant/LEN_DEF
!!  NAME
!!   LEN_DEF
!!  DESCRIPTION
!!   Default size of a string
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_DEF            = LEN_W
!!****d* m_qtlmap_constant/LEN_BUFFER_WORD
!!  NAME
!!   LEN_BUFFER_WORD
!!  DESCRIPTION
!!   Default size of a buffer string
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_BUFFER_WORD    = LEN_L
!!****d* m_qtlmap_constant/FMT_MAX_BUFFER
!!  NAME
!!   FMT_MAX_BUFFER
!!  DESCRIPTION
!!   Default maximum size of a fmt buffer string
!!  NOTES
!!
!!***
   integer, public, parameter             :: FMT_MAX_BUFFER     = 4048
!!****d* m_qtlmap_constant/LEN_BUFFER_LINE
!!  NAME
!!   LEN_BUFFER_LINE
!!  DESCRIPTION
!!   Default maximum size of a buffer string line
!!  NOTES
!!
!!***
   integer, public, parameter             :: LEN_BUFFER_LINE    = 200

!!****d* m_qtlmap_constant/LABEL_NAME_TRAIT
!!  NAME
!!   LABEL_NAME_TRAIT
!!  DESCRIPTION
!!   Name used for each trait when a model is defined by the keyword 'all'
!!  NOTES
!!
!!***
   character (len=10) , public, parameter :: LABEL_NAME_TRAIT   = 'TRAIT_ALL_'

!!****d* m_qtlmap_constant/FMT_REAL
!!  NAME
!!   FMT_REAL
!!  DESCRIPTION
!!   Default Real FMT to print in a console
!!  NOTES
!!
!!***
   character (len=10) , public, parameter :: FMT_REAL           = '(F8.5)'
!!****d* m_qtlmap_constant/FMT_INT
!!  NAME
!!   FMT_INT
!!  DESCRIPTION
!!   Default Integer FMT to print in a console
!!  NOTES
!!
!!***
   character (len=10) , public, parameter :: FMT_INT            = '(I8)'
!!****d* m_qtlmap_constant/FMT_INT_LONG
!!  NAME
!!   FMT_INT_LONG
!!  DESCRIPTION
!!   Default Long Integer FMT to print in a console
!!  NOTES
!!
!!***
   character (len=10) , public, parameter :: FMT_INT_LONG       = '(I9)'
!!****d* m_qtlmap_constant/REAL_NOT_DEFINED
!!  NAME
!!   REAL_NOT_DEFINED
!!  DESCRIPTION
!!   Default initialization of real without correct value
!!  NOTES
!!
!!***
   real (kind=dp),     public,parameter   :: REAL_NOT_DEFINED   = 99999.d0
!!****d* m_qtlmap_constant/INIFINY_REAL_VALUE
!!  NAME
!!   INIFINY_REAL_VALUE
!!  DESCRIPTION
!!   The biggest value using by the optimization package
!!  NOTES
!!
!!***
   real (kind=dp),     public,parameter   :: INIFINY_REAL_VALUE = REAL_NOT_DEFINED
!!****d* m_qtlmap_constant/INT_NOT_DEFINED
!!  NAME
!!   INT_NOT_DEFINED
!!  DESCRIPTION
!!    Default initialization of integer without correct value
!!  NOTES
!!
!!***
   integer, public, parameter             :: INT_NOT_DEFINED    = 9999
!!****d* m_qtlmap_constant/STRING_NOT_DEFINED
!!  NAME
!!   STRING_NOT_DEFINED
!!  DESCRIPTION
!!    Default initialization of string without correct value
!!  NOTES
!!
!!***
   character (len=1) , public, parameter  :: STRING_NOT_DEFINED = '*'
!!****d* m_qtlmap_constant/OPT_SIB_HS
!!  NAME
!!   OPT_SIB_HS
!!  DESCRIPTION
!!    Definition of a constante for half sib analysis
!!  NOTES
!!
!!***
   integer     ,parameter , public        :: OPT_SIB_HS         = 1 ! half sib
!!****d* m_qtlmap_constant/OPT_SIB_FS
!!  NAME
!!   OPT_SIB_FS
!!  DESCRIPTION
!!    Definition of a constante for full sib analysis
!!  NOTES
!!
!!***
   integer     ,parameter , public        :: OPT_SIB_FS         = 2 ! full sib
!!****d* m_qtlmap_constant/DEFAULT_PARAM_MIN
!!  NAME
!!   DEFAULT_PARAM_MIN
!!  DESCRIPTION
!!    Minimum value for a generic parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: DEFAULT_PARAM_MIN  = -1.d4
!!****d* m_qtlmap_constant/DEFAULT_PARAM_MAX
!!  NAME
!!   DEFAULT_PARAM_MAX
!!  DESCRIPTION
!!    Maximum value for a generic parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: DEFAULT_PARAM_MAX  = 1.d4
!!****d* m_qtlmap_constant/SIG_MIN
!!  NAME
!!   SIG_MIN
!!  DESCRIPTION
!!    Minimum value for a standart deviation parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: SIG_MIN            = 1.d-4
!!****d* m_qtlmap_constant/SIG_MAX
!!  NAME
!!   SIG_MAX
!!  DESCRIPTION
!!    Maximum value for a standart deviation parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: SIG_MAX            = DEFAULT_PARAM_MAX
!!****d* m_qtlmap_constant/XMU_MIN
!!  NAME
!!   XMU_MIN
!!  DESCRIPTION
!!    Minimum value for a mean parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: XMU_MIN            = DEFAULT_PARAM_MIN
!!****d* m_qtlmap_constant/XMU_MAX
!!  NAME
!!   XMU_MAX
!!  DESCRIPTION
!!    Maximum value for a mean parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: XMU_MAX            = DEFAULT_PARAM_MAX
!!****d* m_qtlmap_constant/AP_MIN
!!  NAME
!!   AP_MIN
!!  DESCRIPTION
!!    Minimum value for a sire qtl effect parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: AP_MIN             = DEFAULT_PARAM_MIN
!!****d* m_qtlmap_constant/AP_MAX
!!  NAME
!!   AP_MIN
!!  DESCRIPTION
!!    Maximum value for a sire qtl effect parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: AP_MAX             = DEFAULT_PARAM_MAX
!!****d* m_qtlmap_constant/AM_MIN
!!  NAME
!!   AM_MIN
!!  DESCRIPTION
!!    Minimum value for a dam qtl effect parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: AM_MIN             = DEFAULT_PARAM_MIN
!!****d* m_qtlmap_constant/AM_MAX
!!  NAME
!!   AM_MAX
!!  DESCRIPTION
!!    Maximum value for a dam qtl effect parameter for a optimization context
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: AM_MAX             = DEFAULT_PARAM_MAX

!!****d* m_qtlmap_constant/EPS_EM
!!  NAME
!!   EPS_EM
!!  DESCRIPTION
!!    maximum difference between two successive estimation of the haplotype proability
!!    in prob_haplo_complet
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: EPS_EM             = 1.d-4
!!****d* m_qtlmap_constant/ITER_EM_MAX
!!  NAME
!!   ITER_EM_MAX
!!  DESCRIPTION
!!    Maximum number of iteration in the EM in prob_haplo_complet
!!  NOTES
!!
!!***
   real (kind=dp)     ,parameter , public        :: ITER_EM_MAX             = 1000

!!!!*********************************************************** ID ANALYSIS

!!****d* m_qtlmap_constant/TYPE_DATA_CONTINUE
!!  NAME
!!   TYPE_DATA_CONTINUE
!!  DESCRIPTION
!!   Constant for real data
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: TYPE_DATA_CONTINUE                  = 1
!!****d* m_qtlmap_constant/TYPE_DATA_COX
!!  NAME
!!   TYPE_DATA_COX
!!  DESCRIPTION
!!   Constant for model cox data
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: TYPE_DATA_COX                       = 2

!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT
!!  NAME
!!   ANALYSE_UNITRAIT
!!  DESCRIPTION
!!   Constant id LA single trait with pre-corrected data
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT                    = 1
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_MODLIN
!!  NAME
!!   ANALYSE_UNITRAIT_MODLIN
!!  DESCRIPTION
!!   Constant id LA single trait with model description
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_MODLIN             = 2
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LA_HOMO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LA_HOMO
!!  DESCRIPTION
!!   Constant id LA for a single data with a model description (likelihood linearised ­ homoscedatic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LA_HOMO     = 3
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LA_HETERO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LA_HETERO
!!  DESCRIPTION
!!   Constant id LA for a single data with a model description (likelihood linearised ­ heteroscedastic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LA_HETERO   = 4
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LD_HOMO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LD_HOMO
!!  DESCRIPTION
!!   Constant id LD  for a single data with a model description   (likelihood linearised ­ homoscedatic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LD_HOMO     = 25
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LD_HETERO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LD_HETERO
!!  DESCRIPTION
!!   Constant id LD for a single data with a model description (likelihood linearised ­ heteroscedastic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LD_HETERO   = 26
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO
!!  DESCRIPTION
!!   Constant id LDLA  for a single data with a model description (likelihood linearised ­ homoscedatic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO   = 27
!!****d* m_qtlmap_constant/ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO
!!  NAME
!!   ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO
!!  DESCRIPTION
!!   Constant id LDLA  for a single data with a model description (likelihood linearised ­ homoscedatic)
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO = 28

!!   Constant id LD for a single data with a model description (likelihood linearised ­homoscedatic)
!!   A inbreed matrix is computed and used to estimate the model

   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO  = 29

!!   Constant id LA for a single data with a model description (likelihood linearised ­homoscedatic)
!!   A inbreed matrix is computed and used to estimate the model

   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO  = 30


   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LINEAR_RACE_HETERO  = 31

!!   Constant id LA for a set of traits with a multivariate analysis (based on a multi­normal penetrance function)
   integer        ,parameter,public              :: ANALYSE_MULTITRAIT                   = 5

!!    Constant id  LA for a set of traits (without missing data) with a discriminante analysis
   integer        ,parameter,public              :: ANALYSE_MULTITRAIT_DA               = 6

!!    Constant id LA for a single survey trait with the cox model
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_MODLIN_COX    = 7

!!    Constant id LD for a single data with a model description
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LD            = 8

!!    Constant id LDLA  for a single data with a model description
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LDLA          = 9

!!    Constant id LA for a single survey trait with the model description
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LA            = 10

!!    Constant id LA for a single survey trait with a model description ( DEV model LDJH )
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_LDJH          = 11

!!    Constant id 2 qtl interaction ( DEV )
   integer        ,parameter,public              :: ANALYSE_2QTL_INTERACTION       = 19

!!    Constant id LA for a single survey trait with a model description without interaction fixed effect - qtl (multiqtl version).
!!    this analysis used the structured data of m_qtlmap_incidence.F95
   integer        ,parameter,public              :: ANALYSE_UNITRAIT_CONTINGENCE   = 21

!!    Constant id LA for a single survey trait with a model description (other performance are treat as covariates (xiaoquiang.wang@rennes.inra.fr))
!!    this analysis used the structured data of m_qtlmap_incidence.F95
   integer        ,parameter,public              :: ANALYSE_TRAIT_COV_CONTINGENCE  = 22

!!    Constant id LA for a set of traits with a multivariate analysis with a model description (based on a multi­normal penetrance function) (multi-qtl version)
!!    this analysis used the structured data of m_qtlmap_incidence.F95
   integer        ,parameter,public              :: ANALYSE_MULTITRAIT_INCIDENCE   = 23

!!    Constant id LA for a set of traits with a multivariate analysis with a model description (based on a multi­normal penetrance function) (multi-qtl version)
!!    this analysis used the structured data of m_qtlmap_incidence.F95
   integer        ,parameter,public              :: ANALYSE_MULTITRAIT_INCIDENCE_LU   = 24
   

   integer        ,parameter,public              :: ANALYSE_TRAIT_BIALL_FARNIR        = 40

!!   Internal ID dev of method
   integer        ,parameter,public              :: ANALYSE_DEV_1   = 100

!!   Internal ID dev of method
   integer        ,parameter,public              :: ANALYSE_DEV_2   = 101

!!    Constant id for context analysis (analysis case)
   integer, public, parameter                    :: COMMON_ANALYSE         = 0
!!    Constant id for context analysis (simulation case)
   integer, public, parameter                    :: SIMULATION             = 1
!!    Constant id for context analysis (transcriptom analysis case)
   integer, public, parameter                    :: TRANSCIPTOME_ANALYSE   = 2
!!   "F2" label in the simulation parameter file.
   character(len=LEN_BUFFER_WORD), public, parameter :: F2_KEYWORD      = 'F2'
!!   Backcross label in the simulation parameter file.
   character(len=LEN_BUFFER_WORD), public, parameter :: BC_KEYWORD      = 'BC'
!!   Outbred label in the simulation parameter file.
   character(len=LEN_BUFFER_WORD), public, parameter :: OUTBRED_KEYWORD = 'OUTBRED'
!!   Label markers description
  character(len=20)                                     ,parameter    :: LABEL_MARKERS       = "MARKERS"
!!   Label genealogy description
  character(len=20)                                     ,parameter    :: LABEL_GENEALOGY     = "GENEALOGY"
!!   Label QTL description
  character(len=20)                                     ,parameter    :: LABEL_QTL           = "QTL"
!!   Label phenotypic data description
  character(len=20)                                     ,parameter    :: LABEL_TRAITS        = "TRAITS"
!!   Label simulation phenotypic data description
character(len=20)                                       ,parameter    :: LABEL_SIMULTRAITS   = "SIMULTRAITS"
!!   Keyword to apply a generic model for all files
character(len=10)                                       ,parameter    :: ALL_LABEL_MODEL     = 'ALL'
!!   Keyword to apply a generic model for all files
character(len=10)                                       ,parameter  :: FILTER_TRAITS_MODEL   = 'TRAITS'

!!   Maximum number of generic key (opt_key_dev1,opt_key_dev2,...)
   integer                    ,parameter    :: MAXNB_KEY_DEV_GEN               = 100

!!======================================================================================================
!! CONFIDENCE Intervals Methods

!!****d* m_qtlmap_constant/DROP_OFF_CI
!!  NAME
!!   DROP_OFF_CI
!!  DESCRIPTION
!!
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: DROP_OFF_CI      = 1

!!****d* m_qtlmap_constant/BOOTSTRAP_FULL_CI
!!  NAME
!!   BOOTSTRAP_FULL_CI
!!  DESCRIPTION
!!
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: BOOTSTRAP_FULL_CI = 2

!!****d* m_qtlmap_constant/BOOTSTRAP_FULL_CI
!!  NAME
!!   BOOTSTRAP_FULL_CI
!!  DESCRIPTION
!!
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: BOOTSTRAP_SIB_CI  = 3

!!****d* m_qtlmap_constant/HENGDE_LI_CI
!!  NAME
!!   HENGDE_LI_CI
!!  DESCRIPTION
!!   Ref article : a quick method to calculate QTL confidence interval
!!  NOTES
!!
!!***
   integer        ,parameter,public              :: HENGDE_LI_CI  = 4


!!****d* m_qtlmap_constant/OPTI_LBFGS
!! NAME
!!   OPTI_LBFGS
!! DESCRIPTION
!!   constant id L-BFGS (Nocedal) optimization
!! NOTES
!!   http://users.eecs.northwestern.edu/~nocedal/lbfgs.html
!! SOURCE
      integer , parameter , public            :: OPTI_LBFGS                = 2
!!***

!!****d* m_qtlmap_constant/OPTI_LAST
!! NAME
!!   OPTI_LAST
!! DESCRIPTION
!!   Last index of optimization listed
!! SOURCE
      integer , parameter , public            :: OPTI_LAST                 = OPTI_LBFGS



!!   Constant when phases are provided by a external software
   integer, public, parameter             :: VERSION_HAPLOTYPE_PARENTAL_EXTERNAL = 0
!!   Constant for the first version of the haplotypes routines
   integer, public, parameter             :: VERSION_HAPLOTYPE_V1                = 1  ! Original version

!!   Constant for the second version of the haplotypes routines
   integer, public, parameter             :: VERSION_HAPLOTYPE_V2                = 2  !

!!   Constant for the third version of the haplotypes routines
   integer, public, parameter             :: VERSION_HAPLOTYPE_V3                = 3  !
!
!!   Constant for the 4th version of the haplotypes routines
!!   Elsen JM. A fast algorithm for estimating transmission probabilities in QTL detection designs with dense maps.
   integer, public, parameter             :: VERSION_HAPLOTYPE_SNP               = 4  !

!!   Constant for the 5th version of the haplotypes routines
!!   Phase construction (Favier et al. 2010), Elsen JM. A fast algorithm for estimating transmission probabilities in QTL detection designs with dense maps.
   integer, public, parameter             :: VERSION_HAPLOTYPE_SYMMAX2SAT_SNP    = 5  !

!!***

end module m_qtlmap_constant
