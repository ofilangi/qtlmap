!!  NAME
!!    m_qtlmap_cli
!!  SYNOPSIS
!!   A command-line interface (CLI) is a mechanism for interacting with th QTLMap software by typing commands to perform specific tasks.
!!
!!  DESCRIPTION
!!   This module manage the options give as argument to qtlmap.
!!

module m_qtlmap_type_cli
    use m_qtlmap_constant
    use m_qtlmap_base

    implicit none

    !! Option maximum length
    integer                            ,parameter,private   :: LEN_OPT         = 44
    !! blank to print before the help output
    integer                            ,parameter,private   :: MARG_HELP       = 10
    character(len=2)                   ,parameter,private   :: CMARG_HELP      = "10"

    private :: key_exist

    type , public :: QTLMAP_CLI

         !!// default simulation number if is not specified by the user
        integer                                      ,private   :: NSIM_DEFAULT    = 1000
        !!// default hypothesis number (qtl) to test if is not specified by the user
        integer                                      ,private   :: NQTL_DEFAULT    = 1
        !!// default value of the family option (full sib family is the default).
        integer                                      ,private   :: SIB_DEFAULT     = 2
        !!// default optimization method using by qtlmap ( LA optimized analysis )
        integer                                      ,private   :: OPTIM_DEFAULT   = 1

        !!// /* OPTION NAME */
        !!// calcul option name, see m_qtlmap_analyse module to have the available analysis
        character(len=LEN_OPT)                       ,private :: OPT_CALCUL              ='--calcul'
        !!// hypothesis option name, some analysis have limitation with this value.
        character(len=LEN_OPT)                       ,private :: OPT_QTL                 ='--qtl'
        !!// haplotype option name, see m_qtlmap_haplotype/haplotype to see the availables haplotypes versions
        character(len=LEN_OPT)                       ,private :: OPT_HAPLOTYPE           ='--haplotype'
        !!// family option name
        character(len=LEN_OPT)                       ,private :: OPT_SIBF                ='--family'
        !!// optimization option name
        character(len=LEN_OPT)                       ,private :: OPT_OPTIM               ='--optim'
        !!// random simulation option name
        character(len=LEN_OPT)                       ,private :: OPT_SIM                 ='--simulation'
        !!// number of simulation option name
        character(len=LEN_OPT)                       ,private :: OPT_NSIM                ='--nsim'
        !!// expression phenotypics values option name
        character(len=LEN_OPT)                       ,private :: OPT_EQTL                ='--data-transcriptomic'
        !!// permutation simulation option name
        character(len=LEN_OPT)                       ,private :: OPT_PERMUT              ='--permute'
        !!// snp option name
        character(len=LEN_OPT)                       ,private :: OPT_SNP                 ='--snp'
        !!// haplotype dam option name
        character(len=LEN_OPT)                       ,private :: OPT_NO_HDAM             ='--no-hdam'
        !!// biallelic option
        character(len=LEN_OPT)                       ,private :: OPT_BIQ                 ='--biq'
        !!// interaction between 2 qtls option
        character(len=LEN_OPT)                       ,private :: OPT_INTERACTION         ='--interaction'
        !!// censored data option name, Used to take care about a third generation with trait
        character(len=LEN_OPT)                       ,private :: OPT_CENSORED_CALCUL     ='--calcul-cd'
        !!// threshold estimation option name
        character(len=LEN_OPT)                       ,private :: OPT_ESTIME_THRESHOLDS   = '--estimate-thresholds'
        !!// threshold estimation short option name
        character(len=LEN_OPT)                       ,private :: OPT_ESTIME_THRESHOLDS_SHORT = '-s'
        !!// quiet output console option name
        character(len=LEN_OPT)                       ,private :: OPT_QUIET              ='--quiet'
        !!// verbose output console option name
        character(len=LEN_OPT)                       ,private :: OPT_VERBOSE            ='--verbose'
        !!// debug output console option name
        character(len=LEN_OPT)                       ,private :: OPT_DEBUG              ='--debug'
        !!// quiet output console short option name
        character(len=LEN_OPT)                       ,private :: OPT_QUIET_SHORT        ='-q'
        !!// verbose output console short option name
        character(len=LEN_OPT)                       ,private :: OPT_VERBOSE_SHORT      ='-v'
        !!// debug output console short option name
        character(len=LEN_OPT)                       ,private :: OPT_DEBUG_SHORT        ='-d'
        !!// used to print a classical output analysis in a simulation case or expressions phenotypics values analysis.
        character(len=LEN_OPT)                       ,private :: OPT_PRINT_CLASSIC      ='--print-all'
        !!// help parameter file option name
        character(len=LEN_OPT)                       ,private :: OPT_HELP_PANALYSE      = "--help-panalyse"
        !!// help simulation parameter file option name
        character(len=LEN_OPT)                       ,private :: OPT_HELP_PARAMSIM      = "--help-paramsim"
        !!// help option name
        character(len=LEN_OPT)                       ,private :: OPT_HELP               ='--help'
        !!// help option name
        character(len=LEN_OPT)                       ,private :: OPT_HELP_ALL           ='--help-dev'
        !!// suppprime l'estimation des effets polygeniques
        character(len=LEN_OPT)                       ,private :: OPT_NOPOLY             ='--nopoly'
        !!// overload keys defined in the parameter file
        character(len=LEN_OPT)                       ,private :: OPT_OVERLOADED_KEYS    ='--overloaded-keys'
        !!// do not estime Qtl sire effect (useful for analyse of female hetero-chromosme)
        character(len=LEN_OPT)                       ,public  :: DISABLE_SIRE_QTL       ='--disable-sire-qtl'
        !!// confidence intervalle option name option name
        character(len=LEN_OPT)                       ,public  :: OPT_CI                 ='--ci'
        !!// number of simulation for the bootstrap option to calculate Confidence Intervalle
        character(len=LEN_OPT)                       ,public  :: OPT_CI_NSIM            ='--ci-nsim'

    contains

        !// Print QtlMAP Version
        procedure,public :: cli_print_qtlmap_version

        !// Print help
        procedure,public :: cli_print_help

        !// True if the user ask help message
        procedure,public :: cli_is_help

        !// True if the user ask help-dev message
        procedure,public :: cli_is_help_all

        !// Get the analyse parameters file
        procedure,public :: cli_get_parameters_file

        !// Get the analyse type
        procedure,public :: cli_get_analyse

        !// Get the Haplotype version
        procedure,public :: cli_get_haplotype

        !// Get the Family version
        procedure,public :: cli_get_family

        !// Get the optimization
        procedure,public :: cli_get_optimization

        !// True if the user want estimate rejection threholds of the current analysis
        procedure,public :: cli_is_estimate_thresholds

        !// True if the user ask a simulation
        procedure,public :: cli_is_simulation

        !// True if the user ask a transcriptomic analysis/simulation
        procedure,public :: cli_is_transcriptomic_data

        !// True if the user ask a permutation method for the simulation
        procedure,public :: cli_is_permute

        !// True if the user ask won't estimate haplotype dam in the model LD and derivated LD analysis
        procedure,public :: cli_is_no_hdam

        !// True if the user ask to estmate a biallelic qtl
        procedure,public :: cli_is_biq

        !// True if the user ask to estmate an qtl in interaction
        procedure,public :: cli_is_interaction

        !// True if the user ask an iteration and so a simulation
        procedure,public :: cli_is_nsim

        !// Get the number of the simulation
        procedure,public :: cli_get_nsim

        !// Get the number of the qtl to detect
        procedure,public :: cli_get_nqtl

        !// True if the user ask none information
        procedure,public :: cli_is_quiet

        !// True if the user ask verbose information
        procedure,public :: cli_is_verbose

        !// True if the user ask to print all information
        procedure,public :: cli_is_print_all

        !// True if the user ask debug information
        procedure,public :: cli_is_debug

        !// True if the user ask help p_analyse information
        procedure,public :: cli_is_help_panalyse

        !// True if the user ask help paramsim information
        procedure,public :: cli_is_help_paramsim

        !// True if the user want used the F3 phenotypes
        procedure,public  :: cli_is_calcul_cd

        procedure,public  :: cli_get_overloaded_keys

        procedure,public  :: release => release_qtlmap_cli
        procedure,public  :: key_exist
        procedure,public  :: get_key_list_values
        procedure,public  :: get_key_value
        procedure ,public :: copy    => copy_qtlmap_cli

    end type QTLMAP_CLI

contains

    !!****f* m_qtlmap_type_cli/copy_qtlmap_cli
    !!  NAME
    !!    copy_qtlmap_cli
    !!  DESCRIPTION
    !!    copy
    !!  NOTES
    !!  SOURCE
    subroutine copy_qtlmap_cli(this,copy)
        class(QTLMAP_CLI)  , intent(in)     :: this
        type(QTLMAP_CLI)   , intent(inout)  :: copy

    end subroutine copy_qtlmap_cli
    !!***


    !!****f* m_qtlmap_type_cli/release_qtlmap_cli
    !!  NAME
    !!    release_qtlmap_cli
    !!  DESCRIPTION
    !!    release memory
    !!  NOTES
    !!  SOURCE
    subroutine release_qtlmap_cli(this)
        class(QTLMAP_CLI)  , intent(inout)  :: this

    end subroutine release_qtlmap_cli
    !!***

    !!****f* m_qtlmap_type_cli/cli_print_qtlmap_version
    !!  NAME
    !!    cli_print_qtlmap_version
    !!  DESCRIPTION
    !!    Print the version of qtlmap and the date of the binary build.
    !!  NOTES
    !!  SOURCE
    subroutine cli_print_qtlmap_version(this)
        class(QTLMAP_CLI)  , intent(in)  :: this
#ifdef QTLMAP_VERSION
#ifdef DATE_BUILD
            write (*,fmt="(a)")'Release-Build:',QTLMAP_VERSION,'-',DATE_BUILD
#endif
#endif
    end subroutine cli_print_qtlmap_version
    !!***

    !!****f* m_qtlmap_type_cli/cli_print_help
    !!  NAME
    !!    cli_print_help
    !!  DESCRIPTION
    !!    Print command line interface help
    !!  NOTES
    !!  SOURCE
    subroutine cli_print_help(this,complete)
        class(QTLMAP_CLI)  , intent(in)  :: this

        logical, optional :: complete
        integer        :: i,nlopt,iopt
        character(len=200) , dimension(40)  :: info,test
        external :: get_nb_opti_nlopt_qtlmap,info_nlopt_qtlmap
        logical :: compl = .false.

        if ( present(complete) ) compl = complete

        write (*,fmt="(a,a,a)") 'Usage: ',trim(get_binary_name()),' [FILE_ANALYSE] [OPTIONS]...'
        write (*,fmt="(a)")&
            'QTLmap is a linkage analysis software dedicated to the analysis of outbred population in an experimental design.'
        write (*,fmt="(a)")''
        write (*,fmt="(a)") &
            'The QTLMAP software is developed at INRA (french National Institute for Agronomical Research) to map QTL'
        write (*,fmt="(a)") &
            'following linkage analyses approaches in various experimental populations. Specific models were implemented to'
        write (*,fmt="(a)")&
            'take maximum advantage of large data sets now available for QTL detection and get fine locus descriptions'
        write (*,fmt="(a)")&
            '(indications of pleiotropy, of linked loci…) to facilitate further explorations. In QTLMAP, '
        write (*,fmt="(a)")&
            'the basis of the QTL detection technique is an interval mapping applied to half-sib families or a mixture '
        write (*,fmt="(a)")&
            'of full- and half- sib families where the exact computation of phase probability is undertaken from a 2 or 3'
        write (*,fmt="(a)")&
            'generation pedigree. The test statistics are computed either using exact calculation of likelihood ratios '
        write (*,fmt="(a)")&
            'or using partial linear approximations of the likelihood or complete regression on marker genotypes. '
        write (*,fmt="(a)")''
        write (*,fmt="(a)")'   ***  support : send email to qtlmap-users@listes.inra.fr. *** '
        write (*,fmt="(a)")'                              *                                  '

        !----------- ANALYSIS ----------
        write (*,fmt="(a)")'Analysis options'
        write (*,fmt="(a)")'----------------'
        write (*,fmt="(a,a)") trim(this%OPT_CALCUL)   ,'='
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT)),&
            ' : Analyse Trait by trait (precorrected data) - qtl={1,2} '
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_MODLIN)) &
            ,' : Analyse Trait by trait with fixed effects, covariates and qtl interaction with fixed effect',&
            ' - qtl={1} '
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LA_HOMO)),&
            ' : Linear homoscedastic LA with fixed effects and covariates (using the most probable genotype dam)',&
            '- qtl={1,..,n}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LA_HETERO)),&
            ' : Linear heteroscedastic LA with fixed effects and covariates (using the most probable genotype dam)',&
            ' - qtl={1,..,n}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_MULTITRAIT)), &
            ' : Join traits analysis (precorrected data) - qtl={1} '
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_MULTITRAIT_DA)),' : Joint analysis of the traits',&
            ' using a discriminant function (precorrected data) - qtl={1}.'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a,a)")trim(str(ANALYSE_UNITRAIT_MODLIN_COX)),' :',&
            ' survival analysis with fixed ','effects and covariates (Cox model) - qtl={1}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)")trim(str(ANALYSE_UNITRAIT_CONTINGENCE)) &
            ,' : Analyse Trait by trait with fixed effects and covariates  - qtl={1,..,n} '
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LD_HOMO)),&
            ' : Linear homoscedastic LD with fixed effects and covariates (using the most probable genotype dam)',&
            '  - qtl={1,..,n}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LD_HETERO)),&
            ' : Linear heteroscedastic LD with fixed effects and covariates (using the most probable genotype dam)',&
            ' - qtl={1,..,n}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO)),&
            ' : Linear homoscedastic LDLA with fixed effects and covariates (using the most probable genotype dam)',&
            '  - qtl={1,..,n}'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)")trim(str(ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO)),&
            ' : Linear heteroscedastic LDLA with fixed effects and covariates (using the most probable genotype dam)',&
            ' - qtl={1,..,n}'

        ! Mettre ici Toutes les analyses en developpements ou à ne pas montrer pour une utyilisation courante de QTLMAP

        if (compl) then
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LD)),&
                ' :(STABLE/REF) LD analysis with fixed effects covariates and qtl interaction with fixed effect - qtl={1} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LDLA)),&
                ' : (STABLE/REF)  LDLA analysis with fixed effects covariates and qtl interaction with fixed effect - qtl={1} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LA)),&
                ' : (STABLE/REF)  LA analysis with fixed effects covariates and qtl interaction with fixed effect - qtl={1} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LDJH)),&
                ' : (DEV) LDJH analysis(dev jm) with fixed effects , covariates and qtl interaction with fixed effect - qtl={1} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_TRAIT_BIALL_FARNIR)) &
                ,' : (DEV) Analyse BiAllelic Farnir              - with fixed effects and covariates  - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_TRAIT_COV_CONTINGENCE)) &
                ,' : (DEV) Analyse with other trait as covariate - with fixed effects and covariates  - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_2QTL_INTERACTION)) &
                ,' : (DEV) Analyse 2 QTL with an interaction  - with fixed effects and covariates  - qtl=2 '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_MULTITRAIT_INCIDENCE)), &
                ' : (DEV) Join traits analysis with fixed effects, covariates - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_MULTITRAIT_INCIDENCE_LU)), &
                ' : (DEV) Join traits analysis (LU decomposition) with fixed effects, covariates - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO)) &
                ,' : (DEV) LD - A inbreed matrix is computed and used to estimate the model - ',&
                'with fixed effects and covariates  - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO)) &
                ,' : (DEV) LA - A inbreed matrix is computed and used to estimate the model - ',&
                'with fixed effects and covariates  - qtl={1,..,n} '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)) &
                ,' : '

        end if

        write (*,fmt="(a,a)")trim(this%OPT_QTL)      ,'='
        write (*,fmt="("//CMARG_HELP//"(' '),a)")'<integer> : Number of QTL to include in the model (default=1)'
        write (*,*)
        write (*,fmt="(a)")trim(this%DISABLE_SIRE_QTL)
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") 'Disable the Estimation of QTL sire family '//&
            '(analyse of female hetero chromosome).'
        write (*,*)
        write (*,fmt="(a)")"Confidence Intervals options"
        write (*,fmt="(a)")"----------------------------"
        write (*,fmt="(a,a)")trim(this%OPT_CI)       ,'='
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(DROP_OFF_CI)),' : Drop-off method'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(BOOTSTRAP_FULL_CI)),' : Boostrap resampling method'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a)") trim(str(BOOTSTRAP_SIB_CI)),' : Boostrap resampling method '//&
            '(keep exactly the number of progeny within a family)'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(HENGDE_LI_CI)),&
            ' : Hengde Li method using the Relative Frequency Ratio'
        write (*,*)
        write (*,fmt="(a,a)")trim(this%OPT_CI_NSIM)     ,'='
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a,a)") 'Number of simulation for Confidence Intervals/boostrap methods',&
            '[default=',trim(str(this%NSIM_DEFAULT)),'].'

        !        write (*,fmt="(a)")trim(OPT_BIQ)
        !        print *                    ,(' ',i=1,MARG_HELP),'** Devel ** To evaluate a biallelic qtl effect.'
        write (*,fmt="(a)")trim(this%OPT_INTERACTION)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'To evaluate the effect of an interaction between 2 qtls.'
        write (*,fmt="(a)")trim(this%OPT_NO_HDAM)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'remove haplotype dam in the LD/LDLA model analysis.'
        write (*,fmt="(a)")trim(this%OPT_CENSORED_CALCUL)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Calcul the cd of each kd from the 3th generation.'
        write (*,fmt="(a)")trim(this%OPT_EQTL)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Perform an analysis with transcriptomic datas.'
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") 'The Traits file is reading by line and the result',&
            ' log is significantly different than a classic analysis'
        write (*,fmt="(a)")trim(this%OPT_SNP)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Perform an analysis with SNP markers.'
        write (*,fmt="(a,a,a)")trim(this%OPT_ESTIME_THRESHOLDS),',',trim(this%OPT_ESTIME_THRESHOLDS_SHORT)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Estimate the rejection thresholds.'

        !----------- SIMUL ----------
        write (*,fmt="(a)")''
        write (*,fmt="(a)")'Simulation options'
        write (*,fmt="(a)")'------------------'
        write (*,fmt="(a)")trim(this%OPT_SIM)
        !        print *                    ,(' ',i=1,MARG_HELP),'To perform a simulation for a threshold estimation or',&
        !        ' looking for a design optimization.'
        write (*,fmt="(a,a)")trim(this%OPT_NSIM)     ,'='
        write (*,fmt="("//CMARG_HELP//"(' '),a,a,a,a)") 'Number of simulation ','[ default=',trim(str(this%NSIM_DEFAULT)),'].'

        write (*,fmt="(a)")trim(this%OPT_PERMUT)
        write (*,fmt="("//CMARG_HELP//"(' '),a,a)") 'To permut performances and associated effects between genotyped',&
            ' animals within half or full-sib families..'

        !-----------OUTPUT ----------
        write (*,fmt="(a)")''
        write (*,fmt="(a)")'Output options'
        write (*,fmt="(a)")'------------------'
        write (*,fmt="(a,a,a)")trim(this%OPT_QUIET),',',trim(this%OPT_QUIET_SHORT)

        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Display only warning and error messages.'
        write (*,fmt="(a,a,a)")trim(this%OPT_VERBOSE),',',trim(this%OPT_VERBOSE_SHORT)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Display verbose message.'

        write (*,fmt="(a,a,a)")trim(this%OPT_DEBUG) ,',',trim(this%OPT_DEBUG_SHORT)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Display debug and warning message.'
        write (*,fmt="(a)")trim(this%OPT_PRINT_CLASSIC)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Print all information in the result file (useful for eqtl analysis).'
        if ( compl ) then
            !----------- ADNVANCED OPTION ----------
            write (*,fmt="(a)")''
            write (*,fmt="(a)")'Advanced options'
            write (*,fmt="(a)")'------------------'
            write (*,fmt="(a,a)")trim(this%OPT_SIBF)      ,'='
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '1 : Half Sib family'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '2 : Mixture of full and half sib family (default)'
            write (*,fmt="(a,a)")trim(this%OPT_OPTIM)      ,'='
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '1 : * NAG * E04JYF '
            write (*,fmt="("//CMARG_HELP//"(' '),a,a)") '2 : * L-BFGS-B * Broyden–Fletcher–Goldfarb–Shanno update',&
                ' to approximate the Hessian matrix (L-BFGS stands for (limited memory BFGS))'
#ifdef HAVE_NLOPT
        call get_nb_opti_nlopt_qtlmap(nlopt)
        info=''
        call info_nlopt_qtlmap(info)

        do iopt=1,nlopt
           write (*,fmt="("//CMARG_HELP//"(' '),a,a)") trim(str(OPTI_LAST+iopt))//' : * NLopt * ',trim(info(iopt))
        end do
#endif
            write (*,fmt="(a)")''
            write (*,fmt="(a)")trim(this%OPT_OVERLOADED_KEYS)
            write (*,fmt="("//CMARG_HELP//"(' '),a)") "='<keys1>=<value1>,<keys2>=<value2>...'"//&
                " overload keys defined in the parameter file."
            write (*,fmt="(a)")''


            write (*,fmt="(a,a)")trim(this%OPT_HAPLOTYPE),'='
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '0 : parents phases are given by the user',&
                ' [reading order in the genotype file]'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '1 : old implementation'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '2 : PDEGM_V4,PDEGP_V4,PDED_V4'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '3 : PDEGM_V4,PDEGP_V4,PDED_V5'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '4 : PHASE FAST, PDED_V5'
            write (*,fmt="("//CMARG_HELP//"(' '),a)") '5 : PHASE SYM2SAT,PDED_V5 (--snp)'
        end if

        !-----------  HELP ---------
        write (*,fmt="(a)")'Help'
        write (*,fmt="(a)")'----'
        write (*,fmt="(a)")trim(this%OPT_HELP_PANALYSE)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Print help to build p_analyse file.'
        write (*,fmt="(a)")trim(this%OPT_HELP_PARAMSIM)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Print help to build simulation parameter file.'
        write (*,fmt="(a)")trim(this%OPT_HELP)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Display the help.'
        write (*,fmt="(a)")trim(this%OPT_HELP_ALL)
        write (*,fmt="("//CMARG_HELP//"(' '),a)") 'Display the help for advanced users.'

    end subroutine cli_print_help
    !!***

    !!    True if the user ask help message
    function cli_is_help(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_HELP)
    end function cli_is_help

    !!    True if the user ask help message
    function cli_is_help_all(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_HELP_ALL)
    end function cli_is_help_all

    !!    Get the analyse parameters file
    function cli_get_parameters_file(this) result(f)
        class(QTLMAP_CLI)  , intent(in)  :: this
        character(len=LENGTH_MAX_FILE) :: f
        if (get_nb_arg() <=0) then
            write (*,fmt="(a)")trim(get_binary_name()),' can not find parameters file analyse (first argument).'
            call this%cli_print_help
            stop
        end if

        call get_arg(1,f)
    end function cli_get_parameters_file

    !!    Get the analyse parameters file
    subroutine cli_get_overloaded_keys(this,listKeys,listValues,sizeList,nkeys)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer , intent(in)  :: sizeList
        character(len=LEN_DEF) ,dimension(sizeList),intent(out):: listKeys,listValues
        integer , intent(out) :: nkeys
        logical :: val
        character(len=LEN_DEF)   :: buf,subbuf
        integer :: i, ikey

        nkeys=0
        val = this%key_exist(this%OPT_OVERLOADED_KEYS)
        if (val) then
            call this%get_key_value(this%OPT_OVERLOADED_KEYS,buf)
            nkeys=1
            do while (trim(buf) /= '')
                i = index (buf,',')
                if ( i< len(buf)) then
                    subbuf = buf(1:i)
                else
                    subbuf = buf
                end if
                buf = buf(i+1:)
                i = index (subbuf,'=')
                listKeys(nkeys) = subbuf(1:i-1)
                listValues(nkeys) = subbuf(i+1:)
                nkeys = nkeys + 1
            end do

        end if

    end subroutine cli_get_overloaded_keys
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_analyse
    !!  NAME
    !!    cli_get_analyse
    !!  DESCRIPTION
    !!   Get the analyse type
    !!  NOTES
    !!  SOURCE
    function cli_get_analyse(this) result(type)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer                :: type
        character(len=LEN_DEF)   :: value

        type = 0
        call this%get_key_value(this%OPT_CALCUL,value)
        type = get_int(value)
        if (type == 0) type = 1  ! Default first option analyse

    end function cli_get_analyse
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_haplotype
    !!  NAME
    !!    cli_get_haplotype
    !!  DESCRIPTION
    !!   Get the Haplotype version
    !!  NOTES
    !!  SOURCE
    function cli_get_haplotype(this) result(type)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer                :: type
        character(len=LEN_DEF)   :: value

        type = 0
        if (.not. this%key_exist(this%OPT_HAPLOTYPE)) then
            if (this%key_exist(this%OPT_SNP)) then
                type = VERSION_HAPLOTYPE_SYMMAX2SAT_SNP
            else
                type = VERSION_HAPLOTYPE_V3  ! Default first option haplotype
            end if
        else
            call this%get_key_value(this%OPT_HAPLOTYPE,value)
            type = get_int(value)
        end if

    end function cli_get_haplotype
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_simulation
    !!  NAME
    !!    cli_is_simulation
    !!  DESCRIPTION
    !!   True if the user ask a simulation
    !!  NOTES
    !!  SOURCE
    function cli_is_simulation(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_SIM)
    end function cli_is_simulation
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_estimate_thresholds
    !!  NAME
    !!    cli_is_estimate_thresholds
    !!  DESCRIPTION
    !!   True if the user want estimate rejection threholds of the current analysis
    !!
    !!  OUTPUT
    !!   nsim : the number of simulation
    !!  NOTES
    !!  SOURCE
    function cli_is_estimate_thresholds(this,nsim) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        integer, intent(out)   :: nsim
        character(len=LEN_DEF)   :: value

        val = this%key_exist(this%OPT_ESTIME_THRESHOLDS)
        if ( val ) then
            call this%get_key_value(this%OPT_ESTIME_THRESHOLDS,value)
        else
            val = this%key_exist(this%OPT_ESTIME_THRESHOLDS_SHORT)
            if ( val ) then
                call this%get_key_value(this%OPT_ESTIME_THRESHOLDS_SHORT,value)
            end if
        end if

        nsim = 0
        if ( .not. val ) return

        nsim = this%NSIM_DEFAULT ! default value
        if (trim(value)/='') then
            nsim = get_int(value)
            if (nsim == huge(nsim)) then
                write (*,fmt="(a)")'bad value defined :'// trim(value) //'['//trim(this%OPT_NSIM)//']'
                call this%cli_print_help
                stop
            end if
        end if

    end function cli_is_estimate_thresholds
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_transcriptomic_data
    !!  NAME
    !!    cli_is_transcriptomic_data
    !!  DESCRIPTION
    !!   True if the user ask a transcriptomic analysis/simulation
    !!  NOTES
    !!  SOURCE
    function cli_is_transcriptomic_data(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_EQTL)
    end function cli_is_transcriptomic_data
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_permute
    !!  NAME
    !!    cli_is_permute
    !!  DESCRIPTION
    !!   True if the user ask permutation for the current simulation
    !!  NOTES
    !!  SOURCE
    function cli_is_permute(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_PERMUT)
    end function cli_is_permute
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_calcul_cd
    !!  NAME
    !!    cli_is_calcul_cd
    !!  DESCRIPTION
    !!    True if the user ask to compute a CD from a third generation with trait value
    !!  NOTES
    !!  SOURCE
    function cli_is_calcul_cd(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_CENSORED_CALCUL)
    end function cli_is_calcul_cd
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_hdam
    !!  NAME
    !!    cli_is_hdam
    !!  DESCRIPTION
    !!    True if the user ask to take care about haplotype dam in the model LD and derivated LD analysis
    !!  NOTES
    !!  SOURCE
    function cli_is_no_hdam(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_NO_HDAM)
    end function cli_is_no_hdam
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_biq
    !!  NAME
    !!    cli_is_biq
    !!  DESCRIPTION
    !!    True if the user ask to take care about a biallelic qtl effect
    !!  NOTES
    !!  SOURCE
    function cli_is_biq(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_BIQ)
    end function cli_is_biq
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_interaction
    !!  NAME
    !!    cli_is_interaction
    !!  DESCRIPTION
    !!    True if the user ask to take care about an interaction betwwen 2 qtl
    !!  NOTES
    !!  SOURCE
    function cli_is_interaction(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_INTERACTION)
    end function cli_is_interaction
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_nopoly
    !!  NAME
    !!    cli_is_nopoly
    !!  DESCRIPTION
    !!    True if the user ask to avoid the estimation of polygenic effect
    !!  NOTES
    !!  SOURCE
    function cli_is_nopoly(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_NOPOLY)
    end function cli_is_nopoly
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_nsim
    !!  NAME
    !!    cli_is_nsim
    !!  DESCRIPTION
    !!    True if the user ask an iteration and so a simulation
    !!  NOTES
    !!  SOURCE
    function cli_is_nsim(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_NSIM)
    end function cli_is_nsim
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_nsim
    !!  NAME
    !!    cli_get_nsim
    !!  DESCRIPTION
    !!    Get the number of the simulation
    !!  NOTES
    !!  SOURCE
    function cli_get_nsim(this) result(nsim)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer     :: nsim
        character(len=LEN_DEF)   :: value
        nsim = this%NSIM_DEFAULT ! default value
        call this%get_key_value(this%OPT_NSIM,value)
        if (trim(value)/='') then
            nsim = get_int(value)
            if (nsim == huge(nsim)) then
                write (*,fmt="(a)")'bad value defined :'// trim(value) //'['//trim(this%OPT_NSIM)//']'
                call this%cli_print_help()
                stop
            end if
        end if
    end function cli_get_nsim
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_nqtl
    !!  NAME
    !!    cli_get_nqtl
    !!  DESCRIPTION
    !!    Get the number of the qtl hypothesis to test
    !!  NOTES
    !!  SOURCE
    function cli_get_nqtl(this) result(nqtl)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer     :: nqtl
        character(len=LEN_DEF)   :: value

        nqtl = this%NQTL_DEFAULT ! default value
        call this%get_key_value(this%OPT_QTL,value)
        if (trim(value)/='') then
            nqtl = get_int(value)
            if (nqtl == huge(nqtl)) then
                write (*,fmt="(a)")'bad value defined :'// trim(value) //'['//trim(this%OPT_QTL)//']'
                call this%cli_print_help
                stop
            end if
        end if
    end function cli_get_nqtl
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_family
    !!  NAME
    !!    cli_get_family
    !!  DESCRIPTION
    !!    Get the type of family to analyse
    !!  NOTES
    !!  SOURCE
    function cli_get_family(this) result(o_sib)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer     :: o_sib
        character(len=LEN_DEF)   :: value

        o_sib = this%SIB_DEFAULT ! default value
        call this%get_key_value(this%OPT_SIBF,value)
        if (trim(value)/='') then
            o_sib = get_int(value)
            if (o_sib /= 1 .and. o_sib /= 2) then
                write (*,fmt="(a)")'bad value defined :'// trim(value) //'['//trim(this%OPT_SIBF)//']'
                call this%cli_print_help
                stop
            end if
        end if
    end function cli_get_family
    !!***

    !!****f* m_qtlmap_type_cli/cli_get_optimization
    !!  NAME
    !!    cli_get_optimization
    !!  DESCRIPTION
    !!    Get the optimization number id for LA analysis with optimized method
    !!  NOTES
    !!    see m_qtlmap_optimization
    !!  SOURCE
    function cli_get_optimization(this) result(o_optim)
        class(QTLMAP_CLI)  , intent(in)  :: this
        integer     :: o_optim
        character(len=LEN_DEF)   :: value

        o_optim = this%OPTIM_DEFAULT ! default value
        call this%get_key_value(this%OPT_OPTIM,value)
        if (trim(value)/='') then
            o_optim = get_int(value)
        end if
    end function cli_get_optimization
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_quiet
    !!  NAME
    !!    cli_is_quiet
    !!  DESCRIPTION
    !!    True if the user ask none information
    !!  SOURCE
    function cli_is_quiet(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_QUIET)
        if (.not. val) val = this%key_exist(this%OPT_QUIET_SHORT)
    end function cli_is_quiet
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_verbose
    !!  NAME
    !!    cli_is_verbose
    !!  DESCRIPTION
    !!    True if the user ask verbose information
    !!  SOURCE
    function cli_is_verbose(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_VERBOSE)
        if (.not. val) then
            val = this%key_exist(this%OPT_VERBOSE_SHORT)
        end if
    end function cli_is_verbose
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_print_all
    !!  NAME
    !!    cli_is_print_all
    !!  DESCRIPTION
    !!    True if the user ask to print all information
    !!  SOURCE
    function cli_is_print_all(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_PRINT_CLASSIC)
    end function cli_is_print_all
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_debug
    !!  NAME
    !!    cli_is_debug
    !!  DESCRIPTION
    !!    True if the user ask debug information
    !!  SOURCE
    function cli_is_debug(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_DEBUG)
        if (.not. val) val = this%key_exist(this%OPT_DEBUG_SHORT)
    end function cli_is_debug
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_help_panalyse
    !!  NAME
    !!    cli_is_help_panalyse
    !!  DESCRIPTION
    !!    True if the user ask help p_analyse information
    !!  SOURCE
    function  cli_is_help_panalyse(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_HELP_PANALYSE)

    end function  cli_is_help_panalyse
    !!***

    !!****f* m_qtlmap_type_cli/cli_is_help_paramsim
    !!  NAME
    !!    cli_is_help_paramsim
    !!  DESCRIPTION
    !!    True if the user ask help paramsim information
    !!  SOURCE
    function cli_is_help_paramsim(this) result(val)
        class(QTLMAP_CLI)  , intent(in)  :: this
        logical :: val
        val = this%key_exist(this%OPT_HELP_PARAMSIM)
    end function  cli_is_help_paramsim
    !!***

    !*******************************************************************************************************************************
    !*
    !*                               INTERNAL PROCEDURES


    !!****f* m_qtlmap_type_cli/get_nb_arg
    !!  NAME
    !!    get_nb_arg
    !!  DESCRIPTION
    !!    get the number of option give in the command line
    !!  SOURCE
    function get_nb_arg() result(nb)
        integer           :: nb

#ifdef AIX_SYS
         character(len=LEN_OPT) :: v
         v=''
         nb = 0
         call GETARG(nb+1,v)
         v = trim(v)
         do while (v/='')
           nb = nb+1
           call getarg(nb+1,v)
           v = trim(v)
         end do
#else
        nb = COMMAND_ARGUMENT_COUNT()
#endif
        return

    end function get_nb_arg
    !!***

    !!****f* m_qtlmap_type_cli/get_arg
    !!  NAME
    !!    get_arg
    !!  DESCRIPTION
    !!    get the argument value of the ith argument
    !!  SOURCE
    subroutine get_arg(it,value)
        integer , intent(in)                    :: it
        character(len=*),intent(out)      :: value

#ifdef AIX_SYS
        call GETARG(it,value)
#else
        call GET_COMMAND_ARGUMENT(it,value)
#endif
    end subroutine get_arg
    !!***

    !!****f* m_qtlmap_type_cli/get_binary_name
    !!  NAME
    !!    get_binary_name
    !!  DESCRIPTION
    !!    Get the name of the binary
    !!  SOURCE
    function get_binary_name() result(name)
        character(len=LEN_OPT) :: v
        character(len=LEN_DEF)   :: name
        integer                :: id

        call get_arg(0,v)
        name=trim(v)
        id=SCAN(name,'/')
        do while(id>0 .and. (id+1)<=LEN(name))
            name = trim(EXTRACT(name,id+1,LEN(name)))
            id=SCAN(name,'/')
        end do
        return
    end function get_binary_name
    !!***

    !!****f* m_qtlmap_type_cli/get_key_value
    !!  NAME
    !!    get_key_value
    !!  DESCRIPTION
    !!    Get value according to the key give in parameter
    !!  SOURCE
    subroutine get_key_value(this,keysearch,value)
        class(QTLMAP_CLI)  , intent(in)  :: this
        character(len=LEN_OPT) ,intent(in)   :: keysearch
        character(len=LEN_DEF) ,intent(out)    :: value

        integer                :: i,id
        character(len=LEN_OPT) :: v
        character(len=LEN_DEF)   :: option,key

        value=''

        do i=1,get_nb_arg()
            call get_arg(i,v)
            if (trim(v) == '') cycle
            option = v
            id=SCAN(option,'=')
            if (id<=0) cycle
            key = trim(EXTRACT(option,1,id-1))

            if (key == keysearch) then
                value = trim(EXTRACT(option,id+1,LEN(option)))
                return
            end if
        end do

    end subroutine get_key_value
    !!***

    !!****f* m_qtlmap_type_cli/get_key_list_values
    !!  NAME
    !!    get_key_list_values
    !!  DESCRIPTION
    !!    Get value according to the key give in parameter
    !!  SOURCE
    subroutine get_key_list_values(this,keysearch,n,values,narg)
        class(QTLMAP_CLI)  , intent(in)  :: this
        character(len=LEN_OPT) ,intent(in)   :: keysearch
        integer                ,intent(in)   :: n
        character(len=LEN_DEF) ,dimension(n),intent(out) :: values
        integer                ,intent(out)  :: narg

        integer                :: i,id
        character(len=LEN_OPT) :: v
        character(len=LEN_DEF)   :: option,key

        character(len=LEN_DEF) :: value

        call this%get_key_value(keysearch,value)
        id=1
        narg=0
        do while (id>0 .or. n == narg)
            if ( trim(value) == "") exit
            id=SCAN(value,',')
            narg=narg+1
            if ( id <= 0 ) then
                values(narg) = value
                exit
            else
                values(narg) = trim(EXTRACT(value,1,id-1))
                value = trim(EXTRACT(value,id+1,LEN(value)))
            end if
        end do

    end subroutine get_key_list_values
    !!***

    !!****f* m_qtlmap_type_cli/key_exist
    !!  NAME
    !!    key_exist
    !!  DESCRIPTION
    !!    True if the key exist in the argument command line
    !!  SOURCE
    function key_exist(this,keysearch) result(vexist)
        class(QTLMAP_CLI)  , intent(in)  :: this
        character(len=LEN_OPT) ,intent(in)   :: keysearch

        logical                :: vexist
        integer                :: i,id
        character(len=LEN_OPT) :: v
        character(len=LEN_DEF)   :: option,key

        vexist=.false.
        do i=1,get_nb_arg()
            call get_arg(i,v)
            if (trim(v) == '') cycle
            if (trim(v) == keysearch) then
                vexist = .true.
                return
            end if
            option = v
            id=SCAN(option,'=')
            if (id<=0) cycle
            key = trim(EXTRACT(option,1,id-1))
            if (key == keysearch) then
                vexist = .true.
                return
            end if
        end do
    end function key_exist
!!***


end module m_qtlmap_type_cli
