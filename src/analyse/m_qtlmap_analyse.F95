!!  NAME
!!    ANALYSE
!!  DESCRIPTION
!!  Package data :
!!     This package implements all analysis availabale on QTLMap.
!!  CREATION DATE
!!  01/01/2009
!!  COPYRIGHT
!!***

module m_qtlmap_analyse
    use m_qtlmap_log
    use m_qtlmap_math
    use m_qtlmap_output_handler
    use m_qtlmap_analyse_gen
    use m_qtlmap_incidence
    use m_qtlmap_types

    implicit none

    !USER INTERFACE
    !--------------
    public :: get_name_analyse
    public :: get_type_data_analyse
    public :: is_multitrait_analysis
    public :: need_normalize_data
    public :: analyse
    public :: quantile

contains

    !!   Returns the string name of a analysis identified by the calcul id.
    function get_name_analyse(opt_calcul) result (name)
        integer        , intent(in)                    :: opt_calcul
        character(len=300)                             :: name
        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT)
                name ="UNITRAIT ANALYSIS"
            CASE (ANALYSE_UNITRAIT_MODLIN)
                name = "MODLIN ANALYSIS"
            CASE (ANALYSE_MULTITRAIT)
                name = "MULTI ANALYSIS"
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HOMO)
                name = "LINEAR LA HOMOSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO)
                name = "LINEAR LA HETEROSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LD_HOMO)
                name = "LINEAR LD HOMOSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LD_HETERO)
                name = "LINEAR LD HETEROSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO)
                name = "LINEAR LDLA HOMOSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO)
                name = "LINEAR LDLA HETEROSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO)
                name = "LINEAR LA RELATIONSHIP MATRIX HOMOSCEDASTIC"
            CASE (ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO)
                name = "LINEAR LD RELATIONSHIP MATRIX HOMOSCEDASTIC"
            CASE( ANALYSE_UNITRAIT_LDLA)
                name = "LDLA ANALYSIS"
            CASE( ANALYSE_UNITRAIT_LDJH)
                name = "LDJH ANALYSIS"
            CASE( ANALYSE_UNITRAIT_LD)
                name = "LD ANALYSIS"
            CASE( ANALYSE_UNITRAIT_LA)
                name = "LA ANALYSIS (MODLIN)"
            CASE (ANALYSE_UNITRAIT_MODLIN_COX)
                name = "COX ANALYSIS"
            CASE (ANALYSE_MULTITRAIT_DA)
                name = "DISCRIMINANT ANALYSIS"
            CASE (ANALYSE_2QTL_INTERACTION)
                name = "2QTL INTERACTION ANALYSIS"
            CASE (ANALYSE_UNITRAIT_CONTINGENCE)
                name = "NEW MODLIN"
            CASE (ANALYSE_TRAIT_COV_CONTINGENCE)
                name = "TRAITS AS COVARIATES"
            CASE (ANALYSE_MULTITRAIT_INCIDENCE_LU)
                name = "ANALYSE_MULTITRAIT_INCIDENCE_LU"
            CASE (ANALYSE_MULTITRAIT_INCIDENCE)
                name = "MULTI-TRAITS ANALYSIS WITH COVARIATE/FIXED EFFECT"
            CASE (ANALYSE_DEV_1)
                name = "ANALYSE_DEV_1"
            CASE (ANALYSE_DEV_2)
                name = "ANALYSE_DEV_2"
            CASE (ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)
                name = "ANALYSE_UNITRAIT_LINEAR_RACE_HETERO"
            CASE (ANALYSE_TRAIT_BIALL_FARNIR)
                name = "ANALYSE_TRAIT_BIALL_FARNIR"
            CASE DEFAULT
                call stop_application("Unknown name description for analysis [OPT_CALCUL="//trim(str(opt_calcul))//"]")
        END SELECT
    end function get_name_analyse

    !!   Returns the data type which is manageg by the analysis identified by the calcul id
    function get_type_data_analyse(opt_calcul) result (typeData)
        integer        , intent(in)                    :: opt_calcul
        integer                                        :: typeData
        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT,ANALYSE_UNITRAIT_MODLIN,ANALYSE_MULTITRAIT,ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HETERO,ANALYSE_UNITRAIT_LDLA,ANALYSE_UNITRAIT_LA,&
                ANALYSE_UNITRAIT_LD,ANALYSE_UNITRAIT_LDJH, ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO, &
                ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_MULTITRAIT_DA,ANALYSE_2QTL_INTERACTION,ANALYSE_UNITRAIT_CONTINGENCE,&
                ANALYSE_TRAIT_COV_CONTINGENCE,ANALYSE_MULTITRAIT_INCIDENCE,ANALYSE_MULTITRAIT_INCIDENCE_LU,&
                ANALYSE_DEV_1,ANALYSE_DEV_2,ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO,ANALYSE_TRAIT_BIALL_FARNIR)
                typeData = TYPE_DATA_CONTINUE
            CASE (ANALYSE_UNITRAIT_MODLIN_COX)
                typeData = TYPE_DATA_COX
            CASE DEFAULT
                call stop_application("Unknown data type description for analysis [OPT_CALCUL="//trim(str(opt_calcul))//"]")
        END SELECT
    end function get_type_data_analyse


    !!   Returns true if the analysis is a multitrait method, false otherwise
    function is_multitrait_analysis(opt_calcul) result (is_multi)
        integer        , intent(in)                    :: opt_calcul
        logical                                        :: is_multi
        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT,ANALYSE_UNITRAIT_MODLIN,ANALYSE_UNITRAIT_MODLIN_COX,ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HETERO,ANALYSE_UNITRAIT_LDLA,ANALYSE_UNITRAIT_LA,ANALYSE_UNITRAIT_LD,&
                ANALYSE_UNITRAIT_LDJH,ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO, &
                ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_2QTL_INTERACTION,ANALYSE_UNITRAIT_CONTINGENCE,&
                ANALYSE_TRAIT_COV_CONTINGENCE,ANALYSE_DEV_1,ANALYSE_DEV_2,&
                ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO,ANALYSE_TRAIT_BIALL_FARNIR)
                is_multi=.false.
                return
            CASE (ANALYSE_MULTITRAIT,ANALYSE_MULTITRAIT_INCIDENCE,ANALYSE_MULTITRAIT_DA,ANALYSE_MULTITRAIT_INCIDENCE_LU)
                is_multi=.true.
                return
            CASE DEFAULT
                call stop_application("Unknown data type description for analysis [OPT_CALCUL="//trim(str(opt_calcul))//"]")
        END SELECT

    end function is_multitrait_analysis

    !!   Returns true if the analysis need to normalize data

    function need_normalize_data(opt_calcul) result (normalize)
        integer        , intent(in)                    :: opt_calcul
        logical                                        :: normalize
        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT,ANALYSE_UNITRAIT_MODLIN,ANALYSE_MULTITRAIT,ANALYSE_MULTITRAIT_DA,&
                ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HETERO,&
                ANALYSE_UNITRAIT_LDLA,ANALYSE_UNITRAIT_LA,ANALYSE_UNITRAIT_LD,&
                ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO, &
                ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,ANALYSE_UNITRAIT_LDJH,&
                ANALYSE_2QTL_INTERACTION,ANALYSE_UNITRAIT_CONTINGENCE,&
                ANALYSE_MULTITRAIT_INCIDENCE,ANALYSE_TRAIT_COV_CONTINGENCE,ANALYSE_MULTITRAIT_INCIDENCE_LU,ANALYSE_DEV_1,&
                ANALYSE_DEV_2,ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO,ANALYSE_TRAIT_BIALL_FARNIR)
                normalize=.true.
                return
            CASE (ANALYSE_UNITRAIT_MODLIN_COX)
                normalize=.false.
                return
            CASE DEFAULT
                call stop_application("need_normalize_data : Unknown data type description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT

    end function need_normalize_data

    subroutine set_kind_variance_analysis(opt_calcul,type,force)
        use m_qtlmap_incidence
        integer        , intent(in)                    :: opt_calcul
        integer        , intent(inout)                 :: type
        integer        , intent(in)   ,optional        :: force

        if (present(force)) then
            call set_parameter(force,type)
            return
        end if

        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,&
                ANALYSE_UNITRAIT_CONTINGENCE,ANALYSE_MULTITRAIT_INCIDENCE,ANALYSE_TRAIT_COV_CONTINGENCE,&
                ANALYSE_MULTITRAIT_INCIDENCE_LU,ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)
                call set_parameter(KIND_HETEROSCEDASTIC,type)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,ANALYSE_2QTL_INTERACTION)
                call set_parameter(KIND_HOMOSCEDASTIC,type)
            CASE DEFAULT
                call stop_application("set_kind_variance_analysis : Unknown type variance description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT
    end subroutine set_kind_variance_analysis

    subroutine set_model(opt_calcul,type,force)
        use m_qtlmap_incidence
        integer        , intent(in)                       :: opt_calcul
        integer        , intent(inout)                    :: type
        integer        , intent(in)   ,optional           :: force

        if (present(force)) then
            call set_parameter(force,type)
            return
        end if

        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,&
                ANALYSE_UNITRAIT_CONTINGENCE,ANALYSE_TRAIT_COV_CONTINGENCE,ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,ANALYSE_2QTL_INTERACTION,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)
                call set_parameter(MODEL_UNITRAIT_POLYGENIC,type)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO)
                call set_parameter(MODEL_UNITRAIT_RELATIONSHIP,type)
            CASE (ANALYSE_MULTITRAIT_INCIDENCE,ANALYSE_MULTITRAIT_INCIDENCE_LU)
                call set_parameter(MODEL_MULTITRAIT_POLYGENIC,type)
            CASE DEFAULT
                call stop_application("set_model : Unknown model type description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT
    end subroutine set_model

    subroutine set_resolution(opt_calcul,type,force)
        use m_qtlmap_incidence
        integer        , intent(in)                       :: opt_calcul
        integer        , intent(inout)                    :: type
        integer        , intent(in)   ,optional           :: force

        if (present(force)) then
            call set_parameter(force,type)
            return
        end if


        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO,ANALYSE_UNITRAIT_LINEAR_LD_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LD_HOMO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_2QTL_INTERACTION,ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)
                call set_parameter(RESOLUTION_LINEAR,type)
            CASE (ANALYSE_UNITRAIT_CONTINGENCE,ANALYSE_MULTITRAIT_INCIDENCE_LU,ANALYSE_MULTITRAIT_INCIDENCE,&
                ANALYSE_TRAIT_COV_CONTINGENCE)
                call set_parameter(RESOLUTION_OPTIMISATION,type)
            CASE DEFAULT
                call stop_application("set_resolution : Unknown resolution type description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT
    end subroutine set_resolution


    subroutine set_contingence_effect(opt_calcul,cont_eff)
        use m_qtlmap_incidence
        integer        , intent(in)                          :: opt_calcul
        type(type_effect_contingence)        , intent(inout) :: cont_eff

        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_MULTITRAIT_INCIDENCE_LU,&
                ANALYSE_MULTITRAIT_INCIDENCE,&
                ANALYSE_UNITRAIT_CONTINGENCE,&
                ANALYSE_2QTL_INTERACTION)

                cont_eff = CONTINGENCE_QTL

            CASE (ANALYSE_UNITRAIT_LINEAR_LD_HETERO,ANALYSE_UNITRAIT_LINEAR_LD_HOMO)

                cont_eff = CONTINGENCE_HAPLO

            CASE (ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO)

                cont_eff = CONTINGENCE_HAPLO_QTL

            CASE (ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO)

                cont_eff = CONTINGENCE_QTL_WITHOUT_POLY

            CASE (ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO)

                cont_eff = CONTINGENCE_HAPLO_WITHOUT_POLY

            CASE (ANALYSE_TRAIT_COV_CONTINGENCE)

                cont_eff = CONTINGENCE_QTL_TRAITS

            CASE(ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)

                cont_eff = CONTINGENCE_RACE

            CASE DEFAULT
                call stop_application("set_contingence_effect : Unknown contingence type description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT
    end subroutine set_contingence_effect


    subroutine check_qtl_compatibility_incidence(dataset)
        use m_qtlmap_incidence
        type(QTLMAP_DATASET)       ,intent(in) :: dataset
        integer :: opt_qtl,opt_calcul
        logical :: check

        opt_qtl = dataset%cli%cli_get_nqtl()
        opt_calcul = dataset%cli%cli_get_analyse()

        SELECT CASE (opt_calcul)
            CASE (ANALYSE_UNITRAIT_LINEAR_LA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_MULTITRAIT_INCIDENCE_LU,&
                ANALYSE_MULTITRAIT_INCIDENCE,&
                ANALYSE_UNITRAIT_CONTINGENCE,&
                ANALYSE_UNITRAIT_LINEAR_LD_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LD_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,&
                ANALYSE_TRAIT_COV_CONTINGENCE,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)

                check =  opt_qtl <=0

            CASE (ANALYSE_2QTL_INTERACTION)

                check =  opt_qtl /= 2

            CASE DEFAULT
                call stop_application("set_contingence_effect : Unknown contingence type description for analysis [OPT_CALCUL="//&
                    trim(str(opt_calcul))//"]")
        END SELECT

        if ( check ) then
            call stop_application("check_qtl_compatibility : Incompatibility between [OPT_CALCUL="//&
                trim(str(opt_calcul))//"] and [OPT_QTL="//trim(str(opt_qtl))//"]")
        end if

    end subroutine check_qtl_compatibility_incidence


    subroutine check_qtl_compatibility(dataset)
        use m_qtlmap_incidence
        type(QTLMAP_DATASET)       ,intent(in) :: dataset
        integer :: opt_calcul,opt_qtl,ic
        logical :: car_real

        opt_calcul = dataset%cli%cli_get_analyse()
        opt_qtl = dataset%cli%cli_get_nqtl()

        ! TEST SI ANALYSE EST VALABLE POUR L OPTION DISABLE_SIRE_QTL

        SELECT CASE (opt_calcul)
            case(ANALYSE_UNITRAIT_LINEAR_LA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LD_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LD_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO,&
                ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LA_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO,&
                ANALYSE_UNITRAIT_LINEAR_RACE_HETERO)

            continue
        case default
            if (dataset%cli%key_exist(dataset%cli%DISABLE_SIRE_QTL) ) then
                call stop_application("option "//trim(dataset%cli%DISABLE_SIRE_QTL)//&
                    " is not implemented for this analysis !")
            end if
    END SELECT

    if ( opt_qtl <= 0) then
        call stop_application("bad definition of opt_qtl ["//trim(str(opt_qtl))//"]")
    end if

    ! arret si l'utilisateur demande d'analyser que des variables entieres ou categorielles selon
    ! une option non disponible (opt_qtl=2, opt_calcul=ANALYSE_MULTITRAIT  ou ANALYSE_UNITRAIT_MODLIN_cox

    car_real=.false.

    do ic=1,dataset%phenoModel%ncar
        car_real= car_real .or. &
            ((dataset%phenoModel%natureY(ic) == 'r') .or.( dataset%phenoModel%natureY(ic) == 'a') )
    end do

    if(.not.car_real) then
        if (opt_calcul == ANALYSE_UNITRAIT_MODLIN_COX .or.&
            opt_calcul == ANALYSE_MULTITRAIT) then
            call stop_application("Analysis option is not available [OPT_CALCUL="//trim(str(opt_calcul))//"]"&
                //"[OPT_QTL="//trim(str(opt_qtl))//"] for theses kind of data.")
        end if
    end if

    if ( opt_qtl>1) then
        if (opt_calcul == ANALYSE_UNITRAIT_MODLIN_COX .or.&
            opt_calcul == ANALYSE_UNITRAIT_MODLIN .or. &
            opt_calcul == ANALYSE_MULTITRAIT  .or. &
            opt_calcul == ANALYSE_UNITRAIT_LD .or. &
            opt_calcul == ANALYSE_UNITRAIT_LDLA .or. &
            opt_calcul == ANALYSE_UNITRAIT_LA .or. &
            opt_calcul == ANALYSE_UNITRAIT_LDJH .or. &
            opt_calcul == ANALYSE_MULTITRAIT_DA ) then
            call stop_application("Analysis QTL option is not available [OPT_CALCUL="//trim(str(opt_calcul))//"]"&
                //"[OPT_QTL="//trim(str(opt_qtl))//"] ")
        end if
    end if

end subroutine check_qtl_compatibility


subroutine check_simulation_compatibility(dataset)
    use m_qtlmap_incidence
    type(QTLMAP_DATASET)       ,intent(in) :: dataset

    integer :: opt_calcul
    opt_calcul = dataset%cli%cli_get_analyse()

    if ( dataset%cli%cli_is_permute() .and. ( &
        opt_calcul == ANALYSE_UNITRAIT_LD .or. &
        opt_calcul == ANALYSE_UNITRAIT_LDLA .or. &
        opt_calcul == ANALYSE_UNITRAIT_LDJH .or. &
        opt_calcul == ANALYSE_UNITRAIT_LINEAR_LD_HETERO .or. &
        opt_calcul == ANALYSE_UNITRAIT_LINEAR_LD_HOMO .or. &
        opt_calcul == ANALYSE_UNITRAIT_LINEAR_LDLA_HETERO .or. &
        opt_calcul == ANALYSE_UNITRAIT_LINEAR_LDLA_HOMO .or. &
        opt_calcul == ANALYSE_UNITRAIT_LINEAR_LD_MA_HOMO) ) then

        call stop_application("You can not use the permutation option with a LD or LALDA analysis.")
    end if

end subroutine check_simulation_compatibility

!!   Main function to analyse a user dataset with a specific analysis.
!!
!!
!!  INPUTS
!!    opt_calcul : analysis id
!!    opt_qtl    : number of qtl hypothesis (analysis can be not compatible with this option)
!!    hdam       : LD option, take dam haplotypes in considration
!!    opt_sim    : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!!    rhoi       : residual correlation matrix finded for each hypothesis

subroutine analyse(dataset,spt,opt_calcul,opt_qtl,lrtsol,listincsol,rhoi,opt_sim)
    type(QTLMAP_DATASET)       ,intent(inout) :: dataset
    type(PDD_BUILD)            ,intent(inout) :: spt
    integer        , intent(in)                    :: opt_calcul,opt_qtl
    integer ,intent(in)                            :: opt_sim
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,opt_qtl)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,opt_qtl+1)   :: listincsol
    real (kind=dp) ,dimension(dataset%phenoModel%ncar,dataset%phenoModel%ncar), intent(out)  :: rhoi
    logical                        :: biq,interaction
    logical                        :: hdam
    logical                        :: car_real
    integer                        :: type_model = MASK_INIT_MOD_INC
    type(type_effect_contingence)  :: cont_eff

    hdam        = .not. dataset%cli%cli_is_no_hdam()
    biq         = dataset%cli%cli_is_biq()
    interaction = dataset%cli%cli_is_interaction()

    rhoi=0.d0

    call check_qtl_compatibility(dataset)

    if ( opt_calcul == ANALYSE_UNITRAIT .and. opt_qtl == 1) then
        call opti_unitrait(dataset,spt,lrtsol,listincsol,opt_sim,1)
    else if ( opt_calcul == ANALYSE_UNITRAIT .and. opt_qtl == 2 ) then
        call opti_unitrait(dataset,spt,lrtsol,listincsol,opt_sim,2)
    else if ( opt_calcul == ANALYSE_UNITRAIT_MODLIN .and. opt_qtl == 1) then
        call opti_unitrait_modlin(dataset,spt,lrtsol,listincsol,opt_sim)
        if (opt_sim /= TRANSCIPTOME_ANALYSE) call opti_unitrait_discret(dataset,spt,lrtsol,listincsol,opt_sim)
    else if ( opt_calcul == ANALYSE_UNITRAIT_LDLA .and. opt_qtl == 1) then
        call opti_unitrait_modlin_new(dataset,spt,lrtsol(:,1),listincsol,"hete","LDLA",hdam,opt_sim)
    else if ( opt_calcul == ANALYSE_UNITRAIT_LD .and. opt_qtl == 1) then
        call opti_unitrait_modlin_new(dataset,spt,lrtsol(:,1),listincsol,"hete","LD  ",hdam,opt_sim)
    else if ( opt_calcul == ANALYSE_UNITRAIT_LA .and. opt_qtl == 1) then
        call opti_unitrait_modlin_new(dataset,spt,lrtsol(:,1),listincsol,"hete","LA  ",hdam,opt_sim)
    else if ( opt_calcul == ANALYSE_UNITRAIT_LDJH .and. opt_qtl == 1) then
        call opti_unitrait_modlin_new(dataset,spt,lrtsol(:,1),listincsol,"hete","LDJH",hdam,opt_sim)
    else if ( opt_calcul == ANALYSE_MULTITRAIT_DA .and. opt_qtl == 1) then
        call opti_multitrait_DA(dataset,spt,lrtsol(1,1),opt_sim)
    else if ( opt_calcul == ANALYSE_UNITRAIT_MODLIN_cox .and. opt_qtl == 1) then
        call opti_unitrait_modlin_cox(dataset,spt,lrtsol,listincsol,opt_sim)
    else if ( opt_calcul == ANALYSE_MULTITRAIT .and. opt_qtl == 1) then
        call opti_multitraits(dataset,spt,lrtsol(1,1),listincsol,rhoi,opt_sim)
    else if ( opt_calcul == ANALYSE_TRAIT_BIALL_FARNIR ) then
        call opti_unitrait_biall_farnir(dataset,spt,lrtsol(:,1),listincsol,opt_sim,opt_qtl)
    !            else if ( opt_calcul == ANALYSE_2QTL_INTERACTION .and. opt_qtl == 2) then
    !               call opti_unitrait_interaction(lrtsol,opt_sim,biq)
    else if ( opt_calcul == ANALYSE_MULTITRAIT_INCIDENCE .and. opt_qtl >= 1) then
        call set_kind_variance_analysis(opt_calcul,type_model)
        call set_model(opt_calcul,type_model)
        call set_resolution(opt_calcul,type_model)
        call opti_multitrait_modlin(dataset,spt,lrtsol(1,:),listincsol,opt_sim,opt_qtl,type_model)
    else if ( opt_calcul == ANALYSE_MULTITRAIT_INCIDENCE_LU .and. opt_qtl >= 1) then

        call set_kind_variance_analysis(opt_calcul,type_model)
        call set_model(opt_calcul,type_model)
        call set_resolution(opt_calcul,type_model)
        call opti_multitrait_modlin(dataset,spt,lrtsol(1,:),listincsol,opt_sim,opt_qtl,type_model)
    else if ( opt_calcul == ANALYSE_DEV_1 .and. opt_qtl == 1 ) then

        call set_kind_variance_analysis(opt_calcul,type_model)
        call set_model(opt_calcul,type_model)
        call set_resolution(opt_calcul,type_model)
        call opti_unitrait_dev1(dataset,spt,lrtsol,listincsol,opt_sim,type_model,opt_qtl,hdam,biq)
    else
        ! qtl compatibility
        call check_qtl_compatibility_incidence(dataset)
        ! homoscedastic/heteroscedastic
        call set_kind_variance_analysis(opt_calcul,type_model)
        ! Polygenic / Multitrait Polygenic / Animal
        call set_model(opt_calcul,type_model)
        ! Linear / optimization
        call set_resolution(opt_calcul,type_model)
        ! kind of build depending the calculus
        call set_contingence_effect(opt_calcul,cont_eff)
        ! Computation
        call opti_unitrait_incidence&
            (dataset,spt,1,dataset%phenoModel%ncar,lrtsol,listincsol,opt_sim,cont_eff,type_model,interaction,opt_qtl,hdam,biq)
    end if

    if (opt_sim /= SIMULATION) call print_allelic_origin(dataset,spt)


end subroutine analyse
!!***

!!****f* m_qtlmap_analyse/analyse_cuda
!!  NAME
!!    analyse_cuda
!!  DESCRIPTION
!!   Main function to analyse a user dataset with a specific analysis.
!!
!!  INPUTS
!!    opt_calcul : analysis id
!!    opt_qtl    : number of qtl hypothesis (analysis can be not compatible with this option)
!!    hdam       : LD option, take dam haplotypes in considration
!!    opt_sim    : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!!    rhoi       : residual correlation matrix finded for each hypothesis
!! SOURCE

#ifdef MANAGE_CUDA
    subroutine analyse_cuda(dataset,spt,nsim,opt_sim,opt_calcul,opt_qtl,lrtsol,listincsol,listLrtSolSim,ySIMUL)
          type(QTLMAP_DATASET)       ,intent(inout) :: dataset
          type(PDD_BUILD)            ,intent(inout) :: spt
          integer        , intent(in)                    :: nsim,opt_calcul,opt_sim,opt_qtl
          type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,opt_qtl)         :: lrtsol
          type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,opt_qtl+1) :: listincsol
          type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,opt_qtl,nsim)    :: listLrtSolSim
          real (kind=dp), dimension (:,:,:),allocatable  , intent(in) :: ySIMUL
          logical                        :: hdam
          logical                        :: biq,interaction
          integer                        :: type_model = MASK_INIT_MOD_INC
          type(type_effect_contingence)  :: cont_eff

          hdam        = dataset%cli%cli_is_no_hdam()
          biq         = dataset%cli%cli_is_biq()
          interaction = dataset%cli%cli_is_interaction()

          if ( opt_qtl <= 0) then
             call stop_application("bad definition of opt_qtl ["//trim(str(opt_qtl))//"]")
          end if

    ! qtl compatibility
          call check_qtl_compatibility_incidence(dataset)
    ! homoscedastic/heteroscedastic
          call set_kind_variance_analysis(opt_calcul,type_model)
    ! Polygenic / Multitrait Polygenic / Animal
          call set_model(opt_calcul,type_model)
    ! Linear / optimization
          call set_resolution(opt_calcul,type_model)
    ! kind of build depending the calculus
          call set_contingence_effect(opt_calcul,cont_eff)
    ! Computation
          call opti_unitrait_incidence_cuda&
               (dataset,spt,nsim,lrtsol,listincsol,listLrtSolSim,opt_sim,cont_eff,interaction,&
                type_model,opt_qtl,ySIMUL,hdam,biq)

        call print_allelic_origin(dataset,spt)

    end subroutine analyse_cuda
#endif

!!***


!!****f* m_qtlmap_analyse/opti_unitrait
!!  NAME
!!    opti_unitrait
!!  DESCRIPTION
!!   LA single trait with pre-corrected data (1 or 2 qtl to test)
!!  INPUTS
!!    nbqtl      : number of qtl to test. 1 <= nbqtl <= 2
!!    opt_sim    : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_unitrait(dataset,spt,lrtsol,listincsol,opt_sim,nbqtl)
    use m_qtlmap_analyse_unitrait
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    integer                , intent(in)                                       :: nbqtl
    type(TYPE_LRT_SOLUTION)         , intent(inout) ,dimension(dataset%phenoModel%ncar,nbqtl)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,  intent(inout) ,dimension(dataset%phenoModel%ncar,nbqtl+1)  :: listincsol
    integer                        ,  intent(in)                              :: opt_sim

    integer                                           :: i,ic

    if ( nbqtl > 2 .or. nbqtl < 1 ) then
        call stop_application("Devel error : Module UNITRAIT can not perform analysis with nqtl:"//trim(str(nbqtl)))
    end if

    !$OMP PARALLEL DEFAULT(SHARED) IF (dataset%phenoModel%ncar>1)
    call init_analyse_gen(dataset)
    call init_analyse_unitrait_1QTL(dataset)
    if ( nbqtl == 2) call init_analyse_unitrait_2QTL(dataset)

    !$OMP DO
    do ic=1,dataset%phenoModel%ncar
        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) then
            call stop_application("This analyse does not manage this nature of data! TRAIT["//trim(dataset%phenoModel%carac(ic))&
                //"] nature:"//trim(dataset%phenoModel%natureY(ic)))
        end if

        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call optinit(dataset,ic)
        call opti_0qtl(dataset,ic)
        call opti_1qtl(dataset,spt,ic,lrtsol(ic,1))

        if ( opt_sim /= SIMULATION) then
            call set_solution_hypothesis0(dataset,ic,listincsol(ic,1))
            call set_solution_hypothesis1(dataset,spt,ic,listincsol(ic,2))
        end if

        if ( nbqtl == 2 ) then
            call opti_1car_2qtl(dataset,spt,ic,lrtsol(ic,2))
            if ( opt_sim /= SIMULATION) then
                call set_solution_hypothesis2(dataset,spt,ic,listincsol(ic,3))
            end if
        end if
    end do
    !$OMP END DO
    call end_analyse_unitrait_1QTL
    if ( nbqtl == 2) call end_analyse_unitrait_2QTL
    call end_analyse_gen
    !$OMP END PARALLEL

    if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
        do ic=1,dataset%phenoModel%ncar
            if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            call print_incidence_solution(dataset,listincsol(ic,1))
            call print_lrt_solution(dataset,lrtsol(ic,1))
            call print_incidence_solution(dataset,listincsol(ic,2))
            call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
            call print_LRT(dataset,ic,lrtsol(ic,1))
            call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))

            if ( nbqtl == 2 ) then
                call print_lrt_solution(dataset,lrtsol(ic,2))
                call print_incidence_solution(dataset,listincsol(ic,3))
              !  call create_grid_file_2QTL(dataset,lrtsol(ic,2))
                call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
                call print_LRT(dataset,ic,lrtsol(ic,2))
                call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,2))
            end if
        end do
    end if


    ! print txt transcriptome format
    if (opt_sim == TRANSCIPTOME_ANALYSE) then
        call  print_transcriptome_H0(dataset,listincsol(:,1))
        call  print_transcriptome(dataset,lrtsol(:,1),listincsol(:,2))
        if ( nbqtl == 2) call  print_transcriptome(dataset,lrtsol(:,2),listincsol(:,3))
    end if


    if ( opt_sim /= SIMULATION) then
        call print_summary_analyse(dataset,lrtsol(:,1),listincsol(:,2),1)
        if ( nbqtl == 2) call print_summary_analyse(dataset,lrtsol(:,2),listincsol(:,3),2)
    end if
end subroutine opti_unitrait
!!***

!!****f* m_qtlmap_analyse/opti_unitrait_incidence
!!  NAME
!!    opti_unitrait_incidence
!!  DESCRIPTION
!!   LA / LD / LDLA single trait with model description with optimised or linearised likelihood.
!!  INPUTS
!!    starticar      : start index
!!    endicar        : end index
!!    type_incidence : type of incidence matrix construction (INCIDENCE_TYPE_CLASSIC, INCIDENCE_TYPE_LD, INCIDENCE_TYPE_LDLA, INCIDENCE_TRAIT_COV)
!!    type_resolution: type of resolution (LINEAR_HOMOSCEDASTIC, LINEAR_HETEROSCEDASTIC, OPTIM_HETEROSCEDASTIC )
!!    interaction    : the last test , test if there are interaction betwwen qtl (qtl >= 2)
!!    nqtl           : number of qtl to test. nqtl >= 1
!!    hdam           : LD option, take dam haplotypes in considration
!!    opt_sim        : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
!   subroutine opti_unitrait_incidence(lrtsol,listincsol,opt_sim,type_incidence,type_resolution,nqtl,hdam,nopoly)
subroutine opti_unitrait_incidence(dataset,spt,starticar,endicar,lrtsol,listincsol,opt_sim,contingence_effects,&
    type_resolution,interaction,nqtl,hdam,biq,startTraitCov,endTraitCov)
    use m_qtlmap_incidence
    use m_qtlmap_incidence_optim
    use m_qtlmap_incidence_linear
    type(QTLMAP_DATASET)       ,intent(inout) :: dataset
    type(PDD_BUILD)            ,intent(inout) :: spt
    integer                  ,intent(in)                                   :: nqtl,starticar,endicar
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,nqtl)         :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,nqtl+1) :: listincsol
    !    logical                        ,intent(in)                             :: hdam,nopoly
    logical                        ,intent(in)                             :: hdam,biq,interaction
    type(type_effect_contingence) , intent(in)                             :: contingence_effects
    integer                  ,intent(in)                                   :: opt_sim,type_resolution
    integer                  ,optional,intent(in)                          :: startTraitCov,endTraitCov
    type(TYPE_INCIDENCE_SOLUTION)  ,dimension(dataset%phenoModel%ncar)                        :: incsolInteraction
    type(TYPE_LRT_SOLUTION)  , dimension(dataset%phenoModel%ncar)                             :: lrtsolInteraction
    integer                               :: ic,i,j,ip,jm,stc,etc

    type(INCIDENCE_GEN_STRUCT)  ,dimension(:) ,allocatable :: workstruct
    logical                                                :: islinear,check
    integer :: MaxHypothesis

    !check if interaction=> can not estime dam
    if (count(dataset%phenoAnimal%estime(:,:))> 0 .and. interaction ) then
        call stop_application("Interaction is only available for half-sib analysis. increase ndmin key.")
    end if

    !start and end trait as cov
    stc=0
    etc=0
    if ( present(startTraitCov) )  stc = startTraitCov
    if ( present(endTraitCov) )  etc = endTraitCov

    if ( nqtl < 2 .and. interaction) then
        call stop_application("can not test an interaction with one qtl.")
    end if

    if ( starticar < 1 ) then
        call stop_application("Dev error [opti_unitrait_incidence] bad def of starticar :"//trim(str(starticar)))
    end if

    if ( endicar > dataset%phenoModel%ncar ) then
        call stop_application("Dev error [opti_unitrait_incidence] bad def of endicar :"//trim(str(endicar)))
    end if

    allocate (workstruct(starticar:endicar))

    islinear = .false.
    if ( is_parameter(RESOLUTION_LINEAR,type_resolution)) then
        islinear = .true.
        call init_analyse_linear(dataset,spt)
    end if

    MaxHypothesis = nqtl
    if ( interaction ) MaxHypothesis = MaxHypothesis + 1

    !$OMP PARALLEL DEFAULT(SHARED) IF (dataset%phenoModel%ncar>100)
    call init_analyse_gen(dataset)
    !$OMP DO
    do ic=starticar,endicar
        if ( .not. is_parameter(RESOLUTION_OPTIMISATION,type_resolution)) then
            check = dataset%phenoModel%natureY(ic) == 'i'.or. dataset%phenoModel%natureY(ic) == 'a'
        else
            check = dataset%phenoModel%natureY(ic) == 'a'
        end if

        if ( check ) then
            call stop_application("This analyse does not manage this nature of data! TRAIT["//trim(dataset%phenoModel%carac(ic))&
                //"] nature:"//trim(dataset%phenoModel%natureY(ic)))
        end if

        workstruct(ic)%hdam = hdam
        workstruct(ic)%biq  = biq
        workstruct(ic)%startTraitAsCov=stc
        workstruct(ic)%endTraitAsCov=etc

        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call init_workstruct(dataset,ic,workstruct(ic),contingence_effects,type_resolution,nqtl,MaxHypothesis,&
            (opt_sim == COMMON_ANALYSE).and.(contingence_effects%qtl >0) ,(opt_sim == COMMON_ANALYSE),.false.)

        if (islinear) then
            call opti_0qtl( dataset, ic, listincsol(ic,1),workstruct(ic),model_lin_h0)
        else
            call opti_0qtl( dataset, ic, listincsol(ic,1),workstruct(ic),model_optim_h0)
        end if

        do i=1,nqtl
            workstruct(ic)%nqtl=i
            workstruct(ic)%hypothesis=i
            if (islinear) then
                call opti_nqtl(dataset,spt,i,i,ic,workstruct(ic),listincsol(ic,1+i),lrtsol(ic,i),model_lin_hn)
            else
                call opti_nqtl(dataset,spt,i,i,ic,workstruct(ic),listincsol(ic,1+i),lrtsol(ic,i),model_optim_hn)
            end if
        end do ! iqtl

        if ( interaction) then
            !the last test is the interaction
            workstruct(ic)%hypothesis = workstruct(ic)%hypothesis + 1
            workstruct(ic)%effects    = CONTINGENCE_QTL_INTER
            if (islinear) then
                call opti_nqtl(dataset,spt,nqtl,nqtl+1,ic,workstruct(ic),&
                    incsolInteraction(ic),lrtsolInteraction(ic),model_lin_hn)
            else
                call opti_nqtl(dataset,spt,nqtl,nqtl+1,ic,workstruct(ic),&
                    incsolInteraction(ic),lrtsolInteraction(ic),model_optim_hn)
            end if
            print *,lrtsolInteraction(ic)%lrtmax
        end if
    end do ! ic
    !$OMP END DO
    call end_analyse_gen
    !$OMP END PARALLEL

    if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
        do ic=starticar,endicar
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            call print_incidence_solution(dataset,listincsol(ic,1))
            do i=1,nqtl
                call print_lrt_solution(dataset,lrtsol(ic,i))
                call print_incidence_solution(dataset,listincsol(ic,1+i))
                call print_confusion(dataset,i,workstruct(ic)%nalert(i),workstruct(ic)%alertQtl,workstruct(ic)%corrmax(i))
                if ( workstruct(ic)%ntest /= 0 )  then
                    call print_test_nuisances(workstruct(ic)%ntest,workstruct(ic)%listtestnuis)
                end if
                call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,i))
                if ( i==1 ) then
                    call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
                end if

                call print_LRT(dataset,ic,lrtsol(ic,i))

                if ( i==2 ) then
                   ! call create_grid_file_2QTL(dataset,lrtsol(ic,2))
                    call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
                   ! call print_LRT_2QTL(dataset,lrtsol(ic,2))
                end if
            end do ! iqtl
            if ( interaction ) then
                call print_lrt_solution(dataset,lrtsolInteraction(ic))
                call print_incidence_solution(dataset,incsolInteraction(ic))
                call print_confusion(dataset,nqtl+1,workstruct(ic)%nalert(nqtl+1),&
                    workstruct(ic)%alertQtl,workstruct(ic)%corrmax(nqtl+1))
                call lrtsolInteraction(ic)%release()
                call incsolInteraction(ic)%release()
            end if
        end do
    end if

    if ( opt_sim /= SIMULATION .and. (contingence_effects%qtl >0) ) then

        do i=1,nqtl
            call print_summary_analyse(dataset,lrtsol(:,i),listincsol(:,i+1),i,starticar,endicar)
        end do
    end if

     ! print txt transcriptome format
    if (opt_sim == TRANSCIPTOME_ANALYSE) then
        call  print_transcriptome_H0(dataset,listincsol(:,1))
        do i=1,nqtl
            call  print_transcriptome(dataset,lrtsol(:,i),listincsol(:,i+1))
        end do
    end if

    if ( islinear ) call end_analyse_linear

    do ic=starticar,endicar
        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
        call release_ws(workstruct(ic))
    end do

    deallocate (workstruct)

end subroutine opti_unitrait_incidence
!!***


subroutine opti_unitrait_dev1(dataset,spt,lrtsol,listincsol,opt_sim,type_resolution,nqtl,hdam,biq)
    use m_qtlmap_incidence
    type(QTLMAP_DATASET)       ,intent(inout) :: dataset
    type(PDD_BUILD)            ,intent(inout) :: spt
    integer                  ,intent(in)                                   :: nqtl
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,nqtl)         :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,nqtl+1) :: listincsol
    logical                        ,intent(in)                             :: hdam,biq
    integer                   ,intent(in)                                  :: opt_sim,type_resolution

    !local
    integer                  :: i,icar,ncarY,ncarX,bufNcar,opt_sim2,iq,ic,ip,ieff,ios,ic2,kd,jm,chr,n,ig
    real(kind=dp)                                       :: pptt,numerator,denominator
    real(kind=dp) ,dimension(:,:),allocatable           :: ahat,bhat
    real(kind=dp) ,dimension(:),allocatable             :: ysave,myprbp,tposA,LRT
    logical                                :: ok
    character(len=LEN_L)     :: out
    real(kind=dp)    ,dimension(:,:),pointer :: y

    y => dataset%phenoAnimal%y

    call log_mess(" ****************** XIAO DEV **************************************")

    !Number of performance Y (quantitative trait)
    ncarY=get_real(dataset%params%tabDevKey(1),ok)
    if ( .not. ok ) then
        call stop_application("bad init number of Y : ["//trim(dataset%params%tabDevKey(1))//"]")
    end if

    call log_mess("number of Y :"//trim(str(ncarY)),INFO_DEF)

    allocate (ahat(ncarY,dataset%genea%np))
    allocate (bhat(ncarY:dataset%phenoModel%ncar,dataset%genea%np))
    allocate (tposA(ncarY))
    allocate (LRT(dataset%phenoModel%ncar))
    ! output file
    out = dataset%params%tabDevKey(2)
    open(unit=998,file=out)

    i=3
    do icar=1,ncarY
        call log_mess(" =======>TRAIT "//trim(str(icar)),INFO_DEF)
        !ahat(icar) = get_real(tabDevKey(i),ok) fmt="("//trim(str(dataset%genea%np))//"a)"
        read (dataset%params%tabDevKey(i),*,iostat=ios) ahat(icar,:)
        read (dataset%params%tabDevKey(i),*,iostat=ios) tposA(icar)
        do ip=1,dataset%genea%np
            call log_mess("      =======>   SIRE "//trim(str(ip))//&
                " QTL effect :"//trim(str(ahat(icar,ip))),INFO_DEF)
        end do
        i=i+1
    end do

    opt_sim2 = COMMON_ANALYSE
    !step 1 ** Analyse unitrait for each trait **
    call log_mess(" ****************** STEP1 Analyse Y **************************************")
    call opti_unitrait_incidence&
        (dataset,spt,1,ncarY,lrtsol,listincsol,opt_sim2,&
        CONTINGENCE_QTL,type_resolution,.false.,nqtl,hdam,biq)

    do ic=1,ncarY
        LRT(ic) = lrtsol(ic,1)%lrtmax(0)
    end do

    call log_mess(" ****************** STEP2 Analyse X **************************************")
    opt_sim2 = TRANSCIPTOME_ANALYSE
    call opti_unitrait_incidence&
        (dataset,spt,ncarY+1,dataset%phenoModel%ncar,lrtsol,listincsol,&
        opt_sim2,CONTINGENCE_QTL,type_resolution,.false.,nqtl,hdam,biq)

    iq=1

    !on recupere tout les effets qtl des genes etudie pour chaque famille de pere
    !----------------------------------------------------------------------------
    write ( 998, * ) 'car , QTL effect Pour chaque famille de pere (non standardise)'
    do ic=ncarY+1,dataset%phenoModel%ncar
        LRT(ic) = lrtsol(ic,1)%lrtmax(0)
        if ( .not. associated(listincsol(ic,2)%qtl_groupeName) ) then
            call stop_application("Devel error : [opti_unitrait_dev1]**"//&
            " qtl_groupeName in TYPE_INCIDENCE_SOLUTION is not initialized **")
        end if
        ieff=listincsol(ic,2)%qtl_groupeName(1,iq)
        !         print *,"ieff:",ieff,listincsol(ic,2)%nbParameterGroup(ieff),listincsol(ic,2)%paramaterValue(ieff,:)
        do ip=1,dataset%genea%np
            bhat(ic,ip) = listincsol(ic,2)%paramaterValue(ieff,ip)
        end do
        write ( 998, * ) ic,bhat(ic,:)
    end do
    write ( 998, * )

    ! pour chaque caractere Y et pour chaque X on analyse avec X en covariable
    !---------------------------------------------------------------------------
    call log_mess(" ****************** STEP3 Analyse Y avec covariable X **************************************")
    write (998,*) "                    ****                          "
    write (998,*) " Analyse Y avec covariable X  "
    opt_sim2=COMMON_ANALYSE!SIMULATION
    do ic=1,ncarY
        do ic2=ncarY+1,dataset%phenoModel%ncar
            call opti_unitrait_incidence&
                (dataset,spt,ic,ic,lrtsol,listincsol,opt_sim2,&
                CONTINGENCE_QTL_TRAITS,type_resolution,.false.,nqtl,hdam,biq,ic2,ic2)
            write (998,*) lrtsol(ic,1)%chrmax(0) , dataset%map%absi(lrtsol(ic,1)%chrmax(0),lrtsol(ic,1)%nxmax(0)), &
                lrtsol(ic,1)%lrtmax(0)
        end do
    end do


    allocate (ysave(size(y,2)))

    call log_mess(" ****************** STEP4 Analyse Y avec Y - COV(Y,X)/VAR(X) X **************************************")
    write ( 998, * )
    write (998,*) "                    ****                          "
    write (998,*) " Analyse Y avec covariable Y - COV(Y,X)/VAR(X) X  "
    do ic=1,ncarY
        ysave = y(ic,:)
        !
        !       COV(X,Y) = [ SUM i ( X(i)*Y(i) ) / N ]  - [ [ SUM i X(i)  * SUM i Y(i)  ] / N^2 ]
        !
        !
        do ic2=ncarY+1,dataset%phenoModel%ncar
            !cov(X,Y)
            numerator = sum ( (/((y(ic,i)*y(ic2,i)),i=1,size(y,2))/) )/real(size(y,2))
            numerator = numerator - ( sum(y(ic,:))/real(size(y,2))) * ( sum(y(ic2,:))/real(size(y,2)))
            !var(X)
            denominator = sum ( (/((y(ic2,i)*y(ic2,i)),i=1,size(y,2))/) )/real(size(y,2))
            denominator = denominator - ( (sum(y(ic2,:))**2 )/real(size(y,2))**2)

            y(ic,:) = y(ic,:) - (numerator/denominator)*y(ic2,:)
            call opti_unitrait_incidence&
                (dataset,spt,ic,ic,lrtsol,listincsol,opt_sim,CONTINGENCE_QTL,type_resolution,.false.,nqtl,hdam,biq)
            write (998,*) lrtsol(ic,1)%chrmax(0) , dataset%map%absi(lrtsol(ic,1)%chrmax(0),lrtsol(ic,1)%nxmax(0)), &
                lrtsol(ic,1)%lrtmax(0)
            y(ic,:) = ysave
        end do !ic2
    end do ! ic


    call log_mess(" ****************** STEP5 Analyse Y avec Yi - Ai/Bi Xi **************************************")
    write ( 998, * )
    write (998,*) "                    ****                          "
    write (998,*) " Analyse Y avec covariable Yi - Ai/Bi Xi  - ",&
        "Ai correspond a l effet qtl (pour Y) du pere i et Bi effet QTL du pere i (pour X)   "
    do ic=1,ncarY
        ysave = y(ic,:)

        do ic2=ncarY+1,dataset%phenoModel%ncar
            !correction
            do ip=1,dataset%genea%np
                do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
                    do kd=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
                        y(ic,kd) = y(ic,kd) - (ahat(ic,ip)/bhat(ic2,ip))*y(ic2,kd)
                     ! print *,ahat(ic,ip)/bhat(ic2,ip),ahat(ic,ip),bhat(ic2,ip),y(ic,kd),ysave(kd)
                    end do
                end do
            end do

            call opti_unitrait_incidence&
                (dataset,spt,ic,ic,lrtsol,listincsol,opt_sim,CONTINGENCE_QTL,type_resolution,.false.,nqtl,hdam,biq)
            write (998,*) lrtsol(ic,1)%chrmax(0) , dataset%map%absi(lrtsol(ic,1)%chrmax(0),lrtsol(ic,1)%nxmax(0)),&
                lrtsol(ic,1)%lrtmax(0)
            y(ic,:) = ysave
        end do !ic2
    end do ! ic
    !stop
    call log_mess(" ****************** STEP6 Analyse Y avec  Y - B * pdd(t0) **************************************")

    write ( 998, * )
    write (998,*) "                    ****                          "
    write (998,*) " Analyse Y avec Y - B * pdd(t0)     ,",&
        "t0 position du maximum de vraissemblance pour Y (donné par l utilisateur)  "

    allocate (myprbp(size(y,2)))
    do ic=1,ncarY
        ysave = y(ic,:)
        chr=1
        n=1
        do while ( tposA(ic) > dataset%map%absi(chr,n) )
            n=n+1
        end do
        print *,"pos :",n

        do ic2=ncarY+1,dataset%phenoModel%ncar
            !correction
            do ip=1,dataset%genea%np
                do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
                    myprbp=0.d0
                    do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
                        do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
                            pptt=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                            !    print *,pptt,pdd(chr,kd,1,n),pdd(chr,kd,2,n),pdd(chr,kd,3,n),pdd(chr,kd,4,n)
                            myprbp(spt%ndesc(chr,kd))=myprbp(spt%ndesc(chr,kd))+spt%probg(chr,ig)*pptt
                        end do
                    end do
                    do kd=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
                        y(ic,kd) = y(ic,kd) - (bhat(ic2,ip))*myprbp(kd)
                    !    print *,y(ic,kd),ysave(kd),myprbp(kd)
                    end do
                 !  stop
                end do
            end do

            call opti_unitrait_incidence&
                (dataset,spt,ic,ic,lrtsol,listincsol,opt_sim,CONTINGENCE_QTL,type_resolution,.false.,nqtl,hdam,biq)
            write (998,*) lrtsol(ic,1)%chrmax(0) , dataset%map%absi(lrtsol(ic,1)%chrmax(0),lrtsol(ic,1)%nxmax(0)),&
                lrtsol(ic,1)%lrtmax(0)
            y(ic,:) = ysave
        end do !ic2
    end do ! ic

    deallocate (myprbp)
    deallocate (ysave)

    write (998,*)
    write (998,*) " DELTA Y, Y-Beta*X"
    write (998,*)
    write (998,*) "indexY indexX Delta"
    do ic=1,ncarY
        do ic2=ncarY+1,dataset%phenoModel%ncar
            write (998,*) ic,ic2,(LRT(ic2)-LRT(ic))
        end do
    end do


    close (998)
    deallocate (ahat)
    deallocate (bhat)
    deallocate (tposA)
    deallocate (LRT)

end subroutine opti_unitrait_dev1

!!****f* m_qtlmap_analyse/opti_unitrait_incidence_cuda
!!  NAME
!!    opti_unitrait_incidence_cuda
!!  DESCRIPTION
!!   LA / LD / LDLA single trait with model description with optimised or linearised likelihood.
!!  INPUTS
!!    type_incidence : type of incidence matrix construction (INCIDENCE_TYPE_CLASSIC, INCIDENCE_TYPE_LD, INCIDENCE_TYPE_LDLA, INCIDENCE_TRAIT_COV)
!!    type_resolution: type of resolution (LINEAR_HOMOSCEDASTIC, LINEAR_HETEROSCEDASTIC, OPTIM_HETEROSCEDASTIC )
!!    nqtl           : number of qtl to test. nqtl >= 1
!!    hdam           : LD option, take dam haplotypes in considration
!!    opt_sim        : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE

#ifdef MANAGE_CUDA
   subroutine opti_unitrait_incidence_cuda(dataset,spt,nsim,lrtsol,listincsol,listLrtSolSim,opt_sim,&
           contingence_effects,interaction,&
           type_resolution,nqtl,ySIMUL,hdam,biq)
       use m_qtlmap_incidence_linear
       type(QTLMAP_DATASET)       ,intent(inout) :: dataset
       type(PDD_BUILD)            ,intent(inout) :: spt
       integer                  ,intent(in)                                   :: nqtl,opt_sim,nsim
       type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,nqtl)         :: lrtsol
       type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,nqtl+1) :: listincsol
       type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,nqtl,nsim)    :: listLrtSolSim
       real (kind=dp), dimension (nsim,dataset%phenoModel%ncar,dataset%genea%nd)    , intent(in) :: ySIMUL
       logical                        ,intent(in)                             :: hdam
       logical                        ,intent(in)                             :: biq,interaction
       type(type_effect_contingence) , intent(in)                             :: contingence_effects
       integer                  ,intent(in)                                   :: type_resolution

       integer                               :: ic,i,j,ip,jm,kd,jj
       real (kind=dp), dimension (:,:)    ,allocatable                    :: AllsigQtlMoinUn

       type(INCIDENCE_GEN_STRUCT)  ,dimension(:) ,allocatable :: workstruct
       logical                                                :: islinear
       integer                                                :: nkd,sizeF(dataset%genea%np),MaxHypothesis
       real (kind=dp), dimension (:,:),allocatable :: YCOMPACT
       type(TYPE_INCIDENCE_SOLUTION)  ,dimension(dataset%phenoModel%ncar)                        :: incsolInteraction
       type(TYPE_LRT_SOLUTION)  , dimension(dataset%phenoModel%ncar)                             :: lrtsolInteraction

   !check if interaction=> can not estime dam
       if (count(dataset%phenoAnimal%estime(:,:))> 0 .and. interaction ) then
          call stop_application("Interaction is only available for half-sib analysis. increase ndmin key.")
       end if

       if ( nqtl < 2 .and. interaction) then
         call stop_application("can not test an interaction with one qtl.")
       end if

       allocate (YCOMPACT(nsim+1,dataset%genea%nd))
       allocate (workstruct(dataset%phenoModel%ncar))
       allocate (AllsigQtlMoinUn(nsim+1,dataset%genea%np))

       do ic=1,dataset%phenoModel%ncar
         if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) then
             call stop_application("This analyse does not manage this nature of data! TRAIT["//trim(dataset%phenoModel%carac(ic))&
              //"] nature:"//trim(dataset%phenoModel%natureY(ic)))
         end if
       end do

       MaxHypothesis = nqtl
       if ( interaction ) MaxHypothesis = MaxHypothesis + 1

       call init_analyse_linear(dataset,spt)
       call init_analyse_gen(dataset)
       do ic=1,dataset%phenoModel%ncar
            YCOMPACT=0.d0
            AllsigQtlMoinUn=0.d0
            nkd=0
            do ip=1,dataset%genea%np
             do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
              do kd=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
                if (dataset%phenoAnimal%presentc(ic,kd)) then
                  nkd=nkd+1
                  YCOMPACT(1,nkd) = dataset%phenoAnimal%y(ic,kd)
                  YCOMPACT(2:nsim+1,nkd)= ySIMUL(:,ic,kd)
             !YCOMPACT(2:nsim,nkd)= y(ic,kd) ! pour les test
                end if
              end do !kd
             end do ! jm
            end do ! ip

           workstruct(ic)%hdam = hdam
           workstruct(ic)%biq  = biq
           call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
           call init_workstruct(dataset,ic,workstruct(ic),contingence_effects,type_resolution,nqtl,MaxHypothesis,&
            (contingence_effects%qtl > 0) ,.true.,.false.)

           workstruct(ic)%hypothesis=0
           call opti_0qtl_cuda( dataset,nsim,ic, listincsol(ic,1),workstruct(ic),&
            YCOMPACT,sizeF,AllsigQtlMoinUn,model_lin_h0)

           do i=1,nqtl
             workstruct(ic)%nqtl=i
             workstruct(ic)%hypothesis=i
             call opti_nqtl_cuda(dataset,spt,nsim,i,ic,workstruct(ic),listincsol(ic,1+i),lrtsol(ic,i),&
                                 YCOMPACT,nqtl,listLrtSolSim,sizeF,AllsigQtlMoinUn,model_lin_hn)
           end do ! iqtl

           if ( interaction) then
         !the last test is the interaction
             workstruct(ic)%hypothesis = workstruct(ic)%hypothesis + 1
             workstruct(ic)%effects    = CONTINGENCE_QTL_INTER
             call opti_nqtl_cuda(dataset,spt,nsim,nqtl,ic,workstruct(ic),&
                                 incsolInteraction(ic),lrtsolInteraction(ic),&
                                 YCOMPACT,nqtl+1,listLrtSolSim,sizeF,AllsigQtlMoinUn,model_lin_hn)
           end if

       end do ! ic
       call end_analyse_gen

       call release_allocated_internal_sig

       if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
        do ic=1,dataset%phenoModel%ncar
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            call print_incidence_solution(dataset,listincsol(ic,1))
            do i=1,nqtl
              call print_lrt_solution(dataset,lrtsol(ic,i))
              call print_incidence_solution(dataset,listincsol(ic,1+i))
              call print_confusion(dataset,i,workstruct(ic)%nalert(i),workstruct(ic)%alertQtl,workstruct(ic)%corrmax(i))
              if ( workstruct(ic)%ntest /= 0 )  then
               call print_test_nuisances(workstruct(ic)%ntest,workstruct(ic)%listtestnuis)
              end if
              if ( i==1 ) then
                call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
                call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))
              end if
              call print_LRT(dataset,ic,lrtsol(ic,i))
              if ( i==2 ) then
!                call create_grid_file_2QTL(dataset,lrtsol(ic,2))
                call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
                call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,2))
              end if
              if ( interaction ) then
                call print_lrt_solution(dataset,lrtsolInteraction(ic))
                call print_incidence_solution(dataset,incsolInteraction(ic))
                call print_confusion(dataset,nqtl+1,workstruct(ic)%nalert(nqtl+1),&
                 workstruct(ic)%alertQtl,workstruct(ic)%corrmax(nqtl+1))
                call lrtsolInteraction(ic)%release()
                call incsolInteraction(ic)%release()
              end if
           end do ! iqtl
        end do

        if ( (contingence_effects%qtl > 0) ) then
         do i=1,nqtl
          call print_summary_analyse(dataset,lrtsol(:,i),listincsol(:,i+1),i)
         end do
        end if

       end if

    ! print txt transcriptome format
       if (opt_sim == TRANSCIPTOME_ANALYSE) then
          call  print_transcriptome_H0(dataset,listincsol(:,1))
          do i=1,nqtl
            call  print_transcriptome(dataset,lrtsol(:,i),listincsol(:,i+1))
          end do
       end if

       call end_analyse_linear

       do ic=1,dataset%phenoModel%ncar
        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
        call release_ws(workstruct(ic))
       end do

       deallocate (workstruct)
       deallocate (YCOMPACT)
       deallocate (AllsigQtlMoinUn)

   ! appel fonction CUDA pour liber de la memoire sur le gpu
       call release_allocated_internal_sig()

    end subroutine opti_unitrait_incidence_cuda
#endif
!!***


!!****f* m_qtlmap_analyse/opti_multitrait_modlin
!!  NAME
!!    opti_multitrait_modlin
!!  DESCRIPTION
!!   LA multi trait with model description (multi qtl version, interaction with fixed effect and qtl are not managed)
!!  INPUTS
!!    nqtl       : number of qtl to test. nqtl >= 1
!!    opt_sim    : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_multitrait_modlin(dataset,spt,lrtsol,listincsol,opt_sim,nqtl,type_resol)
    use m_qtlmap_incidence_multi
    type(QTLMAP_DATASET)       ,intent(inout) :: dataset
    type(PDD_BUILD)            ,intent(inout) :: spt
    integer                  ,intent(in)                                   :: nqtl,type_resol
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(nqtl)              :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,nqtl+1) :: listincsol
    integer                  ,intent(in)                                   :: opt_sim

    integer :: i,j,iq
    type(INCIDENCE_GEN_STRUCT)   :: workstruct

    call init_analyse_gen(dataset)
    call init_workstruct(dataset,-1,workstruct,CONTINGENCE_QTL,type_resol,nqtl,nqtl,&
        (opt_sim == COMMON_ANALYSE),(opt_sim == COMMON_ANALYSE),.true.)

    if ( is_parameter(MODEL_MULTITRAIT_POLYGENIC,type_resol) ) then
        call opti_0qtl_multi( dataset, listincsol(1,1) ,workstruct, nqtl,model_optim_multi_h0,.false.)
    !      call opti_0qtl_multi( listincsol(1,1) ,workstruct, nqtl,.false.)
    else if ( is_parameter(MODEL_MULTITRAIT_POLYGENIC_LU,type_resol)  ) then
        call opti_0qtl_multi(dataset,listincsol(1,1) ,workstruct, nqtl,model_optim_multi_h0_LU,.true.)
    else
        call stop_application("ERROR DEV opti_multitrait_modlin");
    end if

    !workstruct%sigsquare(1,:)=workstruct%sigsquareEstime

    if (opt_sim == COMMON_ANALYSE) then
        call print_start_multitraits()
        call print_incidence_solution(dataset,listincsol(1,1))
    end if

    do i=1,nqtl
        workstruct%nqtl=i
        workstruct%hypothesis=i
        if ( is_parameter(MODEL_MULTITRAIT_POLYGENIC,type_resol) ) then
            call opti_nqtl_multi(dataset,spt,i,i,workstruct,listincsol(1,1+i),lrtsol(i),model_optim_multi_hn,nqtl,.false.)
            ! call opti_nqtl_multi(i,workstruct,listincsol(1,1+i),lrtsol(i),model_optim_multi_hn,nqtl)
        else if ( is_parameter(MODEL_MULTITRAIT_POLYGENIC_LU,type_resol)  ) then
            call opti_nqtl_multi(dataset,spt,i,i,workstruct,listincsol(1,1+i),lrtsol(i),model_optim_multi_hn_LU,nqtl,.true.)
        else
            call stop_application("ERROR DEV opti_multitrait_modlin");
        end if

        if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
            call print_lrt_solution(dataset,lrtsol(i))
            call print_incidence_solution(dataset,listincsol(1,1+i))
            call print_marker_information_at_the_max(dataset,spt,lrtsol(i))

            !call print_confusion(workstruct%nalert,workstruct%alertQtl,workstruct%corrmax)
            !              if ( workstruct(ic)%ntest /= 0 )  then
            !               call print_test_nuisances(workstruct(ic)%ntest,workstruct(ic)%listtestnuis)
            !              end if
            !if ( i==1 ) then
                !call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))

            !end if
             call print_LRT(dataset,0,lrtsol(i))


            if ( i==2 ) then
                !call create_grid_file_2QTL(dataset,lrtsol(2))
                !call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
                !call print_LRT_2QTL(dataset,lrtsol(2))
            end if

        end if
    end do ! iqtl

    call end_analyse_gen


end subroutine opti_multitrait_modlin
!!***


subroutine opti_unitrait_biall_farnir(dataset,spt,lrtsol,listincsol,opt_sim,nbqtl)
    use m_qtlmap_analyse_biallelic_farnir
    type(QTLMAP_DATASET)       ,intent(in)      :: dataset
    type(PDD_BUILD)            ,intent(in)      :: spt
    integer                  ,intent(in)                                   :: opt_sim,nbqtl
    type(TYPE_LRT_SOLUTION)  , intent(inout)    ,dimension(dataset%phenoModel%ncar,nbqtl)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,nbqtl+1) :: listincsol


    type(work_bicfarnir) :: dataw
    integer :: ic

    do ic=1,dataset%phenoModel%ncar
        call init_work(dataset,spt,ic,dataw)
        call opti_0qtl(dataset,dataw)
        !call opti_1qtl(dataset,spt,dataw,lrtsol(ic))
        call opti_1qtl(dataset,spt,dataw,lrtsol(ic,1))
        call set_solution_hypothesis0(dataset,dataw,listincsol(ic,1))
        call set_solution_hypothesis1(dataset,dataw,listincsol(ic,2))
        call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
        call print_incidence_solution(dataset,listincsol(ic,1))
        call print_lrt_solution(dataset,lrtsol(ic,1))
        call print_incidence_solution(dataset,listincsol(ic,2))
    end do

end subroutine opti_unitrait_biall_farnir

!****f* m_qtlmap_analyse/opti_unitrait_interaction
!  NAME
!    opti_unitrait_interaction
!  DESCRIPTION
!    *** DEV interaction qtl analysis ****
! SOURCE
!     subroutine opti_unitrait_interaction(lrtsol,opt_sim,biq)
!       use m_qtlmap_incidence
!       use m_qtlmap_incidence_linear
!       use m_qtlmap_incidence_optim
!       use m_qtlmap_analyse_lin_gen
!       type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,2)  :: lrtsol
!       integer                  ,intent(in)                                :: opt_sim
!
!       integer                                                             :: type_incidence
!       logical                                                             :: biq
!       !type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,3)   :: listincsol
!
!       integer                               :: ic,supnbnivest,typeModel,i,j,ip
!       integer :: ma_model
!
!       real (kind=dp)   ,dimension(0:2,dataset%genea%np)                                 :: sigsquare
!       real (kind=dp)          ,dimension(dataset%genea%np)                                 :: sigsquareEstime
!
!       !1 ->H0, 2->H1, 3->H2, 4->H3 (avec interaction)
!       type(TYPE_INCIDENCE_SOLUTION)    ,dimension(dataset%phenoModel%ncar,4)        :: listincsol
!       type(TYPE_LRT_SOLUTION)                  ,dimension(dataset%phenoModel%ncar)  :: lrtsolInteraction
!
!       type(INCIDENCE_GEN_STRUCT)  ,dimension(:) ,allocatable :: workstruct
!       logical                                                :: islinear
!
!
!       call init_analyse_gen(dataset)
!       allocate (workstruct(dataset%phenoModel%ncar))
!
!       ma_model = MASK_INIT_MOD_INC
!       call set_parameter(KIND_HOMOSCEDASTIC,ma_model)
!       call set_parameter(MODEL_UNITRAIT_POLYGENIC,ma_model)
!       call set_parameter(RESOLUTION_LINEAR,ma_model)
!
!       if ( is_parameter(RESOLUTION_LINEAR,ma_model)) then
!            call init_analyse_linear()
!            islinear = .true.
!       else
!            islinear = .false.
!       end if
!
!       do ic=1,dataset%phenoModel%ncar
!           if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or. dataset%phenoModel%natureY(ic) == 'a')) cycle
!           sigsquare=0.d0
!
!           call init_workstruct(ic,workstruct(ic),CONTINGENCE_QTL,ma_model,2,3,&
!           (opt_sim == COMMON_ANALYSE),(opt_sim == COMMON_ANALYSE),.false.)
!           if (islinear) then
!            call opti_0qtl( ic, listincsol(ic,1),workstruct(ic),model_lin_h0)
!           else
!            call opti_0qtl( ic, listincsol(ic,1),workstruct(ic),model_optim_h0)
!           end if
!
!           if (opt_sim == COMMON_ANALYSE) then
!              call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
!              call print_incidence_solution(dataset,listincsol(ic,1))
!           end if
!
!           do i=1,2
!             workstruct(ic)%nqtl=i
!             workstruct(ic)%hypothesis=i
!             workstruct(ic)%biq=biq
!             if (islinear) then
!               call opti_nqtl(i,i,ic,workstruct(ic),listincsol(ic,1+i),lrtsol(ic,i),model_lin_hn)
!             else
!               call opti_nqtl(i,i,ic,workstruct(ic),listincsol(ic,1+i),lrtsol(ic,i),model_optim_hn)
!             end if
!
!             if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
!              call print_lrt_solution(dataset,lrtsol(ic,i))
!              call print_incidence_solution(dataset,listincsol(ic,1+i))
!              call print_confusion(i,workstruct(ic)%nalert(i),workstruct(ic)%alertQtl,workstruct(ic)%corrmax(i))
!              if ( workstruct(ic)%ntest /= 0 )  then
!               call print_test_nuisances(workstruct(ic)%ntest,workstruct(ic)%listtestnuis)
!              end if
!              if ( i==1 ) then
!                call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
!                call print_LRT(dataset,lrtsol(ic,1))
!              end if
!
!              if ( i==2 ) then
!                call create_grid_file_2QTL(dataset,lrtsol(ic,2))
!                call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
!                call print_LRT_2QTL(dataset,lrtsol(ic,2))
!              end if
!
!             end if
!           end do ! i
!
!           !Test avec interaction seulement lineaire la resolution.....a voir par la suite...
!           workstruct(ic)%effects=CONTINGENCE_QTL_INTER
!           call opti_2qtl_linear_interaction(ic,workstruct(ic),listincsol(ic,4),lrtsolInteraction(ic),typeModel)
!
!           if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
!              call print_lrt_solution(dataset,lrtsol(ic,2))
!              call print_incidence_solution(dataset,listincsol(ic,4))
!              call create_grid_file_2QTL(dataset,lrtsol(ic,2))
!              call print_pat_mat_effect_2QTL(dataset,ic,lrtsol(ic,2))
!              call print_LRT_2QTL(dataset,lrtsol(ic,2))
!           end if
!
!            ! end if
!       end do ! ic
!       call end_analyse_linear
!
!       do ic=1,dataset%phenoModel%ncar
!        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
!        do i=1,3
!         call release(listincsol(ic,i))
!        end do
!       end do
!
!       do ic=1,dataset%phenoModel%ncar
!        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
!        call release_ws(workstruct(ic))
!       end do
!
!       deallocate (workstruct)
!       call end_analyse_gen
!
!
!    end subroutine opti_unitrait_interaction
!***

!!****f* m_qtlmap_analyse/opti_unitrait_modlin
!!  NAME
!!    opti_unitrait_modlin
!!  DESCRIPTION
!!   LA single trait with model description (1 qtl)
!!  INPUTS
!!    opt_sim    : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_unitrait_modlin(dataset,spt,lrtsol,listincsol,opt_sim)
    use m_qtlmap_analyse_modlin
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,1)    :: lrtsol
    integer          ,intent(in)                                          :: opt_sim
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,2)     :: listincsol
    real (kind=dp)                        :: fmax

    integer                               :: ic,supnbnivest
    !modif carole : l'estimationd e la moyenne est parametrable
    logical                               :: est_moy = .true.
    logical                               :: est_var = .true.
    integer                                           :: i

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(fmax,supnbnivest) IF ( opt_sim /= COMMON_ANALYSE .and. (dataset%phenoModel%ncar>1) )
    call init_analyse_gen(dataset)
    call init_analyse_modlin(dataset,spt)

    !$OMP DO
    do ic=1,dataset%phenoModel%ncar
        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call init_analyse_lin_gen(dataset,spt,ic,1)
        call optinit(dataset,ic)
        call opti_0qtl_modlin(dataset,spt,ic)
        call opti_1qtl_modlin(dataset,spt,ic,lrtsol(ic,1),fmax,supnbnivest)


        if ( opt_sim /= SIMULATION) then
            ! Formattage des solutions sous H0 et H1
            call set_solution_hypothesis0(dataset,ic,listincsol(ic,1))
 
            call set_solution_hypothesis1(dataset,ic,listincsol(ic,2))
        end if

        !a cause de test_lin ,on ne peut pas sortir l affichage de la boucle de traitement
        !cette region devient donc critique pour afficher correctement et succesivement par chaque thread les resultats....
        if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
            !$OMP CRITICAL
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            !  colinearite donne des indications sur la colinearite entre chacun
            !  des effets QTL a estimer et tous les autres effets (avant elimination
            !  des effets non estimables)
            !
            call confusion(dataset,spt,'avant',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')
            ! Impression du resultat sous H0
            call print_incidence_solution(dataset,listincsol(ic,1))
            ! Impression de la courbe LRT sous H1
            call print_lrt_solution(dataset,lrtsol(ic,1))
            ! Impression des resultats sous H1
            call print_incidence_solution(dataset,listincsol(ic,2))
            ! Effet paternel et maternel sous H1 dans un fichier dedie
            call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
            ! Affichage de la valeur LRT a chaque position dans un fichier dedie
            call print_LRT(dataset,ic,lrtsol(ic,1))
            !
            ! test_modlin donne la signification des effets du modele
            !
            call test_lin(dataset,spt,lrtsol(ic,1)%chrmax(0),ic,est_moy,supnbnivest,fmax,lrtsol(ic,1)%nxmax(0))
            !
            !  confusion indique si'l y a une confusion entre les effets finalement
            !  retenus dans le modele
            call confusion(dataset,spt,'apres',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')

            call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))
            !$OMP END CRITICAL
        end if
        call end_analyse_lin_gen
    end do
    !$OMP END DO
    call end_analyse_modlin
    call end_analyse_gen
    !$OMP END PARALLEL

    if ( opt_sim /= SIMULATION) then
        call  print_summary_analyse(dataset,lrtsol(:,1),listincsol(:,2),1)
    end if

     ! print txt transcriptome format
    if (opt_sim == TRANSCIPTOME_ANALYSE) then
        call  print_transcriptome_H0(dataset,listincsol(:,1))
        call  print_transcriptome(dataset,lrtsol(:,1),listincsol(:,2))
    end if

end subroutine opti_unitrait_modlin
!!***

!!****f* m_qtlmap_analyse/opti_unitrait_modlin_new
!!  NAME
!!    opti_unitrait_modlin_new
!!  DESCRIPTION
!!   LA / LD / LDLA /LDJH single trait with model description (1 qtl) optimisation method
!!  INPUTS
!!    option_var  : "hete" : heteroscedastic, "homo" : homoscedastic
!!    option_anal : "LDLA", "LD  ", "LA  ","LDJH"
!!    hdam        : LD option, take dam haplotypes in considration
!!    opt_sim     : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_unitrait_modlin_new(dataset,spt,lrtsol,listincsol,option_var,option_anal,hdam,opt_sim)
    use m_qtlmap_analyse_modlin_ldla
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,1)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,2)     :: listincsol
    character(len=4)         , intent(in)                                 :: option_var,option_anal
    integer          ,intent(in)                                          :: opt_sim
    logical          ,intent(in)                                          :: hdam

    real (kind=dp)       , dimension(dataset%phenoModel%ncar,dataset%genea%np)  :: std_all,ap_all
    logical              , dimension(dataset%genea%nm)              :: subphasm
    real(kind=dp)                                     :: fmax
    integer                               :: ic,supnbnivest,nposx,CHR
    !modif carole : l'estimationd e la moyenne est parametrable
    logical                               :: est_moy = .true.
    logical                               :: est_var = .true.
    logical                               :: qtl,hsire
    type(HAPLOTYPE_POSITION_BUILD)        :: shp_opt

    hsire=( (option_anal == 'LD  ') .or. (option_anal == 'LDLA') .or. (option_anal == 'LDJH') )
    qtl  = ( (option_anal == 'LA  ') .or. (option_anal == 'LDLA') .or. (option_anal == 'LDJH') )

    call init_analyse_gen(dataset)
    call init_analyse_modlin_ldla(dataset)

    do ic=1,dataset%phenoModel%ncar
        if (.not.(dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a')) cycle
        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call shp_opt%set(dataset,spt)

        call init_analyse_lin_gen(dataset,spt,ic,1)
        call optinit(dataset,ic)
        call opti_0qtl_modlin_ldla(dataset,ic,.true.,option_var)
        call opti_1qtl_modlin_ldla(dataset,spt,shp_opt,ic,lrtsol(ic,1),fmax,supnbnivest,.true.,option_var,option_anal,hdam)

        if ( opt_sim /= SIMULATION) then
            ! Formattage des solutions sous H0 et H1
            call set_solution_hypothesis0(dataset,ic,(option_var=="hete"),listincsol(ic,1))
            call set_solution_hypothesis1(dataset,shp_opt,ic,(option_var=="hete"),qtl,hsire,hdam,listincsol(ic,2))
        end if

        if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            !  colinearite donne des indications sur la colinearite entre chacun
            !  des effets QTL a estimer et tous les autres effets (avant elimination
            !  des effets non estimables)
            !
            call confusion(dataset,spt,'avant',lrtsol(ic,1)%chrmax(0),ic,est_moy,option_anal)
            ! Impression du resultat sous H0
            call print_incidence_solution(dataset,listincsol(ic,1))
            ! Impression de la courbe LRT sous H1
            call print_lrt_solution(dataset,lrtsol(ic,1))
            ! Impression des resultats sous H1
            call print_incidence_solution(dataset,listincsol(ic,2))
            ! Effet paternel et maternel sous H1 dans un fichier dedie
            call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
            ! Affichage de la valeur LRT a chaque position dans un fichier dedie
            call print_LRT(dataset,ic,lrtsol(ic,1))
            !
            ! test_modlin donne la signification des effets du modele
            !
            ! call test_lin(lrtsol(ic,1)%chrmax(0),ic,est_moy,supnbnivest,fmax,lrtsol(ic,1)%nxmax(0))
            !
            !  confusion indique si'l y a une confusion entre les effets finalement
            !  retenus dans le modele
            call confusion(dataset,spt,'apres',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')
            call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))
        end if
        call shp_opt%free()

        call end_analyse_lin_gen
    end do

    !       if ( opt_sim /= OPT_SIMULATION) then
    !           call get_eff_paternal_and_total(effp,efft)
    !           call print_summary_0versus1(dataset%phenoModel%ncar,carac,sigt,size(pere),pere,efft,effp,std_all,ap_all,xlrmax0_1,dxma0_1)
    !       end if

    call end_analyse_modlin_ldla
    call end_analyse_gen

end subroutine opti_unitrait_modlin_new
!!***

!!****f* m_qtlmap_analyse/opti_unitrait_modlin_cox
!!  NAME
!!    opti_unitrait_modlin_cox
!!  DESCRIPTION
!!   LA analysis with cox model
!!  INPUTS
!!    opt_sim     : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_unitrait_modlin_cox(dataset,spt,lrtsol,listincsol,opt_sim)
    use m_qtlmap_analyse_modlin_cox
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,1)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,2)     :: listincsol
    integer          ,intent(in)                        :: opt_sim

    real (kind=dp)                        :: fmax
    integer                               :: ic,supnbnivest,npar,i
    logical                               :: est_moy = .false.
    logical                               :: est_var = .false.

    call init_analyse_gen(dataset)
    do ic=1,dataset%phenoModel%ncar
        if ( dataset%phenoModel%natureY(ic) /= 'r' ) then
            call stop_application("This analyse does not manage this nature of data! TRAIT["//trim(dataset%phenoModel%carac(ic))&
                //"] nature:"//trim(dataset%phenoModel%natureY(ic)))
        end if

        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call init_analyse_lin_gen(dataset,spt,ic,1)
        !npar is modified because it is equal to 3*np + 2nm in m_qtlmap_analyse_lin_gen.f95
        npar=dataset%genea%np+dataset%genea%nm
        call init_analyse_modlin_cox(dataset)
        call optinit(dataset,ic)
        call calcul_rang(dataset,ic)
        call opti_0qtl_modlin_cox(dataset,spt,ic,est_moy)
        call opti_1qtl_modlin_cox(dataset,spt,ic,lrtsol(ic,1),fmax,supnbnivest)
        !call stop_application("Survival Analysis with a Cox model (option=5) is in development ")

        if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
             !  colinearite donne des indications sur la colinearite entre chacun
             !  des effets QTL a estimer et tous les autres effets (avant elimination
             !  des effets non estimables)
            call confusion(dataset,spt,'avant',lrtsol(ic,1)%chrmax(0),ic,.false.,'LA  ')
             ! Impression de la courbe LRT sous H1
            call print_lrt_solution(dataset,lrtsol(ic,1))

            ! Formattage des solutions sous H0 et H1
            call set_solution_hypothesis0(dataset,ic,listincsol(ic,1))
            call set_solution_hypothesis1(dataset,ic,listincsol(ic,2))
            ! Impression du resultat sous H0
            call print_incidence_solution_risk_factor(dataset,listincsol(ic,1))
            ! Impression des resultats sous H1
            call print_incidence_solution_risk_factor(dataset,listincsol(ic,2))
             ! Effet paternel et maternel sous H1 dans un fichier dedie
            call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
            ! Affichage de la valeur LRT a chaque position dans un fichier dedie
            call print_LRT(dataset,ic,lrtsol(ic,1))
            !
            ! test_modlin donne la signification des effets du modele
            !
            call test_lin_cox(dataset,spt,lrtsol(ic,1)%chrmax(0),ic,est_moy,supnbnivest,fmax,lrtsol(ic,1)%nxmax(0))
            !
            !  confusion indique si'l y a une confusion entre les effets finalement
            !  retenus dans le modele
            call confusion(dataset,spt,'apres',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')
            call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))
        end if
        call end_analyse_lin_gen
        call end_analyse_modlin_cox
    enddo

    !TODO: Olivier il faut mettre le call print_summary_analyse(dataset,efft,effpp,lrtsol(:,1),listincsol(:,2),1) je n'y arrive pas!!!!! (carole)
     if ( opt_sim /= SIMULATION) then
      call print_summary_analyse(dataset,lrtsol(:,1),listincsol(:,2),1)
     end if


    do ic=1,dataset%phenoModel%ncar
        if (dataset%phenoModel%natureY(ic) /= 'r') cycle
        do i=1,2
            call listincsol(ic,i)%release()
        end do
    end do


    call end_analyse_gen
end subroutine opti_unitrait_modlin_cox
!!***


!!****f* m_qtlmap_analyse/opti_unitrait_discret
!!  NAME
!!    opti_unitrait_discret
!!  DESCRIPTION
!!   LA analysis for discrete data with model description (1 qtl)
!!  INPUTS
!!    opt_sim     : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!! SOURCE
subroutine opti_unitrait_discret(dataset,spt,lrtsol,listincsol,opt_sim)
    use m_qtlmap_analyse_discret_unitrait
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(inout) ,dimension(dataset%phenoModel%ncar,1)    :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(dataset%phenoModel%ncar,2)     :: listincsol

    integer          ,intent(in)                      :: opt_sim

    real (kind=dp)                        :: fmax

    logical              , dimension(dataset%genea%nm)              :: subphasm

    integer                               :: ic,supnbnivest,i
    logical                               :: est_moy = .false.
    logical                               :: est_var = .true.

    call init_analyse_gen(dataset)
    do ic=1,dataset%phenoModel%ncar
        if (dataset%phenoModel%natureY(ic) /= 'i') cycle
        call log_mess("TRAIT ["//trim(dataset%phenoModel%carac(ic))//"]",INFO_DEF)
        call init_analyse_lin_gen(dataset,spt,ic,1)
        call init_analyse_discret_unitrait(dataset,spt)
        call optinit(dataset,ic)
        call opti_0qtl_discret_unitrait(dataset,spt,ic)
        call opti_1qtl_discret_unitrait(dataset,spt,ic,lrtsol(ic,1),fmax,supnbnivest)

        if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
            call print_start_unitrait(dataset,dataset%phenoModel%carac(ic))
            !  colinearite donne des indications sur la colinearite entre chacun
            !  des effets QTL a estimer et tous les autres effets (avant elimination
            !  des effets non estimables)
            !
            call confusion(dataset,spt,'avant',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')

            ! Formattage des solutions sous H0 et H1
            call set_solution_hypothesis0(dataset,ic,listincsol(ic,1))
            call set_solution_hypothesis1(dataset,ic,listincsol(ic,2))
            ! Impression du resultat sous H0
            call print_incidence_solution(dataset,listincsol(ic,1))
            ! Impression de la courbe LRT sous H1
            call print_lrt_solution(dataset,lrtsol(ic,1))
            ! Impression des resultats sous H1
            call print_incidence_solution(dataset,listincsol(ic,2))
            ! Effet paternel et maternel sous H1 dans un fichier dedie
            call print_paternal_maternal_effect(dataset,ic,lrtsol(ic,1))
            ! Affichage de la valeur LRT a chaque position dans un fichier dedie
            call print_LRT(dataset,ic,lrtsol(ic,1))
            !
            ! test_modlin donne la signification des effets du modele
            !
            call test_lin(dataset,spt,lrtsol(ic,1)%chrmax(0),ic,est_moy,supnbnivest,fmax,lrtsol(ic,1)%nxmax(0))
            !
            !  confusion indique si'l y a une confusion entre les effets finalement
            !  retenus dans le modele
            call confusion(dataset,spt,'apres',lrtsol(ic,1)%chrmax(0),ic,est_moy,'LA  ')
            call print_marker_information_at_the_max(dataset,spt,lrtsol(ic,1))
        end if

        call end_analyse_lin_gen
        call end_analyse_discret_unitrait
    end do

    call end_analyse_gen

end subroutine opti_unitrait_discret
!!***

!!****f* m_qtlmap_analyse/opti_multitraits
!!  NAME
!!    opti_multitraits
!!  DESCRIPTION
!!    LA analysis for a set of traits with a multivariate analysis
!!  INPUTS
!!    opt_sim     : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    listincsol : solution for each hypothesis finded by the analysis (see TYPE_INCIDENCE_SOLUTION)
!!    rhoi       : residual correlation matrix finded for each hypothesis
!! SOURCE
subroutine opti_multitraits(dataset,spt,lrtsol,listincsol,rhoi,opt_sim)
    use m_qtlmap_analyse_multitrait
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(inout)                        :: lrtsol
    type(TYPE_INCIDENCE_SOLUTION)  ,intent(inout)  ,dimension(2)    :: listincsol
    real (kind=dp) ,dimension(dataset%phenoModel%ncar,dataset%phenoModel%ncar), intent(out) :: rhoi

    integer , intent(in)                                            :: opt_sim

    logical              , dimension(dataset%genea%nm)              :: subphasm
    character(len=LEN_DEF) , dimension(dataset%genea%nm)            :: submere
    logical                                           :: impfem
    real (kind=dp)       , dimension(dataset%genea%nm)              :: xmoym1,am1
    integer                                           :: nmsub,efft
    integer            , dimension(dataset%genea%np)                :: effpp

    integer    :: ic,i

    call init_analyse_gen(dataset)

    if (dataset%phenoModel%ncar <= 1) then
        call stop_application("can not perform a multitrait analysis on "//trim(str(dataset%phenoModel%ncar))//" traits")
    end if

    do ic=1,dataset%phenoModel%ncar
        if( dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a') cycle
        call stop_application("Multitrait analysis is only available for continuous traits")
    end do

    call init_analyse_multitrait(dataset)
    call optinit_mcar(dataset)
    call opti_mcar_0qtl(dataset)
    call opti_mcar_1qtl(dataset,spt,lrtsol)

    if (opt_sim == COMMON_ANALYSE .or. opt_sim == TRANSCIPTOME_ANALYSE) then ! Representation graphique

        call print_start_multitraits()
        call set_solution_hypothesis0(dataset,listincsol(1))
        call print_incidence_solution(dataset,listincsol(1))
        call print_lrt_solution(dataset,lrtsol)
        call set_solution_hypothesis1(dataset,spt,listincsol(2))
        call print_incidence_solution(dataset,listincsol(2))

        call print_paternal_maternal_effect(dataset,1,lrtsol)
        call print_LRT(dataset,0,lrtsol)
        call print_marker_information_at_the_max(dataset,spt,lrtsol)

    end if

    if ( opt_sim /= SIMULATION) then
      ! call get_eff_paternal_and_total(effpp,efft)
       ! TODO : summary pour multi
    end if

    call end_analyse_multitrait
    call end_analyse_gen

end subroutine opti_multitraits
!!***

!!****f* m_qtlmap_analyse/opti_multitrait_DA
!!  NAME
!!    opti_multitrait_DA
!!  DESCRIPTION
!!    LA for a set of traits (without missing data) with a discriminante analysis
!!  INPUTS
!!    opt_sim     : contexte of the analysis ( COMMON_ANALYSE | SIMULATION | TRANSCIPTOME_ANALYSE )
!!  OUTPUTS
!!    lrtsol     : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    rhoi       : residual correlation matrix finded for each hypothesis
!! SOURCE
subroutine opti_multitrait_DA(dataset,spt,lrtsol,opt_sim)
    use m_qtlmap_analyse_multitrait_DA

    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    type(TYPE_LRT_SOLUTION)  , intent(out)          :: lrtsol
    integer          ,intent(in)                    :: opt_sim

    real (kind=dp) , dimension(:,:) ,allocatable    :: coeff
    real (kind=dp) , dimension(dataset%phenoModel%ncar)     :: coeff_max
    real (kind=dp)       , dimension(dataset%genea%np)  :: std_all,ap_all
    real (kind=dp)       , dimension(dataset%genea%nm)              :: xmoym1,am1
    logical              , dimension(dataset%genea%nm)              :: subphasmm
    character(len=LEN_DEF) , dimension(dataset%genea%nm)              :: submere
    logical                                           :: impfem

    integer            , dimension(dataset%genea%np)                :: effpp
    integer                                           :: i,ic,nmsub,efft,n,ilong,ix,kd,ip,chr
    real (kind=dp), dimension (dataset%genea%nd)                    :: yda
    integer                ,dimension(dataset%phenoModel%ncar)    :: lic
    real (kind=dp)                                    :: dxmax
    type(TYPE_LRT_SOLUTION) ,dimension(1) :: listLrtSol
    type(TYPE_INCIDENCE_SOLUTION)   :: incsol

    call init_analyse_gen(dataset)
    allocate (coeff(dataset%phenoModel%ncar,dataset%map%get_maxnpo()))

    if (dataset%phenoModel%ncar <= 1) then
        call stop_application("can not perform a multitrait analysis on "//trim(str(dataset%phenoModel%ncar))//" traits")
    end if

    do ic=1,dataset%phenoModel%ncar
        if( dataset%phenoModel%natureY(ic) == 'r'.or.dataset%phenoModel%natureY(ic) == 'a') cycle
        call stop_application("Discriminant multitrait analysis is only available for continuous traits")
    end do

    call init_analyse_DA_1QTL(dataset,spt)

    coeff=0.d0

    call lrtsol%new(dataset,1)
    !       ncar=1
    lrtsol%lrtmax(0)=-1.d75
    dxmax=-1.d6
    lrtsol%nxmax=0
    lrtsol%chrmax=1
    do chr=1,dataset%map%nchr
        n=0
        ilong=dataset%map%get_ilong(chr)
        do ix=0,ilong,dataset%map%pas
            n=n+1
            call perform_DA_1QTL(dataset,spt,chr,n,yda,coeff,dataset%map%get_maxnpo())
            call optinit_da(dataset,yda)
            call opti_DA_0qtl(dataset)
            call opti_DA_1qtl(dataset,spt,chr,ix,n,yda,lrtsol)
            if ( lrtsol%chrmax(0) > 0 .and. lrtsol%nxmax(0) > 0 ) then
                if(dataset%map%absi(lrtsol%chrmax(0),lrtsol%nxmax(0)) >  dxmax) then
                    dxmax=dataset%map%absi(lrtsol%chrmax(0),lrtsol%nxmax(0))
                    coeff_max=coeff(:,n)
                end if
                call log_mess('at position ix:'//str(ix)//' LRT max:'//str(lrtsol%lrtmax(0))//&
                    " POS:"//str(dataset%map%absi(lrtsol%chrmax(0),lrtsol%nxmax(0))),DEBUG_DEF)
            else
                dxmax=0.d0
                coeff_max=0.d0
            end if

        end do
    end do

    if (opt_sim == COMMON_ANALYSE) then ! Representation graphique
        std_all(:)=std
        ap_all(:)=ap
        call print_start_multitrait_DA
        call print_lrt_solution(dataset,lrtsol)
        call set_solution_hypothesis1(dataset,incsol)
        call print_incidence_solution(dataset,incsol)
        call incsol%release()

        call print_coeff_linear_combination_max(dataset,coeff_max)
        call print_paternal_maternal_effect(dataset,1,lrtsol)
        call print_LRT(dataset,0,lrtsol)
        call print_coeff_linear_combination(dataset,lrtsol%chrmax(0),&
            dataset%map%get_ilong(lrtsol%chrmax(0)),dataset%map%get_maxnpo(),coeff)
        call print_marker_information_at_the_max(dataset,spt,lrtsol)
    end if

    if ( opt_sim /= SIMULATION) then
    ! call get_eff_paternal_and_total(effpp,efft)
     ! TODO : summary pour multi
    end if

    !    ncar=dataset%phenoModel%ncar
    deallocate (coeff)

    call end_analyse_DA_1QTL
    call end_analyse_gen

end subroutine opti_multitrait_DA
!!***

!!****f* m_qtlmap_analyse/quantile
!!  NAME
!!    quantile
!!  DESCRIPTION
!!   computes quantiles (Harrell & Davis 1982) and write a result file.
!!  INPUTS
!!    simula_file : name of file to print the output
!!    nsim        : number of simulation
!!    opt_calcul  : analysis id
!!    opt_qtl     : number of qtl hypothesis (analysis can be not compatible with this option)
!!    nqtl        : number of qtl simulated
!!    lrtsol      : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!  NOTES
!!    use : MATH_QTLMAP_M01CAF, MATH_QTLMAP_G01AAF, QUANTILE1, print_resume_simulation, print_resume_simulation_2
!! SOURCE
subroutine quantile(dataset,simula_file,nsim,opt_calcul,opt_qtl,nqtl,lrtsol)
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    character(len=LENGTH_MAX_FILE) , intent(in)   :: simula_file
    integer                        , intent(in)   :: nsim,opt_calcul,opt_qtl,nqtl ! nqtl => le nombre de qtl definit dans la simulation
                                                                                  ! opt_qtl => nombre de qtl a detecter pour le calcul
    type(TYPE_LRT_SOLUTION)  , dimension(dataset%phenoModel%ncar,opt_qtl,nsim),intent(in)      :: lrtsol

    real (kind=dp) :: prob(16)
    real (kind=dp) :: x(nsim),wt(nsim)
    real (kind=dp)  ,dimension (16,dataset%phenoModel%ncar)     :: z
    real (kind=dp)  ,dimension (16,1)     :: z2

    CHARACTER(len=LEN_BUFFER_WORD) :: temp,nameT,carac_t(1)
    integer :: nc,ifail,i,im,j,ii,iqtl,iwt
    real (kind=dp) :: xmin1,xmax1,ymu1,sig1,s21,s31,pg,wtsum,p_t,d

    data prob /.1d0,.2d0,.3d0,.4d0,.5d0,.6d0,.7d0,.8d0,.9d0,.95d0,&
        .99d0,.995d0,.9973d0,.999d0,.99947d0,.999947d0/

    if ( nsim < 2 ) return

    nc=dataset%phenoModel%ncar
    nameT=""
    if ( is_multitrait_analysis(opt_calcul) ) then
        nc = 1
        nameT="* All traits *"
    end if

    ifail=0
    iwt = 0
    ! pour tous les lrts trouve sous avec qtl==opt_qtl
    do iqtl=1,opt_qtl   ! ? / opt_qtl
        do im=0,iqtl-1 ! on ecrit tous les test  H: im / opt_qtl
            temp=trim(str(im))//'vs'//trim(str(iqtl))//'Q'
            do i=1,nc
                do j=1,nsim
                    x(j)=lrtsol(i,iqtl,j)%lrtmax(im)
                end do
                if (nameT=="") nameT=trim(dataset%phenoModel%carac(i))

                !Sort X from index 1 to ns ('a' ascendant)....
                call MATH_QTLMAP_M01CAF(X,1,nsim,'A',IFAIL)
                !Computing  Mean, variance, skewness, kurtosis Minimum Maximum
                call MATH_QTLMAP_G01AAF(nsim,x,IWT,WT,YMU1,SIG1,S21,S31,XMIN1,XMAX1,WTSUM,IFAIL)
                do ii=9,16
                    z(ii,i)=QUANTILE1(nsim,prob(ii),X)
                end do

                call print_resume_simulation(nameT,temp,nsim,&
                    YMU1,SIG1,S21,S31,XMIN1,XMAX1,prob,z(:,i),9,16)
                nameT=""
            end do
            if ( .not. is_multitrait_analysis(opt_calcul)  ) then
                call print_resume_simulation_2(temp,dataset%phenoModel%ncar,dataset%phenoModel%carac,z)
            else
                carac_t(1)=trim(nameT)
                z2(:,1)=z(:,1)
                call print_resume_simulation_2(temp,1,carac_t,z2)
            end if
        end do
    end do

end subroutine quantile
!!***

!!****f* m_qtlmap_analyse/QUANTILE1
!!  NAME
!!    QUANTILE1
!!  DESCRIPTION
!!    CALCUL DES QUANTILES D'UNE DISTRIBUTION ORDONNEE
!!  INPUTS
!!    n :
!!    p :
!!    x :
!!  RETURN
!!
!!  NOTES
!!    use : W
!! SOURCE
FUNCTION QUANTILE1(n,p,X) result (Q)
    integer , intent(in)                     :: n
    real (kind=dp) ,intent(in)               :: p
    real (kind=dp) , dimension(:),intent(in) :: X

    real (kind=dp)              :: Q
    integer                     :: i

    Q=0.D0
    IF (N>SIZE(X)) THEN
        call stop_application("QUANTILE1 [m_qtlmap_analyse] N ["//trim(str(N))//&
            '] to high [max='//trim(str(SIZE(X)))//']')
    END IF
    6010 DO 6011 I=1,N
        Q=Q+W(n,I,p)*X(I)
6011 Continue
     RETURN
 END FUNCTION QUANTILE1
 !!***

 !!****f* m_qtlmap_analyse/W
 !!  NAME
 !!     W
 !!  DESCRIPTION
 !!   CALCUL DE W
 !!  INPUTS
 !!    n   :
 !!    i   :
 !!    p   :
 !!  NOTES
 !!    use : MATH_QTLMAP_LOWERTAIL_BETA, MATH_QTLMAP_LOWERTAIL_BETA
 !! SOURCE
 FUNCTION W(n,i,p) result(Wr)
     integer , intent(in)        :: n
     integer , intent(in)        :: i
     real (kind=dp) , intent(in) :: p
     real (kind=dp)              :: Wr,TOL,X1,X2,a,b
     real (kind=dp)              :: p1,p2
     integer                     :: ifail

     TOL=1.D-5
     X1=dble(i)/dble(n)
     X2=dble(i-1)/dble(n)
     a=p*dble(n+1)
     b=(1.D0-p)*dble(n+1)
     ifail=1
     !CALL MATH_QTLMAP_G01EEF(X1,a,b,tol,p1,q1,pdf1,ifail)
     call MATH_QTLMAP_LOWERTAIL_BETA(X1, a, b, TOL, p1, IFAIL)
     ifail=1
     !CALL MATH_QTLMAP_G01EEF(X2,a,b,tol,p2,q2,pdf2,ifail)
     call MATH_QTLMAP_LOWERTAIL_BETA(X2,a,b,tol,p2,ifail)
     Wr=p1-p2
     RETURN
 END FUNCTION W
 !!***


 end module m_qtlmap_analyse
