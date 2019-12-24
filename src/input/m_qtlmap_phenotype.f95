module m_qtlmap_phenotype
    use m_qtlmap_base
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_output_handler
    use m_qtlmap_genealogy

    implicit none
    save

    !!   the unit fortran assigned to the phenotypic file
    integer                                            ,parameter,private  :: unit_perf = 777
    !!   the unit fortran assigned to the model file
    integer                                            ,parameter,private  :: unit_mod  = 3

    !!   nbre minimum de descendant a considerer dans une famille (pere, pere-mere) pour faire des permutations
    integer                                            ,parameter,public   :: NB_DES_MIN=10


    ! INTERNAL VAR
    ! ----------------------------------------------------------------------
    logical                                            ,private             :: GLOBAL_all_mode = .false.
    logical                                            ,private             :: GLOBAL_filter_mode = .false.

    !!   indicate if fixed effect / covariate are take in care (1 otherwise 0) for a specific trait ic
    !!  DIMENSIONS
    !!   nc,nfix+ncov
    integer           , dimension (:,:),allocatable                ,private :: nuis

    !!   buffer array for nuis
    !!  DIMENSIONS
    !!   nc,nfix+ncov
    integer           , dimension (:,:),allocatable                ,private :: nuis_t

    !!   indicate if a qtl is in interaction with a fixed effect in the model trait (1 otherwise 0)
    !!  DIMENSIONS
    !!   nc,nqtl,nfix
    integer         , dimension (:,:,:),allocatable                ,private :: int_qtl

    !!   buffer array for int_qtl
    !!  DIMENSIONS
    !!   nc,nqtl,nfix
    integer         , dimension (:,:,:),allocatable                ,private :: int_qtl_t

    !!   number of qtl - fixed effect interaction defined for each trait
    !!   nombre de qtl en interaction definit pour le caractere
    !!  DIMENSIONS
    !!   ncar
    integer          , dimension (:),allocatable                   ,private :: nb_qtl_def

    integer          , dimension (:),allocatable                   ,private :: nb_qtl_def_t

    !!  DIMENSIONS
    !!   ncar
    !!   nombre de qtl en interaction definit pour le caractere
    !!***
    character(len=LEN_DEF)          , dimension (:),allocatable     ,public :: carac_t
    !!   buffer array for natureY
    !!  DIMENSIONS
    !!   ncar
    character(len=1)               ,  dimension (:),allocatable     ,public :: natureY_t
    !!****d* m_qtlmap_traits/MAX_QTL
    !!  NAME
    !!   MAX_QTL
    !!  DESCRIPTION
    !!   Maximum number of interaction
    !!  NOTES
    !!   nombre maximum possible de qtl en interaction definit dans le fichier model
    !!***
    integer                    , parameter, public                          :: MAX_QTL=10
    !!   buffer array for h2
    !!  DIMENSIONS
    !!   ncar
    real (kind=dp)                   , dimension (:), allocatable   ,public :: h2_t

    !!   buffer array for RhoP
    !!  DIMENSIONS
    !!   ncar,ncar
    real (kind=dp)                 , dimension (:,:), allocatable   ,public :: RhoP_t

    !!   buffer array for RhoG
    !!  DIMENSIONS
    !!   ncar,ncar
    real (kind=dp)                 , dimension (:,:), allocatable   ,public :: RhoG_t

    type, private :: WORK_PHENOTYPE

        integer                                                       ,public  :: nbete = 0

        character(len=LEN_DEF) , dimension (:),pointer                ,public  :: bete     => null()

        !!   level of the fixed effect fe for a progeny kd
        !!  DIMENSIONS
        !!   nd,nfix
        character(len=LEN_DEF)         , dimension (:,:),pointer   ,public :: niv => null()

        !!   value of the covariate cov for a progeny kd
        !!  DIMENSIONS
        !!   nd,ncov
        real (kind=dp)          , dimension (:,:), pointer         ,public :: cov => null()

        !!   value of censured data for a trait ic/progeny kd
        !!  DIMENSIONS
        !!   nc,nbete

        real (kind=dp)          , dimension (:,:), pointer         ,public :: cdt => null()

        !!   phenotypic value for a progeny kd
        !!  DIMENSIONS
        !!   nc,nbete
        real (kind=dp)          , dimension (:,:), pointer         ,public :: valeur => null()

        !!   validity of a phenotypic value for a progeny kd
        !!  DIMENSIONS
        !!   nc,nbete
        integer           , dimension (:,:),pointer                ,public :: ndelt => null()

    contains
        procedure, public :: release => release_work_phenotype

    end type WORK_PHENOTYPE

    !END DATA PERF
    !-------------

    public :: read_traits
    public :: read_model
    public :: fixe_structure_model
    public :: write_perf
    public :: set_estime
    public :: normalize_data
    public :: set_count_discrete
    public :: set_proportion_discrete
    public :: manage_data

contains

    !!  DESCRIPTION
    !!    Read the user phenotypic files and the model file
    !!  INPUTS
    !!   calculCd            : infer the censored data
    subroutine read_traits(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        logical  :: calculCd

        integer :: i,j

        character(len=LEN_DEF) :: temp
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        type(WORK_PHENOTYPE) :: work

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal
        calculCd = dataset%cli%cli_is_calcul_cd()

        call log_mess('SUBROUTINE read_traits',DEBUG_DEF)

        if (dataset%cli%cli_is_transcriptomic_data()) then ! transcriptome information ...
            call read_perf_by_line(dataset,work)
        else
            if ( GLOBAL_all_mode ) then
              call stop_application("You can not used keyword 'all' in the model file without the option : --data-transcriptomic")
            end if
            call read_perf_by_column(dataset,work)
        end if

        !************* NCAR ACCORDING TO FILTER *********
        if ( dpm%n_filter_car/=0 ) dpm%ncar = dpm%n_filter_car

        if ( calculCd ) then
            do i=1,work%nbete
                !print *,work%bete(i),work%valeur(:,i),work%cdt(:,i)
                call dataset%datasetUser%init_perf_animal(dpm,work%bete(i),work%valeur(:,i),work%cdt(:,i))
            end do
           !stop
        end if

        !peut etre mettre cette methode dans la routine de lecture du model....
        call init_model_struct(dataset)
        !!==
        call initialise_struct_internal(dataset,work)
        call check_traits_and_fathers(dataset)
        call set_estime(dataset)

        call work%release()

        deallocate (nuis)
        deallocate (int_qtl)

        call log_mess('END SUBROUTINE read_traits',DEBUG_DEF)
    end subroutine read_traits

    !!    Read the model file
    !!  INPUTS
    !!   file            : model file name
    !!   filter_car      : vector of index to select a subset of trait defined in the model file
    subroutine read_model(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)       :: dataset
        integer                     :: ios,alloc_stat
        character(len=LEN_LINE)            ::line_read,line_orig,saveLineRead
        character(len=LEN_DEF)        :: word,word2
        integer                     :: i,j,l,k,ii,iqtl,ic
        logical                     :: is_ok,interQtl
        character(len=LEN_DEF),dimension(:),allocatable :: lcar,carac_t2,natureY_t2
        type(DATAMODEL_BASE) , pointer :: dpm
        character(len=LENGTH_MAX_FILE) :: in_param_ef

        call log_mess('reading model file...',INFO_DEF)
        allocate (dataset%phenoModel)
        dpm => dataset%phenoModel

        in_param_ef = dataset%params%get_file_val(K_MODEL)

        open(unit_mod,FILE=in_param_ef,IOSTAT=ios,action="read",status="old")

        if ( ios /= 0 ) then
            call stop_application('Cannot open model file ['//trim(in_param_ef)//']')
        endif

        call log_mess('Model description....',VERBOSE_DEF)

        read(unit_mod,*,IOSTAT=ios) dpm%ncar
        if (ios /= 0) then
            call stop_application('Cannot read number of trait ['//trim(in_param_ef)//']' &
                //' line:1')
        end if

        call log_mess('Number of traits            :'//trim(str(dpm%ncar)),VERBOSE_DEF)

        read(unit_mod,*,IOSTAT=ios) dpm%nfix,dpm%ncov
        if (ios /= 0) then
            call stop_application('Cannot read number of fixed effect and number of covariate ['//&
                trim(in_param_ef)//']' &
                //' line:2')
        end if
        call log_mess('Number of fixed effects     :'//trim(str(dpm%nfix)),VERBOSE_DEF)
        call log_mess('Number of covariates        :'//trim(str(dpm%ncov)),VERBOSE_DEF)

        !Lecture des cov et effets fixe
        CALL GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
        if ( ios /= 0) then
            call stop_application('give name of fixed effect and covariate in the file ['//&
                trim(in_param_ef)//']' &
                //' line:3'//' by default write "none"')
        endif

        allocate (dpm%namefix(dpm%nfix))
        allocate (dpm%namecov(dpm%ncov))

        l = 3
        do i=1,dpm%nfix
            dpm%namefix(i) = trim(next_word(line_read,is_ok))
            if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,'fixed effect name ['//trim(str(i))//']')
            call log_mess('Fixed effect   :'//trim(trim(dpm%namefix(i))),VERBOSE_DEF)
        end do
        do i=1,dpm%ncov
            dpm%namecov(i) = next_word(line_read,is_ok)
            if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,'covariate name ['//trim(str(i))//']')
            call log_mess('Covariate      :'//trim(trim(dpm%namecov(i))),VERBOSE_DEF)
        end do

        !lecture des caracteres
        !----------------------
        ALLOCATE (carac_t(dpm%ncar), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'carac')
        ALLOCATE (nuis_t(dpm%ncar,dpm%nfix+dpm%ncov), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'nuis')
        ALLOCATE (int_qtl_t(dpm%ncar,MAX_QTL,dpm%nfix), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'int_qtl')
        int_qtl_t=0

        ALLOCATE (nb_qtl_def_t(dpm%ncar))
        nb_qtl_def_t=0
        allocate(natureY_t(dpm%ncar), stat = alloc_stat)
        call check_allocate(alloc_stat,'natureY')

        GLOBAL_all_mode = .false.
        dpm%ncarcat = 0
        do i=1,dpm%ncar
            !for each trait we read the model description
            line_read=""
            do while(trim(line_read)=="" .and. ios==0)
             l = l+1
             call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
            end do

            if ( ios /= 0) then
                call stop_application('none description trait find for the'// trim(str(i))//' trait ' &
                    //' line:'//trim(str(l))//'model file:['//trim(in_param_ef)//']')
            end if

            word = next_word(line_read,is_ok)
            if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,'*missing name of trait*')

            !! ALL is a generic name for apply model on all traits
            if (word == ALL_LABEL_MODEL) then
                call log_mess('Keywords ['//ALL_LABEL_MODEL//'] detected. apply the description for each trait',INFO_DEF)
                GLOBAL_all_mode = .true.
            end if

            if (.not. GLOBAL_all_mode) then
                carac_t(i) = trim(word)
            end if

            word = trim(next_word(line_read,is_ok))

            if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,'*missing nature of trait*')

            if ( word /= 'r' .and. word /= 'i' .and. word /= 'c' .and.word /= 'a') then
                call stop_application('Trait nature have to be r(real), i(integer), c(categorial) or a(average)'//&
                    ' token:'//trim(word))
            end if

            if ( .not. GLOBAL_all_mode ) then
                natureY_t(i)=trim(word)
            else
                do j=i,dpm%ncar
                    natureY_t(j)=trim(word)
                end do
            end if

            if ( .not. GLOBAL_all_mode ) then
                if ( natureY_t(i) == 'c') dpm%ncarcat = dpm%ncarcat +1
            else

                do j=i,dpm%ncar
                    if ( natureY_t(j) == 'c') dpm%ncarcat = dpm%ncarcat +1
                end do
            end if

            do j=1,(dpm%nfix+dpm%ncov)
                word = trim(next_word(line_read,is_ok))
                if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,'model for trait ['//&
                    trim(carac_t(i))//']')

                nuis_t (i,j) = get_int(word,is_ok)

                if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,&
                    'possible value for model : 0,1 [val:'//trim(word)//']')
                !ALL mode
                if (GLOBAL_all_mode) then
                    do k=i+1,dpm%ncar
                        nuis_t (k,j) = nuis_t(i,j)
                    end do
                end if
            end do

            iqtl=0
            interQtl=(dpm%nfix/=0)

            do while (interQtl)
                saveLineRead=line_read
                word = trim(next_word(saveLineRead,is_ok))

                if ( is_ok ) then
                    iqtl=iqtl+1
                    if ( MAX_QTL < iqtl ) call stop_on_error(1,in_param_ef,l,&
                        "You can not define more than ["//trim(str(MAX_QTL))//"] qtl interaction.")
                    do j=1,dpm%nfix
                        word = trim(next_word(line_read,is_ok))
                        if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,&
                            'possible value for model : 0,1 [val:'//trim(word)//'] numQTL '//&
                            '['// trim(str(iqtl))//']')

                        int_qtl_t(i,iqtl,j) = get_int(word,is_ok)
                        if ( .not. is_ok ) call stop_on_error (1,in_param_ef,l,&
                            'possible value for model : 0,1 [val:'//trim(word)//'] numQTL '//&
                            '['// trim(str(iqtl))//']')
                        !ALL mode
                        if (GLOBAL_all_mode) then
                            do k=i+1,dpm%ncar
                                int_qtl_t (k,iqtl,j) = int_qtl_t(i,iqtl,j)
                            end do
                        end if
                    end do
                else
                    interQtl = .false.
                end if
            end do

            if (GLOBAL_all_mode) then
                nb_qtl_def_t = iqtl
            else
                nb_qtl_def_t(i) = iqtl
            end if

            call log_mess('Description of ['//trim(carac_t(i))//']... ok',VERBOSE_DEF)

            if (GLOBAL_all_mode) then
                exit ! go out
            end if
        end do

        !Lecture des heritabilite et des matrices de correlations genetique et phenotypique
        ! on lit cette matrice si la ligne "matrice_correlation" existe
        allocate(h2_t(dpm%ncar))
        allocate(RhoP_t(dpm%ncar,dpm%ncar))
        allocate(RhoG_t(dpm%ncar,dpm%ncar))

        h2_t=0.5d0
        RhoP_t=0.5d0
        RhoG_t=0.5d0

        ios=0
        line_read=''
        do while( ios == 0 .and. (trim(line_read)=='') )
            call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
        end do
        saveLineRead=line_read
        word = trim(next_word(line_read,is_ok))

        if ( trim(word) == "CORRELATION_MATRIX" ) then
            ios=0
            ic=0
            do while( ios == 0 .and. ic < dpm%ncar )
                call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
                if (ios /= 0) cycle
                if (trim(line_read) =='' ) cycle
                saveLineRead=line_read
                ic=ic+1
                do i=1,dpm%ncar
                    !Heritability
                    if ( ic == i ) then
                        word = trim(next_word(line_read,is_ok))
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of heritability :"//trim(saveLineRead))
                        end if
                        h2_t(ic) = get_real(word,is_ok)
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of heritability :"//trim(saveLineRead))
                        end if

                         !checking heritability
                        if (h2_t(ic)>1) then
                            call stop_application("Heritability for the "//trim(str(ic))//" is greater than 1")
                        end if
                        if (h2_t(ic)<0) then
                            call stop_application("Heritability for the "//trim(str(ic))//" is less than 0")
                        end if
                    end if
                    !Phenotype correlation
                    if ( i < ic ) then
                        word = trim(next_word(line_read,is_ok))
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of Phenotypique correlation :"//trim(saveLineRead))
                        end if
                        RhoP_t(ic,i) = get_real(word,is_ok)
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of Phenotypique correlation :"//trim(saveLineRead))
                        end if
                        RhoP_t(i,ic) = RhoP_t(ic,i)
                    end if

                    !Genetic correlation
                    if ( i > ic ) then
                        word = trim(next_word(line_read,is_ok))
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of Phenotypique correlation :"//trim(saveLineRead))
                        end if
                        RhoG_t(ic,i) = get_real(word,is_ok)
                        if ( .not. is_ok ) then
                            call stop_application("Bad definition of Phenotypique correlation :"//trim(saveLineRead))
                        end if
                        RhoG_t(i,ic) = RhoG_t(ic,i)
                    end if
                end do
            end do

            if ( ic /= dpm%ncar ) then
                print *,' Format : '
                print *,'CORRELATION_MATRIX'
                print *,'h2[C1]             RhoG[C1,C2]          RhoG[C1,C3] ....'
                print *,'RhoP[C1,C2]        h2[C2]               RhoG[C2,C3] ....'
                print *,'RhoP[C1,C3]        RhoP[C2,C3]          h2[C3] ....'
                print *,' ** '
                call stop_application(" ** Bad definition of CORRELATION_MATRIX **")
            end if

            saveLineRead=''
        else
            call log_mess("CORRELATION_MATRIX keyword are not found. Default value is 0.5 for "//&
                "heritability and phenotypic,genotypic correlations",WARNING_DEF);

        end if

        GLOBAL_filter_mode = .false.

        !AJOUT OFI
        ! SI L UTILISATEUR MET UNE LISTE, ON CONSIDERE CETTE LISTE COMME ETANT LA LISTE DES CARACTERES A PRENDRE EN COMPTE DANS L ANALYSE
        ! PAR DEFAUT, donc,ON A ANALYSE TOUS LES CARACTERES DEFINIS

        dpm%n_filter_car=0
        line_read=saveLineRead
        ios=0
        do while( ios == 0 .and. (trim(line_read)=='') )
            call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
        end do

        if ( ios /= 0 ) return

        if ( line_read /= FILTER_TRAITS_MODEL ) then
           call stop_application("[model] expected Keyword "//trim(FILTER_TRAITS_MODEL))
        end if

        GLOBAL_filter_mode = .true.
        allocate (lcar(dpm%ncar))
        lcar=''
        line_read=''
        ! Selection of particular trait
        ios=0
        do while ( trim(line_read)=='' .and. ios == 0 )
            call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
        end do

        do while( ios == 0)
            ! liste des caracteres a prendre en consideration lors de l analyse
            line_orig=trim(line_read)
            do while ( trim(line_read)/='' )
                dpm%n_filter_car = dpm%n_filter_car + 1
                if (dpm%n_filter_car > dpm%ncar ) then
                    call log_mess("Filter trait:"//trim(line_orig),WARNING_DEF)
                    ios=-1
                    exit
                end if
                lcar(dpm%n_filter_car)=trim(next_word(line_read,is_ok))
            end do
            line_read=''
            do while ( trim(line_read)=='' .and. ios == 0 )
                call GET(unit_mod,line_read,maxlen=LEN_LINE,IOSTAT=ios)
            end do
        end do

        if (associated(dataset%phenoModel%filter_car)) deallocate(dataset%phenoModel%filter_car)
        allocate (dataset%phenoModel%filter_car(dpm%n_filter_car))
        allocate (dataset%phenoModel%filter_car_name(dpm%n_filter_car))

        if (dpm%n_filter_car == 0 ) then
          call stop_application("[model] missing a list traits to analyse.")
        end if

        if (dpm%n_filter_car > dpm%ncar) then
            call stop_application ("Too many traits defines in the list trait filter [Model file:"//&
                trim(in_param_ef)//"]")
        end if

        close(unit_mod)

        do ii=1,dpm%n_filter_car
            do i=1,dpm%ncar
                if ( trim(carac_t(i)) == trim(lcar(ii))) exit
            end do

            dataset%phenoModel%filter_car_name(ii)=trim(lcar(ii))

            if ( i <= dpm%ncar ) then
                dataset%phenoModel%filter_car(ii)=i
            else
                !eqtl...
                if ( GLOBAL_all_mode ) then
                 dataset%phenoModel%filter_car(ii)=-1
                else
                 call stop_application("Missing definition of trait ["//trim(lcar(ii))//"]")
                end if
            end if
        end do

        deallocate (lcar)

        call log_mess('END SUBROUTINE read_model',DEBUG_DEF)
    end subroutine  read_model

    !!    Initialize the general array : carac,nuis,int_qtl,natureY,h2,RhoP,RhoG and release buffer array associated
    !!  INPUTS
    !!   filter      : vector of index to select a subset of trait defined in the model file
    subroutine fixe_structure_model(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        integer :: i,j
        type(DATAMODEL_BASE) , pointer :: dpm
        call log_mess("fixe_structure_model",DEBUG_DEF)
        dpm => dataset%phenoModel

        if (GLOBAL_all_mode) then
           if (dpm%n_filter_car /= 0) then
             allocate (dpm%carac(dpm%n_filter_car))
             do i=1,dpm%n_filter_car
              dpm%carac(i)=dataset%phenoModel%filter_car_name(i)
             end do

             allocate (nuis(dpm%n_filter_car,dpm%nfix+dpm%ncov))
             allocate (int_qtl(dpm%n_filter_car,MAX_QTL,dpm%nfix))
             allocate (nb_qtl_def(dpm%n_filter_car))
             allocate (dpm%natureY(dpm%n_filter_car))
             allocate (dpm%h2(dpm%n_filter_car))
             allocate (dpm%RhoP(0,0))
             allocate (dpm%RhoG(0,0))

             do i=1,dpm%n_filter_car
              nuis(i,:)=nuis_t(1,:)
              int_qtl(i,:,:)=int_qtl_t(1,:,:)
              nb_qtl_def(i)=nb_qtl_def_t(1)
              dpm%natureY(i)=natureY_t(1)
              dpm%h2(i) = h2_t(1)
             end do
           end if
        else ! GLOBAL_MODE
         if (  dpm%n_filter_car /= 0 ) then
            allocate (dpm%carac(dpm%n_filter_car))
            allocate (nuis(dpm%n_filter_car,dpm%nfix+dpm%ncov))
            allocate (int_qtl(dpm%n_filter_car,MAX_QTL,dpm%nfix))
            allocate (nb_qtl_def(dpm%n_filter_car))
            allocate (dpm%natureY(dpm%n_filter_car))

            allocate (dpm%h2(dpm%n_filter_car))
            allocate (dpm%RhoP(dpm%n_filter_car,dpm%n_filter_car))
            allocate (dpm%RhoG(dpm%n_filter_car,dpm%n_filter_car))

            do i=1,dpm%n_filter_car
                dpm%carac(i)=carac_t(dpm%filter_car(i))
                nuis(i,:)=nuis_t(dpm%filter_car(i),:)
                int_qtl(i,:,:)=int_qtl_t(dpm%filter_car(i),:,:)
                nb_qtl_def(i)=nb_qtl_def_t(dpm%filter_car(i))
                dpm%natureY(i)=natureY_t(dpm%filter_car(i))
                dpm%h2(i) = h2_t(dpm%filter_car(i))

                do j=i+1,dpm%n_filter_car
                    dpm%RhoP(i,j) = RhoP_t(dpm%filter_car(i),dpm%filter_car(j))
                    dpm%RhoP(j,i) = dpm%RhoP(i,j)
                    dpm%RhoG(i,j) = RhoG_t(dpm%filter_car(i),dpm%filter_car(j))
                    dpm%RhoG(j,i) = dpm%RhoG(i,j)
                end do
            end do
         end if
         end if ! GLOBAL MODE

        if (  dpm%n_filter_car == 0 ) then
            allocate (dpm%carac(dpm%ncar))
            allocate (nuis(dpm%ncar,dpm%nfix+dpm%ncov))
            allocate (int_qtl(dpm%ncar,MAX_QTL,dpm%nfix))
            allocate (nb_qtl_def(dpm%ncar))
            allocate (dpm%natureY(dpm%ncar))
            allocate (dpm%h2(dpm%ncar))
            allocate (dpm%RhoP(dpm%ncar,dpm%ncar))
            allocate (dpm%RhoG(dpm%ncar,dpm%ncar))

            if (.not. GLOBAL_all_mode) then
              dpm%carac=carac_t
            else
              ! Il faut initialiser dpm%carac lors de la lecture du fichier des phenotypes
              dpm%carac="<NOT_INIT>"
            end if

            nuis=nuis_t
            int_qtl=int_qtl_t
            nb_qtl_def=nb_qtl_def_t
            dpm%natureY=natureY_t
            dpm%h2=h2_t
            dpm%RhoP=RhoP_t
            dpm%RhoG=RhoG_t
        end if

        deallocate(carac_t)
        deallocate(nuis_t)
        deallocate(int_qtl_t)
        deallocate(natureY_t)
        deallocate(nb_qtl_def_t)
        deallocate (h2_t)
        deallocate (RhoP_t)
        deallocate (RhoG_t)

    end subroutine fixe_structure_model
    !!***


    !!****f* m_qtlmap_traits/init_model_struct
    !!  NAME
    !!    init_model_struct
    !!  DESCRIPTION
    !!    fill the array model and listModelTrait variable
    !!  NOTES
    !!   * the array model is used by the m_qtlmap_modlin routines and derived
    !!   * the listModelTrait variable is used by the m_qtlmap_incidence routines and derived
    !!  SOURCE
    subroutine init_model_struct(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        integer                     :: alloc_stat
        integer                     :: i,j,ic,k,nf,nc,ifix,ico,nqf,nl,il,kc,iqtl
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        call log_mess('init_model_struct',INFO_DEF)

        allocate (dpm%modele(dpm%ncar,3+(3*dpm%nfix)+dpm%ncov), stat = alloc_stat)
        call check_allocate(alloc_stat,'model')
        allocate (dpm%listModelTrait(dpm%ncar))

        do ic=1,dpm%ncar
            call log_mess("Model - trait:"//dpm%carac(ic),DEBUG_DEF)


            allocate ( dpm%listModelTrait(ic)%indexFixedEffect(dpm%nfix))
            nf=0
            ! Pour les effets fixes et covariables
            do ifix=1,dpm%nfix
                if(nuis(ic,ifix) == 1) then
                    call log_mess("fix :"//trim(str(ifix)),DEBUG_DEF)
                    nf=nf+1
                    dpm%modele(ic,3+nf)=ifix
                    dpm%listModelTrait(ic)%indexFixedEffect(nf)=ifix
                end if
            end do
            dpm%listModelTrait(ic)%nbfe = nf
            dpm%modele(ic,1)=nf

            allocate ( dpm%listModelTrait(ic)%indexCovariate(dpm%ncov))
            nc=0
            do ico=1,dpm%ncov
                if(nuis(ic,dpm%nfix+ico) == 1) then
                    call log_mess("cov :"//trim(str(ico)),DEBUG_DEF)
                    nc=nc+1
                    dpm%modele(ic,3+nf+nc)=ico
                    dpm%listModelTrait(ic)%indexCovariate(nc)=ico
                end if
            end do
            dpm%modele(ic,2)=nc
            dpm%listModelTrait(ic)%nbco = nc
            ! Pour les interactions effets fixes, QTL
            nqf=0

            if ( nb_qtl_def(ic) >= 1 ) then
                do ifix=1,dpm%nfix
                    if(int_qtl(ic,1,ifix) == 1) then
                        call log_mess("interaction fix :"//trim(str(ifix)),DEBUG_DEF)
                        nqf=nqf+1
                        dpm%modele(ic,3+nf+nc+nqf)=ifix
                    end if
                end do
            end if
            dpm%modele(ic,3)=nqf

            allocate (dpm%listModelTrait(ic)%nbint(MAX_QTL))
            dpm%listModelTrait(ic)%nbint=0
            allocate ( dpm%listModelTrait(ic)%indexFixedEffectWithInteraction(nb_qtl_def(ic),dpm%nfix))
            do iqtl=1,nb_qtl_def(ic)
                nqf = 0
                do ifix=1,dpm%nfix
                    if(int_qtl(ic,iqtl,ifix) == 1) then
                        call log_mess("new model:interaction fix :"//trim(str(ifix)),DEBUG_DEF)
                        nqf=nqf+1
                        dpm%listModelTrait(ic)%indexFixedEffectWithInteraction(iqtl,nqf)=ifix
                    end if
                end do
                dpm%listModelTrait(ic)%nbint(iqtl)=nqf
            end do

        end do

        deallocate (nb_qtl_def)

    end subroutine init_model_struct
    !!***


    !!****f* m_qtlmap_traits/initialise_struct_internal
    !!  NAME
    !!    initialise_struct_internal
    !!  DESCRIPTION
    !!    allocation of main array of dataset user : corperf,niveau,covar,y,ycategorial,presentc,ndelta,cd
    !!  NOTES
    !!
    !!  SOURCE
    subroutine initialise_struct_internal(dataset,work)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        type(WORK_PHENOTYPE)       ,intent(inout)            :: work
        integer                     :: alloc_stat
        integer                     :: i,j,ic,k,nf,nc,ifix,ico,nqf,nl,il,icCategorial,kc,iqtl
        !dim : nd
        integer , dimension (:)  ,allocatable    :: nbvx
        !dim: ndx,nfix
        character(len=LEN_DEF) , dimension (:,:),allocatable   :: niveau
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg

        dpa => dataset%phenoAnimal
        dpm => dataset%phenoModel
        dg => dataset%genea

        call log_mess('initialise_struct_internal model and performance',INFO_DEF)

        allocate (dpa%bete(work%nbete))
        dpa%bete = work%bete(:work%nbete)

        allocate (dpa%corperf(dg%nd))
        allocate (dpa%corperf(dg%nd))
        allocate (niveau(dg%nd,dpm%nfix))
        allocate (dpa%covar(dg%nd,dpm%ncov))
        allocate (dpa%y(dpm%ncar,dg%nd))

        if ( dpm%ncarcat /=0 ) then
            allocate (dpa%ycategorial(dpm%ncarcat,dg%nd))
        end if

        allocate (dpa%presentc(dpm%ncar,dg%nd))
        allocate (dpa%ndelta(dpm%ncar,dg%nd))
        allocate (dpa%cd(dpm%ncar,dg%nd))


        dpa%y=REAL_NOT_DEFINED
        dpa%presentc=.false.
        dpa%covar=REAL_NOT_DEFINED
        niveau=STRING_NOT_DEFINED
        dpa%corperf=INT_NOT_DEFINED
        dpa%cd=0.d0
        !---------------------- corperf et presentc
        !for each descendant
        do i=1,dg%nd
            !find the animal associated
            do j=1,size(dpa%bete)
                if ( dg%animal(i) == dpa%bete(j) ) then
                    exit ! exit from this do
                endif
            end do

            if ( j >= size(dpa%bete)+1 ) then  !we does not find the animal
                cycle ! next descendant
            else ! animal find....
                dpa%corperf(i)=j
                icCategorial = 0
                do ic=1,dpm%ncar
                    if (dpm%natureY(ic) == 'c') icCategorial = icCategorial + 1
                    dpa%cd(ic,i)=work%cdt(ic,j)
                    dpa%ndelta(ic,i)=work%ndelt(ic,j)
                    if(work%cdt(ic,j) /= 0.d0) then
                        if (dpm%natureY(ic) == 'c') then
                            dpa%ycategorial(icCategorial,i)=str(work%valeur(ic,j))
                        else
                            dpa%y(ic,i)=work%valeur(ic,j)
                        end if
                        dpa%presentc(ic,i)=.true.
                    endif
                end do

                do ifix=1,dpm%nfix
                    niveau(i,ifix)=work%niv(j,ifix)

                    ! un phenotype manquant est considere lorsque l'effet fixe,
                    ! considere dans le modele, est non renseigne

                    if (niveau(i,ifix) == STRING_NOT_DEFINED) then
                        do ic=1, dpm%ncar
                            do k=1, dpm%modele(ic,1)
                                if(dpm%modele(ic,3+k) == ifix) then
                                    dpa%presentc(ic,i)=.false.
                                end if
                            enddo
                            do k=1, dpm%modele(ic,3)
                                if(dpm%modele(ic,3+dpm%modele(ic,1)+dpm%modele(ic,2)+k)==ifix) then
                                    dpa%presentc(ic,i)=.false.
                                end if
                            enddo
                        enddo
                    endif
                end do

                do ico=1,dpm%ncov
                    dpa%covar(i,ico)=work%cov(j,ico)

                    !  un phenotype manquant est considere lorsque la covariable
                    !  est non renseignee

                    if (dpa%covar(i,ico) == REAL_NOT_DEFINED) then
                        do ic=1, dpm%ncar
                            do k=1, dpm%modele(ic,2)
                                if (dpm%modele(ic,3+dpm%modele(ic,1)+k) == ico) then
                                    dpa%presentc(ic,i)=.false.
                                endif
                            enddo
                        enddo
                    endif
                end do
            end if
        end do
          !-----------------------------------------------------------------------------
          ! DET DU NOMBRE DE NIVEAUX DES EFFETS FIXES
          !-----------------------------------------------------------------------------
        allocate (dpm%nlev(dpm%nfix), stat = alloc_stat)
        call check_allocate(alloc_stat,'nlev')

        allocate (nbvx(dg%nd), stat = alloc_stat)
        call check_allocate(alloc_stat,'nbvx')

        allocate (dpa%nivx(dg%nd,dpm%nfix), stat = alloc_stat)
        call check_allocate(alloc_stat,'nivx')

        allocate (dpm%listelev(dpm%nfix,dg%nd), stat = alloc_stat)
        call check_allocate(alloc_stat,'listelev')

        do ifix=1,dpm%nfix
            nl=0
            dpm%nlev(ifix)=0
            do i=1,dg%nd
                nbvx(i)=0
            enddo
            !Comptage du nombre de fois qu'un niveau est rencontre pour le premier i et mise a
            ! 0 pour les autres
            do i=1,dg%nd
                if (niveau(i,ifix) /= STRING_NOT_DEFINED) then
                    if (nbvx(i) == 2) then
                        nbvx(i)=0
                        cycle
                    else
                        nbvx(i)=1
                        do j=i+1,dg%nd
                            if (niveau(i,ifix) == niveau(j,ifix)) then
                                nbvx(i)=nbvx(i)+1
                                nbvx(j)=2
                            end if
                        end do
                    end if
                end if
            end do

            do i=1,dg%nd
                if ( nbvx(i) /= 0 ) then
                    nl=nl+1
                    dpm%listelev(ifix,nl)= niveau(i,ifix)
                endif
            enddo

            dpm%nlev(ifix)=nl

            ! affectation des differents niveaux d'effets fixes renumerotes de 1 a nl
            do i=1,dg%nd
                if (niveau(i,ifix) == STRING_NOT_DEFINED) dpa%nivx(i,ifix)=0
                do il=1,dpm%nlev(ifix)
                    if (dpm%listelev(ifix,il) == niveau(i,ifix)) dpa%nivx(i,ifix)=il
                enddo
            enddo
        enddo

        deallocate(niveau)
        deallocate(nbvx)
        call log_mess('END SUBROUTINE initialise_struct_internal',DEBUG_DEF)
    end subroutine  initialise_struct_internal
    !!***


    !!****f* m_qtlmap_traits/read_perf_by_line
    !!  NAME
    !!    read_perf_by_line
    !!  DESCRIPTION
    !!    read the phenotypic file with the expression quantitative trait value format :
    !!    The header line is the list of animals phenotyped. The following line are the fixed effects, covariates and finally the phenotype.
    !!    The format of the nuisances effects and phenotype line is  :
    !!    <IDANIMAL> <FIXED_EFFECT1> <FIXED_EFFECT2>...<COV1> <COV2>...<VALUE1> <VALUE2>...
    !!    For missing data, insert a character string which is not interpretable as a numeric(e.g. n/a).
    !!    The following array are filling : bete, niv, cov, valeur, cdt, ndelt
    !!
    !!  INPUTS
    !!   files  : list of phenotypic files
    !!   filter : vector of index to select a subset of trait defined in the model file
    !!  NOTES
    !!
    !!  SOURCE
    subroutine read_perf_by_line(dataset,work)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        type(WORK_PHENOTYPE)       ,intent(inout)            :: work
        integer                     :: ios,alloc_stat
        integer                     :: nbete,current_number,l,i,j,futureNcar,k,ii
        character(len=LEN_DEF)      :: token
        character(len=LEN_LINE)     :: line_read
        logical                     :: is_ok
        character(len=1)            :: buff_read
        character(len=LENGTH_MAX_FILE) :: file
        character(len=LEN_BUFFER_WORD) :: name
        character(len=LEN_BUFFER_WORD),dimension(:),allocatable :: bete_char
        character(len=LEN_DEF) ,dimension(:),allocatable :: niv_t
        character(len=LEN_DEF) ,dimension(:),allocatable ::locvaleur
        real(kind=dp) ,dimension(:),allocatable :: cov_t
        real(kind=dp)                  :: v

        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        call log_mess('SUBROUTINE read_perf_by_line',DEBUG_DEF)
        call log_mess('reading transcriptom data file...',INFO_DEF)
        nbete = 0
        futureNcar = dpm%ncar
        if (dpm%n_filter_car /=0 ) futureNcar = dpm%n_filter_car

        call log_mess(' ** transcriptomic file : do not write number of animal at the first line ** ',WARNING_DEF)

        file = dataset%params%get_file_val(K_TRAITS)
        if  (dpm%ncar <= 0 .or. dpm%nfix < 0 .or. dpm%ncov < 0 ) then
            call stop_application('Devel error: call read_model before read_traits')
        end if

        open(unit_perf,FILE=file,IOSTAT=ios,FORM="formatted",recl=2**20,action="read",status="old")

        if ( ios /= 0 ) then
            call stop_application('The application can not open the traits file ['//trim(file)//']')
        endif

        line_read = ''
        l = 1
        ios = 0

        ! blank line....
        do while ( (ios == 0) .and. trim(line_read)=='')
            l = l + 1
            call GET(unit_perf,line_read,maxlen=LEN_LINE,IOSTAT=ios)
        end do

        if ( ios /= 0 ) then
            call stop_application('The application can not read the first line ['//trim(file)//'].'// &
                ' Split the file to resolve this error')
        endif

        is_ok = .true.

        do while ( is_ok )
            token = next_word(line_read,is_ok)
            nbete = nbete + 1
        end do

        nbete = nbete -1

        rewind(unit_perf)

        call log_mess('Number of animals detected:'//trim(str(nbete)),VERBOSE_DEF)

        allocate (work%bete(nbete), stat = alloc_stat)
        call check_allocate(alloc_stat,'bete')
        allocate (bete_char(nbete), stat = alloc_stat)
        call check_allocate(alloc_stat,'bete_char')
        allocate (work%niv(nbete,dpm%nfix), stat = alloc_stat)
        call check_allocate(alloc_stat,'niv')

        allocate(niv_t(nbete))

        allocate (work%cov(nbete,dpm%ncov), stat = alloc_stat)
        call check_allocate(alloc_stat,'cov')

        allocate (cov_t(nbete), stat = alloc_stat)

        allocate (locvaleur(nbete), stat = alloc_stat)
        allocate (work%valeur(futureNcar,nbete), stat = alloc_stat)
        call check_allocate(alloc_stat,'locvaleur')
        allocate (work%cdt(futureNcar,nbete), stat = alloc_stat)
        call check_allocate(alloc_stat,'cdt')
        allocate (work%ndelt(futureNcar,nbete), stat = alloc_stat)
        call check_allocate(alloc_stat,'ndelt')

        l = 1
        ! animal name header
        read(unit_perf,*,iostat=ios) (bete_char(i),i=1,nbete)
        work%bete = bete_char
        if ( ios /= 0 ) then
            call log_mess('Problem detected at the line:'//trim(str(l)),ERROR_DEF)
            call stop_application('No header (animals id) are finding in traits file:'//file)
        end if

        i = 1
        do while ( i <= dpm%nfix )
            l = l + 1
            read(unit_perf,*,iostat=ios) name,(niv_t(j),j=1,nbete)
            if ( trim(name) == '') cycle

            if ( ios /= 0 ) then
                call log_mess('Problem detected at the line:'//trim(str(l)),ERROR_DEF)
                call stop_application('The application can not initialize the fixed effect ['//&
                    trim(str(i))//']  file:['//trim(file)//'].')
            end if

            !on cherche l effet fixe correspondant
            do j=1,dpm%nfix
                if ( trim(dpm%namefix(j))==trim(name)) then
                    exit
                end if
            end do

            if ( j<=dpm%nfix) then
                work%niv(:,j)=niv_t(:)
            else
                call log_mess("QTLMap do not use fixed effect :"//trim(name),WARNING_DEF)
            end if

            i=i+1
        end do

        i = 1
        do while ( i <= dpm%ncov )
            l = l + 1
            read(unit_perf,*,iostat=ios) name,(cov_t(j),j=1,nbete)
            if ( trim(name) == '') cycle
            if ( ios /= 0 ) then
                call log_mess('Problem detected at the line:'//trim(str(l)),ERROR_DEF)
                call stop_application('The application can not initialize the covariate ['//&
                    trim(str(i))//']  file:['//trim(file)//'].')
            end if

             !on cherche l effet fixe correspondant
            do j=1,dpm%ncov
                if ( trim(dpm%namecov(j))==trim(name)) then
                    exit
                end if
            end do

            if ( j<=dpm%ncov) then
                work%cov(:,j)=cov_t(:)
            else
                call log_mess("QTLMap do not use Covariate :"//trim(name),WARNING_DEF)
            end if

            i=i+1
        end do

        i = 1
        ii=0
        ios=0

        do while ( i <= futureNcar )
            name=""
            ios=0
            do while ( ios == 0 .and. trim(name)=="")
                !reading line
                read(unit_perf,*,iostat=ios) name,(locvaleur(j),j=1,nbete)
            end do

            !No enought line defined in the file
            if ( ios /= 0 ) then
                call log_mess('Problem detected at the line:'//trim(str(l)),ERROR_DEF)
                call stop_application('The application can not initialize the trait ['//&
                    trim(str(i))//']  file:['//trim(file)//'].')
            end if

            if ( GLOBAL_filter_mode ) then
              ! recherche du caractere
              do k=1,futureNcar
                if ( trim(name) == trim(dpm%filter_car_name(k)) ) exit
              end do

              if (k>futureNcar) cycle
            else
              k = i
            end if

            !first colomn : name
            dpm%carac(k) = trim(name)

            do j=1,nbete
                work%valeur(k,j) = get_real(locvaleur(j),is_ok)
                if ( .not. is_ok ) then
                    call log_mess('** value unknown for trait ['//trim(dpm%carac(i))//&
                        '] of animal ['//trim(work%bete(j))//']',WARNING_DEF)
                    work%cdt(k,j)       = 0
                else
                    work%cdt(k,j)       = 1
                end if
            end do

            work%ndelt(k,:) = 1

            i = i + 1
        end do

        work%nbete = nbete

        close(unit_perf)
        deallocate (locvaleur)
        deallocate(bete_char)
        deallocate (niv_t)
        deallocate (cov_t)

        call log_mess('END SUBROUTINE read_perf_by_line',DEBUG_DEF)
    end subroutine read_perf_by_line
    !!***


    !!****f* m_qtlmap_traits/read_perf_by_column
    !!  NAME
    !!    read_perf_by_column
    !!  DESCRIPTION
    !!    read the phenotypic file with the format :
    !!    For each animal, its ID (identical to the ID given in the pedigree file) is followed by information
    !!    about nuisance effects (fixed effect levels, covariable value) and then by three information for each trait :
    !!      * the performance
    !!      * an 0/1 variable IP which indicates if (IP=1) or not (IP=0) the trait was measured for this animal and must be included in the analysis
    !!      *  and 0/1 variable (IC) which indicates if (IC=0) it was censored or not (IC=1), this IC information being needed for survival analysis (by default IC=1).
    !!
    !!    Grammar : <animal_id> <fixed_effect1> <fixed_effect2>... <cov1>..<cov2> <IP> <IC> <VALUE>
    !!
    !!    The following array are filling : bete, niv, cov, valeur, cdt, ndelt
    !!
    !!  INPUTS
    !!   files  : list of phenotypic files
    !!   filter : vector of index to select a subset of trait defined in the model file
    !!  NOTES
    !!
    !!  SOURCE
    subroutine read_perf_by_column(dataset,work)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        type(WORK_PHENOTYPE)       ,intent(inout)            :: work

        integer                        :: ios,alloc_stat
        integer                        :: nbete
        character(len=LEN_DEF)         :: vs,temp_vs
        character(len=LEN_LINE)        ::line_read
        character(len=LEN_BUFFER_WORD) :: word_token
        integer                     :: l,i,j,k
        logical                     :: is_ok,filterActive
        character(len=LEN_BUFFER_WORD) :: file
        integer :: futureNcar
        real (kind=dp), dimension (:,:), allocatable         :: cdt1,valeur1
        !dim : nc,nbete
        integer , dimension (:,:),allocatable                :: ndelt1

        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        call log_mess('SUBROUTINE read_perf_by_column',DEBUG_DEF)
        file = dataset%params%get_file_val(K_TRAITS)
        open(unit_perf,FILE=file,IOSTAT=ios,action="read",status="old")
        call log_mess('reading traits file...',INFO_DEF)
        if ( ios /= 0 ) then
            call stop_application('Cannot open the traits file ['//trim(file)//']')
        endif

        futureNcar = dpm%ncar
        filterActive=.false.
        if (dpm%n_filter_car /=0 ) then
            filterActive=.true.
            futureNcar= dpm%n_filter_car
        end if

        nbete=MAX_ANIMAL
        allocate (work%bete(nbete))
        allocate (work%niv(nbete,dpm%nfix))
        allocate (work%cov(nbete,dpm%ncov))
        allocate (valeur1(dpm%ncar,nbete))
        allocate (cdt1(dpm%ncar,nbete))
        allocate (ndelt1(dpm%ncar,nbete))

        ios = 0
        nbete = 0
        i=1
        rewind(unit_perf)

        do while ( ios >= 0 )
            read (unit_perf,*,iostat=ios) work%bete(i),(work%niv(i,j),j=1,dpm%nfix),(work%cov(i,j),j=1,dpm%ncov),&
                (valeur1(j,i),cdt1(j,i),ndelt1(j,i),j=1,dpm%ncar)
!            print *,ios
!            print *,trim(work%bete(i)),(work%niv(i,j),j=1,dpm%nfix),(work%cov(i,j),j=1,dpm%ncov),&
!                           (valeur1(j,i),cdt1(j,i),ndelt1(j,i),j=1,dpm%ncar),"**"
            if ( trim(work%bete(i))/='' .and. ( ios == 0) ) then
                nbete = nbete+1

                i = i + 1
            else if ( ios > 0 ) then
                call log_mess('Traits file - Bad definition of animal ['//trim(work%bete(i))//"]",ERROR_DEF)
                if ( dpm%nfix> 0 )&
                call log_mess("-"//trim(str(dpm%nfix))//" fixed effect(s)",ERROR_DEF)
                if ( dpm%ncov> 0 )&
                call log_mess("-"//trim(str(dpm%ncov))//" covariate(s)",ERROR_DEF)
                call log_mess("-"//trim(str(dpm%ncar))//" trait(s)",ERROR_DEF)
                call stop_application("check yours phenotype and model files.")
            endif
        end do

        close(unit_perf)

        if ( nbete == 0 ) then
            call stop_application("Can not read in the trait file:["//trim(file)//"] animal:["//trim(work%bete(i))//"]")
        end if

        allocate (work%valeur(futureNcar,nbete))
        allocate (work%cdt(futureNcar,nbete))
        allocate (work%ndelt(futureNcar,nbete))

        do i=1,nbete

            !           do j=1,(i-1)
            !              if ( bete(i) == bete(j)) then
            !                call stop_on_error (1,file,l,'animal ['//trim(bete(i))//'] have two lines definitions of traits !!')
            !              end if
            !           end do
            do j=1,dpm%ncar
                !             if ( futureNcar /= ncar ) then
                if ( filterActive ) then
                    do k=1,dpm%n_filter_car
                        if ( dpm%filter_car(k) == j ) exit
                    end do
                    ! the index car is not used
                    if ( k > dpm%n_filter_car) cycle
                else
                    k = j
                end if
                work%cdt(k,i) = cdt1(j,i)
                work%valeur(k,i) = valeur1(j,i)
                work%ndelt(k,i) = ndelt1(j,i)
            end do
        end do
        deallocate (valeur1)
        deallocate (cdt1)
        deallocate (ndelt1)

        work%nbete = nbete

        call log_mess('END SUBROUTINE read_perf_by_column',DEBUG_DEF)
    end subroutine read_perf_by_column
    !!***


    !!****f* m_qtlmap_traits/check_traits_and_fathers
    !!  NAME
    !!    check_traits_and_fathers
    !!  DESCRIPTION
    !!    check if number of progenies is greater than 1 otheriwse (stop the program with an error message)
    !!
    !!  NOTES
    !!
    !!  SOURCE
    subroutine check_traits_and_fathers(dataset)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer :: i,ic,nm1,nm2,nd1,nd2,jm,kd,effp

        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpa => dataset%phenoAnimal
        dpm => dataset%phenoModel
        dg => dataset%genea

        call log_mess("check traits of sire family...",INFO_DEF)

        if ( .not. associated(dg%nmp) ) then
            call stop_application("Dev error : try to read traits file before genealogy initialization...")
        end if

        do ic=1,dpm%ncar
            do i=1,dg%np
                nm1=dg%nmp(i)+1
                nm2=dg%nmp(i+1)
                effp=0
                do jm=nm1,nm2
                    nd1=dg%ndm(jm)+1
                    nd2=dg%ndm(jm+1)
                    effp = effp+count(dpa%presentc(ic,nd1:nd2))
                end do
                if (effp == 0) then
                    call stop_application('Father ['//trim(dg%pere(i))// &
                        '] have no child with the trait ['//trim(dpm%carac(ic))//']')
                end if
                if (effp == 1) then
                    call stop_application('Father ['//trim(dg%pere(i))// &
                        '] have only one child with the trait ['//trim(dpm%carac(ic))//']')
                end if
            end do
        end do
        call log_mess("end check traits...",INFO_DEF)
    end subroutine check_traits_and_fathers
    !!***

    !!****f* m_qtlmap_traits/sim_perf_tirage
    !!  NAME
    !!    sim_perf_tirage
    !!  DESCRIPTION
    !!
    !!  NOTES
    !!   Tirage des performances pour deux caracteres correles , avec ou sans qtl
    !!  SOURCE
    subroutine sim_perf_tirage(dataset,ncar,rho,h2,ys)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer , intent(in)                                     :: ncar
        !heritabilities of traits , dim : ncar
        real (kind=dp),dimension (ncar), intent(in)  :: h2
        ! correlation matrix , dim : ncar,ncar
        real (kind=dp),dimension (ncar,ncar) ,intent(in):: rho
        !dim : ncar
        real           ,dimension(ncar)           :: u0
        real (kind=dp) ,dimension(:),allocatable  :: upv,umv,udv,ugpv,ugmv

        !dim ncar,ncar
        real (kind=dp) ,dimension(:,:),allocatable :: varcovg,varcovgd
        real (kind=dp) ,dimension(:,:),allocatable :: varcove
        !dim ncar*ncar
        real           ,dimension(ncar,ncar)       :: ccovvar
        !dim : ncar,nr
        real (kind=dp) ,dimension(:,:),allocatable :: urv
        !dim : nd
        !integer ,dimension(:),allocatable  :: itest

        !dim : ncar,nd
        real (kind=dp) ,dimension(ncar,dataset%genea%nd),intent(out) :: ys
        !
        ! dim : nrsx=((ncar+1)*(ncar+2))/2)
        real  ,dimension(:),allocatable :: zr,zd,zh,zg

        real  ,dimension(ncar) :: WORK

        ! dim : ncar*(ncar+3)/2 + 1
        real  ,dimension(:),allocatable :: refgg,refg,refe

        integer       :: i,j,ngm1,ngm2,igm,nr1,nr2,ir,nm2
        integer       :: jm,nd1,nd2,kd,ifail,igp,ic,ip,nm1
        real (kind=dp) :: vare,vargd,varg,varej
        character(len=LEN_DEF) :: tempvs
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal


        call log_mess('simulation of traits...',VERBOSE_DEF)

        if (size(rho) <= 0 ) then
            call stop_application("DEVEL ERROR: sim_perf_tirage RHO is not initialized")
        end if

        if (size(h2) <= 0 ) then
            call stop_application("DEVEL ERROR: sim_perf_tirage H2 is not initialized")
        end if

        allocate (varcove(ncar,ncar))
        allocate (varcovg(ncar,ncar))
        allocate (varcovgd(ncar,ncar))
        allocate (upv(ncar))
        allocate (umv(ncar))
        allocate (udv(ncar))
        allocate (ugpv(ncar))
        allocate (ugmv(ncar))
        allocate (urv(ncar,dg%nr))


        allocate (zr( ((ncar+1)*(ncar+2))/2) )
        allocate (zd( ((ncar+1)*(ncar+2))/2) )
        allocate (zh( ((ncar+1)*(ncar+2))/2) )
        allocate (zg( ((ncar+1)*(ncar+2))/2) )

        allocate (refgg(ncar*(ncar+3)/2 + 1))
        allocate (refg(ncar*(ncar+3)/2 + 1))
        allocate (refe(ncar*(ncar+3)/2 + 1))



        ! call log_mess("Simulation trait only for real value",WARNING_DEF);
        ! natureY='r'
        !
        !*****************************************************************************
        !
        ! Initialisation des matrices de variance covariance residuelles et polygen
        !
        do i=1,ncar
            u0(i)=0.d0
            varg = sqrt(h2(i))
            vargd =sqrt(0.5d0*h2(i))
            vare=sqrt(1.d0-h2(i))
            varcovg(i,i)=varg*varg
            varcovgd(i,i)=vargd*vargd
            varcove(i,i)=vare*vare
            do j=i+1,ncar
                varej = sqrt(1.d0-h2(j))
                !           varcovg(i,j)=varg(i)*varg(j)
                varcovg(i,j)=0.0d0
                varcovg(j,i)=varcovg(i,j)
                !           varcovgd(i,j)=vargd(i)*vargd(j)
                varcovgd(i,j)=0.0d0
                varcovgd(j,i)=varcovgd(i,j)
                varcove(i,j)=rho(i,j)*vare*varej
                varcove(j,i)=varcove(i,j)
            end do
        end do
        ! 12   format(1x,f10.5,3x,10(f10.5,1x))
        !
        !
        !*****************************************************************************
        ! Initialisation des vecteurs de tirages des performances dans une binormale
        !
        do i=1,ncar
            do j=1,ncar
                ccovvar(i,j) = varcovg(i,j)
            end do
        end do

        call setgmn(u0,ccovvar,ncar,ncar,refgg)
        !
        do i=1,ncar
            do j=1,ncar
                ccovvar(i,j) = varcovgd(i,j)
            end do
        end do

        call setgmn(u0,ccovvar,ncar,ncar,refg)
        !
        do i=1,ncar
            do j=1,ncar
                ccovvar(i,j) = varcove(i,j)
            end do
        end do


        call setgmn(u0,ccovvar,ncar,ncar,refe)
        !
        !
        !*****************************************************************************
        ! Tirage des valeurs g�n�tiques des gparents et des reproducteurs
        !
        ifail=0
        ! Grd peres
        do igp=1,dg%ngp
            call genmn(refgg,zg,WORK)
            do ic=1,ncar
                ugpv(ic)=0.5d0*dble(zg(ic))
            end do
            !
            ! Grd meres
            ngm1=dg%ngmgp(igp)+1
            ngm2=dg%ngmgp(igp+1)
            do igm=ngm1,ngm2
                call genmn(refgg,zg,WORK)
                do ic=1,ncar
                    ugmv(ic)=0.5d0*dble(zg(ic))
                end do
                !
                ! Reproducteurs
                nr1=dg%nrgm(igm)+1
                nr2=dg%nrgm(igm+1)

                do ir=nr1,nr2
                    call genmn(refg,zr,WORK)
                    do ic=1,ncar
                        urv(ic,ir)=dble(zr(ic))+ugpv(ic)+ugmv(ic)
                    end do
                end do
            end do
        end do
        !*****************************************************************************
        ! Tirage des performances des descendants
        !
        do ip=1,dg%np
            if (dg%reppere(ip).eq.9999) call genmn(refgg,zg,WORK)
            do ic=1,ncar
                if (dg%reppere(ip).eq.9999) then
                    upv(ic)=0.5d0*dble(zg(ic))
                else
                    upv(ic)=0.5d0*urv(ic,dg%reppere(ip))
                end if
            end do
            nm1=dg%nmp(ip)+1
            nm2=dg%nmp(ip+1)
            do jm=nm1,nm2
                if (dg%repmere(jm).eq.9999) call genmn(refgg,zg,WORK)
                do ic=1,ncar
                    if (dg%repmere(jm).eq.9999) then
                        umv(ic)=0.5d0*dble(zg(ic))
                    else
                        umv(ic)=0.5d0*urv(ic,dg%repmere(jm))
                    end if
                end do
                nd1=dg%ndm(jm)+1
                nd2=dg%ndm(jm+1)
                do kd=nd1,nd2
                    call genmn(refg,zd,WORK)
                    do ic=1,ncar
                        udv(ic)=upv(ic)+umv(ic)+dble(zd(ic))
                    end do
                    !
                    ! Residuelle
                    call genmn(refe,zh,WORK)

                    do ic=1,ncar
                        ys(ic,kd)=udv(ic)+dble(zh(ic))
                    end do
                end do
            end do
        end do

        deallocate (urv)
        deallocate (varcove)
        deallocate (varcovg)
        deallocate (varcovgd)
        deallocate (upv)
        deallocate (umv)
        deallocate (udv)
        deallocate (ugpv)
        deallocate (ugmv)

        deallocate (zr)
        deallocate (zd)
        deallocate (zh)
        deallocate (zg)

        deallocate (refgg)
        deallocate (refg)
        deallocate (refe)

    end subroutine sim_perf_tirage

    !!   Tirage des performances pour deux caracteres correles , avec ou sans qtl
    subroutine sim_perf_interaction(dataset,spt,varp,ch1,snp1,ch2,snp2,effects_tab,outfile)
        type(QTLMAP_DATASET)       ,intent(in)                 :: dataset
        type(PDD_BUILD)            ,intent(in)                 :: spt
        real(kind=dp)              ,intent(in)                 :: varp
        integer                    ,intent(in)                 :: ch1,ch2,snp1,snp2
        real(kind=dp)              ,intent(in), dimension(4,4) :: effects_tab
        character(len=1024)        ,intent(in)                 :: outfile

        integer(kind=KIND_PHENO) :: p1,p2,snp1p1refs,snp1p2refs
        integer(kind=KIND_PHENO) :: snp2p1refs,snp2p2refs,ref1
        integer :: i1,i2,k
        logical :: isok
        real  ,external :: snorm
        !Effects_tab : description of effects according combinaison snp1,snp2

        ! SNP1/SNP2   11    12    21    22
        !      11     1,1   1,2   1,3   1,4
        !      12     2,1
        !      21     3,1
        !      22     4,1


        integer :: ip,jm,kd
        type(GENEALOGY_BASE) , pointer :: dg
        type(GENOTYPE_BASE)  , pointer :: geno

        real(kind=dp)  :: vare,vari,r,y

        dg => dataset%genea
        geno => dataset%genoAnimal

        open(345,file=trim(outfile))

        print *,"manque:",geno%nmanque

        print *,"Freq SNP1:",(geno%pc_all(ch1,snp1,k), k=1,geno%nall(ch1,snp1))
        print *,"Freq SNP2:",(geno%pc_all(ch2,snp2,k), k=1,geno%nall(ch2,snp2))

        ref1=0
        do ip=1,dg%np
            ! s
            !        snp1p1refs = spt%genotyp(ch1,snp1,geno%correp(ip),1)
            !        snp1p2refs = spt%genotyp(ch1,snp1,geno%correp(ip),2)
            !        snp2p1refs = spt%genotyp(ch2,snp2,geno%correp(ip),1)
            !        snp2p2refs = spt%genotyp(ch2,snp2,geno%correp(ip),2)

            print *,"pere ",dg%pere(ip),snp1p1refs,snp1p2refs,snp2p1refs,snp2p2refs

            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)

                !          snp1p1refs = spt%genotyp(ch1,snp1,geno%correm(jm),1)
                !          snp1p2refs = spt%genotyp(ch1,snp1,geno%correm(jm),2)
                !          snp2p1refs = spt%genotyp(ch2,snp2,geno%correm(jm),1)
                !          snp2p2refs = spt%genotyp(ch2,snp2,geno%correm(jm),2)
                !
                !          print *,"mere ",dg%mere(jm),snp1p1refs,snp1p2refs,snp2p1refs,snp2p2refs
                !          print *,"     ",geno%pheno(ch1,snp1,geno%correm(jm),1),geno%pheno(ch1,snp1,geno%correm(jm),2)

                do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
                    isok=.true.

                    y = snorm()*varp

                    if ( geno%presentg(ch1,kd) ) then
                        p1 = spt%genotyp(ch1,snp1,geno%corred(kd),1)
                        p2 = spt%genotyp(ch1,snp1,geno%corred(kd),2)

                        if ( p1 == geno%nmanque .or. p2 == geno%nmanque ) isok=.false.
                        if (isok) then
                         if ( ref1 == 0 ) ref1=p1!init
                         print *,p1,p2
                         if ( p1 == ref1 .and. p2 == ref1 ) then
                            i1 = 1
                         else if ( p1 == ref1 .and. p2 /= ref1 ) then
                            i1 = 2
                         else if ( p1 /= ref1 .and. p2 == ref1 ) then
                            i1 = 3
                         else
                            i1 = 4
                         end if
                        end if
                    else
                        isok=.false.
                    end if

                    if ( geno%presentg(ch2,kd) .and. isok ) then
                        p1 = spt%genotyp(ch2,snp2,geno%corred(kd),1)
                        p2 = spt%genotyp(ch2,snp2,geno%corred(kd),2)

                        if ( p1 == geno%nmanque .or. p2 == geno%nmanque ) isok=.false.

                        if (isok) then
                         if ( p1 == ref1 .and. p2 == ref1 ) then
                            i2 = 1
                         else if ( p1 == ref1 .and. p2 /= ref1 ) then
                            i2 = 2
                         else if ( p1 /= ref1 .and. p2 == ref1 ) then
                            i2 = 3
                         else
                            i2 = 4
                         end if
                        end if
                    else
                        isok=.false.
                    end if

                    if (isok) then
                       y = y + effects_tab(i1,i2)
                      print *,"kids ",kd,i1,i2
                    end if

                    write(345,*) dg%animal(kd),y,' 1 1 '

                end do
            end do
        end do

        close(345)

        stop


    end subroutine sim_perf_interaction


    !!****f* m_qtlmap_traits/sim_perf_shuffling
    !!  NAME
    !!    sim_perf_shuffling
    !!  DESCRIPTION
    !!
    !!  NOTES
    !!     realise un suffling intrafamille:
    !!      - de pere et mere lorsque le couple a plus de ndmin descendants
    !!      !-de pere, pour tous les descendants des meres ayant moins de ndmin
    !!      !descendants avec ce pere (l'ensemble etant > a ndmin)
    !!  SOURCE
    subroutine sim_perf_shuffling(dataset,analyse_is_multi_traits)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        logical                                     :: analyse_is_multi_traits
        real (kind=dp) ,dimension(:,:),allocatable   :: zy,zcd, zcovar
        integer ,dimension(:,:),allocatable          ::  zndelta,znivx
        integer ,dimension(:),allocatable          ::  effpr,iv
        logical ,dimension(:),allocatable          ::  permut_pere
        logical ,dimension(:,:),allocatable          ::  permut_pere_mere, zpresentc,zpresentg
        integer ,dimension(:,:),allocatable        ::   effpm
        integer                                    :: ifail,ip,jm,kd,nm1,nm2,nd1,nd2,i,j,id,nc, npresp
        integer                                   :: alloc_stat
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(GENEALOGY_BASE) , pointer :: dg
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENOTYPE_BASE) , pointer :: dga
        type(MAP_BASE) , pointer :: map

        dga => dataset%genoAnimal
        dpa => dataset%phenoAnimal
        dg => dataset%genea
        dpm => dataset%phenoModel
        map => dataset%map

        ALLOCATE (zy(dpm%ncar,dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zy')
        ALLOCATE (znivx(dg%nd,dpm%nfix), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'znivx')
        ALLOCATE (zcovar(dg%nd,dpm%ncov), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zcovar')
        ALLOCATE (zndelta(dpm%ncar,dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zndelta')
        ALLOCATE (zcd(dpm%ncar,dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zcd')
        ALLOCATE (zpresentc(dpm%ncar,dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zpresentc')
        ALLOCATE (permut_pere(dg%np), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'permut_pere')
        ALLOCATE (zpresentg(map%nchr,dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'zpresentg')
        ALLOCATE (permut_pere_mere(dg%np,dg%nm), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'permut_pere_mere')
        ALLOCATE (effpr(dg%np), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'effpr')
        ALLOCATE (iv(dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'iv')

        ALLOCATE (effpm(dg%np,dg%nm), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'effpm')

        call log_mess('simulation of traits by permutation (multitrait mode)...',VERBOSE_DEF)
        !*****************************************************************************
        !1�boucle sur les peres mere descendants
        !   --> comptage du nombre de descendants par couple pere-mere: effpm(ip,jm)
        !   --> comptage du nombre minimum de descendants d'un pere parmi
        !   les couple ayant plus de nb_des_min descendants: si ce nbre est inf�rieur � ndim
        !   alors ont donne un message d'erreur � l'execution


        zy=0.d0
        zndelta=0
        do ip=1,dg%np
            permut_pere(ip)=.false.
            effpr(ip)=0
            nm1=dg%nmp(ip)+1
            nm2=dg%nmp(ip+1)
            do jm=nm1,nm2
                effpm(ip,jm)=0
                nd1=dg%ndm(jm)+1
                nd2=dg%ndm(jm+1)
                permut_pere_mere(ip,jm)=.true.
                do kd=nd1,nd2
                    zy(:,kd)=dpa%y(:,kd)
                    zndelta(:,kd)= dpa%ndelta(:,kd)
                    zcd(:,kd)=dpa%cd(:,kd)
                    zpresentc(:,kd)=dpa%presentc(:,kd)
                    zpresentg(:,kd)=dga%presentg(:,kd)
                    znivx(kd,:)=dpa%nivx(kd,:)
                    zcovar(kd,:)=dpa%covar(kd,:)
                    if ((analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd))==dpm%ncar)) &
                        .or.(.not.analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd)).ne.0))) then
                        effpm(ip,jm)=effpm(ip,jm)+1
                        effpr(ip)=effpr(ip)+1
                    endif
                end do
                ! si une famille de pere-mere est trop petite les permutations se font en utilisant les performances de toute la famille de ce pere
                if (effpm(ip,jm)<dataset%params%ndmin.or.effpm(ip,jm)<nb_des_min) then
                    permut_pere(ip)=.true.
                    permut_pere_mere(ip,jm)=.false.
                    ! impression d'un warning lorsque l'utilisateur a permut� dans la famille de p�re au lieu de la famille pere mere
                    call log_mess ('WARNING: permutation will be performed within sire family for sire ['//trim(dg%pere(ip))//&
                        '], and dam ['//trim(dg%mere(jm))//']',VERBOSE_DEF)
                endif
            end do

            ! impression d'un warning lorsque des famille de pere sont trop petites pour �tre permut�es
            if (effpr(ip).LT.nb_des_min) then
                call log_mess ('family  size of sire ['//trim(dg%pere(ip))//'] is too small (['//&
                    str(effpr(ip))//'] animals) to use CONFIDENTLY permutation method, try the simulation method ',WARNING_DEF)
            endif
          !print *, 'EFFECTIF PERMUTE',effpr(ip)
        end do
        !

        !*****************************************************************************
                   !PERMUTATION DES PERFORMANCES
        !*****************************************************************************
        iv=0
        do ip=1,dg%np
            npresp=0
            nm1=dg%nmp(ip)+1
            nm2=dg%nmp(ip+1)
            !*****************************************************************************
            !   1--> r�attribution des performances par permutation au hasard intrafamille de pere( issus
            !   de MERE AVEC - de nb_des_min DESCENDANTS CHACUNE )
            IF (permut_pere(ip)) then
                do jm=nm1,nm2
                    nd1=dg%ndm(jm)+1
                    nd2=dg%ndm(jm+1)
                    do kd=nd1,nd2
                        if ((analyse_is_multi_traits.and.&
                        (count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd))==dpm%ncar)) &
                         .or.(.not.analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.&
                                                               (count(zpresentc(:,kd)).ne.0))) then
                            npresp=npresp+1
                            iv(npresp)=kd
                        endif
                    end do
                end do
                ifail=1
                !permutation au hasard dans iv
                call MATH_QTLMAP_G05EHF(iv,npresp,ifail)
                !permutation du vecteur de performances complet (sans donn�es manquantes) lorsqu'on fait une analyse multicaract�re
                !permutation du vecteur de performances incluant les donn�es manquantes lorsqu'on fait une analyse unicaract�re
                !(seuls les animaux avec aucune performance sont �cart�s
                id=0
                do jm=nm1,nm2
                    nd1=dg%ndm(jm)+1
                    nd2=dg%ndm(jm+1)
                    do kd=nd1,nd2
                        if (((analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd))==dpm%ncar)) &
                            .or.(.not.analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd)).ne.0))) &
                            .and..not.permut_pere_mere(ip,jm)) then
                            id=id+1
                            !  print *, 'kd', kd, 'Y', y(:,kd)
                            dpa%y(:,kd)=zy(:,iv(id))
                            ! print *, 'kd', kd, 'Y', y(:,kd)
                            dpa%ndelta(:,kd)=zndelta(:,iv(id))
                            dpa%cd(:,kd)=zcd(:,iv(id))
                            dpa%nivx(kd,:)=znivx(iv(id),:)
                            dpa%covar(kd,:)=zcovar(iv(id),:)
                            dpa%presentc(:,kd)=zpresentc(:,iv(id))
                            dga%presentg(:,kd)=zpresentg(:,iv(id))
                        endif
                    end do
                end do
               !   print*, 'EFFECTIF FINAL', id
            ENDIF
            !*****************************************************************************
            !   2--> permutation par famille de pere-mere
            do jm=nm1,nm2
                npresp=0
                id=0
                nd1=dg%ndm(jm)+1
                nd2=dg%ndm(jm+1)
                IF (permut_pere_mere(ip,jm)) THEN
                    do kd=nd1,nd2
                        if ((analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd))==dpm%ncar)) &
                            .or.(.not.analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.&
                                                                  (count(zpresentc(:,kd)).ne.0))) then
                            npresp=npresp+1
                            iv(npresp)=kd
                        endif
                    end do
                    ifail=1
                    !permutation au hasard dans iv
                    call MATH_QTLMAP_G05EHF(iv,npresp,ifail)
                    !permutation du vecteur de performances complet (sans donn�es manquantes) lorsqu'on fait une analyse multicaract�re
                    !permutation du vecteur de performances incluant les donn�es manquantes lorsqu'on fait une analyse unicaract�re
                    !(seuls les animaux avec aucune performance sont �cart�s
                    do kd=nd1,nd2
                        if ((analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.(count(zpresentc(:,kd))==dpm%ncar)) &
                            .or.(.not.analyse_is_multi_traits.and.(count(zpresentg(:,kd))>0).and.&
                                                                  (count(zpresentc(:,kd)).ne.0))) then
                            id=id+1
                            dpa%y(:,kd)=zy(:,iv(id))
                            dpa%ndelta(:,kd)=zndelta(:,iv(id))
                            dpa%cd(:,kd)=zcd(:,iv(id))
                            dpa%nivx(kd,:)=znivx(iv(id),:)
                            dpa%covar(kd,:)=zcovar(iv(id),:)
                            dpa%presentc(:,kd)=zpresentc(:,iv(id))
                            dga%presentg(:,kd)=zpresentg(:,iv(id))
                        endif
                    end do
                ENDIF
            end do
        !*****************************************************************************
        end do
        !*****************************************************************************
        ! IMPRESSION DES DONNEES PERMUTEES
        !*****************************************************************************
        ! do ic=1,ncar
        !      call log_mess('Trait ['//CHAR(carac(ic))//'] qtlmap normalize data with means:'//&
        !       str(xmut(ic))//'+-'//str(sigt(ic)),INFO_DEF)
        ! enddo
        !    do ip=1,np
        !	nm1=nmp(ip)+1
        !nm2=nmp(ip+1)
        !do jm=nm1,nm2
        !   nd1=ndm(jm)+1
        !   nd2=ndm(jm+1)
        !	   print *, 'permute pere',permut_pere(ip)
        !	   print *, 'permute pere-mere',permut_pere_mere(ip,jm), effpm(ip,jm), ndmin,nb_des_min
        !    do kd=nd1,nd2
        !        print *, trim(animal(kd)), ' ', trim(pere(ip)),' ',  trim(mere(jm)),&
        !	' data ',(trim(carac(ic)),' ',zy(ic,kd), zpresentc(ic,kd),ic=1,ncar),(znivx(kd,i), i=1,nfix), (zcovar(kd, j), j=1, ncov), &
        !        ' sim1 ', (trim(carac(ic)),' ',y(ic,kd), presentc(ic,kd),ic=1,ncar), (nivx(kd,i), i=1,nfix), (covar(kd, j), j=1, ncov)
        !     enddo
        !enddo
        !  enddo
        DEALLOCATE (zy)
        DEALLOCATE (znivx)
        DEALLOCATE (zcovar)
        DEALLOCATE (zndelta)
        DEALLOCATE (zpresentc)
        DEALLOCATE (effpr)
        DEALLOCATE (iv)
        DEALLOCATE (permut_pere)
        DEALLOCATE (permut_pere_mere)
        DEALLOCATE (zcd)
        DEALLOCATE (effpm)
        return

    end subroutine sim_perf_shuffling
    !!***


    !!****f* m_qtlmap_traits/set_estime
    !!  NAME
    !!    set_estime
    !!  DESCRIPTION
    !!    Initialize the following data :
    !!      * estime  (dim : ncar,nm)   : indicate the estimability of a full sib family (condition ndmin involve full sib family are used in the qtl analysis)
    !!      * nmumest (dim ncar)        : number of full sib family with enough progenies (ndmin>=)
    !!      * namest  (dim ncar)        : number of unique female (several full sib family are building with the same dam)
    !!      * iam     (dim : ncar,nfem) : get the index of female (estimable)
    !!
    !!  NOTES
    !!     anciennement dans optinit, on initialisse estime,nmumest,namest,iam dans cette routine
    !!  SOURCE

    subroutine set_estime(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        logical         :: estfem(dataset%genea%nfem)
        integer :: ic,ip,jm,kd,eff,ifem
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(GENEALOGY_BASE) , pointer :: dg
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        call log_mess('SUBROUTINE set_estime',DEBUG_DEF)

        if ( associated (dpa%estime) ) then
            deallocate( dpa%estime )
            deallocate( dpa%namest )
            deallocate( dpa%nmumest )
        end if

        if ( associated (dpa%iam) )  deallocate( dpa%iam )

        allocate (dpa%estime(dpm%ncar,dg%nm) )
        allocate (dpa%namest(dpm%ncar) )
        allocate (dpa%nmumest(dpm%ncar))
        allocate (dpa%iam(dpm%ncar,dg%nfem) )

        dpa%estime=.false.
        dpa%nmumest=0
        dpa%namest=0
        dpa%iam=0.d0

        do ic=1,dpm%ncar
            estfem=.false.
            do ip=1,dg%np
                do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                    eff=0
                    do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
                        if(dpa%presentc(ic,kd)) eff=eff+1
                    end do
                    ifem=dg%repfem(jm)
                    if(eff.ge.dataset%params%ndmin) then
                        dpa%estime(ic,jm)=.true.
                        estfem(ifem)=.true.
                        dpa%nmumest(ic)=dpa%nmumest(ic)+1
                    end if
                end do
            end do
            do ifem=1,dg%nfem
                dpa%iam(ic,ifem)=0
                if(estfem(ifem)) then
                    dpa%namest(ic)=dpa%namest(ic)+1
                    dpa%iam(ic,ifem)=dpa%namest(ic)
                end if
            end do
        end do

        call log_mess('END SUBROUTINE set_estime',DEBUG_DEF)

    end subroutine set_estime
    !!***


    !!****f* m_qtlmap_traits/check_cd
    !!  NAME
    !!    check_cd
    !!  DESCRIPTION
    !!    set cd to 1.0 if the trait is not 'average' and cd /= 0 for each animal
    !!
    !!  NOTES
    !!
    !!  SOURCE
    subroutine check_cd(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        integer :: i,j
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dpa => dataset%phenoAnimal

        call log_mess("check cd...",INFO_DEF)
        do i=1,dpm%ncar
            if ( dpm%natureY(i) /= 'a' ) then
                do j=1,size(dpa%cd,2)
                    if (dpa%cd(i,j)/=0.d0) dpa%cd(i,j)=1.d0
                end do
            end if
        end do ! i

    end subroutine check_cd
    !!***

    !!****f*  m_qtlmap_traits/normalize_data
    !!  NAME
    !!    normalize_data
    !!  SYNOPSIS
    !!    Normalizing data for continue traits
    !!  FUNCTION
    !!    Compute xmut (Mean) and sigt (Standart deviation) and normalize Y array.
    !!  INPUTS
    !!    * is_transcriptom    -- boolean to set for no printing information whith transcriptomic data
    !!    * is_simul           -- boolean to set for not printing information in simulation case
    !!  RESULT
    !!    none.
    !!  EXAMPLE
    !!
    !!  NOTES
    !!  If none traits are found found for one founder, the analyse stop
    !!  BUGS
    !!
    !!  SEE ALSO
    !!  xmut
    !!  sigt
    !!  natureY
    !!  MATH_QTLMAP_G01FAF
    !!
    !!  SOURCE
    subroutine normalize_data(dataset,is_transcriptom,is_simul)
        type(QTLMAP_DATASET)       ,intent(inout)      :: dataset
        logical  , intent(in)                          :: is_transcriptom
        logical  , intent(in)                          :: is_simul
        integer                                        :: kd,i
        real (kind=dp)                                 :: sy,sy2,eff
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        call log_mess('SUBROUTINE normalize_data',DEBUG_DEF)

        dpm%xmut = 0.d0
        dpm%sigt = 1.d0


        do i=1,dpm%ncar
            if (.not.(dpm%natureY(i) == 'r' .or. dpm%natureY(i) == 'a')) cycle
            sy=0.d0
            sy2=0.d0
            eff=0.d0
            dpm%xmut(i)=0.d0
            do kd=1,dg%nd
                if(dpa%presentc(i,kd)) then
                    sy2=sy2+dpa%y(i,kd)*dpa%y(i,kd)
                    sy=sy+dpa%y(i,kd)
                    eff=eff+1.d0
                endif
            enddo
            if (eff == 0.d0) then
                call stop_application("Can not find animal with performance for trait :"//trim(dpm%carac(i)))
            end if

            dpm%xmut(i)=sy/eff
            dpm%sigt(i)=dble(sqrt((sy2/eff)-dpm%xmut(i)*dpm%xmut(i)))

            !  normalisation des donn�es
            if ( .not. is_transcriptom .and. .not. is_simul ) then
                call log_mess('Trait ['//trim(dpm%carac(i))//'] qtlmap normalize data with means:'//&
                    str(dpm%xmut(i))//'+-'//str(dpm%sigt(i)),INFO_DEF)
            end if

            do kd=1,dg%nd
                if(dpa%presentc(i,kd)) then
                    dpa%y(i,kd)=(dpa%y(i,kd)-dpm%xmut(i))/dpm%sigt(i)
                end if
            end do
        end do

        call log_mess('END SUBROUTINE normalize_data',DEBUG_DEF)
    end subroutine normalize_data
    !!***

    !!****f*  m_qtlmap_traits/set_count_discrete
    !!  NAME
    !!    set_count_discrete
    !!  SYNOPSIS
    !!    Manage discrete data to initialized associated structure
    !!  FUNCTION
    !!    Initialize structure ydiscretord , nmod and indicemod for discrete data
    !!  INPUTS
    !!
    !!  RESULT
    !!    none.
    !!  EXAMPLE
    !!
    !!  NOTES
    !!
    !!  BUGS
    !!
    !!  SEE ALSO
    !!  ydiscretord
    !!  nmod
    !!  indicemod
    !!
    !!  SOURCE
    subroutine set_count_discrete(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        logical  :: found
        integer  :: i,ic,ideb,m,m1,j,temp
        integer , dimension(dataset%phenoModel%ncar,dataset%genea%nd) :: indicemod_t
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea
        dpm => dataset%phenoModel
        dpa => dataset%phenoAnimal

        if (.not. associated(dpa%ydiscretord)) allocate (dpa%ydiscretord(dpm%ncar,dg%nd))
        if (.not. associated(dpm%nmod)) allocate (dpm%nmod(dpm%ncar))
        dpa%ydiscretord = 0

        !******************************************************************************
        !******************************************************************************
        !
        ! transformation des donn�es discret � des donn�e utilisable par QTLMAP
        !
        !
        !*****************************************************************************
        !******************************************************************************
        !


        !  Parametres de maximisation

        !  nmod est le nombre des modalit�s du caract�re �tudi�
        !  on reserve donc la partie du vecteur param�tre qui contiendra le lambda � estimer
        dpm%nmod=0

        do ic=1,dpm%ncar
            if (dpm%natureY(ic) /= 'i') cycle
            dpm%nmod(ic)=1
            i=1
            found=.false.
            do while(i.le.dg%nd.and..not.found)
                if(dpa%presentc(ic,i)) then
                    found=.true.
                    ideb=i+1
                    indicemod_t(ic,1) = int(dpa%y(ic,i))
                end if
            end do

            do i=ideb,dg%nd
                if(dpa%presentc(ic,i))then
                    m=1
                    found=.false.
                    do while (m<=dpm%nmod(ic) .and. .not. found)
                        if (int(dpa%y(ic,i))==indicemod_t(ic,m)) then
                            found=.true.
                        else
                            m=m+1
                        end if
                    enddo

                    if (.not. found) then
                        dpm%nmod(ic)=dpm%nmod(ic)+1
                        indicemod_t(ic,dpm%nmod(ic))=dpa%y(ic,i)
                    end if
                end if
            enddo


            ! maintenant on tri le tableau des modalit�s du plus petit au plus grand
            !
            !*************************************************************************
            !**************************tri a bulle************************************
            !*************************************************************************
            do j=1, dpm%nmod(ic)
                do m1=1, j-1
                    if (indicemod_t(ic,j)<indicemod_t(ic,m1)) then
                        temp=indicemod_t(ic,m1)
                        indicemod_t(ic,m1)=indicemod_t(ic,j)
                        indicemod_t(ic,j)=temp
                    end if
                end do
            end do
            !**************************************************************************
            !
            ! nmod est le nombre de modalit� possible pour le caract�re �tudi�
            !
            ! maintenent il faut donc cr�er le nouveau tableau de performances qui sera
            !utilis� dans QTLMap, ydiscretord
            !
            do i=1,dg%nd
                do j=1,dpm%nmod(ic)
                    if (dpa%presentc(ic,i)) then
                        if ( int(dpa%y(ic,i))==indicemod_t(ic,j) ) dpa%ydiscretord(ic,i)=j
                    end if
                end do
            end do

        end do ! ic



        if (.not.associated(dpm%indicemod)) allocate (dpm%indicemod(dpm%ncar,maxval(dpm%nmod)))
        dpm%indicemod = 0
        do ic=1,dpm%ncar
            do j=1,dpm%nmod(ic)
                dpm%indicemod(ic,j) = indicemod_t(ic,j)
            end do
        end do

    end subroutine set_count_discrete
    !!***

    !!****f*  m_qtlmap_traits/set_proportion_discrete
    !!  NAME
    !!    set_proportion_discrete
    !!  SYNOPSIS
    !!    Mange process according to the type of the data (continue, discrete or categorial)
    !!  FUNCTION
    !!    Initialize structure Prop (proportion) and Seuil (Threshold) for discrete data
    !!  INPUTS
    !!
    !!  RESULT
    !!    none.
    !!  EXAMPLE
    !!
    !!  NOTES
    !!
    !!  BUGS
    !!
    !!  SEE ALSO
    !!    prop
    !!    seuil
    !!    MATH_QTLMAP_G01FAF
    !!
    !!  SOURCE
    subroutine set_proportion_discrete(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        integer       :: ic,i,ifail,ii
        real(kind=dp) :: cumul,eff
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        ifail = 0

        !
        !  comptage pour la cr�ation du point de d�part
        !
        if(.not. associated(dpm%seuil)) allocate (dpm%seuil(dpm%ncar,maxval(dpm%nmod)))
        if (.not. associated(dpm%prop)) allocate (dpm%prop(dpm%ncar,maxval(dpm%nmod)))

        dpm%prop=0.d0
        dpm%seuil = 0.d0

        do ic=1,dpm%ncar
            if (dpm%natureY(ic) /= 'i') cycle
            eff=0.d0
            do ii=1,dg%nd
                if(dpa%presentc(ic,ii)) then
                    dpm%prop(ic,dpa%ydiscretord(ic,ii))=dpm%prop(ic,dpa%ydiscretord(ic,ii))+1.d0
                    eff=eff+1.d0
                end if
            end do
            do ii=1,dpm%nmod(ic)
                dpm%prop(ic,ii) = dpm%prop(ic,ii)/eff
            end do

            cumul=0.d0
            do i=1,dpm%nmod(ic)-1
                cumul=cumul+dpm%prop(ic,i)
                dpm%seuil(ic,i)=MATH_QTLMAP_G01FAF('L',cumul,ifail)
            end do
        !
        !  fin du comptage
        end do ! ic
    end subroutine set_proportion_discrete
    !!***

    !!****f*  m_qtlmap_traits/manage_data
    !!  NAME
    !!    manage_data
    !!  SYNOPSIS
    !!    Mange process according to the type of the data (continue, discrete or categorial)
    !!  FUNCTION
    !!    Normalize data for continue data and initilialized data structure sigt, xmut
    !!    Count discrete data for discrete data and compute Proportions (prop structure) and Threshold (seuil structure)
    !!  INPUTS
    !!    is_transcriptom    -- boolean to set for no printing information whith transcriptomic data
    !!    is_simul           -- boolean to set for not printing information in simulation case
    !!  RESULT
    !!    none.
    !!  EXAMPLE
    !!
    !!  NOTES
    !!
    !!  BUGS
    !!
    !!  SEE ALSO
    !!    is_transcriptom,
    !!    set_count_discrete,
    !!    set_proportion_discrete
    !!
    !!  SOURCE
    subroutine manage_data(dataset,is_transcriptom,is_simul,normalize)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        logical  , intent(in)                          :: is_transcriptom
        logical  , intent(in)                          :: is_simul
        logical  , intent(in)                          :: normalize

        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        if ( .not. associated(dpm%xmut) )allocate(dpm%xmut(dpm%ncar))
        if ( .not. associated(dpm%sigt) )allocate(dpm%sigt(dpm%ncar))

        if ( normalize ) then
            call normalize_data(dataset,is_transcriptom,is_simul)
        else
            dpm%xmut=0.d0
            dpm%sigt=1.d0
        end if

        call set_count_discrete(dataset)
        call set_proportion_discrete(dataset)
    end subroutine manage_data
    !!***





    !!****f*  m_qtlmap_traits/init_random
    !!  NAME
    !!    init_random
    !!  DESCRIPTION
    !!
    !!
    !!  NOTES
    !!   faire un module random
    !!  SOURCE
    subroutine init_random()
        integer :: n,clock,i
        integer, dimension(:), allocatable :: seed

        !initialisation pour les routines intrasec de fortran
        call random_seed(size=n)
        allocate(seed(n))
        call system_clock(count=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)
        deallocate(seed)

    end subroutine init_random
    !!***


    !!****f*  m_qtlmap_traits/normal
    !!  NAME
    !!    normal
    !!  DESCRIPTION
    !!
    !!
    !!  NOTES
    !!   faire un module random
    !!  SOURCE
    function normal(mean,sigma) !returns a normal distribution
        real(kind=dp)  , intent(in) :: mean,sigma
        real(kind=dp) normal,tmp
        integer flag
        real(kind=dp) fac,gsave,rsq,r1,r2,u
        save flag,gsave
        data flag /0/
        !$omp threadprivate (flag)
        if (flag.eq.0) then
            rsq=2.0
            do while(rsq.ge.1.0.or.rsq.eq.0.0) ! new from for do
                call random_number(u)
                r1=2.0*u-1.0
                call random_number(u)
                r2=2.0*u-1.0
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal
    !!***

      !! Fait par jm en australie.....

    subroutine create_trait(dataset,genin,perfin,perfout,nqtl,hdeux,position_qtl,effet_qtl)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        character(len=LENGTH_MAX_FILE)                         :: genin,perfin,perfout
        integer      ,            intent(in)                   :: nqtl
        real(kind=dp)       , intent(in)                       :: hdeux
        integer      ,    dimension(nqtl),        intent(in)   :: position_qtl ! chaque qtl est position sur un indice de marqueur
        real (kind=dp) ,dimension(nqtl)   ,        intent(in)  :: effet_qtl ! effet par qtl

        real   :: vh2(2,2),u0(2),refgg(6),z(6),work(2)
        real   :: gennor


        real (kind=dp)   :: up,um,ud,varE,varP
        real (kind=dp) , dimension(8000)   :: perf_d
        integer code1(8000),code2(8000),a1(8000)
        integer nkd,ip,jm,kd,ld,i,nh2,ii
        integer(kind=KIND_PHENO)   ,dimension (:,:),allocatable  ::   gene_o_qtl

        character*4  num(8000)
        character a2(8000)
        real val
        logical found

        character*7 snp(5111)
        integer m1(5111),m2(5111),i_ind
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        call init_random()

        open(111,file=genin)
        read(111,*)(snp(i),i=1,100)
        print*,'position_qtl',position_qtl
        i_ind=1

        allocate (gene_o_qtl(nqtl,8000))

13      read(111,*,end=14)num(i_ind),(m1(i),m2(i), i=1,100)
        do ii=1,nqtl
            i=position_qtl(ii)
            gene_o_qtl(ii,i_ind)= m1(i) + m2(i) -1
            i_ind=i_ind+1
        end do
        go to 13
14  continue
    close(111)

    !      vh2=0.d0
    !      vh2(1,1)=0.15d0
    !      vh2(2,2)=0.85d0
    !      u0=0.d0
    !      nh2=2
    !      call setgmn(u0,vh2,nh2,refgg)
    !      do ip=1,iseed
    !        call genmn(refgg,z,work)
    !      end do

    varP = 0.25d0*hdeux
    varE = 1 - varP


    open(111,file=perfin)
    open(222,file=perfout, form="formatted",recl=26)

    nkd=1
1   read(111,*,end =2)num(nkd),a1(nkd),a2(nkd),val,code1(nkd),code2(nkd)
    nkd=nkd+1
    go to 1
2 continue

  do ip=1,dg%np
      up=normal(0.d0,sqrt(varP))
      do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          um=normal(0.d0,sqrt(varP))
          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
              ud=normal(0.d0,sqrt(varE))
              !       print*,kd,corred(kd),gene_o_qtl(corred(kd))
              perf_d(kd)=up+um+ud
              do ii=1,nqtl
                  print *,( ( real(gene_o_qtl(ii,dga%corred(kd)) - 1.d0 ) / 2.d0 ) )
                  perf_d(kd)= perf_d(kd) + effet_qtl(ii)*( ( real(gene_o_qtl(ii,dga%corred(kd)) - 1.d0 ) / 2.d0 ) )
              end do
          end do
      end do
  end do

  do kd=1,nkd-1
      found=.false.
      do ld=1,dg%nd
          if(dg%animal(ld) ==num(kd)) then
              found=.true.
              exit
          end if
      end do
      !       print*,kd,ld
      write(222,2222)num(kd),a1(kd),a2(kd),perf_d(ld),code1(kd),code2(kd)
2222  format(a4,1x,i4,1x,a1,1x,F10.5,1x,i1,1x,i1)
  end do
  close(111)
  rewind(222)
  close(222)
  deallocate (gene_o_qtl)
  !    stop
  end subroutine create_trait

  subroutine release_work_phenotype(work)
      class(WORK_PHENOTYPE) , intent(inout) :: work

      deallocate (work%bete)
      deallocate (work%ndelt)
      deallocate (work%cdt)
      deallocate (work%valeur)
      deallocate (work%niv)
      deallocate (work%cov)
  end subroutine


  subroutine create_dataset_phenotype(dataset,occurences,array_sample,newdataset)
      type(QTLMAP_DATASET)             ,intent(in)            :: dataset
      integer     ,intent(in)  , dimension(dataset%genea%nd)  :: occurences
      integer             ,dimension(:),intent(in)            :: array_sample
      type(QTLMAP_DATASET)             ,intent(inout)         :: newdataset

      type(WORK_PHENOTYPE) :: work
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm
      integer :: i,j,idd,io,ic

      dpm => dataset%phenoModel
      dpa => dataset%phenoAnimal

      call dataset%phenoModel%copy(newdataset%phenoModel)

      work%nbete = size(dpa%bete)

      !initialisation pour les redondances d'individus
      do i=1,dataset%genea%nd
          if ( occurences(i) > 1 ) then
              work%nbete = work%nbete + occurences(i) - 1
          end if
      end do

      allocate (work%bete(work%nbete))
      work%bete(:size(dpa%bete)) = dpa%bete
      allocate (work%cdt(dpm%ncar,work%nbete))
      work%cdt(:,:size(dpa%cd,2)) = dpa%cd
      allocate (work%ndelt(dpm%ncar,work%nbete))
      work%ndelt(:,:size(dpa%ndelta,2)) = dpa%ndelta
      allocate (work%valeur(dpm%ncar,work%nbete))
      work%valeur(:,:size(dpa%y,2)) = dpa%y

      allocate (work%niv(work%nbete,dpm%nfix))
      do i=1,size(dpa%nivx,1)
          do j=1,size(dpa%nivx,2)
              work%niv(i,j) = trim(str(dpa%nivx(i,j)))
          end do
      end do

      allocate (work%cov(work%nbete,dpm%ncov))
      work%cov(:size(dpa%covar,1),:) = dpa%covar

      idd=dataset%genea%nd
      do i=1,newdataset%genea%nd
          if ( occurences(i) > 1 ) then
              do io=2,occurences(i)
                  idd=idd+1
                  work%bete(idd)=trim(dataset%genea%animal(i))//"_"//trim(str(io))
                  work%cdt(:,idd)=dpa%cd(:,i)
                  work%ndelt(:,idd) = dpa%ndelta(:,i)
                  work%valeur(:,idd) = dpa%y(:,i)
                  do j=1,size(dpa%nivx,2)
                      work%niv(idd,j) = trim(str(dpa%nivx(idd,j)))
                  end do
                  work%cov(idd,:) = dpa%covar(i,:)
              end do
          end if
      end do

      !        if ( present(calculCd) ) then
      !         if ( calculCd ) then
      !          do i=1,size(dpa%bete)
      !           call init_perf_animal(dataset,datasetUser,work%bete(i),work%valeur(:,i),work%cdt(:,i))
      !          end do
      !        end if
      !       end if

      call initialise_struct_internal(newdataset,work)
      call check_traits_and_fathers(newdataset)
      call set_estime(newdataset)

      call work%release()

  end subroutine create_dataset_phenotype

  !!****f* m_qtlmap_traits/write_perf
  !!  NAME
  !!    write_perf
  !!  DESCRIPTION
  !!
  !!  NOTES
  !!
  !!  SOURCE
  subroutine write_perf(dpa,dataset,file_name)
      class(PHENOTYPE_BASE) , intent(in) :: dpa
      type(QTLMAP_DATASET)   ,intent(in) :: dataset
      character(len=*),intent(in)     :: file_name
      integer :: kd,ic
      integer :: myunit=99999
      integer ,dimension(dataset%phenoModel%ncar)    :: cd
      real(kind=dp) ,dimension(dataset%phenoModel%ncar)    :: y
      character(len=10) :: nc
      type(GENEALOGY_BASE) ,pointer :: dg
      type(DATAMODEL_BASE) ,pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea

      write (nc,fmt='(i7)') dpm%ncar

      open (myunit,file=file_name)

      do kd=1,dg%nd
          cd = 1
          y = dpa%y(:,kd)
          do ic=1,dpm%ncar
              if(.not. dpa%presentc(ic,kd)) then
                  y(ic)=-99.d0
                  cd(ic) = 0
              end if
          end do

          write(myunit,FMT='(a12,'// trim(nc) //'(1x,f9.5,1x,i1,1x,"1"))')&
              trim(dg%animal(kd)),(dpa%y(ic,kd),cd(ic),ic=1,dpm%ncar)
      end do
      close(myunit)

  end subroutine write_perf
  !!***



  end module m_qtlmap_phenotype
