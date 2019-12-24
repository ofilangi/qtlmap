!!****m* INPUT/m_qtlmap_simulation
!!  NAME
!!    m_qtlmap_simulation -- Simulation routines.
!!  SYNOPSIS

!!  DESCRIPTION
!!
!!  NOTES
!!
!!  BUGS
!!
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/
module m_qtlmap_simulation

  use m_qtlmap_base
  use m_qtlmap_output_handler
  use m_qtlmap_log
 ! use m_qtlmap_map
 ! use m_qtlmap_genealogy
 ! use m_qtlmap_genotype
  use m_qtlmap_phenotype
  use m_qtlmap_types

  implicit none

  !PARAMETERS

!!****d* m_qtlmap_simulation/ios_simul_file
!!  NAME
!!   ios_simul_file
!!  DESCRIPTION
!!   Unit fortran record to read the simulation file
!!***
  integer                     ,parameter    :: ios_simul_file  = 40
!!****v* m_qtlmap_simulation/current_line
!!  NAME
!!   ios_simul_file
!!  DESCRIPTION
!!   number of the current line while the parse
!!***
  integer                                   :: current_line    = 0


!!****t* QTLMAP_TYPES/SIMULATION_INFO
!!  NAME
!!     SIMULATION_INFO
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type ,public :: SIMULATION_INFO
     logical                      ,public      :: simTyp                 = .false. !!   flag to activate the simulation genotype
     logical                      ,public      :: simulMap               = .false. !!   flag to activate a simulation map
     logical                      ,public      :: simulTraits            = .false. !!   flag to activate a simulation trait
     logical                      ,public      :: simulGenea             = .false. !!   flag to activate a simulation genealogy
     logical                      ,public      :: simulQtl               = .false. !!   flag to activate a simulation with a QTL

     integer                      ,public      :: nqtlsimul              = 0
     !(np+nfem+nd,nqtl,2)
     integer      ,dimension(:,:,:)  ,pointer  :: qtl                    => null()

     character(len=LEN_DEF)       ,public      :: croisement

     ! marker description
     !-------------------
     integer                      ,public      :: nalle
     real (kind=dp)               ,public      :: dens
     real (kind=dp)               ,public      :: taille

     ! Genealogy description
     !-------------------
     integer                      ,public      :: nbpere
     integer                      ,public      :: inmp
     integer                      ,public      :: indm

       ! qtl effect ,dim  ncar,nqtl
     real (kind=dp), dimension (:,:), pointer ,public :: ue                => null()
     ! dim : 1er trait, 2nd modality
     real (kind=dp),dimension(:,:)  , pointer ,public :: soglia             => null()! user threhold
     ! dim : 1er trait,2nd modality
     real (kind=dp),dimension(:,:)  , pointer ,public :: freqmodsimultrait  => null()
     ! dim : trait
     integer,dimension(:)         , pointer ,public   :: nbmodsimultrait    => null()
     ! Simulation H1 - QTL description
     !-------------------
     real (kind=dp), dimension (:), pointer,   public  :: posiqtl           => null()
     character(len=LEN_S),  dimension (:), pointer,   public  :: chrqtl     => null()
     real (kind=dp), dimension (:), pointer   ,public  :: xlim              => null()

    !!****v* m_qtlmap_genealogy/genealogy_outbred_gen
    !!  NAME
    !!   genealogy_outbred_gen
    !!  DESCRIPTION
    !!   Indicates the outbred generation in simulation case.
    !!***
    logical                             ,private      :: genealogy_outbred_gen = .false.

     contains

      procedure ,public :: copy     => copy_simulation_info
      procedure ,public :: release  => release_simulation_info
      procedure ,public :: read_simulation_file
      procedure , public :: sim_genea

      procedure , private :: set_simulation_marker
      procedure , private :: read_correlation_matrix
      procedure , private :: read_attributes_discrete_trait
      procedure , private :: set_soglia
      procedure , private :: set_simulation_genealogy
      procedure , private :: read_qtl_effect_on_trait
      procedure , private :: sim_QTL
      procedure , private :: sim_transform_discret
      procedure , private :: init_dg
      procedure , private :: sim_genea_outbread
      procedure , private :: sim_genea_F2_BC
      procedure , private :: set_simulation_traits
      procedure , private :: set_simulation_qtl
      procedure , private :: set_simulation_simultraits

      procedure ,public :: init_simul_marker
      procedure ,public :: init_genotype_simul
      procedure ,public :: init_traits_simul
      procedure ,public :: init_permutation
      procedure ,public :: log_infoqtldefinedbyuser

      procedure ,public :: manage_simulator_traits

   end type SIMULATION_INFO
!!***

   public :: write_perf_qtl
   public :: create_simulation_file
   public :: help_paramsim


  contains

    subroutine copy_simulation_info(simulInfo,copySimulInfo)
      class(SIMULATION_INFO)            ,intent(in)    :: simulInfo
      type(SIMULATION_INFO)             ,intent(inout) :: copySimulInfo


    end subroutine copy_simulation_info

    subroutine release_simulation_info(simulInfo)
      class(SIMULATION_INFO)            ,intent(in)    :: simulInfo


    end subroutine release_simulation_info

!!****f* m_qtlmap_simulation/read_simulation_file
!!  NAME
!!    read_simulation_file
!!  DESCRIPTION
!!   read simulation user file and initialize internal array
!!  INPUT
!!    dataset        :
!!  OUTPUTS
!!    dataset        : the dataset simulated
!!    simul_info     : information about simulation read in simulation parameter file
!!  NOTES
!!  SOURCE
      subroutine read_simulation_file(simul_info,dataset)
         class(SIMULATION_INFO) ,intent(out)       :: simul_info
         type(QTLMAP_DATASET) ,intent(inout)       :: dataset

         !local
         character(len=LEN_BUFFER_WORD)            :: token
         integer                                   :: ios
         logical                                   :: label_traits_present,label_simultraits_present

         call log_mess('SUBROUTINE : read_simulation_file',DEBUG_DEF)
         call log_mess('reading simulation file...',INFO_DEF)
         !dataset%files%in_parsim = parsim_file
         simul_info%simulQtl   = .false.
         simul_info%simulMap   = .false.
         simul_info%nqtlsimul  = 0

         label_traits_present = .false.
         label_simultraits_present = .false.

         call file_exist(dataset%params%get_file_val(K_PARAMSIM))

         open(unit=ios_simul_file,file=dataset%params%get_file_val(K_PARAMSIM),action="read",form="formatted")
         ios = 0
         do while ( ios == 0 )
            current_line = current_line + 1
            read(ios_simul_file,*,iostat=ios) token
            if ( ios == 0 ) then
               select case (trim(token))
                 case (LABEL_MARKERS)
                    call log_mess('setting simulation of markers...',INFO_DEF)
                    call simul_info%set_simulation_marker(dataset)
                 case (LABEL_GENEALOGY)
                    call log_mess('setting simulation of genealogy...',INFO_DEF)
                    call simul_info%set_simulation_genealogy(dataset)
                    call simul_info%sim_genea(dataset)
                 case (LABEL_QTL)
                    call log_mess('setting simulation of qtl...',INFO_DEF)
                    call simul_info%set_simulation_qtl(dataset)
                 case (LABEL_TRAITS) ! trait are described in the model and unknown data are taken
                    if ( dataset%params%get_file_val(K_TRAITS) == "") then
                       call stop_application("None traits file are defined but a LABEL "//trim(LABEL_TRAITS)&
                       //" is detected. Please use a LABEL "//trim(LABEL_SIMULTRAITS)//" to not use real traits.")
                    end if
                    if ( label_simultraits_present ) then
                         call stop_application(trim(LABEL_SIMULTRAITS)//" have been defined. You cannot define "//&
                         trim(LABEL_TRAITS)//"." )
                    end if
                    call log_mess('setting simulation of traits...',INFO_DEF)
                    if (simul_info%simulGenea) then
                      call stop_application("You cannot use real data traits with a simulated genealogy.")
                    end if
                    call log_mess('setting simul with real traits...',INFO_DEF)
                    simul_info%croisement=F2_KEYWORD

                    if ( dataset%params%get_file_val(K_GENEA) == "" ) then
                      call stop_application("Genealogy file is not defined and none genealogy are simulated.")
                    end if

                    call read_genealogy(dataset)
                    call read_model(dataset)
                    call simul_info%set_simulation_traits(dataset) !! initilialize dpm%filter_car private integer array
                    call fixe_structure_model(dataset)
                    call read_traits(dataset)
                    label_traits_present=.true.

                  case (LABEL_SIMULTRAITS) ! traits are not real (not defined in the model file)
                       if ( label_traits_present ) then
                         call stop_application(LABEL_TRAITS//" have been defined. You cannot define "//LABEL_SIMULTRAITS//"." )
                       end if
                       if (.not. simul_info%simulGenea) then
                          if ( dataset%params%get_file_val(K_GENEA) == "" ) then
                             call stop_application("Genealogy file is not defined and none genealogy are simulated.")
                          end if
                          simul_info%croisement=F2_KEYWORD
                          call read_genealogy(dataset)
                       end if
                       call log_mess('setting simul traits...',INFO_DEF)
                       call simul_info%set_simulation_simultraits(dataset)
                       label_simultraits_present=.true.
                 case default
                    call stop_on_error(-1,dataset%params%get_file_val(K_PARAMSIM), current_line,"label unknown:"//trim(token))
               end select
            end if
         end do

         ! checking at least simulation traits
         if ( .not. simul_info%simulTraits ) then
             call stop_application('Check simulation file. '//trim(LABEL_SIMULTRAITS)//' Label is not defined.')
         end if

         if ( .not. simul_info%simulMap .and. dataset%params%get_file_val(K_GENOTYPE) == "" ) then
           call stop_application("Genotype file is not defined and none genotype are simulated.")
         end if
         if ( .not. simul_info%simulMap .and. dataset%params%get_file_val(K_MAP) == "") then
           call stop_application("Map file is not defined and none map are simulated.")
         end if
        !
         close(ios_simul_file)

         call log_mess('END SUBROUTINE : read_simulation_file',DEBUG_DEF)
      end subroutine read_simulation_file
!!***

!!****f* m_qtlmap_simulation/init_permutation
!!  NAME
!!    init_permutation
!!  DESCRIPTION
!!   initialise a permutation (none simulation user file is readed)
!!  INPUTS
!!   is_transcriptome : data are transcriptomic ?
!!  NOTES
!!  SOURCE
      subroutine init_permutation(simul_info,dataset)
         class(SIMULATION_INFO)      ,intent(inout)           :: simul_info
         type(QTLMAP_DATASET)       ,intent(inout)            :: dataset


          simul_info%simulMap=.false.
          simul_info%simulGenea=.false.
          call read_genealogy(dataset)
          call read_model(dataset)
          call simul_info%set_simulation_traits(dataset) !! initilialize dpm%filter_car private integer array
          call fixe_structure_model(dataset)
          call read_traits(dataset)

      end subroutine init_permutation
!!***

!!****f* m_qtlmap_simulation/set_simulation_marker
!!  NAME
!!    set_simulation_marker
!!  DESCRIPTION
!!   read the MARKER part
!!  OUTPUTS
!!   simMap : simulation map ?
!!  NOTES
!!  SOURCE
      subroutine set_simulation_marker(simul_info,dataset)
       class(SIMULATION_INFO) , intent(inout)   :: simul_info
       type(QTLMAP_DATASET)  , intent(inout)   :: dataset

       !local
       character(len=LEN_DEF)        :: message
       integer                       :: ios,nb_marker_used,alloc_stat,iq,ic
       real                          :: temp

       simul_info%simulMap = .true.
       current_line = current_line + 1
       read(ios_simul_file,*,iostat=ios) simul_info%dens,&
                                         simul_info%nalle,&
                                         simul_info%taille

       if ( simul_info%dens <= 0.d0 ) then
          call stop_application("Density of marker can not be equal or less than 0.")
       end if

       if ( simul_info%taille <= 0.d0 ) then
          call stop_application("Size of map can not be equal or less than 0.")
       end if

       if ( simul_info%nalle <= 0.d0 ) then
          call stop_application("Number allele per marker can not be equal or less than 0.")
       end if

       if ( simul_info%dens > simul_info%taille ) then
           call stop_application("Density can not be greater than the size map.")
       end if


       message = "Problem detected at the markers description:" // &
       " Allowed: density(float) allele number(int) chromosome size(float)" // &
       " unity(int) unknwon char(char)"

       call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,trim(message))
       ! size marker...
       !allocate (nmk(1),stat=ios)
       !call check_allocate(ios,'nmk (m_qtlmap_simulation)')
        !print *,'HOLA'
       temp = simul_info%taille
       dataset%map%nmk(1) = 0
       do while ( temp >= 0 )
          dataset%map%nmk(1) = dataset%map%nmk(1)+1
          temp = temp - simul_info%dens
       end do

       nb_marker_used = dataset%map%nmk(1)

       do iq=1,simul_info%nqtlsimul
         do ic=1,dataset%phenoModel%ncar
            call log_mess('Trait :'//trim(str(ic))//' Qtl:'//trim(str(iq))//&
             ' val:'//trim(str(simul_info%ue(ic,iq)))//' M')
         end do                              !! fin ic
       end do                              !! fin iq

      end subroutine set_simulation_marker
!!***

!!****f* m_qtlmap_simulation/init_simul_marker
!!  NAME
!!    init_simul_marker
!!  DESCRIPTION
!!   initialize the array of marker coming from m_qtlmap_data
!!  NOTES
!!  SOURCE
      subroutine init_simul_marker(simul_info,dataset)
        class(SIMULATION_INFO) , intent(inout)   :: simul_info
        type(QTLMAP_DATASET)   , intent(inout)   :: dataset
        integer               :: alloc_stat,nballelepossible,i,maxnmk
        type(GENOTYPE_BASE)  , pointer :: dga

        dga => dataset%genoAnimal
        maxnmk = maxval(dataset%map%nmk)

        allocate (dataset%map%mark(size(dataset%map%nmk,1),maxnmk), stat = alloc_stat)
	    call check_allocate(alloc_stat,'mark')

        do i=1,dataset%map%nmk(1)
          dataset%map%mark(1,i)='gen_mark_'//str(i)
        end do

	    allocate (dataset%map%posi(size(dataset%map%nmk,1),maxnmk), stat = alloc_stat)
	    call check_allocate(alloc_stat,'posi')

	    allocate (dataset%map%posif(size(dataset%map%nmk,1),maxnmk), stat = alloc_stat)
	    call check_allocate(alloc_stat,'posif')

	    allocate (dataset%map%posim(size(dataset%map%nmk,1),maxnmk), stat = alloc_stat)
	    call check_allocate(alloc_stat,'posim')

        allocate (dataset%map%rm(size(dataset%map%nmk,1),maxnmk,maxnmk), stat = alloc_stat)
        call check_allocate(alloc_stat,'rm')

        allocate (dataset%map%rf(size(dataset%map%nmk,1),maxnmk,maxnmk), stat = alloc_stat)
        call check_allocate(alloc_stat,'rf')

        allocate ( dga%nall(size(dataset%map%nmk,1),maxnmk) , stat = alloc_stat)
        call check_allocate(alloc_stat,'nall')

        nballelepossible = int(VAL_MAX_INDEX_PHENO) - int(VAL_MIN_INDEX_PHENO) + 1

        allocate ( dga%alleles(size(dataset%map%nmk,1),maxnmk,nballelepossible) , stat = alloc_stat)
        call check_allocate(alloc_stat,'alleles')

        allocate ( dga%pc_all(size(dataset%map%nmk,1),maxnmk,nballelepossible) , stat = alloc_stat)
        call check_allocate(alloc_stat,'pc_all')

      end subroutine init_simul_marker
!!***

!!****f* m_qtlmap_simulation/set_simulation_genealogy
!!  NAME
!!    set_simulation_genealogy
!!  DESCRIPTION
!!   read the GENEALOGY part
!!
!!  OUTPUTS
!!   croisement : OUTBRED_KEYWORD, F2_KEYWORD, BC_KEYWORD
!!  NOTES
!!  SOURCE
      subroutine set_simulation_genealogy(simul_info,dataset)
        class(SIMULATION_INFO)      ,intent(inout)            :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        character(len=LEN_DEF)      :: message
        integer                   :: ios
	    type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        simul_info%simulGenea = .true.
        current_line = current_line + 1
        read(ios_simul_file,*,iostat=ios) simul_info%croisement
        message = "None genealogy description is finding."
        call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,trim(message))

        if ( simul_info%croisement /= F2_KEYWORD .and. &
             simul_info%croisement /= BC_KEYWORD .and. &
             simul_info%croisement /= OUTBRED_KEYWORD ) then
           message = dataset%params%get_file_val(K_PARAMSIM) // "-Problem detected at the genealogy description:"//&
       ". SIMULATION FOR NEW PEDIGREE IS ONLY AVAILABLE FOR F2 TYPE CROSSES : value ["//trim(simul_info%croisement)//"]"
           call stop_application(trim(message))
        endif


        current_line = current_line + 1
        read(ios_simul_file,*,iostat=ios) simul_info%nbpere, simul_info%inmp, simul_info%indm
        call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
         current_line,"number of sires , number of dams per sire" // &
        ", number of progeny per dam")

      end subroutine set_simulation_genealogy
!!***



!!****f* m_qtlmap_simulation/init_genotype_simul
!!  NAME
!!    init_genotype_simul
!!  DESCRIPTION
!!    initialize the array of marker coming from m_qtlmap_data
!!  NOTES
!!  SOURCE
    subroutine  init_genotype_simul(simul_info,dataset)
         class(SIMULATION_INFO)      ,intent(in)      :: simul_info
         type(QTLMAP_DATASET)    ,intent(inout)       :: dataset
         integer                    :: ios
         type(GENEALOGY_BASE) , pointer :: dg
         type(GENOTYPE_BASE) , pointer :: dga

         dga => dataset%genoAnimal
         dg => dataset%genea

         dga%nmes = dg%ngp+dg%ngm+dg%np+dg%nm+dg%nd
         allocate ( dga%numero(dga%nmes) , stat=ios)
         call check_allocate(ios,'numero')

         allocate ( dga%pheno(size(dataset%map%nmk),maxval(dataset%map%nmk),dga%nmes,2) , stat=ios)
         call check_allocate(ios,'pheno')

         allocate (dga%corregp(dg%ngp), stat=ios)
         call check_allocate(ios,'corregp')

         allocate (dga%corregm(dg%ngm), stat=ios)
         call check_allocate(ios,'corregm')

         allocate (dga%correr(dg%nr), stat=ios)
         call check_allocate(ios,'correr')

         allocate (dga%correp(dg%np), stat=ios)
         call check_allocate(ios,'correp')

         allocate (dga%correm(dg%nm), stat=ios)
         call check_allocate(ios,'correm')

         allocate (dga%corred(dg%nd), stat=ios)
         call check_allocate(ios,'corred')

         allocate (dga%presentg(dataset%map%nchr,dg%nd), stat=ios)
         call check_allocate(ios,'presentg')

    end subroutine init_genotype_simul
!!***

!!****f* m_qtlmap_simulation/init_traits_simul
!!  NAME
!!    init_traits_simul
!!  DESCRIPTION
!!    initialize the array of marker coming from m_qtlmap_data
!!  NOTES
!!  SOURCE
    subroutine  init_traits_simul(simulInfo,dataset)
         class(SIMULATION_INFO)  , intent(in)         :: simulInfo
         type(QTLMAP_DATASET)    ,intent(inout)       :: dataset
         integer                    :: alloc_stat,ic
         type(GENEALOGY_BASE) , pointer :: dg
         type(PHENOTYPE_BASE) , pointer :: dpa
         type(DATAMODEL_BASE) , pointer :: dpm

         dpm => dataset%phenoModel
         dg => dataset%genea
         dpa => dataset%phenoAnimal

         allocate (dpm%carac(dpm%ncar) , stat=alloc_stat)
         call check_allocate(alloc_stat,'carac')

         allocate (dpa%corperf(dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'corperf')

         allocate (dpa%y(dpm%ncar,dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'y')

         allocate (dpa%presentc(dpm%ncar,dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'presentc')

         allocate(dpm%natureY(dpm%ncar), stat = alloc_stat)
         call check_allocate(alloc_stat,'natureY')

         allocate (dpa%cd(dpm%ncar,dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'cd')
         dpa%cd=1.d0
         dpa%presentc = .true.

         allocate (dpa%ndelta(dpm%ncar,dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'ndelta')

         allocate (dpa%covar(dg%nd,dpm%ncov), stat = alloc_stat)
         call check_allocate(alloc_stat,'covar')

         allocate (dpm%namefix(dpm%nfix), stat = alloc_stat)
         call check_allocate(alloc_stat,'namefix')

         allocate (dpm%namecov(dpm%ncov), stat = alloc_stat)
         call check_allocate(alloc_stat,'namecov')

         allocate (dpm%listelev(dpm%nfix,dg%nd), stat = alloc_stat)
         call check_allocate(alloc_stat,'listelev')

         allocate (dpm%modele(dpm%ncar,3+(2*dpm%nfix)+dpm%ncov), stat = alloc_stat)
         call check_allocate(alloc_stat,'modele')
         dpm%modele=0

         allocate (dpm%listModelTrait(dpm%ncar))
         do ic=1,dpm%ncar
            allocate ( dpm%listModelTrait(ic)%indexFixedEffect(0))
            dpm%listModelTrait(ic)%nbfe=0
            allocate ( dpm%listModelTrait(ic)%indexCovariate(0))
            dpm%listModelTrait(ic)%nbco=0
            allocate (dpm%listModelTrait(ic)%nbint(MAX_QTL))
            dpm%listModelTrait(ic)%nbint=0
         end do

         allocate (dpm%nlev(dpm%nfix), stat = alloc_stat)
         call check_allocate(alloc_stat,'nlev')

         allocate (dpa%nivx(dg%nd,dpm%nfix), stat = alloc_stat)
         call check_allocate(alloc_stat,'nivx')

    end subroutine init_traits_simul
!!***

!!****f* m_qtlmap_simulation/set_simulation_qtl
!!  NAME
!!    set_simulation_qtl
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
      subroutine set_simulation_qtl(simul_info,dataset)
          class(SIMULATION_INFO)      ,intent(inout)       :: simul_info
          type(QTLMAP_DATASET)       ,intent(inout)        :: dataset

          integer                   :: ios,i,iq,jq,chr
          real :: l
          character(len=LEN_DEF)                 :: word
          character(len=LEN_LINE)                :: line_read
          logical                                :: is_ok

          simul_info%simtyp   = .true.

          current_line = current_line + 1
          read(ios_simul_file,*,iostat=ios) simul_info%nqtlsimul

           call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,"QTL number")

           if (simul_info%nqtlsimul.gt.0) then
               simul_info%simulQtl=.true. ! More than one qtl are simulated
           end if

           allocate(simul_info%posiqtl(simul_info%nqtlsimul))
           allocate(simul_info%chrqtl(simul_info%nqtlsimul))
           allocate (simul_info%xlim(simul_info%nqtlsimul))


           !*** POSITION  ***
           !*****************
           current_line = current_line + 1
           call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
           call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,&
            "QTL : waiting for keyword 'position' and values' ")

           word = trim(next_word(line_read,is_ok))
           if (word /= 'position') then
             call stop_application("'position' keyword expected! **"//trim(word)//"**")
           end if

           do i=1,simul_info%nqtlsimul
             word = trim(next_word(line_read,is_ok))
             if (.not. is_ok ) then
               call stop_application("qtl "// trim(str(i))//" position  expected!")
             end if
             simul_info%posiqtl(i)=get_real(word,is_ok)
             if ( .not. is_ok ) then
                call stop_application("qtl "// trim(str(i))//" position  expected! **"//trim(word)//"**")
             end if
           end do

           !*** CHROMOSOME  ***
           !*******************
           current_line = current_line + 1
           call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
           call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
            current_line,"QTL : waiting for keyword 'chromosome' and values' ")

           word = trim(next_word(line_read,is_ok))
           if (word /= 'chromosome') then
             call stop_application("'chromosome' keyword expected! **"//trim(word)//"**")
           end if

           do i=1,simul_info%nqtlsimul
             simul_info%chrqtl(i) = trim(next_word(line_read,is_ok))
             if (.not. is_ok ) then
               call stop_application("qtl "// trim(str(i))//" chromosome  expected!")
             end if
           end do

           !*** FREQUENCIES  ***
           !********************
           current_line = current_line + 1
           call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
           call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,&
            "QTL : waiting for keyword 'frequency' and values' ")

           word = trim(next_word(line_read,is_ok))
           if (word /= 'frequency') then
             call stop_application("'frequency' keyword expected! **"//trim(word)//"**")
           end if

           do i=1,simul_info%nqtlsimul
             word = trim(next_word(line_read,is_ok))
             if (.not. is_ok ) then
               call stop_application("qtl "// trim(str(i))//" frequency  expected!")
             end if
             simul_info%xlim(i)=get_real(word,is_ok)
             if ( .not. is_ok ) then
                call stop_application("qtl "// trim(str(i))//" frequency  expected! **"//trim(word)//"**")
             end if
           end do


           ! check : * position have to be ordered and not at the same position
           do iq=1,simul_info%nqtlsimul-1
              do jq=iq+1,simul_info%nqtlsimul
               if (simul_info%posiqtl(iq) == simul_info%posiqtl(jq) .and. &
                simul_info%chrqtl(iq) == simul_info%chrqtl(jq) ) then
                call stop_application("QTLs"//trim(str(iq))//" and "//trim(str(jq))//"are simulated on the same position")
               end if

               if (simul_info%posiqtl(iq) > simul_info%posiqtl(jq) .and. &
                simul_info%chrqtl(iq) == simul_info%chrqtl(jq) ) then
                call stop_application("QTLs have to be ordered : qtl "//trim(str(iq))//"> qtl "//trim(str(jq)))
               end if
             end do
           end do

           call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
            current_line,"QTL frequencies")

      end subroutine set_simulation_qtl
!!***

!!****f* m_qtlmap_simulation/set_simulation_traits
!!  NAME
!!    set_simulation_traits
!!  DESCRIPTION
!!    set the array dpm%filter_car. This filter contains index trait referenced in the model.
!!  NOTES
!!  SOURCE
      subroutine set_simulation_traits(simul_info,dataset)
        class(SIMULATION_INFO)      ,intent(inout)            :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        integer                                :: ic,jc,iq,j
        integer                                :: ios
        character(len=LEN_DEF)                 :: word
        character(len=LEN_LINE)                :: line_read
        logical                                :: is_ok
        integer                                :: nc
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        simul_info%simulTraits = .true.
        current_line = current_line + 1
        read (ios_simul_file,*,iostat=ios) nc

         if ( nc > dpm%ncar ) then
             call log_mess("               Number of trait in the model file:"//trim(str(dpm%ncar)),ERROR_DEF)
             call log_mess("Number of trait in the simulation parameter file:"//trim(str(nc)),ERROR_DEF)
             call stop_application("Your model file do not match with the simulation parameter file."//&
             " Too many traits defined in the parameter file.")
         end if

        call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
         current_line,"Number of traits")

        ! allocation of filter
        allocate(dpm%filter_car(nc))
        dpm%filter_car=0

        if ( (.not. allocated(carac_t)) .or. size(carac_t) <= 0 ) then
           call stop_application("DEVEL ERROR- set_simulation_traits :: carac_t is empty!")
        end if
        ic=1
!       FORMAT: TRAIT_NAME H2
        do while (ic<=nc)
          current_line = current_line + 1
          call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
          if ( (ios /=0) ) then
             call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
              current_line,'traits ['//trim(str(ic))//'] is not defined.')
          end if
          if (trim(line_read)=='' ) cycle ! empty line

          word = trim(next_word(line_read,is_ok))
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
           current_line,'traits ['//&
           trim(str(ic))//'] is not defined.')

          ! carac_t is a buffer array read from the model file
          do j=1,size(carac_t)
            ! first case : real name, second case : name generated by the keyword all
            if ( carac_t(j) == word .or. (LABEL_NAME_TRAIT//trim(str(j)) == carac_t(j)) ) exit
          end do

          if ( j>size(carac_t) ) then
             call stop_application("Trait["//trim(word)//"] is not defined in the model!")
          end if
          dpm%filter_car(ic)=j

          if (natureY_t(ic) == 'i') call simul_info%read_attributes_discrete_trait(dataset,ic,line_read)
          ic=ic+1
        end do

        call simul_info%set_soglia(nc,natureY_t)

        if (nc == 1) call log_mess('Only one trait --> no correlation matrix to define',VERBOSE_DEF)

        !call simul_info%read_correlation_matrix(dataset)

        call simul_info%read_qtl_effect_on_trait(dataset,nc)

      end subroutine set_simulation_traits
!!***

!!****f* m_qtlmap_simulation/read_correlation_matrix
!!  NAME
!!    read_correlation_matrix
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
     subroutine read_correlation_matrix(simul_info,dataset)
          class(SIMULATION_INFO)      ,intent(inout)          :: simul_info
          type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

          integer                                :: nc
          character(len=LEN_DEF)                 :: word
          character(len=LEN_LINE)                :: line_read
          character(len=1000)                    :: buf_line
          logical                                :: is_ok
          integer                                :: j,ic,ios,jc
          type(DATAMODEL_BASE) , pointer :: dpm

          dpm => dataset%phenoModel

          if ( dpm%ncar <= 1 ) return
          !correlation matrix definition
          line_read=''

          do while ( trim(line_read)=='')
            current_line = current_line + 1
            call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
            if ( (ios /=0) ) then
              !read(ios_simul_file,*,iostat=ios) (rho(ic,jc),jc=1,ic-1)
              call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
               current_line,"The correlation matrix is not defined")
            end if
          end do

          !FORMAT [ [corr(1)] [ corr(2,1) corr(2,2)] ... ]

          !!keyword correlation
          word = trim(next_word(line_read,is_ok))
          if (word /= 'correlation') then
             call stop_application("'correlation' keyword expected! **"//trim(word)//"**")
          end if
          !! first "["
          buf_line = trim(line_read)

          j=index(buf_line,"[")
          if ( j == 0 ) then
              call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), &
               current_line,"Expecting '[' :**"//trim(line_read)//"**")
          end if
          buf_line=buf_line(j+1:)

          do ic=2,dpm%ncar
            call parse_real_array(buf_line,ic-1,dpm%RhoP(ic,:),is_ok)
            if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
             current_line,'Correlation matrix definition for ['//&
             trim(str(ic))//'] unknown.')
            do jc=1,ic-1

             !Add test : a correlation is between -1 and 1
             if ( dpm%RhoP(ic,jc) > 1 .or. dpm%RhoP(ic,jc) < -1 ) then
                call stop_application("Correlation have to be defined betwwen -1 and 1 bad value :["//&
                 trim(str(dpm%RhoP(ic,jc)))//"]")
             end if
             call log_mess ('Correlation matrix definition ['//trim(str(ic))//']:'//trim(str(dpm%RhoP(ic,jc))),VERBOSE_DEF)
            end do
          end do

          j=index(buf_line,"]")
          if ( j == 0 ) then
              call stop_on_error(1,dataset%params%get_file_val(K_PARAMSIM),&
                current_line,"Expecting ']' to close the matrix definition")
          end if
         buf_line=buf_line(:j-1)
         if ( trim(buf_line) /= '' ) then
           call stop_on_error(1,dataset%params%get_file_val(K_PARAMSIM),&
             current_line,"Too many information in correlation matrix :"//&
            trim(buf_line)//']')
         end if

          !fill the matrix
          do ic=2,dpm%ncar
           do jc=1,ic-1
            dpm%RhoP(jc,ic)=dpm%RhoP(ic,jc)
           end do
          end do

     end subroutine read_correlation_matrix
!!***

!!****f* m_qtlmap_simulation/read_qtl_effect_on_trait
!!  NAME
!!    read_qtl_effect_on_trait
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
     subroutine read_qtl_effect_on_trait(simul_info,dataset,nc)
        class(SIMULATION_INFO)      ,intent(inout)           :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        !local
        integer, intent(in)    :: nc
        character(len=LEN_DEF)                 :: word
        character(len=LEN_LINE)                :: line_read
        logical                                :: is_ok
        integer                                :: ic,iq,ios

        allocate (simul_info%ue(nc,simul_info%nqtlsimul))
        simul_info%ue=0.d0
        if (simul_info%nqtlsimul > 0) then
           !qtl effect

           current_line = current_line + 1
           line_read=''
           do while ( trim(line_read) == '' )
             call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
             if ( (ios /= 0) ) then
              call stop_on_error(1,dataset%params%get_file_val(K_PARAMSIM),&
                current_line,"QTL effects on traits")
             end if
           end do
           word = trim(next_word(line_read,is_ok))
           if (word /= 'qtleffect') then
             call stop_application("'qtleffect' keyword expected! **"//trim(word)//"**")
           end if


           do ic=1,nc
             do iq=1,simul_info%nqtlsimul
               word = next_word(line_read,is_ok)
               if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
                current_line,&
               'QTL ['//trim(str(iq))//'] effects on traits ['//trim(str(ic))//'] unknown.')
               simul_info%ue(ic,iq) = get_real(word,is_ok)
               if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
                current_line,&
                   'QTL ['//trim(str(iq))//'] effects on traits ['//trim(str(ic))//'] have to be a real.')
               call log_mess ( 'QTL ['//trim(str(iq))//'] effects on traits ['//trim(str(ic))//']:'&
               //trim(str(simul_info%ue(ic,iq))),VERBOSE_DEF)
             end do
           end do

        else
        call log_mess('NQTL=0 --> No qtl effect to define',VERBOSE_DEF)
       end if
     end subroutine read_qtl_effect_on_trait
!!***

!!****f* m_qtlmap_simulation/set_soglia
!!  NAME
!!    set_soglia
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
     subroutine set_soglia(simul_info,nc,nature)
       class(SIMULATION_INFO)     ,intent(inout)  :: simul_info
       integer                   ,intent(in)     :: nc
       character(len=*),dimension(nc),intent(in) :: nature
       double precision         :: ppf,freqtot
       integer   :: ic,j
    
       if ( .not. associated(simul_info%freqmodsimultrait) ) return

       allocate (simul_info%soglia(nc,size(simul_info%freqmodsimultrait,2)))

       do ic=1,nc

        if ( nature(ic) /='i') cycle

        freqtot=0.d0
        do j=1,simul_info%nbmodsimultrait(ic)-1
          freqtot=freqtot+simul_info%freqmodsimultrait(ic,j)
          call NORPPF(freqtot,PPF)
          simul_info%soglia(ic,j)=PPF
        end do
       end do

     end  subroutine set_soglia
!!***
!**********************************************************************

           SUBROUTINE NORPPF(P,PPF)
!
!     PURPOSE--THIS SUBROUTINE COMPUTES THE PERCENT POINT
!              FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
!              DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1.
!              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
!              THE PROBABILITY DENSITY FUNCTION
!              F(X) = (1/SQRT(2*PI))*EXP(-X*X/2).
!              NOTE THAT THE PERCENT POINT FUNCTION OF A DISTRIBUTION
!              IS IDENTICALLY THE SAME AS THE INVERSE CUMULATIVE
!              DISTRIBUTION FUNCTION OF THE DISTRIBUTION.
!     INPUT  ARGUMENTS--P      = THE SINGLE PRECISION VALUE
!                                (BETWEEN 0.0 AND 1.0)
!                                AT WHICH THE PERCENT POINT
!                                FUNCTION IS TO BE EVALUATED.
!     OUTPUT ARGUMENTS--PPF    = THE SINGLE PRECISION PERCENT
!                                POINT FUNCTION VALUE.
!     OUTPUT--THE SINGLE PRECISION PERCENT POINT
!             FUNCTION VALUE PPF.
!     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS.
!     RESTRICTIONS--P SHOULD BE BETWEEN 0.0 AND 1.0, EXCLUSIVELY.
!     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
!     FORTRAN LIBRARY SUBROUTINES NEEDED--SQRT, ALOG.
!     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
!     LANGUAGE--ANSI FORTRAN.
!     REFERENCES--ODEH AND EVANS, THE PERCENTAGE POINTS
!                 OF THE NORMAL DISTRIBUTION, ALGORTIHM 70,
!                 APPLIED STATISTICS, 1974, PAGES 96-97.
!               --EVANS, ALGORITHMS FOR MINIMAL DEGREE
!                 POLYNOMIAL AND RATIONAL APPROXIMATION,
!                 M. SC. THESIS, 1972, UNIVERSITY
!                 OF VICTORIA, B. C., CANADA.
!               --HASTINGS, APPROXIMATIONS FOR DIGITAL
!                 COMPUTERS, 1955, PAGES 113, 191, 192.
!               --NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS
!                 SERIES 55, 1964, PAGE 933, FORMULA 26.2.23.
!               --FILLIBEN, SIMPLE AND ROBUST LINEAR ESTIMATION
!                 OF THE LOCATION PARAMETER OF A SYMMETRIC
!                 DISTRIBUTION (UNPUBLISHED PH.D. DISSERTATION,
!                 PRINCETON UNIVERSITY), 1969, PAGES 21-44, 229-231.
!               --FILLIBEN, 'THE PERCENT POINT FUNCTION',
!                 (UNPUBLISHED MANUSCRIPT), 1970, PAGES 28-31.
!               --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
!                 DISTRIBUTIONS--1, 1970, PAGES 40-111.
!               --THE KELLEY STATISTICAL TABLES, 1948.
!               --OWEN, HANDBOOK OF STATISTICAL TABLES,
!                 1962, PAGES 3-16.
!               --PEARSON AND HARTLEY, BIOMETRIKA TABLES
!                 FOR STATISTICIANS, VOLUME 1, 1954,
!                 PAGES 104-113.
!     COMMENTS--THE CODING AS PRESENTED BELOW
!               IS ESSENTIALLY IDENTICAL TO THAT
!               PRESENTED BY ODEH AND EVANS
!               AS ALGORTIHM 70 OF APPLIED STATISTICS.
!               THE PRESENT AUTHOR HAS MODIFIED THE
!               ORIGINAL ODEH AND EVANS CODE WITH ONLY
!               MINOR STYLISTIC CHANGES.
!             --AS POINTED OUT BY ODEH AND EVANS
!               IN APPLIED STATISTICS,
!               THEIR ALGORITHM REPRESENTES A
!               SUBSTANTIAL IMPROVEMENT OVER THE
!               PREVIOUSLY EMPLOYED
!               HASTINGS APPROXIMATION FOR THE
!               NORMAL PERCENT POINT FUNCTION--
!               THE ACCURACY OF APPROXIMATION
!               BEING IMPROVED FROM 4.5*(10**-4)
!               TO 1.5*(10**-8).
!     WRITTEN BY--JAMES J. FILLIBEN
!                 STATISTICAL ENGINEERING LABORATORY (205.03)
!                 NATIONAL BUREAU OF STANDARDS
!                 WASHINGTON, D. C. 20234
!                 PHONE:  301-921-2315
!     ORIGINAL VERSION--JUNE      1972.
!     UPDATED         --SEPTEMBER 1975.
!     UPDATED         --NOVEMBER  1975.
!     UPDATED         --OCTOBER   1976.
!
!---------------------------------------------------------------------
!
     double precision  :: P,PPF,P0,P1,P2,P3,P4,Q0,Q1,Q2,Q3,Q4,ADEN,T,ANUM,R
     integer           :: IPR

      DATA P0,P1,P2,P3,P4   &
     /-.322232431088,-1.0,  &
      -.342242088547,-.204231210245E-1, &
      -.453642210148E-4/
      DATA Q0,Q1,Q2,Q3,Q4       &
     /.993484626060E-1,.588581570495, &
      .531103462366,.103537752850,  &
      .38560700634E-2/
!
      IPR=6
!
!     CHECK THE INPUT ARGUMENTS FOR ERRORS
!
      IF(P.LE.0.0.OR.P.GE.1.0) GOTO 50
      GOTO 90
   50 WRITE(IPR,1)
 !     WRITE(IPR,46)P
      RETURN
   90 CONTINUE
    1 FORMAT('115H***** FATAL ERROR--THE FIRST  INPUT ARGUMENT TO THE'//  &
      'NORPPF SUBROUTINE IS OUTSIDE THE ALLOWABLE (0,1) INTERVAL *****')
 !  46 FORMAT('35H***** THE VALUE OF THE ARGUMENT IS ,F15.8,6H *****')
!
!-----START POINT-----------------------------------------------------
!
      IF(P.NE.0.5)GOTO 150
      PPF=0.0
      RETURN
!
  150 R=P
      IF(P.GT.0.5)R=1.0-R
      T=SQRT(-2.0*ALOG(real(R)))
      ANUM=((((T*P4+P3)*T+P2)*T+P1)*T+P0)
      ADEN=((((T*Q4+Q3)*T+Q2)*T+Q1)*T+Q0)
      PPF=T+(ANUM/ADEN)
      IF(P.LT.0.5)PPF=-PPF
      RETURN
!
      END  SUBROUTINE

!!****f* m_qtlmap_simulation/read_attributes_discrete_trait
!!  NAME
!!    read_attributes_discrete_trait
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
      subroutine read_attributes_discrete_trait(simul_info,dataset,ic,line_read)
        class(SIMULATION_INFO)      ,intent(inout)            :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        integer              , intent(in)        :: ic
        character(len=LEN_DEF) ,intent(inout)    :: line_read
        character(len=LEN_DEF)                   :: word
        integer                                  :: nbmod,i
        logical                                  :: is_ok
        real                                     :: summod
        real (kind=dp),dimension(:,:)  , allocatable :: freqmodsimultrait_t
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        if ( .not. associated(simul_info%nbmodsimultrait)) then
           allocate(simul_info%nbmodsimultrait(dpm%ncar))
           simul_info%nbmodsimultrait=0
        end if
        !number of modality
        word = trim(next_word(line_read,is_ok))

        if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
         current_line,&
         'number of modality ['//trim(str(ic))//'] is not defined.')
        nbmod=get_int(word,is_ok)

        if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
         current_line,&
          'modality of traits ['//trim(str(ic))//'] have to be an integer.')
        simul_info%nbmodsimultrait(ic)=nbmod

        if ( .not. associated(simul_info%freqmodsimultrait)) then
            allocate(simul_info%freqmodsimultrait(dpm%ncar,nbmod))
            simul_info%freqmodsimultrait=0.d0
        else
          if (size(simul_info%freqmodsimultrait,2)<nbmod) then
            allocate(freqmodsimultrait_t(size(simul_info%freqmodsimultrait,1),size(simul_info%freqmodsimultrait,2)))
            freqmodsimultrait_t=simul_info%freqmodsimultrait
            deallocate(simul_info%freqmodsimultrait)
            allocate(simul_info%freqmodsimultrait(size(freqmodsimultrait_t,1),nbmod))
            simul_info%freqmodsimultrait=0.0
            simul_info%freqmodsimultrait(:,1:size(freqmodsimultrait_t,2))=freqmodsimultrait_t
            deallocate(freqmodsimultrait_t)
          end if
        end if

        summod=0.d0
        do i=1,nbmod
          word = trim(next_word(line_read,is_ok))
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,'freq [trait:'//&
          trim(str(ic))//' mod:"//trim(str(i))//"] is not defined.')
          simul_info%freqmodsimultrait(ic,i)=get_real(word,is_ok)
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
          'frequency of modality :'//trim(str(i))//'for trait ['//trim(str(ic))//'] have to be a real.')
          if ( (simul_info%freqmodsimultrait(ic,i) > 1.d0) .or. (simul_info%freqmodsimultrait(ic,i) < 0.d0) ) then
            call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
          'frequency of modality :'//trim(str(i))//'for trait ['//trim(str(ic))//'] have to be 0<=mod<=1.')
          end if
          summod = summod + simul_info%freqmodsimultrait(ic,i)
        end do

        if ( summod /= 1.d0 ) then
             call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
             'sum of frequency have to be equal to 1. Trait ['//trim(str(ic))//' : sum:'//trim(str(summod))//']')
        end if

      end subroutine read_attributes_discrete_trait
!!***

!!****f* m_qtlmap_simulation/set_simulation_simultraits
!!  NAME
!!    set_simulation_simultraits
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
      subroutine set_simulation_simultraits(simul_info,dataset)
        class(SIMULATION_INFO)      ,intent(inout)           :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        integer                                :: ic,jc,iq,j
        integer                                :: ios
        character(len=LEN_DEF)                 :: word
        character(len=LEN_LINE)                :: line_read
        logical                                :: is_ok
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg  => dataset%genea
        dpa => dataset%phenoAnimal

        simul_info%simulTraits = .true.
        current_line = current_line + 1
        ios=0
        word=""

        do while (ios == 0 .and. word == "" )
           read (ios_simul_file,*,iostat=ios) word
           dpm%ncar = get_int(word,is_ok)
           if ( .not. is_ok ) then
             call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
              current_line,'bad definition of number trait: ['//trim(word)//'].')
           end if
        end do

        call stop_on_error(ios,dataset%params%get_file_val(K_PARAMSIM), current_line,"Number of traits")

        call simul_info%init_traits_simul(dataset)
        allocate(dpm%h2(dpm%ncar))
        dpm%h2=0.d0
        allocate(dpm%RhoP(dpm%ncar,dpm%ncar))
        dpm%RhoP=0.d0
        allocate(dpm%RhoG(dpm%ncar,dpm%ncar))
        dpm%RhoG=0.d0

        ! FORMAT TRAITNAME NATURE HERITABILITY N=NUMBER_OF_MOD FREQ1 FREQ2..FREQN
        ! Constraint : FREQ1+FREQ2+..+FREQN=1
        do ic=1,dpm%ncar
          current_line = current_line + 1
          call GET(ios_simul_file,line_read,maxlen=LEN_LINE,IOSTAT=ios)
          if (trim(line_read)=='') cycle ! empty line

          if ( (ios /=0) ) then
             call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
              current_line,'simul traits ['//trim(str(ic))//'] is not defined.')
          end if

          word = trim(next_word(line_read,is_ok))
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
           current_line,'traits ['//trim(str(ic))//'] is not defined.')
          dpm%carac(ic)=word

          word = trim(next_word(line_read,is_ok))
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
           'nature of traits ['//trim(str(ic))//'] is not defined.')

          if ( word /= 'i' .and. word /= 'r') then
             call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
             "Availbale value for nature of traits is 'i' or 'r' **"//trim(word)//'** .')
          end if

          dpm%natureY(ic)=word

          word = trim(next_word(line_read,is_ok))
          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),&
           current_line,'Heritability of traits ['//&
          trim(dpm%carac(j))//'] is not defined.')
          dpm%h2(ic) = get_real(word,is_ok)

          if ( .not. is_ok ) call stop_on_error (1,dataset%params%get_file_val(K_PARAMSIM),current_line,&
          'Heritability of traits ['//trim(dpm%carac(j))//'] have to be a real.')
          call log_mess ('Heritability of traits ['//trim(str(ic))//']:'//trim(str(dpm%h2(ic))),VERBOSE_DEF)

          if (dpm%natureY(ic) == 'i') call simul_info%read_attributes_discrete_trait(dataset,ic,line_read)

        end do
        call simul_info%set_soglia(dpm%ncar,dpm%natureY)
         !checking heritability
        do ic=1,dpm%ncar
           if (dpm%h2(ic)>1) then
              call stop_application("Heritability for the "//trim(str(ic))//" is greater than 1")
           end if
           if (dpm%h2(ic)<0) then
              call stop_application("Heritability for the "//trim(str(ic))//" is less than 0")
           end if
        end do

        call simul_info%read_correlation_matrix(dataset)
        call simul_info%read_qtl_effect_on_trait(dataset,dpm%ncar)

      end subroutine set_simulation_simultraits
!!***

      subroutine write_perf_qtl(dataset,file_name,qtl)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        character(len=*),intent(in)     :: file_name
        integer , dimension(:,:,:)   ,intent(in)   :: qtl! (dataset%genea%np+dataset%genea%nfem+dataset%genea%nd,simul_info%nqtlsimul,2)
        integer ,dimension(dataset%phenoModel%ncar)        :: ipresent
        integer                         :: kd,ic,iph,iq
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        call log_mess("writing performance qtl in file :["//file_name//']',INFO_DEF)
        open(1,file=file_name)
        do kd=1,dg%nd
         do ic=1,dpm%ncar
           if(dpa%presentc(ic,kd)) then
              ipresent(ic)=1
           else
              ipresent(ic)=0
           end if
        end do

        write(1,*) trim(dg%animal(kd)),          &
       (dpa%y(ic,kd),ipresent(ic),'1',ic=1,dpm%ncar), &
       ((qtl(dg%np+dg%nfem+kd,iq,iph),iph=1,2),iq=1,2)
       end do
       close(1)

      end subroutine write_perf_qtl


      subroutine create_simulation_file(dataset,parsim_file,irho,ih2)
         type(QTLMAP_DATASET)       ,intent(in)            :: dataset
         character(len=*),intent(in)     :: parsim_file
         real (kind=dp), dimension (:,:) :: irho
         real (kind=dp), dimension (:)   :: ih2

         integer :: i,j

         type(DATAMODEL_BASE) , pointer :: dpm

         dpm => dataset%phenoModel

         call log_mess("creating simulation file:["//trim(parsim_file)//']',INFO_DEF)

         open(unit=ios_simul_file,file=parsim_file,action="write",form="formatted")

         write (ios_simul_file ,* ) LABEL_TRAITS
         write (ios_simul_file ,fmt="(i5)" ) dpm%ncar
         write (ios_simul_file ,fmt="("//trim(str(dpm%ncar))//"(f7.3,1x))" ) (ih2(i),i=1,dpm%ncar)

         if ( dpm%ncar > 1 ) then
           do i=2,dpm%ncar
             write (ios_simul_file ,fmt="(30(f6.3,1x))" ) (irho(i,j),j=1,i-1)
           end do
         end if

         !write (ios_simul_file ,* ) ('0.0 ',i=1,ncar)
         close (ios_simul_file)

!         do i=1,nqtl
!            do j=1,ncar
!                write (ios_simul_file ,* )
!            end do
!         end do

      end subroutine create_simulation_file

      subroutine help_paramsim

      end subroutine help_paramsim

          !!****f* m_qtlmap_genealogy/sim_genea
    !!  NAME
    !!    sim_genea
    !!  DESCRIPTION
    !!    Simulates a population according :
    !!      * number of dam by sire
    !!      * number of progenies by dam
    !!      * croisement type
    !!  INPUTS
    !!   simul_info%inmp         : number of dam by sire
    !!   simul_info%indm         : number of progenies by dam
    !!   croisement   : OUTBRED_KEYWORD, F2_KEYWORD, BC_KEYWORD
    !!
    !!  NOTES
    !!  SOURCE
    subroutine sim_genea(simul_info,dataset)
        class(SIMULATION_INFO)     ,intent(inout)            :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset


        if ( simul_info%croisement == OUTBRED_KEYWORD ) then
            !ne necessite pas de regneration si une genealogie a deja ete creer
            if ( .not. simul_info%genealogy_outbred_gen ) then
                call simul_info%sim_genea_outbread(dataset)
            end if
        else
            call simul_info%sim_genea_F2_BC(dataset)
        end if

    end subroutine sim_genea
    !!***

    !!****f* m_qtlmap_genealogy/sim_genea_outbread
    !!  NAME
    !!    sim_genea_outbread
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!   simul_info%inmp         : number of dam by sire
    !!   simul_info%indm         : number of progenies by dam
    !!
    !!  NOTES
    !!  SOURCE
    subroutine sim_genea_outbread(simul_info,dataset)
        class(SIMULATION_INFO)     ,intent(inout)            :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset


        integer   :: ind,i,j,kd
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        dg%np = simul_info%nbpere
        dg%ngp = simul_info%nbpere*simul_info%inmp + simul_info%nbpere
        dg%ngm = simul_info%nbpere*simul_info%inmp + simul_info%nbpere
        dg%nm = dg%np*simul_info%inmp
        dg%nd=  dg%nm*simul_info%indm
        dg%nr=dg%np+dg%nm

        ! les donnees ngp,ngm,np,nm,nd,nr ont ete initialisé a la lecture du fichier de parametre de simulation
        call simul_info%init_dg(dg)

        ! creation of sires and parent o them
        dg%ngmgp(1)=0
        dg%nrgm(1)=0
        ind = 1
        do i=1,dg%np
            dg%pere(i)=str(ind)
            dg%repro(i)=dg%pere(i)
            ind = ind + 1
            dg%gpere(i)=str(ind)
            ind = ind +1
            dg%gmere(i)=str(ind)
            ind = ind + 1
            dg%nrgm(i+1)=dg%nrgm(i)+1
            dg%ngmgp(i+1)=dg%ngmgp(i)+1
        end do
        ! dams
        do i=1,dg%nm
            dg%mere(i)=str(ind)
            dg%repro(dg%np+i)=dg%mere(i)
            ind = ind + 1
            dg%gpere(dg%np+i)=str(ind)
            ind = ind +1
            dg%gmere(dg%np+i)=str(ind)
            ind = ind + 1
            dg%nrgm(dg%np+i+1)=dg%nrgm(dg%np+i)+1
            dg%ngmgp(dg%np+i+1)=dg%ngmgp(dg%np+i)+1
        end do
        !progeny
        dg%nmp(1)=0
        dg%ndm(1)=0
        kd=1

        do i=1,dg%np
            dg%nmp(i+1)=dg%nmp(i)+simul_info%inmp
            do j=dg%nmp(i)+1,dg%nmp(i+1)
                dg%ndm(j+1)=dg%ndm(j)+simul_info%indm
                do kd=dg%ndm(j)+1,dg%ndm(j+1)
                    dg%animal(kd) = str(ind)
                    ind = ind+1
                end do
            end do
        end do

        call CREATE_STRUCT_DERIVED_GENEALOGY(dg)

        simul_info%genealogy_outbred_gen = .true.
     ! call log_debug_genea()
     ! stop

    end subroutine sim_genea_outbread
    !!***


    !!****f* m_qtlmap_simulation/init_simul_genealogy
    !!  NAME
    !!    init_simul_genealogy
    !!  DESCRIPTION
    !!    initialize the array of marker coming from m_qtlmap_data
    !!  NOTES
    !!  SOURCE
    subroutine init_dg(simul_info,dg)
        class(SIMULATION_INFO)      ,intent(in)   :: simul_info
        type(GENEALOGY_BASE), intent(inout)  :: dg
        integer                 :: stat

        dg%nfem = 0

        allocate (dg%gmere(dg%ngm), STAT = stat)
        call check_allocate(stat,'gmere')
        allocate (dg%gpere(dg%ngp), STAT = stat)
        call check_allocate(stat,'gpere')
        allocate (dg%repro(dg%nr), STAT = stat)
        call check_allocate(stat,'repro')
        allocate (dg%animal(dg%nd), STAT = stat)
        call check_allocate(stat,'animal')
        allocate (dg%pere(dg%np), STAT = stat)
        call check_allocate(stat,'pere')
        allocate (dg%mere(dg%nm), STAT = stat)
        call check_allocate(stat,'mere')
        allocate (dg%femelle(dg%nm), STAT = stat)
        call check_allocate(stat,'femelle')
        allocate (dg%ngmgp(dg%ngp+1), STAT = stat)
        call check_allocate(stat,'ngmgp')
        allocate (dg%nrgm(dg%ngm+1), STAT = stat)
        call check_allocate(stat,'nrgm')
        allocate (dg%nmp(dg%np+1), STAT = stat)
        call check_allocate(stat,'nmp')
        allocate (dg%ndm(dg%nm+1), STAT = stat)
        call check_allocate(stat,'ndm')
        allocate (dg%reppere(dg%np), STAT = stat)
        call check_allocate(stat,'reppere')
        allocate (dg%repmere(dg%nm), STAT = stat)
        call check_allocate(stat,'repmere')
        allocate (dg%repfem(dg%nm), STAT = stat)
        call check_allocate(stat,'repfem')
    end subroutine init_dg
    !!***


    !!****f* m_qtlmap_simulation/sim_genea_F2_BC
    !!  NAME
    !!    sim_genea_F2_BC
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!   simul_info%inmp         : number of dam by sire
    !!   simul_info%indm         : number of progenies by dam
    !!
    !!  NOTES
    !!    Sous programme de simulation de la genealogie de la population: familles de tailles equilibrees np peres, nmp meres par pere et ndm descendants par mere
    !!
    !!  SOURCE
    subroutine sim_genea_F2_BC(simul_info,dataset)
        class(SIMULATION_INFO)      ,intent(in)              :: simul_info
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset

        !Divers
        integer ind,igp,jgm,ir,ip,jm,kr,ijm,nm1,nm2,kd,i,j
        integer , dimension(:,:),allocatable :: ncr
        integer , dimension(:),allocatable   :: meres
        real                                 :: xcr
        real,external                        :: ranf
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        dg%np  = simul_info%nbpere
        dg%ngp = dg%np
        dg%ngm = dg%np
        dg%nm  = dg%np*simul_info%inmp
        dg%nd  = dg%nm*simul_info%indm
        dg%nr  = dg%np+dg%nm

        ! les donnees ngp,ngm,np,nm,nd,nr ont ete initialisé a la lecture du fichier de parametre de simulation
        call simul_info%init_dg(dg)

        !
        !******************************************************************************
        !******************************************************************************
        !    Dispositif F2 �quilibre
        !         ngp=ngm=np, 1 male et simul_info%inmp femelles par famille
        !                      nm*simul_info%indm descendants
        !******************************************************************************
        !******************************************************************************
        !
        !
        !**************************************************************************
        !        Construction des numeros d'animaux : ind
        !
        !  GRAND-PARENTS
        !
        !    - 1 � ngp => gd peres
        !    - ngp+1 a ngp+ngm => gd meres
        !**************************************************************************
        !
        ind=1
        ! Initialisation des numeros des gd peres
        !
        dg%gpere(1)=str(ind)
        dg%ngmgp(1)=0
        do igp=2,dg%ngp
            ind=ind+1
            dg%gpere(igp)=str(ind)
            !       sexe(ind)=1
            !       gener(ind)=0
            dg%ngmgp(igp)=dg%ngmgp(igp-1)+1             !! 1 par defaut
        end do
        dg%ngmgp(dg%ngp+1)=dg%ngmgp(dg%ngp)+1
        !
        ! Initialisation des numeros des gd peres
        !
        dg%nrgm(1)=0
        do jgm=1,dg%ngm
            ind=ind+1
            dg%gmere(jgm)=str(ind)
            !       sexe(ind)=2
            !       gener(ind)=0
            if(jgm.gt.1)dg%nrgm(jgm)=dg%nrgm(jgm-1)+1+simul_info%inmp    !! chaq couple F0 => 1pere+nmp meres
        end do
        dg%nrgm(dg%ngm+1)=dg%nrgm(dg%ngm)+1+simul_info%inmp
        !
        !**************************************************************************
        !  REPRODUCTEURS
        !
        !    -  ngp+ngm +1 a ngp+ngm+nr => repro
        !**************************************************************************
        ! Initialisation
        !
        allocate(meres(size(dg%mere)))
        ir=0
        ip=0
        jm=0
        dg%nmp(1)=0
        !
        do igp=1,dg%ngp
            jgm=igp                      !! 1 couple F0,
            ir=ir+1
            ip=ip+1
            ind=ind+1
            !
            ! Creation des males F1 et tables de correspondaces peres
            dg%repro(ir)=str(ind)
            dg%pere(ip)=dg%repro(ir)
            !        reppere(ip)=ir
            !        gener(ind)=1
            !        sexe(ind)=1
            if(ip.gt.1) dg%nmp(ip)=dg%nmp(ip-1)+simul_info%inmp
            ! write (1,1000) trim(repro(ir)),trim(gpere(igp)),trim(gmere(jgm)),' 1'
            !
            ! Creation des femelles F1
            do kr=1,simul_info%inmp
                ind=ind+1
                jm=jm+1
                ir=ir+1
                dg%repro(ir)=str(ind)
                meres(jm)=ind
            !          gener(ind)=1
            !          sexe(ind)=2
             ! write (1,1000) trim(repro(ir)),trim(gpere(igp)),trim(gmere(jgm)),' 1'
            end do
        end do
        !
        ! Affectations des croisements F1 (BOURRIN)
        allocate (ncr(dg%np,simul_info%inmp))
        do ijm=1,simul_info%inmp
            ip=1
            !XXX       xcr=g05caf(xcr)
            xcr=ranf()
            xcr=xcr*(dg%np+1)+dg%np*(ijm-1)
            ncr(ip,ijm)=xcr
            ! OFI: modif, sinon ncr(ip,ijm) peut valoir 0,vu que c est un index de tableaux.....
            if (ncr(ip,ijm)==0) ncr(ip,ijm) = 1
            do ip=2,dg%np
                if(ncr(ip-1,ijm).eq.dg%np*ijm)then
                    ncr(ip,ijm)=1+dg%np*(ijm-1)
                else
                    ncr(ip,ijm)=ncr(ip-1,ijm)+1
                end if
            end do
        end do

        !
        ! Tables de correspondances meres
        dg%nmp(1)=0
        do ip=1,dg%np
            if(ip.gt.1)dg%nmp(ip)=dg%nmp(ip-1)+simul_info%inmp
            do ijm=1,simul_info%inmp
                jm=dg%nmp(ip)+ijm
                dg%mere(jm)=str(meres(ncr(ip,ijm)))
                dg%repfem(jm)=jm
            !         do ir=1,nr
            !           if(mere(jm).eq.repro(ir)) repmere(jm)=ir
            !         end do
            end do
        end do



        dg%nfem=dg%nm
        dg%nmp(dg%np+1)=dg%nmp(dg%np)+simul_info%inmp
        !
        ! Descendants: genealogie et numeros
        dg%ndm(1)=0
        do ip=1,dg%np
            nm1=dg%nmp(ip)+1
            nm2=dg%nmp(ip+1)
            do jm=nm1,nm2
                if(jm.gt.1)dg%ndm(jm)=dg%ndm(jm-1)+simul_info%indm
                do kd=dg%ndm(jm)+1,dg%ndm(jm)+simul_info%indm
                    ind=ind+1
                    dg%animal(kd)=str(ind)
                !           gener(ind)=2
                !           sexe(ind)=2
                !            xcr=g05caf(xcr)
                !           if(xcr.gt.0.5)sexe(ind)=1
                !  write (1,1000) trim(animal(kd)),trim(pere(ip)),trim(mere(jm)),' 2'
                end do
            end do
        end do
        dg%ndm(dg%nm+1)=dg%ndm(dg%nm)+simul_info%indm
        !1000  format(1x,a,1x,a,1x,a,1x,a3)
        !   close(1)

        deallocate (ncr)
        deallocate(meres)
        call CREATE_STRUCT_DERIVED_GENEALOGY(dg)

    ! call log_debug_genea()
    ! stop

    end subroutine sim_genea_F2_BC
    !!***

    !!****f* m_qtlmap_traits/sim_QTL
!!  NAME
!!    sim_QTL
!!  DESCRIPTION
!!
!!
!!  NOTES
!!   Simulation des typages sur trois generations , reperage des alleles au qtl des descendants.
!!  SOURCE
      subroutine sim_QTL(simul_info,dataset,spt,icar,pas)
      class(SIMULATION_INFO)      ,intent(in)              :: simul_info
      type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
      type(PDD_BUILD)            ,intent(in)               :: spt
      integer                    ,intent(in)               :: icar,pas
      !
! Tableaux dimensionnes selon n, le nombre d'individus dans le pedigree
  !    integer ,intent(out)    :: qtl(dataset%genea%np+dataset%genea%nm+dataset%genea%nd,nqtl,2)
!
!C Divers
      real              :: x_rand
      integer :: geno,iq,ip,ifem,ipos,nm1,nm2,jm,ngeno1,ngeno2,nd1,nd2,chr
      integer :: kd,kkd,iph,n
      real (kind=dp) :: pp,pm
      real,external                        :: ranf
      logical :: ok

      integer ,dimension(:,:,:)   ,allocatable :: qtl

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm


      allocate (qtl(dataset%genea%np+dataset%genea%nm+dataset%genea%nd,simul_info%nqtlsimul,2))

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal


!
!
!***********************************************************************
!                Simulation des alleles au QTL des parents
!
!***********************************************************************
      do iq=1,simul_info%nqtlsimul
       do ip=1,dg%np
! simul_info%xlim(iq)=frequence de allele Q2 chez les GP et de Q1 chez les GM
! si simul_info%xlim(iq)=1, tous GP Q1Q1 et toutes GM Q2Q2
        qtl(ip,iq,1)=2
        qtl(ip,iq,2)=1
        x_rand=ranf()
        if (dble(x_rand).lt.simul_info%xlim(iq))  qtl(ip,iq,1)=1
        x_rand=ranf()
        if (dble(x_rand).lt.simul_info%xlim(iq))  qtl(ip,iq,2)=2
       end do
!
       do ifem=1,dg%nfem
! simul_info%xlim(iq)=frequence de allele Q2 chez les GP et de Q1 chez les GM
! si simul_info%xlim(iq)=1, tous GP Q1Q1 et toutes GM Q2Q2
        qtl(dg%np+ifem,iq,1)=2
        qtl(dg%np+ifem,iq,2)=1
        if (simul_info%croisement /= BC_KEYWORD) then
          x_rand=ranf()
          if (dble(x_rand).gt.simul_info%xlim(iq))  qtl(dg%np+ifem,iq,1)=1
          x_rand=ranf()
          if (dble(x_rand).gt.simul_info%xlim(iq))  qtl(dg%np+ifem,iq,2)=2
        else
          qtl(dg%np+ifem,iq,2)=2
        end if

       end do
      end do


!
!
!**********************************************************************
!           Simulation de la transmission des alleles QTL des parents
! aux descendants : les estimations de phase des meres sont correctes
! la phase la plus probable est systematiquement transmise
!  ==> simulation de la deviation a la moyenne nulle
!**********************************************************************
!
      do iq=1,simul_info%nqtlsimul
       chr=0
       ok=.false.
       do while ( .not. ok )
         chr = chr + 1
         if ( chr > dataset%map%nchr ) call stop_application("none chromosome "//trim(simul_info%chrqtl(iq))//" are defined !");
         ok = simul_info%chrqtl(iq) == dataset%map%chromo(chr)
       end do

       n=dataset%map%get_pos(chr,(simul_info%posiqtl(iq)-dataset%map%posi(chr,1)))
       if ( n > dataset%map%get_ilong(chr)) n = dataset%map%get_ilong(chr)
       if ( n <= 0 ) n = 1
       do ip=1,dg%np
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
         ifem=dg%repfem(jm)
         ngeno1=spt%ngenom(chr,jm)+1
         ngeno2=spt%ngenom(chr,jm+1)
         geno=ngeno1
         nd1=spt%ngend(chr,geno)+1
         nd2=spt%ngend(chr,geno+1)
         do kd=nd1,nd2
            kkd=spt%ndesc(chr,kd)
! Transmission de l'allele paternel: a priori re�oit allele 1 (grand p�re) du p�re
! si pt� transmission pp non d�pass�e, re�oit all�le 2
          qtl(dg%np+dg%nfem+kkd,iq,1)=qtl(ip,iq,1)
        !  pp=(-pdd(chr,kd,1,n)-pdd(chr,kd,2,n)+pdd(chr,kd,3,n)+pdd(chr,kd,4,n))/2.d0+0.5d0
          ! MODIF OFI :  proba qu'on est recu le 2eme allele du pere
          pp=spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
          x_rand=ranf()
          if(dble(x_rand).lt.pp)qtl(dg%np+dg%nfem+kkd,iq,1)=qtl(ip,iq,2)
!
! Transmission de l'allele maternel: a priori re�oit allele 1 (grand p�re) de la m�re
          qtl(dg%np+dg%nfem+kkd,iq,2)=qtl(dg%np+ifem,iq,1)
           ! pm=(-pdd(chr,kd,1,n)+pdd(chr,kd,2,n)-pdd(chr,kd,3,n)+pdd(chr,kd,4,n))/2.d0+0.5d0
          ! MODIF OFI :  proba qu'on est recu le 2eme allele de la mere
          pm=spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,4,n)
          x_rand=ranf()
          if(dble(x_rand).lt.pm)qtl(dg%np+dg%nfem+kkd,iq,2)=qtl(dg%np+ifem,iq,2)
!
! Mise a jour des perf simulees sous H0
! UE correspond a l effet indique par l utilisateur dans le fichier d entree de simulation : (The mean of absolute value of substitution effect)
! On eneleve ou ajoute donc la moitie de l'effet au descendants selon l allele recu au QTL
           if (dpa%presentc(icar,kkd)) then
            if (qtl(dg%np+dg%nfem+kkd,iq,1).eq.qtl(dg%np+dg%nfem+kkd,iq,2)) then
              if (qtl(dg%np+dg%nfem+kkd,iq,1).eq.1) then
                dpa%y(icar,kkd)=dpa%y(icar,kkd)+simul_info%ue(icar,iq)/2.d0
              else
                dpa%y(icar,kkd)=dpa%y(icar,kkd)-simul_info%ue(icar,iq)/2.d0
              end if
            end if
           end if
         end do
        end do
       end do
      end do

      deallocate (qtl)
      end subroutine sim_QTL
!!***

!!****f* m_qtlmap_traits/manage_simulator_traits
!!  NAME
!!    manage_simulator_traits
!!  DESCRIPTION
!!      Manage the type of simulator (depends opt_cal/nature of traits)
!!
!!      4 cases :
!!       1) permutation are switch on : -> permuted data according to the type analysis
!!                          a) MULTITRAIT ANALYSIS -> PERMUTATION BY LINE
!!                          b) UNITRAIT ANALYSIS   -> PERMUTATION on each trait
!!       2) permutation are switch off and analysis is a multitrait
!!                         --> simulation of traits with heritability and correlation matrix (according to the nature)
!!       3) permutation are switch off and analysis is unitrait
!!                         --> simulation of traits with heritability according to the nature of the trait
!!
!!  NOTES
!!   Simulation des typages sur trois generations , reperage des alleles au qtl des descendants.
!!  SOURCE
    subroutine manage_simulator_traits(simulInfo,dataset,spt,analyse_is_multi_traits,permute_mode,step)
        class(SIMULATION_INFO)      ,intent(inout)            :: simulInfo
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        type(PDD_BUILD)            ,intent(in)               :: spt

        logical ,intent(in)                          :: analyse_is_multi_traits ! generaly the simulation
        logical ,intent(in)                          :: permute_mode
        integer,intent(in)                           :: step

 !       integer      ,dimension(dataset%genea%np+dataset%genea%nfem+dataset%genea%nd,nqtl,2)  ,intent(inout) :: qtl

        real (kind=dp),dimension (1)     :: h2_t

        real (kind=dp),dimension (1,dataset%genea%nd) :: ys

        real (kind=dp),dimension (1,1)          :: bidonrho
        integer :: i,j,ic
        character(len=1) :: nat
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        if (dpm%ncar <= 0 ) then
           call stop_application("DEVEL ERROR : NCAR is not initialized")
        end if

        if ( .not. associated(dpm%natureY) ) then
            call stop_application("DEVEL ERROR : NATUREY is not initialized")
        endif
        if ( size(dpm%natureY) <= 0  ) then
            call stop_application("DEVEL ERROR : NATUREY is empty")
        end if

        !the permutation mode is active, none nimulator are called
        if ( permute_mode ) then
              call set_estime(dataset)
              if ( .not. analyse_is_multi_traits ) then
                  call log_mess("** Permutation for unitrait analysis **",VERBOSE_DEF)
                  call sim_perf_shuffling(dataset,.false.)
              else !traits line are permuted between animals
                  call log_mess("** Permutation for multitraits analysis **",VERBOSE_DEF)
                  call sim_perf_shuffling(dataset,.true.)
              end if
              return
        end if


       if ( analyse_is_multi_traits .or. all(dpm%natureY=='r')) then
          !check the nature of all traits (have to bo the same)
          nat=dpm%natureY(1)
          do i=1,size(dpm%natureY)-1
                if (dpm%natureY(i) /= dpm%natureY(i+1)) then
                  call stop_application("You have to defined a model value with the same nature of trait"//&
                  " for a multitrait analysis - [trait "//trim(str(i))//" nature:"//dpm%natureY(i)//"] [trait "//&
                  trim(str(i+1))//" nature:"//dpm%natureY(i+1)//"]")
                end if
          end do

          if ( nat == 'r' ) then
             call log_mess("** Simulation for real multitrait analysis ** ",VERBOSE_DEF)
             call sim_perf_tirage(dataset,dpm%ncar,dpm%RhoP,dpm%h2,dpa%y)
             if (simulInfo%simulQtl) then
              do ic=1,dpm%ncar
               call simulInfo%sim_QTL(dataset,spt,ic,step)
              end do
             end if
          end if

          if ( nat == 'i' ) then
              call log_mess(" ** Simulation for discrete multitrait analysis ** ",VERBOSE_DEF)
              call sim_perf_tirage(dataset,dpm%ncar,dpm%RhoP,dpm%h2,dpa%y)
              if (simulInfo%simulQtl) then
               do ic=1,dpm%ncar
                 call simulInfo%sim_QTL(dataset,spt,ic,step)
               end do
              end if
             do i=1,dpm%ncar
                call simulInfo%sim_transform_discret(dataset,i)
             end do
          end if

          if ( nat == 'c' ) then
              call stop_application("None simulator are defined for categorial data")
          end if

          if ( nat == 'a' ) then
             call log_mess("** Simulation for real multitrait analysis ** ",VERBOSE_DEF)
             call sim_perf_tirage(dataset,dpm%ncar,dpm%RhoP,dpm%h2,dpa%y)
             do i=1,dpm%ncar
              do j=1,dg%nd
                 if (dpa%cd(i,j)/=0) dpa%y(i,j)=ys(1,j)/sqrt(dpa%cd(i,j))
              end do
             end do

             if (simulInfo%simulQtl) then
              do ic=1,dpm%ncar
               call simulInfo%sim_QTL(dataset,spt,ic,step)
              end do
            end if
          end if

       else !otherwise we call for each trait the simulator corresponding to the nature of the traits
          bidonrho=0.d0
          do i=1,dpm%ncar
            h2_t(1) = dpm%h2(i)
            call log_mess("** Simulator for TRAIT ["//trim(str(i))//"]->Nature:"//trim(dpm%natureY(i)//"] **"),VERBOSE_DEF)
            select case (dpm%natureY(i))
            case ('r')  !! real / continue value
             call sim_perf_tirage(dataset,1,bidonrho,h2_t,ys)
             dpa%y(i,:)=ys(1,:)
             if (simulInfo%simulQtl) then
              call simulInfo%sim_QTL(dataset,spt,i,step)
             end if
            case ('i')  !! discrete value
             call sim_perf_tirage(dataset,1,bidonrho,h2_t,ys)
             dpa%y(i,:)=ys(1,:)
             if (simulInfo%simulQtl) then
                call simulInfo%sim_QTL(dataset,spt,i,step)
             end if
             call simulInfo%sim_transform_discret(dataset,i)
            case ('c')
              call stop_application("None simulator are defined for categorial data")
            case ('a')
              call sim_perf_tirage(dataset,1,bidonrho,h2_t,ys)

               do j=1,size(dpa%y,2)
                 if (dpa%cd(i,j)/=0) dpa%y(i,j)=ys(1,j)/sqrt(dpa%cd(i,j))
               end do

               if (simulInfo%simulQtl) then
                call simulInfo%sim_QTL(dataset,spt,i,step)
               end if

               !call stop_application("None simulator are defined for average data")
             case default
              call stop_application("Unknown nature of data ["//dpm%natureY(i)//"]")
            end select
          end do
        end if

!
        call set_estime(dataset)

    end subroutine manage_simulator_traits
!!***

    !!****f*  m_qtlmap_traits/log_infoqtldefinedbyuser
!!  NAME
!!    log_infoqtldefinedbyuser
!!  DESCRIPTION
!!   print in the console information about localisation of qtl simulated.
!!
!!  NOTES
!!
!!  SOURCE
     subroutine log_infoqtldefinedbyuser(simul_info,dataset,nqtl)
      class(SIMULATION_INFO)      ,intent(in)            :: simul_info
      type(QTLMAP_DATASET)       ,intent(inout)       :: dataset
      integer , intent(in) :: nqtl
      integer :: i,chr
      real(kind=dp) :: l
      logical :: ok

      if ( nqtl <=0 ) return
      ! LOG*****
      call log_mess("--------------- QTL DEFINED --------------------",INFO_DEF)

      do i=1,nqtl
          chr=0
          ok=.false.
          do while ( .not. ok )
            chr = chr + 1
            if ( chr > dataset%map%nchr ) call stop_application("none chromosome "//&
             trim(simul_info%chrqtl(i))//" are defined !");
            ok = simul_info%chrqtl(i) == dataset%map%chromo(chr)
          end do

          call log_mess(trim(str(i))//":Position defined by the user :"//trim(str(simul_info%posiqtl(i))),INFO_DEF)
          l=(dataset%map%get_pos(chr,&
          (simul_info%posiqtl(i)-dataset%map%posi(chr,1)))*&
          ((dataset%map%posi(chr,dataset%map%nmk(chr))-dataset%map%posi(chr,1))/dataset%map%get_npo(chr)))+dataset%map%posi(chr,1)
          call log_mess(trim(str(i))//":Position according the sampling :"//trim(str(l)),INFO_DEF)
       end do

       call log_mess("--------------------------------------------------------",INFO_DEF)
     end subroutine log_infoqtldefinedbyuser
!!***


!!****f* m_qtlmap_traits/sim_transform_discret
!!  NAME
!!    sim_transform_discret
!!  DESCRIPTION
!!
!!  NOTES
!!   dans cette subroutine, on genere les donnees discretes (a mettre dans le tableau y) a partir des valeurs de la sous jacente (trouvees dans le tableau y)
!!  SOURCE
  subroutine sim_transform_discret(simul_info,dataset,ic)
      class(SIMULATION_INFO)     ,intent(inout)            :: simul_info
      type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
      integer , intent(in)                                 :: ic

      integer     :: kd,m
      type(PHENOTYPE_BASE) , pointer :: dpa
      dpa => dataset%phenoAnimal

       do kd=1,size(dpa%y,2)
         if(dpa%presentc(ic,kd))then
            m=1
            do while (m < simul_info%nbmodsimultrait(ic) .and. dpa%y(ic,kd) > simul_info%soglia(ic,m))
             m=m+1
            enddo
            dpa%y(ic,kd)=m
         end if
      end do

  end subroutine sim_transform_discret
!!***


end module m_qtlmap_simulation
