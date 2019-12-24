!!****m* INPUT/m_qtlmap_genealogy
!!  NAME
!!    m_qtlmap_genealogy -- Genealogy routines.
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
module m_qtlmap_genealogy
    !Internal parameter
    use m_qtlmap_base
    use m_qtlmap_log
    use m_qtlmap_types


    implicit none
    save

    integer , private,                           parameter                :: GENERATION_MAX=40


    public :: read_genealogy
    public :: check_and_print_correc_genealogy
    public :: check_permutation_context
    public :: write_genea

CONTAINS

!!  DESCRIPTION
!!    Read the genealogy user file. In a computation of censure data case, the genealogy number (4th record) can be greater than 3 :
!!    new traits values (and a censured data) are computated for each progeny which have progenies.
    subroutine read_genealogy(dataset)
        type(QTLMAP_DATASET)      ,intent(inout)            :: dataset
        integer                            :: nb_max_indiv = 0
        integer                            :: ios,eof,i,alloc_stat,j,niv, k,k1,k2,ip,jm
        integer                            :: max_i
        character(len=LEN_DEF)        , dimension(:,:),allocatable   :: genea_list
        integer                       , dimension(:),allocatable     :: genea_niv
        character(len=LEN_DEF)        , dimension(:,:),allocatable   :: rac
        logical :: change
        character(len=LEN_DEF)     :: rac1,an
        type(GENEALOGY_BASE)               :: dg

        call log_mess('SUBROUTINE : read_genealogy',DEBUG_DEF)
        call log_mess('reading genealogy file...',INFO_DEF)

        allocate (genea_list(MAX_ANIMAL,4))
        allocate (genea_niv(-GENERATION_MAX:GENERATION_MAX))
        allocate (rac(MAX_ANIMAL,2))
        genea_niv=0
        rac=''
        rac1=''
        ios = 57

        ! compte max indiv and check line
        open(ios,file=dataset%params%get_file_val(K_GENEA),action="read",iostat=eof,status="old")
        if ( eof /= 0 ) then
         call stop_application("Can not find the genealogy file:"//trim(dataset%params%get_file_val(K_GENEA)))
        end if

        eof = 0
        i=1
        do while ( eof == 0 )
            read(ios,*,iostat=eof) (genea_list(i,j),j=1,4)
            if ( trim(genea_list(i,1)) /= '' .and. eof == 0 ) i=i+1
        end do

        close(ios)
        max_i=i-1
        call log_mess(trim(str(max_i))//" entries in the genealogy file");

        !on compte le nombre d'individu par generation
        do i=1,max_i
            niv = get_int(genea_list(i,4))
            genea_niv(niv) = genea_niv(niv)+1
        end do

        dataset%geneaRace%NB_RACES=1
        rac='UNKNOWN'

        if ( dataset%params%get_file_val(K_RACE) /= '' ) then
            ! lecture du fichier race
            eof=0
            open(ios,file=dataset%params%get_file_val(K_RACE))
            do while ( eof == 0 )
                read(ios,*,iostat=eof) an, rac1
                do j=1,  genea_niv(1)
                    if (trim(an)==trim(genea_list(j,2))) rac(j,1)=rac1
                    if (trim(an)==trim(genea_list(j,3))) rac(j,2)=rac1
                enddo
                i=i+1
            end do
            close(ios)
!            do j=1,  genea_niv(1)
!              print *,rac(j,1),rac(j,2)
!            end do
!            stop
            call set_info_race(dataset%geneaRace,genea_list,max_i,genea_niv,rac)
!            print *,dataset%geneaRace%NB_RACES
!            print *,dataset%geneaRace%nom_race
!            stop
          else
            allocate (dataset%geneaRace%nom_race(dataset%geneaRace%NB_RACES))
            dataset%geneaRace%nom_race='UNKNOWN'
          end if

        call create_genealogy(dataset,&
                              genea_list,&
                              max_i,&
                              genea_niv,&
                              rac)

        deallocate (genea_list)
        deallocate (genea_niv)
        deallocate(rac)

    end subroutine read_genealogy
    !!***

    !!****f* m_qtlmap_genealogy/set_info_race
    !!  NAME
    !!    set_info_race
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!
    !!
    !!  NOTES
    !!  SOURCE
    subroutine set_info_race(dgr,genea_list,max_i,genea_niv,rac)
        type(GENEALOGY_RACE)       ,intent(inout)                       :: dgr
        character(len=LEN_DEF)     , dimension(:,:)      ,intent(inout) :: genea_list
        integer                                          ,intent(in)    :: max_i
        integer                       , dimension(-GENERATION_MAX:GENERATION_MAX)     ,intent(in)    :: genea_niv
        character(len=LEN_DEF)        , dimension(:,:)   ,intent(in)    :: rac

        character(len=LEN_DEF)     :: nom_race_t(30)
        integer                    :: k,j,i,k2,k1,alloc_stat

        ! identification du nombre de races dans le fichier
        k=0
        dgr%NB_RACES=0
        nom_race_t=''
        do j=1,genea_niv(1)
            if (rac(j,1)==''.or.rac(j,2)=='') then
                print *, 'Breed origin has to be given for all or no parents. There is a missing breed origin for parent ', &
                trim(genea_list(j,1))
                stop
            endif

            do k1=1,2
             do k2=1,dgr%NB_RACES
              if ( nom_race_t(k2) == trim(rac(j,k1))) exit
             end do

             if ( k2 > dgr%NB_RACES) then !new race
               dgr%NB_RACES=dgr%NB_RACES+1
               nom_race_t(dgr%NB_RACES)=trim(rac(j,k1))
             end if
            end do
        enddo

        allocate (dgr%nom_race(dgr%NB_RACES), stat = alloc_stat)
        call check_allocate(alloc_stat,'geneaRace%nom_race')

        dgr%nom_race='UNKNOWN'
        do k=1,dgr%NB_RACES
            dgr%nom_race(k)=nom_race_t(k)
        enddo

    end subroutine set_info_race
    !!***


    !!****f* m_qtlmap_genealogy/create_genealogy
    !!  NAME
    !!    create_genealogy
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!
    !!
    !!  NOTES
    !!  SOURCE
    subroutine create_genealogy(dataset,genea_list,max_i,genea_niv,rac)
        type(QTLMAP_DATASET)  ,intent(inout)                       :: dataset
        character(len=LEN_DEF)        , dimension(:,:)   ,intent(inout) :: genea_list
        integer                                          ,intent(in)    :: max_i
        integer                       , dimension(-GENERATION_MAX:GENERATION_MAX),intent(in)    :: genea_niv
        character(len=LEN_DEF)        , dimension(:,:)   ,intent(in)    :: rac
        integer                            :: nb_max_indiv = 0
        integer                            :: alloc_stat,i,j,niv, k,k1,k2,ip,jm

        logical :: change

        type(DATASET_QTLMAP_TRIO)  ,pointer :: datasetUser
        type(GENEALOGY_BASE)       ,pointer :: dg
        type(GENEALOGY_RACE)       ,pointer :: dgr

        datasetUser => dataset%datasetUser
        dg => dataset%genea
        dgr => dataset%geneaRace

        call sort_genea_list(max_i,genea_list,change)
        !call check_stop_genealogy(max_i,genea_list)

        nb_max_indiv=0
        i=1
        dg%OldGenealogySize=0
        allocate (dg%OldGenealogy(MAX_ANIMAL,3), stat = alloc_stat)
        CALL check_allocate(alloc_stat,'OldGenealogy')
        dg%OldGenealogy='0'

        do while ( i<=max_i )
            niv = get_int(genea_list(i,4))
            nb_max_indiv = nb_max_indiv+1
            if (dataset%cli%cli_is_calcul_cd()) then
                if ( niv == 2 ) then
                    call datasetUser%add_animal_genea(genea_list(i,2),genea_list(i,3),genea_list(i,1),&
                      nb_max_indiv-sum(genea_niv(:1)))
                else
                    call datasetUser%add_animal_genea(genea_list(i,2),genea_list(i,3),genea_list(i,1))
                end if
            end if

            ! pour le modele animal on peut avoir une genealogie complete....
            dg%OldGenealogySize=dg%OldGenealogySize+1
            dg%OldGenealogy(dg%OldGenealogySize,1)=genea_list(i,1)
            dg%OldGenealogy(dg%OldGenealogySize,2)=genea_list(i,2)
            dg%OldGenealogy(dg%OldGenealogySize,3)=genea_list(i,3)

            if ( niv > GENERATION_MAX .or. niv < -GENERATION_MAX ) then
                call stop_application("Number of generation is too low or too high ["//&
                trim(genea_list(i,4))//"] ind["//trim(genea_list(i,1))//"]")
            end if

            i = i+1
        end do

        !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(2)
        !$OMP SECTIONS
        !$OMP SECTION
        CALL CREATE_STRUCT_GRAND_PARENT(dg,dgr,nb_max_indiv,genea_list,genea_niv,rac)
        !$OMP SECTION
        CALL CREATE_STRUCT_PARENT(nb_max_indiv,genea_list,dg,genea_niv)
        !$OMP END SECTIONS NOWAIT
        !$OMP END PARALLEL

        ! create repere
        ALLOCATE (dg%reppere(size(dg%pere)))
        ALLOCATE (dg%repmere(size(dg%mere)))
        ALLOCATE (dg%femelle(size(dg%mere)))
        ALLOCATE (dg%repfem(size(dg%mere)))

        dg%reppere = 0
        dg%repmere = 0
        dg%repfem  = 0

        CALL CREATE_STRUCT_DERIVED_GENEALOGY(dg)
        call log_mess('NP='//trim(str(dg%np)),VERBOSE_DEF)
        call log_mess('NM='//trim(str(dg%nm)),VERBOSE_DEF)
        call log_mess('ND='//trim(str(dg%nd)),VERBOSE_DEF)
        call log_mess('END SUBROUTINE : read_genealogy',DEBUG_DEF)

    end subroutine create_genealogy
    !!***


    !!****f* m_qtlmap_genealogy/check_and_print_correc_genealogy
    !!  NAME
    !!    check_and_print_correc_genealogy
    !!  NOTES
    !!  SOURCE
    subroutine check_and_print_correc_genealogy(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        integer                            :: nb_max_indiv = 0
        integer                            :: ios,eof,i,alloc_stat
        integer , parameter                :: GENERATION_MAX=10
        integer                            :: max_i,maxsize,j
        logical                            :: change
        character(len=LEN_DEF)        , dimension(:,:),allocatable   :: genea_list

        call log_mess('check genealogy file...',INFO_DEF)

        allocate (genea_list(MAX_ANIMAL,4))

        ios = 57
        ! compte max indiv and check line
        open(ios,file=dataset%params%get_file_val(K_GENEA))
        eof = 0
        i=1
        maxsize=0
        do while ( eof == 0 )
            read(ios,*,iostat=eof) (genea_list(i,j),j=1,4)
            if ( trim(genea_list(i,1)) /= '' .and. eof == 0 ) then
                i=i+1
            end if
        end do

        max_i=i-1
        call log_mess(trim(str(max_i))//" entries in the genealogy file",INFO_DEF);
        change=.false.
        call sort_genea_list(max_i,genea_list,change)

        close(ios)

        if ( change ) then
            call log_mess("***> Print new pedigree file : "//trim(dataset%params%get_file_val(K_GENEA))//".new",INFO_DEF);
            open(ios,file=trim(dataset%params%get_file_val(K_GENEA))//".new")

            do i=1,max_i
                write ( ios,fmt="( 3(a25,1x),1x,a5)" ) genea_list(i,1),genea_list(i,2),&
                genea_list(i,3),genea_list(i,4)
            end do
            close(ios)
        end if

        deallocate (genea_list)
        call read_genealogy(dataset)

    end subroutine check_and_print_correc_genealogy
    !!***


    subroutine check_stop_genealogy(nlist,genea_list)
        integer                       , intent(in)   :: nlist
        character(len=LEN_DEF)        , dimension(:,:),intent(inout)  :: genea_list
        integer, dimension(2) , parameter :: indexv = (/4,2/)
        integer :: v,i,j,vi
        logical :: endvalue
        character(len=LEN_DEF) :: value_check,generation,sire

        generation='Z'
        do v=1,2
            vi=indexv(v);
            do i=1,nlist
                value_check=genea_list(i,vi)
                if ( vi == 2 ) generation = genea_list(i,4)
                endvalue=.false.
                do j=i+1,nlist
                    if ( generation ==  genea_list(j,4) ) then
                        if ( value_check /= genea_list(j,vi) ) then
                            endvalue=.true.
                            cycle
                        end if
                        if ( endvalue .and. value_check == genea_list(j,vi)) then
                            call stop_application("bad definition of genealogy [values, line:"//&
                            trim(str(i))//" and line:"//trim(str(j))//" column:"//trim(str(vi))&
                            //"]. build a new genealogy file with qtlmap-check.")
                        end if
                    else
                        exit
                    end if
                end do
            end do
        end do

        do i=1,nlist
            vi=3
            value_check=genea_list(i,vi)
            sire = genea_list(i,2)
            endvalue=.false.
            do j=i+1,nlist
                if ( sire == genea_list(j,2) ) then ! intrafamille de pere pour la verif des meres
                    if ( value_check /= genea_list(j,vi) ) then
                        endvalue=.true.
                        cycle
                    end if
                    if ( endvalue .and. value_check == genea_list(j,vi)) then
                        call stop_application("bad definition of genealogy [values, line:"//&
                        trim(str(i))//" and line:"//trim(str(j))//" column:"//trim(str(vi))&
                        //"]. build a new genealogy file with qtlmap-check.")
                    end if
                else
                    exit ! on sort sinon
                end if
            end do
        end do

    end subroutine check_stop_genealogy

    !!****f* m_qtlmap_genealogy/sort_genea_list
    !!  NAME
    !!    sort_genea_list
    !!  DESCRIPTION
    !!    sort the genealogy by generation,sire,dam and progeny
    !!  INPUTS
    !!   genea_list   : the list of entry from the genealogy file
    !!
    !!  NOTES
    !!  SOURCE
    subroutine sort_genea_list(nlist,genea_list,changeSort)
        integer                       , intent(in)   :: nlist
        character(len=LEN_DEF)        , dimension(:,:),intent(inout)  :: genea_list
        logical                       ,intent(out)      :: changeSort

        character(len=LEN_DEF) , dimension(4) :: buf
        character(len=LEN_DEF) :: nextgen,nextsire
        integer :: i,j,idnextgen,idnextsire
        logical :: change
        integer :: countchange

        changeSort=.false.
        countchange=0
        ! generation
        change=.true.
        do while (change)
            change=.false.
            do i=1,nlist-1
                if ( genea_list(i,4) > genea_list(i+1,4)  ) then
                    !              if (stopOnChange) call stop_application("bad definition of genealogy [generation, line:"//&
                    !                 trim(str(i))//"]. build a new genealogy file with qtlmap-check.")
                    buf = genea_list(i,:)
                    genea_list(i,:) = genea_list(i+1,:)
                    genea_list(i+1,:) = buf
                    change = .true.
                    countchange=countchange+1
                end if
            end do
        end do

        !Sire
        idnextgen = 1
        nextgen=genea_list(idnextgen,4)

        do while (idnextgen /= 0)
            change=.true.
            do while (change)
                change=.false.
                do i=1,nlist-1
                    ! intra generation sort sire
                    if ( nextgen /= genea_list(i,4) .or. nextgen /= genea_list(i+1,4) ) cycle
                    if ( genea_list(i,2) > genea_list(i+1,2)  ) then
                        !                if (stopOnChange) call stop_application("bad definition of genealogy [sire, line:"//&
                        !                 trim(str(i))//"]. build a new genealogy file with qtlmap-check.")
                        buf = genea_list(i,:)
                        genea_list(i,:) = genea_list(i+1,:)
                        genea_list(i+1,:) = buf
                        change = .true.
                        countchange=countchange+1
                    end if
                end do
            end do

            j=idnextgen+1
            idnextgen=0
            do i=j,nlist
                if ( nextgen /= genea_list(i,4) ) then
                    idnextgen = i
                    nextgen = genea_list(idnextgen,4)
                    exit
                end if
            end do
        end do

        ! Dams
        idnextsire = 1
        nextsire=genea_list(idnextsire,2)

        do while (idnextsire /= 0)
            change=.true.
            do while (change)
                change=.false.
                do i=1,nlist-1
                    ! intra sire sort dam
                    if ( nextsire /= genea_list(i,2) .or. nextsire /= genea_list(i+1,2)) cycle
                    if ( genea_list(i,3) > genea_list(i+1,3)  ) then
                        !                   if (stopOnChange) call stop_application("bad definition of genealogy [dam, line:"//&
                        !                 trim(str(i))//"]. build a new genealogy file with qtlmap-check.")
                        buf = genea_list(i,:)
                        genea_list(i,:) = genea_list(i+1,:)
                        genea_list(i+1,:) = buf
                        change = .true.
                        countchange=countchange+1
                    end if
                end do
            end do

            j=idnextsire+1
            idnextsire=0

            do i=j,nlist
                if ( nextsire /= genea_list(i,2) ) then
                    idnextsire = i
                    nextsire = genea_list(idnextsire,2)
                    exit
                end if
            end do
        end do

        if ( countchange > 0 ) then
            call log_mess(" ****> pedigree  [change:"//trim(str(countchange))//"]",WARNING_DEF)
            changeSort=.true.
        end if

    end subroutine sort_genea_list
    !!***

    !!****f* m_qtlmap_genealogy/CREATE_STRUCT_GRAND_PARENT
    !!  NAME
    !!    CREATE_STRUCT_GRAND_PARENT
    !!  DESCRIPTION
    !!    Fill ngp,ngm,ngmgp,nrgm,repro,gpere,gmere arrays from the information readed (genea_list,genea_niv).
    !!  INPUTS
    !!   nb_max_indiv   : get the number of animal defined in the genalogy file
    !!
    !!  NOTES
    !!  SOURCE
    SUBROUTINE CREATE_STRUCT_GRAND_PARENT(dg,dgr,nb_max_indiv,genea_list,genea_niv,rac)
        integer,intent(in)                                          :: nb_max_indiv
        character(len=LEN_DEF)        , dimension(:,:),intent(in)   :: genea_list
        type(GENEALOGY_BASE)                         ,intent(inout) :: dg
        type(GENEALOGY_RACE)                         ,intent(inout) :: dgr
        integer   , dimension(-GENERATION_MAX:GENERATION_MAX)     ,intent(in):: genea_niv
        character(len=LEN_DEF)        , dimension(:,:)   ,intent(in)    :: rac

        !local
        character(len=LEN_DEF)            :: ind
        character(len=LEN_DEF)            :: father
        character(len=LEN_DEF)            :: mother
        integer                           :: gen,alloc_stat,l,i,ir1
        character(len=LEN_DEF)            :: word_token
        character(len=LEN_LINE)           :: line_read

        ! nombre de gd meme par gd pere : dim nbre gp
        integer , dimension (:),allocatable :: ngmgp_t
        ! nombre de parent par gd mere
        integer , dimension (:),allocatable :: nrgm_t

        character(len=LEN_DEF) , dimension (:),allocatable :: gmere_t
        character(len=LEN_DEF) , dimension (:),allocatable :: gpere_t
        character(len=LEN_DEF) , dimension (:),allocatable :: repro_t
        character(len=LEN_DEF) , dimension (:),allocatable :: reprop_t
        character(len=LEN_DEF) , dimension (:),allocatable :: reprom_t
        character(len=LEN_DEF) , dimension (:),allocatable :: racep_t
        character(len=LEN_DEF) , dimension (:),allocatable :: racem_t
        logical                                            :: is_ok
        integer                                            :: startGp
        call log_mess('SUBROUTINE : CREATE_STRUCT_GRAND_PARENT',DEBUG_DEF)
        ! Initialize Buffer with the indiv max
        ALLOCATE (ngmgp_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (nrgm_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (gmere_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (gpere_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (reprop_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (reprom_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (repro_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (racep_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (racem_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        dg%nr = 0 ; dg%ngp = 0 ; dg%ngm = 0
        ngmgp_t(1)=0 ; ngmgp_t(2)=0 ; nrgm_t(1)=0 ; nrgm_t(2)=0
         !********------- GEN=1 ------*********
        !lire old..
        !si gen==1
        ind=""

        startGp=sum(genea_niv(:0))+1

        if (genea_niv(1)>0) then
            ind = trim(genea_list(startGp,1))
            father = trim(genea_list(startGp,2))
            mother = trim(genea_list(startGp,3))
            dg%nr = 1
            gpere_t(dg%nr) = trim(genea_list(startGp,2))
            gmere_t(dg%nr) = trim(genea_list(startGp,3))
            repro_t(dg%nr) = trim(genea_list(startGp,1))
            reprop_t(dg%nr) = trim(genea_list(startGp,2))
            reprom_t(dg%nr) = trim(genea_list(startGp,3))
            racep_t(dg%nr) =  trim(rac(startGp,1))
            racem_t(dg%nr) =  trim(rac(startGp,2))
            dg%ngp = 1 ; dg%ngm = 1
            nrgm_t(dg%ngm+1)=nrgm_t(dg%ngm+1)+1
            ngmgp_t(dg%ngp+1)=ngmgp_t(dg%ngp+1)+1
        end if

        do l=startGp+1,startGp+genea_niv(1)-1
            ind = trim(genea_list(l,1))
            father = trim(genea_list(l,2))
            mother = trim(genea_list(l,3))
            dg%nr=dg%nr+1
            repro_t(dg%nr) = ind
            reprop_t(dg%nr) = ind
            reprom_t(dg%nr) = ind
            racep_t(dg%nr) =  trim(rac(l,1))
            racem_t(dg%nr) =  trim(rac(l,2))

            !  print *,'last gp:',trim(gpere_t(ngp)),' current:',trim(father)
            !   print *,'last gm:',trim(gmere_t(ngm)),' current:',trim(mother)
            !New grandfather and grandmother
            if ( gpere_t(dg%ngp) /= father) then
                !print *,'nvx gp:',trim(father)
                dg%ngm=dg%ngm+1
                nrgm_t(dg%ngm+1)=nrgm_t(dg%ngm)+1
                gmere_t(dg%ngm)=mother
                dg%ngp=dg%ngp+1
                gpere_t(dg%ngp)=father
                ! print *,'ngmgp_t[',ngp,']:',ngmgp_t(ngp)
                ngmgp_t(dg%ngp+1)=ngmgp_t(dg%ngp)+1

            !New grandmother
            ELSE IF ( gmere_t(dg%ngm) /= mother) THEN
                 !print *,'nvx gm:',trim(mother)
                dg%ngm=dg%ngm+1
                nrgm_t(dg%ngm+1)=nrgm_t(dg%ngm)+1
                gmere_t(dg%ngm)=mother
                ngmgp_t(dg%ngp+1)=ngmgp_t(dg%ngp+1)+1
             ! print *,ngm,gmere_t(ngm)
            ELSE
                nrgm_t(dg%ngm+1)=nrgm_t(dg%ngm+1)+1
            ENDIF
        END DO

        !ALLOCATES TABLE AND DESALLOCATE BUFFER TAB
        !-------------------------------------------
        ALLOCATE (dg%ngmgp(dg%ngp+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        DO i=1,dg%ngp+1
            dg%ngmgp(i) = ngmgp_t(i)
        END DO
        DEALLOCATE(ngmgp_t)

        ALLOCATE (dg%nrgm(dg%ngm+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        DO i=1,dg%ngm+1
            dg%nrgm(i) = nrgm_t(i)
        END DO
        DEALLOCATE(nrgm_t)

        ALLOCATE (dg%gmere(dg%ngm), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        DO i=1,dg%ngm
            dg%gmere(i) = gmere_t(i)
        END DO
        DEALLOCATE(gmere_t)

        ALLOCATE (dg%gpere(dg%ngp), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        DO i=1,dg%ngp
            dg%gpere(i) = gpere_t(i)
        END DO
        DEALLOCATE(gpere_t)

        ALLOCATE (dg%repro(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dg%reprop(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dg%reprom(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dgr%racep(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dgr%racem(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dg%rep_reprop(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        ALLOCATE (dg%rep_reprom(dg%nr), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%rep_reprop=0; dg%rep_reprom=0
        DO i=1,dg%nr
            dg%repro(i)  = repro_t(i)
            dg%reprop(i) = reprop_t(i)
            dg%reprom(i) = reprom_t(i)
            dgr%racep(i)  = racep_t(i)
            dgr%racem(i)  = racem_t(i)
            !print *, 'RACE',i, repro(i), racep(i), racem(i)
            do ir1=1,dg%nr
                if (dg%reprop(i)==dg%repro(ir1)) dg%rep_reprop(i)=ir1
                if (dg%reprom(i)==dg%repro(ir1)) dg%rep_reprom(i)=ir1
                if (dg%rep_reprom(i).ne.0.and.dg%rep_reprop(i).ne.0)   exit
            end do ! ir1label
        END DO
        DEALLOCATE(repro_t)
        DEALLOCATE(reprop_t)
        DEALLOCATE(reprom_t)

        call log_mess('END SUBROUTINE : CREATE_STRUCT_GRAND_PARENT',DEBUG_DEF)
    END SUBROUTINE CREATE_STRUCT_GRAND_PARENT
    !!***

    !!****f* m_qtlmap_genealogy/CREATE_STRUCT_PARENT
    !!  NAME
    !!    CREATE_STRUCT_PARENT
    !!  DESCRIPTION
    !!    Fill ndm,nmp,pere,mere,animal,nd,nm,np arrays from the information readed (genea_list,genea_niv).
    !!  INPUTS
    !!   nb_max_indiv   : get the number of animal defined in the genalogy file
    !!
    !!  NOTES
    !!  SOURCE
    SUBROUTINE CREATE_STRUCT_PARENT(nb_max_indiv,genea_list,dg,genea_niv)
        integer,intent(in)              ::nb_max_indiv
        character(len=LEN_DEF)        , dimension(:,:),intent(in)   :: genea_list
        type(GENEALOGY_BASE)                         ,intent(inout) :: dg
        integer     , dimension(-GENERATION_MAX:GENERATION_MAX)     ,intent(in)    :: genea_niv


        character(len=LEN_DEF)            ::ind
        character(len=LEN_DEF)            ::father
        character(len=LEN_DEF)            ::mother
        integer                         ::gen,alloc_stat,err,eof,l,i,start
        character(len=LEN_DEF)            ::word_token
        character(len=LEN_LINE)            ::line_read
        integer , dimension (:),allocatable :: ndm_t
        integer , dimension (:),allocatable :: nmp_t

        character(len=LEN_DEF) , dimension (:),allocatable ::mere_t
        character(len=LEN_DEF) , dimension (:),allocatable ::pere_t
        character(len=LEN_DEF) , dimension (:),allocatable ::animal_t
        logical                                          ::is_ok

        call log_mess('SUBROUTINE : CREATE_STRUCT_PARENT',DEBUG_DEF)

        if ( genea_niv(2) <= 0 ) then
            call stop_application("none animals with generation 2 is detected");
        end if

         !***
        ! Initialize Buffer with the indiv max
        ALLOCATE (ndm_t(nb_max_indiv+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (nmp_t(nb_max_indiv+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (mere_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (pere_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        ALLOCATE (animal_t(nb_max_indiv), stat = alloc_stat)
        CALL check_allocate(alloc_stat)

        dg%nd=1
        dg%nm=1
        dg%np=1

        ndm_t(1)=0
        ndm_t(2)=1

        nmp_t(1)=0
        nmp_t(2)=1

        start = sum(genea_niv(:1)) + 1
        ind = trim(genea_list(start,1))
        father = trim(genea_list(start,2))
        mother = trim(genea_list(start,3))

        animal_t(dg%nd) = ind
        pere_t(dg%np) = father
        mere_t(dg%nm) = mother

        DO l=start+1,start+genea_niv(2)-1
            ind = trim(genea_list(l,1))
            father = trim(genea_list(l,2))
            mother = trim(genea_list(l,3))
            dg%nd = dg%nd+1
            animal_t(dg%nd) = ind
            IF ( trim(pere_t(dg%np)) /= trim(father) ) THEN
                dg%nm = dg%nm + 1
                mere_t(dg%nm)= mother
                ndm_t(dg%nm+1)=ndm_t(dg%nm)+1
                dg%np=dg%np+1
                pere_t(dg%np)=father
                nmp_t(dg%np+1)=nmp_t(dg%np)+1
            ELSE IF ( trim(mere_t(dg%nm)) /= trim(mother) ) THEN
                dg%nm = dg%nm + 1
                mere_t(dg%nm)= mother
                ndm_t(dg%nm+1)= ndm_t(dg%nm)+1
                nmp_t(dg%np+1)= nmp_t(dg%np+1)+1
            ELSE
                ndm_t(dg%nm+1)=ndm_t(dg%nm+1)+1
            END IF
        END DO

        !ALLOCATES TABLE AND DESALLOCATE BUFFER TAB
        !-------------------------------------------
        ALLOCATE (dg%ndm(dg%nm+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%ndm = ndm_t(:dg%nm+1)
        DEALLOCATE(ndm_t)

        ALLOCATE (dg%nmp(dg%np+1), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%nmp = nmp_t(:dg%np+1)
        DEALLOCATE(nmp_t)

        ALLOCATE (dg%mere(dg%nm), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%mere = mere_t(:dg%nm)
        DEALLOCATE(mere_t)

        ALLOCATE (dg%pere(dg%np), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%pere = pere_t(:dg%np)
        DEALLOCATE(pere_t)
        ALLOCATE (dg%animal(dg%nd), stat = alloc_stat)
        CALL check_allocate(alloc_stat)
        dg%animal = animal_t(:dg%nd)
        DEALLOCATE(animal_t)


        call log_mess('END SUBROUTINE : CREATE_STRUCT_PARENT',DEBUG_DEF)
    END SUBROUTINE CREATE_STRUCT_PARENT

     !**********************************************
     ! SUBROUTINE : CREATE_STRUCT_DERIVED_GENEALOGY
     !**********************************************
    SUBROUTINE CREATE_STRUCT_DERIVED_GENEALOGY (dg)
        type(GENEALOGY_BASE)                         ,intent(inout) :: dg
        integer                 :: alloc_stat
        integer                 :: ip,im,ir,ifem,i
        call log_mess('SUBROUTINE CREATE_STRUCT_DERIVED_GENEALOGY',DEBUG_DEF)

        do  ip=1,size(dg%pere)
            dg%reppere(ip)=INT_NOT_DEFINED
            do ir=1,size(dg%repro)
                if ( dg%pere(ip) == dg%repro(ir) ) then
                    dg%reppere(ip)= ir
                    exit
                endif
            end do ! irlabel
        end do ! iplabel

        dg%nfem=1
        dg%femelle(dg%nfem)=dg%mere(dg%nfem)
        do  im=1,size(dg%mere)
            dg%repmere(im)=INT_NOT_DEFINED
            do ir=1,size(dg%repro)
                if ( dg%mere(im) == dg%repro(ir) ) then
                    dg%repmere(im)= ir
                    exit
                endif
            end do ! irlabel
            do ifem=1,dg%nfem
                if(dg%mere(im).eq.dg%femelle(ifem)) then
                    dg%repfem(im)=ifem
                    exit
                end if
            end do
            if (ifem > dg%nfem) then
                dg%nfem=dg%nfem+1
                dg%femelle(dg%nfem)=dg%mere(im)
                dg%repfem(im)=dg%nfem
            end if
        end do ! iplabel

        call log_mess('END SUBROUTINE CREATE_STRUCT_DERIVED_GENEALOGY',DEBUG_DEF)
    END SUBROUTINE CREATE_STRUCT_DERIVED_GENEALOGY
    !!***

    !!****f* m_qtlmap_genealogy/write_genea
    !!  NAME
    !!    write_genea
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!   file_name         : path name of the output file
    !!
    !!  NOTES
    !!
    !!  SOURCE
    subroutine write_genea(file_name,dg)
        character(len=*),intent(in)        :: file_name
        type(GENEALOGY_BASE)   ,intent(in) :: dg
        integer :: ip,jm,kr,jgm,igp,id

        open(1,file=file_name)

        call log_mess('TODO:write generation file for genealogy...')
        do igp=1,dg%ngp
            do jgm=dg%ngmgp(igp)+1,dg%ngmgp(igp+1)
                do kr=dg%nrgm(jgm)+1,dg%nrgm(jgm+1)
                    write (1,*) trim(dg%repro(kr)),' ',trim(dg%gpere(igp))&
                    ,' ',trim(dg%gmere(jgm)),' 1'
                end do
            end do
        end do
        do ip=1,dg%np
            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                do id=dg%ndm(jm)+1,dg%ndm(jm+1)
                    write (1,*) trim(dg%animal(id)),' ',&
                    trim(dg%pere(ip)),' ',trim(dg%mere(jm)),' 2'
                end do
            end do
        end do
        close(1)

    end subroutine write_genea
    !!***

    !!****f* m_qtlmap_genealogy/log_debug_genea
    !!  NAME
    !!    log_debug_genea
    !!  DESCRIPTION
    !!
    !!  NOTES
    !!
    !!  SOURCE
    !     subroutine log_debug_genea()
    !       integer  :: i,j,k
    !
    !       do i=1,ngp
    !         print *,'------------------------------------------------'
    !         print *,'index gp:',i,' ngmgp(',i,')=',ngmgp(i),' ngmgp(',(i+1),')=',ngmgp(i+1)
    !         do j=ngmgp(i)+1,ngmgp(i+1)
    !
    !            do k=nrgm(j)+1,nrgm(j+1)
    !              print *,k
    !              print *,trim(repro(k)),' ',trim(gpere(i)),' ',trim(gmere(j)),' 1'
    !            end do
    !         end do
    !       end do
    !     end subroutine log_debug_genea
    !!***

    !put in newdg, the new genealogy
    subroutine create_dataset_genealogy(dataset,array_sample,newdataset,occurences)
        type(QTLMAP_DATASET)             ,intent(in)            :: dataset
        integer             ,dimension(:),intent(in)            :: array_sample
        type(QTLMAP_DATASET)             ,intent(inout)         :: newdataset
        integer                   ,intent(out)     , dimension(dataset%genea%nd)  :: occurences

        character(len=LEN_DEF)        , dimension(:,:),allocatable   :: genea_list
        character(len=LEN_DEF)        , dimension(:,:),allocatable   :: genea_list_enfants,genea_list_parents
        logical                                                      :: find
        integer                                                      :: id,ip,jm,ik,idd
        integer                       , dimension(:),allocatable     :: genea_niv
        type(GENEALOGY_BASE)     ,pointer                   :: dg
        character(len=LEN_DEF)        , dimension(:,:) ,allocatable  :: rac
        logical                       , dimension(dataset%genea%np)  :: pereOk
        logical                       , dimension(dataset%genea%nm)  :: mereOk


        dg => dataset%genea

        allocate (rac(dg%nd,2))
        rac='UNKNOWN'
        allocate (genea_list_enfants(dg%nd,4),genea_list_parents(dg%np+dg%nm,4))
        allocate (genea_list(dg%nd+dg%np+dg%nm,4))
        allocate (genea_niv(-GENERATION_MAX:GENERATION_MAX))
        !allocate (rac(dg%nd,2))

        genea_niv=0
        genea_niv(2)=dg%nd
        pereOk=.false.
        mereOk=.false.
        occurences=0
        genea_list=''

        do id=1,dg%nd
            find = .false.
            do ip=1,dg%np
                do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                    do ik=dg%ndm(jm)+1,dg%ndm(jm+1)
                        if ( ik == array_sample(id)) then
                         occurences(ik)=occurences(ik)+1
                         if ( occurences(ik) == 1) then
                            genea_list_enfants(id,1)=dg%animal(array_sample(id))
                         else
                            genea_list_enfants(id,1)=trim(dg%animal(array_sample(id)))//&
                             "_"//trim(str(occurences(ik)))
                         end if
                            genea_list_enfants(id,2)=dg%pere(ip)
                            genea_list_enfants(id,3)=dg%mere(jm)
                            genea_list_enfants(id,4)='2'
                            find = .true.
                            if (.not. pereOk(ip) ) pereOk(ip)=.true.
                            if (.not. mereOk(jm) ) mereOk(jm)=.true.
                        end if
                        if ( find ) exit
                    end do ! kd
                    if ( find ) exit
                end do !jm
                if ( find ) exit
            end do !ip
        end do !id
        idd=0

        do id=1,dg%np+dg%nm

         if ( id <= dg%np ) then
            if ( .not. pereOk(id) ) cycle
         else
            if ( .not. mereOk(id-dg%np) ) cycle
         end if

         find = .false.

         do ip=1,dg%ngp
          do jm=dg%ngmgp(ip)+1,dg%ngmgp(ip+1)
            do ik=dg%nrgm(jm)+1,dg%nrgm(jm+1)
             if ( id <= dg%np ) then
              if ( dg%repro(ik) == dg%pere(id) ) then
               idd=idd+1
               genea_list_parents(idd,1)=dg%pere(id)
               genea_list_parents(idd,2)=dg%gpere(ip)
               genea_list_parents(idd,3)=dg%gmere(jm)
               genea_list_parents(idd,4)='1'
               find = .true.
              end if
             else
              if ( dg%repro(ik) == dg%mere(id-dg%np) ) then
               idd=idd+1
               genea_list_parents(idd,1)=dg%mere(id-dg%np)
               genea_list_parents(idd,2)=dg%gpere(ip)
               genea_list_parents(idd,3)=dg%gmere(jm)
               genea_list_parents(idd,4)='1'
               find = .true.
              end if
             end if
             if ( find ) exit
            end do ! kd
            if ( find ) exit
          end do ! jm
          if ( find ) exit
         end do !ip
        end do

        genea_niv(1)=idd

        genea_list(:idd,:)=genea_list_parents(:idd,:)
        genea_list(idd+1:idd+dg%nd,:)=genea_list_enfants(:,:)

        deallocate (genea_list_parents,genea_list_enfants)

!        do ik=1,idd+dg%nd
!          print *,genea_list(ik,:4)
!        end do

        call dataset%geneaRace%copy(newdataset%geneaRace)
        call create_genealogy(newdataset,&
                              genea_list,&
                              idd+dg%nd,&
                              genea_niv,&
                              rac)

        deallocate (genea_list)
        deallocate (genea_niv)
        deallocate (rac)

    end subroutine create_dataset_genealogy

    ! Les phenotypes doivent etre existant. On verifie qu'il y a suffisamment de descendants
    ! dans un cas ou l'utilisateur lance une simule avec permutation (>=10)
    subroutine check_permutation_context(dataset)
       type(QTLMAP_DATASET)  ,intent(in)  :: dataset

       type(GENEALOGY_BASE)     ,pointer                   :: dg
       integer :: ic,ip,jm,kd
       integer :: cfs
       call log_mess("genealogy:check_permutation_context",INFO_DEF)
       dg => dataset%genea

       if (dataset%cli%cli_is_permute()) then

        do ic=1,dataset%phenoModel%ncar
         do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if (dataset%phenoAnimal%estime(ic,jm)) then
              cfs=0
              do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
               if ( dataset%phenoAnimal%presentc(ic,kd)) cfs=cfs+1
              end do
              if (cfs < MIN_NDMIN_PERMUTATION ) then
                call stop_application("Number of progeny of the Full-Sib Family (sire "//trim(dg%pere(ip))//&
                ", dam "//trim(dg%mere(jm))//") = "//trim(str(cfs))&
                //". You can not used permutation option with a value below to "//trim(str(MIN_NDMIN_PERMUTATION)))
              end if
             end if
           end do
          end do
         end do

       end if
       call log_mess("genealogy:end check_permutation_context",DEBUG_DEF)

    end subroutine check_permutation_context

END MODULE m_qtlmap_genealogy




