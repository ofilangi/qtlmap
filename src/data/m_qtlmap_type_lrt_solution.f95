module m_qtlmap_type_lrt_solution
    use m_qtlmap_constant
    use m_qtlmap_type_dataset
    use omp_lib

    implicit none

    type CONFIDENCE_INTERVALS_SOLUTION
        !// Name method
        character(len=LEN_DEF)                   ,public :: method = ""
        !// number of average
        integer                                  ,public :: nci      = 0
        !// average 95%,90%,...
        real (kind=dp)   ,dimension(:),pointer           :: ci_seuil => NULL()
        !// lower and upper for each qtl/each threshold
        !// 1:threshold, 2: hypothesys tested, 3: QTL, 4: 1/2
        real (kind=dp)   ,dimension(:,:,:,:),pointer     :: ci_intervals => NULL()

    contains
        procedure, public ::  release => release_confidence_intervals_solution

    end type CONFIDENCE_INTERVALS_SOLUTION

    ! Structure generique pour enregistrer la vraisemblance/effets paternel,... quelque soit la valeur de nqtl
    ! On fixe une valeur maximum au nombre de point teste NQTL <=10
    type LIST_TESTED_POSTION_VALUE
        integer                                         :: sizevaluehyp = 0   ! allocation size
        integer                                         :: sizevaluechr = 0   ! allocation size
        integer                                         :: sizevaluepos = 0   ! allocation size
        ! * Le premier indice est l'hypothese qu'on test (pour nqtl=1 c'est H0 (hyp=1),pour nqtl=2 c'est H0(hyp=1) ou H1(hyp=2)
        ! * deuxieme indice la valeur de l'index calculé sur la combinaison des chromosomes testés
        ! * troisieme indice la valeur de l'index calculé sur la combinaison des positions testés
        real(kind=dp)       ,dimension(:,:,:),pointer   :: values  => null()
    contains

        procedure, public :: release   => release_list_tested_position

        procedure, public :: add       => add_tested_position
        procedure, public :: add1      => add_lrt_tested_position_qtl1
        procedure, public :: add1p     => add_pos_tested_position_qtl1
        procedure, public :: add2      => add_tested_position_qtl2

        procedure, public :: get       => get_tested_position
        procedure, public :: get1      => get_tested_position_qtl1
        procedure, public :: get2      => get_tested_position_qtl2

        procedure, public :: list_available  => list_is_available   ! check the list status
        procedure, public :: pos_available => position_is_available ! check if the value associated with a position existe

        procedure, public :: get_iterator_available_index

    end type LIST_TESTED_POSTION_VALUE


        public :: get_iterator_index

    !!//    nqtl       : hypothesis to test
    !!//    lrtmax     : maximum reached
    !!//    dxmax      : position where are the maximum
    !!//    nxmax      : position in the data structure where the maximum are reached
    !!//    chrmax     : chromosome where the maximum are reached
    !!//
    !!//    lrt1       : LRT curve under Hypothesis One (Only this hypothesis case)
    !!//    pater_eff  : paternal effect curve under Hypothesis One
    !!//    mater_eff  : maternal effect curve under Hypothesis One
    !!//    xlrp       : LRT sires surve under Hypothesis One
    !!//    xlrm       : LRT dams surve under Hypothesis One
    !!//
    !!//    lrt0_2     : LRT curve under Hypothesis Two against Zero (Only this hypothesis case)
    !!//    lrt1_2     : LRT curve under Hypothesis Two against One (Only this hypothesis case)
    !!//    pater_eff2 : paternal effect curve under Hypothesis Two
    !!//    mater_eff2 : maternal effect curve under Hypothesis Two
    !!//    xlrp2      : LRT sires surve under Hypothesis Two against One
    !!//    xlrm2      : LRT dams surve under Hypothesis Two against One
    type TYPE_LRT_SOLUTION
        !// The hypothesis (usually equivalent to nqtl)
        integer                                                   :: hypothesis = -1
        !// LRT Maximum reached
        real (kind=dp)   ,dimension(:),pointer                    :: lrtmax => NULL()
        !// Centimorgan positions of Qtls
        !// real (kind=dp)   ,dimension(:),pointer         :: dxmax => NULL()
        !// Position of Qtl in data structure
        integer          ,dimension(:),pointer                    :: nxmax => NULL()
        !// Chromosome where the position is localised
        integer          ,dimension(:),pointer                    :: chrmax => NULL()

        type(LIST_TESTED_POSTION_VALUE)                           :: LRT
        ! Size = le nombre de pere
        type(LIST_TESTED_POSTION_VALUE)  ,dimension(:),pointer    :: LRT_SIRES => NULL()
        ! Size = le nombre de mere
        type(LIST_TESTED_POSTION_VALUE)  ,dimension(:),pointer    :: LRT_DAMS => NULL()

        !//--------------------- QTL = 1 -----------------------------------
        real (kind=dp)  ,dimension(:,:,:),pointer    :: pater_eff => NULL()
        real (kind=dp)   ,dimension(:,:,:),pointer   :: mater_eff => NULL()

        !//--------------------- QTL = 2 ---------------------------------

        real (kind=dp)  ,dimension(:,:,:,:,:,:),pointer    :: pater_eff2 => NULL()
        real (kind=dp)  ,dimension(:,:,:,:,:,:),pointer    :: mater_eff2 => NULL()

        type (CONFIDENCE_INTERVALS_SOLUTION)  , dimension(:) , pointer :: list_ci => NULL()

    contains

        procedure ,public :: new              => new_lrt_solution
        procedure ,public :: release          => release_lrt_solution

    end type TYPE_LRT_SOLUTION

    !!//    Description of a alert of a correlation to high between classical effect and a qtl effect
    !!//
    !!//    attributes members :
    !!//    --------------------
    !!//    qtl    : number of qtl to test
    !!//    ip     : sire index where the correlation are finding
    !!//    jm     : dams index where the correlation are finding
    !!//    corr   : value of the correlation
    !!//    ntlev  : level of the qtl
    !!//    name_effect  : name of the classical effect
    !!//    name_level   : name of qtl level

    ! Type pour l impression des alerts de correlation trop eleve entre les effets polygenic, de nuisance et les effets qtls
    type CORR_ALERT_TYPE
        integer              :: qtl=0
        integer              :: ip = -1         ! indice du pere sinon -1
        integer              :: jm = -1         ! indice de la mere concerne sinon -1
        real(kind=dp)        :: corr            ! la correlation
        integer              :: ntlev =-1       ! le niveau du qtl concerne
        character(len=LEN_W) :: name_effect     ! le nom de l effet
        character(len=LEN_W) :: name_level      ! le nom du niveau
    end type CORR_ALERT_TYPE


    !!//    Description of a nuisance test
    !!//
    !!//    attributes members :
    !!//    --------------------
    !!//    directeffect  : false if the qtl are in interaction with the effect
    !!//    name          : name without the nuisance
    !!//    df            : freedom degree
    !!//    lrt           : likelihood ratio test
    !!//    pvalue        :

    !type pour l impression de test lin
    type TEST_NUISANCES_TYPE
        logical              :: directeffect    ! intra qtl (fixed effect) or direct effect
        character(len=LEN_W) :: name            ! name without the nuisances
        integer              :: df
        real(kind=dp)        :: lrt
        real(kind=dp)        :: pvalue
    end type TEST_NUISANCES_TYPE

contains


    subroutine add_lrt_tested_position_qtl1(list,dataset,nchr,nposdim,lrtval)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout) :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: nchr
        integer                          , intent(in)    :: nposdim
        real(kind=dp) ,dimension(nchr,nposdim) , intent(in)    :: lrtval

        integer ,dimension(1) :: chr_int
        integer ,dimension(1) :: pos_int
        integer :: ich,ipos

        do ich=1,nchr
            chr_int(1) = ich
            do ipos=1,dataset%map%get_npo(ich)
                pos_int(1) = ipos
                call add_tested_position(list,dataset,1,chr_int,pos_int,lrtval(ich,ipos),1)
            end do
        end do

    end subroutine add_lrt_tested_position_qtl1

    subroutine add_pos_tested_position_qtl1(list,dataset,chr,pos,val)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout) :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: chr
        integer                          , intent(in)    :: pos
        real(kind=dp)                    , intent(in)    :: val

        integer ,dimension(1) :: chr_int
        integer ,dimension(1) :: pos_int

        chr_int(1) = chr
        pos_int(1) = pos
        call add_tested_position(list,dataset,1,chr_int,pos_int,val,1)
    end subroutine add_pos_tested_position_qtl1



    subroutine add_tested_position_qtl2(list,dataset,nchr,npos,lrtval,hyp)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout) :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: nchr
        integer                          , intent(in)    :: npos
        real(kind=dp) ,dimension(nchr,nchr,npos,npos)    , intent(in)    :: lrtval
        integer                          , intent(in)    :: hyp !contre H0 hyp=1, H1=>hyp=2,etc,....

        integer ,dimension(2) :: chr_int
        integer ,dimension(2) :: pos_int

        integer :: ich,ipos,ich2,ipos2,spos

        if ( hyp <= 0 ) then
            call stop_application("LIST_TESTED_POSTION_VALUE::add_tested_position_qtl2"//&
                " hyp value > 0 hyp=["//trim(str(hyp))//"]")
        end if

        do ich=1,nchr
            chr_int(1) = ich
            do ipos=1,dataset%map%get_npo(ich)
                pos_int(1) = ipos
                do ich2=ich,nchr
                    chr_int(2) = ich2
                    spos=1
                    if ( ich == ich2 ) spos=ipos+1
                    do ipos2=spos,dataset%map%get_npo(ich2)
                        pos_int(2) = ipos2
                        call add_tested_position(list,dataset,2,chr_int,pos_int,lrtval(ich,ich2,ipos,ipos2),hyp)
                    end do
                end do
            end do
        end do

    end subroutine add_tested_position_qtl2

    subroutine add_tested_position(list,dataset,nqtl,chr,pos,value,hyp)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout) :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: nqtl
        integer      ,dimension(nqtl)    , intent(in)    :: chr
        integer      ,dimension(nqtl)    , intent(in)    :: pos
        real(kind=dp)                    , intent(in)    :: value
        integer                          , intent(in)    :: hyp !contre H0 hyp=1, H1=>hyp=2,etc,....

        integer       ,parameter :: BLOCK_ALLOCATE_POS = 200

        real(kind=dp)  ,dimension(:,:,:) ,allocatable :: bufvalues

        integer                         :: i,block_pos,chrid,posid

        !$OMP CRITICAL
        !!(add_tested_position)
        call dataset%map%getChridPosId(nqtl,chr,pos,chrid,posid)
        !  print *,hyp,chrid,posid
        ! print *,list%sizevaluechr,list%sizevaluepos
        ! Resize Values
        if ( (list%sizevaluechr < chrid) .or. (list%sizevaluepos < posid)&
            .or. ( list%sizevaluehyp < hyp ) ) then
            !Attention une allocation peut s effectuer dans un environnement multithreadé
            !        print *,"=============================================================================================="
            !        print *,"ALLOCATE sizechr=",list%sizevaluechr," sizepos=",list%sizevaluepos,' chr=',chrid," pos=",posid
            if ( chrid > list%sizevaluechr .and. list%sizevaluechr /= 0 ) then
                block_pos = max(list%sizevaluepos,BLOCK_ALLOCATE_POS)
                block_pos = max(block_pos,posid)
            else if ( list%sizevaluehyp < hyp .and. list%sizevaluehyp /= 0  ) then
                block_pos = max(list%sizevaluepos,BLOCK_ALLOCATE_POS)
                block_pos = max(block_pos,posid)
            else
                block_pos = max(list%sizevaluepos + BLOCK_ALLOCATE_POS,posid)
            end if

            if ( list%sizevaluechr > 0 ) then
                allocate (bufvalues(list%sizevaluehyp,list%sizevaluechr,list%sizevaluepos))
                bufvalues = list%values
                deallocate(list%values)
            end if

            allocate (list%values(nqtl,chrid,list%sizevaluepos+block_pos))

            if ( list%sizevaluechr > 0 ) then
                list%values(:size(bufvalues,1),:size(bufvalues,2),:size(bufvalues,3)) = bufvalues
                deallocate(bufvalues)
            end if

            list%sizevaluehyp = hyp
            list%sizevaluechr = chrid
            list%sizevaluepos = list%sizevaluepos + block_pos

        end if

        list%values(hyp,chrid,posid) = value
     !$OMP END CRITICAL
     !!(add_tested_position)
    end subroutine add_tested_position



    function get_tested_position(list,dataset,n,chr,pos,hyp) result(value)
        class(LIST_TESTED_POSTION_VALUE) , intent(in)    :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          ,intent(in)     :: n
        integer     ,dimension(n)        ,intent(in)     :: chr
        integer     ,dimension(n)        ,intent(in)     :: pos
        integer                          , intent(in)    :: hyp !contre H0 hyp=1, H1=>hyp=2,etc,....

        real(kind=dp)   :: value
        integer         :: chrid,dimid,posid,i

        call dataset%map%getChridPosId(n,chr,pos,chrid,posid)

        if ( (chrid > list%sizevaluechr) .or.(posid > list%sizevaluepos) ) then
            call log_mess("Position ( {Chr} {Pos} ) does not exit ! :",DEBUG_DEF)
            do i=1,n
             call log_mess("{"//trim(str(chr(i)))//"} {"//trim(str(pos(i)))//"}",DEBUG_DEF)
            end do
!            call stop_application("Devel error: LIST_TESTED_POSTION_VALUE::get_tested_position("//&
!                trim(str(chrid))//","//trim(str(posid))//&
!                ")  - DIM_LIST=("//&
!                trim(str(list%sizevaluechr))//","//trim(str(list%sizevaluepos))//")")
            value = 0.d0
        else
         value = list%values(hyp,chrid,posid)
        end if

    end function get_tested_position


    function position_is_available(list,dataset,n,chr,pos) result(value)
        class(LIST_TESTED_POSTION_VALUE) , intent(in)    :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          ,intent(in)     :: n
        integer     ,dimension(n)        ,intent(in)     :: chr
        integer     ,dimension(n)        ,intent(in)     :: pos

        logical :: value
        integer :: chrid,posid

        call dataset%map%getChridPosId(n,chr,pos,chrid,posid)

        value = .true.

        if ( (chrid > list%sizevaluechr) .or.(posid > list%sizevaluepos) ) then
          value = .false.
        end if

    end function position_is_available

    function get_tested_position_qtl1(list,dataset,chr,pos) result(value)
        class(LIST_TESTED_POSTION_VALUE) , intent(in)    :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: chr
        integer                          , intent(in)    :: pos

        real(kind=dp)     :: value

        integer ,dimension(1) :: chr_int
        integer ,dimension(1) :: pos_int

        chr_int(1) = chr
        pos_int(1) = pos

        value = get_tested_position(list,dataset,1,chr_int,pos_int,1)

    end function get_tested_position_qtl1


    function get_tested_position_qtl2(list,dataset,chr1,chr2,pos1,pos2,hyp) result(value)
        class(LIST_TESTED_POSTION_VALUE) , intent(in)    :: list
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          , intent(in)    :: chr1,chr2
        integer                          , intent(in)    :: pos1,pos2
        integer                          , intent(in)    :: hyp !contre H0 hyp=1, H1=>hyp=2,etc,....
        real(kind=dp)  :: value

        integer ,dimension(2) :: chr_int
        integer ,dimension(2) :: pos_int

        chr_int(1) = chr1
        pos_int(1) = pos1
        chr_int(2) = chr2
        pos_int(2) = pos2

        value = get_tested_position(list,dataset,2,chr_int,pos_int,hyp)

    end function get_tested_position_qtl2


    ! Donne la liste des positions et des chromosomes sous forme de tableau à 1 dimension pour parcourir
    ! l'ensemble des position possible lineairement (une unique boucle for)
    subroutine get_iterator_index(dataset,n,listChrRes,listPosRes,nGLTotal)
        type(QTLMAP_DATASET)             , intent(in)    :: dataset
        integer                          ,intent(in)     :: n
        integer     ,dimension(:,:),pointer ,intent(inout)  :: listChrRes
        integer     ,dimension(:,:),pointer ,intent(inout)  :: listPosRes
        integer                          ,intent(inout)     :: nGLTotal

        integer :: nGL,i,iqtl,nlinValide,chr,nlin
        integer ,dimension(n) :: listN,listChr
        integer ,dimension(:,:) ,pointer :: lchr,lnpos
        logical :: ok

        call log_mess("m_qtlmap_type_lrt_solution : get_iterator_index",DEBUG_DEF)

        !Nombre de point sur le groupe de liaison
        !on cherche le nombre de point dependant de l hyp nqtl du numbre de point sur le groupe de liason
        nGL=0
        do chr=1,dataset%map%nchr
            nGL=nGL + dataset%map%get_npo(chr)
        end do

        nGLTotal=1
        do iqtl=1,n
            nGLTotal=nGLTotal*nGL
        end do

        listN(n)=0
        listN(1:n-1)=1
        listChr=1

        allocate (lchr(nGLTotal,n))
        allocate (lnpos(nGLTotal,n))

        nlinValide=0
        do nlin=1,nGLTotal
            i=n
            ok=.false.
            do while ( (.not. ok) .and. (i>=1))
                listN(i) = listN(i)+1
                if ( listN(i) > dataset%map%get_npo(listChr(i)) ) then
                    if ( dataset%map%nchr > listChr(i) ) then
                        listN(i)=1
                        listChr(i)=listChr(i)+1
                        ok=.true.
                    else if ( i>1 ) then
                        listN(i)=listN(i-1)
                        listChr(i)=listChr(i-1)
                        i=i-1
                    else
                        listN(i)=1
                        listChr(i)=1
                        i=i-1
                    end if
                else
                    ok = .true.
                end if
            end do
            ok=.true.
            do iqtl=2,n
                ok = ok .and. ( listN(iqtl-1)<listN(iqtl) )
            end do
            ! Afichage des indexes générés
            !print *,(listChr(i),listN(i),i=1,n),ok

            if ( ok ) then
                nlinValide=nlinValide+1
                lchr(nlinValide,:)=listChr
                lnpos(nlinValide,:)=listN
            end if
            !condition d arret
            do iqtl=n,1,-1
             ok = ( listN(iqtl) == (dataset%map%get_npo(dataset%map%nchr) - (iqtl-1)) &
             .and. ( listChr(iqtl) == dataset%map%nchr ))
            end do
            if ( ok ) exit ! on sort, c'est le dernier point
        end do
        nGLTotal=nlinValide

        allocate (listChrRes(nlinValide,n))
        allocate (listPosRes(nlinValide,n))

        listChrRes = lchr(:nlinValide,:)
        listPosRes = lnpos(:nlinValide,:)

        deallocate (lchr)
        deallocate (lnpos)


    end subroutine get_iterator_index



    subroutine get_iterator_available_index(list,dataset,n,listChrRes,listPosRes,nGLTotal)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout)    :: list
        type(QTLMAP_DATASET)             , intent(in)       :: dataset
        integer                          ,intent(in)        :: n
        integer     ,dimension(:,:),pointer ,intent(inout)  :: listChrRes
        integer     ,dimension(:,:),pointer ,intent(inout)  :: listPosRes
        integer                          ,intent(inout)     :: nGLTotal

        integer     ,dimension(:,:),pointer  :: listChrRes_tp
        integer     ,dimension(:,:),pointer  :: listPosRes_tp
        integer                              :: nGLTotal_tp

        integer :: i

        call get_iterator_index(dataset,n,listChrRes_tp,listPosRes_tp,nGLTotal_tp)

        allocate (listChrRes(nGLTotal_tp,n),listPosRes(nGLTotal_tp,n))
        nGLTotal = 0
        do i=1,nGLTotal_tp
           if ( list%pos_available(dataset,n,listChrRes_tp(i,:),listPosRes_tp(i,:)) ) then
             nGLTotal = nGLTotal + 1
             listChrRes(nGLTotal,:)=listChrRes_tp(i,:)
             listPosRes(nGLTotal,:)=listPosRes_tp(i,:)
           end if
        end do

        deallocate (listChrRes_tp,listPosRes_tp)

    end subroutine get_iterator_available_index


    subroutine release_list_tested_position(list)
        class(LIST_TESTED_POSTION_VALUE) , intent(inout) :: list

        if (associated(list%values)) deallocate (list%values)
        list%sizevaluechr = 0
        list%sizevaluepos = 0

    end subroutine release_list_tested_position


    function list_is_available(list) result(avail)
        class(LIST_TESTED_POSTION_VALUE) , intent(in) :: list

        logical :: avail

        avail = (list%sizevaluechr>0) .and. (list%sizevaluepos>0)

    end function list_is_available

    subroutine release_confidence_intervals_solution(this)
        class(CONFIDENCE_INTERVALS_SOLUTION) , intent(inout) :: this

        if ( associated(this%ci_seuil)) deallocate (this%ci_seuil)
        if ( associated(this%ci_intervals)) deallocate (this%ci_intervals)

    end subroutine release_confidence_intervals_solution

    !!   initialize a variable of type TYPE_LRT_SOLUTION
    subroutine new_lrt_solution(lrtsol,dataset,hypothesis,nbqtl)
        class(TYPE_LRT_SOLUTION) , intent(inout) :: lrtsol
        type(QTLMAP_DATASET) , intent(in)       :: dataset
        integer   , intent(in) :: hypothesis ! le test d'hypothese
        integer   , intent(in), optional :: nbqtl ! le nombre de qtl dans le model

        integer  :: npo,i,iq,ip,jm,nbq

        type(GENEALOGY_BASE) ,pointer :: dg
        type(MAP_BASE)       ,pointer :: map

        dg => dataset%genea
        map => dataset%map

        nbq = hypothesis

        if ( present(nbqtl) ) nbq = nbqtl

        ! solution under hypothesis nqtl
        lrtsol%hypothesis=hypothesis

        ! begin index <--0
        allocate (lrtsol%lrtmax(0:hypothesis-1))
        allocate (lrtsol%nxmax(0:nbq-1))
        allocate (lrtsol%chrmax(0:nbq-1))

        lrtsol%lrtmax = 0.d0
        lrtsol%nxmax = 0
        lrtsol%chrmax = 0

        allocate (lrtsol%LRT_SIRES(dataset%genea%np))
        allocate (lrtsol%LRT_DAMS(dataset%genea%nm))

        if ( hypothesis == 1 .or. hypothesis == 2 ) then
            npo = map%get_maxnpo()
            if ( hypothesis == 1) then
                !           allocate (lrtsol%lrt1(map%nchr,npo))
                allocate (lrtsol%pater_eff(map%nchr,dg%np,npo))
                allocate (lrtsol%mater_eff(map%nchr,dg%nm,npo))
                !           allocate (lrtsol%xlrp(map%nchr,dg%np,npo))
                !           allocate (lrtsol%xlrm(map%nchr,dg%nm,npo))
                lrtsol%pater_eff=0.d0
                lrtsol%mater_eff=0.d0
            !           lrtsol%xlrp=0.d0
            !           lrtsol%xlrm=0.d0
            end if

            if ( hypothesis == 2 ) then
                !            allocate (lrtsol%lrt0_2(map%nchr,map%nchr,npo,npo))
                !            allocate (lrtsol%lrt1_2(map%nchr,map%nchr,npo,npo))
                allocate (lrtsol%pater_eff2(map%nchr,map%nchr,dg%np,npo,npo,2))
                allocate (lrtsol%mater_eff2(map%nchr,map%nchr,dg%nm,npo,npo,2))
                !            allocate (lrtsol%xlrp2(map%nchr,map%nchr,dg%np,npo,npo))
                !            allocate (lrtsol%xlrm2(map%nchr,map%nchr,dg%nm,npo,npo))
                !            lrtsol%lrt0_2=0.d0
                !            lrtsol%lrt1_2=0.d0
                lrtsol%pater_eff2=0.d0
                lrtsol%mater_eff2=0.d0
            !            lrtsol%xlrp2=0.d0
            !            lrtsol%xlrm2=0.d0
            end if
        end if
    end subroutine new_lrt_solution

    !!   release a variable of type TYPE_LRT_SOLUTION
    subroutine release_lrt_solution(lrtsol)
        class(TYPE_LRT_SOLUTION) , intent(inout) :: lrtsol
        integer :: i,ip,jm

        if ( lrtsol%hypothesis < 0 ) return

        if ( associated (lrtsol%lrtmax) ) deallocate (lrtsol%lrtmax)
        if ( associated (lrtsol%chrmax) ) deallocate (lrtsol%chrmax)
        if ( associated (lrtsol%nxmax) ) deallocate (lrtsol%nxmax)

        lrtsol%lrtmax=>null()
        lrtsol%chrmax=>null()
        lrtsol%nxmax=>null()

        do ip=1,size(lrtsol%LRT_SIRES,1)
            call lrtsol%LRT_SIRES(ip)%release()
        end do
        do jm=1,size(lrtsol%LRT_DAMS,1)
            call lrtsol%LRT_DAMS(jm)%release()
        end do

        call lrtsol%LRT%release()

        deallocate(lrtsol%LRT_SIRES)
        deallocate(lrtsol%LRT_DAMS)

        if ( lrtsol%hypothesis == 1 ) then
            !        if ( associated (lrtsol%lrt1) ) deallocate (lrtsol%lrt1)
            if ( associated (lrtsol%pater_eff) ) deallocate (lrtsol%pater_eff)
            if ( associated (lrtsol%mater_eff) ) deallocate (lrtsol%mater_eff)
        !        if ( associated (lrtsol%xlrp) ) deallocate (lrtsol%xlrp)
        !        if ( associated (lrtsol%xlrm) ) deallocate (lrtsol%xlrm)
        end if

        !       lrtsol%lrt1=>null()
        lrtsol%pater_eff=>null()
        lrtsol%mater_eff=>null()
        !       lrtsol%xlrp=>null()
        !       lrtsol%xlrm=>null()

        if ( lrtsol%hypothesis == 2 ) then
            !         if ( associated (lrtsol%lrt0_2) ) deallocate (lrtsol%lrt0_2)
            !         if ( associated (lrtsol%lrt1_2) ) deallocate (lrtsol%lrt1_2)
            if ( associated (lrtsol%pater_eff2) ) deallocate (lrtsol%pater_eff2)
            if ( associated (lrtsol%mater_eff2) ) deallocate (lrtsol%mater_eff2)
        !         if ( associated (lrtsol%xlrp2) ) deallocate (lrtsol%xlrp2)
        !         if ( associated (lrtsol%xlrm2) ) deallocate (lrtsol%xlrm2)
        end if

        !       lrtsol%lrt0_2=>null()
        !       lrtsol%lrt1_2=>null()
        lrtsol%pater_eff2=>null()
        lrtsol%mater_eff2=>null()
        !      lrtsol%xlrp2=>null()
        !      lrtsol%xlrm2=>null()

        if ( associated (lrtsol%list_ci) ) then
            do i=1,size(lrtsol%list_ci)
                call lrtsol%list_ci(i)%release()
            end do
            deallocate (lrtsol%list_ci)
        end if

    end subroutine release_lrt_solution

end module m_qtlmap_type_lrt_solution
