!!****m* ANALYSE/m_qtlmap_incidence_multi
!!  NAME
!!    m_qtlmap_incidence_multi
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
module m_qtlmap_incidence_multi
    use m_qtlmap_types
    use m_qtlmap_math
    use m_qtlmap_log
    use m_qtlmap_incidence

    implicit none
    save

    type(DATASET_QTLMAP_TRIO) , pointer :: p_dataTrio
    type(GENEALOGY_BASE) , pointer :: p_dg
    type(PDD_BUILD)      , pointer :: p_spt
    type(PHENOTYPE_BASE) , pointer :: p_dpa
    type(DATAMODEL_BASE) , pointer :: p_dpm

    !!****v* m_qtlmap_incidence_multi/estime_multi
    !!  NAME
    !!   estime_multi
    !!  DESCRIPTION
    !!   Estimability per progeny phenotyped for all trait combined
    !!  NOTES
    !!
    !!***
    logical , dimension(:), allocatable  ,private   :: estime_multi
    !!****v* m_qtlmap_incidence_multi/estfem_multi
    !!  NAME
    !!   estfem_multi
    !!  DESCRIPTION
    !!   Estimability per progeny genotyped for all trait combined
    !!  NOTES
    !!
    !!***
    logical , dimension(:), allocatable  ,private   :: estfem_multi
    !!****v* m_qtlmap_incidence_multi/iam_multi
    !!  NAME
    !!   iam_multi
    !!  DESCRIPTION
    !!   Indexe of the femal for all trait combined
    !!  NOTES
    !!
    !!***
    integer , dimension(:), allocatable  ,private   :: iam_multi
    !!****v* m_qtlmap_incidence_multi/namest_multi
    !!  NAME
    !!   namest_multi
    !!  DESCRIPTION
    !!   Number of female with estimability for all trait combined
    !!  NOTES
    !!
    !!***
    integer                              ,private   :: namest_multi
    !!****v* m_qtlmap_incidence_multi/dataset
    !!  NAME
    !!   dataset
    !!  DESCRIPTION
    !!   The current dataset using bi the likelihood function
    !!  NOTES
    !!
    !!***
    type(DATASET_TYPE)                    ,pointer    ,private         :: dataset     => null()
    !!****v* m_qtlmap_incidence_multi/ntnivmaxtotal
    !!  NAME
    !!   ntnivmaxtotal
    !!  DESCRIPTION
    !!   maximum level finded from all contingence matrix
    !!  NOTES
    !!
    !!***
    integer                                       ,private             :: ntnivmaxtotal
    !!****v* m_qtlmap_incidence_multi/my_listdesc
    !!  NAME
    !!   my_listdesc
    !!  DESCRIPTION
    !!   description of each contingence matrix
    !!  NOTES
    !!
    !!***
    type(INCIDENCE_TYPE)   ,dimension(:)  ,pointer    ,private         :: my_listdesc => null()
    !!****v* m_qtlmap_incidence_multi/my_xincreduitmul
    !!  NAME
    !!   my_xincreduitmul
    !!  DESCRIPTION
    !!   contingence matrix for each trait
    !!  NOTES
    !!
    !!***
    real (kind=dp),dimension(:,:,:),allocatable,private                :: my_xincreduitmul
    !!****v* m_qtlmap_incidence_multi/current_chr
    !!  NAME
    !!   current_chr
    !!  DESCRIPTION
    !!   current linkage group tested
    !!  NOTES
    !!
    !!***
    integer , dimension(:)   ,pointer    ,private                      :: current_chr

    !$omp threadprivate (dataset,my_listdesc,my_xincreduitmul,current_chr)

    public :: init_incidence_multi
    public :: end_incidence_multi

contains
    !!****f* m_qtlmap_incidence_multi/init_incidence_multi
    !! NAME
    !!   init_incidence_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine init_incidence_multi(data,spt,nqtl,xinc,listDesc,effects,resolution_LU)
        type(QTLMAP_DATASET)                  ,intent(inout)      :: data
        type(PDD_BUILD)                       ,intent(inout)      :: spt
        integer                               , intent(in)        :: nqtl
        type(type_effect_contingence)         , intent(in)        :: effects
        real (kind=dp) , dimension(:,:,:), pointer, intent(inout) :: xinc
        type(INCIDENCE_TYPE)     ,dimension(data%phenoModel%ncar), intent(inout)  :: listDesc
        logical, intent(in)                                       :: resolution_LU

        integer :: ip,jm,kd,ifem,eff,nteffmaxtotal,ic,nparmax,k,jc,s
      
        ! General array to get index information about femal estimable
        !                      ***
        allocate (estime_multi(data%genea%nm))
        allocate (estfem_multi(data%genea%nfem))
        allocate (iam_multi(data%genea%nfem))
        estime_multi=.false.
        estfem_multi=.false.
        namest_multi=0
        iam_multi=0

        do ip=1,data%genea%np
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                eff=0
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if(count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                        eff=eff+1
                    end if
                end do
                ifem=data%genea%repfem(jm)
                if( eff>=data%params%NDMIN ) then
                    estime_multi(jm)=.true.
                    estfem_multi(ifem)=.true.
                end if
            end do !jm
        end do !ip

        do ifem=1,data%genea%nfem
            if(estfem_multi(ifem)) then
                namest_multi=namest_multi+1
                iam_multi(ifem)=namest_multi
            end if
        end do

        ntnivmaxtotal=0
        nteffmaxtotal=0
        do ic=1,data%phenoModel%ncar
            listDesc(ic)%ic=ic
            listDesc(ic)%nqtl=nqtl
            call set_ntnivmax(data,ic,nqtl,listDesc(ic)%ntnivmax,&
                listDesc(ic)%nteffmax,listDesc(ic)%ntlev,listDesc(ic)%nbniv)
            ntnivmaxtotal=max(ntnivmaxtotal,listDesc(ic)%ntnivmax)
            nteffmaxtotal=max(nteffmaxtotal,listDesc(ic)%nteffmax)

            if ( effects%traits > 0 ) then
                listDesc(ic)%ntnivmax=listDesc(ic)%ntnivmax+data%phenoModel%ncar
                listDesc(ic)%nteffmax=listDesc(ic)%nteffmax+1
            end if


            listDesc(ic)%ntniv    = 0
            listDesc(ic)%nbnivest = 0
            listDesc(ic)%nteff    = 0

            call log_mess("********************* INIT CARAC="// trim(str(ic))//" *****************************",DEBUG_DEF)
            call log_mess("NTNIVMAX:"//str(listDesc(ic)%ntnivmax),DEBUG_DEF)
            call log_mess("NTEFFMAX:"//str(listDesc(ic)%nteffmax),DEBUG_DEF)
            call log_mess("NTLEVQTL / Reproductor:"//str(listDesc(ic)%ntlev),DEBUG_DEF)
            call log_mess("NBLEV FIXED EFFECT:"//str(listDesc(ic)%nbniv),DEBUG_DEF)
            call log_mess("*********************************************************",DEBUG_DEF)

            allocate(listDesc(ic)%desc(listDesc(ic)%nteffmax))
            allocate(listDesc(ic)%vecsol(listDesc(ic)%ntnivmax))
            listDesc(ic)%vecsol=.false.
            !  allocate(listDesc(ic)%corniv(ntnivmax))
            allocate(listDesc(ic)%precis(listDesc(ic)%ntnivmax))
            listDesc(ic)%precis=0.d0

            if ( nqtl > 0) then
                allocate (listDesc(ic)%ntniv_qtlsires(nqtl))
                allocate (listDesc(ic)%ntniv_qtldams(nqtl))
            end if
            allocate (listDesc(ic)%corr_niv_nivb(listDesc(ic)%ntnivmax))
        end do

        allocate (xinc(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal))
        xinc=0.d0

        allocate (listDesc(1)%dataset)
        call init_dataset_ncar(data,spt,nqtl,listDesc(1)%dataset)

        do ic=1,data%phenoModel%ncar
            allocate (listDesc(ic)%nonull(listDesc(1)%dataset%nkd,ntnivmaxtotal))
            allocate (listDesc(ic)%nbnonull(listDesc(1)%dataset%nkd))
        end do

        nparmax=data%phenoModel%ncar*data%genea%np+&
            data%phenoModel%ncar*(data%phenoModel%ncar-1)/2+data%phenoModel%ncar*ntnivmaxtotal

        allocate (listDesc(1)%borni(nparmax))
        allocate (listDesc(1)%borns(nparmax))
        allocate (listDesc(1)%par(nparmax))
      
        if ( resolution_LU ) then
            listDesc(1)%par(1:(data%genea%np-1)*data%phenoModel%ncar)=1.d0 ! matH
            listDesc(1)%par((data%genea%np-1)*data%phenoModel%ncar+1:(data%genea%np-1)*&
                data%phenoModel%ncar+data%phenoModel%ncar*(data%phenoModel%ncar+1)/2)=0.d0 ! matL

            s=(data%genea%np-1)*data%phenoModel%ncar+1
            listDesc(1)%par(s)=1.d0

            do ic=2,data%phenoModel%ncar
                listDesc(1)%par(s+ic)=1.d0
                s = s + ic
            end do

            listDesc(1)%borni(:)=DEFAULT_PARAM_MIN
            listDesc(1)%borns(:)=DEFAULT_PARAM_MAX

            listDesc(1)%borni(1:(data%genea%np-1)*data%phenoModel%ncar)=0.1d0

            listDesc(1)%par(data%genea%np*data%phenoModel%ncar+data%phenoModel%ncar*(data%phenoModel%ncar-1)/2+1:)=0.d0

           !	print *,"**********>PAR:",listDesc(1)%par(:(np-1)*ncar)
        else

            listDesc(1)%borni(1:data%phenoModel%ncar*data%genea%np)=SIG_MIN
            listDesc(1)%borns(1:data%phenoModel%ncar*data%genea%np)=SIG_MAX

            k=data%phenoModel%ncar*data%genea%np
            !bornes des correlations residuelles
            do ic=2,data%phenoModel%ncar
                do jc=1,ic-1
                    k=k+1
                    ! ** MODIF RESCALE ** OFI =>
                    !listDesc(1)%borni(k)=-1.d0
                    !listDesc(1)%borns(k)=1.d0
                    listDesc(1)%borni(k)=DEFAULT_PARAM_MIN
                    listDesc(1)%borns(k)=DEFAULT_PARAM_MAX
                end do
            end do

            listDesc(1)%borni(k+1:)=DEFAULT_PARAM_MIN
            listDesc(1)%borns(k+1:)=DEFAULT_PARAM_MAX

            !initialisation des variances
            do ic=1,data%phenoModel%ncar
                do ip=1,data%genea%np
                    listDesc(1)%par((ic-1)*data%genea%np+ip)=listDesc(1)%dataset%lSires(ip)%sig0(ic)
                end do
            end do

            listDesc(1)%par(data%phenoModel%ncar*data%genea%np+1:)=0.d0

        end if
      
        ! Add Avril 2013 => Computation of sub LRT (each full and half family)
        allocate (listDesc(1)%fperemax(data%genea%np))
        allocate (listDesc(1)%fmeremax(data%genea%nm))

        do ic=2,data%phenoModel%ncar
            listDesc(ic)%dataset=>listDesc(1)%dataset
            listDesc(ic)%borni=>listDesc(1)%borni
            listDesc(ic)%borns=>listDesc(1)%borns
            listDesc(ic)%par=>listDesc(1)%par
        end do

    end subroutine init_incidence_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/end_incidence_multi
    !! NAME
    !!   end_incidence_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine end_incidence_multi(data,listDesc)
        type(QTLMAP_DATASET)       ,intent(in) :: data
        type(INCIDENCE_TYPE)     ,dimension(data%phenoModel%ncar), intent(inout)  :: listDesc
        integer :: ic

        deallocate (estime_multi)
        deallocate (iam_multi)
        deallocate (estfem_multi)

        call end_dataset(listDesc(1)%dataset)
        deallocate (listDesc(1)%dataset)
        do ic=2,data%phenoModel%ncar
            listDesc(ic)%dataset=>null()
            listDesc(ic)%borni=>null()
            listDesc(ic)%borns=>null()
            listDesc(ic)%par=>null()
        end do

        do ic=1,data%phenoModel%ncar
            call end_incidence(listDesc(ic))
        end do

    end subroutine end_incidence_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/init_dataset_ncar
    !! NAME
    !!   init_dataset_ncar
    !! DESCRIPTION
    !!
    !! NOTE
    !!  constition des famille de plein et demi frere via des indexes
    !!  pas de notion d estimable...on comptabilise seulement les effectifs et
    !!  on repere les lignes de la matrice de contingence
    !! SOURCE
    subroutine init_dataset_ncar(data,spt,nqtl,dataset)
        type(QTLMAP_DATASET)       ,target     ,intent(inout) :: data
        type(PDD_BUILD)            ,target     ,intent(inout) :: spt
        type (DATASET_TYPE)                    ,intent(inout) :: dataset
        integer                                ,intent(in)    :: nqtl
        logical ,dimension(data%genea%nd) :: pres
        integer :: ip,jm,kd,ifem,kdd,s1,s2,chr,iq,ic,kkd,effp,efft,jc,eff
        real(kind=dp) :: somyp(data%phenoModel%ncar),somy(data%phenoModel%ncar)

        call log_mess("build_family for incidence multitrait analysis...",DEBUG_DEF)

        dataset%data=>data
        dataset%spt=>spt

        estime_multi=.false.
        allocate(dataset%lSires(data%genea%np))
        kdd=0
        do ip=1,data%genea%np
            allocate(dataset%lSires(ip)%full_sib(data%genea%nmp(ip+1)-data%genea%nmp(ip)))
            ifem=0
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                ifem=ifem+1
                kd=data%genea%ndm(jm)+1
                do while ( kd<=data%genea%ndm(jm+1) .and. .not. (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) )
                    kd=kd+1
                end do
                if (kd>data%genea%ndm(jm+1)) then
                    dataset%lSires(ip)%full_sib(ifem)%firstKd=-1
                    dataset%lSires(ip)%full_sib(ifem)%lastKd=-2
                    cycle
                end if

                kdd=kdd+1
                ! on a trouve le premier kd valide du pere ip dans la matrice d incidence
                dataset%lSires(ip)%full_sib(ifem)%firstKd=kdd
                do kd=kd+1,data%genea%ndm(jm+1)
                    if ( (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) ) then
                        kdd=kdd + 1
                    end if
                end do
                dataset%lSires(ip)%full_sib(ifem)%lastKd=kdd

                !NOUVEAU ESTIME....
                estime_multi(jm)=&
                    (dataset%lSires(ip)%full_sib(ifem)%lastKd-dataset%lSires(ip)%full_sib(ifem)%firstKd)>=data%params%NDMIN

                call log_mess("Full Sib ["//trim(data%genea%pere(ip))//"-"//trim(data%genea%mere(jm))//"] :"//&
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
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                ifem=ifem+1
                if (dataset%lSires(ip)%full_sib(ifem)%firstKd<0) cycle
                dataset%lSires(ip)%half_sib%firstKd= dataset%lSires(ip)%full_sib(ifem)%firstKd
                exit
            end do
            ifem=data%genea%nmp(ip+1)-data%genea%nmp(ip)+1
            do jm=data%genea%nmp(ip+1),data%genea%nmp(ip)+1,-1
                ifem=ifem-1
                if (dataset%lSires(ip)%full_sib(ifem)%lastKd<0) cycle
                dataset%lSires(ip)%half_sib%lastKd=dataset%lSires(ip)%full_sib(ifem)%lastKd
                exit
            end do

            call log_mess("Half Sib ["//trim(data%genea%pere(ip))//"] :"//&
                trim(str(dataset%lSires(ip)%half_sib%firstKd))//"->"//trim(str(dataset%lSires(ip)%half_sib%lastKd)),&
                VERBOSE_DEF)
            !s1=lSires(ip)%half_sib%lastKd-lSires(ip)%half_sib%firstKd+1
            if ( nqtl >= 1) then
                allocate(dataset%lSires(ip)%ppt(nqtl,&
                    dataset%lSires(ip)%half_sib%firstKd:dataset%lSires(ip)%half_sib%lastKd))
            end if
        end do ! ip

        do ic=2,data%phenoModel%ncar
            dataset%lSires = dataset%lSires
        end do

        dataset%nkd=0
        dataset%nkd_max_by_fam=0

        do ip=1,data%genea%np
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                eff=0
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                        dataset%nkd=dataset%nkd+1
                        eff=eff+1
                    end if
                end do !kd
                dataset%nkd_max_by_fam=max(dataset%nkd_max_by_fam,eff)
            end do ! jm
        end do ! ip

        allocate (dataset%Y(data%phenoModel%ncar,dataset%nkd))
        allocate (dataset%CD(data%phenoModel%ncar,dataset%nkd))

        do ic=1,data%phenoModel%ncar
            kkd=0
            do ip=1,data%genea%np
                do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                    do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                        if (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                            kkd=kkd+1
                            if ( data%phenoModel%natureY(ic) /= 'i' ) then
                                dataset%Y(ic,kkd)=data%phenoAnimal%y(ic,kd)
                                dataset%CD(ic,kkd)=data%phenoAnimal%cd(ic,kd)
                            else
                                call stop_application("Type of trait 'i' are not managed in multi-trait analysis!")
                                dataset%YDISCRETORD(ic,kkd)=data%phenoAnimal%ydiscretord(ic,kd)
                            end if
                        end if
                    end do !kd
                end do ! jm
            end do ! ip
        end do !ic

        somy=0
        efft=0
        do ip=1,data%genea%np
            allocate(dataset%lSires(ip)%sig0(data%phenoModel%ncar))
            allocate(dataset%lSires(ip)%xmu0(data%phenoModel%ncar))
            effp=0
            somyp=0
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                        effp=effp+1
                        do ic=1,data%phenoModel%ncar
                            somyp(ic)=somyp(ic)+data%phenoAnimal%y(ic,kd)!*cd(ic,kd)
                        end do
                    end if
                end do
            end do

            !somme total sur le caractere
            do ic=1,data%phenoModel%ncar
                somy(ic)=somy(ic)+somyp(ic)
            end do
            !effectif total pour le caractere
            efft=efft+effp

            if (effp == 0.d0) then
                call stop_application('Father ['//trim(data%genea%pere(ip))//'] has got no child with trait value ['&
                    //trim(data%phenoModel%carac(ic))//"]")
            end if

            !moyenne du caractere dans la famille du pere ip
            do ic=1,data%phenoModel%ncar
                dataset%lSires(ip)%xmu0(ic)=somyp(ic)/dble(effp)
                call log_mess("Mean for Sire/trait ["//trim(data%genea%pere(ip))//","//trim(data%phenoModel%carac(ic))&
                    //"]="//str(dataset%lSires(ip)%xmu0(ic)),INFO_DEF)
            end do
            !ecart type du caractere dans la famille du pere ip
            dataset%lSires(ip)%sig0(:)=0.d0
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                        do ic=1,data%phenoModel%ncar
                            dataset%lSires(ip)%sig0(ic)=dataset%lSires(ip)%sig0(ic)+(data%phenoAnimal%y(ic,kd)-&
                                dataset%lSires(ip)%xmu0(ic))*(data%phenoAnimal%y(ic,kd)-dataset%lSires(ip)%xmu0(ic))
                        end do
                    end if
                end do !kd
            end do !jm

            do ic=1,data%phenoModel%ncar
                dataset%lSires(ip)%sig0(ic)=sqrt(dataset%lSires(ip)%sig0(ic)/dble(effp-1))
                call log_mess("Standart deviation for Sire/trait ["//&
                 trim(data%genea%pere(ip))//","//trim(data%phenoModel%carac(ic))&
                    //"]="//str(dataset%lSires(ip)%sig0(ic)),INFO_DEF)
            end do
        end do !ip

        !moyenne du caractere sur l ensemble de la population
        allocate(dataset%xmu(data%phenoModel%ncar))
        do ic=1,data%phenoModel%ncar
            dataset%xmu(ic)=somy(ic)/dble(efft)
            call log_mess("Mean for trait ["//trim(data%phenoModel%carac(ic))//"]="//str(dataset%xmu(ic)),INFO_DEF)
        end do

        allocate(dataset%sig(data%phenoModel%ncar))
        allocate(dataset%cov(data%phenoModel%ncar,data%phenoModel%ncar))
        allocate(dataset%rhoi(data%phenoModel%ncar,data%phenoModel%ncar))
        dataset%cov=0.d0
        dataset%sig=0.d0

        do ip=1,data%genea%np
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if (count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar) then
                        do ic=1,data%phenoModel%ncar
                            dataset%sig(ic)=dataset%sig(ic)+&
                                (data%phenoAnimal%y(ic,kd)-dataset%xmu(ic))*(data%phenoAnimal%y(ic,kd)-dataset%xmu(ic))
                            do jc=1,ic-1
                                dataset%cov(jc,ic)=dataset%cov(jc,ic)+&
                                    (data%phenoAnimal%y(jc,kd)-dataset%xmu(jc))*(data%phenoAnimal%y(ic,kd)-dataset%xmu(ic))
                            end do
                        end do
                    end if
                end do
            end do
        end do

        !variance sur l ensemble de la population
        do ic=1,data%phenoModel%ncar
            dataset%sig(ic)=sqrt(dataset%sig(ic)/dble(efft-1))
            call log_mess("Standart deviation for trait ["//trim(data%phenoModel%carac(ic))//"]="//str(dataset%sig(ic)),INFO_DEF)
        end do
        !calcul des covariances
        do ic=1,data%phenoModel%ncar
            do jc=1,ic-1
                dataset%cov(jc,ic)=dataset%cov(jc,ic)/dble(efft)
                dataset%cov(ic,jc)=dataset%cov(jc,ic)
                dataset%rhoi(jc,ic)=dataset%cov(jc,ic)/(dataset%sig(ic)*dataset%sig(jc))
                dataset%rhoi(ic,jc)=dataset%rhoi(jc,ic)
            end do
        end do

    !   stop
    end subroutine init_dataset_ncar
    !!***

    !!****f* m_qtlmap_incidence_multi/add_general_mean_multi
    !! NAME
    !!   add_general_mean_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine add_general_mean_multi(data,xinc,listDesc)
        type(QTLMAP_DATASET)       ,intent(in) :: data
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal) , intent(inout) :: xinc
        type(INCIDENCE_TYPE) , dimension(data%phenoModel%ncar) , intent(inout) :: listDesc
        integer :: ic,ntniv
        call log_mess("multi incidence : ** Add General mean ** ",DEBUG_DEF)

        do ic=1,data%phenoModel%ncar
            listDesc(ic)%ntniv = listDesc(ic)%ntniv+1
            ntniv = listDesc(ic)%ntniv
            listDesc(ic)%nteff = listDesc(ic)%nteff+1
            listDesc(ic)%desc(listDesc(ic)%nteff)%name="General Mean"
            listDesc(ic)%desc(listDesc(ic)%nteff)%start=ntniv
            listDesc(ic)%desc(listDesc(ic)%nteff)%end=ntniv
            listDesc(ic)%desc(listDesc(ic)%nteff)%haveSubDesc=.true.
            allocate (listDesc(ic)%desc(listDesc(ic)%nteff)%listSubDesc(1))
            listDesc(ic)%desc(listDesc(ic)%nteff)%listSubDesc(1)%name="General Mean"!//"["//trim(carac(ic))//"]"
            listDesc(ic)%desc(listDesc(ic)%nteff)%listSubDesc(1)%start=ntniv
            listDesc(ic)%desc(listDesc(ic)%nteff)%listSubDesc(1)%end=ntniv
            xinc(ic,:,ntniv) = 1
        end do
       !incidenceDesc%nivdir(:,incidenceDesc%nteff)   = 1
    end subroutine add_general_mean_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/add_polygenic_effect_multi
    !! NAME
    !!   add_polygenic_effect_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine add_polygenic_effect_multi(data,xinc,listDesc)
        type(QTLMAP_DATASET)       ,intent(in) :: data
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal) , intent(inout) :: xinc
        type(INCIDENCE_TYPE) , dimension(data%phenoModel%ncar) , intent(inout) :: listDesc
        integer :: nkd
        integer :: ip,jm,kd,nt,j,fem,nteff,ic
        call log_mess("incidence : ** Multi : Add Polygenic effect to estimate **",DEBUG_DEF)

        do ic=1,data%phenoModel%ncar
            nt=listDesc(ic)%ntniv
            nteff = listDesc(ic)%nteff
            listDesc(ic)%desc(listDesc(ic)%nteff+1)%name="Sire polygenic effects"
            listDesc(ic)%desc(listDesc(ic)%nteff+1)%start=nt+1
            listDesc(ic)%desc(listDesc(ic)%nteff+1)%end=nt+data%genea%np
            listDesc(ic)%desc(listDesc(ic)%nteff+1)%haveSubDesc=.true.
            allocate (listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(data%genea%np))
            do ip=1,data%genea%np
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(ip)%name="Sire "//trim(data%genea%pere(ip))!//"["//trim(carac(ic))//"]"
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(ip)%start=nt+ip
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(ip)%end=nt+ip
            end do
            listDesc(ic)%nteff = listDesc(ic)%nteff+1
            listDesc(ic)%ntniv = listDesc(ic)%ntniv+data%genea%np

            if ( count(estime_multi(:))>0) then
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%name="Dam polygenic effects"
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%start=listDesc(ic)%desc(listDesc(ic)%nteff)%end+1
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%end=nt+data%genea%np+count(estime_multi(:))!nbfemLoc(ic)
                listDesc(ic)%desc(listDesc(ic)%nteff+1)%haveSubDesc=.true.
                allocate (listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(count(estime_multi(:))))
                fem=0
                do ip=1,data%genea%np
                    do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                        if ( estime_multi(jm) ) then
                            fem=fem+1
                            listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(fem)%name=&
                                "Dam "//trim(data%genea%mere(jm))//" [Sire "//trim(data%genea%pere(ip))//"]"!&//"["//trim(carac(ic))//"]"
                            listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(fem)%start=nt+data%genea%np+fem
                            listDesc(ic)%desc(listDesc(ic)%nteff+1)%listSubDesc(fem)%end=nt+data%genea%np+fem
                        end if
                    end do
                end do
                listDesc(ic)%nteff = listDesc(ic)%nteff+1
                listDesc(ic)%ntniv = listDesc(ic)%ntniv+count(estime_multi(:))!nbfemLoc(ic)
            else

            end if
            j=0

            !initialisation
            xinc(ic,:,nt+1:nt+data%genea%np+count(estime_multi(:data%genea%nm))) = 0

            nkd=0
            do ip=1,data%genea%np
                xinc(ic,:,nt+ip) = 0
                do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                    do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                        if ( count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar ) then
                            nkd =nkd +1
                            xinc(ic,nkd,nt+ip) = 1
                            !    listDesc(ic)%nivdir(kd,nteff+1)= nt+ip
                            if (estime_multi(jm)) then
                                fem=count(estime_multi(:jm)) ! on compte le nombre estimable dans la famille de pere ip
                                xinc(ic,nkd,nt+data%genea%np+fem) = 1
                            !    listDesc(ic)%nivdir(kd,nteff+2)= nt+np+fem
                            end if
                        end if
                    end do
                end do
            end do
        end do !ic
    end subroutine add_polygenic_effect_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/add_qtleffect_multi
    !! NAME
    !!   add_qtleffect_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine add_qtleffect_multi(data,spt,nteff,numqtl,chr,n,xinc,listDesc,mint,linear)
        type(QTLMAP_DATASET)                        ,intent(in)      :: data
        type(PDD_BUILD)                             ,intent(in)      :: spt
        integer                                      ,intent(in)     :: nteff,numqtl,chr,n
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal) ,intent(inout)  :: xinc
        type(INCIDENCE_TYPE) , dimension(data%phenoModel%ncar)      , intent(inout)  :: listDesc
        integer                                     , intent(in)     :: mint ! index of interaction to remove from incidence
        logical                                     , intent(in)     :: linear

        integer :: ic,ntniv,nteffstart,ip,fem,jm,l,s,u,nbint,nbt,jef
        logical ,dimension(data%genea%nfem) :: femInit

        call log_mess("multi incidence : ** Add Qtl effect to estimate **",DEBUG_DEF)

        do ic=1,data%phenoModel%ncar

            nteffstart = listDesc(ic)%nteff
            listDesc(ic)%nteff = listDesc(ic)%nteff+1
            if ( count(estime_multi(:)) > 0 ) listDesc(ic)%nteff = listDesc(ic)%nteff+1

            nbint=data%phenoModel%modele(ic,3)
            nbt=3+data%phenoModel%modele(ic,1)+data%phenoModel%modele(ic,2)
            listDesc(ic)%ntlev=1
            do jef=1,nbint
                if ( mint == jef) cycle ! this interaction is not adding
                listDesc(ic)%ntlev=listDesc(ic)%ntlev*data%phenoModel%nlev(data%phenoModel%modele(ic,nbt+jef))
            end do

            listDesc(ic)%desc(nteffstart+1)%name="Sire QTL effects ["//trim(str(numqtl))//"]"
            listDesc(ic)%desc(nteffstart+1)%start=listDesc(ic)%ntniv+1
            listDesc(ic)%ntniv_qtlsires(numqtl)=listDesc(ic)%desc(nteffstart+1)%start
            listDesc(ic)%desc(nteffstart+1)%end=listDesc(ic)%ntniv+listDesc(ic)%ntlev*data%genea%np

            listDesc(ic)%desc(nteffstart+1)%haveSubDesc=.true.
            allocate (listDesc(ic)%desc(nteffstart+1)%listSubDesc(data%genea%np*listDesc(ic)%ntlev))
            u=0
            do ip=1,data%genea%np
                l=0
                do s=listDesc(ic)%ntniv+(ip-1)*listDesc(ic)%ntlev+1,listDesc(ic)%ntniv+ip*listDesc(ic)%ntlev
                    l=l+1
                    u=u+1
                    listDesc(ic)%desc(nteffstart+1)%listSubDesc(u)%name="Sire "//trim(data%genea%pere(ip))//" Level "//trim(str(l))!&//"["//trim(carac(ic))//"]"
                    listDesc(ic)%desc(nteffstart+1)%listSubDesc(u)%start=s
                    listDesc(ic)%desc(nteffstart+1)%listSubDesc(u)%end=s
                end do
            end do

            if ( count(estime_multi(:)) > 0 ) then
                listDesc(ic)%desc(nteffstart+2)%name="Dam QTL effects ["//trim(str(numqtl))//"]"
                listDesc(ic)%desc(nteffstart+2)%start=listDesc(ic)%desc(nteffstart+1)%end+1
                listDesc(ic)%ntniv_qtldams(numqtl)=listDesc(ic)%desc(nteffstart+2)%start
                femInit=.false.
                l=0
                do jm=1,data%genea%nm
                    if (estime_multi(jm) ) then
                        fem=iam_multi(data%genea%repfem(jm))
                        if (femInit(fem)) then
                            cycle
                        end if
                        femInit(fem)=.true.
                        l=l+1
                    end if
                end do

                allocate (listDesc(ic)%desc(nteffstart+2)%listSubDesc(l*listDesc(ic)%ntlev))
                listDesc(ic)%desc(nteffstart+2)%end=listDesc(ic)%ntniv+listDesc(ic)%ntlev*(data%genea%np+l)
                listDesc(ic)%desc(nteffstart+2)%haveSubDesc=.true.
                !fem=0
                femInit=.false.
                u=0
                do jm=1,data%genea%nm
                    if (estime_multi(jm) ) then
                        fem=iam_multi(data%genea%repfem(jm))
                        if (femInit(fem)) then
                            cycle
                        end if
                        femInit(fem)=.true.
                        l=0
                        do s=listDesc(ic)%desc(nteffstart+2)%start+(fem-1)*listDesc(ic)%ntlev,&
                            listDesc(ic)%desc(nteffstart+2)%start-1+fem*listDesc(ic)%ntlev
                            l=l+1
                            u=u+1
                            !fem=fem+1
                            listDesc(ic)%desc(nteffstart+2)%listSubDesc(u)%name=&
                             "Dam "//trim(data%genea%mere(jm))//" Level "//trim(str(l))!&//"["//trim(carac(ic))//"]"
                            listDesc(ic)%desc(nteffstart+2)%listSubDesc(u)%start=s
                            listDesc(ic)%desc(nteffstart+2)%listSubDesc(u)%end=s
                        end do
                    end if
                end do

                listDesc(ic)%ntniv = listDesc(ic)%desc(nteffstart+2)%end
            else
                listDesc(ic)%ntniv = listDesc(ic)%desc(nteffstart+1)%end
            end if
        end do !ic

        call change_qtleffect_multi(data,spt,nteff,numqtl,chr,n,xinc,listDesc,mint)

    end subroutine add_qtleffect_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/change_qtleffect_multi
    !! NAME
    !!   change_qtleffect_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine change_qtleffect_multi(data,spt,nteff,iq,chr,n,xinc,listDesc,mint)
        type(QTLMAP_DATASET)                    ,intent(in)     :: data
        type(PDD_BUILD)                         ,intent(in)     :: spt
        integer                                 ,intent(in)     :: nteff,iq,chr,n
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal) ,intent(inout)  :: xinc
        type(INCIDENCE_TYPE) , dimension(data%phenoModel%ncar)      , intent(inout)  :: listDesc
        integer                                 ,intent(in)     :: mint ! index of interaction to remove from incidence

        !local
        integer :: nbef,nbco,nbint,nbt,jef,nivint,ntlev,ip,jm,kd,ntt,km,ntniv,nkd,ic
        real (kind=dp) ,dimension(:,:),allocatable:: pp,pm

        allocate (pp(data%genea%nd,data%phenoModel%ncar))
        allocate (pm(data%genea%nd,data%phenoModel%ncar))

        nkd=0

        do ip=1,data%genea%np
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                !if ( linear ) then
                 !call getprob_lin(chr,n,ic,ip,jm,ndm(jm)+1,ndm(jm+1),pp,pm)
                !else
                do ic=1,data%phenoModel%ncar
                    call getglobalprob(data,spt,chr,n,ic,ip,jm,data%genea%ndm(jm)+1,&
                        data%genea%ndm(jm+1),pp(:,ic),pm(:,ic))
                end do
                !end if
                do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)
                    if ( count(data%phenoAnimal%presentc(:,kd)) == data%phenoModel%ncar ) then
                        nkd=nkd+1

                        !  si des effets fixes sont en interaction avec l'allele paternel au qtl
                        !  il faut estimer les effets qtl pour chacun des niveaux concernes
                        !
                        !  on commence donc par creer un effet composites regroupant l'ensemble de
                        !  ces effets
                        !
                        do ic=1,data%phenoModel%ncar
                            ntniv = listDesc(ic)%desc(nteff)%start -1
                            nivint=1
                            ntlev=1
                            nbt=3+data%phenoModel%modele(ic,1)+data%phenoModel%modele(ic,2)
                            !nt=ntniv
                            do jef=1,data%phenoModel%modele(ic,3)
                                if ( mint == jef) cycle ! this interaction is not adding
                                nivint=nivint+ntlev*(data%phenoAnimal%nivx(kd,data%phenoModel%modele(ic,nbt+jef))-1)
                                ntlev=ntlev*data%phenoModel%nlev(data%phenoModel%modele(ic,nbt+jef))
                            end do
                            ntt=ntniv+nivint+ntlev*(ip-1)
                            xinc(ic,nkd,ntt)=pp(kd,ic)

                            ! incidenceDesc%nivdir(kd,nteff+1)=ntt
                            ! print *,'nt:',nivdir(kd,nteff+1),' vp:',pp(kd),' nteff:',nteff,' kd:',kd

                            if(estime_multi(jm)) then
                                km=iam_multi(data%genea%repfem(jm))
                                ntt= ntniv+ntlev*data%genea%np+nivint+ntlev*(km-1)
                                xinc(ic,nkd,ntt)=pm(kd,ic)
                            !   incidenceDesc%nivdir(kd,nteff+2)=ntt
                             ! print *,'nt:',nivdir(kd,nteff+2),' vp:',pm(kd)
                            end if
                        end do !ic
                    end if !presentc
                end do
            end do
        end do


        deallocate (pp)
        deallocate (pm)
        ! a modifier
        ! --> on devrait integerer l intialisation des pdds dans la boucle precedente
        ! attention les interaction qtl ne sont pas encore pris en compte
        ! if (.not. linear) then
        call pdd_at_position_multi(data,spt,listDesc(1)%dataset,iq,chr,n)
     !  end if
    end subroutine change_qtleffect_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/pdd_at_position_multi
    !! NAME
    !!   pdd_at_position_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine pdd_at_position_multi(data,spt,dataset,iq,chr,n)
        type(QTLMAP_DATASET)       ,intent(in) :: data
        type(PDD_BUILD)            ,intent(in) :: spt
        type(DATASET_TYPE)             , intent(inout)  :: dataset
        integer   ,intent(in)   :: iq,chr,n

        integer :: ip,jm,kd,ig,kkd,igg,ifem,kkkd,kds
        real (kind=dp) :: ppt,pmt

        do ip=1,data%genea%np
            dataset%lSires(ip)%ppt(iq,:)=0.d0
            ifem=0
            kkkd=0
            do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                ifem=ifem+1
                if (dataset%lSires(ip)%full_sib(ifem)%firstkd>=0) then
                    dataset%lSires(ip)%full_sib(ifem)%ppt(iq,:,:)=0.d0
                    dataset%lSires(ip)%full_sib(ifem)%pmt(iq,:,:)=0.d0
                    igg=0
                    kds=kkkd
                    ! print *,"ip,jm,firstkd:",ip,jm,lSires(ip)%half_sib%firstKd
                    do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
                        igg=igg+1
                        kkd=0
                        kkkd=kds
                        do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
                            if(count(data%phenoAnimal%presentc(:,spt%ndesc(chr,kd)))==data%phenoModel%ncar) then
                                kkd=kkd+1
                                kkkd=kkkd+1
                                if( estime_multi(jm) )then
                                    dataset%lSires(ip)%full_sib(ifem)%ppt(iq,igg,kkd)=&
                                        -spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                                    dataset%lSires(ip)%full_sib(ifem)%pmt(iq,igg,kkd)=&
                                        -spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                                else
                                    dataset%lSires(ip)%ppt(iq,dataset%lSires(ip)%half_sib%firstKd+kkkd-1)=&
                                        -spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n)
                                end if
                            end if
                        end do !! kd
                    end do !! ig
                end if
            end do ! jm
        end do ! ip
    end subroutine pdd_at_position_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/add_effcov_multi
    !! NAME
    !!   add_effcov_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine add_effcov_multi(data,xinc,listDesc,meff,mcov)
        type(QTLMAP_DATASET)                                  ,intent(in) :: data
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%nd,ntnivmaxtotal) , intent(inout) :: xinc
        type(INCIDENCE_TYPE) , dimension(data%phenoModel%ncar)      , intent(inout)  :: listDesc
        integer                                 , intent(in)    :: meff   ! indice on effet fixed list to remove (<=0 otherwise)
        integer                                 , intent(in)    :: mcov   ! indice on covariate list to remove (<=0 otherwise)
        integer :: ic,ntniv,ip,jm,kd,ntnifix,nbt,nbef,nbco,nteff,jef,jcov,ieff,nbniv,nivd,i,j,nbef2,nbco2,nkd

        call log_mess("incidence : ** Add eff and cov ** ",DEBUG_DEF)

        do ic=1,data%phenoModel%ncar

            nbef=data%phenoModel%modele(ic,1)
            nbco=data%phenoModel%modele(ic,2)

            nbef2=nbef
            if ( meff /= 0 ) nbef2=nbef-1
            nbco2=nbco
            if ( mcov /=0 ) nbco2=nbco-1


            if (nbef2+nbco2 <=0)return

            nbniv=0
            ieff=listDesc(ic)%nteff

            if ( nbef2 > 0 ) then
                ieff=ieff+1
                listDesc(ic)%eqtl_print(ieff)=.false.
                do i=1,nbef
                    if ( i == meff ) cycle
                    nbniv=data%phenoModel%nlev(data%phenoModel%modele(ic,3+i))+nbniv
                end do
                listDesc(ic)%desc(ieff)%name="Fixed effects"
                listDesc(ic)%desc(ieff)%start=listDesc(ic)%ntniv+1
                listDesc(ic)%desc(ieff)%end=listDesc(ic)%ntniv+nbniv
                listDesc(ic)%desc(ieff)%haveSubDesc=.true.

                !init of xinc
                xinc(ic,:,listDesc(ic)%desc(ieff)%start:listDesc(ic)%desc(ieff)%end)=0.d0

                allocate ( listDesc(ic)%desc(ieff)%listSubDesc(nbniv) )

                j=0
                do jef=1,nbef
                    if ( jef == meff ) cycle
                    do i=1,data%phenoModel%nlev(data%phenoModel%modele(ic,3+jef))
                        j=j+1
                        listDesc(ic)%desc(ieff)%listSubDesc(j)%name=trim(data%phenoModel%namefix(data%phenoModel%modele(ic,3+jef)))&
                            //' level '//trim(str(i))!&//"["//trim(carac(ic))//"]"
                        listDesc(ic)%desc(ieff)%listSubDesc(j)%start=listDesc(ic)%ntniv+j
                        listDesc(ic)%desc(ieff)%listSubDesc(j)%end=listDesc(ic)%ntniv+j
                    end do
                end do
                listDesc(ic)%nteff = listDesc(ic)%nteff+1
            end if

            if ( nbco2 > 0 ) then
                ieff=ieff+1
                listDesc(ic)%eqtl_print(ieff)=.false.
                listDesc(ic)%desc(ieff)%name="Covariates"
                listDesc(ic)%desc(ieff)%start=listDesc(ic)%ntniv+nbniv+1
                listDesc(ic)%desc(ieff)%end=listDesc(ic)%ntniv+nbniv+nbco2
                listDesc(ic)%desc(ieff)%haveSubDesc=.true.

                !init of xinc
                xinc(ic,:,listDesc(ic)%desc(ieff)%start:listDesc(ic)%desc(ieff)%end)=0.d0

                allocate ( listDesc(ic)%desc(ieff)%listSubDesc(nbco2) )

                j=0
                do jcov=1,nbco
                    if ( jcov == mcov ) cycle
                    j=j+1
                    listDesc(ic)%desc(ieff)%listSubDesc(j)%name=&
                        trim(data%phenoModel%namecov(data%phenoModel%modele(ic,3+data%phenoModel%modele(ic,1)+jcov)))!//"["//trim(carac(ic))//"]"
                    listDesc(ic)%desc(ieff)%listSubDesc(j)%start=listDesc(ic)%ntniv+nbniv+j
                    listDesc(ic)%desc(ieff)%listSubDesc(j)%end=listDesc(ic)%ntniv+nbniv+j
                end do
                listDesc(ic)%nteff = listDesc(ic)%nteff+1
            end if

            nkd = 0
            nbniv=0
            do ip=1,data%genea%np
                do jm=data%genea%nmp(ip)+1,data%genea%nmp(ip+1)
                    do kd=data%genea%ndm(jm)+1,data%genea%ndm(jm+1)

                        if ( count(data%phenoAnimal%presentc(:,kd))==data%phenoModel%ncar ) then
                            nkd=nkd+1
                            ntniv=listDesc(ic)%ntniv

                            !
                            !  autres effets fixes
                            !
                            nbniv=0
                            do jef=1,nbef
                                if ( jef == meff ) cycle ! this effet fixed is not adding
                                ! listDesc(ic)%nivdir(kd,nteff+jef)=ntniv+nbniv+nivx(kd,modele(ic,3+jef))
                                nivd=ntniv+nbniv+data%phenoAnimal%nivx(kd,data%phenoModel%modele(ic,3+jef))
                                xinc(ic,nkd,nivd)=1
                                nbniv=nbniv+data%phenoModel%nlev(data%phenoModel%modele(ic,3+jef))
                            end do

                            ntnifix=ntniv+nbniv
                            nbt=3+nbef
                            !covariate
                            j=0
                            do jcov=1,nbco
                                if ( jcov == mcov ) cycle ! this covariate is not adding
                                j=j+1
                                ! listDesc(ic)%covdir(kd,jcov)=covar(kd,modele(ic,nbt+jcov))
                                xinc(ic,nkd,ntnifix+j)=data%phenoAnimal%covar(kd,data%phenoModel%modele(ic,nbt+jcov))
                            end do
                          !ntniv=ntnifix+nbco
                          !nbt=nbt+nbco
                        end if
                    end do
                end do
            end do

            listDesc(ic)%ntniv=listDesc(ic)%ntniv+nbniv+nbco2
        end do !ic

        call log_mess("FIN ** Add eff and cov ** ",DEBUG_DEF)

    end subroutine add_effcov_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/confusion_multi_type1
    !! NAME
    !!   confusion_multi_type1
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine confusion_multi_type1(data,hypothesis,listDesc,xxx,workstruct)
        type(QTLMAP_DATASET)                          ,intent(in)      :: data
        integer                                       , intent(in)     :: hypothesis
        type(INCIDENCE_TYPE)   , dimension(data%phenoModel%ncar)      , intent(inout)  :: listDesc
        real (kind=dp)  ,dimension(data%phenoModel%ncar,ntnivmaxtotal,ntnivmaxtotal) ,intent(in)      :: xxx
        type(INCIDENCE_GEN_STRUCT)            ,intent(inout)       :: workstruct

        real (kind=dp)  ,dimension(:,:,:),allocatable    :: xcor
        type(CORR_ALERT_TYPE) ,dimension(:),allocatable  :: alerts
        type(CORR_ALERT_TYPE) ,dimension(:),allocatable  :: as
        integer :: ic,ialert,sizealert,eff,i,j,starteff

        allocate (xcor(data%phenoModel%ncar,ntnivmaxtotal,ntnivmaxtotal))

        sizealert=0
        do ic=1,data%phenoModel%ncar
            call get_correlations(listDesc(ic),xxx(ic,:,:),xcor(ic,:,:))
            sizealert=sizealert+data%genea%np*namest_multi*listDesc(ic)%ntlev*listDesc(ic)%nbnivest
        end do

        workstruct%corrmax(hypothesis)=0.d0

        allocate (alerts(sizealert))
        allocate (as(sizealert*workstruct%nqtl))

        workstruct%nalert(hypothesis)=0
        starteff = workstruct%listnteff(workstruct%nqtl)+1
        if (namest_multi>0) starteff = starteff+1

        do ic=1,data%phenoModel%ncar
            do eff=starteff,listDesc(ic)%nteff
                call log_mess("computing confusion beetween qtls and:"//listDesc(ic)%desc(eff)%name,VERBOSE_DEF)
                call confusion_between_qtl_and_effect(listDesc(ic),xcor(ic,:,:),&
                    eff,sizealert,alerts,ialert,workstruct%corrmax(hypothesis))
                if (ialert>0) then
                    !       print *,"ialert:",alerts(1)%corr,alerts(2)%corr,alerts(3)%corr
                    as(workstruct%nalert(hypothesis)+1:workstruct%nalert(hypothesis)+ialert)=alerts(:ialert)
                    workstruct%nalert(hypothesis)=workstruct%nalert(hypothesis)+ialert
                end if
            end do
        end do

        workstruct%alertQtl(hypothesis,:workstruct%nalert(hypothesis))=as(:workstruct%nalert(hypothesis))

        deallocate (alerts)
        deallocate (as)
        deallocate (xcor)


    end subroutine confusion_multi_type1
    !!***

    !!****f* m_qtlmap_incidence_multi/opti_0qtl_multi
    !! NAME
    !!   opti_0qtl_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine opti_0qtl_multi(  data,         &
        incsol   ,       &  !
        workstruct, &
        maxqtl,FUNCT_MODEL,resol_lu)          ! maximum qtl to find

        type(QTLMAP_DATASET)                    ,intent(inout)    :: data
        integer                                 ,intent(in)       :: maxqtl
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)    :: workstruct
        type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)      :: incsol
        logical                                 ,intent(in)       :: resol_lu
        external :: FUNCT_MODEL   ! function with a model specific

        integer :: nkd,ip,i,listnteff(1)

        real (kind=dp) ,dimension(:,:,:),pointer      :: xinc
        real (kind=dp) , dimension(:,:) ,pointer      :: Bestim,rhoiestim
        type(INCIDENCE_TYPE)  ,dimension(data%phenoModel%ncar)        :: listDesc
        type(POSITION_LRT_INCIDENCE)                  :: curPos
        real (kind=dp) :: sigsquareEstime(data%phenoModel%ncar,data%genea%np),f
        integer :: k,ic,jc
        type(PDD_BUILD)           , pointer :: spt => null()

        p_dpm => data%phenoModel
        p_dg => data%genea
        p_dpa => data%phenoAnimal
        p_dataTrio => data%datasetUser

        !initialisation of incidence matrix
        call init_incidence_multi(data,spt,0,xinc,listDesc,workstruct%effects,resol_lu)
        dataset=>listDesc(1)%dataset
        allocate (Bestim(p_dpm%ncar,ntnivmaxtotal))
        allocate (rhoiestim(p_dpm%ncar,p_dpm%ncar))
        Bestim=0.d0
        rhoiestim=0.d0
        sigsquareEstime=0.d0
        workstruct%nqtl=0
        Bestim=0.d0
        call add_general_mean_multi(data,xinc,listDesc)
        workstruct%effects%general=1
        !add polygenic effect
        call add_polygenic_effect_multi(data,xinc,listDesc)
        !add fixed effect and covariate
        call add_effcov_multi(data,xinc,listDesc,0,0)

        !      call model_optim_multi_h0(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,.true.)
        call FUNCT_MODEL(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,.true.)
        call set_solution_multi(data,0,workstruct,sigsquareEstime,rhoiestim,Bestim,listDesc,incsol,1,listnteff)

        if (size(workstruct%sigsquare,1)>=1) then
            do ip=1,p_dg%np
                do ic=1,p_dpm%ncar
                    workstruct%sigsquare(1,ip,ic)=sigsquareEstime(ic,ip)
                end do
            end do
        end if

        if (size(workstruct%fnqtlsires,1)>=1) then
            workstruct%fnqtlsires(1,:)=listDesc(1)%fperemax(:)
        end if

        if (size(workstruct%fnqtldams,1)>=1) then
            workstruct%fnqtldams(1,:)=listDesc(1)%fmeremax(:)
        end if

        if (associated(workstruct%rhoi)) then
            workstruct%rhoi(1,:,:) = rhoiestim
        end if

        !call model analysis
        !call FUNCT_MODEL(xinc,incidenceDesc,workstruct,Bestim,.true.)
        ! call set_solution(ic,0,workstruct,Bestim,incidenceDesc,incsol,1,listnteff)

        call end_incidence_multi(data,listDesc)
        deallocate (xinc)
        deallocate (Bestim)
        deallocate (rhoiestim)
    end subroutine opti_0qtl_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/init_startpoint
    !! NAME
    !!   init_startpoint
    !! DESCRIPTION
    !!
    !! NOTE
    !!  initialise le point de depart pour la minimisation
    !! SOURCE
    subroutine init_startpoint(workStruct,incDesc)
        type(INCIDENCE_TYPE)              ,intent(inout)       :: incDesc
        type(INCIDENCE_GEN_STRUCT)        ,intent(in)          :: workstruct

        real(kind=dp) :: r1,tempLin(p_dpm%ncar*p_dpm%ncar)
        integer :: k,ic,jc,ip,i

        !initialisation d'un point de depart pour les variances
        k=0
        do ic=1,p_dpm%ncar
            do ip=1,p_dg%np
                k=k+1
                incDesc%par(k)=sqrt(workstruct%sigsquare(workstruct%nqtl,ip,ic))
            end do
        end do

        i=0
        do jc=1,p_dpm%ncar-1
            do ic=jc+1,p_dpm%ncar
                i = i + 1
                tempLin(i) = workstruct%rhoi(workstruct%nqtl,ic,jc)
            end do
        end do

        i=0
        !initialisation d'un point de depart pour les correlations residuelles
        do ic=2,p_dpm%ncar
            do jc=1,ic-1
                k=k+1
                i=i+1
                ! *** RESCALE OFI ***
                r1=tempLin(i)
                incDesc%par(k)=dlog((1.d0+r1)/(1.d0-r1))
             !incDesc%par(k)=workstruct%rhoi(workstruct%nqtl,jc,ic)
            end do
        end do





        incDesc%par(p_dpm%ncar*p_dg%np+1+p_dpm%ncar*(p_dpm%ncar-1)/2:)=0.d0

    end subroutine init_startpoint
    !!***

    !!****f* m_qtlmap_incidence_multi/opti_nqtl_multi
    !! NAME
    !!   opti_nqtl_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE

    subroutine opti_nqtl_multi(     data,      &
        spt,       &
        nqtl,      &  ! number of qtl
        hyp,      &  ! under hypothesis hyp
        workStruct,       &  ! variance estimated under NQTL-1 hypothesis
        incsol   ,       &  ! incidence solution (estimation of each effect)
        lrtsol   ,       &  ! maximum lrt finding at a position
        FUNCT_MODEL,       &  ! function with a model specific
        maxqtl,resol_lu)
        type(QTLMAP_DATASET)                    ,intent(inout)       :: data
        type(PDD_BUILD)   ,target               ,intent(inout)       :: spt
        integer                                 ,intent(in)          :: nqtl,hyp,maxqtl
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)       :: workstruct
        type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)         :: incsol
        type(TYPE_LRT_SOLUTION)                 ,intent(out)         :: lrtsol
        logical                                 ,intent(in)       :: resol_lu
        external                                                     :: FUNCT_MODEL


        integer :: i

        real (kind=dp) ,dimension(:,:,:),pointer      :: xxx,xinc   ! incidence matrix and temp for testing nuisance effect
        type(INCIDENCE_TYPE)  ,dimension(data%phenoModel%ncar)        :: listDesc ! description incidence matrix
        type(POSITION_LRT_INCIDENCE)                  :: curPos        ! the position with additional info for the model
        integer :: n,j,chr,k,ic,jc,ip

         ! solution of the estimation
        real (kind=dp) , dimension(:,:) ,pointer      :: Bestim,rhoiestim
        real(kind=dp)                                 :: sigsquareEstime(data%phenoModel%ncar,data%genea%np)

        p_dpm => data%phenoModel
        p_dg => data%genea
        p_dpa => data%phenoAnimal
        p_spt => spt
        p_dataTrio => data%datasetUser

        if ( nqtl <= 0 ) then
            call stop_application("Error dev: can not call opti_nqtl_linear with nqtl:"//trim(str(nqtl)))
        end if

        !Position Allocation
        call init_position(data,hyp,nqtl,curPos)

        !initialisation of incidence matrix
        call init_incidence_multi(data,spt,nqtl,xinc,listDesc,workstruct%effects,resol_lu)

        dataset=>listDesc(1)%dataset

        allocate (Bestim(p_dpm%ncar,ntnivmaxtotal))
        allocate (rhoiestim(p_dpm%ncar,p_dpm%ncar))
        allocate (xxx(p_dpm%ncar,ntnivmaxtotal,ntnivmaxtotal))

        if ( .not. resol_lu ) call init_startpoint(workstruct,listDesc(1))

        Bestim=0.d0
        rhoiestim=0.d0
        sigsquareEstime=0.d0

        do i=1,nqtl
            n=0
            j=i
            do chr=1,data%map%nchr
                n=n+data%map%get_ilong(chr)
                if (i<=n) exit
                j=1
            end do

            if (chr>data%map%nchr) call stop_application("Not enough sampling point to detect QTLs")

            curPos%listChr(i)=1 ! attention si pas assez d echantillonage on doit passer au chromosome suivant....
            curPos%listN(i)=i
        end do

        !      call build_incidence_matrix(workstruct,curPos,xinc,incidenceDesc,0,0,0)
        Bestim=0.d0
        call add_general_mean_multi(data,xinc,listDesc)
        workstruct%effects%general=1

        workstruct%listnteff(1)=2
        do i=1,workstruct%nqtl
            !add qtl effect/interaction at position n to estim
            call add_qtleffect_multi(data,spt,workstruct%listnteff(i),i,curPos%listChr(i),curPos%listN(i),&
                xinc,listDesc,0,workstruct%linear)
            if (i<workstruct%nqtl) workstruct%listnteff(i+1)=listDesc(1)%nteff+1
        end do
        !add polygenic effect
        call add_polygenic_effect_multi(data,xinc,listDesc)
        !add fixed effect and covariate
        call add_effcov_multi(data,xinc,listDesc,0,0)

        !initalizing lrtsolution
        call lrtsol%new(data,nqtl)

        lrtsol%lrtmax=-INIFINY_REAL_VALUE

        !  call gen_loop_opti_nqtl(1,nqtl,curPos,workstruct,xinc,&
        ! incidenceDesc,Bestim,lrtsol,.true.,FUNCT_MODEL)
        !call gen_loop_opti_nqtl_multi(1,nqtl,curPos,workstruct,&
        !         xinc,listDesc,sigsquareEstime,rhoiestim,Bestim,lrtsol,.false.,FUNCT_MODEL)


        call gen_opti_nqtl_multi(data,spt,nqtl,hyp,curPos,workstruct,xinc,&
            listDesc,sigsquareEstime,rhoiestim,Bestim,lrtsol,FUNCT_MODEL)

        ! Compute precision at the maximum LRT in position posx
        do i=1,nqtl
            call change_qtleffect_multi(data,spt,workstruct%listnteff(i),i,lrtsol%chrmax(i-1),lrtsol%nxmax(i-1),&
                xinc,listDesc,0)
            curPos%listChr(i)=lrtsol%chrmax(i-1)
            curPos%listN(i)=lrtsol%nxmax(i-1)
        end do

        !    call debug_write_incidence(xinc,incidenceDesc)
        call FUNCT_MODEL(xinc,listDesc,curPos,workstruct,sigsquareEstime,rhoiestim,Bestim,&
            .true.,workstruct%performConfusion,xxx,.false.)

        !    print *,"****************************************************************"
        if ( workstruct%performConfusion) then
            print *," ** devel message : no call to confusion function  ** "
        ! call confusion_multi_type1(hyp,listDesc,xxx,workstruct)
        end if
        if ( workstruct%performTestNuis ) then
        !    call test_nuisances(ic,incidenceDesc%nbnivest,lrtsol,curPos,workstruct,FUNCT_MODEL)
        end if
        call set_solution_multi(data,nqtl,workstruct,sigsquareEstime,&
            rhoiestim,Bestim,listDesc,incsol,nqtl,workstruct%listnteff)

        if (size(workstruct%sigsquare,1)>nqtl) then
            do ip=1,p_dg%np
                do ic=1,p_dpm%ncar
                    workstruct%sigsquare(nqtl,ip,ic)=sigsquareEstime(ic,ip)
                end do
            end do
        end if

        if (size(workstruct%fnqtlsires,1)>nqtl) then
            workstruct%fnqtlsires(nqtl,:)=listDesc(1)%fperemax
        end if

        if (size(workstruct%fnqtldams,1)>nqtl) then
            workstruct%fnqtldams(nqtl,:)=listDesc(1)%fmeremax
        end if

        if (size(workstruct%rhoi,1)>nqtl) then
            workstruct%rhoi(nqtl,:,:) = rhoiestim
        end if

        call end_incidence_multi(data,listDesc)

        deallocate (xxx)
        deallocate (xinc)
        deallocate (Bestim)
        deallocate (rhoiestim)
        call end_position(curPos)
    end subroutine opti_nqtl_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/model_optim_multi_h0
    !! NAME
    !!   model_optim_multi_h0
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE

    subroutine model_optim_multi_h0(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,performPrecision)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)                   ,intent(inout)    :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)      , intent(inout)  :: listDesc
        logical                                    ,intent(in)         :: performPrecision

        real(kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) :: xxx


        call model_optim_multi_family(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
            performPrecision,likelihood_ncar_h0_family,likelihood_ncar_h0_family_withcd,.true.,xxx)

        !stop
       !workstruct%sigsquareEstime=osig*osig

    end subroutine model_optim_multi_h0
    !!***

    !!****f* m_qtlmap_incidence_multi/model_optim_multi_hn
    !! NAME
    !!   model_optim_multi_hn
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine model_optim_multi_hn(xinc,listDesc,curPos,workstruct,sigsquareEstime,rhoiestim,Bestim,&
        performPrecision,tConf,tempForConfusion,invlrt)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)                   ,intent(inout)    :: workstruct
        type(POSITION_LRT_INCIDENCE)                 ,intent(inout)    :: curPos
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)      , intent(inout)  :: listDesc
        logical                                    ,intent(in)         :: performPrecision,tConf
        real (kind=dp),dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) ,intent(out)  :: tempForConfusion
        logical                                          ,intent(in)   :: invlrt

        integer :: i,hypothesis
        character(len=LEN_L) :: dx

        allocate (current_chr(workstruct%nqtl))
        current_chr=curPos%listChr

        call model_optim_multi_family(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
            performPrecision,likelihood_ncar_hn_family,likelihood_ncar_hn_family_withcd,tConf,tempForConfusion)
        !         print *,"F:",workstruct%fmax


        if ((any(listDesc(1)%fperemax == INIFINY_REAL_VALUE)) ) then
            dx=""
            do i=1,workstruct%nqtl-1
                dx=trim(dx)//trim(str(dataset%data%map%absi(curPos%listChr(i),curPos%listN(i))))//","
            end do

            dx=trim(dx)//str(dataset%data%map%absi(curPos%listChr(workstruct%nqtl),curPos%listN(workstruct%nqtl)))
            call log_mess("dx ["//trim(dx)// "]. Can not optimize likelihood....The start point is reinitializing",WARNING_DEF)
            call init_startpoint(workStruct,listDesc(1))
            curPos%lrtSires=0.d0
            listDesc(1)%fperemax=0.d0
            listDesc(1)%fmeremax=0.d0
            return
        end if

        !compute LRT
        if (invLrt) then
            do hypothesis=1,workstruct%nqtl
                curPos%lrtSires(hypothesis,:)=-2.d0*(workstruct%fnqtlsires(hypothesis,:)-listDesc(1)%fperemax)
                curPos%lrtDams(hypothesis,:)=-2.d0*(workstruct%fnqtldams(hypothesis,:)-listDesc(1)%fmeremax)
            end do
        else
            do hypothesis=1,workstruct%nqtl
                curPos%lrtSires(hypothesis,:)=-2.d0*(listDesc(1)%fperemax-workstruct%fnqtlsires(hypothesis,:))
                curPos%lrtDams(hypothesis,:)=-2.d0*(listDesc(1)%fmeremax-workstruct%fnqtldams(hypothesis,:))
            end do
        end if

        deallocate (current_chr)

    end subroutine model_optim_multi_hn
    !!***

    !!****f* m_qtlmap_incidence_multi/model_optim_multi_family
    !! NAME
    !!   model_optim_multi_family
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine model_optim_multi_family(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
        performPrecision,FUNCT_PART,FUNCT_PART_CENSURE,tConf,tempForConfusion)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)         :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)   , intent(inout),target  :: listDesc
        logical                                    ,intent(in)         :: performPrecision,tConf
        real (kind=dp),dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) ,intent(out)  :: tempForConfusion
        external                                                       :: FUNCT_PART,FUNCT_PART_CENSURE

        integer :: j,i,ip,ifail,ibound,npar
        !        real(kind=dp) ,dimension(np+ntnivmax) :: par
        real (kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal)   :: XX
        real (kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal)   :: triang
        integer                                    :: iuser(1),ix,kd1,kd2,jm,jjm,ic,s,jc,iic,np,ncar,nm
        real (kind=dp)   ,dimension(1)    :: user
        logical       ,dimension(p_dpm%ncar,ntnivmaxtotal)    :: lastvecsol
        real(kind=dp) :: vci_t(p_dpm%ncar,p_dpm%ncar,p_dg%np)
        real(kind=dp) :: determ_t(p_dg%np),f
        real(kind=dp) :: r,tempLin(p_dpm%ncar*p_dpm%ncar)

        !Filter for the mimimization - dim (np,nm,nd)
        logical              , dimension(:,:,:) ,pointer      :: filter_inc
        !Filter for calcul on the vci matrix : dim np,
        logical              , dimension(:,:) , pointer       :: filter_vci
        logical              , dimension(:,:) , pointer       :: p_filter_vci

        np = p_dg%np
        nm = p_dg%nm
        ncar = p_dpm%ncar

        iuser=0
        user=0
        my_listDesc=>listDesc
        dataset => listDesc(1)%dataset
        allocate (my_xincreduitmul(p_dg%nd,ntnivmaxtotal,ncar))
        allocate (filter_inc(np,nm,ncar*np+ncar*(ncar-1)/2+ncar*ntnivmaxtotal))
        allocate (filter_vci(np,ncar*np+ncar*(ncar-1)/2+ncar*ntnivmaxtotal))


        !filter car initialisation
        filter_vci=.false.
        do ip=1,np
            do ic=1,ncar
                filter_vci(ip,(ic-1)*np+ip)=.true.
            end do
        end do

        filter_vci(:,ncar*np+1:ncar*np+ncar*(ncar-1)/2)=.true.


        my_xincreduitmul=0.d0
        do ic=1,ncar
            lastvecsol(ic,:my_listDesc(ic)%ntniv)=my_listDesc(ic)%vecsol(:my_listDesc(ic)%ntniv)
            ! create X'.X matrix from incidence matrix
            call model_XT_X(xinc(ic,:,:),my_listDesc(ic),XX)
            ! Check all parameters to remove from the estimation
            call estim_cholesky(XX,my_listDesc(ic),ntnivmaxtotal,triang)
            ! compute the precision of each parameter
            if (performPrecision) call get_precision(XX,tempForConfusion(:,:,ic),my_listDesc(ic))
            call set_corrxinc(xinc(ic,:,:),my_listDesc(ic),my_xincreduitmul(:,:,ic))
            !optimisation : on sauvegarde les index des elements non nul de la matrice reduite
            call fill_nonull_elements(my_listDesc(ic),size(my_xincreduitmul,1),size(my_xincreduitmul,2),my_xincreduitmul(:,:,ic))
          !call debug_write_incidence(xinc(ic,:,:),my_listDesc(ic))
        end do

          ! Optimisation de la vraisemblance a la position dx
        ifail=1
        ibound=0
        npar=ncar*np+ncar*(ncar-1)/2

        j=ncar*np+ncar*(ncar-1)/2
        do ic=1,ncar
            npar=npar+my_listDesc(ic)%nbnivest

            if (count(lastvecsol(ic,:))>0) then ! autre que initialisation
                do i=1,my_listDesc(ic)%ntniv
                    if(my_listDesc(ic)%vecsol(i))then
                        j=j+1
                        if(.not. lastvecsol(ic,i)) my_listDesc(ic)%par(j)=0.d0
                    end if
                end do
            else  ! initialisation des correlations phenotypiques
                s=ncar*np
                do iic=2,ncar
                    do jc=1,iic-1
                        s=s+1
                        my_listDesc(ic)%par(s)=log( (1+p_dpm%RhoP(iic,jc)) / (1-(p_dpm%RhoP(iic,jc))))
                    !        print *, my_listDesc(ic)%par(s)
                    end do
                end do
            end if
        end do
        !  allocate (filter_inc(np,npar))
        s=np*ncar+ncar*(ncar-1)/2
        filter_inc=.true.
        filter_inc(:,:,1:np*ncar)=.false.
        filter_inc(:,:,np*ncar+1:s)=.true. !rho
        do ic=1,ncar
            do ip=1,np
                filter_inc(ip,p_dg%nmp(ip)+1:p_dg%nmp(ip+1),(ic-1)*p_dg%np+ip)=.true.
                jjm=0
                do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
                    jjm=jjm+1
                    if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
                        kd1=dataset%lSires(ip)%full_sib(jjm)%firstkd
                        kd2=dataset%lSires(ip)%full_sib(jjm)%lastkd
                        filter_inc(ip,jm,s+1:s+my_listDesc(ic)%nbnivest)=&
                         any(my_xincreduitmul(kd1:kd2,:my_listDesc(ic)%nbnivest,ic)/=0.d0,dim=1)
                    else
                        !     print *,'pas de desc....'
                        filter_inc(ip,jm,s+1:s+my_listDesc(ic)%nbnivest)=.false.
                    end if
                end do
            end do
            s=s+my_listDesc(ic)%nbnivest
        end do

        if ( dataset%data%datasetUser%na>0 ) then
            !prise en compte des donnees censures
            call minimizing_funct_family(dataset%data,npar,ibound,FUNCT_PART_CENSURE,filter_inc,&
                listDesc(1)%fmeremax,listDesc(1)%fperemax,&
                my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,f,iuser,user,ifail)
        else
            p_filter_vci=>filter_vci
            !pas de donnee censure=> optimisation pour le calcul de l'inverse de matrice de covariance
            call minimizing_funct_family_multi(dataset%data,p_dpm,npar,ibound,FUNCT_PART,filter_inc,&
             p_filter_vci,vci_t,determ_t,listDesc(1)%fmeremax,listDesc(1)%fperemax,&
                my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,f,iuser,user,ifail)
        end if

        Bestim=0.d0
        s=np*ncar+ncar*(ncar-1)/2
        !getting standart deviation
        do ic=1,ncar
            do ip=1,np
                sigsquareEstime(ic,ip)=my_listDesc(1)%par((ic-1)*np+ip)*my_listDesc(1)%par((ic-1)*np+ip)
            end do
             !The solution
            Bestim(ic,:my_listDesc(ic)%nbnivest)=my_listDesc(1)%par(s+1:s+my_listDesc(ic)%nbnivest)
            s=s+my_listDesc(ic)%nbnivest
        end do

        rhoiestim=0.d0
        s=ncar*np
        do ic=2,ncar
            do jc=1,ic-1
                s=s+1
                ! *** RESCALE OFI ***
                !rhoiestim(jc,ic)=my_listDesc(1)%par(s)
                r = my_listDesc(1)%par(s)
                rhoiestim(jc,ic)=(dexp(r)-1.d0)/(dexp(r)+1.d0)
                rhoiestim(ic,jc)=rhoiestim(jc,ic)
            end do
        end do

        i=0
        do ic=2,ncar
            do jc=1,ic-1
                i=i+1
                tempLin(i)=rhoiestim(ic,jc)
            end do
        end do

        i=0
        do jc=1,ncar-1
            do ic=jc+1,ncar
                i = i + 1
                rhoiestim(ic,jc) = tempLin(i)
                rhoiestim(jc,ic) = rhoiestim(ic,jc)
            end do
        end do

        deallocate (my_xincreduitmul)
        deallocate (filter_inc)
        deallocate (filter_vci)
    end subroutine model_optim_multi_family
    !!***

    !!****f* m_qtlmap_incidence_multi/get_inv_residual_covariance_matrix
    !! NAME
    !!   get_inv_residual_covariance_matrix
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine get_inv_residual_covariance_matrix(ip,n,x,vci,determ,ok)
        integer         , intent(in)                   :: ip,n
        real (kind=dp)  ,dimension(n)   , intent(in)   :: x
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar),intent(out) :: VCI
        real(kind=dp),                     intent(out) :: determ
        logical                           ,intent(out) :: ok
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: rhoc
        real(kind=dp),dimension(p_dpm%ncar+1,p_dpm%ncar)    :: VCIInverse

        real(kind=dp) :: rh
        integer :: ic,jc,k,irank

        ok=.true.
        !RHO(I,J) avec i/=j => x((i-1)*(ncar-i+1))
        k=0
        do ic=2,p_dpm%ncar
            do jc=1,ic-1
                k=k+1
                ! *** RESCALE OFI ***
                !rhoc(jc,ic)=x(ncar*np+k)
                rh = x(p_dpm%ncar*p_dg%np+k)
                rhoc(jc,ic) = ((dexp(rh)-1.d0)/(dexp(rh)+1.d0))
                rhoc(ic,jc)=rhoc(jc,ic)
            end do
        end do

        !VCi : residual covariance matrix from sire I
        !         | SIGi1^2  SIGi1*SIGi2*RHO(1,2) ....  SIGi1*SIGip*RHO(1,p)
        !         | ..       SIGi2^2
        !    VCi  | ..                ...
        !         | ..                       ...
        !         | ..
        !         | SIGi1*SIGip*RHO(1,p)        SIGip^2

        do ic=1,p_dpm%ncar
            do jc=ic,p_dpm%ncar
                VCI(ic,jc)=x((ic-1)*p_dg%np+ip)*x((jc-1)*p_dg%np+ip)
                if ( ic /= jc ) then
                    VCI(ic,jc)=VCI(ic,jc)*rhoc(ic,jc)
                    VCI(jc,ic)=VCI(ic,jc)
                end if
            end do
        end do

        !Calcul Determinant de la matrice de covariance residuelle du pere ip
        !--------------------------------------------------------------------
        irank=0
        determ=0
        ! call MATH_QTLMAP_F03ABF(VCI,ncar,ncar,determ,irank)
        !      if (irank /= 0) then
        !        print *,"mauvais determinant matrice"
        !        ok=.false.
        !        return
        !      end if
        !Calcul Inverse de la matrice de covariance residuelle du pere ip
        !-----------------------------------------------------------------
        VCIInverse(:p_dpm%ncar,:p_dpm%ncar)=VCI
        CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,VCIInverse,p_dpm%ncar+1,determ,irank)
        !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
        if (irank /= 0) then
            print *,"mauvaise inversion matrice"
            ok=.false.
            return
        end if

        do ic=1,p_dpm%ncar
            do jc=1,ic
                VCI(ic,jc)=VCIInverse(ic+1,jc)
                VCI(jc,ic)=VCI(ic,jc)
            end do
        end do



    end subroutine get_inv_residual_covariance_matrix
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_h0_family
    !! NAME
    !!   likelihood_ncar_h0_family
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine likelihood_ncar_h0_family(ip,jm,n,x,vci,determ,f,iuser,user)
        integer         , intent(in)                  :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in) :: x
        real (kind=dp)  , intent(inout)               :: f
        integer ,       dimension(1), intent(inout)      :: iuser
        real (kind=dp)      ,dimension(1), intent(inout) :: user
        real (kind=dp)                   ,intent(in) :: determ
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar) ,intent(in)    :: VCI
        real (kind=dp) :: vpf
        integer        :: ifem
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V

        !
        integer :: ic,jc,s,nbnivest,kd1,kd2,kd,na
        logical :: ok,valide

        ! print *,x
        f = 0.d0
        V = 0.d0
        s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2

        ifem=jm-p_dg%nmp(ip)
        kd1=dataset%lSires(ip)%full_sib(ifem)%firstkd
        kd2=dataset%lSires(ip)%full_sib(ifem)%lastkd
        na = kd2-kd1+1
        if (kd2<kd1) then
            f=0
            return
        end if

        ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
        ! pour chaque caractere
        do ic=1,p_dpm%ncar
            nbnivest=my_listdesc(ic)%nbnivest
            call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
            V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
            s=s+my_listdesc(ic)%nbnivest
        end do


        vpf=0
        !calcul de la vraissemblance de plein frere
        do kd=1,na
            vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
        end do
        !finallement, la vraissemblance de la famille ip/jm:
        f=+0.5*vpf+dble(kd2-kd1+1)*0.5*log(determ)
    end subroutine likelihood_ncar_h0_family
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_hn_family
    !! NAME
    !!   likelihood_ncar_hn_family
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_hn_family(ip,jm,n,x,VCI,determ,f,iuser,user)
        integer         , intent(in)                                             :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in)                            :: x
        real (kind=dp)  , intent(inout)                                          :: f
        integer ,       dimension(1), intent(inout)                              :: iuser
        real (kind=dp) ,dimension(1), intent(inout)                              :: user
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)                           ,intent(in) :: VCI
        real(kind=dp)                                                ,intent(in) :: determ

        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V


        !
        integer        :: jjm,ig(my_listdesc(1)%nqtl),ifem,indf,indm,iq,ngg,iig,z,irank,ic,jc
        real(kind=dp)  :: effm,pbr,vmere,vpf
        integer        :: kd1,kd2,jj,nqtl,nbnivest,s,kd,na
        logical        :: ok,valide

        f=0.d0
        nqtl=my_listdesc(1)%nqtl

        jjm=jm-p_dg%nmp(ip)

        if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
            kd1=dataset%lSires(ip)%full_sib(jjm)%firstKd
            kd2=dataset%lSires(ip)%full_sib(jjm)%lastKd
            na = kd2-kd1+1
            effm=dble(na)
        else
            return
        end if

        vmere=0.d0
        if ( estime_multi(jm) ) ifem=iam_multi(p_dg%repfem(jm))
        !ngg : le nombre de genotype possible sur tous les qtls...
        ngg=1
        do iq=1,nqtl
            ig(iq)=p_spt%ngenom(current_chr(iq),jm)+1
            ngg=ngg*(p_spt%ngenom(current_chr(iq),jm+1)-p_spt%ngenom(current_chr(iq),jm))
        end do
        !pour toutes les combinaisons possibles des genotypes
        do iig=1,ngg
            pbr=1
            !on modifie la matrice d incidence pour les n qtl
            do iq=1,nqtl
                do ic=1,p_dpm%ncar
                    indm=my_listdesc(ic)%ntniv_qtlsires(iq)+ip-1
                    if ( estime_multi(jm) ) then
                        !Si la mere est estimable, on place dans la matrice les
                        !pdds dam et pdds sires (si ces effets sont estimables)
                        indf=my_listdesc(ic)%ntniv_qtldams(iq)+ifem-1
                        if ( my_listdesc(ic)%vecsol(indf)) then
                            indf=my_listdesc(ic)%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indf,ic)=dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                    else
                        !la mere n est pas estimable, on place seulement les pdds males
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%ppt(iq,kd1:kd2)
                        end if
                    end if
                end do ! ic
                !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
                pbr=pbr*p_spt%probg(current_chr(iq),ig(iq))
            end do ! iq
            ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
            ! pour chaque caractere
            s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2
            do ic=1,p_dpm%ncar
                nbnivest=my_listdesc(ic)%nbnivest
                call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                    My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                s=s+my_listdesc(ic)%nbnivest
            end do
            vpf=0
            !calcul de la vraissemblance de plein frere
            do kd=1,na
                vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
            end do
            vmere=vmere+pbr*dexp(-0.5d0*vpf)

            ! on increment
            ok=.true.
            do iq=1,nqtl
                if (ok) then
                    if ((ig(iq) < p_spt%ngenom(current_chr(iq),jm+1))) then
                        ig(iq)=ig(iq)+1
                        ok=.false.
                    end if
                end if
            end do
        end do ! iig


        if (vmere == 0) then
            f=INIFINY_REAL_VALUE
        else
            !finallement, la vraissemblance de la famille ip/jm:
            f=f-dlog(vmere)+dble(kd2-kd1+1)*0.5*log(determ)
        end if
        ! end do
          !  print *,ip,jm,f
          !  stop

    end subroutine likelihood_ncar_hn_family
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_h0_family_withcd
    !! NAME
    !!   likelihood_ncar_h0_family_withcd
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_h0_family_withcd(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                  :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in) :: x
        real (kind=dp)  , intent(inout)               :: f
        integer ,       dimension(1), intent(inout)      :: iuser
        real (kind=dp)      ,dimension(1), intent(inout) :: user


        real (kind=dp),dimension(dataset%nkd_max_by_fam) :: determ
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar,dataset%nkd_max_by_fam)    :: VCI
        real (kind=dp) :: vpf
        integer        :: ifem,kkk
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V

        !
        integer :: ic,jc,s,nbnivest,kd1,kd2,kd,na
        logical :: ok,valide

        determ=0.d0
        kd=1
        do kkk=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
            if ( count(p_dpa%presentc(:,kkk)) == p_dpm%ncar ) then
                call get_inv_residual_covariance_matrix_cd(ip,kkk,n,x,vci(:,:,kd),determ(kd),ok)
                if ( .not. ok ) then
                    f=INIFINY_REAL_VALUE
                    return
                end if
                kd = kd + 1
            end if
        end do

        f = 0.d0
        V = 0.d0
        s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2

        ifem=jm-p_dg%nmp(ip)

        kd1=dataset%lSires(ip)%full_sib(ifem)%firstkd
        kd2=dataset%lSires(ip)%full_sib(ifem)%lastkd
        na = kd2-kd1+1

        if (kd2<kd1) then
            f=0
            return
        end if

        ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
        ! pour chaque caractere
        do ic=1,p_dpm%ncar
            nbnivest=my_listdesc(ic)%nbnivest
            call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
            V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
            s=s+my_listdesc(ic)%nbnivest
        end do

        vpf=0
        !calcul de la vraissemblance de plein frere
        do kd=1,na
            vpf=vpf+dot_product(matmul(V(kd,:),VCI(:,:,kd)),V(kd,:))
        end do

        !finallement, la vraissemblance de la famille ip/jm:
        f=+0.5*vpf+0.5*sum(log(determ(:na)))
        !print *,f

    end subroutine likelihood_ncar_h0_family_withcd
    !!***

    !!****f* m_qtlmap_incidence_multi/get_inv_residual_covariance_matrix_cd
    !! NAME
    !!   get_inv_residual_covariance_matrix_cd
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine get_inv_residual_covariance_matrix_cd(ip,kd,n,x,vci,determ,ok)
        integer         , intent(in)                   :: ip,kd,n
        real (kind=dp)  ,dimension(n)   , intent(in)   :: x
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar),intent(out) :: VCI
        real(kind=dp),                     intent(out) :: determ
        logical                           ,intent(out) :: ok
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: rhoc
        real(kind=dp),dimension(p_dpm%ncar+1,p_dpm%ncar)    :: VCIInverse

        integer , parameter :: MAX_COUNT = 200

        real(kind=dp) :: rh
        integer :: ic,jc,k,irank,countt

        ok=.true.
        !RHO(I,J) avec i/=j => x((i-1)*(ncar-i+1))
        k=0
        do ic=2,p_dpm%ncar
            do jc=1,ic-1
                k=k+1
                ! *** RESCALE OFI ***
                !rhoc(jc,ic)=x(ncar*np+k)
                rh = x(p_dpm%ncar*p_dg%np+k)
                rhoc(jc,ic) = ((dexp(rh)-1.d0)/(dexp(rh)+1.d0))
                rhoc(ic,jc)=rhoc(jc,ic)
            end do
        end do

        !VCi : residual covariance matrix from sire I
        !         | SIGi1^2  SIGi1*SIGi2*RHO(1,2) ....  SIGi1*SIGip*RHO(1,p)
        !         | ..       SIGi2^2
        !    VCi  | ..                ...
        !         | ..                       ...
        !         | ..
        !         | SIGi1*SIGip*RHO(1,p)        SIGip^2
        !print *,'kd:',kd,cd(1,kd),my_corcd(1,2,kd)
        do ic=1,p_dpm%ncar
            do jc=ic,p_dpm%ncar
                VCI(ic,jc)=x((ic-1)*p_dg%np+ip)*x((jc-1)*p_dg%np+ip)
                if ( ic /= jc ) then
                    VCI(ic,jc)=VCI(ic,jc)*rhoc(ic,jc)*p_dataTrio%corcd(ic,jc,kd)
                    VCI(jc,ic)=VCI(ic,jc)
                else
                    VCI(ic,jc)=VCI(ic,jc)*p_dpa%cd(ic,kd)
                end if
            end do
        end do

        !Calcul Determinant de la matrice de covariance residuelle du pere ip
        !--------------------------------------------------------------------
        irank=0
        determ=0
        ! call MATH_QTLMAP_F03ABF(VCI,ncar,ncar,determ,irank)
        !      if (irank /= 0) then
        !        print *,"mauvais determinant matrice"
        !        ok=.false.
        !        return
        !      end if
        !Calcul Inverse de la matrice de covariance residuelle du pere ip
        !-----------------------------------------------------------------
        VCIInverse(:p_dpm%ncar,:p_dpm%ncar)=VCI
        irank = 1
        countt=0

        ! do while (irank /= 0 .and. countt <= MAX_COUNT)
        do while (irank /= 0)
            CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,VCIInverse,p_dpm%ncar+1,determ,irank)
            !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
            if (irank /= 0) then
                countt=countt+1
                !        print *,"mauvaise inversion matrice:",irank,' determ:',determ
                !        ok=.false.
                !        print *," *************** "
                !        print *," ** START VCI ** "
                !        do ic=1,ncar
                !         do jc=ic,ncar
                !          print *,x((ic-1)*np+ip)*x((jc-1)*np+ip)
                !         end do
                !        end do
                !        print *," ** RHOC ** "
                !        do ic=1,ncar
                !        print *,rhoc(ic,:)
                !        end do
                !        print *," ** VCI **"
                !        do ic=1,ncar
                !        print *,VCI(ic,:)
                !        end do
                !        print *," ** CORCD pour KD:",kd," **"
                !        do ic=1,ncar
                !        do jc=ic,ncar
                !            print *,ic,jc,datasetUser%corcd(ic,jc,kd)
                !        end do
                !        end do
                !        print *," ** CD pour KD:",kd," **"
                !        do ic=1,ncar
                !         print *,ic,cd(ic,kd)
                !        end do
                VCIInverse(:p_dpm%ncar,:p_dpm%ncar)=VCI
                do ic=1,p_dpm%ncar
                    VCIInverse(ic,ic)=VCIInverse(ic,ic)+0.1d0*real(countt)
                end do
            !        stop
            !  return
            end if
        end do

        !    if ( countt > MAX_COUNT ) then
        !        print *,"mauvaise inversion matrice countt:",countt
        !        return
        !    end if

        do ic=1,p_dpm%ncar
            do jc=1,ic
                VCI(ic,jc)=VCIInverse(ic+1,jc)
                VCI(jc,ic)=VCI(ic,jc)
            end do
        end do



    end subroutine get_inv_residual_covariance_matrix_cd
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_hn_family_withcd
    !! NAME
    !!   likelihood_ncar_hn_family_withcd
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_hn_family_withcd(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                                             :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in)                            :: x
        real (kind=dp)  , intent(inout)                                          :: f
        integer ,       dimension(1), intent(inout)                              :: iuser
        real (kind=dp) ,dimension(1), intent(inout)                              :: user


        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar,dataset%nkd_max_by_fam)  :: VCI
        real(kind=dp),dimension(dataset%nkd_max_by_fam)            :: determ
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar)       :: V
        !
        integer        :: jjm,ig(my_listdesc(1)%nqtl),ifem,indf,indm,iq,ngg,iig,z,irank,ic,jc
        real(kind=dp)  :: effm,pbr,vmere,vpf
        integer        :: kd1,kd2,jj,nqtl,nbnivest,s,kd,kkk,na
        logical        :: ok,valide

        determ=0.d0
        kd=1
        do kkk=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
            if ( count(p_dpa%presentc(:,kkk)) == p_dpm%ncar ) then
                call get_inv_residual_covariance_matrix_cd(ip,kkk,n,x,vci(:,:,kd),determ(kd),ok)
                if ( .not. ok ) then
                    f=INIFINY_REAL_VALUE
                    return
                end if
                kd = kd + 1
            end if
        end do

        f=0.d0
        nqtl=my_listdesc(1)%nqtl

        jjm=jm-p_dg%nmp(ip)

        if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
            kd1=dataset%lSires(ip)%full_sib(jjm)%firstKd
            kd2=dataset%lSires(ip)%full_sib(jjm)%lastKd
            na = kd2-kd1+1
            effm=dble(na)
        else
            return
        end if

        vmere=0.d0
        if ( estime_multi(jm) ) ifem=iam_multi(p_dg%repfem(jm))
        !ngg : le nombre de genotype possible sur tous les qtls...
        ngg=1
        do iq=1,nqtl
            ig(iq)=p_spt%ngenom(current_chr(iq),jm)+1
            ngg=ngg*(p_spt%ngenom(current_chr(iq),jm+1)-p_spt%ngenom(current_chr(iq),jm))
        end do
        !pour toutes les combinaisons possibles des genotypes
        do iig=1,ngg
            pbr=1
            !on modifie la matrice d incidence pour les n qtl
            do iq=1,nqtl
                do ic=1,p_dpm%ncar
                    indm=my_listdesc(ic)%ntniv_qtlsires(iq)+ip-1
                    if ( estime_multi(jm) ) then
                        !Si la mere est estimable, on place dans la matrice les
                        !pdds dam et pdds sires (si ces effets sont estimables)
                        indf=my_listdesc(ic)%ntniv_qtldams(iq)+ifem-1
                        if ( my_listdesc(ic)%vecsol(indf)) then
                            indf=my_listdesc(ic)%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indf,ic)=dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                    else
                        !la mere n est pas estimable, on place seulement les pdds males
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%ppt(iq,kd1:kd2)
                        end if
                    end if
                end do ! ic
                !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
                pbr=pbr*p_spt%probg(current_chr(iq),ig(iq))
            end do ! iq
            ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
            ! pour chaque caractere
            s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2
            do ic=1,p_dpm%ncar
                nbnivest=my_listdesc(ic)%nbnivest
                call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                    My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                s=s+my_listdesc(ic)%nbnivest
            end do
            vpf=0
            !calcul de la vraissemblance de plein frere
            do kd=1,na
                vpf=vpf+dot_product(matmul(V(kd,:),vci(:,:,kd)),V(kd,:))
            end do
            vmere=vmere+pbr*dexp(-0.5d0*vpf)
            ! on increment
            ok=.true.
            do iq=1,nqtl
                if (ok) then
                    if ((ig(iq) < p_spt%ngenom(current_chr(iq),jm+1))) then
                        ig(iq)=ig(iq)+1
                        ok=.false.
                    end if
                end if
            end do
        end do ! iig


        if (vmere == 0) then
            f=INIFINY_REAL_VALUE
        else
            !finallement, la vraissemblance de la famille ip/jm:
            f=f-dlog(vmere)+0.5*sum(log(determ(:na)))
        end if
    !        print *,f

    end subroutine likelihood_ncar_hn_family_withcd
    !!***



    !************************************************************* IMPLEMENTATION LU *****************************************************************************************************************************


    !!****f* m_qtlmap_incidence_multi/model_optim_multi_h0_LU
    !! NAME
    !!   model_optim_multi_h0_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE

    subroutine model_optim_multi_h0_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,performPrecision)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)                   ,intent(inout)    :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)      , intent(inout)  :: listDesc
        logical                                    ,intent(in)         :: performPrecision

        real(kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) :: xxx


        call model_optim_multi_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
            performPrecision,likelihood_ncar_h0_family_LU,likelihood_ncar_h0_LU_family_withcd,.true.,xxx)

    !         call model_optim_multi_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
    !          performPrecision,likelihood_ncar_h0_LU,likelihood_ncar_h0_LU,.true.,xxx)


    end subroutine model_optim_multi_h0_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/model_optim_multi_hn_LU
    !! NAME
    !!   model_optim_multi_hn_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine model_optim_multi_hn_LU(xinc,listDesc,curPos,workstruct,sigsquareEstime,rhoiestim,Bestim,&
        performPrecision,tConf,tempForConfusion,invlrt)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)                   ,intent(inout)    :: workstruct
        type(POSITION_LRT_INCIDENCE)                 ,intent(inout)    :: curPos
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)      , intent(inout)  :: listDesc
        logical                                    ,intent(in)         :: performPrecision,tConf
        real (kind=dp),dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) ,intent(out)  :: tempForConfusion
        logical                                          ,intent(in)   :: invlrt

        integer :: i,hypothesis
        character(len=LEN_L) :: dx

        allocate (current_chr(workstruct%nqtl))
        current_chr=curPos%listChr

        call model_optim_multi_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
            performPrecision,likelihood_ncar_hn_family_LU,likelihood_ncar_hn_LU_family_withcd,tConf,tempForConfusion)

        !         call model_optim_multi_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
        !           performPrecision,likelihood_ncar_hn_LU,likelihood_ncar_hn_LU,tConf,tempForConfusion)
        !         print *,"F:",workstruct%fmax


        if ((any(listDesc(1)%fperemax == INIFINY_REAL_VALUE)) ) then
            dx=""
            do i=1,workstruct%nqtl-1
                dx=trim(dx)//trim(str(dataset%data%map%absi(curPos%listChr(i),curPos%listN(i))))//","
            end do

            dx=trim(dx)//str(dataset%data%map%absi(curPos%listChr(workstruct%nqtl),curPos%listN(workstruct%nqtl)))
            call log_mess("dx ["//trim(dx)// "]. Can not optimize likelihood....The start point is reinitializing",WARNING_DEF)
            call init_startpoint(workStruct,listDesc(1))
            curPos%lrtSires=0.d0
            listDesc(1)%fperemax=0.d0
            listDesc(1)%fmeremax=0.d0
            return
        end if

        !compute LRT
        if (invLrt) then
            do hypothesis=1,workstruct%nqtl
                curPos%lrtSires(hypothesis,:)=-2.d0*(workstruct%fnqtlsires(hypothesis,:)-listDesc(1)%fperemax)
                curPos%lrtDams(hypothesis,:)=-2.d0*(workstruct%fnqtldams(hypothesis,:)-listDesc(1)%fmeremax)
            end do
        else
            do hypothesis=1,workstruct%nqtl
                curPos%lrtSires(hypothesis,:)=-2.d0*(listDesc(1)%fperemax-workstruct%fnqtlsires(hypothesis,:))
                curPos%lrtDams(hypothesis,:)=-2.d0*(listDesc(1)%fmeremax-workstruct%fnqtldams(hypothesis,:))
            end do
        end if

        deallocate (current_chr)

    end subroutine model_optim_multi_hn_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/model_optim_multi_LU
    !! NAME
    !!   model_optim_multi_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine model_optim_multi_LU(xinc,listDesc,workstruct,sigsquareEstime,rhoiestim,Bestim,&
        performPrecision,FUNCT_PART,FUNCT_PART_CENSURE,tConf,tempForConfusion)
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal), intent(in)       :: xinc
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)         :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)        ,intent(out)        :: sigsquareEstime
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)  ,intent(out)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)        :: rhoiestim
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)   , intent(inout),target  :: listDesc
        logical                                    ,intent(in)         :: performPrecision,tConf
        real (kind=dp),dimension(ntnivmaxtotal,ntnivmaxtotal,p_dpm%ncar) ,intent(out)  :: tempForConfusion
        external                                                       :: FUNCT_PART,FUNCT_PART_CENSURE

        integer :: j,i,ip,ifail,ibound,npar
        !         real(kind=dp) ,dimension(np+ntnivmax) :: par
        real (kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal)   :: XX
        real (kind=dp) , dimension(ntnivmaxtotal,ntnivmaxtotal)   :: triang
        integer                                    :: iuser(1),ix,kd1,kd2,jm,k,jjm,ic,s,jc,ncar,np
        real (kind=dp)   ,dimension(1)    :: user
        logical       ,dimension(p_dpm%ncar,ntnivmaxtotal)    :: lastvecsol
        real(kind=dp) :: vci_t(p_dpm%ncar,p_dpm%ncar,p_dg%np),determ_t(p_dg%np),f,r
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: mat_L,V
        real(kind=dp),dimension(p_dg%np,p_dpm%ncar)        :: mat_H
        real(kind=dp)                           :: sigtemp

        !Filter for the mimimization - dim (np,nm,nd)
        logical, dimension(p_dg%np,p_dg%nm,p_dpm%ncar*p_dg%np+p_dpm%ncar*(p_dpm%ncar-1)/2+p_dpm%ncar*ntnivmaxtotal):: filter_inc
        ncar = p_dpm%ncar
        np = p_dg%np
        iuser=0
        user=0
        my_listDesc=>listDesc
        dataset => listDesc(1)%dataset
        allocate (my_xincreduitmul(p_dg%nd,ntnivmaxtotal,ncar))

        my_xincreduitmul=0.d0
        do ic=1,p_dpm%ncar
            lastvecsol(ic,:my_listDesc(ic)%ntniv)=my_listDesc(ic)%vecsol(:my_listDesc(ic)%ntniv)
            ! create X'.X matrix from incidence matrix
            call model_XT_X(xinc(ic,:,:),my_listDesc(ic),XX)
            ! Check all parameters to remove from the estimation
            call estim_cholesky(XX,my_listDesc(ic),ntnivmaxtotal,triang)
            ! compute the precision of each parameter
            if (performPrecision) call get_precision(XX,tempForConfusion(:,:,ic),my_listDesc(ic))
            call set_corrxinc(xinc(ic,:,:),my_listDesc(ic),my_xincreduitmul(:,:,ic))
            !optimisation : on sauvegarde les index des elements non nul de la matrice reduite
            call fill_nonull_elements(my_listDesc(ic),size(my_xincreduitmul,1),size(my_xincreduitmul,2),my_xincreduitmul(:,:,ic))
          !call debug_write_incidence(xinc(ic,:,:),my_listDesc(ic))
        end do

          ! Optimisation de la vraisemblance a la position dx
        ifail=1
        ibound=0
        npar=ncar*np+ncar*(ncar-1)/2

        j=ncar*np+ncar*(ncar-1)/2
        do ic=1,ncar
            npar=npar+my_listDesc(ic)%nbnivest

            if (count(lastvecsol(ic,:))>0) then ! autre que initialisation
                do i=1,my_listDesc(ic)%ntniv
                    if(my_listDesc(ic)%vecsol(i))then
                        j=j+1
                        if(.not. lastvecsol(ic,i)) my_listDesc(ic)%par(j)=0.d0
                    end if
                end do
            end if
        end do

        filter_inc=.true.

        !         allocate (filter_inc(np,npar))
        s=np*ncar+ncar*(ncar-1)/2
        !         filter_inc=.true.
        !         filter_inc(:,:,1:np*ncar)=.false.
        !         filter_inc(:,:,np*ncar+1:s)=.true. !rho
        filter_inc(:,:,1:(np-1)*ncar)=.false.
        do ic=1,ncar
            do ip=1,np
                if ( ip>1)  then
                    filter_inc(ip,:,(ip-2)*ncar+1:(ip-1)*ncar)=.true. !rho
                endif
                filter_inc(ip,p_dg%nmp(ip)+1:p_dg%nmp(ip+1),(ic-1)*p_dg%np+ip)=.true.
                jjm=0
                do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
                    jjm=jjm+1
                    if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
                        kd1=dataset%lSires(ip)%full_sib(jjm)%firstkd
                        kd2=dataset%lSires(ip)%full_sib(jjm)%lastkd
                        filter_inc(ip,jm,s+1:s+my_listDesc(ic)%nbnivest)=&
                         any(my_xincreduitmul(kd1:kd2,:my_listDesc(ic)%nbnivest,ic)/=0.d0,dim=1)
                    else
                        !     print *,'pas de desc....'
                        filter_inc(ip,jm,s+1:s+my_listDesc(ic)%nbnivest)=.false.
                    end if
                end do
            end do
            s=s+my_listDesc(ic)%nbnivest
        end do

        !         print *,"PAR*** H"
        !         print *,my_listDesc(1)%par(1:(np-1)*ncar)
        !         print *,"PAR*** L"
        !         print *,my_listDesc(1)%par((np-1)*ncar+1:(np-1)*ncar+ncar*(ncar+1)/2)
        !         print *,"PAR*** B"
        !         print *,my_listDesc(1)%par((np-1)*ncar+ncar*(ncar+1)/2+1:(np-1)*ncar+ncar*(ncar+1)/2+my_listDesc(1)%ntniv)

        if ( dataset%data%datasetUser%na>0 ) then
            !prise en compte des donnees censures
            call minimizing_funct_family(dataset%data,npar,ibound,FUNCT_PART_CENSURE,filter_inc,&
                my_listDesc(1)%fmeremax,my_listDesc(1)%fperemax,&
                my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,f,iuser,user,ifail)
        !        call minimizing_funct(npar,ibound,FUNCT_PART_CENSURE,&
        !          my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,listDesc(1)%fmax,iuser,user,ifail)
        else
            !pas de donnee censure=> optimisation pour le calcul de l'inverse de matrice de covariance
            call minimizing_funct_family(dataset%data,npar,ibound,FUNCT_PART,filter_inc,&
                my_listDesc(1)%fmeremax,my_listDesc(1)%fperemax,&
                my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,f,iuser,user,ifail)
        !          call minimizing_funct(npar,ibound,FUNCT_PART,&
        !          my_listDesc(1)%borni,my_listDesc(1)%borns,my_listDesc(1)%par,listDesc(1)%fmax,iuser,user,ifail)
        end if

        Bestim=0.d0


        ! on fixe pour ip=1 H
        mat_H(1,:)=1.d0

        do ip=2,np
            do jc=1,ncar
                mat_H(ip,jc)=my_listDesc(1)%par((ip-2)*ncar+jc)
            end do
        end do

        mat_L=0.d0
        k=0
        do ic=1,ncar
            do jc=1,ic
                k=k+1
                mat_L(ic,jc) = my_listDesc(1)%par(ncar*(np-1)+k)
            end do
        end do

        V = matmul(mat_L,transpose(mat_L))

        s=np*ncar+ncar*(ncar-1)/2
        !getting standart deviation
        do ic=1,ncar
            do ip=1,np
                sigsquareEstime(ic,ip)=V(ic,ic)*mat_H(ip,ic)*mat_H(ip,ic)
            end do
             !The solution
            Bestim(ic,:my_listDesc(ic)%nbnivest)=my_listDesc(1)%par(s+1:s+my_listDesc(ic)%nbnivest)
            s=s+my_listDesc(ic)%nbnivest
        end do

        rhoiestim=0.d0

        do ic=1,ncar
            do jc=1,ic
                rhoiestim(ic,jc)=V(ic,jc) / ( dsqrt( V(ic,ic)*V(jc,jc)) )
                rhoiestim(jc,ic)=rhoiestim(ic,jc)
            end do
        end do

        deallocate (my_xincreduitmul)

    end subroutine model_optim_multi_LU
    !!***


    !!****f* m_qtlmap_incidence_multi/get_inv_residual_covariance_matrix_LU
    !! NAME
    !!   get_inv_residual_covariance_matrix_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine get_inv_residual_covariance_matrix_LU(ip,n,x,vci,determ,ok)
        integer         , intent(in)                   :: ip,n
        real (kind=dp)  ,dimension(n)   , intent(in)   :: x
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar),intent(out) :: VCI
        real(kind=dp)       ,         intent(out) :: determ
        logical                           ,intent(out) :: ok
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: mat_L,V
        real(kind=dp),dimension(p_dpm%ncar)           :: mat_H
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: VCITemp
        real(kind=dp),dimension(p_dpm%ncar+1,p_dpm%ncar)    :: VCIInverse

        integer :: ic,jc,k,irank

        ok=.true.
        !RHO(I,J) avec i/=j => x((i-1)*(ncar-i+1))
        k=0
        mat_L=0.d0
        do ic=1,p_dpm%ncar
            do jc=1,ic
                k=k+1
                mat_L(ic,jc) = x(p_dpm%ncar*(p_dg%np-1)+k)
            end do
        end do

        V = matmul(mat_L,transpose(mat_L))
        mat_H=0.d0
         ! on fixe pour ip=1 H
        if ( ip == 1) then
            mat_H(:)=1.d0
        else
            !  do ip=2,np
            do jc=1,p_dpm%ncar
                mat_H(jc)=x((ip-2)*p_dpm%ncar+jc)
            end do
        !   end do
        end if

        determ=0

        do ic=1,p_dpm%ncar
            do jc=1,p_dpm%ncar
                VCIInverse(ic,jc) = V(ic,jc)*mat_H(ic)*mat_H(jc)
            end do
        end do

        !Calcul Determinant de la matrice de covariance residuelle du pere ip
        !--------------------------------------------------------------------
        irank=0

        ! call MATH_QTLMAP_F03ABF(VCI,ncar,ncar,determ,irank)
        !      if (irank /= 0) then
        !        print *,"mauvais determinant matrice"
        !        ok=.false.
        !        return
        !      end if
        !Calcul Inverse de la matrice de covariance residuelle du pere ip
        !-----------------------------------------------------------------
        !VCIInverse(:ncar,:ncar)=VCI
        CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,VCIInverse,p_dpm%ncar+1,determ,irank)
        !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
        if (irank /= 0) then
            print *,"mat_L"
            do ic=1,p_dpm%ncar
                print *,VCIInverse(ic,:p_dpm%ncar)
            end do
            print *,"mauvaise inversion matrice"
            ok=.false.
            stop
            return
        end if

        do ic=1,p_dpm%ncar
            do jc=1,ic
                VCI(ic,jc)=VCIInverse(ic+1,jc)
                VCI(jc,ic)=VCI(ic,jc)
            end do
        end do


    end subroutine get_inv_residual_covariance_matrix_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/get_inv_residual_covariance_matrix_LU_cd
    !! NAME
    !!   get_inv_residual_covariance_matrix_LU_cd
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine get_inv_residual_covariance_matrix_LU_cd(ip,kd,n,x,vci,determ,ok)
        integer         , intent(in)                   :: ip,kd,n
        real (kind=dp)  ,dimension(n)   , intent(in)   :: x
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar),intent(out) :: VCI
        real(kind=dp),                     intent(out) :: determ
        logical                           ,intent(out) :: ok

        real(kind=dp),dimension(p_dpm%ncar+1,p_dpm%ncar)    :: VCIInverse

        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar)      :: mat_L,V
        real(kind=dp),dimension(p_dpm%ncar)           :: mat_H

        real(kind=dp) :: rh
        integer :: ic,jc,k,irank

        ok=.true.
        !RHO(I,J) avec i/=j => x((i-1)*(ncar-i+1))
        k=0
        mat_L=0.d0
        do ic=1,p_dpm%ncar
            do jc=1,ic
                k=k+1
                mat_L(ic,jc) = x(p_dpm%ncar*(p_dg%np-1)+k)
            end do
        end do

        V = matmul(mat_L,transpose(mat_L))

        mat_H=0.d0
         ! on fixe pour ip=1 H
        if ( ip == 1) then
            mat_H(:)=1.d0
        else
            !  do ip=2,np
            do jc=1,p_dpm%ncar
                mat_H(jc)=x((ip-2)*p_dpm%ncar+jc)
            end do
        !   end do
        end if


        determ=0
        VCIInverse=0.d0
        do ic=1,p_dpm%ncar
            do jc=1,p_dpm%ncar
                VCIInverse(ic,jc) = V(ic,jc)*mat_H(ic)*mat_H(jc)
                if ( ic /= jc ) then
                    VCIInverse(ic,jc) = VCIInverse(ic,jc)*p_dataTrio%corcd(ic,jc,kd)
                else
                    VCIInverse(ic,jc) = VCIInverse(ic,jc)*p_dpa%cd(ic,kd)
                end if
            !          VCIInverse(jc,ic) = VCIInverse(ic,jc)
            end do
        end do

        !Calcul Determinant de la matrice de covariance residuelle du pere ip
        !--------------------------------------------------------------------
        irank=0

        ! call MATH_QTLMAP_F03ABF(VCI,ncar,ncar,determ,irank)
        !      if (irank /= 0) then
        !        print *,"mauvais determinant matrice"
        !        ok=.false.
        !        return
        !      end if
        !Calcul Inverse de la matrice de covariance residuelle du pere ip
        !-----------------------------------------------------------------
        !  VCIInverse(:ncar,:ncar)=VCI

        !      CALL MATH_QTLMAP_INVDETMAT(ncar,VCIInverse,determ,irank)
        !
        !    !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
        !       if (irank /= 0) then
        !        do ic=1,ncar
        !          print *,VCIInverse(ic,:ncar)
        !       end do
        !        print *,"mauvaise inversion matrice"
        !        ok=.false.
        !        stop
        !        return
        !      end if
        !
        !      do ic=1,ncar
        !        do jc=ic,ncar
        !            VCI(ic,jc)=VCIInverse(ic,jc)
        !            VCI(jc,ic)=VCI(ic,jc)
        !        end do
        !      end do

        CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,VCIInverse,p_dpm%ncar+1,determ,irank)
        !  call MATH_QTLMAP_F01ADF(ncar,VCIInverse,ncar+1,irank)
        if (irank /= 0) then
            print *,"mauvaise inversion matrice"
            ok=.false.
            return
        end if

        do ic=1,p_dpm%ncar
            do jc=1,ic
                VCI(ic,jc)=VCIInverse(ic+1,jc)
                VCI(jc,ic)=VCI(ic,jc)
            end do
        end do

    end subroutine get_inv_residual_covariance_matrix_LU_cd
    !!***


    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_h0_family_LU
    !! NAME
    !!   likelihood_ncar_h0_family_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine likelihood_ncar_h0_family_LU(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                     :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in)    :: x
        real (kind=dp)  , intent(inout)                  :: f
        integer ,       dimension(1), intent(inout)      :: iuser
        real (kind=dp)      ,dimension(1), intent(inout) :: user
        real (kind=dp)                                   :: determ
        real(kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar)              :: VCI
        ! real(kind=dp),dimension(ncar,ncar)               :: VCITemp
        real (kind=dp) :: vpf
        integer        :: ifem
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V

        !
        integer :: ic,jc,s,nbnivest,kd1,kd2,kd,na
        logical :: ok,valide

        call get_inv_residual_covariance_matrix_LU(ip,n,x,vci,determ,ok)
        f = 0.d0

        !  do ip=1,np
        !   VCITemp = VCI(:,:)
        !  do jm=nmp(ip)+1,nmp(ip+1)
        V = 0.d0
        s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2

        ifem=jm-p_dg%nmp(ip)
        kd1=dataset%lSires(ip)%full_sib(ifem)%firstkd
        kd2=dataset%lSires(ip)%full_sib(ifem)%lastkd
        na = kd2-kd1+1
        if (kd2<kd1) then
            f=0
            return
        end if

        ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
        ! pour chaque caractere
        do ic=1,p_dpm%ncar
            nbnivest=my_listdesc(ic)%nbnivest
            call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
            V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
            s=s+my_listdesc(ic)%nbnivest
        end do


        vpf=0
        !calcul de la vraissemblance de plein frere
        do kd=1,na
            vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
        end do
        !finallement, la vraissemblance de la famille ip/jm:
        f=0.5*vpf+dble(kd2-kd1+1)*0.5*log(determ)
      !  print *,f
    !      end do
    !     end do
        !stop

    end subroutine likelihood_ncar_h0_family_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_h0_LU
    !! NAME
    !!   likelihood_ncar_h0_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine likelihood_ncar_h0_LU(n,x,f,iuser,user)
        integer         , intent(in)                     :: n
        real (kind=dp)      ,dimension(n), intent(in)    :: x
        real (kind=dp)  , intent(inout)                  :: f
        integer ,       dimension(1), intent(inout)      :: iuser
        real (kind=dp)      ,dimension(1), intent(inout) :: user
        real (kind=dp)                                   :: determ
        real(kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar)              :: VCI
        ! real(kind=dp),dimension(ncar,ncar)               :: VCITemp
        real (kind=dp) :: vpf
        integer        :: ifem,ip,jm
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V

        !
        integer :: ic,jc,s,nbnivest,kd1,kd2,kd,na
        logical :: ok,valide

        f = 0.d0

        do ip=1,p_dg%np
            call get_inv_residual_covariance_matrix_LU(ip,n,x,vci,determ,ok)
            !       VCITemp = VCI(:,:)
            do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
                V = 0.d0
                s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2

                ifem=jm-p_dg%nmp(ip)
                kd1=dataset%lSires(ip)%full_sib(ifem)%firstkd
                kd2=dataset%lSires(ip)%full_sib(ifem)%lastkd
                na = kd2-kd1+1
                if (kd2<kd1) then
                    cycle
                end if

                ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
                ! pour chaque caractere
                do ic=1,p_dpm%ncar
                    nbnivest=my_listdesc(ic)%nbnivest
                    call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                        My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                    V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                    s=s+my_listdesc(ic)%nbnivest
                end do


                vpf=0
                !calcul de la vraissemblance de plein frere
                do kd=1,na
                    vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
                end do
                !finallement, la vraissemblance de la famille ip/jm:
                f=f+0.5*vpf+dble(kd2-kd1+1)*0.5*log(determ)
            !  print *,f
            end do
        end do

      ! print *,f
        !stop

    end subroutine likelihood_ncar_h0_LU
    !!***


    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_h0_LU_family_withcd
    !! NAME
    !!   likelihood_ncar_h0_LU_family_withcd
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_h0_LU_family_withcd(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                  :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in) :: x
        real (kind=dp)  , intent(inout)               :: f
        integer ,       dimension(1), intent(inout)      :: iuser
        real (kind=dp)      ,dimension(1), intent(inout) :: user


        real (kind=dp),dimension(dataset%nkd_max_by_fam) :: determ
        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar,dataset%nkd_max_by_fam)    :: VCI
        real (kind=dp) :: vpf
        integer        :: ifem,kkk
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V

        !
        integer :: ic,jc,s,nbnivest,kd1,kd2,kd,na
        logical :: ok,valide

        determ=0.d0
        kd=1
        do kkk=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
            if ( count(p_dpa%presentc(:,kkk)) == p_dpm%ncar ) then
                call get_inv_residual_covariance_matrix_LU_cd(ip,kkk,n,x,vci(:,:,kd),determ(kd),ok)
                if ( .not. ok ) then
                    f=INIFINY_REAL_VALUE
                    return
                end if
                kd = kd + 1
            end if
        end do

        f = 0.d0
        V = 0.d0
        s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2

        ifem=jm-p_dg%nmp(ip)

        kd1=dataset%lSires(ip)%full_sib(ifem)%firstkd
        kd2=dataset%lSires(ip)%full_sib(ifem)%lastkd
        na = kd2-kd1+1

        if (kd2<kd1) then
            f=0
            return
        end if

        ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
        ! pour chaque caractere
        do ic=1,p_dpm%ncar
            nbnivest=my_listdesc(ic)%nbnivest
            call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
            V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
            s=s+my_listdesc(ic)%nbnivest
        end do

        vpf=0
        !calcul de la vraissemblance de plein frere
        do kd=1,na
            vpf=vpf+dot_product(matmul(V(kd,:),VCI(:,:,kd)),V(kd,:))
        end do

        !finallement, la vraissemblance de la famille ip/jm:
        f=+0.5*vpf+0.5*sum(log(determ(:na)))
        !print *,f

    end subroutine likelihood_ncar_h0_LU_family_withcd
    !!***


    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_hn_family_LU
    !! NAME
    !!   likelihood_ncar_hn_family_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_hn_family_LU(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                                             :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in)                            :: x
        real (kind=dp)  , intent(inout)                                          :: f
        integer ,       dimension(1), intent(inout)                              :: iuser
        real (kind=dp) ,dimension(1), intent(inout)                              :: user

        !local
        real(kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar)               :: VCI
        real(kind=dp)                                     :: determ
        !   real(kind=dp),dimension(ncar,ncar)                   :: VCITemp
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V


        !
        integer        :: jjm,ig(my_listdesc(1)%nqtl),ifem,indf,indm,iq,ngg,iig,z,irank,ic,jc
        real(kind=dp)  :: effm,pbr,vmere,vpf
        integer        :: kd1,kd2,jj,nqtl,nbnivest,s,kd,na
        logical        :: ok,valide

        f=0.d0
        nqtl=my_listdesc(1)%nqtl

        call get_inv_residual_covariance_matrix_LU(ip,n,x,vci,determ,ok)

        !      do ip=1,np
        !       VCITemp = VCI(ip,:,:)
        !      do jm=nmp(ip)+1,nmp(ip+1)

        jjm=jm-p_dg%nmp(ip)

        if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
            kd1=dataset%lSires(ip)%full_sib(jjm)%firstKd
            kd2=dataset%lSires(ip)%full_sib(jjm)%lastKd
            na = kd2-kd1+1
            effm=dble(na)
        else
            return
        end if

        vmere=0.d0
        if ( estime_multi(jm) ) ifem=iam_multi(p_dg%repfem(jm))
        !ngg : le nombre de genotype possible sur tous les qtls...
        ngg=1
        do iq=1,nqtl
            ig(iq)=p_spt%ngenom(current_chr(iq),jm)+1
            ngg=ngg*(p_spt%ngenom(current_chr(iq),jm+1)-p_spt%ngenom(current_chr(iq),jm))
        end do
        !pour toutes les combinaisons possibles des genotypes
        do iig=1,ngg
            pbr=1
            !on modifie la matrice d incidence pour les n qtl
            do iq=1,nqtl
                do ic=1,p_dpm%ncar
                    indm=my_listdesc(ic)%ntniv_qtlsires(iq)+ip-1
                    if ( estime_multi(jm) ) then
                        !Si la mere est estimable, on place dans la matrice les
                        !pdds dam et pdds sires (si ces effets sont estimables)
                        indf=my_listdesc(ic)%ntniv_qtldams(iq)+ifem-1
                        if ( my_listdesc(ic)%vecsol(indf)) then
                            indf=my_listdesc(ic)%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indf,ic)=dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                    else
                        !la mere n est pas estimable, on place seulement les pdds males
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%ppt(iq,kd1:kd2)
                        end if
                    end if
                end do ! ic
                !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
                pbr=pbr*p_spt%probg(current_chr(iq),ig(iq))
            end do ! iq
            ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
            ! pour chaque caractere
            s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2
            do ic=1,p_dpm%ncar
                nbnivest=my_listdesc(ic)%nbnivest
                call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                    My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                s=s+my_listdesc(ic)%nbnivest
            end do
            vpf=0
            !calcul de la vraissemblance de plein frere
            do kd=1,na
                vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
            end do
            vmere=vmere+pbr*dexp(-0.5d0*vpf)

            ! on increment
            ok=.true.
            do iq=1,nqtl
                if (ok) then
                    if ((ig(iq) < p_spt%ngenom(current_chr(iq),jm+1))) then
                        ig(iq)=ig(iq)+1
                        ok=.false.
                    end if
                end if
            end do
        end do ! iig


        if (vmere == 0) then
            f=INIFINY_REAL_VALUE
        else
            !finallement, la vraissemblance de la famille ip/jm:
            f=f-dlog(vmere)+dble(kd2-kd1+1)*0.5*log(determ)
        end if
    !     end do
    !     end do

    end subroutine likelihood_ncar_hn_family_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_hn_LU
    !! NAME
    !!   likelihood_ncar_hn_LU
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_hn_LU(n,x,f,iuser,user)
        integer         , intent(in)                                             :: n
        real (kind=dp)      ,dimension(n), intent(in)                            :: x
        real (kind=dp)  , intent(inout)                                          :: f
        integer ,       dimension(1), intent(inout)                              :: iuser
        real (kind=dp) ,dimension(1), intent(inout)                              :: user

        !local
        real(kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar)               :: VCI
        real(kind=dp)                                     :: determ
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar) :: V


        !
        integer        :: jjm,ig(my_listdesc(1)%nqtl),ifem,indf,indm,iq,ngg,iig,z,irank,ic,jc
        real(kind=dp)  :: effm,pbr,vmere,vpf
        integer        :: kd1,kd2,jj,nqtl,nbnivest,s,kd,na,ip,jm
        logical        :: ok,valide

        f=0.d0
        nqtl=my_listdesc(1)%nqtl

        do ip=1,p_dg%np
            call get_inv_residual_covariance_matrix_LU(ip,n,x,vci,determ,ok)
            do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)

                jjm=jm-p_dg%nmp(ip)

                if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
                    kd1=dataset%lSires(ip)%full_sib(jjm)%firstKd
                    kd2=dataset%lSires(ip)%full_sib(jjm)%lastKd
                    na = kd2-kd1+1
                    effm=dble(na)
                else
                    cycle
                end if

                vmere=0.d0
                if ( estime_multi(jm) ) ifem=iam_multi(p_dg%repfem(jm))
                !ngg : le nombre de genotype possible sur tous les qtls...
                ngg=1
                do iq=1,nqtl
                    ig(iq)=p_spt%ngenom(current_chr(iq),jm)+1
                    ngg=ngg*(p_spt%ngenom(current_chr(iq),jm+1)-p_spt%ngenom(current_chr(iq),jm))
                end do
                !pour toutes les combinaisons possibles des genotypes
                do iig=1,ngg
                    pbr=1
                    !on modifie la matrice d incidence pour les n qtl
                    do iq=1,nqtl
                        do ic=1,p_dpm%ncar
                            indm=my_listdesc(ic)%ntniv_qtlsires(iq)+ip-1
                            if ( estime_multi(jm) ) then
                                !Si la mere est estimable, on place dans la matrice les
                                !pdds dam et pdds sires (si ces effets sont estimables)
                                indf=my_listdesc(ic)%ntniv_qtldams(iq)+ifem-1
                                if ( my_listdesc(ic)%vecsol(indf)) then
                                    indf=my_listdesc(ic)%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                                    My_XINCREDUITMUL(kd1:kd2,indf,ic)=dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                                        p_spt%ngenom(current_chr(iq),jm)),:)
                                end if
                                if ( my_listdesc(ic)%vecsol(indm) ) then
                                    indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                                    My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                                        p_spt%ngenom(current_chr(iq),jm)),:)
                                end if
                            else
                                !la mere n est pas estimable, on place seulement les pdds males
                                if ( my_listdesc(ic)%vecsol(indm) ) then
                                    indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                                    My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%ppt(iq,kd1:kd2)
                                end if
                            end if
                        end do ! ic
                        !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
                        pbr=pbr*p_spt%probg(current_chr(iq),ig(iq))
                    end do ! iq
                    ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
                    ! pour chaque caractere
                    s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2
                    do ic=1,p_dpm%ncar
                        nbnivest=my_listdesc(ic)%nbnivest
                        call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                            My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                        V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                        s=s+my_listdesc(ic)%nbnivest
                    end do
                    vpf=0
                    !calcul de la vraissemblance de plein frere
                    do kd=1,na
                        vpf=vpf+dot_product(matmul(V(kd,:),VCI),V(kd,:))
                    end do
                    vmere=vmere+pbr*dexp(-0.5d0*vpf)

                    ! on increment
                    ok=.true.
                    do iq=1,nqtl
                        if (ok) then
                            if ((ig(iq) < p_spt%ngenom(current_chr(iq),jm+1))) then
                                ig(iq)=ig(iq)+1
                                ok=.false.
                            end if
                        end if
                    end do
                end do ! iig


                if (vmere == 0) then
                    f=INIFINY_REAL_VALUE
                else
                    !finallement, la vraissemblance de la famille ip/jm:
                    f=f-dlog(vmere)+dble(kd2-kd1+1)*0.5*log(determ)
                end if
            end do
        end do

    end subroutine likelihood_ncar_hn_LU
    !!***

    !!****f* m_qtlmap_incidence_multi/likelihood_ncar_hn_LU_family_withcd
    !! NAME
    !!   likelihood_ncar_hn_LU_family_withcd
    !! DESCRIPTION
    !!
    !! NOTE
    !!    multivariate analysis
    !!     (indicecar-1)*np +ip  : SIG ip,indice_caracetere
    !!     ncar*np + (ic-1)*ncar+jc        : RHO (ic,jc) residual correlation between two traits....
    !!     ncar*np + ncar*ncar + nbnivest
    !! SOURCE
    subroutine likelihood_ncar_hn_LU_family_withcd(ip,jm,n,x,f,iuser,user)
        integer         , intent(in)                                             :: ip,jm,n
        real (kind=dp)      ,dimension(n), intent(in)                            :: x
        real (kind=dp)  , intent(inout)                                          :: f
        integer ,       dimension(1), intent(inout)                              :: iuser
        real (kind=dp) ,dimension(1), intent(inout)                              :: user


        real(kind=dp),dimension(p_dpm%ncar,p_dpm%ncar,dataset%nkd_max_by_fam)  :: VCI
        real(kind=dp),dimension(dataset%nkd_max_by_fam)            :: determ
        real(kind=dp),dimension(dataset%nkd_max_by_fam,p_dpm%ncar)       :: V
        !
        integer        :: jjm,ig(my_listdesc(1)%nqtl),ifem,indf,indm,iq,ngg,iig,z,irank,ic,jc
        real(kind=dp)  :: effm,pbr,vmere,vpf
        integer        :: kd1,kd2,jj,nqtl,nbnivest,s,kd,kkk,na
        logical        :: ok,valide

        determ=0.d0
        kd=1
        do kkk=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
            if ( count(p_dpa%presentc(:,kkk)) == p_dpm%ncar ) then
                call get_inv_residual_covariance_matrix_LU_cd(ip,kkk,n,x,vci(:,:,kd),determ(kd),ok)
                !  print *,kd,determ(kd)
                if ( .not. ok ) then
                    f=INIFINY_REAL_VALUE
                    return
                end if
                kd = kd + 1
            end if
        end do


        f=0.d0
        nqtl=my_listdesc(1)%nqtl

        jjm=jm-p_dg%nmp(ip)

        if (dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
            kd1=dataset%lSires(ip)%full_sib(jjm)%firstKd
            kd2=dataset%lSires(ip)%full_sib(jjm)%lastKd
            na = kd2-kd1+1
            effm=dble(na)
        else
            return
        end if

        vmere=0.d0
        if ( estime_multi(jm) ) ifem=iam_multi(p_dg%repfem(jm))
        !ngg : le nombre de genotype possible sur tous les qtls...
        ngg=1
        do iq=1,nqtl
            ig(iq)=p_spt%ngenom(current_chr(iq),jm)+1
            ngg=ngg*(p_spt%ngenom(current_chr(iq),jm+1)-p_spt%ngenom(current_chr(iq),jm))
        end do
        !pour toutes les combinaisons possibles des genotypes
        do iig=1,ngg
            pbr=1
            !on modifie la matrice d incidence pour les n qtl
            do iq=1,nqtl
                do ic=1,p_dpm%ncar
                    indm=my_listdesc(ic)%ntniv_qtlsires(iq)+ip-1
                    if ( estime_multi(jm) ) then
                        !Si la mere est estimable, on place dans la matrice les
                        !pdds dam et pdds sires (si ces effets sont estimables)
                        indf=my_listdesc(ic)%ntniv_qtldams(iq)+ifem-1
                        if ( my_listdesc(ic)%vecsol(indf)) then
                            indf=my_listdesc(ic)%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indf,ic)=dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                                p_spt%ngenom(current_chr(iq),jm)),:)
                        end if
                    else
                        !la mere n est pas estimable, on place seulement les pdds males
                        if ( my_listdesc(ic)%vecsol(indm) ) then
                            indm=my_listdesc(ic)%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                            My_XINCREDUITMUL(kd1:kd2,indm,ic)=dataset%lSires(ip)%ppt(iq,kd1:kd2)
                        end if
                    end if
                end do ! ic
                !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
                pbr=pbr*p_spt%probg(current_chr(iq),ig(iq))
            end do ! iq
            ! Dans la famille du pere ip, on calcul le Y (performance auquel on soustrait l'ensemble des effets)
            ! pour chaque caractere
            s=p_dg%np*p_dpm%ncar+p_dpm%ncar*(p_dpm%ncar-1)/2
            do ic=1,p_dpm%ncar
                nbnivest=my_listdesc(ic)%nbnivest
                call matmul_incidence(kd1,kd2,my_listdesc(ic),p_dg%nd,ntnivmaxtotal,&
                    My_XINCREDUITMUL(:,:,ic),X(s+1:s+nbnivest),dataset%nkd_max_by_fam,1,V(:,ic))
                V(:na,ic) = my_listdesc(ic)%dataset%Y(ic,kd1:kd2)-V(:na,ic)
                s=s+my_listdesc(ic)%nbnivest
            end do
            vpf=0
            !calcul de la vraissemblance de plein frere
            do kd=1,na
                vpf=vpf+dot_product(matmul(V(kd,:),vci(:,:,kd)),V(kd,:))
            end do
            vmere=vmere+pbr*dexp(-0.5d0*vpf)
            ! on increment
            ok=.true.
            do iq=1,nqtl
                if (ok) then
                    if ((ig(iq) < p_spt%ngenom(current_chr(iq),jm+1))) then
                        ig(iq)=ig(iq)+1
                        ok=.false.
                    end if
                end if
            end do
        end do ! iig


        if (vmere == 0) then
            f=INIFINY_REAL_VALUE
        else
            !finallement, la vraissemblance de la famille ip/jm:
            f=f-dlog(vmere)+0.5*sum(log(determ(:na)))
        end if
        !    print *,f

    end subroutine likelihood_ncar_hn_LU_family_withcd
    !!***

    !********************************************************************************* FIN LU *********************************************************************




    !!****f* m_qtlmap_incidence_multi/gen_loop_opti_nqtl_multi
    !! NAME
    !!   gen_loop_opti_nqtl_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    recursive subroutine gen_loop_opti_nqtl_multi(data, &
        spt,       &
        iqtl,      &
        nqtl,      &
        curPos,    &
        workstruct,&
        xinc,      &
        listDesc  , &
        sigsquareEstimeMax,&
        rhoiestimMax,      &
        BestimMax,     &
        lrtsol,        &
        printpercent,  &
        FUNCT_MODEL)
        type(QTLMAP_DATASET)                        ,intent(in) :: data
        type(PDD_BUILD)                             ,intent(in) :: spt
        integer                                 , intent(in)    :: nqtl  ! under hypothesis nqtl
        integer                                 , intent(in)    :: iqtl  ! iteration of inside loop
        type(POSITION_LRT_INCIDENCE)             ,intent(inout) :: curPos
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)  :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal) , intent(inout) :: xinc  ! incidence matrix
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)   , intent(inout) :: listDesc ! description of the incidence matrix
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal) , intent(inout)  :: BestimMax      ! estimation of parameter on the maximum finded
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)    :: rhoiestimMax
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np) , intent(inout)     :: sigsquareEstimeMax      ! estimation of parameter on the maximum finded
        type(TYPE_LRT_SOLUTION)               ,intent(inout)    :: lrtsol
        logical                              , intent(in)       :: printpercent
        external                                                :: FUNCT_MODEL


        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)              :: sigsquareEstime
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)            :: rhoiestim
        integer :: n,hypothesis,chr,chr_start,iip,ifem,ii,lnpo,chr1,n1,ip,i
        logical hypotheseOne, hypotheseTwo
        real (kind=dp) ,dimension(:,:,:),pointer      :: xxx

        allocate (xxx(p_dpm%ncar,ntnivmaxtotal,ntnivmaxtotal))
        ! print *,"       **************    loop_opti_nqtl ***********************"
        !   incidenceDescSave = incidenceDesc
        if ( iqtl == nqtl ) then
            hypotheseOne = ( lrtsol%hypothesis == 1 )
            hypotheseTwo = ( lrtsol%hypothesis == 2 )
            if (nqtl >= 2 ) then
                chr_start = curPos%listChr(nqtl-1)
            else
                chr_start = 1
            end if

            do chr=chr_start,data%map%nchr
                ! Step above the chromosome
                curPos%listChr(nqtl)=chr
                if ( (nqtl >= 2) .and. (chr == chr_start ) ) then
                    n = curPos%listN(nqtl-1)
                else
                    n=0
                end if

                do while (n < data%map%get_npo(chr))
                    n=n+1
                    curPos%listN(nqtl)=n
                    ! if (nqtl == 3)print *,"chr:",chr," n:",n
                     !
                    !add qtl effect/interaction at position n to estim
                    call change_qtleffect_multi(data,spt,workstruct%listnteff(iqtl),iqtl,chr,n,xinc,listDesc,0)
                    !call debug_write_incidence(xinc,incidenceDesc)
                    !                do ii=1,ncar
                    !             call debug_write_incidence(xinc(ii,:,:),listDesc(ii))
                    !           end do
                    !call model
                    call FUNCT_MODEL(xinc,listDesc,curPos,workstruct,sigsquareEstime,rhoiestim,Bestim,.false.,.false.,xxx,.false.)

                     do i=1,nqtl
                       call lrtsol%LRT%add(data,nqtl,curPos%listChr,curPos%listN,sum(curPos%lrtSires(i,:)),i)
                       do ip=1,data%genea%np
                        call lrtsol%LRT_SIRES(ip)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtSires(i,ip),i)
                       end do
                       do ii=1,data%genea%nm
                        call lrtsol%LRT_DAMS(ii)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtDams(i,ii),i)
                       end do
                     end do

                    if ( hypotheseOne ) then
!                        lrtsol%lrt1(chr,n)=sum(curPos%lrtSires(1,:))
!                        lrtsol%xlrp(curPos%listChr(1),:,curPos%listN(1))=curPos%lrtSires(1,:)
!                        lrtsol%xlrm(curPos%listChr(1),:,curPos%listN(1))=curPos%lrtDams(1,:)
                        iip=1
                        do ip=1,p_dg%np
                            if ( listDesc(1)%vecsol(1+ip)) then
                                iip=iip+1
                                ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                                lrtsol%pater_eff(chr,ip,n)=Bestim(1,iip)
                            end if
                        end do
                        ifem=0
                        do ii=1,p_dg%nm
                            if ( estime_multi(ii) ) then
                                ifem=ifem+1
                                if ( listDesc(1)%vecsol(p_dg%np+1+ifem) ) then
                                    iip=iip+1
                                    ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                                    lrtsol%mater_eff(chr,ifem,n)=Bestim(1,iip)
                                end if
                            end if
                        end do
                    end if

                    !keep maximum
                    if(lrtsol%lrtmax(nqtl-1) < sum(curPos%lrtSires(nqtl,:))) then
                        lrtsol%nxmax=curPos%listN
                        lrtsol%chrmax=curPos%listChr
                        BestimMax = Bestim
                        sigsquareEstimeMax=sigsquareEstime
                        rhoiestimMax=rhoiestim
                        do i=0,nqtl-1
                          lrtsol%lrtmax(i)=sum(curPos%lrtSires(i+1,:))
                        end do
                      !print *,"max sigsquareestime:",workstruct%sigsquareEstime
                    end if
                end do
            end do
        else

            if ( iqtl-2 >= 0  ) then
                chr_start=curPos%listChr(iqtl-1)
            else
                chr_start=1
            end if

            do chr=chr_start,data%map%nchr
                ! Step above the chromosome

                if ( iqtl-2 >= 0  .and. chr == chr_start) then
                    n=curPos%listN(iqtl-1)
                else
                    n=0
                end if

                curPos%listChr(iqtl)=chr
                lnpo=data%map%get_npo(chr)
                do while (n < lnpo) !ix=startx,ilong,pas
                    if (printpercent) call log_mess( trim(str((float(n)/float(lnpo))*100.d0))//"%", VERBOSE_DEF )
                    n=n+1
                    curPos%listN(iqtl)=n
                    !add qtl effect/interaction at position n to estim
                    call change_qtleffect_multi(data,spt,workstruct%listnteff(iqtl),iqtl,chr,n,xinc,listDesc,0)

                    ! ****recursivite ****
                    call gen_loop_opti_nqtl_multi(data,spt,iqtl+1,nqtl,curPos,workstruct,&
                        xinc,listDesc,sigsquareEstimeMax,rhoiestimMax,BestimMax,lrtsol,.false.,FUNCT_MODEL)

                end do
            end do
        end if
        deallocate (xxx)


    end subroutine gen_loop_opti_nqtl_multi
    !!***

    !!****f* m_qtlmap_incidence_multi/gen_opti_nqtl_multi
    !! NAME
    !!   gen_opti_nqtl_multi
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine gen_opti_nqtl_multi( data,     &
        spt,       &
        nqtl,      &
        hyp,         &
        curPosMax,    &
        workstruct,&
        xinc,      &
        listDescMax  , &
        sigsquareEstimeMax,&
        rhoiestimMax,      &
        BestimMax,     &
        lrtsol,        &
        FUNCT_MODEL)
        type(QTLMAP_DATASET)                    , intent(in)    :: data
        type(PDD_BUILD)                         , intent(in)    :: spt
        integer                                 , intent(in)    :: nqtl,hyp  ! under hypothesis nqtl
        type(POSITION_LRT_INCIDENCE)             ,intent(inout) :: curPosMax
        type(INCIDENCE_GEN_STRUCT)              ,intent(inout)  :: workstruct
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%nd,ntnivmaxtotal) , intent(inout) :: xinc  ! incidence matrix
        type(INCIDENCE_TYPE)   , dimension(p_dpm%ncar)   , intent(inout) :: listDescMax ! description of the incidence matrix
        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal) , intent(inout)  :: BestimMax      ! estimation of parameter on the maximum finded
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)      ,intent(out)    :: rhoiestimMax
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np) , intent(inout)     :: sigsquareEstimeMax      ! estimation of parameter on the maximum finded
        type(TYPE_LRT_SOLUTION)               ,intent(inout)    :: lrtsol
        external                                                :: FUNCT_MODEL


        real (kind=dp) , dimension(p_dpm%ncar,ntnivmaxtotal)        :: Bestim
        real (kind=dp) , dimension(p_dpm%ncar,p_dg%np)              :: sigsquareEstime
        real (kind=dp) , dimension(p_dpm%ncar,p_dpm%ncar)            :: rhoiestim

        real (kind=dp) ,dimension(:,:,:),pointer         :: xxx
        real (kind=dp) , dimension(:,:,:),allocatable :: xinc2

        logical :: hypotheseOne,hypotheseTwo,ok
        integer :: i,chr,nGL,nGLTotal,nGLFact,InGL,nlinMax,nlin,nmax(nqtl),chmax(nqtl)
        integer :: nlinValide,iqtl,ic,ii,ifem,iip,ip
        integer , dimension(:,:,:) ,allocatable :: bestim_save,sigsq_save,rhoi_save
        integer , dimension(:,:)   ,allocatable :: lchr,lnpos
        real(kind=dp) , dimension(:,:) ,allocatable :: lrt_save
        type(POSITION_LRT_INCIDENCE)             :: curPos
        type(INCIDENCE_TYPE) , dimension(p_dpm%ncar)   :: listDesc

        hypotheseOne = ( lrtsol%hypothesis == 1 )
        hypotheseTwo = ( lrtsol%hypothesis == 2 )

        nGL = 0

        !Nombre de point sur le groupe de liaison
        !on cherche le nombre de point dependant de l hyp nqtl du numbre de point sur le groupe de liason

        do chr=1,data%map%nchr
            nGL=nGL + data%map%get_npo(chr)
        end do

        nGLTotal=1
        do iqtl=1,nqtl
            nGLTotal=nGLTotal*nGL
        end do

        curPosMax%listN(nqtl)=0
        curPosMax%listN(1:nqtl-1)=1
        curPosMax%listChr=1

        allocate (lchr(nGLTotal,nqtl))
        allocate (lnpos(nGLTotal,nqtl))

        nlinValide=0
        do nlin=1,nGLTotal
            i=nqtl
            ok=.false.
            do while ( (.not. ok) .and. (i>=1))
                curPosMax%listN(i) = curPosMax%listN(i)+1
                if ( curPosMax%listN(i) > data%map%get_npo(curPosMax%listChr(i)) ) then
                    if ( data%map%nchr > curPosMax%listChr(i) ) then
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

            if ( ok ) then
                nlinValide=nlinValide+1
                lchr(nlinValide,:)=curPosMax%listChr
                lnpos(nlinValide,:)=curPosMax%listN
            end if
        end do
        nGLTotal=nlinValide

        allocate (bestim_save(nGLTotal,p_dpm%ncar,ntnivmaxtotal))
        allocate (sigsq_save(nGLTotal,p_dpm%ncar,p_dg%np))
        allocate (rhoi_save(nGLTotal,p_dpm%ncar,p_dpm%ncar))
        allocate (lrt_save(nGLTotal,nqtl))
        lrt_save=0.d0
        !$OMP PARALLEL DEFAULT(SHARED)  &
        !$OMP PRIVATE(curPos,i,ic,iqtl,ok,rhoiestim,sigsquareEstime,Bestim,xxx,ip,iip,ifem,listDesc,xinc2)
        call init_position (data,hyp,nqtl,curPos)
        allocate (xxx(p_dpm%ncar,ntnivmaxtotal,ntnivmaxtotal))
        allocate (xinc2(p_dpm%ncar,p_dg%nd,ntnivmaxtotal))
        xinc2=xinc
        !initialisation of incidence matrix
        do ic=1,p_dpm%ncar
            call copy_incidence_desc (listDescMax(ic),listDesc(ic))
        end do

        curPos%lrtSires(nqtl,:)=0
        curPos%listN(nqtl)=0
        curPos%listN(1:nqtl-1)=1
        curPos%listChr=1

        !$OMP DO
        do nlin=1,nGLTotal
            call log_mess( trim(str((float(nlin)/float(nGLTotal))*100.d0))//"%", INFO_DEF )

            curPos%listN=lnpos(nlin,:)
            curPos%listChr=lchr(nlin,:)

            do iqtl=1,nqtl
                call change_qtleffect_multi(data,spt,workstruct%listnteff(iqtl),&
                    iqtl,curPos%listChr(iqtl),curPos%listN(iqtl),xinc2,listDesc,0)
            end do
            !           do ic=1,ncar
            !             call debug_write_incidence(xinc2(ic,:,:),listDesc(ic))
            !           end do

            call FUNCT_MODEL(xinc2,listDesc,curPos,workstruct,sigsquareEstime,rhoiestim,Bestim,.false.,.false.,xxx,.false.)

            !save value
            bestim_save(nlin,:,:)=Bestim
            sigsq_save(nlin,:,:)=sigsquareEstime
            rhoi_save(nlin,:,:)=rhoiestim
            do i=1,nqtl
              lrt_save(nlin,i)=sum(curPos%lrtSires(i,:))
              call lrtsol%LRT%add(data,nqtl,curPos%listChr,curPos%listN,lrt_save(nlin,i),i)

              do ip=1,data%genea%np
                 call lrtsol%LRT_SIRES(ip)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtSires(i,ip),i)
              end do

              do ii=1,data%genea%nm
                 call lrtsol%LRT_DAMS(ii)%add(data,nqtl,curPos%listChr,curPos%listN,curPos%lrtDams(i,ii),i)
              end do
            end do

            if ( hypotheseOne ) then
                iip=1
                do ip=1,p_dg%np
                    if ( listDesc(1)%vecsol(1+ip)) then
                        iip=iip+1
                        ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                        lrtsol%pater_eff(curPos%listChr(1),ip,curPos%listN(1))=Bestim(1,iip)
                    end if
                end do
                ifem=0
                do ii=1,p_dg%nm
                    if ( estime_multi(ii)) then
                        ifem=ifem+1
                        if ( listDesc(1)%vecsol(p_dg%np+1+ifem) ) then
                            iip=iip+1
                            ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
                            lrtsol%mater_eff(curPos%listChr(1),ifem,curPos%listN(1))=Bestim(1,iip)
                        end if
                    end if
                end do
            end if
        !    end if !ok
        end do ! nlin
        !$OMP END DO
        call end_position (curPos)
        deallocate (xxx)
        deallocate (xinc2)
        do ic=1,p_dpm%ncar
            call release_copy_incidence_desc(listDesc(ic))
        end do

        !$OMP END PARALLEL

        nlinMax=1
        lrtsol%lrtmax(0)=lrt_save(nlinMax,nqtl)
        nmax=1
        chmax=1

        call init_position (data,hyp,nqtl,curPos)

        curPos%listN(nqtl)=0
        curPos%listN(1:nqtl-1)=1
        curPos%listChr=1

        do nlin=1,nGLTotal
            curPos%listN=lnpos(nlin,:)
            curPos%listChr=lchr(nlin,:)

            if(lrtsol%lrtmax(nqtl-1) < lrt_save(nlin,nqtl)) then
                nmax=curPos%listN
                chmax=curPos%listChr
                nlinMax=nlin
                lrtsol%lrtmax=lrt_save(nlin,:)
            end if
        end do

        do iqtl=1,nqtl
            lrtsol%nxmax(iqtl-1)=nmax(iqtl)
            lrtsol%chrmax(iqtl-1)=chmax(iqtl)
            BestimMax = bestim_save(nlinMax,:,:)
            sigsquareEstimeMax=sigsq_save(nlinMax,:,:)
            rhoiestimMax=rhoi_save(nlinMax,:,:)
            curPosMax%listN=curPos%listN
            curPosMax%listChr=curPos%listChr
        end do

        deallocate (lchr)
        deallocate (lnpos)

        call end_position (curPos)

        deallocate (bestim_save)
        deallocate (sigsq_save)
        deallocate (rhoi_save)
        deallocate (lrt_save)

    end subroutine gen_opti_nqtl_multi
    !!***

    !!****f*ANALYSE/m_qtlmap_incidence_multi/set_solution_multi
    !! NAME
    !!    fill incsol the solution from the estimimation give in parameters
    !! DESCRIPTION
    !!
    !! HISTORY
    !!    09/09/2010 * modification de l affichage
    !! NOTES
    !!
    !! SOURCE
    subroutine set_solution_multi(data,hypothesis,workstruct,sigsquareEstime,rhoiestim,Bestim,listDesc,incsol,lennteff,listnteff)
        type(QTLMAP_DATASET)                     ,intent(in) :: data
        integer                        , intent(in)          :: hypothesis
        type(INCIDENCE_GEN_STRUCT)            ,intent(in)    :: workstruct
        real (kind=dp) , dimension(data%phenoModel%ncar,ntnivmaxtotal) , intent(in)    :: Bestim      ! estimation of parameter on the maximum finded
        real (kind=dp) , dimension(data%phenoModel%ncar,data%genea%np) , intent(in)     :: sigsquareEstime
        real (kind=dp) , dimension(data%phenoModel%ncar,data%phenoModel%ncar),intent(in)     :: rhoiestim
        type(INCIDENCE_TYPE)   ,dimension(data%phenoModel%ncar),intent(in)   :: listDesc
        type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol
        integer                               , intent(in)   :: lennteff
        integer        , dimension(lennteff)  , intent(in)   :: listnteff
        integer :: ic,jc,i,j,k,ip,maxNbPar,g,nb,INDEX,ntniv,nt(data%phenoModel%ncar),nteff,ii,ncar,np,nm
        logical :: addxmut
        type(GENEALOGY_BASE) , pointer :: dg
        type(DATAMODEL_BASE) , pointer :: dpm
        !      real(kind=dp) :: tempLin(ncar*ncar)

        dg => data%genea
        dpm => data%phenoModel

        ncar = dpm%ncar
        np = dg%np
        nm =dg%nm


        incsol%hypothesis=hypothesis
        call log_mess( "" , VERBOSE_DEF)
        call log_mess( "-------------------------------", VERBOSE_DEF)
        call log_mess( "Estimation of parameters under "//trim(str(incsol%hypothesis)), VERBOSE_DEF)
        call log_mess( "-------------------------------", VERBOSE_DEF)
        call log_mess( "rank max to resolve :"//trim(str(ntnivmaxtotal)), VERBOSE_DEF)

        allocate (incsol%sig(ncar,np))

        do ic=1,ncar-1
            if (listDesc(ic)%nteff /= listDesc(ic+1)%nteff) then
                print *,"Dev error: number of effect in the incidence matrix have to be equal..."
                stop
            end if
        end do

        nteff=listDesc(1)%nteff

        allocate (incsol%groupeName(nteff*ncar))
        allocate (incsol%eqtl_print(nteff*ncar))
        allocate (incsol%nbParameterGroup(nteff*ncar))

        if (hypothesis > 0) then
            allocate (incsol%qtl_groupeName(ncar,hypothesis))

            if ( lennteff < hypothesis ) then
                call stop_application("Devel error : bad init of array ** listnteff **")
            end if

            do i=1,hypothesis
                do ic=1,ncar
                    incsol%qtl_groupeName(ic,i)=ncar -1 + listnteff(i)+(ic-1)
                !    print *,"NUM EFFET:",incsol%qtl_groupeName(ic,i)
                end do
            end do
        end if

        maxNbPar=1
        do i=1,nteff
            do ic=1,ncar
                nb=0
                if ( listDesc(ic)%desc(i)%haveSubDesc) then
                    do j=1,size(listDesc(ic)%desc(i)%listSubDesc)
                        nb = nb + listDesc(ic)%desc(i)%listSubDesc(j)%end - listDesc(ic)%desc(i)%listSubDesc(j)%start + 1
                    end do
                end if
                maxNbPar = max(maxNbPar,nb)
            end do
        end do

        maxNbPar = max(maxNbPar,ncar) ! pour le general mean

        allocate (incsol%parameterName(nteff*ncar,maxNbPar))
        allocate (incsol%paramaterValue(nteff*ncar,maxNbPar))
        allocate (incsol%parameterVecsol(nteff*ncar,maxNbPar))
        allocate (incsol%parameterPrecis(nteff*ncar,maxNbPar))

        do ic=1,ncar
            do ip=1,np
                incsol%sig(ic,ip)=sqrt(sigsquareEstime(ic,ip))*dpm%sigt(ic)
            end do
        end do

        ii=0
        nt=0
        do i=1,nteff
            addxmut = .false.
            if ( workstruct%effects%general == i ) addxmut = .true.
            do ic=1,ncar
                ii = ii + 1
                incsol%eqtl_print(ii) = listDesc(ic)%eqtl_print(i)
                incsol%groupeName(ii) = trim(listDesc(ic)%desc(i)%name)//" ["//trim(dpm%carac(ic))//"]"
                if ( listDesc(ic)%desc(i)%haveSubDesc) then
                    incsol%nbParameterGroup(ii)=0
                    !do ic=1,ncar
                    do j=1,size(listDesc(ic)%desc(i)%listSubDesc)
                        INDEX=0
                        do k=listDesc(ic)%desc(i)%listSubDesc(j)%start,listDesc(ic)%desc(i)%listSubDesc(j)%end
                            INDEX=INDEX+1
                            incsol%nbParameterGroup(ii) = incsol%nbParameterGroup(ii) +1
                            g=incsol%nbParameterGroup(ii)
                            incsol%parameterName(ii,g)=listDesc(ic)%desc(i)%listSubDesc(j)%name//" "//adjustl(STR(INDEX))
                            ntniv=listDesc(ic)%desc(i)%listSubDesc(j)%start+INDEX -1
                            incsol%parameterVecsol(ii,g)=listDesc(ic)%vecsol(ntniv)
                            if ( incsol%parameterVecsol(ii,g) ) then
                                nt(ic) = nt(ic) + 1
                                incsol%paramaterValue(ii,g)= Bestim(ic,nt(ic))*dpm%sigt(ic)
                                if ( addxmut ) incsol%paramaterValue(ii,g) = incsol%paramaterValue(ii,g)+dpm%xmut(ic)
                                incsol%parameterPrecis(ii,g) = listDesc(ic)%precis(ntniv)
                             !  print *,"----"
                             !  print *,ic,nt(ic),Bestim(ic,nt(ic)),sigt(ic)
                             !  print *,incsol%groupeName(i),incsol%parameterName(i,g),incsol%paramaterValue(i,g),incsol%parameterPrecis(i,g)
                            end if
                        end do
                    end do
                else
                    !do ic=1,ncar
                    incsol%nbParameterGroup(ii)=1
                    incsol%parameterName(ii,ic)=trim(listDesc(ic)%desc(i)%name)//"*["//trim(dpm%carac(ic))//"]*"
                    incsol%parameterVecsol(ii,ic)=listDesc(ic)%vecsol(listDesc(ic)%desc(i)%start)
                    if ( incsol%parameterVecsol(ii,ic) ) then
                        nt(ic) = nt(ic) + 1
                        incsol%paramaterValue(ii,ic)  = Bestim(ic,nt(ic))*dpm%sigt(ic)
                        if ( addxmut ) incsol%paramaterValue(ii,g) = incsol%paramaterValue(ii,g)+dpm%xmut(ic)
                        incsol%parameterPrecis(ii,ic) = listDesc(ic)%precis(listDesc(ic)%desc(i)%start)
                    end if
                 !end do
                end if
            end do
        end do

        allocate (incsol%rhoi(ncar,ncar))
        incsol%rhoi=rhoiestim

    !          i=0
    !       do ic=2,ncar
    !         do jc=1,ic-1
    !           i=i+1
    !           tempLin(i)=rhoiestim(ic,jc)
    !         end do
    !       end do
    !
    !       i=0
    !       do jc=1,ncar-1
    !          do ic=jc+1,ncar
    !           i = i + 1
    !           incsol%rhoi(ic,jc) = tempLin(i)
    !           incsol%rhoi(jc,ic) = incsol%rhoi(ic,jc)
    !          end do
    !       end do

    end subroutine set_solution_multi
!!***


end module m_qtlmap_incidence_multi
