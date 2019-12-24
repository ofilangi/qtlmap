!!****m* ANALYSE/m_qtlmap_calcul_ic
!!  NAME
!!    m_qtlmap_calcul_ic -- Interface analysis to calcul confidence Intervals in QTL Mapping by Bootstrapping
!!  SYNOPSIS
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
module m_qtlmap_calcul_ic
    use m_qtlmap_log
    use m_qtlmap_types
    use m_qtlmap_output_handler
    use m_qtlmap_analyse
    use m_qtlmap_haplotype
    use m_qtlmap_genealogy
    use m_qtlmap_genotype
    use m_qtlmap_phenotype

    implicit none

    !USER INTERFACE
    !--------------
    public :: computingIC
   

contains

    !!****f* m_qtlmap_calcul_ic/boostrap
    !!  NAME
    !!    boostrap
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!
    !!  RETURN
    !!
    !!  SOURCE
    subroutine computingIC(dataset,lrtsol,listincsol,opt_qtl)
        type(QTLMAP_DATASET)                                      , intent(in)  :: dataset
        integer                                                   ,intent(in)   :: opt_qtl
        type(TYPE_LRT_SOLUTION)        ,dimension(dataset%phenoModel%ncar,opt_qtl)   , intent(inout)  :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION)  ,dimension(dataset%phenoModel%ncar,opt_qtl+1) , intent(in)  :: listincsol

        real (kind=dp) ,dimension(dataset%phenoModel%ncar,dataset%phenoModel%ncar)          :: rhoi
        character(len=LEN_DEF) ,dimension(10) :: values
        character(len=LEN_DEF) :: v,nameMethod
        integer :: nargs,opt_calcul,nsim,cli_haplo,i,id_ci,ic,iq
        logical :: hdam,biq,interaction

        call log_mess(" = Confidence Intervals =",INFO_DEF)

        opt_calcul = dataset%cli%cli_get_analyse()
        cli_haplo = dataset%cli%cli_get_haplotype()
        if ( dataset%cli%key_exist(dataset%cli%OPT_CI_NSIM) ) then
            call dataset%cli%get_key_value(dataset%cli%OPT_CI_NSIM,v)
            nsim = get_int(v)
        else
            nsim = 1000
        end if
        biq = dataset%cli%cli_is_biq()
        hdam = .not. dataset%cli%cli_is_no_hdam()
        interaction = dataset%cli%cli_is_interaction()

        !Liste des options pour le calcul Interval de Confiance
        if ( dataset%cli%key_exist(dataset%cli%OPT_CI) ) then
            call dataset%cli%get_key_list_values(dataset%cli%OPT_CI,10,values,nargs)
            if (nargs==0) then
                nargs=1
                values(nargs) = trim(str(DROP_OFF_CI))
                nargs=2
                values(nargs) = trim(str(HENGDE_LI_CI))
            end if
        else
            nargs=1
            values(nargs) = trim(str(DROP_OFF_CI))
            nargs=2
            values(nargs) = trim(str(HENGDE_LI_CI))
        end if

        !allocation
        do ic=1,dataset%phenoModel%ncar
            do iq=1,opt_qtl
                allocate (lrtsol(ic,iq)%list_ci(nargs))
            end do
        end do

        do i=1,nargs
            id_ci = get_int(values(i))
            if ( id_ci == DROP_OFF_CI ) then
                call drop_off(dataset,i,opt_calcul,opt_qtl,lrtsol)
            else if ( id_ci == BOOTSTRAP_FULL_CI ) then
                nameMethod = "boostrap 1"
                call boostrap(dataset,i,nsim,opt_calcul,opt_qtl,cli_haplo,lrtsol,listincsol,&
                    hdam,biq,interaction,create_dataset_full_sample,nameMethod)
            else if ( id_ci == BOOTSTRAP_SIB_CI ) then
                nameMethod = "boostrap 2"
                call boostrap(dataset,i,nsim,opt_calcul,opt_qtl,cli_haplo,lrtsol,listincsol,&
                    hdam,biq,interaction,create_dataset_sib_family_sample,nameMethod)
            else if ( id_ci == HENGDE_LI_CI ) then
                call hengde_li(dataset,i,nsim,opt_calcul,opt_qtl,cli_haplo,lrtsol,listincsol,&
                    hdam,biq,interaction)
            else
                call log_mess("confidence intervals ["//trim(str(id_ci))//"] does not exist!",WARNING_DEF)
            end if
        end do

    end subroutine
    !!***

    !!****f* m_qtlmap_calcul_ic/drop_off
    !!  NAME
    !!    drop_off
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!
    !!  RETURN
    !!
    !!  SOURCE
    subroutine drop_off(dataset,iarg,opt_calcul,opt_qtl,lrtsol)
        type(QTLMAP_DATASET)                                      , intent(in)  :: dataset
        integer                                                   , intent(in)  :: iarg,opt_calcul,opt_qtl
        type(TYPE_LRT_SOLUTION)     ,target   ,dimension(dataset%phenoModel%ncar,opt_qtl)   , intent(inout)  :: lrtsol

        !automatic/local
        type(TYPE_LRT_SOLUTION), pointer :: plrtsol
        integer                     :: i,posl,posr,iThres,icar,u,qtl,iq,hyp
        logical                     :: left,right
        real(kind=dp)               :: lrtThreshold,val
        real(kind=dp),dimension(3)  :: chisqArrayDfOne
        integer      ,dimension(3)  :: Average
        integer      ,dimension(opt_qtl) :: idx

                            ! p=0.05, p=0.01, p=0.001
        ! data (chisqArrayDfOne(i), i = 1, 3) /3.84, 6.64, 10.83/
                             ! p=0.05, p=0.025
        data (chisqArrayDfOne(i), i = 1, 3) /3.84, 5.02, 6.64 /
        data (Average(i), i = 1, 3) / 90, 95, 98 /

        call log_mess (" == drop off ==",INFO_DEF)

        do qtl=1,opt_qtl
            do icar=1,dataset%phenoModel%ncar
                plrtsol=>lrtsol(icar,qtl)
                if ( .not. associated(plrtsol%nxmax) ) cycle
                if ( .not. associated(plrtsol%chrmax) ) cycle
                if ( plrtsol%nxmax(0)>0 .and. plrtsol%chrmax(0)>0) then
                    allocate (plrtsol%list_ci(iarg)%ci_seuil(size(Average)))
                    allocate (plrtsol%list_ci(iarg)%ci_intervals(qtl,size(Average),qtl,2))
                    plrtsol%list_ci(iarg)%method="Drop off"


                    do iThres=1,size(Average)
                        plrtsol%list_ci(iarg)%nci=lrtsol(icar,qtl)%list_ci(iarg)%nci+1
                        u = plrtsol%list_ci(iarg)%nci
                        do hyp=1,qtl
                            do iq=0,qtl-1

                                lrtThreshold = plrtsol%lrtmax(iq) - chisqArrayDfOne(iThres)
                                left=.false.;right=.false.
                                idx(:qtl) = plrtsol%nxmax(0:qtl-1)

                                !find position uderthreshold on left
                                do i=plrtsol%nxmax(iq)-1,1,-1
                                    idx(iq+1) = i
                                    val = plrtsol%LRT%get(dataset,qtl,plrtsol%chrmax,idx,hyp)
                                    !on ne prend pas en compte les vraisemblances negatives
                                    if ( (val>0.d0) .and. (val < lrtThreshold) ) then
                                        left=.true.
                                        exit ! go outside the loop....
                                    end if
                                end do
                                posl=i
                                !find position uderthreshold on right
                                do i=plrtsol%nxmax(iq)+1,dataset%map%get_npo(plrtsol%chrmax(iq))
                                    idx(iq+1) = i
                                    val = plrtsol%LRT%get(dataset,qtl,plrtsol%chrmax,idx,hyp)
                                    if ( (val>0.d0) .and. (val < lrtThreshold) ) then
                                        right=.true.
                                        exit ! go outside the loop....
                                    end if
                                end do
                                posr=i

                                if (posl <= 0) posl=1
                                if (posr>size(dataset%map%absi,2)) posr=size(dataset%map%absi,2)

                                !print *,qtl,hyp,dataset%map%absi(plrtsol%chrmax(iq),posl)
                                plrtsol%list_ci(iarg)%ci_seuil(u) = Average(iThres)
                                plrtsol%list_ci(iarg)%ci_intervals(hyp,u,iq+1,1) = &
                                    dataset%map%absi(plrtsol%chrmax(iq),posl)
                                plrtsol%list_ci(iarg)%ci_intervals(hyp,u,iq+1,2) = &
                                    dataset%map%absi(plrtsol%chrmax(iq),posr)
                            !                print *,icar,chisqArrayDfOne(iThres),&
                            !                 dataset%map%absi(plrtsol%chrmax(0),posl),dataset%map%absi(plrtsol%chrmax(0),posr)
                            end do !iq
                        end do !hyp
                    end do !thres

                end if
            end do !ncar
        end do

    end subroutine drop_off
    !!***

    !!****f* m_qtlmap_calcul_ic/boostrap
    !!  NAME
    !!    boostrap
    !!  DESCRIPTION
    !!
    !!  INPUTS
    !!
    !!  RETURN
    !!
    !!  SOURCE
    subroutine boostrap(dataset,iarg,nsim,opt_calcul,opt_qtl,cli_haplo,lrtsolFromAnalyse,&
        listincsolFromAnalyse,hdam,biq,interaction,FUNCT_BOOSTRAP,nameMethod)
        !      use m_mrgrnk

        type(QTLMAP_DATASET)                                      , intent(in)  :: dataset
        integer                                                   , intent(in)  :: iarg,nsim,opt_calcul
        integer                                                   , intent(in)  :: opt_qtl,cli_haplo
        logical                                                   , intent(in)  :: biq,interaction
        type(TYPE_LRT_SOLUTION)        ,dimension(dataset%phenoModel%ncar,opt_qtl)   , intent(inout)  :: lrtsolFromAnalyse
        type(TYPE_INCIDENCE_SOLUTION)  ,dimension(dataset%phenoModel%ncar,opt_qtl+1) , intent(in)  :: listincsolFromAnalyse
        logical                                                   , intent(in)  :: hdam
        external                                                     :: FUNCT_BOOSTRAP
        character(len=LEN_DEF)                                    ,intent(in)   :: nameMethod

        integer                                      :: i,isim,iqtl,j,save_level,ncar
        type(QTLMAP_DATASET)                         :: newdataset
        type(PDD_BUILD)                              :: spt
        type(TYPE_LRT_SOLUTION)  , dimension(:,:),allocatable       :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION) ,dimension(:,:), allocatable  :: listincsol

        integer , parameter                                    :: sTab=1000
        integer , dimension(opt_qtl,opt_qtl,opt_qtl,dataset%phenoModel%ncar,sTab)  :: count_pos_max
        integer , dimension(sTab)                              :: RANK
        real    , dimension(sTab)      :: pos_max
        real (kind=dp),dimension (:,:,:),allocatable           :: rhoi
#ifdef MANAGE_CUDA
        type(TYPE_LRT_SOLUTION)  , dimension(:,:,:),allocatable      :: listLrtSolSimul
        real (kind=dp), dimension (:,:,:),allocatable                :: ySIMUL
        integer                                                      :: nsim_seuil
#endif
                                                   !10% (5 a gauche + 5 a droite), 5% et 1%
        real (kind=dp) ,dimension(3) :: quantile = (/ 0.05,0.025, 0.01 /)
        real :: sizeMap
        integer :: c,ic,xd,xg,id,q,prop,sxd,sxg,u,iq,qtl,hyp

        call log_mess (" == "//trim(nameMethod)//" ===",INFO_DEF)

        if (is_multitrait_analysis(opt_calcul)) then
            ncar = 1
        else
            ncar = dataset%phenoModel%ncar
        end if

        c=1
        sizeMap = dataset%map%posi(c,dataset%map%nmk(c))

        count_pos_max=0

        !calcul de la discretisation de la map...
        do i=1,sTab
            pos_max(i)=(real(i)*sizeMap/real(sTab))
        end do

        print *," == bootstrap / Confidence Intervals =="

        allocate(lrtsol(dataset%phenoModel%ncar,opt_qtl))
        allocate (listincsol(dataset%phenoModel%ncar,opt_qtl+1))
        allocate (rhoi(dataset%map%nchr,dataset%phenoModel%ncar,dataset%phenoModel%ncar))

        ! ON ENLEVE LES MESSAGES POUR LES SIMUL....
        save_level=get_log_level()
        call init_log(ERROR_DEF)

        do isim=1,nsim

            call newdataset%set()
            ! Bootstrap sample
            ! call create_dataset_full_sample(dataset,newdataset)
            ! call create_dataset_sib_family_sample(dataset,newdataset)
            call FUNCT_BOOSTRAP(dataset,newdataset)
            !Apply dataset
            call haplotype(newdataset,spt,cli_haplo)
#ifdef MANAGE_CUDA
            nsim_seuil=0
            !aucune simulation pour les seuil !
            allocate(listLrtSolSimul(newdataset%phenoModel%ncar,opt_qtl,nsim_seuil))
            allocate (ySIMUL(nsim_seuil,dataset%phenoModel%ncar,dataset%genea%nd))

            call analyse_cuda(newdataset,spt,nsim_seuil,SIMULATION,opt_calcul,opt_qtl,lrtSol,listincsol,&
                listLrtSolSimul,ySIMUL)

            deallocate (listLrtSolSimul,ySIMUL)
#else
            call analyse(newdataset,spt,opt_calcul,opt_qtl,lrtsol,listincsol,rhoi,SIMULATION)
#endif
            do qtl=1,opt_qtl
                do hyp=1,qtl
                    do iq=0,qtl-1
                        do i=1,ncar
                            id=nint((dataset%map%absi(lrtsol(i,qtl)%chrmax(iq),lrtsol(i,qtl)%nxmax(iq))*real(sTab))/sizeMap)
                            if (id == 0 ) cycle
                            !   if ( id > )
                            print *,qtl,iq,"isim:",isim,"icar:",i,&
                                " pos:",dataset%map%absi(lrtsol(i,qtl)%chrmax(iq),lrtsol(i,qtl)%nxmax(iq)),&
                                'lrt:',lrtsol(i,qtl)%lrtmax(iq),id,&
                                dataset%map%absi(lrtsol(i,qtl)%chrmax(iq),lrtsol(i,qtl)%nxmax(iq))

                            count_pos_max(qtl,hyp,iq+1,i,id)=count_pos_max(qtl,hyp,iq+1,i,id)+1
                        end do
                    end do
                end do
            end do


            !release data
            call release_dataset_boostrap(newdataset)

            call spt%release()
            do i=1,size(lrtsol,1)
                do j=1,size(lrtsol,2)
                    call lrtsol(i,j)%release()
                end do
            end do

            do i=1,size(listincsol,1)
                do j=1,size(listincsol,2)
                    call listincsol(i,j)%release()
                end do
            end do

        end do !isim
        do qtl=1,opt_qtl
            do ic=1,dataset%phenoModel%ncar
                allocate (lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%ci_seuil(size(quantile)))
                allocate (lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%ci_intervals(qtl,size(quantile),qtl,2))
                lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%nci = size(quantile)
                lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%method=nameMethod

                do q=1,size(quantile)
                    ! proportion de l'echantillong pour le quantile q
                    prop = nint(real(nsim)*quantile(q))
                    do hyp=1,qtl
                        do iq=0,qtl-1
                            ! on cherche la limite à droite
                            sxd=0
                            do xd=sTab,1,-1
                                sxd = sxd +  count_pos_max(qtl,hyp,iq+1,ic,xd)
                                if (sxd >= prop) exit
                            end do
                            sxg=0
                            do xg=1,sTab
                                sxg = sxg +  count_pos_max(qtl,hyp,iq+1,ic,xg)
                                if (sxg >= prop) exit
                            end do
                            ! corrections
                            if ( sxd > prop .and. xd < sTab ) xd = xd+1
                            if ( sxg > prop .and. xg > 1 ) xg = xg-1

                            lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%ci_seuil(q) = nint((1.d0 - 2*quantile(q))*100.d0)

                            lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%ci_intervals(hyp,q,iq+1,1) = pos_max(xg)
                            lrtsolFromAnalyse(ic,qtl)%list_ci(iarg)%ci_intervals(hyp,q,iq+1,2) = pos_max(xd)
                        end do !iq
                    ! print *,nint((1.d0 - 2*quantile(q))*100.d0),"%"," proportion:",prop,pos_max(ic,xg),pos_max(ic,xd)
                    end do !hyp
                end do !q
            end do !ic
        end do!qtl
        !  print *,count_pos_max(1,:)

        deallocate(lrtsol)
        deallocate (listincsol)
        deallocate (rhoi)

        call init_log(save_level)

    end subroutine boostrap
    !!***

    subroutine create_dataset_full_sample(dataset,newdataset)
        type(QTLMAP_DATASET)              ,intent(in)              :: dataset
        type(QTLMAP_DATASET)              ,intent(inout)           :: newdataset

        !local
        integer                :: id,i
        real                   :: r
        integer, dimension(dataset%genea%nd) :: array_sample

        do id=1,dataset%genea%nd
            call random_number(r)
            array_sample(id)=int(r*dataset%genea%nd)+1
        end do

        call create_dataset(dataset,array_sample,newdataset)

    end subroutine create_dataset_full_sample


    subroutine create_dataset_sib_family_sample(dataset,newdataset)
        type(QTLMAP_DATASET)              ,intent(in)              :: dataset
        type(QTLMAP_DATASET)              ,intent(inout)           :: newdataset

        !local
        integer                :: ic,ip,jm,kd,id,i
        real                   :: r
        integer, dimension(dataset%genea%nd) :: array_sample
        integer, dimension(dataset%genea%nd) :: half_sib
        integer                 :: nb_half_sib

        array_sample=0
        ic=1
        !do ic=1,dataset%phenoModel%ncar
        do ip=1,dataset%genea%np
            nb_half_sib=0
            do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
                if (dataset%phenoAnimal%estime(ic,jm)) then
                    !full sib
                    do id=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
                        call random_number(r)
                        array_sample(id)=int(r*(dataset%genea%ndm(jm+1)-dataset%genea%ndm(jm)))+dataset%genea%ndm(jm)+1
                    !    print *,id,r,int(r*(dataset%genea%ndm(jm+1)-dataset%genea%ndm(jm))),dataset%genea%ndm(jm)
                    end do
                else
                    do id=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
                        nb_half_sib=nb_half_sib+1
                        half_sib(nb_half_sib)=id
                    end do !id
                end if
            end do !jm
            !  print *,array_sample
            !   stop
            !sample for half sib
            do id=1,nb_half_sib
                call random_number(r)
                array_sample(half_sib(id))=half_sib(int(r*(nb_half_sib))+1)
            end do
        end do !ip

        !   print *,array_sample
        !  stop
        !end do !ic

        call create_dataset(dataset,array_sample,newdataset)

    end subroutine create_dataset_sib_family_sample

    subroutine create_dataset(dataset,array_sample,newdataset)
        type(QTLMAP_DATASET)              ,intent(in)              :: dataset
        integer ,dimension(dataset%genea%nd),intent(in)            :: array_sample
        type(QTLMAP_DATASET)              ,intent(inout)           :: newdataset

        integer , dimension(dataset%genea%nd) :: occurences

        call dataset%params%link(newdataset%params)
        call dataset%map%link(newdataset%map)
        !genealogy
        call create_dataset_genealogy(dataset,array_sample,newdataset,occurences)
        !genotype
        call create_dataset_genotype(dataset,occurences,array_sample,newdataset)
        !phenotype
        call create_dataset_phenotype(dataset,occurences,array_sample,newdataset)
           !call newdataset%phenoAnimal%write_file(newdataset,"test_perf")
           !call write_typ(newdataset,"test_typ")

    end subroutine create_dataset


    subroutine release_dataset_boostrap(dataset)
        type(QTLMAP_DATASET)              ,intent(inout)              :: dataset

        call dataset%genea%release()
        call dataset%phenoAnimal%release()
        call dataset%genoAnimal%release()

           !! la carte etait un link sur l ancienne carte donc pas de release !!

    end subroutine release_dataset_boostrap


    !!  DESCRIPTION
    !!    Implementation the Hengde Li method provide by the article "A quick method to calculate QTL confidence interval"
    subroutine hengde_li(dataset,iarg,nsim,opt_calcul,opt_qtl,cli_haplo,lrtsol,&
        listincsolFromAnalyse,hdam,biq,interaction)
        type(QTLMAP_DATASET)                                      , intent(in)  :: dataset
        integer                                                   , intent(in)  :: iarg,nsim,opt_calcul
        integer                                                   , intent(in)  :: opt_qtl,cli_haplo
        logical                                                   , intent(in)  :: biq,interaction
        type(TYPE_LRT_SOLUTION)        ,dimension(dataset%phenoModel%ncar,opt_qtl)   , intent(inout)  :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION)  ,dimension(dataset%phenoModel%ncar,opt_qtl+1) , intent(in)  :: listincsolFromAnalyse
        logical                                                   , intent(in)  :: hdam

        real (kind=dp) , dimension(:) , allocatable :: f,area
        real (kind=dp) :: thr,tol,alpha
        logical    :: loop
        real (kind=dp) :: lrtq,sArea,fd,tArea,lArea,rArea,dl,dr,val,denom
        integer :: ic,npos,ipos,iq,qtl,ncar,hyp,iqu
        real (kind=dp) ,dimension(3) :: quantile = (/ 0.05,0.025,0.01 /)
        integer ,dimension(0:opt_qtl-1) :: lpos,rpos,lchr

        real (kind=dp)   , parameter :: THRES_NUM = 1.d-5 ! seuil pour le calcul de surface et eviter les probleme numerique
        real (kind=dp)   , parameter :: STEP_DOWN = 0.01  ! pour calculer la surface en dessous de la courbre on descend un seuil (=> pas du seuil)

        call log_mess (" == Hengde Li ===",INFO_DEF)
        lpos=0
        rpos=0
        lchr=0
        if (is_multitrait_analysis(opt_calcul)) then
            ncar = 1
        else
            ncar = dataset%phenoModel%ncar
        end if

        do qtl=1,opt_qtl
            do ic=1,ncar
                allocate (lrtsol(ic,qtl)%list_ci(iarg)%ci_seuil(size(quantile)))
                allocate (lrtsol(ic,qtl)%list_ci(iarg)%ci_intervals(qtl,size(quantile),qtl,2))
                lrtsol(ic,qtl)%list_ci(iarg)%method="Hengde Li"
                lrtsol(ic,qtl)%list_ci(iarg)%nci=lrtsol(ic,qtl)%list_ci(iarg)%nci+size(quantile)
                do iqu=1,size(quantile)
                    lrtsol(ic,qtl)%list_ci(iarg)%ci_seuil(iqu) = nint((1.d0 - 2*quantile(iqu))*100)
                end do

                do hyp=1,qtl
                    do iq=0,qtl-1
                        lchr(0:qtl-1)=lrtsol(ic,qtl)%chrmax
                        npos = dataset%map%get_npo(lchr(iq))
                        allocate (f(0:npos),area(0:npos))
                        lpos(0:qtl-1) = lrtsol(ic,qtl)%nxmax
                        !Calculate RFR (Relative Frequence Ratio) at each scanned position
                        lrtq = lrtsol(ic,qtl)%lrtmax(iq)
                        f=0.d0
                        do ipos=1,npos
                            lpos(iq)=ipos
                            val = lrtsol(ic,qtl)%LRT%get(dataset,qtl,lchr,lpos,hyp)
                            if ( val > 0.d0 ) then
                                f(ipos) = exp(((val - lrtq)) * 0.5d0)
                            end if
                        end do
                        sArea = 0
                        area=0
                        do ipos=1,npos-1
                            !Calculate the Area Ai for position i to i+1
                            area(ipos) = dataset%map%absi(lchr(iq),ipos+1) - dataset%map%absi(lchr(iq),ipos)
                            area(ipos) = area(ipos) * ( f(ipos) + f(ipos+1) )
                            area(ipos) = area(ipos) * 0.5d0
                            if(area(ipos) /= area(ipos)) area(ipos) = 0.d0 !isnan...
                            sArea = sArea + area(ipos)
                        end do

                        !Rescale
                        area = area / sArea

                        do iqu=1,size(quantile)
                            alpha = quantile(iqu)
                            !Calculate CI
                            thr=lrtq
                            loop=.true.
                            tol = alpha / 100.d0
                            do while (loop)
                                !on augmente la surface, si la quantité deja calculé n'est pas suffisante pour atteindre le seuil iqu
                                thr = thr-STEP_DOWN
                                fd = exp((thr -lrtq)*0.5d0)
                                tArea=0;lArea=0;rArea=0;dl=0;dr=0

                                lpos(0:qtl-1) = lrtsol(ic,qtl)%nxmax
                                val = lrtsol(ic,qtl)%LRT%get(dataset,qtl,lchr,lpos,hyp)
                                do while ( (val < 0.d0) .or. (val >=thr))
                                    lpos(iq)=lpos(iq)-1
                                    if ( lpos(iq) == 0 ) exit
                                    val = lrtsol(ic,qtl)%LRT%get(dataset,qtl,lchr,lpos,hyp)
                                end do

                                rpos(0:qtl-1) = lrtsol(ic,qtl)%nxmax
                                val = lrtsol(ic,qtl)%LRT%get(dataset,qtl,lchr,rpos,hyp)
                                do while ( (val < 0.d0) .or. (val>=thr) )
                                    rpos(iq)=rpos(iq)+1
                                    if ( rpos(iq) > npos ) exit
                                    val = lrtsol(ic,qtl)%LRT%get(dataset,qtl,lchr,rpos,hyp)
                                end do
                                if ( lpos(iq)>0 ) then
                                    denom = f(lpos(iq)+1)-f(lpos(iq))
                                    if ( denom < THRES_NUM ) then
                                        dl = (dataset%map%absi(lchr(iq),lpos(iq)+1) - &
                                            dataset%map%absi(lchr(iq),lpos(iq)))*(fd - f(lpos(iq)))
                                    ! print *,"LEFT...."
                                    else
                                        dl = (dataset%map%absi(lchr(iq),lpos(iq)+1) - &
                                            dataset%map%absi(lchr(iq),lpos(iq)))*&
                                             (fd - f(lpos(iq)))/(f(lpos(iq)+1)-f(lpos(iq)))
                                    end if

                                    ! print *,"dl:",dl,lpos,fd,f(lpos(iq)),f(lpos(iq)+1),(fd - f(lpos(iq)))
                                    lArea = 0.5*dl*(f(lpos(iq))+fd)/SArea
                                end if
                                if ( lpos(iq) > 1 ) then
                                    !print *,"l:",sum(area(1:lpos(iq)-1)),lpos(iq)
                                    lArea = lArea + sum(area(1:lpos(iq)-1))
                                end if

                                if ( rpos(iq)<=npos .and. rpos(iq)>1 ) then
                                    denom = f(rpos(iq)-1)-f(rpos(iq))
                                    if ( denom < THRES_NUM ) then
                                        dr = (dataset%map%absi(lchr(iq),rpos(iq)) - &
                                            dataset%map%absi(lchr(iq),rpos(iq)-1))*(fd - f(rpos(iq)))
                                    !    print *,"RIGHT...."
                                    else
                                        dr = (dataset%map%absi(lchr(iq),rpos(iq)) - &
                                            dataset%map%absi(lchr(iq),rpos(iq)-1))*&
                                            (fd - f(rpos(iq)))/(f(rpos(iq)-1)-f(rpos(iq)))
                                    end if

                                    !print *,"dr:",dr,rpos," total:",npos,val
                                    rArea=0.5d0*dr*(f(rpos(iq))+fd)/Sarea
                                end if

                                if ( rpos(iq)>1 ) then
                                    !print *,"r:",sum(area(npos:rpos(iq)+1:-1)),rpos(iq)
                                    rArea = rArea + sum(area(npos:rpos(iq)+1:-1))
                                end if

                                tArea= lArea + rArea
                               ! print *,iqu," A:",alpha,tArea," B:", abs(tArea - alpha),tol,&
                               ! " C:",lArea ,rArea,lpos(iq),rpos(iq)
                                if ( (tArea <= alpha) .or. (abs(tArea - alpha)<= tol ) .or. thr < 0.d0 ) &
                                 loop = .false.
                            end do

                            !if ( qtl == 2 ) stop

                            if (lpos(iq)<1)lpos(iq)=1
                            if (rpos(iq)>npos)rpos(iq)=npos

                            lrtsol(ic,qtl)%list_ci(iarg)%ci_intervals(hyp,iqu,iq+1,1) = dataset%map%absi(lchr(iq),lpos(iq))+dl
                            lrtsol(ic,qtl)%list_ci(iarg)%ci_intervals(hyp,iqu,iq+1,2) = dataset%map%absi(lchr(iq),rpos(iq))-dr
                        end do !iqu
                        deallocate (f,area)
                    end do!iq
                end do !hyp
            end do
        end do !qtl

    end subroutine hengde_li
!!***

end module m_qtlmap_calcul_ic
