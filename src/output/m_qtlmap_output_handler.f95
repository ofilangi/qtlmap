!!****p* QTLMap/OUTPUT
!!  NAME
!!    INPUT
!!  DESCRIPTION
!!  Package output :
!!
!!  CREATION DATE
!!  01/01/2009
!!  COPYRIGHT
!!***

!!****m* OUTPUT/m_qtlmap_output_handler
!!  NAME
!!    m_qtlmap_output_handler -- Print analysisresult information on a ASCII format
!!  SYNOPSIS
!!
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
module m_qtlmap_output_handler
    use m_qtlmap_constant
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_base
    use m_qtlmap_types
    use m_qtlmap_math
    use m_qtlmap_haplotype_ldla
    use m_qtlmap_haplotype

    implicit none
    save
#ifdef GNU_COMP
    integer ,private , parameter              :: BUF_ALLOC_FILE=2**12
#elif INTEL_COMP
    integer ,private , parameter              :: BUF_ALLOC_FILE=2**10   ! Intel plante si >
#else
    ERREUR_DEFINED_KIND_OF_COMPILER => GNU_COMP or INTEL_COMP
    stop
#endif
    ! doit disparaitre....
    ! logical                      ,public      :: t_imp

    ! METTRE CES CONSTANTES EN PRIVE POUR INTERDIRE L UTILISATION EN DEHORS DE CE MODULE
    integer, public, parameter             :: nficout                = 16
    integer, public, parameter             :: nficerr                = 7
    integer, public, parameter             :: nficopti               = 18
    integer, public, parameter             :: nficpere               = 19
    integer, public, parameter             :: nficmere               = 20

    integer, public, parameter             :: unit_pded              = 1222
    integer, public, parameter             :: unit_pdedjoin          = 1223
    integer, public, parameter             :: unit_summary           = 9
    integer, public, parameter             :: unit_summary_2qtl      = 90
    integer, public, parameter             :: unit_coeff             = 80
    integer, public, parameter             :: unit_haplotypes        = 81
    integer, public, parameter             :: unit_phases_offspring  = 84
    integer, public, parameter             :: unit_freqall           = 83

    integer, public, parameter             :: MAX_GENOTYP_PRINT  = 15

    logical , private        :: xml       = .false.
    logical , private        :: file      = .true.


    public :: init_output_handler
    public :: end_output_handler
    public :: set_xml_output
    public :: set_file_output

    !General information
    !-------------------
    public :: log_descriptif_genealogy
    public :: log_marker_description
    public :: log_descriptif_traits
    public :: log_simulation_message

    !haplotype print
    !---------------
    public :: print_pded
    public :: print_phases
    public :: print_offspring_phase
    public :: print_freqall
    ! OFI - sept 2012 - informativity of each sire, of each marker
    public :: print_informativity_markers

    !Analyse print
    !-------------
    public :: print_start_multitraits
    public :: print_start_multitrait_DA
    public :: print_start_unitrait
    public :: print_end_multitraits
    public :: print_end_unitrait
    public :: print_courbe_LRT
    public :: print_maximum_LRT

    public :: print_transcriptome
    public :: print_transcriptome_Struct_famille


    public :: print_paternal_maternal_effect
    public :: print_pat_mat_effect_2QTL
    public :: print_LRT
    public :: print_allelic_origin
    public :: print_incidence_solution
    public :: print_incidence_solution_risk_factor
    public :: print_lrt_solution
    public :: print_confidence_intervals_solution

    !Simulation
    !-----------
    public :: print_resume_simulation
    public :: print_resume_simulation_2




contains

    !!****f* m_qtlmap_output_handler/init_output_handler
    !! NAME
    !!   init_output_handler
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine init_output_handler(dataset,opt_qtl,opt_calcul,name_funct,nb_thread)
        type(QTLMAP_DATASET)       ,intent(in)        :: dataset
        integer, intent(in)                           :: opt_qtl,opt_calcul,nb_thread
        character(len=300)         , intent(in)       :: name_funct
        integer :: ios,i
        integer :: val_recl = BUF_ALLOC_FILE ! anciennement 2048
        character(8)         :: date
        character(10)        :: time
        character(len=50)    :: valArg
        character(len=10000) :: com
        character(len=LENGTH_MAX_FILE)    :: result

        open(UNIT=nficout,file=dataset%params%get_file_val(K_OUTPUT), form="formatted",recl=val_recl,iostat=ios)
        if (ios/=0) then
            call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_OUTPUT)))
        end if

        call date_and_time(DATE=date,TIME=time)
        write (nficout,*) "   ***  "
        write (nficout,*) '       DATE      = ',date(1:4),'/',date(5:6),'/',date(7:8),'-',time(1:2),':',time(3:4),':',time(5:6)
#ifdef QTLMAP_VERSION
#ifdef DATE_BUILD
        write (nficout,*) '       Release-build      = ',QTLMAP_VERSION,'-',DATE_BUILD
#endif
#endif

        com="       ARGUMENTS   = "
        do i=1,COMMAND_ARGUMENT_COUNT()
            call GET_COMMAND_ARGUMENT(i,valArg)
            com=trim(com)//" "//trim(valArg)
        end do

        write (nficout,*) trim(com)
        !write (nficout,*) '       --QTL     = ',opt_qtl
        write (nficout,*) '       --CALCUL  = ',opt_calcul,' (',trim(name_funct),')'
        write (nficout,*) 'OMP_NUM_THREADS  = ',nb_thread

        open(UNIT=unit_summary,file=dataset%params%get_file_val(K_SUMM), form="formatted",recl=val_recl,iostat=ios)
        if (ios/=0) then
            call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_SUMM)))
        end if

        if ( trim(dataset%params%get_file_val(K_FREQALL))/='' ) then
            open(UNIT=unit_freqall,file=dataset%params%get_file_val(K_FREQALL),form="formatted",recl=val_recl,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_FREQALL)))
            end if
        end if

        if ( trim(dataset%params%get_file_val(K_HAPLOTYPES))/='' ) then
            open(UNIT=unit_haplotypes,file=dataset%params%get_file_val(K_HAPLOTYPES),form="formatted",recl=val_recl,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_HAPLOTYPES)))
            end if
        end if

    end subroutine init_output_handler
    !!***

    !!****f* m_qtlmap_output_handler/end_output_handler
    !! NAME
    !!   end_output_handler
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine end_output_handler(dataset,opt_qtl)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer, intent(in)     :: opt_qtl
        character(8)  :: date
        character(10) :: time

        call date_and_time(DATE=date,TIME=time)

        write (nficout,*) '    ***    '
        write (nficout,*) '    DATE      = ',date(1:4),'/',date(5:6),'/',date(7:8),'-',time(1:2),':',time(3:4),':',time(5:6)
        write (nficout,*) '    ***    '
        close(nficout)

    end subroutine end_output_handler
    !!***

    !!****f* m_qtlmap_output_handler/set_file_output
    !! NAME
    !!   set_file_output
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine set_file_output(active)
        logical, intent(in)     :: active
        file = active
    end subroutine set_file_output
    !!***

    !!****f* m_qtlmap_output_handler/set_xml_output
    !! NAME
    !!   set_xml_output
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine set_xml_output(active)
        logical, intent(in)     :: active
        xml = active
    end subroutine set_xml_output
    !!***

    !!****f* m_qtlmap_output_handler/print_summary_panalyse
    !! NAME
    !!   print_summary_panalyse
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_summary_panalyse(n,index_key,values)
        integer                           ,intent(in) :: n
        character(len=LEN_L) ,dimension(:),intent(in) :: index_key
        character(len=LEN_L) ,dimension(:),intent(in) :: values
        integer :: i
        call log_mess("",VERBOSE_DEF)
        call log_mess("-------------------------------",VERBOSE_DEF)
        call log_mess("********* P_ANALYSE KEYS ******",VERBOSE_DEF)
        call log_mess("-------------------------------",VERBOSE_DEF)

        write(nficout,FMT='(/,/,5x,'                 // &
            '"*****************  PARAMETERS ANALYSE FILE SUMMARY *****************",/)')
        do i=1,n
            write (nficout,FMT="(a40,'=',a)") index_key(i),values(i)
            call log_mess("["//trim(index_key(i))//"]-->"//trim(values(i)),VERBOSE_DEF)
        end do

    end subroutine print_summary_panalyse
    !!***

    !!****f* m_qtlmap_output_handler/log_descriptif_genealogy
    !! NAME
    !!   log_descriptif_genealogy
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine log_descriptif_genealogy(dataset)
        type(QTLMAP_DATASET)             , intent(in)  :: dataset
        integer                                        :: ip,ir,im

        character(len=LEN_DEF)                           :: unknown_r
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        write(nficout,FMT='(/,/,5x,'                 // &
            '"*****************  GENEALOGY DESCRIPTION *****************",/)')


        write (nficout,FMT='(1x,"The pedigree file includes ",i7," parents",'    //  &
            '" born from",i7," grand sires and ",i7," grand dams")') &
            size(dg%repro),size(dg%gpere),size(dg%gmere)

        write(nficout,FMT='(1x,"and ",i7," progeny born from",'           // &
            'i7," sires and ",i7," dams")')size(dg%animal),size(dg%pere),size(dg%mere)

        unknown_r=''
        do ip=1,size(dg%pere)
            do ir=1,size(dg%repro)
                if(dg%pere(ip) == dg%repro(ir)) then
                    exit
                endif
            end do
            if (ir > size(dg%repro) ) then
                unknown_r = unknown_r//' '//trim(dg%pere(ip))
            end if
        end do
        if (unknown_r /= '' ) then
            write(nficout,FMT='(1x,"Sires",a10," have no known ancestor")') trim(unknown_r)
        end if

        unknown_r=''
        do im=1,size(dg%mere)
            do ir=1,size(dg%repro)
                if(dg%mere(im) == dg%repro(ir)) then
                    exit
                endif
            end do
            if (ir > size(dg%repro) ) then
                unknown_r = unknown_r//' '//trim(dg%mere(im))
            end if
        end do

        if (unknown_r /= '' ) then
            write(nficout,FMT='(1x,"Dams ",a10," have no known ancestor")') trim(unknown_r)
        end if

    end subroutine log_descriptif_genealogy
    !!***

    !!****f* m_qtlmap_output_handler/log_marker_description
    !! NAME
    !!   log_marker_description
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine  log_marker_description(dataset)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer                                 :: i,j,k,jm,ip,igp,jgm
        integer                                 :: nmes0,l
        character(len=LEN_DEF) , dimension(:),allocatable    :: num_t
        type(GENEALOGY_BASE) , pointer :: dg
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dg => dataset%genea

        allocate( num_t (size(dga%numero)) )
        num_t=''
        !! compute number of animal which are not informative
        nmes0 = 0
        do i=1,size(dg%animal)
            if ( count(dga%presentg(:,i)) == 0 ) nmes0 = nmes0 + 1
        end do
        l = 0
        do i=1,dga%nmes
            do j=1,size(dg%animal)
                if (dga%numero(i) == dg%animal(j)) then
                    exit
                endif
            end do

            if ( j <= size(dg%animal) ) then
                cycle
            end if

            do jm=1,dg%nm
                if (dg%mere(jm) == dga%numero(i)) exit
            enddo
            if ( jm <= dg%nm ) cycle
            do ip=1,dg%np
                if (dg%pere(ip) == dga%numero(i)) exit
            enddo
            if ( ip <= dg%np ) cycle
            do  jgm=1,dg%ngm
                if(dg%gmere(jgm) == dga%numero(i)) exit
            enddo
            if ( jgm <= dg%ngm ) cycle
            do igp=1,dg%ngp
                if(dg%gpere(igp) == dga%numero(i)) exit
            enddo
            if ( igp <= dg%ngp ) cycle
            l = l + 1
            num_t(l)=dga%numero(i)
        end do

        write(nficout,FMT='(/,/,5x,'                 // &
            '"***************** MARKER DESCRIPTION *****************",/)')

        write (nficout,FMT='(i5, " animals are present in the genotype file ")') size(dg%animal)

        if (l > 0) then
            if ( l<50) then
                write(nficout,*) 'animal',(trim(num_t(i))//' ',i=1,l),'of genotype file ', &
                    'are not in the pedigree file'
            else
                write(nficout,*) l,' animals of genotype file are not in the pedigree file'
            end if

        end if

        deallocate( num_t )



        if (nmes0.eq.0) then
            write(nficout,*) 'where all animals are genotyped for at least ', &
                'one marker.'
        else
            write(nficout,*) 'where',nmes0,'animals have no genotyped marker .'
        endif

        write(nficout,FMT='("markers were selected among ",i4," markers")') sum(dataset%map%nmk)
        write(nficout,FMT='("There are ",i4," genotyped animals")') (size(dg%animal)-nmes0)

    end subroutine log_marker_description
    !!***

    !!****f* m_qtlmap_output_handler/print_freqall
    !! NAME
    !!   print_freqall
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_freqall(dataset)
        type(QTLMAP_DATASET)       ,intent(in) :: dataset
        integer :: i,j,k
        type(GENOTYPE_BASE) , pointer :: dga

        if ( trim(dataset%params%get_file_val(K_FREQALL))/='' ) then

            dga => dataset%genoAnimal

            do i=1,dataset%map%nchr ! chr
                do j=1,dataset%map%nmk(i) ! mark/chr
                    write(unit_freqall, &
                        FMT='(a," (Chr ",a,")"," has",i3," alleles:",90(" ",a,"(",f5.1,"%) ",a,"(",f5.1,"%)"))') &
                        trim(dataset%map%mark(i,j)),trim(dataset%map%chromo(i)),dga%nall(i,j), ( trim(dga%alleles(i,j,k)),&
                        dga%pc_all(i,j,k), k=1,dga%nall(i,j))
                end do
            end do
        end if

    end subroutine print_freqall
    !!***

    !!****f* m_qtlmap_output_handler/log_descriptif_traits
    !! NAME
    !!   log_descriptif_traits
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine log_descriptif_traits(dataset)
        type(QTLMAP_DATASET)                           :: dataset
        integer                                        :: ngeno,kd,ic,ifx,ilev,ico,alloc_stat
        integer                                        :: iperf,i,nnn,jfx,l
        real (kind=dp)                                 :: sy,sy2
        real (kind=dp), dimension (:), allocatable     :: covmu,covsig,covmin,covmax,ymax,ymin
        integer, dimension (:), allocatable            :: ncens,nmanq
        character(len=LEN_DEF) , dimension(size(dataset%phenoAnimal%bete)) :: manq
        type(PHENOTYPE_BASE)  ,pointer                 :: dpa
        type(GENEALOGY_BASE)  , pointer                :: dg
        type(DATAMODEL_BASE) , pointer                 :: dpm
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dpm => dataset%phenoModel
        dpa => dataset%phenoAnimal
        dg => dataset%genea

        write(nficout,FMT='(/,/,'    // &
            '"    ***************** TRAITS DESCRIPTION *****************",/)')!  // &
        !         '"The performance(s) file(s) is ",15(a10),/)') perfs

        l = 0

        do iperf=1,size(dpa%bete)
            do kd=1,size(dg%animal)
                if(dg%animal(kd) == dpa%bete(iperf)) exit
            enddo
            if ( kd <= size(dg%animal)) cycle
            l = l + 1
            manq(l)=dpa%bete(iperf)
        end do
        if (l > 0) then
            write(nficout,*) l,' animals with performance are not in the pedigree file'
        end if
        !! ---------  BUFFER ARRAY FOR LOG ----------------------------------------------------------------
        allocate ( covmu(dpm%ncov), STAT = alloc_stat )
        call check_allocate(alloc_stat,'covmu')
        allocate ( covsig(dpm%ncov), STAT = alloc_stat )
        call check_allocate(alloc_stat,'covsig')
        allocate ( covmin(dpm%ncov), STAT = alloc_stat )
        call check_allocate(alloc_stat,'covmin')
        allocate ( covmax(dpm%ncov), STAT = alloc_stat )
        call check_allocate(alloc_stat,'covmax')
        allocate ( ncens(dpm%ncar), STAT = alloc_stat )
        call check_allocate(alloc_stat,'ncens')
        allocate ( nmanq(dpm%ncar), STAT = alloc_stat )
        call check_allocate(alloc_stat,'nmanq')
        allocate ( ymax(dpm%ncar), STAT = alloc_stat )
        call check_allocate(alloc_stat,'ymax')
        allocate ( ymin(dpm%ncar), STAT = alloc_stat )
        call check_allocate(alloc_stat,'ymin')

        do i=1, dpm%ncar
            ncens(i)=0
            nmanq(i)=0
            ymax(i)=0.d0
            ymin(i)=INIFINY_REAL_VALUE
            do kd=1,dg%nd
                if((dpa%presentc(i,kd)) .and. dpa%ndelta(i,kd) == 0) ncens(i)=ncens(i)+1
                if(.not. dpa%presentc(i,kd)) nmanq(i)=nmanq(i)+1
                if(dpa%presentc(i,kd)) then
                    if(dpa%y(i,kd).GT.ymax(i)) ymax(i)=dpa%y(i,kd)
                    if(dpa%y(i,kd).LT.ymin(i)) ymin(i)=dpa%y(i,kd)
                endif
            enddo
        enddo

        do ico=1,dpm%ncov
            covmax(ico)=0.d0
            covmin(ico)=INIFINY_REAL_VALUE
            sy=0.d0
            nnn=0
            sy2=0.d0
            covmu(ico)=0.d0
            covsig(ico)=0.d0
            do kd=1,dg%nd
                if(dpa%covar(kd,ico) .ne. INIFINY_REAL_VALUE) then
                    if(dpa%covar(kd,ico).GT.covmax(ico)) covmax(ico)=dpa%covar(kd,ico)
                    if(dpa%covar(kd,ico).LT.covmin(ico)) covmin(ico)=dpa%covar(kd,ico)
                    sy=sy+dpa%covar(kd,ico)
                    nnn=nnn+1
                endif
            enddo
            covmu(ico)=sy/(dble(nnn))
            do kd=1,dg%nd
                sy2=sy2+((dpa%covar(kd,ico)-covmu(ico))*(dpa%covar(kd,ico)-covmu(ico)))
            enddo
            covsig(ico)=sqrt(sy2/(dble(nnn-1.d0)))
        enddo

        !!--------------------------------------------------------------------------

        ! IMPRESSION DU DESCRIPTIF DES DONNEES

        ngeno=0
        do kd=1,dg%nd
            do ic=1, size(dpa%presentc,1)
                !si il existe au moin un chromosome genotype pour l animal kd
                if ((count(dga%presentg(:,kd)) >= 1) .and. dpa%presentc(ic,kd))then
                    ngeno=ngeno+1
                    exit
                end if
            end do
        end do


        write(nficout,FMT='(/,5x,"NUMBER OF PHENOTYPED ANIMALS   : ",i5/,'               // &
            '5x,"NUMBER OF PHENOTYPED AND GENOTYPED ANIMALS : ",i5/,'              // &
            '5x,"NUMBER OF TRAITS               : ",i5/,'                          // &
            '5x,"NUMBER OF FIXED EFFECTS        : ",i5/,'                          // &
            '5x,"NUMBER OF COVARIABLES          : ",i5)') size(dpa%bete),ngeno, dpm%ncar, dpm%nfix, dpm%ncov


        !POUR LES COVARIABLES

        if( dpm%nfix /= 0 ) then
            do ifx=1,dpm%nfix
                write(nficout,*) 'FIXED EFFECT Num:', ifx,trim(dpm%namefix(ifx)),         &
                    'HAS', dpm%nlev(ifx),'LEVELS: ',                         &
                    (trim(dpm%listelev(ifx,ilev))//" ",ilev=1,dpm%nlev(ifx))
            end do
        end if

        if(dpm%ncov /= 0) then
            do ico=1,dpm%ncov
                write(nficout,FMT='(" COVARIABLE Num",i2,", ", a15,'                             // &
                    '" MEAN = ",f8.3,"+-",f8.3," (MIN=",f8.3,",MAX=",f8.3,")")')&
                    ico,trim(trim(dpm%namecov(ico))),covmu(ico),covsig(ico), covmin(ico),covmax(ico)
            end do
        end if

        !
        ! POUR CHAQUE CARACTERE
        !
        do ic=1, dpm%ncar
            if ((dpm%natureY(ic) == 'r' ).or.( dpm%natureY(ic) == 'a')) then
                write (nficout,FMT= '(/,1x,"TRAIT :",a6/, 5x                                    '  // &
                    ',"NUMBER OF PHENOTYPED PROGENY               : ",i5/, 5x   '  // &
                    ',"MEANS                                      : ",f8.3, "+-"'  // &
                    ',f8.3/, 5x,"MINIMUM                                    : " '  // &
                    ',f8.3/, 5x,"MAXIMUM                                    :",'  // &
                    'f8.3/ , 5x,"NUMBER OF MISSING PHENOTYPES               : ",'  // &
                    'i5)')&
                    trim(trim(dpm%carac(ic))),(dg%nd-nmanq(ic)),dpm%xmut(ic),dpm%sigt(ic),&
                    (ymin(ic)*dpm%sigt(ic)+dpm%xmut(ic)),&
                    (ymax(ic)*dpm%sigt(ic)+dpm%xmut(ic)),nmanq(ic)
            endif
            if ( dpm%natureY(ic) == 'i' ) then
                write (nficout,FMT= '(/,1x,"TRAIT :",a6/, 5x                                    '  // &
                    ',"NUMBER OF PHENOTYPED PROGENY               : ",i5/   '  // &
                    ',5x,"MINIMUM                                    : " '  // &
                    ',f8.3/, 5x,"MAXIMUM                                    :",'  // &
                    'f8.3/ , 5x,"NUMBER OF MISSING PHENOTYPES               : ",'  // &
                    'i5)')&
                    trim(trim(dpm%carac(ic))),(dg%nd-nmanq(ic)),ymin(ic),ymax(ic),nmanq(ic)

                write (nficout,FMT= '(/,1x,"  Class  ",3x," Frequency ")')

                do i=1,dpm%nmod(ic)
                    write (nficout,FMT= '(1x,i9,3x,f7.3)') dpm%indicemod(ic,i),dpm%prop(ic,i)
                end do

            end if

            !    if (opt_calcul == ANALYSE_UNITRAIT_MODLIN_COX) then
            write(nficout,FMT='(  5x,"NUMBER OF CENSORED PHENOTYPES              :",I5)') ncens(ic)
            !         end if

             !POUR LES COVARIABLES
            if  (dpm%modele(ic,1).ne.0.AND.dpm%modele(ic,2).eq.0.AND.dpm%modele(ic,3).eq.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namefix(dpm%modele(ic,3+ifx))),ifx=1,dpm%modele(ic,1))

            else if  (dpm%modele(ic,1).eq.0.AND.dpm%modele(ic,2).ne.0.AND.dpm%modele(ic,3).eq.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namecov(dpm%modele(ic,3+ico))),ico=1,dpm%modele(ic,2))

            else if  (dpm%modele(ic,1).ne.0.AND.dpm%modele(ic,2).ne.0.AND.dpm%modele(ic,3).eq.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namefix(dpm%modele(ic,3+ifx))),ifx=1,dpm%modele(ic,1)),      &
                    (' + ',trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ico))),ico=1,dpm%modele(ic,2))

            else if  (dpm%modele(ic,1).ne.0.AND.dpm%modele(ic,2).eq.0.AND.dpm%modele(ic,3).ne.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namefix(dpm%modele(ic,3+ifx))),ifx=1,dpm%modele(ic,1)),     &
                    (' + QTL*',trim(dpm%namefix(dpm%modele(ic,3+dpm%modele(ic,1)+dpm%modele(ic,2)+ifx))),ifx=1,dpm%modele(ic,3))

            else if  (dpm%modele(ic,1).ne.0.AND.dpm%modele(ic,2).ne.0.AND.dpm%modele(ic,3).ne.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namefix(dpm%modele(ic,3+ifx))),ifx=1,dpm%modele(ic,1)),     &
                    (' + ',trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ico))),ico=1,dpm%modele(ic,2)),(' + QTL*',            &
                    trim(dpm%namefix(dpm%modele(ic,3+dpm%modele(ic,1)+dpm%modele(ic,2)+jfx))),jfx=1,dpm%modele(ic,3))

            else if  (dpm%modele(ic,1).eq.0.AND.dpm%modele(ic,2).eq.0.AND.dpm%modele(ic,3).ne.0) then
                write(nficout,*) '     MODEL = mu',(' + QTL*',&
                    trim(dpm%namefix(dpm%modele(ic,3+dpm%modele(ic,1)+dpm%modele(ic,2)+jfx))), &
                    jfx=1,dpm%modele(ic,3))

            else if  (dpm%modele(ic,1).eq.0.AND.dpm%modele(ic,2).ne.0.AND.dpm%modele(ic,3).ne.0) then
                write(nficout,*) '     MODEL = mu',(' + ',trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ico))),         &
                    ico=1,dpm%modele(ic,2)),(' + QTL*',trim(dpm%namefix(dpm%modele(ic,3+dpm%modele(ic,1)+dpm%modele(ic,2)+ifx))), &
                    ifx=1,dpm%modele(ic,3))
            else
                write(nficout,*)'WITHOUT MODEL for fixed effects and covariables'
            end if
        end do

        deallocate ( covmu )
        deallocate ( covsig )
        deallocate ( covmin )
        deallocate ( covmax )
        deallocate ( ncens )
        deallocate ( nmanq )
        deallocate ( ymax )
        deallocate ( ymin )

    end subroutine log_descriptif_traits
    !!***

    !!****f* m_qtlmap_output_handler/log_simulation_message
    !! NAME
    !!   log_simulation_message
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine log_simulation_message(dataset,simulMap,h2,ue,dens,nalle,taille)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        logical , intent(in)          :: simulMap
        integer         , intent(in)  :: nalle
        real (kind=dp)  ,intent(in)   :: dens,taille
        real (kind=dp)  ,dimension(:),intent(in)   :: h2
        real (kind=dp)  ,dimension(:,:),intent(in)   :: ue

        integer                              :: ic,iq
        type(GENEALOGY_BASE) , pointer :: dg
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel
        dg => dataset%genea

        write(nficout,*)'SIMULATION FOR NEW MOLECULAR DATA IS ONLY ',&
            'AVAILABLE FOR REGULARLY SPACED MARKERS AND EQUAL',          &
            ' ALLELE FREQUENCIES'

        if ( simulMap ) then
            write(nficout,FMT="(1x,a8,'|',1x,f4.3,' M')") 'Densite ',dens
            write(nficout,FMT="(1x,a8,'|',2x,i3,' alleles')") 'Nb all  ',nalle
            write(nficout,FMT="(1x,a8,'|',2x,f3.1,' M')") 'LongChr ',taille
            write(nficout,FMT="(1x,a8,'|',2x,i3,' marq')") 'NbMarq  ',sum(dataset%map%nmk)
        end if

        if ( size(dpm%h2) <= 0 ) return

        if ( size(dpm%h2) == size(dpm%carac)) then
            write(nficout,*)
            write(nficout,*) 'Trait    heritabil effetsQTL'
            write(nficout,*)
            do ic=1,dpm%ncar
                write(nficout,'(1x,a5,1x,i2,1x,40(f8.3,2x))')'trait  ',ic,h2(ic), &
                    (ue(ic,iq), iq=1,size(ue,2))
            end do
        else
            write(nficout,*)
            write(nficout,*) 'Trait  effetsQTL'
            write(nficout,*)
            do ic=1,dpm%ncar
                write(nficout,'(1x,a5,1x,i2,1x,40(f8.3,2x))')'trait  ',ic, &
                    (ue(ic,iq), iq=1,size(ue,2))
            end do
        end if

    end subroutine log_simulation_message
    !!***

    !!****f* m_qtlmap_output_handler/print_start_multitraits
    !! NAME
    !!   print_start_multitraits
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_start_multitraits()

        write(nficout,1100)
1100    format(//9x,36('*')/9x,'*',34x,'*'/9x,            &
            '*  Joint analysis of the traits    *'/9x,    &
            '*',34x,'*'/9x,36('*')//)

    end subroutine print_start_multitraits
    !!***

    !!****f* m_qtlmap_output_handler/print_start_multitrait_DA
    !! NAME
    !!   print_start_multitrait_DA
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_start_multitrait_DA

        write(nficout,1100)
1100    format(//9x,36('*')/9x,'*',34x,'*'/9x,            &
            '*  Joint analysis of the traits    *'/9x,    &
            '*  using a discriminant  function  *'/9x,    &
            '*',34x,'*'/9x,36('*')//)

    end subroutine print_start_multitrait_DA
    !!***

    !!****f* m_qtlmap_output_handler/print_start_unitrait
    !! NAME
    !!   print_start_unitrait
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_start_unitrait(dataset,name_trait)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        character(len=LEN_DEF) ,intent(in)   :: name_trait
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        write(nficout,1000) trim(name_trait)

1000    format(//9x,36('*')/9x,'*',34x,'*'/9x,'*  Analysis of trait      ',a8,' *'/9x,'*',&
            34x,'*'/9x,'*',34x,'*'/9x,36('*')//)
        write(nficout,1001)dg%np,dg%nfem
1001    format(1x,'LRT profile on the linkage group :'/1x,' position, test statistic  , '/3x,i3,&
            ' sire QTL effects , '/3x,i3,' dam QTL effects')

    end subroutine print_start_unitrait
    !!***

    !!****f* m_qtlmap_output_handler/print_courbe_LRT
    !! NAME
    !!   print_courbe_LRT
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_courbe_LRT(dataset,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(TYPE_LRT_SOLUTION)       ,intent(in)    :: lrtsol

        integer                       :: ipos,chr,nlong,s,t
        real       ,dimension(:), allocatable    :: x_p,y_p
        real :: v

        if ( .not. lrtsol%LRT%list_available() ) then
            call log_mess("Devel error - LRT curves is not set !")
            return
        end if

        do chr=1,dataset%map%nchr
            nlong=dataset%map%get_npo(chr)
            s=0
            allocate (x_p(nlong),y_p(nlong))
            t=0
            do ipos=1,nlong
                s=s+1
                v = lrtsol%LRT%get1(dataset,chr,s)
                if (  v > 0.0d0 ) then
                    !      if(lrtsol%lrt1(chr,s) > 0.0d0) then
                    t=t+1
                    x_p(t) = dataset%map%absi(chr,s)
                    ! y_p(t) = lrtsol%lrt1(chr,s)
                    y_p(t) =  v
                end if
            end do
            write (nficout,*) "  ** chromosome  "//trim(dataset%map%chromo(chr))//" ** "
            call plott(y_p,x_p,t,nficout)
            deallocate (x_p,y_p)
        end do
    end subroutine print_courbe_LRT
    !!***

    !!****f* m_qtlmap_output_handler/print_transcriptome_H0
    !! NAME
    !!   print_transcriptome_H0
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_transcriptome_H0(dataset,listincsol)
        type(QTLMAP_DATASET)                           :: dataset
        type(TYPE_INCIDENCE_SOLUTION) , intent(in) ,dimension(size(dataset%phenoModel%carac))    :: listincsol

        type(TYPE_INCIDENCE_SOLUTION)      :: incsol
        character(len=LEN_LINE)                                 ::  title,likelyhood
        integer                                            :: ip,j,i,nbvalue,g,icar,ind,nbProfil,iprofile
        real (kind=dp)       , dimension(:)  ,allocatable         :: values
        integer , dimension(:),allocatable                        :: profilCar,profileNbValue
        character(len=LEN_LINE) , dimension(:),allocatable        :: profileTitle,profileFmt
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg

        dpm => dataset%phenoModel
        dg => dataset%genea

        allocate (profilCar(dpm%ncar))
        allocate (profileTitle(20)) ! maximum profile....
        allocate (profileFmt(20))
        allocate (profileNbValue(20))

        profilCar=0
        nbProfil = 0

        do icar=1,dpm%ncar
            nbvalue = 0
            title=""
            do i=1,size(listincsol(icar)%groupeName)
                if ( .not. listincsol(icar)%eqtl_print(i) ) cycle
                if ( listincsol(icar)%nbParameterGroup(i) == 1 ) then
                    title=trim(title)//"["//trim(listincsol(1)%groupeName(i))//"] "
                    nbvalue = nbvalue +1
                else
                    title=trim(title)//"[ *"//trim(listincsol(1)%groupeName(i))//"* "
                    do g=1,listincsol(icar)%nbParameterGroup(i)
                        title=trim(title)//trim(listincsol(1)%parameterName(i,g))
                        if ( g <  listincsol(icar)%nbParameterGroup(i) )title=trim(title)//","
                        nbvalue = nbvalue +1
                    end do
                    title=trim(title)//"] "
                end if
            end do

            !parcours des profil pour savoir si ce profil existe deja
            i=1
            do while ( i <= nbProfil )
                if ( trim(title) ==  profileTitle(i) ) then
                    exit
                end if
                i=i+1
            end do

            if ( i > nbProfil ) then ! ajout d un nouveau profil
                profileTitle(i)   = trim(title)
                profileFmt(i)     = "(1x,a20,"// trim(str(nbvalue+dg%np))//"f7.3)"
                profileNbValue(i) = nbvalue
                nbProfil = nbProfil + 1
            end if

            profilCar(icar)=i
        end do

        likelyhood=""
        likelyhood=trim(likelyhood)//trim('[ *std dev *')

        if ( size(listincsol)>0 ) then
            if ( associated(listincsol(1)%siga)) then
                likelyhood=trim(likelyhood)//trim(' *siga *')
            end if
        end if

        do i=1,dg%np
            likelyhood=trim(likelyhood)//trim(dg%pere(i))
            if ( i< dg%np) likelyhood=trim(likelyhood)//','
        end do
        likelyhood=trim(likelyhood)//']'

        do iprofile=1,nbProfil

            allocate(values(dg%np+1+profileNbValue(iprofile)))

            write(nficout,*) 'Profile    :',iprofile
            write(nficout,*) 'Hypothesis :0'
            write(nficout,*) 'Given parameters are respectively :'
            write(nficout,*) 'Gene position on the array, '//trim(likelyhood)//trim(profileTitle(iprofile))
            write(nficout,*)
            write(nficout,*) 'note : 0.0 means not estimable '
            write(nficout,*)


            values=0.0

            do icar=1,size(listincsol)
                if (profilCar(icar) /= iprofile ) cycle
                incsol=listincsol(icar)
                ind=0
                do i=1,dg%np
                    ind=ind+1
                    values(ind)=incsol%sig(1,i)
                end do

                ! variance animal feb 2012
                if ( associated(incsol%siga)) then
                    ind=ind+1
                    values(ind)=incsol%siga(1)
                end if

                do i=1,size(listincsol(icar)%groupeName)
                    if ( .not. listincsol(icar)%eqtl_print(i) ) cycle
                    if ( listincsol(icar)%nbParameterGroup(i) == 1 ) then
                        ind=ind+1
                        if ( .not. associated(incsol%parameterVecsol) ) then
                            values(ind)=0.d0
                        else
                            if ( incsol%parameterVecsol(i,1) ) then
                                values(ind)=incsol%paramaterValue(i,1)
                            else
                                values(ind)=0.d0
                            end if
                        end if
                    else
                        do g=1,listincsol(icar)%nbParameterGroup(i)
                            ind=ind+1
                            if ( .not. associated(incsol%parameterVecsol) ) then
                                values(ind)=0.d0
                            else
                                if ( incsol%parameterVecsol(i,g) ) then
                                    values(ind)=incsol%paramaterValue(i,g)
                                else
                                    values(ind)=0.d0
                                end if
                            end if
                        end do
                    end if
                end do
                write(nficout,FMT=profileFmt(iprofile))trim(dpm%carac(icar)),(values(i),i=1,ind)
            end do

            deallocate(values)
        end do


        deallocate (profilCar)
        deallocate (profileTitle)
        deallocate (profileFmt)
        deallocate (profileNbValue)

    end subroutine print_transcriptome_H0
    !!***

    !!****f* m_qtlmap_output_handler/print_transcriptome
    !! NAME
    !!   print_transcriptome
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_transcriptome(dataset,listlrtsol,listincsol)
        type(QTLMAP_DATASET)                           :: dataset
        type(TYPE_LRT_SOLUTION)       , intent(in) ,dimension(dataset%phenoModel%ncar)    :: listlrtsol
        type(TYPE_INCIDENCE_SOLUTION) , intent(in) ,dimension(dataset%phenoModel%ncar)    :: listincsol

        type(TYPE_LRT_SOLUTION)            :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION)      :: incsol
        character(len=LEN_LINE)                             :: fmt,title,likelyhood
        integer                                             :: ip,j,i,nbvalue,g,icar,ind,nbpos,nbProfil,iprofile
        real (kind=dp)       , dimension(:)  ,allocatable   :: values
        integer , dimension(:),allocatable                        :: profilCar,profileNbValue
        character(len=LEN_LINE) , dimension(:),allocatable        :: profileTitle,profileFmt
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg

        dpm => dataset%phenoModel
        dg => dataset%genea


        allocate (profilCar(dpm%ncar))
        allocate (profileTitle(20)) ! maximum profile....
        allocate (profileFmt(20))
        allocate (profileNbValue(20))

        profilCar=0
        nbProfil = 0

        do icar=1,dpm%ncar
            title=""

            ! on compte le nombre de valeur que doit contenir le tableau de transcrit
            ! + les position et les LRT selon l hypothese
            nbvalue = 0

            !construction du titre
            ! on prend comme reference le premier caractere pour l'estimabilite des effets
            do i=1,size(listincsol(icar)%groupeName)
                if ( .not. listincsol(icar)%eqtl_print(i) ) cycle
                if ( listincsol(icar)%nbParameterGroup(i) == 1 ) then
                    title=trim(title)//"["//trim(listincsol(icar)%groupeName(i))//"] "
                    nbvalue = nbvalue +1
                else
                    title=trim(title)//"[ *"//trim(listincsol(icar)%groupeName(i))//"* "
                    do g=1,listincsol(icar)%nbParameterGroup(i)
                        title=trim(title)//trim(listincsol(icar)%parameterName(i,g))
                        if ( g <  listincsol(icar)%nbParameterGroup(i) )title=trim(title)//","
                        nbvalue = nbvalue +1
                    end do
                    title=trim(title)//"] "
                end if
            end do

            !parcours des profil pour savoir si ce profil existe deja
            i=1
            do while ( i <= nbProfil )
                if ( trim(title) ==  profileTitle(i) ) then
                    exit
                end if
                i=i+1
            end do

            if ( i > nbProfil ) then ! ajout d un nouveau profil
                profileTitle(i) = trim(title)
                profileFmt(i)   = "(1x,a20,"//trim(str(listlrtsol(1)%hypothesis))//"(1x,a5,1x,f7.3)"//&
                    trim(str(listlrtsol(1)%hypothesis+nbvalue+dg%np))//"(f7.3,1x))"
                profileNbValue(i) = nbvalue
                nbProfil = nbProfil + 1
            end if

            profilCar(icar)=i
        end do

        likelyhood=""

        do i=1,listlrtsol(1)%hypothesis
            likelyhood=trim(likelyhood)//&
                'Chromosome '//trim(str(i))//', QTL Position '//trim(str(i))//','
        end do
        do i=1,listlrtsol(1)%hypothesis
            likelyhood=trim(likelyhood)//&
                'H'//trim(str(i-1))//"/H"//trim(str(listlrtsol(1)%hypothesis))//","
        end do


        likelyhood=trim(likelyhood)//'[ *std dev *'

        if ( size(listincsol)>0 ) then
            if ( associated(listincsol(1)%siga)) then
                likelyhood=trim(likelyhood)//trim(' *siga *')
            end if
        end if

        do i=1,dg%np
            likelyhood=trim(likelyhood)//dg%pere(i)
            if ( i< dg%np) likelyhood=trim(likelyhood)//','
        end do
        likelyhood=trim(likelyhood)//']'


        do iprofile=1,nbProfil

            allocate(values(dg%np+1+profileNbValue(iprofile)+3*size(dataset%map%absi,2)))

            write(nficout,*) 'Profile    :',iprofile
            write(nficout,*) 'Hypothesis :'//trim(str(listlrtsol(1)%hypothesis))
            write(nficout,*) 'Given parameters are respectively :'
            write(nficout,*) 'Gene position on the array, ',trim(likelyhood),trim(profileTitle(iprofile))
            write(nficout,*)
            write(nficout,*) 'note : 0.0 means not estimable '
            write(nficout,*)

            values=0.0

            do icar=1,size(listincsol)
                if (profilCar(icar) /= iprofile ) cycle
                incsol=listincsol(icar)
                lrtsol=listlrtsol(icar)
                ind=0

                do i=1,dg%np
                    ind=ind+1
                    values(ind)=incsol%sig(1,i)
                end do

                ! variance animal feb 2012
                if ( associated(incsol%siga)) then
                    ind=ind+1
                    values(ind)=incsol%siga(1)
                end if

                do i=1,size(listincsol(icar)%groupeName)
                    if ( .not. listincsol(icar)%eqtl_print(i) ) cycle
                    if ( listincsol(icar)%nbParameterGroup(i) == 1 ) then
                        ind=ind+1
                        if ( incsol%parameterVecsol(i,1) ) then
                            values(ind)=incsol%paramaterValue(i,1)
                        else
                            values(ind)=0.d0
                        end if
                    else
                        do g=1,listincsol(icar)%nbParameterGroup(i)
                            ind=ind+1
                            if ( incsol%parameterVecsol(i,g) ) then
                                values(ind)=incsol%paramaterValue(i,g)
                            else
                                values(ind)=0.d0
                            end if
                        end do
                    end if
                end do

                write(nficout,FMT=profileFmt(iprofile))trim(dpm%carac(icar)),&
                    (dataset%map%chromo(lrtsol%chrmax(i)),dataset%map%absi(lrtsol%chrmax(i),lrtsol%nxmax(i)),i=0,&
                    lrtsol%hypothesis-1),(lrtsol%lrtmax(i),i=0,lrtsol%hypothesis-1),(values(i),i=1,ind)
            end do

            deallocate(values)
        end do


        deallocate (profilCar)
        deallocate (profileTitle)
        deallocate (profileFmt)
        deallocate (profileNbValue)

        close(nficout)

        open(UNIT=nficout,file=dataset%params%get_file_val(K_OUTPUT), form="formatted",recl=BUF_ALLOC_FILE,position="append")

    end subroutine print_transcriptome
    !!***

    !!****f* m_qtlmap_output_handler/print_transcriptome_Struct_famille
    !! NAME
    !!   print_transcriptome_Struct_famille
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_transcriptome_Struct_famille(dataset)
        type(QTLMAP_DATASET)                    :: dataset
        integer  , dimension(:,:) , allocatable :: lestypes      ! correspondance type avec lmes caracteres associe a ce type
        integer  , dimension(:)   , allocatable :: nbcarbytype   ! nombre de carac dans le type
        logical  , dimension(:,:) , allocatable :: struct_fam    ! pour chaque type on stocke le profil c.a.d : presence du phenotype pour l individu kd
        integer                                 :: ntype = 0
        character(len=LEN_L) :: listUnknown
        integer :: kd,ic,typ
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        allocate (lestypes (dpm%ncar,dpm%ncar))
        allocate (nbcarbytype(dpm%ncar))
        allocate (struct_fam (dpm%ncar,dg%nd))
        nbcarbytype=0

        ntype=1

        struct_fam(ntype,:)=dpa%presentc(1,:)
        nbcarbytype(ntype)=1
        lestypes(ntype,nbcarbytype(ntype))=1


        do ic=2,dpm%ncar
            do typ=1,ntype
                if ( all(dpa%presentc(ic,:) .EQV. struct_fam(typ,:))) then
                    !on a trouve le meme profil
                    nbcarbytype(typ)=nbcarbytype(typ)+1
                    lestypes(typ,nbcarbytype(ntype))=ic
                    exit
                end if
            end do

            !new kind....
            if (typ>ntype) then
                ntype=ntype+1
                struct_fam(ntype,:)=dpa%presentc(ic,:)
                nbcarbytype(ntype)=1
                lestypes(ntype,nbcarbytype(ntype))=ic
            end if
        end do

        write (nficout,fmt="(2x,'TYPE',2x,'Transcript',2x,'Unknown',2x)")

        do typ=1,ntype
            listUnknown=''
            do kd=1,dg%nd
                if (.not. struct_fam(typ,kd)) then
                    if ( trim(listUnknown) == '') then
                        listUnknown=adjustl(dg%animal(kd))
                    else
                        listUnknown=trim(listUnknown)//','//adjustl(dg%animal(kd))
                    end if
                end if
            end do
            write (nficout,*) typ,lestypes(typ,:nbcarbytype(ntype)),trim(listUnknown)
        end do


        deallocate (nbcarbytype)
        deallocate (struct_fam)
        deallocate (lestypes)

    end subroutine print_transcriptome_Struct_famille
    !!***

    !!****f* m_qtlmap_output_handler/print_maximum_LRT
    !! NAME
    !!   print_maximum_LRT
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_maximum_LRT(dxmax,xlrmax)
        real (kind=dp)  ,intent(in)        :: dxmax,xlrmax

        write(nficout,3003)
3003    format(//1x,'Maximum likelihood ratio test :'/)
        write(nficout,3000)dxmax,xlrmax
3000    format(1x,'The maximum is reached at position ',f8.4,' M, with value ',f8.3/)

    end subroutine print_maximum_LRT
    !!***

    !!****f* m_qtlmap_output_handler/print_residual_correlation
    !! NAME
    !!   print_residual_correlation
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_residual_correlation(rhoi)
        real (kind=dp)  ,intent(in)  ,dimension(:,:) :: rhoi
        integer   :: ic,jc

        if (size(rhoi,1) /= size(rhoi,2)) then
            call stop_application("Devel error: print_residual_correlation")
        end if

        write(nficout,*)
        write(nficout,*) 'Residual Correlations '
        write(nficout,3015) (ic,ic=1,size(rhoi,1)-1)

        do ic=2,size(rhoi,1)
            write(nficout,3014) ic,(rhoi(ic,jc),jc= 1,ic-1)
        end do
3014    format(7x,i3,1x,30(f6.3,4x))
3015    format(10x,30(1x,i3,9x))

    end subroutine print_residual_correlation
    !!***

    !!****f* m_qtlmap_output_handler/print_coeff_linear_combination_max
    !! NAME
    !!   print_coeff_linear_combination_max
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_coeff_linear_combination_max(dataset,coeff_max)
        type(QTLMAP_DATASET)               ,intent(in)                      :: dataset
        real (kind=dp) ,intent(in) , dimension(dataset%phenoModel%ncar)     :: coeff_max
        integer   :: ic
        type(DATAMODEL_BASE) , pointer :: dpm

        dpm => dataset%phenoModel

        write(nficout,*)
        write(nficout,*) 'Coefficients of the linear  combination'
        write(nficout,3016) 'Trait ', (trim(dpm%carac(ic)),ic=1,size(dpm%carac))
        write(nficout,3017) (coeff_max(ic),ic=1,size(dpm%carac))

3016    format(1x,a6,1x,30(a10,4x))
3017    format(7x,30(4x,f8.4,2x))

    end subroutine print_coeff_linear_combination_max
    !!***

    !!****f* m_qtlmap_output_handler/print_coeff_linear_combination
    !! NAME
    !!   print_coeff_linear_combination
    !! DESCRIPTION
    !!
    !! SOURCE
    subroutine print_coeff_linear_combination(dataset,chr,ilong,npo,coeff)
        type(QTLMAP_DATASET)                    :: dataset
        integer, intent(in)         :: chr,npo,ilong
        real (kind=dp) ,intent(in) , dimension(size(dataset%phenoModel%carac),npo)     :: coeff
        integer   :: ic,ix,n,ios

        if (trim(dataset%params%get_file_val(K_COEFFDA)) == '' ) then
            return
        end if

        open(UNIT=unit_coeff, file=trim(dataset%params%get_file_val(K_COEFFDA))//trim(str(chr)),&
            form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
        if (ios/=0) then
            call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_COEFFDA)))
        end if

        write(unit_coeff,*) '# Coefficients of the linear  combination'
        write(unit_coeff,3018) '# Position ', (trim(dataset%phenoModel%carac(ic)),ic=1,size(dataset%phenoModel%carac))
        n=0
        do ix=0,ilong,dataset%map%pas
            n=n+1
            write(unit_coeff,3019) ix,(coeff(ic,n),ic=1,size(dataset%phenoModel%carac))
        end do

3018    format(1x,a9,1x,30(a10,4x))
3019    format(2x,i4,3x,30(4x,f8.4,2x))
        close(unit_coeff)

    end subroutine print_coeff_linear_combination


    subroutine print_confusion(dataset,hypothesis,nalert,alertCorrQtl,corrmax)
        type(QTLMAP_DATASET)             ,intent(in)            :: dataset
        integer                          ,intent(in)            :: hypothesis,nalert
        type(CORR_ALERT_TYPE) ,dimension(:,:)  ,intent(in)      :: alertCorrQtl
        real (kind=dp)                   ,intent(in)            :: corrmax
        character(len=LEN_L)  :: fmts,fmtd
        character(LEN_W)      :: last
        integer :: i

        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        write (nficout,fmt="(//,80('*')/'Test of confusion between QTL ','and other effects in the final constained model',/,"//&
            "'(test based on the correlation between columns of the ','incidence matrix)',//)")

        if (nalert == 0) then
            write(nficout,fmt="(/,' No confusion detected',/,' the highest correlation is : ', F7.3,/,80('*')/)") corrmax
            return
        end if

        fmts="('Risk for sire ',a,' of confusion between the QTL ',i3,' level',i3,' and :')"
        fmtd="('Risk for dam ',a,' of confusion between the QTL ',i3,' level',i3,' and :')"
        last=""
        do i=1,nalert
            if (alertCorrQtl(hypothesis,i)%ip > 0 ) then
                if (trim(last)/=trim(dg%pere(alertCorrQtl(hypothesis,i)%ip))) then
                    write(nficout,*)
                    write(nficout,fmt=fmts) trim(dg%pere(alertCorrQtl(hypothesis,i)%ip)),&
                        alertCorrQtl(hypothesis,i)%qtl,alertCorrQtl(hypothesis,i)%ntlev
                    last=trim(dg%pere(alertCorrQtl(hypothesis,i)%ip))
                end if
            end if
            if (alertCorrQtl(hypothesis,i)%jm > 0 ) then
                if (trim(last)/=trim(dg%mere(alertCorrQtl(hypothesis,i)%jm))) then
                    write(nficout,fmt=fmtd) trim(dg%mere(alertCorrQtl(hypothesis,i)%jm)),&
                        alertCorrQtl(hypothesis,i)%qtl,alertCorrQtl(hypothesis,i)%ntlev
                    last=trim(dg%mere(alertCorrQtl(hypothesis,i)%jm))
                end if
            end if
            write(nficout,*) alertCorrQtl(hypothesis,i)%name_effect,&
                alertCorrQtl(hypothesis,i)%name_level,alertCorrQtl(hypothesis,i)%corr
        end do
    end subroutine print_confusion
    !!**

    subroutine print_test_nuisances(ntest,listtestnuis)
        integer                                    ,intent(in)  :: ntest
        type(TEST_NUISANCES_TYPE) ,dimension(ntest),intent(in)  :: listtestnuis
        character(len=LEN_LINE)  :: fmt1,fmt2,myfmt,warn
        integer :: i
        if ( ntest == 0 ) return

        myfmt="(//,80('*')/'testing model effects',//)"
        write(nficout,fmt=myfmt)
        myfmt="('Tested effect     df.    Likelihood     p-value'/'                         ratio                 '/)"
        !header
        write(nficout,fmt=myfmt)

        fmt1="(a15,'(direct effect)',1x,i5,3x,f8.3,3x,f7.3)"
        fmt2="(a15,'  (intra qtl)  ',1x,i5,3x,f8.3,3x,f7.3)"

        do i=1,ntest
            if ( listtestnuis(i)%df == 0 ) then
                write(nficout,*) "The effect ["//trim(listtestnuis(i)%name)//"] might be confused with another effect !"
            end if
            if ( listtestnuis(i)%directeffect ) then
                write(nficout,fmt=fmt1) listtestnuis(i)%name,listtestnuis(i)%df,listtestnuis(i)%lrt,listtestnuis(i)%pvalue
            else
                write(nficout,fmt=fmt2) listtestnuis(i)%name,listtestnuis(i)%df,listtestnuis(i)%lrt,listtestnuis(i)%pvalue
            end if
        end do

        write(nficout,*)
        write(nficout,*) "When this probability exceeds the standard threshold corresponding to the 5, 1 or 0.1 Pent level",&
            ", you might consider removing this effect from the model"

    end subroutine print_test_nuisances
    !!***

    !!****f* m_qtlmap_output_handler/print_summary_analyse
    !! NAME
    !!   print_summary_analyse
    !! DESCRIPTION
    !!
    !! NOTE
    !!  A MODIFIER print_summary_analyse pour prendre en compte les effets d interaction QTL
    !! SOURCE
    subroutine print_summary_analyse(dataset,listlrtsol,listincsol,nqtl,starticar,endicar)
        type(QTLMAP_DATASET)              ,intent(in)                       :: dataset
        integer              , intent(in)                                   :: nqtl
        type(TYPE_LRT_SOLUTION)  , intent(in) ,dimension(dataset%phenoModel%ncar)       :: listlrtsol
        type(TYPE_INCIDENCE_SOLUTION) , intent(in) ,dimension(dataset%phenoModel%ncar)  :: listincsol
        integer                                 , optional                  :: starticar,endicar


        integer   :: ip,jm,ic,ifail,iq,i,ieff,ntlev,jef,nbtp,c,lp,effp(dataset%genea%np),efft,s,e
        real (kind=dp)    :: deffp,tst,prob
        character(len=4) :: ctest(size(dataset%genea%pere)*20,nqtl),tail,nqtlc
        character(len=LEN_LINE) :: fmt1,summary,sum2,sum3
        type(GENEALOGY_BASE) , pointer :: dg
        type(DATAMODEL_BASE) , pointer :: dpm
        type(PHENOTYPE_BASE) , pointer :: dpa

        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        s=1
        e=dpm%ncar
        if ( present(starticar) ) s = starticar
        if ( present(endicar) ) e = endicar

        summary=""
        sum2=""
        sum3=""
        nqtlc=trim(str(nqtl))
        do iq=1,nqtl
            summary=trim(summary)//trim(str(iq-1))//" QTL versus "//trim(nqtlc)//" QTL"
            if (iq < nqtl) summary=trim(summary)//","
            sum2=trim(sum2)//" "//trim(str(iq-1))//"/"//trim(str(nqtl))//"QTL"
            sum3=trim(sum3)//"   Chr              Pos"//trim(str(iq))
        end do

        write(unit_summary,*)
        write(unit_summary,*) '*************************************************************************************'
        write(unit_summary,*) 'Summary '//trim(summary)

        write(unit_summary,3017)(trim(dg%pere(ip)),ip=1,dg%np)
        write(unit_summary,3018) trim(sum2)//" "//trim(sum3),             &
            ( (' eff'//trim(str(iq)),iq=1,nqtl),' SD ', (' Student-Test'//trim(str(iq))//"(*)",iq=1,nqtl),ip=1,dg%np)
3017    format('Variable  N        Max Lik        Pos (M)    Sire',        &
            '          ',             50(1x,a12,25x))
3018    format(14x,26a,9x,50(28a,2x))

        fmt1="(a8,1x,i3,"//trim(nqtlc)//"(2x,f8.2),"//trim(nqtlc)//"(2x,a5,2x,f8.2),7x"//&
            trim(str(dg%np))//"( "//trim(nqtlc)//"(1x,f8.2),2x,f8.2,1x,"//trim(nqtlc)//"(1x,a4)))"

        do ic=s,e
            if (dpm%natureY(ic) /= 'r') cycle

            !OFI 02/09/2010 - le calcul des effectifs est dependant du caracteres, on ajoute dans la boucle interne de ic, le calcul
            !et on le retire des parametres de la procedure
            efft = 0
            do ip=1,dg%np
                effp(ip) = 0
                do jm=dg%nmp(ip)+1,dg%nmp(ip+1)

                    if ( dg%ndm(jm)+1 > dg%ndm(jm+1) ) cycle
                    if ( dg%ndm(jm)+1 <= 0 ) cycle

                    effp(ip) = effp(ip) + count(dpa%presentc(ic,dg%ndm(jm)+1:dg%ndm(jm+1)))
                end do
                efft = efft + effp(ip)
            end do

            ntlev=1
            nbtp = 3 + dpm%modele(ic,1)+dpm%modele(ic,2)!+modele(ic,3)
            do jef=1,dpm%modele(ic,3)
                ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
            end do
            if ( ntlev > 1 ) then
                write(unit_summary,*) " .....Summary do not describe interaction*qtl effect..... : trait "//trim(dpm%carac(ic))
                cycle
            end if

            ctest=' ns '
            c=0
            do ip=1,dg%np
                do lp=1,ntlev
                    c=c+1
                    ifail=0
                    deffp=dble(effp(ip)/2 -1)
                    do iq=1,nqtl
                        if ( .not. associated(listincsol(ic)%qtl_groupeName) ) then
                            call stop_application("Devel error : ** qtl_groupeName"//&
                            " in TYPE_INCIDENCE_SOLUTION is not initialized **")
                        end if
                        ieff=listincsol(ic)%qtl_groupeName(1,iq)
                        tst=sqrt(deffp)*dabs(listincsol(ic)%paramaterValue(ieff,c)*2.d0)/listincsol(ic)%sig(1,ip)
                        tail='U'
                        if (deffp.ge.1.d0) then
                            prob=MATH_QTLMAP_G01EBF(tail,tst,deffp, ifail)
                            if(0.05.gt.prob) ctest(c,iq)='sign'
                        else
                            ctest(ip,iq)='na'
                        end if
                    end do
                end do
            end do
            write(unit_summary,fmt=fmt1)trim(dpm%carac(ic)),efft,(listlrtsol(ic)%lrtmax(i),i=0,nqtl-1), &
                (dataset%map%chromo(listlrtsol(ic)%chrmax(i)),&
                dataset%map%absi(listlrtsol(ic)%chrmax(i),listlrtsol(ic)%nxmax(i)),i=0,nqtl-1), &
                !(dataset%map%absi(listlrtsol(ic)%chrmax(i),listlrtsol(ic)%nxmax(i)),i=0,nqtl-1),&
                (((listincsol(ic)%paramaterValue(listincsol(ic)%qtl_groupeName(1,iq),lp),&
                lp=(ntlev*(ip-1))+1,ip*ntlev),iq=1,nqtl),&
                listincsol(ic)%sig(1,ip),           &
                ((ctest(lp,iq),iq=1,nqtl),lp=(ntlev*(ip-1))+1,ip*ntlev),ip=1,dg%np)

        end do

        write(unit_summary,*)
        write(unit_summary,*)
        write(unit_summary,*) "   (*) Approximate test of the significance of the QTL effect within sire."
     !  close(unit_summary)

    end subroutine print_summary_analyse

    !!**
    subroutine print_end_multitraits

    end subroutine print_end_multitraits
    !!**
    subroutine print_end_unitrait

    end subroutine print_end_unitrait
    !!**

    subroutine print_paternal_maternal_effect(dataset,ic,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)    :: dataset
        integer                    ,intent(in)    :: ic
        type(TYPE_LRT_SOLUTION)    ,intent(in)    :: lrtsol

        integer    :: i,ii,npo,chr,ios
        character(len=100) :: FMT1,FMT2
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa

        integer              :: unit_paternal_effects  = 61
        integer              :: unit_maternal_effects  = 62


        if ( trim(dataset%params%get_file_val(K_PATEFF))/='' ) then
            open(UNIT=unit_paternal_effects,file=dataset%params%get_file_val(K_PATEFF),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_PATEFF)))
            end if

            write(unit_paternal_effects,*) ' ********************************************* '
            write(unit_paternal_effects,*) ' This file is invalid if interaction qtl case '
            write(unit_paternal_effects,*) ' ********************************************* '

        else
            unit_paternal_effects = 0
        end if

        if ( trim(dataset%params%get_file_val(K_MATEFF))/='' ) then
            open(UNIT=unit_maternal_effects,file=dataset%params%get_file_val(K_MATEFF),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_MATEFF)))
            end if

            write(unit_maternal_effects,*) ' ********************************************* '
            write(unit_maternal_effects,*) ' This file is invalid if interaction qtl case '
            write(unit_maternal_effects,*) ' ********************************************* '
        else
            unit_maternal_effects = 0
        end if


        dpa => dataset%phenoAnimal
        dg => dataset%genea

        if ( unit_maternal_effects == 0 .and. unit_paternal_effects == 0 ) then
            return
        end if

        FMT1="(1x,i3,f8.3,"// trim(str(dg%np)) //"(1x,f7.2))"
        if ( dpa%namest(ic) /= 0 ) then
            FMT2="(1x,i3,f8.3,"// trim(str(dpa%namest(ic))) //"(1x,f7.2))"
        else
            FMT2="(1x,i3,f8.3,1x,f7.2)"
        end if

        if ( unit_paternal_effects /= 0 ) &
            write(unit_paternal_effects,*)"  Chr  Pos  ",("   "// trim(dg%pere(ii)),ii=1,dg%np )
        if ( unit_maternal_effects /= 0 ) &
            write(unit_maternal_effects,*)"  Chr  Pos  ",("   "// trim(dg%mere(ii)),ii=1,dpa%namest(ic) )

        if (file) then
            do chr=1,dataset%map%nchr
                npo = dataset%map%get_npo(chr)
                do i=1,npo
                    if ( unit_paternal_effects /= 0 ) &
                        write(unit_paternal_effects,FMT=FMT1)&
                        chr, dataset%map%absi(chr,i),(lrtsol%pater_eff(chr,ii,i),ii=1,dg%np)
                    if ( unit_maternal_effects /= 0 ) &
                        write(unit_maternal_effects,FMT=FMT2)&
                        chr,dataset%map%absi(chr,i),(lrtsol%mater_eff(chr,ii,i),ii=1,dpa%namest(ic))
                end do
            end do
        end if
    end subroutine print_paternal_maternal_effect


    subroutine print_LRT(dataset,ic,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)         :: dataset
        integer                    ,intent(in)         :: ic
        type(TYPE_LRT_SOLUTION)    ,intent(inout)      :: lrtsol   ! attention lrtsol est modifié dans la procedure lrt<0 => init à 0

        integer    :: i,ii,j,npo,chr,ios
        character(len=100) :: FMT1,FMT2
        logical :: ok1,ok2
        type(GENEALOGY_BASE) , pointer :: dg
        integer     :: unit_sire_res          = 7
        integer     :: unit_dam_res           = 8
        integer :: nhyp
        integer,dimension(:,:),pointer :: listChr,listPos
        integer :: nGLTotal,ihyp,nestime
        integer ,dimension(:),pointer :: listm

        ok1=.true.;ok2=.true.

        if ( trim(dataset%params%get_file_val(K_LRTSIRE))/='' ) then
            open(UNIT=unit_sire_res,file=dataset%params%get_file_val(K_LRTSIRE),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_LRTSIRE)))
            end if
        else
            ok1=.false.
        end if

        if ( trim(dataset%params%get_file_val(K_LRTDAM))/='' ) then
            open(UNIT=unit_dam_res,file=dataset%params%get_file_val(K_LRTDAM), &
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_LRTDAM)))
            end if
        else
            ok2 = .false.
        end if

        dg => dataset%genea

        if ( .not. (ok1 .or. ok2)) then
            return
        end if

        nhyp=size(lrtsol%LRT%values,1)

        if ( ok1 ) then

            FMT1="(1x,i7,1x,"//trim(str(nhyp))//"(a,1x,f9.4,1x),"//trim(str(nhyp))&
            //"(f8.3,1x),"//trim(str(dg%np))//"(f8.3,1x))"

            write(unit_sire_res,&
             fmt="('IC',5x,"//trim(str(nhyp))//"(' Chr ',1x,' Pos ',1x),"//trim(str(nhyp))//"('LRT(H"&
              //trim(str(nhyp))//"/H',i1,')',2x),"//trim(str(dg%np))//&
               "(a,'(H"//trim(str(nhyp))//"/H"//trim(str(nhyp-1))//")',1x))")&
               (ihyp,ihyp=0,nhyp-1),(trim(dg%pere(ii)),ii=1,dg%np)
        end if

        if ( ok2 ) then

            call dataset%get_list_dam_estime(ic,listm,nestime)

            FMT2="(1x,i7,1x,"//trim(str(nhyp))//"(a,1x,f9.4,1x),"//trim(str(nestime))//"(1x,f8.3))"
            write(unit_dam_res,fmt="('IC',5x,"//trim(str(nhyp))//"(' Chr ',1x,' Pos ',1x),"//&
            trim(str(nestime))//"(1x,a))")&
             (trim(dg%mere(listm(ii))),ii=1,nestime )
        end if

         if (file) then
           if ( ok1 ) then
             ! On met à 0 les LRT qui n'ont pas de signification ! ticket #2178
             where ( lrtsol%LRT%values < 0)
               lrtsol%LRT%values = 0.d0
             end where
             do ii=1,dg%np
               where ( lrtsol%LRT_SIRES(ii)%values < 0)
               lrtsol%LRT_SIRES(ii)%values = 0.d0
             end where
             end do
           end if

           if ( ok2 ) then
            do ii=1,nestime
               where ( lrtsol%LRT_DAMS(listm(ii))%values < 0)
               lrtsol%LRT_DAMS(listm(ii))%values = 0.d0
             end where
             end do
           end if

           ! On va chercher la liste des positions calculées pour parcourir l ensembles des solutions
           call lrtsol%LRT%get_iterator_available_index(dataset,nhyp,listChr,listPos,nGLTotal)

           do i=1,nGLTotal
            if ( ok1 ) then
             write (unit_sire_res,fmt=FMT1) ic,&
             (trim(dataset%map%chromo(listChr(i,j))),dataset%map%absi(listChr(i,j),listPos(i,j)),j=1,nhyp),&
              (lrtsol%LRT%get(dataset,nhyp,listChr(i,:),listPos(i,:),ihyp),ihyp=1,nhyp),&
              (lrtsol%LRT_SIRES(ii)%get(dataset,nhyp,listChr(i,:),listPos(i,:),nhyp),ii=1,dg%np)
            end if
            if ( ok2 ) then
               write (unit_dam_res,fmt=FMT2) ic,&
               (trim(dataset%map%chromo(listChr(i,j))),dataset%map%absi(listChr(i,j),listPos(i,j)),j=1,nhyp),&
              (lrtsol%LRT_DAMS(listm(ii))%get(dataset,nhyp,listChr(i,:),listPos(i,:),nhyp),ii=1,nestime)
            end if

           end do
           deallocate (listChr,listPos)
           if ( ok2 ) deallocate (listm)
         end if

    end subroutine print_LRT
    !!***

    subroutine print_pat_mat_effect_2QTL(dataset,ic,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer , intent(in)                              :: ic
        type(TYPE_LRT_SOLUTION)         ,intent(in)       :: lrtsol

        integer    :: i,j,ip,jm,iq,init,chr,chr2,ios
        character(len=100) :: FMT1,FMT2
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa

        integer              :: unit_paternal_effects  = 61
        integer              :: unit_maternal_effects  = 62


        if ( trim(dataset%params%get_file_val(K_PATEFF))/='' ) then
            open(UNIT=unit_paternal_effects,file=dataset%params%get_file_val(K_PATEFF),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_PATEFF)))
            end if

            write(unit_paternal_effects,*) ' ********************************************* '
            write(unit_paternal_effects,*) ' This file is invalid if interaction qtl case '
            write(unit_paternal_effects,*) ' ********************************************* '

        else
            unit_paternal_effects = 0
        end if

        if ( trim(dataset%params%get_file_val(K_MATEFF))/='' ) then
            open(UNIT=unit_maternal_effects,file=dataset%params%get_file_val(K_MATEFF),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_MATEFF)))
            end if

            write(unit_maternal_effects,*) ' ********************************************* '
            write(unit_maternal_effects,*) ' This file is invalid if interaction qtl case '
            write(unit_maternal_effects,*) ' ********************************************* '
        else
            unit_maternal_effects = 0
        end if


        dpa => dataset%phenoAnimal
        dg => dataset%genea

        if (  unit_paternal_effects == 0  .and. unit_maternal_effects == 0) then
            return
        end if

        FMT1="(1x,i5,1x,i5,1x,f8.3,1x,f8.3,"// trim(str(dg%np*2)) //"(1x,f7.2))"
        if ( dpa%namest(ic) /= 0 ) then
            FMT2="(1x,i5,1x,i5,1x,f8.3,1x,f8.3,"// trim(str(dpa%namest(ic)*2)) //"(1x,f7.2))"
        else
            FMT2="(1x,i5,1x,i5,1x,f8.3,1x,f8.3,1x,f7.2)"
        end if

        if ( unit_paternal_effects /= 0 ) write(unit_paternal_effects,*)&
            "   Chr1    Chr2   Pos1  ","   Pos2  ",(( "   "// trim(dg%pere(ip))//"/Qtl["//&
            trim(str(iq))//"] ",iq=1,2),ip=1,dg%np)

        !il faut afficher que les mere estimable.....
        !  if ( trim(mateff)/='' ) write(unit_maternal_effects,*)&
        !  "   Chr1    Chr2   Pos1  ","   Pos2  ",(( "   "// trim(dg%mere(jm))//"/Qtl["//&
        !  trim(str(iq))//"] ",iq=1,2),jm=1,namest(ic))

        if (file) then
            do chr=1,dataset%map%nchr
                do i=1,dataset%map%get_npo(chr)-1
                    do chr2=chr,dataset%map%nchr
                        init=1
                        if ( chr2 == chr ) init=i+1
                        do j=init,dataset%map%get_npo(chr2)
                            if ( unit_paternal_effects /= 0 ) &
                                write(unit_paternal_effects,FMT=FMT1)&
                                chr,chr2,dataset%map%absi(chr,i),dataset%map%absi(chr2,j),&
                                ((lrtsol%pater_eff2(chr,chr2,ip,i,j,iq),iq=1,2),ip=1,dg%np)
                            if ( unit_maternal_effects /= 0 ) &
                                write(unit_maternal_effects,FMT=FMT2)&
                                chr,chr2,dataset%map%absi(chr,i),dataset%map%absi(chr2,j),&
                                ((lrtsol%mater_eff2(chr,chr2,jm,i,j,iq),iq=1,2),jm=1,dpa%namest(ic))
                        end do
                    end do
                end do
            end do
        end if
    end subroutine print_pat_mat_effect_2QTL


    subroutine print_phases(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)          :: dataset
        type(PDD_BUILD)            ,intent(inout)       :: spt

        integer :: ip,i,maxp,chr,ios,nombfem,jm,nd1,nd2,ngeno1,ngeno2,ifem,j,ii
        real(kind=dp) :: pb
        character(len=LEN_DEF)  :: sep,sep2,endline
        type(GENEALOGY_BASE) , pointer :: dg
        type(GENOTYPE_BASE) , pointer :: dga
        integer, parameter       :: unit_phases            = 82

        if ( trim(dataset%params%get_file_val(K_PHASES))/='' ) then
            open(UNIT=unit_phases,file=dataset%params%get_file_val(K_PHASES),&
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_PHASES)))
            end if
        else
            ! Pas de fichier definit pour mettre les phases, on sort...
            return
        end if

        dga => dataset%genoAnimal
        dg => dataset%genea

        do chr=1,dataset%map%nchr
            sep=' '
            sep2=' / '

            write(nficout,*)
            write(nficout,*) '    ***************** PARENTAL PHASES *****************'
            write(nficout,*) '                      FILE :',trim(dataset%params%get_file_val(K_PHASES))
            write(nficout,*) '    ***************************************************'

            write(unit_phases,*)
            write(unit_phases,*) '    ***************** SIRE PARENTAL PHASES *****************'
            write(unit_phases,*) '                      CHROMOSOME :',dataset%map%chromo(chr)
            write(unit_phases,*)

            do ip=1,dg%np
                if(spt%phasp(chr,ip)) then
                    write(unit_phases,*) trim(dg%pere(ip)),' s ',(trim(get_pheno(dga,chr,spt%genotyp(chr,i,dga%correp(ip),1)))&
                        //' ',i=1,dataset%map%nmk(chr))!,' / ',&
                    write(unit_phases,*) trim(dg%pere(ip)),' d ',(trim(get_pheno(dga,chr,spt%genotyp(chr,i,dga%correp(ip),2)))&
                        //' ',i=1,dataset%map%nmk(chr))
                else
                    write(unit_phases,*) trim(dg%pere(ip)),' ? ' ,(trim(get_pheno(dga,chr,spt%genotyp(chr,i,dga%correp(ip),1)))&
                        //' ',i=1,dataset%map%nmk(chr))!,' / ',&
                    write(unit_phases,*) trim(dg%pere(ip)),' ? ' ,(trim(get_pheno(dga,chr,spt%genotyp(chr,i,dga%correp(ip),2)))&
                        //' ',i=1,dataset%map%nmk(chr))
                end if

1000            format(/1x,'Sire ',a,' genotype (paternal / maternal phases): ')
1001            format(/1x,'Sire ',a,' genotype (unknown phase origins): ')

            end do
        end do

        nombfem=0
        do chr=1,dataset%map%nchr

            if ( MAX_GENOTYP_PRINT < dataset%map%nmk(chr) ) then
                maxp = MAX_GENOTYP_PRINT
                endline="..."
            else
                maxp = dataset%map%nmk(chr)
                endline=""
            end if

            endline=""

            write(unit_phases,*)
            write(unit_phases,*) '    ***************** DAM PARENTAL PHASES *****************'
            write(unit_phases,*) '                      CHROMOSOME :',dataset%map%chromo(chr)
            write(unit_phases,*)

            do jm=1,dg%nm
                nd1=dg%ndm(jm)+1
                nd2=dg%ndm(jm+1)
                ngeno1=spt%ngenom(chr,jm)+1
                ngeno2=spt%ngenom(chr,jm+1)
                ifem=dg%repfem(jm)

                if ( dga%estfem(ifem) ) then
                    nombfem=nombfem+1

                    if (spt%phasm(chr,jm)) then
                        write(unit_phases,2000)trim(dg%femelle(ifem)),(ngeno2-ngeno1+1)
                    else
                        write(unit_phases,1010)trim(dg%femelle(ifem)),(ngeno2-ngeno1+1)
                    end if

2000                format(/1x,'Dam ',a,' has ',i3,' likely genotypes',' (paternal / maternal phases):')
1010                format(/1x,'Dam ',a,' has ',i3,' likely genotypes',' (unknown phase origins) :')

                    do j=ngeno1,ngeno2
                        pb = spt%probg(chr,j)
                        write(unit_phases,*)(trim(get_pheno(dga,chr,spt%genotypm(chr,ii,j,1)))&
                            //' ',ii=1,dataset%map%nmk(chr)),&
                            trim(endline),' / ',       &
                            (trim(get_pheno(dga,chr,spt%genotypm(chr,ii,j,2)))//' ',ii=1,dataset%map%nmk(chr)),trim(endline),&
                            ' proba : ',pb
                    end do
                end if

            end do
        end do

        if(nombfem.eq.0) write(unit_phases,1020)

1020    format(/,'None of the females had more than the minimum number of progeny needed to estimate its possible phases')



    end subroutine  print_phases
    !!***

    !!****f* m_qtlmap_output_handler/print_informativity_markers
    !! NAME
    !!   print_informativity_markers
    !! DESCRIPTION
    !!   print in a file informativity Halfsib-Family / Linkage Group with different format file
    !!
    !!   By Sire - Sire - rate informativity (1>x>0)
    !!   By Linkage Group - nb sire informative (heterozygote) - rate informativity (1>x>0)
    !! NOTE
    !!
    !! SOURCE
    subroutine print_informativity_markers(dataset,file_informativity,countH,markH)
        type(QTLMAP_DATASET)                           , intent(in)   :: dataset
        character(len=LENGTH_MAX_FILE)                 , intent(in)   :: file_informativity
        real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr)            , intent(in)   :: countH
        real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)), intent(in)   :: markH

        integer :: ip,uni,ik,k,c,ios
        integer :: val_recl = BUF_ALLOC_FILE
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea

        uni=10000
        open(UNIT=uni,file=file_informativity, form="formatted",recl=val_recl,iostat=ios)
        if (ios/=0) then
            call stop_application("Can not open the file :"//trim(file_informativity))
        end if

        !Affichage de l'informativité global d'un pere pour un groupe de liaison donné

        do ip=1,dg%np
            do c=1,dataset%map%nchr
                write(uni,'(i7,i7,f7.3)') ip,c,countH(ip,c)
            end do
            write(uni,*) !obligatoire pour gnuplot
        end do
        ! 1 - Affichage du nombre de pere informatif le long des groupes de liaison
        ! 2 - Affichage de l'informativité de charque marqueur le long du groupe de liaison
        k=0
        do c=1,dataset%map%nchr
            do ik=1,dataset%map%nmk(c)
                k=k+1
                write (uni,'(i7,i7,i7,i7,f7.3)') c,ik,k,&
                    count( markH(:,c,ik)>0.d0),sum(markH(:,c,ik))/dg%np
            end do
        end do

        close (uni)

    end subroutine print_informativity_markers

    !!   print_offspring_phase
    subroutine print_offspring_phase(dataset,spt,c,mktot1,mktot2,lrtsol,opt_qtl,namefile,multi)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
        integer                               ,intent(in) :: opt_qtl
        type(TYPE_LRT_SOLUTION), dimension(dataset%phenoModel%ncar,opt_qtl),intent(in) :: lrtsol
        character(len=LENGTH_MAX_FILE) , intent(in)       :: namefile
        integer                               ,intent(in) :: c,mktot1, mktot2
        logical                               ,intent(in) :: multi
        !local
        integer                 :: ip,jm,kd,k,ios,nx

        Character(len=LEN_DEF)  :: hap_print(maxval(dataset%map%nmk),2)
        integer :: unit_p=888,i,ic,nqtl,iq,sm,em,c2,scar,flleft,flright
        character(len=LEN_LINE)  :: phrase
        Character(len=LEN_DEF)   :: h1,h2
        type(GENEALOGY_BASE) , pointer :: dg
        type(PHENOTYPE_BASE) , pointer :: dpa
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENOTYPE_BASE) , pointer :: dga
        Character(len=LEN_LINE)  :: string

        dga => dataset%genoAnimal
        dpm => dataset%phenoModel
        dg => dataset%genea
        dpa => dataset%phenoAnimal

        scar = dpm%ncar
        if ( multi) scar = 1

        open(UNIT=unit_p,file=namefile, form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
        if ( ios /= 0 ) then
            call stop_application("Can not create file:"//trim(namefile))
        end if

        do i=1,scar
            open(UNIT=unit_p+i,file=trim(namefile)//"_"//trim(dpm%carac(i)), form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
        end do

        write(unit_p,*)
        write(unit_p,*) '    ***************** OFFSPRING PHASES *****************'
        write(unit_p,*) '                      CHROMOSOME :',dataset%map%chromo(c)
        write(unit_p,*)
        write(unit_p,fmt="('    MARKER:',a,' (',f9.5,' M)',' ==> ','MARKER:',a,' (',f9.5,' M)')") &
            trim(dataset%map%mark(c,mktot1)),dataset%map%posi(c,mktot1),&
            trim(dataset%map%mark(c,mktot2)),dataset%map%posi(c,mktot2)

        do ic=1,scar
            write(unit_p+ic,*) '   ***************** OFFSPRING HAPLOTYPES ',&
                ' where the maximum likelihood was reached ***************** '
            string="ID Y"
            do nqtl=1,opt_qtl
              do iq=nqtl-1,0,-1
              string=trim(string)//" HAPSIRE(H"//trim(str(nqtl))//"/H"//trim(str(iq))//") "
              string=trim(string)//" HAPDAM(H"//trim(str(nqtl))//"/H"//trim(str(iq))//") "
              end do
            end do
            write(unit_p+ic,*) trim(string)
        end do

        do ip=1,dg%np
            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
                    hap_print='.'
                    IF(count(dga%presentg(:,kd))>0) then !condition pour sélectionner les individus génotypés
                        DO k=mktot1,mktot2 !boucle sur les marqueurs

                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            ! 3) CREATION D'un haplotype contenant des lettres en minuscule lorsqu'il est prédit pour les génotypes en charactère
                            !                                      le genotype suivit d'un p lorsqu'il est prédit pour les génotypes en chiffres
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            hap_print(k,1)=get_pheno(dga,c,spt%genotyp(c,k,dga%corred(kd),1))
                            hap_print(k,2)=get_pheno(dga,c,spt%genotyp(c,k,dga%corred(kd),2))
                            if (spt%reconstructed(c,dga%corred(kd),k)) then
                                hap_print(k,1)=trim(get_pheno(dga,c,spt%genotyp(c,k,dga%corred(kd),1)))//'p'
                                hap_print(k,2)=trim(get_pheno(dga,c,spt%genotyp(c,k,dga%corred(kd),2)))//'p'
                            !write(unit_p,*) 'marqueur',k,' ', ph,' ', ph1,' ', ph2
                            endif
                        !! A mettre dans m_qtlmap_output_handler par olivier -->CMO
                        ENDDO ! boucle marker mktot1 mktot2


                        ! Ajout OFI
                        ! On ecrit dans NCAR Fichiers : ID progeny, PERF, HAPLOTYPES
                        do ic=1,scar
                            if (.not. dpa%presentc(ic,kd)) cycle
                            phrase=trim(dg%animal(kd))//" "//trim(str(dpa%y(ic,kd)))
                            do nqtl=1,opt_qtl
                                do iq=nqtl-1,0,-1
                                    c2=lrtsol(ic,nqtl)%chrmax(iq)
                                    nx=lrtsol(ic,nqtl)%nxmax(iq)
                                    call dataset%map%get_flanking_marker(c2,nx,flleft,flright)
                                    sm=flleft-(int(dataset%params%LONGHAP/2)-1)
                                    em=flright+(int(dataset%params%LONGHAP/2)-1)
                                    if (sm<1) sm=1
                                    if (em>dataset%map%nmk(c2)) em = dataset%map%nmk(c2)
                                    h1='';h2=''
                                    do k=sm,em
                                        h1=trim(h1)//trim(get_pheno(dga,c,spt%genotyp(c2,k,dga%corred(kd),1)))
                                        h2=trim(h2)//trim(get_pheno(dga,c,spt%genotyp(c2,k,dga%corred(kd),2)))
                                        if (spt%reconstructed(c2,dga%corred(kd),k)) then
                                            h1=trim(h1)//'p'
                                            h2=trim(h2)//'p'
                                        end if
                                    end do
                                    phrase=trim(phrase)//" "//trim(h1)//"/"//trim(h2)
                                end do !iq
                            end do !nqtl

                            write(unit_p+ic,*) trim(phrase)
                        end do !ic
                        !! A mettre dans m_qtlmap_output_handler par olivier -->CMO
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        ! 4) IMPRESSION DES PHASES COMPLETES
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        write(unit_p,*) &
                            trim(dg%animal(kd)),' ',trim(dg%pere(ip)),' ',trim(dg%mere(jm)), ' ','s',' ',  &
                            (trim(hap_print(k,1))//' ',k=mktot1, mktot2) !, &
                        write(unit_p,*) &
                            trim(dg%animal(kd)),' ',trim(dg%pere(ip)),' ',trim(dg%mere(jm)), ' ','d',' ',  &
                            ( trim(hap_print(k,2))//' ',k=mktot1, mktot2) !, &
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        do i=0,scar
            write(unit_p+i,*)
            write(unit_p+i,*) 'When the phase is not found, it is predicted using the closest',&
                ' flanking markers with known phase. If the probability of the phases',&
                ' is upper the threshold, the most likely marker genotype coming from',&
                ' the sire and the dam are noted followed by p (for predicted) '
            close(unit_p+i)
        end do

    end subroutine  print_offspring_phase

    !!***


    !!****f* m_qtlmap_output_handler/print_pded
    !! NAME
    !!   print_pded
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine print_pded(dataset,spt,chr)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)      ,intent(in)                  :: spt
        integer, intent(in)                               :: chr

        real (kind=dp)    :: app,apm,pp,pm,dx
        integer           :: ip,nm1,nm2,ngeno1,ngeno2,jm,geno,igeno,kd,kkd,n,nd1,nd2,ios
        type(GENEALOGY_BASE) , pointer :: dg
        type(GENOTYPE_BASE) , pointer :: dga
        logical :: ok1,ok2

        ok1=.false.;ok2=.false.
        dga => dataset%genoAnimal
        dg => dataset%genea

        if ( (trim(dataset%params%get_file_val(K_PDED)) == '' ) .and. ( trim(dataset%params%get_file_val(K_PDECPLE)) =='') ) then
            return
        end if

        if ( trim(dataset%params%get_file_val(K_PDED)) /= '' ) then
            open(UNIT=unit_pded, file=trim(dataset%params%get_file_val(K_PDED))//&
                '_Chr_'//trim(dataset%map%chromo(chr)), form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_PDED)))
            end if
            ok1=.true.
        end if

        if ( trim(dataset%params%get_file_val(K_PDECPLE)) /= '' ) then
            open(UNIT=unit_pdedjoin, file=trim(dataset%params%get_file_val(K_PDECPLE))//&
                '_Chr_'//trim(dataset%map%chromo(chr)), &
                form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
            if (ios/=0) then
                call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_PDECPLE)))
            end if
            ok2=.true.
        end if

        if ( ok1 ) &
            write(unit_pded,*) 'Position   Sire    Dam Dam_Phase Animal   p(2nd sire allele)   p(2nd dam allele) '

        if ( ok2 ) &
            write(unit_pdedjoin,*) 'Position   Sire    Dam  Dam_Phase Animal   p(Hs1/Hd1 )  p(Hs1/Hd2 )  p(Hs2/Hd1 )  p(Hs2/Hd2 ) '

        do ip=1,dg%np
            nm1=dg%nmp(ip)+1
            nm2=dg%nmp(ip+1)
            do jm=nm1,nm2
                ngeno1=spt%ngenom(chr,jm)+1
                ngeno2=spt%ngenom(chr,jm+1)
                do geno=ngeno1,ngeno2
                    igeno=ngeno2-geno+1
                    nd1=spt%ngend(chr,geno)+1
                    nd2=spt%ngend(chr,geno+1)
                    do kd=nd1,nd2
                        kkd=spt%ndesc(chr,kd)
                        do n=1,(dataset%map%get_npo(chr)-1)
                            dx=dataset%map%absi(chr,n)

                            if(.not.(dga%estfem(dg%repfem(jm))).or.dataset%params%opt_sib.eq.OPT_SIB_HS) then
                                pp=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n)
                                app=(pp+1.d0)/2.d0
                                if ( ok1 ) then
                                    write(unit_pded,333) dx,trim(dg%pere(ip)),trim(dg%mere(jm)),igeno,trim(dg%animal(kkd)),app,0.5d0
                                end if
                                !
                                ! Impression des probabilites des couples d'haplo parentaux
                                if ( ok2 ) then
                                    write(unit_pdedjoin,334) dx,trim(dg%pere(ip)),trim(dg%mere(jm)),igeno,trim(dg%animal(kkd)),&
                                        0.5d0*spt%pdd(chr,kd,1,n),0.5d0*spt%pdd(chr,kd,1,n),&
                                        0.5d0*spt%pdd(chr,kd,3,n),0.5d0*spt%pdd(chr,kd,3,n)
                                !332           format(2x,f8.3,1x,a12,1x,a12,1x,2(5x,f8.3))
                                end if

                            else    !OPT_SIB_FS

                                pp=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                                pm=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                                app=(pp+1.d0)/2.d0
                                apm=(pm+1.d0)/2.d0


                                !Impression des probabilites des haplo parentaux
                                if ( ok1 ) then
                                    write(unit_pded,333) dx,trim(dg%pere(ip)),trim(dg%mere(jm)),igeno,trim(dg%animal(kkd)),app,apm
                                end if

333                             format(2x,f8.3,1x,2(a12,1x),i5,1x,a12,1x,2(12x,f8.3))
                                !Impression des probabilites des couples d'haplo parentaux
                                if ( ok2 ) then
                                    write(unit_pdedjoin,334) dx,trim(dg%pere(ip)),trim(dg%mere(jm)),igeno,&
                                        trim(dg%animal(kkd)),spt%pdd(chr,kd,1,n),spt%pdd(chr,kd,2,n),&
                                        spt%pdd(chr,kd,3,n),spt%pdd(chr,kd,4,n)
                                end if

334                             format(2x,f8.3,1x,2(a12,1x),i5,1x,a12,1x,4(5x,f8.3))

                            end if

                        end do ! N
                    end do ! KD
                end do ! GENO
            end do  ! JM
        end do ! IP

        if ( ok1 ) close(unit_pded)
        if ( ok2 ) close(unit_pdedjoin)

    end subroutine print_pded
    !!***

    !!****f* m_qtlmap_output_handler/print_resume_simulation
    !! NAME
    !!   print_resume_simulation
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine print_resume_simulation(trait,test,ns,YMU1,SIG1,S21,S31,XMIN1,XMAX1,prob,z,minP,maxP)
        character(len=*),intent(in)              :: trait,test
        integer         ,intent(in)              :: ns
        real(kind=dp)   ,intent(in)              :: YMU1,SIG1,S21,S31,XMIN1,XMAX1
        real(kind=dp) ,dimension(16),intent(in)  :: prob,z
        integer         ,intent(in)              :: minP,maxP
        integer        :: ii
        write(nficout,*)
        write(nficout,4004) trim(trait)
4004    format(                                             &
            1X,'*---------------------------------------*'        &
            ,/,                                                  &
            1x,'           Variable ',a)

        write(nficout,4015) test
4015    format(                                      &
            1X,'*---------------------------------------*' &
            ,/,                                           &
            1x,'          Test ',a   &
            ,/,                                              &
            &1X,'*---------------------------------------*',/)

        WRITE(nficout,6103)ns,YMU1,SIG1,S21,S31,XMIN1,XMAX1
6103    FORMAT(                                      &
            1X,' Test statistic distribution  :',/,       &
            1X,'     Number of simulations : ',i6,/,      &
            1X,'     Mean                  : ',F12.5,/,   &
            1X,'     Standard deviation    : ',F12.5,/,   &
            1X,'     Skewness              : ',F12.5,/,   &
            1X,'     Kurtosis              : ',F12.5,/,   &
            1X,'     Minimum               : ',F12.5,/,   &
            1X,'     Maximum               : ',F12.5,/)

        write(nficout,4005)
4005    format(                                              &
            1X,'*--------------------------------------*'         &
            ,/,                                                   &
            1x,'| chromosome | genome     |  Threshold |'         &
            ,/,                                                   &
            1x,'|          level          |            |'         &
            ,/,                                                   &
            1X,'|--------------------------------------|')
        do ii=minP,maxP
            if (ii.eq.11)write(nficout,4007) 1-prob(ii),'chrom_level',z(ii)
            if (ii.eq.12)write(nficout,4007) 1-prob(ii),'     *     ',z(ii)
            if (ii.eq.13)write(nficout,4007) 1-prob(ii),'  nb_chrom ',z(ii)
            if (ii.lt.11.or.ii.gt.13)write(nficout,4007) 1-prob(ii),'           ',z(ii)
4007        format(1x,'|',2x,F7.4,3x,'|',a12,'|',2x,F8.2,2x,'|')
        end do
        write(nficout,*)'*--------------------------------------*'

    end subroutine print_resume_simulation
    !!***

    !!****f* m_qtlmap_output_handler/print_resume_simulation_2
    !! NAME
    !!   print_resume_simulation_2
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine print_resume_simulation_2(test,ncar,carac,z)
        character(len=*),intent(in)                   :: test
        integer         ,intent(in)                   :: ncar
        character(len=*),dimension(ncar),intent(in)   :: carac
        real(kind=dp) ,dimension(16,ncar),intent(in)  :: z

        integer    :: i


        write(nficout,4009) test
4009    format(/,                                                  &
            '    ',a8,'               p_value at                 '        &
            ,/,                                                           &
            'Trait                 chromosome level                ',/    &
            '              5%        1%        0.1%       ',/)

        do i=1,ncar
            write(nficout,4010)carac(i),z(10,i),z(11,i),z(14,i)
        end do
4010    format(a8,4x,F6.2,4x,F6.2,4x,F6.2,8x)

    end subroutine print_resume_simulation_2
    !!***

    !!****f* m_qtlmap_output_handler/print_allelic_origin
    !! NAME
    !!   print_allelic_origin
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine print_allelic_origin(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)      ,intent(in)                  :: spt

        type(GENEALOGY_BASE) , pointer :: dg

        integer :: i,g,ip,chr

        dg => dataset%genea

        write (nficout,*)
        do ip=1,dg%np
            write (nficout,fmt="(' Allelic origin for ',a16)") dg%pere(ip)
            do chr=1,dataset%map%nchr
                if ( spt%phasp(chr,ip)) then
                    write (nficout,fmt="('Chromosome ',a4,' : known')") dataset%map%chromo(chr)
                else
                    write (nficout,fmt="('Chromosome ',a4,' : unknown')") dataset%map%chromo(chr)
                end if
            end do
        end do
        write (nficout,*)
        write (nficout,*) 'NOTE: known allelic origin means QTL effect =  maternal - paternal allele effects'


    end subroutine print_allelic_origin
    !!***

    !!****f* m_qtlmap_analyse_unitrait/print_incidence_solution
    !! NAME
    !!    print_incidence_solution
    !! DESCRIPTION
    !!
    !! HISTORY
    !!  09/09/2010 . add information :
    !!               - mean of absolute value of substitution effect (within sires family)
    !!               - mean of part of s.d (within sire family)
    !! HISTORY
    !!  14/09/2010 : correction bug : pas de calcul wq si pas d effet qtl (LD analyse)
    !! SOURCE
    subroutine print_incidence_solution(dataset,incsol)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(TYPE_INCIDENCE_SOLUTION)     , intent(in)   :: incsol
        integer :: i,g,ip,chr,ic,l,iq,hy
        real(kind=dp) :: wq
        character(LEN=400) ,dimension(2) :: nhr
        type(DATAMODEL_BASE) , pointer :: dpm
        type(GENEALOGY_BASE) , pointer :: dg

        dpm => dataset%phenoModel
        dg => dataset%genea

        write(nficout,*) "---------------------------------------------------------------"
        write(nficout,*) 'Estimation of parameters under H'//trim(str(incsol%hypothesis))
        write(nficout,*) "---------------------------------------------------------------"
        write(nficout,*)
        write(nficout,*) 'Within sire standard deviation'
        do ic=1,size(incsol%sig,1)
            if ( size(incsol%sig,1) > 1 ) then
                write(nficout,*) " ** Trait ",trim(dpm%carac(ic))," **"
            end if
            do ip = 1,dg%np
                write(nficout,fmt="(' sire ',a, '  s.d. :',f10.3)") trim(dg%pere(ip)),incsol%sig(ic,ip)
                if ( associated (incsol%haplotypes) ) then
                    do hy=1,incsol%hypothesis
                        write(nficout,*)  'haplotype:',hy,'[',&
                            trim(incsol%haplotypes(ip,hy,1)), ', breed origin=',trim(incsol%races_haplotypes(ip,hy,1)),'/', &
                            trim(incsol%haplotypes(ip,hy,2)), ', breed origin=',trim(incsol%races_haplotypes(ip,hy,2)),']'
                        write(nficout,*)
                    end do
                end if
            end do
        end do

        write(nficout,*)

        ! variance animal feb 2012
        if ( associated (incsol%siga)  ) then
            write(nficout,fmt="(' genetic s.d ',f10.3)") incsol%siga(1)
        end if

        write(nficout,"(//,'  parameter    ','        estimable ?    value     ','precision'/)")

        do i=1,size(incsol%groupeName)
            write (nficout,*) incsol%groupeName(i)
            write (nficout,*)

            if ( incsol%nbParameterGroup(i) == 1 ) then
                if ( incsol%parameterVecsol(i,1) ) then
                    write (nficout,fmt="(a50,'  yes ',2f10.3,1x)") " " , incsol%paramaterValue(i,1),&
                        incsol%parameterPrecis(i,1)
                else
                    write (nficout,fmt="(a50,'  no ')") incsol%groupeName(i)
                end if
            else

                do g=1,incsol%nbParameterGroup(i)
                    if ( incsol%parameterVecsol(i,g) ) then
                        write (nficout,fmt="(a50,'  yes ',2f10.3,1x)") incsol%parameterName(i,g)&
                            , incsol%paramaterValue(i,g) ,  incsol%parameterPrecis(i,g)
                    else
                        write (nficout,fmt="(a50,'  no ')") incsol%parameterName(i,g)
                    end if
                end do
            end if
            write (nficout,*)
        end do
        !  write(nficout,*)' NOTE: known allelic origin means QTL effect =  maternal - paternal allele effects'

        if (associated(incsol%rhoi)) call print_residual_correlation(incsol%rhoi)

        if (incsol%hypothesis>0) then
            write (nficout,*) '                        ***                          '
            write (nficout,*) ' The mean of absolute value of substitution effect WQ (in std unit) ='
            write (nficout,*) ' -------------------------- '
            do iq=1,incsol%hypothesis
                if ( .not. associated(incsol%qtl_groupeName) ) then
                    call stop_application("Devel error*** incsol%qtl_groupeName is not allocated"//&
                        " [m_qtlmap_output_handler:print_incidence_solution] ")
                end if

                ! Moyenne des valeurs absolues des effets de substitution
                do ic=1,size(incsol%sig,1)
                    !          print *,iq,'QTLEFFECT:',incsol%qtl_groupeName(ic,iq)
                    g = incsol%qtl_groupeName(ic,iq) ! get the index of qtl position effect
                    if (g <= 0 .or. g > size(incsol%parameterVecsol,1)) cycle
                    wq = 0
                    l = 0
                    do ip=1,dg%np
                        !  print *,incsol%parameterName(g,ip),incsol%groupeName(g)
                        if (incsol%parameterVecsol(g,ip)) then
                            wq = wq + abs( incsol%paramaterValue(g,ip) / incsol%sig(ic,ip))
                            l = l + 1
                        end if
                    end do
                    wq = wq / l
                    if ( size(incsol%sig,1) > 1 ) then
                        write (nficout,fmt="(' | qtl ',i5, ' | wq :',f10.3,' |',a)") iq,(wq*2.d0),trim(dpm%carac(ic))
                    else
                        write (nficout,fmt="(' | qtl ',i5, ' | wq :',f10.3,' |')") iq,(wq*2.d0)
                    end if
                end do
            end do
            write (nficout,*) ' -------------------------- '
        end if

        close(nficout)
        open(UNIT=nficout,file=dataset%params%get_file_val(K_OUTPUT), form="formatted",recl=BUF_ALLOC_FILE,position="append")

    end subroutine print_incidence_solution


    subroutine print_marker_information_at_the_max(dataset,spt,lrtsol)

        type(QTLMAP_DATASET)       ,intent(in)       :: dataset
        type(PDD_BUILD)            ,intent(in)       :: spt
        type(TYPE_LRT_SOLUTION)    ,intent(in)       :: lrtsol

        integer :: ip,iq,chr,n,jhr,totaln,flleft,flright

        type(HAPLOTYPE_POSITION_BUILD)  :: shp
        character(LEN=100)  :: nhr(2)

        write(nficout,*)
        write(nficout,*) "       ****        "
        write(nficout,*) "Marker informativity at the maximum likelihood estimation"
        write(nficout,*)
        write(nficout,*) " 0 < informativity < 1 "

        totaln=0
        do chr=1,dataset%map%nchr
            totaln = totaln + dataset%map%get_ilong(chr)
        end do

        do iq=1,lrtsol%hypothesis
            chr = lrtsol%chrmax(iq-1)
            n = lrtsol%nxmax(iq-1)
            write(nficout,*)
            write(nficout,*) " -------------------------------------"
            write(nficout,fmt="(a,i3,' position=',f7.4)") " QTL ",iq,dataset%map%absi(chr,n)
            write(nficout,*)
            write(nficout,fmt='(a,i3,a,a,a)') "Chromosome number tested = " , chr, " [Chromosome=",trim(dataset%map%chromo(chr)),"]"
            write(nficout,fmt='(a,i3,a,i3)') "Position number tested = ",n,"/",totaln
            write(nficout,*)

            call dataset%map%get_flanking_marker(chr,n,flleft,flright)
            write(nficout,*)
            write(nficout,fmt="(a,a,a,f7.4)") "Left marker = ",trim(dataset%map%mark(chr,flleft)) &
                ,' position = ',dataset%map%posi(chr,flleft)
            write(nficout,fmt="(a,a,a,f7.4)") "Right marker = ",trim(dataset%map%mark(chr,flright))&
                ,' position = ',dataset%map%posi(chr,flright)

            call shp%set(dataset,spt)
            call shp%set_haplo_for_ldla(chr,dataset%map%absi(chr,n),n,.true.,.false.)

            do ip=1,dataset%genea%np
                write (nficout,*) "                   *** "
                write (nficout,*) "                   Sire  ",trim(dataset%genea%pere(ip))
                write (nficout,fmt="(a,f7.3)") "  - Informativity = ",get_informativity_position_sire(dataset,spt,ip,chr,n)
                nhr= 'unknown'
                do jhr=1,2
                    if(shp%num_haplo_pere(ip,jhr,1) /= 0) nhr(jhr) =shp%name_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
                end do
                write (nficout,*) " - Haplotype     = [",trim(nhr(1)) ,"/", trim(nhr(2)), ']'
            end do
        end do
        write(nficout,*)


    end subroutine print_marker_information_at_the_max

    subroutine print_incidence_solution_risk_factor(dataset,incsol)
        type(QTLMAP_DATASET)       ,intent(in)           :: dataset
        type(TYPE_INCIDENCE_SOLUTION)     , intent(in)   :: incsol

        integer :: i,g,ip,chr,ic

        write(nficout,*) "---------------------------------------------------------------"
        write(nficout,*) 'Estimation of parameters under H'//trim(str(incsol%hypothesis))
        write(nficout,*) "---------------------------------------------------------------"
        write(nficout,*)
        write(nficout,*)

        write(nficout,"(//,'  parameter    ','        estimable ?     risk factor '/)")

        do i=1,size(incsol%groupeName)
            write (nficout,*) incsol%groupeName(i)
            write (nficout,*)

            if ( incsol%nbParameterGroup(i) == 1 ) then
                if ( incsol%parameterVecsol(i,1) ) then
                    write (nficout,fmt="(a25,'  yes ',2f10.3,1x)") " " , incsol%paramaterValue(i,1)
                else
                    write (nficout,fmt="(a25,'  no ')") incsol%groupeName(i)
                end if
            else

                do g=1,incsol%nbParameterGroup(i)
                    if ( incsol%parameterVecsol(i,g) ) then
                        write (nficout,fmt="(a25,'  yes ',f10.3,1x)") incsol%parameterName(i,g), incsol%paramaterValue(i,g)
                    else
                        write (nficout,fmt="(a25,'  no ')") incsol%parameterName(i,g)
                    end if
                end do
            end if
            write (nficout,*)
        end do

        !write(nficout,*)' NOTE: known allelic origin means QTL effect =  maternal - paternal allele effects'
        close(nficout)
        open(UNIT=nficout,file=dataset%params%get_file_val(K_OUTPUT), form="formatted",recl=BUF_ALLOC_FILE,position="append")

    end subroutine print_incidence_solution_risk_factor
    !!***

    !!****f* m_qtlmap_output_handler/print_lrt_solution
    !! NAME
    !!   print_lrt_solution
    !! DESCRIPTION
    !!
    !! NOTE
    !!
    !! SOURCE
    subroutine print_lrt_solution(dataset,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)       :: dataset
        type(TYPE_LRT_SOLUTION)       ,intent(in)    :: lrtsol
        integer :: iq,flleft,flright,i,j
        character(len=1) :: l

        if ( lrtsol%hypothesis <= 0 ) return
        write(nficout,*) "                      *****                    "
        if ( lrtsol%hypothesis == 1 ) then
            call print_courbe_LRT(dataset,lrtsol)
        end if

        write(nficout,fmt="(//1x,'Maximum likelihood ratio test :'/)")
        l='1'
        if (lrtsol%hypothesis >= 10 ) l='2'

        do iq=1,lrtsol%hypothesis
            write(nficout,fmt="('Test H',i"//l//",' / H',i"//l//",' : ',f9.5,/)") iq-1, lrtsol%hypothesis, lrtsol%lrtmax(iq-1)
        end do

        write(nficout,fmt="(1x,'The maximum is reached at position(s) ',"//trim(str(lrtsol%hypothesis))&
            //"(f9.4,'(Chr :',a4,') '))") &
            (dataset%map%absi(lrtsol%chrmax(iq),lrtsol%nxmax(iq)),&
            dataset%map%chromo(lrtsol%chrmax(iq)),iq=0,size(lrtsol%nxmax)-1)

        do iq=1,size(lrtsol%chrmax)
            call dataset%map%get_flanking_marker(lrtsol%chrmax(iq-1),lrtsol%nxmax(iq-1),flleft,flright)
            write (nficout,fmt="(1x,'flanking marker (qtl:',i5,') ,',a20,',',a20)") iq,&
                dataset%map%mark(lrtsol%chrmax(iq-1),flleft),dataset%map%mark(lrtsol%chrmax(iq-1),flright)
        end do

    end subroutine print_lrt_solution
    !!***
    subroutine print_confidence_intervals_solution(dataset,opt_qtl,lrtsol,multitraitanalysis)
        type(QTLMAP_DATASET)       ,intent(in)       :: dataset
        integer                    ,intent(in)       :: opt_qtl
        type(TYPE_LRT_SOLUTION)  ,dimension(dataset%phenoModel%ncar,opt_qtl)    ,intent(in)    :: lrtsol
        logical , intent(in)                         :: multitraitanalysis
        integer :: i,j,ic,qtl,chr,iq,iml,imr,hyp,scar
        real :: posi

        character(len=1) :: l

        scar = dataset%phenoModel%ncar

        if (multitraitanalysis) scar = 1

        write(nficout,*)
        write(nficout,*) " === Confidence Intervals === "

        do qtl=1,opt_qtl
            do ic=1,scar
                if ( associated(lrtsol(ic,qtl)%list_ci)) then
                    write(nficout,fmt="(a,i1,a)") "--------------------------- QTL = ",qtl,"------------------------------"
                    if ( .not. multitraitanalysis) then
                        write(nficout,*) " Trait ["//trim(dataset%phenoModel%carac(ic))//"]"
                    end if
                    write(nficout,fmt="(a14,a11,a11,a11,a11,a11,a20,a8,a20,a11)") "Name","Position","Method","Average",&
                        "Pos Left","Pos Right","Left flank marker","Pos","Right flank marker","Pos"
                    ! Pouchaque hypothese testé (ex : H2/H1, H2/H0)
                    do hyp=qtl,1,-1
                        write(nficout,fmt="(a,i1,a1,i1,a)") "                       ==  H ",qtl,"/",(hyp-1)," == "
                        ! Pouch chaque QTL du model, on affiche les infos
                        do iq=0,qtl-1
                            chr=lrtsol(ic,qtl)%chrmax(iq)
                            posi = dataset%map%absi(chr,lrtsol(ic,qtl)%nxmax(iq))
                            do i=1,size(lrtsol(ic,qtl)%list_ci)
                                do j=1,lrtsol(ic,qtl)%list_ci(i)%nci
                                    !search flanking marker at the left position
                                    do iml=2,dataset%map%nmk(chr)
                                        if ( dataset%map%posi(chr,iml) > lrtsol(ic,qtl)%list_ci(i)%ci_intervals(hyp,j,iq+1,1) ) then
                                            exit
                                        end if
                                    end do

                                    iml = iml-1

                                    do imr=dataset%map%nmk(chr)-1,1,-1
                                        if ( dataset%map%posi(chr,imr) < lrtsol(ic,qtl)%list_ci(i)%ci_intervals(hyp,j,iq+1,2) ) then
                                            exit
                                        end if
                                    end do

                                    imr = imr+1
                                    write(nficout,fmt="(a14,f9.3,a14,f9.1,'% ',f9.3,' ',f9.3,' ',a18,' ',f11.3,' ',a15,' ',f9.3)")&
                                        "QTL_"//trim(str(iq+1)),posi,&
                                        trim(lrtsol(ic,qtl)%list_ci(i)%method),&
                                        lrtsol(ic,qtl)%list_ci(i)%ci_seuil(j),lrtsol(ic,qtl)%list_ci(i)%ci_intervals(hyp,j,iq+1,1),&
                                        lrtsol(ic,qtl)%list_ci(i)%ci_intervals(hyp,j,iq+1,2),&
                                        trim(dataset%map%mark(chr,iml)),dataset%map%posi(chr,iml),&
                                        trim(dataset%map%mark(chr,imr)),dataset%map%posi(chr,imr)

                                end do
                                write(nficout,fmt="('                       ***                                      ')")
                            end do
                        end do ! iq
                    end do ! hyp
                end if

            end do
        end do

    end subroutine print_confidence_intervals_solution

    subroutine write_simulation_file(dataset,opt_qtl,is_multi,nsim,lrtsol)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer       , intent(in)     :: opt_qtl,nsim
        logical       , intent(in)     :: is_multi

        type(TYPE_LRT_SOLUTION)  , dimension(dataset%phenoModel%ncar,opt_qtl,nsim),intent(in)      :: lrtsol
        integer, parameter             :: unit_simula        = 14
        integer :: isim,i,ic,iiq,iq,ios,j,k
        character(len=LEN_LINE)                                 :: myfmt,mystr
        type(DATAMODEL_BASE) , pointer :: dpm

        mystr=""
        dpm => dataset%phenoModel

        if ( trim(dataset%params%get_file_val(K_OUTSIM))=='' ) then
            call log_mess(" * "//&
                " Can not write (Maximum LRT/Position) of each simulation "//&
                "(define the out_maxlrt key in the parameter file)*",WARNING_DEF)
            return
        end if

        call log_mess(" Write LRT/Position of each simulation in file:"//trim(dataset%params%get_file_val(K_OUTSIM)),INFO_DEF)

        open(unit_simula,file=dataset%params%get_file_val(K_OUTSIM),form="formatted",recl=BUF_ALLOC_FILE,iostat=ios)
        if (ios/=0) then
            call stop_application("Can not open the file :"//trim(dataset%params%get_file_val(K_OUTSIM)))
        end if

       myfmt="("
         do iq=1,opt_qtl
           myfmt=trim(myfmt)//"(f10.4"
           do iiq=1,iq
              myfmt=trim(myfmt)//",1x,a5,1x,f10.4"
           end do
           myfmt=trim(myfmt)//")"
        end do
        myfmt=trim(myfmt)//")"

        if (.not. is_multi ) then


            do ic=1,dpm%ncar
                write(unit_simula,*) "Trait ["//trim(dpm%carac(ic))//"]"

                do i=1,opt_qtl
                 !  do j=i-1,0,-1
                     mystr=trim(mystr)//" LRTMAX(H"//trim(str(i))//"/H"//trim(str(i-1))//")"
                     do k=1,i
                      mystr=trim(mystr)//" CHR"//trim(str(k))//" POS"//trim(str(k))
                     end do
                 !  end do
                end do
                write(unit_simula,fmt="(a)") trim(mystr)

                do isim=1,nsim
                    if (.not. associated(lrtsol(ic,1,isim)%chrmax) ) cycle
                    write(unit_simula,fmt=trim(myfmt)) &
                        (lrtsol(ic,iq,isim)%lrtmax(iq-1),(dataset%map%chromo(lrtsol(ic,iq,isim)%chrmax(iiq)),&
                        dataset%map%absi(lrtsol(ic,iq,isim)%chrmax(iiq),lrtsol(ic,iq,isim)%nxmax(iiq)),iiq=0,iq-1)&
                        ,iq=1,opt_qtl)
                end do
            end do
        else
            write(unit_simula,*) "* All Traits *"
            mystr=""
            do i=1,opt_qtl
               mystr=trim(mystr)//" LRTMAX(H"//trim(str(i))//"/H"//trim(str(i-1))//")"
               do k=1,i
                 mystr=trim(mystr)//" CHR"//trim(str(k))//" POS"//trim(str(k))
               end do
            end do

            write(unit_simula,fmt="(a)") trim(mystr)

            do isim=1,nsim
                write(unit_simula,fmt=trim(myfmt)) &
                    (lrtsol(1,iq,isim)%lrtmax(iq-1),(dataset%map%chromo(lrtsol(1,iq,isim)%chrmax(iiq)),&
                    dataset%map%absi(lrtsol(1,iq,isim)%chrmax(iiq),lrtsol(1,iq,isim)%nxmax(iiq)),iiq=0,iq-1)&
                    ,iq=1,opt_qtl)
            end do

        end if

        close(unit_simula)

    end subroutine write_simulation_file

    !! bilan_shared_haplo  informe sur le nombre et les carateristiques des haplotypes paternels
    !! retrouves chez les meres
    !! Utilise nmk,np,posi,correp,nm                         de m_qtlmap_data
    !! Utilise genotyp                                       de m_qtlmap_haplotype_data
    !! SOURCE
    subroutine bilan_shared_haplo(dataset,shp,c,dx,n,ip)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
        !
        integer       , intent(in) :: c,ip
        real(kind=dp) , intent(in) :: dx
        integer, intent (in) :: n
        integer itab,jtab,ktab,stat,imax1,imax2
        integer  , dimension (:,:)  , allocatable   :: tab
        type(MAP_BASE) , pointer    :: map

        map => dataset%map

        ! allocation
        !
        allocate (tab(2,maxval(map%nmk)))

        if(n ==1 .and. ip ==1)then
            write(nficout,20) ip
20          format(' Dam carried Sire ',i6,' haplotypes length distribution (grouped by 20 SNP)')
            write(nficout,30)
30          format(' position *  sire  * 1st haplotype / 2nd haplotype'//72('*')//&
                '          *        *  10  20  30  40  50  60  70  80  90 100 110 120',&
                ' 130 140 150 160 170 180 190 200.....'// 72('*'))
        end if

        tab=0
        jtab=0
        ktab=1
        do itab =1,map%nmk(c)
            jtab=jtab+1
            tab(1,ktab)=tab(1,ktab)+ shp%tab_shared_haplo(ip,1,1+itab)
            tab(2,ktab)=tab(2,ktab)+ shp%tab_shared_haplo(ip,2,1+itab)
            if(jtab == 10) then
                jtab=0
                ktab=ktab+1
            end if
        end do
        do imax1 =ktab,1,-1
            if(tab(1,imax1) /=0) exit
        end do
        do imax2 =ktab,1,-1
            if(tab(2,imax2) /=0) exit
        end do

        write(nficout,31)  dx,ip,(tab(1,itab),itab=1,imax1)
31      format(3x,f7.3,'  *',i8,'*',30(1x,i3))
        write(nficout,32)      (tab(2,itab),itab=1,imax2)
32      format(3x,5x,11x,'*',30(1x,i3))
        deallocate (tab)

    end subroutine bilan_shared_haplo
!!***


end module m_qtlmap_output_handler
