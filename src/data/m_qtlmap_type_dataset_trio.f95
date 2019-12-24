module m_qtlmap_type_dataset_trio
    use m_qtlmap_log
    use m_qtlmap_type_genealogy
    use m_qtlmap_type_phenotype

    implicit none

    integer           , parameter      ,private                 :: MAXCAR = 20

     type QTLMAP_TRIO
          ! max trait for QTLMAP_TRIO type description
          integer                                               :: idSire    = -1   ! sire
          integer                                               :: idDam     = -1   ! dam
          integer                                               :: idKd      = -1   ! Kd
          real(kind=dp)    ,dimension(MAXCAR)                   :: perfKd    = 0.d0
          logical          ,dimension(MAXCAR)                   :: presentc  =.false.

     end type QTLMAP_TRIO

    type DATASET_QTLMAP_TRIO
         integer                                                 :: na        = 0       ! number of animal (the population)
         character(LEN=LEN_DEF) ,dimension(:) ,pointer           :: listNA    => null() ! the key ID => Name for each animal of the population
         type(QTLMAP_TRIO)  ,dimension(:) ,pointer               :: trioKd    => null() ! List of all trio parent link
         integer           ,dimension(:) ,pointer                :: corranim  => null()

         integer          ,dimension(:,:),pointer                :: NCK       => null() ! nombre de descendant ayant une perf par KD pour le Caracteres i
         integer          ,dimension(:,:,:),pointer              :: NCiCjK    => null() ! nombre de descendant ayant une perf par KD pour le Caracteres i et j
         real(kind=dp)    ,dimension(:,:,:),pointer              :: corcd     => null()

         contains

         !Ajoute un trio Pere,Mere,Desc dans la structure generale dataset
         procedure, public :: add_animal_genea
         !Ajoute les performance et l info donne manquante pour chaque Kd
         procedure, public :: init_perf_animal
         !Donne la liste des descendants d'un Individu
         procedure, public :: get_listProgenies

         !Construit l'esperances des performances pour chaque KD de generation 2 ainsi que les CD
         !initilise Nck et Nc1c2k
         procedure,public :: calcul_y_cd

         procedure,public :: calcul_corcd

         procedure ,public :: copy    => copy_qtlmap_dataset_trio
         procedure ,public :: release => release_qtlmap_dataset_trio

    end type DATASET_QTLMAP_TRIO

    interface extend
        module procedure extendListNa
        module procedure extendListTrio
        module procedure extendListCorrAnim
    end interface extend

  contains

  subroutine copy_qtlmap_dataset_trio(this,copy)
    class(DATASET_QTLMAP_TRIO) , intent(in)    :: this
    type(DATASET_QTLMAP_TRIO) , intent(inout) :: copy

    copy%na = this%na

    if ( associated(this%listNA)) then
      allocate (copy%listNA(size(this%listNA)))
      copy%listNA = this%listNA
    end if

    if ( associated(this%trioKd)) then
      allocate (copy%trioKd(size(this%trioKd)))
      copy%trioKd = this%trioKd
    end if

    if ( associated(this%corranim)) then
      allocate (copy%corranim(size(this%corranim)))
      copy%corranim = this%corranim
    end if

    if ( associated(this%NCK)) then
      allocate (copy%NCK(size(this%NCK,1),size(this%NCK,2)))
      copy%NCK = this%NCK
    end if

    if ( associated(this%NCiCjK)) then
      allocate (copy%NCiCjK(size(this%NCiCjK,1),size(this%NCiCjK,2),size(this%NCiCjK,3)))
      copy%NCiCjK = this%NCiCjK
    end if

    if ( associated(this%corcd)) then
      allocate (copy%corcd(size(this%corcd,1),size(this%corcd,2),size(this%corcd,3)))
      copy%corcd = this%corcd
    end if


  end subroutine copy_qtlmap_dataset_trio


  subroutine release_qtlmap_dataset_trio(this)
     class(DATASET_QTLMAP_TRIO) , intent(inout)    :: this

     if ( associated(this%listNA)) then
      deallocate (this%listNA)
     end if

  end subroutine release_qtlmap_dataset_trio


    !!****f* m_qtlmap_cd/add_animal_genea
!! NAME
!!    add_animal_genea
!! DESCRIPTION
!!   Add a animal-triplet (sire-dam-progeny) inside a list of type DATASET_QTLMAP_TRIO
!!
!! INPUTS
!!   sire   : name of the sire
!!   dam    : name of the dam
!! progeny  : name of the progeny
!!  nd      : the index corresponding to the animal array
!!
!! INPUTS/OUTPUTS
!!  dataset : anima-triplet list
!!
!! SOURCE
     subroutine add_animal_genea(dataset,sire,dam,progeny,nd)
          class(DATASET_QTLMAP_TRIO)  , intent(inout)  :: dataset
          character(LEN=LEN_DEF) , intent(in)     :: sire,dam,progeny
          integer  ,optional     , intent(in)     :: nd             ! correspond a l index du tableau animal =>progeny

          integer  ,  parameter :: BLOCK_ALLOC = 1000
          logical :: findSire,findDam,findProg
          integer :: i,idSire,idDam,a

          if ( .not. associated(dataset%listNA) ) then
             call extend(dataset%listNA,BLOCK_ALLOC)
             call extend(dataset%trioKd,BLOCK_ALLOC)
             call extend(dataset%corranim,BLOCK_ALLOC)
          end if

          ! Le pere et la mere existe ?
          findSire=.false.
          findDam=.false.

         ! print *,sire,dam,progeny

          do i=1,dataset%na
            if ( trim(dataset%listNA(i)) == trim(sire)) then
               findSire=.true.
               idSire = i
            else if ( trim(dataset%listNA(i)) == trim(dam)) then
               findDam=.true.
               idDam = i
            end if
            if ( findSire .and. findDam ) exit
          end do

          if ( .not. findSire ) then
              !il faut creer une entree pour le pere
              if ( size(dataset%listNA) <= dataset%na ) then
                  call extend(dataset%listNA,BLOCK_ALLOC)
                  call extend(dataset%trioKd,BLOCK_ALLOC)
              end if
              dataset%na = dataset%na + 1
              dataset%listNA(dataset%na) = trim(sire)
              idSire = dataset%na
              !dataset%trioKd(dataset%na) => null()
          end if

          if ( .not. findDam ) then
              !il faut creer une entree pour la mere
              if ( size(dataset%listNA) <= dataset%na ) then
                  call extend(dataset%listNA,BLOCK_ALLOC)
                  call extend(dataset%trioKd,BLOCK_ALLOC)
              end if
              dataset%na = dataset%na + 1
              dataset%listNA(dataset%na) = trim(dam)
              idDam = dataset%na
              !dataset%trioKd(dataset%na) => null()
          end if

          !check pour savoir si la progeniture n existe pas deja !
          findProg = .false.

!          do i=1,dataset%na
!             if ( trim(dataset%listNA(i)) == trim(progeny)) then
!                print *,"double entry for kd:",progeny
!                stop
!             end if
!          end do

          if ( size(dataset%listNA) <= dataset%na ) then
            call extend(dataset%listNA,BLOCK_ALLOC)
            call extend(dataset%trioKd,BLOCK_ALLOC)
          end if

          dataset%na = dataset%na + 1

          dataset%listNA(dataset%na) = trim(progeny)

          dataset%trioKd(dataset%na)%idSire = idSire
          dataset%trioKd(dataset%na)%idDam  = idDam
          dataset%trioKd(dataset%na)%idKd   = dataset%na

          if (present(nd)) then
             if ( nd > size(dataset%corrAnim) ) then
               call extend(dataset%corranim,BLOCK_ALLOC)
             end if
             dataset%corranim(nd) = dataset%trioKd(dataset%na)%idKd
          end if

     end subroutine add_animal_genea
!!***

!!****f* m_qtlmap_cd/init_perf_animal
!! NAME
!!    init_perf_animal
!! DESCRIPTION
!!   Add performance and cd information about a progeniture contents in the animal-triplet list
!!
!! INPUTS
!! progeny   : name of the progeny
!! ncar      : number of trait (dimension of listPerf,listCdt)
!! listPerf  : phenotypic values
!! listCdt   : cd values
!!
!! INPUTS/OUTPUTS
!!  dataset : anima-triplet list
!!
!! SOURCE
     subroutine init_perf_animal(this,dpm,progeny,listPerf,listCdt)
       class(DATASET_QTLMAP_TRIO)  , intent(inout)   :: this
       type(DATAMODEL_BASE)       ,intent(inout)     :: dpm
       character(len=LEN_DEF)        , intent(in)    :: progeny
       real(kind=dp)  ,dimension(dpm%ncar) , intent(in)  :: listPerf,listCdt

       integer :: i,j
       call log_mess("init_perf_animal "//trim(progeny),DEBUG_DEF)

       do i=1,this%na
         if (trim(this%listNA(i))  == trim(progeny)) then
             this%trioKd(i)%presentc=.false.
             this%trioKd(i)%perfKd(:dpm%ncar) = listPerf(:dpm%ncar)

             do j=1,dpm%ncar
                if ( listCdt(j) /= 0 ) this%trioKd(i)%presentc(j)=.true.
             end do
             return
         end if
       end do
     end subroutine init_perf_animal
!!***

!!****f* m_qtlmap_cd/get_listProgenies
!! NAME
!!    get_listProgenies
!! DESCRIPTION
!!   get the list of progeniture which the sire or the dam is the animal indexed by ijkd
!!
!! INPUTS
!! dataset     : anima-triplet list
!! ijkd        : index of the parental progeny
!! nameAnimal  : name of the parental progeny
!!
!! OUTPUTS
!! sizeList    : size of the result list
!! list        : list of progenies
!! SOURCE
     subroutine get_listProgenies(dataset,ijkd,nameAnimal,sizeList,list)
        class (DATASET_QTLMAP_TRIO)         ,intent(in)     :: dataset
        integer                       , intent(in)    :: ijkd
        character(len=LEN_DEF)        , intent(in)    :: nameAnimal
        integer                      , intent(out)    :: sizeList
        type(QTLMAP_TRIO) ,dimension(:) , intent(out) , pointer  :: list

        integer  ,  parameter :: BLOCK_ALLOC = 100
        integer :: i,id

        list =>null()

        id = dataset%corranim(ijkd)

        if ( id <=0 .or. id > size(dataset%corranim)) then
            print *,'(get_listProgenies) can not get corresponding index of animal :',ijkd,' name:',nameAnimal
            stop
        end if

        call extendListTrio(list,BLOCK_ALLOC)
        sizeList = 0
        do i=1,dataset%na
            if ( dataset%trioKd(i)%idSire == id .or.  dataset%trioKd(i)%idDam == id ) then
              if ( sizeList >= size(list) ) then
                call extendListTrio(list,BLOCK_ALLOC)
              end if
              sizeList = sizeList + 1
              list(sizeList)=dataset%trioKd(i)
            end if
        end do

     end subroutine get_listProgenies
!!***

!!****f* m_qtlmap_cd/calcul_y_cd
!! NAME
!!    calcul_y_cd
!! DESCRIPTION
!!   Compute the new value phenotypic and a cd value of progenies who have a progeniture list. The old phenotypic value is removed if exist.
!!   depends the heritability of trait h2
!!
!!   NCK (ic,kd)      : number of progenies from sire/dam kd which theses progenies have a value for the trait ic
!!   NCiCjK(kd,ic,jc) : number of progenies from sire/dam kd which theses progenies have a value for the trait ic and jc jointly
!!
!!   SUM Yprog(kd)    : sum of trait values of progenies of sire/dam kd
!!   h2(ic)           : heritability of the trait ic
!!
!!   Y(ic,kd)  = SYM Yprog(kd) / NCK (ic,kd)
!!   CD(ic,kd) = 1.d0 + [h2(ic) * (NCK(kd,ic) - 1.d0) / 4] /  NCK (ic,kd)
!!
!! INPUTS/OUTPUTS
!! dataset     : anima-triplet list
!!
!! SOURCE
     subroutine calcul_y_cd(this,dg,dpa,dpm)
       class(DATASET_QTLMAP_TRIO)  ,intent(inout)       :: this

       type(GENEALOGY_BASE) ,intent(in)                 :: dg
       type(PHENOTYPE_BASE) ,intent(inout)              :: dpa
       type(DATAMODEL_BASE) ,intent(in)                 :: dpm

       type(QTLMAP_TRIO) ,dimension(:) ,pointer  :: list
       integer :: ip,jm,kd,s,i,j,ic,k,jc
       real(kind=dp) :: sumy(dpm%ncar)

       call log_mess("Calcul CD and Y with generation 3",INFO_DEF)

       allocate (this%NCK(dg%nd,dpm%ncar))
       this%NCK=0
       allocate (this%NCiCjK(dg%nd,dpm%ncar,dpm%ncar))
       this%NCiCjK=0

       do ip=1,dg%np
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
             call this%get_listProgenies(kd,dg%animal(kd),s,list)
             call log_mess(trim(dg%animal(kd))//":Number of progeny :"//str(s),DEBUG_DEF)
             sumy=0.d0
             do i=1,s
              do j=1,dpm%ncar
               print *,list(i)%presentc(j)
               if ( list(i)%presentc(j) ) then
                 this%NCK(kd,j) = this%NCK(kd,j) + 1
                 sumy(j) = sumy(j) + list(i)%perfkd(j)
                 do k=j+1,dpm%ncar
                    if ( list(i)%presentc(k) ) then
                     this%NCiCjK(kd,j,k) = this%NCiCjK(kd,j,k) + 1
                     this%NCiCjK(kd,k,j) = this%NCiCjK(kd,j,k)
                    end if
                 end do
               end if
              end do
             end do

             do ic=1,dpm%ncar
               dpa%presentc(ic,kd) = .false. !on invalide la performave de l individu kd
               if ( this%NCK(kd,ic) > 0 ) then
                print *,kd
                dpa%presentc(ic,kd) = .true.
                dpa%y(ic,kd)  = sumy(ic) / real(this%NCK(kd,ic))
                dpa%cd(ic,kd) = 1.d0 + 0.25d0*dpm%h2(ic) * (real(this%NCK(kd,ic)) - 1.d0)
                dpa%cd(ic,kd) = dpa%cd(ic,kd) / real(this%NCK(kd,ic))
                !cd = 0.25d0*h2(ic)
       !         print *,ic,kd,cd(ic,kd)
               end if
             end do


             if (s>0) deallocate(list)
           end do
         end do
       end do

       !cd=0.0d0

     end subroutine calcul_y_cd
!!***

!!****f* m_qtlmap_cd/calcul_corcd
!! NAME
!!    calcul_corcd
!! DESCRIPTION
!!   Compute the cd correlation (used by the multi trait analysis)
!!
!!   ic,jc 2 traits
!!   kd the sire/dam progenies
!!   NCK (ic,kd)      : number of progenies from sire/dam kd which theses progenies have a value for the trait ic
!!   NCiCjK(kd,ic,jc) : number of progenies from sire/dam kd which theses progenies have a value for the trait ic and jc jointly
!!   RhoG(ic,jc)      : genetic correlation between ic and jc
!!   RhoP(ic,jc)      : phenotypic correlation between ic and jc
!!
!!
!!                          RhoG(ic,jc) * SQRT( h2(ic)*h2(jc) )
!!   v(ic,jc)          =    -----------------------------------
!!                                     RhoP(ic,jc)
!!
!!                       NCiCjK(kd,ic,jc) + 0.25d0*v(ic,jc)* [ NCK(kd,ic)*NCK(kd,jc) - NCiCjK(kd,ic,jc) ]
!!   CORR_CD(ic,jc,kd) = ---------------------------------------------------------------------------------
!!                                                     NCK(kd,ic) * NCK(kd,jc)
!!
!! INPUTS/OUTPUTS
!! dataset     : anima-triplet list
!!
!! SOURCE
     subroutine calcul_corcd(this,dg,dpa,dpm)
       class(DATASET_QTLMAP_TRIO)         ,intent(inout)       :: this

       type(GENEALOGY_BASE) ,intent(in)                :: dg
       type(PHENOTYPE_BASE) ,intent(in)                :: dpa
       type(DATAMODEL_BASE) ,intent(in)                :: dpm


       integer :: ip,jm,kd,s,i,j,ic,k,jc

       real(kind=dp) :: r,v(dpm%ncar,dpm%ncar)

       allocate (this%corcd(dpm%ncar,dpm%ncar,dg%nd))

       do ic=1,dpm%ncar
        do jc=ic+1,dpm%ncar
             v(ic,jc) = (dpm%RhoG(ic,jc)/dpm%RhoP(ic,jc))*sqrt(dpm%h2(ic)*dpm%h2(jc))
        end do
       end do

       this%corcd = 1.d0
       do ip=1,dg%np
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
             do ic=1,dpm%ncar
               if ( this%NCK(kd,ic) > 0 ) then
                 do jc=ic+1,dpm%ncar
                   if ( this%NCK(kd,jc) > 0 ) then
                     r = real( this%NCK(kd,ic)*this%NCK(kd,jc) - this%NCiCjK(kd,ic,jc))
                     this%corcd(ic,jc,kd) = real(this%NCiCjK(kd,ic,jc)) + 0.25d0*v(ic,jc)*r
                     this%corcd(ic,jc,kd) = this%corcd(ic,jc,kd) / &
                      ( real( this%NCK(kd,ic)) * real( this%NCK(kd,jc)))
                     this%corcd(jc,ic,kd) = this%corcd(ic,jc,kd)
                   end if
                 end do
               end if
             end do
           end do
         end do
       end do

     end subroutine calcul_corcd
!!***

!**********************************************************************************************************************************
!     subroutine copy_trio_link(in,out)
!       type(QTLMAP_TRIO) , intent(in) :: in
!       type(QTLMAP_TRIO) , intent(out) :: out
!       integer :: i
!
!       out%idSire = in%idSire
!       out%idDam = in%idDam
!       out%idKd = in%idKd
!       allocate (out%perfKd(size(in%perfKd)))
!       allocate (out%presentc(size(in%presentc)))
!       out%perfKd = in%perfKd
!       out%presentc = in%presentc
!
!    end subroutine copy_trio_link

      !    subroutine release_trio_link(inout)
!      type(QTLMAP_TRIO) , intent(in) :: inout
!
!      deallocate (inout%perfKd)
!      deallocate (inout%presentc)
!
!    end subroutine release_trio_link



!!****f* m_qtlmap_cd/extendListNa
!! NAME
!!    extendListNa
!! DESCRIPTION
!!   increase the list of type character(LEN=LEN_DEF) with BLOCK_ALLOC values
!!
!! INPUTS
!!  BLOCK_ALLOC : the size to extends
!!
!! INPUTS/OUTPUTS
!!  list        : the character list
!!
!! SOURCE
      subroutine extendListNa(list,BLOCK_ALLOC)
        character(LEN=LEN_DEF) ,dimension(:) ,pointer  ,intent(inout) :: list
        integer                                        ,intent(in)    :: BLOCK_ALLOC

        character(LEN=LEN_DEF) ,dimension(:) ,pointer  :: listTemp => null()

        !start
        if ( .not. associated(list) ) then
           allocate (list(BLOCK_ALLOC))
           return
        end if

        allocate ( listTemp(size(list)) )
        listTemp = list
        deallocate ( list )
        allocate ( list(size(listTemp)+BLOCK_ALLOC) )
        list(:size(listTemp)) = listTemp
        deallocate (listTemp)

        return
     end subroutine extendListNa
!!***

!!****f* m_qtlmap_cd/extendListTrio
!! NAME
!!    extendListTrio
!! DESCRIPTION
!!   increase the list of type type(QTLMAP_TRIO) with BLOCK_ALLOC values
!!
!! INPUTS
!!  BLOCK_ALLOC : the size to extends
!!
!! INPUTS/OUTPUTS
!!  list        : the type(QTLMAP_TRIO) list
!!
!! SOURCE
      subroutine extendListTrio(list,BLOCK_ALLOC)
        type(QTLMAP_TRIO) ,dimension(:) ,pointer  ,intent(inout) :: list
        integer                                        ,intent(in)    :: BLOCK_ALLOC

        type(QTLMAP_TRIO) ,dimension(:) ,pointer  :: listTemp => null()
        integer :: i

        !start
        if ( .not. associated(list) ) then
           allocate (list(BLOCK_ALLOC))
           return
        end if

        allocate ( listTemp(size(list)) )
        listTemp = list
        deallocate ( list )
        allocate ( list(size(listTemp)+BLOCK_ALLOC) )
        list(:size(listTemp)) = listTemp
        deallocate (listTemp)

        return
     end subroutine extendListTrio
!!***

!!****f* m_qtlmap_cd/extendListCorrAnim
!! NAME
!!    extendListCorrAnim
!! DESCRIPTION
!!   increase the list of type integer with BLOCK_ALLOC values
!!
!! INPUTS
!!  BLOCK_ALLOC : the size to extends
!!
!! INPUTS/OUTPUTS
!!  list        : the integer list
!!
!! SOURCE
     subroutine extendListCorrAnim(list,BLOCK_ALLOC)
        integer                ,dimension(:) ,pointer  ,intent(inout) :: list
        integer                                        ,intent(in)    :: BLOCK_ALLOC

        integer ,dimension(:) ,pointer  :: listTemp => null()
        integer :: i

        !start
        if ( .not. associated(list) ) then
           allocate (list(BLOCK_ALLOC))
           list=0
           return
        end if

        allocate ( listTemp(size(list)) )
        listTemp = list
        deallocate ( list )
        allocate ( list(size(listTemp)+BLOCK_ALLOC) )
        list=0
        list(:size(listTemp)) = listTemp
        deallocate (listTemp)

        return
     end subroutine extendListCorrAnim
!!***


end module m_qtlmap_type_dataset_trio
