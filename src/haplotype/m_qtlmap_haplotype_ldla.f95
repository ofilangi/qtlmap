!!****m* HAPLOTYPE/m_qtlmap_haplotype_ldla
!!  NAME
!!    m_qtlmap_haplotype_ldla
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
module m_qtlmap_haplotype_ldla
    use m_qtlmap_base
    use m_qtlmap_types
    use m_qtlmap_log

    implicit none

!!****v* m_qtlmap_haplotype_ldla/haplo
!! NAME
!!   haplo
!! DESCRIPTION
!!   list of the possible haplotype (with unknwon allele) at the current position
!! DIMENSIONS
!!   1: index of the haplotype
!!   2: index of the allele (1 < index allele < longap)
!! NOTES
!! SOURCE
!    integer(kind=KIND_PHENO), dimension (:,:)    , allocatable   :: haplo
!!***

!!****v* m_qtlmap_haplotype_ldla/race_h
!! NAME
!!   race_h
!! DESCRIPTION
!!   breed origin of haplotype haplo
!! DIMENSIONS
!!   1: index of the haplotype
!! NOTES
!! SOURCE
!     character(len=LEN_DEF), dimension (:)    , allocatable      :: race_h
!!***

!!****v* m_qtlmap_haplotype_ldla/haplo_complet
!! NAME
!!   haplo_complet
!! DESCRIPTION
!!   list of the possible haplotype (without unknwon allele) at the current position
!! DIMENSIONS
!!   1: index of the complete haplotype
!!   2: index of the allele (1 < index allele < dataset%params%longhap)
!! NOTES
!!  see liste_haplo_complet
!! SOURCE
!    integer(kind=KIND_PHENO), dimension (:,:)    , allocatable   :: haplo_complet
!!***

!!****v* m_qtlmap_haplotype_ldla/race_haplo_complet
!! NAME
!!   race_haplo_complet
!! DESCRIPTION
!!   breed origin of haplotype haplo_complet
!! DIMENSIONS
!!   1: index of the complete haplotype
!! NOTES
!! SOURCE
!     character(len=LEN_DEF), dimension (:)    , allocatable      :: race_haplo_complet
!!***

!!****v* m_qtlmap_haplotype_ldla/gamete
!! NAME
!!   gamete
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1: 1,2
!!   2: marker (1 <= <= dataset%params%longhap)
!!   3: 2**dataset%params%longhap
!! NOTES
!!  gamete(1,lk,j_gam) l allele paternel au lk eme marqueur pour le j_gam eme gamete
!!	gamete(2,lk,j_gam) l allele maternel au lk eme marqueur pour le j_gam eme gamete
!! SOURCE
!    integer(kind=KIND_PHENO), dimension (:,:,:) ,allocatable :: gamete
!!***

!!****v* m_qtlmap_haplotype_ldla/loc_haplo
!! NAME
!!   loc_haplo
!! DESCRIPTION
!!
!!                !-----!
!!  A A A T T T A A A A T T T T A T A T A T A T T T T
!!                | | | |
!!                | | | ---------------
!!                | | -------------   |
!!                | -------       |   |
!!         loc_haplo(1)  |       |   |
!!                 loc_haplo(2)  |   |
!!                      loc_haplo(3) |
!!                           loc_haplo(4)
!!
!! DIMENSIONS
!!   1: index of the marker (size:dataset%params%longhap)
!! NOTES
!!  loc_haplo, le vecteur de dimension dataset%params%longhap qui contient les numeros des marqueurs
!!  pour l haplotype entourant la position testee dx, et lk le premier marqueur flanquant a droite
!! SOURCE
!    integer                 , dimension (:)      , allocatable   :: loc_haplo
!!***

!!****v* m_qtlmap_haplotype_ldla/count_haplo_complet
!! NAME
!!   count_haplo_complet
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1:
!! NOTES
!! SOURCE
!    real(kind=dp)        , dimension (:)      , allocatable   :: count_haplo_complet
!!***

!!****v* m_qtlmap_haplotype_ldla/pb_haplo_complet
!! NAME
!!   pb_haplo_complet
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1:
!!   2:
!!   3:
!! NOTES
!! SOURCE
!    real(kind=dp)        , dimension (:)      , allocatable   :: pb_haplo_complet
!!***

!!****v* m_qtlmap_haplotype_ldla/tab_haplo_complet
!! NAME
!!   tab_haplo_complet
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1:
!!   2:
!! NOTES
!! SOURCE
!    integer(kind=KIND_PHENO), dimension (:,:)    , allocatable   :: tab_haplo_complet
!!***

!!****v* m_qtlmap_haplotype_ldla/nb_haplo_possible
!! NAME
!!   nb_haplo_possible
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1:
!! NOTES
!! SOURCE
!    integer              , dimension (:)    , allocatable     :: nb_haplo_possible
!!***
!!****v* m_qtlmap_haplotype_ldla/comp_haplo
!! NAME
!!   comp_haplo
!! DESCRIPTION
!!
!! DIMENSIONS
!!   1:
!!   2:
!! NOTES
!! SOURCE
!    logical                 , dimension (:,:)    , allocatable   :: comp_haplo
!!***

    !!!!$omp threadprivate (haplo,haplo_complet,race_h,race_haplo_complet,gamete,loc_haplo,count_haplo_complet)
    !!!!$omp threadprivate (pb_haplo_complet,tab_haplo_complet,comp_haplo,nb_haplo_possible)

!!****t* QTLMAP_TYPES/HAPLOTYPE_POSITION_BUILD
!!  NAME
!!     HAPLOTYPE_POSITION_BUILD
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type ,public :: HAPLOTYPE_POSITION_BUILD

      ! Internal
      type(QTLMAP_DATASET)             ,private , pointer :: dataset

      type(PDD_BUILD)                  ,private , pointer :: spt

      integer(kind=KIND_PHENO)    ,dimension(:,:,:),pointer           ,public :: tab_IBS        => null()
      integer        ,dimension(:),pointer                            ,public :: shared_haplo   => null()
      integer        ,dimension(:,:,:),pointer                        ,public :: tab_shared_haplo => null()
      integer                   , dimension (:,:,:)    , pointer    ,  public :: num_haplo_pere => null()
      integer                   , dimension (:,:,:)    , pointer    ,  public :: num_haplo_mere => null()
      integer                   , dimension (:,:)      , pointer    ,  public :: num_haplo_desc => null()

      integer                   , dimension (:,:)      , pointer    ,  public :: nb_gam_pere => null()
      integer                   , dimension (:,:)      , pointer    ,  public :: nb_gam_mere => null()
      integer                   , dimension (:)        , pointer    ,  public :: nb_gam_desc => null()
      integer                   , dimension (:)        , pointer    , private :: nb_gam      => null()
      real   (kind=dp)          , dimension (:,:)      , pointer    , private :: prob_gam    => null()
      double precision          , dimension (:,:)      , pointer    ,  public :: pb_haplo_desc => null()
      integer                                                        , public :: nb_max_haplo = 0

      integer ,                                              public      :: nb_haplo           = 0
      integer ,                                              public      :: nb_haplo_complet   = 0
      integer ,                                              public      :: nb_haplo_reduit    = 0
      integer                   , dimension (:) , pointer , private      :: nb_haplo_possible => null()

      integer                 , dimension (:)      , pointer             :: liste_haplo_reduit => null()
      real(kind=dp)                 , dimension (:)          , pointer   :: pb_haplo_reduit    => null()
      real(kind=dp)                 , dimension (:)          , pointer   :: pb_haplo_reduit_r  => null()
      logical                 , dimension (:,:)              , pointer   :: comp_haplo_reduit  => null()

      character(LEN=400)                  , dimension (:)    , pointer    ,  public :: name_haplo         => null()
      character(LEN=400)                  , dimension (:)    , pointer    ,  public :: name_haplo_complet => null()
      character(LEN=400)                  , dimension (:)    , pointer    ,  public :: name_haplo_reduit  => null()
      character(LEN=400)                  , dimension (:)    , pointer    ,  public :: race_haplo_reduit  => null()
      real(kind=dp)                  , dimension (:)         , pointer    ,  public :: count_haplo => null()



      !provient du module haplotype_ldla
      integer(kind=KIND_PHENO), dimension (:,:)    , pointer   :: haplo => null()
      character(len=LEN_DEF), dimension (:)    , pointer       :: race_h => null()
  !    integer(kind=KIND_PHENO), dimension (:,:)    , pointer   :: haplo_complet => null()
      character(len=LEN_DEF), dimension (:)    , pointer       :: race_haplo_complet => null()
      integer(kind=KIND_PHENO), dimension (:,:,:) ,pointer     :: gamete  => null()
      integer                 , dimension (:)      , pointer   :: loc_haplo => null()
      real(kind=dp)        , dimension (:)      , pointer      :: count_haplo_complet => null()
      real(kind=dp)        , dimension (:)      , pointer , private   :: pb_haplo_complet => null()
      integer(kind=KIND_PHENO), dimension (:,:)    , pointer   :: tab_haplo_complet => null()

      logical                 , dimension (:,:)    , pointer   :: comp_haplo => null()

      contains
        !initialisation of structure
        procedure ,public :: set
        ! Computation of a list of haplotype for each desc
        procedure ,public :: set_haplo_for_ldla
        !destructor
        procedure ,public :: free

        procedure ,public :: sort_haplo
        procedure ,public :: liste_haplo

        procedure ,public :: liste_shared_haplo

        ! Internal procedure
        procedure , private :: check_change_marker_windows

        procedure ,private  :: liste_haplo_complet
        procedure ,private  :: proba_haplo_complet
        procedure ,private  :: tri_haplo_complet
        procedure ,private  :: set_haplo_final

        procedure ,private :: local_gamete
        procedure ,private :: point_gamete
        procedure ,private :: liste_gamete
        procedure ,private :: set_tab_ibs

        procedure ,private :: global_gamete

   end type HAPLOTYPE_POSITION_BUILD
!!***

   contains

!!****f* m_qtlmap_haplotype_ldla/set_allocation
!! NAME
!!    set_allocation
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
      subroutine set(shp,dataset,spt)
        class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)        :: shp
        type(QTLMAP_DATASET) ,target      ,intent(in)           :: dataset
        type(PDD_BUILD)     ,target       ,intent(in)           :: spt


        integer  :: stat,value_alloc_1,value_alloc_2,c
        type(GENEALOGY_BASE) , pointer :: dg
        type(GENEALOGY_RACE) , pointer :: dgr

        shp%dataset => dataset
        shp%spt => spt

        dg => shp%dataset%genea
        dgr => shp%dataset%geneaRace

      allocate (shp%loc_haplo(dataset%params%longhap))
      allocate (shp%nb_gam(dg%nd))
      allocate (shp%nb_gam_pere(dg%np,2))
      allocate (shp%nb_gam_mere(maxval(spt%ngenom(:,dg%nm+1)),2))
      allocate (shp%tab_shared_haplo(dg%np,2,1+maxval(dataset%map%nmk)))

      value_alloc_1= (2*(dg%np+dg%nm)+dg%nd*2**shp%dataset%params%longhap)*dgr%nb_races
      ! OFI Janv 2013 => si plus de deux alleles ca peut planter car pas assez de memoire...a corriger
      value_alloc_2=2**(2*dataset%params%longhap)+1
      allocate (shp%haplo(value_alloc_1,dataset%params%longhap))
      allocate (shp%nb_haplo_possible(value_alloc_1))
      !add race info
      allocate (shp%race_h(value_alloc_1))
      allocate (shp%race_haplo_complet(value_alloc_1))
      allocate (shp%count_haplo_complet(value_alloc_1))
      allocate (shp%prob_gam(dg%nd,value_alloc_2))
      allocate (shp%nb_gam_desc(dg%nd))
      allocate (shp%gamete(2,dataset%params%longhap,value_alloc_2))
      allocate (shp%tab_IBS(dg%nd,maxval(dataset%map%nmk),2))
      allocate (shp%name_haplo(value_alloc_2))
      allocate (shp%name_haplo_reduit(value_alloc_2))
      allocate (shp%race_haplo_reduit(value_alloc_2))
      allocate (shp%name_haplo_complet(value_alloc_2))
      allocate (shp%count_haplo(value_alloc_2))
      allocate (shp%pb_haplo_complet(10*dgr%nb_races))
      allocate (shp%comp_haplo(10,10*dgr%nb_races))
      allocate (shp%liste_haplo_reduit(10*dgr%nb_races))
      allocate (shp%tab_haplo_complet(10*dgr%nb_races ,dataset%params%longhap))
      allocate (shp%comp_haplo_reduit(10,50))
      shp%comp_haplo_reduit=.false.
      allocate (shp%pb_haplo_reduit(50))
      allocate (shp%pb_haplo_reduit_r(50))
      allocate (shp%num_haplo_pere(dg%np,2,10))
      allocate (shp%num_haplo_mere(maxval(spt%ngenom(:,dg%nm+1)),2,10))
      allocate (shp%num_haplo_desc(dg%nd,100))
      allocate (shp%pb_haplo_desc(dg%nd,100))

    end subroutine set
!!***

!!****f* m_qtlmap_haplotype_ldla/free
!! NAME
!!    free
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine free(shp)
        class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

        deallocate (shp%haplo)
    !    deallocate (shp%haplo_complet)
        deallocate (shp%race_h)
        deallocate (shp%loc_haplo)
        deallocate (shp%race_haplo_complet)
        deallocate (shp%gamete)
        deallocate (shp%count_haplo_complet)
        deallocate (shp%nb_gam)
        deallocate (shp%prob_gam)
        deallocate (shp%num_haplo_pere)
        deallocate (shp%num_haplo_mere)
        deallocate (shp%num_haplo_desc)
        deallocate (shp%nb_gam_pere)
        deallocate (shp%nb_gam_mere)
        deallocate (shp%nb_gam_desc)
        deallocate (shp%tab_IBS)
        deallocate (shp%name_haplo)
        deallocate (shp%name_haplo_reduit)
        deallocate (shp%race_haplo_reduit)
        deallocate (shp%name_haplo_complet)
        deallocate (shp%count_haplo)
        deallocate (shp%tab_shared_haplo)
        deallocate (shp%comp_haplo)
        deallocate (shp%comp_haplo_reduit)
        deallocate (shp%liste_haplo_reduit)
        deallocate (shp%pb_haplo_reduit)
        deallocate (shp%pb_haplo_reduit_r)
        deallocate (shp%pb_haplo_complet)
        deallocate (shp%tab_haplo_complet)
        deallocate (shp%nb_haplo_possible)
        deallocate (shp%pb_haplo_desc)


    end subroutine free
!!***

!!****f* m_qtlmap_haplotype_ldla/set_haplo_for_ldla
!! NAME
!!    set_haplo_for_ldla
!! DESCRIPTION
!!
!! NOTES
!!  set_haplo construit la liste des haplotypes de longueur dataset%params%longhap
!!
!! Utilise nchr,posi,pas,nmk de m_qtlmap_data
!! utilise dataset%params%longhap de m_qtlmap_data
!! SOURCE
      subroutine set_haplo_for_ldla(shp,c,dx,n,hsire,hdam)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
      integer c,i,j,k
      double precision dx
      integer, intent (in) :: n
      logical,intent(in)   :: hsire,hdam

      !print *,"ldla:",c,dx,shp%check_change_marker_windows(c,dx)
      !bug OFI Janvier 2013 => check_change_marker_windows ne fonctionne pas...
      !if (.not. shp%check_change_marker_windows(c,dx)) return
!
!  dataset%params%longhap est la longueur de l haplotype dont on estime l effet
!  dataset%params%longhap est defini en nombre de marqueurs ( multiple de 2)
 !       PROB_HAPLO_MIN=0.0d0
        !call free
        !call set
        call shp%liste_haplo(c,dx,n,hsire,hdam)
        call shp%liste_haplo_complet(c)
        call shp%proba_haplo_complet(c,dx)
        call shp%tri_haplo_complet(c,dx)
        call shp%set_haplo_final(c,n,hsire,hdam)
      end subroutine set_haplo_for_ldla
!!***

!!****f* m_qtlmap_haplotype_ldla/set_tab_ibs
!! NAME
!!    set_tab_ibs
!! DESCRIPTION
!!
!! NOTES
!! set_tab_ibs construit la liste des haplotypes correspondant au chromosome entier
!! SOURCE
      subroutine set_tab_ibs(shp,c)
        class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
!
! set_tab_ibs construit la liste des haplotypes correspondant au chromosome entier
!
! Utilise nchr,posi,pas,nmk de m_qtlmap_data
! utilise dataset%params%longhap de m_qtlmap_data
!
      integer, intent (in) :: c
!
!  dataset%params%longhap est la longueur de l haplotype dont on estime l effet
!  dataset%params%longhap est defini en nombre de marqueurs ( multiple de 2)
        !call free
        !call set
        call shp%liste_gamete(c)
       return
 
     end subroutine set_tab_ibs
!!***

!!****f* m_qtlmap_haplotype_ldla/sort_haplo
!! NAME
!!    sort_haplo
!! DESCRIPTION
!!
!! NOTES
!!  sort_haplo trie les haplotypes dans la liste haplo et les recodifie
!!  chaque reproducteur est affecte de ses haplotype recodifies
!!  la liste des noms d'haplotype (name_haplo) est cree
!!
!! Utilise np,nm,ndm,estfem,repfem de m_qtlmap_data
!! utilise ngenom de m_qtlmap_haplotpe_data
!! SOURCE
     subroutine sort_haplo(shp,c,n,hsire,hdam)
       class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

        logical :: hsire,hdam
! dim allocatable
        integer(kind=KIND_PHENO)  , dimension (:,:)    , allocatable    :: kept_haplo
        integer                   , dimension (:)      , allocatable    :: corr_haplo
        integer                   , dimension (:)      , allocatable    :: corr_corr_haplo
        integer                   , dimension (:)      , allocatable    :: inv_corr_haplo
! divers
      integer c
      integer i,i_haplo,i_kept,ip,j,jm,j_gam,kd,geno,previous,j_haplo
      integer nb_kept,value_alloc
      integer  :: stat,lstr1,lstr2
      integer, intent (in) :: n
      character(len=LEN_DEF) :: buf
      double precision dencount
      logical inconnu, new
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga

      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea

! allocation
      value_alloc=shp%nb_haplo
      allocate (kept_haplo(value_alloc,shp%dataset%params%longhap))
      allocate (corr_haplo(value_alloc))
      corr_haplo=0

      allocate (corr_corr_haplo(value_alloc))
      allocate (inv_corr_haplo(value_alloc))

! on trie les haplotypes trouves en eliminant les redondances
!
! la liste reduite sera dans kept_haplo
! de longueur nb_kept
! le i eme ancien haplotype est a la corr_haplo(i) eme place dans kept_haplo
!
      shp%name_haplo=''
      kept_haplo(1,1)=shp%haplo(1,1)
      lstr1=len(trim(get_pheno(dga,c,shp%haplo(1,1))))
      shp%name_haplo(1)(1:lstr1)=trim(get_pheno(dga,c,shp%haplo(1,1)))

      do i=2,shp%dataset%params%longhap
         kept_haplo(1,i)=shp%haplo(1,i)
         lstr2=len(trim(get_pheno(dga,c,shp%haplo(1,i))))
         shp%name_haplo(1)(lstr1+1:lstr1+lstr2)=trim(get_pheno(dga,c,shp%haplo(1,i)))
         lstr1=lstr1+lstr2
      end do !i
      corr_haplo(1)=1
      nb_kept=1

      do i_haplo=2, shp%nb_haplo
        inconnu=.true.
        do i_kept =1,nb_kept
          new=.false.         
          do i=1,shp%dataset%params%longhap
            if(kept_haplo(i_kept,i) /= shp%haplo(i_haplo,i))new=.true.
            if(new)exit 
          end do
          if(.not.new) previous=i_kept
          if(.not.(inconnu.and.new))inconnu=.false.
        end do ! i_kept
        if(inconnu)then
          nb_kept=nb_kept+1
          kept_haplo(nb_kept,1)=shp%haplo(i_haplo,1)
          lstr1=len(trim(get_pheno(dga,c,shp%haplo(i_haplo,1))))
          shp%name_haplo(nb_kept)(1:lstr1)=trim(get_pheno(dga,c,shp%haplo(i_haplo,1)))
          do i=2,shp%dataset%params%longhap
            kept_haplo(nb_kept,i)=shp%haplo(i_haplo,i)
            lstr2=len(trim(get_pheno(dga,c,shp%haplo(i_haplo,i))))
            shp%name_haplo(nb_kept)(lstr1+1:lstr1+lstr2)=trim(get_pheno(dga,c,shp%haplo(i_haplo,i)))
            lstr1=lstr1+lstr2
            buf=get_pheno(dga,c,shp%haplo(i_haplo,i))
          end do !i
          corr_haplo(i_haplo)=nb_kept
        else
          corr_haplo(i_haplo)=previous
        end if
      end do !i_haplo
!
!  Compteur d'haplotypes
!
      dencount=1.d0/dble(shp%nb_haplo)
      shp%count_haplo(:)=0.d0
      do i_haplo=1,shp%nb_haplo
        shp%count_haplo(corr_haplo(i_haplo))=shp%count_haplo(corr_haplo(i_haplo))+dencount
      end do

      shp%nb_max_haplo=nb_kept

      if(shp%nb_max_haplo > shp%dataset%params%NB_HAPLO_PRIOR)then
       call stop_application('more haplotypes found than expected ['//trim(str(shp%nb_haplo))//&
                             '] opt_nb_haplo_prior must be increased')
      end if

!
!  regroupement des haplotypes rares dans le paquet des donnees manquantes
!
       nb_kept=1
       corr_corr_haplo(1)=1
       inv_corr_haplo(1)=1
       do i_haplo=2,shp%nb_max_haplo
        if(shp%count_haplo(i_haplo) < shp%dataset%params%PROB_HAPLO_MIN) then
          corr_corr_haplo(i_haplo)=1
          inv_corr_haplo(1)=1!peut etre que c est i_haplo....
        else
          nb_kept=nb_kept+1
          corr_corr_haplo(i_haplo)=nb_kept
          inv_corr_haplo(nb_kept)=i_haplo
        end if
      end do ! i_haplo
      shp%count_haplo(:)=0.d0
      do i_haplo=1,shp%nb_haplo
        j_haplo=corr_haplo(i_haplo)
        shp%count_haplo(corr_corr_haplo(j_haplo))=shp%count_haplo(corr_corr_haplo(j_haplo))+dencount
      end do
      shp%nb_max_haplo=nb_kept

!  creation d un indice d haplotype par reproducteur
!
!  pour les peres
       j_haplo=1
       if(hsire)then
       do ip=1,dg%np
         do j=1,2
           j_haplo=j_haplo+1
           shp%num_haplo_pere(ip,j,1)=corr_corr_haplo(corr_haplo(j_haplo))
         end do !j
       end do !ip
       end if ! opt_model

       if(hdam)then
!
!  pour les meres 
!
       do jm=1,dg%nm
         if(dga%estfem(dg%repfem(jm)))then
!
!  meres phasees
!
           do geno=shp%spt%ngenom(c,jm)+1,shp%spt%ngenom(c,jm+1)
             do j=1,2
               j_haplo=j_haplo+1
               shp%num_haplo_mere(geno,j,1)=corr_corr_haplo(corr_haplo(j_haplo))
             end do !j
           end do !geno
         else
! 
!   non phasees
           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
             if(dga%presentg(c,kd))then
             do j_gam=1,shp%nb_gam(kd)
               j_haplo=j_haplo+1
!               shp%num_haplo_desc(kd,j_gam,1)=corr_corr_haplo(corr_haplo(j_haplo))
               shp%num_haplo_desc(kd,j_gam)=corr_corr_haplo(corr_haplo(j_haplo))
             end do !j_gam
             end if
           end do !kd
         end if
       end do !jm
       end if ! opt_model

      deallocate(kept_haplo)
      deallocate(corr_haplo)
      deallocate(corr_corr_haplo)
      deallocate(inv_corr_haplo)

     end  subroutine sort_haplo
!!***

!!****f* m_qtlmap_haplotype_ldla/check_change_maker_windows
!! NAME
!!    liste_haplo
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
     function check_change_marker_windows(shp,c,dx) result(change)
       class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
       integer       , intent(in) :: c
       real(kind=dp) , intent(in) :: dx

       logical :: change
       integer :: lkmin,lkmax,lk,dx2
       type(MAP_BASE) , pointer    :: map

       map => shp%dataset%map

       change = .true.
       if (.not. associated(shp%loc_haplo) ) then
      !   print *,'not allocated...'
         return
       end if
      ! lkmax est la limite  pour le marqueur flanquant l haplotype sur sa droite
!
       lkmin=ceiling(shp%dataset%params%longhap/2.d0)
       if ( lkmin == 0 ) lkmin = 1
       lkmax=map%nmk(c)-lkmin+1
!  la position testee doit etre compatible avec la longueur de l'haplotype
!

!       if((dx < dataset%map%posi(c,lkmin)) .or. (dx > dataset%map%posi(c,lkmax))) then
!         call stop_application("check_change_marker_windows : "//&
!         "incompatibilty between haplotype length and tested position."//&
!         " position:"//trim(str(dx)))
!       end if

       dx2=dx

       if ( dx2 < map%posi(c,lkmin) ) dx2 = map%posi(c,lkmin)
       if ( dx2 > map%posi(c,lkmax) ) dx2 = map%posi(c,lkmax)

       lk=lkmin +1

       do while (dx2 > map%posi(c,lk) .and. lk < lkmax)
         lk=lk+1
       end do

       if ( shp%loc_haplo(1) == lk-lkmin ) change=.false.

     end function check_change_marker_windows
!!***


!!****f* m_qtlmap_haplotype_ldla/liste_haplo
!! NAME
!!    liste_haplo
!! DESCRIPTION
!!
!! NOTES
!!  liste_haplo cherche  les haplotypes de longeur dataset%params%longhap possibles
!!  autour de la position dx du chromosome c
!! Utilise nmk,np,posi,correp,nm,repfem,correm,estfem    de m_qtlmap_data
!! Utilise genotyp,ngenom,genotypm                       de m_qtlmap_haplotype_data
!! SOURCE
    subroutine liste_haplo(shp,c,dx_user,n,hsire,hdam)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
      integer       , intent(in) :: c
      real(kind=dp) , intent(in) :: dx_user
      integer, intent (in)       :: n
      logical       , intent(in) :: hsire,hdam

      integer lkmax,lkmin,i,ip,j,jm,geno,j_gam,kd,lk
      real(kind=dp) :: dx
      logical :: raceFileDefined

      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_RACE) , pointer    :: dgr
      type(MAP_BASE) , pointer    :: map

      map => shp%dataset%map
      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea
      dgr => shp%dataset%geneaRace

      raceFileDefined = shp%dataset%params%get_file_val(K_RACE) /= ""

      shp%nb_gam=0

      dx = dx_user

! lkmax est la limite  pour le marqueur flanquant l haplotype sur sa droite
!
       lkmin=ceiling(shp%dataset%params%longhap/2.d0)
       if ( lkmin == 0 ) lkmin = 1
       lkmax=map%nmk(c)-lkmin+1
!  la position testee doit etre compatible avec la longueur de l'haplotype
!
       if ((dx < map%posi(c,lkmin))) dx = map%posi(c,lkmin)
       if (dx > map%posi(c,lkmax)) dx = map%posi(c,lkmax)

!       if((dx < dataset%map%posi(c,lkmin)) .or. (dx > dataset%map%posi(c,lkmax))) then
!         call stop_application("liste_haplo : incompatibilty between haplotype length and tested position.")
!       end if

!
! initialisation de loc_haplo, le vecteur de dimension dataset%params%longhap qui contient les numeros des marqueurs
! pour l haplotype entourant la position testee dx, et lk le premier marqueur flanquant a droite
       do i=1,shp%dataset%params%longhap
         shp%loc_haplo(i)=i
       end do !i

       lk=shp%loc_haplo(lkmin) +1

       do while (dx > map%posi(c,lk) .and. lk < lkmax)
         lk=lk+1
       end do

!
! lk est le premier marqueur a droite de dx
       do i=1,shp%dataset%params%longhap
         shp%loc_haplo(i)=lk-lkmin-1+i
       end do  !i 

       !print *,shp%loc_haplo
!
!  le premier haplotype dans la liste sera celui des donnees manqantes, il y en aura un par race quand il y a plusieurs races
!
       shp%nb_haplo=0
       shp%race_h='UNKNOWN'

       do j=1,dgr%nb_races
          shp%nb_haplo=shp%nb_haplo+1
          do i=1,shp%dataset%params%longhap
           shp%haplo(shp%nb_haplo,i)=dga%nmanque
         end do !i
           shp%race_h(shp%nb_haplo)=dgr%nom_race(j)
       enddo
!
!  haplotypes trouves chez les peres
!
       if(hsire)then
       do ip=1,dg%np
         do j=1,2
           shp%nb_haplo=shp%nb_haplo+1
           do i=1,shp%dataset%params%longhap
             shp%haplo(shp%nb_haplo,i)=shp%spt%genotyp(c,shp%loc_haplo(i),dga%correp(ip),j)
           end do !i
!! Si les parents du père sont connus (reppere.ne.INT_NOT_DEFINED) et sont 
!! soit de meme race soit génotypés sur au moins un marqueur
           if (raceFileDefined) then
               if (dg%reppere(ip).eq.INT_NOT_DEFINED) then
                    call stop_application('Both parents of sires and dams have to be defined in the ped file '//&
           'to use breed origin of haplotypes. The parents of sire  '//trim(dg%pere(ip))//&
                    ' are not in the ped file! ')
               else          
                 if (( dgr%racep(dg%reppere(ip))==dgr%racem(dg%reppere(ip)) )  &
                .or. (shp%spt%phasp(c,ip)                             )) then
                    if (j==1) shp%race_h(shp%nb_haplo)= dgr%racep(dg%reppere(ip))
                    if (j==2) shp%race_h(shp%nb_haplo)= dgr%racem(dg%reppere(ip))
                 endif
              endif 
            endif
         end do !j

       end do !ip
       end if ! hsire

       if(hdam)then

!
!  haplotypes trouves chez les meres phasees
!
  
       do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
         if(dga%estfem(dg%repfem(jm)))then
!
!  cas des femelles dont on estime la phase
!
           do geno=shp%spt%ngenom(c,jm)+1,shp%spt%ngenom(c,jm+1)
             do j=1,2
               shp%nb_haplo=shp%nb_haplo+1
               do i=1,shp%dataset%params%longhap
                  shp%haplo(shp%nb_haplo,i)=shp%spt%genotypm(c,shp%loc_haplo(i),geno,j)
               end do !i
!! Si les parents de la mère sont connus (repmere.ne.INT_NOT_DEFINED) et sont 
!! soit de meme race soitgénotypés sur au moins un marqueur
           if (raceFileDefined) then
               if (dg%repmere(jm).eq.INT_NOT_DEFINED) then
                   call stop_application ("Both parents of sires and dams have to be defined in the ped file "//&
           "to use breed origin of haplotypes. the parents of dam  "//trim(dg%mere(jm))//&
                    " are not in the ped file!")
               else          
                  if ( ( dgr%racep(dg%repmere(jm))==dgr%racem(dg%repmere(jm)) ) &
                  .or. ( shp%spt%phasm(c,jm))) then
                    if (j==1)  shp%race_h(shp%nb_haplo)= dgr%racep(dg%repmere(jm))
                    if (j==2)  shp%race_h(shp%nb_haplo)= dgr%racem(dg%repmere(jm))
                end if
               end if
             end if
             end do !j
           end do !geno
! 
!  cas des femelles non phasees (issues d'une seule race lorsqu'elle est renseignée)
!
         else 
            geno=0
            do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
              if(dga%presentg(c,kd)) then
                call shp%local_gamete(ip,jm,geno,kd,c,dx)
                 do j_gam=1,shp%nb_gam(kd)
                    shp%nb_haplo=shp%nb_haplo+1
                    do i=1,shp%dataset%params%longhap
                      shp%haplo(shp%nb_haplo,i)=shp%gamete(2,i,j_gam)
                    end do !i
                    if (raceFileDefined) then
                       if (dg%repmere(jm).eq.INT_NOT_DEFINED) then
                          call stop_application ("Both parents of sires and dams have to be defined in the ped file "//&
                          "to use breed origin of haplotypes. the parents of dam  "//trim(dg%mere(jm))//&
                          " are not in the ped file!")
                          STOP
                       else          
                        if (dgr%racep(dg%repmere(jm))==dgr%racem(dg%repmere(jm))) &
                          shp%race_h(shp%nb_haplo)= dgr%racep(dg%repmere(jm))
                       endif
                    endif
                end do !j_gam
              end if !presentg
            end do !kd
         end if

         end do !jm
        end do !ip
         end if ! hdam
!ATTENTION: lorsqu'un haplotype n'a pas une race connue -->le programme ne continue pas
! normalement cette erreur a été detectée avant
        
         do j=1,shp%nb_haplo
            if (raceFileDefined.and.shp%race_h(j)=='UNKNOWN') then
              call stop_application("An haplotype has an unknown breed origin, PLEASE check if you give"//&
              " the breeds of all funder parents!")
            endif
           ! print *, 'HAPLO_pop:',j, race_h(j), (shp%haplo(j,i),i=1,dataset%params%longhap)
         enddo
         end subroutine liste_haplo
!!***

!!****f* m_qtlmap_haplotype_ldla/liste_haplo_complet
!! NAME
!!    liste_haplo_complet
!! DESCRIPTION
!!
!! NOTES
!! liste_haplo_complet etablit la liste des haplotypes complets possible
!!
!!  en entree, il lui faut
!!	loc_haplo , dataset%params%longhap  et loc_haplo  créés par la routine liste_haplo du module m_qtlmap_haplotype_ldla
!!       nall et alleles créés par la routine set_allele_info_vector du module m_qtlmap_genotype
!!  en sortie il fournit
!!	tab_haplo_complet et shp%nb_haplo_complet
!! SOURCE
     subroutine liste_haplo_complet(shp,c)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
      integer       , intent(in) :: c
      integer i1,i2,j,j1,j2
      integer nb_row,nb_col
      integer  :: stat
     ! type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_RACE) , pointer    :: dgr

      dgr => shp%dataset%geneaRace
      dga => shp%dataset%genoAnimal
      !dg => dataset%genea

! allocation

      shp%nb_haplo_complet = 1
      do j=1,shp%dataset%params%longhap
        shp%nb_haplo_complet = shp%nb_haplo_complet * dga%nall(c,shp%loc_haplo(j))
      end do !j


      if ( (size(shp%tab_haplo_complet,1)<shp%nb_haplo_complet*dgr%nb_races) ) then
        deallocate (shp%tab_haplo_complet)
        allocate (shp%tab_haplo_complet(shp%nb_haplo_complet*dgr%nb_races ,shp%dataset%params%longhap))
      end if



      shp%tab_haplo_complet = dga%nmanque


      do i1 = 1,dga%nall(c,shp%loc_haplo(1))
        shp%tab_haplo_complet(i1,1) = set_pheno(dga,shp%dataset%map,c,dga%alleles(c,shp%loc_haplo(1),i1))
      end do !i1


      nb_row=dga%nall(c,shp%loc_haplo(1))
      nb_col=1
 
      do j=2,shp%dataset%params%longhap

        do i1 = 1,nb_row
          shp%tab_haplo_complet(i1,nb_col+1) = set_pheno(dga,shp%dataset%map,c,dga%alleles(c,shp%loc_haplo(j),1))
        end do !i1

        do j1= 2,dga%nall(c,shp%loc_haplo(j))

          do i2 = 1,nb_row
            do j2= 1,nb_col
              shp%tab_haplo_complet(i2+(j1-1)*nb_row,j2) = shp%tab_haplo_complet(i2,j2)
            end do !j2
              shp%tab_haplo_complet(i2+(j1-1)*nb_row,nb_col+1) = set_pheno(dga,shp%dataset%map,c,dga%alleles(c,shp%loc_haplo(j),j1))
          end do !i2

        end do !j1
          nb_row=nb_row * dga%nall(c,shp%loc_haplo(j))
          nb_col= nb_col+1
      end do !j

 
      
      end subroutine liste_haplo_complet
!!***

!!****f* m_qtlmap_haplotype_ldla/proba_haplo_complet
!! NAME
!!    proba_haplo_complet
!! DESCRIPTION
!!
!! NOTES
!!  proba_haplo_complet etablit les probabilites des haplotypes complets
!!
!! en entree, il lui faut
!!	shp%nb_haplo_complet et tab_haplo_complet créés par la routine liste_haplo_complet du module m_qtlmap_haplotype_ldla
!!       haplo et shp%nb_haplo créés par la routine set_allele_info_vector du module m_qtlmap_genotype
!! en sortie il fournit pb_haplo_complet
!!
!! HISTORY
!!   correction BUG : leak memory 05/11/2010
!! SOURCE
     subroutine proba_haplo_complet(shp,c,dx)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
      integer       , intent(in) :: c
      real(kind=dp) , intent(in) :: dx
      double precision erreur
      integer iht,jhc,i,iter,ri
      logical compatible
      integer  :: stat
      real(kind=dp)  :: S
! dim allocatable
 
      real(kind=dp)     , dimension (:)      , allocatable    :: phc
      real(kind=dp)     , dimension (:)      , allocatable    :: pht
      type(GENEALOGY_RACE) , pointer    :: dgr
      type(GENOTYPE_BASE) , pointer :: dga

      dga => shp%dataset%genoAnimal
      dgr => shp%dataset%geneaRace

! local allocation

      allocate (phc(shp%nb_haplo_complet*dgr%nb_races))
      allocate (pht(shp%nb_haplo))

      !resize if too small
      if ( size(shp%comp_haplo,1) < shp%nb_haplo .or. size(shp%comp_haplo,2) < shp%nb_haplo_complet*dgr%nb_races ) then
         deallocate (shp%comp_haplo)
         allocate (shp%comp_haplo(shp%nb_haplo,shp%nb_haplo_complet*dgr%nb_races))
      end if

     if ( size(shp%pb_haplo_complet) < shp%nb_haplo_complet*dgr%nb_races ) then
         deallocate (shp%pb_haplo_complet)
         !deallocate (shp%tab_haplo_complet)
         deallocate (shp%race_haplo_complet)
         allocate (shp%pb_haplo_complet(shp%nb_haplo_complet*dgr%nb_races))
        ! allocate (shp%tab_haplo_complet(shp%nb_haplo_complet*dgr%nb_races,shp%dataset%params%longhap))
         allocate (shp%race_haplo_complet(shp%nb_haplo_complet*dgr%nb_races))
     end if

 
      shp%pb_haplo_complet = 0.d0
      shp%comp_haplo=.false.
      shp%race_haplo_complet='UNKNOWN'

      do ri=1,dgr%nb_races
         do jhc = 1,shp%nb_haplo_complet
            do i=1,shp%dataset%params%longhap
              shp%tab_haplo_complet(jhc+shp%nb_haplo_complet*(ri-1),i)=shp%tab_haplo_complet(jhc,i)
              shp%race_haplo_complet(jhc+shp%nb_haplo_complet*(ri-1))=dgr%nom_race(ri)
            enddo
          enddo
       enddo
       shp%nb_haplo_complet= shp%nb_haplo_complet*dgr%nb_races

!
!  creation de la table de compatibilites entre les haplotypes observes et les haplotypes complets possibles
!
      do iht=1,shp%nb_haplo
        do jhc = 1,shp%nb_haplo_complet
          compatible =.true.
         !haplotypes de races différentes
          if (shp%race_haplo_complet(jhc) /= shp%race_h(iht)) then
              compatible = .false.
              !exit
          else
         !haplotypes différents de meme race
            do i=1,shp%dataset%params%longhap
              if( (shp%haplo(iht,i) /= shp%tab_haplo_complet(jhc,i)) .and. (shp%haplo(iht,i) /= dga%nmanque) )then
                compatible = .false.
                !exit
              end if
            end do !i
          endif
          if(compatible) shp%comp_haplo(iht,jhc)=.true.
        end do !jhc
      end do !iht


!
!  Algorithme EM pour l'estimation des probabilites des haplotypes complets
!
      shp%pb_haplo_complet =1.d0 / dble(shp%nb_haplo_complet)

      erreur=1.d0
      iter=0

      do while( erreur > EPS_EM)

        pht=0.d0
        do iht = 1, shp%nb_haplo
          do jhc=1,shp%nb_haplo_complet
           if(shp%comp_haplo(iht,jhc)) pht(iht)=pht(iht)+shp%pb_haplo_complet(jhc)
          end do !jhc
        end do !iht
 

        erreur=0.d0
        phc=0.d0 
        do jhc=1,shp%nb_haplo_complet
          do iht = 1, shp%nb_haplo
            if(shp%comp_haplo(iht,jhc)) phc(jhc)=phc(jhc)+shp%pb_haplo_complet(jhc)/pht(iht)
          end do !iht
          phc(jhc)=phc(jhc)/dble(shp%nb_haplo)
          erreur=erreur+dabs(phc(jhc)-shp%pb_haplo_complet(jhc))
        end do !jhc

        iter=iter+1
        if ( iter > ITER_EM_MAX) then
          call stop_application( ' the EM algorithm for the estimation of haplotypes frequency does not converge ')
        end if
       do jhc=1,shp%nb_haplo_complet
          shp%pb_haplo_complet(jhc)=phc(jhc)
       end do !jhc
      end do !while
        S=0.D0
      ! print *,'RACE       PB_HAPLO' 
       do jhc=1,shp%nb_haplo_complet
          S=S+shp%pb_haplo_complet(jhc)
          !print *,trim(shp%race_haplo_complet(jhc)), shp%pb_haplo_complet(jhc)
       end do !jhc
       ! print *, 'Somme PROB_haplo_complet', S, shp%nb_haplo_complet

      deallocate (phc)
      deallocate (pht)

      end subroutine proba_haplo_complet
!!***

!!****f* m_qtlmap_haplotype_ldla/tri_haplo_complet
!! NAME
!!    tri_haplo_complet
!! DESCRIPTION
!!
!! NOTES
!! tri_haplo_complet elimine de la liste les haplotypes complets peu probables
!!
!! en entree, il lui faut
!!	shp%nb_haplo_complet et tab_haplo_complet créés par la routine liste_haplo_complet du module m_qtlmap_haplotype_ldla
!!       haplo et shp%nb_haplo créés par la routine set_allele_info_vector du module m_qtlmap_genotype
!! en sortie il fournit
!!	pb_haplo_complet recalcule en, eliminant les haplotypes peu probables
!!        name_haplo_complet le nom compacte des haplotypes complets
!! HISTORY
!!  correction BUG : leak memory 05/11/2010
!! SOURCE
     subroutine tri_haplo_complet(shp,c,dx)
       class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
       integer       , intent(in) :: c
       real(kind=dp) , intent(in) :: dx
       integer  :: stat,lstr1,lstr2
       integer  jhc,i,idx,ipb,khr,iht,ir
       integer nb_possible
       double precision tphc
       logical une_race_rare
       character(kind=1) :: nm_race_rare
       real(kind=dp) :: spb_haplo_reduit
       real(kind=dp) :: pb,pbtot,spb, spb_r !cmo
! dim allocatable
     double precision          , dimension (:)    , allocatable    :: phc
     double precision          , dimension (:)    , allocatable    :: pbtot_haplo_complet !cmo
     logical :: raceFileDefined
     type(GENEALOGY_RACE) , pointer    :: dgr
     type(GENEALOGY_BASE) , pointer :: dg
     type(GENOTYPE_BASE) , pointer :: dga

     dga => shp%dataset%genoAnimal
     dg => shp%dataset%genea
     dgr => shp%dataset%geneaRace

     raceFileDefined = shp%dataset%params%get_file_val(K_RACE) /= ""

! allocation
      allocate (phc(shp%nb_haplo_complet*dgr%NB_RACES))

      if ( size(shp%liste_haplo_reduit,1) < shp%nb_haplo_complet*dgr%NB_RACES ) then
        deallocate (shp%liste_haplo_reduit)
        allocate (shp%liste_haplo_reduit(shp%nb_haplo_complet*dgr%NB_RACES))
      end if

      allocate (pbtot_haplo_complet(dgr%NB_RACES)) !CMO

      !correction BUG OFI : 14/09/2010
      if ( size (shp%name_haplo_complet)<shp%nb_haplo_complet*dgr%NB_RACES ) then
        deallocate(shp%name_haplo_complet)
        allocate (shp%name_haplo_complet(shp%nb_haplo_complet*dgr%NB_RACES))
      end if

!
!  nom compacte des haplotypes mis dans shp%name_haplo_complet
!
      shp%name_haplo_complet=''

      do jhc=1,shp%nb_haplo_complet
        lstr1=len(trim(get_pheno(dga,c,shp%tab_haplo_complet(jhc,1))))
        shp%name_haplo_complet(jhc)(1:lstr1)=trim(get_pheno(dga,c,shp%tab_haplo_complet(jhc,1)))

        do i=2,shp%dataset%params%longhap
          lstr2=len(trim(get_pheno(dga,c,shp%tab_haplo_complet(jhc,i))))
          shp%name_haplo_complet(jhc)(lstr1+1:lstr1+lstr2)=trim(get_pheno(dga,c,shp%tab_haplo_complet(jhc,i)))
          lstr1=lstr1+lstr2
        end do !i

      end do ! jhc
!
!  comptage des haplotypes complets poossibles par haplotype trouvé
!
        do iht =1,shp%nb_haplo
         shp%nb_haplo_possible(iht)=0
         do jhc=1,shp%nb_haplo_complet
           if(shp%comp_haplo(iht,jhc)) shp%nb_haplo_possible(iht)=shp%nb_haplo_possible(iht)+1
         end do ! jhc
      end do !iht 

!!!CMO
   ! Recalcul de pb_haplo_complet(jhc): Calcul des freq haplotypique intrarace
        pbtot_haplo_complet=0.d0
        pbtot=0.d0
            do ir=1, dgr%nb_races
              do jhc=1,shp%nb_haplo_complet
                    if (shp%race_haplo_complet(jhc)==dgr%nom_race(ir)) &
                    pbtot_haplo_complet(ir)=pbtot_haplo_complet(ir)+shp%pb_haplo_complet(jhc)
                enddo
                pbtot=pbtot+pbtot_haplo_complet(ir)
                !print *, 'PB_race',ir, pbtot_haplo_complet(ir)
             enddo
             !print *, 'PBtot', pbtot
!!!CMO

!  regroupement des haplotypes rares dans le paquet des donnees manquantes
!  on prepare le calcul des probabilites des haplotypes non rares (phc(jhc))
!  et on les compte (shp%nb_haplo_reduit)
!
        phc=0.d0
        tphc=0.d0
        i=0
        shp%nb_haplo_reduit=0
        shp%liste_haplo_reduit=0
        spb_haplo_reduit=0.D0
        une_race_rare=.true.
        nm_race_rare='UNKNOWN'
        !print *,PROB_HAPLO_MIN
        do jhc=1,shp%nb_haplo_complet
!!!CMO
             pb=0.d0
              do ir=1, dgr%nb_races
                if (shp%race_haplo_complet(jhc)==dgr%nom_race(ir)) then
                    pb=shp%pb_haplo_complet(jhc)
                    shp%pb_haplo_complet(jhc)=(pb*pbtot)/(pbtot_haplo_complet(ir))
                    !print *, shp%race_haplo_complet(jhc), 'AVANT ',pb, 'APRES ', shp%pb_haplo_complet(jhc)
                endif
              enddo
!!!CMO
          if(shp%pb_haplo_complet(jhc) > shp%dataset%params%PROB_HAPLO_MIN) then
            phc(jhc)=shp%pb_haplo_complet(jhc)
            tphc=tphc+phc(jhc)
            shp%nb_haplo_reduit= shp%nb_haplo_reduit+1
            shp%liste_haplo_reduit(shp%nb_haplo_reduit)=jhc
            spb_haplo_reduit=spb_haplo_reduit+phc(jhc)
          else
              
              if (shp%race_haplo_complet(jhc).ne.'UNKNOWN') i=i+1
              if (i==1) nm_race_rare=shp%race_haplo_complet(jhc)
              if(nm_race_rare.ne.'UNKNOWN'.and.shp%race_haplo_complet(jhc).ne.nm_race_rare) &
                une_race_rare=.false.   
          end if
        end do ! jhc 
 !print *,'somme des p_haplo=',    spb_haplo_reduit , shp%nb_haplo_complet
! 
!  arret si aucun haplotype possible
!
        if ( tphc == 0.d0) then
          return
          idx=int(100.*dx)
          ipb=int(100.*shp%dataset%params%PROB_HAPLO_MIN)
          call stop_application(' No possible haplotypes found at location '//trim(str(idx))//&
                                ' cM with a threshold '//trim(str(ipb))//' %')
        end if
!
! allocations selon le nobre d'haplotypes non rare
! 
      if ( size(shp%comp_haplo_reduit,1) < shp%nb_haplo .or. size(shp%comp_haplo_reduit,2) < shp%nb_haplo_reduit+dgr%nb_races ) then
       deallocate (shp%comp_haplo_reduit)
       allocate (shp%comp_haplo_reduit(shp%nb_haplo,shp%nb_haplo_reduit+dgr%nb_races))
       shp%comp_haplo_reduit=.false.
      end if
 
      if ( size(shp%pb_haplo_reduit,1) < shp%nb_haplo_reduit+dgr%nb_races ) then
        deallocate (shp%pb_haplo_reduit)
        allocate (shp%pb_haplo_reduit(shp%nb_haplo_reduit+dgr%nb_races))
      end if

      if ( size(shp%pb_haplo_reduit_r,1) < shp%nb_haplo_reduit+dgr%nb_races ) then
        deallocate (shp%pb_haplo_reduit_r)
        allocate (shp%pb_haplo_reduit_r(shp%nb_haplo_reduit+dgr%nb_races))
      end if

      if ( size(shp%race_haplo_reduit,1) < shp%nb_haplo_reduit+dgr%nb_races ) then
        deallocate (shp%race_haplo_reduit)
        allocate (shp%race_haplo_reduit(shp%nb_haplo_reduit+dgr%nb_races))
      end if

      if ( size(shp%name_haplo_reduit,1) <  shp%nb_haplo_reduit+dgr%nb_races ) then
        deallocate (shp%name_haplo_reduit)
        allocate (shp%name_haplo_reduit(shp%nb_haplo_reduit+dgr%nb_races))
      end if

!
!  creation des probabilites et correspondances
!
      shp%pb_haplo_reduit =0.d0
      shp%pb_haplo_reduit_r =0.d0
      shp%race_haplo_reduit='UNKNOWN'

          do khr=1,shp%nb_haplo_reduit
!         shp%pb_haplo_reduit(khr)=phc(shp%liste_haplo_reduit(khr))/tphc
          shp%name_haplo_reduit(khr)=shp%name_haplo_complet(shp%liste_haplo_reduit(khr))
          shp%race_haplo_reduit(khr)=shp%race_haplo_complet(shp%liste_haplo_reduit(khr))
! Calcul of probability of frequent haplotypes 
          if (raceFileDefined) then
            do ir=1, dgr%nb_races
              if (shp%race_haplo_reduit(khr)==dgr%nom_race(ir)) &
! probability of haplotype wathever the breed origin
              shp%pb_haplo_reduit(khr)=phc(shp%liste_haplo_reduit(khr))*pbtot_haplo_complet(ir)
! probability of haplotype given the breed origin
              shp%pb_haplo_reduit_r(khr)=phc(shp%liste_haplo_reduit(khr))
            enddo
          else
           shp%pb_haplo_reduit(khr)=phc(shp%liste_haplo_reduit(khr))
           shp%pb_haplo_reduit_r(khr)=phc(shp%liste_haplo_reduit(khr))
          endif
      end do ! khr

      do iht =1,shp%nb_haplo
         do ir=1,dgr%nb_races
         if (shp%race_h(iht)==dgr%nom_race(ir)) then
            nb_possible=0
            do khr=1,shp%nb_haplo_reduit
               shp%comp_haplo_reduit(iht,khr)=shp%comp_haplo(iht,shp%liste_haplo_reduit(khr))
               if(shp%comp_haplo_reduit(iht,khr) &
               .and.shp%race_h(iht)==shp%race_haplo_complet(shp%liste_haplo_reduit(khr))) &
               nb_possible=nb_possible+1
            end do ! khr
! mise à true de shp%comp_haplo_reduit pour les haplotypes rares
            if(nb_possible < shp%nb_haplo_possible(iht)) &
             shp%comp_haplo_reduit(iht,shp%nb_haplo_reduit+ir)=.true.
         endif
         enddo !ir
        !shp%name_haplo_reduit(shp%nb_haplo_reduit+1)='race1?'
        !shp%name_haplo_reduit(shp%nb_haplo_reduit+2)='race2?'
        !do khr=1,shp%nb_haplo_reduit+dgr%nb_races
            !print *, iht, race_h(iht), khr,trim(shp%name_haplo_reduit(khr)),&
           ! shp%comp_haplo_reduit(iht,khr)
         !enddo
      end do !iht    

! Calcul of probability of rare haplotypes 
      if(shp%nb_haplo_reduit < shp%nb_haplo_complet) then
         if (raceFileDefined) then
             do ir=1, dgr%nb_races
               spb_haplo_reduit=0.d0
               shp%nb_haplo_reduit= shp%nb_haplo_reduit+1
               do jhc=1,shp%nb_haplo_complet
                 if ((shp%race_haplo_complet(jhc)==dgr%nom_race(ir)) &
                 .and.(shp%pb_haplo_complet(jhc).LT.shp%dataset%params%PROB_HAPLO_MIN)) &
                  spb_haplo_reduit=spb_haplo_reduit+shp%pb_haplo_complet(jhc)
                end do !jhc
! probability of haplotype given the breed origin
               shp%pb_haplo_reduit_r(shp%nb_haplo_reduit)= spb_haplo_reduit
! probability of haplotype wathever the breed origin
               shp%pb_haplo_reduit(shp%nb_haplo_reduit)= spb_haplo_reduit*pbtot_haplo_complet(ir)
              shp%name_haplo_reduit(shp%nb_haplo_reduit)='RARE'
               shp%race_haplo_reduit(shp%nb_haplo_reduit)=dgr%nom_race(ir)

            enddo
         else
          shp%nb_haplo_reduit= shp%nb_haplo_reduit+1
          shp%pb_haplo_reduit(shp%nb_haplo_reduit)=1.d0-tphc
          shp%pb_haplo_reduit_r(shp%nb_haplo_reduit)=1.d0-tphc
          shp%name_haplo_reduit(shp%nb_haplo_reduit)='RARE'
          if (une_race_rare) shp%race_haplo_reduit(shp%nb_haplo_reduit)=nm_race_rare
         end if
       end if

       shp%nb_max_haplo=shp%nb_haplo_reduit
 
!      write(nficout,*)  'Num_HAPLO   RACE    HAPLO    pb_HAPLO'
!      do i= 1, shp%nb_haplo_reduit
!         write(nficout,*) i,trim(shp%race_haplo_reduit(i)),' ',&
 !                   trim(shp%name_haplo_reduit(i)), &
!                    shp%pb_haplo_reduit(i)
 !     enddo

       !  correctiopn BUG : leak memory 05/11/2010
       deallocate (phc)
       deallocate (pbtot_haplo_complet)

      end subroutine tri_haplo_complet
!!***

!!****f* m_qtlmap_haplotype_ldla/set_haplo_final
!! NAME
!!    set_haplo_final
!! DESCRIPTION
!!
!! NOTES
!! set_haplo_final met au propre les listes de correspondances des haplotypes paternels et maternels
!!  shp%num_haplo_pere(ip,j,l)  pour le pere ip, chromosome j et haplotype complet de la liste reduite possible l
!!   (avec = 1,... nb_gam_pere(ip,j))
!!  shp%num_haplo_mere(geno,j,l)  pour la combinaison mere x phase geno, chromosome j et haplotype complet de la liste
!!    reduite possible l  (avec = 1,... nb_gam_mere(geno,j))
!!  shp%num_haplo_desc(kd,j,l)  pour le descendant kd, chromosome j et haplotype complet de la liste
!!    reduite possible l  (avec = 1,... nb_gam_desc(kd))
!!
!!  Utilise np,nm,ndm,estfem,repfem de m_qtlmap_data
!!  utilise shp%nb_haplo_reduit et shp%comp_haplo_reduit  de m_qtlmap_haplotype_ldla (tri_haplo_complet)
!!  utilise ngenom de m_qtlmap_haplotpe_data
!!  utilise presentg
!!  utilise nb_gam
!! SOURCE
     subroutine set_haplo_final(shp,c,n,hsire,hdam)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

      logical, intent (in) :: hsire,hdam
      integer, intent (in) :: n,c
      integer  :: stat
! divers
      integer ip,j,j_haplo,k_haplo,jm,geno,kd,j_gam, ir
      integer  :: nb_gam_t,i_gam,value_alloc,value_alloc_2
      real (kind=dp)  :: su,r
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_RACE) , pointer    :: dgr

      dgr => shp%dataset%geneaRace
      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea

      if ( size(shp%num_haplo_pere,3) <  shp%nb_haplo_complet ) then
        deallocate (shp%num_haplo_pere)
        allocate (shp%num_haplo_pere(dg%np,2,shp%nb_haplo_complet))
        end if
      shp%num_haplo_pere=0

      if ( size(shp%num_haplo_mere,3) <  shp%nb_haplo_complet ) then
        deallocate (shp%num_haplo_mere)
        allocate (shp%num_haplo_mere(maxval(shp%spt%ngenom(:,dg%nm+1)),2,shp%nb_haplo_complet))
      end if

      shp%num_haplo_mere = 0
      value_alloc=shp%nb_haplo_reduit*maxval(shp%nb_gam)
      ! OFI Janv 2013 => si plus de deux alleles ca peut planter car pas assez de memoire...a corriger
      value_alloc_2=2**(2*shp%dataset%params%longhap)+10

      if ( size(shp%num_haplo_desc,2) <  value_alloc ) then
        deallocate (shp%num_haplo_desc)
        allocate (shp%num_haplo_desc(dg%nd,value_alloc))
      end if

      shp%num_haplo_desc = 0

      if ( size(shp%pb_haplo_desc,2) <  value_alloc_2) then
        deallocate (shp%pb_haplo_desc)
        allocate (shp%pb_haplo_desc(dg%nd,value_alloc_2))
      end if

      shp%pb_haplo_desc = 0.d0

! On saute les haplotypes manquants
       j_haplo=0
       do ir=1,dgr%nb_races
          j_haplo=j_haplo+1
       enddo

!
!  pour les peres
!
       shp%nb_gam_pere=0
       if(hsire)then

       do ip=1,dg%np
         do j=1,2
           j_haplo=j_haplo+1
           do k_haplo=1,shp%nb_haplo_reduit
             if(shp%comp_haplo_reduit(j_haplo,k_haplo)) then
               shp%nb_gam_pere(ip,j)=shp%nb_gam_pere(ip,j)+1
               shp%num_haplo_pere(ip,j,shp%nb_gam_pere(ip,j))= k_haplo
             end if
           end do !k_haplo
          end do !j
        end do !ip

       end if ! opt_model(hsire)
!
!  pour les meres 
!
       shp%nb_gam_mere=0
       shp%nb_gam_desc=0
       if(hdam)then

       do jm=1,dg%nm
         if(dga%estfem(dg%repfem(jm)))then
!
!  meres phasees
!
           do geno=shp%spt%ngenom(c,jm)+1,shp%spt%ngenom(c,jm+1)
             do j=1,2
               j_haplo=j_haplo+1
                 do k_haplo=1,shp%nb_haplo_reduit
                   if(shp%comp_haplo_reduit(j_haplo,k_haplo)) then
                     shp%nb_gam_mere(geno,j)=shp%nb_gam_mere(geno,j)+1
                     shp%num_haplo_mere(geno,j,shp%nb_gam_mere(geno,j))= k_haplo
                   end if
                 end do !k_haplo
             end do !j
           end do !geno
         else
! 
!   non phasees
!           do kd=ndm(jm)+1,ndm(jm+1)
!             if(presentg(c,kd))then
!             do j_gam=1,shp%nb_gam(kd)
!               j_haplo=j_haplo+1
!                 do k_haplo=1,shp%nb_haplo_reduit
!                   if(shp%comp_haplo_reduit(j_haplo,k_haplo)) then
!                     nb_gam_desc(kd,j_gam)=nb_gam_desc(kd,j_gam)+1
!                     shp%num_haplo_desc(kd,j_gam,nb_gam_desc(kd,j_gam))= k_haplo
!                   end if
!                 end do !k_haplo
!             end do !j_gam
!             end if
!            end do !kd

           do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
             if(dga%presentg(c,kd))then
               do j_gam=1,shp%nb_gam(kd)
                 j_haplo=j_haplo+1
                 nb_gam_t=shp%nb_gam_desc(kd)+1
                 do k_haplo=1,shp%nb_haplo_reduit
                   if(shp%comp_haplo_reduit(j_haplo,k_haplo)) then
                     shp%nb_gam_desc(kd)=shp%nb_gam_desc(kd)+1
                     shp%num_haplo_desc(kd,shp%nb_gam_desc(kd))=k_haplo
                   end if
                 end do  !k
                 su = dble(sum(shp%pb_haplo_reduit(shp%num_haplo_desc(kd,nb_gam_t:shp%nb_gam_desc(kd)))))
                 if ( su == 0 ) cycle
                 do i_gam=nb_gam_t,shp%nb_gam_desc(kd)
 !                   shp%pb_haplo_desc(kd,i_gam)  = shp%pb_haplo_desc(kd,i_gam)  &
 !                       +shp%prob_gam(kd,j_gam) * (shp%pb_haplo_reduit(shp%num_haplo_desc(kd,i_gam)) /su)
                   shp%pb_haplo_desc(kd,shp%num_haplo_desc(kd,i_gam))  =  &
                         shp%pb_haplo_desc(kd,shp%num_haplo_desc(kd,i_gam))  &
                        +shp%prob_gam(kd,j_gam) * (shp%pb_haplo_reduit(shp%num_haplo_desc(kd,i_gam)) /su)
                 end do !i_gam
               end do !j_gam
             end if !presentg(c,kd)
           end do !kd
         end if ! estfem
       end do !jm
       end if ! opt_model (hdam)

     end  subroutine set_haplo_final
!!***

!!****f* m_qtlmap_haplotype_ldla/local_gamete
!! NAME
!!    local_gamete
!! DESCRIPTION
!!
!! NOTES
!!  local_gamete trouve les gametes parentaux transmis par ses parents (ip et jm)
!!  pour une phase geno de la mere
!!  a un individu particulier (kd)
!!  autour d'un point precis defini par le chromosome (c) et la position (dx)
!!  ainsi que le nombre de marqueurs l'entourant (dataset%params%longhap)
!!
!!  en sortie, local_gamete fournit
!!	nb_gam = le nombre de couples de gametes parentaux (2 ** nb indeterminations)
!!	gamete(1,lk,j_gam) l allele paternel au lk eme marqueur pour le j_gam eme gamete
!!	gamete(2,lk,j_gam) l allele maternel au lk eme marqueur pour le j_gam eme gamete
!!   prob_gam(j_gam) la probabilite de cet allele
!!
!!  Utilise posi,nmk,nmanque,correp,corred,correm,pheno     de m_qtlmap_data
!!  Utilise genotyp,genotypm                                de m_qtlmap_haplotype_data
!!
!!  HISTORY
!!    correction BUG : leak memory 05/11/2010
!! SOURCE
   subroutine local_gamete(shp,ip,jm,geno,kd,c,dx)
      class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)  :: shp
      integer , intent(in) ::ip,jm,geno,kd,c
      real (kind=dp), intent(in) :: dx

! divers
      integer lk,lkleft,lkright,l,llk
      integer left_known_sire,left_known_dam,right_known_sire,right_known_dam
      integer lks,lkd,rks,rkd
      integer isex,igam
      integer i,j,k,idup
      integer ngamkept,nb_inc
      integer value_alloc
      logical known_sire_left,known_dam_left,known_sire_right,known_dam_right,manquant
      logical ksr,ksl,kdr,kdl
      integer (kind=KIND_PHENO) :: allele(2),unknown
      integer :: orig(2)
      real(kind=dp) recm,recf,probtot
      integer  :: stat

! dimension allocatable
      integer                   , dimension (:,:)    , allocatable   :: base
      integer(kind=KIND_PHENO)  , dimension (:,:,:)  , allocatable   :: origine
      integer                   , dimension (:,:)    , allocatable   :: tabdup
      type(GENEALOGY_RACE) , pointer    :: dgr
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(MAP_BASE) , pointer    :: map

      map => shp%dataset%map
      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea
      dgr => shp%dataset%geneaRace
! allocation
!
      allocate (base(2,map%nmk(c)))
      value_alloc= 2**(2*shp%dataset%params%longhap)

      allocate (origine(2,map%nmk(c),value_alloc))
      allocate (tabdup(value_alloc,2*shp%dataset%params%longhap))
      unknown=0

      shp%gamete=dga%nmanque

!  la position testee doit etre compatible avec la longueur de l'haplotype
!

      if ( ceiling(shp%dataset%params%longhap/2.d0) > 0 ) then
       if((dx < map%posi(c,ceiling(shp%dataset%params%longhap/2.d0))) .or. &
         (dx > map%posi(c,map%nmk(c)-ceiling(shp%dataset%params%longhap/2.d0)+1)))then
       call stop_application("local_gamete : incompatibilty between haplotype length and tested position.")
       end if
      end if
!
! on cherche le marqueur flanquant a gauche la position dx
! 
       lk=1
       do while (dx > map%posi(c,lk))
         lk=lk+1
       end do
      lk=lk-1
!
! lkleft et lkright seront les numeros des marqueurs bornant l haplotype a determiner
!
       lkleft  =lk-ceiling(shp%dataset%params%longhap/2.d0)+1
       if(lk <= ceiling(shp%dataset%params%longhap/2.d0)) lkleft=1
       lkright =lk+ceiling(shp%dataset%params%longhap/2.d0)
       if(lk >= (map%nmk(c)-ceiling(shp%dataset%params%longhap/2.d0))) lkright=map%nmk(c)
!
!  du fait des incertitudes 
!      (pere et descendant heterozygotes identiques si la mere n est pas marquee
!       deux parents et descendants heterozygotes identiques si la mere est marquee)
!     plusieurs gametes peuvent avoir forme un descendant
!
!  nb_gam = 2**(nb d'incertitudes) est le nombre de gametes possibles
!  gamete(1,lk,j_gam) est le nom de l'allele au i eme marqueur transmis par le pere 
!    (par la mere pour gamete(2,lk,j_gam))
!     pour la j_gam eme configuration gametique
!  origine(1,lk,j_gam) en est l origine grand paternel 
!     (0 pour inconnue, 1 si grand pere, 2 si grand mere, 3 le pere est homozygote, 4 si manquant)
!   (grand maternelle pour origine(2,lk))
!
!  construction des gametes correspondant a l haplotype
!  pour chaque marqueur de l'intervalle, on cherche (avec point_gamete) le nom des alleles recus
!  de chacun des parents, ainsi que leur origine grand parentale (mis dans base)
! base vaut 1 si origine = grand pere; 2 = grand mere ; 
!           0 = indetermine (duo ou triplet heterozygote identique); 3= parent homozygote, 4= manquant
!
       origine =0
       base=0

      manquant=.false.
       do lk = lkleft,lkright
         call shp%point_gamete(ip,jm,geno,kd,c,lk,allele,base(:,lk))
         if(allele(1) == dga%nmanque .or. allele(2) == dga%nmanque) manquant=.true.
          origine(1,lk,1)=base(1,lk)
          origine(2,lk,1)=base(2,lk)
       end do !lk
 
!
!  je mets a manque total les haplotypes compremnant une donnee manquante
!
     if(manquant) then
       shp%nb_gam(kd)=1
       shp%prob_gam(kd,1)=1.d0

       llk=0
       do lk = lkleft,lkright
         llk=llk+1
         do isex=1,2
           shp%gamete(isex,llk,1)=dga%nmanque
         end do ! isex
       end do !lk

       deallocate (base)
       deallocate (origine)
       deallocate (tabdup)

       return

      end if
!
! il y a 2**(nombre d'origines inconnues) gametes possibles dont on doit calculer les probabilites
!
! quand les marqueurs flanquants sont d'origine inconnue, il faut rechercher une information 
! a l exterieur de l'haplotype
!
!  vers la gauche
!
       left_known_sire=0
       left_known_dam=0
       known_sire_left=.false.
       if(base(1,lkleft) ==1 .or.base(1,lkleft) == 2) then
         known_sire_left=.true.
         left_known_sire=lkleft
       end if
       known_dam_left=.false.
       if(base(2,lkleft) ==1 .or.base(2,lkleft) == 2) then
         known_dam_left=.true.
         left_known_dam=lkleft
       end if

       if(.not.manquant)then
       if(dga%estfem(dg%repfem(jm))) then

         if(.not. (known_sire_left .and. known_dam_left))  then
           do lk=lkleft-1,1,-1 
             call shp%point_gamete(ip,jm,geno,kd,c,lk,allele,base(:,lk))

          if(.not. known_sire_left .and. (base(1,lk) == 1 .or. base(1,lk) == 2)) then
               known_sire_left=.true.
               left_known_sire=lk 
               origine(1,lk,1)=base(1,lk)
             end if
             if(.not. known_dam_left .and. (base(2,lk) == 1 .or. base(2,lk) == 2)) then
               known_dam_left=.true.
               left_known_dam=lk 
               origine(2,lk,1)=base(2,lk)
             end if

             if(known_sire_left .and. known_dam_left) exit
           end do !lk
         end if  
  
       else  

         if(.not.known_sire_left)  then
           do lk=lkleft-1,1,-1 
            call shp%point_gamete(ip,jm,geno,kd,c,lk,allele,base(:,lk))

             if(.not. known_sire_left .and. (base(1,lk) == 1 .or. base(1,lk) == 2)) then
               known_sire_left=.true.
               left_known_sire=lk 
               origine(1,lk,1)=base(1,lk)
             end if          
             if(known_sire_left) exit
           end do !lk
         end if

       end if
       end if
!
!  vers la droite
!
       right_known_sire=0
       right_known_dam=0
       known_sire_right=.false.
       if(base(1,lkright) ==1 .or.base(1,lkright) == 2) then
         known_sire_right=.true.
         right_known_sire=lkright
       end if
       known_dam_right=.false.
       if(base(2,lkright) ==1 .or.base(2,lkright) == 2) then
         known_dam_right=.true.
         right_known_dam=lkright
       end if

       if(.not.manquant) then
       if(dga%estfem(dg%repfem(jm))) then

         if(.not. (known_sire_right .and. known_dam_right))  then

           do lk=lkright+1,map%nmk(c)

             call shp%point_gamete(ip,jm,geno,kd,c,lk,allele,base(:,lk))


             if(.not. known_sire_right .and. (base(1,lk) == 1 .or. base(1,lk) == 2)) then
               known_sire_right=.true.
               right_known_sire=lk
               origine(1,lk,1)=base(1,lk)
             end if
             if(.not. known_dam_right .and. (base(2,lk) == 1 .or. base(2,lk) == 2)) then
               known_dam_right=.true.
               right_known_dam=lk
               origine(2,lk,1)=base(2,lk)
             end if

             if(known_sire_right .and. known_dam_right) exit
           end do !lk
         end if
       else

         if(.not.known_sire_right)  then
           do lk=lkright+1,map%nmk(c)
             call shp%point_gamete(ip,jm,geno,kd,c,lk,allele, base(:,lk))

             if(.not. known_sire_right .and. ( base(1,lk) == 1 .or.  base(1,lk) == 2)) then
               known_sire_right=.true.
               right_known_sire=lk
               origine(1,lk,1)= base(1,lk)
             end if
             if(known_sire_right) exit
           end do !lk
         end if
       end if
       end if
!
! Construction de la liste des gametes possibles
! en parcourant les marqueurs le long de l haplotype et en dedoublant la liste 
! a chaque fois qu'une incertitude (base =0) est rencontree
! nb_inc est le nombre d'incertitudes
! shp%nb_gam(kd) le nombre de gametes possibles pour kd (shp%nb_gam(kd)=2**nb_inc dans un premier temps)
!
        nb_inc=0

        do lk = lkleft,lkright
          do isex=1,2
            if(base(isex,lk) == 0) nb_inc=nb_inc + 1
          end do !isex
        end do !lk
        shp%nb_gam(kd)=2**nb_inc
!
!tabdup contient la liste des vecteurs d'origines possibles pour les shp%nb_gam(kd) incertitudes
!
       if(nb_inc > 0) then
         tabdup(1,nb_inc)=1
         tabdup(2,nb_inc)=2
         do i=2,nb_inc
           do j=1,2**(i-1)
             tabdup(j,nb_inc-i+1)=1
             tabdup(2**(i-1)+j,nb_inc-i+1)=2
             do k=1,i-1
               tabdup(2**(i-1)+j,nb_inc-k+1)=tabdup(j,nb_inc-k+1)
             end do !k
           end do !j
         end do !i

         do igam=1,shp%nb_gam(kd)
           idup=0
           do lk = lkleft,lkright
             do isex=1,2
               origine(isex,lk,igam)=origine(isex,lk,1)
               if(base(isex,lk) == 0) then
                 idup=idup+1
                 origine(isex,lk,igam)=tabdup(igam,idup)
               end if
             end do ! isex
           end do ! lk
           if(known_sire_left)  origine(1,left_known_sire,igam) =origine(1,left_known_sire,1)
           if(known_sire_right) origine(1,right_known_sire,igam)=origine(1,right_known_sire,1)
           if(known_dam_left)   origine(2,left_known_dam,igam) = origine(2,left_known_dam,1)
           if(known_dam_right)  origine(2,right_known_dam,igam)= origine(2,right_known_dam,1)
         end do ! igam

       end if !nb_gam


     ksr=known_sire_right
     ksl=known_sire_left
     kdr=known_dam_right
     kdl=known_dam_left
     lks=left_known_sire 
     rks=right_known_sire 
     lkd=left_known_dam
     rkd=right_known_dam 
  
!
!  calcul des probabilites des gametes
!
           do igam=1,shp%nb_gam(kd)
            shp%prob_gam(kd,igam)=1.d0
            known_sire_right=ksr
            known_sire_left=ksl
            known_dam_right=kdr
            known_dam_left=kdl
            left_known_sire=lks
            right_known_sire=rks 
            left_known_dam=lkd
            right_known_dam=rkd 

            do lk=lkleft,lkright

             if(origine(1,lk,igam) < 3) then
               if(known_sire_left) then
                 recm=xaldane(map%posim(c,lk)-map%posim(c,left_known_sire))
                 if(origine(1,lk,igam) ==  origine(1,left_known_sire,igam)) then
                   shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * (1.d0 -recm)
                 else
                   shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * recm
                 end if
               else
                 known_sire_left=.true.
               end if
               left_known_sire=lk
             end if

             if(geno /=0 .and. origine(2,lk,igam) < 3) then
               if(known_dam_left) then
                 recf=xaldane(map%posif(c,lk)-map%posif(c,left_known_dam))
                 if(origine(2,lk,igam) ==  origine(2,left_known_dam,igam)) then
                   shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * (1.d0 -recf)
                 else
                   shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * recf
                 end if
               else
                 known_dam_left=.true.
               end if
               left_known_dam=lk
             end if

            end do ! lk
!
! bout du segment
             if(known_sire_right .and. right_known_sire > lkright) then
              recm=xaldane(map%posim(c,right_known_sire)-map%posim(c,lkright))
               if(origine(1,lkright,igam) ==  origine(1,right_known_sire,igam)) then
                 shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * (1.d0 -recm)
               else
                 shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * recm
               end if
              end if
             if(geno /= 0 .and. known_dam_right .and. right_known_dam > lkright) then
              recf=xaldane(map%posif(c,right_known_dam)-map%posif(c,lkright))
               if(origine(2,lkright,igam) ==  origine(2,right_known_dam,igam)) then
                 shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * (1.d0 -recf)
               else
                 shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam) * recf
               end if
              end if

          end do !igam
 
!  au cas ou, standardisation des probabilites
!
       probtot=0.d0
       do igam=1,shp%nb_gam(kd)
        probtot=probtot+shp%prob_gam(kd,igam)
       end do
       do igam=1,shp%nb_gam(kd)
        shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam)/probtot
       end do
!
!  on enrichir les informations sur les gametes pour chacune des configurations possibles
!
      shp%gamete=dga%nmanque
      do igam=1,shp%nb_gam(kd)
        llk=0
        do lk=lkleft,lkright
          llk=llk+1
          if(origine(1,lk,igam) == 4) shp%gamete(1,llk,igam) = dga%nmanque
          if(origine(1,lk,igam) == 3) shp%gamete(1,llk,igam) = shp%spt%genotyp(c,lk,dga%correp(ip),1)
          if(origine(1,lk,igam) == 2) shp%gamete(1,llk,igam) = shp%spt%genotyp(c,lk,dga%correp(ip),2)
          if(origine(1,lk,igam) == 1) shp%gamete(1,llk,igam) = shp%spt%genotyp(c,lk,dga%correp(ip),1)


          if(origine(1,lk,igam) <= 3) then
       
             if(shp%gamete(1,llk,igam) == dga%pheno(c,lk,dga%corred(kd),1)) then
                  shp%gamete(2,llk,igam) = dga%pheno(c,lk,dga%corred(kd),2)
             else
                  shp%gamete(2,llk,igam) = dga%pheno(c,lk,dga%corred(kd),1)
             end if
          end if
        end do ! lk
      end do !igam
!
!
! elimination des cas improbables
!
     ngamkept=0
     do igam=1,shp%nb_gam(kd)
       if(shp%prob_gam(kd,igam) >= shp%dataset%params%prob_seuil_gam)then
          ngamkept=ngamkept+1
          shp%prob_gam(kd,ngamkept)=shp%prob_gam(kd,igam)
          do llk=1,shp%dataset%params%longhap
            shp%gamete(1,llk,ngamkept)=shp%gamete(1,llk,igam)
            shp%gamete(2,llk,ngamkept)=shp%gamete(2,llk,igam)
          end do !i
       end if
     end do !igam
!
! standardisation
!
       probtot=0.d0
       do igam=1,ngamkept
        probtot=probtot+shp%prob_gam(kd,igam)
       end do
       do igam=1,ngamkept
        shp%prob_gam(kd,igam)=shp%prob_gam(kd,igam)/probtot
       end do

     shp%nb_gam(kd)=ngamkept

       ! correctiopn BUG : leak memory 05/11/2010
       deallocate (base)
       deallocate (origine)
       deallocate (tabdup)
       end subroutine local_gamete
!!***

!!****f* m_qtlmap_haplotype_ldla/point_gamete
!! NAME
!!    point_gamete
!! DESCRIPTION
!!
!! NOTES
!!
!! point_gamete determine le nom des alleles (allele) au marqueur lk du chromosome c
!!  recus par les descendant kd du pere ip et de la mere jm  affectee de la phase geno
!!  ainsi que leur origine grand parentale (orig)
!! allele(1) est l'allele recu du pere, allele(2) de la mere
!! allele(1)=allele(2)=0 si le duo ou le triplet est heterozygte identique
!! allele(1)=allele(2)=nmanque si le phenotype du descednant et ou du pere sont manquants
!! orig(1) vaut 1 si origine = grand pere paternel ; 2 = grand mere paternelle;
!!           0 = indetermine (duo ou triplet heterozygote identique); 3= pere homozygote
!! orig(2) est le symetrique pour l'origne maternelle
!!
!! Utilise nmanque,correp,pheno,corred,correm		de m_qtlmap_data
!! Utilise genotyp,genotypm				de m_qtlmap_genotype_data
!! SOURCE
       subroutine point_gamete(shp,ip,jm,geno,kd,c,lk,allele,orig)
           class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

     integer(kind=KIND_PHENO) :: allele(2),unknown

     integer :: orig(2)
     integer , intent(in) :: ip,jm,geno,kd,c,lk
     logical found
     type(GENEALOGY_BASE) , pointer    :: dg
     type(GENOTYPE_BASE) , pointer     :: dga
     type(GENEALOGY_RACE) , pointer    :: dgr

     dgr => shp%dataset%geneaRace
     dga => shp%dataset%genoAnimal
     dg => shp%dataset%genea

        allele(1)=dga%nmanque
        allele(2)=dga%nmanque
        unknown=0
        found=.false.
!  le genotype du pere et le phenotype du descendant ne doivent pas etre manquants
!
         if(shp%spt%genotyp(c,lk,dga%correp(ip),1) /= dga%nmanque .and. &
          dga%pheno(c,lk,dga%corred(kd),1) /= dga%nmanque) then  !1
!
!  cas du pere homozygote
!
         if(shp%spt%genotyp(c,lk,dga%correp(ip),1) == shp%spt%genotyp(c,lk,dga%correp(ip),2)) then  !2

           allele(1)=shp%spt%genotyp(c,lk,dga%correp(ip),1)
           if(dga%pheno(c,lk,dga%corred(kd),1) == allele(1)) then
             allele(2)=dga%pheno(c,lk,dga%corred(kd),2)
           else
             allele(2)=dga%pheno(c,lk,dga%corred(kd),1)
           end if
           found=.true.
!
! pere heterozygote

         else
!    descendant homozygote
!
           if(dga%pheno(c,lk,dga%corred(kd),1) == dga%pheno(c,lk,dga%corred(kd),2)) then !3
             allele(1)=dga%pheno(c,lk,dga%corred(kd),1)
             allele(2)=dga%pheno(c,lk,dga%corred(kd),2)
             found=.true.
           end if

! descendant heterozygote
!  son premier allele ne pouvant pas venir du pere
!
           if(dga%pheno(c,lk,dga%corred(kd),1) /=   shp%spt%genotyp(c,lk,dga%correp(ip),1) .and. &
                  dga%pheno(c,lk,dga%corred(kd),1) /=   shp%spt%genotyp(c,lk,dga%correp(ip),2)) then    !4
              allele(1)=dga%pheno(c,lk,dga%corred(kd),2)
              allele(2)=dga%pheno(c,lk,dga%corred(kd),1)
              found=.true.
           end if
!
! descendant heterozygote
!  son deuxieme allele ne pouvant pas venir du pere
!
           if(dga%pheno(c,lk,dga%corred(kd),2) /=   shp%spt%genotyp(c,lk,dga%correp(ip),1) .and. &
                  dga%pheno(c,lk,dga%corred(kd),2) /=   shp%spt%genotyp(c,lk,dga%correp(ip),2)) then      !5
               allele(1)=dga%pheno(c,lk,dga%corred(kd),1)
               allele(2)=dga%pheno(c,lk,dga%corred(kd),2)
               found=.true.
           end if
!
! incertitude pere - descendant, la mere etant connue
!
           if (.not. found ) then   !6
            allele(1)=unknown
            allele(2)=unknown
            if ( dga%pheno(c,lk,dga%correm(jm),1) /= dga%nmanque ) then   !6
!
!  cas de la mere homozygote
!
             if(dga%pheno(c,lk,dga%correm(jm),1) == dga%pheno(c,lk,dga%correm(jm),2)) then  !7

               allele(2)=dga%pheno(c,lk,dga%correm(jm),1)
               if(dga%pheno(c,lk,dga%corred(kd),1) == allele(2)) then
                allele(1)=dga%pheno(c,lk,dga%corred(kd),2)
               else
                 allele(1)=dga%pheno(c,lk,dga%corred(kd),1)
               end if
               found=.true.
             end if
!
! mere heterozygote
! descendant heterozygote
!  son premier allele ne pouvant pas venir de la mere
!
             if(dga%pheno(c,lk,dga%corred(kd),1) /=   dga%pheno(c,lk,dga%correm(jm),1) .and. &
                  dga%pheno(c,lk,dga%corred(kd),1) /=   dga%pheno(c,lk,dga%correm(jm),2)) then  !8
              allele(1)=dga%pheno(c,lk,dga%corred(kd),1)
              allele(2)=dga%pheno(c,lk,dga%corred(kd),2)
              found=.true.
             end if
!
! descendant heterozygote
!  son deuxieme allele ne pouvant pas venir de la mere
!
             if(dga%pheno(c,lk,dga%corred(kd),2) /=   dga%pheno(c,lk,dga%correm(jm),1) .and. &
                  dga%pheno(c,lk,dga%corred(kd),2) /=   dga%pheno(c,lk,dga%correm(jm),2)) then  !9
               allele(1)=dga%pheno(c,lk,dga%corred(kd),2)
               allele(2)=dga%pheno(c,lk,dga%corred(kd),1)
               found=.true.
             end if

         end if
         end if
         end if
         end if

!
!  on peut alors determiner les origines grands parentales des gametes recus par kd
!
!  0 = indetermine (le descendant est heterozygote comme ses parents)
!  3 = indeterminable (les parents sont homozygotes )
!  4 =  donnees manquantes)
!
         orig(1)=0
         if(shp%spt%genotyp(c,lk,dga%correp(ip),1) == dga%nmanque .or. dga%pheno(c,lk,dga%corred(kd),1) == dga%nmanque) orig(1) =4
         if(orig(1) == 0 .and. (shp%spt%genotyp(c,lk,dga%correp(ip),1) == shp%spt%genotyp(c,lk,dga%correp(ip),2))) orig(1)=3
         if(orig(1) == 0 .and. allele(1) /= unknown) then
             if(allele(1) == shp%spt%genotyp(c,lk,dga%correp(ip),1)) then
               orig(1)=1
             else
               orig(1)=2
             end if
         end if

         orig(2)=0
         if(geno /=  0) then
!
! cas ou la mere est phasee (orig correspond aux origines grand parentales)
!
          if(shp%spt%genotypm(c,lk,geno,1) == dga%nmanque .or. dga%pheno(c,lk,dga%corred(kd),1) == dga%nmanque) orig(2)= 4
          if(orig(2) == 0 .and. (shp%spt%genotypm(c,lk,geno,1) == shp%spt%genotypm(c,lk,geno,2))) orig(2)=3
          if(orig(2) == 0 .and. allele(2) /= unknown) then
                if(allele(2) == shp%spt%genotypm(c,lk,geno,1)) then
                  orig(2)=1
                else
                  orig(2)=2
                end if
           end if
         else
!
! cas ou la mere n'est pas phasee (orig correspond a l'ordre de lecture des phenotypes)
!
          if(dga%pheno(c,lk,dga%correm(jm),1) == dga%nmanque .or. dga%pheno(c,lk,dga%corred(kd),1) == dga%nmanque) orig(2)= 4
          if(orig(2) == 0 .and. (dga%pheno(c,lk,dga%correm(jm),1) == dga%pheno(c,lk,dga%correm(jm),2))) orig(2)=3
          if(orig(2) == 0 .and. allele(2) /= unknown) then
                if(allele(2) == dga%pheno(c,lk,dga%correm(jm),1)) then
                 orig(2)=1
                  else
                orig(2)=2
                end if
              end if
          end if
      end subroutine point_gamete
!!***

!!****f* m_qtlmap_haplotype_ldla/liste_gamete
!! NAME
!!    liste_gamete
!! DESCRIPTION
!!
!! NOTES
!!
!! liste_gamete trouve les gametes du chromosome c ayant forme
!! chacun des descendant
!!
!! Utilise nmk,np,posi,correp,nm,repfem,correm,estfem    de m_qtlmap_data
!! Utilise genotyp,ngenom,genotypm                       de m_qtlmap_haplotype_data
!! SOURCE
       subroutine liste_gamete(shp,c)
         class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp
!
      integer       , intent(in) :: c

      integer ip,jm,kd,geno,i,lk!,indtab(10000)
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga

      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea

      geno=0
!
        do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            if (dga%presentg(c,kd)) call shp%global_gamete(ip,jm,geno,kd,c)

!            do i=1,2
!              do lk=1,map%nmk(c)
!                indtab(lk)=0
!                if(shp%tab_IBS(kd,lk,i) /= 0)indtab(lk)=shp%tab_IBS(kd,lk,i)- VAL_MIN_INDEX_PHENO
!              end do ! lk
     !       end do ! i
          end do !kd
        end do !jm
      end do !ip
!

         end subroutine liste_gamete
!!***

!!****f* m_qtlmap_haplotype_ldla/global_gamete
!! NAME
!!    global_gamete
!! DESCRIPTION
!!
!! NOTES
!!  global_gamete trouve les gametes parentaux transmis par ses parents (ip et jm)
!!  pour une phase geno de la mere
!!  a un individu particulier (kd)
!!  sur l'ensemble du chromosome (c)
!!  les doubles recombinaisons entre marqueurs informatifs flanquants sont supposees absentes
!!
!!  en sortie, global_gamete fournit
!!	shp%tab_IBS(kd,lk,i) pour lk=1,dataset%map%nmk(c) et i=1 (gamete paternel) et 2 gamete maternel)
!!
!!
!!  Utilise nmk,nmanque,correp,corred,correm,pheno     de m_qtlmap_data
!!  Utilise genotyp                                    de m_qtlmap_haplotype_data
!!
!! HYSTORY
!!   correction BUG : leak memory 05/11/2010
!! SOURCE
         subroutine global_gamete(shp,ip,jm,geno,kd,c)
           class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

      integer , intent(in) ::ip,jm,geno,kd,c

! divers
      integer lk,llk
      integer (kind=KIND_PHENO) orig_left,orig_right
      logical found
      integer (kind=KIND_PHENO) :: allele(2)
      integer :: orig(2)
      integer  :: stat
      integer (kind=KIND_PHENO) :: unknown=0
! dimension allocatable
      integer(kind=KIND_PHENO)  , dimension (:)  , allocatable   :: origine
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_RACE) , pointer    :: dgr
      type(MAP_BASE) , pointer    :: map

      dgr => shp%dataset%geneaRace
      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea
      map => shp%dataset%map

! allocation
!
      allocate (origine(maxval(map%nmk)))

       shp%tab_IBS(kd,:,:)=0

       do lk = 1,map%nmk(c)
!
!  pour chaque marqueur de l'intervalle, on cherche (avec point_gamete) le nom des alleles recus
!  de chacun des parents (allele), ainsi que leur origine grand parentale
!  orig vaut 1 si origine = grand pere; 2 = grand mere ;
!            0 = indetermine (duo ou triplet heterozygote identique); 3= parent homozygote, 4= manquant
!
         call shp%point_gamete(ip,jm,geno,kd,c,lk,allele,orig)
         origine(lk)=orig(1)
         if(origine(lk) ==1 .or. origine(lk) ==2 .or. origine(lk) == 3) shp%tab_IBS(kd,lk,1)=allele(1)
!
! si l'allele paternel est indermine on regarde ce qui s'est passe a gauche et a droite
         if(origine(lk) == 0) then

             orig_left=0
             found=.false.
             llk=lk
             do while((.not. found ).and. llk > 1)
               llk=llk-1
               if(origine(llk) ==1 .or. origine(llk) ==2 )  then
                  found=.true.
                  orig_left=origine(llk)
               end if
             end do ! found
             orig_right=0
             found=.false.
             llk=lk
             do while((.not. found ) .and. llk < map%nmk(c))
               llk=llk+1
               call shp%point_gamete(ip,jm,geno,kd,c,llk,allele,orig)
               if(orig(1) ==1 .or. orig(1) ==2 )  then
                  found=.true.
                  orig_right=orig(1)
               end if
             end do ! found

             if(orig_left == orig_right) shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),orig_left)
             !if(orig_left == 0)          shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),orig_right)    !!ATTENTION orig_right peut valoir =0
             !if(orig_right == 0)         shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),orig_left)
             if( ((orig_left == 0).or.(orig_left == 4) ) .and. ( (orig_right /= 0) .and. (orig_right /= 4)) )  then
               if ( orig_right /= 3 ) then
                 shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),orig_right)
               else
                 shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),1)
               end if
             end if

             if( ((orig_right == 0).or.(orig_right == 4) ) .and. ( (orig_left /= 0) .and. (orig_left /= 4)) )  then
               if ( orig_left /= 3 ) then
                 shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),orig_left)
               else
                 shp%tab_IBS(kd,lk,1)=shp%spt%genotyp(c,lk,dga%correp(ip),1)
               end if
             end if

             !!*** pas de calcul de proba....a voir ***
             !! si orig_right et orig_left inconnu (0 ou 4) alors shp%tab_IBS pas initialise...

         end if
!
!  connaissant l'allele recu du pre on en deduit celui recu de la mere
!
!         if(dga%pheno(c,lk,dga%corred(kd),1) == dga%pheno(c,lk,dga%corred(kd),2)) then
!           shp%tab_IBS(kd,lk,2)=dga%pheno(c,lk,dga%corred(kd),1)
!         else
!           if(dga%pheno(c,lk,dga%corred(kd),1) == shp%tab_IBS(kd,lk,1)) shp%tab_IBS(kd,lk,2)=dga%pheno(c,lk,dga%corred(kd),2)
!           if(dga%pheno(c,lk,dga%corred(kd),2) == shp%tab_IBS(kd,lk,1)) shp%tab_IBS(kd,lk,2)=dga%pheno(c,lk,dga%corred(kd),1)
!         end if

           if(dga%pheno(c,lk,dga%corred(kd),1) == shp%tab_IBS(kd,lk,1)) then
             shp%tab_IBS(kd,lk,2)=dga%pheno(c,lk,dga%corred(kd),2)
           else
             shp%tab_IBS(kd,lk,2)=dga%pheno(c,lk,dga%corred(kd),1)
           end if

       end do !lk

       ! correctiopn BUG : leak memory 05/11/2010
       deallocate (origine)
       end subroutine global_gamete
!!***

!!****f* m_qtlmap_haplotype_ldla/liste_shared_haplo
!! NAME
!!    liste_shared_haplo
!! DESCRIPTION
!!
!! NOTES
!!
!! liste_shared_haplo  determine si un descendant a recu de sa mere , autour de la positon dx,
!! un allele ibs a un de ceux de son pere, et sur quelle longueur
!! l'effet haplotypique sera ajoute au modle de description des performances si cette longueur
!! est superieure a LONG_MIN_IBS
!!
!! Utilise nmk,np,posi,correp,nm                         de m_qtlmap_data
!! Utilise genotyp                                       de m_qtlmap_haplotype_data
!! SOURCE
       subroutine liste_shared_haplo(shp,c,dx,n)
         class(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)         :: shp

!
      integer       , intent(in) :: c
      real(kind=dp) , intent(in) :: dx
      integer, intent (in) :: n

      integer lkleft,ip,jm,kd,jp,iall,icount,icount_max,lk,stat
      type(GENEALOGY_BASE) , pointer :: dg
      type(GENOTYPE_BASE) , pointer :: dga
      type(MAP_BASE) , pointer    :: map

      map => shp%dataset%map
      dga => shp%dataset%genoAnimal
      dg => shp%dataset%genea

      if (associated(shp%shared_haplo)) deallocate(shp%shared_haplo)
      allocate (shp%shared_haplo(dg%nd))


!  la position testee doit etre entre les bornes du chromsome
!
      if((dx < map%posi(c,1)) .or. (dx > map%posi(c,map%nmk(c)))) then
        call stop_application("liste_shared_haplo : incompatibilty between haplotype length and tested position.")
      end if

      shp%tab_shared_haplo=0
!
!  recherche du marqueur a gauche de dx
      do lkleft=1,map%nmk(c)-1
        if(dx <= map%posi(c,lkleft)) exit
      end do
!
      do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
            icount_max=0
            shp%shared_haplo(kd)=0
!
! icount est la longueur du segment de l'haplotype maternel ibs au chromosome iall du pere jp
            do jp=1,dg%np
              do iall=1,2
                icount=0
                do lk=lkleft,1,-1
                  if(shp%tab_IBS(kd,lk,2) /= shp%spt%genotyp(c,lk,dga%correp(jp),iall)) exit
                  if(shp%spt%genotyp(c,lk,dga%correp(jp),iall)/= dga%nmanque) icount=icount+1
                end do
                do lk=lkleft+1,map%nmk(c)
                  if(shp%tab_IBS(kd,lk,2) /= shp%spt%genotyp(c,lk,dga%correp(jp),iall)) exit
                  if(shp%spt%genotyp(c,lk,dga%correp(jp),iall)/= dga%nmanque) icount=icount+1
                end do

                if(icount > icount_max .and. icount >= shp%dataset%params%LONG_MIN_IBS) then
                  icount_max=icount
                  shp%shared_haplo(kd)=iall+2*(jp-1)
                end if
              end do ! iall
            end do ! jp

            if( shp%shared_haplo(kd) /= 0) then
              iall=2-modulo(shp%shared_haplo(kd),2)
              jp=(shp%shared_haplo(kd)+modulo(shp%shared_haplo(kd),2))/2
              shp%tab_shared_haplo(jp,iall,1+icount_max)=shp%tab_shared_haplo(jp,iall,1+icount_max)+1 !nombre d'haplotyped e longueur count_max
              shp%tab_shared_haplo(jp,iall,1)=shp%tab_shared_haplo(jp,iall,1)+1 ! nombre d haplotype partage
            end if
          end do !kd
        end do !jm
      end do !ip
!
      end subroutine liste_shared_haplo
!!***

      end module m_qtlmap_haplotype_ldla
