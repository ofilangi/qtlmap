module m_qtlmap_type_phenotype
    use m_qtlmap_base
    use m_qtlmap_constant
    implicit none


!!****t* QTLMAP_TYPES/model_trait
!! NAME
!!   DESC_EFFECT_TYPE
!! DESCRIPTION
!!
!! NOTES
!!    This type is used by all analysis to make homogene the print of a estimation parameter solution
!! SOURCE
  type , public :: model_trait
     ! number of fixed effect
     integer                                      :: nbfe
     ! number of covariate
     integer                                      :: nbco
     ! number of qtl with interaction
     integer                                      :: nbqtlinter
     ! number of interaction with the qtl for each qtl with interaction defined
     !size : nbqtlinter
     integer     ,dimension(:) ,pointer           :: nbint                           => null()
     ! size : nbfe : reference
     integer     ,dimension(:) ,pointer           :: indexFixedEffect                => null()
     ! size : nbco : reference
     integer     ,dimension(:) ,pointer           :: indexCovariate                  => null()
     !size : nbqtlinter,nbint
     integer     ,dimension(:,:) ,pointer         :: indexFixedEffectWithInteraction => null()

   contains

     procedure, public :: release => release_model_trait
     procedure, public :: copy => copy_model_trait

  end type model_trait
!!***

!!****t* QTLMAP_TYPES/PHENOTYPE_BASE
!!  NAME
!!     PHENOTYPE_BASE
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type ,public :: PHENOTYPE_BASE
     character(len=LEN_DEF) , dimension (:),pointer        :: bete     => null()
     real(kind=dp)    ,dimension(:,:),pointer              :: y        => null()
     real(kind=dp)    ,dimension(:,:),pointer              :: cd       => null()
     integer , dimension (:,:),pointer                     :: ndelta   => null()
     logical          ,dimension (:,:),pointer             :: presentc => null()
     integer          ,dimension (:,:),pointer             :: nivx     => null()
     real (kind=dp)   ,dimension (:,:),pointer             :: covar    => null()

     integer        ,dimension(:,:)    ,pointer            :: ydiscretord => null()
     character(len=LEN_DEF)  ,dimension(:,:)  ,pointer     :: ycategorial => null()
     integer , dimension (:),pointer                       :: corperf => null() !!Contains the corresponding index of y from an index animal
     logical , dimension (:,:),pointer                     :: estime => null()
     integer  ,    dimension (:),pointer                   :: nmumest => null() ! Number of dam with estimability
     integer  ,    dimension (:),pointer                   :: namest => null() ! Number of female with estimability
     integer              , dimension (:,:)      ,pointer  :: iam => null() !!   Indexe of the female

     contains

      procedure ,public :: copy    => copy_phenotype_base
      procedure ,public :: release => release_phenotype_base

   end type PHENOTYPE_BASE
!!***

   type ,public :: DATAMODEL_BASE
     integer , dimension (:,:),pointer                     :: modele   => null()
     integer                                               :: ncar     = 0
     integer                                               :: ncarcat  = 0
     character(len=LEN_DEF) , dimension (:),pointer        :: carac    => null()
     character(len=1)    ,  dimension (:),pointer          :: natureY  => null()

     character(len=LEN_DEF) ,dimension(:),pointer ,public  :: filter_car_name => null()
     integer  , dimension(:),pointer             ,public   :: filter_car => null()
     integer                                     ,public   :: n_filter_car = 0

     integer        ,dimension(:,:)    ,pointer            :: indicemod   => null() ! Correspond value from a original discrete value and the array ydiscretord
     integer          ,dimension(:)   ,pointer             :: nmod        => null() ! Number of classes corresponding with discrete value
     real(kind=dp)          ,dimension(:,:) ,pointer       :: seuil       => null() ! Thresholds on the underlying scale
     real(kind=dp)          ,dimension(:,:) ,pointer       :: prop        => null() ! Proportion of the discrete classes

     integer                                               :: ncov     = 0
     character(len=LEN_DEF) , dimension (:),pointer        :: namecov  => null()

     integer                                               :: nfix     = 0
     character(len=LEN_DEF) , dimension (:),pointer        :: namefix  => null()
     integer , dimension (:),  pointer                     :: nlev     => null()
     character(len=LEN_DEF) , dimension (:,:)  ,pointer    :: listelev => null()

     real (kind=dp), dimension (:)   , pointer             :: h2       => null()
     real (kind=dp), dimension (:,:) , pointer             :: RhoP     => null()
     real (kind=dp), dimension (:,:) , pointer             :: RhoG     => null()
     real (kind=dp), dimension (:)   , pointer             :: xmut     => null()
     real (kind=dp), dimension (:)   , pointer             :: sigt     => null()

     type(model_trait) , dimension(:),pointer             :: listModelTrait => null()


     contains

      procedure ,public :: copy    => copy_datamodel_base
      procedure ,public :: release => release_datamodel_base
      procedure ,public :: write_file  => write_datamodel_base

   end type DATAMODEL_BASE
!!***


   contains


!!****f* m_qtlmap_types/delete_traits
!!  NAME
!!    delete_traits
!!  DESCRIPTION
!!    release memory store in the PHENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
    subroutine release_phenotype_base(mapTraits)
        class(PHENOTYPE_BASE), intent(inout)  :: mapTraits
        if (associated(mapTraits%bete)) deallocate (mapTraits%bete)
        if (associated(mapTraits%y)) deallocate (mapTraits%y)
        if (associated(mapTraits%cd)) deallocate (mapTraits%cd)
        if (associated(mapTraits%presentc)) deallocate (mapTraits%presentc)
        if (associated(mapTraits%nivx)) deallocate (mapTraits%nivx)
        if (associated(mapTraits%covar)) deallocate (mapTraits%covar)
        if (associated(mapTraits%ndelta)) deallocate (mapTraits%ndelta)
        if (associated(mapTraits%ydiscretord)) deallocate (mapTraits%ydiscretord)
        if (associated(mapTraits%ycategorial)) deallocate (mapTraits%ycategorial)
        if (associated(mapTraits%corperf)) deallocate (mapTraits%corperf)
        if (associated(mapTraits%estime)) deallocate (mapTraits%estime)
        if (associated(mapTraits%nmumest)) deallocate (mapTraits%nmumest)
        if (associated(mapTraits%namest)) deallocate (mapTraits%namest)
        if (associated(mapTraits%iam)) deallocate (mapTraits%iam)

    end subroutine release_phenotype_base

!!****f* m_qtlmap_types/copy_traits
!!  NAME
!!    copy_traits
!!  DESCRIPTION
!!    copy all structure in the type PHENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
   subroutine copy_phenotype_base(mapTraits,copyMapTraits)
        class(PHENOTYPE_BASE), intent(in)     :: mapTraits
        type(PHENOTYPE_BASE), intent(inout)  :: copyMapTraits

        allocate (copyMapTraits%bete(size(mapTraits%bete)))
        copyMapTraits%bete     = mapTraits%bete

        allocate (copyMapTraits%y(size(mapTraits%y,1),size(mapTraits%y,2)))
        copyMapTraits%y        = mapTraits%y

        allocate (copyMapTraits%cd(size(mapTraits%cd,1),size(mapTraits%cd,2)))
        copyMapTraits%cd       = mapTraits%cd

        allocate (copyMapTraits%presentc(size(mapTraits%presentc,1),size(mapTraits%presentc,2)))
        copyMapTraits%presentc = mapTraits%presentc

        allocate (copyMapTraits%ndelta(size(mapTraits%ndelta,1),size(mapTraits%ndelta,2)))
        copyMapTraits%ndelta   = mapTraits%ndelta

        allocate (copyMapTraits%nivx(size(mapTraits%nivx,1),size(mapTraits%nivx,2)))
        copyMapTraits%nivx     = mapTraits%nivx

        allocate (copyMapTraits%covar(size(mapTraits%covar,1),size(mapTraits%covar,2)))
        copyMapTraits%covar    = mapTraits%covar

        allocate (copyMapTraits%ydiscretord(size(mapTraits%ydiscretord,1),size(mapTraits%ydiscretord,2)))
        copyMapTraits%ydiscretord  = mapTraits%ydiscretord

        allocate (copyMapTraits%ycategorial(size(mapTraits%ycategorial,1),size(mapTraits%ycategorial,2)))
        copyMapTraits%ycategorial  = mapTraits%ycategorial

        allocate (copyMapTraits%corperf(size(mapTraits%corperf)))
        copyMapTraits%corperf      = mapTraits%corperf

        allocate (copyMapTraits%estime(size(mapTraits%estime,1),size(mapTraits%estime,2)))
        copyMapTraits%estime  = mapTraits%estime

        allocate (copyMapTraits%nmumest(size(mapTraits%nmumest)))
        copyMapTraits%nmumest = mapTraits%nmumest

        allocate (copyMapTraits%namest(size(mapTraits%namest)))
        copyMapTraits%namest  = mapTraits%namest

        allocate (copyMapTraits%iam(size(mapTraits%iam,1),size(mapTraits%iam,2)))
        copyMapTraits%iam     = mapTraits%iam

    end subroutine copy_phenotype_base
!!***

!!****f* m_qtlmap_types/release_datamodel_base
!!  NAME
!!    delete_traits
!!  DESCRIPTION
!!    release memory store in the PHENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
    subroutine release_datamodel_base(mapModel)
        class(DATAMODEL_BASE), intent(inout)  :: mapModel

        integer :: i

        if (associated(mapModel%modele)) deallocate (mapModel%modele)
        if (associated(mapModel%carac)) deallocate (mapModel%carac)
        if (associated(mapModel%natureY)) deallocate (mapModel%natureY)
        if (associated(mapModel%filter_car)) deallocate (mapModel%filter_car)
        if (associated(mapModel%indicemod)) deallocate (mapModel%indicemod)
        if (associated(mapModel%nmod)) deallocate (mapModel%nmod)
        if (associated(mapModel%seuil)) deallocate (mapModel%seuil)
        if (associated(mapModel%prop)) deallocate (mapModel%prop)
        if (associated(mapModel%namecov)) deallocate (mapModel%namecov)
        if (associated(mapModel%namefix)) deallocate (mapModel%namefix)
        if (associated(mapModel%nlev)) deallocate (mapModel%nlev)
        if (associated(mapModel%listelev)) deallocate (mapModel%listelev)
        if (associated(mapModel%h2)) deallocate (mapModel%h2)
        if (associated(mapModel%RhoP)) deallocate (mapModel%RhoP)
        if (associated(mapModel%RhoP)) deallocate (mapModel%RhoP)
        if (associated(mapModel%xmut)) deallocate (mapModel%xmut)
        if (associated(mapModel%sigt)) deallocate (mapModel%sigt)
        if (associated(mapModel%listModelTrait)) then
         do i=1,size(mapModel%listModelTrait)
          call mapModel%listModelTrait(i)%release()
         end do
         deallocate (mapModel%listModelTrait)
        end if

    end subroutine release_datamodel_base

!!****f* m_qtlmap_types/copy_datamodel_base
!!  NAME
!!    copy_traits
!!  DESCRIPTION
!!    copy all structure in the type PHENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
   subroutine copy_datamodel_base(mapModel,copyMapModel)
        class(DATAMODEL_BASE), intent(in)     :: mapModel
        type(DATAMODEL_BASE), intent(inout)  :: copyMapModel

        integer :: i

        if (associated(mapModel%modele)) then
         allocate (copyMapModel%modele(size(mapModel%modele,1),size(mapModel%modele,2)))
         copyMapModel%modele = mapModel%modele
        end if

        copyMapModel%ncar = mapModel%ncar
        copyMapModel%ncarcat = mapModel%ncarcat

        if (associated(mapModel%carac)) then
         allocate (copyMapModel%carac(size(mapModel%carac)))
         copyMapModel%carac = mapModel%carac
        end if

        if (associated(mapModel%natureY)) then
         allocate (copyMapModel%natureY(size(mapModel%natureY)))
         copyMapModel%natureY = mapModel%natureY
        end if

        if (associated(mapModel%filter_car)) then
         allocate (copyMapModel%filter_car(size(mapModel%filter_car)))
         copyMapModel%filter_car = mapModel%filter_car
        end if

        if (associated(mapModel%indicemod)) then
         allocate (copyMapModel%indicemod(size(mapModel%indicemod,1),size(mapModel%indicemod,2)))
         copyMapModel%indicemod = mapModel%indicemod
        end if

        if (associated(mapModel%nmod)) then
         allocate (copyMapModel%nmod(size(mapModel%nmod)))
         copyMapModel%nmod = mapModel%nmod
        end if

        if (associated(mapModel%seuil)) then
         allocate (copyMapModel%seuil(size(mapModel%seuil,1),size(mapModel%seuil,2)))
         copyMapModel%seuil = mapModel%seuil
        end if

        if (associated(mapModel%prop)) then
         allocate (copyMapModel%prop(size(mapModel%prop,1),size(mapModel%prop,2)))
         copyMapModel%prop = mapModel%prop
        end if

        copyMapModel%ncov = mapModel%ncov

        if (associated(mapModel%namecov)) then
         allocate (copyMapModel%namecov(size(mapModel%namecov)))
         copyMapModel%namecov = mapModel%namecov
        end if

        copyMapModel%nfix = mapModel%nfix

        if (associated(mapModel%namefix)) then
         allocate (copyMapModel%namefix(size(mapModel%namefix)))
         copyMapModel%namefix = mapModel%namefix
        end if

        if (associated(mapModel%nlev)) then
         allocate (copyMapModel%nlev(size(mapModel%nlev)))
         copyMapModel%nlev = mapModel%nlev
        end if

        if (associated(mapModel%listelev)) then
         allocate (copyMapModel%listelev(size(mapModel%listelev,1),size(mapModel%listelev,2)))
         copyMapModel%listelev = mapModel%listelev
        end if

        if (associated(mapModel%h2)) then
         allocate (copyMapModel%h2(size(mapModel%h2)))
         copyMapModel%h2 = mapModel%h2
        end if

        if (associated(mapModel%RhoP)) then
         allocate (copyMapModel%RhoP(size(mapModel%RhoP,1),size(mapModel%RhoP,2)))
         copyMapModel%RhoP = mapModel%RhoP
        end if

        if (associated(mapModel%RhoG)) then
         allocate (copyMapModel%RhoG(size(mapModel%RhoG,1),size(mapModel%RhoG,2)))
         copyMapModel%RhoP = mapModel%RhoG
        end if

        if (associated(mapModel%xmut)) then
         allocate (copyMapModel%xmut(size(mapModel%xmut)))
         copyMapModel%xmut = mapModel%xmut
        end if

        if (associated(mapModel%xmut)) then
         allocate (copyMapModel%sigt(size(mapModel%sigt)))
         copyMapModel%sigt = mapModel%sigt
        end if

        if (associated(mapModel%listModelTrait)) then
         allocate (copyMapModel%listModelTrait(size(mapModel%listModelTrait)))

         do i=1,size(mapModel%listModelTrait)
          call mapModel%listModelTrait(i)%copy(copyMapModel%listModelTrait(i))
         end do
        end if

    end subroutine copy_datamodel_base
!!***

!!****f* m_qtlmap_tools/release_model_trait
!!  NAME
!!    release_model_trait
!!  DESCRIPTION
!!   release a variable of type model_trait
!!  INPUTS/OUTPUTS
!!    inout_desc  : the variable to release
!! SOURCE
    subroutine release_model_trait(modelTrait)
      class(model_trait) , intent(inout)  :: modelTrait

     if ( associated(modelTrait%nbint) ) deallocate(modelTrait%nbint)
     if ( associated(modelTrait%indexFixedEffect) ) deallocate(modelTrait%indexFixedEffect)
     if ( associated(modelTrait%indexCovariate) ) deallocate(modelTrait%indexCovariate)
     if ( associated(modelTrait%indexFixedEffectWithInteraction) ) deallocate(modelTrait%indexFixedEffectWithInteraction)

    end subroutine release_model_trait
!!***

    subroutine copy_model_trait(modelTrait,copy)
      class(model_trait) , intent(in)    :: modelTrait
      type(model_trait) , intent(inout)  :: copy

     copy%nbfe = modelTrait%nbfe
     copy%nbco = modelTrait%nbco
     copy%nbqtlinter = modelTrait%nbqtlinter

     allocate (copy%nbint(size(modelTrait%nbint)))
     copy%nbint                           = modelTrait%nbint
     allocate (copy%indexFixedEffect(size(modelTrait%indexFixedEffect)))
     copy%indexFixedEffect                = modelTrait%indexFixedEffect
     allocate (copy%indexCovariate(size(modelTrait%indexCovariate)))
     copy%indexCovariate                  = modelTrait%indexCovariate
     allocate (copy%indexFixedEffectWithInteraction(size(modelTrait%indexFixedEffectWithInteraction,1),&
                                                    size(modelTrait%indexFixedEffectWithInteraction,2)))
     copy%indexFixedEffectWithInteraction = modelTrait%indexFixedEffectWithInteraction
    end subroutine copy_model_trait



  subroutine write_datamodel_base(dpm,file)
    class(DATAMODEL_BASE)            ,intent(in) :: dpm
    character(len=LENGTH_MAX_FILE)   ,intent(in) :: file
    integer  :: u
    u=789

    open (u,file=file)

    write (u,*) dpm%ncar
    write (u,*) dpm%nfix,dpm%ncov
    write (u,*) dpm%namefix,dpm%namecov

    close (u)

  end subroutine write_datamodel_base

end module m_qtlmap_type_phenotype
