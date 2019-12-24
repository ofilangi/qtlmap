module m_qtlmap_type_incidence_solution
    use m_qtlmap_constant
    use m_qtlmap_type_dataset

    implicit none

!!****t* QTLMAP_TYPES/DESC_EFFECT_TYPE
!!  NAME
!!   DESC_EFFECT_TYPE
!!  DESCRIPTION
!!    Description of an effect enumered in the INCIDENCE_TYPE
!!
!!    attributes members :
!!    --------------------
!!    name        : name of the effect
!!    start       : level start of the effect in the incidence matrix
!!    end         : level end of the effect in the incidence matrix
!!    haveSubDesc : true if the effect have contains sub effect :
!!                  sample :'Fixed effect A' is a effect that contains several level (value 1, value 2,...)
!!    listSubDesc : recursive structure that contains the sub effect
!!  NOTES
!!    This type is used by all analysis to make homogene the print of a estimation parameter solution
!! SOURCE
    type DESC_EFFECT_TYPE
        ! Name of effect
        character(len=LEN_W)                                   :: name
        ! start ntniv
        integer                                                :: start   = -1
        ! end ntniv
        integer                                                :: end     = -1
        !If sub description :
        logical                                                :: haveSubDesc = .false.
        ! sub effect name
        type(DESC_EFFECT_TYPE),dimension(:),pointer            :: listSubDesc => NULL()
        ! true if the effect depends from the position in the linkage group
        logical                                                :: isVar = .false.

        contains

        procedure, public :: release       => release_description_type
        procedure, public :: copy          => copy_description_type

    end type DESC_EFFECT_TYPE
!!***
!!****t* QTLMAP_TYPES/TYPE_INCIDENCE_SOLUTION
!!  NAME
!!   TYPE_INCIDENCE_SOLUTION
!!  DESCRIPTION
!!    Description of a solution (estimation) of a analysis. All analysis fill this structure and give return the results.
!!
!!    attributes members :
!!    --------------------
!!    hypothesis       : Test hypothesis of the analysis
!!    sig              : solution of the standart deviation of each sire
!!    eqtl_print       : logical array that referenced the effect to print in the context of a eqtl printing
!!                       eqtl solution are printing in a array, we give only the following information :
!!                       st.dev, mean, qtl effect. All estimation of fixed effect and covariate are skipped
!!    haplotype        : Haplotypes
!!    groupeName       : All effect are categorized : General Mean, Polygenic effect, Fixed effect, Covariate, Qtl Effect
!!    nbParameterGroup : Number of parameter inside a group. sample, inside Polygenic effect : Sire A, Sire B,...
!!    parameterName    : Name of parameters by groupname
!!    paramaterValue   : Value of parameters by groupname
!!    parameterVecsol  : Estimability of parameters by groupname
!!    parameterPrecis  : Precision value of parameters by groupname
!!    qtl_groupeName   : Referenced the group of the QTL effect estimation
!!    rhoi             : residual correlation of traits
!! HISTORY
!!  09/09/2010    * add a dimension to qtl_groupeName (for detected the qtl effect with a specific trait on multi-trait analysis)
!! SOURCE
    type TYPE_INCIDENCE_SOLUTION
      ! Test based on the number of Qtls
      integer                                      :: hypothesis
      ! formated to print
      real (kind=dp)     ,dimension(:,:),pointer   :: sig => NULL()
      ! animal inbreeded s.d
      real (kind=dp)     ,dimension(:)  ,pointer   :: siga => NULL()
      ! Standart deviation for unknown haplotype dam
      real (kind=dp)     ,dimension(:,:),pointer   :: unknown_dam_sig => NULL()
      ! LD case => haplotype by sire , np, hypothesis,2
      character(len=LEN_W) ,dimension(:,:,:),pointer :: haplotypes       => NULL()
      character(len=LEN_W) ,dimension(:,:,:),pointer :: races_haplotypes => NULL()
      !
      logical            ,dimension(:) ,pointer    :: eqtl_print => NULL()
      !
      character(len=LEN_W) ,dimension(:),pointer   :: groupeName => NULL()
      !
      integer              ,dimension(:),pointer   :: nbParameterGroup => NULL()
      ! List of parameter name
      character(len=LEN_W) ,dimension(:,:),pointer :: parameterName => NULL()
      ! List of value
      real (kind=dp)     ,dimension(:,:),pointer   :: paramaterValue => NULL()

      logical        , dimension(:,:) , pointer    :: parameterVecsol => NULL()

      real (kind=dp) ,dimension(:,:),pointer       :: parameterPrecis => NULL()

      ! Index postion in groupeName of qtl effect : ncar,nqtl
      integer       ,dimension(:,:),pointer        :: qtl_groupeName => NULL()

      real (kind=dp)     ,dimension(:,:),pointer   :: rhoi => NULL()

      contains

      procedure, public :: release       => release_incidence_solution
      procedure, public :: debug_solution


    end type TYPE_INCIDENCE_SOLUTION
!!***


    contains
!!****f* m_qtlmap_tools/release_incidence_solution
!!  NAME
!!    release_incidence_solution
!!  DESCRIPTION
!!   release a variable of type TYPE_INCIDENCE_SOLUTION
!!  INPUTS/OUTPUTS
!!    obj   : the variable to release
!! SOURCE
      subroutine release_incidence_solution(obj)
        class(TYPE_INCIDENCE_SOLUTION) , intent(inout) :: obj

        if ( associated (obj%sig) ) deallocate (obj%sig)
        if ( associated (obj%unknown_dam_sig) ) deallocate (obj%unknown_dam_sig)
        if ( associated (obj%haplotypes) ) deallocate (obj%haplotypes)
        if ( associated (obj%races_haplotypes) ) deallocate (obj%races_haplotypes)
        if ( associated (obj%paramaterValue) ) deallocate (obj%paramaterValue)
        if ( associated (obj%parameterVecsol) ) deallocate (obj%parameterVecsol)
        if ( associated (obj%eqtl_print) ) deallocate (obj%eqtl_print)
        if ( associated (obj%parameterPrecis) ) deallocate (obj%parameterPrecis)
        if ( associated(obj%rhoi) ) deallocate (obj%rhoi)

        if ( associated (obj%groupeName) ) deallocate (obj%groupeName)
        if ( associated (obj%nbParameterGroup) ) deallocate (obj%nbParameterGroup)
        if ( associated (obj%parameterName) ) deallocate (obj%parameterName)
        if (associated(obj%qtl_groupeName)) deallocate (obj%qtl_groupeName)
        if (associated(obj%siga)) deallocate (obj%siga)

      end subroutine release_incidence_solution
!!***

!!****f* m_qtlmap_tools/copy_description_type
!!  NAME
!!    copy_description_type
!!  DESCRIPTION
!!   copy the variable in_desc of type TYPE_INCIDENCE_SOLUTION to out_desc
!!  INPUTS
!!    in_desc  : the variable to copy
!!    out_desc : the copy
!! SOURCE
      recursive subroutine copy_description_type(in_desc,out_desc)
       class(DESC_EFFECT_TYPE) , intent(in) :: in_desc
       type(DESC_EFFECT_TYPE) , intent(out) :: out_desc
       integer :: i

        out_desc%name =  in_desc%name
        out_desc%isVar =  in_desc%isVar
        out_desc%start =  in_desc%start
        out_desc%end =  in_desc%end
        out_desc%haveSubDesc =  in_desc%haveSubDesc

        if ( in_desc%haveSubDesc ) then
          allocate ( out_desc%listSubDesc(size(in_desc%listSubDesc)) )
          do i=1,size(in_desc%listSubDesc)
            call copy_description_type(in_desc%listSubDesc(i),out_desc%listSubDesc(i))
          end do
        end if

    end subroutine copy_description_type
!!***

!!****f* m_qtlmap_tools/release_description_type
!!  NAME
!!    release_description_type
!!  DESCRIPTION
!!   release a variable of type DESC_EFFECT_TYPE
!!  INPUTS/OUTPUTS
!!    inout_desc  : the variable to release
!! SOURCE
    recursive subroutine release_description_type(inout_desc)
      class(DESC_EFFECT_TYPE) , intent(inout) :: inout_desc
      integer :: i

        if ( inout_desc%haveSubDesc ) then
          do i=1,size(inout_desc%listSubDesc)
            call release_description_type(inout_desc%listSubDesc(i))
          end do
          deallocate (inout_desc%listSubDesc)
        end if

    end subroutine release_description_type
!!***

!!****f* m_qtlmap_tools/debug_solution
!! NAME
!!    debug_solution
!! DESCRIPTION
!!
!! HISTORY
!!
!! SOURCE
    subroutine debug_solution(incsol,dataset,nq,listChr,listN)
         class(TYPE_INCIDENCE_SOLUTION)     , intent(in)   :: incsol
         type(QTLMAP_DATASET)          , intent(in)       :: dataset
         integer                           , intent(in)   :: nq
         integer  ,dimension(nq)           , intent(in)   :: listChr,listN

         integer :: i,g,ip,chr,ic,l,iq
         real(kind=dp) :: wq
         type(DATAMODEL_BASE) , pointer :: dpm

         dpm => dataset%phenoModel

         write(*,*) "---------------------------------------------------------------"
         write(*,*) 'Estimation of parameters under H'//trim(str(incsol%hypothesis))
         do i=1,nq
           write(*,*) "POSITION ",i,":",dataset%map%absi(listChr(i),listN(i))
         end do
         write(*,*) "---------------------------------------------------------------"
         write(*,*)
         write(*,*) 'Within sire standard deviation'
         do ic=1,size(incsol%sig,1)
           if ( size(incsol%sig,1) > 1 ) then
             write(*,*) " ** Trait ",trim(dpm%carac(ic))," **"
           end if
          do ip = 1,dataset%genea%np
           write(*,fmt="(' sire ',a, '  s.d. :',f10.3)") trim(dataset%genea%pere(ip)),incsol%sig(ic,ip)
          end do
         end do

         write(*,*)

         if ( associated(incsol%siga)) then
            write(*,*) " SIGA = ",incsol%siga(:)," **"
         end if

         write(*,"(//,'  parameter    ','        estimable ?    value     ','precision'/)")

         do i=1,size(incsol%groupeName)
           write (*,*) incsol%groupeName(i)
           write (*,*)

           if ( incsol%nbParameterGroup(i) == 1 ) then
             if ( incsol%parameterVecsol(i,1) ) then
               write (*,fmt="(a50,'  yes ',2f10.3,1x)") " " , incsol%paramaterValue(i,1),&
                        incsol%parameterPrecis(i,1)
             else
               write (*,fmt="(a50,'  no ')") incsol%groupeName(i)
             end if
           else

              do g=1,incsol%nbParameterGroup(i)
                if ( incsol%parameterVecsol(i,g) ) then
                  write (*,fmt="(a50,'  yes ',2f10.3,1x)") incsol%parameterName(i,g)&
                   , incsol%paramaterValue(i,g) ,  incsol%parameterPrecis(i,g)
                else
                  write (*,fmt="(a50,'  no ')") incsol%parameterName(i,g)
                end if
              end do
           end if
         end do

    end subroutine debug_solution


end module m_qtlmap_type_incidence_solution
