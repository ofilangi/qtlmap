module m_qtlmap_type_pdd
    use m_qtlmap_constant
    use m_qtlmap_type_dataset
    implicit none

!Description
!https://forge-dga.jouy.inra.fr/projects/qtlmap/wiki/HaploOrga
   type ,public :: PDD_BUILD
       real (kind=dp)   ,dimension(:,:,:,:),pointer   ,public :: pdd    => null()
       integer          ,dimension(:,:)    ,pointer   ,public :: ngenom => null()
       integer          ,dimension(:,:)    ,pointer   ,public :: ngend  => null()
       integer          ,dimension(:,:)    ,pointer   ,public :: ndesc  => null()
       real (kind=dp)   ,dimension(:,:)    ,pointer   ,public :: probg  => null()
       logical          ,dimension(:,:)    ,pointer   ,public :: phasp  => null()
       logical          ,dimension(:,:)    ,pointer   ,public :: phasm  => null()
       real (kind=dp)   ,dimension(:,:,:,:),pointer   ,public :: prot   => null()
       integer(kind=KIND_PHENO),dimension(:,:,:,:),pointer   ,public :: genotyp => null()
       integer(kind=KIND_PHENO),dimension(:,:,:,:),pointer   ,public :: genotypm => null()
       logical  , dimension (:,:,:), pointer    ,public     :: reconstructed => null()

     contains

      procedure ,public :: copy    => copy_pdd_build
      procedure ,public :: release => release_pdd_build
      procedure ,public :: get_maxnbgenotypedam

   end type PDD_BUILD

  contains

      subroutine copy_pdd_build(spt,copySpt)
      class(PDD_BUILD)            ,intent(in)    :: spt
      type(PDD_BUILD)             ,intent(inout) :: copySpt


    end subroutine copy_pdd_build


   subroutine release_pdd_build(spt)
     class(PDD_BUILD)            ,intent(inout) :: spt

     if (associated(spt%pdd)) deallocate (spt%pdd)
     if (associated(spt%ngenom)) deallocate (spt%ngenom)
     if (associated(spt%ngend)) deallocate (spt%ngend)
     if (associated(spt%ngend)) deallocate (spt%ngend)
     if (associated(spt%probg)) deallocate (spt%probg)
     if (associated(spt%phasp)) deallocate (spt%phasp)
     if (associated(spt%phasm)) deallocate (spt%phasm)
     if (associated(spt%prot)) deallocate (spt%prot)
     if (associated(spt%genotyp)) deallocate (spt%genotyp)
     if (associated(spt%genotypm)) deallocate (spt%genotypm)
     if (associated(spt%reconstructed)) deallocate (spt%reconstructed)

   end subroutine release_pdd_build


!!****f* m_qtlmap_tools/get_maxnbgenotypedam
!!  NAME
!!    get_maxnbgenotypedam
!!  DESCRIPTION
!!   get the maximum number of genotype (all dams confused) found by the haplotype on the genome wide
!!  RESULTS
!!    maxng  : the variable to release
!! SOURCE
      function get_maxnbgenotypedam(prob_trans,dataset) result(maxng)
         class(PDD_BUILD)     , intent(in)       :: prob_trans
         type(QTLMAP_DATASET) , intent(in)       :: dataset

         integer             :: maxng

         maxng = maxval(prob_trans%ngenom(:,dataset%genea%nmp(dataset%genea%np+1)+1))

      end function get_maxnbgenotypedam
!!***


end module m_qtlmap_type_pdd
