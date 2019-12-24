!!****m* HAPLOTYPE/m_qtlmap_haplotype
!!  NAME
!!    m_qtlmap_haplotype
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
module m_qtlmap_haplotype
   use m_qtlmap_types
   use m_qtlmap_base
   use m_qtlmap_log
   use m_qtlmap_haplotype_V1
   use m_qtlmap_haplotype_V2
   use m_qtlmap_phase_offspring
   use m_qtlmap_haplotype_external

   implicit none

   public :: haplotype
   public :: check_genotype
   public :: get_information_informative_marker
   public :: get_informativity_position_sire

   contains
!!****f* m_qtlmap_haplotype/haplotype
!!  NAME
!!    m_qtlmap_haplotype
!!  DESCRIPTION
!!   Apply the algorithm identify by opt_version
!!    - Initialisation of the following structures : pdd, genotyp, genotypm, ngenom, ngend, ndesc, probg, phasp, phasm
!!  INPUTS
!!    opt_version : the version of the haplotype to used
!!  SOURCE
      subroutine haplotype(dataset,spt,opt_version)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
        integer                    ,intent(in)            :: opt_version

        !free structure if allocated
        call spt%release()

        select case (opt_version)

        case (VERSION_HAPLOTYPE_V1)
         call log_mess('Computation of transmission probabilities V1',INFO_DEF)
         call haplotype_V1(dataset,spt)

        case (VERSION_HAPLOTYPE_V2)
          call log_mess('Computation of transmission probabilities V2',INFO_DEF)
          call haplotype_V2(dataset,spt)

        case (VERSION_HAPLOTYPE_V3)
            call log_mess('Computation of transmission probabilities V3',INFO_DEF)
            call haplotype_V3(dataset,spt)

        case (VERSION_HAPLOTYPE_SNP)
          call log_mess('Computation of transmission probabilities SNP',INFO_DEF)
          call haplotype_SNP(dataset,spt)

        case (VERSION_HAPLOTYPE_SYMMAX2SAT_SNP)
          call log_mess('Phases SYMMAX2SAT + Computation of transmission probabilities SNP',INFO_DEF)
          call haplotype_SYMMAX2SAT_SNP(dataset,spt)

        case (VERSION_HAPLOTYPE_PARENTAL_EXTERNAL)
          call log_mess('Parents phases are given by the user + Computation of transmission probabilities',INFO_DEF)
          call haplotype_external(dataset,spt)


        case default
         call stop_application('bad value of opt_version['//trim(str(opt_version))//'].')
       end select

       call log_mess("** END module haplotype ** ",DEBUG_DEF)
      end subroutine haplotype
!!***


!!****f* m_qtlmap_haplotype/check_genotype
!!  NAME
!!    check_genotype
!!  DESCRIPTION
!!
!!  SOURCE
      subroutine check_genotype(dataset,spt)
         type(QTLMAP_DATASET)       ,intent(in)            :: dataset
         type(PDD_BUILD)            ,intent(inout)         :: spt
         call check_genotype_mod(dataset,spt)

      end subroutine check_genotype
!!***


      subroutine get_information_informative_marker(dataset,spt,countH,markH)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
        real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr), intent(out)               :: countH
        real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)),intent(out)    :: markH

        call internal_get_information_informative_marker(dataset,spt,countH,markH)

      end subroutine get_information_informative_marker

!!****f* m_qtlmap_haplotype/haplotype_offspring
!!  NAME
!!    haplotype_offspring
!!  DESCRIPTION
!!   Apply the algorithm identify by opt_haplotype_offspring
!!    - Initialisation of the following structures : genotyp,reconstructed
!!  INPUTS
!!    opt_version : the version of the haplotype to used
!!  SOURCE
      subroutine haplotype_offspring(dataset,spt,opt_version)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
        integer                         , intent(in)                           :: opt_version

        !free structure if allocated
        if ( associated(spt%reconstructed) ) deallocate( spt%reconstructed )
        call haplotype_offspring_v1(dataset,spt)
         

      end subroutine haplotype_offspring
!!***


!!****f* m_qtlmap_haplotype/get_pdd_at_mark
!!  NAME
!!    get_pdd_at_mark
!!  DESCRIPTION
!!    Get the index position of the array pdd (prob of transmission) at the position : marker lk on the chromosome ch.
!!  INPUTS
!!    ch : index chromosome
!!    lk : index marker
!!  OUTPUTS
!!    pos : index position
!!  SOURCE
      subroutine get_pdd_at_mark(ch,lk,pos)
          integer , intent(in)  :: ch
          integer , intent(in)  :: lk
          integer , intent(out) :: pos

          real(kind=dp) ,parameter :: epsil = 0.001

          integer :: istart,iend



      end subroutine get_pdd_at_mark
!!***

!!****f* m_qtlmap_haplotype/test_print_haplotype
!!  NAME
!!    test_print_haplotype
!!  DESCRIPTION
!!
!!  SOURCE
     subroutine test_print_haplotype

!      real (kind=dp)   ,dimension(:,:,:,:),allocatable :: prot
!      allocate (prot(nchr,maxval(nmk),nd,4),STAT=stat)
!
!      call gammapf(prot)
!
!      do ip=1,np
!	   do jm=nmp(ip)+1,nmp(ip+1)
!        do geno=ngenom(jm)+1,ngenom(jm+1)
!          do kd=ngend(geno)+1,ngend(geno+1)
!             kkd=ndesc(kd)
!             !! animal s ALL_M1 ALL_M2 ....  origine : 1/2
!          end do
!        end do
!       end do
!      end do
!
!      deallocate (prot)

     end subroutine test_print_haplotype

     ! Get informativity at the position pos (GL chr).
     ! let probp = pdd(1)+ pdd(2)  , the probability to receive the chr1 of sire for a progeny
     ! inf = 2 * Sum | probp - 0.5 | / nd_ip
     ! 0 < inf < 1
     ! Evolution #1579
     function get_informativity_position_sire(dataset,spt,ip,chr,pos) result(inf)
       type(QTLMAP_DATASET)       ,intent(in)         :: dataset
       type(PDD_BUILD)            ,intent(in)         :: spt
       integer                    ,intent(in)         :: ip  ! index sire
       integer                    ,intent(in)         :: chr ! chromosome
       integer                    ,intent(in)         :: pos ! position

       integer :: jm,kd,ig

       real(kind=dp) :: inf,subinf
       integer       :: ndk

       inf = 0.d0
       ndk = 0

       do jm=dataset%genea%nmp(ip)+1,dataset%genea%nmp(ip+1)
        ig = spt%ngenom(chr,jm)+1
        ndk = ndk + spt%ngend(chr,ig+1) - spt%ngend(chr,ig)
        do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
          subinf=0.d0
          do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
           subinf = subinf + abs( spt%pdd(chr,kd,1,pos)+spt%pdd(chr,kd,2,pos) - 0.5d0)
          end do
          inf = inf + spt%probg(chr,ig)*subinf
        end do
       end do

       inf = ( 2 * inf ) / real(ndk)

     end function get_informativity_position_sire

end module m_qtlmap_haplotype
