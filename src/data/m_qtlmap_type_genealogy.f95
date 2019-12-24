module m_qtlmap_type_genealogy
    use m_qtlmap_base
    use m_qtlmap_constant

    implicit none

!!****t* QTLMAP_TYPES/GENEALOGY_BASE
!!  NAME
!!     GENEALOGY_BASE
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type , public :: GENEALOGY_BASE

    integer                                        :: ngp=0
    integer                                        :: ngm=0
    integer , dimension (:), pointer               :: ngmgp => null()
    integer , dimension (:), pointer               :: nrgm => null()
    character(len=LEN_DEF) , dimension (:),pointer :: gmere => null()
    character(len=LEN_DEF) , dimension (:),pointer :: gpere => null()
    integer                                        :: nr=0

    character(len=LEN_DEF) , dimension (:),pointer :: repro => null()
    character(len=LEN_DEF) , dimension (:),pointer :: reprop => null()
    character(len=LEN_DEF) , dimension (:),pointer :: reprom => null()
    integer , dimension (:), pointer               :: rep_reprop => null()
    integer , dimension (:), pointer               :: rep_reprom => null()

    integer , dimension (:), pointer               :: ndm => null()
    integer , dimension (:), pointer               :: nmp => null()
    character(len=LEN_DEF) , dimension (:),pointer :: mere => null()
    character(len=LEN_DEF) , dimension (:),pointer :: pere => null()
    character(len=LEN_DEF) , dimension (:),pointer :: animal => null()

    integer                                        :: nd=0
    integer                                        :: nm=0
    integer                                        :: np=0
    integer                                        :: nfem = 0
    integer                , dimension (:),pointer :: reppere => null()
    integer                , dimension (:),pointer :: repmere => null()
    character(len=LEN_DEF) , dimension (:),pointer :: femelle => null()
    integer                , dimension (:),pointer :: repfem => null()
    character(len=LEN_DEF),dimension(:,:),pointer  :: OldGenealogy => null()
    integer                                        :: OldGenealogySize = 0

    contains
     procedure, public :: copy    => copy_genealogy_base
     procedure, public :: release => release_genealogy_base
     procedure ,public :: print   => print_genealogy_base

   end type GENEALOGY_BASE
!!***

!!****t* QTLMAP_TYPES/GENEALOGY_RACE
!!  NAME
!!     GENEALOGY_RACE
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type , public :: GENEALOGY_RACE
    character(len=LEN_DEF) , dimension (:),pointer , public :: racep => null()
    character(len=LEN_DEF) , dimension (:),pointer , public :: racem => null()
    character(len=LEN_DEF) , dimension (:),pointer , public :: nom_race => null()
    integer                                        , public :: NB_RACES=1

    contains

     procedure ,public :: copy    => copy_genealogy_race
     procedure ,public :: release  => release_genealogy_race

   end type GENEALOGY_RACE
!!***

   contains


       !!****f* m_qtlmap_types/release_genealogy_base
    !!  NAME
    !!    release_genealogy_base
    !!  DESCRIPTION
    !!    release memory store in the GENEALOGY_BASE
    !!  NOTES
    !!
    !!  SOURCE
    subroutine release_genealogy_base(dg)
        class(GENEALOGY_BASE), intent(inout)  :: dg

        if (associated(dg%ngmgp)) deallocate (dg%ngmgp)
        if (associated(dg%nrgm)) deallocate (dg%nrgm)
        if (associated(dg%gmere)) deallocate (dg%gmere)
        if (associated(dg%gpere)) deallocate (dg%gpere)
        if (associated(dg%reprop)) deallocate (dg%reprop)
        if (associated(dg%reprom)) deallocate (dg%reprom)
        if (associated(dg%repro)) deallocate (dg%repro)
        if (associated(dg%rep_reprop)) deallocate (dg%rep_reprop)
        if (associated(dg%rep_reprom)) deallocate (dg%rep_reprom)
        if (associated(dg%ndm)) deallocate (dg%ndm)
        if (associated(dg%nmp)) deallocate (dg%nmp)
        if (associated(dg%mere)) deallocate (dg%mere)
        if (associated(dg%pere)) deallocate (dg%pere)
        if (associated(dg%animal)) deallocate (dg%animal)
        if (associated(dg%reppere)) deallocate (dg%reppere)
        if (associated(dg%repmere)) deallocate (dg%repmere)
        if (associated(dg%femelle)) deallocate (dg%femelle)
        if (associated(dg%repfem)) deallocate (dg%repfem)
        if (associated(dg%OldGenealogy)) deallocate (dg%OldGenealogy)

    end subroutine release_genealogy_base
    !!***
    !!****f* m_qtlmap_types/delete_genealogy
    !!  NAME
    !!    delete_genealogy
    !!  DESCRIPTION
    !!    release memory store in the GENEALOGY_BASE
    !!  NOTES
    !!
    !!  SOURCE
    subroutine release_genealogy_race(dgr)
        class(GENEALOGY_RACE), intent(inout)  :: dgr
        if (associated(dgr%racep)) deallocate (dgr%racep)
        if (associated(dgr%racem)) deallocate (dgr%racem)
        if (associated(dgr%nom_race)) deallocate (dgr%nom_race)
    end subroutine release_genealogy_race
    !!***






    subroutine copy_genealogy_race(mapGeneaRace,copyMapGeneaRace)
       class(GENEALOGY_RACE), intent(in)  :: mapGeneaRace
       type(GENEALOGY_RACE), intent(inout)  :: copyMapGeneaRace

       if ( associated(mapGeneaRace%racep)) then
        allocate(copyMapGeneaRace%racep(size(mapGeneaRace%racep)))
        copyMapGeneaRace%racep = mapGeneaRace%racep
       end if

       if ( associated(mapGeneaRace%racem)) then
        allocate(copyMapGeneaRace%racem(size(mapGeneaRace%racem)))
        copyMapGeneaRace%racem = mapGeneaRace%racem
       end if

       if ( associated(mapGeneaRace%nom_race)) then
        allocate(copyMapGeneaRace%nom_race(size(mapGeneaRace%nom_race)))
        copyMapGeneaRace%nom_race = mapGeneaRace%nom_race
       end if

       copyMapGeneaRace%NB_RACES=mapGeneaRace%NB_RACES

    end subroutine copy_genealogy_race

    !!****f* m_qtlmap_types/copy_genealogy
    !!  NAME
    !!    copy_genealogy
    !!  DESCRIPTION
    !!    copy all structure in the type GENEALOGY_BASE
    !!  NOTES
    !!
    !!  SOURCE
    subroutine copy_genealogy_base(genea,geneaCopy)
    class(GENEALOGY_BASE) , intent(in)  :: genea
    type(GENEALOGY_BASE) , intent(inout) :: geneaCopy

    geneaCopy%ngp=genea%ngp
    geneaCopy%ngm=genea%ngm
    allocate (geneaCopy%ngmgp(size(genea%ngmgp)))
    geneaCopy%ngmgp = genea%ngmgp
    allocate (geneaCopy%nrgm(size(genea%nrgm)))
    geneaCopy%nrgm = genea%nrgm
    allocate (geneaCopy%gmere(size(genea%gmere)))
    geneaCopy%gmere = genea%gmere
    allocate (geneaCopy%gpere(size(genea%gpere)))
    geneaCopy%gpere = genea%gpere
    geneaCopy%nr=genea%nr

    allocate (geneaCopy%repro(size(genea%repro)))
    geneaCopy%repro = genea%repro
    allocate (geneaCopy%reprop(size(genea%reprop)))
    geneaCopy%reprop = genea%reprop
    allocate (geneaCopy%reprom(size(genea%reprom)))
    geneaCopy%reprom = genea%reprom
    allocate (geneaCopy%rep_reprop(size(genea%rep_reprop)))
    geneaCopy%rep_reprop = genea%rep_reprop
    allocate (geneaCopy%rep_reprom(size(genea%rep_reprom)))
    geneaCopy%rep_reprom = genea%rep_reprom

    allocate (geneaCopy%ndm(size(genea%ndm)))
    geneaCopy%ndm = genea%ndm
    allocate (geneaCopy%nmp(size(genea%nmp)))
    geneaCopy%nmp = genea%nmp
    allocate (geneaCopy%mere(size(genea%mere)))
    geneaCopy%mere = genea%mere
    allocate (geneaCopy%pere(size(genea%pere)))
    geneaCopy%pere = genea%pere
    allocate (geneaCopy%animal(size(genea%animal)))
    geneaCopy%animal = genea%animal

    geneaCopy%nd=genea%nd
    geneaCopy%nm=genea%nm
    geneaCopy%np=genea%np
    geneaCopy%nfem=genea%nfem

    allocate (geneaCopy%reppere(size(genea%reppere)))
    geneaCopy%reppere = genea%reppere
    allocate (geneaCopy%repmere(size(genea%repmere)))
    geneaCopy%repmere = genea%repmere
    allocate (geneaCopy%femelle(size(genea%femelle)))
    geneaCopy%femelle = genea%femelle
    allocate (geneaCopy%repfem(size(genea%repfem)))
    geneaCopy%repfem = genea%repfem
    allocate (geneaCopy%OldGenealogy(size(genea%OldGenealogy,1),size(genea%OldGenealogy,2)))
    geneaCopy%OldGenealogy = genea%OldGenealogy

    geneaCopy%OldGenealogySize=genea%OldGenealogySize

   end subroutine copy_genealogy_base
!!***

   subroutine print_genealogy_base(genea,unit)
    class(GENEALOGY_BASE) , intent(in)  :: genea
    integer             ,intent(in)   :: unit ! 6 stdout

    write (unit=6,fmt=*) "======================="
    write (unit=6,fmt=*) "NGP            :",genea%ngp
    write (unit=6,fmt=*) "NGM            :",genea%ngm
    write (unit=6,fmt=*) "NR             :",genea%nr
    write (unit=6,fmt=*) "size(ngmgp)              :",size(genea%ngmgp)
    write (unit=6,fmt=*) "size(nrgm)               :",size(genea%nrgm)
    write (unit=6,fmt=*) "size(gmere)              :",size(genea%gmere)
    write (unit=6,fmt=*) "size(gpere)              :",size(genea%gpere)
    write (unit=6,fmt=*)
    write (unit=6,fmt=*) "size(repro)              :",size(genea%repro)
    write (unit=6,fmt=*) "size(reprop)             :",size(genea%reprop)
    write (unit=6,fmt=*) "size(reprom)             :",size(genea%reprom)
    write (unit=6,fmt=*) "size(rep_reprop)         :",size(genea%rep_reprop)
    write (unit=6,fmt=*) "size(rep_reprom)         :",size(genea%rep_reprom)
    write (unit=6,fmt=*)
    write (unit=6,fmt=*) "NP                       :",genea%np
    write (unit=6,fmt=*) "NM                       :",genea%nm
    write (unit=6,fmt=*) "ND                       :",genea%nd
    write (unit=6,fmt=*) "NFEM                     :",genea%nfem
    write (unit=6,fmt=*) "size(ndm)                :",size(genea%ndm)
    write (unit=6,fmt=*) "size(nmp)                :",size(genea%nmp)
    write (unit=6,fmt=*) "size(mere)               :",size(genea%mere)
    write (unit=6,fmt=*) "size(pere)               :",size(genea%pere)
    write (unit=6,fmt=*) "size(animal)             :",size(genea%animal)
    write (unit=6,fmt=*)
    write (unit=6,fmt=*) "size(reppere)            :",size(genea%reppere)
    write (unit=6,fmt=*) "size(repmere)            :",size(genea%repmere)
    write (unit=6,fmt=*) "size(femelle)            :",size(genea%femelle)
    write (unit=6,fmt=*) "size(repfem)             :",size(genea%repfem)
    write (unit=6,fmt=*)
    write (unit=6,fmt=*) "OldGenealogySize         :",genea%OldGenealogySize
    write (unit=6,fmt=*) "size(OldGenealogy)       :",size(genea%OldGenealogy,1),size(genea%OldGenealogy,2)

   end subroutine print_genealogy_base

  !!*******************************************************************************************************************

end module m_qtlmap_type_genealogy
