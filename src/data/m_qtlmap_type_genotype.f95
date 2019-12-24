module m_qtlmap_type_genotype
    use m_qtlmap_base
    use m_qtlmap_constant
    use m_qtlmap_type_map
    implicit none

    !!****t* QTLMAP_TYPES/GENOTYPE_BASE
!!  NAME
!!     GENOTYPE_BASE
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
   type ,public :: GENOTYPE_BASE
!!****v* m_qtlmap_data/nall
!!  NAME
!!    nall
!!  DESCRIPTION
!!   Contains the number of allele by chromosome and marker
!!  DIMENSIONS
!!   1 : chromosome index
!!   2 : marker index
!!
!!***
   integer            , dimension (:,:),allocatable        ,public :: nall
   ! dim nd : true if genotype for animal (genealogy) is defined

   !============== VERIF ALLELES NECESSAIRE ====================
!!****v* m_qtlmap_data/alleles
!!  NAME
!!    alleles
!!  DESCRIPTION
!!    Get the string value of the allele referenced.
!!  DIMENSIONS
!!   1: chromosome index
!!   2: marker index
!!   3: allele index ( 1 < allele index < nall(ch,m) )
!!***
   character(len=LEN_DEF) , dimension (:,:,:),allocatable    ,public :: alleles
!!****v* m_qtlmap_data/pc_all
!!  NAME
!!    pc_all
!!  DESCRIPTION
!!    Frequencies alleles.
!!  DIMENSIONS
!!   1: chromosome index (ch)
!!   2: marker index     (m)
!!   3: allele index ( 1 < allele index < nall(ch,m) )
!!***
   real (kind=dp), dimension (:,:,:), allocatable          ,public :: pc_all
     integer                                                 :: nmes = 0
     character(len=LEN_DEF)  , dimension (:),pointer         :: numero => null()
     integer(kind=KIND_PHENO),dimension (:,:,:,:),pointer    :: pheno => null()
     character(len=LEN_DEF) , dimension (:,:) ,pointer       :: value_pheno => null()
     integer                , dimension (:)    ,pointer      :: nb_value_pheno => null()
     integer            , dimension (:),pointer              :: corregp => null()
     integer            , dimension (:),pointer              :: corregm => null()
     integer            , dimension (:),pointer              :: correr => null()
     integer            , dimension (:),pointer              :: correm => null()
     integer            , dimension (:),pointer              :: correp => null()
     integer            , dimension (:),pointer              :: corred => null()
     logical            , dimension (:,:), pointer           :: presentg => null()
     logical            , dimension (:),pointer              :: estfem => null() ! Estimabilité / nb de descendant genotypé    ( 1ere dim : indice female)
     ! Genotype code for missing value
     integer(kind=KIND_PHENO)             ,      public      :: NMANQUE

     contains

      procedure ,public :: copy    => copy_genotype_base
      procedure ,public :: release => release_genotype_base
      procedure ,public :: print   => print_genotype_base
      procedure ,public :: init_pheno
      procedure ,public :: end_pheno
      procedure ,public :: set_pheno
      procedure ,public :: get_pheno
      procedure ,public :: get_ind_pheno
      !procedure ,public :: write_file => write_typ
   end type GENOTYPE_BASE
!!***

    contains

!!****f* m_qtlmap_types/delete_genotype
!!  NAME
!!    delete_genotype
!!  DESCRIPTION
!!    release memory store in the GENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
    subroutine release_genotype_base(mapGenotype)
      class(GENOTYPE_BASE), intent(inout)  :: mapGenotype

     deallocate(mapGenotype%numero)
     deallocate(mapGenotype%pheno)
     deallocate(mapGenotype%value_pheno)
     deallocate(mapGenotype%nb_value_pheno)
     deallocate(mapGenotype%corregp)
     deallocate(mapGenotype%corregm)
     deallocate(mapGenotype%correr)
     deallocate(mapGenotype%correm)
     deallocate(mapGenotype%correp)
     deallocate(mapGenotype%corred)
     deallocate(mapGenotype%presentg)
     deallocate(mapGenotype%nall)
     deallocate(mapGenotype%alleles)
     deallocate(mapGenotype%pc_all)

    end subroutine release_genotype_base
!!***

    subroutine print_genotype_base(mapGenotype,unit)
      class(GENOTYPE_BASE), intent(in)  :: mapGenotype
      integer             ,intent(in)   :: unit ! 6 stdout

      write (unit=6,fmt=*) "======================="
      write (unit=6,fmt=*) "NMES            :",mapGenotype%nmes
      write (unit=6,fmt=*) "SIZE(NUMERO)    :",size(mapGenotype%numero)
      write (unit=6,fmt=*) "NUMERO(1)       :",mapGenotype%numero(1)
      write (unit=6,fmt=*) "NUMERO(",size(mapGenotype%numero),")       :",&
        mapGenotype%numero(size(mapGenotype%numero))
      write (unit=6,fmt=*) "SIZE(PHENO)     :",size(mapGenotype%pheno,1),size(mapGenotype%pheno,2),&
      size(mapGenotype%pheno,3),size(mapGenotype%pheno,4)
      write (unit=6,fmt=*) "SIZE(value_pheno):",size(mapGenotype%value_pheno)
      write (unit=6,fmt=*) "SIZE(nb_value_pheno):",size(mapGenotype%nb_value_pheno)
      write (unit=6,fmt=*) "SIZE(corregp):",size(mapGenotype%corregp)
      write (unit=6,fmt=*) "SIZE(corregm):",size(mapGenotype%corregm)
      write (unit=6,fmt=*) "SIZE(correr):",size(mapGenotype%correr)
      write (unit=6,fmt=*) "SIZE(correm):",size(mapGenotype%correm)
      write (unit=6,fmt=*) "SIZE(correp):",size(mapGenotype%correp)
      write (unit=6,fmt=*) "SIZE(corred):",size(mapGenotype%corred)
      write (unit=6,fmt=*) "SIZE(presentg):",size(mapGenotype%presentg,1),size(mapGenotype%presentg,2)
      write (unit=6,fmt=*) "SIZE(nall):",size(mapGenotype%nall,1),size(mapGenotype%nall,2)
      write (unit=6,fmt=*) "SIZE(alleles):",size(mapGenotype%alleles,1),&
          size(mapGenotype%alleles,2),size(mapGenotype%alleles,3)
      write (unit=6,fmt=*) "SIZE(pc_all):",size(mapGenotype%pc_all,1),&
          size(mapGenotype%pc_all,2),size(mapGenotype%pc_all,3)

    end subroutine print_genotype_base

!!****f* m_qtlmap_types/copy_genotype_base
!!  NAME
!!    copy_genotype_base
!!  DESCRIPTION
!!    copy all structure in the type GENOTYPE_BASE
!!  NOTES
!!
!!  SOURCE
    subroutine copy_genotype_base(mapGenotype,copyMapGenotype)
        class(GENOTYPE_BASE), intent(in)     :: mapGenotype
        type(GENOTYPE_BASE), intent(inout)  :: copyMapGenotype

     copyMapGenotype%nmes =mapGenotype%nmes
     allocate (copyMapGenotype%numero(size(mapGenotype%numero)))
     copyMapGenotype%numero=mapGenotype%numero

     allocate (copyMapGenotype%pheno(size(mapGenotype%pheno,1),size(mapGenotype%pheno,2),&
      size(mapGenotype%pheno,3),size(mapGenotype%pheno,4)))
     copyMapGenotype%pheno=mapGenotype%pheno

     allocate (copyMapGenotype%value_pheno(size(mapGenotype%value_pheno,1),size(mapGenotype%value_pheno,2)))
     copyMapGenotype%value_pheno=mapGenotype%value_pheno

     allocate (copyMapGenotype%nb_value_pheno(size(mapGenotype%nb_value_pheno)))
     copyMapGenotype%nb_value_pheno=mapGenotype%nb_value_pheno

     allocate (copyMapGenotype%corregp(size(mapGenotype%corregp)))
     copyMapGenotype%corregp=mapGenotype%corregp

     allocate (copyMapGenotype%corregm(size(mapGenotype%corregm)))
     copyMapGenotype%corregm=mapGenotype%corregm

     allocate (copyMapGenotype%correr(size(mapGenotype%correr)))
     copyMapGenotype%correr=mapGenotype%correr

     allocate (copyMapGenotype%correm(size(mapGenotype%correm)))
     copyMapGenotype%correm=mapGenotype%correm

     allocate (copyMapGenotype%correp(size(mapGenotype%correp)))
     copyMapGenotype%correp=mapGenotype%correp

     allocate (copyMapGenotype%corred(size(mapGenotype%corred)))
     copyMapGenotype%corred=mapGenotype%corred

     allocate (copyMapGenotype%presentg(size(mapGenotype%presentg,1),size(mapGenotype%presentg,2)))
     copyMapGenotype%presentg=mapGenotype%presentg

     allocate (copyMapGenotype%nall(size(mapGenotype%nall,1),size(mapGenotype%nall,2)))
     copyMapGenotype%nall=mapGenotype%nall

     allocate (copyMapGenotype%alleles(size(mapGenotype%alleles,1),&
      size(mapGenotype%alleles,2),size(mapGenotype%alleles,3)))
     copyMapGenotype%alleles=mapGenotype%alleles

     allocate (copyMapGenotype%pc_all(size(mapGenotype%pc_all,1),&
      size(mapGenotype%pc_all,2),size(mapGenotype%pc_all,3)))
     copyMapGenotype%pc_all=mapGenotype%pc_all


    end subroutine copy_genotype_base
!!***

!!****f* m_qtlmap_tools/init_pheno
!!  NAME
!!    init_pheno
!!  DESCRIPTION
!!     allocate arrays to manage the internal encoded format of allele : value_pheno,nb_value_pheno
!!  INPUTS
!!    nch : number of chromosome to manage
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!! SOURCE
       subroutine init_pheno(dga,nch)
         class(GENOTYPE_BASE), intent(inout) :: dga
         integer ,intent(in) :: nch
         allocate (dga%value_pheno(nch,256))
         allocate (dga%nb_value_pheno(nch))
         dga%value_pheno=""
         dga%nb_value_pheno=0

       end subroutine init_pheno
!!***


       subroutine end_pheno(dga)
         class(GENOTYPE_BASE), intent(inout) :: dga
         if (associated (dga%value_pheno)) deallocate (dga%value_pheno)
         if (associated (dga%nb_value_pheno)) deallocate (dga%nb_value_pheno)

       end subroutine end_pheno

!!****f* m_qtlmap_tools/set_pheno
!!  NAME
!!    set_pheno
!!  DESCRIPTION
!!    Setting a new identifiant number allele with allele string value or get the value of the identifiant
!!    if this value exist.
!!    * Reallocation or Allocation (first call) of the value_pheno for each call
!!  INPUTS
!!    * value     -- Value of the allele string
!!  RESULT
!!    The idenfiant allele number
!!  EXAMPLE
!!
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!  m_qtlmap_data/value_pheno
!!  m_qtlmap_data/VAL_MIN_INDEX_PHENO
!!  m_qtlmap_data/VAL_MAX_INDEX_PHENO
!! SOURCE
       function set_pheno(dga,map,c,value) result(res)
           class(GENOTYPE_BASE), intent(inout) :: dga
           type (MAP_BASE)     ,intent(in)     :: map
           integer                 ,intent(in)  :: c
           character(len=LEN_DEF)  ,intent(in)  :: value
           integer(kind=KIND_PHENO)           :: res
           integer                            :: nbvalue
           integer                            :: i,ch
           character(len=LEN_DEF) , dimension (:,:)    ,allocatable  :: buf_value_pheno
           res = VAL_MIN_INDEX_PHENO

           do i=1,dga%nb_value_pheno(c)
              if ( dga%value_pheno(c,i) == trim(value) ) then
                  return
              end if
              res = res + 1
           end do

           !manage error
           if ( i > VAL_MAX_INDEX_PHENO ) then
             nbvalue = VAL_MAX_INDEX_PHENO
             nbvalue = nbvalue - VAL_MIN_INDEX_PHENO

             write (0,*) "* Too many allele detected for internal data structure "
             write (0,*) "number max [",nbvalue,']'
             do ch=1,map%nchr
             print *,"------chromo : "//trim(map%chromo(ch))//"---------- size: ",dga%nb_value_pheno(ch)
             print *,(trim(dga%value_pheno(ch,i))//",",i=1,dga%nb_value_pheno(ch))
             end do
             stop
           end if

           if ( size(dga%value_pheno,2) <= dga%nb_value_pheno(c) ) then
            allocate (buf_value_pheno(map%nchr,size(dga%value_pheno,2)))
            buf_value_pheno = dga%value_pheno
            deallocate(dga%value_pheno)
            allocate(dga%value_pheno(map%nchr,size(buf_value_pheno,2)+256) )
            dga%value_pheno=""
            do ch=1,map%nchr
             do i=1,dga%nb_value_pheno(ch)
              dga%value_pheno(ch,i)=trim(buf_value_pheno(ch,i))
             end do
            end do
            deallocate(buf_value_pheno)
           end if

           dga%nb_value_pheno(c)=dga%nb_value_pheno(c)+1
           dga%value_pheno(c,dga%nb_value_pheno(c))=trim(value)
           return

       end function set_pheno
!!***

!!****f* m_qtlmap_tools/get_pheno
!!  NAME
!!    get_pheno
!!  DESCRIPTION
!!    Get the corresponding identifiant number of a coded initial allele string
!!  INPUTS
!!    * value    -- Value of the allele string
!!  RESULT
!!    The string value of the allele
!!  EXAMPLE
!!
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!  m_qtlmap_data/value_pheno
!!  m_qtlmap_data/VAL_MIN_INDEX_PHENO
!! SOURCE
       function get_pheno(dga,c,value) result(res)
           class(GENOTYPE_BASE), intent(in) :: dga
           integer                   ,intent(in)  :: c
           integer(kind=KIND_PHENO)  ,intent(in)  :: value
           character(len=LEN_DEF)                 :: res
           integer                                :: ind,i

           ind = value - VAL_MIN_INDEX_PHENO + 1
           if ( ind < 1 .or. ind > dga%nb_value_pheno(c)) then
              write (0,*) "Devel error : get_pheno * none index correspond [",value,"]"
              print *, (trim(dga%value_pheno(c,i))//" ",i=1,dga%nb_value_pheno(c))
              stop
           end if
           res = dga%value_pheno(c,ind)
       end function get_pheno
!!***

!!****f* m_qtlmap_tools/get_ind_pheno
!!  NAME
!!    get_ind_pheno
!!  DESCRIPTION
!!    Get the corresponding internal identifiant number of a coded initial allele string
!!  INPUTS
!!    * c      : the chromosome
!!    * value  : Value of the string allele
!!  RESULT
!!    The identifiant value of the value allele
!!  EXAMPLE
!!
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!  m_qtlmap_data/value_pheno
!!  m_qtlmap_data/VAL_MIN_INDEX_PHENO
!! SOURCE
       function get_ind_pheno(dga,c,value) result(res)
           class(GENOTYPE_BASE), intent(in) :: dga
           integer                 ,intent(in)  :: c
           character(len=LEN_DEF)  ,intent(in)  :: value
           integer(kind=KIND_PHENO)             :: res
           integer :: i

            do i=1,dga%nb_value_pheno(c)
              if ( dga%value_pheno(c,i) == value ) then
                  res = VAL_MIN_INDEX_PHENO + i + 1
                  return
              end if
           end do
           write (0,*) "Devel error : get_ind_pheno * none value correspond [",trim(value),"]"
           stop
       end function get_ind_pheno
!!***


end module m_qtlmap_type_genotype
