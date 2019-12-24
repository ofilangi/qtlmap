!!****m* OUTPUT/m_qtlmap_format_carthagene
!!  NAME
!!    m_qtlmap_format_carthagene
!!  SYNOPSIS
!!    Manage allele information and print the information on a ASCII carthagne format
!!  DESCRIPTION
!!    Print information phases for each animal in a encoded carthagene format according the following array :
!!
!!   Notation Synonym      Possible Genotypes
!!
!!   1        A                 F0 |M0
!!   2                          F0 |M1
!!   3                       F0 |M0 , F0 |M1
!!   4                           F1 |M0
!!   5                       F0 |M0 , F1 |M0
!!   6        H              F0 |M1 , F1 |M0
!!   7        D          F0 |M0 , F0 |M1 , F1 |M0
!!   8        B                   F1 |M1
!!   9                       F0 |M0 , F1 |M1
!!   a                       F0 |M1 , F1 |M1
!!   b                  F0 |M0 , F0 |M1 , F1 |M1
!!   c                       F1 |M0 , F1 |M1
!!   d                  F0 |M0 , F1 |M0 , F1 |M1
!!   e                 F0 |M1 , F1 |M0 , F1 |M1
!!   f              F0 |M0 , F0 |M1 , F1 |M0 , F1 |M1
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
module m_qtlmap_format_carthagene
    use m_qtlmap_log
    use m_qtlmap_base
    use m_qtlmap_types

    implicit none

    public :: print_transmission_allele

    contains

!!****f* m_qtlmap_format_carthagene/print_transmission_allele
!!  NAME
!!    print_transmission_allele
!!  DESCRIPTION
!!    use the prot array (probabilities to receive a haplotype from parents) and encode in the carthagene format.
!!
!!  NOTES
!!     prot ( *,*,*,1) --> F0/M0 => 1
!!     prot ( *,*,*,2) --> F0/M1 => 2
!!     prot ( *,*,*,3) --> F1/M0 => 4
!!     prot ( *,*,*,4) --> F1/M1 => 8
!!  SOURCE
     subroutine print_transmission_allele(dataset,prot,outfile,type_out)
       type(QTLMAP_DATASET)           ,intent(in) :: dataset
       character(len=LENGTH_MAX_FILE) ,intent(in) :: outfile
       real (kind=dp)   ,dimension(:,:,:,:), intent(in) :: prot
       integer                        ,intent(in) :: type_out
       integer :: val_recl = 2**14
       integer :: ios,unitf,valhexa(dataset%genea%nd),ll,kd,i,chr,stat
       character(len=1) ,dimension(0:15) :: hexa,synonym

       !si on tombe sur 0, cas impossible, on considere que toutes les transmissions sont possibles
       data(hexa(i),i=0,15) /'f','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'/
    data(synonym(i),i=0,15) /' ','A',' ',' ',' ',' ','H','D','B',' ',' ',' ',' ',' ',' ','-'/

       if (trim(outfile) == '' ) then
         call stop_application("Dev error : print_transmission_allele need a input file")
       end if

       do chr=1,dataset%map%nchr
           open(UNIT=unitf,file=trim(outfile)//trim(str(chr)), form="formatted",recl=val_recl,iostat=ios)
           if (ios/=0) then
            call stop_application("Can not open the file :"//trim(outfile)//trim(str(chr)))
           end if
           if ( type_out == 1 ) then
            write (unitf,fmt='(a)') "data type f2 intercross"
           else
            write (unitf,fmt='(a)') "data type f2 backcross"
           end if

           write (unitf,fmt='(i5,i5," 0 0 ")') dataset%genea%nd,dataset%map%nmk(chr)

           do ll=1,dataset%map%nmk(chr)
              valhexa = 0
              do kd=1,dataset%genea%nd

                 if ( prot (chr,ll,kd,1) /= 0 ) valhexa(kd) = valhexa(kd) + 1
                 if ( prot (chr,ll,kd,2) /= 0 ) valhexa(kd) = valhexa(kd) + 2
                 if ( prot (chr,ll,kd,3) /= 0 ) valhexa(kd) = valhexa(kd) + 4
                 if ( prot (chr,ll,kd,4) /= 0 ) valhexa(kd) = valhexa(kd) + 8
              end do
              select case (type_out)
              case (1)
                write (unitf,fmt='(a," ",'//trim(str(dataset%genea%nd))//'(a1))') '*'//&
                trim(dataset%map%mark(chr,ll)),&
                 ( hexa(valhexa(kd)),kd=1,dataset%genea%nd )
              case (2)
                write (unitf,fmt='(a," ",'//trim(str(dataset%genea%nd))//'(a1))') '*'//&
                trim(dataset%map%mark(chr,ll)),&
                 ( synonym(valhexa(kd)),kd=1,dataset%genea%nd )
              case default

              end select
           end do
           close(unitf)
           call log_mess(" ** Generate ["//trim(outfile)//trim(str(chr))//"] ** ",INFO_DEF)
       end do

     end subroutine
!!***


end module m_qtlmap_format_carthagene
