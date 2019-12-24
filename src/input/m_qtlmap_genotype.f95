!!****m* INPUT/m_qtlmap_genotype
!!  NAME
!!    m_qtlmap_genotype -- Genotype routines module
!!  SYNOPSIS

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
module m_qtlmap_genotype

  use m_qtlmap_base
  use m_qtlmap_types
  use m_qtlmap_types
  use m_qtlmap_log
  use m_qtlmap_output_handler

  implicit none
  save

!!****d* m_qtlmap_genotype/unit_genotype
!!  NAME
!!   unit_genotype
!!  DESCRIPTION
!!
!!***
  integer , parameter, private                              :: unit_genotype = 56

  type , public :: WORK_GENOTYPE
!!****v* m_qtlmap_genotype/index_position_marker
!!  NAME
!!   index_position_marker
!!  DESCRIPTION
!!   give the position (column index) according a marker index
!!  NOTES
!!   each i element contains the position in the file of for the marker mark(i)
!!***
    integer , dimension (:), allocatable,private              :: index_position_marker,index_chr
  !
!!****v* m_qtlmap_genotype/nmarker_in_genotype_file
!!  NAME
!!   nmarker_in_genotype_file
!!  DESCRIPTION
!!    number of marker in the genotype file
!!  NOTES
!!***
     integer                   ,private                        :: nmarker_in_genotype_file = 0

   contains
     procedure, private :: check_marker_name
     procedure, private :: set_pheno_structure
     procedure, public  :: release => release_work_genotype

  end type WORK_GENOTYPE

  public :: read_genotype
  public :: sim_typ
  public :: recup
  public :: write_typ
  public :: check_typage
  public :: set_estfem
  public :: check_HWE

  contains

!!****f* m_qtlmap_genotype/read_genotype
!!  NAME
!!    read_genotype
!!  DESCRIPTION
!!    Read the user genotypic file
!!  INPUTS
!!   ndmin   : minimum number of progenies to build full sib family
!!
!!  NOTES
!!    ndmin used by set_estfem
!!  SOURCE
     subroutine read_genotype(dataset)
       type(QTLMAP_DATASET)      ,intent(inout)   :: dataset
       integer                                    :: ios

       type(WORK_GENOTYPE) :: work

       call log_mess('START SUBROUTINE : read_genotype',DEBUG_DEF)

       open (unit_genotype,FILE=dataset%params%get_file_val(K_GENOTYPE),IOSTAT=ios,&
        recl=2**25,action="read",status="old")

       if ( ios /= 0 ) then
          call stop_application('Cannot open genotype file ['//trim(dataset%params%get_file_val(K_GENOTYPE))//']')
       endif

       call log_mess('reading genotype file...',INFO_DEF)
       call allocate_vector_(dataset)
       call work%check_marker_name(dataset)
       call internal_read_file(dataset,work)
       close(unit_genotype)
       call set_corresponding_vector(dataset)
       call set_allele_info_vector(dataset)
       call set_estfem(dataset)

       call work%release()

       call log_mess('END SUBROUTINE : read_genotype',INFO_DEF)

     end subroutine read_genotype

!!  DESCRIPTION
!!   read and check the header genotypic file. fill internal arrays index_chr,index_position_marker,
!!  INPUTS
!!   file   : genotypic file
     subroutine check_marker_name(work,dataset)
     class(WORK_GENOTYPE)      ,intent(inout)   :: work
     type(QTLMAP_DATASET)      ,intent(inout)   :: dataset
     character(len=LEN_DEF)                  :: token
     integer                                 :: ios,i,j,k,c,long,s1
     logical                                 :: is_ok
     integer, parameter                      :: LENBUF=LEN_LINE*10
     character(len=LENBUF) ,pointer          :: line_read

     character(len=LEN_DEF) , dimension (:),allocatable ::mark00
     character(len=1)                        :: OneCar
     integer                                 :: nmk_nearly

     call log_mess('SUBROUTINE : check_marker_name',DEBUG_DEF)

     rewind (unit_genotype)

     ios = 0
     allocate (line_read)
     line_read = ''

     do while ( trim(line_read) == '' .and. (ios == 0)  )
       call GET(unit_genotype,line_read,maxlen=LENBUF,IOSTAT=ios)
     end do

     if (ios /= 0) then
		  call stop_application('No marker name are found in the genotype file. [file:'&
		    //trim(dataset%params%get_file_val(K_GENOTYPE))//']')
     end if

     nmk_nearly=1
     do i=1,len(trim(line_read))
        if ( line_read(i:i) == " ") nmk_nearly = nmk_nearly +1
     end do

     allocate (mark00(nmk_nearly))

     line_read=trim(line_read)
     long = len(trim(line_read))
     work%nmarker_in_genotype_file=1
     is_ok=.false.
     s1=1
     do i=1,long
      if ( line_read(i:i) == ' ') then
         if (is_ok) then
            mark00(work%nmarker_in_genotype_file)=adjustl(line_read(s1:i))
            s1=i
            work%nmarker_in_genotype_file=work%nmarker_in_genotype_file+1
            is_ok=.false.
         end if
      else if ( .not. is_ok ) then
         is_ok=.true.
      end if
     end do
     mark00(work%nmarker_in_genotype_file)=adjustl(line_read(s1:long))
     !stop
     call log_mess('number of marker in the header genotype file:'// &
      trim(str(work%nmarker_in_genotype_file)) , VERBOSE_DEF)

     allocate (work%index_position_marker(work%nmarker_in_genotype_file))
     allocate (work%index_chr(work%nmarker_in_genotype_file))

     work%index_position_marker = INT_NOT_DEFINED
     work%index_chr             = 0

     do c=1,dataset%map%nchr
      do i=1,dataset%map%nmk(c)
        ! we search the marker associated
        j = 1
        do while ( j <= size(mark00) .and. dataset%map%mark(c,i) /= mark00(j) )
            !print *,mark(c,i), mark00(j)
            j = j+1
            if ( j > size(mark00) ) exit
        end do

        if ( j > size(mark00) ) then
           	call log_mess('bad construction of genotype file. a marker from map file' &
           	//' is not finding on the genotype file ['//trim(dataset%map%mark(c,i))//']',ERROR_DEF)
       	    stop 1
        end if

        work%index_position_marker(j) = i
        work%index_chr(j) = c
      end do
     end do

     deallocate (mark00)
     deallocate (line_read)

     call log_mess('END SUBROUTINE : check_marker_name',DEBUG_DEF)
     end subroutine check_marker_name
!!***

!!****f* m_qtlmap_genotype/internal_read_file
!!  NAME
!!    internal_read_file
!!  DESCRIPTION
!!   read the genotypic file.
!!   This file contains the animals phenotypes at the markers.
!!   The first line gives the marker names, the markers must belong to the marker map file.
!!   For each animal, a line gives :
!!    * its ID (as decribed in the pedigree file)
!!    * the markers phenotypes, ranked following in the first line order .
!!   Each phenotype is made of 2 alleles, unordered. When an animal has no phenotype for a marker, both alleles must be given the missing value code as
!!   given in the parametrisation of the analysis.
!!  INPUTS
!!   file   : genotypic file
!!
!!  SOURCE
     subroutine internal_read_file(dataset,work)
         type(QTLMAP_DATASET)      ,intent(inout)     :: dataset
         type(WORK_GENOTYPE)       ,intent(inout)     :: work
         character(len=LEN_DEF)          :: buffer_char,name
         character(len=LEN_DEF)                  :: val_nmanque,w
         integer                                 :: eof,ios,l,i_ind,i,j,k
         logical                                 :: is_ok
         character(len=LEN_DEF), dimension(:,:,:) , allocatable     :: temp_marker_list
         type(GENOTYPE_BASE)  , pointer :: dga
         type(GENEALOGY_BASE) , pointer :: dg

         call log_mess('SUBROUTINE : internal_read_file',DEBUG_DEF)

         dga => dataset%genoAnimal
         dg => dataset%genea

         rewind (unit_genotype)

         eof = 0
         dga%nmes = 0

         do while ( eof == 0 )
           read(unit_genotype,*,iostat=eof) buffer_char
           if ( trim(buffer_char) /= '' .and. eof == 0 ) then
                dga%nmes = dga%nmes+1
           end if
         end do

         ! the first line is the header (marker name)
         dga%nmes = dga%nmes - 1
         call log_mess('find ['//trim(str(dga%nmes))//'] animals in the genotype file.',INFO_DEF)

         ! fix bug : parents and grand parent can be add later....
         allocate ( dga%numero(dga%nmes+dg%ngp+dg%ngm+dg%np+dg%nm))
         dga%numero=STRING_NOT_DEFINED

         allocate (temp_marker_list(dga%nmes,work%nmarker_in_genotype_file,2))

         rewind (unit_genotype)
         read(unit_genotype,*,iostat=eof) buffer_char
         !Recupere les nom de tous les individus pour l affichage d erreur
         i=1
         do while ( eof == 0 .and. i <= dga%nmes )
           read (unit_genotype,*,iostat=eof) dga%numero(i),&
            ((temp_marker_list(i,k,j),j=1,2),k=1,work%nmarker_in_genotype_file)
           if ( trim(dga%numero(i)) /= '' .and. eof == 0 ) then
                i = i+1
           end if
         end do

         call work%set_pheno_structure(dataset,dga%nmes,temp_marker_list)

         deallocate (temp_marker_list)

         call log_mess('number of animals defined in genotype file:'// trim(str(dga%nmes)),VERBOSE_DEF)
         call log_mess('END SUBROUTINE : internal_read_file',DEBUG_DEF)
     end subroutine internal_read_file
!!***

     subroutine set_pheno_structure(work,dataset,nmes,temp_marker_list)
         class(WORK_GENOTYPE)      ,intent(inout)            :: work
         type(QTLMAP_DATASET)      ,intent(inout)            :: dataset
         integer                   , intent(in)              :: nmes
         character(len=LEN_DEF), dimension(nmes,work%nmarker_in_genotype_file,2),intent(in) :: temp_marker_list

         character(len=LEN_DEF)                  :: val_nmanque
         integer                                 :: i_ind,i
         type(GENOTYPE_BASE)  , pointer :: dga

         call log_mess('SUBROUTINE : set_pheno_structure',DEBUG_DEF)

         dga => dataset%genoAnimal
         allocate ( dga%pheno(dataset%map%nchr,maxval(dataset%map%nmk),size(dga%numero),2))
         dga%pheno=dataset%genoAnimal%nmanque

         !first encode, unknown char
         val_nmanque = get_pheno(dga,1,dataset%genoAnimal%nmanque)
         !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i)
         !$OMP DO
         do i_ind=1,dga%nmes
             !******************************************************
             ! 05/07/2010 : Test si l individu a deja ete lu (plusieurs definition de typage existe)
             !              test beaucoup trop long
!             do i=1,(i_ind-1)
!               if ( numero(i_ind) == numero(i)) then
!                 call stop_on_error (1,file,l,'animal ['//trim(numero(i_ind))//&
!                 '] have two lines definitions of phenotype marker !!')
!               end if
!             end do

             ! save now in the pheno structure....
             do i=1,work%nmarker_in_genotype_file ! for all markers in the genotype file
               if ( work%index_position_marker(i) == INT_NOT_DEFINED) cycle
              if ( trim(temp_marker_list(i_ind,i,1)) == val_nmanque .or. trim(temp_marker_list(i_ind,i,2)) == val_nmanque) then
        if ( .not. (trim(temp_marker_list(i_ind,i,1)) == val_nmanque .and. trim(temp_marker_list(i_ind,i,2)) == val_nmanque) ) then
                        call log_mess("*** Only one phenotype marker for [ animal:"//trim(dga%numero(i_ind))//&
                                   "] [marker:"// trim(dataset%map%mark(work%index_chr(i),work%index_position_marker(i)))//&
                                   "] is not allowed *** [setting -->" &
                                  //trim(val_nmanque)//" "//trim(val_nmanque)//"]",WARNING_DEF)
                     end if
                     dga%pheno(work%index_chr(i),work%index_position_marker(i),i_ind,1) = dataset%genoAnimal%nmanque
                     dga%pheno(work%index_chr(i),work%index_position_marker(i),i_ind,2) = dataset%genoAnimal%nmanque
                 else
                   !$OMP CRITICAL
                   !w=trim(temp_marker_list(i_ind,i,1))
                   dga%pheno(work%index_chr(i),work%index_position_marker(i),i_ind,1) = &
                     set_pheno(dga,dataset%map,work%index_chr(i),temp_marker_list(i_ind,i,1))
                   !w=trim(temp_marker_list(i_ind,i,2))
                   dga%pheno(work%index_chr(i),work%index_position_marker(i),i_ind,2) = &
                     set_pheno(dga,dataset%map,work%index_chr(i),temp_marker_list(i_ind,i,2))
                   !w=''
                   !$OMP END CRITICAL
                 end if
             end do
         end do
         !$OMP END DO
         !$OMP END PARALLEL

    end subroutine set_pheno_structure

    subroutine release_work_genotype(work)
     class(WORK_GENOTYPE) , intent(inout) :: work

     deallocate (work%index_position_marker,work%index_chr)

    end subroutine release_work_genotype

!!****f* m_qtlmap_genotype/allocate_vector_
!!  NAME
!!    allocate_vector_
!!  DESCRIPTION
!!   Allocation of main arrays coming from m_qtlmap_data :  corregp,corregm,correr,correp,correm,corred,presentg
!!
!!  SOURCE
    subroutine allocate_vector_(dataset)
         type(QTLMAP_DATASET)      ,intent(inout)            :: dataset
         integer                            :: ios
         type(GENOTYPE_BASE) , pointer :: dga
         type(GENEALOGY_BASE) , pointer :: dg

         dga => dataset%genoAnimal
         dg => dataset%genea

         allocate (dga%corregp(dg%ngp))
         allocate (dga%corregm(dg%ngm))
         allocate (dga%correr(dg%nr))
         allocate (dga%correp(dg%np))
         allocate (dga%correm(dg%nm))
         allocate (dga%corred(dg%nd))
         allocate (dga%presentg(dataset%map%nchr,dg%nd))

         call log_mess('END SUBROUTINE : allocate_vector_correspond',DEBUG_DEF)
    end subroutine allocate_vector_
!!***

!!****f* m_qtlmap_genotype/set_corresponding_vector
!!  NAME
!!    set_corresponding_vector
!!  DESCRIPTION
!!   * Build the correspondence vector between  genealogy array and the array numero (name of genotyped animal) : corregp,corregm,correp,correm,correr,corred
!!   * initialize presentg array
!!
!!  SOURCE
    subroutine  set_corresponding_vector(dataset)
        type(QTLMAP_DATASET)      ,intent(inout)            :: dataset
        integer                               :: igp,jgm,kr,ip,jm,kd,mes,ich
        integer :: origin_nmes

        type(GENOTYPE_BASE) , pointer :: dga
        type(GENEALOGY_BASE) , pointer :: dg

        dga => dataset%genoAnimal
        dg => dataset%genea

        origin_nmes = dga%nmes

      ! Creation d'un vecteur de correspondance pour les grands peres (corregp)
         do igp=1,dg%ngp
          dga%corregp(igp)=INT_NOT_DEFINED
          do mes=1,dga%nmes
           if(dg%gpere(igp) == dga%numero(mes)) then
             dga%corregp(igp)=mes
             exit
           endif
          end do
          !Ajout OFI Juin 2009 --> Add unknown genotype to avoid error with array acces
          if ( dga%corregp(igp) == INT_NOT_DEFINED ) then
             call log_mess('The grand sire '// trim(dg%gpere(igp))// ' has no genotype',WARNING_DEF)
             dga%nmes = dga%nmes + 1
             dga%numero(dga%nmes) = dg%gpere(igp)
             dga%corregp(igp)=dga%nmes
             dga%pheno(:,:,dga%corregp(igp),:)=dataset%genoAnimal%nmanque
          end if

         end do
!         if (t_imp) write(16,4000)dg%gpere(igp)
 !4000    format(1x,'The grand sire ',i6,' has no genotype')
  ! 42   continue

! Creation d'un vecteur de correspondance pour les grands meres (corregm)
         do jgm=1,dg%ngm
           dga%corregm(jgm)=INT_NOT_DEFINED
           do mes=1,dga%nmes
            if(dg%gmere(jgm) == dga%numero(mes)) then
             dga%corregm(jgm)=mes
             exit
            endif
          end do

           !Ajout OFI Juin 2009 --> Add unknown genotype to avoid error with array acces
          if ( dga%corregm(jgm) == INT_NOT_DEFINED ) then
             call log_mess('The grand dam '// trim(dg%gmere(jgm))// ' has no genotype',WARNING_DEF)
             dga%nmes = dga%nmes + 1
             dga%numero(dga%nmes) = dg%gmere(jgm)
             dga%corregm(jgm)=dga%nmes
             dga%pheno(:,:,dga%corregm(jgm),:)=dataset%genoAnimal%nmanque
          end if
         end do

! Creation d'un vecteur de correspondance pour les peres (correp)
        do ip=1,dg%np
         dga%correp(ip)=INT_NOT_DEFINED
         do mes=1,dga%nmes
          if(dg%pere(ip) == dga%numero(mes)) then
            dga%correp(ip)=mes
            exit
          endif
         end do

         !Ajout OFI Juin 2009 --> Add unknown genotype to avoid error with array acces
         if ( dga%correp(ip) == INT_NOT_DEFINED ) then
             call log_mess('The sire '// trim(dg%pere(ip))// ' has no genotype',WARNING_DEF)
             dga%nmes = dga%nmes + 1
             dga%numero(dga%nmes) = dg%pere(ip)
             dga%correp(ip)=dga%nmes
             dga%pheno(:,:,dga%correp(ip),:)=dataset%genoAnimal%nmanque
          end if

        end do
 !       if (t_imp) write(16,4002)pere(ip)
 !4002   format(1x,'The sire ',i6,' has no genotype')
 !  44   continue

 !Creation d'un vecteur de correspondance pour les meres (correm)
        do jm=1,dg%nm
         dga%correm(jm)=INT_NOT_DEFINED
         do mes=1,dga%nmes
          if(dg%mere(jm) == dga%numero(mes)) then
            dga%correm(jm)=mes
            exit
          endif
         end do

         !Ajout OFI Juin 2009 --> Add unknown genotype to avoid error with array acces
         if ( dga%correm(jm) == INT_NOT_DEFINED ) then
             call log_mess('The dam '// trim(dg%mere(jm))// ' has no genotype',WARNING_DEF)
             dga%nmes = dga%nmes + 1
             dga%numero(dga%nmes) = dg%mere(jm)
             dga%correm(jm)=dga%nmes
             dga%pheno(:,:,dga%correm(jm),:)=dataset%genoAnimal%nmanque
          end if
!         if (t_imp) write(16,4003)mere(jm)
! 4003    format(1x,'The dam ',i6,' has no genotype')
!   45   continue
        end do

        ! Creation d'un vecteur de correspondance pour les parents (correr)
        do kr=1,dg%nr
         dga%correr(kr)=INT_NOT_DEFINED
         do mes=1,dga%nmes
          if(dg%repro(kr) == dga%numero(mes)) then
            dga%correr(kr)=mes
              exit
          endif
         end do

        end do

! Creation d'un vecteur de correspondance pour les descendants (corred)
        do kd=1,dg%nd
            dga%presentg(:,kd)=.false.
            dga%corred(kd)=INT_NOT_DEFINED
         do mes=1,dga%nmes
          if(dg%animal(kd) == dga%numero(mes)) then
            dga%corred(kd)=mes
            do ich=1,dataset%map%nchr
               dga%presentg(ich,kd)=.not. (all(dga%pheno(ich,:,mes,:)==dataset%genoAnimal%nmanque))
            end do
            exit
          endif
         end do
        end do

        dga%nmes = origin_nmes
    end subroutine set_corresponding_vector
!!***

!!****f* m_qtlmap_genotype/set_allele_info_vector
!!  NAME
!!    set_allele_info_vector
!!  DESCRIPTION
!!   Set basic statistic and information about genotypic information (frequency, number of allele by marker,...)
!!
!!  SOURCE
    subroutine   set_allele_info_vector(dataset)
         type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
         integer                          :: i,j,kd,imes,k,nb,kmes,ios,c
         logical  , dimension(:) , allocatable       :: adress
         logical                                     :: tri_termine
         character(len=LEN_DEF)                        :: temp
         ! dim : nb marker,2 * nmes
         integer(kind=KIND_PHENO) , dimension (:,:),allocatable    :: mm

         type(GENOTYPE_BASE) , pointer :: dga
         type(GENEALOGY_BASE) , pointer :: dg

         dga => dataset%genoAnimal
         dg => dataset%genea

         allocate (mm(maxval(dataset%map%nmk),size(dga%numero)*2))
         allocate (dga%alleles(dataset%map%nchr,maxval(dataset%map%nmk),size(dga%numero)*2))
         allocate (dga%pc_all(dataset%map%nchr,maxval(dataset%map%nmk),size(dga%numero)*2))
         allocate (dga%nall(dataset%map%nchr,maxval(dataset%map%nmk)))

         do c=1,dataset%map%nchr
            mm = dataset%genoAnimal%nmanque
            call log_mess('['//trim(dataset%map%chromo(c))//&
             '] count number of allele by marker,sort and computing frequencies...',VERBOSE_DEF)
            do i=1,dataset%map%nmk(c)
              do j=1,dga%nmes
               mm(i,j) = dga%pheno(c,i,j,1)
               mm(i,dga%nmes+j) = dga%pheno(c,i,j,2)
              end do
             end do
         !if none information are found for a
          do i=1,dga%nmes
           do j=1,dataset%map%nmk(c)
              if ( dga%pheno(c,j,i,1) /= dataset%genoAnimal%nmanque .and. &
               dga%pheno(c,j,i,2) /= dataset%genoAnimal%nmanque) then
                exit ! this animal is ok!!
              end if
           end do

           if ( j > dataset%map%nmk(c) ) then
             do  kd=1,size(dg%animal)
                if (dg%animal(kd) == dga%numero(i)) then
                   dga%presentg(c,kd)=.false.
                   exit
                end if
             end do
           end if
         end do


 ! COMPTAGE DES ALLELES PAR MARQUEUR ET CREATION D'UN VECTEUR D'ALLELES par MARQUEUR
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(adress,imes,kmes,tri_termine,temp,k,nb)
      allocate (adress(size(dga%numero)*2))

      !$OMP DO
      do i=1,dataset%map%nmk(c)
        dga%nall(c,i)=0
        adress= .true.
        do imes=1,dga%nmes*2
           if (mm(i,imes)==dataset%genoAnimal%nmanque) then
              adress(imes)=.false.
           else
              do kmes=1,imes-1
!
!  debut modif JME le 11/03/2010

                 if (mm(i,imes)==mm(i,kmes)) then
                    adress(imes)=.false.
                    exit
                 end if
!  fin modif JME le 11/03/2010
       !          if (mm(i,imes)==mm(i,kmes)) adress(imes)=.false.
              enddo
           endif
        enddo
        do imes=1,dga%nmes*2
          if (adress(imes))then
            dga%nall(c,i)=dga%nall(c,i)+1
            dga%alleles(c,i,dga%nall(c,i))=get_pheno(dga,c,mm(i,imes))
          endif
 ! TRIE DES VECTEUR D'ALLELES  PAR MARQUEUR
          tri_termine = .true.
          do while ( tri_termine )
            tri_termine = .false.
            do k=2,dga%nall(c,i)
              if ( dga%alleles(c,i,k).LT.dga%alleles(c,i,k-1) ) then
                 tri_termine = .true.
                 temp = dga%alleles(c,i,k-1)
                 dga%alleles(c,i,k-1) = dga%alleles(c,i,k)
                 dga%alleles(c,i,k) = temp
             endif
            enddo
          enddo
        enddo

! CALCUL DE LA FREQUENCE DE CHAQUE ALLELE
        do k=1,dga%nall(c,i)
          nb=0
          do imes=1,dga%nmes*2
             if (get_pheno(dga,c,mm(i,imes))==dga%alleles(c,i,k)) nb=nb+1
          enddo
          dga%pc_all(c,i,k)=dble(nb)*100.d0/(dble(dga%nmes)*2.d0)
        enddo
      enddo
      !$OMP END DO
      deallocate ( adress )
      !$OMP END PARALLEL
     end do ! end chr

      deallocate ( mm )

    end subroutine set_allele_info_vector
!!***

!!****f* m_qtlmap_genotype/check_typage
!!  NAME
!!    check_typage
!!  DESCRIPTION
!!   check_typage
!!     check phenotype marker of the progeny
!!     p1,p2 phenotype marker of the kd progeny
!!     sp1,sp2  phenotype marker of kd's sire
!!     dp1,dp2  phenotype marker of kd's dam
!!
!!     if p1==sp1 .or. sp2 then p2==dp1 .or. dp2
!!     if p1==dp1 .or. dp2 then p2==sp1 .or. sp2
!!
!!  SOURCE
  subroutine check_typage(dataset)
       type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
       integer :: ll,ip,jm,kd,igp,jgm,kr
       logical :: typp,typm,typgp,typgm
       integer :: i,c,nberr,typerr
       integer ,parameter :: MAX_ERR=100
       integer (kind=KIND_PHENO) :: error1(MAX_ERR,2),error2(MAX_ERR,2),error3(MAX_ERR,2)
       integer  :: kderror(MAX_ERR),iperror(MAX_ERR),jmerror(MAX_ERR),chrerror(MAX_ERR)
       character(len=LEN_DEF) :: merror(MAX_ERR)
       type(GENOTYPE_BASE) , pointer :: dga
       type(GENEALOGY_BASE) , pointer :: dg

       dga => dataset%genoAnimal
       dg => dataset%genea

      nberr=0
      do c=1,dataset%map%nchr
       !for each marker
       do ll=1,dataset%map%nmk(c)
          do ip=1,dg%np ! for each sire
            typp=.false.
            if(dga%correp(ip) /= INT_NOT_DEFINED)then
              if( dga%pheno(c,ll,dga%correp(ip),1) /= dataset%genoAnimal%nmanque ) then
                 typp=.true.
              endif
	        end if

            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
              typm=.false.
              if(dga%correm(jm) /= INT_NOT_DEFINED)then
                if(dga%pheno(c,ll,dga%correm(jm),1) /= dataset%genoAnimal%nmanque) then
                  typm=.true.
                endif
	          end if

              ! MODIF OFI 21 janv 2010 : (.not. typp) .or. (.not. typm) => (.not. typp) .and. (.not. typm)
	          !no information about parents
	          if ( (.not. typp) .and. (.not. typm) ) then
	        !      call log_mess('none information find about sire ['//trim(dg%pere(ip))//&
	        !      '] dam ['//trim(dg%mere(jm))//'] for marker ['//trim(dataset%map%mark(c,ll))//']',DEBUG_DEF)
	              cycle
	          end if

	          !checking progenies
	          do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
                if ( dga%corred(kd) == INT_NOT_DEFINED ) then
                   cycle
                end if

                if ( dga%pheno(c,ll,dga%corred(kd),1) == dataset%genoAnimal%nmanque .or.&
                  dga%pheno(c,ll,dga%corred(kd),2) == dataset%genoAnimal%nmanque) then
                   cycle
                end if

                  if (.not. typp .or. dga%pheno(c,ll,dga%corred(kd),1) == dga%pheno(c,ll,dga%correp(ip),1) .or.  &
                      dga%pheno(c,ll,dga%corred(kd),1) == dga%pheno(c,ll,dga%correp(ip),2)) then
                       !if ( correm(jm) == INT_NOT_DEFINED ) cycle
                       if (.not. typm .or. (dga%pheno(c,ll,dga%corred(kd),2) == dga%pheno(c,ll,dga%correm(jm),1)) .or. &
                          (dga%pheno(c,ll,dga%corred(kd),2) == dga%pheno(c,ll,dga%correm(jm),2)) ) then
                           cycle
                       end if
                  end if

                if ( .not. typm                                      .or. &
                    dga%pheno(c,ll,dga%corred(kd),1) == dga%pheno(c,ll,dga%correm(jm),1) .or.  &
                    dga%pheno(c,ll,dga%corred(kd),1) == dga%pheno(c,ll,dga%correm(jm),2)) then

                  if ( .not. typp                                    .or.  &
                       (dga%pheno(c,ll,dga%corred(kd),2) == dga%pheno(c,ll,dga%correp(ip),1)) .or. &
                       (dga%pheno(c,ll,dga%corred(kd),2) == dga%pheno(c,ll,dga%correp(ip),2)) ) then
                      cycle
                  end if
                end if

                if (nberr >= MAX_ERR) exit
                nberr=nberr+1
                error1(nberr,1)=dga%pheno(c,ll,dga%correp(ip),1)
                error2(nberr,1)=dga%pheno(c,ll,dga%correm(jm),1)
                error3(nberr,1)=dga%pheno(c,ll,dga%corred(kd),1)
                error1(nberr,2)=dga%pheno(c,ll,dga%correp(ip),2)
                error2(nberr,2)=dga%pheno(c,ll,dga%correm(jm),2)
                error3(nberr,2)=dga%pheno(c,ll,dga%corred(kd),2)
                kderror(nberr)=kd
                iperror(nberr)=ip
                jmerror(nberr)=jm
                merror(nberr)=dataset%map%mark(c,ll)
                chrerror(nberr)=c

            end do ! end progeny
           end do ! end dam
         end do ! end sire
       end do ! end marker
     end do ! end chr
     typerr=nberr
     if ( nberr < MAX_ERR ) then
!! CHECKING Reproducter..........
      do c=1,dataset%map%nchr
       do ll=1,dataset%map%nmk(c)
          do igp=1,dg%ngp
           typgp=.false.
           if(dga%corregp(igp) /= INT_NOT_DEFINED) then
             if(dga%pheno(c,ll,dga%corregp(igp),1) /= dataset%genoAnimal%nmanque) then
              typgp=.true.
             end if
           end if

           do jgm=dg%ngmgp(igp)+1,dg%ngmgp(igp+1)
             typgm = .false.
             if(dga%corregm(jgm) /= INT_NOT_DEFINED) then
               if(dga%pheno(c,ll,dga%corregm(jgm),1) /= dataset%genoAnimal%nmanque) then
                typgm=.true.
               end if
             end if
            !no information about parents
	        if ( (.not. typgp) .or. (.not. typgm) ) then
	        !      call log_mess('none information find about grand sire ['//trim(dg%gpere(igp))//&
	        !      '] grand dam ['//trim(dg%gmere(jgm))//'] for marker ['//trim(dataset%map%mark(c,ll))//']',DEBUG_DEF)
	              cycle
	        end if

            do kr=dg%nrgm(jgm)+1,dg%nrgm(jgm+1)
               if ( dga%correr(kr) == INT_NOT_DEFINED ) then
                   cycle
                end if

               if ( dga%pheno(c,ll,dga%correr(kr),1) == dataset%genoAnimal%nmanque .or.&
                 dga%pheno(c,ll,dga%correr(kr),2) == dataset%genoAnimal%nmanque) then
                   cycle
               end if
                 if ( .not. typgp                                      .or.  &
                    dga%pheno(c,ll,dga%correr(kr),1) == dga%pheno(c,ll,dga%corregp(igp),1) .or.  &
                    dga%pheno(c,ll,dga%correr(kr),1) == dga%pheno(c,ll,dga%corregp(igp),2)) then
                  if ( .not. typgm                                    .or.  &
                       (dga%pheno(c,ll,dga%correr(kr),2) == dga%pheno(c,ll,dga%corregm(jgm),1)) .or. &
                       (dga%pheno(c,ll,dga%correr(kr),2) == dga%pheno(c,ll,dga%corregm(jgm),2)) ) then
                       cycle
                  end if
                end if

                if ( .not. typgm                                      .or.  &
                    dga%pheno(c,ll,dga%correr(kr),1) == dga%pheno(c,ll,dga%corregm(jgm),1) .or.  &
                    dga%pheno(c,ll,dga%correr(kr),1) == dga%pheno(c,ll,dga%corregm(jgm),2)) then

                  if ( .not. typgp                                    .or.  &
                       (dga%pheno(c,ll,dga%correr(kr),2) == dga%pheno(c,ll,dga%corregp(igp),1)) .or. &
                       (dga%pheno(c,ll,dga%correr(kr),2) == dga%pheno(c,ll,dga%corregp(igp),2)) ) then
                      cycle
                  end if
                end if
                if (nberr >= MAX_ERR) exit
                nberr=nberr+1
                error1(nberr,1)=dga%pheno(c,ll,dga%corregp(igp),1)
                error2(nberr,1)=dga%pheno(c,ll,dga%corregm(jgm),1)
                error3(nberr,1)=dga%pheno(c,ll,dga%correr(kr),1)
                error1(nberr,2)=dga%pheno(c,ll,dga%corregp(igp),2)
                error2(nberr,2)=dga%pheno(c,ll,dga%corregm(jgm),2)
                error3(nberr,2)=dga%pheno(c,ll,dga%correr(kr),2)
                kderror(nberr)=kr
                iperror(nberr)=igp
                jmerror(nberr)=jgm
                merror(nberr)=dataset%map%mark(c,ll)
                chrerror(nberr)=c
            end do ! end repro
          end do ! end grand dam
        end do ! end grand sire
      end do ! end marker
     end do ! end chr
    end if

     if ( nberr > 0 ) then
        do i=1,typerr
             c=chrerror(i)
             call log_mess('     ----------   ',ERROR_DEF)
             call log_mess('Phenotype sire    :'//trim(get_pheno(dga,c,error1(i,1)))//' '//&
                    trim(get_pheno(dga,c,error1(i,2))),ERROR_DEF)
             call log_mess('Phenotype dam     :'//trim(get_pheno(dga,c,error2(i,1)))//' '//&
                    trim(get_pheno(dga,c,error2(i,2))),ERROR_DEF)
             call log_mess('Phenotype progeny :'//trim(get_pheno(dga,c,error3(i,1)))//' '//&
                    trim(get_pheno(dga,c,error3(i,2))),ERROR_DEF)
             call log_mess('progeny ['//trim(dg%animal(kderror(i)))//'] from sire ['//trim(dg%pere(iperror(i)))//&
                      '] dam ['//trim(dg%mere(jmerror(i)))//'] for marker ['//trim(merror(i))//']',ERROR_DEF)

        end do

        do i=typerr+1,nberr
             c=chrerror(i)
             call log_mess('     ----------   ',ERROR_DEF)
             call log_mess('Phenotype grand sire    :'//trim(get_pheno(dga,c,error1(i,1)))//' '//&
                    trim(get_pheno(dga,c,error1(i,2))),ERROR_DEF)
             call log_mess('Phenotype grand dam     :'//trim(get_pheno(dga,c,error2(i,1)))//' '//&
                    trim(get_pheno(dga,c,error2(i,2))),ERROR_DEF)
             call log_mess('Phenotype reproductor   :'//trim(get_pheno(dga,c,error3(i,1)))//' '//&
                    trim(get_pheno(dga,c,error3(i,2))),ERROR_DEF)
             call log_mess('reproductor ['//trim(dg%repro(kderror(i)))//'] from grand sire ['//trim(dg%gpere(iperror(i)))//&
                    '] grand dam ['//trim(dg%gmere(jmerror(i)))//'] for marker ['//trim(merror(i))//']',ERROR_DEF)
        end do
       call stop_application(' ** NB ERROR :'// trim(str(nberr))//'** ')
     end if
    end subroutine check_typage
!!***

!!****f* m_qtlmap_genotype/sim_typ
!!  NAME
!!    sim_typ
!!  DESCRIPTION
!!
!!  NOTES
!!   Simulation des typages sur trois generations , reperage des alleles au qtl des descendants.
!!  SOURCE
      subroutine sim_typ(dataset)
!
      type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
! Divers
      logical :: typ2
      integer :: c,p
      double precision liminf,limsup
      integer gpare,k,pare,i,j,iq,jgm,ngm1,ngm2,igp
      integer l1,l2,nr1,nr2,ir,ip,nm1,nm2,jm,id,ic,mes,nd1,nd2
      real (kind=dp) :: r
      real           :: x_rand

      real,  external :: ranf
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_BASE) , pointer :: dg

      dga => dataset%genoAnimal
      dg => dataset%genea
!
!
!*******************************************************************
!         Lecture parametres de simulation du QTL
!*******************************************************************
!
      dga%pheno=dataset%genoAnimal%nmanque
!
!***********************************************************************
!                Simulation des typages des grands-parents
!
!***********************************************************************
! Alleles des Grands peres
!
      if (.not. allocated (dga%pc_all) ) then
        call stop_application('Devel error: pc_all array is not allocated. Can not simulate type.')
      end if

      if (.not. allocated (dga%alleles) ) then
        call stop_application('Devel error: all array is not allocated. Can not simulate type.')
      end if


      mes=0
      do igp=1,dg%ngp
       mes=mes+1
       dga%numero(mes)=dg%gpere(igp)
       dga%corregp(igp)=mes
! Alleles aux marqueurs
       do c=1,dataset%map%nchr
        do j=1,2
         do i=1,dataset%map%nmk(c)
          x_rand=ranf()
          limsup=0.d0
          do k=1,dga%nall(c,i)
           liminf=limsup
           limsup=liminf+dga%pc_all(c,i,k)/1.d2
           if(dble(x_rand).gt.liminf.and.dble(x_rand).le.limsup) then
            dga%pheno(c,i,mes,j)=set_pheno(dga,dataset%map,c,dga%alleles(c,i,k))
           end if
          end do
         end do
        end do
       end do
      end do
!
! Alleles des Grands meres
      do jgm=1,dg%ngm
       mes=mes+1
       dga%numero(mes)=dg%gmere(jgm)
       dga%corregm(jgm)=mes
!
! Alleles aux marqueurs
      do c=1,dataset%map%nchr
       do j=1,2
        do i=1,dataset%map%nmk(c)
         x_rand=ranf()
         limsup=0.d0
         do k=1,dga%nall(c,i)
          liminf=limsup
          limsup=liminf+dga%pc_all(c,i,k)/1.d2
          if(dble(x_rand).gt.liminf.and.dble(x_rand).le.limsup) then
            dga%pheno(c,i,mes,j)=set_pheno(dga,dataset%map,c,dga%alleles(c,i,k))
          end if
         end do
        end do
       end do
      end do
     end do !! jgm
!
!
!**********************************************************************
!           Simulation de la transmission des phases aux parents
! Les parents qui sont aussi gparents doivent etre enregistres les premiers
!
!**********************************************************************
!
      do igp=1,dg%ngp
        ngm1=dg%ngmgp(igp)+1
        ngm2=dg%ngmgp(igp+1)
        do jgm=ngm1,ngm2
          nr1=dg%nrgm(jgm)+1
          nr2=dg%nrgm(jgm+1)
          do ir=nr1,nr2
             mes=mes+1
             dga%numero(mes)=dg%repro(ir)
             dga%correr(ir)=mes
!
! Deux chromosomes a la suite : p=1 gd pere ; p=2 gd mere
             do p=1,2
               if (p.eq.1) gpare=dga%corregp(igp)
               if (p.eq.2) gpare=dga%corregm(jgm)
               l1=1
               k=0
               iq=1
! Premier marqueur
! Choix de la phase heritee au 1er marqueur : k
               x_rand=ranf()
               if (dble(x_rand).le.0.5d0) then
                 dga%pheno(:,l1,mes,p)=dga%pheno(:,l1,gpare,1)
                 k=1
                else
                 dga%pheno(:,l1,mes,p)=dga%pheno(:,l1,gpare,2)
                 k=2
               end if
!
! Marqueurs suivants
             do c=1,dataset%map%nchr
               do l2=2,dataset%map%nmk(c)
                 if (p.eq.1) r=(dataset%map%rm(c,l1,l2))
                 if (p.eq.2) r=(dataset%map%rf(c,l1,l2))
                 x_rand=ranf()
                 typ2=.true.
!
! Il y a recombinaison
                 if (dble(x_rand).lt.r) then
! Phase 1 heritee -> echange
                   if (k .eq. 1 .and. typ2) then
                    dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,gpare,3-k)
!
                    k=2
                    typ2=.false.
                   end if
!
! Phase 2 heritee -> echange si pas deja fait de phase 1 a 2
                   if (k .eq. 2 .and. typ2) then
                    dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,gpare,3-k)
                    k=1
!
                   end if
                 end if

!
! Pas de recombinaison
                 if (dble(x_rand).ge.r)dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,gpare,k)
!
                 l1=l1+1
               end do            !! fin chr
              end do
             end do                !!fin deux gparents
!
         end do                !! repro
       end do                  !! gmere
      end do                    !! gpere
!
!
!**********************************************************************
!       Simulation de la transmission des phases aux descendants
!
!**********************************************************************
!
      do ip=1,dg%np
        dga%correp(ip)=dga%correr(dg%reppere(ip))
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
          dga%correm(jm)=dga%correr(dg%repmere(jm))
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          do id=nd1,nd2
            mes=mes+1
            dga%numero(mes)=dg%animal(id)
            dga%corred(id)=mes
            do p=1,2
             if (p.eq.1) pare=dga%correp(ip)
             if (p.eq.2) pare=dga%correm(jm)
             l1=1
             k=0
             iq=1
! Premier marqueur
             x_rand=ranf()
             if (dble(x_rand).le.0.5d0) then
               dga%pheno(:,l1,mes,p)=dga%pheno(:,l1,pare,1)
               k=1
             else
               dga%pheno(:,l1,mes,p)=dga%pheno(:,l1,pare,2)
               k=2
             end if
! Marqueurs suivants
           do c=1,dataset%map%nchr
             do l2=2,dataset%map%nmk(c)
               if (p.eq.1) r=(dataset%map%rm(c,l1,l2))
               if (p.eq.2) r=(dataset%map%rf(c,l1,l2))
              x_rand=ranf()
               typ2=.true.
!
               if (dble(x_rand).lt.r) then
                 if (k .eq. 1 .and. typ2) then
                   dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,pare,2)
                   k=2
                   typ2=.false.
                   goto 10
                 end if
!
                 if (k .eq. 2 .and. typ2) then
                   dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,pare,1)
                   k=1
                   goto 10
                 end if
               end if
!
               if (dble(x_rand).ge.r) then
                 if (k.eq.1) dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,pare,1)
                 if (k.eq.2) dga%pheno(c,l2,mes,p)=dga%pheno(c,l2,pare,2)
                 goto 10
               end if
!
 10            continue
               l1=l1+1
              end do            !! fin chr
             end do
            end do                !!fin deux gparents

         end do                !!  animal
        end do                  !! mere
      end do                    !! pere
!
!**********************************************************************
!           Ecriture des genotypes dans le fichier typ_sim
!**********************************************************************
!
       dga%nmes=mes

      call set_estfem(dataset)
      end subroutine sim_typ
!!***

!!****f* m_qtlmap_genotype/recup
!!  NAME
!!    recup
!!  DESCRIPTION
!!
!!  NOTES
!!   Deduction des phenotypes des reproducteurs non types d apres leur descendance
!!  SOURCE
      subroutine recup(dataset)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset

! Divers
        logical typp,typm,typgp,typgm
        integer :: ll,igp,nagp,ngm1,ngm2,jgm,nagm,nr1,nr2,kr,il
        integer :: jmes,imes,ip,nap,nm1,nm2,jm,nam,nd1,nd2,kd,c
        integer(kind=KIND_PHENO) :: m1,m2,mc1,mc2,md1,md2,mx
        type(GENOTYPE_BASE) , pointer :: dga
        type(GENEALOGY_BASE) , pointer :: dg

        dga => dataset%genoAnimal
        dg => dataset%genea

!
!******************************************************************************
! Reconstitution des phenotypes des grands parents
!******************************************************************************
!

        call log_mess('finding phenotype marker...',INFO_DEF)
! Recherche des phenotypes des grands parents
      do c=1,dataset%map%nchr
       do ll=1,dataset%map%nmk(c)
        do 10 igp=1,dg%ngp
           if(dga%corregp(igp).eq.9999) then
             typgp=.false.
              nagp=0
           else if(dga%pheno(c,ll,dga%corregp(igp),1).eq.dataset%genoAnimal%nmanque) then
              typgp=.false.
              nagp=0
           else
            typgp=.true.
              nagp=2
            m1=dga%pheno(c,ll,dga%corregp(igp),1)
            m2=dga%pheno(c,ll,dga%corregp(igp),2)
            end if
          ngm1=dg%ngmgp(igp)+1
          ngm2=dg%ngmgp(igp+1)
          do 11 jgm=ngm1,ngm2
             if(dga%corregm(jgm).eq.9999) then
                typgm=.false.
                nagm=0
             else if(dga%pheno(c,ll,dga%corregm(jgm),1).eq.dataset%genoAnimal%nmanque) then
                typgm=.false.
                nagm=0
             else
              typgm=.true.
              nagm=2
              mc1=dga%pheno(c,ll,dga%corregm(jgm),1)
              mc2=dga%pheno(c,ll,dga%corregm(jgm),2)
            end if
            if(typgp.and.typgm)go to 11

!            if ( dga%corregp(igp) /= 9999 ) then
!              call log_mess('Grand sire ['//trim(dg%gpere(igp))//'] -->'//trim(get_dga%pheno(c,dga%pheno(c,ll,dga%corregp(igp),1)))//&
!              ' '//trim(get_dga%pheno(c,dga%pheno(c,ll,dga%corregp(igp),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!            end if
!
!            if ( dga%corregm(jgm) /= 9999 ) then
!            call log_mess('Grand dam ['//trim(gmere(jgm))//'] -->'//trim(get_dga%pheno(c,dga%pheno(c,ll,dga%corregm(jgm),1)))//' '//&
!            trim(get_dga%pheno(c,dga%pheno(c,ll,dga%corregm(jgm),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!            end if
!
! Passage en revue des parents
            nr1=dg%nrgm(jgm)+1
            nr2=dg%nrgm(jgm+1)
            do 12 kr=nr1,nr2
               if(dga%correr(kr).eq.9999)  go to 12
               if(dga%pheno(c,ll,dga%correr(kr),1).eq.dataset%genoAnimal%nmanque) go to 12
              md1=dga%pheno(c,ll,dga%correr(kr),1)
              md2=dga%pheno(c,ll,dga%correr(kr),2)
              mx=dataset%genoAnimal%nmanque

!              if ( dga%correr(kr) /= 9999 ) then
!                 call log_mess('Repro ['//trim(repro(kr))//'] -->'//trim(get_pheno(dga,c,pheno(c,ll,dga%correr(kr),1)))//' '//&
!                 trim(get_pheno(dga,c,pheno(c,ll,dga%correr(kr),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!              end if
!
! Parent homozygote
! Grand pere inconnu
              if(md1.eq.md2)then
             !   call log_mess('Reproductor homo ['//trim(get_pheno(dga,c,md1))//' '//trim(get_pheno(dga,c,md2))//']',DEBUG_DEF)
                if(nagp.ne.2) then
                    if(nagp.eq.0) then
                      nagp=1
                      m1=md1
                  !    call log_mess('finding one allele for grand sire ['//trim(get_pheno(dga,c,m1))//']',DEBUG_DEF)
                    else
                      if(md1.ne.m1) then
                        nagp=2
                        m2=md1
                    !    call log_mess('finding one allele for grand sire ['//trim(get_pheno(dga,c,m2))//']',DEBUG_DEF)
                      end if
                    end if
               end if
!
! Grand mere inconnue
                if(nagm.ne.2) then
                    if(nagm.eq.0) then
                    nagm=1
                    mc1=md1
                !    call log_mess('finding one allele for grand dam ['//trim(get_pheno(dga,c,mc1))//']',DEBUG_DEF)
                  else
                      if(md1.ne.mc1) then
                        nagm=2
                        mc2=md1
                    !    call log_mess('finding one allele for grand dam ['//trim(get_pheno(dga,c,mc2))//']',DEBUG_DEF)
                        end if
                      end if
                   end if
!
! Parent heterozygote
! Grand mere connue et grand pere inconnu
              else
              !    call log_mess('Reproductor hetero ['//trim(get_pheno(dga,c,md1))//','//trim(get_pheno(dga,c,md2))//']',DEBUG_DEF)
                  if(nagm.eq.2.and.nagp.ne.2) then
                  if((md1.eq.mc1.and.md2.ne.mc2).or.(md1.eq.mc2.and.md2.ne.mc1)) mx=md2
                  if((md2.eq.mc1.and.md1.ne.mc2).or.(md2.eq.mc2.and.md1.ne.mc1)) mx=md1
                  if(mx.ne.dataset%genoAnimal%nmanque)then
                  if(nagp.eq.0) then
                    nagp=1
                    m1=mx
                !    call log_mess('finding one allele for grand sire ['//trim(get_pheno(dga,c,m1))//']',DEBUG_DEF)
                  else
                      if(mx.ne.m1) then
                      nagp=2
                      m2=mx
                  !    call log_mess('finding one allele for grand sire ['//trim(get_pheno(dga,c,m2))//']',DEBUG_DEF)
                    end if
                  end if
                 end if
                end if
!
! Grand pere connu et grand mere inconnue
                  if(nagp.eq.2.and.nagm.ne.2) then
                  if((md1.eq.m1.and.md2.ne.m2).or.(md1.eq.m2.and.md2.ne.m1)) mx=md2
                  if((md2.eq.m1.and.md1.ne.m2).or.(md2.eq.m2.and.md1.ne.m1)) mx=md1
                  if(mx.ne.dataset%genoAnimal%nmanque)then
                  if(nagm.eq.0) then
                      nagm=1
                      mc1=mx
                    !  call log_mess('finding one allele for grand dam ['//trim(get_pheno(dga,c,mc1))//']',DEBUG_DEF)
                  else
                      if(mx.ne.mc1) then
                      nagm=2
                      mc2=mx
                  !    call log_mess('finding one allele for grand dam ['//trim(get_pheno(dga,c,mc2))//']',DEBUG_DEF)
                    end if
                  end if
                 end if
                end if
              end if
   12       continue
!
! Stockage du phenotype reconstitue de la grand mere
            if(typgm.or.nagm.lt.2) go to 11
            if(dga%corregm(jgm).eq.9999) then
              dga%nmes=dga%nmes+1
              dga%numero(dga%nmes)=dg%gmere(jgm)
              dga%corregm(jgm)=dga%nmes
              dga%pheno(c,:,dga%corregm(jgm),1)=dataset%genoAnimal%nmanque
              dga%pheno(c,:,dga%corregm(jgm),2)=dataset%genoAnimal%nmanque
            end if
            jmes=dga%corregm(jgm)
            dga%pheno(c,ll,jmes,1)=mc1
            dga%pheno(c,ll,jmes,2)=mc2

            call log_mess('setting phenotype marker for grand dam ['//trim(dg%gmere(jgm))//'] :'//&
            trim(get_pheno(dga,c,dga%pheno(c,ll,dga%corregm(jgm),1)))//' '//&
            trim(get_pheno(dga,c,dga%pheno(c,ll,dga%corregm(jgm),2)))&
            //' marker ['//trim(dataset%map%mark(c,ll))//']',DEBUG_DEF)

   11     continue
!
! Stockage du phenotype reconstitue du grand pere
          if(typgp.or.nagp.lt.2) go to 10
          if(dga%corregp(igp).eq.9999) then
            dga%nmes=dga%nmes+1
            dga%numero(dga%nmes)=dg%gpere(igp)
            dga%corregp(igp)=dga%nmes
            dga%pheno(c,:,dga%corregp(igp),1)=dataset%genoAnimal%nmanque
            dga%pheno(c,:,dga%corregp(igp),2)=dataset%genoAnimal%nmanque
          end if
          imes=dga%corregp(igp)
          dga%pheno(c,ll,imes,1)=m1
          dga%pheno(c,ll,imes,2)=m2
!          call log_mess('setting phenotype marker for grand sire ['//trim(dg%gpere(igp))//'] :'//&
!            trim(get_pheno(dga,c,pheno(c,ll,dga%corregp(igp),1)))//' '//trim(get_pheno(dga,c,pheno(c,ll,dga%corregp(igp),2)))//&
!            ' marker ['//trim(mark(c,ll))//']',DEBUG_DEF)
   10   continue
        end do
!
!******************************************************************************
! Reconstitution des phenotypes des parents
!******************************************************************************
!
! Recherche des phenotypes des parents
        do ll=1,dataset%map%nmk(c)
        do 20 ip=1,dg%np
          if(dga%correp(ip).eq.9999) then
              typp=.false.
              nap=0
          else if (dga%pheno(c,ll,dga%correp(ip),1).eq.dataset%genoAnimal%nmanque) then
              typp=.false.
              nap=0
            else
            typp=.true.
            nap=2
            m1=dga%pheno(c,ll,dga%correp(ip),1)
            m2=dga%pheno(c,ll,dga%correp(ip),2)
            end if
          nm1=dg%nmp(ip)+1
          nm2=dg%nmp(ip+1)
          do 21 jm=nm1,nm2
            if(dga%correm(jm).eq.9999) then
                typm=.false.
                nam=0
            else if(dga%pheno(c,ll,dga%correm(jm),1).eq.dataset%genoAnimal%nmanque) then
                typm=.false.
                nam=0
              else
              typm=.true.
              nam=2
              mc1=dga%pheno(c,ll,dga%correm(jm),1)
              mc2=dga%pheno(c,ll,dga%correm(jm),2)
              end if
              if(typp.and.typm)go to 21

!              if ( dga%correp(ip) /= 9999 ) then
!                  call log_mess('Sire ['//trim(pere(ip))//'] -->'//trim(get_pheno(dga,c,pheno(c,ll,dga%correp(ip),1)))//&
!                  ' '//trim(get_pheno(dga,c,pheno(c,ll,dga%correp(ip),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!              end if
!
!              if ( dga%correm(jm) /= 9999 ) then
!                  call log_mess('Dam ['//trim(mere(jm))//'] -->'//trim(get_pheno(dga,c,pheno(c,ll,dga%correm(jm),1)))//' '&
!                  //trim(get_pheno(dga,c,pheno(c,ll,dga%correm(jm),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!              end if
! Passage en revue des descendants
            nd1=dg%ndm(jm)+1
            nd2=dg%ndm(jm+1)
            do 22 kd=nd1,nd2
              if(dga%corred(kd).eq.9999)go to 22
              if(dga%pheno(c,ll,dga%corred(kd),1).eq.dataset%genoAnimal%nmanque)go to 22
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              mx=dataset%genoAnimal%nmanque

!              if ( dga%corred(kd) /= 9999 ) then
!                  call log_mess('Progeny ['//trim(animal(kd))//'] -->'//trim(get_pheno(dga,c,pheno(c,ll,dga%corred(kd),1)))//&
!                  ' '//trim(get_pheno(dga,c,pheno(c,ll,dga%corred(kd),2)))//' marker:'//trim(mark(c,ll)),DEBUG_DEF)
!              end if
!
! Descendant homozygote
! Pere inconnu
              if(md1.eq.md2)then
             !   call log_mess('Progeny homo',DEBUG_DEF)
                if(nap.ne.2) then
                    if(nap.eq.0) then
                    nap=1
                    m1=md1
                 !   call log_mess('finding one allele for sire ['//trim(get_pheno(dga,c,m1))//']',DEBUG_DEF)
                  else
                      if(md1.ne.m1) then
                      nap=2
                      m2=md1
                   !   call log_mess('finding one allele for sire ['//trim(get_pheno(dga,c,m2))//']',DEBUG_DEF)
                      end if
                    end if
                  end if
!
! Mere inconnue
                if(nam.ne.2) then
                    if(nam.eq.0) then
                    nam=1
                    mc1=md1
                 !   call log_mess('finding one allele for dam ['//trim(get_pheno(dga,c,mc1))//']',DEBUG_DEF)
                  else
                      if(md1.ne.mc1) then
                      nam=2
                      mc2=md1
                   !   call log_mess('finding one allele for dam ['//trim(get_pheno(dga,c,mc2))//']',DEBUG_DEF)
                      end if
                    end if
                  end if
!
! Descendant heterozygote
! Mere connue et pere inconnu
              else
                !  call log_mess('Progeny hetero',DEBUG_DEF)
                  if(nam.eq.2.and.nap.ne.2) then
                  if((md1.eq.mc1.and.md2.ne.mc2).or.       &
                    (md1.eq.mc2.and.md2.ne.mc1)) mx=md2
                  if((md2.eq.mc1.and.md1.ne.mc2).or.       &
                    (md2.eq.mc2.and.md1.ne.mc1)) mx=md1
                  if(mx.ne.dataset%genoAnimal%nmanque) then
                  if(nap.eq.0) then
                    nap=1
                    m1=mx
                  !  call log_mess('finding one allele for sire ['//trim(get_pheno(dga,c,m1))//']',DEBUG_DEF)
                  else
                      if(mx.ne.m1) then
                      nap=2
                      m2=mx
                   !   call log_mess('finding one allele for sire ['//trim(get_pheno(dga,c,m2))//']',DEBUG_DEF)
                    end if
                  end if
                 end if
                end if
!
! Pere connu et mere inconnue
                if(nap.eq.2.and.nam.ne.2) then
                  if((md1.eq.m1.and.md2.ne.m2).or.   &
                    (md1.eq.m2.and.md2.ne.m1)) mx=md2
                  if((md2.eq.m1.and.md1.ne.m2).or.   &
                    (md2.eq.m2.and.md1.ne.m1)) mx=md1
                  if(mx.ne.dataset%genoAnimal%nmanque) then
                  if(nam.eq.0) then
                    nam=1
                    mc1=mx
                  !  call log_mess('finding one allele for dam ['//trim(get_pheno(dga,c,mc1))//']',DEBUG_DEF)
                  else
                      if(mx.ne.mc1) then
                      nam=2
                      mc2=mx
                    !  call log_mess('finding one allele for dam ['//trim(get_pheno(dga,c,mc2))//']',DEBUG_DEF)
                    end if
                  end if
                end if
               end if
              end if
   22       continue
!
! Stockage du phenotype reconstitue de la mere
            if(typm.or.nam.lt.2) go to 21
            if(dga%correm(jm).eq.9999)then
              dga%nmes=dga%nmes+1
              dga%numero(dga%nmes)= dg%mere(jm)
              dga%correm(jm)=dga%nmes
              dga%pheno(c,:,dga%correm(jm),1)=dataset%genoAnimal%nmanque
              dga%pheno(c,:,dga%correm(jm),2)=dataset%genoAnimal%nmanque
            end if



            jmes=dga%correm(jm)
            dga%pheno(c,ll,jmes,1)=mc1
            dga%pheno(c,ll,jmes,2)=mc2
!            call log_mess('setting phenotype marker for dam ['//trim(mere(jm))//'] :'//&
!            trim(get_pheno(dga,c,pheno(c,ll,dga%correm(jm),1)))//' '//trim(get_pheno(dga,c,pheno(c,ll,dga%correm(jm),2)))//&
!            ' marker ['//trim(mark(c,ll))//']',DEBUG_DEF)
   21       continue
!
! Stockage du phenotype reconstitue du pere
          if(typp.or.nap.lt.2) go to 20
          if(dga%correp(ip).eq.9999) then
            dga%nmes=dga%nmes+1
            dga%numero(dga%nmes)=dg%pere(ip)
            dga%correp(ip)=dga%nmes
            dga%pheno(c,:,dga%correp(ip),1)=dataset%genoAnimal%nmanque
            dga%pheno(c,:,dga%correp(ip),2)=dataset%genoAnimal%nmanque
          end if
          imes=dga%correp(ip)
          dga%pheno(c,ll,imes,1)=m1
          dga%pheno(c,ll,imes,2)=m2
!          call log_mess('setting phenotype marker for sire ['//trim(pere(ip))//'] :'//&
!            trim(get_pheno(dga,c,pheno(c,ll,dga%correp(ip),1)))//' '//trim(get_pheno(dga,c,pheno(c,ll,dga%correp(ip),2)))//&
!            ' marker ['//trim(mark(c,ll))//']',DEBUG_DEF)
   20   continue

        end do
       end do ! end nchr
     end subroutine recup
!!***

!!****f* m_qtlmap_genotype/write_typ
!!  NAME
!!    write_typ
!!  DESCRIPTION
!!
!!  INPUTS
!!   file_name : name of the ouput file
!!  SOURCE
      subroutine write_typ(dataset,file_name)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        character(len=*),intent(in)     :: file_name
        integer          :: mes,i,j,c,lk
        type(GENOTYPE_BASE) , pointer :: dga
        dga => dataset%genoAnimal

        open(1111,file=file_name)
        write (1111,*) ((trim(dataset%map%mark(c,lk))//' ',lk=1,dataset%map%nmk(c)),c=1,dataset%map%nchr)

        do mes=1,dga%nmes
         write (1111,*) trim(dga%numero(mes))//' ',&
         (((trim(get_pheno(dga,c,dga%pheno(c,i,mes,j)))//' ',j=1,2),i=1,dataset%map%nmk(c)),c=1,dataset%map%nchr)
        end do

        close (1111)
      end subroutine write_typ
!!***

    !!****f* m_qtlmap_genotype/set_estfem
!!  NAME
!!    set_estfem
!!  DESCRIPTION
!!   initialize the array estfem.
!!  INPUTS
!!   ndmin : minimum number of progenies to build full sib family
!!  SOURCE
      subroutine set_estfem(dataset)
        type(QTLMAP_DATASET)       ,intent(inout)            :: dataset
        integer         :: ip,jm
        type(GENOTYPE_BASE) , pointer :: dga
        type(GENEALOGY_BASE) , pointer :: dg

        dg => dataset%genea
        dga => dataset%genoAnimal

        if ( associated(dga%estfem)) deallocate ( dga%estfem )
        allocate (dga%estfem(dg%nm) )

        dga%estfem=.false.

        do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          if((dg%ndm(jm+1)-dg%ndm(jm)).ge.dataset%params%ndmin) dga%estfem(dg%repfem(jm))=.true.
        end do
      end do
      end subroutine set_estfem
!!***

!!****f* m_qtlmap_genotype/check_HWE
!!  NAME
!!    check_HWE
!!  DESCRIPTION
!!   check_Hardy Weinberg : check the equilibrium of marker transmission within each family
!!  INPUTS
!!   ndmin : minimum number of progenies to build full sib family
!!  SOURCE
    subroutine check_HWE(dataset)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset

       integer :: ll,ip,jm,kd,nb_hete_prog,nb_known_prog
       integer :: c,borne(2)
       type(GENEALOGY_BASE) , pointer :: dg
       type(GENOTYPE_BASE) , pointer :: dga
       dga => dataset%genoAnimal
       dg => dataset%genea

       write (nficout,*)
       write (nficout,*) "** Check the equilibrium of marker transmission within each family **"
       write (nficout,*)

       do c=1,dataset%map%nchr
       !for each marker
        do ll=1,dataset%map%nmk(c)
         if ( dga%nall(c,ll) > 2 ) cycle
         do ip=1,dg%np ! for each sire

           if( dga%pheno(c,ll,dga%correp(ip),1) /= dga%pheno(c,ll,dga%correp(ip),2)) then
            nb_known_prog=0
            nb_hete_prog=0

            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
                      do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
                if ( dga%corred(kd) /= INT_NOT_DEFINED ) then
                  if ( dga%pheno(c,ll,dga%corred(kd),1) /= dataset%genoAnimal%nmanque .and.&
                   dga%pheno(c,ll,dga%corred(kd),2) /= dataset%genoAnimal%nmanque) then
                    nb_known_prog=nb_known_prog+1
                    if (dga%pheno(c,ll,dga%corred(kd),1) /= dga%pheno(c,ll,dga%corred(kd),2)) nb_hete_prog=nb_hete_prog+1
                  end if
                end if
              end do ! end progeny
            end do  ! end dam

            call confbin(nb_known_prog,dataset%params%pseuilHWE,borne)

            if(nb_hete_prog < borne(1) .or. nb_hete_prog > borne(2)) then
              write(nficout,*)'Marker ['//trim(dataset%map%mark(c,ll))//'] for sire :'//trim(dg%pere(ip))//&
                 ' not in HWE : ',nb_hete_prog, ' heterozygous progeny amongst ',nb_known_prog
!              call log_mess('Marker ['//trim(mark(c,ll))//'] for sire :'//trim(pere(ip))//&
!                 ' not in HWE ... should be excluded of the analysis' ,WARNING_DEF)

            end if
          end if
         end do ! end sire
       end do ! end marker
      end do ! end chr

      call log_mess('check_HWE ok',INFO_DEF)
     end subroutine
!!***

!!****f* m_qtlmap_genotype/confbin
!!  NAME
!!    confbin
!!  DESCRIPTION
!!
!!  INPUTS
!!   n     :
!!   p     :
!!   borne :
!!
!!  SOURCE
     subroutine confbin(n,p,borne)
       real(kind=dp) :: p,c,seuil,tc
       integer       :: n,m,l,borne(2)

       seuil=log(p/2.d0)+dble(n) *log(2.d0)
       tc=1.d0
       c=1.d0
       m=n
       l=1

       do while (log(tc) < seuil)
         c=c *dble(m) / dble(l)
         tc=tc+c
         l=l+1.d0
         m=m-1.d0
       end do

       borne(1)=n-m
       borne(2)=m

       end subroutine confbin
!!***

       subroutine create_dataset_genotype(dataset,occurences,array_sample,newdataset)
        type(QTLMAP_DATASET)             ,intent(in)            :: dataset
        integer     ,intent(in)  , dimension(dataset%genea%nd)  :: occurences
        integer             ,dimension(:),intent(in)            :: array_sample
        type(QTLMAP_DATASET)             ,intent(inout)         :: newdataset

        integer             :: c,i,ik,ll,idd,io
        character(len=LEN_DEF), dimension(:,:,:) , allocatable   :: temp_marker_list
        character(len=LEN_DEF), dimension(:,:) , allocatable     :: temp_marker

        type(WORK_GENOTYPE)            :: work
        type(GENOTYPE_BASE) , pointer :: dga
        character(len=LEN_W)   :: buf

        dga => dataset%genoAnimal
        call allocate_vector_(newdataset)

        !on recupere les genotypes lu dans le fichiers de base
        newdataset%genoAnimal%nmes = dataset%genoAnimal%nmes

        !initialisation pour les redondances d'individus
        do i=1,dataset%genea%nd
         if ( occurences(i) > 1 ) then
          newdataset%genoAnimal%nmes = newdataset%genoAnimal%nmes + occurences(i) - 1
         end if
        end do

        !on recupere toute la map
        work%nmarker_in_genotype_file = sum(dataset%map%nmk(:))

        allocate (temp_marker_list(newdataset%genoAnimal%nmes, &
                                   work%nmarker_in_genotype_file,2))

        allocate (newdataset%genoAnimal%numero(newdataset%genoAnimal%nmes))

        allocate (work%index_position_marker(work%nmarker_in_genotype_file))
        allocate (work%index_chr(work%nmarker_in_genotype_file))

        !On prend les marqueurs de l analyse (dans l 'ordre)
        ik=0
        do c=1,dataset%map%nchr
          do ll=1,dataset%map%nmk(c)
             ik=ik+1
             work%index_chr(ik)=c
             work%index_position_marker(ik)=ik
          end do
        end do

        buf=trim(get_pheno(dga,1,dga%nmanque))
        call init_pheno(newdataset%genoAnimal,dataset%map%nchr)
        do c=1,dataset%map%nchr
         newdataset%genoAnimal%nmanque = set_pheno(newdataset%genoAnimal,dataset%map,c,buf)
     !    print *,"nmanque:",dataset%genoAnimal%nmanque
        end do
        !idd=0

        !les repro
        do i=1,dataset%genoAnimal%nmes
         !idd=idd+1
         newdataset%genoAnimal%numero(i)=dga%numero(i)
         ik=0
         do c=1,dataset%map%nchr
          do ll=1,dataset%map%nmk(c)
            ik=ik+1
            temp_marker_list(i,ik,1)=trim(get_pheno(dga,c,dga%pheno(c,ll,i,1)))
            temp_marker_list(i,ik,2)=trim(get_pheno(dga,c,dga%pheno(c,ll,i,2)))
          end do
         end do
        end do


        allocate (temp_marker(work%nmarker_in_genotype_file,2))
        idd=dataset%genoAnimal%nmes
        !les F2 en plus....
        do i=1,newdataset%genea%nd
         if ( occurences(i) > 1 ) then
          ik=0
          do c=1,dataset%map%nchr
            do ll=1,dataset%map%nmk(c)
             ik=ik+1
             temp_marker(ik,1)=trim(get_pheno(dga,c,dga%pheno(c,ll,dga%corred(i),1)))
             temp_marker(ik,2)=trim(get_pheno(dga,c,dga%pheno(c,ll,dga%corred(i),2)))
            end do ! c
           end do ! ll
          do io=2,occurences(i)
           idd=idd+1
           !correspondance
           newdataset%genoAnimal%numero(idd)=trim(dga%numero(dga%corred(i)))//"_"//trim(str(io))
           temp_marker_list(idd,:,1)=temp_marker(:,1)
           temp_marker_list(idd,:,2)=temp_marker(:,2)
           end do !io
          end if
        end do !i
        deallocate (temp_marker)

      !  print *,"NMES:",newdataset%genoAnimal%nmes," IDD:",IDD
        call work%set_pheno_structure(newdataset,newdataset%genoAnimal%nmes,temp_marker_list)

        deallocate (temp_marker_list)

        call set_corresponding_vector(newdataset)
        call set_allele_info_vector(newdataset)
        call set_estfem(newdataset)

!        print *,"OLD"
!        print *,"==="
!        call dataset%genea%print(6)
!        call dataset%genoAnimal%print(6)
!        print *,"NEW"
!        print *,"==="
!        call newdataset%genea%print(6)
!        call newdataset%genoAnimal%print(6)

!        print *, trim(get_pheno(dga,1,dga%pheno(1,5,dga%corred(8),1))),",",&
!        trim(get_pheno(dga,1,dga%pheno(1,5,dga%corred(8),2))),dga%pheno(1,5,dga%corred(8),2)
!        dga=>newdataset%genoAnimal
!        print *, trim(get_pheno(dga,1,dga%pheno(1,5,dga%corred(8),1))),",",&
!         trim(get_pheno(dga,1,dga%pheno(1,5,dga%corred(8),2))),dga%pheno(1,5,dga%corred(8),2)
!        stop



       end subroutine create_dataset_genotype

end module m_qtlmap_genotype
