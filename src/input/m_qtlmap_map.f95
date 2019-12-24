!!****m* INPUT/m_qtlmap_map
!!  NAME
!!    m_qtlmap_map -- Map routines.
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
module m_qtlmap_map
  !Internal parameter
  use m_qtlmap_base
  use m_qtlmap_log
  use m_qtlmap_types
  use m_qtlmap_output_handler

  implicit none
  save

!!****d* m_qtlmap_map/MAX_SIZE_MAP
!!  NAME
!!   MAX_SIZE_MAP
!!  DESCRIPTION
!!   Maximum size in Morgan allowed
!!***
      real (kind=dp) ,   parameter, private                       :: MAX_SIZE_MAP = 20.d0


      ! Name of the chromosome
   character(len=LEN_DEF) , dimension (:), pointer,public        :: ch    => null()
!!****v* m_qtlmap_map/mark0
!!  NAME
!!   mark0
!!  DIMENSIONS
!!   number of marker read in the map file
!!  DESCRIPTION
!!   maker name list read from the map file. buffer array
!!***
      character(len=LEN_DEF) , dimension (:), allocatable,private :: mark0
!!****v* m_qtlmap_map/posi0
!!  NAME
!!   posi0
!!  DIMENSIONS
!!   number of marker read in the map file
!!  DESCRIPTION
!!   average map list read from the map file. buffer array
!!***
      real (kind=dp), dimension (:), allocatable,private          :: posi0
!!****v* m_qtlmap_map/posim0
!!  NAME
!!   posim0
!!  DIMENSIONS
!!   number of marker read in the map file
!!  DESCRIPTION
!!   male map list read from the map file. buffer array
!!***
      real (kind=dp), dimension (:), allocatable,private          :: posim0
!!****v* m_qtlmap_map/posif0
!!  NAME
!!   posif0
!!  DIMENSIONS
!!   number of marker read in the map file
!!  DESCRIPTION
!!   female map list read from the map file. buffer array
!!***
      real (kind=dp), dimension (:), allocatable,public           :: posif0
!!****v* m_qtlmap_map/mselec
!!  NAME
!!   mselec
!!  DIMENSIONS
!!   number of marker read in the map file
!!  DESCRIPTION
!!   select column value list read from the map file. buffer array
!!***
      integer , dimension (:), allocatable,private                :: mselec

  public :: read_map
  public :: sim_carte
  public :: write_map

  CONTAINS
!!****f* m_qtlmap_map/read_map
!!  NAME
!!    read_map
!!  DESCRIPTION
!!   read the map file and fill buffer arrays mark0,posi0,posim0,posif0,mselec.
!!   This file gives the locations of the markers on the chromosome(s). Each line corresponds to a single marker, and gives (order to be followed) :
!!     * marker name (alphanumerique) ;
!!     * name of the chromosome carrying the marker  (alphanumerique) ;
!!     * marker position of the marker on the average map (in Morgan) ;
!!     * marker position of the marker on the male map (in Morgan) ;
!!     * marker position of the marker on the female map (in Morgan) ;
!!     * inclusion key (=1 if the marker has to be included in the analysis, 0 if not)
!!
!!  NOTES
!!  SOURCE
           subroutine read_map(dataset,allChr)
             type(QTLMAP_DATASET)       ,intent(inout)       :: dataset
             logical                    ,intent(in),optional :: allChr

             integer                         :: ios = 56
             integer                         :: eof,err
             character(len=LEN_BUFFER_WORD)  :: word_token_char
             character(len=LEN_DEF)          :: word_token
             character(len=LEN_LINE)         :: line_read
             integer                         :: dimArray = 0
             integer                         :: alloc_stat,l0,nmark,i,ich,j
             logical                         :: is_ok
             character(len=LEN_DEF) , dimension(:),allocatable  :: chromo_temp
             logical :: check
             character(len=LEN_DEF) :: buf

             call log_mess('SUBROUTINE read_map',DEBUG_DEF)
             call log_mess('reading map file...',INFO_DEF)

             open(ios,file=dataset%params%get_file_val(K_MAP),action="read",iostat=eof,status="old")
             if ( eof /= 0 ) then
              call stop_application("Can not find the map file:"//trim(dataset%params%get_file_val(K_MAP)))
             end if


             ! pour savoir le nombre de ligne
             ! read(ios,*,iostat=eof) line_read
             eof = 0
             do while ( eof == 0 )
                read(ios,*,iostat=eof) word_token_char

                if ( trim(word_token_char) /= "" .and. eof == 0 ) then
                    dimArray = dimArray+1
               endif
             end do
             call log_mess('number of line detected in the map file:'//trim(str(dimArray)),DEBUG_DEF)

             ! Allocate array
             ALLOCATE (mark0(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             ALLOCATE (ch(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             ALLOCATE (posi0(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             ALLOCATE (posim0(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             ALLOCATE (posif0(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             ALLOCATE (mselec(dimArray), stat = alloc_stat)
             CALL check_allocate(alloc_stat)

             !!go to the 1rst line
             rewind(ios)

             nmark = 0
             l0 = 0
             eof = 0
             do while (nmark < dimArray .and. (eof == 0) )
              l0 = l0+1
              read(ios,*,IOSTAT=eof) mark0(nmark+1),ch(nmark+1),posi0(nmark+1),posim0(nmark+1),posif0(nmark+1),mselec(nmark+1)
              if ( eof /= 0 ) THEN
                   cycle
              end if
              nmark = nmark + 1
             end do

             !force a garder tous les chromosomes (utiliser pour la verification des genotypages)
             if (present(allChr)) then
               if (allChr .and. nmark>0) then
                   allocate (chromo_temp(100))
                   ich=0
                   do i=1,nmark
                     check=.false.
                     do j=1,ich
                       if ( ch(i) == chromo_temp(j)) then
                         check = .true.
                         exit
                       end if
                     end do
                     if ( .not. check ) then
                       ich=ich+1
                       chromo_temp(ich)=ch(i)
                     end if
                   end do

                   dataset%map%nchr=ich
                   deallocate (dataset%map%chromo)
                   allocate (dataset%map%chromo(dataset%map%nchr))
                   dataset%map%chromo = chromo_temp(:dataset%map%nchr)
                   call log_mess("manage all chromosomes number:"//trim(str(dataset%map%nchr)),INFO_DEF)
                   deallocate(chromo_temp)
                   deallocate (dataset%map%nmk)
                   allocate (dataset%map%nmk(dataset%map%nchr))

                   buf = get_pheno(dataset%genoAnimal,1,dataset%genoAnimal%nmanque)
                   call end_pheno(dataset%genoAnimal)
                   call init_pheno(dataset%genoAnimal,dataset%map%nchr)
                   dataset%genoAnimal%nmanque = set_pheno(dataset%genoAnimal,dataset%map,1,buf)
               end if
             end if


             CALL INIT_INTERNAL_MAP_STRUCTURE(dataset,size(mark0))
             call check_map(dataset)
             call dataset%map%set_absi()
             call log_mess('END SUBROUTINE read_map',DEBUG_DEF)
           end subroutine read_map
!!***

!!****f* m_qtlmap_map/INIT_INTERNAL_MAP_STRUCTURE
!!  NAME
!!    INIT_INTERNAL_MAP_STRUCTURE
!!  DESCRIPTION
!!    initialize the persitents map data array :mark,posi,posim,posif,rm,rf
!!     * select only marker selected
!!     * select only marker positioned on the chromosome list : chromo
!!  INPUTS
!!      nb_marker : number of marker defined in the map file
!!
!!  NOTES
!!  SOURCE
          subroutine INIT_INTERNAL_MAP_STRUCTURE(dataset,nb_marker)
            type(QTLMAP_DATASET)       ,intent(inout)       :: dataset
            integer, intent(in)              :: nb_marker
            integer                          :: i,j,k,c,l1,l2,nb_marker_used,im(dataset%map%nchr)
            integer                          :: alloc_stat,max_mark
            real (kind=dp)                   :: dm,df,p,pf,pm
            character(len=LEN_DEF)             :: t
            logical                          :: sort
            integer,dimension(dataset%map%nchr)                                :: nb_marker_by_ch
            integer,dimension(nb_marker)                           :: ind_chromo
            call log_mess('SUBROUTINE INIT_INTERNAL_MAP_STRUCTURE',DEBUG_DEF)

            nb_marker_by_ch=0
            ind_chromo=0
            do i=1,nb_marker
              !manage only marker selectionned
              if ( mselec(i) /= 1 ) cycle
              ind_chromo(i)=chromo_is_select(dataset,ch(i))
              if ( ind_chromo(i) <= 0 ) cycle
              nb_marker_by_ch(ind_chromo(i))=nb_marker_by_ch(ind_chromo(i))+1
            end do

            !find the greatest to allocate array to the maximum
            max_mark=0
            max_mark=maxval(nb_marker_by_ch)

            ALLOCATE (dataset%map%mark(dataset%map%nchr,max_mark), stat = alloc_stat)
            CALL check_allocate(alloc_stat)

            ALLOCATE (dataset%map%posi(dataset%map%nchr,max_mark), stat = alloc_stat)
            CALL check_allocate(alloc_stat)

            ALLOCATE (dataset%map%posif(dataset%map%nchr,max_mark), stat = alloc_stat)
            CALL check_allocate(alloc_stat)

            ALLOCATE (dataset%map%posim(dataset%map%nchr,max_mark), stat = alloc_stat)
            CALL check_allocate(alloc_stat)
            dataset%map%nmk(:)=0
            do c=1,dataset%map%nchr
              dataset%map%nmk(c)=nb_marker_by_ch(c)
            end do

            im=0
            do i=1,nb_marker

               !manage only marker selectionned
              if ( mselec(i) /= 1 ) cycle
              ind_chromo(i)=chromo_is_select(dataset,ch(i))
              if ( ind_chromo(i) <= 0 ) cycle
              im(ind_chromo(i))=im(ind_chromo(i))+1
              !call log_mess('Marker '//trim(mark0(i))//' is selectionned...',VERBOSE_DEF)

              dataset%map%mark(ind_chromo(i),im(ind_chromo(i)))=mark0(i)
              dataset%map%posi(ind_chromo(i),im(ind_chromo(i)))=posi0(i)
              dataset%map%posif(ind_chromo(i),im(ind_chromo(i)))=posif0(i)
              dataset%map%posim(ind_chromo(i),im(ind_chromo(i)))=posim0(i)
            end do

            deallocate(mark0)
            deallocate(posi0)
            deallocate(posif0)
            deallocate(posim0)
            deallocate(mselec)
!*******************************
! Sort marker
           do c=1,dataset%map%nchr
            sort = .true.
            do while ( sort )
              sort = .false.

              if ( dataset%map%nmk(c) <= 0 ) then
                call stop_application("None marker for the chromosome ["//trim(dataset%map%chromo(c))//"] are selectionned.")
              end if

              do i=1,dataset%map%nmk(c)-1
                if ( dataset%map%posi(c,i) > dataset%map%posi(c,i+1)) then
                    sort = .true.
                    t = dataset%map%mark(c,i)
                    p = dataset%map%posi(c,i)
                    pf= dataset%map%posif(c,i)
                    pm= dataset%map%posim(c,i)
                    dataset%map%mark(c,i) = dataset%map%mark(c,i+1)
                    dataset%map%posi(c,i) = dataset%map%posi(c,i+1)
                    dataset%map%posif(c,i)= dataset%map%posif(c,i+1)
                    dataset%map%posim(c,i)= dataset%map%posim(c,i+1)
                    dataset%map%mark(c,i+1) = t
                    dataset%map%posi(c,i+1) = p
                    dataset%map%posif(c,i+1)= pf
                    dataset%map%posim(c,i+1)= pm
                end if
              end do
             end do
           end do


             if (size(dataset%map%mark) == 0) then
               call stop_application('None marker is selectionned. '// &
               'chromosome selectionned in analyse file:')
             endif



            allocate (dataset%map%rm(dataset%map%nchr,max_mark,max_mark), stat = alloc_stat)
            call check_allocate(alloc_stat)

            allocate (dataset%map%rf(dataset%map%nchr,max_mark,max_mark), stat = alloc_stat)
            call check_allocate(alloc_stat)

           do c=1,dataset%map%nchr
            do l1=1,dataset%map%nmk(c)-1
               do l2=l1+1,dataset%map%nmk(c)
                 dm=dataset%map%posim(c,l2)-dataset%map%posim(c,l1)
                 df=dataset%map%posif(c,l2)-dataset%map%posif(c,l1)
                 dataset%map%rm(c,l1,l2)=xaldane(dm)
                 dataset%map%rf(c,l1,l2)=xaldane(df)
               end do
             end do
           end do
           call log_mess(' END SUBROUTINE INIT_INTERNAL_MAP_STRUCTURE',DEBUG_DEF)
       end subroutine INIT_INTERNAL_MAP_STRUCTURE
!!***

!!****f* m_qtlmap_map/check_map
!!  NAME
!!    check_map
!!  DESCRIPTION
!!    check integrity of information
!!     * The markers can not be overlap
!!     * Marker have to ordered
!!     * size map have a maximum size (MAX_SIZE_MAP)
!!
!!  NOTES
!!  SOURCE
       subroutine check_map(dataset)
            type(QTLMAP_DATASET)       ,intent(in)       :: dataset
            integer         :: i,c
            real (kind=dp)  :: temp
            call log_mess('checking map......',VERBOSE_DEF)

            do c=1,dataset%map%nchr
               if ( (dataset%map%posi(c,dataset%map%nmk(c))-dataset%map%posi(c,1)) < dataset%map%get_long_step_morgan() ) then
                 call stop_application("CHR ["//trim(dataset%map%chromo(c))//&
                 "]The step is biggest than the consensus map [sizemap:"//&
                trim(str(dataset%map%posi(c,dataset%map%nmk(c))-dataset%map%posi(c,1)))//&
                "] [step:"//trim(str(dataset%map%get_long_step_morgan()))//"]"//&
                ". Bad definition of average map : Fist marker :"//trim(str(dataset%map%posi(c,1)))//&
                  " Last marker :"//trim(str(dataset%map%posi(c,dataset%map%nmk(c)))))
               end if

               if ( (dataset%map%posim(c,dataset%map%nmk(c))-dataset%map%posim(c,1)) < dataset%map%get_long_step_morgan()  ) then
                call stop_application("CHR ["//trim(dataset%map%chromo(c))//&
                 "]The step is biggest than the male map [sizemap:"//&
                trim(str(dataset%map%posim(c,dataset%map%nmk(c))-dataset%map%posim(c,1)))//&
                "] [step:"//trim(str(dataset%map%get_long_step_morgan()))//"]"//&
                ". Bad definition of male map : Fist marker :"//trim(str(dataset%map%posim(c,1)))//&
                 " Last marker :"//trim(str(dataset%map%posim(c,dataset%map%nmk(c)))))
               end if

               if ( (dataset%map%posif(c,dataset%map%nmk(c))-dataset%map%posif(c,1)) < dataset%map%get_long_step_morgan()  ) then
                call stop_application("CHR ["//trim(dataset%map%chromo(c))//&
                "]The step is biggest than the female map [sizemap:"//&
                trim(str(dataset%map%posif(c,dataset%map%nmk(c))-dataset%map%posif(c,1)))//"] [step:"&
                //trim(str(dataset%map%get_long_step_morgan()))//"]"//&
                ". Bad definition of female map : Fist marker :"//trim(str(dataset%map%posif(c,1)))//&
                 " Last marker :"//trim(str(dataset%map%posif(c,dataset%map%nmk(c)))))
               end if


               do i=2,dataset%map%nmk(c)
                 if (dataset%map%posi(c,i-1)>=dataset%map%posi(c,i)) then
                   call stop_application('CHR ['//trim(dataset%map%chromo(c))//'] Marker ['//trim(dataset%map%mark(c,i-1))// &
                   '] defined in map file is greater or equal than Marker ['&
                   //trim(dataset%map%mark(c,i))//'] for average map')
                 end if
               end do

               do i=2,dataset%map%nmk(c)
                 if (dataset%map%posim(c,i-1)>=dataset%map%posim(c,i)) then
                   call stop_application('CHR ['//trim(dataset%map%chromo(c))//'] Marker ['//trim(dataset%map%mark(c,i-1))// &
                   '] defined in map file is greater or equal than Marker ['&
                   //trim(dataset%map%mark(c,i))//'] for male map')
                 end if
               end do

               do i=2,dataset%map%nmk(c)
                 if (dataset%map%posif(c,i-1)>=dataset%map%posif(c,i)) then
                   call stop_application('CHR ['//dataset%map%chromo(c)//'] Marker ['//trim(dataset%map%mark(c,i-1))// &
                   '] defined in map file is greater or equal than Marker ['&
                   //trim(dataset%map%mark(c,i))//'] for female map')
                 end if
               end do
             !!checkin the map size
             if ( dataset%map%posif(c,dataset%map%nmk(c))-dataset%map%posif(c,1) > MAX_SIZE_MAP ) then
                  temp = dataset%map%posif(c,dataset%map%nmk(c))-dataset%map%posif(c,1)
                  call log_mess('The map is very large : ['//trim(str(temp))//&
                  ']. The map have to be defined in Morgan!',WARNING_DEF)
             end if

            end do

            call log_mess('map is checked.',VERBOSE_DEF)
       end subroutine check_map
!!***

!!****f* m_qtlmap_map/sim_carte
!!  NAME
!!    sim_carte
!!  DESCRIPTION
!!
!!  INPUTS
!!   c     : index of chromosome to fill
!!  dens   : density, number of marker / Morgan
!!  nalle  : number of allele by marker
!!
!!  NOTES
!!   Simulation d'une carte genetique et des caracteristiques
!!   des marqueurs genetiques Appele par lect_carte
!!  SOURCE
      subroutine sim_carte(dataset,c,dens,nalle)
      type(QTLMAP_DATASET)       ,intent(inout) :: dataset
      integer        , intent(in)      :: c ! chromosome
      real (kind=dp) , intent(in)      :: dens
      integer        , intent(in)      :: nalle
!
      integer        :: il,l1,l2,i
      real (kind=dp) :: dm,df
      type(GENOTYPE_BASE)  , pointer :: dga

      dga => dataset%genoAnimal

!***********************************************************************
!                      Construction de la carte
!  - deduit nb marqueurs de densite et taille chr
!  - attribue meme position des marqueurs sur cartes male et femelle
!***********************************************************************
!
!
      call log_mess('simulation of map...',VERBOSE_DEF)

      dataset%map%posi(c,1)=0
      dataset%map%posim(c,1)=0
      dataset%map%posif(c,1)=0
      do il=2,dataset%map%nmk(c)
        dataset%map%posi(c,il)=dataset%map%posi(c,il-1)+dens
        dataset%map%posim(c,il)=dataset%map%posim(c,il-1)+dens
        dataset%map%posif(c,il)=dataset%map%posif(c,il-1)+dens
      end do

!
! Initialisation des taux de recombinaison entre marqueurs
        do l1=1,dataset%map%nmk(c)-1
         do l2=l1+1,dataset%map%nmk(c)
          dm=dataset%map%posim(c,l2)-dataset%map%posim(c,l1)
          df=dataset%map%posif(c,l2)-dataset%map%posif(c,l1)
          dataset%map%rm(c,l1,l2)=xaldane(dm)
          dataset%map%rf(c,l1,l2)=xaldane(df)
        end do
       end do
!


!************************************************************************
!                    Caract�ristiques des marqueurs: tous identiques
!  - nom des marqueurs : Mark l
!  - nom des all�les : 1,2...
!
!***********************************************************************
!
      dga%nall=nalle
      !manque=manq
      !typ=1

      do il=1,dataset%map%nmk(c)
        do i=1,dga%nall(c,il)
            dga%alleles(c,il,i)=trim(str(i))
            dga%pc_all(c,il,i)=1.d2/dble(dga%nall(c,il))
        end do
      end do

      end subroutine sim_carte
!!***

!!****f* m_qtlmap_map/write_map
!!  NAME
!!    write_map
!!  DESCRIPTION
!!    write a map file in the qtlmap format
!!  INPUTS
!!   file_name : name of the output file
!!
!!  NOTES
!!  SOURCE
      subroutine write_map(map,file_name)
       type (MAP_BASE) ,intent(in)     :: map
       character(len=*),intent(in)     :: file_name
       integer         :: il,c
       call log_mess('writing map file :['//trim(file_name)//']', INFO_DEF)
       open(1,file=file_name)
       do c=1,map%nchr
         do il=1,map%nmk(c)
           write(1,*) map%mark(c,il),c,' ',map%posi(c,il),map%posim(c,il),map%posif(c,il),' 1'
         end do
       end do
       close(1)

      end subroutine write_map
!!***

!!****f* m_qtlmap_map/chromo_is_select
!!  NAME
!!    chromo_is_select
!!  DESCRIPTION
!!    compare the name ch with all name records in the chromo array.
!!  INPUTS
!!   ch : index of chromosome
!!
!!  RETURN
!!   true if the chromosome is selected
!!  NOTES
!!  SOURCE
      function chromo_is_select(dataset,ch) result(ind)
        type(QTLMAP_DATASET)       ,intent(inout) :: dataset
        character(len=LEN_DEF)          :: ch
        integer :: ind,i

        ind = 0
        do i=1,dataset%map%nchr
           if ( trim(dataset%map%chromo(i)) == trim(ch) ) then
              ind = i
              return
           end if
        end do
        return
      end function chromo_is_select
!!***

!!    Get index of markers (start and end) to print offspring haplotypes
!!   the binary qtlmap2cartha.
      subroutine set_haplotype_offspring_context(dataset,C,M1,M2,ok,namefile)
          type(QTLMAP_DATASET)       ,intent(in)    :: dataset
          integer   , intent(out) :: M1,M2,C
          logical   , intent(out) :: ok
          character(len=LENGTH_MAX_FILE) , intent(out) :: namefile
          character(len=LEN_W)   :: buf
          integer :: C1=0,C2=0,ik
          logical :: found

          ok=.false.
          namefile=""
          if ( dataset%params%key_exist(K_PHASES_OFFSPRING) ) then
           call dataset%params%get_string_val(K_PHASES_OFFSPRING,namefile)
           if ( dataset%params%key_exist(K_PHASES_OFFSPRING_MARK_START) ) then
             call dataset%params%get_string_val(K_PHASES_OFFSPRING_MARK_START,buf)
             found=.false.
             do c=1,dataset%map%nchr
               do ik=1,dataset%map%nmk(c)
                if ( trim(buf) == dataset%map%mark(c,ik) ) then
                   M1 = ik
                   C1 = c
                   found=.true.
                   exit
                end if
                end do
                if (found) exit
             end do
             if ( .not. found ) then
               call stop_application("Key ["//trim(K_PHASES_OFFSPRING_MARK_START)//&
                "] can not found the marker:"//trim(buf))
             end if
          else
            M1 = 1
            C1 = 1
          end if
          if ( dataset%params%key_exist(K_PHASES_OFFSPRING_MARK_END) ) then
             call dataset%params%get_string_val(K_PHASES_OFFSPRING_MARK_END,buf)
             found=.false.
             do c=1,dataset%map%nchr
               do ik=1,dataset%map%nmk(c)
                if ( trim(buf) == dataset%map%mark(c,ik) ) then
                   M2 = ik
                   C2 = c
                   found=.true.
                   exit
                end if
                end do
                if (found) exit
             end do
             if ( .not. found ) then
               call stop_application("Key ["//K_PHASES_OFFSPRING_MARK_END//"] can not found the marker:"//trim(buf))
             end if
           else
               M2 = dataset%map%nmk(C1)
               C2 = C1
           end if

          if ( C1 /= C2 ) then
             call stop_application("Key ["//K_PHASES_OFFSPRING_MARK_START//"] and key ["//&
               K_PHASES_OFFSPRING_MARK_END//"] are defined on different chromosomes.")
          end if
          C = C1
          ok = .true.
         end if

      end subroutine set_haplotype_offspring_context


END MODULE m_qtlmap_map
