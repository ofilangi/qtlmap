module m_qtlmap_type_dataset
    use m_qtlmap_constant
    use m_qtlmap_type_genealogy
    use m_qtlmap_type_map
    use m_qtlmap_type_phenotype
    use m_qtlmap_type_genotype
    use m_qtlmap_type_cli
    use m_qtlmap_type_parameter
    use m_qtlmap_type_dataset_trio
    implicit none

   type ,public :: QTLMAP_DATASET

     type(MAP_BASE)       , pointer    :: map          => null()
     type(GENEALOGY_BASE) , pointer    :: genea        => null()
     type(GENEALOGY_RACE) , pointer    :: geneaRace    => null()
     type(GENOTYPE_BASE)  , pointer    :: genoAnimal   => null()
     type(PHENOTYPE_BASE) , pointer    :: phenoAnimal  => null()
     type(DATAMODEL_BASE) , pointer    :: phenoModel   => null()

     type(PARAMETER_BASE) ,pointer     :: params       => null()

     type(QTLMAP_CLI),pointer           :: cli => null()

     !New structure used for multitraits analysis (These Mohamed Kileh wais)
     type (DATASET_QTLMAP_TRIO), public , pointer :: datasetUser

    contains

      procedure ,public :: set     => set_qtlmap_dataset
      procedure ,public :: set_parameters_qtlmap_dataset !
      procedure ,public :: copy    => copy_qtlmap_dataset
      procedure ,public :: release => release_qtlmap_dataset
      procedure ,public :: get_list_dam_estime ! donne la liste des indexes des meres estimables
   end type QTLMAP_DATASET

   contains

   subroutine set_qtlmap_dataset(dataset)
         class(QTLMAP_DATASET)       ,intent(inout) :: dataset

         allocate (dataset%genea)
         allocate (dataset%map)
         allocate (dataset%geneaRace)
         allocate (dataset%genoAnimal)
         allocate (dataset%phenoAnimal)
         allocate (dataset%phenoModel)
         allocate (dataset%params)
         allocate (dataset%cli)
         allocate (dataset%datasetUser)

     end subroutine set_qtlmap_dataset


     subroutine set_parameters_qtlmap_dataset(dataset,p_analyse)
       class(QTLMAP_DATASET)       ,intent(inout) :: dataset
       character(len=LENGTH_MAX_FILE),intent(in) :: p_analyse

       integer :: i,chr
       integer ,parameter :: NSIZE = 30
       character(len=LEN_DEF) ,dimension(NSIZE) :: listKeys,listValues
       integer :: nkeys
       character(len=LEN_W)   :: buf,buf2
       character(len=LEN_S),dimension(100)   :: chromo_t
       logical :: ok
       real(kind=dp) :: pas_temp

       listKeys=''
       listValues=''

       call dataset%cli%cli_get_overloaded_keys(listKeys,listValues,NSIZE,nkeys)
       call dataset%params%set(p_analyse,listKeys,listValues,NSIZE,nkeys)


        call dataset%params%get_string_val(K_STEP,buf)
        pas_temp    =  get_real(buf,ok)
        if ( .not. ok ) call stop_application("bad definition of key ["//trim(K_STEP)//"] expecting int value.")
        call dataset%map%set_base_and_step(buf)

        call dataset%params%get_string_val(K_CHROM,buf2)
        i=1;
        dataset%map%nchr=0

        do while ( i /= 0)
          i=index(buf2,",")
          if ( i /= 0 ) then
            dataset%map%nchr = dataset%map%nchr+1
            chromo_t(dataset%map%nchr) = trim(buf2(:i-1))
            buf2 = buf2(i+1:)
          else
            dataset%map%nchr = dataset%map%nchr+1
            chromo_t(dataset%map%nchr) = trim(buf2)
          end if
        end do

        allocate (dataset%map%chromo(dataset%map%nchr))
        dataset%map%chromo(:dataset%map%nchr)=chromo_t(:dataset%map%nchr)
        ALLOCATE (dataset%map%nmk(dataset%map%nchr))

        call dataset%params%get_string_val(K_UNKNOWN_GENO,buf)
        buf=trim(buf)
        call init_pheno(dataset%genoAnimal,dataset%map%nchr)
        do chr=1,dataset%map%nchr
         dataset%genoAnimal%nmanque =  set_pheno(dataset%genoAnimal,dataset%map,chr,buf)
        end do

        if ( dataset%cli%cli_is_permute() .and. .not. dataset%params%key_exist(K_MODEL) ) then
           call stop_application("key ["//trim(K_MODEL)//&
           "] is not defined. To perform a permutation, you have to defined a model.")
        end if


     end subroutine set_parameters_qtlmap_dataset



     subroutine release_qtlmap_dataset(dataset)
         class(QTLMAP_DATASET)       ,intent(inout) :: dataset

         if (associated(dataset%genea)) then
           call dataset%genea%release()
           deallocate (dataset%genea)
         end if

         if (associated(dataset%map)) then
           call dataset%map%release()
           deallocate (dataset%map)
         end if

         if (associated(dataset%geneaRace)) then
          call dataset%geneaRace%release()
          deallocate (dataset%geneaRace)
         end if

         if (associated(dataset%genoAnimal)) then
          call dataset%genoAnimal%release()
          deallocate (dataset%genoAnimal)
         end if

         if (associated(dataset%phenoAnimal)) then
          call dataset%phenoAnimal%release()
          deallocate (dataset%phenoAnimal)
         end if

         if (associated(dataset%phenoModel)) then
          call dataset%phenoModel%release()
          deallocate (dataset%phenoModel)
         end if

         if (associated(dataset%params)) then
           call dataset%params%release()
           deallocate(dataset%params)
         end if

         if (associated(dataset%cli)) then
           call dataset%cli%release()
           deallocate(dataset%cli)
         end if

         if (associated(dataset%datasetUser)) then
           call dataset%datasetUser%release()
           deallocate(dataset%datasetUser)
         end if

     end subroutine release_qtlmap_dataset

    subroutine copy_qtlmap_dataset(mapDataset,copyDataset)
        class(QTLMAP_DATASET), intent(in)     :: mapDataset
        type(QTLMAP_DATASET), intent(inout)  :: copyDataset

        call mapDataset%map%copy(copyDataset%map)
        call mapDataset%genea%copy(copyDataset%genea)
        call mapDataset%geneaRace%copy(copyDataset%geneaRace)
        call mapDataset%genoAnimal%copy(copyDataset%genoAnimal)
        call mapDataset%phenoAnimal%copy(copyDataset%phenoAnimal)
        call mapDataset%phenoModel%copy(copyDataset%phenoModel)
        call mapDataset%params%copy(copyDataset%params)
        call mapDataset%cli%copy(copyDataset%cli)
        call mapDataset%datasetUser%copy(copyDataset%datasetUser)

    end subroutine copy_qtlmap_dataset


    subroutine get_list_dam_estime(dataset,ic,listmere,nestime)
       class(QTLMAP_DATASET)         ,intent(in)  :: dataset
       integer                       ,intent(in)  :: ic        ! le caractere etudié, si 0 => multicaractere
       integer, dimension(:),pointer ,intent(out) :: listmere
       integer                       ,intent(out) :: nestime

       integer :: jm,countnd,kd
       integer, dimension(:),pointer :: listmere_tps

       nestime = 0
       allocate (listmere_tps(dataset%genea%nm))

       do jm=1,dataset%genea%nm
         countnd=0
         do kd=dataset%genea%ndm(jm)+1,dataset%genea%ndm(jm+1)
          if (ic>0) then
           if ( dataset%phenoAnimal%presentc(ic,kd)) then
            countnd = countnd +1
           end if
          else
           if ( count(dataset%phenoAnimal%presentc(:,kd)) == dataset%phenoModel%ncar) then
            countnd = countnd +1
           end if
          end if

         end do

         !if they are enough progenies with perf so...
         if ( countnd >= dataset%params%ndmin ) then
           nestime=nestime+1
           listmere_tps(nestime) = jm
         end if

       end do

       allocate (listmere(nestime))
       listmere = listmere_tps(:nestime)
       deallocate (listmere_tps)

    end subroutine get_list_dam_estime

end module m_qtlmap_type_dataset
