!Copyright or © or Copr. INRA (2010)
!contributors :
!
!  Pascale Le Roy     <Pascale.LeRoy@rennes.inra.fr>
!  Jean-Michel Elsen  <Jean-Michel.Elsen@toulouse.inra.fr>
!  Helene Gilbert     <Helene.Gilbert@jouy.inra.fr>
!  Carole Moreno      <Carole.Moreno@toulouse.inra.fr>
!  Andres Legarra     <Andres.Legarra@toulouse.inra.fr>
!  Olivier Filangi    <Olivier.Filangi@rennes.inra.fr>
!
!
!QTLMap is a computer program dedicated to the detection
!of QTL from experimental designs in outbred population.
!The statistical techniques used are linkage analysis (LA) and
!linkage disequilibrium linkage analysis (LDLA) using interval mapping.
!Different versions of the LA are proposed from a  quasi Maximum Likelihood
!approach to a fully linear (regression) model.
!The LDLA is a regression approach (Legarra and Fernando, 2009).
!The population may be sets of half-sib families or mixture of
!full- and half- sib families.
!The computations of Phase and Transmission probabilities are optimized to be rapid and as exact as possible.
!QTLMap is able to deal with large numbers of markers (SNP) and traits (eQTL).

!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software.  You can  use,
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info".
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software's author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability.
!
!In this respect, the user's attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software's suitability as regards their
!requirements in conditions enabling the security of their systems and/or
!data to be ensured and,  more generally, to use and operate it in the
!same conditions as regards security.
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.

program qtlmap
#ifdef MPI_ACTIVE
      use mpi
#endif
#ifdef MANAGE_OMP
       use omp_lib
#endif
      use m_qtlmap_log
      use m_qtlmap_constant
      use m_qtlmap_haplotype
      use m_qtlmap_analyse
      use m_qtlmap_genealogy
      use m_qtlmap_genotype
      use m_qtlmap_phenotype
      use m_qtlmap_map
      use m_qtlmap_simulation
      use m_qtlmap_output_handler
      use m_qtlmap_optimization
      use m_qtlmap_math
      use m_qtlmap_calcul_ic

      implicit none

      ! Dec  2012  - QTLMap dataset
      type(QTLMAP_DATASET)       :: dataset
      type(PDD_BUILD)            :: spt


! Divers
      real (kind=dp), dimension (:),allocatable :: xmutanalyse,sigtanalyse
      integer :: nsim,thread_rang = 0 ! by default / identifiant du processus fils mpi
      integer :: insim=0,iseed1,iseed2

      character(len=100)             :: phrase
      character(len=8)               :: date
      character(len=10)              :: time
      character(len=5)               :: zone
      integer, dimension(8)          :: value

#ifdef MPI_ACTIVE
        integer :: code          ! identifiant du processus global mpi
        integer :: nb_procs      ! nombre de processus alloue par l utilisateur
        integer :: nb_slaves=1   ! nombre de processus alloue
        logical :: drapeau=.true.
        logical :: logique
        integer :: intra_comm,inter_comm
        integer(kind = MPI_ADDRESS_KIND) :: nb_procs_total

        call MPI_INIT (code)
        call MPI_COMM_SIZE ( MPI_COMM_WORLD , nb_procs , code)
        call MPI_COMM_RANK (MPI_COMM_WORLD, thread_rang , code)
        if ( thread_rang == 0 ) then
          call log_mess(" ",INFO_DEF)
          call log_mess("   ****************************** ",INFO_DEF)
          call log_mess("   ****           MPI         *** ",INFO_DEF)
          call log_mess("   NB_PROC:        :"//trim(strint(nb_procs)),INFO_DEF )
          call log_mess("   ****************************** ",INFO_DEF)
          call log_mess(" ",INFO_DEF)
        end if
#endif
#ifdef MANAGE_OMP
     !parallelisation sur le parcours du groupe de liaison
     !call omp_set_nested(.true.)

     call log_mess(" ",INFO_DEF)
     call log_mess("   ****************************** ",INFO_DEF)
     call log_mess("   **** Use OpenMP directives *** ",INFO_DEF)
     call log_mess("   OMP_NUM_THREADS   :"//trim(strint(omp_get_max_threads())),INFO_DEF )
     call log_mess("   TO CHANGE  EXPORT OMP_NUM_THREADS=<NB> or SETENV OMP_NUM_THREADS <NB>",INFO_DEF )
     !call log_mess("   OMP_NESTED        :"//trim(strbool(omp_get_nested())),INFO_DEF )
     call log_mess("   ****************************** ",INFO_DEF)
     call log_mess(" ",INFO_DEF)
#endif
      call dataset%set()

      ! Printing Version
      call dataset%cli%cli_print_qtlmap_version()

      !----- PRINT HELP IF OPTION --HELP DETECTED
      if (dataset%cli%cli_is_help_panalyse()) then
         call dataset%params%help_panalyse()
         stop
      end if

       if (dataset%cli%cli_is_help_paramsim()) then
         call help_paramsim
         stop
      end if

      if (dataset%cli%cli_is_help()) then
          call dataset%cli%cli_print_help(.false.)
          stop
      end if

      if (dataset%cli%cli_is_help_all()) then
          call dataset%cli%cli_print_help(.true.)
          stop
      end if

      ! ----- INIT LOG LEVEL
      call init_log(INFO_DEF)

      if (dataset%cli%cli_is_quiet()) call init_log(ERROR_DEF)
      if (dataset%cli%cli_is_verbose()) call init_log(VERBOSE_DEF)
      if (dataset%cli%cli_is_debug()) call init_log(DEBUG_DEF)

      call inrgcm()
      ! Init the random generator
      call DATE_AND_TIME(date, time, zone, value)
      phrase=date//time//zone
      call phrtsd(phrase,iseed1,iseed2)
      call setall(iseed1,iseed2)

      !keep information about command line interface

     !Getting the optimization implementation for minimizing likelihood
      call set_optimization(dataset%cli%cli_get_optimization())

       !----- Give information about library used by the program
      call MATH_QTLMAP_INFO

      ! ------ READING PARAMETER FILE
      call dataset%set_parameters_qtlmap_dataset(dataset%cli%cli_get_parameters_file())

      ! --> PARAMETER PART OF QTLMAP_DATA IS INITIALIZED
#ifdef MANAGE_CUDA
      insim=2
#else
      if ( dataset%cli%cli_is_simulation() .or. dataset%cli%cli_is_nsim() ) then
       insim = 1
      else
       insim = 0
      end if
#endif

      call check_simulation_compatibility(dataset)

      ! ------ PRINTING ANALAYSE AND HAPLOTYPE VERSION DETECTED
      call log_mess('opt_calcul:'//trim(str(dataset%cli%cli_get_analyse())),INFO_DEF)
      call log_mess('opt_version:'//trim(str(dataset%cli%cli_get_haplotype())),INFO_DEF)
      dataset%params%opt_sib = dataset%cli%cli_get_family()
      call log_mess('opt_sib:'//trim(str(dataset%params%opt_sib)),INFO_DEF)

#ifdef MANAGE_CUDA
         call qtlmap_cuda_simul_and_analyse(dataset,dataset%cli%cli_get_nsim(),&
          dataset%cli%cli_get_analyse(),dataset%cli%cli_get_nqtl())
#else
      !------- IF SIMULATION OPTION OR NSIM OPTION IS DETECTED --> SIMULATION
     !print *, 'SIMUL', cli_is_simulation(),cli_is_nsim() 
      if(dataset%cli%cli_is_simulation() .or. dataset%cli%cli_is_nsim())then
         nsim = dataset%cli%cli_get_nsim()
#ifdef MPI_ACTIVE
         nsim = nsim - mod(nsim,nb_procs)
#endif
         ! ----RUNNING SIMULATION
        call qtlmap_simul(dataset,.false.,nsim,dataset%cli%cli_get_analyse(),dataset%cli%cli_get_nqtl())
      else
         ! ----RUNNING ANALYZE
         call qtlmap_analyze(dataset,dataset%cli%cli_get_analyse(),dataset%cli%cli_get_nqtl())
      end if
#endif

#ifdef MPI_ACTIVE
      print *,"PROC ["//trim(str(thread_rang))//"] END ** "
      call MPI_Finalize (code);
#endif

      call dataset%release()

      contains

!***********************************  SUBROUTINE *********************************************************
      subroutine qtlmap_analyze(dataset,opt_calcul,opt_qtl)
        type(QTLMAP_DATASET) ,         intent(inout)      :: dataset
        integer         , intent(in)              :: opt_calcul,opt_qtl

        type(TYPE_LRT_SOLUTION)  , dimension(:,:),allocatable       :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION) ,dimension(:,:), allocatable  :: listincsol

        real (kind=dp), dimension (:),allocatable   :: ih2
        real (kind=dp),dimension (:,:,:),allocatable :: rhoi

        real(kind=dp) :: hdeux
        integer ,dimension(10) :: position_qtl
        real(kind=dp) ,dimension(10) ::effet_qtl

        integer                  :: c,stat,opt_sim,jm,igeno,nsim,i,j,nnnqtl,sizeList,nb_thread
        integer                 :: M1,M2
        character(len=LEN_DEF)   :: simula_t,cc
        character(len=LEN_L) ,dimension(200) :: index_key,values
        character(len=LENGTH_MAX_FILE) :: genin,perfin,perfout,phases_off
        character(len=300) :: name_funct
        
        real (kind=dp) , dimension(:,:),   allocatable   :: countH
        real (kind=dp) , dimension(:,:,:), allocatable   :: markH

        real (kind=dp) , dimension(4,4)   :: effects_tab
        character(len=1024)               :: outf

        logical :: ok
        !! Information of generated by haplotype module
        type(PDD_BUILD)     :: spt

        if (dataset%cli%cli_is_transcriptomic_data() .and. .not. dataset%cli%cli_is_print_all() ) then
           opt_sim = TRANSCIPTOME_ANALYSE
        else
           opt_sim = COMMON_ANALYSE
        end if

        name_funct = get_name_analyse(opt_calcul)
#ifdef MANAGE_OMP
        nb_thread = omp_get_max_threads()
#else
        nb_thread = 1
#endif
        call init_output_handler(dataset,opt_qtl,opt_calcul,name_funct,nb_thread)

        call dataset%params%get_summary_panalyse(index_key,values,c)
        call print_summary_panalyse(c,index_key,values)

        !******** Lecture de la map et de la genealogie en parallele ****

        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP SECTIONS
        !$OMP SECTION
        call read_map(dataset)
        !$OMP SECTION
        call read_genealogy(dataset)
        !$OMP END SECTIONS NOWAIT
        !$OMP END PARALLEL

        !******** Lecture des genotypes et de la genealogie en parallele ****

        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP SECTIONS
        !$OMP SECTION
        call read_genotype(dataset)
        !$OMP SECTION
        call read_model(dataset)
        call fixe_structure_model(dataset)
        call read_traits(dataset)
        !$OMP END SECTIONS NOWAIT
        !$OMP END PARALLEL

        call recup(dataset)
        call check_typage(dataset)
        call log_descriptif_genealogy(dataset)
        call log_marker_description(dataset)
        call print_freqall(dataset)
        call check_HWE(dataset)

        if ( dataset%cli%cli_is_calcul_cd() ) then
          call dataset%datasetUser%calcul_y_cd(dataset%genea,dataset%phenoAnimal,dataset%phenoModel)
          call dataset%datasetUser%calcul_corcd(dataset%genea,dataset%phenoAnimal,dataset%phenoModel)
        end if

        ! Normalizing data and print information about performance trait
        call manage_data(dataset,dataset%cli%cli_is_transcriptomic_data(),.false.,need_normalize_data(opt_calcul))

        if (.not. dataset%cli%cli_is_transcriptomic_data()) then
           call log_descriptif_traits(dataset)
        end if

! **************************************************
! Information sire Heterozygote

      if ( dataset%params%get_file_val(K_INFORMATIVITY) /= '' ) then
       allocate (countH(dataset%genea%np,dataset%map%nchr))
       allocate (markH(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)))
       ! print information in a file : heterozygote marker from given sire/marker/linkage group
       call get_information_informative_marker(dataset,spt,countH,markH)
       call print_informativity_markers(dataset,dataset%params%get_file_val(K_INFORMATIVITY),countH,markH)
       deallocate (countH)
       deallocate (markH)
      end if
! **************************************************


        allocate(lrtsol(dataset%phenoModel%ncar,opt_qtl))
        allocate (listincsol(dataset%phenoModel%ncar,opt_qtl+1))
        allocate (rhoi(dataset%map%nchr,dataset%phenoModel%ncar,dataset%phenoModel%ncar))
        allocate (ih2(dataset%phenoModel%ncar))

        !computing phase,haplotype
        call haplotype(dataset,spt,dataset%cli%cli_get_haplotype())

       ! get the offspring haplotype around the Maximum Likelihood Position
        call haplotype_offspring(dataset,spt,0)


        do c=1,dataset%map%nchr
         do jm=1,dataset%genea%nm
           call log_mess('Mother :'//trim(dataset%genea%mere(jm))//&
            '-->genotypes found [Chr '//trim(dataset%map%chromo(c))//']:'//&
       trim(str(spt%ngenom(c,jm+1)-spt%ngenom(c,jm))),VERBOSE_DEF)
           do igeno=spt%ngenom(c,jm)+1,spt%ngenom(c,jm+1)
                call log_mess('genotype '//trim(str((igeno)))//' : '//trim(str(spt%probg(c,igeno)*100.d0))//'%',VERBOSE_DEF)
           end do
          end do
        end do


        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP SECTIONS
        !$OMP SECTION
        do c=1,dataset%map%nchr
          ! Printing probabilities of transmission in associated file
          call print_pded(dataset,spt,c)
        end do
        !$OMP SECTION
        ! Printing genotype finded for sire/ all genotyp available for dam
        call print_phases(dataset,spt)
        !$OMP END SECTIONS NOWAIT
        !$OMP END PARALLEL

!        effects_tab=0.d0
!        !biq additif
!        effects_tab(1,1)=-30.d0
!        effects_tab(2,1)=-15d0
!        effects_tab(3,1)=15d0
!        effects_tab(4,1)=30.d0
!
!        effects_tab(2,1)=-2d0
!        effects_tab(2,2)=-1.5d0
!        effects_tab(2,3)=-0.5d0
!        effects_tab(2,4)=0d0
!
!        effects_tab(3,1)=0d0
!        effects_tab(3,2)=0.5d0
!        effects_tab(3,3)=1.5d0
!        effects_tab(3,4)=2.d0
!
!        effects_tab(4,1)=1.d0
!        effects_tab(4,2)=1.5d0
!        effects_tab(4,3)=2.5d0
!        effects_tab(4,4)=3d0
!
!        outf="pheno_2qtl_inter2"
!        call sim_perf_interaction(dataset,spt,0.5d0,1,65,1,380,effects_tab,outf)
!        stop

        call log_mess("-- ANALYSIS --  CALCUL="//get_name_analyse(opt_calcul),INFO_DEF)
        !analysis

        call analyse(dataset,spt,opt_calcul,opt_qtl,lrtsol,listincsol,rhoi,opt_sim)

      !  if ( .not. is_multitrait_analysis(opt_calcul) ) then
         ! Dev Confidence Intervals/QTL Mapping
         call computingIC(dataset,lrtsol,listincsol,dataset%cli%cli_get_nqtl())
         call print_confidence_intervals_solution(dataset,dataset%cli%cli_get_nqtl(),&
          lrtsol,is_multitrait_analysis(opt_calcul))
      !  end if

        ! Printing genotype finded for offspring
        call set_haplotype_offspring_context(dataset,C,M1,M2,ok,phases_off)
        if (ok) call print_offspring_phase(dataset,spt,C,M1,M2,lrtsol,opt_qtl,phases_off,&
         is_multitrait_analysis(opt_calcul))


        if (.not. dataset%cli%cli_is_transcriptomic_data()) then
           if (.not. is_multitrait_analysis(opt_calcul)) then
             !print to the console max LRT
              call print_information_max_lrt(dataset,opt_qtl,INFO_DEF,lrtsol)
          else
             !print to the console max LRT
              call print_information_max_lrt_multi(dataset,opt_qtl,INFO_DEF,lrtsol(1,:))
          end if
        end if

       do i=1,size(lrtsol,1)
          do j=1,size(lrtsol,2)
            call lrtsol(i,j)%release()
          end do
        end do

         do i=1,size(listincsol,1)
          do j=1,size(listincsol,2)
            call listincsol(i,j)%release()
          end do
        end do


        deallocate(lrtsol)
        deallocate(rhoi)
        deallocate (ih2)


        call end_output_handler(dataset,opt_qtl)
        call spt%release

      end subroutine qtlmap_analyze




#ifdef MANAGE_CUDA

  subroutine qtlmap_cuda_simul_and_analyse(dataset,nsim,opt_calcul,opt_qtl)
         type(QTLMAP_DATASET)  ,intent(inout)      :: dataset
         integer         , intent(in)              :: nsim
         integer         , intent(in)              :: opt_calcul
         integer         , intent(in)              :: opt_qtl

         !local
         logical                        :: permut_active=.false.
         integer                        :: stat,i,ic,isim,opt_sim,opt_version,j,iqtl,nqtlsimul,k
         integer, dimension(8)          :: value

         integer                        :: nb_thread,c,igeno,jm,M1,M2
         character(len=LEN_BUFFER_WORD) :: croisement
         character(len=LENGTH_MAX_FILE) :: phases_off

         integer      ,dimension(:,:,:)  ,allocatable         :: qtl !(np+nfem+nd,nqtl,2)
         type(TYPE_LRT_SOLUTION)  , dimension(:,:,:),allocatable      :: listLrtSolSimul
         type(TYPE_LRT_SOLUTION)  ,dimension(:,:) ,allocatable        :: lrtSol
         type(TYPE_INCIDENCE_SOLUTION) ,dimension(:,:), allocatable   :: listincsol

         real (kind=dp), dimension (:,:),allocatable     :: ySAVE
         real (kind=dp), dimension (:,:,:),allocatable   :: ySIMUL
         logical                                         :: existe,ok

         real (kind=dp) , dimension(:,:),   allocatable   :: countH
         real (kind=dp) , dimension(:,:,:), allocatable   :: markH

         type(SIMULATION_INFO)                     :: simul_info

         opt_sim = COMMON_ANALYSE

         if ( dataset%cli%cli_is_transcriptomic_data() ) then
           opt_sim = TRANSCIPTOME_ANALYSE
         end if

         if ( dataset%cli%cli_is_print_all() ) then
           opt_sim = COMMON_ANALYSE
         end if

         opt_version = dataset%cli%cli_get_haplotype()

         permut_active = ( (  get_type_data_analyse(opt_calcul) == TYPE_DATA_COX ) .or. dataset%cli%cli_is_permute() )
         if ( .not. permut_active ) then
          call log_mess("*******************************",INFO_DEF)
          call log_mess("******* SIMULATION TRAITS *****",INFO_DEF)
          call log_mess("*******************************",INFO_DEF)
         else
          call log_mess("*******************************",INFO_DEF)
          call log_mess("****** PERMUTATION TRAITS *****",INFO_DEF)
          call log_mess("*******************************",INFO_DEF)
         endif

#ifdef MANAGE_OMP
        nb_thread = omp_get_max_threads()
#else
        nb_thread = 1
#endif
         !a enlever pour enlever les logs...
         call init_output_handler(dataset,opt_qtl,opt_calcul,get_name_analyse(opt_calcul),nb_thread)
         nqtlsimul=0

         ! si pas de fichier qui decrit la simulation, on simule les caracteres avec prise en comptes des donnees manquantes
         inquire( file=dataset%params%get_file_val(K_PARAMSIM), exist=existe)

         if ( .not. permut_active .and. existe  ) then
          call simul_info%read_simulation_file(dataset)
         else
          call read_genealogy(dataset)
          call read_model(dataset)
          call fixe_structure_model(dataset)
          call read_traits(dataset)
        end if

        call read_map(dataset)
        call read_genotype(dataset)
        call recup(dataset)
        call check_typage(dataset)
        call log_descriptif_genealogy(dataset)
        call log_marker_description(dataset)
        call print_freqall(dataset)
        call check_HWE(dataset)

        if ( dataset%cli%cli_is_calcul_cd() ) then
          call dataset%datasetUser%calcul_y_cd(dataset%genea,dataset%phenoAnimal,dataset%phenoModel)
          call dataset%datasetUser%calcul_corcd(dataset%genea,dataset%phenoAnimal,dataset%phenoModel)
        end if

        ! Normalizing data and print information about performance trait
        call manage_data(dataset,dataset%cli%cli_is_transcriptomic_data(),.false.,need_normalize_data(opt_calcul))

        if (.not. dataset%cli%cli_is_transcriptomic_data()) then
           call log_descriptif_traits(dataset)
        end if

        call check_permutation_context(dataset)

! **************************************************
! Information sire Heterozygote

     if ( dataset%params%get_file_val(K_INFORMATIVITY) /= '' ) then
       allocate (countH(dataset%genea%np,dataset%map%nchr))
       allocate (markH(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)))
       ! print information in a file : heterozygote marker from given sire/marker/linkage group
       call get_information_informative_marker(dataset,spt,countH,markH)
       call print_informativity_markers(dataset,dataset%params%get_file_val(K_INFORMATIVITY),countH,markH)
       deallocate (countH)
       deallocate (markH)
      end if
! **************************************************

        call haplotype(dataset,spt,opt_version)

        ! get the offspring haplotype around the Maximum Likelihood Position
        call haplotype_offspring(dataset,spt,0)


        do c=1,dataset%map%nchr
         do jm=1,dataset%genea%nm
           call log_mess('Mother :'//trim(dataset%genea%mere(jm))//'-->genotypes found [Chr '//trim(dataset%map%chromo(c))//']:'//&
       trim(str(spt%ngenom(c,jm+1)-spt%ngenom(c,jm))),VERBOSE_DEF)
           do igeno=spt%ngenom(c,jm)+1,spt%ngenom(c,jm+1)
                call log_mess('genotype '//trim(str((igeno)))//' : '//trim(str(spt%probg(c,igeno)*100.d0))//'%',VERBOSE_DEF)
           end do
          end do
        end do

        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP SECTIONS
        !$OMP SECTION
        do c=1,dataset%map%nchr
          ! Printing probabilities of transmission in associated file
          call print_pded(dataset,spt,c)
        end do
        !$OMP SECTION
        ! Printing parents genotype
        call print_phases(dataset,spt)
        !$OMP END SECTIONS NOWAIT
        !$OMP END PARALLEL

        call simul_info%log_infoqtldefinedbyuser(dataset,nqtlsimul)

        allocate (qtl(dataset%genea%np+dataset%genea%nm+dataset%genea%nd,nqtlsimul,2),STAT=stat)
        call check_allocate(stat,'qtl')

        allocate (lrtSol(dataset%phenoModel%ncar,opt_qtl))
        allocate (listincsol(dataset%phenoModel%ncar,opt_qtl+1))
        allocate(listLrtSolSimul(dataset%phenoModel%ncar,opt_qtl,nsim))

        allocate (ySAVE(dataset%phenoModel%ncar,dataset%genea%nd))
        ySAVE = dataset%phenoAnimal%y
        allocate (ySIMUL(nsim,dataset%phenoModel%ncar,dataset%genea%nd))
        call log_mess("generating simul traits...");
        do isim=1,nsim
        	call simul_info%manage_simulator_traits(dataset,spt,is_multitrait_analysis(opt_calcul),&
        	 permut_active,dataset%map%pas)
            ySIMUL(isim,:,:)=dataset%phenoAnimal%y
        end do

        dataset%phenoAnimal%y = ySAVE
        deallocate (ySAVE)

        ! Analyse according to the option opt_calcul
        call analyse_cuda(dataset,spt,nsim,opt_sim,opt_calcul,opt_qtl,lrtSol,listincsol,&
          listLrtSolSimul,ySIMUL)

       !if ( .not. is_multitrait_analysis(opt_calcul) ) then
       ! Dev Confidence Intervals/QTL Mapping
        call computingIC(dataset,lrtsol,listincsol,dataset%cli%cli_get_nqtl())
        call print_confidence_intervals_solution(dataset,dataset%cli%cli_get_nqtl(),lrtsol,is_multitrait_analysis(opt_calcul))
       !end if

        ! Printing genotype's offspring found
        call set_haplotype_offspring_context(dataset,C,M1,M2,ok,phases_off)
        if (ok) call print_offspring_phase(dataset,spt,C,M1,M2,lrtSol,opt_qtl,phases_off,&
         is_multitrait_analysis(opt_calcul))
        
        if (.not. dataset%cli%cli_is_transcriptomic_data()) then
           if (.not. is_multitrait_analysis(opt_calcul)) then
             !print to the console max LRT
              call print_information_max_lrt(dataset,opt_qtl,INFO_DEF,lrtsol)
          else
             !print to the console max LRT
              call print_information_max_lrt_multi(dataset,opt_qtl,INFO_DEF,lrtsol(1,:))
          end if
        end if
        
        call write_simulation_file(dataset,opt_qtl,is_multitrait_analysis(opt_calcul),nsim,listLrtSolSimul)
        call log_mess('computing quantiles...',VERBOSE_DEF)

        dataset%phenoModel%sigt=1.d0
        dataset%phenoModel%xmut=0.d0

        call quantile(dataset,dataset%params%get_file_val(K_OUTSIM),nsim,opt_calcul,opt_qtl,nqtlsimul,listLrtSolSimul)

        call end_output_handler(dataset,opt_qtl)

        do i=1,size(listLrtSolSimul,1)
          do j=1,size(listLrtSolSimul,2)
           do k=1,size(listLrtSolSimul,3)
            call listLrtSolSimul(i,j,k)%release()
           end do
          end do
        end do

         do i=1,size(listincsol,1)
          do j=1,size(listincsol,2)
            call listincsol(i,j)%release()
          end do
        end do


        deallocate (lrtSol)
        deallocate (listincsol)
        deallocate(listLrtSolSimul)
        deallocate (ySIMUL)



  end subroutine qtlmap_cuda_simul_and_analyse

#endif

       subroutine qtlmap_simul(dataset,isThresholdEstim,nsim,opt_calcul,opt_qtl)
        type(QTLMAP_DATASET)             ,intent(inout)            :: dataset
        logical,         intent(in)               :: isThresholdEstim
        integer         , intent(in)              :: nsim,opt_calcul,opt_qtl
        type(TYPE_LRT_SOLUTION)  , dimension(:,:,:),allocatable      :: lrtsol
        type(TYPE_INCIDENCE_SOLUTION) ,dimension(:,:), allocatable   :: listincsol
        real (kind=dp),dimension (:,:,:),allocatable :: rhoi

        integer                        :: stat,i,ic,isim,opt_sim,opt_version,j,iqtl

        logical                        :: permut_active=.false.
        character(len=LEN_DEF)           :: tempfile
        character(len=LEN_BUFFER_WORD) :: croisement

        integer                        :: nbAff,nAff,o,nqtlsimul,startSim,endSim,iq,s1,s2,sp1,sp2,nb_thread,k
        real                           :: pasAff
        logical                                   :: existe
        type(SIMULATION_INFO)                     :: simul_info
        type(GENEALOGY_BASE) ,target              :: dg

#ifdef MPI_ACTIVE
        integer                                  :: code              ! code erreur pour les routine MPI
        integer                                  :: nbsimulByThread   ! nombre de simulation par thread
        integer                                  :: dim3double        ! identifiant de type MPI pour le vecteur de double passer par message
        integer                                  :: dim3integer       ! identifiant de type MPI pour le vecteur d'entier passer par message
        integer                                  :: deb_indice
        logical           ,allocatable ,dimension(:,:,:,:) :: unpackAll ! utiliser pour decompacter un vecteur passer en message MPI

        double precision  ,allocatable ,dimension(:,:,:,:) :: lrtmaxTab
        double precision  ,allocatable ,dimension(:)       :: lrtmaxTab_lin   ! compactage des valeurs du LRT pour chaque processus
        double precision  ,allocatable ,dimension(:)       :: lrtmaxTab_final ! agregation de l ensemble des resultats et gerer par le processus maitre

        integer           ,allocatable ,dimension(:,:,:,:) :: chrmaxTab
        integer  ,allocatable ,dimension(:)                :: chrmaxTab_lin   ! compactage des valeurs du LRT pour chaque processus
        integer  ,allocatable ,dimension(:)                :: chrmaxTab_final ! agregation de l ensemble des resultats et gerer par le processus maitre

        integer           ,allocatable ,dimension(:,:,:,:) :: nxmaxTab
        integer  ,allocatable ,dimension(:)                :: nxmaxTab_lin   ! compactage des valeurs du LRT pour chaque processus
        integer  ,allocatable ,dimension(:)                :: nxmaxTab_final ! agregation de l ensemble des resultats et gerer par le processus maitre
#endif
        !INITIALISATION DE NQTLSIMUL car sinon il etait egal à 14000 (CAROLE)
        nqtlsimul=0

        if ( dataset%cli%cli_is_print_all() ) then
           opt_sim = COMMON_ANALYSE
        else
           opt_sim = SIMULATION
        end if
        opt_version = dataset%cli%cli_get_haplotype()

        permut_active = ( (  get_type_data_analyse(opt_calcul) == TYPE_DATA_COX ) .or. dataset%cli%cli_is_permute() )

         if ( .not. permut_active ) then
          call log_mess("*******************************",INFO_DEF)
          call log_mess("******* SIMULATION TRAITS *****",INFO_DEF)
          call log_mess("*******************************",INFO_DEF)
         else
          call log_mess("*******************************",INFO_DEF)
          call log_mess("****** PERMUTATION TRAITS *****",INFO_DEF)
          call log_mess("*******************************",INFO_DEF)
         endif

#ifdef MANAGE_OMP
        nb_thread = omp_get_max_threads()
#else
        nb_thread = 1
#endif

        if (.not. isThresholdEstim ) then
         call init_output_handler(dataset,opt_qtl,opt_calcul,get_name_analyse(opt_calcul),nb_thread)
        end if
        ! Initialize Random number Generator CoMmon (for Randlib)
        call inrgcm()
        ! Init the random generator
        call DATE_AND_TIME(date, time, zone, value)
        phrase=date//time//zone
        call phrtsd(phrase,iseed1,iseed2)
        call setall(iseed1,iseed2)

        ! si pas de fichier qui decrit la simulation, on simule les caracteres avec prise en comptes des donnees manquantes
        inquire( file=dataset%params%get_file_val(K_PARAMSIM), exist=existe)

        if ( .not. permut_active .and. existe  ) then
          call simul_info%read_simulation_file(dataset)
        else
          call read_genealogy(dataset)
          call read_model(dataset)
          call fixe_structure_model(dataset)
          call read_traits(dataset)

        end if
        if (.not. isThresholdEstim ) then

          if (simul_info%simulGenea .and. permut_active ) then
             call stop_application("Permutation only are available for threshold estimation.")
          end if
          if (simul_info%simulMap) then
            call simul_info%init_simul_marker(dataset)
            call sim_carte(dataset,1,simul_info%dens,simul_info%nalle)
            call dataset%map%set_absi()
            call simul_info%init_genotype_simul(dataset)
            call sim_typ(dataset)
          else
            call read_map(dataset)
            call read_genotype(dataset)
          end if
        end if

        call simul_info%log_infoqtldefinedbyuser(dataset,nqtlsimul)

        allocate (simul_info%qtl(dataset%genea%np+dataset%genea%nm+dataset%genea%nd,simul_info%nqtlsimul,2),STAT=stat)
        call check_allocate(stat,'qtl')
 
        allocate(lrtsol(dataset%phenoModel%ncar,opt_qtl,nsim))
        allocate(listincsol(dataset%phenoModel%ncar,opt_qtl+1))

        if ( isThresholdEstim ) then
          dataset%phenoModel%xmut = xmutanalyse
          dataset%phenoModel%sigt = sigtanalyse
        end if

        allocate (rhoi(dataset%map%nchr,dataset%phenoModel%ncar,dataset%phenoModel%ncar))

        if (.not. isThresholdEstim .and. .not. simul_info%simulMap) then
           call recup(dataset)
           call check_typage(dataset)
           call haplotype(dataset,spt,opt_version)
        end if

        call check_permutation_context(dataset)

        ! ON ENLEVE LES MESSAGES POUR LES SIMUL....
        call init_log(ERROR_DEF)

        !**************************************
        ! MAIN LOPP
        !**************************************
#ifdef MPI_ACTIVE

       ! on divise l ensemble des simulation en bloc
       ! chaque processus a a sa charge NSIM / NB_PROCESSUS simulation
       nbsimulByThread = nsim / (nb_procs)

       ! Ce qui revient pour nos structure de donnees a faire une boucle partielle pour chaque processus
       startSim = 1!nbsimulByThread * (thread_rang)+1
       endSim   = nbsimulByThread !* (thread_rang+1)

       allocate (lrtmaxTab_final(nsim*ncar*opt_qtl*opt_qtl))
       allocate (chrmaxTab_final(nsim*ncar*opt_qtl*opt_qtl))
       allocate (nxmaxTab_final(nsim*ncar*opt_qtl*opt_qtl))
       allocate (unpackAll(nbsimulByThread,ncar,opt_qtl,opt_qtl))
       unpackAll = .true.

       allocate (lrtmaxTab(nbsimulByThread,ncar,opt_qtl,opt_qtl))
       allocate (chrmaxTab(nbsimulByThread,ncar,opt_qtl,opt_qtl))
       allocate (nxmaxTab(nbsimulByThread,ncar,opt_qtl,opt_qtl))

       allocate (lrtmaxTab_lin(nbsimulByThread*ncar*opt_qtl*opt_qtl))
       allocate (chrmaxTab_lin(nbsimulByThread*ncar*opt_qtl*opt_qtl))
       allocate (nxmaxTab_lin(nbsimulByThread*ncar*opt_qtl*opt_qtl))

       lrtmaxTab=0.d0
       chrmaxTab=0
       nxmaxTab=0

#else
       startSim = 1
       endSim   = nsim
        if (thread_rang == 0  ) then
          nbAff=20
          if ( endSim < nbAff ) nbAff = endSim
          pasAff = real(endSim) / real(nbAff)
          nAff = 1
          write (*,fmt='(a3,'//trim(str(nbAff))//'(a4),a5)') '0%[',(' ',i=1,nbAff),']100%'
          write (*,fmt='(a3)',advance='no') ' '
        end if

#endif
        do isim=startSim,endSim
#ifndef MPI_ACTIVE
         if ( isim > int(real(nAff)*pasAff) .and. (isim < 10000) ) then
            write (*,fmt='(i4)',advance='no') isim
            nAff = nAff + 1
         end if
#endif
          if (simul_info%simulMap) then
            if (simul_info%simulGenea) then
              call dataset%genea%release()
              call simul_info%sim_genea(dataset)
            end if

            call sim_typ(dataset)
            call recup(dataset)
            call check_typage(dataset)
            call haplotype(dataset,spt,opt_version)
          end if

          call simul_info%manage_simulator_traits(dataset,spt,is_multitrait_analysis(opt_calcul),&
           permut_active,dataset%map%pas)
          call manage_data(dataset,dataset%cli%cli_is_transcriptomic_data() ,.true.,need_normalize_data(opt_calcul))

          ! Analyse according to the option opt_calcul
          call analyse(dataset,spt,opt_calcul,opt_qtl, lrtsol(:,:,isim),listincsol(:,:),rhoi,opt_sim)

          ! release incsol
          do i=1,size(listincsol,1)
           do j=1,size(listincsol,2)
            call listincsol(i,j)%release()
           end do
          end do

#ifdef MPI_ACTIVE
          do ic=1,ncar
           do iqtl=1,opt_qtl
            do iq=1,iqtl
             lrtmaxTab(isim,ic,iqtl,iq)=lrtsol(ic,iqtl,isim)%lrtmax(iq-1)
             chrmaxTab(isim,ic,iqtl,iq)=lrtsol(ic,iqtl,isim)%chrmax(iq-1)
             nxmaxTab(isim,ic,iqtl,iq)=lrtsol(ic,iqtl,isim)%nxmax(iq-1)
            end do
           end do
          end do
#endif
          !print *,'thread_rang:',thread_rang,' isim=',isim,lrtsol(1,1,isim)%lrtmax(0),lrtsol(1,1,isim)%nxmax(0)
        end do

        !On recupere à nouveau les messages
        call init_log(INFO_DEF)

#ifdef MPI_ACTIVE
!        print *,"proc:",thread_rang,"car1",lrtmaxTab(:,1,1,1),chrmaxTab(:,1,1,1),nxmaxTab(:,1,1,1)
!        print *,"proc:",thread_rang,"car2",lrtmaxTab(:,2,1,1),chrmaxTab(:,2,1,1),nxmaxTab(:,2,1,1)
!        print *,"proc:",thread_rang,"car3",lrtmaxTab(:,3,1,1),chrmaxTab(:,3,1,1),nxmaxTab(:,3,1,1)

        !compactage des infos a envoyer au maitre
        lrtmaxTab_lin = pack(lrtmaxTab,mask=.true.)
        chrmaxTab_lin = pack(chrmaxTab,mask=.true.)
        nxmaxTab_lin = pack(nxmaxTab,mask=.true.)

        !Declaration des types a passer par message
        call MPI_TYPE_CONTIGUOUS (nbsimulByThread*ncar*opt_qtl*opt_qtl, MPI_DOUBLE_PRECISION, dim3double,code)
        call MPI_TYPE_COMMIT(dim3double,code)
        call MPI_TYPE_CONTIGUOUS (nbsimulByThread*ncar*opt_qtl*opt_qtl, MPI_INTEGER, dim3integer,code)
        call MPI_TYPE_COMMIT(dim3integer,code)


        ! Ici le processus Maitre (0) recupere l'ensemble des resultats de la simulation....
        deb_indice=nbsimulByThread * (thread_rang)+1
        call MPI_GATHER (lrtmaxTab_lin(1),1, dim3double ,lrtmaxTab_final(deb_indice),1,dim3double, 0 ,MPI_COMM_WORLD ,code)
        call MPI_GATHER (chrmaxTab_lin(1),1, dim3integer ,chrmaxTab_final(deb_indice),1,dim3integer, 0 ,MPI_COMM_WORLD ,code)
        call MPI_GATHER (nxmaxTab_lin(1),1, dim3integer ,nxmaxTab_final(deb_indice),1,dim3integer, 0 ,MPI_COMM_WORLD ,code)

        call MPI_TYPE_FREE(dim3double,code)
        call MPI_TYPE_FREE(dim3integer,code)

        !Gestion seulement par le processus maitre
        if (thread_rang == 0  ) then
           !on libere l allocation
           deallocate (lrtmaxTab)
           deallocate (chrmaxTab)
           deallocate (nxmaxTab)

           !on realloue pour l ensemble de la simulation
           allocate (lrtmaxTab(nsim,ncar,opt_qtl,opt_qtl))
           allocate (chrmaxTab(nsim,ncar,opt_qtl,opt_qtl))
           allocate (nxmaxTab(nsim,ncar,opt_qtl,opt_qtl))

           !agregation des infos
           do i=1,nb_procs
             s1 = nbsimulByThread*ncar*opt_qtl*opt_qtl * (i-1)+1
             s2 = nbsimulByThread*ncar*opt_qtl*opt_qtl * (i)

             sp1 = nbsimulByThread * (i-1)+1
             sp2 = nbsimulByThread * (i)

             lrtmaxTab(sp1:sp2,:,:,:) = unpack(lrtmaxTab_final(s1:s2),mask=unpackAll,field=lrtmaxTab(sp1:sp2,:,:,:))
             chrmaxTab(sp1:sp2,:,:,:) = unpack(chrmaxTab_final(s1:s2),mask=unpackAll,field=chrmaxTab(sp1:sp2,:,:,:))
             nxmaxTab(sp1:sp2,:,:,:) = unpack(nxmaxTab_final(s1:s2),mask=unpackAll,field=nxmaxTab(sp1:sp2,:,:,:))
           end do
!             print *,""
!             print *,"Processus maitre...",nsim,ncar,opt_qtl
!             print *,""
!             print *,"car1",lrtmaxTab(:,1,1,1),chrmaxTab(:,1,1,1),nxmaxTab(:,1,1,1)
!             print *,"car2",lrtmaxTab(:,2,1,1),chrmaxTab(:,2,1,1),nxmaxTab(:,2,1,1)
!             print *,"car3",lrtmaxTab(:,3,1,1),chrmaxTab(:,3,1,1),nxmaxTab(:,3,1,1)

             ! on libere les lrtsols effectuer par le maitre
             do isim=startSim,endSim
               do i=1,size(lrtsol,1)
                do j=1,size(lrtsol,2)
                 call release(lrtsol(i,j,isim))
                end do
               end do
             end do

             !remplissage de la structure lrtsol pour calculer les quantile et ecrire par le processus maitre le fichier de simulation
             do isim=1,nsim
                do ic=1,ncar
                 do iqtl=1,opt_qtl
                  call new(opt_qtl,lrtsol(ic,iqtl,isim))
                  do iq=1,iqtl
                   lrtsol(ic,iqtl,isim)%lrtmax(iq-1) = lrtmaxTab(isim,ic,iqtl,iq)
                   lrtsol(ic,iqtl,isim)%chrmax(iq-1) = chrmaxTab(isim,ic,iqtl,iq)
                   lrtsol(ic,iqtl,isim)%nxmax(iq-1) = nxmaxTab(isim,ic,iqtl,iq)
                  end do
                 end do
                end do
             end do

             ! Pour le maitre on doit liberer tous les lrtsol
             startSim = 1
             endSim  = nsim
        end if
#endif
        !retour chariot pour l affichage de progression
        if (thread_rang == 0  ) then
         write(*,*)''

         call write_simulation_file(dataset,opt_qtl,is_multitrait_analysis(opt_calcul),nsim,lrtsol)
         call log_mess('computing quantiles...',VERBOSE_DEF)
         call quantile(dataset,dataset%params%get_file_val(K_OUTSIM),nsim,opt_calcul,opt_qtl,nqtlsimul,lrtsol)

         if (.not. isThresholdEstim ) call end_output_handler(dataset,opt_qtl)
        end if

         do isim=startSim,endSim
          do i=1,size(lrtsol,1)
           do j=1,size(lrtsol,2)
            call lrtsol(i,j,isim)%release()
           end do
          end do
         end do

        deallocate(lrtsol)
        deallocate (simul_info%qtl)
        deallocate (rhoi)

#ifdef MPI_ACTIVE
        deallocate (lrtmaxTab)
        deallocate (chrmaxTab)
        deallocate (nxmaxTab)
        deallocate (lrtmaxTab_lin)
        deallocate(lrtmaxTab_final)
        deallocate (chrmaxTab_lin)
        deallocate(chrmaxTab_final)
        deallocate (nxmaxTab_lin)
        deallocate(nxmaxTab_final)
        deallocate (unpackAll)
#endif
       end subroutine qtlmap_simul


       subroutine print_information_max_lrt(dataset,opt_qtl,level,lrtsol)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         integer            , intent(in)                      :: opt_qtl,level
         type(TYPE_LRT_SOLUTION)  , dimension(dataset%phenoModel%ncar,opt_qtl), intent(in)      :: lrtsol

         integer :: i,c,iqtl,nqtl,iqtl2
         character(len=LEN_LINE) :: dxlist

         do nqtl=1,opt_qtl
          do iqtl=0,nqtl-1

           call log_mess('    ***************************************',level)
           call log_mess('    *******   Test '//trim(str(iqtl))//' QTL / '//trim(str(nqtl))//' QTL   ********',level)
           call log_mess('    ***************************************',level)

           do i=1,dataset%phenoModel%ncar
            if ( .not. associated(lrtsol(i,nqtl)%lrtmax) ) cycle
            call log_mess('Trait ['//trim(trim(dataset%phenoModel%carac(i)))//']',level)
            call log_mess('Maximum likelihood ratio test :',level)

            dxlist=adjustl("[")
          !  vallist="["

            do iqtl2=0,nqtl-1
              dxlist=trim(dxlist)//trim(str(dataset%map%absi(lrtsol(i,nqtl)%chrmax(iqtl2),lrtsol(i,nqtl)%nxmax(iqtl2))))//&
              ", Ch "//trim(dataset%map%chromo(lrtsol(i,nqtl)%chrmax(iqtl2)))
              if ( iqtl2 < nqtl-1) dxlist=trim(dxlist)//"] ["
            !  vallist=vallist//" "//trim(str(lrtsol(i,nqtl)%lrtmax(iqtl2)))
            end do
            dxlist=trim(dxlist)//"]"

            call log_mess('The maximum is reached at position '//trim(dxlist)&
            //' Morgan, with value '//trim(str(lrtsol(i,nqtl)%lrtmax(iqtl))),level)
          end do
         end do
        end do

       end subroutine


       subroutine print_information_max_lrt_multi(dataset,opt_qtl,level,lrtsol)
         type(QTLMAP_DATASET)                     ,intent(in) :: dataset
         integer            , intent(in)                      :: opt_qtl,level
         type(TYPE_LRT_SOLUTION)  , dimension(opt_qtl), intent(in)      :: lrtsol
         character(len=LEN_L)         :: list_traits,sep

         integer :: i,c
         list_traits = ''
         sep=' '
         do i=1,dataset%phenoModel%ncar
            list_traits=list_traits//trim(dataset%phenoModel%carac(i))//sep
         end do
         call log_mess('    ****   Join traits [ '//trim(list_traits)//']  ****',level)
          call log_mess('    ***************************************',level)
          call log_mess('    ***** Test 0 QTL / 1 QTL   ************',level)
          call log_mess('    ***************************************',level)
          call log_mess('Maximum likelihood ratio test :',level)
          call log_mess('The maximum is reached at position ['//trim(str(dataset%map%absi(lrtsol(1)%chrmax(0),lrtsol(1)%nxmax(0))))&
          //", Ch "//trim(dataset%map%chromo(lrtsol(1)%chrmax(0)))//&
          '] M, with value '//trim(str(lrtsol(1)%lrtmax(0))),level)

!          if (opt_qtl == 2 ) then
!           call log_mess('Test 0 QTL / 2 QTL',level)
!           call log_mess('Maximum likelihood ratio test :',level)
!           call log_mess('The maximum is reached at position ['//trim(str(absi(lrtsol(1)%chrmax(0),lrtsol(1)%nxmax(0))))//&
!             ',Ch '//trim(chromo(lrtsol(1)%chrmax(0)))//'] and ['//trim(str(absi(lrtsol(1)%chrmax(1),lrtsol(1)%nxmax(1))))//&
!             ',Ch '//trim(chromo(lrtsol(1)%chrmax(1)))//'] M, with values [0/2] '//trim(str(lrtsol(1)%lrtmax(0)))//&
!             ' and [1/2] '//trim(str(lrtsol(1)%lrtmax(1))),level)
!          end if
       end subroutine


end program qtlmap

