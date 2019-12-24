! Ce module a à sa charge la fabrication d'un jeu de donnees QTLMAP_DATATYPE
! La fabrication peut etre de plusieurs nature :
! - construction à partir d'un jeu de donnees utilisateur
! - constrcution partiel a partir des fichiers de donnees utilisateur puis simulation des autres de donnees
! - construction par simulation
!
!
! auteur : olivier.filangi@rennes.inra.fr
module m_qtlmap_factory_dataset
    use m_qtlmap_types
    use m_qtlmap_genealogy
    use m_qtlmap_map
    use m_qtlmap_genotype
    use m_qtlmap_phenotype
    use m_qtlmap_simulation

    implicit none

    type, public :: FACTORY_QTLMAP_DATASET


    contains
      procedure, public :: buildQTLMapDatasetWithParameterAnalyse

      procedure, private :: buildQTLMapDatasetUserFiles
      procedure, private :: check_file

    end type FACTORY_QTLMAP_DATASET

    contains

      ! - Par defaut si une option de simulation est demandé, si les parametres de simulation sont initailisé on genere le jeu de donnees
      ! - Si analyse seul, aucune generation par simulation est possible
      subroutine buildQTLMapDatasetWithParameterAnalyse(this,dataset,rebuild_pdd)
       class(FACTORY_QTLMAP_DATASET), intent(in)   :: this
       type(QTLMAP_DATASET)         ,intent(inout) :: dataset
       logical                      ,intent(inout) :: rebuild_pdd ! need to recompute the PDD_BUILD structure


       ! Si l'utilisateur ne demande pas de simulation nous analysons seulement
       ! sans prendre en compte de fichier de parametre des simulations

       if ( dataset%cli%cli_is_simulation()) then
         call this%check_file(dataset,K_GENEA)
         call this%check_file(dataset,K_GENOTYPE)
         call this%check_file(dataset,K_MAP)
         call this%check_file(dataset,K_TRAITS)
         call this%check_file(dataset,K_MODEL)

         call this%buildQTLMapDatasetUserFiles(dataset)
         rebuild_pdd = .false.
       else ! simulation CASE

       end if


      end subroutine buildQTLMapDatasetWithParameterAnalyse


      subroutine buildQTLMapDatasetUserFiles(this,dataset)
        class(FACTORY_QTLMAP_DATASET), intent(in) :: this
        type(QTLMAP_DATASET)       ,intent(inout) :: dataset
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

      end subroutine buildQTLMapDatasetUserFiles


      subroutine check_file(this,dataset,key)
        class(FACTORY_QTLMAP_DATASET), intent(in) :: this
        type(QTLMAP_DATASET)       ,intent(inout) :: dataset
        character(len=LEN_L)       ,intent(in)    :: key

        if ( .not. dataset%params%key_exist(key) ) then
           call stop_application(trim(key)//" is not defined in the parameter user file.")
         end if

      end subroutine check_file

end module m_qtlmap_factory_dataset
