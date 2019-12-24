module m_qtlmap_type_parameter_simulation
    use m_qtlmap_constant
    use m_qtlmap_log
    use m_qtlmap_base
    implicit none

      type, private :: PARAMETER_GENEALOGY
         character(len=LEN_DEF)    , public :: croisement  !// F2_KEYWORD,BC_KEYWORD,OUTBRED_KEYWORD
         integer                   , public :: nbpere      !// number of sire
         integer                   , public :: inmp        !// number of dam by sire
         integer                   , public :: indm        !// number of prog by dams

      end type PARAMETER_GENEALOGY

      type, private :: PARAMETER_MAP

         integer                   , public :: nalle !// allele number by marker
         real (kind=dp)            , public :: dens  !// proportion of marker
         real (kind=dp)            , public :: taille
      end type



    type, public :: PARAMETER_SIMULATION

       type (PARAMETER_GENEALOGY) , pointer, public :: genealogy => null()
       type (PARAMETER_MAP)       , pointer, public :: map       => null()

    end type PARAMETER_SIMULATION

end module m_qtlmap_type_parameter_simulation
