module test_type_lrt_solution
  use fruit
  use m_qtlmap_type_map
  use m_qtlmap_type_lrt_solution

  implicit none

  public :: all_test_type_lrt_solution

  contains

    subroutine all_test_type_lrt_solution

        type(MAP_BASE):: map
        real(kind=dp) :: posi = 0.01
        !initialisation d'une carte
        posi=0.01
        call map%add_marker("Z","m1",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m2",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m3",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m4",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m5",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m6",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m7",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m8",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m9",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m10",posi,posi,posi);posi=posi+0.01
        call map%add_marker("Z","m11",posi,posi,posi);posi=posi+0.01

        call tested_position(map)

    end subroutine all_test_type_lrt_solution


    subroutine tested_position(map)
      type(MAP_BASE)      ,intent(in) :: map
      type(LIST_TESTED_POSTION_VALUE) :: list
      type(QTLMAP_DATASET)            :: dataset
      integer :: nqtl,hyp,chr(10),pos(10)
      real(kind=dp) :: val

      call assert_true(associated(list%values))
      call dataset%set()

      nqtl=1;chr=1;pos=1;val=0.0;hyp=1
      call list%add(dataset,nqtl,chr,pos,val,hyp)

      call dataset%release()
      !call assertEquals(list%values == null())
      !call assertEquals(list%values /= null())

    end subroutine tested_position

end module

