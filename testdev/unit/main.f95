program test_main
  use fruit
  use test_type_lrt_solution  
  use test_type_map
  implicit none

   print *,"data_test_package"
    call init_fruit
    call all_test_type_map
    !call tested_position
    call fruit_summary


end program test_main
