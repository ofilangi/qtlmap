module test_type_map
    use fruit
    use m_qtlmap_type_map
    implicit none


    contains

       subroutine all_test_type_map
        type(MAP_BASE) :: map
        real(kind=dp) :: posi = 0.01

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

         call test_creation_map
         call test_getChridPosId(map)
       end subroutine all_test_type_map

       subroutine test_creation_map
        type(MAP_BASE) :: map
        real(kind=dp) :: posi = 0.01

        call map%add_marker("X","m1",0.01_dp,0.01_dp,0.01_dp)
        call assertEquals(map%nchr,1,"nombre de chromosome (nchr)")
        call assertEquals(map%chromo(1),"X","nom chromosome")
        call map%add_marker("Y","m1",0.01_dp,0.01_dp,0.01_dp)
        call assertEquals(map%nchr,2,"nombre de chromosome (nchr)")
        call assertEquals(map%chromo(1),"X","le nom du chromosome 1 existe tjs")
        call assertEquals(map%chromo(2),"Y","le nom du chromosome 2")
        call assertEquals(map%nmk(1),1,"Il y a toujours 1 seul marker sur le chromosome 1")
        call map%add_marker("X","m1",0.01_dp,0.01_dp,0.01_dp)
        call assertEquals(map%nmk(1),1,"Ajout d'un marqueur avec les memes position (le marqueur ne doit pas etre enregistré).")
        call map%add_marker("X","m1",0.02_dp,0.02_dp,0.02_dp)
        call assertEquals(map%nmk(1),1,"Ajout d'un marqueur avec le meme nom (le marqueur ne doit pas etre enregistré).")
        call map%add_marker("X","m2",0.02_dp,0.02_dp,0.02_dp)
        call assertEquals(map%nmk(1),2,"Ajout d'un deuxieme marqueur sur le choromose 1. test sur le nombre de marquer.")

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
        call assertEquals(map%nmk(3),11,"Ajout de 11 marqueurs sur le chromosome Z.")
        call assertEquals(map%chromo(1),"X","le nom du chromosome 1")
        call assertEquals(map%chromo(2),"Y","le nom du chromosome 2")
        call assertEquals(map%chromo(3),"Z","le nom du chromosome 3")
        call assertEquals(map%nchr,3,"nombre de chromosome (nchr)")
        call assertEquals(map%nmk(1),2,"2 marqueurs sur X.")
        call assertEquals(map%nmk(2),1,"1 marqueur sur Y.")

       end subroutine test_creation_map

     subroutine test_getChridPosId(map)
      type(MAP_BASE)      ,intent(in) :: map
      integer, dimension(3) :: chr,pos
      integer :: chrid,posid,chrid2,posid2,chrid3,posid3,n
      chr(1)=1
      pos(1)=1
      n=1
      call map%getChridPosId(n,chr,pos,chrid,posid)
      call assertEquals(chr(1),chrid,"Test basic d'equivalence chr sur la routine map%getChridPosId")
      call assertEquals(pos(1),posid,"Test basic d'equivalence chr sur la routine map%getChridPosId")
      chr(1)=1
      pos(1)=10
      call map%getChridPosId(n,chr,pos,chrid,posid)
      call assertEquals(chr(1),chrid,"Test basic d'equivalence chr sur la routine map%getChridPosId")
      call assertEquals(pos(1),posid,"Test basic d'equivalence chr sur la routine map%getChridPosId")

    end subroutine test_getChridPosId


end module test_type_map
