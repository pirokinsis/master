      program main     

      Implicit None

      Real, Dimension (:,:), allocatable :: a    

      Integer :: i, j

      ! one line per read
C      Write( *, * ) 'Line at a time'
      Open( 10, file = 'example.txt' )
      Read( 10 , *) j

      allocate ( a(j,j) )

      Do i = 1, j
            Read ( 10, * ) a( i, : )
C           Write(  *, * ) a( i, : )
      End Do

      Write (*,*) a
      Close( 10 )

      end program main 