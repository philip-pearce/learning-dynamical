!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   OUTPUT ROUTINES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! These are functions to output the data in a neat format for displaying     !!
!! with a python script.                                                      !!                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module inoutput

    contains

subroutine output_array(X, output, NDIMENSION, N)
        double precision, intent(in) :: X(:)
        integer, intent(in) :: NDIMENSION, N
        double precision :: G(NDIMENSION)
        CHARACTER (LEN=*), intent(in) :: output
        integer :: i, j
        double precision :: E

        open (unit=41, file=output, access='append', action='write')

        do j=1, NDIMENSION
            do i=1, N+2*NDIMENSION, NDIMENSION
                write(41, fmt='(A1 E13.5E3)', advance='no') ' ', X(i+j-1)
            end do
                write(41, *)
        end do

        write(41, *)
        write(41, *)
        write(41, *)

        close (unit=41)

    end subroutine output_array

    subroutine output_minima(minima, output, NDIMENSION, N)
        ! Variables passed in
        double precision, intent(in) :: minima(:)
        integer, intent(in) :: NDIMENSION, N
        CHARACTER (LEN=*), intent(in) :: output

        ! Local variables
        integer :: i, j

        open (unit=120, file=output, status='replace', action='write')

        j=1
        do i=1, size(minima), NDIMENSION+1
            write(120, *) "Minima: ", j
            write(120, *) "Coord: ", minima(i:i+NDIMENSION-1)
            write(120, *) "Energy: ", minima(i+NDIMENSION)
            write(120, fmt='()', advance='yes')             
                j=j+1
        end do

        close (unit=120)

    end subroutine output_minima

    subroutine output_barriers(X, barrier, barrier_pos, saddle, min_num1, min_num2, NDIMENSION, N)
        ! Variables passed in
        double precision, intent(in) :: X(:)
        double precision, intent(in) :: barrier, saddle
        integer, intent(in) :: min_num1, min_num2
        integer, intent(in) :: NDIMENSION, N
        double precision, intent(in) :: barrier_pos(NDIMENSION)

        ! Local variables
        integer :: i, j

        ! Writes to the fine storing important values
        open (unit=41, file='output/barriers.txt', access='append', action='write')
            write(41, *) "Barrier between minima ", min_num1, " and ", min_num2
            write(41, *) "Coords of min ", min_num1, ":", X(1:NDIMENSION)
            write(41, *) "Coords of min ", min_num2, ":", X(N+NDIMENSION+1:N+2*NDIMENSION)
            write(41, *) "Coords of barrier: ", barrier_pos(1:NDIMENSION)
            write(41, *) "Barrier energy: ", barrier
            write(41, *) "Saddle energy: ", saddle
            write(41, fmt='()', advance='yes')
        close(unit=41)

    end subroutine output_barriers


    !! These read and write in a binary format for storing midway through a run and then continuing on.
    subroutine write_minima(minima_array, NDIMENSION)

        double precision, intent(in) :: minima_array(:)
        integer, intent(in) :: NDIMENSION

        integer :: i
        integer :: inu = 20, outu = 21

        open(unit=outu,form="unformatted",file='output/minima_snapshot',action="write")
        write (outu) size(minima_array)
        do i=1,size(minima_array),NDIMENSION+1
            write (outu) minima_array(i:i+NDIMENSION) ! write one row at a time
        end do
        close (outu)

    end subroutine    

    subroutine read_minima(minima_array, NDIMENSION)

        double precision, allocatable, intent(out) :: minima_array(:)
        integer, intent(in) :: NDIMENSION

        integer :: i, minima_array_size
        integer :: inu = 20, outu = 21

        open (unit=inu,form="unformatted",file='output/minima_snapshot',action="read")
        read (inu) minima_array_size
        allocate(minima_array(minima_array_size))
        do i=1,size(minima_array),NDIMENSION+1
            read (inu) minima_array(i:i+NDIMENSION) ! write one row at a time
        end do
        close (inu)

    end subroutine    


end module inoutput
