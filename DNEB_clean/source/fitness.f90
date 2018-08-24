!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   FITNESS ROUTINES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! These are functions to load the data from the GMM and then to calc         !!
!! the gradient and the energy at these positions.                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module fitness
    integer, parameter :: NDIMENSION = 7, NCOMPONENT = 8
    integer :: ORIG_DIM = 7
    integer :: SCALING_METHOD = 2
    ! 0 is no scaling, 1 is the symmetric
    ! 2 is the general scaling
    ! ONLY NEED ORIG DIM IF SCALING METHOD IS 1
    ! TODO: This definitely needs neatening up at some point!

    double precision, allocatable :: WEIGHT(:), DETERMINANTS(:), MEANVECTOR(:,:)
    double precision, allocatable :: MATRIX(:,:,:), SIG(:), SCALING(:)
    real, parameter :: PI = 3.1415927

    CONTAINS

    subroutine fitness_init()

        INTEGER J1, J2, J3
        DOUBLE PRECISION PI
        character(len=50) :: filename
        INTEGER TEST(2)

        ALLOCATE(SIG(NCOMPONENT))
        ALLOCATE(WEIGHT(NCOMPONENT))
        ALLOCATE(DETERMINANTS(NCOMPONENT))
        ALLOCATE(MEANVECTOR(NDIMENSION,NCOMPONENT))
        ALLOCATE(MATRIX(NDIMENSION, NDIMENSION,NCOMPONENT))
        ALLOCATE(SCALING(NCOMPONENT))

        OPEN(UNIT=222,FILE='data/fortran/weight.out',status='old')
        READ(222,*) WEIGHT(1:NCOMPONENT)
        CLOSE(222)

        OPEN(UNIT=222,FILE='data/fortran/determinants.out',status='old')
        READ(222,*) DETERMINANTS(1:NCOMPONENT)
        CLOSE(222)

        OPEN(UNIT=222,FILE='data/fortran/meanvector.out',status='old')
        DO J1 = 1, NCOMPONENT
            READ(222,*) MEANVECTOR(1:NDIMENSION,J1)
        END DO
        CLOSE(222)

        DO J2=1, NCOMPONENT

            write (filename, fmt='(A I8 A)') "data/fortran/matrix", J2, ".out"
            call strip_spaces(filename)
            filename = trim(filename)
            OPEN(UNIT=222,FILE=filename,status='old')

            DO J1 = 1, NDIMENSION
                READ(222,*) MATRIX(1:NDIMENSION,J1,J2)
            END DO
            CLOSE(222)
        END DO

        ! ONLY LOAD SIGMA VALUES IF WE NEED THEM
        IF(SCALING_METHOD == 1) THEN
            OPEN(UNIT=222,FILE='data/fortran/sig.out',status='old')
            READ(222,*) SIG(1:NCOMPONENT)
            CLOSE(222)
        END IF

        ! ONLY LOAD THE SCALINGS IF WE NEED THEM
        IF(SCALING_METHOD == 2) THEN
            OPEN(UNIT=222,FILE='data/fortran/scalings.out',status='old')
            READ(222,*) SCALING(1:NCOMPONENT)
            CLOSE(222)
        END IF

    end subroutine fitness_init

    ! Calculate the gradient of the fitness surface
    subroutine fitness_gradient(COORDS, E, GRAD)

        double precision, intent(in) ::  COORDS(NDIMENSION)
        double precision, intent(out) :: GRAD(NDIMENSION)
        double precision, intent(out) :: E
        INTEGER J1, J2, J3
        DOUBLE PRECISION PROD(NCOMPONENT), PROB(NCOMPONENT), VECT(NDIMENSION)

        E = 0.0D0
        PROD(:) = 0.0D0
        PROB(:) = 0.0D0
        GRAD = 0d0

        DO J1 = 1, NCOMPONENT
             ! Starting vector = 0 for this component
         VECT(:) = 0.0D0

         DO J2 = 1, NDIMENSION
            DO J3 = 1, NDIMENSION
                         ! This is doing the sigma * (x - mu)
               VECT(J2) = VECT(J2) + MATRIX(J2,J3, J1) * ( COORDS(J3)-MEANVECTOR(J3, J1) )
            ENDDO
                    ! This is doing th (x - mu) * that
            PROD(J1) = PROD(J1) + VECT(J2) * (COORDS(J2)-MEANVECTOR(J2,J1))
         ENDDO
             ! This calculates the probability and the energy

         IF (SCALING_METHOD == 0) THEN
                PROB(J1) = EXP( -PROD(J1)/2.0 ) / SQRT( DETERMINANTS(J1) )
             ELSE IF (SCALING_METHOD == 1) THEN
                PROB(J1) = EXP( -PROD(J1)/2.0 ) / SQRT( DETERMINANTS(J1) )*(1/( SQRT( 2*PI*SIG(J1) ) ))**(ORIG_DIM-NDIMENSION)
             ELSE
                  PROB(J1) = EXP( -PROD(J1)/2.0 ) / SQRT( DETERMINANTS(J1) )*SCALING(J1)
         END IF

         E = E + WEIGHT(J1)*PROB(J1)

             ! Once the probability has been calculated we can calculate the gradient for each direction
         DO J2 = 1, NDIMENSION
            GRAD(J2) = GRAD(J2) - PROB(J1) * VECT(J2) * WEIGHT(J1)
         END DO

        END  DO

        ! Logging the energy
        GRAD = -GRAD/E
        E = -log(E)

    end subroutine fitness_gradient

    ! This is a dirty subroutine for removing spaces from a string
    ! this is neccesary for the file names matrix1.out
    ! otherwise we end up with matrix   1.out
    subroutine strip_spaces(string)
        character(len=*) :: string
        integer :: string_len
        integer :: last, actual

        string_len = len (string)
        last = 1
        actual = 1

        do while (actual < string_len)
            if (string(last:last) == ' ') then
                actual = actual + 1
                string(last:last) = string(actual:actual)
                string(actual:actual) = ' '
            else
                last = last + 1
                if (actual < last) &
                    actual = last
            endif
        end do
    end subroutine strip_spaces

end module fitness
