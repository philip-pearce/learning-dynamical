!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   DNEB ROUTINES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This contains the sub-routines required to find the minimum points using   !!
!! steepest descent and then the neb method to find the MEP                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module neb
    use fitness
    use inoutput

    integer, parameter:: IMAGES = 100            !! This is how many IMAGES are used in the neb method
    integer, parameter :: N = IMAGES*NDIMENSION   !! How many co-ordinates will be in the path between the minima
   
    contains
    !! This subroute finds all the minima in the system via jumping around max_jumps number of times,
    !! it then minimises using steepest descent with dx and finds the minimum when the gradient is less 
    !! than max_rms. If the minima is unique (defined by similarity) it is then stored in minima_array
    subroutine find_minima(max_jumps, jump_amp, similarity, iter_max, dx, max_rms, minima_array)
        !! Variables which are passed in and out
        integer, intent(in) :: max_jumps, iter_max
        double precision, intent(in) :: dx, similarity, max_rms, jump_amp
        double precision, allocatable, intent(out) :: minima_array(:)

        !! Local variables
        double precision :: g(NDIMENSION), x(NDIMENSION), e, rms
        integer :: i, iter, converged, jump, accept

        !! Initial values
        x = 0d0
        jump = 0

        !! Main loop to find minima
        print *, "########## MINIMA ##########"
        do while ( jump<max_jumps )

            jump = jump+1
            iter = 1
            converged = 0

            !! Jump to new set of co-ordinates
            call basin_hop(jump_amp, x)

            !! Little cheat for finding gMM minima
            if (max_jumps <= 40 .and. jump <= NCOMPONENT) then
                x = meanvector(:, jump)
            end if

            !! Steepest descent until minimised
            do while ( (iter<iter_max) .and. converged .eq. 0)
                call fitness_gradient(x,e,g)
                x = x - g*dx
                rms = sum( g**2 )

                if (rms < max_rms) then
                    converged = 1
                end if

                !! This checks if we're at one the infinite energy regions
                if ( rms .ne. rms ) then
                    iter = iter_max
                    converged = 0
                end if

                iter = iter+1
            end do

            !! Check similarity before storing the minima
            if ( (converged .eq. 1) ) then

                !! Intialise the arrays if not already
                if ( .not. allocated(minima_array) ) then 
                    minima_array = [x, e]
                    print *, "Located minima 1"
                end if    

                if( allocated(minima_array) ) then
                    accept = 1
                    do i=1, size(minima_array), NDIMENSION+1
                        if( (abs(e - minima_array(i+NDIMENSION)) < similarity) .and. &
                            (  sum( (minima_array(i:i+NDIMENSION-1)-x)**2 )/NDIMENSION < similarity ) ) then
                            accept = 0
                        end if
                    end do

                    if (accept == 1) then
                        minima_array = [minima_array, x, e]
                        print '((A15,I2))' ,'Located minima', size(minima_array)/(NDIMENSION+1)
                    end if
                end if    
            end if

        end do

    call output_minima(minima_array, 'output/minima.txt', NDIMENSION, N)

    end subroutine find_minima

    subroutine find_paths(minima_array, iter_max, dx, max_rms, spring_const, perturb)
        double precision, allocatable, intent(in) :: minima_array(:)
        double precision, intent(in) :: dx, spring_const, max_rms, perturb
        integer, intent(in) :: iter_max

        ! Local variables
        double precision :: min_1(NDIMENSION), min_2(NDIMENSION), path(N+2*NDIMENSION), saddle, barrier, barrier_pos(NDIMENSION)
        integer :: i, j, k
        integer, allocatable :: start_point(:)
        integer, allocatable :: end_point(:)

        !! Clear the barrier from a previous run
        open (unit=41, file='output/barriers.txt', status='replace', action='write')
        write(41, fmt='()', advance='yes')
        close(unit=41)
        open (unit=41, file='output/path.txt', status='replace', action='write')
        write(41, fmt='()', advance='yes')
        close(unit=41)


        !! Make arrays of staring points and end points
        if ( (size(minima_array) > NDIMENSION+1) .and. (allocated(minima_array)) ) then
            do i=1, size(minima_array)-( 2*(NDIMENSION) ), NDIMENSION+1
                do j=i+NDIMENSION+1, size(minima_array), NDIMENSION+1
                    if ( .not. allocated(start_point) ) then 
                        start_point = [i]
                        end_point = [j]
                    else
                        start_point = [start_point, i]
                        end_point = [end_point, j]
                    end if 
                end do
            end do
        end if

        !! Loop through each pair of start points and end points and find the barrier and saddle
        print *, "######### BARRIERS #########"

        !$OMP PARALLEL DO PRIVATE(i, j, path, saddle, barrier, barrier_pos)
        do k=1, size(end_point), 1
            i = start_point(k)
            j = end_point(k)

            ! Input: min_1, min_2, iter_max, dx, max_rms, spring_const
            ! Ouput: path, saddle, barrier
            call minimise_path(minima_array(i:i+NDIMENSION-1), minima_array(j:j+NDIMENSION-1), &
                               iter_max, dx, max_rms, dx, path, saddle, barrier, barrier_pos, &
                                (i/(NDIMENSION+1)+1), (j/(NDIMENSION+1)+1), perturb)

            ! Outout the barrier, we put this in a critical wrapper to make sure it's only done by 1 thread at one
            ! time, luckily as this is such a quick step the effect of this is negligable
            !$OMP CRITICAL 
            call output_barriers(path, barrier, barrier_pos, saddle, (i/(NDIMENSION+1)+1), (j/(NDIMENSION+1)+1), NDIMENSION, N)
            call output_array(path, 'output/path.txt', NDIMENSION, N)
            !$OMP END CRITICAL 

        end do 
        !$OMP END PARALLEL DO

    end subroutine find_paths

    !! This determines the new perturbed coordinates
    subroutine basin_hop(jump_amp, x)
        !! A random x is outputted
        double precision, intent(in) :: jump_amp
        double precision, intent(out) :: x(NDIMENSION)

        double precision :: delta_x(NDIMENSION)

        call random_number(delta_x)
        x = delta_x*jump_amp

    end subroutine

    !! This subroutine calculates the minimum path using the neb routine
    !! based on steepest descent and outputs the path, saddle and barrier
    subroutine minimise_path(min_1, min_2, iter_max, dx, max_rms, spring_const, path, &
        saddle, barrier, barrier_pos, loc1, loc2, perturb)
        !! Variables which are passed in and out
        double precision, intent(in) :: dx, min_1(NDIMENSION), min_2(NDIMENSION), spring_const, max_rms, perturb
        integer, intent(in) :: iter_max, loc1, loc2
        double precision, intent(out) :: path(N+2*NDIMENSION), barrier, saddle, barrier_pos(NDIMENSION)
        double precision :: dx_normalised

        double precision :: dx_array(N)

        !! Local variables
        double precision :: g(N), x(N), e, rms
        integer :: i, iter, converged, climbing_start_iter, climbing_start, max_image

        ! Linearly interpolated coordinates with the minima at either end as our first guess for the MEP
        call linear_interp(min_1, min_2, x)

        iter = 1
        converged = 0
        climbing_start = 0

        !! this makes sure that the step doesn't scale with dimension
        dx_normalised = dx / (NDIMENSION**0.5)

        call gradient_neb(x, min_1, min_2, spring_const, climbing_start, e, g, saddle, barrier, barrier_pos, max_image)
        print *, e, rms, iter, barrier, saddle

        do while ( (iter<=iter_max) .and. converged .eq. 0)

            ! Calculating the gradient and updating the path
            call gradient_neb(x, min_1, min_2, spring_const, climbing_start, e, g, saddle, barrier, barrier_pos, max_image)
            dx_array = abs(g)*dx

            x = x - g*dx_normalised
            !x = x - g*dx_array

            ! checking for full convergence
            rms = sum( g**2 )/(IMAGES*NDIMENSION)

            if ( rms .ne. rms ) then
                iter = iter_max
                converged = 0
            end if    

            if (rms < max_rms) then

                if( (climbing_start .eq. 1) .and. (iter > climbing_start_iter+10000) ) then
                    converged = 1
                end if

                if( (climbing_start .eq. 0) ) then
                    print *, '   ******* Beginning to climb *******'
                    climbing_start = 1
                    climbing_start_iter = iter
                end if   

            end if

            if ( mod(iter,1000000) == 0) then
                print *, e, rms, iter, barrier, saddle
            end if

            if ( mod(iter,100000000) == 0) then
                print *, 'Resetting the initial guess to perturb linear path'
                call linear_interp_perturb(min_1, min_2, x, perturb)
                !! [reset x to a random perturbation of a linear interpolation]
            end if

            iter = iter+1
        end do

        print '((A31,I2,A4,I2))', " ## Summary for MEP from minima ", loc1, " to ", loc2

        print *, 'Iter: ', iter
        print *, 'RMS: ', rms

        if(converged .eq. 0) then
            print *, "DID NOT CONVERGE, DO NOT TRUST MEP: TRY DIFFERENT SETTINGS"
            !! dirty way of generating nan's, some compilers may not like this
            barrier = barrier / 0.
            saddle = saddle / 0.
        end if    

        if ( rms .ne. rms ) then
            print *, "NaN RMS suggests a discontinuous surface"
        end if

        ! Storing the minima and the path in one array
        path(1:NDIMENSION) = min_1
        path(N+NDIMENSION+1 : N+NDIMENSION*2) = min_2
        path(NDIMENSION+1:(NDIMENSION)*(IMAGES+1)) = x

    end subroutine minimise_path

    subroutine gradient_neb(x, min_1, min_2, spring_const, climbing_start, e, g, saddle, barrier, barrier_pos, max_image)
        !! Variables passed in and out
        double precision, intent(in) :: x(N), min_1(NDIMENSION), min_2(NDIMENSION), spring_const
        integer, intent(in) :: climbing_start ! whether climbing is on
        integer, intent(inout) :: max_image
        double precision, intent(out) :: g(N), e, barrier, saddle, barrier_pos(NDIMENSION)

        !! Local variables
        integer :: i, image
        double precision :: pos_image(NDIMENSION), pos_image_ahead(NDIMENSION)
        double precision :: pos_image_behind(NDIMENSION), tangent(NDIMENSION), energy(IMAGES)
        double precision :: sep_a, sep_b, surface(NDIMENSION), spring_tot(NDIMENSION), e1, e2
        double precision :: e_image, e_image_behind, e_image_ahead, e_max, e_min, e_surf

        ! Loops through calcualting the gradient for each image
        e = 0d0
        do i=1, N-NDIMENSION+1, NDIMENSION
            ! grabs the position of the current image
            pos_image = x(i:i+NDIMENSION-1)

            ! Finding the position directly behind us and directly ahead
            ! if we are the first/last image then grabbing the minima
            if (i .eq. 1) then
                pos_image_behind = min_1
            else
                pos_image_behind = x(i-NDIMENSION:i-1)
            end if

            ! Finding the points directly ahead
            if (i .eq. N-NDIMENSION+1) then
                pos_image_ahead = min_2
            else
                pos_image_ahead = x(i+NDIMENSION:i+2*NDIMENSION-1)
            end if

            ! calculating the energy of the current image from the surface and 
            ! the ones directly adjacent setting up the energies behind and ahead
            if (i .eq. 1) then
                call fitness_gradient(pos_image, e_image, g)
                call fitness_gradient(pos_image_ahead, e_image_ahead, g)
                call fitness_gradient(min_1, e_image_behind, g)
            else if (i .eq. N-NDIMENSION+1) then
                e_image_behind = e_image
                e_image = e_image_ahead
                call fitness_gradient(min_2, e_image_ahead, g)
            else
                e_image_behind = e_image
                e_image = e_image_ahead
                call fitness_gradient(pos_image_ahead, e_image_ahead, g)
            end if

            ! calculates the seperations
            sep_a = norm2( pos_image-pos_image_ahead )
            sep_b = norm2( pos_image-pos_image_behind )

            ! There are multiple ways to calculate the tangent and the spring force
            ! as given by http://aip.scitation.org/doi/pdf/10.1063/1.1323224
            ! We use the improved method with less kinks and more equal spacing between the IMAGES
            if (( (e_image_ahead < e_image) .AND. (e_image_behind < e_image) ) .or. &
             ( (e_image_ahead > e_image) .AND. (e_image_behind > e_image) )) then
                
                 e_max = max(abs(e_image_behind-e_image), abs(e_image_ahead-e_image))
                 e_min = min(abs(e_image_behind-e_image), abs(e_image_ahead-e_image))

                if( e_image_behind < e_image_ahead ) then
                    tangent = (pos_image-pos_image_behind)*e_max &
                    + (pos_image_ahead - pos_image)*e_min
                else
                    tangent = (pos_image-pos_image_behind)*e_min &
                    + (pos_image_ahead - pos_image)*e_max
                end if

            else if ( e_image_ahead > e_image ) then
                tangent = (pos_image_ahead - pos_image)

            else if ( e_image_behind > e_image ) then
                tangent = (pos_image-pos_image_behind)
            end if

            tangent = tangent / norm2(tangent) ! normalising
            spring_tot = -spring_const*(sep_a-sep_b)*tangent
    
            ! Calculates the gradient of the surface perpendicular to the tangents
            call fitness_gradient(pos_image, e_surf, surface)
            surface = surface - tangent * dot_product(surface,tangent)

            ! Total gradient due to the surface and the spring
            ! add the energy of this image to the total
            g(i:i+NDIMENSION-1) = surface + spring_tot
            e = e + e_surf

            ! Calculating the image number and storing the energy for use in the climbing code
            ! For the climbing method we want to check if the energy is the max (i.e close to saddle)
            ! if so we need to calculate a different gradient and over-write the gradient already stored 
            image = (i+NDIMENSION-1) / (NDIMENSION)
            energy( image ) = e_surf

	    if ( (climbing_start .eq. 1) .and. (image .eq. max_image) ) then
            call fitness_gradient(pos_image, e_surf, surface)
            ! this is the gradient we use to climb up the surface to the top of the saddle point
            barrier_pos = pos_image
            g(i:i+NDIMENSION-1) = surface - 2*dot_product( surface, tangent )*tangent
	    end if
        end do

        call fitness_gradient(min_1, e1, g)
        call fitness_gradient(min_2, e2, g)

        barrier = maxval(energy)-min(e1, e2)
        saddle = maxval(energy)
        max_image = maxloc(energy, 1)

    end subroutine gradient_neb

    ! Calculates the initial MEP guess between minimum point 1 and minimum point 2
    ! This is just linear interpolation
    subroutine linear_interp(min_1, min_2, x)
        !! Variables passed in and out
        double precision, intent(out) :: x(N)
        double precision, intent(in) :: min_1(NDIMENSION), min_2(NDIMENSION)

        !! Local variables
        double precision :: delta(NDIMENSION), behind(NDIMENSION)
        integer :: i

        delta = (min_2-min_1) / (IMAGES+1)

        do i=1,N,NDIMENSION
            if (i .eq. 1) then
                behind = min_1
            else
                behind = x(i-NDIMENSION:i-1)
            end if

            x(i:i+NDIMENSION-1) = behind + delta
        end do

    end subroutine linear_interp

    ! Calculates the initial MEP guess between minimum point 1 and minimum point 2
    ! This is just linear interpolation
    subroutine linear_interp_perturb(min_1, min_2, x, perturb)
        !!Variables passed in and out
        double precision, intent(out) :: x(N)
        double precision, intent(in) :: perturb
        double precision, intent(in) :: min_1(NDIMENSION), min_2(NDIMENSION)

        double precision :: rand_num, eps(NDIMENSION)

        !! Local variables
        double precision :: delta(NDIMENSION), behind(NDIMENSION)
        integer :: i, p

        delta = (min_2-min_1) / (IMAGES+1)

        do i=1,N,NDIMENSION
            if (i .eq. 1) then
                behind = min_1
            else
                behind = x(i-NDIMENSION:i-1)
            end if

            !! this generates a random perturbatsion to the linearly interpolated path
            do p=1,NDIMENSION
                call random_number(rand_num)
                eps(p) = 0.5-rand_num
            end do    

            eps = eps / norm2(eps)

            x(i:i+NDIMENSION-1) = behind + delta + eps*perturb
        end do

    end subroutine linear_interp_perturb

 ! Calculates the initial MEP guess between minimum point 1 and minimum point 2
    ! This is just linear interpolation
    subroutine linear_interp_midpoint(min_1, min_2, x, perturb)
        !!Variables passed in and out
        double precision, intent(out) :: x(N)
        double precision, intent(in) :: perturb
        double precision, intent(in) :: min_1(NDIMENSION), min_2(NDIMENSION)

        double precision :: delta(NDIMENSION), behind(NDIMENSION)
        double precision :: rand_num, eps(NDIMENSION), midpoint(NDIMENSION)
        integer :: i, p

        ! Generating the perturbed midpoint as the actual midpoint
        ! plus a perturbation
        do p=1,NDIMENSION
            call random_number(rand_num)
            eps(p) = 0.5-rand_num
        end do    
        eps = eps / norm2(eps) * perturb

        midpoint = min_1 + (min_2 -  min_1)/2
        midpoint = midpoint + eps

        ! interpolating min_1 to the midpoint
        delta = (midpoint-min_1) / ((IMAGES/2)+1)
        do i=1,N/2,NDIMENSION
            if (i .eq. 1) then
                behind = min_1
            else
                behind = x(i-NDIMENSION:i-1)
            end if

            x(i:i+NDIMENSION-1) = behind + delta
        end do

        ! interpolating from midpoint to min2
        delta = (min_2-midpoint) / ( (IMAGES/2) +1)
        do i=N/2+1,N,NDIMENSION
            if (i .eq. 1) then
                behind = min_1
            else
                behind = x(i-NDIMENSION:i-1)
            end if

            x(i:i+NDIMENSION-1) = behind + delta
        end do

    end subroutine linear_interp_midpoint

    subroutine seed_rng()
        ! Seed the rng from /dev/urandom
        integer :: n
        integer, dimension(:), allocatable :: i
        call random_seed(size=n)
        allocate(i(n))
        open(11, file='/dev/urandom', access='stream', form='unformatted')
        read(11) i
        close(11)
        call random_seed(put=i)
    end subroutine seed_rng

end module neb
