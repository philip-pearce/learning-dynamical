program main
    use neb

    double precision, allocatable :: minima_array(:) ! This stores all possible minimumm
    double precision jump_amp, similarity, dx, max_rms, spring_const
    integer max_jumps, iter_max
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: Ideally these will be read in from a .yml style file
    !! at some point, and the images in the neb will be set here
    !! currently that is set in neb.f90
    !!
    !! Parameters used by the alogrithm (flag - add 5 zeros to this!!)
    iter_max = 100000000
    max_jumps = 39

    !! Finding all the minima 1D+0
    jump_amp = 1D+0
    similarity = 1D-2

    !! Steepest descent parameters
    dx = 1D-4
    max_rms = 1D-9

    !! Spring const for neb method
    spring_const = 10D+0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Seed the random number generator
    call seed_rng()

    !! Load the files for the GMM
    call fitness_init()

    !! This calculates all the minima and their corresponding energy
    !! Input: max_jumps, jump_amp, similarity, iter_max, dx, max_rms
    !! Output: minima array
    call find_minima(max_jumps, jump_amp, similarity, iter_max, dx, max_rms, minima_array)


    call write_minima(minima_array, NDIMENSION)
    !call read_minima(minima_array, NDIMENSION)

    !! Find all the barriers
    !! Input: minima_array, iter_max, dx, max_rms, spring_const
    !call find_paths(minima_array, iter_max, dx, max_rms, spring_const)
    
end program main
