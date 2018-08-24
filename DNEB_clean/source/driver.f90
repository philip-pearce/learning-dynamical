program main
    use neb

    double precision, allocatable :: minima_array(:) ! This stores all possible minimumm
    double precision jump_amp, similarity, dx, max_rms, spring_const, perturb
    integer max_jumps, iter_max
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: Ideally these will be read in from a .yml style file
    !! at some point, and the images in the neb will be set here
    !! currently that is set in neb.f90
    !!
    !! Parameters used by the alogrithm
    iter_max = 500000
    max_jumps = 10

    !! Finding all the minima
    jump_amp = 2D+0
    similarity = 1D-2

    !! Steepest descent parameters
    dx = 1D-5
    max_rms = 1D-5

    !! Spring const for neb method
    spring_const = 1D+1

    !! perturb value for the linear path
    perturb = 1d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Seed the random number generator
    call seed_rng()

    !! Load the files for the GMM
    call fitness_init()

    !! This calculates all the minima and their corresponding energy
    !! Input: max_jumps, jump_amp, similarity, iter_max, dx, max_rms
    !! Output: minima array
    call find_minima(max_jumps, jump_amp, similarity, iter_max, dx, max_rms, minima_array)

    !! Find all the barriers
    !! Input: minima_array, iter_max, dx, max_rms, spring_const
    call find_paths(minima_array, iter_max, 1D-5, 1D-6, spring_const, perturb)
    
end program main
