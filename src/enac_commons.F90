module enac_commons
    !! TODO: Add description of module here

    use, intrinsic :: iso_fortran_env, dp => real64

    implicit none

    ! Parameters for dimensioning arrays
    ! TODO: Add description of each parameter using double !!
    integer, parameter :: maxspec = 3 !! Maximum number of species
    integer, parameter :: maxage = 20
    integer, parameter :: maxyear = 1000
    integer, parameter :: maxtimestep = 12
    integer, parameter :: maxsiperyear = 100
    integer, parameter :: maxsi = (maxyear + maxage + 2) * maxsiperyear
    integer, parameter :: maxfunc = 7
    integer, parameter :: maxfuncparam = 6
    integer, parameter :: maxiter = 100
    integer, parameter :: maxlkl = 50
    integer, parameter :: maxlinkfunc = 3
    integer, parameter :: maxlinkfuncparam = 12
    integer, parameter :: maxenvvariables = 3
    integer, parameter :: maxregimes = 3
    integer, parameter :: maxyear2 = 1

    ! Species names
    character(len=50) :: species_name(maxspec)

    ! Actual bounds for loops, must not exceed the max-bounds
    integer :: nspec
    integer :: nyear
    integer :: ntimestep
    integer :: niter
    integer :: nage(maxspec)
    integer :: nsiperyear(maxspec)
    real :: biom(maxspec)
    real :: mlengths(maxspec, 14)
    real :: mweights(maxspec, 14)

    ! Fixed input numbers:
    integer :: recruit_season(maxspec)
    integer :: spawning_season(maxspec)

    ! Current species, year, season etc
    integer :: species
    integer :: year
    integer :: season
    integer :: iter

    ! Book-keeping of SI numbers
    integer :: nsistart(maxspec)
    integer :: nsi1(maxspec)
    integer :: nsi2(maxspec)

    ! Link functions specification by species
    integer :: nlinkfunc(maxspec)
    real :: link_parameters(maxspec, 0:maxlinkfunc, maxlinkfuncparam)

    ! Environment variables as input
    real :: varenv(maxyear, maxtimestep, maxenvvariables, 2)
    real :: varenv1(2, maxyear, 2)

    ! Matrices for the population and its parameters
    ! sistate is the state vector for the SIs. There are 8 state variables
    real :: sistate(maxspec, maxsi, 8)
    real :: sistatekeep(maxspec, maxsi, 8)

    ! For keeping n-values as inout for decmodul, sistate is dumped to 
    ! sistatekeep when requested
    ! siparam are the individul function parameters for each SI
    ! There are 8 functions, each with max 6 function parameters:
    !   1: is the basic value
    !   2: is the current value
    real :: siparam(maxspec, maxsi, maxfunc, maxfuncparam, 2)

    ! stdparameters are the input standard values for the function parameters
    ! the 3rd dimension is 1: Value, 2: sigma
    real :: stdparameters(maxspec, maxfunc, maxfuncparam, 2)

    ! Parameters caused by environment as seen by managers
    real :: decparam(maxspec, maxfunc, maxfuncparam)

    ! The recruit_parameters are by species, regime, parameter number (1-16) and source
    ! The source dimension is: 
    !   1. Standard values
    !   2: modified after exposure to link functions,
    !   3: As assumed by managers
    real :: recruit_parameters(maxspec, maxregimes, 20, 3)

    ! The prime factor is a multiplier for the a-parameter, 
    ! to avoid assuming optimal conditions in the priming phase.
    real :: prime_factor(maxspec)

    ! Regimeshifts
    integer :: nregimes(1:maxspec)
    integer :: icodeshift(1:maxspec)
    integer :: lastreg(1:maxspec)
    integer :: reg_out(1:maxspec)
    integer :: nregact(1:maxspec, -maxage:maxyear)
    integer :: listreg(1:maxspec, 0:maxregimes)
    real :: regintv(1:maxspec)
    real :: r1_out(1:maxspec)

    ! spasmodic recruitments
    integer :: ispasm(1:maxspec, -maxage:maxyear)
    integer :: spasm(maxspec, 3, -maxage:maxyear)

    ! recruitment residuals (for calculating autocorrelated residuals)
    real :: rec_res(1:maxspec, -maxage:maxyear)

    ! F-levels (decided from the HCR implementation)
    real :: f_level(maxspec,-maxage:maxyear,maxtimestep)

    ! Priming F
    real :: fprime(maxspec)

    ! Biomasses
    real :: tsb(maxspec, 0:maxyear)
    real :: ssb(maxspec, 0:maxyear)

    ! Catch
    real  :: tac(maxspec,maxyear)

    ! Management options
    real :: manage_opt(maxspec, 22)
    
    ! Predation from mackerel
    real :: effect(maxspec)

    ! Climate effect
    integer :: recregime(maxspec)

    ! Random seed
    real(dp) :: rintvar

    ! Flags
    integer :: flag_prime

    ! Testing with fix recruitment
    real :: rec_in(45, 3)

    ! For reading in age-specific values for maturity, mortality, and selection
    real :: bio_opt(maxspec, 3, 3)
    real :: mat_pattern(maxspec,maxyear2, 14, 2)
    real :: mor_pattern(maxspec,maxyear2, 14, 2)
    real :: sel_pattern(maxspec,maxyear2, 14, 2)

    ! Natage matrices (species,year,iteration,ages)
    real :: natage(maxspec, maxyear, maxage+1)
    real :: catage(maxspec, maxyear, maxage+1)

    ! Abundance of Mackerel ages 2+ - for plotting effects...
    real :: macnum

    ! Other by Kjell. TODO: make descriptive comment
    real :: actdev(maxspec, maxyear)

end module enac_commons