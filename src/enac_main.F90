program enac
    use enac_commons
    use enac_functions
    use enac_inout

    use enac_get_calcs
    use enac_links
    use enac_superindividuals
    use enac_update_pop
    use enac_decision_model

    implicit none

    !! Global and internal procedure variables
    real :: catch, f, r1, ts, y
    integer :: i, k, icflag, ntacyr, assesserror
    real, dimension(1:maxspec, 0:maxtacyr) :: dectac

    call setframe
    call popinput
    call recinput
    call linkinput
    call agebasedinput
    call envinput
    call ruleinput
    call statein

    !! Create files for output
    !! JTT: A better way to do this? More descriptive file names?
    open (13,file='out/'//'SSB.txt')
    open (11,file='out/'//'TAC.txt')
    open (14,file='out/'//'NAA.txt')
    write(14,*) 'iter',char(9),'year',char(9),'species',char(9),'age',char(9),'number',char(9),'catch'
    open (15,file='out/'//'Sistate_bw.txt')
    open (16,file='out/'//'Sistate_mac.txt')
    open (17,file='out/'//'Sistate_her.txt')
    open (23,file='out/'//'Mort_Mat.txt')
    open (30,file='out/'//'SSB_decmodul.txt')
    open (31,file='out/'//'Recruitment_r1.txt')
    write(31,*) 'iter',char(9),'year',char(9),'blu_rec',char(9),'mac_rec',char(9), &
                'her_rec',char(9),'blu_regime',char(9),'mac_regime',char(9),'her_regime'
    open (32,file='out/'//'Rec_residuals.txt')
    write(32,*) 'iter',char(9),'year',char(9),'blu_rec_resid',char(9),'mac_rec_resid', &
                char(9),'her_rec_resid',char(9),'blu_regime',char(9),'mac_regime',char(9),'her_regime'
    open (41,file='out/'//'Weight_kol.txt')
    open (42,file='out/'//'Weight_mac.txt')
    open (43,file='out/'//'Weight_her.txt')
    open (44,file='out/'//'Length_kol.txt')
    open (45,file='out/'//'Length_mac.txt')
    open (46,file='out/'//'Length_her.txt')
    ! open (68,file='out/'//'rng_sequence.txt')

    ! Open the file with the vectors of random assessment error from a multivariate normal distribution
    open (47,file='in/blue_whiting_rndmvnormval_with_R.txt')
    read (47,*) rndmvnorm_kol
    open (48,file='in/mackerel_rndmvnormval_with_R.txt')
    read (48,*) rndmvnorm_mac
    open (49,file='in/herring_rndmvnormval_with_R.txt')
    read (49,*) rndmvnorm_her
  
    do species=1,maxspec
      if(species.eq.1) then
        rndmvnorm(species,1:1000,1:maxage) =rndmvnorm_kol
      else if(species.eq.2) then
          rndmvnorm(species,1:1000,1:maxage) =rndmvnorm_mac
        else
          rndmvnorm(species,1:1000,1:maxage) =rndmvnorm_her
        endif
    enddo

    !! Set random seed. JTT: Maybe done outside main program, to be read in?
    rintvar = 15349668.0_8
    ! rng_count = 0

    ! Set 1  or 0 (on and off respectively) the variable assesserror to apply or not random error to assessment (to the observations from the OM)
    assesserror = 0 

    !! Simulation loop of MSE
    do iter = 1, niter

        write(*,*) 'Iteration ',iter,'of',niter 
    
        !! Initialization.
        !! Start in the year -maxage.
        !! nsistart and nsi2 are the initial values for Superindividual (SI) numbers. 
        !! nsistart is the index for the first SI in the population, it is later on cut at twice max age (JTT: ?).
        !! Actual numbers are set in initialize_sis where new SIs get numbers
        do species = 1, nspec
            nsistart = 1
            nsi2(species) = 0
            year = -maxage
        end do
    
        !! Initialize f_level  
        f_level(:,:,:) = -1
    
        !! Initialize numbers at age
        natage(:,:,:) = 0.0
        catage(:,:,:) = 0.0
    
        !! setup_recruits outputs table of when abnormal recruitment years occur for each species
        call setup_recruits
        
        !! Build population to starting year (1) for MSE (i.e. priming)
        !! Years -maxage:0 are priming, and turned on by flag_prime
        !! The starting reference TAC for year 1 is derived from the catch at priming F
        !! and set in year 0
        flag_prime = 1
   
        !! Year loop, including priming
        do year = -maxage, nyear
    
            write(*,*) 'Year', year
            
            !! Once normal years of MSE start,
            !! Turn off priming and increment the indices of each year's cohort of new SI (by 100). 
            if (year .gt. 0) then
                nsistart = nsistart + 100  
                flag_prime = 0
            end if
        
            !! Season loop (Current monthly at 12 steps)
            do season = 1, ntimestep
                
                !! Species loop (3) 
                do species = 1, nspec
                    !! Step 1. Get the F-level and catch (for y>0) in each season and update parameters (JTT: ?)
                    if (year .le. 0) then
                        !! In priming years, fixed annual F is divided by number of seasons
                        !! fprime is set in In/framework file and read in by enac_inout
                        f_level(species,year,season) = fprime(species) / float(ntimestep)
                    else 
                        !! Ordinary years
                        ! Note: update_parameters and decis_parameters only are effective when year >0
                        if (season .eq. 1) call update_parameters ! in links.f90:  itself calls linkfunc
            
                        !! Set F and catch unless it is already in place
                        if (f_level(species,year,season) .gt. -0.00001 .and. &
                        f_level(species,year,season) .lt. 0.00001) then
                            !! Case where we have run out of fish. F-level is already defined (see below) and set to 0 
                            y = 0.0
                        else if (f_level(species,year,season) .lt. -0.999) then
                            !! Case where no F defined yet, but there should be a TAC
                            y = tac(species,year) / float(ntimestep) 
                            
                            ! write(*,*) 'TAC',tac(species,year)    
                            call get_f_from_yield(y,icflag)         
                            
                            !! icflag signals that there were not enough fish.
                            if (icflag .eq. 1) then
                                !! In this case, f_level is set to 0 for the rest of the year 
                                !! and the TAC is adjusted to the fraction of original TAC up until current season
                                tac(species,year) = tac(species,year) * float((season - 1)) / float(ntimestep)            
                
                                !! Fill in rest of year with 0 fishing  
                                do k = season, ntimestep
                                    f_level(species,year,k) = 0.0
                                end do      
                            end if
                        end if
                    end if
            
                    !! Step 2. Recruitment
                    if (season .eq. recruit_season(species)) then      
                        !! superindividuals module: newsis calls getnrecuits and initialize_sis
                        call newsis(r1)
                        ! write(*,*) 'Year= ', year, ' Species= ', species, ' Recruits= ', r1
                    end if
            
                    !! Step 3. Get TSB & SSB if start of the year 
                    if (season .eq. 1 .and. year .ge. 0) then
                        !! enac_get_calcs module
                        call get_tsb(ts)
                        tsb(species,year) = ts
                        ! write(*,*) 'Year= ', year, ' Species= ', species, ' tsb= ', ts
                    end if
                    !! Get SSB 
                    if (season .eq. spawning_season(species) .and. year .ge. 0) then
                        !! enac_get_calcs module
                        call get_ssb(ts)
                        ssb(species,year) = ts
                    end if
            
                    !! Step 4. Reference season for decision, dump seastate to seastatekeep
                    ! Note that this season must be before spawning, therefore it is forced to be 1.
                    if (season .eq. manage_opt(species,1)) then
                        !! From enac_links module
                        call decis_param ! links.f90
                        do i = nsistart(species), nsi2(species) 
                            do k = 1, 8 
                                sistatekeep(species,i,k) = sistate(species,i,k)
                            end do
                        end do
                    end if
                end do
        
                !! Step 5. Decide TAC for next year (done in season 1)
                if (year .eq. 0) then
                    !! The TAC for year 1 is set according to fprime at decision time in year 0
                    if (season .eq. 1) then
                        do species = 1, maxspec
                            f = fprime(species) / float(ntimestep)
                            
                            !! enac_get_calcs
                            call get_yield(f,catch)

                            !! Set the TAC for year 1 from yield in year 0 (according to priming F) and a tacfactor (JTT: ?)
                            tac(species,1) = catch * float(ntimestep)
                        end do     
                    end if
                else if (year .gt. 0 .and. year .lt. maxyear) then
                    !! TAC for next year is set according to decision model
                    if(season .eq. 1) then
                        !! enac_decision_model 
                        !! Changed by Alfonso: the rndmvnorm had to be included also in the arguments for decision_model(), 
                        !! it contains the vectors with the random multivariate error
                        call decision_model(maxtacyr,maxspec,maxtimestep,maxlkl,dectac,ntacyr,rndmvnorm, assesserror)
                        !! call decision_model(maxtacyr,maxspec,maxtimestep,maxlkl,dectac,ntacyr)
            
                        !! ntacyr is the number of years filled in dectac - it is an internal parameter in decision_model
                        do species = 1, maxspec
                            tac(species,year + 1) = dectac(species,1)
                        end do
                    end if
                end if
        
                !! Step 6. Update population (grow, mature and die)
                !! These processes have annual values and are calculated for each SI
                do species=1,maxspec
                    call update_population ! update_population.f90
                end do !species loop
            end do ! season loop    

            call outstate 
            if (year .gt. 0) then
                call setstate 
            end if

            call sistateout
            !! Proceed to next year
        end do
        
        !! Write output
        call mort_mat_out
        call output
        !! Proceed to next iteration
    end do

    !! Close file connections
    close (11)
    close (13)
    close (14)
    close (15)
    close (16)
    close (17)
    ! close (20)
    ! close (21)
    ! close (22)
    
    write(*,*) 'Finished'

end program enac