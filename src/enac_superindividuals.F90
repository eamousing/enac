module enac_superindividuals

    use enac_commons

    implicit none

contains

    subroutine newsis()
        !! 1. Provides the total number of a year class
        !! 2. Initializes new SIs.
        real :: r1
    
        ! Get this year's total recruitment: r1
        call getnrecruits(r1)
        ! distribute on SIs
        call initialize_sis(r1)
    end subroutine newsis

    subroutine initialize_sis(r1)
        !! Fills in, for each SI that is generated, the sistandard table, with
        !! standard parameters. 
        integer :: i, j, k, regime
        real :: rx, r1, xsi, rancum, rnx, snrn
      
        ! identify SIs that recruit this year, incl. update nsi2: the highest SI number
        nsi2(species) = nsi2(species) + nsiperyear(species)
        nsi1(species) = nsi2(species) - nsiperyear(species) + 1

        ! find values for siparam
        do i = nsi1(species), nsi2(species)
            do j = 1, maxfunc
                do k = 1, maxfuncparam
    
                    if (stdparameters(species,j,k,1) .gt. -999.0) then
                        if (stdparameters(species,j,k,2) .gt. 0.0) then
                            ! rnx is rng that returns random # between 0 and 1
                            rx = rnx(rintvar)
                            ! snrn returns probability of 0-1 number in normal distribuion
                            rancum = snrn(rx)
                            xsi = stdparameters(species,j,k,2)*rancum
                        else
                            xsi = 0.0
                        endif
                    else
                        ! if the stdparam is unused (-999.0) this carries over to siparam
                        xsi = 0.0
                    endif   
                    siparam(species,i,j,k,1) = stdparameters(species,j,k,1) !*exp(xsi)   
                    siparam(species,i,j,k,2) = siparam(species,i,j,k,1) 
          
                enddo
            enddo
      
            ! fill in initial state variables in the table sistate for the new SIs
            ! numbers in SI
            sistate(species,i,1) = r1/float(nsiperyear(species))
          
            ! Age
            sistate(species,i,2) = 0.0

            ! Initial Length
            if (icodeshift(species) .lt. 3 .or. year .le. 0) regime = nregact(species,year)
            if (icodeshift(species) .eq. 3 .and. year .gt. 0) regime = recregime(species)
            rx = rnx(rintvar)
            rancum = snrn(rx)
            xsi = rancum*recruit_parameters(species,regime,18,2)
            sistate(species,i,3) = recruit_parameters(species,regime,17,2)*exp(xsi) 
            
            ! Weight
            sistate(species,i,4) = siparam(species,i,2,1,2)*sistate(species,i,3)**siparam(species,i,2,2,2)
            
            ! Maturity: recruit as immature
            sistate(species,i,5) = 0.0
            sistate(species,i,6) = 0.0
            
            ! Positions, set to 0 for the time being
            sistate(species,i,7) = 0.0
            sistate(species,i,8) = 0.0
        enddo
    end subroutine initialize_sis

    subroutine getnrecruits(r1)
        ! Get the total number of recruits for the current species in the current year in r1
        real :: r0, r1, a, b, sigma, rho, rx, xsi, rnx, snrn
        real :: amplitude, period, phase, xm, trunchigh, trunclow
        integer :: isrfunc, distribution, regime, trunc_type
      
        if (icodeshift(species) .eq. 3 .and. year .gt. 0) then
            regime= recregime(species)
        else
            regime= nregact(species,year) 
        endif
      
        isrfunc = nint(recruit_parameters(species,regime,1,2))
        a = recruit_parameters(species,regime,4,2)
        b = recruit_parameters(species,regime,5,2)
        sigma = recruit_parameters(species,regime,6,2)
        rho = recruit_parameters(species,regime,20,2)
      
        ! Spasmodic on a-parameter
        if (icodeshift(species) .ne. 3 .and. ispasm(species,year) .eq. 1) then
            a = a*recruit_parameters(species,regime,15,2)
        end if
        if (icodeshift(species) .eq. 3 .and. spasm(species,regime,year) .eq. 1) then
            a = a*recruit_parameters(species,regime,15,2)
        end if
        
        ! Periodic:
        amplitude = recruit_parameters(species,regime,9,2)
        period = recruit_parameters(species,regime,10,2)
        phase = recruit_parameters(species,regime,11,2)
        a = a*(1.0 + amplitude*cos(6.28*(float(year) - phase)/period))
       
        if (sigma .le. -999.0) sigma=0.0
      
        if (flag_prime .eq. 1) then
            ! In priming phase, there is no SSB, so we use the a-parameter and a 'priming multiplier'
            ! The a-parameter has been modified by periods and spasms
            ! Periods always apply, but spasms only from the year started as the last spasmodic year, 
            ! which may before year 0
            r0 = a*prime_factor(species)
        else
            ! not priming, use SSB to get r0
            if(isrfunc .eq. 1) then
                ! Hockey stick (R = a above b)
                r0 = a*min(ssb(species,year)/b,1.0)
            else if (isrfunc.eq.2) then
                ! Beverton-Holt
                r0 = a*ssb(species,year)/(b + ssb(species,year))
            else if (isrfunc .eq. 3) then
                ! Ricker - private parametrization
                ! r0=a*ssb(species,year)/b*exp(1.0-ssb(species,year)/b)
      
                ! Ricker - conventional form
                r0 = a*ssb(species,year)*exp(-b*ssb(species,year))
            end if
        end if ! flag_prime yes or no
        ! get random noise - always
10      rx = rnx(rintvar)
        xsi = snrn(rx)
        distribution = nint(recruit_parameters(species,regime,2,2))
      
        ! calculate effect of mackerel predation (J. Trochta)
        ! This version modifies sigma and is based on quantile regression fit (90th)
        ! i.e. key assumption is that mackerel 'limits' recruitment by dampening recruitment variance, not the mean
        if (link_parameters(species,3,1) .eq. 1 .and. year .gt. 0 .and. species .eq. 3) then 
            sigma = effect(species)
        end if
      
        rec_res(species,year) = xsi*sigma
      
        ! Autocorrelation in residuals if rho!=0.0
        if (rho .ne. 0.0 .and. year .gt. -maxage) then
            rec_res(species,year) = rho*rec_res(species,year-1) + &
                sqrt(1 - rho**2)*rec_res(species,year)
        end if
      
        if (distribution .eq. 1) then
            xm = 1.0 + rec_res(species,year)
        else if (distribution .eq. 2) then
            xm = exp(rec_res(species,year) - sigma**2/2.0)
        end if
      
        ! Truncation
        trunc_type = nint(recruit_parameters(species,regime,3,2))
        trunclow = recruit_parameters(species,regime,7,2)
        trunchigh = recruit_parameters(species,regime,8,2)
           
        ! Variant 1: Trunc is lower and upper fraction of sigma (lower negative).
        if (trunc_type .eq. 1) then
            if (xsi .lt. trunclow .or. xsi .gt. trunchigh) then
                goto 10
            endif
        else
        ! Variant 2:  upper and lower bounds on the multiplier xm,
            if (xm .lt. trunclow .or. xm .gt. trunchigh) then
                !if(species.eq.3) print *, xsi,sigma*xsi,sigma**2/2.0,xm,trunclow,trunchigh
                goto 10
            endif
        endif
        ! Final nummber, never negative
      
        r1=max(xm*r0,0.0)       
      
        if (link_parameters(species,3,1) .eq. 1 .and. year .gt. 0 .and. species .eq. 3) then 
            if (year.eq.1) open(67, file='Out/'//'Mac_vs_Her.txt')
                write(67,*) iter,year,macnum,r1,effect(3)
            if(year.eq.nyear) close (67)
        endif
      
        ! calculate mac predation - Kjell's old code commented out 14.11.2021
        ! if (link_parameters(species,3,1).eq.1.and.year.gt.0) then 
        !     r1 = r1 * effect(species)
        ! endif
      
        ! variables for output  
        r1_out(species) = r1
        reg_out(species) = regime
             
        if (species.eq.maxspec.and.year.gt.0) then
            !write(*,*) 'year',year, 'R',r1_out
            write(31,*) iter,year,r1_out,reg_out !,ssb(:,year)
            write(32,*) iter,year,rec_res(:,year),reg_out
        endif
        
1       format(181(f10.1,x)) ! EAM: Unused, remove?
      
    end subroutine getnrecruits

    


end module enac_superindividuals