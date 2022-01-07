module enac_update_pop

    use enac_commons
    use enac_functions

    implicit none

contains

    subroutine update_population()
        !! Updates all SIs of the current species in the current season
        real :: s, t
    
        if (season .eq. 1) call setstate()

        call grow_and_mature()

        if (season .eq. 1 .and. year .ge. 0) then
            call get_ssb(s)
            ssb(species,year) = s
            call get_tsb(t)
            tsb(species,year) = t  
        endif

        call die()
    end subroutine update_population

    subroutine die()
        !! Reduces all SIs for the species in the year and season,
        !! according to mortalities
        !! It also updates age if the next season is the first in a year
        real  :: length, selection, l50, slope, natmor, cc
        real :: f, z
        integer :: si, age
       
        do si = nsistart(species), nsi2(species)
            ! This is calculating selection by calling logistic directly (i.e. same as get_selection)
            ! length=sistate(species,si,3)
            ! l50=siparam(species,si,4,1,2)
            ! slope=siparam(species,si,4,2,2)
            ! call logistic(l50,slope,length,selection)
          
            if (bio_opt(species,3,1) .eq. 1) then
                call get_selection(si,selection)
            else
                call get_selection_2(si,selection)
            end if
      
            if (bio_opt(species,2,1) .eq. 1) then
                call get_natmor(si,natmor)
            else
                call get_natmor_2(si,natmor)
            end if
          
            ! John:  So he fixes mortality by age anyways for herring..
            ! if(species.eq.3.and.sistate(species,si,2).lt.3) natmor=0.9
            ! if(species.eq.3.and.sistate(species,si,2).gt.2) natmor=0.15 
            
            ! Natmor is on annual scale from the outset
            natmor = natmor / float(ntimestep)
            ! f-level is seasonal from the outset
            f = f_level(species,year,season)*selection
            z = f + natmor
        
            ! The stock rec function for blue whiting is based on 1 year old fish,
            ! therefore we assume no mortality for blue whiting as age 0
            if (species .eq. 1 .and. sistate(species,si,2) .eq. 0) z = 0
      
            sistate(species,si,1) = sistate(species,si,1) * exp(-z)
          
            if (year .gt. 0) then
                age = nint(sistate(species,si,2)) + 1
                natage(species,year,age) = sistate(species,si,1) + natage(species,year,age)
            
                cc = sistate(species,si,1)/exp(-z) * (1.0-exp(-z)) * f/z
                catage(species,year,age) = cc + catage(species,year,age)
            end if
            if (season .eq. ntimestep) sistate(species,si,2) = sistate(species,si,2)+1.0
        enddo ! end SI loop
    end subroutine die

    subroutine grow_and_mature()
        !! Increases length of all SIs for the species in the year and season
        !! Updates weight
        real :: deltat, deltal, lastlen, mod1, mod2, pm, rx, rnx, biomass
        integer :: si, a1, a2, a3
      
        deltat = 1.0 / float(ntimestep)
        do si = nsistart(species), nsi2(species)
            !if(sistate(species,si,2).ge.nage(species)) sistate(species,si,1)=0.0
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
            if (species .eq. 3 .and. sistate(species,si,3) .gt. 32) then
                call vonbert(siparam(species,si,1,1,2),siparam(species,si,1,2,2),sistate(species,si,3),deltat,deltal)
            else
                call vonbert(siparam(species,si,1,1,2),siparam(species,si,1,2,2),sistate(species,si,3),deltat,deltal)
                !call vonbert2(siparam(species,si,1,1,2),siparam(species,si,1,2,2),sistate(species,si,3),deltat,siparam(species,si,1,3,2),deltal) 
            endif
            !if(year.gt.0.and.species.eq.3) print *, species,deltal,siparam(species,si,1,2,2),sistate(species,si,2)
      
            lastlen = sistate(species,si,3) 
            sistate(species,si,3) = lastlen + deltal
            sistate(species,si,4) = siparam(species,si,2,1,2)*sistate(species,si,3)**siparam(species,si,2,2,2)
      
            if (sistate(species,si,5) .lt. 1.0 .and. sistate(species,si,2) .gt. 0.0) then
                if (bio_opt(species,1,1) .eq. 1) then
                    ! If length-based, use logistic to calculate maturity prob at lengths in current and last time step
                    call logistic(siparam(species,si,3,1,2),siparam(species,si,3,2,2),sistate(species,si,3),mod1)
                    call logistic(siparam(species,si,3,1,2),siparam(species,si,3,2,2),lastlen,mod2)
                else if (sistate(species,si,2) .lt. 14.0) then
                    ! If age-based (non-parametric), use linear interpolation to compute mature prob at incremental ages in current and last time step
                    ! i.e. mature prob incrementally increases with monthly progression from this year's to next year's age (age increments by 1/12)
                    a1 = nint(sistate(species,si,2))
                    ! NOTE: potential type conversion problems for a2 and a3 below?
                    a2 = sistate(species,si,2) + float(season)/float(ntimestep)
                    a3 = a2 - 1.0 / float(ntimestep)
                    ! NOTE: introduced type conversion for a2 and a3 below. 
                    ! Check if a1-a3 should be defined as reals instead of integers?
                    call linint(real(a2),float(a1),mat_pattern(species,1,a1,2),float(a1+1),mat_pattern(species,1,a1+1,2),mod1)
                    call linint(real(a3),float(a1),mat_pattern(species,1,a1,2),float(a1+1),mat_pattern(species,1,a1+1,2),mod2)
                else
                    mod1 = 1.0
                    mod2 = 0.0
                endif
            
                pm = (mod1 - mod2) / (1.0 - mod2) 
                rx = rnx(rintvar) 
                if (rx .le. pm) then
                    sistate(species,si,5) = 1.0
                    sistate(species,si,6) = sistate(species,si,2)
                endif
            endif
        enddo
    end subroutine grow_and_mature

    subroutine setstate
        integer :: i,j 
        ! real :: age
        integer :: age
    
        do i = 1, maxspec
            do j = nsistart(i), nsi2(i)
                age = nint(sistate(i,j,2))
                if (age .gt. 0 .and. age .le. 14) then
                    sistate(i,j,3) = mlengths(i,age)
                    sistate(i,j,4) = mweights(i,age)
                endif 
            enddo
        enddo
    end subroutine setstate

end module enac_update_pop