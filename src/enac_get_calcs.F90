module enac_get_calcs

    use enac_commons
    use enac_functions

    implicit none

contains

    subroutine get_yield(f,y)
        !! Calculates the yield corresponding to a level of F, for the species and year and season in the calling context
        !! Basically, this is a seasonal F, and the catch is the seasonal catch
        real, intent(in) :: f
        real, intent(out) :: y
        real :: selection,rn,m,c,sel,tot
        integer :: si

        y = 0.0
        sel = 0.0
        tot = 0.0

        !! Yield from each SI
        do si = nsistart(species), nsi2(species)
            if (bio_opt(species,3,1) .eq. 1) then
                call get_selection(si,selection)
            else
                call get_selection_2(si,selection)
            endif

            rn = sistate(species,si,1)

            if (bio_opt(species,2,1) .eq. 1) then
                call get_natmor(si,m)
            else
                call get_natmor_2(si,m)
            endif
        
            ! if (species.eq.1.and.si.eq.nsi2(1)) then  
            !   call get_selection_2(si,selection)
            !   write(*,*) 'age-based selection: species 1  SI',si, 'sel',selection
            !   call get_selection(si,selection)
            !   write(*,*) 'length-based selection: species 1  SI',si, 'sel',selection
            ! endif

            sel = sel + selection
            tot = tot + 1.0
            m = m / float(ntimestep)
            c = rn * f * selection * (1.0 - exp( -m - f * selection)) / (f * selection + m)
            y = y + c * sistate(species,si,4)
        enddo
        !write(*,*) 'full selections',sel, 'tot SI',tot,'prop_sel',sel/tot
        !write(*,*) 'nsistart',nsistart, 'nsi2',nsi2
        !write(*,*) 'season',season, 'species',species,'yield',y
    end subroutine get_yield

    subroutine get_selection(si,selection)
        integer, intent(in) :: si
        real, intent(out) :: selection
        real :: length,l50,slope

        length = sistate(species,si,3)
        l50 = siparam(species,si,4,1,2)
        slope = siparam(species,si,4,2,2)
        
        call logistic(l50, slope, length, selection)
    end subroutine get_selection

    subroutine get_selection_2(si,selection)
        integer, intent(in) :: si
        real, intent(out) :: selection
        integer :: a1,a2
        real :: length,rx,rnx
        real, dimension(14) :: si_sel_len,si_sel_val
      
        length = sistate(species,si,3)
        !! If options are not set correctly, no selection so no fishing
        selection = 0.
       
        if (bio_opt(species,3,1) .eq. 2) then
            a1 = nint(sistate(species,si,2))
            !! Year is indexed at 1 - will eventually allow option for resampling...
            if (a1 .gt. 0 .and. a1 .le. 14) selection = sel_pattern(species,1,a1,2)
            if (a1 .gt. 14) selection = sel_pattern(species,1,14,2)
        else if (bio_opt(species,3,1) .eq. 3) then
            do a1 = 1, 14
                si_sel_len(a1) = sel_pattern(species,1,a1,1)
                si_sel_val(a1) = sel_pattern(species,1,a1,2)
            enddo
        
            !! Linearly interpolate selection for length
            !! Loop thru each mean length in si_len, check to see if between mean lengths, then linearly interpolate
            do a2 = 2, 14
                if (length .le. si_sel_len(1)) then
                    selection = 0.0
                else if (length .gt. si_sel_len(14)) then
                    selection = si_sel_val(14)
                else if (length .gt. si_sel_len(a2 - 1) .and. length .le. si_sel_len(a2)) then
                    ! Use continuation character '&' because this line is too long for FORTRAN (spits out error)
                    call linint(length, si_sel_len(a2 - 1), si_sel_val(a2 - 1), si_sel_len(a2), si_sel_val(a2), selection)
                end if
            enddo
        
            ! If linear interpolation spits out >1., fix to 1.
            if (selection .gt. 1.0) selection = 1.0
        endif
    end subroutine get_selection_2

    subroutine get_natmor(si,m)
        !! It is assumed that the natural mortality is on an annual scale…
        integer, intent(in) :: si
        real, intent(out) :: m
        real :: length,mlow,mhigh,slope1,slope2,l1_50,l2_50

        length = sistate(species,si,3)
        mlow = siparam(species,si,5,1,1)
        mhigh = siparam(species,si,5,2,1)
        l1_50 = siparam(species,si,5,3,1)
        l2_50 = siparam(species,si,5,4,1)
        slope1 = siparam(species,si,5,5,1)
        slope2 = siparam(species,si,5,6,1)

        call comblog(mlow,mhigh,l1_50,slope1,l2_50,slope2,length,m)
    end subroutine get_natmor

    subroutine get_natmor_2(si,m)
        !! This calculates natural mortality at age, rather than length
        integer, intent(in) :: si
        real, intent(out) :: m
        integer :: a1
        
        !! Initialize m (m for age 0 fish, which is equal to age 1)
        m = mor_pattern(species,1,1,2)

        !! Takes age from state array
        a1 = nint(sistate(species,si,2))

        !!  Year is indexed at 1 - hope to allow option for resampling...
        if (a1 .gt. 0 .and. a1 .lt. 14) m = mor_pattern(species,1,a1,2)

        if (a1 .ge. 14) m = mor_pattern(species,1,14,2)
    end subroutine get_natmor_2

    subroutine get_ssb(s)
        real, intent(out) :: s
        integer :: si
       
        s = 0.0
        do si = nsistart(species), nsi2(species)
            s = s + sistate(species,si,1) * sistate(species,si,4) * sistate(species,si,5) 
        end do
    end subroutine get_ssb

    subroutine get_tsb(ts)
        !! Get total biomass in the current state
        !! Note: In the decision model, the speclen holds lower limits for lengths to be included in TSB
        !! So far, they are hard coded.
        !! At some stage, they should be part of a common information base
        !! The routine here gives the total over all lengths for the time being, so it is not quite comparable
        real, intent(out) :: ts
        integer :: si
        
        ts = 0.0
        do si = nsistart(species), nsi2(species)
            ts = ts + sistate(species,si,1) * sistate(species,si,4) 
        enddo
    end subroutine get_tsb

    subroutine get_f_from_yield(y,icflag)
        !! Converts TAC to F-values for the true stock in the calling season
        !! The annual tac is split into seasonal tacs, for the time being assuming
        !! equal quantities every time step, but this may be made seasonal later on.
        !! The transferred y is supposed to be a seasonal tac.
        !! The selection is so far the fixed selection, calculated in subroutine get_yield from paramin.
        !! Flevel is the result, gets filled in here and is used by die in each time step
        real, intent(in) :: y
        !real, intent(out) :: f
        integer, intent(out) :: icflag
        integer :: iround
        real :: tacnow, fnow, ytest, ynow, crit
        icflag = 0     
        tacnow = y
      
        if (tacnow .gt. 0.0) then
            !! Check that there is enough for the tac, i.e. the tac is lower than the catch 
            !! corresponding to an annual F = 2.0
            fnow = 2.0 / float(ntimestep)

            call get_yield(fnow, ytest)
            
            if (ytest .lt. tacnow) then
                !! Not enough, so the fishery is limited by F=2.0 per year
                icflag = 1
                f_level(species,year,season) = 0.0
            else
                !! Enough fish, so implement searching routine to get the flevel that gives tacnow
                !! The flevels are all scaled as annual F
                !! start with fprime just to start at a sensible level
                fnow = fprime(species) / float(ntimestep)
                iround = 0
                
                !! JTT: search_loop replaces GOTO statement below - remove GOTO once tested
                search_loop: do iround = 0, 20
                    call get_yield(fnow, ynow)
                    crit = abs(1.0 - ynow / tacnow) 
                    if (crit .gt. 0.0001) then
                        if (iround .eq. 20) then
                            write(*,*) 'No convergence in search for F'
                            fnow = (fnow + fnow * tacnow / ynow) / 2.0
                            write(*,*) 'Search broken'
                        else
                            fnow = fnow * tacnow / ynow
                        end if
                    else
                        exit search_loop
                    end if
                end do search_loop

                !! 10 call get_yield(fnow, ynow)
                !! crit = abs(1.0 - ynow / tacnow) 
                !! if (crit .gt. 0.0001) then
                !!     iround = iround + 1
                !!     if (iround .gt. 20) then
                !!         write(*,*) 'No convergence in search for F'
                !!         fnow = (fnow + fnow * tacnow / ynow) / 2.0
                !!         write(*,*) 'Search broken'
                !!     else
                !!         fnow = fnow * tacnow / ynow
                !!         goto 10
                !!     end if  
                !! end if

                f_level(species,year,season) = fnow
            end if ! end if there is enough fish

        else if (tac(species,year) .eq. 0.0) then
          !! if tac=0, just set f to 0
          f_level(species,year,season) = 0.0
        else
          !! tac set negative: undefined (in priming)
          f_level(species,year,season) = fprime(species) / float(ntimestep)
        end if
    end subroutine get_f_from_yield

end module enac_get_calcs