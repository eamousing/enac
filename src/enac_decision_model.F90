module enac_decision_model
    !! This module contains the procedures for running the decision model,
    !! which includes the assessment model (incl. observation and estimation error),
    !! forecast model, and management model (harvest control rule)

    !! The variables in this module are stand-alone, without the use of those in enac_commons module,
    !! All subroutines have names starting with d, to keep them apart form similar
    !! subroutines in the main part.
    
    ! The subroutine decision_model sets up what was the input to the decomdul.
    ! It extracts the necessary information from common-variables, and adds some information
    ! that at present is hard-coded.
    ! Then it calls the decmodul itself.
    ! The end product from the decomodul is a TAC in dectac(species,year).

    use enac_commons
    use enac_functions
    
    implicit none

contains

  subroutine decision_model(maxspec, maxtimestep, maxlkl, dectac, ntacyr)
    ! Note: the parameter maxyr sets the outer dimension of the years.
    integer, parameter :: maxyr = 2

    integer, intent(in) :: maxspec, maxtimestep, maxlkl
    integer, intent(out) :: ntacyr
    real, dimension (1:maxspec, 0:maxyr), intent(out) :: dectac

    ! transferred variables
    ! These were transferred to the original decmodul. Here, they are
    ! derived from the common variables, and then tansferred to the present decmodul
    ! Data shall be available for the decision process for all species, years, seasons
    ! not only the species, year and season for which the decision is to be made
    
    ! Note that the years (iy) here refer to the internal years in decompodul. The calling
    ! year is year 0 and the tac year is year 1. Normally, maxyr is 2, since it sometimes
    ! is needed to project the stock into the year after the TAC is taken
  
    integer :: mart, mseas, maxl
    real, dimension (1:maxspec, 1:4) :: speclen  !cols: interval, min size class, max size class, min size for TSB
    real, dimension (1:maxspec, 0:maxyr,10) :: bioparam1  !3rd dim: for each parameter defined below 
    real, dimension (1:maxspec, 0:maxyr, 1:maxtimestep, 1:4, 2) :: bioparam2 !4th dim: parameter, 5th dim: parameter specific to something (juvenile and adult K)
    real, dimension (1:maxspec, 0:maxyr, 1:maxtimestep, maxlkl, 3) :: bioparam3
    real, dimension (1:maxspec, 1:maxlkl, 2) :: rninit ! rninit(,,1) is initial stock numbers, rinit(,,2) is sigma for obs noise on these numbers
    real, dimension (1:maxspec, 1:maxtimestep) :: seastac
    real, dimension (1:maxspec, 22) :: options
    integer, dimension (1:maxspec) :: regim
  
    ! Local variables
    integer :: i, iart, iy, iseas, lkl, lklminus, lklplus, si, regime
    real :: length, mlow, mhigh, slope1, slope2, l1_50, l2_50, natmor

    ! Species currently defined in framework:
    ! 1: kolmule
    ! 2: makrell
    ! 3: sild  
  
    ! to avoid confusion with common-variables:
    ! decmodul only sees the year from the current year until maxyr, normally maxyr =2, 
    ! to allow rules that look inot the year after the tac-year
    mart = maxspec
    mseas = maxtimestep
    maxl = maxlkl ! from enac_commons, only used here
  
    ! The standard values used so far are hard-coded. Should be coordinated with the main program.
    ! Interval for length classes - 1 cm everywhere
    speclen(1,1) = 1.0      
    speclen(2,1) = 1.0      
    speclen(3,1) = 1.0
    ! Smallest length class (minus-class)
    speclen(1,2) = 13.0
    speclen(2,2) = 21.0
    speclen(3,2) = 10.0
    ! Largest length class (plus-class)
    speclen(1,3) = 37.0
    speclen(2,3) = 41.0
    speclen(3,3) = 38.0
    ! Lowest length class to be included in TSB
    speclen(1,4) = 15.0
    speclen(2,4) = 15.0
    speclen(3,4) = 15.0
  
    ! bioparam1
    do iart = 1, mart
        do iy = 0, maxyr
            regime = nregact(iart,iy) ! GLOBAL
            ! Linf
            bioparam1(iart,iy,1) = stdparameters(iart,1,1,1) ! GLOBAL
            ! Recruiting length
            bioparam1(iart,iy,2) = recruit_parameters(iart,regime,17,3) ! GLOBAL
            ! CV recruiting length
            bioparam1(iart,iy,3) = recruit_parameters(iart,regime,18,3) ! GLOBAL
            ! Spawning season
            bioparam1(iart,iy,4) = spawning_season(iart) ! GLOBAL
            ! Recruiting season
            bioparam1(iart,iy,5) = recruit_season(iart) ! GLOBAL
            ! L50 in maturation logistic function
            bioparam1(iart,iy,6) = decparam(iart,3,1) ! GLOBAL
            ! Slope in maturation logistic function
            bioparam1(iart,iy,7) = decparam(iart,3,2) ! GLOBAL
            ! Code for recruitment function
            bioparam1(iart,iy,8) = recruit_parameters(iart,regime,1,1) ! GLOBAL
            ! Note on recruit-parameters: Decmodul will only see the assumed a and b-parameters,
            ! not the spasmodic or periodic variants that go into the real recruitment
            ! But it will know the regime in the calling year    
            ! Recruitment a-parameter
            bioparam1(iart,iy,9) = recruit_parameters(iart,regime,4,3) ! GLOBAL
            ! Recruitment b-parameter
            bioparam1(iart,iy,10) = recruit_parameters(iart,regime,5,3) ! GLOBAL
        end do
    end do
  
    ! bioparam2
    do iart=1, mart
      do iy=0, maxyr
        do iseas=1, mseas
          ! k-values
          bioparam2(iart,iy,iseas,1,1) = decparam(iart,1,2) !Adult K ! GLOBAL
          !if(iart.eq.3) print *, bioparam2(iart,iy,iseas,1,1)
          bioparam2(iart,iy,iseas,1,2) = decparam(iart,1,4) !Juvenile K ! GLOBAL
          ! CV growth increment (arbitrary value - not used in the main program
          bioparam2(iart,iy,iseas,2,1) = decparam(iart,1,3) ! GLOBAL
          ! L50 in selection
          bioparam2(iart,iy,iseas,3,1) = decparam(iart,4,1) ! GLOBAL
          ! slope in selection
          bioparam2(iart,iy,iseas,4,1) = decparam(iart,4,2) ! GLOBAL
        end do
      end do
    end do
  
    ! bioparam3
    do iart = 1, mart
      do iy = 0, maxyr
        do iseas = 1, mseas      
          do lkl = 1, maxl
            ! Natural mortality per year
            ! One cm length classes tacitly assumed
            length = float(lkl)+0.5
            mlow = stdparameters(iart,5,1,1) ! GLOBAL
            mhigh = stdparameters(iart,5,2,1) ! GLOBAL
            l1_50 = stdparameters(iart,5,3,1) ! GLOBAL
            l2_50 = stdparameters(iart,5,4,1) ! GLOBAL
            slope1 = stdparameters(iart,5,5,1) ! GLOBAL
            slope2 = stdparameters(iart,5,6,1) ! GLOBAL
            
            call comblog(mlow, mhigh, l1_50,slope1,l2_50,slope2,length,natmor)
            
            ! Natmor is scaled as annual mortalities
            ! Decmodul assumes natural mortality evenly distributed on seasons
            bioparam3(iart,iy,iseas,lkl,1) = natmor
            ! Condition according to environment 1, a      
            bioparam3(iart,iy,iseas,lkl,2) = decparam(iart,2,1) ! GLOBAL
            ! Condition according to environment 2, b      
            bioparam3(iart,iy,iseas,lkl,3) = decparam(iart,2,2) ! GLOBAL
          end do
        end do
      end do
    end do
  
    ! Initial numbers from sistate  
    do iart = 1, mart
      ! Initialize      
      do lkl = 1, maxl
        rninit(iart,lkl,1) = 0.0
        ! rninit(.,.,2) is used to make noise to the N-values in the observation model.
        ! Here there is a need for extensions, e.g. length dependence, annual effects, etc.
        ! and it should not be hardcoded.
        rninit(iart,lkl,2) = 0.30
      end do
      
      lklminus = int(speclen(iart,2))
      lklplus = int(speclen(iart,3))
      
      ! Loop through each SI 
      do si = nsistart(iart), nsi2(iart)   ! GLOBAL 
        lkl = int(sistatekeep(iart,si,3)) ! GLOBAL
   
        if (lkl .lt. lklminus) then
          rninit(iart,lklminus,1) = rninit(iart,lklminus,1) + sistatekeep(iart,si,1) ! GLOBAL
        else if (lkl .gt. lklplus) then
          rninit(iart,lklplus,1) = rninit(iart,lklplus,1) + sistatekeep(iart,si,1) ! GLOBAL
        else
          rninit(iart,lkl,1) = rninit(iart,lkl,1) + sistatekeep(iart,si,1) ! GLOBAL
        end if
      end do
    end do
  
  ! seastac - hard-coded (seasonal distribution of TAC)
    do iart = 1, mart ! Loop through each species
      do iseas = 1, mseas ! Loop through each season
        seastac(iart,iseas) = 1.0 / float(mseas)
      end do
    end do
  
    do iart = 1, mart
      do i = 1, 22
        options(iart,i) = manage_opt(iart,i) ! col 3 from optinn ! GLOBAL
      end do

      ! Options 1 and 2 are exceptions:
      ! Option (.,1) is set to 1 in manage_opt in the main program, and it is transferred here.   
      ! Option (.,2) tac in the call year is set here to what is already determined as tac in year 0
      options(iart,2) = tac(iart,year) ! GLOBAL
    end do  
  
    regim = recregime  ! GLOBAL
    
    call decmodul(mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2,bioparam3,rninit,options,seastac,rintvar,dectac,regim) ! GLOBAL (rintvar)
    ntacyr = maxyr
  
  end subroutine decision_model

  subroutine decmodul(mart, maxyr, mseas, maxl, speclen, bioparam1, bioparam2, & 
    bioparam3, rninit, options, seastac, rintvar, dectac, regim)
    ! All subroutines have the prefix d to avoid confusion with similar routines in the population model  
    ! Some variables have names similar to enac_commons, but commons are not visible here    
    ! These are transferred variables
   
    integer, intent(in) :: mart, maxyr, mseas, maxl
    real, dimension (1:mart, 1:4), intent(in) :: speclen
    real, dimension (1:mart, 0:maxyr, 10), intent(in) :: bioparam1
    real, dimension (1:mart, 0:maxyr, 1:mseas, 1:4, 2), intent(in) :: bioparam2
    real, dimension (1:mart, 0:maxyr, 1:mseas, maxl, 3), intent(in) :: bioparam3
    real, dimension (1:mart, 1:maxl, 2), intent(in) :: rninit
    real, dimension (1:mart, 22), intent(in) :: options
    integer, dimension (1:mart), intent(in) :: regim

    real, dimension (1:mart, 1:mseas), intent(inout) :: seastac

    real, dimension (mart, 0:maxyr), intent(out) :: dectac
   
    double precision, intent(in) :: rintvar

    ! Local variables som utveksles
    real, dimension(mart, 0:maxyr, mseas) ::  tactemp
    real, dimension (mart, 0:maxyr, mseas, maxl) :: rn
    real, dimension (mart, 0:maxyr, mseas) :: tacprev
    real, dimension (mart, 2) :: decbasis
    real, dimension (mart) :: v
   
    ! Local temporary variables
    integer :: iart, iy, lkl, iseas, istart, nriter
    real :: summ, ssq
   
    ! Initialize rn and tactemp, and build bioparam3 (natural mortality)
    do iart = 1, mart
      do iy = 0, maxyr
        do iseas = 1, mseas
          tactemp(iart,iy,iseas) = -1.0
          do lkl = 1, maxl
            rn(iart,iy,iseas,lkl) = 0.0
          end do
        end do
      end do
    end do
   
    ! This part of code is useful for when seasonal distribution of TAC is not uniform
    ! Currently every season has same proportion of annual TAC, but this could not be the case
    ! E.g. a bell curve (that does not have to integrate to 1)
    do iart = 1, mart
      ! Compute denominator for normalization
      summ = 0.0
      do iseas = 1, mseas
        summ = summ + seastac(iart,iseas)
      end do
      
      ! Normalize seastac (relative tac in season).
      if (summ .le. 0.0) write(*,*) 'Warning: Missing seasonal specification of TAC'
      do iseas=1,mseas
        seastac(iart,iseas) = seastac(iart,iseas) / summ
      end do
    end do
   
    ! Fill in tactemp (copy to tacprev) with the TAC for year 0 as preliminary
    do iart = 1, mart
      ! dectac is updated with the previously decided tac for year 0
      dectac(iart,0) = options(iart,2)
      do iy = 0, maxyr
        if (iy .eq. 0) then
          istart = int(options(iart,1)) ! Point in time for initial stock
        else
          istart = 1
        end if
        
        do iseas = istart, mseas
          tactemp(iart,iy,iseas) = options(iart,2) * seastac(iart,iseas) ! Scale the annual TAC to season now
          tacprev(iart,iy,iseas) = tactemp(iart,iy,iseas)
        end do
      end do
    end do
   
    ! Fill in rn with obs-model values for first year and season, add noise to rn-values
    call dobsmodel(rn,mart,maxyr,mseas,maxl,speclen,rninit,options,rintvar)
   
    ! Iteration loop
    do nriter = 0, 20
      ! Project the stock forwards with the current TAC
      ! Will set tactemp=0.0 if not enough fish for catching; simulates within season closure (depletflag==1)
      ! depletflag assigned in ddie based on high F criteria
      call dproject(rn,tactemp,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2,bioparam3,options)
     
      ! Derive the decision basis
      call dgetdecbasis(rn,decbasis,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam3,options)
     
      ! Apply the rule, get the exploitation measure v from the decision basis for all species
      call drule(mart,decbasis,v,options,regim)
     
      ! Translate the v to TAC for all years and species
      call dtranslate(rn,tactemp,v,mart,maxyr,mseas,maxl,speclen,bioparam2,bioparam3,options,seastac)
     
      ! Modify the TACs (filters, constraints)
      ! Fill in for future years
      call dmodifytac(tactemp,mart,maxyr,mseas,options)

      ! Compare with previous (tactemp vs. tacprev), in case the decided TAC changes the decision basis
      ssq = 0.0

      do iart = 1, mart
        do iy = 1, maxyr
          do iseas = 1, mseas
            if (tacprev(iart,iy,iseas) .gt. 0.0 .and. tactemp(iart,iy,iseas) .gt. 0.0) then
              ssq = ssq + log(tactemp(iart,iy,iseas) / tacprev(iart,iy,iseas)) ** 2  
            else
              ! This means that if the tac is set to 0 in the decision process, the loop is broken, and the 
              ! zero tac is returned
              ssq = 0.0
            end if
          end do
        end do
      end do
     
      ! If not similar, iterate
      if (ssq .gt. 0.0001) then
        if (nriter .eq. 20) then
          ! This happens if the iteration becomes cyclic.
          ! The simple solution is to set the value mid-way between current and previous.
          ! Normally, that should be quite close, but a warning is given.
          write(*,*) 'Warning: Line search in decmodul broken at iteration nr. ', nriter, ': SSQ: ', ssq 
          do iart = 1, mart
            do iy = 1, maxyr
              do iseas = 1, mseas
                !print *, 'first', tacprev(iart,iy,iseas),tactemp(iart,iy,iseas)
                tactemp(iart,iy,iseas) = (tacprev(iart,iy,iseas) + tactemp(iart,iy,iseas)) / 2.0
                !print *, 'sec',tactemp(iart,iy,iseas)
              end do
            end do
          end do   
        else
          do iart = 1,mart
            do iy = 1,maxyr
              do iseas = 1,mseas
                tacprev(iart,iy,iseas) = tactemp(iart,iy,iseas)
              end do
            end do
          end do
        end if
      else
        ! If ssq is satisfactory, just sum up and leave
        do iart = 1, mart
          !if(iart.eq.3) write(30,*) (ssb(i,1),i=1,mart) 
          do iy = 1, maxyr
            dectac(iart,iy) = 0.0
            do iseas = 1, mseas
              dectac(iart,iy) = dectac(iart,iy) + tactemp(iart,iy,iseas)
            end do
          end do
        end do
        exit

      end if
    end do
  end subroutine decmodul

  subroutine dobsmodel(rn,mart,maxyr,mseas,maxl,speclen,rninit,options,rintvar)
    !! dobsmodel fills in rn with obs-model values for first year and season,
    !! by adding lognormal random error to the true numbers in rninit
    
    !! Output and input arguments to dobsmodel
    real, dimension(mart,0:maxyr,mseas,maxl), intent(out) :: rn
    integer, intent(in) :: mart, maxyr, mseas, maxl
    real, dimension(mart,1:4), intent(in) :: speclen
    real, dimension(mart,maxl,2), intent(in) :: rninit
    real, dimension(mart,22), intent(in) :: options
    double precision, intent(in) :: rintvar
        
    !! Local variables
    integer :: iart, iseas, lklstart, lklend, lkl
    real :: sigma, rx, xsi, rnx
  
    ! Calculate noisy initial stock numbers
    do iart = 1, mart
      iseas = int(options(iart,1))
      lklstart = int(speclen(iart,2))
      lklend = int(speclen(iart,3))
      do lkl = lklstart, lklend
        sigma = rninit(iart,lkl,2)
        rx = rnx(rintvar)
        xsi = snrn(rx) * sigma
        rn(iart,0,iseas,lkl) = rninit(iart,lkl,1) * exp(xsi) 
      end do
    end do
  end subroutine dobsmodel

  subroutine dproject(rn,tactemp,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2,bioparam3,options)
    !! dproject generates the stock number table rn by projecting the stock numbers at length forwards in time
    ! with given parameters and loss due to given quotas in tactemp. Adds recruits in the recruiting season.
   
    !! Local variables that are exchanged
    real, dimension(mart,0:maxyr,mseas), intent(inout) ::  tactemp
    real, dimension (mart,0:maxyr,mseas,maxl), intent(inout) :: rn

    ! General variables
    integer, intent(in) :: mart, maxyr, mseas, maxl
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,10), intent(in) :: bioparam1
    real, dimension(1:mart,0:maxyr,1:mseas,1:4,2), intent(in) :: bioparam2
    real, dimension(1:mart,0:maxyr,1:mseas,maxl,3), intent(in) :: bioparam3
    real, dimension(1:mart,22), intent(in) :: options

    ! Local variables
    integer :: depletflag, iart, iy, iseas, initseas, iseasstart, i
    real, dimension (mart,0:maxyr) :: ssb(mart,0:maxyr)
   
    do iart = 1, mart
      initseas = int(options(iart,1))
      do iy = 0, maxyr
        if (iy .eq. 0) then
          iseasstart = initseas
        else
          iseasstart = 1
        end if
      
        do iseas = iseasstart, mseas
          if (.not. (iy .eq. 0 .and. iseas .eq. initseas))  then
            call dgrow(iart,iy,iseas,rn,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2)
            !print *,'dproject dgrow'
          end if

          if (iseas .eq. (int(bioparam1(iart,iy,4)))) then
            call dgetssb(iart,iy,rn,ssb,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam3)
            !print *,'dproject dgetssb'
          end if

          if (iseas .eq. (int(bioparam1(iart,iy,5)))) then
            call drecruit(iart,iy,rn,ssb,mart,maxyr,mseas,maxl,speclen,bioparam1)
            !print *,'dproject drecruit'
          end if

          call ddie(iart,iy,iseas,rn,tactemp,mart,maxyr,mseas,maxl,speclen,bioparam2,bioparam3,depletflag)
          !print *,'dproject ddie'
          
          if (depletflag .eq. 1) then
            ! If not enought fish, set the TAC to 0 for the rest of the year
            do i = iseas, mseas
              tactemp(iart,iy,i) = 0.0
            end do ! iseas to mseas
          end if
          !print *,'dproject complete:',iart,'species',iy,'year',iseas,'season'
        end do
      end do
    end do
  end subroutine dproject

  subroutine drecruit(iart,iy,rn,ssb,mart,maxyr,mseas,maxl,speclen,bioparam1)
    ! In and out arguments
    integer, intent(in) :: iart, iy, mart, maxyr, mseas, maxl
    real, dimension(mart,0:maxyr,mseas,maxl), intent(inout) :: rn
    real, dimension(mart,0:maxyr), intent(in) :: ssb
    real, dimension(mart,1:4), intent(in) :: speclen
    real, dimension(mart,0:maxyr,10), intent(in) :: bioparam1
    
    ! Local
    real, dimension (maxl) :: rnew
    real :: a, b, cumnorm, fract, float, x, fractprev, rall, sigma
    integer :: ifunc, lkl, i, iseas, nlend, nlstart, lklend, lklstart, grlen
    
    ifunc = int(bioparam1(iart,iy,8))
    a = bioparam1(iart,iy,9)
    b = bioparam1(iart,iy,10)
    
    ! Find recruitment as function of SSB
    if (ifunc .eq. 1) then
      ! Hockey stick
      rall = min(a, a*ssb(iart,iy) / b)
    else if (ifunc .eq. 2) then
      ! Beverton-Holt
      rall = a * ssb(iart,iy) / (b + ssb(iart,iy))
    else if (ifunc .eq. 3) then
      ! Ricker
      rall = a * ssb(iart,iy) / b * exp(1.0 - ssb(iart,iy) / b)
    else
      write(*,*) 'Warning: Invalid recruitment function'
      rall = 0.0
    end if
    
    ! spread on length classes
    ! initialize rnew
    do lkl = 1, maxl
      rnew(lkl) = 0.0
    end do
    
    grlen = int(bioparam1(iart,iy,2))
    sigma = bioparam1(iart,iy,3)
    lklstart = int(speclen(iart,2))
    lklend = int(speclen(iart,3))
    nlstart = int(grlen - 3.0 * sigma)
    nlend = int(grlen + 3.0 * sigma) + 1
    iseas = int(bioparam1(iart,iy,5))
    fractprev = 0.0
    
    do i = nlstart, nlend
      x = (float(i) - grlen) / sigma
      fract = cumnorm(x)
      
      if (i .gt. lklend) then
        rnew(lklend) = rnew(lklend) + rall * (fract - fractprev)
      else if (i .lt. lklstart) then
        rnew(lklstart) = rnew(lklstart) + rall * (fract - fractprev)
      else
        rnew(i) = rall * (fract - fractprev)
      end if
      fractprev = fract
    end do
    
    ! Put in place in rn
    do lkl = lklstart, lklend
      rn(iart,iy,iseas,lkl) = rnew(lkl) + rn(iart,iy,iseas,lkl)
    end do
  end subroutine drecruit

  subroutine dgetssb(iart,iy,rn,ssb,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam3)
    !! In and out arguments
    integer, intent(in) :: iart, iy, mart, maxyr, mseas, maxl
    real, dimension(mart,0:maxyr,mseas,maxl), intent(in) :: rn
    real, dimension(mart,0:maxyr), intent(inout) :: ssb
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,10), intent(in) :: bioparam1
    real, dimension(1:mart,0:maxyr,1:mseas,maxl,3), intent(in) :: bioparam3
    
    !! Local variables
    integer :: lklstart, lklend, iseas, lkl
  
    !! Local variables - outputs of dgetweight and dgetmatur
    real, dimension(maxl) :: weight, rmat
  
    lklstart = int(speclen(iart,2))
    lklend = int(speclen(iart,3))
    iseas = int(bioparam1(iart,iy,4))
    call dgetweight(iart,iy,iseas,weight,mart,maxyr,mseas,maxl,speclen,bioparam3)
    call dgetmatur(iart,iy,rmat,mart,maxyr,maxl,speclen,bioparam1)
    ssb(iart,iy) = 0.0
  
    do lkl = lklstart, lklend
      ! print *, lkl,rn(iart,iy,iseas,lkl),weight(lkl),rmat(lkl)
      ssb(iart,iy) = ssb(iart,iy) + rn(iart,iy,iseas,lkl) * weight(lkl) * rmat(lkl)
    end do
  end subroutine dgetssb

  subroutine dgettsb(iart,iy,iseas,rn,tsb,mart,maxyr,mseas,maxl,speclen,bioparam3)
    !! TSB is calculated for each season as the biomass of fish above a certain length
    
    !! In and out arguments
    integer, intent(in) :: iart, iy, iseas, mart, maxyr, mseas, maxl
    real, dimension(mart,0:maxyr,mseas,maxl), intent(in) :: rn
    real, dimension(mart,0:maxyr,mseas), intent(out) :: tsb
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,1:mseas,maxl,3), intent(in) :: bioparam3

    !! Local variables
    integer :: lklstart, lklend, lkl
    real, dimension (maxl) :: weight
    
    lklstart = nint(speclen(iart,4))
    lklend = nint(speclen(iart,3))
    tsb(iart,iy,iseas) = 0.0
   
    call dgetweight(iart,iy,iseas,weight,mart,maxyr,mseas,maxl,speclen,bioparam3)
    
    do lkl = lklstart, lklend 
      tsb(iart,iy,iseas) = tsb(iart,iy,iseas) + rn(iart,iy,iseas,lkl) * weight(lkl)
    end do
  end subroutine dgettsb

  subroutine dgrow(iart,iy,iseas,rn,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2)
    !! dgrow calculates new numbers at length due to growth and redistribution on length classes at
    !! the start of the current time step. This is entirely different from the main program, where 
    !! the length of each super individual can be recorded on a continous scale.
    
    !! This procedure assumes lengths are evenly distributed within each length class, and
    !! all lengths are incrementded equally. Those passing length class limits are put in new length
    !! classes. This avoids bias and avoids stochastic growth in the short term projection (JTT: ??).
   
    !! In and out arguments
    integer, intent(in) :: iart, iy, iseas, mart, maxyr, mseas, maxl
    real, dimension(mart,0:maxyr,mseas,maxl), intent(inout) :: rn
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,10), intent(in) :: bioparam1
    real, dimension(1:mart,0:maxyr,1:mseas,1:4,2), intent(in) :: bioparam2
    
    !! Local variables
    real, dimension(maxl) :: rnew
    real :: deltat, dl, sigma, actlen, fract, deltal, k_value
    integer :: lklstart, lklend, lkl, ndelta
   
    deltat = 1.0 / float(mseas)
    dl = speclen(iart,1)
    lklstart = int(speclen(iart,2))
    lklend = int(speclen(iart,3))
    sigma = bioparam2(iart,iy,iseas,2,1)
    k_value = bioparam2(iart,iy,iseas,1,1) 
   
    do lkl = lklstart, lklend
      rnew(lkl) = 0.0
    end do
   
    do lkl = lklstart, lklend
      !! JTT: If herring, set separate k in Von Bertalanffy for juveniles (CHECK)
      if(iart .eq. 3) then
        k_value = bioparam2(iart,iy,iseas,1,1)
        if(lkl .lt. bioparam1(iart,iy,6)) k_value = bioparam2(iart,iy,iseas,1,2)
      end if
    
      actlen = (float(lkl) + 0.5) / dl
      
      !! Growth increment occurs only if Linf > L
      if (bioparam1(iart,iy,1) .gt. actlen) then
        deltal = (bioparam1(iart,iy,1) - actlen) * (1.0 - exp(-k_value * deltat))
      else
        deltal = 0.0
      end if
    
      ndelta = int(deltal)
      fract = deltal - float(ndelta)
    
      if ((lkl + ndelta) .ge. (lklend)) then
        rnew(lklend) = rnew(lklend) + rn(iart,iy,iseas,lkl) * (1.0 - fract)
      else
        rnew(lkl + ndelta) = rnew(lkl + ndelta) + rn(iart,iy,iseas,lkl) * (1.0 - fract)
      end if

      if((lkl + ndelta + 1) .ge. lklend) then
        rnew(lklend) = rnew(lklend) + rn(iart,iy,iseas,lkl) * fract
      else
        rnew(lkl + ndelta + 1) = rnew(lkl + ndelta + 1) + rn(iart,iy,iseas,lkl) * fract
      end if
      rn(iart,iy,iseas,lkl) = rnew(lkl)
    end do
  end subroutine dgrow

  subroutine ddie(iart,iy,iseas,rn,tactemp,mart,maxyr,mseas,maxl,speclen,bioparam2,bioparam3,depletflag)
    !! ddie calculates loss in each length class in the current time step, due to natural mortality and the catch
    !! in the time step according to the TAC.
   
    !! In and out arguments
    integer, intent(in) :: iart, iy, iseas, mart, maxyr, mseas, maxl
    integer, intent(out) :: depletflag
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,1:mseas,1:4,2), intent(in) :: bioparam2
    real, dimension(1:mart,0:maxyr,1:mseas,maxl,3), intent(in) :: bioparam3
    real, dimension (mart,0:maxyr,mseas,maxl), intent(inout) :: rn
    real, dimension (mart,0:maxyr,mseas), intent(inout) :: tactemp
    
    !! Local variables
    integer :: lklstart, lklend, lkl, iynext, iseasnext, ii
    real, dimension (maxl) :: sel, weight, rnact, rmact, ctemp
    real :: ctempsum, fact, sqdiff, fcheck
    
    lklstart = int(speclen(iart,2))
    lklend = int(speclen(iart,3))
    depletflag = 0
    do lkl = lklstart, lklend
      rnact(lkl) = rn(iart,iy,iseas,lkl)
      !! Natural mortality form bioparam3 is on annual scale, convert here to seasonal
      rmact(lkl) = bioparam3(iart,iy,iseas,lkl,1) / float(mseas)  
    end do
   
    call dgetselect(iart,iy,iseas,sel,mart,maxyr,mseas,maxl,speclen,bioparam2)
    call dgetweight(iart,iy,iseas,weight,mart,maxyr,mseas,maxl,speclen,bioparam3)

    !! Must check that there is enough fish to take the TAC (tactemp)
    !! JTT: Need to detail how this checker works and what the criterion are.
    !! For example, makes an assumption about max F at 2.0 (?)
    if (tactemp(iart,iy,iseas) .gt. 0.0) then
      ctempsum = 0.0
      fcheck = 2.0 / float(mseas)
      do lkl = lklstart, lklend
        ctempsum = ctempsum + rnact(lkl) * fcheck * sel(lkl) / (fcheck * sel(lkl) + rmact(lkl)) * (1.0 - exp(-rmact(lkl) - &
        fcheck * sel(lkl))) * weight(lkl)
      end do

      if (ctempsum .lt. tactemp(iart,iy,iseas)) then
        depletflag = 1
        do lkl = lklstart, lklend
          ! There is not enough fish to take the proposed TAC even with an annual F of 2. 
          ! Set the TAC according to fcheck for all lengths, and set the depletflag
          ctemp(lkl) = rnact(lkl) * fcheck / (fcheck * sel(lkl) + rmact(lkl)) * (1.0 - exp(-rmact(lkl) - fcheck * sel(lkl)))&
           * weight(lkl)
          rnact(lkl) = rnact(lkl) * exp(-rmact(lkl) - fcheck * sel(lkl))
        end do
      else 
        ! Enough fish, line search for the f, start somewhere = 0.02 seasonal
        fact = 0.02
        ctempsum = 0.0
        
        !! JTT: Following do loop replaces previous goto - check to see if this produces same result!
        do ii = 1, 20
          do lkl = lklstart, lklend
            ctempsum = ctempsum + rnact(lkl) * fact * sel(lkl) / (fact * sel(lkl) + rmact(lkl)) * (1.0 - exp(-rmact(lkl) - &
            fact * sel(lkl))) * weight(lkl)
          end do
      
          if (ctempsum .gt. 0.0) then
            sqdiff = (log(tactemp(iart,iy,iseas) / ctempsum)) ** 2
          else
            sqdiff = 0.0
          end if
      
          if (sqdiff .gt. 0.00001 .and. ctempsum .gt. 0.0) then
            ! Sikring mot cyclic iterations (JTT: This is old comment - what does it mean?)
            if (ii .lt. 20) then
              fact = fact * tactemp(iart,iy,iseas) / ctempsum
            else
              write(*,*) 'Warning: Line search in ddie within enac_decision_model module broken'
              fact = (fact + fact * tactemp(iart,iy,iseas) / ctempsum) / 2.0
              tactemp(iart,iy,iseas) = (tactemp(iart,iy,iseas) + ctempsum) / 2.0
              exit
            end if
          else
            exit
          end if ! if sqdiff>0.000001
        end do
      end if ! if enough fish
    else ! case tactemp=0
      fact = 0.0
    end if
    
    do lkl = lklstart, lklend
      rnact(lkl) = rnact(lkl) * exp(-fact * sel(lkl) - rmact(lkl))  
    end do
   
    if (iseas .eq. mseas) then
      iynext = iy + 1
      iseasnext = 1
    else
      iynext = iy
      iseasnext = iseas + 1
    end if
   
    if (iynext .le. maxyr) then
      do lkl = lklstart, lklend
        rn(iart,iynext,iseasnext,lkl) = rnact(lkl)
      end do
    end if
  end subroutine ddie

  subroutine dgetdecbasis(rn,decbasis,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam3,options)
    !! dgetdecbasis generates the basis for decisions (JTT: ?). The decbasis array is by species and has two dimensions, which are for the two trigger points.
    !! Values that refer to all species are equal for all the species.
    
    !! In and out arguments
    integer, intent(in) :: mart, maxyr, mseas, maxl
    real, dimension(mart,0:maxyr,mseas,maxl), intent(in) :: rn
    real, dimension(mart,2), intent(out) :: decbasis
    real, dimension(1:mart,1:4), intent(in) :: speclen
    real, dimension(1:mart,0:maxyr,10), intent(in) :: bioparam1
    real, dimension(1:mart,0:maxyr,1:mseas,maxl,3), intent(in) :: bioparam3
    real, dimension(1:mart,22), intent(in) :: options
    
    !! Local variables
    integer :: iref, iart, idesc, iydesc1, iydesc2, idescart, idesccode, iy, iseas, nart
    real :: summ, sumxy, sumy, sumx, sumxx, xmean
    real, dimension(mart,0:maxyr) :: ssb
    real, dimension(mart,0:maxyr,mseas) :: tsb
    real, dimension(mart,0:maxyr) :: tempdesc
   
    ! Loop over reference point 1 and 2
    do iref = 1, 2
      do iart = 1, mart
        ! Get the options for each reference point
        if (iref .eq. 1) then
          idesc = int(options(iart,8))
          iydesc1 = int(options(iart,4))
          iydesc2 = int(options(iart,5))
          idescart = int(options(iart,10))
          idesccode = int(options(iart,12))
        else
          idesc = int(options(iart,9))
          iydesc1 = int(options(iart,6))
          iydesc2 = int(options(iart,7))
          idescart = int(options(iart,11))
          idesccode = int(options(iart,13))
        end if
   
        ! Get the SSBs and TSBs
        do iy = iydesc1, iydesc2
          if (idesc .eq. 1) then
            call dgetssb(iart,iy,rn,ssb,mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam3)
            tempdesc(iart,iy) = ssb(iart,iy)
          else if (idesc .eq. 2) then
            summ = 0.0
            do iseas = 1, mseas
              call dgettsb(iart,iy,iseas,rn,tsb,mart,maxyr,mseas,maxl,speclen,bioparam3)
              summ = summ + tsb(iart,iy,iseas)
            end do
            tempdesc(iart,iy) = summ / float(mseas)
          end if
        end do
      end do
      
      ! Convert to decbasis
      ! If all species together, place the sum in tempdesc(1,.)
      if (idescart .eq. 2) then
        do iy = iydesc1, iydesc2
          do iart = 2, mart
            tempdesc(1,iy) = tempdesc(1,iy) + tempdesc(iart,iy)
          end do
        end do
        nart = 1
      else
        nart = mart
      end if
   
      ! If average over years  
      if (idesccode .eq. 1) then
        do iart = 1, nart
          summ = 0.0
          do iy = iydesc1, iydesc2
            summ = summ + tempdesc(iart,iy)
          end do   
          decbasis(iart,iref) = summ / float(iydesc2 - iydesc1 + 1)
        end do
      else
        ! trend (only if the two years are different)
        if (iydesc2 .le. iydesc1) then
          write(*,*) 'Warning: You ask for a trend over only one year'
          stop
        end if
        
        do iart = 1, mart
          sumxy = 0.0
          sumy = 0.0
          sumx = 0.0
          sumxx = 0.0 
       
          do iy = iydesc1, iydesc2
            sumxy = sumxy + float(iy) * tempdesc(iart,iy)
            sumy = sumy + tempdesc(iart,iy)
            sumx = sumx + float(iy)
            sumxx = sumxx + float(iy * iy)
          end do
       
          xmean = sumx / float(iydesc2 - iydesc1)
          decbasis(iart,iref) = (sumxy - xmean * sumy) / (sumxx - xmean * sumx)
        end do   
      end if ! end if trend
      
      ! When the decision basis is for all species, it is copied to them all
      if (mart .gt. 1 .and. idescart .eq. 2) then
        do iart = 2, mart
          decbasis(iart,iref) = decbasis(1,iref)
        end do
      end if
    end do  ! end loop over reference points
  end subroutine dgetdecbasis

end module enac_decision_model