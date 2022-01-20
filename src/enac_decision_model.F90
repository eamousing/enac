module enac_decision_model
    ! This is the original decmodul block, which was written as a stand-alone
    ! program that should be called as a .dll file. This block should later on be
    ! incorporated in the general enac program. 
    !
    ! The variables in decmodul are stand-alone, without the use of commons,
    ! All subroutines have names starting with d, to keep them apart form similar
    ! subroutines in the main part.
    !
    ! The subroutine decision_model sets up what was the input to the decomdul.
    ! It extracts the necessary information from common-variables, and adds some information
    ! that at present is hard-coded.
    ! Then it calls the decmodul itself.
    ! The end product from the decomodul is a TAC in dectac(species,year).

    implicit none

contains

subroutine decision_model(dectac,ntacyr)
    ! inputs: maxspec,maxyr,maxtimestep,maxlkl,
  
    ! transferred variables
    ! These were transferred to the original decmodul. Here, they are
    ! derived from the common variables, and then tansferred to the present decmodul
    ! Data shall be available for the decision process for all species, years, seasons
    ! not only the species, year and season for which the decision is to be made
    
    ! Note that the years (iy) here refer to the internal years in decompodul. The calling
    ! year is year 0 and the tac year is year 1. Normally, maxyr is 2, since it sometimes
    ! is needed to project the stock into the year after the TAC is taken
  
    integer,parameter :: maxyr=2
    ! Note: the parameter maxyr sets the outer dimension of the years. 
    integer :: mart,mseas,maxl,ntacyr
    !,j
    real, dimension (1:maxspec,1:4) :: speclen  !cols: interval, min size class, max size class, min size for TSB
    real, dimension (1:maxspec,0:maxyr,10) :: bioparam1  !3rd dim: for each parameter defined below 
    real, dimension (1:maxspec,0:maxyr,1:maxtimestep,1:4,2) :: bioparam2 !4th dim: parameter, 5th dim: parameter specific to something (juvenile and adult K)
    real, dimension (1:maxspec,0:maxyr,1:maxtimestep,maxlkl,3) :: bioparam3
    real, dimension (1:maxspec,1:maxlkl,2) :: rninit ! rninit(,,1) is initial stock numbers, rinit(,,2) is sigma for obs noise on these numbers
    real, dimension (1:maxspec,1:maxtimestep) :: seastac
    real, dimension (1:maxspec,22) :: options
    real, dimension (1:maxspec,0:maxyr) :: dectac
    integer, dimension (1:maxspec) :: regim
  
    ! Local variables
    integer :: i,iart,iy,iseas,lkl,lklminus,lklplus,si,regime
    real :: length,mlow, mhigh, slope1, slope2, l1_50, l2_50,natmor
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
    speclen(1,1)=1.0      
    speclen(2,1)=1.0      
    speclen(3,1)=1.0
    ! Smallest length class (minus-class)
    speclen(1,2)=13.0
    speclen(2,2)=21.0
    speclen(3,2)=10.0
    ! Largest length class (plus-class)
    speclen(1,3)=37.0
    speclen(2,3)=41.0
    speclen(3,3)=38.0
    ! Lowest length class to be included in TSB
    speclen(1,4)=15.0
    speclen(2,4)=15.0
    speclen(3,4)=15.0
  
    ! bioparam1
    do iart = 1, mart
        do iy = 0, maxyr
            regime = nregact(iart,iy)
            ! Linf
            bioparam1(iart,iy,1)=stdparameters(iart,1,1,1)
            ! Recruiting length
            bioparam1(iart,iy,2) = recruit_parameters(iart,regime,17,3)
            ! CV recruiting length
            bioparam1(iart,iy,3) = recruit_parameters(iart,regime,18,3)
            ! Spawning season
            bioparam1(iart,iy,4) = spawning_season(iart)
            ! Recruiting season
            bioparam1(iart,iy,5) = recruit_season(iart)
            ! L50 in maturation logistic function
            bioparam1(iart,iy,6) = decparam(iart,3,1)
            ! Slope in maturation logistic function
            bioparam1(iart,iy,7) = decparam(iart,3,2)
            ! Code for recruitment function
            bioparam1(iart,iy,8) = recruit_parameters(iart,regime,1,1)
            ! Note on recruit-parameters: Decmodul will only see the assumed a and b-parameters,
            ! not the spasmodic or periodic variants that go into the real recruitment
            ! But it will know the regime in the calling year    
            ! Recruitment a-parameter
            bioparam1(iart,iy,9) = recruit_parameters(iart,regime,4,3)
            ! Recruitment b-parameter
            bioparam1(iart,iy,10) = recruit_parameters(iart,regime,5,3)
        enddo
    enddo
  
    ! bioparam2
    do iart=1,mart
     do iy=0,maxyr
      do iseas=1,mseas
       ! k-values
         bioparam2(iart,iy,iseas,1,1) = decparam(iart,1,2) !Adult K
         !if(iart.eq.3) print *, bioparam2(iart,iy,iseas,1,1)
         bioparam2(iart,iy,iseas,1,2) = decparam(iart,1,4) !Juvenile K
         ! CV growth increment (arbitrary value - not used in the main program
       bioparam2(iart,iy,iseas,2,1)= decparam(iart,1,3)
       ! L50 in selection
       bioparam2(iart,iy,iseas,3,1)=decparam(iart,4,1)
       ! slope in selection
       bioparam2(iart,iy,iseas,4,1)=decparam(iart,4,2)
      enddo
     enddo
    enddo
  
    ! bioparam3
    do iart=1,mart
     do iy=0,maxyr
      do iseas=1,mseas      
       do lkl=1,maxl
        ! Natural mortality per year
        ! One cm length classes tacitly assumed
        length=float(lkl)+0.5
        mlow=stdparameters(iart,5,1,1)
        mhigh=stdparameters(iart,5,2,1)
        l1_50=stdparameters(iart,5,3,1)
        l2_50=stdparameters(iart,5,4,1)
        slope1=stdparameters(iart,5,5,1)
        slope2=stdparameters(iart,5,6,1)
        call comblog(mlow, mhigh, l1_50,slope1,l2_50,slope2,length,natmor)
        ! Natmor is scaled as annual mortalities
        ! Decmodul assumes natural mortality evenly distributed on seasons
        bioparam3(iart,iy,iseas,lkl,1)=natmor
        ! Condition according to environment 1, a      
        bioparam3(iart,iy,iseas,lkl,2)=decparam(iart,2,1)
        ! Condition according to environment 2, b      
        bioparam3(iart,iy,iseas,lkl,3)=decparam(iart,2,2)
       enddo
      enddo
     enddo
    enddo
  
  ! Initial numbers from sistate  
    do iart=1,mart
     ! Initialize      
     do lkl=1,maxl
      rninit(iart,lkl,1)=0.0
      ! rninit(.,.,2) is used to make noise to the N-values in the observation model.
      ! Here there is a need for extensions, e.g. length dependence, annual effects, etc.
      ! and it should not be hardcoded.
      rninit(iart,lkl,2)=0.30
     enddo
  
     lklminus=nint(speclen(iart,2))
     lklplus=nint(speclen(iart,3))
  
     ! Loop through each SI 
     do si=nsistart(iart),nsi2(iart)   
      !lkl=nint(sistatekeep(iart,si,3))
        lkl=int(sistatekeep(iart,si,3))
  
       if (lkl.lt.lklminus) then
        rninit(iart,lklminus,1)=rninit(iart,lklminus,1)+sistatekeep(iart,si,1)
       else if (lkl.gt.lklplus) then
        rninit(iart,lklplus,1)=rninit(iart,lklplus,1)+sistatekeep(iart,si,1)
       else
        rninit(iart,lkl,1)=rninit(iart,lkl,1)+sistatekeep(iart,si,1)
       endif
     enddo
  
    enddo
  
  
  ! seastac - hard-coded (seasonal distribution of TAC)
    do iart=1,mart ! Loop through each species
     do iseas=1,mseas ! Loop through each season
      seastac(iart,iseas)=1.0/float(mseas)
     enddo
    enddo
  
    do iart=1,mart
     do i=1,22
      options(iart,i)=manage_opt(iart,i) ! col 3 from optinn
     enddo
     ! Options 1 and 2 are exceptions:
     ! Option (.,1) is set to 1 in manage_opt in the main program, and it is transferred here.   
     ! Option (.,2) tac in the call year is set here to what is already determined as tac in year 0
     options(iart,2) = tac(iart,year)
    enddo  
  
    regim = recregime !KRU 10 sep 2013 transfer climate regime (recruitment regime) to decmodul 
    
    call decmodul(mart,maxyr,mseas,maxl,speclen,bioparam1,bioparam2,bioparam3,rninit,options,seastac,rintvar,dectac,regim)
    ntacyr = maxyr
  
   end subroutine decision_model

end module enac_decision_model