module enac_inout

    use enac_commons
    use enac_functions

    implicit none

    integer :: io

contains

    subroutine setframe()
        !! Reads file 'framework' and sets the basic framework for the simulation
        character(len=100) :: string
        real :: values(maxspec, 8)
        integer :: i, j

        open(10, file="in/"//"framework", form="formatted")
        call read_line(10, string, io); read(string,*) nyear
        call read_line(10, string, io); read(string,*) ntimestep
        call read_line(10, string, io); read(string,*) niter
        call read_line(10, string, io); read(string,*) nspec
        do i = 1, nspec
            call read_line(10, string, io); read(string,*) species_name(i)
        end do
        do i = 1, nspec
            do j = 1, 6
                call read_line(10, string, io); read(string,*) values(i,j)
            end do
        end do
        close(10)

        do i = 1, nspec
            nage(i) = nint(values(i,1))
            nsiperyear(i) = nint(values(i,2))
            spawning_season(i) = nint(values(i,3))
            recruit_season(i) = nint(values(i,4))
            nlinkfunc(i) = nint(values(i,5))
            fprime(i) = values(i,6)
        end do
    end subroutine setframe

    subroutine read_line(unit, line, io_status)
        !! Reads a line from a file, skipping comments and blank lines
        !! Original code by Bjørn Ådlandsvik, IMR, 1997

        !! In and out arguments
        integer, intent(in) :: unit !! File unit 
        character(len=*), intent(out) :: line !! Read character string
        integer, intent(out) :: io_status

        !! Local variables
        character(len=*), parameter :: comchar = "*!#" ! Comment characters
        integer :: ipos

        do
            read(unit, "(a)", iostat=io_status) line
            if (io_status .ne. 0) exit
            ipos = scan(line, comchar)
            if (ipos .ne. 0) then
                line = line(:ipos-1)
            end if
            if (len_trim(line) .ne. 0) then
                exit
            end if
        end do
    end subroutine read_line

    subroutine popinput()
        !! Reads file 'bioinn.txt' to standard parameters array

        character(len=1) :: dummy
        integer :: i,j,k

        stdparameters = -999.0
        open (10, file='in/'//'bioinn.txt', form='formatted')
        ! Skip header line
        read(10, '(a1)') dummy
        do
            read(10, *, iostat=io) i, j, k, stdparameters(i,j,k,1), stdparameters(i,j,k,2)
            if (io .ne. 0) exit
        end do
        close(10)
    end subroutine popinput

    subroutine recinput()
        !! Reads file 'recruitment' and sets parameters for recruitment model

        character(len=100) :: string
        real :: tmp_value
        integer :: i, j, k

        recruit_parameters = -999.0

        open(10, file="in/"//"recruitment", form="formatted")
        do k = 1, nspec
            call read_line(10, string, io); read(string,*) species, nregimes(species)
            call read_line(10, string, io); read(string,*) species, icodeshift(species)
            if (icodeshift(species) .eq. 1) then
                listreg(species,1) = 1
                call read_line(10, string, io)
                read(string,*) species, (listreg(species,i), i = 2, nregimes(species))
                listreg(species, nregimes(species)+1) = nyear
            else
                call read_line(10, string, io)
                read(string,*) species, regintv(species), lastreg(species)
            end if
        end do
        do
            call read_line(10, string, io)
            if (io .ne. 0) exit
            read(string,*) species, i, j, tmp_value
            if (j .gt. 0) then
                recruit_parameters(species,i,j,1) = tmp_value
                ! NOTE(EAM): Why is it necessary with 3 copies of this value?
                recruit_parameters(species,i,j,2) = recruit_parameters(species,i,j,1)
                recruit_parameters(species,i,j,3) = recruit_parameters(species,i,j,1)
            else
                prime_factor(species) = tmp_value
            end if
        end do
        close(10)

        ! Check recruitment parameters, there shall be 20 in total
        do species = 1, nspec
            do i = 1, nregimes(species)
                do j = 1, 20
                    if (recruit_parameters(species,i,j,1) .le. -999.0) then
                        print *, "Recruitment specification missing for"
                        print *, "species", species
                        print *, "Regime", i
                        print *, "Parameter number", j
                    end if
                end do
            end do
        end do
    end subroutine recinput

    subroutine linkinput()
        !! Reads 'linkfile' and sets link parameters
        character(len=100) :: string
        integer :: i, j, k

        link_parameters = -999.0
        open(10, file="in/"//"linkfile")
        call read_line(10,string, io); read(string,*) i,j,k,link_parameters(i,j,k)
        close(10)
    end subroutine linkinput

    subroutine agebasedinput()
        !! TODO(EAM): Insert description here!
        integer :: i, j, k, m
        real :: length, tmp_value
        character(len=100) :: string

        mat_pattern = -999.0
        mor_pattern = -999.0
        sel_pattern = -999.0

        open(10, file="in/"//"biobyage.txt")
        do
            call read_line(10, string, io); read(string,*) i, j, k, bio_opt(i,j,k)
            if (i .eq. nspec .and. j .eq. 3 .and. k .eq. 3) exit
        end do
        do
            call read_line(10, string, io)
            if (io .ne. 0) exit
            read(string,*) i,j,k,m,length,tmp_value
            if (j .eq. 1) then
                mat_pattern(i,k,m,1) = length
                mat_pattern(i,k,m,2) = tmp_value
            else if (j .eq. 2) then
                mor_pattern(i,k,m,1) = length
                mor_pattern(i,k,m,2) = tmp_value
            else if (j .eq. 3) then
                sel_pattern(i,k,m,1) = length
                sel_pattern(i,k,m,2) = tmp_value
            end if
        end do
        close(10)
    end subroutine agebasedinput

    subroutine envinput()
        !! Reads in 'envfile' and sets basis for environmental scenarios

        integer :: i, j
        character(len=100) :: string

        varenv1 = -999.0
        open(10, file="in/"//"envfile", form="formatted")
        do 
            call read_line(10, string, io)
            if (io .ne. 0) exit
            read(string,*) i, j, varenv1(i,j,1)
        end do
        close(10)
        
        !! Index 1 of 3rd dimension is true value. Index 2 is value perceived by management 
        !! Below assumes management has perfect knowledge of environmental factor (in envfile)
        !! If environmental factor is incorporated, this varenv1(,,2) should have error
        !! Currently, varenv1 is only used within r_link_test procedure (enac_links) to specify recruitment regimes
        varenv1(i,j,2) = varenv1(i,j,1)
    end subroutine envinput

    subroutine ruleinput()
        !! Reads in 'optinn' which contain input to management rules
        
        integer :: i, j
        character(len=100) :: string
        
        manage_opt = -999.0
        open (10, file='in/'//'optinn', form='formatted')
        do
            call read_line(10, string, io)
            if (io .ne. 0) exit
            read(string,*) i, j, manage_opt(i,j)
        end do
        close(10)
        
        !! Initial value should be 1; if not, set to 1 here
        !! (JTT: Rewrote original comment from Norwegian, but do not know what this means)
        do i = 1, nspec
            if (manage_opt(i,1) .ne. 1) then
                write(*,*) 'Management option 1 should have value 1'
                write(*,*) 'Setting manage_opt(', i, ',1) = 1 ...'
                manage_opt(i,1) = 1
            end if
        end do
    end subroutine ruleinput

    subroutine statein()
        !! TODO: Write description
        integer :: o, p, i 
 
        open (10, file='in/'//'lenvekt.txt')

        do o = 1, 3
            read(10,*) (mlengths(o,i), i=1, 14)
        end do
        
        do p = 1, 3
            read(10,*) (mweights(p,i), i=1, 14)
        end do
    end subroutine statein

    subroutine output()
        integer :: i, j, k

        do j = 1, nyear
            write(11,*) iter, j, (tac(i,j),i=1,nspec)
            write(13,*) iter, j, (ssb(i,j),i=1,nspec)
        end do

        do i = 1, nspec
            do j = 1, nyear
                do k = 1, (maxage+1)
                    write(14,*) iter, j, i, k-1, natage(i,j,k), catage(i,j,k)
                end do
            end do
        end do
    end subroutine output

    subroutine sistateout()
        integer :: species, i, j, k, nsitemp, nsitem
        real :: a(maxspec,8,nyear + 20 + 1), su2, su1(100) 
       
        1 format(181(f10.2,x))
        do species = 1, maxspec
            nsitemp = int(nsi2(species) * 0.01)
            do i = 1, nsitemp 
                do j = 1, 8
                    do k = 1, 100
                        su1(k) = sistate(species,((i - 1) * 100) + k,j)
                        if (j == 4) su1(k) = su1(k) * 1000 
                    end do
                    
                    su2 = sum(su1)
                    
                    if (su2 .gt. 0) then
                        a(species,j,i) = su2 / nsiperyear(species) 
                    else
                        a(species,j,i) = 0.
                    end if
                end do  
            end do
        
            do j = 1, 8
                if(species .eq. 1) write(15,1) (a(species,j,i), i = 1, nsitemp)
                if(species .eq. 2) write(16,1) (a(species,j,i), i = 1, nsitemp)
                if(species .eq. 3) write(17,1) (a(species,j,i), i = 1, nsitemp)
            end do
        
            if(species .eq. 1) write(15, *) 
            if(species .eq. 2) write(16, *)
            if(species .eq. 3) write(17, *)
       end do
    end subroutine sistateout

    subroutine mort_mat_out()
        integer	:: i, j
        real :: length, mlow, mhigh, l1_50, l2_50, slope1, slope2, m
        real :: ma(5:45)
        
        do i = 1, maxspec
            do j = 5, 45
                length = j * 1.
                mlow = siparam(i,1,5,1,1)
                mhigh = siparam(i,1,5,2,1)
                l1_50 = siparam(i,1,5,3,1)
                l2_50 = siparam(i,1,5,4,1)
                slope1 = siparam(i,1,5,5,1)
                slope2 = siparam(i,1,5,6,1)
                call comblog(mlow,mhigh,l1_50,slope1,l2_50,slope2,length,m)
                ma(j) = m
            end do
            write(23,'(2I4,181F10.2)') iter, i, ma
        end do   
    end subroutine mort_mat_out

    subroutine outstate()
        integer :: i, age
        real :: mweight0, mlength0
        real, dimension(14) :: outweight, outlength
       
        if (year .gt. 0) then
            do species = 1, maxspec
                do i = (nsi2(species) / 100 - 14), (nsi2(species) / 100 - 1)
                    age = int(sistate(species,1 + (i * 100), 2))
                    mweight0 = sum(sistate(species,(1 + (i * 100)):(99 + (i * 100)),4)) / 100
                    mlength0 = sum(sistate(species,(1+(i*100)):(99+(i*100)),3))/100
                    outweight(age) = mweight0
                    outlength(age) = mlength0  
                end do

                if (species .eq. 1) then
                    write(41, '(2I4,14F6.3)') iter, year, outweight
                    write(44, '(2I4,14F6.2)') iter, year, outlength
                end if

                if (species .eq. 2) then
                    write(42, '(2I4,14F6.3)') iter, year, outweight
                    write(45, '(2I4,14F6.2)') iter, year, outlength
                end if

                if(species .eq. 3) then
                    write(43, '(2I4,14F6.3)') iter, year, outweight
                    write(46, '(2I4,14F6.2)') iter, year, outlength
                end if
            end do
        end if
    end subroutine outstate

    subroutine setup_recruits()
        !! setup_recruits specifically sets up regimes and spasmodic recruitment.
        !! This includes randomization routines to determine when these occur, and to what extent.
        !! Outputs nregact(species,year), the regime index in year, and ispasm(species,year), years with spasmodic recruitment.

        !! TODO: Write intents!
        real, dimension(0:maxyear) :: geodist
        real :: x, rinterval, rintvcv, ri, flii, pr, px, rilast, spasmfct, xsi
        integer :: ishift, i, j, k, ii, iy, nn, ifirst, ilast, lastspasm, ioptspasm, p
         
        do species = 1, nspec
            ! Set up regime shifts      
            if (nregimes(species) .le. 1) then
                ! case no regimeshifts
                do year = -maxage, nyear
                    nregact(species,year) = 1
                end do
            else
                ! nregimes >1
                ! Set up sequence of regimes for recruitment
                ! If icodeshift=1, listreg is ready from input
                ! If icodeshift=2, listreg is set up here      
                if (icodeshift(species) .eq. 2 .or. icodeshift(species) .eq. 3) then   
                    call creategeo(geodist,regintv(species),nyear)
                    listreg(species,1) = 1
                    ishift = lastreg(species)
                    nn = 1

                    ! Find next shift
                    do
                        x = rnx(rintvar)
                        ii = 1
                        do
                            ii = ii + 1
                            if (geodist(ii) .ge. x) exit
                            
                            if (ii .lt. nyear) then 
                                cycle
                            else
                                exit
                            end if
                        end do
                        
                        ishift = ishift + ii

                        if (ishift .gt. 0 .and. ishift .le. nyear .and. nn .le. nyear) then
                            nn = nn + 1
                            listreg(species,nn) = ishift
                        end if

                        if (ishift .lt. nyear) then 
                            cycle
                        else
                            exit
                        end if
                    end do

                    listreg(species,nn+1) = nyear 
                else 
                    ! icodeshift not 1        
                    if (icodeshift(species) .ne. 1) then
                        write(*,*) 'Error: Incomprehensible code for regime shift', species, icodeshift(species)
                        stop
                    end if 
                end if
            end if
        
            ! Known shift years are now in listreg (fixed or derived from flex model)
            ! Set up nregact: Regime by year from listreg
            ! Always start with regime 1, so prime nregact with that
            do i = -maxage, nyear
                nregact(species,i) = 1 
            end do  
          
            if (nregimes(species) .gt. 1) then
                j = 1
                ! for each regime block, define the years
                ! The last entry on listreg shall be nyear
                do
                    j = j + 1
                    ifirst = listreg(species,j)
                    ilast = listreg(species,j+1) - 1
                    
                    if (listreg(species,j+1) .eq. nyear) ilast = nyear
            
                    ! Then, find the regime number k for those years
                    ! icodeshift = 2 is for random regimes
                    if (icodeshift(species) .eq. 2) then
                        ! Find which regime k in this interval, each is equally likely
                        x = rnx(rintvar)
                        x = x * float(nregimes(species))
                        k = int(x) + 1

                    ! icodeshift = 1 or 3 
                    else
                        k = j
                    end if
            
                    do iy = ifirst, ilast
                        nregact(species,iy) = k
                    end do
                    
                    !  end loop over regimes
                    if (ilast .lt. nyear) then
                        cycle
                    else 
                        exit
                    end if
                end do
            end if
        
            ! Set up the sequence of spasmodic recruitment successes.  
            do i = -nage(species), nyear
               ispasm(species,i) = 0
            end do
        
            if(icodeshift(species) .lt. 3) then
                ! set up spasmodic years for each regime
                do k = 1, nregimes(species)
                    rinterval = recruit_parameters(species,k,13,1)
                    rintvcv = recruit_parameters(species,k,14,1)
                    spasmfct = recruit_parameters(species,k,15,1)
                    
                    if (k .eq. 1) then
                        ifirst = 0
                    else
                        ifirst = listreg(species,k)
                    end if
                
                    if (k.eq.nregimes(species)) then
                        ilast=nyear
                    else
                        ilast=listreg(species,k+1)
                    end if
                
                    lastspasm = ifirst + int(recruit_parameters(species,k,16,1))
                
                    if (rinterval .gt. 1.0 .and. spasmfct .gt. 1.0) then
                        ! only if the interval and spasmodic factor is greater than 1, 
                        ! if not, no such effects are assumed
                        ioptspasm = int(recruit_parameters(species,k,12,1))		
            
                        if (ioptspasm .eq. 1) then
                            ! case spasms with random intervals
                            rilast = float(lastspasm)

                            do
                                x = rnx(rintvar)
                                xsi = snrn(x)
                                ri = rinterval * exp(rintvcv * xsi) + rilast
                                j = int(ri + 0.5)
                                
                                if (j .le. ilast .and. j .ge. ifirst) then
                                    ! Note: The spasmodic years start in lastspasm
                                    ! Years prior to year 0 only matter when priming the stock.
                                    ! Else, earlier years are only used to propagate the sequence. 
                                    ispasm(species,j) = 1
                                    rilast = ri
                                    cycle
                                else
                                    exit
                                end if
                            end do     
                        else if (ioptspasm .eq. 2) then
                            pr = 1.0 / rinterval
                            ilast = lastspasm

                            do
                                x = rnx(rintvar)
                                ii = 0
                                px = 0.0

                                do
                                    ii = ii + 1
                                    flii = float(ii) - 1.0
                                    px = px + pr * exp(flii * log(1.0 - pr))
                                    
                                    if (px .ge. x) then
                                        exit
                                    else 
                                        if (ii .lt. nyear) then
                                            cycle
                                        else
                                            exit
                                        end if
                                    end if
                                end do

                                j = ii + ilast
                                ilast = j   
                                
                                if (j .le. ilast .and. j .ge. ifirst) then
                                    ! Note: The spasmodic years start in lastspasm
                                    ! Years prior to year 0 only matter when priming the stock.
                                    ! Else, earlier years are only used to propagate the sequence. 
                                    ispasm(species,j) = 1
                                    cycle
                                else
                                    exit
                                end if
                            end do
                        else if (ioptspasm .eq. 3) then
                            ! Flat distribution between rinterval*(1-cv,1+cv)
                            ilast = lastspasm

                            do
                                x = rnx(rintvar)
                                x = rinterval + ((x - 0.5) * rintvcv * 2)
                                ii = int(x + 0.5) + ilast

                                if (ii .lt. nyear) then
                                    if (ii .ge. -nage(species)) ispasm(species,ii) = 1
                                    ilast = ii
                                    cycle
                                else
                                    exit
                                end if
                            end do
                        end if
                    end if  ! end if spasmodic recruitments
                end do  ! end regimes loop
            else  ! icodeshift = 3 
                do p = 1, nregimes(species) 
                    nregact(species,iy) = p
                    rinterval = recruit_parameters(species,p,13,1)
                    rintvcv = recruit_parameters(species,p,14,1)
                    spasmfct = recruit_parameters(species,p,15,1)
                    ioptspasm = int(recruit_parameters(species,p,12,1))
                    ilast = recruit_parameters(species,p,16,1)
                    
                    do
                        x = rnx(rintvar)
                        x = rinterval + ((x - 0.5) * rintvcv * 2)
                        ii = int(x + 0.5) + ilast
                        if (ii .lt. nyear) then
                            if (ii .ge. -nage(species)) spasm(species,p,ii) = 1
                            ilast = ii
                            cycle
                        else
                            exit
                        end if
                    end do
                end do
            end if  ! end if icodeshift = 1, 2 or 3
        end do  ! end species loop
    end subroutine setup_recruits
        
    !subroutine fixrec
    ! integer	:: o,i,test
    
    ! test=10001
    
    ! open(test, file = 'In/'//'rec.txt')
    
    ! 	do o = 1,45
    !	 read(test,*)(rec_in(o,i),i=1,3)
    !	enddo
    
    !end subroutine fixrec

end module enac_inout