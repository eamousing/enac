module enac_inout

    use enac_commons

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
        !!
        !! Original code by Bjørn Ådlandsvik, IMR, 1997
        integer, intent(in) :: unit !! File unit 
        character(len=*), intent(out) :: line !! Read character string
        integer, intent(out) :: io_status

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
            read(10, *, iostat=io) i,j,k,stdparameters(i,j,k,1), stdparameters(i,j,k,2)
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
            call read_line(10, string, io); read(string,*) i,j,k,bio_opt(i,j,k)
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
            read(string,*) i,j,varenv1(i,j,1)
        end do
        close(10)

        ! NOTE(EAM): I have no idea what the comment below means?
        !!!!!!!!!!!!!
        ! For the time being - managers know the real environment - should be corrupted
        varenv1(i,j,2) = varenv1(i,j,1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end subroutine envinput

end module enac_inout