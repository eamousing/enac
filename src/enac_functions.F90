module enac_functions

    use, intrinsic :: iso_fortran_env, dp => real64
    use enac_commons

    implicit none

contains

    subroutine vonbert2(rlinf, rk, rl, deltat, epsi, deltal)
        !! TODO: Description of procedure and arguments using !!
        real, intent(in) :: rlinf
        real, intent(in) :: rk
        real, intent(in) :: rl
        real, intent(in) :: deltat
        real, intent(in) :: epsi
        real, intent(out) :: deltal

        real :: eps, eps1

        eps1 = exp(epsi)**0.5
        eps = rand_normal(0.0, eps1)

        deltal = (rlinf - rl) * (1.0 - exp(-rk * deltat)) + eps 
        if (deltal .lt. 0.0) deltal=0.0
        ! Maybe wrong to include this, but program crash without it
        ! TODO: Figure out why negative numbers causes a crash?
    end subroutine vonbert2

    subroutine logistic(l50, slope, x, v)
        !! value is a logistic function of x with value in (0,1)
        !! TODO: description of arguments using !!
        real, intent(in) :: l50
        real, intent(in) :: slope
        real, intent(in) :: x
        real, intent(out) :: v
              
        v = 1.0 / (1.0 + exp(-4.0 * slope * (x - l50)))
    end subroutine logistic

    subroutine scaledlog(plow, phigh, l50, slope, x, v)
        !! value is a logistic function of x with value in (plow,phigh)
        !! TODO: description of arguments using !!
        real, intent(in) :: plow
        real, intent(in) :: phigh
        real, intent(in) :: l50
        real, intent(in) :: slope
        real, intent(in) :: x
        real, intent(out) :: v
      
        call logistic(l50, slope, x, v)
        v = plow + (phigh - plow) * v
    end subroutine scaledlog

    subroutine comblog(plow, phigh, l51, slope1, l52, slope2, x, v)
        !! Returns the product of a logistic function and an inverse logistic function, 
        !! scaled to the interval (plow,phigh)
        !! TODO: description of arguments using !!
        real, intent(in) :: plow
        real, intent(in) :: phigh
        real, intent(in) :: l51
        real, intent(in) :: l52
        real, intent(in) :: slope1
        real, intent(in) :: slope2
        real, intent(in) :: x
        real, intent(out) :: v
        
        real :: v1,v2
      
        call logistic(l51, slope1, x, v1)
        call logistic(l52, slope2, x, v2)
        v = plow + (phigh - plow) * v1 * (1.0 - v2)
    end subroutine comblog

    subroutine assymbell(a1, a2, p, x, v)
        !! Returns an assymmetric bell function of x, with maximum = 1 at x=0, 
        !! and value in (0,1)
        !!
        !! parameters a and power: 
        !! a = a1 when x<0 and a2 when x>0.
        !! TODO: description of arguments using !!
        real, intent(in) :: a1
        real, intent(in) :: a2
        real, intent(in) :: p
        real, intent(in) :: x
        real, intent(out) :: v
      
        if (x .lt. 0.0) then
         v = exp(-(-x / a1)**p)
        else
         v = exp(-(x / a2)**p)
        endif
    end subroutine assymbell

    subroutine linear(a,b,x,v)
        !! TODO: procedure and arguments description using !!
        real, intent(in) :: a
        real, intent(in) :: b
        real, intent(in) :: x
        real, intent(out) :: v
        
        v=a+b*x
    end subroutine linear

    subroutine linint(x,x1,y1,x2,y2,y)
        !! Returns linear interpolation of x between (x1,y1) and (x2,y2)
        !! TODO: description of arguments using !!
        real, intent(in) :: x
        real, intent(in) :: x1
        real, intent(in) :: y1
        real, intent(in) :: x2
        real, intent(in) :: y2
        real, intent(out) :: y
              
        y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    end subroutine linint

    function cumnorm(x) result(res)
        !! Returns the cumulated probability of x in a standard normal distribution
        !! Based on an approximation by Bagby 1995, cited in 
        !! http://mathworld.wolfram.com/NormalDistributionFunction.html, formula 14
        !! That formula gives the integral from 0 to x,
        !! rearranged here to give the integral from -inf to x.
        !! TODO: description of arguments using !!
        real, intent(in) :: x
        real :: res
        real :: xx

        xx = -x * x
        res = 7.0 * exp(0.5 * xx) + 16.0 * exp(0.58558 * xx) + (7.0 - 0.7854 * xx) *exp(xx)
        res = 0.5 * sqrt(1.0 - res / 30.0)
        
        if (x .lt. 0.0) then
            res =0.5 - res
        else
            res = res + 0.5
        endif
    end function cumnorm

    function snrn(r) result(res)
        !! When r is a random number between 0 and 1
        !! the function returns the number with that probability
        !! in a standard normal distribution
        !! TODO: description of arguments using !!
        real, intent(in) ::  r
        real :: res
        
        real :: t, sg
        if (r .gt. 0.5) then
            t = sqrt(log((1.0 / (r - 0.5))**2))
            sg = 1.0
        else if (r .lt. 0.5) then
            t = sqrt(log((1.0 / (r - 0.5))**2))
            sg = -1.0
        else
            sg = 0.0
        endif
        res = sg * (t - (2.30753 + 0.27061 * t) / (1.0 + 0.99229 * t + 0.04481 * t * t))
    end function snrn

    real function rnx(rintvar)
        !! returns a random variable between 0 and 1, adjusts the seed rintvar
        !! (EAM): Not sure if this is supposed to overwrite rintvar in enax_commons
        !! and I have therefore not updated the structure. 
        !! TODO: Check the intentions of this function before refactoring.
        real(dp) :: rintvar,drn
        rintvar=rintvar*23.0_8
        rintvar=dmod(rintvar,100000001.0_8)
        drn=rintvar/100000001.0_8
        rnx =sngl(drn)
    end function rnx

    subroutine creategeo(geodist, rintv, maxyr)
        !! TODO: Not sure about the in/out intends. Should look up in procedure calls
        !! and add intentions.
        !! TODO: Description of procedure and arguments using !!
        real :: geodist(0:maxyr)
        real :: rintv
        real :: px
        real :: pr
        real :: flii
        integer :: maxyr
        integer :: i
        px = 0.0
        pr = 1.0 / rintv
        do i = 1, maxyr
            flii = float(i) - 1.0
            px = px + pr * exp(flii * log(1.0 - pr))
            geodist(i) = px
        enddo
    end subroutine creategeo

    function rand_normal(mean, stdev) result(c)
        !! This code was broken in the original code. Added theta and r so it compiles. 
        !! CAUTION: The function is called from vonbert procedure, so we should check 
        !! if the return value is as intended.
        !! TODO: Description of procedure and arguments using !!
        real, intent(in) :: mean
        real, intent(in) :: stdev
        real :: c
        
        real :: temp(2)
        real :: pi=3.141592653589793238462
        real :: r, theta
    
        if (stdev .le. 0.) then
            write(*,*) "Standard Deviation must be +ve"
        else
            call random_number(temp)
            r = (-2.0 * log(temp(1)))**0.5
            theta = 2.0 * pi * temp(2)
            c = mean + stdev * r * sin(theta)
        end if
    end function rand_normal

end module enac_functions