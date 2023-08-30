!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! an implementation of the stochastic SOM training algorithm based on
! ideas from tensor algebra
! written by Lutz Hamel, (c) University of Rhode Island
!
! LICENSE: This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free Software
! Foundation.
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! A copy of the GNU General Public License is available at
! http://www.r-project.org/Licenses/
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOTE: kind(1) and kind(1.0) return the respective bytes needed to represent
!       these constants. we use these values to declare integers and reals of
!       appropriate precision -- single precision!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! vsom !!!!!!
subroutine vsom(neurons,dt,dtrows,dtcols,xdim,ydim,alpha,train)
    implicit none

    !!! Input/Output
    ! NOTE: neurons are assumed to be initialized to small random values and then trained.
    integer(kind(1)),intent(in)   :: dtrows,dtcols,xdim,ydim,train
    real(kind(1.0)),intent(in)    :: alpha
    real(kind(1.0)),intent(inout) :: neurons(1:xdim*ydim,1:dtcols)
    real(kind(1.0)),intent(in)    :: dt(1:dtrows,1:dtcols)

    !!! Locals
    ! Note: the neighborhood cache is only valid as long as cache_counter < nsize_step
    integer(kind(1)) :: step_counter
    integer(kind(1)) :: nsize
    integer(kind(1)) :: nsize_step
    integer(kind(1)) :: epoch
    integer(kind(1)) :: i
    integer(kind(1)) :: ca(1)
    integer(kind(1)) :: c
    real(kind(1.0))  :: cache(1:xdim*ydim,1:xdim*ydim)       ! neighborhood cache
    logical          :: cache_valid(1:xdim*ydim)
    real(kind(1.0))  :: diff(1:xdim*ydim,1:dtcols)
    real(kind(1.0))  :: squ(1:xdim*ydim,1:dtcols)
    real(kind(1.0))  :: s(1:xdim*ydim)
    integer(kind(1)) :: coord_lookup(1:xdim*ydim,1:2)
    integer(kind(1)) :: ix
    real(kind(1.0))  :: ix_random
    integer(kind(1)) :: n

    !!! setup
    nsize = max(xdim,ydim) + 1
    nsize_step = ceiling((train*1.0)/nsize)
    step_counter = 0
    cache_valid = .false.

    ! fill the 2D coordinate lookup table that associates each
    ! 1D neuron coordinate with a 2D map coordinate
    do i=1,xdim*ydim
        call coord2D(coord_lookup(i,:),i,xdim)
    end do

    !!! training !!!
    ! the epochs loop
    do epoch=1,train

        step_counter = step_counter + 1
        if (step_counter == nsize_step) then
            step_counter = 0
            nsize = nsize - 1
            cache_valid = .false.
        endif

        ! select a training observation
        ix = mod(epoch,dtrows)+1

        !!! competetive step
        ! find the best matching neuron
        do i=1,dtcols
           diff(:,i) = neurons(:,i) - dt(ix,i)
        enddo
        squ = diff * diff
        call rowsums(s,squ,xdim*ydim,dtcols)
        ca = minloc(s)
        c = ca(1)

        !!! update step
        ! compute neighborhood vector
        call Gamma(cache(:,c),cache_valid,coord_lookup,nsize,xdim,ydim,c)

        do i=1,dtcols
           where (cache(:,c) > 0.0)
              neurons(:,i) = neurons(:,i) - (1.0 - real(epoch, kind=4)/real(train, kind=4)) * alpha * diff(:,i)
           endwhere
        enddo
    enddo
    return
end subroutine vsom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Gamma !!!!!!
subroutine Gamma(neighborhood,cache_valid,coord_lookup,nsize,xdim,ydim,c)
    implicit none

    ! parameters
    integer(kind(1)),intent(in)    :: nsize,xdim,ydim,c
    real(kind(1.0)),intent(inout)  :: neighborhood(1:xdim*ydim)
    logical,intent(inout)          :: cache_valid(1:xdim*ydim)
    integer(kind(1)),intent(in)    :: coord_lookup(1:xdim*ydim,1:2)

    ! locals
    integer(kind(1)) :: m
    integer(kind(1)) :: c2D(1:2),m2D(1:2)
    real(kind(1.0))  :: d

    ! cache is valid - nothing to do
    if (cache_valid(c)) then
        return
    endif

    ! convert the 1D neuron index into a 2D map index
    call coord2D(c2D,c,xdim)

    ! for each neuron m check if on the grid it is
    ! within the neighborhood.
    do m=1,xdim*ydim
        m2D = coord_lookup(m,:)
        d = sqrt(real(sum((c2D-m2D)**2)))
        if (d < nsize*1.5) then
            neighborhood(m) = 1.0
        else
            neighborhood(m) = 0.0
        end if
    end do

    ! cache it
    cache_valid(c) = .true.

    return
end subroutine Gamma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! rowsums !!!!!!
pure subroutine rowsums(s,v,rows,cols)
    implicit none

    ! parameters
    integer(kind(1)),intent(in) :: rows,cols
    real(kind(1.0)),intent(out) :: s(1:rows)
    real(kind(1.0)),intent(in)  :: v(1:rows,1:cols)

    ! locals
    integer(kind(1)) :: i

    s = 0.0

    do i = 1,cols
       s(:) = s(:) + v(:,i)
    enddo

    return
end subroutine rowsums


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! coord2D - convert a 1D rowindex into a 2D map coordinate !!!
pure subroutine coord2D(coord,ix,xdim)
    implicit none

    integer(kind(1)),intent(out) :: coord(1:2)
    integer(kind(1)),intent(in)  :: ix,xdim

    coord(1) = modulo(ix-1,xdim) + 1
    coord(2) = (ix-1)/xdim + 1

    return
end subroutine coord2D

