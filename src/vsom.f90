!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! vectorized version of the stochastic SOM training algorithm
! written by Lutz Hamel, University of Rhode Island (c) 2016
!
! NOTE: the OPENMP code is experimental with a very corse view of parallelism.
! Support for OPENMP in R is very limited and therefore should not be used when
! compiling for R.
!
! License:
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation.
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! define the following to emit debugging code (only tested in single thread mode)
!#define debug 1

!!!!!! vsom !!!!!!
subroutine vsom(neurons,dt,dtrows,dtcols,xdim,ydim,alpha,train)
    !$ use omp_lib
    implicit none

    !!! Input/Output
    ! neurons are initialized to small random values and then trained.
    real*4,intent(inout) :: neurons(1:xdim*ydim,1:dtcols)

    !!! Input
    real*4,intent(in) :: dt(1:dtrows,1:dtcols)
    integer,intent(in) :: dtrows,dtcols,xdim,ydim,train
    real*4,intent(in) :: alpha

    !!! Locals
    ! Note: the neighborhood cache is only valid as long as cache_counter < nsize_step
    integer :: step_counter
    integer :: nsize
    integer :: nsize_step
    integer :: epoch
    integer :: i
    integer :: ca(1)
    integer :: c
    real*4  :: cache(1:xdim*ydim,1:xdim*ydim)       ! neighborhood cache
    logical :: cache_valid(1:xdim*ydim)
    real*4  :: diff(1:xdim*ydim,1:dtcols)
    real*4  :: squ(1:xdim*ydim,1:dtcols)
    real*4  :: s(1:xdim*ydim)
    integer :: coord_lookup(1:xdim*ydim,1:2)
    integer :: ix
    real*4  :: ix_random

    !$OMP THREADPRIVATE(i,ca,c,diff,squ,s,xi,ix_random)

#ifdef debug
    open(unit=1,file="debug.txt",form="formatted",status="replace",action="write")
    write(1,*) 'dtrows',dtrows,'dtcols',dtcols,'xdim',xdim,'ydim',ydim,'alpha',alpha,'train',train
    write(1,*) 'init neuron matrix'
    call write_array(1,neurons,xdim*ydim,dtcols,'f7.3')
#endif

    !!! setup
    nsize = max(xdim,ydim) + 1
    nsize_step = ceiling((train*1.0)/nsize)
    step_counter = 0
    cache_valid = .false.
    call random_seed()

    !$ call OMP_set_num_threads(2)
    !$OMP PARALLEL
    !$ print *, 'no. of threads: ', OMP_get_num_threads() 

    ! fill the 2D coordinate lookup table that associates each
    ! 1D neuron coordinate with a 2D map coordinate
    !$OMP DO
    do i=1,xdim*ydim
        call coord2D(coord_lookup(i,:),i,xdim)
    end do
    !$OMP END DO


    !!! training !!!
    ! the epochs loop
    !$OMP DO
    do epoch=1,train

#ifdef debug
        write(1,*) 'Epoch',epoch,'Neighborbood',nsize
#endif
        !$OMP CRITICAL step
        ! check if we are at the end of a step
        step_counter = step_counter + 1
        if (step_counter == nsize_step) then
            step_counter = 0
            nsize = nsize - 1
            cache_valid = .false.
        endif
        !$OMP END CRITICAL step

        ! select a training observation
        call random_number(ix_random)
        ix = 1 + int(ix_random*dtrows)

        !!! learn the training observation
        ! neuron local computation
        do i=1,dtcols
           diff(:,i) = neurons(:,i) - dt(ix,i)
        enddo
        squ = diff * diff
        call rowsums(s,squ,xdim*ydim,dtcols)

        ! reduce
        ca = minloc(s)
        c = ca(1)

        !!! update step
        ! compute neighborhood vector
        !$OMP CRITICAL cache
        call Gamma(cache(:,c),cache_valid,coord_lookup,nsize,xdim,ydim,c)
        !$OMP END CRITICAL cache

#ifdef debug
        write(1,*) 'neighborhood cache for',c
        call write_array(1,cache(:,c),xdim,ydim,'f2.0')
#endif

        do i=1,dtcols
           where (cache(:,c) > 0.0) 
              !$OMP CRITICAL neurons
              neurons(:,i) = neurons(:,i) - alpha * diff(:,i)
              !$OMP END CRITICAL neurons
           endwhere
        enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

#ifdef debug
    write(1,*) 'trained neuron matrix'
    call write_array(1,neurons,xdim*ydim,dtcols,'f7.3')
    close(unit=1)
#endif

    return
end subroutine vsom

!!!!!! Gamma !!!!!!
subroutine Gamma(neighborhood,cache_valid,coord_lookup,nsize,xdim,ydim,c)
    implicit none

    ! parameters
    real*4,intent(inout)  :: neighborhood(1:xdim*ydim)
    logical,intent(inout) :: cache_valid(1:xdim*ydim)
    integer,intent(in)    :: coord_lookup(1:xdim*ydim,1:2)
    integer,intent(in)    :: nsize,xdim,ydim,c

    ! locals
    integer :: m
    integer :: c2D(1:2),m2D(1:2)
    real*4  :: d

    ! cache is valid - nothing to do
    if (cache_valid(c)) then
#ifdef debug
        write(1,*) 'cache hit',xc,yc
#endif
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

#ifdef debug
    call write_array(1,neighborhood,xdim,ydim,'f2.0')
#endif

    ! cache it
    cache_valid(c) = .true.

    return
end subroutine Gamma


!!!!!! rowsums !!!!!!
pure subroutine rowsums(s,v,rows,cols)
    implicit none

    ! parameters
    real*4,intent(out) :: s(1:rows)
    real*4,intent(in)  :: v(1:rows,1:cols)
    integer,intent(in) :: rows,cols

    integer :: i

    s = 0.0

    do i = 1,cols
       s(:) = s(:) + v(:,i)
    enddo

    return
end subroutine rowsums

!!! convert a 1D rowindex into a 2D map coordinate
pure subroutine coord2D(coord,ix,xdim)
    implicit none

    integer,intent(out) :: coord(1:2)
    integer,intent(in) :: ix,xdim

    coord(1) = modulo(ix-1,xdim) + 1
    coord(2) = (ix-1)/xdim + 1

    return
end subroutine coord2D

#ifdef debug
!!!!!!!!!!!!!!!!!!!
!!! debug stuff !!!
!!!!!!!!!!!!!!!!!!!

subroutine write_array(unit,a,x,y,efmt)
    implicit none

    integer,intent(in) :: unit,x,y
    real*4,intent(in) :: a(x,y)
    character(len=4) :: efmt

    integer :: i
    character(len=5) :: num
    character(len=80) :: fmt

    write(num,'(i5)') y
    num = trim(num)
    efmt = trim(efmt)
    fmt = '('//num//efmt//')'

    do i=1,x
        write(unit,fmt) a(i,:)
    enddo


    return
end subroutine write_array
#endif

