! data-parallel implementation of SOM training
! to compile this for R on a Mac
! gfortran-4.8  -fPIC -Wall -g -Ofast -cpp  -c  vsom.f90 -o vsom.o
! gfortran-4.8 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o vsom.so vsom.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation

! define the following to emit debugging code
!#define debug 1

! define this for pure vector/matrix based operations
#define vector 1

! define this for standard do-loop computations
!#define iterative 1

!!!!!! vsom !!!!!!
subroutine vsom(neurons,dt,dtrows,dtcols,xdim,ydim,alpha,train)
    implicit none

    !!! Input/Output
    ! neurons are initialized to small random values and then trained.
    real*4,intent(inout) :: neurons(1:xdim*ydim,1:dtcols)

    !!! Input
    real*4,intent(in) :: dt(1:dtrows,1:dtcols)
    integer,intent(in) :: dtrows,dtcols,xdim,ydim,train
    real*4,intent(in) :: alpha

    !!! Locals
    integer :: cache_counter
    integer :: nsize
    integer :: nsize_step
    integer :: epoch
    integer :: k,i
    integer :: ca(1)
    integer :: c
    real*4  :: cache(1:xdim*ydim,1:xdim*ydim)       ! neighborhood cache
    logical :: cache_valid(1:xdim*ydim)
    real*4  :: diff(1:xdim*ydim,1:dtcols)
    real*4  :: squ(1:xdim*ydim,1:dtcols)
    real*4  :: s(1:xdim*ydim)
#ifdef iterative
    integer :: j
#endif

#ifdef debug
    open(unit=1,file="debug.txt",form="formatted",status="replace",action="write")
    write(1,*) 'dtrows',dtrows,'dtcols',dtcols,'xdim',xdim,'ydim',ydim,'alpha',alpha,'train',train
    write(1,*) 'init neuron matrix'
    call write_array(1,neurons,xdim*ydim,dtcols,'f7.3')
#endif

    nsize = max(xdim,ydim) + 1
    nsize_step = ceiling((train*1.0)/nsize)
    cache_valid = .false.
    cache_counter = 0

    !!! training !!!
    ! the epochs loop
    do epoch=1,train

#ifdef debug
        write(1,*) 'Epoch',epoch,'Neighborbood',nsize
#endif

        cache_counter = cache_counter + 1
        if (cache_counter == nsize_step) then
            cache_counter = 0
            nsize = nsize - 1
            cache_valid = .false.
        endif

        !!! run through the training set
        do k=1,dtrows
            ! neuron local computation
#ifdef vector
            do i=1,dtcols
               diff(:,i) = neurons(:,i) - dt(k,i)
            enddo
#endif
#ifdef iterative
            do i=1,dtcols
               do j=1,xdim*ydim
                  diff(j,i) = neurons(j,i) - dt(k,i)
               enddo
            enddo
#endif

#ifdef vector
            squ = diff * diff
#endif
#ifdef iterative
            do i=1,dtcols
               do j=1,xdim*ydim
                 squ(j,i) = diff(j,i) * diff(j,i)
               enddo
            enddo
#endif

            call rowsums(s,squ,xdim*ydim,dtcols)

            ! reduce
            ca = minloc(s)
            c = ca(1)

            !!! update step
            ! compute neighborhood vector
            call Gamma(cache(:,c),cache_valid,nsize,xdim,ydim,c)

#ifdef debug
            write(1,*) 'neighborhood cache for',c
            call write_array(1,cache(:,c),xdim,ydim,'f2.0')
#endif

#ifdef vector
            do i=1,dtcols
               where (cache(:,c) > 0.0) 
                  neurons(:,i) = neurons(:,i) - alpha * diff(:,i)
               endwhere
            enddo
#endif
#ifdef iterative
            do i=1,dtcols
               do j=1,xdim*ydim
                  if (cache(j,c) > 0.0) then
                    neurons(j,i) = neurons(j,i) - alpha * diff(j,i)
                  endif
               enddo
            enddo
#endif
        enddo
    enddo

#ifdef debug
    write(1,*) 'trained neuron matrix'
    call write_array(1,neurons,xdim*ydim,dtcols,'f7.3')
    close(unit=1)
#endif

    return
end subroutine vsom

!!!!!! Gamma !!!!!!
subroutine Gamma(neighborhood,cache_valid,nsize,xdim,ydim,c)
    implicit none

    ! parameters
    ! Note: in the cache a neighborhood is a vector, here we
    ! reshaping the vector into a 2D matrix on the fly
    real*4,intent(inout)  :: neighborhood(1:xdim,1:ydim)
    logical,intent(inout) :: cache_valid(1:xdim*ydim)
    integer,intent(in)    :: nsize,xdim,ydim,c

    integer :: xc,yc
    real*4  :: x_v(1:xdim)
    integer :: i
    integer :: y_lb,y_ub
#ifdef iterative
    integer :: j
#endif

    ! cache is valid - nothing to do
    if (cache_valid(c)) then
#ifdef debug
        write(1,*) 'cache hit',xc,yc
#endif
        return
    endif

    ! convert the 1D neuron index into a 2D map index
    xc = modulo(c-1,xdim) + 1
    yc = (c-1)/xdim + 1

    ! take care of simple cases
    if (nsize >= xdim .and. nsize >= ydim) then
        neighborhood = 1
    else if (nsize >= ydim) then
        call build_xvector(x_v,nsize,xc,xdim)

#ifdef vector
        do i=1,ydim
            neighborhood(:,i) = x_v
        enddo
#endif
#ifdef iterative
        do i=1,ydim
           do j=1,xdim
              neighborhood(j,i) = x_v(j)
           enddo
        enddo
#endif
    else ! full blown neighborhood construction
        neighborhood = 0
        call build_xvector(x_v,nsize,xc,xdim)

        ! compute y_lb
        if (yc - nsize + 1 < 1) then
            y_lb = 1
        else
            y_lb = yc - nsize + 1
        endif

        ! compute y_ub
        if (yc + nsize - 1 > ydim) then
            y_ub = ydim
        else
            y_ub = yc + nsize - 1
        endif
            
        ! construct the neighborhood
#ifdef vector
        do i=y_lb,y_ub
            neighborhood(:,i) = x_v
        enddo
#endif
#ifdef iterative
        do i=y_lb,y_ub
           do j=1,xdim
              neighborhood(j,i) = x_v(j)
           enddo
        enddo
#endif
    endif

#ifdef debug
    write(1,*) 'computing neighborhood matrix for',xc,yc
    call write_array(1,neighborhood,xdim,ydim,'f2.0')
#endif

    ! cache it
    cache_valid(c) = .true.

    return
end subroutine Gamma

!!!!!! build_xvector !!!!!!
pure subroutine build_xvector(x_v,nsize,xc,xdim)
    implicit none

    ! parameters
    real*4,intent(out) :: x_v(1:xdim)
    integer,intent(in) :: nsize,xc,xdim

    integer :: x_lb,x_ub,i

    x_v = 0

    ! compute x_lb
    if (xc - nsize + 1 < 1) then
        x_lb = 1
    else
        x_lb = xc - nsize + 1
    endif

    ! compute x_ub
    if (xc + nsize - 1 > xdim) then
        x_ub = xdim
    else
        x_ub = xc + nsize - 1
    endif

    ! construct the xdim neighborhood
    do i = x_lb,x_ub
        x_v(i) = 1
    enddo

    return
end subroutine build_xvector

!!!!!! rowsums !!!!!!
pure subroutine rowsums(s,v,rows,cols)
    implicit none

    ! parameters
    real*4,intent(out) :: s(1:rows)
    real*4,intent(in)  :: v(1:rows,1:cols)
    integer,intent(in) :: rows,cols

    integer :: i
#ifdef iterative
    integer :: j
#endif

    s = 0.0

#ifdef vector
    do i = 1,cols
       s(:) = s(:) + v(:,i)
    enddo
#endif
#ifdef iterative
    do j = 1,cols
       do i=1,rows
          s(i) = s(i) + v(i,j)
       enddo
    enddo
#endif

    return
end subroutine rowsums


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

