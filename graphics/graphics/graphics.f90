!  graphics.f90 
!
!  FUNCTIONS:
!  graphics - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: graphics
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program graphics
    use portlib
    implicit none
    
    ! Variables
    character*255 ::nomearq
    integer :: I, errnum, N
    INTEGER IO_STATUS, graf, ver
    character*50 ,allocatable ::files(:)
    ! Body of graphics
    I = system('dir /b *.gnu > list.txt')
    If (I .eq. -1) then
        errnum = ierrno( )
        print *, 'Error ', errnum  
    end if
    
    IO_STATUS=0
    n=0
    open(12,file='list.txt')
    do 
    read(12,*, IOSTAT=IO_STATUS) nomearq
        IF (IS_IOSTAT_END (IO_STATUS)) THEN
            !write(*,*) 'acabou'
            exit
        else
            n = n+1
        end if
    end do
    close(12)
 
    allocate(files(n))
    10 format(i4, t6 , '- ',a50)
    open(12,file='list.txt')
    do i=1,n
    read(12,*, IOSTAT=IO_STATUS) files(i)
        IF (IS_IOSTAT_END (IO_STATUS)) THEN
           ! write(*,*) 'acabou'
            exit
        end if
    end do
    close(12)
    
    do
    write(*,*) 'Digite o numero do grafico: '
        do i=1, n
            write(*,10) i, files(i)
        end do
    read(*,*) graf
    
    if (graf == -1) then
        exit
    else 
        if(graf >= 1 .and. graf <= n) then
            ver = system('WGNUPLOT.EXE ' // files(graf))
        end if
    end if
    
    
    end do
    
    
    deallocate(files)
    !ver = system('WGNUPLOT.EXE p_1002.gnu')

    end program graphics

