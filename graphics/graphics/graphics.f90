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
integer :: ver
    ! Body of graphics
    ver = system('WGNUPLOT.EXE p_1002.gnu')

    end program graphics

