module solvers_1D

use variaveis

implicit none

contains


!-------------------------------------------------

subroutine TDMA (N,ap,aw,ae,bp,T)

    implicit none

    integer,intent(in) :: n ! número de volumes de controle

    real*8,dimension(:),intent(in)  :: aw ! coeficiente oeste
    real*8,dimension(:),intent(in)  :: ap ! coeficiente central
    real*8,dimension(:),intent(in)  :: ae ! coeficiente leste
    real*8,dimension(:),intent(in)  :: bp ! termo independente

    real*8,dimension(:),intent(out) :: T  ! incógnita

    integer :: i   ! número do volume de controle
    real*8  :: div ! variável auxiliar

    real*8,dimension(:),allocatable :: p ! coeficiente do tdma
    real*8,dimension(:),allocatable :: q ! coeficiente do tdma

    allocate ( p(n), q(n) )

    p(1) = - ae(1) / ap(1)
    q(1) = bp(1) / ap(1)

    do i = 2, n
       div = ap(i) + aw(i)*p(i-1)
       p(i) = - ae(i) / div
       q(i) = (bp(i) - aw(i)*q(i-1))/div
    end do

    T(n) = q(n)

    do i = n-1, 1, -1
       T(i) = p(i)*T(i+1) + q(i)
    end do

    deallocate ( p, q )

  end subroutine TDMA
  
!-------------------------------------------------
  ! calcula a norma l1 média do resíduo das equações
 subroutine Norma_L1 ( n, aw, ap, ae, bp, T, norma )

    implicit none

    integer,intent(in) :: n ! número de volumes de controle

    real*8,dimension(:),intent(in) :: aw ! coeficiente oeste
    real*8,dimension(:),intent(in) :: ap ! coeficiente central
    real*8,dimension(:),intent(in) :: ae ! coeficiente leste
    real*8,dimension(:),intent(in) :: bp ! termo independente
    real*8,dimension(:),intent(in) :: T  ! incógnita

    real*8,intent(out) :: norma ! normal L1 do sistema linear

    integer :: i ! número do volume de controle

    norma = 0.0d0

    ! volume 1
    i = 1
       norma = norma + dabs (                ap(i)*T(i) &
                            + ae(i)*T(i+1) - bp(i)      )

    ! volumes internos
    do i = 2, n-1
       norma = norma + dabs ( aw(i)*T(i-1) + ap(i)*T(i) &
                            + ae(i)*T(i+1) - bp(i)      )
    end do

    ! volume n
    i = n
       norma = norma + dabs ( aw(i)*T(i-1) + ap(i)*T(i) &
                                           - bp(i)      )

  end subroutine Norma_L1

!-------------------------------------------------

end module solvers_1D
