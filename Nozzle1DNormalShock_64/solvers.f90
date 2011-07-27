module solvers_1D

use variaveis

implicit none

contains


!-------------------------------------------------

subroutine TDMA (N,ap,aw,ae,bp,T)

    implicit none

    integer,intent(in) :: n ! número de volumes de controle

    real*16,dimension(:),intent(in)  :: aw ! coeficiente oeste
    real*16,dimension(:),intent(in)  :: ap ! coeficiente central
    real*16,dimension(:),intent(in)  :: ae ! coeficiente leste
    real*16,dimension(:),intent(in)  :: bp ! termo independente

    real*16,dimension(:),intent(out) :: T  ! incógnita

    integer :: i   ! número do volume de controle
    real*16  :: div ! variável auxiliar

    real*16,dimension(:),allocatable :: p ! coeficiente do tdma
    real*16,dimension(:),allocatable :: q ! coeficiente do tdma

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

    real*16,dimension(:),intent(in) :: aw ! coeficiente oeste
    real*16,dimension(:),intent(in) :: ap ! coeficiente central
    real*16,dimension(:),intent(in) :: ae ! coeficiente leste
    real*16,dimension(:),intent(in) :: bp ! termo independente
    real*16,dimension(:),intent(in) :: T  ! incógnita

    real*16,intent(out) :: norma ! normal L1 do sistema linear

    integer :: i ! número do volume de controle

    norma = 0.0q0

    ! volume 1
    i = 1
       norma = norma + qabs (                ap(i)*T(i) &
                            + ae(i)*T(i+1) - bp(i)      )

    ! volumes internos
    do i = 2, n-1
       norma = norma + qabs ( aw(i)*T(i-1) + ap(i)*T(i) &
                            + ae(i)*T(i+1) - bp(i)      )
    end do

    ! volume n
    i = n
       norma = norma + qabs ( aw(i)*T(i-1) + ap(i)*T(i) &
                                           - bp(i)      )

  end subroutine Norma_L1

!-------------------------------------------------

      SUBROUTINE PENTA(N,E,A,D,C,F,B,X)

!   RESULTS:  matrix has 5 bands, EADCF, with D being the main diagonal,
!   E and A are the lower diagonals, and C and F are the upper diagonals.

!     E is defined for rows i = 3:N, but is defined as E(1) to E(N-2)
!     A is defined for rows i = 2:N, but is defined as A(1) to A(N-1)
!     D is defined for rows i = 1:N
!     C is defined for rows i = 1:N-1, but the last element isn't used
!     F is defined for rows i = 1:N-2, but the last 2 elements aren't used

!   B is the right-hand side
!   X is the solution vector

      implicit none
      integer i,n
      double precision E(N),A(N),D(N),C(N),F(N),B(N),X(N),XMULT
      DO 2 I = 2,N-1
        XMULT = A(I-1)/D(I-1)
        D(I) = D(I) - XMULT*C(I-1)
        C(I) = C(I) - XMULT*F(I-1)
        B(I) = B(I) - XMULT*B(I-1)
        XMULT = E(I-1)/D(I-1)
        A(I) = A(I) - XMULT*C(I-1)
        D(I+1) = D(I+1) - XMULT*F(I-1)
        B(I+1) = B(I+1) - XMULT*B(I-1)
   2  CONTINUE
      XMULT = A(N-1)/D(N-1)
      D(N) = D(N) - XMULT*C(N-1)
      X(N) = (B(N) - XMULT*B(N-1))/D(N)
      X(N-1) = (B(N-1) - C(N-1)*X(N))/D(N-1)
      DO 3 I = N-2,1,-1
        X(I) = (B(I) - F(I)*X(I+2) - C(I)*X(I+1))/D(I)
   3  CONTINUE
      RETURN
      
      END SUBROUTINE PENTA
      
      
end module solvers_1D
