module solvers_1D

! objetivo: resolver sistema linear de equa��es alg�bricas
!           originado de problemas unidimensionais

use variaveis

implicit none

contains


!-------------------------------------------------

  ! m�todo direto Tri-Diagonal Matrix Algorithm (TDMA)

  subroutine TDMA (N,a,b,c,d,T)

    implicit none

    integer :: i   ! n�mero do n�
    real*8  :: div ! vari�vel auxiliar

    integer,intent(in) :: N ! n�mero de n�s

    real*8,dimension(:),allocatable :: P ! coeficiente do tdma
    real*8,dimension(:),allocatable :: Q ! coeficiente do tdma

    real*8,intent(in), dimension(N) :: a ! coeficiente aP
    real*8,intent(in), dimension(N) :: b ! coeficiente aW
    real*8,intent(in), dimension(N) :: c ! coeficiente aE
    real*8,intent(in), dimension(N) :: d ! termo fonte bP

    real*8,intent(out),dimension(N) :: T ! inc�gnita

    allocate(P(N),Q(N))
	
	P(1) = c(1) / a(1)
    Q(1) = d(1) / a(1)
	
    do i = 2, N
      div  = a(i) - b(i)*P(i-1)
      P(i) = c(i) / div
      Q(i) = (d(i) + b(i)*Q(i-1))/div
    end do
	
	T(N) = Q(N)
   	
    do i = N-1, 1, -1
      T(i) = P(i)*T(i+1) + Q(i)
    end do
	
	deallocate(P,Q)

  end subroutine tdma

!-------------------------------------------------
  
  ! calcula a norma l1 m�dia do res�duo das equa��es

  subroutine norma (N,a,b,c,d,T,R)

    implicit none

    integer :: i ! n�mero do n�

    integer,intent(in) :: N  ! n�mero de n�s

    real*8,intent(in), dimension(N) :: a ! coeficiente aP
    real*8,intent(in), dimension(N) :: b ! coeficiente aW
    real*8,intent(in), dimension(N) :: c ! coeficiente aE

	real*8,intent(in), dimension(N) :: d ! termo fonte bP

    real*8,intent(inout),dimension(N) :: T ! inc�gnita

    real*8,intent(inout) :: R ! norma l1 m�dia dos res�duos

    R = 0.0d0

	! Oeste
   	R = R + dabs ( c(1)*T(2) + d(1) - a(1)*T(1) )

	! Reais
    do i = 2, N-1
          R = R + dabs ( b(i)*T(i-1) + c(i)*T(i+1) + d(i) - a(i)*T(i) )
    end do

	! Leste
    R = R + dabs ( b(N)*T(N-1) + d(N) - a(N)*T(N) )

  end subroutine norma

!-------------------------------------------------

end module solvers_1D
