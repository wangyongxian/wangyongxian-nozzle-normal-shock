module solvers_1D

use variaveis

implicit none

contains


!-------------------------------------------------

  ! método direto Tri-Diagonal Matrix Algorithm (TDMA)
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
  

 ! subroutine TDMA (N,ap,aw,ae,bp,Tp)

  !  implicit none

   ! integer :: i   ! número do nó
    !integer,intent(in) :: N ! número de nós

    !real*8,dimension(:),allocatable :: Pp ! coeficiente do tdma
    !real*8,dimension(:),allocatable :: Qp ! coeficiente do tdma

    !real*8,intent(in), dimension(N) :: ap ! coeficiente aP
    !real*8,intent(in), dimension(N) :: aw ! coeficiente aW
    !real*8,intent(in), dimension(N) :: ae ! coeficiente aE
    !real*8,intent(in), dimension(N) :: bp ! termo fonte bP

    !real*8,intent(out),dimension(N) :: Tp ! incógnita

    !allocate(Pp(N),Qp(N))
	
	!Pp(1) = ae(1) / ap(1)
    !Qp(1) = bp(1) / ap(1)
	
    !do i = 2, N
    !  Pp(i) = ae(i) / (ap(i) - aw(i)*Pp(i-1))
    !  Qp(i) = (bp(i) + aw(i)*Qp(i-1))/(ap(i) - aw(i)*Pp(i-1))
    !end do
	
	!Tp(N) = Qp(N)
   	
   ! do i = N-1, 1, -1
  !    Tp(i) = Pp(i)*Tp(i+1) + Qp(i)
 !   end do
	
!	deallocate(Pp,Qp)

!  end subroutine TDMA

!-------------------------------------------------
  
  ! calcula a norma l1 média do resíduo das equações

  subroutine norma (N,a,b,c,d,T,R)

    implicit none

    integer :: i ! número do nó

    integer,intent(in) :: N  ! número de nós

    real*8,intent(in), dimension(N) :: a ! coeficiente aP
    real*8,intent(in), dimension(N) :: b ! coeficiente aW
    real*8,intent(in), dimension(N) :: c ! coeficiente aE

	real*8,intent(in), dimension(N) :: d ! termo fonte bP

    real*8,intent(inout),dimension(N) :: T ! incógnita

    real*8,intent(inout) :: R ! norma l1 média dos resíduos

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
