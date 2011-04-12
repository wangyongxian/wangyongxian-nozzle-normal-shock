module coeficientes
  
use variaveis  
use dados
  
implicit none
  
contains

!-------------------------------------------------

  subroutine lista_coeficientes_velocidade

    write(10,4)
    4 format(//,4x,'COEFICIENTES DA VELOCIDADE',//,  &
             t4,'volume',t13,'X',t34,'oeste',t55,'central', &
                     t76,'leste',t97,'fonte',/)
    do i = 1, N
       write(10,2) i, X(i), awu(i), aPu(i), aeu(i), bPu(i)
    end do

    2 format(i4,4x,5(1pe21.11))

  end subroutine lista_coeficientes_velocidade

!-------------------------------------------------

  subroutine lista_coeficientes_pressao

    write(10,10)
    10 format(//,4x,'COEFICIENTES DA PRESSÃO',//,  &
             t4,'volume',t13,'X',t34,'oeste',t55,'central', &
                     t76,'leste',t97,'fonte',/)
    do i = 1, N
       write(10,2) i, X(i), awplinha(i), aPplinha(i), aeplinha(i), bPplinha(i)
    end do

    2 format(i4,4x,5(1pe21.11))

  end subroutine lista_coeficientes_pressao

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml
	
	!real*8 teste1, teste2, teste3
	real*8 bpUDS, bpB
	!real*8,dimension(:),allocatable :: Mx, Mex, ap, aex, aw, bp, uexx, afux, atux, btux, btrux  ! fluxo de massa na face leste
	!allocate (Mx(N),Mex(N), ap(N), aex(N), aw(N), bp(N), uexx(N), afux(N), atux(N), btux(N), btrux(N))
	! Calculando o fluxo de massa na face leste
	do i = 1, N
	   Me(i) = ro*ue(i)*Ae(i) 
	   !Mex(i) = Me(i)
	   !Mx(i) = M(i)
	end do
   !uexx = ue
	! Fictício esquerdo P = 1
    awu(1) = 0.0d0
    aeu(1) = -1.0d0
    aPu(1) = 1.0d0
    bPu(1) = (2.0d0*x(2)/(x(3)-x(2)))*(u(3)-u(2))
	
	! volumes internos
    do i = 2, N-1
       !afu(i)  = (Pi/8.0d0)*ro*fator*u(i)*(2.0d0*raio(i))*deltax
	   !atu(i)  = M(i)/deltat
	   !btu(i)  = (M(i)*u_o(i))/deltat
	   !bpru(i) = - A(i)*(p(i+1)-p(i-1))/2.0d0

	   awu(i)  = roe(i-1)*ue(i-1)*Ae(i-1)
	   aeu(i)  = 0
	   
       aPu(i)  = -(awu(i) + aeu(i)) + rop(i)*A(i)*deltax/deltat
	   !bpu(i)  = btu(i) + bpru(i)
	   bpUDS = rop(i)*A(i)*u(i)*deltax/deltat+A(i)*(p(i-1)-p(i+1))/2.0d0-Pi*fator*rop(i)*(u(i)**2)*raio(i)*deltax/4.0d0
	   bpB = beta*(roe(i-1)*ue(i-1)*Ae(i-1)*(u(i)-u(i-1))-roe(i)*ue(i)*Ae(i)*(u(i+1)-u(i)))/2.0d0
	   bpu(i)  = bpUDS+bpB
	   !aex(i) = aeu(i)
	   !ap(i) = aPu(i)
	   !aw(i) = awu(i)
	   !bp(i) = bpu(i)
	end do
      
    ! Fictício direito P = N
    awu(N) = -1.0d0
    aeu(N) = 0.0d0
    aPu(N) = 1.0d0
    bPu(N) = (2.0d0*(x(N)-x(N-1))/(x(N-1)-x(N-2)))*(u(N-1)-u(N-2))
    
    !aex(1) = aeu(1)
    !aw(1) = awu(1)
    !ap(1) = aPu(1)
    !bp(1) = bpu(1)
    !afux = afu
    !atux = atu
    !btux = btu
    !btrux = bpru
    !aex(N) = aeu(N)
    !aw(N) = awu(N)
    !ap(N) = aPu(N)
    !bp(N) = bpu(N)
!write(*,*), 'cacareco'

  end subroutine coeficientes_e_fontes_qml

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia
	
	! Calculando o fluxo de massa na face leste
	do i = 1, N
	   Me(i) = ro*ue(i)*Ae(i) 
	end do

	! Fictício esquerdo P = 1
    awu(1) = 0.0d0
    aeu(1) = 1.0d0
    aPu(1) = 1.0d0
    bPu(1) = 2.0d0*T0
	
	! volumes internos
    do i = 2, N-1
	   awu(i)  = -cp*Me(i-1)
	   aeu(i)  = 0 
       aPu(i)  = cp*ro*A(i)*deltax/deltat+awu(i)
	   bPu(i)  = cp*A(i)*(p(i+1)-p(i-1))/2.0d0
	end do
      
    ! Fictício direito P = N
    awu(N) = -1.0d0
    aeu(N) = 0.0d0
    aPu(N) = 1.0d0
    bPu(N) = (2.0d0*(x(N)-x(N-1))/(x(N-1)-x(N-2)))*(T(N-1)-T(N-2))

  end subroutine coeficientes_e_fontes_energia
  
!-------------------------------------------------
  
  subroutine velocidades_ue

      real*8 treta
    real*8  :: sigmap, sigmaE ! auxiliar
   real*8,dimension(:),allocatable :: uex, ue0x      ! solução numéri
   allocate (uex(N), ue0x(N))
   
	! Calculando ue no passo deltat+1
	ue(1) = Uin
	
	do i = 2, N-2
	   sigmap = awu(i)*u(i-1) + aeu(i)*u(i+1)
	   sigmaE = awu(i+1)*u(i) + aeu(i+1)*u(i+2)
	   ue(i) = (sigmap+sigmaE+((M(i)+M(i+1))/deltat)*ue_o(i) - 2.0d0*Ae(i)*(p(i+1)-p(i)))/(apu(i)+apu(i+1))
	end do
	
	
	ue(N-1) = (u(N-1)+u(N))/2.0d0
	treta = Ae(N-1)
	Fat = (Uin*Ae(1))/(ue(N-1)*Ae(N-1))
  	ue(N-1) = Fat*ue(N-1)
	treta = ue(N-1)
	ue0x = ue_o
uex = ue
!write(*,*), 'eta'
  end subroutine velocidades_ue

!-------------------------------------------------

  subroutine coeficientes_simplec

    real*8,dimension(:),allocatable :: dpx, dex      ! solução numéri
   allocate (dpx(N),dex(N))
   
    do i = 2, N-1
	   ds(i) = A(i)/(apu(i)-awu(i)-aeu(i))
    end do

	do i = 2, N-2
	   de(i) = (ds(i) + ds(i+1))/2.0d0
	end do

	! Calculando nos contornos
	de(1) = ds(1+1)
	de(N-1) = ds(N-1)
	dpx = ds
	dex = de
!write(*,*), 'eta'

  end subroutine coeficientes_simplec

!-------------------------------------------------

  subroutine coeficientes_fontes_massa

    real*8 teste1, teste2
   !real*8,dimension(:),allocatable :: appl, awpl, awpl, bppl
   real*8,dimension(:),allocatable :: aepl, bppl, awpl, appl, uexx, Aexx
   !allocate (appl(N) , awpl(N)) , awpl(N) , bppl(N))
   allocate (aepl(N) , bppl(N), awpl(N) , appl(N), uexx(N) , Aexx(N))
   
    ! Calculando fictício P = 1	
    awplinha(1) = 0.0d0
    aeplinha(1) = 1.0d0
    aPplinha(1) = 1.0d0
    bpplinha(1) = 0.0d0

    ! Calculando os volumes internos
    do i = 2, N-1
      awplinha(i) = de(i-1)*Ae(i-1)
	   aeplinha(i) = de(i)*Ae(i)
	   aPplinha(i) = awplinha(i) + aeplinha(i)
	   teste1 = ue(i-1)*Ae(i-1)
	   teste2 = ue(i)*Ae(i)
	   bpplinha(i) = ue(i-1)*Ae(i-1) - ue(i)*Ae(i)
    end do
  
    ! Calculando fictício P = N 
    awplinha(N) = -1.0d0
    aeplinha(N) = 0.0d0
    aPplinha(N) = 1.0d0
    bpplinha(N) = (2.0d0*(x(N)-x(N-1))/(x(N-1)-x(N-2)))*(plinha(N-1)-plinha(N-2))
    uexx = ue
    appl = aPplinha
    awpl = awplinha
    aepl = aeplinha
    bppl = bpplinha
    Aexx = Ae
!write(*,*), 'eta'

  end subroutine coeficientes_fontes_massa

!-------------------------------------------------

  subroutine atualizar_ficticios_massa
real*8 teste1, teste2
	! Calculando fictício P = 1
	plinha(1) = 2.0d0*plinha(1+1) - plinha(1+2)

	! Calculando fictício P = N 
	plinha(N) = 2.0d0*plinha(N-1) - plinha(N-2)    
	teste1 = plinha(1)
	teste2 = plinha(N)
!write(*,*), 'eta'
  end subroutine atualizar_ficticios_massa 

!-------------------------------------------------

  subroutine corrigir_pressao
    
  real*8,dimension(:),allocatable :: opa, ppl
   allocate (opa(N), ppl(N))
    
	do i = 1, N
	   p(i) = p(i) + plinha(i)
	end do 
	opa = p
	ppl=plinha
!write(*,*), 'eta'

  end subroutine corrigir_pressao

!-------------------------------------------------
  
  subroutine corrigir_velocidades
     
   real*8,dimension(:),allocatable :: vel
   real*8 teste1, teste2
   allocate (vel(N))
	! Calculando fictício P = 1
	teste1 = u(2)
	teste2 = 2.0d0*Uin
	u(1) = - u(1+1) + 2.0d0*Uin
   
	! Calculando para volumes internos
	do i = 2, N-1
	   u(i) = u(i) - ds(i)*(plinha(i+1)-plinha(i-1))/2.0d0
	end do

	! Calculando fictício P = N
	u(N) = u(N-1) + (u(N-1)-u(N-2))
	vel = u
!write(*,*), 'eta'

  end subroutine corrigir_velocidades

!-------------------------------------------------

  subroutine corrigir_velocidades_faces

  !real*8 testex, teste2, teste3
   !real*8,dimension(:),allocatable :: faces
   !!allocate (faces(N) )
!faces = ue
	! Calculando ue somente para volumes internos
    do i = 2, N-2
 !   testex = ue(i)
  !  teste2 = plinha(i+1)
   ! teste3 = plinha(i)
	   ue(i) = ue(i) - de(i)*(plinha(i+1)-plinha(i))
	   
	end do 
	!faces = ue
!write(*,*), 'eta'

  end subroutine corrigir_velocidades_faces

!-------------------------------------------------
subroutine calculo_massa_especifica

    do i=2, N-1
        rop(i) = plinha(i)/(Rgases*T(i))
    end do    
    
end subroutine calculo_massa_especifica

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces
    
    do i=2, N-1
        roe(i) = rop(i)+Beta*(rop(i+1)-rop(i))/2.0d0
    end do
    
end subroutine calculo_massa_especifica_nas_faces

subroutine corrigir_massa_especifica
    do i=2, N-1
    rop(i) = rop(i)+P(i)/(Rgases*T(i))
    end do
end subroutine corrigir_massa_especifica

subroutine corrigir_massa_especifica_faces
    do i=2, N-1
    roe(i) = roe(i) + P(i)/(Rgases*T(i))
    end do
end subroutine corrigir_massa_especifica_faces


end module coeficientes

