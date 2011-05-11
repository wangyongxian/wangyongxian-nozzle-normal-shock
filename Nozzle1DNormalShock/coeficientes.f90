module coeficientes
  
use variaveis  
use dados
  
implicit none
  
contains

!-------------------------------------------------

  !subroutine lista_coeficientes_velocidade

   ! write(10,4)
   ! 4 format(//,4x,'COEFICIENTES DA VELOCIDADE',//,  &
    !         t4,'volume',t13,'X',t34,'oeste',t55,'central', &
     !                t76,'leste',t97,'fonte',/)
    !do i = 1, N
     !  write(10,2) i, X(i), awu(i), aPu(i), aeu(i), bPu(i)
    !end do

    !2 format(i4,4x,5(1pe21.11))

  !end subroutine lista_coeficientes_velocidade

!-------------------------------------------------

  !subroutine lista_coeficientes_pressao

   ! write(10,10)
    !10 format(//,4x,'COEFICIENTES DA PRESSÃO',//,  &
     !        t4,'volume',t13,'X',t34,'oeste',t55,'central', &
      !               t76,'leste',t97,'fonte',/)
    !do i = 1, N
     !  write(10,2) i, X(i), awplinha(i), aPplinha(i), aeplinha(i), bPplinha(i)
   ! end do

    !2 format(i4,4x,5(1pe21.11))

  !end subroutine lista_coeficientes_pressao

!-------------------------------------------------


  subroutine coeficientes_e_fontes_qml
	
	real*8 :: bpUDS, bpB
	
	! Fictício esquerdo P = 1
	!fator = 2 * xp(2) / ( xp(3) - xp(2) )
    !aw(1) =  0.0d0
    !ap(1) =  1.0d0
    !ae(1) = -1.0d0
    !bp(1) = -fator * ( u(3) - u(2) )
    !bf(1) = 0.0d0
    
    aw(1) = 0.0d0
    ae(1) = -1.0d0
    ap(1) = 1.0d0
    bp(1) = -(2.0d0*xp(2)/(xp(3)-xp(2)))*(u(3)-u(2))
	
	! volumes internos
    do i = 2, N-1
    !aw(i) = - rom(i-1) * um(i-1) * se(i-1)
	   aw(i)  = -roe(i-1)*ue(i-1)*Se(i-1)
	   ae(i)  = 0.0d0
	   !ap(i) = roa(i) * sp(i) * dx / dt - ( aw(i) + ae(i) )
       ap(i)  = ro_o(i)*Sp(i)*dx/dt - (aw(i) + ae(i))
	   !bf(i) = - pi * f(i) * ro(i) * (u(i)**2) * rp(i) * dx / 4   
	   !bp(i) = roa(i) * sp(i) * dx * ua(i) / dt + bf(i)   &
       !      + 0.5d0  * sp(i) * ( p(i-1) - p(i+1) )
	   bpUDS = ro_o(i)*Sp(i)*u_o(i)*dx/dt &
	                + 0.5d0*Sp(i)*(p(i+1)-p(i-1))  &
	                - 0.25d0*Pi*fator*ro(i)*(u(i)**2)*raio(i)*dx
        
       !oeste = rom(i-1) * um(i-1) * se(i-1) * ( u(i) - u(i-1) )
       !leste = rom(i) * um(i) * se(i) * ( u(i+1) - u(i) )
       !bc(i) = 0.5d0 * beta * ( oeste - leste )
       	                
	   bpB = 0.5d0*beta*(roe(i-1)*ue(i-1)*Se(i-1)*(u(i)-u(i-1)) &
	            -roe(i)*ue(i)*Se(i)*(u(i+1)-u(i)))
	   bp(i) = bpUDS + bpB
	end do
      
    ! Fictício direito P = N
    !fator = 2 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    !aw(n) = -1.0d0
    !ap(n) =  1.0d0
    !ae(n) =  0.0d0
    !bp(n) =  fator * ( u(n-1) - u(n-2) )
    !bf(n) =  0.0d0
    
    aw(N) = -1.0d0
    ae(N) = 0.0d0
    ap(N) = 1.0d0
    bp(N) = (2.0d0*(xp(N)-xp(N-1))/(xp(N-1)-xp(N-2)))*(u(N-1)-u(N-2))
    
  end subroutine coeficientes_e_fontes_qml

!-------------------------------------------------
  
  subroutine calculo_velocidades_face

    real*8  :: sigmap, sigmaE, bcP, bcE, bfP, bfE ! auxiliar
   
	! Calculando ue no passo dt+1
	!um(1) = 0.5d0 * ( u(1) + u(2) )
	ue(1) = 0.5d0*(u(1)+u(2))
	
	do i = 2, N-2
	   !massa_p = roa(i) * sp(i) * (xe(i) - xe(i-1))
       !massa_e = roa(i+1) * sp(i+1) * (xe(i+1) - xe(i))
       !somap = aw(i)*u(i-1) + ae(i)*u(i+1)
       !somae = aw(i+1)*u(i) + ae(i+1)*u(i+2)
       !um(i) = (-somap - somae + bc(i) + bc(i+1) + bf(i) + bf(i+1)       &
       !      + (massa_p+massa_e)*uma(i)/dt - 2.0d0*se(i)*(p(i+1)-p(i)))  &
       !      / (ap(i)+ap(i+1))

	   bcP = beta*((roe(i-1)*Se(i-1)*ue(i-1))*(u(i)-u(i-1))-(roe(i)*Se(i)*ue(i))*(u(i+1)-u(i)))/2.0d0
	   bcE = beta*((roe(i)*Se(i)*ue(i))*(u(i+1)-u(i))-(roe(i+1)*Se(i+1)*ue(i+1))*(u(i+2)-u(i+1)))/2.0d0
	   ! arrumar o raio hidraulico do volume de controle
	   bfP = 0 !-Pi*fator*ro(i)*u(i)*raio(i)*dx/4.0d0
	   bfE = 0 !-Pi*fator*ro(i+1)*u(i+1)*raio(i+1)*dx/4.0d0
	   sigmap = aw(i)*u(i-1) + ae(i)*u(i+1)
	   sigmaE = aw(i+1)*u(i) + ae(i+1)*u(i+2)
	   !ue(i) = (-sigmap-sigmaE+bcP+bcE+bfP+bfE+((ro_o(i)*Sp(i)*dx+ro_o(i+1)*Sp(i+1)*dx)/dt)*ue_o(i) - 2.0d0*Se(i)*(p(i+1)-p(i)))/(apu(i)+apu(i+1))
	   ue(i) = (-sigmap-sigmaE+(ro_o(i)*Sp(i)*dx+ro_o(i+1)*Sp(i+1)*dx)*ue_o(i)/dt - 2.0d0*Se(i)*(p(i+1)-p(i)))/(ap(i)+ap(i+1))
	end do
	
	
    !um(n-1) = 0.5d0 * ( u(n) + u(n-1) )
	ue(N-1) = 2.0d0*(u(N-1)+u(N))
	!arrumar o uin u(1)>uin
	!Fat = (u(1)*Se(1))/(ue(N-1)*Se(N-1))
  	!ue(N-1) = Fat*ue(N-1)
 
  end subroutine calculo_velocidades_face

!-------------------------------------------------


  subroutine coeficientes_e_fontes_energia
	
	real*8 :: bUDS, bBeta

	! Fictício esquerdo P = 1
	!aw(1) = 0.0d0
    !ap(1) = 1.0d0
    !ae(1) = 1.0d0
    !bp(1) = 2.0d0 * T_in
    
    aw(1) = 0.0d0
    ae(1) = 1.0d0
    ap(1) = 1.0d0
    bp(1) = 2.0d0*T_in
	
	! volumes internos
    do i = 2, N-1
	   aw(i)  = -cp*roe(i-1)*ue(i-1)*Se(i-1)
	   ae(i)  = 0 
       ap(i)  = cp*Sp(i)*ro_o(i)*dx/dt - (aw(i)+ae(i))
!       aw(i) = - cp(i) * rom(i-1) * um(i-1) * se(i-1)       
!       ap(i) = cp(i) * roa(i) * sp(i) * dx / dt - ( aw(i) + ae(i) ) &
!             + h(i) * C_rec(i) * sn(i)                              &
!             + epsilon_cte * sigma * sn(i) * ( T(i) ** 3 )    

!       bp(i) = cp(i) * roa(i) * sp(i) * dx * Ta(i) / dt        &
!             + sp(i) * dx * ( p(i) - pa(i) ) / dt              &
!             + sp(i) * u(i) * ( p(i+1) - p(i-1) ) * 0.5d0      &
!             + pi * f(i) * ro(i) * (u(i)**3) * rp(i) * dx / 4  &
!             + h(i) * T_wall(i) * sn(i)                        &
!             + epsilon_cte * sigma * sn(i) * ( T_wall(i) ** 4 )
       
       bUDS = cp*Sp(i)*ro_o(i)*T_o(i)*dx/dt &
                +Sp(i)*(p(i)-p_o(i))*dx/dt &
                +Sp(i)*u(i)*(p(i+1)-p(i-1))*0.5d0 &
                +0.25d0*Pi*fator*ro(i)*raio(i)*dx*u(i)**3
       !oeste = rom(i-1) * um(i-1) * se(i-1) * ( T(i) - T(i-1) )
       !leste = rom(i) * um(i) * se(i) * ( T(i+1) - T(i) )
       !bc(i) = 0.5d0 * beta * cp(i) * ( oeste - leste )
       	   
	   bBeta = 0.5d0*beta*cp*(roe(i-1)*ue(i-1)*Se(i-1)*(T(i)-T(i-1)) &
	                -roe(i)*ue(i)*Se(i)*(T(i+1)-T(i)))
	                
	   bp(i)  = bUDS + bBeta
	end do
     
    !fator = 2 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    !aw(n) = -1.0d0
    !ap(n) =  1.0d0
    !ae(n) =  0.0d0
    !bp(n) =  fator * ( T(n-1) - T(n-2) )
    
    ! Fictício direito P = N
    aw(N) = -1.0d0
    ae(N) = 0.0d0
    ap(N) = 1.0d0
    bp(N) = (2.0d0*(xp(N)-xp(N-1))/(xp(N-1)-xp(N-2)))*(T(N-1)-T(N-2))

  end subroutine coeficientes_e_fontes_energia
  
!-------------------------------------------------

  subroutine coeficientes_fontes_massa

    real*8 :: bUDS, bBeta
    !aw(1) = 0.0d0
    !ap(1) = 1.0d0
    !ae(1) = 1.0d0
    !bp(1) = 2.0d0 * pl_in
    
    ! Calculando fictício P = 1	
    aw(1) = 0.0d0
    ae(1) = 1.0d0
    ap(1) = 1.0d0
    !arrumar
    bp(1) = 2.0d0*pl_in

    ! Calculando os volumes internos
    do i = 2, N-1
    
      ! aw(i) = - rom(i-1) * de(i-1) * se(i-1)          &
      !        - um(i-1)  * se(i-1) / ( Rg(i-1) * T(i-1) )

      ! ae(i) = - rom(i) * de(i) * se(i)

      ! ap(i) = ( sp(i) * dx / dt + um(i) * se(i) ) / ( Rg(i) * T(i) )   &
      !       + rom(i-1) * de(i-1) * se(i-1)                              &
      !       + rom(i)   * de(i)   * se(i)

      ! bp(i) = - ( ( ro(i) - roa(i) ) * sp(i) * dx / dt              &
      !       + ro(i) * um(i) * se(i) - ro(i-1) * um(i-1) * se(i-1) )
             
             
       aw(i) = -roe(i-1)*Se(i-1)*de(i-1) &
                        -ue(i-1)*Se(i-1)/(rgases*T(i-1)) 
	   ae(i) = -roe(i)*Se(i)*de(i)
	   ap(i) = (Sp(i)*dx/dt + ue(i)*Se(i))/(rgases*T(i)) &
	                    +roe(i-1)*de(i-1)*Se(i-1) &
	                    +roe(i)*de(i)*Se(i) 
	                    
	   bUDS = -((ro(i)-ro_o(i))*Sp(i)*dx/dt &
	            +ro(i)*Se(i)*ue(i) &
	            -ro(i-1)*Se(i-1)*ue(i-1))
       !oeste = se(i-1) * ( ro(i) - ro(i-1) ) * um(i-1)
       !leste = se(i)   * ( ro(i+1) - ro(i) ) * um(i)
       !bc(i) = 0.5d0 * beta * ( oeste - leste )
	   bBeta = 0.5d0*beta*(Se(i-1)*(ro(i)-ro(i-1))*ue(i-1)-Se(i)*(ro(i+1)-ro(i))*ue(i))
	   
	   bp(i) = bUDS+bBeta
    end do
  
    !fator = 2 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
    !aw(n) = -1.0d0
    !ap(n) =  1.0d0
    !ae(n) =  0.0d0
    !bp(n) =  fator * ( pl(n-1) - pl(n-2) )
    
    ! Calculando fictício P = N 
    aw(N) = -1.0d0
    ae(N) = 0.0d0
    ap(N) = 1.0d0
    bp(N) = (2.0d0*(xp(N)-xp(N-1))/(xp(N-1)-xp(N-2)))*(pl(N-1)-pl(N-2))

  end subroutine coeficientes_fontes_massa

!-------------------------------------------------

  subroutine coeficientes_simplec

    ds(1) = 0.0d0
    ds(N) = 0.0d0
    
    do i = 2, N-1
	   ds(i) = Sp(i)/(ap(i)+aw(i)+ae(i))
    end do

    ! Calculando nos contornos
	de(1) = ds(2)
	de(N-1) = ds(N-1)
	
	do i = 2, N-2
	   de(i) = 0.5d0*(ds(i) + ds(i+1))
	end do

	

  end subroutine coeficientes_simplec


!-------------------------------------------------

  subroutine atualizar_ficticios_massa
	! Calculando fictício P = 1
	pl(1) = 2.0d0*pl(1+1) - pl(1+2)

	! Calculando fictício P = N 
	pl(N) = 2.0d0*pl(N-1) - pl(N-2)    
  
  end subroutine atualizar_ficticios_massa 

!-------------------------------------------------

  subroutine corrigir_pressao
    
	do i = 1, N
	   p(i) = p(i) + pl(i)
	end do 

  end subroutine corrigir_pressao

!-------------------------------------------------
  
  subroutine corrigir_velocidades
     !arrumar
	u(1) = - u(2) + 2.0d0*u(1)
   
	! Calculando para volumes internos
	do i = 2, N-1
	   u(i) = u(i) - ds(i)*(pl(i+1)-pl(i-1))/2.0d0
	end do

	! Calculando fictício P = N
	u(N) = u(N-1) + (u(N-1)-u(N-2))

  end subroutine corrigir_velocidades

!-------------------------------------------------

  subroutine corrigir_velocidades_faces

    do i = 2, N-2
	   ue(i) = ue(i) - de(i)*(pl(i+1)-pl(i))
	end do 

  end subroutine corrigir_velocidades_faces

!-------------------------------------------------

  subroutine corrigir_massa_especifica

    do i = 1, N
	   ro(i) = ro(i) + pl(i)/(Rgases*T(i))
	end do 

  end subroutine corrigir_massa_especifica
  
!-------------------------------------------------

subroutine calculo_massa_especifica

    do i=1, N
        ro(i) = p(i)/(Rgases*T(i))
    end do    
    
end subroutine calculo_massa_especifica

!-------------------------------------------------

subroutine calculo_massa_especifica_nas_faces
    
    do i=1, N-1
        !rom(i) = ro(i) + beta * 0.5d0 * ( ro(i+1) - ro(i) )
        roe(i) = ro(i) + beta * 0.5d0 * (ro(i+1)-ro(i))
    end do
    
end subroutine calculo_massa_especifica_nas_faces

!-------------------------------------------------

subroutine calcula_fluxo_massa
    do i=1,N
        m(i) = ro(i)*u(i)*Sp(i)
    end do
end subroutine calcula_fluxo_massa

!-------------------------------------------------

subroutine calcula_coeficiente_descarga
    call calcula_fluxo_massa
    do i=1,N
        Cd(i) = m(i)/Ma(i)
    end do
end subroutine calcula_coeficiente_descarga

!-------------------------------------------------

subroutine calcula_empuxo
    do i=1, N
        Empuxo(i) = ro(i)*u(i)*Sp(i)*u(i)
    end do
end subroutine calcula_empuxo

!-------------------------------------------------

subroutine correcoes_com_plinha

    real*8 :: fator ! auxiliar

    ! pressão
    p = p + pl

    ! massa específica nodal
    ro = ro + pl / ( Rgases * T )

    ! velocidade média
    do i = 1, n-1
       ue(i) = ue(i) - de(i) * ( pl(i+1) - pl(i) )
    end do

    ! velocidade na entrada da tubeira
    u_in = ue(1)

    ! velocidade nodal
    do i = 2,n-1
       u(i) = u(i) - 0.5d0 * ds(i) * ( pl(i+1) - pl(i-1) )
    end do

    fator = 2.0d0 * xp(2) / ( xp(3) - xp(2) )
	u(1) = u(2) - fator * ( u(3) - u(2) )

    fator = 2.0d0 * ( xp(n) - xp(n-1) ) / ( xp(n-1) - xp(n-2) )
	u(n) = u(n-1) + fator * ( u(n-1) - u(n-2) )
	
end subroutine correcoes_com_plinha
!-------------------------------------------------

end module coeficientes

