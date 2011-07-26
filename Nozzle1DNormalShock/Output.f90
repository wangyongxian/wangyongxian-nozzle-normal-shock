module output

use variaveis
use solvers_1D
use NormalShock1D


implicit none

contains

!-------------------------------------------------
subroutine gera_txt
    integer ::i
    integer*4    :: var(8)    ! data e hora
    character*20 :: vardate   ! data e hora
    character*20 :: vartime   ! data e hora
    character*20 :: varzone   ! data e hora
    character*2  :: aux1,aux2
    character*4  :: aux3
    character*50 :: aux
    character*12 :: dia       ! data da simulação
    character*8  :: hora      ! horário da simulação
    character*62 :: head      ! título do gráfico + dia

    open(10, file="resultados.txt")
        
    call date_and_time(vardate,vartime,varzone,var)
    
    write(aux,*) var(5)
    aux1 = trim(adjustl(aux))
    write(aux,*) var(6)
    aux2 = trim(adjustl(aux))
    write(aux,*) var(7)
    aux3 = trim(adjustl(aux))
    hora = trim(aux1)//':'//trim(aux2)//':'//aux3
    
    head = trim(title)//" "//trim(dia)//"'"
    
    !write(10,2) trim(adjustl(title)), dia, hora
    !2 format(/,'Título = ', a<comp1>, &
    !        //,5x,'Dia = ',a12,5x,'Hora = ',a8)
    !comp = len(trim(adjustl(caso)))

    write(10,1)  P_cam, P_out, T_cam, N, dt, iteracao, Lt, fDarcy, R, Gama, rin, rg, Lc, Ln, Beta

    1 format(/,2x,'DADOS',//,  &
                ! a<comp>,  ' = caso',/, &	 
				 1pe16.8,  ' = pressao de estagnacao',/, &
				 1pe16.8,  ' = pressao da saida',/, &
				 1pe16.8,  ' = temperatura de estagnacao',/, &
				 8x,i8,    ' = número de volumes de controle',/, &
				 1pe16.8,  ' = número de avanços no tempo',/, &
				 8x,i8,    ' = número de iterações (ciclo total)',/, &
				 1pe16.8,  ' = comprimento do domínio de cálculo',/, &
				 1pe16.8,  ' = fator de atrito de Darcy',/, &
				 1pe16.8,  ' = R constante do gas',/,    &
				 1pe16.8,  ' = Gama',/,  & 
				 1pe16.8,  ' = Rin',/,  & 
				 1pe16.8,  ' = Rg',/,  & 
				 1pe16.8,  ' = Lc',/,  & 
				 1pe16.8,  ' = Ln',/,  & 
				 1pe16.8,  ' = Beta',/)

    9 format(i4,4x,4(1pe21.11))
	
	write(10,13)
    13 format(//,t1,'volume',t13,'U(Analitica)',t34,'U(Numerica)',t55,'Erro')

	do i = 1, N
	  write(10,9) i, ua(i), u(i), (ua(i) - u(i))
	end do
	
	!Temperarura
	write(10,14)
    14 format(//,t1,'volume',t13,'T(Analitica)',t34,'T(Numerica)',t55,'Erro')
	do i = 1, N
	    
	  write(10,9) i, Ta(i), T(i), (Ta(i) - T(i))
	end do
	
	!pressao
	write(10,15)
    15 format(//,t1,'volume',t13,'P(Analitica)',t34,'P(Numerica)',t55,'Erro')
	do i = 1, N
	  write(10,9) i, Pa(i), P(i), (Pa(i) - p(i))
	end do
	
	!Mach
	write(10,16)
    16 format(//,t1,'volume',t13,'Mach(Analitico)',t34,'Mach(Numerico)',t55,'Erro')
	do i = 1, N 
	  write(10,9) i, Mach(i), (u(i)/dsqrt(gama*R*T(i))), (Mach(i) - (u(i)/dsqrt(gama*R*T(i))))
	end do
	
	!massa especifica
	write(10,17)
    17 format(//,t1,'volume',t13,'RO(Analitico)',t34,'RO(Numerico)',t55,'Erro')
	do i = 1, N
	  write(10,9) i, roa(i), ro(i), (roa(i) - ro(i))
	end do

    write(10,10) tcpu
    10 format(/, f14.3, ' = tempo de processamento (segundos)')
    
	close(10)
    
end subroutine gera_txt

subroutine gera_graficos(analitica)
logical, intent(in) ::analitica
    integer ::i
    real*8 ::T0
    
    open(23, file='fm.dat')
    do i = 1, N
	  write(23,*) xp(i), Ma(i), M(i)
	end do
    close(23)
    
    open(23, file='u.dat')
    do i = 1, N
	  write(23,48) xp(i), Ua(i), U(i), raio(i)*10000.d0
	  48 format(4(1pe27.18))
	end do
	
    close(23)
    
    open(23, file='Mach.dat')
    do i = 1, N
	  write(23,48) xp(i), Mach(i), (u(i)/dsqrt(gama*R*T(i))), raio(i)
	end do
	
	close(23)
    open(23, file='ro.dat')
    do i = 1, N
	  write(23,48) xp(i), roa(i), ro(i), raio(i)
	end do
	
	open(23, file='p.dat')
    do i = 1, N
	  write(23,48) xp(i), pa(i), p(i), raio(i)
	end do
	
    close(23)
    open(22,file='dominio.dat')
	do i = 1, N
	  write(22,*) xp(i), raio(i)
	end do
	close(22)
	
    open(22,file='T.dat')
	do i = 1, N
	    call T0Calc(gama,T(i), (u(i)/dsqrt(gama*R*T(i))), T0)
	    !call T0Calc(gama,Ta(i), Mach(i), T0)
	  write(22,48) xp(i), Ta(i), T(i), raio(i)
	end do
	close(22)
	
	call create_gnufile('U.gnu', '''u.dat''', '''solucao_analitica.txt''', 3, analitica, '''x(m)''', '''u(m/s)''', '''teste titulo''')
    call create_gnufile('T.gnu', '''T.dat''', '''solucao_analitica.txt''', 5, analitica, '''x(m)''', '''T(K)''', '''teste titulo''')
    call create_gnufile('RO.gnu', '''RO.dat''', '''solucao_analitica.txt''', 7, analitica, '''x(m)''', '''Ro(km/m^3)''', '''teste titulo''')
    call create_gnufile('P.gnu', '''P.dat''', '''solucao_analitica.txt''', 6, analitica, '''x(m)''', '''P(Pa)''', '''teste titulo''')
    call create_gnufile('Mach.gnu', '''Mach.dat''', '''solucao_analitica.txt''', 8, analitica, '''x(m)''', '''Mach''', '''teste titulo''')
    call create_gnufile('fm.gnu', '''fm.dat''', '''solucao_analitica.txt''',  9    , analitica, '''x(m)''', '''Fluxo de Massa (kg/s)''', '''teste titulo''')
    call create_gnufile('coef_descarga.gnu', '''fm.dat''', '''solucao_analitica.txt''',  10    , analitica, '''x(m)''', '''Coeficiente de Descarga (kg/s)''', '''teste titulo''')
    call create_gnufile('empuxo.gnu', '''fm.dat''', '''solucao_analitica.txt''',  11    , analitica, '''x(m)''', '''Empuxo (kg/s)''', '''teste titulo''')
    
    
end subroutine gera_graficos

subroutine mostra_dados

    if (graf_p) then
        ver = system('wgnuplot p.gnu')
    end if
    if (graf_cdesc) then
        ver = system('wgnuplot coef_descarga.gnu')
    end if
    if (graf_e) then
        ver = system('wgnuplot empuxo.gnu')
    end if
    if (graf_dom) then
        ver = system('wgnuplot dominio.gnu')
    end if
    if (graf_ro) then
        ver = system('wgnuplot ro.gnu')
    end if
    if (graf_m) then
        ver = system('wgnuplot fm.gnu')
    end if
    if (graf_v) then
        ver = system('wgnuplot U.gnu')
    end if
    if (graf_t) then
       ver = system('wgnuplot T.gnu')
    end if
    if(graf_mach)then
       ver = system('wgnuplot Mach.gnu')
    end if
    if(res_coef) then
    ver = system('coeficientes.txt')
    end if 
    if (res_result) then
        ver = system('resultados.txt')
    end if
    
    if (res_iter) then
        ver = system('norma_l1.txt')
    end if
    
end subroutine mostra_dados

subroutine gera_arq_coef(prop)
    integer, intent(in) ::prop
    !prop 2 - PL
    !prop 3 - Energia
    !prop 4 - QML
    integer ::i
    open(14,file='coeficientes.txt')
    2 format('Coeficientes de Pl', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    3 format('Coeficientes Energia', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    4 format('Coeficientes QML', /,t4, 'volume', t13, 'x', t34,'oeste',t55,'central',t76,'leste',t97,'fonte',/)
    5 format(i4, 4x, 5(1pe21.11))
    select case(prop)
        case(2)
            write(14, 2)
        case(3)
            write(14, 3)
        case(4)
            write(14, 4)
    end select
    
    do i=1,N
     !   write(14, 5) i, xp(i), aw(i), ap(i), ae(i), bp(i)
    end do
    
    close(14)
end subroutine gera_arq_coef

subroutine create_gnufile(filename, dat_file, filename_analictic, column, analitica, xlabel, ylabel, title)
character(*), intent(in) ::filename
character(*), intent(in) ::dat_file
character(*), intent(in) ::filename_analictic
logical, intent(in) ::analitica
character(*), intent(in) ::title
character(*), intent(in) ::xlabel
character(*), intent(in) ::ylabel
integer, intent(in) ::column
open(12, file=filename) 
if (analitica) then
12 format ('set data style linespoints',1/,                  &
'set grid',1/,                                               &
'plot ',a50,' using 1:', i2,' title ''analítica 1''', 1/,    &
'replot ',a50,' using 1:2 title ''analítica 2''',1/,         &
'replot ',a50,' using 1:3 title ''numérica''',1/,            &
'replot ',a50,' using 1:4 title ''geometria''',1/,           &
'set xlabel ',a10,1/,                                        &
'set ylabel ',a10,1/,                                        &
'set title ',a50,1/,                                         &
'replot')
write(12,12) filename_analictic, column, dat_file, dat_file, dat_file, xlabel, ylabel, title
else
13 format ('set data style linespoints',1/,                  &
'set grid',1/,                                               &
'plot ',a50,' using 1:2 title ''analítica 2''',1/,         &
'replot ',a50,' using 1:3 title ''numérica''',1/,            &
'replot ',a50,' using 1:4 title ''geometria''',1/,           &
'set xlabel ',a10,1/,                                        &
'set ylabel ',a10,1/,                                        &
'set title ',a50,1/,                                         &
'replot')
write(12,13) dat_file, dat_file, dat_file, xlabel, ylabel, title
end if


close(12)

end subroutine create_gnufile

end module output
