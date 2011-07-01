program Nozzle1DNormalShock
!
use dados

use resultados

use NormalShock1D
use CriaArqRichadson
! -----------------------------------------------

implicit none

!-------------------------------------------------
    integer ::i

    call le_dados
    
    
    do i=1, Niveis
      call init
      call ThroatFinder
        
      call solucao_numerica
      ! escrita da variável primária e sua visualização
      call gera_graficos
      call gera_txt
      
      call solucao_analitica
      call dealloc
      N = RazaoRef*N
    end do
    !ESCREVER OS COEFICIENTES E TERMO FONTE TODOS ok
    !propriedades DA GARGANTA
    !RICHARDSON DO ULTIMO
    !gradiente pra localizar o choque
   
    !call teste
! -----------------------------------------------

contains
subroutine teste

!real*8 ::mach1, mach2, gama, P1,P2,T1,T2, U1, U2, AA, mach
!integer ::ver
!    gama = 1.4d0
!!    mach1 = 3.0d0 
!    mach2 = 0.0d0
!    P1 = 0.5d0
    !T1 = 200.d0
!call PrandtRelation(mach1, mach2)
    !call Mach2Calc(gama, Mach1, Mach2)
!Mach1Calc
    !call P2Calc(gama, Mach1, p1, p2)
!RO2Calc
    !call T2Calc(gama, Mach1, T1, T2)
    !call U2Calc(gama, Mach1, U1, U2)
!P0Calc
!T0Calc
!Ro0Calc
    !Mach = 3.368d0
    
!call ACalc
!UCalc
    !ROOT(gama, P01, P02,m)
    !mach2 = 0.0d0
    !mach = ROOT(gama,1000.0d0,910.4277d0,mach2)
    !AA = 0.0d0
    !call AARatioCalc(gama, Mach, AA)
    !call Mach2Calc(gama,mach,mach2)
!    call WriteConfFile(richardson_1, richardson_2)
!    call CreateMeshFile(richardson_2,richardson_3,richardson_4,caso,5)
!    ver = system('notepad.exe ' // richardson_1)
!write(*,*), 'a'

end subroutine teste

end program

