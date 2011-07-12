program Nozzle1DNormalShock

use main
use output
use NormalShock1D
use CriaArqRichadson
! -----------------------------------------------

implicit none

!-------------------------------------------------
    integer ::i
    
    call conf_file_read
    
    if(gerar_analitico) call GenerateSolutionFile(numeroNos)
    
    do i=1, Niveis
      
      call init_alloc(N)
      call solucao_analitica_init(N)
      
      call solucao_numerica
      ! escrita da variável primária e sua visualização
      call gera_txt
      call gera_graficos
      
      
      call mostra_dados
      
      call teste
      
      call dealloc
      N = RazaoRef*N
    end do
    !ESCREVER OS COEFICIENTES E TERMO FONTE TODOS ok
    !propriedades DA GARGANTA
    !RICHARDSON DO ULTIMO
    !adicionar a entropia entalpia nos graficos
    !gradiente pra localizar o choque
    !call teste
! -----------------------------------------------

contains
subroutine teste

!    call WriteConfFile(richardson_1, richardson_2)
!    call CreateMeshFile(richardson_2,richardson_3,richardson_4,caso,5)
!    ver = system('notepad.exe ' // richardson_1)
!write(*,*), 'a'

end subroutine teste

end program

