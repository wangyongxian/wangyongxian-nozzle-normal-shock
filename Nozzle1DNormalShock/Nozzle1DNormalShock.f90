program Nozzle1DNormalShock

use main
use output
use NormalShock1D
use CriaArqRichadson
! -----------------------------------------------

implicit none

!-------------------------------------------------
    integer ::i
    integer ::local
    
    
    call conf_file_read
    
    
    if(gerar_analitico) call GenerateSolutionFile(numeroNos)
    
    
    call WriteConf1File(richardson_1, richardson_2)
    call WriteConf2File(richardson_2,richardson_3,richardson_4,caso,Niveis)
    
    
    do i=1, Niveis
      
      call init_alloc(N)
      local = ThroatFinder()
      call solucao_analitica_init(N)
      
      !para o richardson
      if(i == Niveis) call WriteAnalitico(richardson_4,local)
      
      
      call solucao_numerica
      
      ! escrita da variável primária e sua visualização
      call gera_txt
      call gera_graficos
      call mostra_dados
      
      !call teste(local)
      
      call WriteMeshFile(files(i), local)
      
      call dealloc
      N = RazaoRef*N
    end do
    !ESCREVER OS COEFICIENTES E TERMO FONTE TODOS ok
    !propriedades DA GARGANTA
    !RICHARDSON DO ULTIMO
    !adicionar a entropia entalpia nos graficos
    !gradiente pra localizar o choque
    !call teste
    deallocate(files)
! -----------------------------------------------

contains
subroutine teste
    call WriteConf1File(richardson_1, richardson_2)
    call WriteConf2File(richardson_2,richardson_3,richardson_4,caso,Niveis)
    ver = system('notepad.exe ' // richardson_1)
    ver = system('notepad.exe ' // richardson_2)
    !ver = system('notepad.exe ' // richardson_3)
!write(*,*), 'a'

end subroutine teste

end program

