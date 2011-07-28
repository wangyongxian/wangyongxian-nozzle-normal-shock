program Nozzle1DNormalShock

use main
use output
use NormalShock1D
use CriaArqRichadson
! -----------------------------------------------
use Solucao_CDS
use Solucao_CDS_UDS
use Solucao_UDS_UDS2
use Solucao_TVD
! -----------------------------------------------

implicit none

!-------------------------------------------------
    integer ::i
    integer ::local, shock
    
    call conf_file_read
    
    
    if(gerar_analitico) call GenerateSolutionFile(N_Sol)
    
    call WriteConf1File(richardson_1, richardson_2)
    call WriteConf2File(richardson_2,richardson_3,richardson_4,caso,Niveis)
    
    do i=1, Niveis
      
      call init_alloc(N)
      local = ThroatFinder()
      call solucao_analitica_init(N)
      
      !para o richardson
      if(i == Niveis)   call WriteAnalitico(richardson_4,local,ShockDistAnalitical)
      
      select case(tipo)
      case (1)
        call solucao_numerica(coeficientes_e_fontes_qml_cds, &
        calculo_velocidades_face_cds,                        &
        coeficientes_e_fontes_energia_cds,                   &
        coeficientes_fontes_massa_cds)
      case (2)
      call solucao_numerica(coeficientes_e_fontes_qml_cds_uds, &
        calculo_velocidades_face_cds_uds,                        &
        coeficientes_e_fontes_energia_cds_uds,                   &
        coeficientes_fontes_massa_cds_uds)
        !call solucao_numerica(coeficientes_e_fontes_qml_cds_uds, &
        !calculo_velocidades_face_cds_uds,                        &
        !coeficientes_e_fontes_energia_cds_uds,                   &
        !coeficientes_fontes_massa_cds_uds)
      case (3)
        call solucao_numerica(coeficientes_e_fontes_qml_tvd,    &
        calculo_velocidades_face_tvd,                           &
        coeficientes_e_fontes_energia_tvd,                      &
        coeficientes_fontes_massa_tvd)
      case  default
        call solucao_numerica(coeficientes_e_fontes_qml_cds,    &
        calculo_velocidades_face_cds,                           &
        coeficientes_e_fontes_energia_cds,                      &
        coeficientes_fontes_massa_cds)      
      end select
      
      call ShockNumFinder(p,shock)
      
      ! escrita da variável primária e sua visualização
      call gera_txt
      call gera_graficos(.false.)
      call mostra_dados
      
      !call teste(local)
      
      call WriteMeshFile(files(i), local,shock, Lt/(N-2))
      
      call dealloc
      N = RazaoRef*(N-2)+2
      dt = dt/RazaoRef
      iteracao = 2*iteracao
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
    !call WriteConf1File(richardson_1, richardson_2)
    !call WriteConf2File(richardson_2,richardson_3,richardson_4,caso,Niveis)
    !ver = system('notepad.exe ' // richardson_1)
    !ver = system('notepad.exe ' // richardson_2)
    !ver = system('notepad.exe ' // richardson_3)
!write(*,*), 'a'

end subroutine teste

end program
