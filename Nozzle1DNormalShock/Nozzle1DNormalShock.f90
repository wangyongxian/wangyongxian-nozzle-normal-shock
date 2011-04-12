program Nozzle1DNormalShock

! -----------------------------------------------

use dados

use resultados

! -----------------------------------------------

implicit none

integer ::  comp1

!-------------------------------------------------

  call date_and_time(vardate,vartime,varzone,var)

  write(aux,*) var(3)
  aux1 = trim(adjustl(aux))
  write(aux,*) var(2)
  aux2 = trim(adjustl(aux))
  write(aux,*) var(1)
  aux3 = trim(adjustl(aux))
  dia = '('//trim(aux1)//'/'//trim(aux2)//'/'//aux3//')'

  write(aux,*) var(5)
  aux1 = trim(adjustl(aux))
  write(aux,*) var(6)
  aux2 = trim(adjustl(aux))
  write(aux,*) var(7)
  aux3 = trim(adjustl(aux))
  hora = trim(aux1)//':'//trim(aux2)//':'//aux3

  call le_dados

  head = trim(title)//" "//trim(dia)//"'"

  open(10,file=caso)

!  comp = len(trim(adjustl(nome)))

  comp1 = len(trim(adjustl(title)))

 ! write(10,18) trim(adjustl(nome)), trim(adjustl(title)), dia, hora
  !18 format(/, 'Aluno  = ', a<comp>,  &
	!        //,'Título = ', a<comp1>, &
    !        //,5x,'Dia = ',a12,5x,'Hora = ',a8)

  call mostra_dados 
  ! calcula a area e outras inicializacoes
  call inicializacao

  call solucao_numerica

  close (10)

  note_caso = 'notepad '//caso
  ver = system(note_caso) ! lista arquivo de resultados

! -----------------------------------------------

end
