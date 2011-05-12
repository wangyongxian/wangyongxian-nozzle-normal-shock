module variaveis

use portlib

implicit none

integer :: N           ! n�mero total de n�s

integer :: i,j           ! n�mero do n� e contador

real*8  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! aux�lio do comando System 

integer :: iteracao	   !	N�mero de itera��es

real*8  :: Beta

real*8  :: fator  ! Fator de atrito e velocidade inicial

real*8  :: Lt, dx   ! Comprimento dom�nio de c�lculo e do volume de controle

real*8  :: dt, Pi  ! N�mero de avan�os no tempo

real*8  :: Fat, R_o, R ! Coeficiente corre��o velocidade na face

real*8,dimension(:),allocatable :: u      ! solu��o num�rica
real*8,dimension(:),allocatable :: p      ! solu��o num�rica
real*8,dimension(:),allocatable :: T ! solu��o num�rica
real*8,dimension(:),allocatable :: ro      ! solu��o num�rica

real*8,dimension(:),allocatable :: ue     ! solu��o num�rica
real*8,dimension(:),allocatable :: roe      ! solu��o num�rica
real*8,dimension(:),allocatable :: pl ! solu��o num�rica

real*8,dimension(:),allocatable :: ue_o   ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: u_o    ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: ro_o    ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: p_o    ! solu��o num�rica inicial
real*8,dimension(:),allocatable :: T_o    ! solu��o num�rica inicial

real*8,dimension(:),allocatable :: Sp, Se  ! solu��o num�rica da �rea
real*8,dimension(:),allocatable :: M, Ma, Me  ! fluxo de massa na face leste fluxo de massa analitico
real*8,dimension(:),allocatable :: Empuxo, Mach, Mache, Ua, Cd, ropA, Ta 

real*8,dimension(:),allocatable :: xp, xe  ! coordenada espacial nodal
real*8,dimension(:),allocatable :: Raio      ! raio do duto
real*8,dimension(:),allocatable :: bf, bc

real*8 :: Rgases      ! constante dos gases
real*8 :: cp      ! cp
real*8 :: gama      ! gama
real*8 :: rin      ! rin
real*8 :: rg      ! rg
real*8 :: Lc      ! Lc
real*8 :: Ln      ! Ln

real*8 :: u_in  ! velocidade na entrada da tubeira (m/s)
real*8 :: p_in  ! press�o na entrada da tubeira que satisfaz QM (Pa)
real*8 :: T_in  ! temperatura na entrada da tubeira (K)
real*8 :: ro_in ! massa espec�fica na entrada da tubeira (m/s)
real*8 :: p_cam   ! press�o na c�mara de combust�o (Pa)
real*8 :: T_cam   ! temperatura na c�mara de combust�o (K)
real*8 :: p_ia  ! press�o antiga na entrada da tubeira (Pa)
real*8 :: pl_in ! varia��o da press�o na entrada da tubeira (Pa)
real*8 :: p_ex  ! press�o na sa�da da tubeira (Pa)
real*8 :: T_ex  ! temperatura na sa�da da tubeira (K)
real*8 :: ro_ex ! massa espec�fica na sa�da da tubeira (m/s)

real*8,dimension(:),allocatable :: ap, aw, ae, bp ! coeficiente central de u e p
!real*8,dimension(:),allocatable :: aWu, aWplinha, aWt ! coeficiente esquerdo de u e p
!real*8,dimension(:),allocatable :: aEu, aEplinha, aEt ! coeficiente direito de u	e p
									   
!real*8,dimension(:),allocatable :: bPu, bPplinha, bPt ! termo fonte de u e p

real*8,dimension(:),allocatable :: afu, atu, btu, bpru ! termos inclusos apu e bpu
real*8,dimension(:),allocatable :: ds, de   ! coeficientes do m�todo SIMPLEC

character*20 :: caso      ! nome do arquivo de sa�da

character*50 :: title     ! t�tulo do gr�fico
character*62 :: head      ! t�tulo do gr�fico + dia

character*12 :: dia       ! data da simula��o
character*8  :: hora      ! hor�rio da simula��o
integer*4    :: var(8)    ! data e hora
character*20 :: vardate   ! data e hora
character*20 :: vartime   ! data e hora
character*20 :: varzone   ! data e hora
character*70 :: note_caso ! notepad + caso
character*2  :: aux1,aux2
character*4  :: aux3
character*50 :: aux


!graficos
integer :: graf_m, graf_t, graf_v, graf_ro, graf_p, graf_e, graf_cdesc, graf_dom


end