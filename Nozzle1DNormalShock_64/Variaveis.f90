module variaveis

use portlib

implicit none

integer :: N           ! n�mero total de n�s

real*16  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! aux�lio do comando System 

integer :: iteracao	   !	N�mero de itera��es

real*16  :: Beta

real*16  :: fDarcy  ! Fator de atrito e velocidade inicial

real*16  :: Lt   ! Comprimento dom�nio de c�lculo e do volume de controle

real*16  :: dt, Pi  ! N�mero de avan�os no tempo

real*16,dimension(:),allocatable :: u      ! solu��o num�rica
real*16,dimension(:),allocatable :: p      ! solu��o num�rica
real*16,dimension(:),allocatable :: T ! solu��o num�rica
real*16,dimension(:),allocatable :: ro      ! solu��o num�rica

real*16,dimension(:),allocatable :: ue     ! solu��o num�rica
real*16,dimension(:),allocatable :: roe      ! solu��o num�rica
real*16,dimension(:),allocatable :: pl ! solu��o num�rica

real*16,dimension(:),allocatable :: ue_o   ! solu��o num�rica inicial
real*16,dimension(:),allocatable :: u_o    ! solu��o num�rica inicial
real*16,dimension(:),allocatable :: ro_o    ! solu��o num�rica inicial
real*16,dimension(:),allocatable :: p_o    ! solu��o num�rica inicial
real*16,dimension(:),allocatable :: T_o    ! solu��o num�rica inicial

real*16,dimension(:),allocatable :: Sp, Se  ! solu��o num�rica da �rea
real*16,dimension(:),allocatable :: M, Ma, Me  ! fluxo de massa na face leste fluxo de massa analitico
real*16,dimension(:),allocatable :: Empuxo, Mach, Mache, Ua, Cd, Ta , Pa, Roa

real*16,dimension(:),allocatable :: xp, xe  ! coordenada espacial nodal
real*16,dimension(:),allocatable :: Raio      ! raio do duto
real*16,dimension(:),allocatable :: bf, bc, f
character*255, dimension(:), allocatable ::files

real*16 :: R      ! constante dos gases
real*16 :: cp      ! cp
real*16 :: gama      ! gama
real*16 :: rin      ! rin
real*16 :: rg      ! rg
real*16 :: Lc      ! Lc
real*16 :: Ln      ! Ln
real*16 :: ShockDistAnalitical
real*16 :: u_in  ! velocidade na entrada da tubeira (m/s)
real*16 :: p_in  ! press�o na entrada da tubeira que satisfaz QM (Pa)
real*16 :: T_in  ! temperatura na entrada da tubeira (K)
real*16 :: ro_in ! massa espec�fica na entrada da tubeira (m/s)
real*16 :: p_cam   ! press�o na c�mara de combust�o (Pa)
real*16 :: T_cam   ! temperatura na c�mara de combust�o (K)
real*16 :: p_ia  ! press�o antiga na entrada da tubeira (Pa)
real*16 :: pl_in ! varia��o da press�o na entrada da tubeira (Pa)
real*16 :: p_ex  ! press�o na sa�da da tubeira (Pa)
real*16 :: T_ex  ! temperatura na sa�da da tubeira (K)
real*16 :: ro_ex ! massa espec�fica na sa�da da tubeira (m/s)
real*16 :: P_out

real*16,dimension(:),allocatable :: ap, aw, ae, bp, aww, AAp ! coeficiente central de u e p
real*16,dimension(:),allocatable :: afu, atu, btu, bpru ! termos inclusos apu e bpu
real*16,dimension(:),allocatable :: ds, de   ! coeficientes do m�todo SIMPLEC
character*255 ::richardson_path
character*255 ::richardson_exe
character*100 :: caso      ! nome do arquivo de sa�da
character*255 :: richardson_1      ! nome do arquivo principal do programa Richardson
character*255 :: richardson_2      ! nome do arquivo de entrada de dados para o programa Richardson
character*255 :: richardson_3      ! nome do arquivo de sa�da
character*255 :: richardson_4      ! nome do arquivo de entrada com a solu��o anal�tica exata
character*92 :: title     ! t�tulo do grafico
!graficos
integer :: graf_m, graf_t, graf_v, graf_ro, graf_p, graf_e, graf_cdesc, graf_dom, graf_mach
integer :: res_iter, res_result, res_coef,freq

integer :: RazaoRef, Niveis
!calculo do residuo
real*16 ::Residuo_T, Residuo_T_o,Residuo_P, Residuo_P_o,Residuo_U, Residuo_U_o
real*16 ::p_outa, p_outN, pl_out, M_out, u_out, T_out

integer ::N_Sol, gerar_analitico, tipo
end