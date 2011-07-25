module variaveis

use portlib

implicit none

integer :: N           ! número total de nós

real*8  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! auxílio do comando System 

integer :: iteracao	   !	Número de iterações

real*8  :: Beta

real*8  :: fDarcy  ! Fator de atrito e velocidade inicial

real*8  :: Lt   ! Comprimento domínio de cálculo e do volume de controle

real*8  :: dt, Pi  ! Número de avanços no tempo

real*8,dimension(:),allocatable :: u      ! solução numérica
real*8,dimension(:),allocatable :: p      ! solução numérica
real*8,dimension(:),allocatable :: T ! solução numérica
real*8,dimension(:),allocatable :: ro      ! solução numérica

real*8,dimension(:),allocatable :: ue     ! solução numérica
real*8,dimension(:),allocatable :: roe      ! solução numérica
real*8,dimension(:),allocatable :: pl ! solução numérica

real*8,dimension(:),allocatable :: ue_o   ! solução numérica inicial
real*8,dimension(:),allocatable :: u_o    ! solução numérica inicial
real*8,dimension(:),allocatable :: ro_o    ! solução numérica inicial
real*8,dimension(:),allocatable :: p_o    ! solução numérica inicial
real*8,dimension(:),allocatable :: T_o    ! solução numérica inicial

real*8,dimension(:),allocatable :: Sp, Se  ! solução numérica da Área
real*8,dimension(:),allocatable :: M, Ma, Me  ! fluxo de massa na face leste fluxo de massa analitico
real*8,dimension(:),allocatable :: Empuxo, Mach, Mache, Ua, Cd, Ta , Pa, Roa

real*8,dimension(:),allocatable :: xp, xe  ! coordenada espacial nodal
real*8,dimension(:),allocatable :: Raio      ! raio do duto
real*8,dimension(:),allocatable :: bf, bc, f
character*255, dimension(:), allocatable ::files

real*8 :: R      ! constante dos gases
real*8 :: cp      ! cp
real*8 :: gama      ! gama
real*8 :: rin      ! rin
real*8 :: rg      ! rg
real*8 :: Lc      ! Lc
real*8 :: Ln      ! Ln

real*8 :: u_in  ! velocidade na entrada da tubeira (m/s)
real*8 :: p_in  ! pressão na entrada da tubeira que satisfaz QM (Pa)
real*8 :: T_in  ! temperatura na entrada da tubeira (K)
real*8 :: ro_in ! massa específica na entrada da tubeira (m/s)
real*8 :: p_cam   ! pressão na câmara de combustão (Pa)
real*8 :: T_cam   ! temperatura na câmara de combustão (K)
real*8 :: p_ia  ! pressão antiga na entrada da tubeira (Pa)
real*8 :: pl_in ! variação da pressão na entrada da tubeira (Pa)
real*8 :: p_ex  ! pressão na saída da tubeira (Pa)
real*8 :: T_ex  ! temperatura na saída da tubeira (K)
real*8 :: ro_ex ! massa específica na saída da tubeira (m/s)
real*8 :: P_out

real*8,dimension(:),allocatable :: ap, aw, ae, bp, aww, AAp ! coeficiente central de u e p
real*8,dimension(:),allocatable :: afu, atu, btu, bpru ! termos inclusos apu e bpu
real*8,dimension(:),allocatable :: ds, de   ! coeficientes do método SIMPLEC
character*255 ::richardson_path
character*255 ::richardson_exe
character*100 :: caso      ! nome do arquivo de saída
character*255 :: richardson_1      ! nome do arquivo principal do programa Richardson
character*255 :: richardson_2      ! nome do arquivo de entrada de dados para o programa Richardson
character*255 :: richardson_3      ! nome do arquivo de saída
character*255 :: richardson_4      ! nome do arquivo de entrada com a solução analítica exata
character*92 :: title     ! título do grafico
!graficos
integer :: graf_m, graf_t, graf_v, graf_ro, graf_p, graf_e, graf_cdesc, graf_dom, graf_mach
integer :: res_iter, res_result, res_coef

integer :: RazaoRef, Niveis
!calculo do residuo
real*8 ::Residuo_T, Residuo_T_o,Residuo_P, Residuo_P_o,Residuo_U, Residuo_U_o
real*8 ::p_outa, p_outN, pl_out, M_out, u_out, T_out

integer ::N_Sol, gerar_analitico, tipo
end