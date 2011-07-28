module variaveis

use portlib

implicit none

integer :: N           ! número total de nós

real*16  :: tcpu        ! tempo de CPU em segundos

integer :: ver         ! auxílio do comando System 

integer :: iteracao	   !	Número de iterações

real*16  :: Beta

real*16  :: fDarcy  ! Fator de atrito e velocidade inicial

real*16  :: Lt   ! Comprimento domínio de cálculo e do volume de controle

real*16  :: dt, Pi  ! Número de avanços no tempo

real*16,dimension(:),allocatable :: u      ! solução numérica
real*16,dimension(:),allocatable :: p      ! solução numérica
real*16,dimension(:),allocatable :: T ! solução numérica
real*16,dimension(:),allocatable :: ro      ! solução numérica

real*16,dimension(:),allocatable :: ue     ! solução numérica
real*16,dimension(:),allocatable :: roe      ! solução numérica
real*16,dimension(:),allocatable :: pl ! solução numérica

real*16,dimension(:),allocatable :: ue_o   ! solução numérica inicial
real*16,dimension(:),allocatable :: u_o    ! solução numérica inicial
real*16,dimension(:),allocatable :: ro_o    ! solução numérica inicial
real*16,dimension(:),allocatable :: p_o    ! solução numérica inicial
real*16,dimension(:),allocatable :: T_o    ! solução numérica inicial

real*16,dimension(:),allocatable :: Sp, Se  ! solução numérica da Área
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
real*16 :: p_in  ! pressão na entrada da tubeira que satisfaz QM (Pa)
real*16 :: T_in  ! temperatura na entrada da tubeira (K)
real*16 :: ro_in ! massa específica na entrada da tubeira (m/s)
real*16 :: p_cam   ! pressão na câmara de combustão (Pa)
real*16 :: T_cam   ! temperatura na câmara de combustão (K)
real*16 :: p_ia  ! pressão antiga na entrada da tubeira (Pa)
real*16 :: pl_in ! variação da pressão na entrada da tubeira (Pa)
real*16 :: p_ex  ! pressão na saída da tubeira (Pa)
real*16 :: T_ex  ! temperatura na saída da tubeira (K)
real*16 :: ro_ex ! massa específica na saída da tubeira (m/s)
real*16 :: P_out

real*16,dimension(:),allocatable :: ap, aw, ae, bp, aww, AAp ! coeficiente central de u e p
real*16,dimension(:),allocatable :: afu, atu, btu, bpru ! termos inclusos apu e bpu
real*16,dimension(:),allocatable :: ds, de   ! coeficientes do método SIMPLEC
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
integer :: res_iter, res_result, res_coef,freq

integer :: RazaoRef, Niveis
!calculo do residuo
real*16 ::Residuo_T, Residuo_T_o,Residuo_P, Residuo_P_o,Residuo_U, Residuo_U_o
real*16 ::p_outa, p_outN, pl_out, M_out, u_out, T_out

integer ::N_Sol, gerar_analitico, tipo
end