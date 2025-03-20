module mod_param
  !
  implicit none
  !
  integer, parameter :: ndims = 2
  real(8), parameter :: pi = acos(-1.)
  integer, dimension(ndims), parameter :: dims = (/1,8/)
  !integer, dimension(ndims), parameter :: dims = (/1,16/)
  !integer, dimension(ndims), parameter :: dims = (/1,32/)
  !integer, dimension(ndims), parameter :: dims = (/1,64/)
  !integer, parameter :: itot = 4, jtot = 128, ktot = 64
  integer, parameter :: itot = 4, jtot = 256, ktot = 128
  !integer, parameter :: itot = 4, jtot = 512, ktot = 256
  !integer, parameter :: itot = 4, jtot = 512, ktot = 512
  !integer, parameter :: itot = 4, jtot = 1024, ktot = 512
  !integer, parameter :: itot = 4, jtot = 1024, ktot = 1024
  integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: n  = (/imax,jmax,kmax/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
  real(8), parameter :: lz = 0.0225d0
  real(8), parameter :: ly = lz*2.d0
  !real(8), parameter :: lx = 0.0001d0
  real(8), parameter :: lx = itot*ly/jtot
  !
  real(8), parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
  !real(8), parameter :: dxi = jtot/ly, dyi = jtot/ly, dzi = ktot/lz
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  real(8), parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
  real(8), parameter, dimension(3) :: dl  = (/dx ,dy ,dz /)
  !
  include 'bc.h90'
  !
  ! simulations checks
  !
  !integer, parameter :: nstep       = 2e+9
  !integer, parameter :: nstep       = 25000
  integer, parameter :: nstep       = 10000
  !integer, parameter :: nstep       = 1
  !
  !integer, parameter :: istep_do_av = 1
  integer, parameter :: istep_do_av = 1e+7
  integer, parameter :: isave       = 50000
  !logical, parameter :: restart     = .true.
  logical, parameter :: restart     = .false.
  logical, parameter :: restart_fir = .true.
  !logical, parameter :: var_dt      = .false.
  logical, parameter :: var_dt      = .true.
  !
  ! initial flow condition
  !
  character(len=3), parameter :: inivel = 'zer' 
  character(len=3), parameter :: inivof = 'khi' 
  character(len=3), parameter :: initmp = 'uni'
  character(len=3), parameter :: inisca = 'std'
  !
  ! dimensionless physical parameters
  !
  real(8), parameter :: ra  = 1.00e+06
  real(8), parameter :: pr  = 1.00d0
  real(8), parameter :: we  = 500.00d0
  real(8), parameter :: eps = epsb           ! fr = 2.0*eps
  real(8), parameter :: sc  = 1.00d0         ! eva
  real(8), parameter :: ste = 0.20d0         ! eva
  real(8), parameter :: lambda_rho0 = 20.0d0
  real(8), parameter :: lambda_mu0  = 1.00d0
  real(8), parameter :: lambda_cp0  = 1.00d0
  real(8), parameter :: lambda_ka0  = 1.00d0
  real(8), parameter :: lambda_mm0  = 0.62d0 ! eva
  real(8), parameter :: lambda_dcp0 = 1.00d0 ! eva
  real(8), parameter :: ratio       = ly/lz  
  real(8), parameter :: deltaT_d    = 1.00d0 ! dimensionless
  real(8), parameter :: phi_p1      = 0.20d0 ! eva: (ru/(mm2*cp2))
  real(8), parameter :: phi_p2      = 1.00d0 ! eva: (pth0/((ru/mm2)*tg0))
  real(8), parameter :: phi_p3      = 0.10d0 ! beta_l_th*deltaT
#ifndef LOW_MACH
  real(8), parameter :: phi_p4      = 0.10d0 ! beta_g_th*deltaT ! only OBB (re-check)
#endif
  real(8), parameter :: phi_p5      = 0.2636d0 ! eva (pth_0/pc) 
  real(8), parameter :: phi_p6      = 0.5961d0 ! eva (tg0/tc) 
  !
  ! reference quantities
  !
  real(8), parameter :: tc     = 647.09600d0 
  real(8), parameter :: pc     = 220.64000*101325 
  real(8), parameter :: t_bot  = tb_b
  real(8), parameter :: t_top  = tb_t
  real(8), parameter :: tl0    = 0.5d0*(tb_b+tb_t)
  real(8), parameter :: tg0    = 0.5d0*(tb_b+tb_t)
  real(8), parameter :: pth_0  = 5.8925e+06
  real(8), parameter :: ru     = 8314.0d0
  real(8), parameter :: m2     = 28.013d0
  real(8), parameter :: rho2_0 = pth_0/(phi_p2*(ru/m2)*tg0)
  real(8), parameter, dimension(3) :: gacc = (/0.0,0.0,-9.81d0/)
  !
  real(8), parameter :: lref   = lz
  real(8), parameter :: deltaT = (t_bot-t_top)*deltaT_d
  real(8), parameter :: uref   = sqrt(abs(gacc(3))*lref*2.d0*eps)
  real(8), parameter :: tmp0   = tl0
  real(8), parameter :: t_ff   = lref/uref
  !
  ! output frequency (so we know t_ff)
  !
  real(8), parameter :: frac = 0.5d0
  !
  integer, parameter :: icheck  = 10
  integer, parameter :: iout0d  = 100
  integer, parameter :: iout1d  = 10000 ! print 
  integer, parameter :: iout_av = 100   ! averaging frequency (not necesseraly the printing one)
  integer, parameter :: iout2d  = 5000
  integer, parameter :: iout3d  = 100
  !
  !integer, parameter :: icheck  = 10
  !integer, parameter :: iout0d  = 50*4!*2
  !integer, parameter :: iout1d  = int(frac*t_ff/dt) ! print
  !integer, parameter :: iout_av = 100   ! averaging frequency (not necesseraly the printing one)
  !integer, parameter :: iout2d  = 1000*20!*5
  !integer, parameter :: iout3d  = 5000*40!*5
  !
  ! dimensional quantities
  !
  real(8), parameter :: m1        = lambda_mm0*m2
  real(8), parameter :: cp2       = ru/(phi_p1*m2)
  real(8), parameter :: nu2talp2  = 2.d0*eps*abs(gacc(3))*(lref**3)/ra
  real(8), parameter :: nu2       = sqrt(nu2talp2*pr)
  real(8), parameter :: mu2       = nu2*rho2_0
  real(8), parameter :: alpha2    = nu2/pr
  real(8), parameter :: rho1      = lambda_rho0*rho2_0
  real(8), parameter :: sigmaca   = rho1*(uref**2.d0)*lref/we ! we use the density of the liquid like Lohse-Verzicco
  real(8), parameter :: mu1       = lambda_mu0*mu2
  real(8), parameter :: kappa2    = alpha2*rho2_0*cp2
  real(8), parameter :: kappa1    = lambda_ka0*kappa2
  real(8), parameter :: cp1       = lambda_cp0*cp2
  real(8), parameter :: beta_l_th = phi_p3/deltaT
#ifndef LOW_MACH
  real(8), parameter :: beta_g_th = phi_p4/deltaT ! check
#endif
  real(8), parameter :: d_m12     = (mu2/rho2_0)*(1.d0/sc)
  real(8), parameter :: lheat     = cp2*deltaT/ste ! be sure to put a good Ste
  real(8), parameter :: delta_cp  = lambda_dcp0*cp2
  !
  integer, parameter :: e         = 0
  real(8), parameter :: sinit     = 0.d0
  real(8), parameter :: gri       = 0.d0
  real(8), parameter :: cfl       = 0.25d0
  real(8), parameter :: cfl_o     = 0.50d0
  real(8), parameter :: small     = 1e-09
  real(8), parameter :: theta_thr = 0.25d0
  real(8), parameter :: mflux0    = 0.0d0
  real(8), parameter :: dt_input  = 2.10E-005
  !real(8), parameter :: dt_input  = (2.50e-4)*t_ff!*2.d0
  !real(8), parameter :: dt_input  = (7.50e-4)*t_ff
  !
  ! disperse phase 
  !
  integer, parameter :: n_dp  = 1 ! initial number of droplets
  real(8), parameter :: alpha = acos(lz/sqrt(lx**2.d0+lz**2.d0))
  real(8), parameter :: xcc(1:n_dp) = (/0.500d0/)*lx
  real(8), parameter :: ycc(1:n_dp) = (/0.400d0/)*ly
  real(8), parameter :: zcc(1:n_dp) = (/0.200d0/)*lz
  real(8), parameter :: rd(1:n_dp)  = 0.25d0*lz
  !
  character(len=*), parameter :: datadir = 'data/'
  character(len=*), parameter :: datapos = 'data/post/'
  character(len=*), parameter :: datapos_mass = 'data/post/mass/'
  character(len=*), parameter :: datapos_vol  = 'data/post/vol/'
  !
  logical, parameter, dimension(2,3) :: no_outflow = & 
      reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                .false.,.false.,   & ! no outflow in y lower,upper bound
                .false.,.false./), & ! no outflow in z lower,upper bound
                 shape(no_outflow))
  !
  ! re-calculation of the dimensionless physical parameters
  !
  real(8), parameter :: ra_c          = 2.d0*eps*abs(gacc(3))*(lref**3)/nu2talp2 
  real(8), parameter :: pr_c          = nu2/alpha2
  real(8), parameter :: we_c          = rho1*(uref**2.d0)*lref/sigmaca
  real(8), parameter :: eps_c         = epsb
  real(8), parameter :: sc_c          = nu2/d_m12                  ! only eva
  real(8), parameter :: ste_c         = cp2*tg0/lheat              ! only eva
  real(8), parameter :: lambda_rho0_c = rho1/rho2_0                
  real(8), parameter :: lambda_mu0_c  = mu1/mu2                    
  real(8), parameter :: lambda_cp0_c  = cp1/cp2                    
  real(8), parameter :: lambda_ka0_c  = kappa1/kappa2              
  real(8), parameter :: lambda_mm0_c  = m1/m2                      ! only eva
  real(8), parameter :: lambda_dcp0_c = delta_cp/cp2               ! only eva
  real(8), parameter :: ratio_c       = ly/lz                      
  real(8), parameter :: deltaT_d_c    = deltaT/(t_bot-t_top)        
  real(8), parameter :: phi_p1_c      = (ru/(m2*cp2))              ! only eva
  real(8), parameter :: phi_p2_c      = (pth_0/(rho2_0*(ru/m2)*tg0))       
  real(8), parameter :: phi_p3_c      = beta_l_th*deltaT
#ifndef LOW_MACH        
  real(8), parameter :: phi_p4_c      = beta_g_th*deltaT          
#endif
  real(8), parameter :: phi_p5_c      = pth_0/pc ! only eva
  real(8), parameter :: phi_p6_c      = tg0/tc   ! only eva
  !
end module mod_param
