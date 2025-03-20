module mod_param
  !
  implicit none
  !
  include 'bc.h90'
  !
  ! parameters for the grid
  !
  integer, parameter :: ndims = 2
  real(8), parameter :: pi = acos(-1.)
  integer, dimension(ndims), parameter :: dims = (/1,16/)
  integer, parameter :: itot = 4, jtot = 256, ktot = 256
  integer, parameter :: imax = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: n  = (/imax,jmax,kmax/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter, dimension(3) :: ng = (/itot,jtot,ktot/)! = itot/dims(1), jmax = jtot/dims(2),kmax = ktot
  integer, parameter :: i1 = imax+1, j1 = jmax+1, k1 = kmax+1
  !real(8), parameter :: lx = 0.0016d0
  !real(8), parameter :: ly = 0.0016d0
  !real(8), parameter :: lz = 0.0016d0
  real(8), parameter :: lx = 0.00001562
  real(8), parameter :: ly = 0.001d0
  real(8), parameter :: lz = 0.001d0
  !
  real(8), parameter :: dxi = itot/lx, dyi = jtot/ly, dzi = ktot/lz
  real(8), parameter, dimension(3) :: dli = (/dxi,dyi,dzi/)
  real(8), parameter :: dx = 1./dxi, dy = 1./dyi, dz = 1./dzi
  real(8), parameter, dimension(3) :: dl  = (/dx ,dy ,dz /)
  !
  ! simulations checks
  !
  integer, parameter :: nstep = 7000000
  !
  integer, parameter :: istep_do_av = 1e+7
  integer, parameter :: isave       = 10000
  !
  !logical, parameter :: start_f_ze      = .true.  ! start from 0
  logical, parameter :: start_f_ze      = .false.  ! start from 0
  logical, parameter :: restart_sp      = .false.  ! re-start from a pre-initialized single phase field
  !logical, parameter :: restart_sp      = .true.  ! re-start from a pre-initialized single phase field
  !
  logical, parameter :: restart_mp_isot = .false.  ! re-start from a pre-initialized isothermal multiphase field (no ht, no mt)
  !logical, parameter :: restart_mp_isot = .true.  ! re-start from a pre-initialized isothermal multiphase field (no ht, no mt)
  
  !logical, parameter :: restart_mp_comp = .false.  ! re-start from the full multiphase (with heat and mass) transfer
  logical, parameter :: restart_mp_comp = .true.   ! re-start from the full multiphase (with heat and mass) transfer
  !
  logical, parameter :: var_dt      = .true.
  real(8), parameter :: dt_input    = 2.10E-005
  !
  ! dimensionless physical parameters
  !
  real(8), parameter :: re  = 5750.d0          !checked
  real(8), parameter :: pr  = 0.71d0           !checked
  real(8), parameter :: we  = 0.10d0           !test
  real(8), parameter :: sc  = 0.87d0           !checked
  real(8), parameter :: ste = 0.748d0          !scaled by gas temperature,checked
  real(8), parameter :: lambda_rho0 = 10.d0  ! 575/18
  real(8), parameter :: lambda_mu0  = 63.68d0 
  real(8), parameter :: lambda_cp0  = 4.161d0  
  real(8), parameter :: lambda_ka0  = 0.2307!23.d0  
  real(8), parameter :: lambda_mm0  = 0.62d0  ! eva
  real(8), parameter :: lambda_dcp0 = 1.00d0   ! eva
  real(8), parameter :: phi_p1      = 0.20d0   ! eva (ru/(mm2*cp2))
  real(8), parameter :: phi_p2      = 1.00d0   ! eva (pth0/((ru/mm2)*tmp0))
  real(8), parameter :: phi_p3      = 1.29d0   ! eva (pth_0/pc) 
  real(8), parameter :: phi_p4      = 1.50d0   ! eva (tmp0/tc) 
  !real(8), parameter :: phi_p4      = 1.00d0  ! eva (tmp0/tc) 
  !real(8), parameter :: phi_p4      = 0.75d0  ! eva (tmp0/tc) 
  !
  ! referencie
  !
  real(8), parameter :: k0_wave = 2.d0*pi/lx
  integer, parameter :: k0_freq = 2
  real(8), parameter, dimension(3) :: gacc = (/0.d0,0.d0,0.d0/)
  real(8), parameter :: tc = 647.d0, pc = 220*101325
  real(8), parameter :: ru = 8314.d0
  real(8), parameter :: m2 = 28.9647d0
  !real(8), parameter :: res    = 50.d0
  real(8), parameter :: res    = 64.d0
  real(8), parameter :: lref   = res*dy

  real(8), parameter :: nu2    = 1.8586E-006
  !real(8), parameter :: nu2   = 1.0/re
  real(8), parameter :: d_m12  = 2.23e-5!nu2/sc
  !real(8), parameter :: uref  = d_m12/lref
  real(8), parameter :: uref   = re*k0_wave*nu2
  real(8), parameter :: u_abc  = uref
  !real(8), parameter :: pth_0  = 40*101325!saturation tmp is 352k
  real(8), parameter :: pth_0  = 1.d0*101325!saturation tmp is 352k
  real(8), parameter :: tmp0   = 298.d0!phi_p4*tc
  real(8), parameter :: rho2_0 = 1.d0!pth_0/(phi_p2*(ru/m2)*tmp0)!18.826
  real(8), parameter :: beta_l_th = 1.d0/tmp0
  real(8), parameter :: beta_g_th = 1.d0/tmp0
  !
  real(8), parameter :: tmpg0 = tmp0
  real(8), parameter :: tmpl0 = 298.d0  ! phi_p4 = 1.50
  !tmpl should always be less than the satuartion temperature at the Pth to avoid boiling 
  !at 40 bar the saturation temperature is 353K 
  ! we set the inital liquid temperature to the wet bulb temp= 315


  !real(8), parameter :: tmpl0  = 387.60d0  ! phi_p4 = 1.00
  !real(8), parameter :: tmpl0  = 333.51d0  ! phi_p4 = 0.75
  !
  ! initial flow condition
  !
  character(len=3), parameter :: inivel = 'zer' 
  character(len=3), parameter :: inivof = 'mub' 
  character(len=3), parameter :: initmp = 'sin'
  character(len=3), parameter :: inisca = 'std'
  !
  ! output frequency (so we know t_ff)
  !
  integer, parameter :: icheck    = 10
  !integer, parameter :: iout0d   = 50
  integer, parameter :: iout0d    = 50
  integer, parameter :: iout0d_ta = 5000
  integer, parameter :: iout_av   = 1000    ! averaging frequency (not necesseraly the printing one)
  integer, parameter :: iout1d    = 1000*20 ! print 
  integer, parameter :: iout2d    = 2000!10!1000*20
  !integer, parameter :: iout3d   = 5000*40
  integer, parameter :: iout3d    = 1000
  !integer, parameter :: iout3d   = 500000
  !integer, parameter :: iout3d   = 1
  !integer, parameter :: iout3d   = 20
  !integer, parameter :: iout3d   = 10
  !integer, parameter :: icheck   = 10
  !integer, parameter :: iout0d   = 50*4!*2
  !integer, parameter :: iout1d   = int(frac*t_ff/dt) ! print
  !integer, parameter :: iout_av  = 100   ! averaging frequency (not necesseraly the printing one)
  !integer, parameter :: iout2d   = 1000*20!*5
  !integer, parameter :: iout3d   = 5000*40!*5
  !
  !dimensional quantities
  !
  real(8), parameter :: m1        = lambda_mm0*m2
  real(8), parameter :: cp2       = 1006.d0!ru/(phi_p1*m2)
  real(8), parameter :: mu2       = 1.79e-5!nu2*rho2_0!3.4993e-05
  real(8), parameter :: rho1      = lambda_rho0*rho2_0 !=575.55 
  real(8), parameter :: sigmaca   = 0.0001!rho2_0*(1**2.d0)*lref/we !rho2_0*(uref**2.d0)*lref/we ! we use the density of the liquid like Lohse-Verzicco
  real(8), parameter :: mu1       = lambda_mu0*mu2!1.1143e-04
  real(8), parameter :: kappa2    = 0.026!mu2*cp2/pr
  real(8), parameter :: kappa1    = lambda_ka0*kappa2!0.4401
  real(8), parameter :: cp1       = lambda_cp0*cp2 !  4.9006e+03
  !real(8), parameter :: d_m12    = (mu2/rho2_0)*(1.d0/sc)
  real(8), parameter :: lheat     = 2.33e6!cp2*tmp0/ste !    1.0902e+06
  real(8), parameter :: delta_cp  = lambda_dcp0*cp2
  !
  integer, parameter :: e         = 0
  real(8), parameter :: sinit     = 0.d0
  real(8), parameter :: gri       = 0.d0
#ifdef USE_VOF
  real(8), parameter :: cfl       = 0.25d0 ! we have also We
#else
  real(8), parameter :: cfl       = 0.25d0 ! only convection 
#endif
  real(8), parameter :: cfl_o     = 0.90d0
  real(8), parameter :: small     = 1e-09
  real(8), parameter :: theta_thr = 0.25d0
  real(8), parameter :: mflux0    = 0.0d0
  real(8), parameter :: ext_liq   = 0.0d0
  real(8), parameter :: ext_gas   = 0.0d0 ! keep 0 
  !
  ! disperse phase 
  !
  integer, parameter :: n_dp  = 1 ! initial number of droplets
  real(8), parameter :: xcc(1:n_dp) = (/0.500d0/)*lx
  real(8), parameter :: ycc(1:n_dp) = (/0.500d0/)*ly
  real(8), parameter :: zcc(1:n_dp) = (/0.500d0/)*lz
  real(8), parameter :: rd(1:n_dp)  = lref/2.d0
  !logical, parameter :: manual_dp   = .false.
  logical, parameter :: manual_dp   = .true.
  !
  character(len=*), parameter :: datadir = 'data/'
  character(len=*), parameter :: datapos = 'data/post/'
  character(len=*), parameter :: datadir_ta   = 'data/post/tagging/'
  !
  logical, parameter, dimension(2,3) :: no_outflow = & 
      reshape((/.false.,.false.,   & ! no outflow in x lower,upper bound
                .false.,.false.,   & ! no outflow in y lower,upper bound
                .false.,.false./), & ! no outflow in z lower,upper bound
                 shape(no_outflow))
  !
  ! re-calculation of the dimensionless physical parameters
  !
  real(8), parameter :: re_c          = rho2_0*uref/(mu2*k0_wave)
  real(8), parameter :: pr_c          = mu2*cp2/kappa2
  real(8), parameter :: we_c          = rho2_0*(uref**2.d0)*lref/sigmaca !Re*Ca
  real(8), parameter :: sc_c          = nu2/d_m12                     ! only eva
  real(8), parameter :: ste_c         = cp2*tmp0/lheat                ! only eva
  real(8), parameter :: lambda_rho0_c = rho1/rho2_0                  
  real(8), parameter :: lambda_mu0_c  = mu1/mu2                      
  real(8), parameter :: lambda_cp0_c  = cp1/cp2                      
  real(8), parameter :: lambda_ka0_c  = kappa1/kappa2                
  real(8), parameter :: lambda_mm0_c  = m1/m2                         ! only eva
  real(8), parameter :: lambda_dcp0_c = delta_cp/cp2                  ! only eva
  real(8), parameter :: phi_p1_c      = (ru/(m2*cp2))                 ! only eva
  real(8), parameter :: phi_p2_c      = (pth_0/(rho2_0*(ru/m2)*tmp0))       
  real(8), parameter :: phi_p3_c      = pth_0/pc                      ! only eva
  real(8), parameter :: phi_p4_c      = tmp0/tc                       ! only eva
  !
  ! forcing
  !
  character(len=3), parameter :: turb_type = 'abc' 
  !character(len=3), parameter :: c_or_t    = 'tdp' 
  character(len=3), parameter :: c_or_t    = 'uni' 
  real(8)         , parameter, dimension(3) :: abc = (/u_abc,u_abc,u_abc/)
  real(8)         , parameter :: f0_t = nu2*(k0_wave)**2
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 1.d0    ! m1 - n1
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 10.d0   ! m1 - n2
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 50.d0   ! m1 - n3
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 2.0d0, amp_n = 100.d0  ! m1 - n3
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 1.d0    ! m2 - n1
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 10.d0   ! m2 - n2
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 50.d0   ! m2 - n3
  !real(8)         , parameter :: amp_a = 1.d0, amp_b = 1.0d0, amp_n = 100.d0  ! m2 - n3
  !real(8)         , parameter :: amp_a = 0.d0, amp_b = 20.0d0, amp_n = 1.d0    ! m3 - n1
  real(8)         , parameter :: amp_a = 0.d0, amp_b = 20.0d0, amp_n = 1.d0    ! m3 - n1
  !
end module mod_param
