module mod_thermo
  !
  !use mod_common_mpi, only: istep
  !
  implicit none
  !
  real(8), parameter :: tmin = 0.d0
  !
  contains
  !
  real(8) function mass_fraction(pth,tmp)
    !
    use mod_param, only: m1,m2,tmpl0,tc,pc,lheat,ru 
    !
    real(8), intent(in) :: pth,tmp
    !
    real(8), parameter :: a = -7.28305, &  ! Ammonia
                          b = +1.59166, &  ! Ammonia    
                          c = -1.95126, &  ! Ammonia
                          d = -2.05623     ! Ammonia
    !
    real(8), parameter :: thr = 0.98d0
    real(8), parameter :: eps = 1.0e-09
    !
    real(8) :: pvap,eta
    !
    if(tmp.le.tmin) then
      eta = 1.d0-tmpl0/tc
    else
      eta = 1.d0-tmp  /tc
    endif
    !pvap = pc*exp((a*eta**1.0+b*eta**1.5+c*eta**2.5+d*eta**5.0)/(1.d0-eta))
    pvap = pth*exp(-(lheat*m1/ru)*(1.d0/tmp-1.d0/373.15))
    mass_fraction = pvap*m1/(pvap*m1+(pth-pvap)*m2)
    mass_fraction = max(eps,min(mass_fraction,thr))
    !
    return
  end function mass_fraction
  !
  real(8) function thermo_rhog(pth,tmp,yl)
    !
#ifdef LOW_MACH
    use mod_param, only: ru,m1,m2,tmpg0
#else
    use mod_param, only: rho2_0
#endif
    use mod_param, only: tmpg0
    !
    real(8), intent(in) :: pth,tmp,yl
    !
#ifdef LOW_MACH
    if(tmp.le.tmin) then
      thermo_rhog = ((yl/m1+(1.d0-yl)/m2)**(-1.d0))*pth/(ru*tmpg0) ! 3rd level of approx
    else 
      thermo_rhog = ((yl/m1+(1.d0-yl)/m2)**(-1.d0))*pth/(ru*tmp  ) ! 3rd level of approx
    endif
#else
    thermo_rhog = rho2_0
#endif
    !
    return
  end function thermo_rhog
  !
  real(8) function thermo_d_lg(pth,tmp,yl)
    !
    ! modified Wilke Lee's law
    !
    use mod_param, only: d_m12,tmpg0,pth_0,ru,m1,m2,mu2
    !
    real(8), intent(in) :: pth,tmp,yl
    !
#ifdef VAR_D_LG
    if(tmp.le.tmin) then
      thermo_d_lg = d_m12*((tmpg0/tmpg0)**(3.d0/2.d0))*(pth_0/pth)
    else
      thermo_d_lg = d_m12*((tmp  /tmpg0)**(3.d0/2.d0))*(pth_0/pth)
    endif
#else
    thermo_d_lg = d_m12
#endif
    ! 
    return
  end function thermo_d_lg
  !
  real(8) function thermo_mug(tmp)
    !
    ! modified Sutherland's law
    ! 
    use mod_param, only: tmpg0,mu2
    !
    real(8), intent(in) :: tmp
    !
#ifdef VAR_MUG
    if(tmp.le.tmin) then
      thermo_mug = mu2*((tmpg0/tmpg0)**(2.d0/3.d0))
    else
      thermo_mug = mu2*((tmp  /tmpg0)**(2.d0/3.d0))
    endif
#else
    thermo_mug = mu2
#endif
    !
    return
  end function thermo_mug
  !
  real(8) function thermo_kag(tmp)
    !
    ! kappa from modified Sutherland's law to ensure the same Pr
    ! 
    use mod_param, only: tmpg0,mu2,cp2,kappa2
    !
    real(8), intent(in) :: tmp
    !
#ifdef VAR_KAG
    if(tmp.le.tmin) then
      thermo_kag = (mu2*((tmpg0/tmpg0)**(2.d0/3.d0)))*cp2/pr
    else
      thermo_kag = (mu2*((tmp  /tmpg0)**(2.d0/3.d0)))*cp2/pr
    endif
#else
    thermo_kag = kappa2
#endif
    !
    return
  end function thermo_kag
  ! 
end module mod_thermo
