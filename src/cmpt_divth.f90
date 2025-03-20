module mod_cmpt_divth
  !
  use mpi
  use mod_param     , only: cp2,m1,m2,ru,small,theta_thr,delta_cp
  use mod_thermo    , only: thermo_rhog,mass_fraction,thermo_d_lg
  use mod_common_mpi, only: ierr,comm_cart,coord
  !
  implicit none
  !
  public  :: cmpt_divth,cmpt_sth
#ifdef VAP_MASS
  public  :: cmpt_normals,extended,cmpt_mflux_m3
#endif
  private 
  !
  contains
  !
  subroutine cmpt_divth(n,dli,ptho,pth,dpthdt_n,vof,mflux,rho1,kappa,rho_gas,d_lg,tmp,tmpge,sca,scae,div_th,divg_th)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in )                         :: ptho,pth,dpthdt_n
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: vof,mflux
    real(8), intent(in )                         :: rho1
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: kappa
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rho_gas,d_lg
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: tmpge
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: sca,scae
    real(8), intent(out), dimension( 0:, 0:, 0:) :: div_th
    real(8), intent(out), dimension( 0:, 0:, 0:) :: divg_th
    !
    real(8), dimension(8) :: mx,my,mz
    real(8) :: dvofdx,dvofdy,dvofdz,diff_s,diff_t,m_avg
    real(8) :: s_in,t_in,rhog_in,rhol_in
    real(8) :: dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm
    real(8) :: dscadxp,dscadxm,dscadyp,dscadym,dscadzp,dscadzm
    real(8) :: rhdlgxp,rhdlgxm,rhdlgyp,rhdlgym,rhdlgzp,rhdlgzm
    integer :: i,j,k,im,jm,km,ip,jp,kp
    real(8) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      do j=1,n(2)
        jp = j+1
        jm = j-1
        do i=1,n(1)
          ip = i+1
          im = i-1
          !
          ! a. calculation of the interfacial area
          !
#ifdef TWOD 
          mx(1:8) = 0.d0
#else
          !i+1/2 j+1/2 k+1/2
          mx(1) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mx(2) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mx(3) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mx(4) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mx(5) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mx(6) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mx(7) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mx(8) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)*0.25d0
#endif
          ! 
          !i+1/2 j+1/2 k+1/2
          my(1) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k+1/2
          my(2) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)*0.25d0
          !i+1/2 j+1/2 k-1/2
          my(3) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k-1/2
          my(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k+1/2
          my(5) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k+1/2
          my(6) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k-1/2
          my(7) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k-1/2
          my(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)*0.25d0
          !
          !i+1/2 j+1/2 k+1/2
          mz(1) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli(3)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mz(2) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli(3)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mz(3) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli(3)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mz(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(3)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mz(5) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli(3)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mz(6) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli(3)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mz(7) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli(3)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mz(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(3)*0.25d0
          !
          dvofdx = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dvofdy = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dvofdz = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          !
#ifdef LOW_MACH
          !
          ! b. vapour mass effects
          !
          dscadxm = (sca(i ,j,k)-sca(im,j,k))*dli(1)
          dscadxp = (sca(ip,j,k)-sca(i ,j,k))*dli(1)
          dscadym = (sca(i,j ,k)-sca(i,jm,k))*dli(2)
          dscadyp = (sca(i,jp,k)-sca(i,j ,k))*dli(2)
          dscadzm = (sca(i,j,k )-sca(i,j,km))*dli(3)
          dscadzp = (sca(i,j,kp)-sca(i,j ,k))*dli(3)
          ! 
          rhdlgxm = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(im,j,k)*d_lg(im,j,k))
          rhdlgxp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(ip,j,k)*d_lg(ip,j,k))
          rhdlgym = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(i,jm,k)*d_lg(i,jm,k))
          rhdlgyp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(i,jp,k)*d_lg(i,jp,k))
          rhdlgzm = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(i,j,km)*d_lg(i,j,km))
          rhdlgzp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)+rho_gas(i,j,kp)*d_lg(i,j,kp))
          !
          m_avg  = (sca(i,j,k)/m1+(1.d0-sca(i,j,k))/m2)**(-1.d0)
          diff_s = ((dscadxp*rhdlgxp-dscadxm*rhdlgxm)*dli(1) + &
                    (dscadyp*rhdlgyp-dscadym*rhdlgym)*dli(2) + &
                    (dscadzp*rhdlgzp-dscadzm*rhdlgzm)*dli(3))*(1.d0/rho_gas(i,j,k))*m_avg*(1.d0/m1-1.d0/m2)
          !
          ! c. thermal effects
          !
          dtmpdxm = (tmp(i ,j,k)-tmp(im,j,k))*dli(1)
          dtmpdxp = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
          dtmpdym = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
          dtmpdyp = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
          dtmpdzm = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
          dtmpdzp = (tmp(i,j,kp)-tmp(i,j ,k))*dli(3)
          !
          kappaxp = 0.5d0*(kappa(ip,j,k)+kappa(i,j,k))
          kappaxm = 0.5d0*(kappa(im,j,k)+kappa(i,j,k))
          kappayp = 0.5d0*(kappa(i,jp,k)+kappa(i,j,k))
          kappaym = 0.5d0*(kappa(i,jm,k)+kappa(i,j,k))
          kappazp = 0.5d0*(kappa(i,j,kp)+kappa(i,j,k))
          kappazm = 0.5d0*(kappa(i,j,km)+kappa(i,j,k))
          ! 
          diff_t  = ( (kappaxp*dtmpdxp-kappaxm*dtmpdxm)*dli(1) + &
                      (kappayp*dtmpdyp-kappaym*dtmpdym)*dli(2) + &
                      (kappazp*dtmpdzp-kappazm*dtmpdzm)*dli(3) )*(ru/(cp2*m_avg)) + &
                    delta_cp*d_lg(i,j,k)*( &
                             (dscadxp+dscadxm)*(dtmpdxp+dtmpdxm) + &
                             (dscadyp+dscadym)*(dtmpdyp+dtmpdym) + &
                             (dscadzp+dscadzm)*(dtmpdzp+dtmpdzm)   &
                                         )*0.25d0*rho_gas(i,j,k)*(ru/(cp2*m_avg)) 
          !
#else
          !
          m_avg  = 1.d0 ! arbitrary number to avoid the division by 0 in absence of low Mach.
          diff_t = 0.d0
          diff_s = 0.d0
          !
#endif
          !
          ! d. we have to compute an estimation of the tmp and s_in at the interface 
          !
#ifdef LOW_MACH
          !
          s_in = scae(i,j,k)
          t_in = tmpge(i,j,k) ! for now
          !
#else
          !
          s_in = 0.d0
          t_in = tmp(i,j,k) 
          !
#endif
          !
          rhog_in = thermo_rhog(ptho,t_in,s_in) ! to be consistent with rho_gas after VoF
          rhol_in = rho1
          !
          ! e. put all together
          !
#ifdef LOW_MACH
          divg_th(i,j,k) = - (1.d0-vof(i,j,k))*(1.d0-(ru/(cp2*m_avg)))*dpthdt_n + &
                             diff_s*(1.d0-vof(i,j,k))*(ptho/ptho) + &
                             diff_t*(1.d0-vof(i,j,k))*(1.d0/ptho)
#else
          divg_th(i,j,k) = 0.d0
#endif
          !
          div_th(i,j,k)  = divg_th(i,j,k) + & 
                           mflux(i,j,k)*(1.d0/rhog_in-1.d0/rhol_in)*sqrt(dvofdx**2+dvofdy**2+dvofdz**2)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_divth
  !
  subroutine cmpt_normals(n,dli,phi,normx,normy,normz)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: phi
    real(8), intent(out), dimension( 0:, 0:, 0:) :: normx,normy,normz
    !
    real(8), dimension(8) :: mx,my,mz
    integer :: i,j,k,im,jm,km,ip,jp,kp
    real(8) :: norm
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      do j=1,n(2)
        jp = j+1
        jm = j-1
        do i=1,n(1)
          ip = i+1
          im = i-1
          !
          mx(1)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mx(2)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1)))*dli(1)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mx(3)=((phi(i+1,j  ,k  )+phi(i+1,j+1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mx(4)=((phi(i+1,j  ,k  )+phi(i+1,j-1,k  )+phi(i+1,j  ,k-1)+phi(i+1,j-1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mx(5)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j+1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j+1,k+1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mx(6)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k+1)+phi(i  ,j-1,k+1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k+1)+phi(i-1,j-1,k+1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mx(7)=((phi(i  ,j  ,k  )+phi(i  ,j+1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j+1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j+1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j+1,k-1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mx(8)=((phi(i  ,j  ,k  )+phi(i  ,j-1,k  )+phi(i  ,j  ,k-1)+phi(i  ,j-1,k-1))-&
                 (phi(i-1,j  ,k  )+phi(i-1,j-1,k  )+phi(i-1,j  ,k-1)+phi(i-1,j-1,k-1)))*dli(1)*0.25d0
          !
          !i+1/2 j+1/2 k+1/2
          my(1)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k+1/2
          my(2)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1)))*dli(2)*0.25d0
          !i+1/2 j+1/2 k-1/2
          my(3)=((phi(i  ,j+1,k  )+phi(i+1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k-1/2
          my(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i+1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k+1/2
          my(5)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k+1/2
          my(6)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k-1/2
          my(7)=((phi(i  ,j+1,k  )+phi(i-1,j+1,k  )+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k-1/2
          my(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1))-&
                 (phi(i  ,j-1,k  )+phi(i-1,j-1,k  )+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(2)*0.25d0
      
          !i+1/2 j+1/2 k+1/2
          mz(1)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i+1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  )))*dli(3)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mz(2)=((phi(i  ,j  ,k+1)+phi(i+1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i+1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  )))*dli(3)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mz(3)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j+1,k  )+phi(i+1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i+1,j+1,k-1)))*dli(3)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mz(4)=((phi(i  ,j  ,k  )+phi(i+1,j  ,k  )+phi(i  ,j-1,k  )+phi(i+1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i+1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i+1,j-1,k-1)))*dli(3)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mz(5)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j+1,k+1)+phi(i-1,j+1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  )))*dli(3)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mz(6)=((phi(i  ,j  ,k+1)+phi(i-1,j  ,k+1)+phi(i  ,j-1,k+1)+phi(i-1,j-1,k+1))-&
                 (phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  )))*dli(3)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mz(7)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j+1,k  )+phi(i-1,j+1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j+1,k-1)+phi(i-1,j+1,k-1)))*dli(3)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mz(8)=((phi(i  ,j  ,k  )+phi(i-1,j  ,k  )+phi(i  ,j-1,k  )+phi(i-1,j-1,k  ))-&
                 (phi(i  ,j  ,k-1)+phi(i-1,j  ,k-1)+phi(i  ,j-1,k-1)+phi(i-1,j-1,k-1)))*dli(3)*0.25d0
          !
          normx(i,j,k) = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          normy(i,j,k) = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          normz(i,j,k) = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          norm         = sqrt(normx(i,j,k)**2+normy(i,j,k)**2+normz(i,j,k)**2+1.0e-16)
          normx(i,j,k) = normx(i,j,k)/norm
          normy(i,j,k) = normy(i,j,k)/norm
          normz(i,j,k) = normz(i,j,k)/norm
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_normals
  !
  subroutine cmpt_mflux_m3(n,dli,qg,ql,pth,phi,tmpg,tmpl,s,scae,nx,ny,nz,mflux,rdir)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)              :: n
    real(8), intent(in ), dimension(3)              :: dli
    integer, intent(in )                            :: qg,ql
    real(8), intent(in )                            :: pth
    real(8), intent(in ), dimension( -2:, -2:, -2:) :: phi
    real(8), intent(in ), dimension(-qg:,-qg:,-qg:) :: tmpg
    real(8), intent(in ), dimension(-ql:,-ql:,-ql:) :: tmpl
    real(8), intent(in ), dimension(  0:,  0:,  0:) :: s
    real(8), intent(in ), dimension(  0:,  0:,  0:) :: scae
    real(8), intent(in ), dimension(  0:,  0:,  0:) :: nx,ny,nz
    real(8), intent(out), dimension(  0:,  0:,  0:) :: mflux
    real(8), intent(in )                            :: rdir
    !
    real(8), dimension(3) :: dl
    real(8) :: dsdx,dsdy,dsdz,ss_ip,ss_im,grad_sn
    real(8) :: dsdxp,dsdyp,dsdzp,dsdxm,dsdym,dsdzm
    real(8) :: sl_in,tl_in,tg_in
    real(8) :: d_lg_in,rhog_in
    real(8) :: tmpl_ip,tmpl_im,tmpg_ip,tmpg_im,thetap,thetam
    real(8) :: tmpg_ix,tmpg_iy,tmpg_iz,tmpl_ix,tmpl_iy,tmpl_iz
    real(8) :: dlp,dlm,xx,yy,zz
    integer :: im,ip,jm,jp,km,kp,i,j,k
    !
    dl = dli**(-1.d0)
    !
    do k=1,n(3)
      kp = k + 1
      km = k - 1
      zz = (k-0.5d0)*dl(3)
      do j=1,n(2)
        jp = j + 1
        jm = j - 1
        yy = (j+coord(2)*n(2)-0.5d0)*dl(2)
        do i=1,n(1)
          ip = i + 1
          im = i - 1
          xx = (i+coord(1)*n(1)-0.5d0)*dl(1)
          !
          if(phi(i,j,k).ge.0.d0) then
            !
            mflux(i,j,k) = 0.d0
            !
          else
            !
            dsdxp = (s(ip,j,k)-s(i ,j,k))*dli(1)
            dsdxm = (s(i ,j,k)-s(im,j,k))*dli(1)
            dsdyp = (s(i,jp,k)-s(i,j ,k))*dli(2)
            dsdym = (s(i,j ,k)-s(i,jm,k))*dli(2)
            dsdzp = (s(i,j,kp)-s(i,j,k ))*dli(3)
            dsdzm = (s(i,j,k )-s(i,j,km))*dli(3)
            !
            dsdx = 0.5d0*(dsdxp+dsdxm)
            dsdy = 0.5d0*(dsdyp+dsdym)
            dsdz = 0.5d0*(dsdzp+dsdzm)
            !
            tmpg_ix = tmpg(i,j,k) ! gas
            tmpg_iy = tmpg(i,j,k) ! gas
            tmpg_iz = tmpg(i,j,k) ! gas
            !
            tmpl_ix = tmpl(i,j,k) ! liq
            tmpl_iy = tmpl(i,j,k) ! liq
            tmpl_iz = tmpl(i,j,k) ! liq
            !
            ! along x
            !  --> i   case: occur positive normal component, nx and dsdx is positive;
            !  --> ii  case: occur negative normal component, nx and dsdx is negative;
            !  --> iii case: zero.
            !
            if(     phi(ip,j,k)*phi(i,j,k).lt.0.d0.and.phi(im,j,k)*phi(i,j,k).gt.0.d0 ) then ! im,i same side, ip not
              tmpl_ip = (tmpl(i,j,k)*abs(phi(ip,j,k))+tmpl(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              thetap  = abs(phi(i,j,k))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              if(thetap.gt.theta_thr) then
                !dsdxp = + (ss_ip   -s(i ,j,k))*dli(1)/thetap
                dsdx = + (ss_ip   -s(i ,j,k))*dli(1)/thetap
              else
                !dsdxp = + (ss_ip   -s(im,j,k))*dli(1)/(1.d0+thetap)
                dsdx = + (ss_ip   -s(im,j,k))*dli(1)/(1.d0+thetap)
              endif
              tmpl_ix = tmpl_ip
              tmpg_ix = (tmpg(i,j,k)*abs(phi(ip,j,k))+tmpg(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
            elseif( phi(ip,j,k)*phi(i,j,k).gt.0.d0.and.phi(im,j,k)*phi(i,j,k).lt.0.d0 ) then ! ip,i same side, im not
              tmpl_im = (tmpl(i,j,k)*abs(phi(im,j,k))+tmpl(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              thetam  = abs(phi(i,j,k))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              if(thetam.gt.theta_thr) then
                !dsdxm = - (ss_im   -s(i ,j,k))*dli(1)/thetam
                dsdx = - (ss_im   -s(i ,j,k))*dli(1)/thetam
              else
                !dsdxm = - (ss_im   -s(ip,j,k))*dli(1)/(1.d0+thetam)
                dsdx = - (ss_im   -s(ip,j,k))*dli(1)/(1.d0+thetam)
              endif
              tmpl_ix = tmpl_im
              tmpg_ix = (tmpg(i,j,k)*abs(phi(im,j,k))+tmpg(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k))) 
            elseif( phi(ip,j,k)*phi(i,j,k).lt.0.d0.and.phi(im,j,k)*phi(i,j,k).lt.0.d0 ) then ! i outside interface, im,ip inside
              tmpl_ip = (tmpl(i,j,k)*abs(phi(ip,j,k))+tmpl(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmpl_im = (tmpl(i,j,k)*abs(phi(im,j,k))+tmpl(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              tmpg_ip = (tmpg(i,j,k)*abs(phi(ip,j,k))+tmpg(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmpg_im = (tmpg(i,j,k)*abs(phi(im,j,k))+tmpg(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              !ss_ip  = mass_fraction(pth,tmp_ip)
              !ss_im  = mass_fraction(pth,tmp_im)
              !dlp    = (xx*abs(phi(ip,j,k))+(xx+dl(1))*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              !dlm    = (xx*abs(phi(im,j,k))+(xx-dl(1))*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              !dsdxp  = 0.d0
              !dsdxm  = 0.d0
              dsdx = 0.d0
              !if(abs(dlp-dlm).gt.2.d0*theta_thr*dl(1)) then
              !  dsdxp = + (ss_ip   -s(i,j,k))/abs(dlp-dlm)
              !  dsdxm = - (ss_im   -s(i,j,k))/abs(dlp-dlm)
              !else
              !  dsdxp = 0.d0
              !  dsdxm = 0.d0
              !endif
              tmpl_ix = 0.5d0*(tmpl_ip+tmpl_im)
              tmpg_ix = 0.5d0*(tmpg_ip+tmpg_im)
            endif
            !
            ! along y
            !  --> i   case: occur positive normal component, ny and dsdy is positive;
            !  --> ii  case: occur negative normal component, ny and dsdy is negative;
            !  --> iii case: zero.
            !
            if(     phi(i,jp,k)*phi(i,j,k).lt.0.d0.and.phi(i,jm,k)*phi(i,j,k).gt.0.d0 ) then ! jm,j same side, jp not
              tmpl_ip = (tmpl(i,j,k)*abs(phi(i,jp,k))+tmpl(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              thetap  = abs(phi(i,j,k))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              if(thetap.gt.theta_thr) then
                !dsdyp = + (ss_ip   -s(i,j ,k))*dli(2)/thetap
                dsdy = + (ss_ip   -s(i,j ,k))*dli(2)/thetap
              else
                !dsdyp = + (ss_ip   -s(i,jm,k))*dli(2)/(1.d0+thetap)
                dsdy = + (ss_ip   -s(i,jm,k))*dli(2)/(1.d0+thetap)
              endif
              tmpl_iy = tmpl_ip
              tmpg_iy = (tmpg(i,j,k)*abs(phi(i,jp,k))+tmpg(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
            elseif( phi(i,jp,k)*phi(i,j,k).gt.0.d0.and.phi(i,jm,k)*phi(i,j,k).lt.0.d0 ) then ! jp,j same side, jm not
              tmpl_im = (tmpl(i,j,k)*abs(phi(i,jm,k))+tmpl(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              thetam  = abs(phi(i,j,k))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              if(thetam.gt.theta_thr) then
                !dsdym = - (ss_im   -s(i,j ,k))*dli(2)/thetam
                dsdy = - (ss_im   -s(i,j ,k))*dli(2)/thetam
              else
                !dsdym = - (ss_im   -s(i,jp,k))*dli(2)/(1.d0+thetam)
                dsdy = - (ss_im   -s(i,jp,k))*dli(2)/(1.d0+thetam)
              endif
              tmpl_iy = tmpl_im
              tmpg_iy = (tmpg(i,j,k)*abs(phi(i,jm,k))+tmpg(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
            elseif( phi(i,jp,k)*phi(i,j,k).lt.0.d0.and.phi(i,jm,k)*phi(i,j,k).lt.0.d0 ) then ! j outside interface, jm,jp inside
              tmpl_ip = (tmpl(i,j,k)*abs(phi(i,jp,k))+tmpl(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmpl_im = (tmpl(i,j,k)*abs(phi(i,jm,k))+tmpl(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              tmpg_ip = (tmpg(i,j,k)*abs(phi(i,jp,k))+tmpg(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmpg_im = (tmpg(i,j,k)*abs(phi(i,jm,k))+tmpg(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              !ss_ip  = mass_fraction(pth,tmp_ip)
              !ss_im  = mass_fraction(pth,tmp_im)
              !dlp    = (yy*abs(phi(i,jp,k))+(yy+dl(2))*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              !dlm    = (yy*abs(phi(i,jm,k))+(yy-dl(2))*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              !dsdyp  = 0.d0
              !dsdym  = 0.d0
              dsdy  = 0.d0
              !if(abs(dlp-dlm).gt.2.d0*theta_thr*dl(2)) then
              !  dsdyp = + (ss_ip   -s(i,j,k))/abs(dlp-dlm)
              !  dsdym = - (ss_im   -s(i,j,k))/abs(dlp-dlm)
              !else
              !  dsdyp = 0.d0
              !  dsdym = 0.d0
              !endif
              tmpl_iy = 0.5d0*(tmpl_ip+tmpl_im)
              tmpg_iy = 0.5d0*(tmpg_ip+tmpg_im)
            endif
            !
            ! along z
            !  --> i   case: occur positive normal component, nz and dsdz is positive;
            !  --> ii  case: occur negative normal component, nz and dsdz is negative;
            !  --> iii case: zero.
            !
            if(     phi(i,j,kp)*phi(i,j,k).lt.0.d0.and.phi(i,j,km)*phi(i,j,k).gt.0.d0 ) then ! km,k same side, kp not
              tmpl_ip = (tmpl(i,j,k)*abs(phi(i,j,kp))+tmpl(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              thetap  = abs(phi(i,j,k))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              if(thetap.gt.theta_thr) then
                !dsdzp = + (ss_ip   -s(i,j,k ))*dli(3)/thetap
                dsdz = + (ss_ip   -s(i,j,k ))*dli(3)/thetap
              else
                !dsdzp = + (ss_ip   -s(i,j,km))*dli(3)/(1.d0+thetap)
                dsdz = + (ss_ip   -s(i,j,km))*dli(3)/(1.d0+thetap)
              endif
              tmpl_iz = tmpl_ip
              tmpg_iz = (tmpg(i,j,k)*abs(phi(i,j,kp))+tmpg(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
            elseif( phi(i,j,kp)*phi(i,j,k).gt.0.d0.and.phi(i,j,km)*phi(i,j,k).lt.0.d0 ) then ! kp,k same side, km not
              tmpl_im = (tmpl(i,j,k)*abs(phi(i,j,km))+tmpl(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              thetam  = abs(phi(i,j,k))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              if(thetam.gt.theta_thr) then
                !dsdzm = - (ss_im   -s(i,j,k ))*dli(3)/thetam
                dsdz = - (ss_im   -s(i,j,k ))*dli(3)/thetam
              else
                !dsdzm = - (ss_im   -s(i,j,kp))*dli(3)/(1.d0+thetam)
                dsdz = - (ss_im   -s(i,j,kp))*dli(3)/(1.d0+thetam)
              endif
              tmpl_iz = tmpl_im
              tmpg_iz = (tmpg(i,j,k)*abs(phi(i,j,km))+tmpg(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
            elseif( phi(i,j,kp)*phi(i,j,k).lt.0.d0.and.phi(i,j,km)*phi(i,j,k).lt.0.d0 ) then ! k outside interface, km,kp inside
              tmpl_ip = (tmpl(i,j,k)*abs(phi(i,j,kp))+tmpl(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmpl_im = (tmpl(i,j,k)*abs(phi(i,j,km))+tmpl(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              tmpg_ip = (tmpg(i,j,k)*abs(phi(i,j,kp))+tmpg(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmpg_im = (tmpg(i,j,k)*abs(phi(i,j,km))+tmpg(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              !ss_ip  = mass_fraction(pth,tmp_ip)
              !ss_im  = mass_fraction(pth,tmp_im)
              !dlp    = (zz*abs(phi(i,j,kp))+(zz+dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              !dlm    = (zz*abs(phi(i,j,km))+(zz-dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              !dsdzp  = 0.d0
              !dsdzm  = 0.d0
              dsdz  = 0.d0
              !if(abs(dlp-dlm).gt.2.d0*theta_thr*dl(3)) then
              !  dsdzp = + (ss_ip   -s(i,j,k))/abs(dlp-dlm)
              !  dsdzm = - (ss_im   -s(i,j,k))/abs(dlp-dlm)
              !else
              !  dsdzp = 0.d0
              !  dsdzm = 0.d0
              !endif
              tmpl_iz = 0.5d0*(tmpl_ip+tmpl_im)
              tmpg_iz = 0.5d0*(tmpg_ip+tmpg_im)
            endif
            !
            !sn(i,j,k) = 0.5d0*(nx(i,j,k)*(dsdxp+dsdxm)+ny(i,j,k)*(dsdyp+dsdym)+nz(i,j,k)*(dsdzp+dsdzm))*rdir
            !sn(i,j,k) = ( 0.5d0*(nx(i,j,k)-abs(nx(i,j,k)))*dsdxm + 0.5d0*(nx(i,j,k)+abs(nx(i,j,k)))*dsdxp + &
            !              0.5d0*(ny(i,j,k)-abs(ny(i,j,k)))*dsdym + 0.5d0*(ny(i,j,k)+abs(ny(i,j,k)))*dsdyp + &
            !              0.5d0*(nz(i,j,k)-abs(nz(i,j,k)))*dsdzm + 0.5d0*(nz(i,j,k)+abs(nz(i,j,k)))*dsdzp )*rdir
            !sn(i,j,k) = 0.5d0*(nx(i,j,k)*(dsdxp+dsdxm)+ny(i,j,k)*(dsdyp+dsdym)+nz(i,j,k)*(dsdzp+dsdzm))*rdir
            grad_sn = (nx(i,j,k)*dsdx+ny(i,j,k)*dsdy+nz(i,j,k)*dsdz)*rdir
            !
            ! Note that to compute the mflux, we need: 
            !   --> sl_in distributed in a band;
            !   --> tg_in and tl_in distributed in a band;
            !   --> rhog_in and d_lg_in distributed in a band.
            !       These properties are computed from the corresponding constitutive laws 
            !       using pth (uniform by definition), sl_in, tg_in and tl_in.
            !
            tl_in   = tmpl_ix*nx(i,j,k)**2.d0+tmpl_iy*ny(i,j,k)**2.d0+tmpl_iz*nz(i,j,k)**2.d0 ! liq temperature
            !sl_in   = scae(i,j,k)     ! to be consistent at the interface
            sl_in   = mass_fraction(pth,tl_in)     ! to be consistent at the interface
            tg_in   = tmpg_ix*nx(i,j,k)**2.d0+tmpg_iy*ny(i,j,k)**2.d0+tmpg_iz*nz(i,j,k)**2.d0 ! gas temperature
            rhog_in = thermo_rhog(pth,tg_in,sl_in) ! you have to use sl_in by convection in the thermo_rhog fun.
            d_lg_in = thermo_d_lg(pth,tg_in,sl_in)
            !
            mflux(i,j,k) = ((rhog_in*d_lg_in)/(1.0d0-sl_in+small))*( grad_sn - 0.0d0)
            !
            !mflux(i,j,k) = max(0.d0,mflux(i,j,k))
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_mflux_m3
  !
  subroutine extended(n,dli,dzc,dzf,qmin,phi,nx,ny,nz,q,qext,rdir,delta)
    !
    use mod_bound, only: boundp
    use mod_param, only: bcvof,cbcvof
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                    :: n
    real(8), intent(in ), dimension(3)                    :: dli
    real(8), intent(in ), dimension(-2:)                  :: dzc,dzf
    integer, intent(in )                                  :: qmin
    real(8), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(8), intent(in ), dimension(    0:,    0:,    0:) :: nx,ny,nz
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: q
    real(8), intent(out), dimension(    0:,    0:,    0:) :: qext
    real(8), intent(in )                                  :: rdir,delta
    !
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: phi1
    real(8), dimension(3) :: dl
    integer :: i,j,k,im,jm,km,ip,jp,kp,ax,ay,az,s,p
    integer :: ps_step
    real(8) :: num,den,eps_mfx
    !
    real(8), parameter :: eps = 1e-08
    !
    if(delta.eq.0.d0) then
      ps_step = 5
      eps_mfx = 3.0d0
    else
      ps_step = 10
      eps_mfx = 6.0d0
    endif
    !
    dl(:) = dli(:)**(-1.d0)
    !
    ! 0. define the level-set curve:
    !    note: --> phi1, level 1;
    !          --> normal defined unless of an constant offset, so no need to be recomputed.
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          !
          phi1(i,j,k) = phi(i,j,k)+rdir*delta*dl(2)
          !
        enddo
      enddo
    enddo
    !
    qext(:,:,:) = 0.d0 ! initialize to zero
    !
    ! 1. first populate the cell cut at least in one direction by
    !    the interface (i.e., the interfacial cells)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ax = nint(sign(1.d0,nx  (i,j,k)))
          ay = nint(sign(1.d0,ny  (i,j,k)))
          az = nint(sign(1.d0,nz  (i,j,k)))
          p  = nint(sign(1.d0,phi1(i,j,k)))
          !
          if(phi1(i,j,k)*rdir.gt.0.d0) then
            !
            ! case 1 (cut along all the direction)
            !
            if(     phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then 
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3) 
              den = ax*nx(i,j,k)*dli(1) + &
                    ay*ny(i,j,k)*dli(2) + &
                    az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps) 
            !
            ! case 2 (cut along only x,y, not z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then 
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2)
              den = ax*nx(i,j,k)*dli(1) + &
                    ay*ny(i,j,k)*dli(2) 
              qext(i,j,k) = num/(den+eps) 
            !
            ! case 3 (cut along only x,z, not y)
            ! 
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then 
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = ax*nx(i,j,k)*dli(1) + &
                    az*nz(i,j,k)*dli(3) 
              qext(i,j,k) = num/(den+eps)
            !
            ! case 4 (cut along only y,z, not x)
            ! 
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then 
              num = ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2) + &
                    az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = ay*ny(i,j,k)*dli(2) + &
                    az*nz(i,j,k)*dli(3) 
              qext(i,j,k) = num/(den+eps)
            !
            ! case 5 (cut along only x, not y and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then 
              num = ax*nx(i,j,k)*q(i-p*ax,j,k)*dli(1) 
              den = ax*nx(i,j,k)*dli(1)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 6 (cut along only y, not x and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).lt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then 
              num = ay*ny(i,j,k)*q(i,j-p*ay,k)*dli(2) 
              den = ay*ny(i,j,k)*dli(2)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 7 (cut along only z, not y and z)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).lt.0.d0 ) then 
              num = az*nz(i,j,k)*q(i,j,k-p*az)*dli(3)
              den = az*nz(i,j,k)*dli(3)
              qext(i,j,k) = num/(den+eps)
            !
            ! case 8 (no cut)
            !
            elseif( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                    phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then 
              num = 0.d0
              den = 1.d0 ! arbitrary value
              qext(i,j,k) = num/(den+eps)
            !
            endif
            !
          else
            !
            ! case 8 (no cut)
            !
            if( phi1(i,j,k)*phi1(i-p*ax,j,k).gt.0.d0.and. &
                phi1(i,j,k)*phi1(i,j-p*ay,k).gt.0.d0.and. &
                phi1(i,j,k)*phi1(i,j,k-p*az).gt.0.d0 ) then 
                qext(i,j,k) = 0.d0 !q(i,j,k) !num/den
            else
            !
            ! all the others (with at least one cut)
            !
                qext(i,j,k) = q(i,j,k)
            endif
            !
            !qext(i,j,k) = q(i,j,k)
            !qext(i,j,k) = 0.d0
            !
          endif
          !
        enddo
      enddo
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,qext)
    !
    ! 2. iterate to fill all the cells at least within 3 grid-cells
    !    away from the interface
    !    Note: --> case 8 of previous list accounted as it is sufficient simply not to update qext;
    !          --> we extend in both directions (-/+ normal) using the interger p;
    !          --> ps_step defined above.
    !
    do s=1,ps_step
      !
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            ax = nint(sign(1.d0,nx  (i,j,k)))
            ay = nint(sign(1.d0,ny  (i,j,k)))
            az = nint(sign(1.d0,nz  (i,j,k)))
            p  = nint(sign(1.d0,phi1(i,j,k))) ! so we account both direction (-/+ along the normal)
            !
            if(abs(phi1(i,j,k))*dli(2).lt.eps_mfx) then
              !
              ! case 1 (available neighbours: x,y,z)
              !
              if(     qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3) 
                den = ax*nx(i,j,k)*dli(1) + &
                      ay*ny(i,j,k)*dli(2) + &
                      az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 2 (available neighbours: x,y, not z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2) 
                den = ax*nx(i,j,k)*dli(1) + &
                      ay*ny(i,j,k)*dli(2) 
                qext(i,j,k) = num/(den+eps)
              !
              ! case 3 (available neighbours: x,z, not y)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3) 
                den = ax*nx(i,j,k)*dli(1) + &
                      az*nz(i,j,k)*dli(3) 
                qext(i,j,k) = num/(den+eps)
              !
              ! case 4 (available neighbours: y,z, not x)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2) + &
                      az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3) 
                den = ay*ny(i,j,k)*dli(2) + &
                      az*nz(i,j,k)*dli(3) 
                qext(i,j,k) = num/(den+eps)
              !
              ! case 5 (available neighbours: x, not y,z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).ne.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ax*nx(i,j,k)*qext(i-p*ax,j,k)*dli(1)
                den = ax*nx(i,j,k)*dli(1)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 6 (available neighbours: y, not x,z)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).ne.0.d0.and. &
                                              qext(i,j,k-p*az).eq.0.d0      ) then
                num = ay*ny(i,j,k)*qext(i,j-p*ay,k)*dli(2)
                den = ay*ny(i,j,k)*dli(2)
                qext(i,j,k) = num/(den+eps)
              !
              ! case 7 (available neighbours: z, not x,y)
              !
              elseif( qext(i,j,k).eq.0.d0.and.qext(i-p*ax,j,k).eq.0.d0.and. &
                                              qext(i,j-p*ay,k).eq.0.d0.and. &
                                              qext(i,j,k-p*az).ne.0.d0      ) then
                num = az*nz(i,j,k)*qext(i,j,k-p*az)*dli(3)
                den = az*nz(i,j,k)*dli(3)
                qext(i,j,k) = num/(den+eps)
              !
              endif
              !
            else
              !
              !qext(i,j,k) = q(i,j,k)
              qext(i,j,k) = 0.d0
              !
            endif
            !
            !qext(i,j,k) = max(0.d0,qext(i,j,k))
            !
          enddo
        enddo
      enddo
      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,qext)
      !
    enddo
    !
    return
  end subroutine extended
  !
  subroutine cmpt_sth(n,dli,vof,rho_gas,mflux,mgas,mfxt) ! not nice name
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: vof
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rho_gas,mflux
    real(8), intent(out)                         :: mgas,mfxt
    !
    real(8), dimension(8) :: mx,my,mz
    real(8), dimension(3) :: dl
    real(8) :: dvofdx,dvofdy,dvofdz
    integer :: i,j,k
    !
    dl(:) = dli(:)**(-1.d0)
    !
    mgas = 0.d0
    mfxt = 0.d0 ! initial value
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          mgas = mgas + rho_gas(i,j,k)*(1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)
          !
#ifdef TWOD 
          mx(1:8) = 0.d0
#else     
          !i+1/2 j+1/2 k+1/2
          mx(1) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mx(2) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mx(3) = ((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mx(4) = ((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mx(5) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mx(6) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mx(7) = ((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mx(8) = ((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)) - &
                   (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)*0.25d0
#endif    
          ! 
          !i+1/2 j+1/2 k+1/2
          my(1) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k+1/2
          my(2) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)*0.25d0
          !i+1/2 j+1/2 k-1/2
          my(3) = ((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k-1/2
          my(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k+1/2
          my(5) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k+1/2
          my(6) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k-1/2
          my(7) = ((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k-1/2
          my(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)) - &
                   (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)*0.25d0
          !
          !i+1/2 j+1/2 k+1/2
          mz(1) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli(3)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mz(2) = ((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli(3)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mz(3) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli(3)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mz(4) = ((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(3)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mz(5) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli(3)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mz(6) = ((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)) - &
                   (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli(3)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mz(7) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli(3)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mz(8) = ((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )) - &
                   (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(3)*0.25d0
          !
          dvofdx = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dvofdy = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dvofdz = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          !
          mfxt = mfxt + mflux(i,j,k)*sqrt(dvofdx**2+dvofdy**2+dvofdz**2)*dl(1)*dl(2)*dl(3)
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,mgas,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mfxt,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    return
  end subroutine cmpt_sth
  !
  subroutine tag_dp(n,dli,vof,vof_index)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)           :: n
    real(8), intent(in   ), dimension(3)           :: dli
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: vof
    real(8), intent(inout), dimension( 0:, 0:, 0:) :: vof_index
    !
    !integer, dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vof_index
    integer :: i,j,k!,im,jm,km,ip,jp,kp
    integer :: ii,jj,kk
    !
    real(8), parameter :: small_tg = 1e-02
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(vof(i,j,k).lt.small_tg) then
            !
            vof_index(i,j,k) = 0.d0
            !
          elseif(vof(i,j,k).ge.small_tg.and.int(vof_index(i,j,k)).eq.0) then
            !
            kloop : do kk = k-1,k+1
              jloop: do jj = j-1,j+1
               iloop: do ii = i-1,i+1
                 !
                 if(int(vof_index(ii,jj,kk)).ne.0) then
                   !
                   vof_index(i,j,k) = vof_index(ii,jj,kk)
                   exit iloop 
                   exit jloop
                   exit kloop
                   !
                 endif
                 !
               end do iloop
              end do jloop
            end do kloop
            !
          endif
          !
          !vof_rindex(i,j,k) = 1.d0*vof_index(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine tag_dp
  !
end module mod_cmpt_divth
