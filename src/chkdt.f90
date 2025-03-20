module mod_chkdt
  !
  use mpi
  use mod_common_mpi, only: comm_cart,ierr,myid
  use mod_param     , only: mu1,rho1,cp1,cp2,kappa1,sigmaca,pi,theta_thr,gacc
  !
  implicit none
  !
  private
  public  :: chkdt
  !
  contains
  !
  subroutine chkdt(n,dl,u,v,w,rho2_min,d_lg_max,mu_g_max,ka_g_max,dtmax,dtmax_o)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dl
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in )                         :: rho2_min,d_lg_max,mu_g_max,ka_g_max
    real(8), intent(out)                         :: dtmax,dtmax_o
    !
    real(8) :: dxi,dyi,dzi
    real(8) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    !real(8) :: dtic_m,dtic_mx,dtic_my,dtic_mz,dtc_m
    real(8) :: dtic_mx,dtic_my,dtic_mz,dtic_m,dtiv,dtik,dtig,dti,dlmin,dlmini
    real(8) :: dt_mu,dt_mas,dt_ene,dtk,dtc_m 
    integer :: i,j,k,ip,jp,kp,im,jm,km
    !
    dtic_m = 0.d0
    dxi    = 1.d0/dl(1)
    dyi    = 1.d0/dl(2)
    dzi    = 1.d0/dl(3)
    !
    ! convection
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
          ! along x
          !
          ux = abs(u(i,j,k))
          vx = 0.25d0*abs( v(i,j,k)+v(i,jm,k)+v(ip,j,k)+v(ip,jm,k) )
          wx = 0.25d0*abs( w(i,j,k)+w(i,j,km)+w(ip,j,k)+w(ip,j,km) )
          dtic_mx = ux*dxi+vx*dyi+wx*dzi
          !
          ! along y
          !
          uy = 0.25d0*abs( u(i,j,k)+u(i,jp,k)+u(im,jp,k)+u(im,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25d0*abs( w(i,j,k)+w(i,jp,k)+w(i,jp,km)+w(i,j,km) )
          dtic_my = uy*dxi+vy*dyi+wy*dzi
          !
          ! along z
          !
          uz = 0.25d0*abs( u(i,j,k)+u(im,j,k)+u(im,j,kp)+u(i,j,kp) )
          vz = 0.25d0*abs( v(i,j,k)+v(i,jm,k)+v(i,jm,kp)+v(i,j,kp) )
          wz = abs(w(i,j,k))
          dtic_mz = uz*dxi+vz*dyi+wz*dzi
          !
          ! overall convective contribution
          ! 
          dtic_m = max(dtic_m,dtic_mx,dtic_my,dtic_mz)
          !
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,dtic_m,1,mpi_real8,mpi_max,comm_cart,ierr)
    if(dtic_m.eq.0.d0) dtic_m = 1.d0
    dtc_m = dtic_m**(-1.0)
    !
    ! momentum diffusion
    !
#ifdef TWOD
    dt_mu = ( max(mu1/rho1,mu_g_max/rho2_min)*(2.d0*(       dyi**2+dzi**2)) )**(-1.0)
#else
    dt_mu = ( max(mu1/rho1,mu_g_max/rho2_min)*(2.d0*(dxi**2+dyi**2+dzi**2)) )**(-1.0)
#endif
    !
    ! thermal diffusion
    !
#ifdef ENERGY
#ifdef TWOD
    dt_ene = ( max(kappa1/(cp1*rho1),ka_g_max/(cp2*rho2_min))*(2.d0*(       dyi**2+dzi**2)) )**(-1.0)
#else
    dt_ene = ( max(kappa1/(cp1*rho1),ka_g_max/(cp2*rho2_min))*(2.d0*(dxi**2+dyi**2+dzi**2)) )**(-1.0)
#endif
#else
    dt_ene = 1.d0
#endif
    !
    ! mass diffusion
    !
#ifdef VAP_MASS
#ifdef TWOD
    dt_mas = (theta_thr**2.d0)*( d_lg_max*(2.d0*(       dyi**2+dzi**2)) )**(-1.0)
#else
    dt_mas = (theta_thr**2.d0)*( d_lg_max*(2.d0*(dxi**2+dyi**2+dzi**2)) )**(-1.0)
#endif
#else
    dt_mas = 1.d0
#endif
    !
    ! surface tension
    !
#ifdef USE_VOF
    dtk = sqrt((rho1+rho2_min)*min(dl(1),dl(2),dl(3))**3/(4.d0*pi*sigmaca))
#else
    dtk = 1000.d0
#endif
    !
    if(myid.eq.0) print*, "Here the time-step restrictions: ", dtc_m,dt_mu,dt_ene,dt_mas,dtk,sigmaca,rho2_min
    !
    ! overall contribution
    !
    !dtmax   = dtc_m
    !dtmax_o = min(dt_mu,dt_ene,dt_mas,dtk)
    !
    dtiv    = 1.d0/dt_mu 
    dtik    = 1.d0/dtk
    dtig    = sqrt(maxval(abs(gacc))/(minval(dl(:))))
    !
    dtmax   = 2.d0*(dtic_m+dtiv+sqrt((dtic_m+dtiv)**2+4.d0*(dtig**2+dtik**2)))**(-1)
    dtmax_o = min(dt_ene,dt_mas)
    !
    return
  end subroutine chkdt
  !
end module mod_chkdt
