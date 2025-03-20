module mod_inte_vel
  !
  use mod_param, only: rho1
  !
  implicit none
  !
  private
  public  :: inte_vel
  !
  contains
  !
  subroutine inte_vel(n,nor,mflux,u,v,w,ug,vg,wg)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)                :: n
    real(8), intent(in   ), dimension( 0:, 0:, 0:, 1: ) :: nor
    real(8), intent(in   ), dimension( 0:, 0:, 0:)      :: mflux
    real(8), intent(in   ), dimension(-2:,-2:,-2:)      :: u,v,w
    real(8), intent(inout), dimension( 0:, 0:, 0:)      :: ug,vg,wg
    !
    integer :: i,j,k,ip,jp,kp
    real(8) :: mf_nx_drhoi_m, mf_nx_drhoi_p, mf_nx_drhoi, &
               mf_ny_drhoi_m, mf_ny_drhoi_p, mf_ny_drhoi, &
               mf_nz_drhoi_m, mf_nz_drhoi_p, mf_nz_drhoi
    !
    real(8), parameter :: rho1i = 1.d0/rho1
    !
    do k=1,n(3)
      kp = k+1
      do j=1,n(2)
        jp = j+1
        do i=1,n(1)
          ip = i+1
          !
          mf_nx_drhoi_m = mflux(i ,j,k)*nor(i ,j,k,1)*rho1i
          mf_nx_drhoi_p = mflux(ip,j,k)*nor(ip,j,k,1)*rho1i
          mf_nx_drhoi   = 0.5d0*(mf_nx_drhoi_m + mf_nx_drhoi_p)
          !
          mf_ny_drhoi_m = mflux(i,j ,k)*nor(i,j ,k,2)*rho1i
          mf_ny_drhoi_p = mflux(i,jp,k)*nor(i,jp,k,2)*rho1i
          mf_ny_drhoi   = 0.5d0*(mf_ny_drhoi_m + mf_ny_drhoi_p)
          !
          mf_nz_drhoi_m = mflux(i,j,k )*nor(i,j,k ,3)*rho1i
          mf_nz_drhoi_p = mflux(i,j,kp)*nor(i,j,kp,3)*rho1i
          mf_nz_drhoi   = 0.5d0*(mf_nz_drhoi_m + mf_nz_drhoi_p)
          !
          ug(i,j,k)     = u(i,j,k) - ug(i,j,k) + mf_nx_drhoi
          vg(i,j,k)     = v(i,j,k) - vg(i,j,k) + mf_ny_drhoi
          wg(i,j,k)     = w(i,j,k) - wg(i,j,k) + mf_nz_drhoi
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine inte_vel
  !
end module mod_inte_vel
