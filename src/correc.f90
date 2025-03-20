module mod_correc
  !
  implicit none
  !
  public  :: correc
  private
  !
  contains
  !
  subroutine correc(n,qmin,dli,dzci,dt,p,pp,up,vp,wp,rho,rho0,u,v,w)
    !
    ! corrects the velocity to impose the constrain on the divergence
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                    :: n
    integer, intent(in )                                  :: qmin
    real(8), intent(in ), dimension(3)                    :: dli
    real(8), intent(in ), dimension(-2:)                  :: dzci
    real(8), intent(in )                                  :: dt
    real(8), intent(in ), dimension(    0:,    0:,    0:) :: p,pp,up,vp,wp
    real(8), intent(in ), dimension(    0:,    0:,    0:) :: rho
    real(8), intent(in )                                  :: rho0
    real(8), intent(out), dimension(-qmin:,-qmin:,-qmin:) :: u,v,w
    !
    real(8), dimension(-2:n(3)+3) :: factork
    real(8) :: factori,factorj
    integer :: i,j,k,ip,jp,kp
#ifdef FFT_DIRECT
    real(8) :: rho0i
#else
    real(8) :: rhox,rhoy,rhoz
#endif
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    factork = dt*dzci!dli(3)
#ifdef FFT_DIRECT
    rho0i = 1.d0/rho0
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,factori,factorj,factork,u,v,w,up,vp,wp,p,pp,rho) &
    !$OMP PRIVATE(i,j,k,ip,jp,kp,rhox,rhoy,rhoz)
    do k=1,n(3)
      kp = k+1
      do j=1,n(2)
        jp = j+1
        do i=1,n(1)
          ip = i+1
          !
#ifdef FFT_DIRECT
          !
          u(i,j,k) = up(i,j,k) - factori*(    (p( ip,j,k)-p( i,j,k))*rho0i )
          v(i,j,k) = vp(i,j,k) - factorj*(    (p( i,jp,k)-p( i,j,k))*rho0i )
          w(i,j,k) = wp(i,j,k) - factork(k)*( (p( i,j,kp)-p( i,j,k))*rho0i )
          !
#endif
#ifdef MULTI_GRID
          !
          rhox = 0.5d0*(rho(ip,j,k)+rho(i,j,k))
          rhoy = 0.5d0*(rho(i,jp,k)+rho(i,j,k))
          rhoz = 0.5d0*(rho(i,j,kp)+rho(i,j,k))
          !
          u(i,j,k) = up(i,j,k) - factori*(    p( ip,j,k)-p( i,j,k) )/rhox
          v(i,j,k) = vp(i,j,k) - factorj*(    p( i,jp,k)-p( i,j,k) )/rhoy
          w(i,j,k) = wp(i,j,k) - factork(k)*( p( i,j,kp)-p( i,j,k) )/rhoz
          !
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine correc
  !
end module mod_correc
