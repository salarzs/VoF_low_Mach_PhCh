module mod_fillps
  !
  implicit none
  !
  private
  public  :: fillps
  !
  contains
  !
  subroutine fillps(n,dli,dzci,dzfi,dti,vof,div_th,up,vp,wp,rho,pp,rho0,p)
    !
    !  fills the right-hand side of the Poisson equation for the correction step.
    !
    implicit none
    !
    integer, intent(in ), dimension(3)        :: n
    real(8), intent(in ), dimension(3)        :: dli
    real(8), intent(in ), dimension(-2:)      :: dzci,dzfi
    real(8), intent(in )                      :: dti
    real(8), intent(in ), dimension(0:,0:,0:) :: vof
    real(8), intent(in ), dimension(0:,0:,0:) :: div_th
    real(8), intent(in ), dimension(0:,0:,0:) :: up,vp,wp,rho,pp
    real(8), intent(in )                      :: rho0
    real(8), intent(out), dimension(0:,0:,0:) :: p
    !
    real(8) :: dtidxi,dtidyi!,dtidzi
    real(8), dimension(-2:n(3)+3) :: dtidzfi
    integer :: i,j,k,ip,jp,kp,im,jm,km
    real(8) :: rhoxm,rhoxp,rhoym,rhoyp,rhozm,rhozp
    real(8), dimension(8) :: mx,my,mz
    real(8) :: dvofdx,dvofdy,dvofdz
    !
    dtidxi = dti*dli(1)
    dtidyi = dti*dli(2)
    !dtidzi = dti*dli(3)
    dtidzfi(:) = dti*dzfi(:)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,up,vp,wp,rho,pp,dli,dzci,dzfi,dtidzfi,dtidyi,dtidxi,mflux,dti) &
    !$OMP SHARED(rhoxm,rhoym,rhozm,rhoxp,rhoyp,rhozp,dvofdx,dvofdy,dvofdz,mx,my,mz,vof) &
    !$OMP PRIVATE(i,j,k,im,jm,km,ip,jp,kp)
    do k=1,n(3)
      km = k-1
      kp = k+1
      do j=1,n(2)
        jm = j-1
        jp = j+1
        do i=1,n(1)
          im = i-1
          ip = i+1
          !
          p(i,j,k) = ( &
                      (wp(i,j,k)-wp(i,j,km))*dtidzfi(k)  + &
                      (vp(i,j,k)-vp(i,jm,k))*dtidyi      + &
                      (up(i,j,k)-up(im,j,k))*dtidxi    ) - &
                      dti*div_th(i,j,k)
          !
#ifdef FFT_DIRECT
          !
          p(i,j,k) = p(i,j,k)*rho0 
          !
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine fillps
  !
end module mod_fillps
