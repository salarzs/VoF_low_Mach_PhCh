module mod_vorticity
  !
  implicit none
  !
  private
  public vorticity
  !
  contains
  !
  subroutine vorticity(n,dli,ux,uy,uz,omega_x,omega_y,omega_z)
    !
    implicit none
    !
    integer, intent(in ), dimension(3) :: n
    real(8), intent(in ), dimension(3) :: dli
    real(8), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(8), intent(out), dimension(0:,0:,0:) :: omega_x,omega_y,omega_z
    !
    integer :: i,j,k
    !
    omega_x(1:n(1),1:n(2),1:n(3)) = 0.0
    omega_y(1:n(1),1:n(2),1:n(3)) = 0.0
    omega_z(1:n(1),1:n(2),1:n(3)) = 0.0
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          omega_x(i,j,k) = 0.25*( &
                                (uz(i,j+1,k  )-uz(i,j  ,k  ))*dli(2) - (uy(i,j  ,k+1)-uy(i,j  ,k  ))*dli(3) + &
                                (uz(i,j+1,k-1)-uz(i,j  ,k-1))*dli(2) - (uy(i,j  ,k  )-uy(i,j  ,k-1))*dli(3) + &
                                (uz(i,j  ,k  )-uz(i,j-1,k  ))*dli(2) - (uy(i,j-1,k+1)-uy(i,j-1,k  ))*dli(3) + &
                                (uz(i,j  ,k-1)-uz(i,j-1,k-1))*dli(2) - (uy(i,j-1,k  )-uy(i,j-1,k-1))*dli(3)   &
                                )
          !
#ifdef TWOD
          !
          omega_y(i,j,k) = 0.0
          omega_z(i,j,k) = 0.0
          !
#else
          !
          omega_y(i,j,k) = 0.25*( &
                                (ux(i  ,j,k+1)-ux(i  ,j,k  ))*dli(3) - (uz(i+1,j,k  )-uz(i  ,j,k  ))*dli(1) + &
                                (ux(i  ,j,k  )-ux(i  ,j,k-1))*dli(3) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))*dli(1) + &
                                (ux(i-1,j,k+1)-ux(i-1,j,k  ))*dli(3) - (uz(i  ,j,k  )-uz(i-1,j,k  ))*dli(1) + &
                                (ux(i-1,j,k  )-ux(i-1,j,k-1))*dli(3) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))*dli(1)   &
                                )
          !
          omega_z(i,j,k) = 0.25*( &
                                (uy(i+1,j  ,k)-uy(i  ,j  ,k))*dli(1) - (ux(i  ,j+1,k)-ux(i  ,j  ,k))*dli(2) + &
                                (uy(i+1,j-1,k)-uy(i  ,j-1,k))*dli(1) - (ux(i  ,j  ,k)-ux(i ,j-1 ,k))*dli(2) + &
                                (uy(i  ,j  ,k)-uy(i-1,j  ,k))*dli(1) - (ux(i-1,j+1,k)-ux(i-1,j  ,k))*dli(2) + &
                                (uy(i  ,j-1,k)-uy(i-1,j-1,k))*dli(1) - (ux(i-1,j  ,k)-ux(i-1,j-1,k))*dli(2)   &
                                )
          !
#endif     
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine vorticity
  !
end module mod_vorticity
