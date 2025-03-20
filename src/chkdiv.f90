module mod_chkdiv
  !
  use mpi
  use mod_common_mpi!, only: myid,coords,ierr
  !
  implicit none
  !
  private
  public chkdiv
  !
  contains
  !
  subroutine chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
  !subroutine chkdiv(n,dli,u,v,w,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension(-2:)         :: dzfi
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(out)                         :: divtot,divmax
    !
    real(8) :: dxi,dyi,dzi,div
    integer :: i,j,k,im,jm,km
    !integer :: ii,jj,kk
    !
    dxi = dli(1)
    dyi = dli(2)
    dzi = dli(3)
    divtot = 0.d0
    divmax = 0.d0
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP  SHARED(n,u,v,w,dxi,dyi,dzfi) &
    !$OMP  PRIVATE(i,j,k,im,jm,km,div) &
    !$OMP  REDUCTION(+:divtot) &
    !$OMP  REDUCTION(max:divmax)
    do k=1,n(3)
       km = k-1
       do j=1,n(2)
          jm = j-1
          do i=1,n(1)
             im = i-1
             div = (w(i,j,k)-w(i,j,km))*dzi + &
                   (v(i,j,k)-v(i,jm,k))*dyi + &
                   (u(i,j,k)-u(im,j,k))*dxi
             divmax = max(divmax,abs(div))
             divtot = divtot + div
    !         ii = coord(1)*n(1)+i
    !         jj = coord(2)*n(2)+j
    !         if(abs(div).ge.1.e+3 .AND. myid.eq.1) print*,'DIV',div,ii,jj,myid
         enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid.eq.0) print*, 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
    !
    return
  end subroutine chkdiv
  !
  !subroutine chkdiv_old(n,dli,u,v,w,div_max,div_tot)
  !  !
  !  implicit none
  !  !
  !  integer, intent(in ), dimension(3)        :: n
  !  real(8), intent(in ), dimension(3)        :: dli
  !  real(8), intent(in ), dimension(0:,0:,0:) :: u,v,w
  !  real(8), intent(out)                      :: div_max,div_tot
  !  !
  !  real(8) :: div,divtot,divtot_all,divmax(2),divmax_all(2)
  !  integer :: i,j,k,im,jm,km
  !  !
  !  divmax(1) = 0.0
  !  divmax(2) = 1.0*myid
  !  divtot    = 0.0 
  !  im = 0
  !  jm = 0
  !  km = 0
  !  !
  !  do k=1,n(3)
  !     do j=1,n(2)
  !        do i=1,n(1)
  !          !
  !          div = (w(i,j,k)-w(i,j,k-1))*dli(3) + &
  !                (v(i,j,k)-v(i,j-1,k))*dli(2) + &
  !                (u(i,j,k)-u(i-1,j,k))*dli(1)
  !          !
  !          divtot = divtot+div
  !          div = abs(div)
  !          !
  !          if(div.gt.divmax(1)) then
  !            divmax(1) = div
  !            im = i
  !            jm = j
  !            km = k
  !          endif
  !          !
  !        enddo
  !     enddo
  !  enddo
  !  !
  !  call mpi_allreduce(divtot,divtot_all,1,mpi_real8,mpi_sum,MPI_COMM_WORLD,ierr)
  !  call mpi_allreduce(divmax,divmax_all,1,mpi_2double_precision,mpi_maxloc,MPI_COMM_WORLD,ierr)
  !  !
  !  if(myid.eq.int(divmax_all(2))) then
  !     write(6,111) zstart(1)-1+im, zstart(2)-1+jm,km
  !     write(6,222) divtot_all,divmax_all(1),int(divmax_all(2))
  !  111    format('Maximal divergence at i = ',I5,' j = ', I5,' k = ',I5)
  !  222    format('Divergence: Tot = ',e13.6,' Max = ',e13.6,' Rank = ',I3)
  !  endif
  !  !
  !  div_tot = divtot_all
  !  div_max = divmax_all(1)
  !  !
  !  return
  !end subroutine chkdiv_old
  ! 
end module mod_chkdiv
