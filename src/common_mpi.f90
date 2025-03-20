module mod_common_mpi
  !
  use mpi
  use mod_param, only: dims
  !
  implicit none
  !
  integer               :: myid
  integer               :: left,right,front,back
  integer, dimension(2) :: coord
  integer               :: comm_cart,ierr
  integer               :: xhalo,yhalo,xhalo_3p,yhalo_3p
  integer               :: status(MPI_STATUS_SIZE)
  !
  integer(8), allocatable, dimension(:,:) :: xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf
  !
end module mod_common_mpi
