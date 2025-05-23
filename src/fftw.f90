module mod_fftw_param
  !
  use iso_c_binding
  !
  type, bind(C) :: fftw_iodim
     integer(C_INT) n, is, os
  end type fftw_iodim
  !
  interface
     type(C_PTR) function fftw_plan_guru_r2r(rank,dims, &
       howmany_rank,howmany_dims,in,out,kind,flags)  &
       bind(C, name='fftw_plan_guru_r2r')
       import
       integer(C_INT), value :: rank
       type(fftw_iodim), dimension(*), intent(in) :: dims
       integer(C_INT), value :: howmany_rank
       type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
       real(C_DOUBLE), dimension(*), intent(out) :: in,out
       integer(C_INT) :: kind
       integer(C_INT), value :: flags
     end function fftw_plan_guru_r2r
  end interface
  !
  integer :: FFTW_PATIENT,FFTW_ESTIMATE
  parameter (FFTW_PATIENT=32)
  parameter (FFTW_ESTIMATE=64)
  integer FFTW_R2HC
  parameter (FFTW_R2HC=0)
  integer FFTW_HC2R
  parameter (FFTW_HC2R=1)
  !
  integer REDFT00
  parameter (FFTW_REDFT00=3)
  integer FFTW_REDFT01
  parameter (FFTW_REDFT01=4)
  integer FFTW_REDFT10
  parameter (FFTW_REDFT10=5)
  integer FFTW_REDFT11
  parameter (FFTW_REDFT11=6)
  integer FFTW_RODFT00
  parameter (FFTW_RODFT00=7)
  integer FFTW_RODFT01
  parameter (FFTW_RODFT01=8)
  integer FFTW_RODFT10
  parameter (FFTW_RODFT10=9)
  integer FFTW_RODFT11
  parameter (FFTW_RODFT11=10)
  !
  type(C_PTR) :: fwd_guruplan_y,bwd_guruplan_y
  type(C_PTR) :: fwd_guruplan_z,bwd_guruplan_z
  logical :: planned=.false.
  !
end module mod_fftw_param
