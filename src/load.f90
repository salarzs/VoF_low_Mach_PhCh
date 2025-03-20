module mod_load
  !
  use mpi
  use mod_common_mpi, only: ierr,dims,myid
  use decomp_2d
  use decomp_2d_io
  use mod_param, only: itot,jtot,ktot
  !
  implicit none
  !
  private
  public  :: load,load_scalar!load_ns,load_sf
  !
  contains
  !
  subroutine load(io,filename,n,fld)
    !
    ! reads/writes a restart file
    !
    implicit none
    !
    character(len=1), intent(in   )                            :: io
    character(len=*), intent(in   )                            :: filename
    integer         , intent(in   ), dimension(3)              :: n
    real(8)         , intent(inout), dimension(n(1),n(2),n(3)) :: fld ! generic field to be read/written
    !
    integer(MPI_OFFSET_KIND) :: filesize,disp,good
    integer :: lenr,fh
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      lenr = storage_size(fld(1,1,1))/8
      !good = int(ng(1)*ng(2)*ng(3),MPI_OFFSET_KIND)*lenr
      !good = int(lenr,MPI_OFFSET_KIND)*ng(1)*ng(2)*ng(3)
      good = int(lenr,MPI_OFFSET_KIND)*itot*jtot*ktot
      !
      if(filesize.ne.good) then
        if(myid.eq.0) print*, ''
        if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)
        call exit
      endif
      !
      ! read
      !
      if(myid.eq.0) print*, 'Armin', filename
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,3,fld )
      call MPI_FILE_CLOSE(fh,ierr)
      if(myid.eq.0) print*, 'Mehdi', filename
      ! 
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_write_var(fh,disp,3,fld)
      call MPI_FILE_CLOSE(fh,ierr)
      !
    end select
    !
    return
  end subroutine load
  !
  subroutine load_scalar(io,filename,pth,dpthdt_n,time,istep)
    !
    implicit none
    !
    character(len=1), intent(in   ) :: io
    character(len=*), intent(in   ) :: filename
    real(8)         , intent(inout) :: pth,dpthdt_n
    integer         , intent(inout) :: istep
    real(8)         , intent(inout) :: time
    !
    integer :: fh
    ! integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    ! integer(8), dimension(3) :: ng
    ! integer(8) :: lenr
    !
    select case(io)
    case('r')
      open(88,file=filename,status='old',action='read')
      read(88,*) pth,dpthdt_n,time,istep
      close(88)
    case('w')
      if(myid.eq.0) then
        open(88,file=filename)
        !write(88,'(3E15.7, 1I9.8)') pth,dpthdt_n,time,istep
        write(88,'(3E15.7, 1I12.11)') pth,dpthdt_n,time,istep
        close(88)
      endif
    end select
    !
    return
  end subroutine load_scalar
  !
  !subroutine load(io,filename,n,u,v,w,p,vof, &
  !                us,vs,ws,tmp,sca,time,istep,dt,pth,dpthdt_n)
  !  !
  !  ! reads/writes a restart file
  !  !
  !  implicit none
  !  !
  !  character(len=1), intent(in   )                            :: io
  !  character(len=*), intent(in   )                            :: filename
  !  integer         , intent(in   ), dimension(3)              :: n
  !  real(8)         , intent(inout), dimension(n(1),n(2),n(3)) :: u,v,w,p,vof, &
  !                                                                us,vs,ws,tmp,sca
  !  !
  !  real(8), intent(inout) :: time,istep,dt,pth,dpthdt_n
  !  real(8), dimension(5)  :: fldinfo
  !  integer :: fh
  !  integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
  !  integer(8), dimension(3) :: ng
  !  integer(8) :: lenr
  !  !
  !  select case(io)
  !  case('r')
  !    !
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
  !    !
  !    ! check file size first
  !    !
  !    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
  !    ng(:)   = n(:)
  !    ng(1:2) = ng(1:2)*dims(:)
  !    lenr = sizeof(time)
  !    good = (product(ng)*(3+1+1+3+2)+5)*lenr
  !    if(filesize.ne.good) then
  !      if(myid.eq.0) print*, ''
  !      if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
  !      if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
  !      call decomp_2d_finalize
  !      call MPI_FINALIZE(ierr)
  !      call exit
  !    endif
  !    !
  !    ! read
  !    !
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_read_var(fh,disp,3,u   )
  !    call decomp_2d_read_var(fh,disp,3,v   )
  !    call decomp_2d_read_var(fh,disp,3,w   )
  !    call decomp_2d_read_var(fh,disp,3,p   )
  !    call decomp_2d_read_var(fh,disp,3,vof )
  !    call decomp_2d_read_var(fh,disp,3,us  )
  !    call decomp_2d_read_var(fh,disp,3,vs  )
  !    call decomp_2d_read_var(fh,disp,3,ws  )
  !    call decomp_2d_read_var(fh,disp,3,tmp )
  !    call decomp_2d_read_var(fh,disp,3,sca )
  !    call decomp_2d_read_scalar(fh,disp,5,fldinfo)
  !    time     = fldinfo(1)
  !    istep    = fldinfo(2)
  !    dt       = fldinfo(3)
  !    pth      = fldinfo(4)
  !    dpthdt_n = fldinfo(5)
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  case('w')
  !    !
  !    ! write
  !    !
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
  !    filesize = 0_MPI_OFFSET_KIND
  !    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_write_var(fh,disp,3,u   )
  !    call decomp_2d_write_var(fh,disp,3,v   )
  !    call decomp_2d_write_var(fh,disp,3,w   )
  !    call decomp_2d_write_var(fh,disp,3,p   )
  !    call decomp_2d_write_var(fh,disp,3,vof )
  !    call decomp_2d_write_var(fh,disp,3,us  )
  !    call decomp_2d_write_var(fh,disp,3,vs  )
  !    call decomp_2d_write_var(fh,disp,3,ws  )
  !    call decomp_2d_write_var(fh,disp,3,tmp )
  !    call decomp_2d_write_var(fh,disp,3,sca )
  !    fldinfo = (/time,istep,dt,pth,dpthdt_n/)
  !    call decomp_2d_write_scalar(fh,disp,5,fldinfo)
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  end select
  !  !
  !  return
  !end subroutine load
  !!
  !subroutine load_ns(io,filename,n,u,v,w,p,time,istep,dt)
  !  !
  !  ! reads/writes a restart file
  !  !
  !  implicit none
  !  !
  !  character(len=1), intent(in   )                            :: io
  !  character(len=*), intent(in   )                            :: filename
  !  integer         , intent(in   ), dimension(3)              :: n
  !  real(8)         , intent(inout), dimension(n(1),n(2),n(3)) :: u,v,w,p
  !  real(8)         , intent(inout)                            :: time,istep,dt
  !  !                                                          
  !  real(8), dimension(3) :: fldinfo                           
  !  integer :: fh
  !  integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
  !  integer(8), dimension(3) :: ng
  !  integer(8) :: lenr
  !  !
  !  select case(io)
  !  case('r')
  !    !
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
  !    !
  !    ! check file size first
  !    !
  !    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
  !    ng(:)   = n(:)
  !    ng(1:2) = ng(1:2)*dims(:)
  !    lenr = sizeof(time)
  !    good = (product(ng)*(3+1)+3)*lenr
  !    if(filesize.ne.good) then
  !      if(myid.eq.0) print*, ''
  !      if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
  !      if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
  !      call decomp_2d_finalize
  !      call MPI_FINALIZE(ierr)
  !      call exit
  !    endif
  !    !
  !    ! read
  !    !
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_read_var(fh,disp,3,u)
  !    call decomp_2d_read_var(fh,disp,3,v)
  !    call decomp_2d_read_var(fh,disp,3,w)
  !    call decomp_2d_read_var(fh,disp,3,p)
  !    call decomp_2d_read_scalar(fh,disp,3,fldinfo)
  !    time  = fldinfo(1)
  !    istep = fldinfo(2)
  !    dt    = fldinfo(3)
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  case('w')
  !    !
  !    ! write
  !    !
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
  !    filesize = 0_MPI_OFFSET_KIND
  !    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_write_var(fh,disp,3,u)
  !    call decomp_2d_write_var(fh,disp,3,v)
  !    call decomp_2d_write_var(fh,disp,3,w)
  !    call decomp_2d_write_var(fh,disp,3,p)
  !    fldinfo = (/time,istep,dt/)
  !    call decomp_2d_write_scalar(fh,disp,3,fldinfo)
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  end select
  !  !
  !  return
  !end subroutine load_ns
  !!
  !subroutine load_sf(io,filename,n,p)
  !  !
  !  ! reads/writes a restart file
  !  !
  !  implicit none
  !  character(len=1)  , intent(in) :: io
  !  character(len=*), intent(in) :: filename
  !  integer, intent(in), dimension(3) :: n
  !  real(8), intent(inout), dimension(n(1),n(2),n(3)) :: p
  !  integer :: fh
  !  integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
  !  integer(8), dimension(3) :: ng
  !  integer(8) :: lenr
  !  !
  !  select case(io)
  !  case('r')
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
  !    !
  !    ! check file size first
  !    !
  !    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
  !    ng(:)   = n(:)
  !    ng(1:2) = ng(1:2)*dims(:)
  !    lenr = sizeof(p(1,1,1))
  !    good = product(ng)*lenr
  !    !
  !    if(filesize.ne.good) then
  !      if(myid.eq.0) print*, ''
  !      if(myid.eq.0) print*, '*** Simulation aborted due a checkpoint file with incorrect size ***'
  !      if(myid.eq.0) print*, '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
  !      call decomp_2d_finalize
  !      call MPI_FINALIZE(ierr)
  !      call exit
  !    endif
  !    !
  !    ! read
  !    !
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_read_var(fh,disp,3,p   )
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  case('w')
  !    !
  !    ! write
  !    !
  !    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
  !         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
  !    filesize = 0_MPI_OFFSET_KIND
  !    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
  !    disp = 0_MPI_OFFSET_KIND
  !    call decomp_2d_write_var(fh,disp,3,p   )
  !    call MPI_FILE_CLOSE(fh,ierr)
  !    !
  !  end select
  !  !
  !  return
  !end subroutine load_sf
  !
end module mod_load
