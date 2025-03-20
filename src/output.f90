module mod_output
  !
  use mpi
  use decomp_2d_io
  use mod_param     , only: dims,dx,dy,dz,dxi,dyi,dzi,lx,ly,lz,datapos
  use mod_common_mpi, only: ierr,myid,coord,comm_cart
  use mod_bound     , only: boundp
  !
  implicit none
  !
  private
  public  :: out0d,out1d,out1d_2,out2d,out3d,velstats,cmpt_total_area, &
             cmpt_infq,int_qtn,time_avg,mixed_variables,compute_vorticity, &
             cmpt_delta,energy_th,density_cont,budget, &
             energy_balance
  !
  contains
  !
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n
    real(8), intent(in), dimension(:) :: var
    integer :: iunit
    character(len=30) :: cfmt
    integer :: i
    !
    write(cfmt,'(A,I3,A)') '(',n,'E15.7)'
    iunit = 10
    if (myid .eq. 0) then
      open(iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) (var(i),i=1,n) 
      close(iunit)
    endif
    return
  end subroutine out0d
  !
  subroutine out1d(fname,n,qmin,idir,z,dzlzi,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! n     -> size of the input array
    ! qmin  -> starting index of the array (e.g., -2: vel, 0: p)
    ! idir  -> direction of the profile
    ! z     -> z coordinate (grid is non-uniform in z)
    ! dzlzi -> dz/lz weight of a grid cell for averaging over z
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: qmin
    integer, intent(in) :: idir
    real(8), intent(in), dimension(-2:) :: z,dzlzi
    real(8), intent(in), dimension(-qmin:,-qmin:,-qmin:) :: p
    real(8), allocatable, dimension(:) :: p1d
    integer :: i,j,k,ii,jj
    integer :: iunit
    integer, dimension(3) :: ng
    !
    ng(:) = n(:)
    ng(1:2) = n(1:2)*dims(1:2)
    iunit = 10
    select case(idir)
    case(3)
      allocate(p1d(n(3)))
      do k=1,ng(3)
        p1d(k) = 0.
        do j=1,n(2)
          do i=1,n(1)
            p1d(k) = p1d(k) + p(i,j,k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1)*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do k=1,n(3)
          write(iunit,'(2E15.7)') z(k),p1d(k)
        enddo
        close(iunit)
      endif
    case(2)
      allocate(p1d(ng(2)))
      p1d(:) = 0.
      do j=1,n(2)
        jj = coord(2)*n(2)+j
        p1d(jj) = 0.
        do k=1,n(3)
          do i=1,n(1)
            p1d(jj) = p1d(jj) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(2),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(1))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do j=1,ng(2)
          write(iunit,'(2E15.7)') (1.d0*j-.5d0)/(1.d0*ng(2)),p1d(j)
        enddo
        close(iunit)
      endif
    case(1)
      allocate(p1d(ng(1)))
      p1d(:) = 0.
      do i=1,n(1)
        ii = coord(1)*n(1)+i
        p1d(i) = 0.
        do k=1,n(3)
          do j=1,n(2)
            p1d(ii) = p1d(ii) + p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,p1d(1),ng(1),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d(:) = p1d(:)/(1.*ng(2))
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,'(2E15.7)') (1.d0*i-.5d0)/(1.d0*n(1)),p1d(j)
        enddo
        close(iunit)
      endif
    end select
    deallocate(p1d)
  end subroutine out1d
  !
  subroutine velstats(fname,n,psi,u,v,w,p)
    !
    ! computes some volume-averaged profiles along z (assumes constant grid spacing)
    !
    implicit none
    real(8), parameter :: small = 1.e-9
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension( 0:, 0:, 0:) :: psi,p
    real(8), intent(in), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), dimension(n(3)) :: u_1,u_2,v_1,v_2,w_1,w_2,p_1,p_2, &
                                uu_1,uu_2,vv_1,vv_2,ww_1,ww_2,pp_1,pp_2,ph
    integer :: i,j,k,ii,jj
    integer :: iunit
    integer, dimension(3) :: ng
    !
    ng(:) = n(:)
    ng(1:2) = n(1:2)*dims(1:2)
    do k=1,n(3)
      u_1(k)  = 0. 
      u_2(k)  = 0. 
      v_1(k)  = 0. 
      v_2(k)  = 0. 
      w_1(k)  = 0. 
      w_2(k)  = 0. 
      p_1(k)  = 0. 
      p_2(k)  = 0. 
      uu_1(k) = 0. 
      uu_2(k) = 0. 
      vv_1(k) = 0. 
      vv_2(k) = 0. 
      ww_1(k) = 0. 
      ww_2(k) = 0. 
      pp_1(k) = 0. 
      pp_2(k) = 0. 
      ph(k)   = 0. 
      do j=1,n(2)
        do i=1,n(1)
          u_1(k)  = u_1(k)  + (   psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))
          u_2(k)  = u_2(k)  + (1.-psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))
          v_1(k)  = v_1(k)  + (   psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))
          v_2(k)  = v_2(k)  + (1.-psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))
          w_1(k)  = w_1(k)  + (   psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))
          w_2(k)  = w_2(k)  + (1.-psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))
          p_1(k)  = p_1(k)  + (   psi(i,j,k))*(0.5*(p(i,j,k)           ))
          p_2(k)  = p_2(k)  + (1.-psi(i,j,k))*(0.5*(p(i,j,k)           ))
          uu_1(k) = uu_1(k) + (   psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))**2
          uu_2(k) = uu_2(k) + (1.-psi(i,j,k))*(0.5*(u(i,j,k)+u(i-1,j,k)))**2
          vv_1(k) = vv_1(k) + (   psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))**2
          vv_2(k) = vv_2(k) + (1.-psi(i,j,k))*(0.5*(v(i,j,k)+v(i,j-1,k)))**2
          ww_1(k) = ww_1(k) + (   psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))**2
          ww_2(k) = ww_2(k) + (1.-psi(i,j,k))*(0.5*(w(i,j,k)+w(i,j,k-1)))**2
          pp_1(k) = pp_1(k) + (   psi(i,j,k))*(0.5*(p(i,j,k)           ))**2
          pp_2(k) = pp_2(k) + (1.-psi(i,j,k))*(0.5*(p(i,j,k)           ))**2
          ph(k)   = ph(k)   + psi(i,j,k)
        enddo
      enddo
    enddo
    call mpi_allreduce(MPI_IN_PLACE,u_1(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,u_2(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v_1(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v_2(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w_1(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w_2(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p_1(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,p_2(1) ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,uu_1(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,uu_2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vv_1(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vv_2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ww_1(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ww_2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,pp_1(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,pp_2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ph(1)  ,n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    do k=1,n(3)
      ph(k)   = ph(k)/(1.*ng(1)*ng(2))
      u_1(k)  = u_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small)
      u_2(k)  = u_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small)
      v_1(k)  = v_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small)
      v_2(k)  = v_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small)
      w_1(k)  = w_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small)
      w_2(k)  = w_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small)
      p_1(k)  = p_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small)
      p_2(k)  = p_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small)
      uu_1(k) = uu_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small) - u_1(k)**2
      uu_2(k) = uu_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small) - u_2(k)**2
      vv_1(k) = vv_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small) - v_1(k)**2
      vv_2(k) = vv_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small) - v_2(k)**2
      ww_1(k) = ww_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small) - w_1(k)**2
      ww_2(k) = ww_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small) - w_2(k)**2
      pp_1(k) = pp_1(k)/(1.*ng(1)*ng(2))/(ph(k)   +small) - p_1(k)**2
      pp_2(k) = pp_2(k)/(1.*ng(1)*ng(2))/(1.-ph(k)+small) - p_2(k)**2
    enddo
    iunit = 10
    if(myid.eq.0) then
      open(unit=iunit,file=fname)
      do k=1,n(3)
        write(iunit,'(18E15.7)') (k-0.5)*dz,ph(k), &
                                            u_1(k) ,v_1(k) ,w_1(k) ,p_1(k) , &
                                            uu_1(k),vv_1(k),ww_1(k),pp_1(k), &
                                            u_2(k) ,v_2(k) ,w_2(k) ,p_2(k) , &
                                            uu_2(k),vv_2(k),ww_2(k),pp_2(k)
      enddo
      close(iunit)
    endif
    return
  end subroutine velstats
  subroutine cmpt_total_area(n,psi,area)
    !
    ! computes some volume-averaged profiles along z (assumes constant grid spacing)
    !
    implicit none
    integer, intent(in), dimension(3) :: n
    real(8), intent(in), dimension(0:,0:,0:) :: psi
    real, intent(out) :: area
    real(8), dimension(8) :: mx,my,mz
    real(8) :: norm,dpsidx,dpsidy,dpsidz
    integer :: i,j,k,ii,jj
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !i+1/2 j+1/2 k+1/2
          mx(1)=((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j+1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1)))*dxi*0.25
          !i+1/2 j-1/2 k+1/2
          mx(2)=((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k+1)+psi(i+1,j-1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1)))*dxi*0.25
          !i+1/2 j+1/2 k-1/2
          mx(3)=((psi(i+1,j  ,k  )+psi(i+1,j+1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j+1,k-1))-&
                 (psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1)))*dxi*0.25
          !i+1/2 j-1/2 k-1/2
          mx(4)=((psi(i+1,j  ,k  )+psi(i+1,j-1,k  )+psi(i+1,j  ,k-1)+psi(i+1,j-1,k-1))-&
                 (psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1)))*dxi*0.25
          !i-1/2 j+1/2 k+1/2
          mx(5)=((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j+1,k+1))-&
                 (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j+1,k+1)))*dxi*0.25
          !i-1/2 j-1/2 k+1/2
          mx(6)=((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k+1)+psi(i  ,j-1,k+1))-&
                 (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k+1)+psi(i-1,j-1,k+1)))*dxi*0.25
          !i-1/2 j+1/2 k-1/2
          mx(7)=((psi(i  ,j  ,k  )+psi(i  ,j+1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j+1,k-1))-&
                 (psi(i-1,j  ,k  )+psi(i-1,j+1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j+1,k-1)))*dxi*0.25
          !i-1/2 j-1/2 k-1/2
          mx(8)=((psi(i  ,j  ,k  )+psi(i  ,j-1,k  )+psi(i  ,j  ,k-1)+psi(i  ,j-1,k-1))-&
                 (psi(i-1,j  ,k  )+psi(i-1,j-1,k  )+psi(i-1,j  ,k-1)+psi(i-1,j-1,k-1)))*dxi*0.25
          ! 
          !i+1/2 j+1/2 k+1/2
          my(1)=((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)))*dyi*0.25
          !i+1/2 j-1/2 k+1/2
          my(2)=((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1))-&
                 (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1)))*dyi*0.25
          !i+1/2 j+1/2 k-1/2
          my(3)=((psi(i  ,j+1,k  )+psi(i+1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1))-&
                 (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)))*dyi*0.25
          !i+1/2 j-1/2 k-1/2
          my(4)=((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1))-&
                 (psi(i  ,j-1,k  )+psi(i+1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dyi*0.25
          !i-1/2 j+1/2 k+1/2
          my(5)=((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)))*dyi*0.25
          !i-1/2 j-1/2 k+1/2
          my(6)=((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1))-&
                 (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1)))*dyi*0.25
          !i-1/2 j+1/2 k-1/2
          my(7)=((psi(i  ,j+1,k  )+psi(i-1,j+1,k  )+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1))-&
                 (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)))*dyi*0.25
          !i-1/2 j-1/2 k-1/2
          my(8)=((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1))-&
                 (psi(i  ,j-1,k  )+psi(i-1,j-1,k  )+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dyi*0.25
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i+1,j+1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  )))*dzi*0.25
          !i+1/2 j-1/2 k+1/2
          mz(2)=((psi(i  ,j  ,k+1)+psi(i+1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i+1,j-1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  )))*dzi*0.25
          !i+1/2 j+1/2 k-1/2
          mz(3)=((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j+1,k  )+psi(i+1,j+1,k  ))-&
                 (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i+1,j+1,k-1)))*dzi*0.25
          !i+1/2 j-1/2 k-1/2
          mz(4)=((psi(i  ,j  ,k  )+psi(i+1,j  ,k  )+psi(i  ,j-1,k  )+psi(i+1,j-1,k  ))-&
                 (psi(i  ,j  ,k-1)+psi(i+1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i+1,j-1,k-1)))*dzi*0.25
          !i-1/2 j+1/2 k+1/2
          mz(5)=((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j+1,k+1)+psi(i-1,j+1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  )))*dzi*0.25
          !i-1/2 j-1/2 k+1/2
          mz(6)=((psi(i  ,j  ,k+1)+psi(i-1,j  ,k+1)+psi(i  ,j-1,k+1)+psi(i-1,j-1,k+1))-&
                 (psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  )))*dzi*0.25
          !i-1/2 j+1/2 k-1/2
          mz(7)=((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j+1,k  )+psi(i-1,j+1,k  ))-&
                 (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j+1,k-1)+psi(i-1,j+1,k-1)))*dzi*0.25
          !i-1/2 j-1/2 k-1/2
          mz(8)=((psi(i  ,j  ,k  )+psi(i-1,j  ,k  )+psi(i  ,j-1,k  )+psi(i-1,j-1,k  ))-&
                 (psi(i  ,j  ,k-1)+psi(i-1,j  ,k-1)+psi(i  ,j-1,k-1)+psi(i-1,j-1,k-1)))*dzi*0.25
          !
          dpsidx=0.125*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dpsidy=0.125*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dpsidz=0.125*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          norm=SQRT(dpsidx**2+dpsidy**2+dpsidz**2)
          area = area + norm*dx*dy*dz
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,area ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    return
  end subroutine cmpt_total_area
  !
  subroutine out2d(fname,inorm,islice,p)
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice 
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in) :: inorm,islice
    real(8),intent(in), dimension(:,:,:) :: p
    !
    select case(inorm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(3,p,inorm,islice,fname)
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(3,p,inorm,islice,fname)
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(3,p,inorm,islice,fname)
    end select
    return
  end subroutine out2d
  !
  subroutine out3d(fname,nskip,p)
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: (/1,1,1/) 
    !           writes the full field 
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: nskip
    real(8),intent(in), dimension(:,:,:) :: p
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_every(3,p,nskip(1),nskip(2),nskip(3),fname,.true.)
    call MPI_FILE_CLOSE(fh,ierr)
    return
  end subroutine out3d
  !
  subroutine out1d_2(fname,n,idir,z,u,v,w) ! e.g. for a channel with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(8), intent(in), dimension(-2:) :: z
    real(8), intent(in), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    integer :: i,j,k
    integer :: iunit
    integer, dimension(3) :: ng
    integer :: q
    !
    ng(:) = n(:)
    ng(1:2) = n(1:2)*dims(1:2)
    iunit = 10
    select case(idir)
    case(3)
      q = n(3)
      allocate(um(0:q+1),vm(0:q+1),wm(0:q+1),u2(0:q+1),v2(0:q+1),w2(0:q+1),uw(0:q+1))
      do k=1,n(3)
        um(k) = 0.
        vm(k) = 0.
        wm(k) = 0.
        u2(k) = 0.
        v2(k) = 0.
        w2(k) = 0.
        uw(k) = 0.
        do j=1,n(2)
          do i=1,n(1)
            um(k) = um(k) + u(i,j,k)
            vm(k) = vm(k) + v(i,j,k)
            wm(k) = wm(k) + 0.50d0*(w(i,j,k-1) + w(i,j,k))
            u2(k) = u2(k) + u(i,j,k)**2
            v2(k) = v2(k) + v(i,j,k)**2
            w2(k) = w2(k) + 0.50d0*(w(i,j,k)**2+w(i,j,k-1)**2)
            uw(k) = uw(k) + 0.25d0*(u(i-1,j,k) + u(i,j,k))* &
                                   (w(i,j,k-1) + w(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uw(1),n(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:) = um(:)/(1.*ng(1)*ng(2))
      vm(:) = vm(:)/(1.*ng(1)*ng(2))
      wm(:) = wm(:)/(1.*ng(1)*ng(2))
      u2(:) = sqrt(u2(:)/(1.*ng(1)*ng(2)) - um(:)**2)
      v2(:) = sqrt(v2(:)/(1.*ng(1)*ng(2)) - vm(:)**2)
      w2(:) = sqrt(w2(:)/(1.*ng(1)*ng(2)) - wm(:)**2)
      uw(:) = uw(:)/(1.*ng(1)*ng(2)) - um(:)*wm(:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do k=1,n(3)
          write(iunit,'(8E15.7)') z(k),um(k),vm(k),wm(k), &
                                       u2(k),v2(k),w2(k), &
                                       uw(k)
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,uw)
    case(2)
    case(1)
    end select
  end subroutine out1d_2
  !
  subroutine out2d_2(fname,n,idir,z,u,v,w) ! e.g. for a duct with streamwise dir in x
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(8), intent(in), dimension(-2:) :: z
    real(8), intent(in), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,vw
    integer :: i,j,k,ii,jj,kk
    integer :: iunit
    integer, dimension(3) :: ng
    integer :: p,q
    real(8) :: y
    !
    ng(:) = n(:)
    ng(1:2) = n(1:2)*dims(1:2)
    iunit = 10
    select case(idir)
    case(3)
      p = ng(1)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),vw(p,q))
      !
      um(:,:) = 0.d0
      vm(:,:) = 0.d0
      wm(:,:) = 0.d0
      u2(:,:) = 0.d0
      v2(:,:) = 0.d0
      w2(:,:) = 0.d0
      uv(:,:) = 0.d0
      vw(:,:) = 0.d0
      do k=1,n(3)
        kk = k
        do i=1,n(1)
          ii = i+coord(1)*n(1)
          um(ii,kk) = 0.d0
          vm(ii,kk) = 0.d0
          wm(ii,kk) = 0.d0
          u2(ii,kk) = 0.d0
          v2(ii,kk) = 0.d0
          w2(ii,kk) = 0.d0
          vw(ii,kk) = 0.d0
          uv(ii,kk) = 0.d0
          do j=1,n(2)
            jj = j+coord(2)*n(2)
            um(ii,kk) = um(ii,kk) + 0.5d0*(u(i-1,j,k)+u(i,j,k))
            vm(ii,kk) = vm(ii,kk) + v(i,j,k)
            wm(ii,kk) = wm(ii,kk) + 0.5d0*(w(i,j,k-1)+w(i,j,k))
            u2(ii,kk) = u2(ii,kk) + 0.5d0*(u(i-1,j,k)**2+u(i,j,k)**2)
            v2(ii,kk) = v2(ii,kk) + v(i,j,k)**2
            w2(ii,kk) = w2(ii,kk) + 0.5d0*(w(i,j,k-1)**2+w(i,j,k)**2)
            vw(ii,kk) = vw(ii,kk) + 0.25d0*(v(i,j-1,k) + v(i,j,k))* &
                                           (w(i,j,k-1) + w(i,j,k))
            uv(ii,kk) = uv(ii,kk) + 0.25d0*(u(i-1,j,k) + u(i,j,k))* &
                                           (v(i,j-1,k) + v(i,j,k))
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,um(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vm(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,wm(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,u2(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,v2(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,w2(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,vw(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,uv(1,1),ng(1)*ng(3),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      um(:,:) =      um(:,:)/(1.*ng(2))
      vm(:,:) =      vm(:,:)/(1.*ng(2))
      wm(:,:) =      wm(:,:)/(1.*ng(2))
      u2(:,:) = sqrt(u2(:,:)/(1.*ng(2)) - um(:,:)**2)
      v2(:,:) = sqrt(v2(:,:)/(1.*ng(2)) - vm(:,:)**2)
      w2(:,:) = sqrt(w2(:,:)/(1.*ng(2)) - wm(:,:)**2)
      vw(:,:) =      vw(:,:)/(1.*ng(2)) - vm(:,:)*wm(:,:)
      uv(:,:) =      uv(:,:)/(1.*ng(2)) - um(:,:)*vm(:,:)
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do k=1,ng(3)
          do i=1,ng(1)
            y = (i-.5d0)*dx
            write(iunit,'(10E15.7)') y,z(k),um(i,k),vm(i,k),wm(i,k), &
                                            u2(i,k),v2(i,k),w2(i,k), &
                                            vw(i,k),uv(i,k)
          enddo
        enddo
        close(iunit)
      endif
      deallocate(um,vm,wm,u2,v2,w2,vw,uv)
    case(2)
    case(1)
    end select
  end subroutine out2d_2
  !
  subroutine cmpt_infq(n,dl,num,rdir,phi,sca,tmp,tinf_v,yinf_v)
    !
    use mod_common_mpi, only: mpi_sum,ierr,comm_cart,mpi_real8,myid
    use mod_param,      only: tmpl0,sinit
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dl
    integer, intent(in )                         :: num
    integer, intent(in )                         :: rdir
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: phi
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: sca
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(out), dimension(num+1)       :: tinf_v,yinf_v
    !
    real(8), dimension(-2:n(1)+3,-2:n(2)+3,-2:n(3)+3,num) :: phi_v
    integer, dimension(1:num) :: coun
    real(8), dimension(1:num) :: rcoun
    real(8) :: tmpg_i,yg_i,d
    integer :: i,j,k,ip,jp,kp,im,jm,km,p
    !
    ! 1. first compute the iso-surfaces away from the interface
    !    in the gas phase
    !
    do k=-2,n(3)+3
      do j=-2,n(2)+3
        do i=-2,n(1)+3
          !
          do p=1,num
            phi_v(i,j,k,p) = phi(i,j,k) + p*rdir*dl(2)
          enddo
          !
        enddo
      enddo
    enddo
    !
    ! 2. compute the values on the iso-surfaces 
    !
    tinf_v(1:num) = 0.d0
    yinf_v(1:num) = 0.d0
    coun(1:num)   = 0
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
          if(phi(i,j,k).lt.0.d0) then
            !
            do p=1,num
              !
              ! along x
              !
              if( phi_v(ip,j,k,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(ip,j,k,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(ip,j,k,p))+tmp(ip,j,k)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(ip,j,k,p))+sca(ip,j,k)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              if( phi_v(im,j,k,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(im,j,k,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(im,j,k,p))+tmp(im,j,k)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(im,j,k,p))+sca(im,j,k)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              !
              ! along y
              !
              if( phi_v(i,jp,k,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(i,jp,k,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,jp,k,p))+tmp(i,jp,k)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(i,jp,k,p))+sca(i,jp,k)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              if( phi_v(i,jm,k,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(i,jm,k,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,jm,k,p))+tmp(i,jm,k)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(i,jm,k,p))+sca(i,jm,k)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              !
              ! along z
              !
              if( phi_v(i,j,kp,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(i,j,kp,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,j,kp,p))+tmp(i,j,kp)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(i,j,kp,p))+sca(i,j,kp)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              if( phi_v(i,j,km,p)*phi_v(i,j,k,p).lt.0. ) then
                d         = (abs(phi_v(i,j,km,p))+abs(phi_v(i,j,k,p)))
                tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,j,km,p))+tmp(i,j,km)*abs(phi_v(i,j,k,p)))/d
                yg_i      = (sca(i,j,k)*abs(phi_v(i,j,km,p))+sca(i,j,km)*abs(phi_v(i,j,k,p)))/d
                tinf_v(p) = tinf_v(p) + tmpg_i
                yinf_v(p) = yinf_v(p) + yg_i
                coun(p)   = coun(p) + 1
              endif
              !
              rcoun(p) = 1.d0*coun(p)
              !
            enddo
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,rcoun ,num,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tinf_v,num,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,yinf_v,num,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    do p=1,num
      !
      if(rcoun(p).eq.0.d0) then
        tinf_v(p) = tmpl0
        yinf_v(p) = 0.d0
      else
        tinf_v(p) = tinf_v(p)/rcoun(p)
        yinf_v(p) = yinf_v(p)/rcoun(p)
      endif
      !
    enddo
    tinf_v(num+1) = tmpl0
    yinf_v(num+1) = sinit
    !
    return
  end subroutine cmpt_infq
  !
  subroutine int_qtn(n,rho_gas,kappa,tmp,tmpge,tmple,vof,mflux,sca,scae,phi,num,tinf_v,yinf_v,pth,ptho,dpthdt,time,mas_g,istep)
    !
    ! quantity at the interface (heat and mass transfer)
    !
    use mod_param , only: dl,dli,lx,ly,lz,datapos,lheat,rho1, &
                          d_m12,rho2_0,kappa2,itot,jtot,ktot,lref, &
                          tmpl0,tmpg0
    !
    use mod_thermo, only: thermo_rhog,thermo_mug,thermo_kag
    use mod_thermo, only: thermo_d_lg
    use mod_thermo, only: mass_fraction
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rho_gas,kappa
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: tmpge
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: tmple
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: vof,mflux
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: sca,scae
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: phi
    integer, intent(in )                         :: num
    real(8), intent(in ), dimension(num+1)       :: tinf_v,yinf_v
    real(8), intent(in )                         :: pth,ptho,dpthdt
    real(8), intent(in )                         :: time,mas_g
    integer, intent(in )                         :: istep
    !
    real(8), dimension(num+1) :: sh_vp,sh_cp,nu_sen,nu_lat
    real(8), dimension(8) :: mx,my,mz
    real(8) :: nu_cp,nu_vp
    integer :: i,j,k, &
               ip,jp,kp, &
               im,jm,km, int_s
    real(8) :: area_d,vol_d,dvofdx,dvofdy,dvofdz,dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm, &
               lap_kap,kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm,vofxp,vofxm,vofyp,vofym, &
               vofzp,vofzm,mfx_tot,qka_tot_sen,qka_tot_lat,norm,tmp_d, &
               kappag_in_t,rhog_in_t,d_lg_in_t,s_in_t,tmpl_in_t,tmpg_in_t,h_mas,tmp_inf,s_inf,h_tmp_lat,h_tmp_sen
    real(8) :: y_g,tmp_g,tmpg_i,tmpl_i,s_in,rhog_av,theta,mug_av,rcounter_t,dvol,volll
    integer :: counter_t
    real(8) :: rhog_in,rhol_in,t_in,ss_in,rhog_ap_in,d_lg_ap_in,dtmpdx,dtmpdy,dtmpdz
    ! 
    volll = 1.d0*itot*jtot*ktot
    dvol  = product(dl(:))
    !
    ! quantities for evaporation (tot.)
    !
    vol_d       = 0.d0
    area_d      = 0.d0
    mfx_tot     = 0.d0
    tmp_d       = 0.d0
    tmp_g       = 0.d0
    y_g         = 0.d0
    qka_tot_sen = 0.d0
    qka_tot_lat = 0.d0
    counter_t   = 0
    rcounter_t  = 0.d0
    rhog_in_t   = 0.d0
    d_lg_in_t   = 0.d0
    kappag_in_t = 0.d0
    s_in_t      = 0.d0
    tmpl_in_t   = 0.d0
    tmpg_in_t   = 0.d0
    ! 
    do k=1,n(3)
      kp = k + 1
      km = k - 1
      do j=1,n(2)
        jp = j + 1
        jm = j - 1
        do i=1,n(1)
          ip = i + 1
          im = i - 1
          !
#ifdef LOW_MACH
          !
          ss_in = scae(i,j,k)
          t_in  = tmpge(i,j,k) ! for now
          !
          rhog_in = thermo_rhog(ptho,t_in,ss_in) ! to be consistent with rho_gas after VoF
          rhol_in = rho1
          !
#else
          !
          ss_in = 0.d0
          t_in  = tmp(i,j,k)
          !
          rhog_in = rho_gas(i,j,k) 
          rhol_in = rho1
          !
#endif
          ! 
          ! a1: total liquid volume
          ! 
          vol_d = vol_d+vof(i,j,k)*dl(1)*dl(2)*dl(3) ! total liquid volume  
          ! 
          ! a2: total liquid area
          ! 
          !i+1/2 j+1/2 k+1/2
          mx(1)=((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)*0.25
          !i+1/2 j-1/2 k+1/2
          mx(2)=((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)*0.25
          !i+1/2 j+1/2 k-1/2
          mx(3)=((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)*0.25
          !i+1/2 j-1/2 k-1/2
          mx(4)=((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)*0.25
          !i-1/2 j+1/2 k+1/2
          mx(5)=((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)*0.25
          !i-1/2 j-1/2 k+1/2
          mx(6)=((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25
          !i-1/2 j+1/2 k-1/2
          mx(7)=((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)*0.25
          !i-1/2 j-1/2 k-1/2
          mx(8)=((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)*0.25
          ! 
          !i+1/2 j+1/2 k+1/2
          my(1)=((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)*0.25
          !i+1/2 j-1/2 k+1/2
          my(2)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1))-&
                 (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)*0.25
          !i+1/2 j+1/2 k-1/2
          my(3)=((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)*0.25
          !i+1/2 j-1/2 k-1/2
          my(4)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1))-&
                 (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)*0.25
          !i-1/2 j+1/2 k+1/2
          my(5)=((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)*0.25
          !i-1/2 j-1/2 k+1/2
          my(6)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1))-&
                 (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(2)*0.25
          !i-1/2 j+1/2 k-1/2
          my(7)=((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)*0.25
          !i-1/2 j-1/2 k-1/2
          my(8)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1))-&
                 (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)*0.25
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli(3)*0.25
          !i+1/2 j-1/2 k+1/2
          mz(2)=((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli(3)*0.25
          !i+1/2 j+1/2 k-1/2
          mz(3)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli(3)*0.25
          !i+1/2 j-1/2 k-1/2
          mz(4)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(3)*0.25
          !i-1/2 j+1/2 k+1/2
          mz(5)=((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli(3)*0.25
          !i-1/2 j-1/2 k+1/2
          mz(6)=((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli(3)*0.25
          !i-1/2 j+1/2 k-1/2
          mz(7)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli(3)*0.25
          !i-1/2 j-1/2 k-1/2
          mz(8)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(3)*0.25
          !
          dvofdx = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dvofdy = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dvofdz = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          norm   = sqrt(dvofdx**2+dvofdy**2+dvofdz**2)
          !
          area_d = area_d + norm*dl(1)*dl(2)*dl(3) ! total liquid area
          ! 
          ! a3: total exchanged mass flux (kg/s)
          ! 
          mfx_tot = mfx_tot + mflux(i,j,k)*norm*dl(1)*dl(2)*dl(3)
          ! 
          ! a4: liquid temperature
          ! 
          tmp_d = tmp_d + vof(i,j,k)*tmp(i,j,k)*dl(1)*dl(2)*dl(3) ! droplet temperature
          ! 
          ! a5: gas temperature
          ! 
          tmp_g = tmp_g + (1.d0-vof(i,j,k))*tmp(i,j,k)*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k) ! gas temperature
          ! 
          ! a6: gas vapor mass
          ! 
          y_g = y_g + (1.d0-vof(i,j,k))*sca(i,j,k)*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k) ! gas vapor mass
          ! 
          ! a7: q tot
          ! 
          kappaxp = 0.5d0*(kappa(ip,j,k)+kappa(i,j,k))
          kappaxm = 0.5d0*(kappa(im,j,k)+kappa(i,j,k))
          kappayp = 0.5d0*(kappa(i,jp,k)+kappa(i,j,k))
          kappaym = 0.5d0*(kappa(i,jm,k)+kappa(i,j,k))
          kappazp = 0.5d0*(kappa(i,j,kp)+kappa(i,j,k))
          kappazm = 0.5d0*(kappa(i,j,km)+kappa(i,j,k))
          !
          dtmpdxp = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
          dtmpdxm = (tmp(i ,j,k)-tmp(im,j,k))*dli(1)
          dtmpdyp = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
          dtmpdym = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
          dtmpdzp = (tmp(i,j,kp)-tmp(i,j,k ))*dli(3)
          dtmpdzm = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
          !compute d/dx_j(k*dT/dx_j)
          lap_kap = (kappaxp*dtmpdxp-kappaxm*dtmpdxm)*dli(1) + &
                    (kappayp*dtmpdyp-kappaym*dtmpdym)*dli(2) + &
                    (kappazp*dtmpdzp-kappazm*dtmpdzm)*dli(3)
          !
          qka_tot_sen = qka_tot_sen + lap_kap*vof(i,j,k)*dl(1)*dl(2)*dl(3) ! we integrate over the liquid vol. 
          qka_tot_lat = qka_tot_lat + mflux(i,j,k)*lheat*norm*dl(1)*dl(2)*dl(3) ! we integrate over the liquid surf.
          !
          ! quantities at the interface
          !
          if( phi(ip,j,k)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(ip,j,k))+tmple(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(ip,j,k))+tmpge(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t   = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          if( phi(im,j,k)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(im,j,k))+tmple(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(im,j,k))+tmpge(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t    = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          !
          if( phi(i,jp,k)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(i,jp,k))+tmple(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(i,jp,k))+tmpge(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t   = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          if( phi(i,jm,k)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(i,jm,k))+tmple(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(i,jm,k))+tmpge(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t   = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          !
          if( phi(i,j,kp)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(i,j,kp))+tmple(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(i,j,kp))+tmpge(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t   = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          if( phi(i,j,km)*phi(i,j,k).lt.0. ) then
            tmpl_i = (tmple(i,j,k)*abs(phi(i,j,km))+tmple(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
            tmpg_i = (tmpge(i,j,k)*abs(phi(i,j,km))+tmpge(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
            s_in   = mass_fraction(pth,tmpl_i)
            s_in_t = s_in_t + s_in
            counter_t   = counter_t + 1
            kappag_in_t = kappag_in_t + thermo_kag(tmpg_i)
            tmpl_in_t   = tmpl_in_t + tmpl_i
            tmpg_in_t   = tmpg_in_t + tmpg_i
            rhog_ap_in  = thermo_rhog(pth,tmpg_i,s_in)
            rhog_in_t   = rhog_in_t + rhog_ap_in
            d_lg_ap_in  = thermo_d_lg(pth,tmpg_i,s_in)
            d_lg_in_t   = d_lg_in_t + d_lg_ap_in
          endif
          !
          rcounter_t = 1.d0*counter_t
          !
        enddo
      enddo
    enddo
    !
    ! a1. EVAPORATION (overall contribution)
    ! 
    call mpi_allreduce(MPI_IN_PLACE,vol_d      ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,area_d     ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mfx_tot    ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tmp_d      ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tmp_g      ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,y_g        ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,qka_tot_sen,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,qka_tot_lat,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    call mpi_allreduce(MPI_IN_PLACE,rhog_in_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,d_lg_in_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,kappag_in_t,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,s_in_t     ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,rcounter_t ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tmpl_in_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,tmpg_in_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    tmp_d = tmp_d/vol_d
    tmp_g = tmp_g/mas_g 
    y_g   = y_g  /mas_g 
    !
    if(myid.eq.0) then
      !
      ! a2. evaporation (overall)
      !
#ifdef VAP_MASS
      rhog_in_t   = rhog_in_t  /rcounter_t
      d_lg_in_t   = d_lg_in_t  /rcounter_t
      kappag_in_t = kappag_in_t/rcounter_t
      s_in_t      = s_in_t     /rcounter_t
      tmpl_in_t   = tmpl_in_t  /rcounter_t
      tmpg_in_t   = tmpg_in_t  /rcounter_t
#else
      rhog_in_t   = rho2_0
      d_lg_in_t   = d_m12
      kappag_in_t = kappa2
      s_in_t      = mass_fraction(pth,tmpl0)
      tmpl_in_t   = tmpl0
      tmpg_in_t   = tmpg0
#endif
      !
      open(94,file=trim(datapos)//'evaporation_tot.out',position='append')
      write(94,'(16E15.7)') 1.d0*istep,time,vol_d,area_d,tmp_d,tmp_g,y_g,qka_tot_sen,qka_tot_lat, &
                            mfx_tot,rhog_in_t,d_lg_in_t,kappag_in_t,s_in_t,tmpl_in_t,tmpg_in_t
      close(94)
      !
      ! a3. low Mach
      !
      open(94,file=trim(datapos)//'low_mach.out',position='append')
      write(94,'(5E15.7)') 1.d0*istep,time,pth,dpthdt,mas_g
      close(94)
      ! 
    endif
    !
    ! a4. sherwood/nusselt
    !
    do int_s=1,num+1
      !
      ! Sherwood
      !
#ifdef VAP_MASS
      h_mas        = mfx_tot /(area_d*rhog_in_t*(s_in_t-yinf_v(int_s)))
      sh_vp(int_s) = h_mas*lref/d_lg_in_t   ! sherwood var. prop.
      h_mas        = mfx_tot /(area_d*rho2_0*(s_in_t-yinf_v(int_s)))
      sh_cp(int_s) = h_mas*lref/d_m12       ! sherwood constant prop.
#else
      sh_vp(int_s) = 0.d0
      sh_cp(int_s) = 0.d0
#endif
      !
      ! Nusselt
      !
      h_tmp_sen     = qka_tot_sen/(area_d*(tmpg_in_t-tinf_v(int_s))) ! sensible
      h_tmp_lat     = qka_tot_lat/(area_d*(tmpg_in_t-tinf_v(int_s))) ! latent
      nu_sen(int_s) = h_tmp_sen*(6.d0*vol_d/area_d)/kappag_in_t   ! nusselt sensible
      nu_lat(int_s) = h_tmp_lat*(6.d0*vol_d/area_d)/kappag_in_t   ! nusselt latent
      !
    enddo
    !
    call out0d(trim(datapos)//'sh_ga.out',2+4*(num+1),(/1.d0*istep,time, &
                                                        sh_vp(1:num+1) ,sh_cp(1:num+1), &
                                                        yinf_v(1:num+1),tinf_v(1:num+1)/))
    call out0d(trim(datapos)//'nu_ga.out',2+2*(num+1),(/1.d0*istep,time, &
                                                        nu_sen(1:num+1),nu_lat(1:num+1)/))
    !
    return
  end subroutine int_qtn
  !
  subroutine compute_vorticity(n,dli,nh_u,v,w,vor)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                       :: n
    real(8), intent(in ), dimension(3)                       :: dli
    integer, intent(in )                                     :: nh_u
    real(8), intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: v,w
    real(8), intent(out), dimension(     1:,     1:,     1:) :: vor
    !
    real(8) :: vp,vm,wp,wm
    integer :: i,j,k,im,jm,km,ip,jp,kp
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
          vp = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,kp)+v(i,jm,kp))
          vm = 0.25*(v(i,j,k)+v(i,jm,k)+v(i,j,km)+v(i,jm,km))
          wp = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jp,k)+w(i,jp,km))
          wm = 0.25*(w(i,j,k)+w(i,j,km)+w(i,jm,k)+w(i,jm,km))
          !
          vor(i,j,k) = (wp-wm)*dli(2)-(vp-vm)*dli(3)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine compute_vorticity
  !
  subroutine mixed_variables(n,dli,nh_u,nh_s1,nh_s2,u,v,w,s1,s2, &
                             us1,vs1,ws1,uv,vw,wu,us2,vs2,ws2)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                          :: n
    real(8), intent(in ), dimension(3)                          :: dli
    integer, intent(in )                                        :: nh_u,nh_s1,nh_s2
    real(8), intent(in ), dimension(1-nh_u :,1-nh_u :,1-nh_u :) :: u,v,w
    real(8), intent(in ), dimension(1-nh_s1:,1-nh_s1:,1-nh_s1:) :: s1 ! generic scalar
    real(8), intent(in ), dimension(1-nh_s2:,1-nh_s2:,1-nh_s2:) :: s2 ! generic scalar
    real(8), intent(out), dimension(     1:,     1:,     1:)    :: us1,vs1,ws1
    real(8), intent(out), dimension(     1:,     1:,     1:)    :: uv ,vw ,wu
    real(8), intent(out), dimension(     1:,     1:,     1:)    :: us2,vs2,ws2
    !
    integer :: i,j,k,im,jm,km
    !
    do k=1,n(3)
      km = k-1
      do j=1,n(2)
        jm = j-1
        do i=1,n(1)
          im = i-1
          !
          uv(i,j,k)  = 0.25d0*(u(i,j,k)+u(im,j,k))*(v(i,j,k)+v(i,jm,k))
          vw(i,j,k)  = 0.25d0*(v(i,j,k)+v(i,jm,k))*(w(i,j,k)+w(i,j,km))
          wu(i,j,k)  = 0.25d0*(w(i,j,k)+w(i,j,km))*(u(i,j,k)+u(im,j,k))
          !
          us1(i,j,k) = 0.5d0*(u(i,j,k)+u(im,j,k))*s1(i,j,k)
          vs1(i,j,k) = 0.5d0*(v(i,j,k)+v(i,jm,k))*s1(i,j,k)
          ws1(i,j,k) = 0.5d0*(w(i,j,k)+w(i,j,km))*s1(i,j,k)
          !
          us2(i,j,k) = 0.5d0*(u(i,j,k)+u(im,j,k))*s2(i,j,k)
          vs2(i,j,k) = 0.5d0*(v(i,j,k)+v(i,jm,k))*s2(i,j,k)
          ws2(i,j,k) = 0.5d0*(w(i,j,k)+w(i,j,km))*s2(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine mixed_variables
  !
  subroutine time_avg_mt(do_avg,do_favre,fname,n,ng,istep,istep_av,iout1d,idir, &
                         nh_d,nh_v,nh_p,z,dzlzi,psi,rho_p,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! do_avg   -> do or not averaging
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! nh_      -> halos point
    ! z        -> z coordinate (grid is non-uniform in z)
    ! dzlzi    -> dz/lz weight of a grid cell for averaging over z
    ! psi      -> 3D vof field (0 --> liquid, 1 --> gas)
    ! rho_p    -> density of the phase (for compressible phase)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pvol2    -> second order time statistics of volume averaged field (rms)
    !
    use mod_param, only: dl
    !
    implicit none
    !
    logical         , intent(in   )                                     :: do_avg
    logical         , intent(in   )                                     :: do_favre
    character(len=*), intent(in   )                                     :: fname
    integer         , intent(in   ), dimension(3)                       :: n,ng
    integer         , intent(in   )                                     :: istep,istep_av,iout1d
    integer         , intent(in   )                                     :: idir
    integer         , intent(in   )                                     :: nh_d,nh_v,nh_p
    real(8)         , intent(in   ), dimension(1-nh_d:)                 :: z,dzlzi
    real(8)         , intent(in   ), dimension(1-nh_v:,1-nh_v:,1-nh_v:) :: psi
    real(8)         , intent(in   ), dimension(0     :,0     :,0     :) :: rho_p
    real(8)         , intent(in   ), dimension(1-nh_p:,1-nh_p:,1-nh_p:) :: p
    real(8)         , intent(inout), dimension(1:)                      :: pout1,pout2
    real(8)         , intent(inout)                                     :: pvol1,pvol2
    !
    real(8), dimension(1:n(3)) :: p1d1,p1d2,rhop
    real(8) :: factor,p_dl12
    integer :: i,j,k,ii,jj
    integer :: iunit
    !
    p_dl12 = dl(1)*dl(2)
    !
    iunit  = 10
    factor = 1.d0*istep_av
    !
    if(do_favre) then
      !
      ! Density-based averaging (Favre)
      !
      do k=1,n(3)
        p1d1(k)  = 0.d0
        p1d2(k)  = 0.d0
        rhop(k)  = 0.d0
        do j=1,n(2)
          do i=1,n(1)
            !
            rhop(k) = rhop(k) + p_dl12*rho_p(i,j,k)
            p1d1(k) = p1d1(k) + p_dl12*rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)     
            p1d2(k) = p1d2(k) + p_dl12*rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)**2 
            !p1d1(k) = p1d1(k) + (1.d0-psi(i,j,k))*p(i,j,k)
            !p1d2(k) = p1d2(k) + (1.d0-psi(i,j,k))*p(i,j,k)**2
            !
          enddo
        enddo
      enddo
      !
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      !p1d1(:) = p1d1(:)/(1.d0*ng(1)*ng(2))
      !p1d2(:) = p1d2(:)/(1.d0*ng(1)*ng(2))
      rhop(:) = rhop(:)
      p1d1(:) = p1d1(:)/(rhop(:))
      p1d2(:) = p1d2(:)/(rhop(:))
      !
    else
      !
      ! Volumetric-based averaging (no Favre)
      !
      do k=1,n(3)
        p1d1(k)  = 0.d0
        p1d2(k)  = 0.d0
        !rhop(k)  = 0.d0
        do j=1,n(2)
          do i=1,n(1)
            !
            !rhop(k) = rhop(k) + p_dl12*rho_p(i,j,k)
            !p1d1(k) = p1d1(k) + p_dl12*rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)     
            !p1d2(k) = p1d2(k) + p_dl12*rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)**2 
            p1d1(k) = p1d1(k) + (1.d0-psi(i,j,k))*p(i,j,k)
            p1d2(k) = p1d2(k) + (1.d0-psi(i,j,k))*p(i,j,k)**2
            !
          enddo
        enddo
      enddo
      !
      !call mpi_allreduce(MPI_IN_PLACE,rhop(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),n(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      p1d1(:) = p1d1(:)/(1.d0*ng(1)*ng(2))
      p1d2(:) = p1d2(:)/(1.d0*ng(1)*ng(2))
      !rhop(:) = rhop(:)
      !p1d1(:) = p1d1(:)/(rhop(:))
      !p1d2(:) = p1d2(:)/(rhop(:))
    endif
    !
    ! decide or not to averaging 
    !
    if(.not.do_avg) then
      if(myid.eq.0) print*, "NO, I am NOT doing time-averaging"
      pout1(:) = p1d1(:)
      pout2(:) = p1d2(:)
    else
      if(myid.eq.0) print*, "YES, I am doing time-averaging"
      pout1(:) = ((factor-1.d0)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.d0)*pout2(:)+p1d2(:))/factor
    endif
    pvol1 = sum(pout1)
    pvol2 = sum(pout2)
    !
    ! print 
    !  note: we put this condition on iout1d in order to ensure that the 
    !        averaging frequency (iout_av) is indenpendent 
    !        of the print frequency of the files (iout1d)
    !
    if(mod(istep,iout1d).eq.0) then
      if(myid.eq.0) print*, "I am also printing the sample at", istep
      if(myid.eq.0) then
        open(unit=iunit,file=fname)
        do k=1,n(3)
          write(iunit,'(3E15.7)') z(k),pout1(k),pout2(k)
        enddo
        close(iunit)
      endif
    endif
    !
    return
  end subroutine time_avg_mt
  !
  subroutine time_avg(fname,n,ng,istep,istep_av,iout1d,idir,z,dzlzi,psi,rho_p,p,pout1,pout2,pvol1,pvol2)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname    -> name of the file
    ! n        -> size of the input array
    ! ng       -> total size of computational domain
    ! istep    -> current time step
    ! istep_av -> size of statistical sample
    ! iout1d   -> print file every iout1d time steps
    ! idir     -> direction of the profile (all other directions are averaged)
    ! z        -> z coordinate (grid is non-uniform in z)
    ! dzlzi    -> dz/lz weight of a grid cell for averaging over z
    ! psi      -> 3D vof field (0 --> liquid, 1 --> gas)
    ! rho_p    -> density of the phase (for compressible phase)
    ! p        -> 3D input scalar field to be averaged
    ! pout1    ->  first order time statistics of plane averaged field (mean)
    ! pout2    -> second order time statistics of plane averaged field (rms)
    ! pvol     ->  first order time statistics of volume averaged field (mean)
    ! pout2    -> second order time statistics of volume averaged field (rms)
    !
    implicit none
    !
    character(len=*), intent(in   )                      :: fname
    integer         , intent(in   ), dimension(3)        :: n,ng
    integer         , intent(in   )                      :: istep,istep_av,iout1d
    integer         , intent(in   )                      :: idir
    real(8)         , intent(in   ), dimension(0:)       :: z,dzlzi
    real(8)         , intent(in   ), dimension(1:,1:,1:) :: psi
    real(8)         , intent(in   ), dimension(1:,1:,1:) :: rho_p
    real(8)         , intent(in   ), dimension(1:,1:,1:) :: p
    real(8)         , intent(inout), dimension(1:)       :: pout1,pout2
    real(8)         , intent(inout)                      :: pvol1,pvol2
    !
    real(8), allocatable, dimension(:) :: p1d1,p1d2,rhop
    real(8) :: factor
    integer :: i,j,k,ii,jj
    integer :: iunit
    !
    iunit  = 10
    factor = 1.d0*istep_av
    !
    select case(idir)
    !
    ! along z
    !
    case(3)
      allocate(p1d1(n(3)))
      allocate(p1d2(n(3)))
      allocate(rhop(n(3)))
      do k=1,ng(3)
        p1d1(k) = 0.
        p1d2(k) = 0.
        do j=1,n(2)
          do i=1,n(1)
            rhop(k) = rhop(k) + rho_p(i,j,k)
            p1d1(k) = p1d1(k) + rho_p(i,j,i)*(1.d0-psi(i,j,k))*p(i,j,k)
            p1d2(k) = p1d2(k) + rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)**2
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(3),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      rhop(:)  = rhop(:)/(1.*ng(1)*ng(2)        )
      p1d1(:)  = p1d1(:)/(1.*ng(1)*ng(2)*rhop(:))
      p1d2(:)  = p1d2(:)/(1.*ng(1)*ng(2)*rhop(:))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do k=1,n(3)
            write(iunit,'(3E15.7)') z(k),pout1(k),pout2(k)
          enddo
          close(iunit)
        endif
      endif
    !
    ! along y
    !
    case(2)
      allocate(p1d1(ng(2)))
      allocate(p1d2(ng(2)))
      allocate(rhop(ng(2)))
      p1d1(:) = 0.
      p1d2(:) = 0.
      do j=1,n(2)
        jj = coord(2)*n(2)+j
        p1d1(jj) = 0.
        p1d2(jj) = 0.
        do k=1,n(3)
          do i=1,n(1)
            rhop(jj) = rhop(jj) + rho_p(i,j,k) 
            p1d1(jj) = p1d1(jj) + rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)*dzlzi(k)
            p1d2(jj) = p1d2(jj) + rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)*p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng(2),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(2),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(2),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      rhop(:)  = rhop(:)/(1.*ng(1)*ng(3)        )
      p1d1(:)  = p1d1(:)/(1.*ng(1)*ng(3)*rhop(:))
      p1d2(:)  = p1d2(:)/(1.*ng(1)*ng(3)*rhop(:))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do j=1,ng(2)
            write(iunit,'(3E15.7)') (1.*j-.5)/(1.*ng(2)),pout1(j),pout2(j)
          enddo
          close(iunit)
        endif
      endif
    !
    ! along z
    !
    case(1)
      allocate(p1d1(ng(1)))
      allocate(p1d2(ng(1)))
      allocate(rhop(ng(1)))
      p1d1(:) = 0.
      p1d2(:) = 0.
      do i=1,n(1)
        ii = coord(1)*n(1)+i
        p1d1(i) = 0.
        p1d2(i) = 0.
        do k=1,n(3)
          do j=1,n(2)
            rhop(ii) = rhop(ii) + rho_p(i,j,k) 
            p1d1(ii) = p1d1(ii) + rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)*dzlzi(k)
            p1d2(ii) = p1d2(ii) + rho_p(i,j,k)*(1.d0-psi(i,j,k))*p(i,j,k)*p(i,j,k)*dzlzi(k)
          enddo
        enddo
      enddo
      call mpi_allreduce(MPI_IN_PLACE,rhop(1),ng(1),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d1(1),ng(1),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call mpi_allreduce(MPI_IN_PLACE,p1d2(1),ng(1),mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      rhop(:)  = rhop(:)/(1.*ng(2)*ng(3)        )
      p1d1(:)  = p1d1(:)/(1.*ng(2)*ng(3)*rhop(:))
      p1d2(:)  = p1d2(:)/(1.*ng(2)*ng(3)*rhop(:))
      pout1(:) = ((factor-1.)*pout1(:)+p1d1(:))/factor
      pout2(:) = ((factor-1.)*pout2(:)+p1d2(:))/factor
      pvol1    = sum(pout1)
      pvol2    = sum(pout2)
      if(mod(istep,iout1d).eq.0) then
        if(myid.eq.0) then
          open(unit=iunit,file=fname)
          do i=1,ng(1)
            write(iunit,'(3E15.7)') (1.*i-.5)/(1.*n(1)),pout1(i),pout2(i)
          enddo
          close(iunit)
        endif
      endif
    end select
    deallocate(p1d1,p1d2,rhop)
    !
    return
  end subroutine time_avg
  !
  subroutine vorticity(n,dli,nh_u,ux,uy,uz,omega_x,omega_y,omega_z)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                       :: n
    real(8), intent(in ), dimension(3)                       :: dli
    integer, intent(in )                                     :: nh_u
    integer, intent(in ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: ux,uy,uz
    real(8), intent(out), dimension(     0:,     0:,     0:) :: omega_x,omega_y,omega_z
    !
    integer :: i,j,k
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
  subroutine density_cont(n,h,dl,pth,dpthdt_n,vof,rho,rho_gas,d_lg,kappa,tmp,sca,istep,time) 
    !
    use mod_param, only: delta_cp,datapos,ru,cp2,m1,m2
    !
    implicit none
    !
    integer, intent(in ), dimension(3)        :: n
    integer, intent(in )                      :: h
    real(8), intent(in ), dimension(3)        :: dl
    real(8), intent(in )                      :: pth,dpthdt_n
    real(8), intent(in ), dimension(0:,0:,0:) :: vof,rho,rho_gas,d_lg,kappa
    real(8), intent(in ), dimension(h:,h:,h:) :: tmp
    real(8), intent(in ), dimension(0:,0:,0:) :: sca
    integer, intent(in )                      :: istep
    real(8), intent(in )                      :: time
    !
    real(8), dimension(3) :: dli
    real(8) :: t_pth,t_vap_v,t_tmp_v,t_vap_m,t_tmp_m,vol_g,mas_g, &
               diff_s,diff_t,m_avg
    real(8) :: dscadxm,dscadxp,dscadym,dscadyp,dscadzm,dscadzp
    real(8) :: dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm
    real(8) :: rhdlgxp,rhdlgxm,rhdlgyp,rhdlgym,rhdlgzp,rhdlgzm
    real(8) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    integer :: i,j,k,im,ip,jm,jp,km,kp
    !
    dli(:) = dl(:)**(-1.d0)
    !
    t_vap_v = 0.d0
    t_tmp_v = 0.d0
    t_vap_m = 0.d0
    t_tmp_m = 0.d0
    vol_g   = 0.d0
    mas_g   = 0.d0
    !
    t_pth = -dpthdt_n
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
          ! 0. vol. of the gas
          !
          vol_g = vol_g + (1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)
          mas_g = mas_g + (1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k)
          !
          ! 1. contribution due to p_th
          !
          !term1 = term1 + (-dpthdt_n) ! not needed we have
          !
          ! 2. contribution due to change in composition
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
          t_vap_v = t_vap_v + (1.d0-vof(i,j,k))*(diff_s/pth)*dl(1)*dl(2)*dl(3)
          t_vap_m = t_vap_m + (1.d0-vof(i,j,k))*(diff_s/pth)*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k) ! mass-based
          !
          ! 3. contribution due to change in temperature
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
                                         )*0.25d0*rho(i,j,k)*(ru/(cp2*m_avg)) +  &
                    dpthdt_n*pth/(rho(i,j,k)*cp2)
          !
          t_tmp_v = t_tmp_v + (1.d0-vof(i,j,k))*(diff_t/pth)*dl(1)*dl(2)*dl(3)
          t_tmp_m = t_tmp_m + (1.d0-vof(i,j,k))*(diff_t/pth)*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k) ! mass-based
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,t_vap_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! vol.-based
    call mpi_allreduce(MPI_IN_PLACE,t_tmp_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! vol.-based
    call mpi_allreduce(MPI_IN_PLACE,t_vap_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! mass-based
    call mpi_allreduce(MPI_IN_PLACE,t_tmp_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! mass-based
    call mpi_allreduce(MPI_IN_PLACE,mas_g  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! mass of the gas
    call mpi_allreduce(MPI_IN_PLACE,vol_g  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! vol. of the gas
    !
    t_vap_v = t_vap_v/vol_g
    t_tmp_v = t_tmp_v/vol_g
    t_vap_m = t_vap_m/mas_g
    t_tmp_m = t_tmp_m/mas_g
    !
    if(myid.eq.0) then
      !
      ! a. contribution to the change in density
      !
      open(92,file=trim(datapos)//'rhog_con.out',position='append')
      !write(92,'(A)') '# istep, tS, vol-vap cont, vol-tmp cont, mas-vap cont, mas-tmp cont, pth cont, vol_g, mass'
      !
      write(92,'(9E15.7)') 1.d0*istep,time, &
                           t_vap_v,t_tmp_v,t_vap_m,t_tmp_m,t_pth, &
                           vol_g,mas_g
      close(92)
      !
    endif
    !
    return
  end subroutine density_cont
  !
  subroutine energy_th(n,h,dl,pth,dpthdt,vof,delta,rho,rho_gas,d_lg,kappa,tmp,sca,mflux,istep,time) 
    !
    ! note:
    !  --> we write in the form of enthalpy balance and then we apply the definition
    !      i = h - (p_th/rho);
    !  --> we do not need to compute the overall enthalpy balance. Since the system is
    !      closed, the total energy should be conserved and, therefore, we should compute
    !      e_t=i+0.5*rho*v^2. So, knowing the MKE and e_T(0), we know immediately i;
    !  --> so, we perform for the two phases separately.
    !
    use mod_param, only: datapos,ru, &
                         cp1,cp2,delta_cp,rho1,lheat,tmpg0,tmpl0,sinit
    !
    implicit none
    !
    integer, intent(in), dimension(3)        :: n
    integer, intent(in)                      :: h
    real(8), intent(in), dimension(3)        :: dl
    real(8), intent(in)                      :: pth,dpthdt
    real(8), intent(in), dimension(0:,0:,0:) :: vof,delta,rho,rho_gas,d_lg,kappa
    real(8), intent(in), dimension(h:,h:,h:) :: tmp
    real(8), intent(in), dimension(0:,0:,0:) :: sca
    real(8), intent(in), dimension(0:,0:,0:) :: mflux
    integer, intent(in)                      :: istep
    real(8), intent(in)                      :: time
    !
    real(8), dimension(3) :: dli
    real(8) :: norm
    real(8) :: q_ga,q_s1_p1,q_s1_p2,q_s2_p2,h_p1,h_p2 
    real(8) :: dscadxm,dscadxp,dscadym,dscadyp,dscadzm,dscadzp
    real(8) :: dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm
    real(8) :: rhdlgxp,rhdlgxm,rhdlgyp,rhdlgym,rhdlgzp,rhdlgzm
    real(8) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    real(8) :: lap_kap,diff_ss,vol_d,mas_g
    integer :: i,j,k,im,ip,jm,jp,km,kp
    !
    dli(:) = dl(:)**(-1.d0)
    !
    q_s1_p1 = 0.d0
    q_s1_p2 = 0.d0
    q_s2_p2 = 0.d0
    q_ga    = 0.d0
    !
    h_p1 = 0.d0
    h_p2 = 0.d0
    !
    vol_d = 0.d0
    mas_g = 0.d0
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
          ! vol of the liquid/mass of the gas
          !
          vol_d = vol_d + (     vof(i,j,k))*dl(1)*dl(2)*dl(3)
          mas_g = mas_g + (1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k)
          !
          norm  = delta(i,j,k)
          !
          kappaxp = 0.5d0*(kappa(ip,j,k)+kappa(i,j,k))
          kappaxm = 0.5d0*(kappa(im,j,k)+kappa(i,j,k))
          kappayp = 0.5d0*(kappa(i,jp,k)+kappa(i,j,k))
          kappaym = 0.5d0*(kappa(i,jm,k)+kappa(i,j,k))
          kappazp = 0.5d0*(kappa(i,j,kp)+kappa(i,j,k))
          kappazm = 0.5d0*(kappa(i,j,km)+kappa(i,j,k))
          !
          dtmpdxp = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
          dtmpdxm = (tmp(i ,j,k)-tmp(im,j,k))*dli(2)
          dtmpdyp = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
          dtmpdym = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
          dtmpdzp = (tmp(i,j,kp)-tmp(i,j,k ))*dli(3)
          dtmpdzm = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
          !
          lap_kap = (kappaxp*dtmpdxp-kappaxm*dtmpdxm)*dli(1) + &
                    (kappayp*dtmpdyp-kappaym*dtmpdym)*dli(2) + &
                    (kappazp*dtmpdzp-kappazm*dtmpdzm)*dli(3)
          !
          dscadxm = (sca(i ,j,k)-sca(im,j,k))*dli(1)
          dscadxp = (sca(ip,j,k)-sca(i ,j,k))*dli(1)
          dscadym = (sca(i,j ,k)-sca(i,jm,k))*dli(2)
          dscadyp = (sca(i,jp,k)-sca(i,j ,k))*dli(2)
          dscadzm = (sca(i,j,k )-sca(i,j,km))*dli(3)
          dscadzp = (sca(i,j,kp)-sca(i,j ,k))*dli(3)
          ! 
          rhdlgxm = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(im,j,k)*d_lg(im,j,k)*(tmp(im,j,k)-tmpg0))
          rhdlgxp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(ip,j,k)*d_lg(ip,j,k)*(tmp(ip,j,k)-tmpg0))
          rhdlgym = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(i,jm,k)*d_lg(i,jm,k)*(tmp(i,jm,k)-tmpg0))
          rhdlgyp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(i,jp,k)*d_lg(i,jp,k)*(tmp(i,jp,k)-tmpg0))
          rhdlgzm = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(i,j,km)*d_lg(i,j,km)*(tmp(i,j,km)-tmpg0))
          rhdlgzp = 0.5d0*(rho_gas(i,j,k)*d_lg(i,j,k)*(tmp(i,j,k)-tmpg0)+rho_gas(i,j,kp)*d_lg(i,j,kp)*(tmp(i,j,kp)-tmpg0))
          !
          diff_ss = ((dscadxp*rhdlgxp-dscadxm*rhdlgxm)*dli(1) + &
                     (dscadyp*rhdlgyp-dscadym*rhdlgym)*dli(2) + &
                     (dscadzp*rhdlgzp-dscadzm*rhdlgzm)*dli(3))*delta_cp
          !
          ! a. liquid contribution (integration over the liquid volume)
          !
          q_s1_p1 = q_s1_p1 + lap_kap*vof(i,j,k)*dl(1)*dl(2)*dl(3) 
          !
          ! b. gas contribution (integration over the mass of the gas)
          !
          q_s1_p2 = q_s1_p2 + lap_kap*(1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k)  
          q_s2_p2 = q_s2_p2 + diff_ss*(1.d0-vof(i,j,k))*dl(1)*dl(2)*dl(3)*rho_gas(i,j,k)  
          !
          ! c. interfacial contribution
          !
          q_ga = q_ga + mflux(i,j,k)*lheat*norm*dl(1)*dl(2)*dl(3) ! we integrate over the interface
          !
          ! d. liquid enthalpy
          !
          h_p1 = h_p1 + rho1*cp1*(tmp(i,j,k)-tmpl0)*vof(i,j,k)*dl(1)*dl(2)*dl(3)
          !
          ! e. gas enthalpy
          !
          h_p2 = h_p2 + rho_gas(i,j,k)*(1.d0-vof(i,j,k))*( &
                        rho_gas(i,j,k)*cp2*(tmp(i,j,k)-tmpg0) + &
                        rho_gas(i,j,k)*(delta_cp*(tmp(i,j,k)-tmpg0))*(sca(i,j,k)-sinit) &
                                                         )
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,vol_d,1,mpi_real8,mpi_sum,comm_cart,ierr) ! liquid
    call mpi_allreduce(MPI_IN_PLACE,mas_g,1,mpi_real8,mpi_sum,comm_cart,ierr) ! liquid
    !
    ! right hand side
    ! 
    call mpi_allreduce(MPI_IN_PLACE,q_s1_p1,1,mpi_real8,mpi_sum,comm_cart,ierr) ! liquid
    !
    call mpi_allreduce(MPI_IN_PLACE,q_s1_p2,1,mpi_real8,mpi_sum,comm_cart,ierr) ! gas
    call mpi_allreduce(MPI_IN_PLACE,q_s2_p2,1,mpi_real8,mpi_sum,comm_cart,ierr) ! gas
    !call mpi_allreduce(MPI_IN_PLACE,q_s2_p3,1,mpi_real8,mpi_sum,comm_cart,ierr) ! gas (dpthdt)
    !
    call mpi_allreduce(MPI_IN_PLACE,q_ga   ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! interface
    !
    ! left hand side
    !
    call mpi_allreduce(MPI_IN_PLACE,h_p1,1,mpi_real8,mpi_sum,comm_cart,ierr) ! liquid
    call mpi_allreduce(MPI_IN_PLACE,h_p2,1,mpi_real8,mpi_sum,comm_cart,ierr) ! gas
    !
    ! volumetric/mass average
    !
    q_s1_p1 = q_s1_p1/vol_d ! liquid
    h_p1    = h_p1/vol_d    ! liquid
    !  
    q_s1_p2 = q_s1_p2/mas_g ! gas
    q_s2_p2 = q_s2_p2/mas_g ! gas 
    h_p2    = h_p2/mas_g    ! gas
    !
    if(myid.eq.0) then
      !
      ! liquid mean internal energy balance
      !
      open(92,file=trim(datapos)//'en_th_p1.out',position='append')
      !
      write(92,'(5E15.7)') 1.d0*istep,time, &
                           h_p1,q_s1_p1,-1.d0*q_ga
      close(92)
      !
      ! gas mean internal energy balance
      !
      open(92,file=trim(datapos)//'en_th_p2.out',position='append')
      !
      write(92,'(7E15.7)') 1.d0*istep,time, &
                           h_p2,q_s1_p2,q_s2_p2,dpthdt,+1.d0*q_ga
      close(92)
      !
    endif
    !
    return
  end subroutine energy_th
  !
  subroutine cmpt_delta(n,dli,vof,delta)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)        :: n
    real(8), intent(in ), dimension(3)        :: dli
    real(8), intent(in ), dimension(0:,0:,0:) :: vof
    real(8), intent(out), dimension(0:,0:,0:) :: delta
    !
    real(8), dimension(8) :: mx,my,mz
    real(8) :: dvofdx,dvofdy,dvofdz
    integer :: i,j,k
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
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
          ! compute a delta Dirac function based on the VoF color indicator
          !
          dvofdx = 0.125d0*sum(mx(1:8))
          dvofdy = 0.125d0*sum(my(1:8))
          dvofdz = 0.125d0*sum(mz(1:8))
          !
          delta(i,j,k) = sqrt(dvofdx**2+dvofdy**2+dvofdz**2)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_delta
  !
  subroutine budget(n,dl,e,u,v,w,p,rho,rho_gas,vof,time,mass_gas,istep)
    !
    use mod_param, only: itot,jtot,ktot,rho1
    !
    implicit none
    !
    integer, intent(in ), dimension(3)        :: n
    real(8), intent(in ), dimension(3)        :: dl
    integer, intent(in )                      :: e
    real(8), intent(in ), dimension(e:,e:,e:) :: u,v,w
    real(8), intent(in ), dimension(0:,0:,0:) :: p,rho,rho_gas
    real(8), intent(in ), dimension(0:,0:,0:) :: vof
    real(8), intent(in )                      :: time,mass_gas
    integer, intent(in )                      :: istep
    !
    integer :: i,j,k,ip,jp,kp,im,jm,km
    real(8) :: ke_all,ke_t,ke_1,ke_2v,ke_2m,dvol,volll,vol_t1,vol_t2
    real(8), parameter :: eps = 1.e-12
    ! 
    volll = 1.d0*itot*jtot*ktot
    dvol  = product(dl(:))
    !
    ! volume
    !
    vol_t1 = 0.d0
    vol_t2 = 0.d0
    !
    ! quantities for the mean kinetic energy balance
    ! 
    ke_t  = 0.d0
    ke_1  = 0.d0
    ke_2v = 0.d0
    ke_2m = 0.d0
    !
    do k=1,n(3)
      kp = k + 1
      km = k - 1
      do j=1,n(2)
        jp = j + 1
        jm = j - 1
        do i=1,n(1)
          ip = i + 1
          im = i - 1
          !
          ! 1. volume
          ! 
          vol_t1 = vol_t1 + vof(i,j,k)
          vol_t2 = vol_t2 + (1.d0-vof(i,j,k))
          !
          ! 2. turbulent kinetic energy
          ! 
          ke_t = ke_t + rho(i,j,k)*&
          0.5d0*(0.25d0*(u(i,j,k)+u(im,j,k))**2 + 0.25d0*(v(i,j,k)+v(i,jm,k))**2 + 0.25d0*(w(i,j,k)+w(i,j,km))**2)
          ke_1 = ke_1 + rho1*vof(i,j,k)*&
          0.5d0*(0.25d0*(u(i,j,k)+u(im,j,k))**2 + 0.25d0*(v(i,j,k)+v(i,jm,k))**2 + 0.25d0*(w(i,j,k)+w(i,j,km))**2)
          ke_2v = ke_2v + (rho_gas(i,j,k)   )*(1.d0-vof(i,j,k))*&
          0.5d0*(0.25d0*(u(i,j,k)+u(im,j,k))**2 + 0.25d0*(v(i,j,k)+v(i,jm,k))**2 + 0.25d0*(w(i,j,k)+w(i,j,km))**2)
          ke_2m = ke_2m + (rho_gas(i,j,k)**2)*(1.d0-vof(i,j,k))*&
          0.5d0*(0.25d0*(u(i,j,k)+u(im,j,k))**2 + 0.25d0*(v(i,j,k)+v(i,jm,k))**2 + 0.25d0*(w(i,j,k)+w(i,j,km))**2)
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,vol_t1,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vol_t2,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_1  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_2v ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,ke_2m ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    ke_t  = ke_t /volll
    ke_1  = ke_1 /(vol_t1+eps)
    ke_2v = ke_2v/(vol_t2+eps)
    ke_2m = ke_2m/(mass_gas+eps)
    !
    if(myid.eq.0) then
      !
      ! a. post-processing one fluid
      !
      open(92,file=trim(datapos)//'ke_t.out',position='append')
      write(92,'(3E15.7)') 1.d0*istep,time,ke_t
      close(92)
      ! 
      ! b1. phase 1
      ! 
      open(93,file=trim(datapos)//'ke_1.out',position='append')
      write(93,'(4E15.7)') 1.d0*istep,time,vol_t1,ke_1
      close(93)
      ! 
      ! b2. phase 2
      !
      open(94,file=trim(datapos)//'ke_2m.out',position='append')
      write(94,'(5E15.7)') 1.d0*istep,time,vol_t2,ke_2v,ke_2m 
      close(94)
      !
    endif
    !
    return
  end subroutine budget
  !
!  subroutine energy_balance(n,dli,nh_d,nh_u,dzc,dzf,time,istep,rho,mu,vof,kappa,p, &
!                            rho_gas,mu_gas,div_th,divg_th,u,v,w)
!    !
!    ! energy balance in physical space
!    !
!    use mod_param     , only: datapos,sigmaca,cbcvof,bcvof,dims,rho1,mu1, &
!                              f0_t,abc,k0_freq,pi!,rho2
!    use mod_common_mpi, only: ierr
!    !
!    implicit none
!    !
!    integer, intent(in  ), dimension(3)                       :: n
!    real(8), intent(in  ), dimension(3)                       :: dli
!    integer, intent(in  )                                     :: nh_d,nh_u
!    real(8), intent(in  ), dimension(1-nh_d:)                 :: dzc,dzf
!    real(8), intent(in  )                                     :: time
!    integer, intent(in  )                                     :: istep
!    real(8), intent(in  ), dimension(     0:,     0:,     0:) :: rho,mu,vof,kappa,p
!    real(8), intent(in  ), dimension(     0:,     0:,     0:) :: rho_gas,mu_gas,div_th,divg_th
!    real(8), intent(in  ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
!    !
!    integer :: i,j,k,ip,im,jp,jm,kp,km,q
!    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)   :: u_p,v_p,w_p
!    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1,6) :: S      
!    real(8), dimension(3) :: dl
!    real(8) :: u_avg,v_avg,w_avg,S_S,E_E, &
!               duTxx,duTxy,duTxz,dvTyx,dvTyy,dvTyz,dwTzx,dwTzy,dwTzz, & 
!               tt_com,tg_com,t_prod
!    real(8) :: ur_t,vr_t,wr_t,ur1_v,vr1_v,wr1_v,ur2_v,vr2_v,wr2_v,ur2_m,vr2_m,wr2_m, &
!               um_t,vm_t,wm_t,um1_v,vm1_v,wm1_v,um2_v,vm2_v,wm2_v,um2_m,vm2_m,wm2_m, &
!               eps_t ,eps1_v,eps2_v,eps2_m, &
!               prd_t ,prd1_v,prd2_v,prd2_m, & 
!               ke_t  ,ke1_v ,ke2_v ,ke2_m , & 
!               ens_t ,ens1_v,ens2_v,ens2_m, &
!               Psi_mf,Psi_nu,Psi_gf_v,Psi_gf_m, &
!               Tp1_v ,Tp2_v ,Tp2_m , &
!               Tnu1_v,Tnu2_v,Tnu2_m, &
!               vol1,vol2,mass_g,vol,obt
!    real(8) :: VGRux , VGRuy , VGRuz , &
!               VGRvx , VGRvy , VGRvz , &
!               VGRwx , VGRwy , VGRwz , &
!               dupx  , dupy  , dupz  , &
!               fx_hit, fy_hit, fz_hit, &
!               xc    , yc    , zc
!    !
!    vol    = 1.d0*(n(1)*dims(1))*(n(2)*dims(2))*n(3)
!    dl(:)  = 1.d0/dli(:)
!    obt    = 1.d0/3.d0
!    !
!    E_E    = 0.d0 
!    ens_t  = 0.d0
!    ens1_v = 0.d0
!    ens2_v = 0.d0
!    ens2_m = 0.d0
!    !
!    ur_t   = 0.d0 
!    vr_t   = 0.d0
!    wr_t   = 0.d0
!    ur1_v  = 0.d0
!    vr1_v  = 0.d0
!    wr1_v  = 0.d0
!    ur2_v  = 0.d0
!    vr2_v  = 0.d0
!    wr2_v  = 0.d0
!    ur2_m  = 0.d0
!    vr2_m  = 0.d0
!    wr2_m  = 0.d0
!    !
!    um_t   = 0.d0 
!    vm_t   = 0.d0
!    wm_t   = 0.d0
!    um1_v  = 0.d0
!    vm1_v  = 0.d0
!    wm1_v  = 0.d0
!    um2_v  = 0.d0
!    vm2_v  = 0.d0
!    wm2_v  = 0.d0
!    um2_m  = 0.d0
!    vm2_m  = 0.d0
!    wm2_m  = 0.d0
!    !
!    eps_t  = 0.d0 
!    eps1_v = 0.d0
!    eps2_v = 0.d0
!    eps2_m = 0.d0
!    !
!    prd_t  = 0.d0
!    prd1_v = 0.d0
!    prd2_v = 0.d0
!    prd2_m = 0.d0
!    !
!    ke_t   = 0.d0
!    ke1_v  = 0.d0
!    ke2_v  = 0.d0
!    ke2_m  = 0.d0
!    !
!    Psi_mf   = 0.d0 
!    Psi_nu   = 0.d0
!    Psi_gf_v = 0.d0
!    Psi_gf_m = 0.d0
!    !
!    Tp1_v  = 0.d0 
!    Tp2_v  = 0.d0
!    Tp2_m  = 0.d0
!    !
!    Tnu1_v = 0.d0 
!    Tnu2_v = 0.d0
!    Tnu2_m = 0.d0
!    ! 
!    ! following Dodd and Ferrante, we compute (u_p,v_p,w_p)
!    !
!    u_avg = 0.d0
!    v_avg = 0.d0
!    w_avg = 0.d0
!    !
!    vol1   = 0.d0
!    vol2   = 0.d0
!    mass_g = 0.d0
!    !
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          !
!          ! average velocity
!          !
!          u_avg = u_avg + u(i,j,k)*dl(1)*dl(2)*dl(3)
!          v_avg = v_avg + v(i,j,k)*dl(1)*dl(2)*dl(3)
!          w_avg = w_avg + w(i,j,k)*dl(1)*dl(2)*dl(3)
!          !
!          ! volume of the two phases and mass of the gas
!          !
!          vol1   = vol1 + vof(i,j,k)
!          vol2   = vol2 + (1.d0-vof(i,j,k))
!          mass_g = mass_g + rho_gas(i,j,k)*(1.d0-vof(i,j,k))
!          !
!        enddo
!      enddo
!    enddo 
!    call mpi_allreduce(MPI_IN_PLACE,u_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    call mpi_allreduce(MPI_IN_PLACE,v_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    call mpi_allreduce(MPI_IN_PLACE,w_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    u_avg = u_avg/(lx*ly*lz)
!    v_avg = v_avg/(lx*ly*lz)
!    w_avg = w_avg/(lx*ly*lz)
!    !
!    call mpi_allreduce(MPI_IN_PLACE,vol1  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    call mpi_allreduce(MPI_IN_PLACE,vol2  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    call mpi_allreduce(MPI_IN_PLACE,mass_g,1,mpi_real8,mpi_sum,comm_cart,ierr)
!    vol1   = vol1   + 1e-16 ! to avoid division by 0
!    vol2   = vol2   + 1e-16 ! to avoid division by 0
!    mass_g = mass_g + 1e-16 ! to avoid division by 0
!    !
!    do k=0,n(3)+1
!      do j=0,n(2)+1
!        do i=0,n(1)+1
!          u_p(i,j,k) = u(i,j,k)-u_avg
!          v_p(i,j,k) = v(i,j,k)-v_avg
!          w_p(i,j,k) = w(i,j,k)-w_avg
!        enddo
!      enddo
!    enddo 
!    !
!    do k=1,n(3)
!      kp = k+1
!      km = k-1
!      do j=1,n(2)
!        jp = j+1
!        jm = j-1
!        do i=1,n(1)
!          ip = i+1
!          im = i-1
!          !
!          ! a. MEAN KINETIC ENERGY BALANCE --> part 1
!          !
!          ! a1. calculation and storage (for the next loop) of the 
!          !     strain rate tensor S_ij
!          !
!          ! along x
!          !
!          VGRux = ( u_p(i,j,k) - u_p(im,j,k) )*dli(1)
!          VGRuy = 0.25*( u_p(i,j,k) + u_p(i,jp,k) + u_p(im,j,k) + u_p(im,jp,k) )*dli(2) - &
!                  0.25*( u_p(i,j,k) + u_p(i,jm,k) + u_p(im,j,k) + u_p(im,jm,k) )*dli(2)  
!          VGRuz = 0.25*( u_p(i,j,k) + u_p(i,j,kp) + u_p(im,j,k) + u_p(im,j,kp) )*dli(3) - &
!                  0.25*( u_p(i,j,k) + u_p(i,j,km) + u_p(im,j,k) + u_p(im,j,km) )*dli(3)
!          !
!          ! along y
!          !  
!          VGRvx = 0.25*( v_p(i,j,k) + v_p(ip,j,k) + v_p(i,jm,k) + v_p(ip,jm,k) )*dli(1) - &
!                  0.25*( v_p(i,j,k) + v_p(im,j,k) + v_p(i,jm,k) + v_p(im,jm,k) )*dli(1)
!          VGRvy = ( v_p(i,j,k) - v_p(i,jm,k) )*dli(2)
!          VGRvz = 0.25*( v_p(i,j,k) + v_p(i,j,kp) + v_p(i,jm,k) + v_p(i,jm,kp) )*dli(3) - &
!                  0.25*( v_p(i,j,k) + v_p(i,j,km) + v_p(i,jm,k) + v_p(i,jm,km) )*dli(3)
!          !
!          ! along z
!          ! 
!          VGRwx = 0.25*( w_p(i,j,k) + w_p(ip,j,k) + w_p(i,j,km) + w_p(ip,j,km) )*dli(1) - &
!                  0.25*( w_p(i,j,k) + w_p(im,j,k) + w_p(i,j,km) + w_p(im,j,km) )*dli(1)
!          VGRwy = 0.25*( w_p(i,j,k) + w_p(i,jp,k) + w_p(i,j,km) + w_p(i,jp,km) )*dli(2) - &
!                  0.25*( w_p(i,j,k) + w_p(i,jm,k) + w_p(i,j,km) + w_p(i,jm,km) )*dli(2)
!          VGRwz = ( w_p(i,j,k) - w_p(i,j,km) )*dli(3)
!          !
!          ! incompressible part
!          ! 
!          !-----S11 = dU/dx
!          S(i,j,k,1) = VGRux
!          !-----S12 = 0.5*(dU/dy+dV/dx) = S21
!          S(i,j,k,2) = 0.5 * (VGRuy + VGRvx)
!          !-----S13 = 0.5*(dU/dz+dW/dx) = S31
!          S(i,j,k,3) = 0.5 * (VGRuz + VGRwx)
!          !-----S22 = dV/dy
!          S(i,j,k,4) = VGRvy
!          !-----S23 = 0.5*(dV/dz+dW/dy) = S32
!          S(i,j,k,5) = 0.5 * (VGRvz + VGRwy)
!          !-----S33 = dW/dz
!          S(i,j,k,6) = VGRwz
!          !
!          ! add the compressible contribution
!          ! 
!          S(i,j,k,1) = S(i,j,k,1) - obt*divg_th(i,j,k)
!          S(i,j,k,4) = S(i,j,k,4) - obt*divg_th(i,j,k)
!          S(i,j,k,6) = S(i,j,k,6) - obt*divg_th(i,j,k)
!          !
!          ! a2. enstrophy
!          !
!          E_E    = (VGRwy-VGRvz)**2. + (VGRuz-VGRwx)**2. + (VGRvx-VGRuy)**2. 
!          ens_t  = ens_t  + E_E
!          ens1_v = ens1_v + E_E*vof(i,j,k)
!          ens2_v = ens2_v + E_E*(1.d0-vof(i,j,k))
!          ens2_m = ens2_m + E_E*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          !
!        enddo
!      enddo
!    enddo
!    ! 
!    do q=1,6
!      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,S(:,:,:,q))
!    enddo
!    !
!    do k=1,n(3)
!      kp = k+1
!      km = k-1
!      zc = (k-0.5d0)*dl(3)/lz*2.d0*pi
!      do j=1,n(2)
!        jp = j+1
!        jm = j-1
!        yc = (j+coord(2)*n(2)-0.5d0)*dl(2)/ly*2.d0*pi
!        do i=1,n(1)
!          ip = i+1
!          im = i-1
!          xc = (i+coord(1)*n(1)-0.5d0)*dl(1)/lx*2.d0*pi
!          !
!          ! b1. RMS for one-fluid and two phases
!          !
!          ur_t = ur_t + (u_p(i,j,k))**2
!          vr_t = vr_t + (v_p(i,j,k))**2
!          wr_t = wr_t + (w_p(i,j,k))**2
!          ! 
!          ur1_v = ur1_v + ((u_p(i,j,k))**2)*vof(i,j,k)
!          vr1_v = vr1_v + ((v_p(i,j,k))**2)*vof(i,j,k)
!          wr1_v = wr1_v + ((w_p(i,j,k))**2)*vof(i,j,k)
!          !
!          ur2_v = ur2_v + ((u_p(i,j,k))**2)*(1.d0-vof(i,j,k))
!          vr2_v = vr2_v + ((v_p(i,j,k))**2)*(1.d0-vof(i,j,k))
!          wr2_v = wr2_v + ((w_p(i,j,k))**2)*(1.d0-vof(i,j,k))
!          ur2_m = ur2_m + ((u_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          vr2_m = vr2_m + ((v_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          wr2_m = wr2_m + ((w_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          !
!          ! b2. mean velocities for one-fluid and two phases
!          !
!          um_t = um_t + u_p(i,j,k)
!          vm_t = vm_t + v_p(i,j,k)
!          wm_t = wm_t + w_p(i,j,k)
!          ! 
!          um1_v = um1_v + u_p(i,j,k)*vof(i,j,k)
!          vm1_v = vm1_v + v_p(i,j,k)*vof(i,j,k)
!          wm1_v = wm1_v + w_p(i,j,k)*vof(i,j,k)
!          !
!          um2_v = um2_v + u_p(i,j,k)*(1.d0-vof(i,j,k))
!          vm2_v = vm2_v + v_p(i,j,k)*(1.d0-vof(i,j,k))
!          wm2_v = wm2_v + w_p(i,j,k)*(1.d0-vof(i,j,k))
!          um2_m = um2_m + u_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          vm2_m = vm2_m + v_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          wm2_m = wm2_m + w_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          !
!          ! c. MEAN KINETIC ENERGY BALANCE --> part 2
!          !
!          ! c1. dissipation
!          !
!          !S_S =     (S(i,j,k,1)+obt*divg_th(i,j,k))*(S(i,j,k,1)+obt*divg_th(i,j,k)) + &
!          !      2.d0*S(i,j,k,2)*S(i,j,k,2) + &
!          !      2.d0*S(i,j,k,3)*S(i,j,k,3) + & 
!          !          (S(i,j,k,4)+obt*divg_th(i,j,k))*(S(i,j,k,4)+obt*divg_th(i,j,k)) + &
!          !      2.d0*S(i,j,k,5)*S(i,j,k,5) + &
!          !          (S(i,j,k,6)+obt*divg_th(i,j,k))*(S(i,j,k,6)+obt*divg_th(i,j,k))
!          S_S =  2.d0*S(i,j,k,2)*S(i,j,k,2) + &
!                 2.d0*S(i,j,k,3)*S(i,j,k,3) + & 
!                 2.d0*S(i,j,k,5)*S(i,j,k,5) + &
!                 2.d0*S(i,j,k,1)*S(i,j,k,1) + &
!                 2.d0*S(i,j,k,4)*S(i,j,k,4) + &
!                 2.d0*S(i,j,k,6)*S(i,j,k,6) 
!          !
!          tt_com = 2.d0*obt*divg_th(i,j,k)*divg_th(i,j,k) + ( & ! mflux doesn't contribute to compressibility
!                   (u_p(i,j,k)+u_p(im,j,k))*dli(1)*(div_th(ip,j,k)-div_th(im,j,k))*0.25d0 + &
!                   (v_p(i,j,k)+v_p(i,jm,k))*dli(2)*(div_th(i,jp,k)-div_th(i,jm,k))*0.25d0 + &
!                   (w_p(i,j,k)+w_p(i,j,km))*dli(3)*(div_th(i,j,kp)-div_th(i,j,km))*0.25d0 )
!          tg_com = 2.d0*obt*divg_th(i,j,k)*divg_th(i,j,k) + ( &
!                   (u_p(i,j,k)+u_p(im,j,k))*dli(1)*(divg_th(ip,j,k)-divg_th(im,j,k))*0.25d0 + &
!                   (v_p(i,j,k)+v_p(i,jm,k))*dli(2)*(divg_th(i,jp,k)-divg_th(i,jm,k))*0.25d0 + &
!                   (w_p(i,j,k)+w_p(i,j,km))*dli(3)*(divg_th(i,j,kp)-divg_th(i,j,km))*0.25d0 )
!          eps_t  = eps_t  + 2.d0*(S_S - tt_com)*mu(i,j,k)
!          eps1_v = eps1_v + 2.d0*(S_S         )*vof(i,j,k)*mu1
!          eps2_v = eps2_v + 2.d0*(S_S - tg_com)*(1.d0-vof(i,j,k))*mu_gas(i,j,k)
!          eps2_m = eps2_m + 2.d0*(S_S - tg_com)*(1.d0-vof(i,j,k))*mu_gas(i,j,k)*rho_gas(i,j,k)
!          !
!          ! c2. production
!          !
!          fx_hit = f0_t*(abc(1)*sin(k0_freq*zc) + abc(3)*cos(k0_freq*yc))!*rhox/rhox ! to revise this rhox/rhox
!          fy_hit = f0_t*(abc(2)*sin(k0_freq*xc) + abc(1)*cos(k0_freq*zc))!*rhox/rhoy ! to revise this rhox/rhox
!          fz_hit = f0_t*(abc(3)*sin(k0_freq*yc) + abc(2)*cos(k0_freq*xc))!*rhox/rhoz ! to revise this rhox/rhox
!          t_prod = 0.5d0*(u_p(i,j,k)+u_p(im,j,k))*fx_hit + &
!                   0.5d0*(v_p(i,j,k)+v_p(i,jm,k))*fy_hit + &
!                   0.5d0*(w_p(i,j,k)+w_p(i,j,km))*fz_hit
!          prd_t  = prd_t  + t_prod*rho(i,j,k)
!          prd1_v = prd1_v + t_prod*vof(i,j,k)*rho1
!          prd2_v = prd2_v + t_prod*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          prd2_m = prd2_m + t_prod*(1.d0-vof(i,j,k))*rho_gas(i,j,k)**2
!          !
!          ! c3. kinetic energy
!          !
!          ke_t  = ke_t  + rho(i,j,k)*&
!          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
!          ke1_v = ke1_v + rho1*vof(i,j,k)*&
!          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
!          ke2_v = ke2_v + (rho_gas(i,j,k)   )*(1.d0-vof(i,j,k))*&
!          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
!          ke2_m = ke2_m + (rho_gas(i,j,k)**2)*(1.d0-vof(i,j,k))*&
!          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
!          !
!          ! c4. external forces
!          !
!          Psi_nu = Psi_nu + sigmaca*kappa(i,j,k)*( &
!                                    dli(1)*(0.5d0*(vof(ip,j,k)-vof(im,j,k)))*0.5d0*(u_p(i,j,k)+u_p(im,j,k)) + & 
!                                    dli(2)*(0.5d0*(vof(i,jp,k)-vof(i,jm,k)))*0.5d0*(v_p(i,j,k)+v_p(i,jm,k)) + & 
!                                    dli(3)*(0.5d0*(vof(i,j,kp)-vof(i,j,km)))*0.5d0*(w_p(i,j,k)+w_p(i,j,km))   &
!                                                 ) 
!          !
!          ! c5. pressure power due to volume expansion 
!          !
!#ifdef VAP_MASS
!          Psi_mf = Psi_mf + p(i,j,k)*div_th(i,j,k)
!#else
!          Psi_mf = 0.d0
!#endif
!          !
!          ! c6. pressure power due to volume expansion in the gas phase
!          !
!#ifdef LOW_MACH
!          Psi_gf_v = Psi_gf_v + p(i,j,k)*divg_th(i,j,k)*(1.d0-vof(i,j,k))
!          Psi_gf_m = Psi_gf_m + p(i,j,k)*divg_th(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!#else
!          Psi_gf_v = 0.d0
!          Psi_gf_m = 0.d0
!#endif
!          !
!          ! c7. pressure power (contribution per phase)
!          !
!          dupx = ( u_p(i,j,k)*0.5d0*(p(ip,j,k)+p(i,j,k))-u_p(im,j,k)*0.5d0*(p(i,j,k)+p(im,j,k)) )*dli(1)
!          dupy = ( v_p(i,j,k)*0.5d0*(p(i,jp,k)+p(i,j,k))-v_p(i,jm,k)*0.5d0*(p(i,j,k)+p(i,jm,k)) )*dli(2)
!          dupz = ( w_p(i,j,k)*0.5d0*(p(i,j,kp)+p(i,j,k))-w_p(i,j,km)*0.5d0*(p(i,j,k)+p(i,j,km)) )*dli(3)
!          !
!          Tp1_v = Tp1_v + (dupx+dupy+dupz)*vof(i,j,k) 
!          Tp2_v = Tp2_v + (dupx+dupy+dupz)*(1.d0-vof(i,j,k))
!          Tp2_m = Tp2_m + (dupx+dupy+dupz)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          !
!          ! c8. viscous power
!          !
!          duTxx = 0.500d0*(u_p(i ,j,k)*(S(i,j,k,1)+S(ip,j,k,1))*(mu(i,j,k)+mu(ip,j,k))-&
!                           u_p(im,j,k)*(S(i,j,k,1)+S(im,j,k,1))*(mu(i,j,k)+mu(im,j,k)))*dli(1)
!          duTxy = 0.125d0*(u_p(i,j,k)+u_p(i,jp,k)+u_p(im,j,k)+u_p(im,jp,k))*&
!                          (S(i,jp,k,2)+S(i,j,k,2))*(mu(i,jp,k)+mu(i,j,k))*dli(2)-&
!                  0.125d0*(u_p(i,j,k)+u_p(i,jm,k)+u_p(im,j,k)+u_p(im,jm,k))*&
!                          (S(i,jm,k,2)+S(i,j,k,2))*(mu(i,j,k)+mu(i,jm,k))*dli(2)  
!          duTxz = 0.125d0*(u_p(i,j,k)+u_p(i,j,kp)+u_p(im,j,k)+u_p(im,j,kp))*&
!                          (S(i,j,kp,3)+S(i,j,k,3))*(mu(i,j,kp)+mu(i,j,k))*dli(3)-&
!                  0.125d0*(u_p(i,j,k)+u_p(i,j,km)+u_p(im,j,k)+u_p(im,j,km))*&
!                          (S(i,j,k,3)+S(i,j,km,3))*(mu(i,j,k)+mu(i,j,km))*dli(3)
!          !
!          dvTyx = 0.125d0*(v_p(i,j,k)+v_p(ip,j,k)+v_p(i,jm,k)+v_p(ip,jm,k))*&
!                          (S(ip,j,k,2)+S(i,j,k,2))*(mu(ip,j,k)+mu(i,j,k))*dli(1)-&
!                  0.125d0*(v_p(i,j,k)+v_p(im,j,k)+v_p(i,jm,k)+v_p(im,jm,k))*&
!                          ((S(i,j,k,2)+S(im,j,k,2))*(mu(i,j,k)+mu(im,j,k)))*dli(1)
!          dvTyy = 0.500d0*(v_p(i,j ,k)*(S(i,j,k,4)+S(i,jp,k,4))*(mu(i,j,k)+mu(i,jp,k))-&
!                           v_p(i,jm,k)*(S(i,j,k,4)+S(i,jm,k,4))*(mu(i,j,k)+mu(i,jm,k)))*dli(2)
!          dvTyz = 0.125d0*(v_p(i,j,k)+v_p(i,j,kp)+v_p(i,jm,k)+v_p(i,jm,kp))*&
!                          (S(i,j,k,5)+S(i,j,kp,5))*(mu(i,j,k)+mu(i,j,kp))*dli(3)-&
!                  0.125d0*(v_p(i,j,k)+v_p(i,j,km)+v_p(i,jm,k)+v_p(i,jm,km))*&
!                          (S(i,j,k,5)+S(i,j,km,5))*(mu(i,j,k)+mu(i,j,km))*dli(3)
!          !
!          dwTzx = 0.125d0*(w_p(i,j,k)+w_p(ip,j,k)+w_p(i,j,km)+w_p(ip,j,km))*&
!                          (S(ip,j,k,3)+S(i,j,k,3))*(mu(ip,j,k)+mu(i,j,k))*dli(1)-&
!                  0.125d0*(w_p(i,j,k)+w_p(im,j,k)+w_p(i,j,km)+w_p(im,j,km))*& 
!                          (S(i,j,k,3)+S(im,j,k,3))*(mu(i,j,k)+mu(im,j,k))*dli(1)
!          dwTzy = 0.125d0*(w_p(i,j,k)+w_p(i,jp,k)+w_p(i,j,km)+w_p(i,jp,km))*&
!                          (S(i,j,k,5)+S(i,jp,k,5))*(mu(i,j,k)+mu(i,jp,k))*dli(2)-&
!                  0.125d0*(w_p(i,j,k)+w_p(i,jm,k)+w_p(i,j,km)+w_p(i,jm,km))*&    
!                          (S(i,j,k,5)+S(i,jm,k,5))*(mu(i,j,k)+mu(i,jm,k))*dli(2)
!          dwTzz = 0.500d0*(w_p(i,j,k )*(S(i,j,k,6)+S(i,j,kp,6))*(mu(i,j,k)+mu(i,j,kp))-&
!                           w_p(i,j,km)*(S(i,j,k,6)+S(i,j,km,6))*(mu(i,j,k)+mu(i,j,km)))*dli(3)
!          !
!          Tnu1_v = Tnu1_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*vof(i,j,k)
!          Tnu2_v = Tnu2_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*(1.d0-vof(i,j,k))
!          Tnu2_m = Tnu2_m + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
!          !
!        enddo
!      enddo
!    enddo
!    !
!    call mpi_allreduce(MPI_IN_PLACE,E_E   ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
!    call mpi_allreduce(MPI_IN_PLACE,ens1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
!    call mpi_allreduce(MPI_IN_PLACE,ens2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
!    call mpi_allreduce(MPI_IN_PLACE,ens2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
!    E_E    = E_E   /vol
!    ens1_v = ens1_v/vol1
!    ens2_v = ens2_v/vol2
!    ens2_m = ens2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,ur_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,vr_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,wr_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,ur1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,vr1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,wr1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,ur2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,vr2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,wr2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,ur2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,vr2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    call mpi_allreduce(MPI_IN_PLACE,wr2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
!    ur_t  = sqrt(ur_t/vol)
!    vr_t  = sqrt(vr_t/vol)
!    wr_t  = sqrt(wr_t/vol)
!    ur1_v = sqrt(ur1_v/vol1)
!    vr1_v = sqrt(vr1_v/vol1)
!    wr1_v = sqrt(wr1_v/vol1)
!    ur2_v = sqrt(ur2_v/vol2)
!    vr2_v = sqrt(vr2_v/vol2)
!    wr2_v = sqrt(wr2_v/vol2)
!    ur2_m = sqrt(ur2_m/mass_g)
!    vr2_m = sqrt(vr2_m/mass_g)
!    wr2_m = sqrt(wr2_m/mass_g)
!    !
!    call mpi_allreduce(MPI_IN_PLACE,um_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,vm_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,wm_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,um1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,vm1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,wm1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,um2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,vm2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,wm2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,um2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,vm2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    call mpi_allreduce(MPI_IN_PLACE,wm2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
!    um_t  = um_t/vol
!    vm_t  = vm_t/vol
!    wm_t  = wm_t/vol
!    um1_v = um1_v/vol1
!    vm1_v = vm1_v/vol1
!    wm1_v = wm1_v/vol1
!    um2_v = um2_v/vol2
!    vm2_v = vm2_v/vol2
!    wm2_v = wm2_v/vol2
!    um2_m = um2_m/mass_g
!    vm2_m = vm2_m/mass_g
!    wm2_m = wm2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,eps_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
!    call mpi_allreduce(MPI_IN_PLACE,eps1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
!    call mpi_allreduce(MPI_IN_PLACE,eps2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
!    call mpi_allreduce(MPI_IN_PLACE,eps2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
!    eps_t  = eps_t /vol1
!    eps2_v = eps2_v/vol2
!    eps2_m = eps2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,prd_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    call mpi_allreduce(MPI_IN_PLACE,prd1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    call mpi_allreduce(MPI_IN_PLACE,prd2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    c
!    eps1_v = eps1_v/vol1
!    eps2_v = eps2_v/vol2
!    eps2_m = eps2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,prd_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    call mpi_allreduce(MPI_IN_PLACE,prd1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    call mpi_allreduce(MPI_IN_PLACE,prd2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    call mpi_allreduce(MPI_IN_PLACE,prd2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
!    prd_t  = prd_t /vol
!    prd1_v = prd1_v/vol1
!    prd2_v = prd2_v/vol2
!    prd2_m = prd2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
!    call mpi_allreduce(MPI_IN_PLACE,ke1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
!    call mpi_allreduce(MPI_IN_PLACE,ke2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
!    call mpi_allreduce(MPI_IN_PLACE,ke2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
!    ke_t  = ke_t /vol
!    ke1_v = ke1_v/vol1
!    ke2_v = ke2_v/vol2
!    ke2_m = ke2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,Psi_mf  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c4
!    call mpi_allreduce(MPI_IN_PLACE,Psi_nu  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c5
!    call mpi_allreduce(MPI_IN_PLACE,Psi_gf_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c6-a
!    call mpi_allreduce(MPI_IN_PLACE,Psi_gf_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c6-b
!    Psi_mf   = Psi_mf/vol
!    Psi_nu   = Psi_nu/vol
!    Psi_gf_v = Psi_gf_v/vol2
!    Psi_gf_m = Psi_gf_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,Tp1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
!    call mpi_allreduce(MPI_IN_PLACE,Tp2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
!    call mpi_allreduce(MPI_IN_PLACE,Tp2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
!    Tp1_v  = Tp1_v/vol1
!    Tp2_v  = Tp2_v/vol2
!    Tp2_m  = Tp2_m/mass_g
!    !
!    call mpi_allreduce(MPI_IN_PLACE,Tnu1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
!    call mpi_allreduce(MPI_IN_PLACE,Tnu2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
!    call mpi_allreduce(MPI_IN_PLACE,Tnu2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
!    Tnu1_v = Tnu1_v/vol1
!    Tnu2_v = Tnu2_v/vol2
!    Tnu2_m = Tnu2_m/mass_g
!    !
!    if(myid.eq.0) then
!      !
!      ! overall budget
!      !
!      open(94,file=trim(datapos)//'budget_overall.out',position='append')
!      write(94,'(16E15.7)') 1.d0*istep,time,um_t,vm_t,wm_t,ur_t,vr_t,wr_t, &
!                            ens_t,ke_t,prd_t,eps_t,psi_nu,psi_mf,psi_gf_v,psi_gf_m 
!      close(94)
!      !
!      ! phase 1 - vol
!      !
!      open(94,file=trim(datapos)//'budget_phase1_v.out',position='append')
!      write(94,'(14E15.7)') 1.d0*istep,time,um1_v,vm1_v,wm1_v,ur1_v,vr1_v,wr1_v, &
!                            ens1_v,ke1_v,prd1_v,eps1_v,Tp1_v,Tnu1_v 
!      close(94)
!      !
!      ! phase 2 - vol
!      !
!      open(94,file=trim(datapos)//'budget_phase2_v.out',position='append')
!      write(94,'(14E15.7)') 1.d0*istep,time,um2_v,vm2_v,wm2_v,ur2_v,vr2_v,wr2_v, &
!                            ens2_v,ke2_v,prd2_v,eps2_v,Tp2_v,Tnu2_v 
!      close(94)
!      !
!      ! phase 2 - mass
!      !
!      open(94,file=trim(datapos)//'budget_phase2_m.out',position='append')
!      write(94,'(14E15.7)') 1.d0*istep,time,um2_m,vm2_m,wm2_m,ur2_m,vr2_m,wr2_m, &
!                            ens2_m,ke2_m,prd2_m,eps2_m,Tp2_m,Tnu2_m 
!      close(94)
!      !
!    endif
!    !
!    return
!  end subroutine energy_balance
! 
!end module mod_output
 
 subroutine energy_balance(n,dli,nh_d,nh_u,dzc,dzf,time,istep,rho,mu,vof,kappa,p, &
                            rho_gas,mu_gas,div_th,divg_th,u,v,w)
    !
    ! energy balance in physical space
    !
    use mod_param     , only: datapos,sigmaca,cbcvof,bcvof,dims,rho1,mu1, &
                              f0_t,abc,k0_freq,pi,ng,mu2!,rho2
    use mod_common_mpi, only: ierr
    !
    implicit none
    !
    integer, intent(in  ), dimension(3)                       :: n
    real(8), intent(in  ), dimension(3)                       :: dli
    integer, intent(in  )                                     :: nh_d,nh_u
    real(8), intent(in  ), dimension(1-nh_d:)                 :: dzc,dzf
    real(8), intent(in  )                                     :: time
    integer, intent(in  )                                     :: istep
    real(8), intent(in  ), dimension(     0:,     0:,     0:) :: rho,mu,vof,kappa,p
    real(8), intent(in  ), dimension(     0:,     0:,     0:) :: rho_gas,mu_gas,div_th,divg_th
    real(8), intent(in  ), dimension(1-nh_u:,1-nh_u:,1-nh_u:) :: u,v,w
    !
    integer :: i,j,k,ip,im,jp,jm,kp,km,q
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)   :: u_p,v_p,w_p
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1,6) :: S      
    real(8), dimension(3) :: dl
    real(8) :: u_avg,v_avg,w_avg,S_S,E_E, &
               duTxx,duTxy,duTxz,dvTyx,dvTyy,dvTyz,dwTzx,dwTzy,dwTzz, & 
               tt_com,tg_com,t_prod
    real(8) :: ur_t,vr_t,wr_t,ur1_v,vr1_v,wr1_v,ur2_v,vr2_v,wr2_v,ur2_m,vr2_m,wr2_m, &
               um_t,vm_t,wm_t,um1_v,vm1_v,wm1_v,um2_v,vm2_v,wm2_v,um2_m,vm2_m,wm2_m, &
               eps_t ,eps1_v,eps2_v,eps2_m, &
               prd_t ,prd1_v,prd2_v,prd2_m, & 
               ke_t  ,ke1_v ,ke2_v ,ke2_m , & 
               ens_t ,ens1_v,ens2_v,ens2_m, &
               Psi_mf,Psi_nu,Psi_gf_v,Psi_gf_m, &
               Tp1_v ,Tp2_v ,Tp2_m , &
               Tnu1_v,Tnu2_v,Tnu2_m, &
               vol1,vol2,mass_g,vol,obt
    real(8) :: VGRux , VGRuy , VGRuz , &
               VGRvx , VGRvy , VGRvz , &
               VGRwx , VGRwy , VGRwz , &
               dupx  , dupy  , dupz  , &
               fx_hit, fy_hit, fz_hit, &
               xc    , yc    , zc
    !real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1)   :: S_S_A,eps1_v_A,ke1_v_A
    real(8) :: S_S_A,eps1_v_A,ke1_v_A
    real(8) :: Re_L 
    !
    vol    = 1.d0*(n(1)*dims(1))*(n(2)*dims(2))*n(3)
    dl(:)  = 1.d0/dli(:)
    obt    = 1.d0/3.d0
    !
    E_E    = 0.d0 
    ens_t  = 0.d0
    ens1_v = 0.d0
    ens2_v = 0.d0
    ens2_m = 0.d0
    !
    ur_t   = 0.d0 
    vr_t   = 0.d0
    wr_t   = 0.d0
    ur1_v  = 0.d0
    vr1_v  = 0.d0
    wr1_v  = 0.d0
    ur2_v  = 0.d0
    vr2_v  = 0.d0
    wr2_v  = 0.d0
    ur2_m  = 0.d0
    vr2_m  = 0.d0
    wr2_m  = 0.d0
    !
    um_t   = 0.d0 
    vm_t   = 0.d0
    wm_t   = 0.d0
    um1_v  = 0.d0
    vm1_v  = 0.d0
    wm1_v  = 0.d0
    um2_v  = 0.d0
    vm2_v  = 0.d0
    wm2_v  = 0.d0
    um2_m  = 0.d0
    vm2_m  = 0.d0
    wm2_m  = 0.d0
    !
    eps_t  = 0.d0 
    eps1_v = 0.d0
    eps2_v = 0.d0
    eps2_m = 0.d0
    !
    prd_t  = 0.d0
    prd1_v = 0.d0
    prd2_v = 0.d0
    prd2_m = 0.d0
    !
    ke_t   = 0.d0
    ke1_v  = 0.d0
    ke2_v  = 0.d0
    ke2_m  = 0.d0
    !
    Psi_mf   = 0.d0 
    Psi_nu   = 0.d0
    Psi_gf_v = 0.d0
    Psi_gf_m = 0.d0
    !
    Tp1_v  = 0.d0 
    Tp2_v  = 0.d0
    Tp2_m  = 0.d0
    !
    Tnu1_v = 0.d0 
    Tnu2_v = 0.d0
    Tnu2_m = 0.d0
    ! 
    ! following Dodd and Ferrante, we compute (u_p,v_p,w_p)
    !
    u_avg = 0.d0
    v_avg = 0.d0
    w_avg = 0.d0
    !
    vol1   = 0.d0
    vol2   = 0.d0
    mass_g = 0.d0
    !
    Re_L   = 0.d0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! average velocity
          !
          u_avg = u_avg + u(i,j,k)*dl(1)*dl(2)*dl(3)
          v_avg = v_avg + v(i,j,k)*dl(1)*dl(2)*dl(3)
          w_avg = w_avg + w(i,j,k)*dl(1)*dl(2)*dl(3)
          !
          ! volume of the two phases and mass of the gas
          !
          vol1   = vol1 + vof(i,j,k)
          vol2   = vol2 + (1.d0-vof(i,j,k))
          mass_g = mass_g + rho_gas(i,j,k)*(1.d0-vof(i,j,k))
          !
        enddo
      enddo
    enddo 
    call mpi_allreduce(MPI_IN_PLACE,u_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,v_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,w_avg ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    u_avg = u_avg/(lx*ly*lz)
    v_avg = v_avg/(lx*ly*lz)
    w_avg = w_avg/(lx*ly*lz)
    !
    call mpi_allreduce(MPI_IN_PLACE,vol1  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,vol2  ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mass_g,1,mpi_real8,mpi_sum,comm_cart,ierr)
    vol1   = vol1   + 1e-16 ! to avoid division by 0
    vol2   = vol2   + 1e-16 ! to avoid division by 0
    mass_g = mass_g + 1e-16 ! to avoid division by 0
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          u_p(i,j,k) = u(i,j,k)-u_avg
          v_p(i,j,k) = v(i,j,k)-v_avg
          w_p(i,j,k) = w(i,j,k)-w_avg
        enddo
      enddo
    enddo 
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
          ! a. MEAN KINETIC ENERGY BALANCE --> part 1
          !
          ! a1. calculation and storage (for the next loop) of the 
          !     strain rate tensor S_ij
          !
          ! along x
          !
          VGRux = ( u_p(i,j,k) - u_p(im,j,k) )*dli(1)
          VGRuy = 0.25*( u_p(i,j,k) + u_p(i,jp,k) + u_p(im,j,k) + u_p(im,jp,k) )*dli(2) - &
                  0.25*( u_p(i,j,k) + u_p(i,jm,k) + u_p(im,j,k) + u_p(im,jm,k) )*dli(2)  
          VGRuz = 0.25*( u_p(i,j,k) + u_p(i,j,kp) + u_p(im,j,k) + u_p(im,j,kp) )*dli(3) - &
                  0.25*( u_p(i,j,k) + u_p(i,j,km) + u_p(im,j,k) + u_p(im,j,km) )*dli(3)
          !
          ! along y
          !  
          VGRvx = 0.25*( v_p(i,j,k) + v_p(ip,j,k) + v_p(i,jm,k) + v_p(ip,jm,k) )*dli(1) - &
                  0.25*( v_p(i,j,k) + v_p(im,j,k) + v_p(i,jm,k) + v_p(im,jm,k) )*dli(1)
          VGRvy = ( v_p(i,j,k) - v_p(i,jm,k) )*dli(2)
          VGRvz = 0.25*( v_p(i,j,k) + v_p(i,j,kp) + v_p(i,jm,k) + v_p(i,jm,kp) )*dli(3) - &
                  0.25*( v_p(i,j,k) + v_p(i,j,km) + v_p(i,jm,k) + v_p(i,jm,km) )*dli(3)
          !
          ! along z
          ! 
          VGRwx = 0.25*( w_p(i,j,k) + w_p(ip,j,k) + w_p(i,j,km) + w_p(ip,j,km) )*dli(1) - &
                  0.25*( w_p(i,j,k) + w_p(im,j,k) + w_p(i,j,km) + w_p(im,j,km) )*dli(1)
          VGRwy = 0.25*( w_p(i,j,k) + w_p(i,jp,k) + w_p(i,j,km) + w_p(i,jp,km) )*dli(2) - &
                  0.25*( w_p(i,j,k) + w_p(i,jm,k) + w_p(i,j,km) + w_p(i,jm,km) )*dli(2)
          VGRwz = ( w_p(i,j,k) - w_p(i,j,km) )*dli(3)
          !
          ! incompressible part
          ! 
          !-----S11 = dU/dx
          S(i,j,k,1) = VGRux
          !-----S12 = 0.5*(dU/dy+dV/dx) = S21
          S(i,j,k,2) = 0.5 * (VGRuy + VGRvx)
          !-----S13 = 0.5*(dU/dz+dW/dx) = S31
          S(i,j,k,3) = 0.5 * (VGRuz + VGRwx)
          !-----S22 = dV/dy
          S(i,j,k,4) = VGRvy
          !-----S23 = 0.5*(dV/dz+dW/dy) = S32
          S(i,j,k,5) = 0.5 * (VGRvz + VGRwy)
          !-----S33 = dW/dz
          S(i,j,k,6) = VGRwz
          !
          ! add the compressible contribution
          ! 
          !S(i,j,k,1) = S(i,j,k,1) - obt*divg_th(i,j,k)
          !S(i,j,k,4) = S(i,j,k,4) - obt*divg_th(i,j,k)
          !S(i,j,k,6) = S(i,j,k,6) - obt*divg_th(i,j,k)
          !
          ! a2. enstrophy
          !
          E_E    = (VGRwy-VGRvz)**2. + (VGRuz-VGRwx)**2. + (VGRvx-VGRuy)**2. 
          ens_t  = ens_t  + E_E
          ens1_v = ens1_v + E_E*vof(i,j,k)
          ens2_v = ens2_v + E_E*(1.d0-vof(i,j,k))
          ens2_m = ens2_m + E_E*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    ! 
    do q=1,6
      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,S(:,:,:,q))
    enddo
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      zc = (k-0.5d0)*dl(3)/lz*2.d0*pi
      do j=1,n(2)
        jp = j+1
        jm = j-1
        yc = (j+coord(2)*n(2)-0.5d0)*dl(2)/ly*2.d0*pi
        do i=1,n(1)
          ip = i+1
          im = i-1
          xc = (i+coord(1)*n(1)-0.5d0)*dl(1)/lx*2.d0*pi
          !
          ! b1. RMS for one-fluid and two phases
          !
          ur_t = ur_t + (u_p(i,j,k))**2
          vr_t = vr_t + (v_p(i,j,k))**2
          wr_t = wr_t + (w_p(i,j,k))**2
          ! 
          ur1_v = ur1_v + ((u_p(i,j,k))**2)*vof(i,j,k)
          vr1_v = vr1_v + ((v_p(i,j,k))**2)*vof(i,j,k)
          wr1_v = wr1_v + ((w_p(i,j,k))**2)*vof(i,j,k)
          !
          ur2_v = ur2_v + ((u_p(i,j,k))**2)*(1.d0-vof(i,j,k))
          vr2_v = vr2_v + ((v_p(i,j,k))**2)*(1.d0-vof(i,j,k))
          wr2_v = wr2_v + ((w_p(i,j,k))**2)*(1.d0-vof(i,j,k))
          ur2_m = ur2_m + ((u_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          vr2_m = vr2_m + ((v_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          wr2_m = wr2_m + ((w_p(i,j,k))**2)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          !
          ! b2. mean velocities for one-fluid and two phases
          !
          um_t = um_t + u_p(i,j,k)
          vm_t = vm_t + v_p(i,j,k)
          wm_t = wm_t + w_p(i,j,k)
          ! 
          um1_v = um1_v + u_p(i,j,k)*vof(i,j,k)
          vm1_v = vm1_v + v_p(i,j,k)*vof(i,j,k)
          wm1_v = wm1_v + w_p(i,j,k)*vof(i,j,k)
          !
          um2_v = um2_v + u_p(i,j,k)*(1.d0-vof(i,j,k))
          vm2_v = vm2_v + v_p(i,j,k)*(1.d0-vof(i,j,k))
          wm2_v = wm2_v + w_p(i,j,k)*(1.d0-vof(i,j,k))
          um2_m = um2_m + u_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          vm2_m = vm2_m + v_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          wm2_m = wm2_m + w_p(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          !
          ! c. MEAN KINETIC ENERGY BALANCE --> part 2
          !
          ! c1. dissipation
          !
          !S_S and S_S_A should be the same for incompressible flow
          !
          S_S =     (S(i,j,k,1)+obt*divg_th(i,j,k))*(S(i,j,k,1)+obt*divg_th(i,j,k)) + &
                2.d0*S(i,j,k,2)*S(i,j,k,2) + &
                2.d0*S(i,j,k,3)*S(i,j,k,3) + & 
                    (S(i,j,k,4)+obt*divg_th(i,j,k))*(S(i,j,k,4)+obt*divg_th(i,j,k)) + &
                2.d0*S(i,j,k,5)*S(i,j,k,5) + &
                    (S(i,j,k,6)+obt*divg_th(i,j,k))*(S(i,j,k,6)+obt*divg_th(i,j,k))
         
          !
          S_S_A =  2.d0*S(i,j,k,2)*S(i,j,k,2) + &
                   2.d0*S(i,j,k,3)*S(i,j,k,3) + & 
                   2.d0*S(i,j,k,5)*S(i,j,k,5) + &
                        S(i,j,k,1)*S(i,j,k,1) + &
                        S(i,j,k,4)*S(i,j,k,4) + &
                        S(i,j,k,6)*S(i,j,k,6) 
          !
          tt_com = 2.d0*obt*divg_th(i,j,k)*divg_th(i,j,k) + ( & ! mflux doesn't contribute to compressibility
                   (u_p(i,j,k)+u_p(im,j,k))*dli(1)*(div_th(ip,j,k)-div_th(im,j,k))*0.25d0 + &
                   (v_p(i,j,k)+v_p(i,jm,k))*dli(2)*(div_th(i,jp,k)-div_th(i,jm,k))*0.25d0 + &
                   (w_p(i,j,k)+w_p(i,j,km))*dli(3)*(div_th(i,j,kp)-div_th(i,j,km))*0.25d0 )
          tg_com = 2.d0*obt*divg_th(i,j,k)*divg_th(i,j,k) + ( &
                   (u_p(i,j,k)+u_p(im,j,k))*dli(1)*(divg_th(ip,j,k)-divg_th(im,j,k))*0.25d0 + &
                   (v_p(i,j,k)+v_p(i,jm,k))*dli(2)*(divg_th(i,jp,k)-divg_th(i,jm,k))*0.25d0 + &
                   (w_p(i,j,k)+w_p(i,j,km))*dli(3)*(divg_th(i,j,kp)-divg_th(i,j,km))*0.25d0 )
          !
          !eps_t and eps1_v_A should lead to the same value for an incompressible flow 
          !eps_t will be used for avergaed k-epsilon method and the corresponding Re_lambda will be printed at the budget_overall
          !eps1_v_A will be used for averaged Re_lambda method
          !
          !eps_t      = eps_t  + 2.d0*(S_S - tt_com)*mu(i,j,k)
          eps_t      = eps_t  + 2.d0*(S_S_A)*mu(i,j,k)
          eps1_v     = eps1_v + 2.d0*(S_S         )*vof(i,j,k)*mu1
          eps1_v_A   = 2.d0*(S_S_A       )*mu_gas(i,j,k)
          eps2_v     = eps2_v + 2.d0*(S_S - tg_com)*(1.d0-vof(i,j,k))*mu_gas(i,j,k)
          eps2_m     = eps2_m + 2.d0*(S_S - tg_com)*(1.d0-vof(i,j,k))*mu_gas(i,j,k)*rho_gas(i,j,k)
          !
          ! c2. production
          !
          fx_hit = f0_t*(abc(1)*sin(k0_freq*zc) + abc(3)*cos(k0_freq*yc))!*rhox/rhox ! to revise this rhox/rhox
          fy_hit = f0_t*(abc(2)*sin(k0_freq*xc) + abc(1)*cos(k0_freq*zc))!*rhox/rhoy ! to revise this rhox/rhox
          fz_hit = f0_t*(abc(3)*sin(k0_freq*yc) + abc(2)*cos(k0_freq*xc))!*rhox/rhoz ! to revise this rhox/rhox
          t_prod = 0.5d0*(u_p(i,j,k)+u_p(im,j,k))*fx_hit + &
                   0.5d0*(v_p(i,j,k)+v_p(i,jm,k))*fy_hit + &
                   0.5d0*(w_p(i,j,k)+w_p(i,j,km))*fz_hit
          prd_t  = prd_t  + t_prod*rho(i,j,k)
          prd1_v = prd1_v + t_prod*vof(i,j,k)*rho1
          prd2_v = prd2_v + t_prod*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          prd2_m = prd2_m + t_prod*(1.d0-vof(i,j,k))*rho_gas(i,j,k)**2
          !
          ! c3. kinetic energy
          !
          !ke_t will be used for averaged k-epsilon method and the correspoding Re_lambda will be printed at the budget_overall
          !ke1_v_A will be used for averaged Re_lambda method
          !
          ke_t           =  ke_t  + rho(i,j,k)*&
          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
          ke1_v          =  ke1_v + rho1*vof(i,j,k)*&
          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
          ke1_v_A        =  rho_gas(i,j,k)*&
          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
          ke2_v          =  ke2_v + (rho_gas(i,j,k)   )*(1.d0-vof(i,j,k))*&
          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
          ke2_m          =  ke2_m + (rho_gas(i,j,k)**2)*(1.d0-vof(i,j,k))*&
          0.5d0*(0.25d0*(u_p(i,j,k)+u_p(im,j,k))**2 + 0.25d0*(v_p(i,j,k)+v_p(i,jm,k))**2 + 0.25d0*(w_p(i,j,k)+w_p(i,j,km))**2)
          !
          !
          !
          !computing Reynolds_lambda locally
          ! 
          Re_L = Re_L + sqrt(20/3*(ke1_v_A**2/(eps1_v_A*mu_gas(i,j,k)+1e-16)))  
          !
          !
          ! c4. external forces
          !
          Psi_nu = Psi_nu + sigmaca*kappa(i,j,k)*( &
                                    dli(1)*(0.5d0*(vof(ip,j,k)-vof(im,j,k)))*0.5d0*(u_p(i,j,k)+u_p(im,j,k)) + & 
                                    dli(2)*(0.5d0*(vof(i,jp,k)-vof(i,jm,k)))*0.5d0*(v_p(i,j,k)+v_p(i,jm,k)) + & 
                                    dli(3)*(0.5d0*(vof(i,j,kp)-vof(i,j,km)))*0.5d0*(w_p(i,j,k)+w_p(i,j,km))   &
                                                 ) 
          !
          ! c5. pressure power due to volume expansion 
          !
#ifdef VAP_MASS
          Psi_mf = Psi_mf + p(i,j,k)*div_th(i,j,k)
#else
          Psi_mf = 0.d0
#endif
          !
          ! c6. pressure power due to volume expansion in the gas phase
          !
#ifdef LOW_MACH
          Psi_gf_v = Psi_gf_v + p(i,j,k)*divg_th(i,j,k)*(1.d0-vof(i,j,k))
          Psi_gf_m = Psi_gf_m + p(i,j,k)*divg_th(i,j,k)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
#else
          Psi_gf_v = 0.d0
          Psi_gf_m = 0.d0
#endif
          !
          ! c7. pressure power (contribution per phase)
          !
          dupx = ( u_p(i,j,k)*0.5d0*(p(ip,j,k)+p(i,j,k))-u_p(im,j,k)*0.5d0*(p(i,j,k)+p(im,j,k)) )*dli(1)
          dupy = ( v_p(i,j,k)*0.5d0*(p(i,jp,k)+p(i,j,k))-v_p(i,jm,k)*0.5d0*(p(i,j,k)+p(i,jm,k)) )*dli(2)
          dupz = ( w_p(i,j,k)*0.5d0*(p(i,j,kp)+p(i,j,k))-w_p(i,j,km)*0.5d0*(p(i,j,k)+p(i,j,km)) )*dli(3)
          !
          Tp1_v = Tp1_v + (dupx+dupy+dupz)*vof(i,j,k) 
          Tp2_v = Tp2_v + (dupx+dupy+dupz)*(1.d0-vof(i,j,k))
          Tp2_m = Tp2_m + (dupx+dupy+dupz)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          !
          ! c8. viscous power
          !
          duTxx = 0.500d0*(u_p(i ,j,k)*(S(i,j,k,1)+S(ip,j,k,1))*(mu(i,j,k)+mu(ip,j,k))-&
                           u_p(im,j,k)*(S(i,j,k,1)+S(im,j,k,1))*(mu(i,j,k)+mu(im,j,k)))*dli(1)
          duTxy = 0.125d0*(u_p(i,j,k)+u_p(i,jp,k)+u_p(im,j,k)+u_p(im,jp,k))*&
                          (S(i,jp,k,2)+S(i,j,k,2))*(mu(i,jp,k)+mu(i,j,k))*dli(2)-&
                  0.125d0*(u_p(i,j,k)+u_p(i,jm,k)+u_p(im,j,k)+u_p(im,jm,k))*&
                          (S(i,jm,k,2)+S(i,j,k,2))*(mu(i,j,k)+mu(i,jm,k))*dli(2)  
          duTxz = 0.125d0*(u_p(i,j,k)+u_p(i,j,kp)+u_p(im,j,k)+u_p(im,j,kp))*&
                          (S(i,j,kp,3)+S(i,j,k,3))*(mu(i,j,kp)+mu(i,j,k))*dli(3)-&
                  0.125d0*(u_p(i,j,k)+u_p(i,j,km)+u_p(im,j,k)+u_p(im,j,km))*&
                          (S(i,j,k,3)+S(i,j,km,3))*(mu(i,j,k)+mu(i,j,km))*dli(3)
          !
          dvTyx = 0.125d0*(v_p(i,j,k)+v_p(ip,j,k)+v_p(i,jm,k)+v_p(ip,jm,k))*&
                          (S(ip,j,k,2)+S(i,j,k,2))*(mu(ip,j,k)+mu(i,j,k))*dli(1)-&
                  0.125d0*(v_p(i,j,k)+v_p(im,j,k)+v_p(i,jm,k)+v_p(im,jm,k))*&
                          ((S(i,j,k,2)+S(im,j,k,2))*(mu(i,j,k)+mu(im,j,k)))*dli(1)
          dvTyy = 0.500d0*(v_p(i,j ,k)*(S(i,j,k,4)+S(i,jp,k,4))*(mu(i,j,k)+mu(i,jp,k))-&
                           v_p(i,jm,k)*(S(i,j,k,4)+S(i,jm,k,4))*(mu(i,j,k)+mu(i,jm,k)))*dli(2)
          dvTyz = 0.125d0*(v_p(i,j,k)+v_p(i,j,kp)+v_p(i,jm,k)+v_p(i,jm,kp))*&
                          (S(i,j,k,5)+S(i,j,kp,5))*(mu(i,j,k)+mu(i,j,kp))*dli(3)-&
                  0.125d0*(v_p(i,j,k)+v_p(i,j,km)+v_p(i,jm,k)+v_p(i,jm,km))*&
                          (S(i,j,k,5)+S(i,j,km,5))*(mu(i,j,k)+mu(i,j,km))*dli(3)
          !
          dwTzx = 0.125d0*(w_p(i,j,k)+w_p(ip,j,k)+w_p(i,j,km)+w_p(ip,j,km))*&
                          (S(ip,j,k,3)+S(i,j,k,3))*(mu(ip,j,k)+mu(i,j,k))*dli(1)-&
                  0.125d0*(w_p(i,j,k)+w_p(im,j,k)+w_p(i,j,km)+w_p(im,j,km))*& 
                          (S(i,j,k,3)+S(im,j,k,3))*(mu(i,j,k)+mu(im,j,k))*dli(1)
          dwTzy = 0.125d0*(w_p(i,j,k)+w_p(i,jp,k)+w_p(i,j,km)+w_p(i,jp,km))*&
                          (S(i,j,k,5)+S(i,jp,k,5))*(mu(i,j,k)+mu(i,jp,k))*dli(2)-&
                  0.125d0*(w_p(i,j,k)+w_p(i,jm,k)+w_p(i,j,km)+w_p(i,jm,km))*&    
                          (S(i,j,k,5)+S(i,jm,k,5))*(mu(i,j,k)+mu(i,jm,k))*dli(2)
          dwTzz = 0.500d0*(w_p(i,j,k )*(S(i,j,k,6)+S(i,j,kp,6))*(mu(i,j,k)+mu(i,j,kp))-&
                           w_p(i,j,km)*(S(i,j,k,6)+S(i,j,km,6))*(mu(i,j,k)+mu(i,j,km)))*dli(3)
          !
          Tnu1_v = Tnu1_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*vof(i,j,k)
          Tnu2_v = Tnu2_v + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*(1.d0-vof(i,j,k))
          Tnu2_m = Tnu2_m + (duTxx+duTxy+duTxz+dvTyx+dvTyy+dvTyz+dwTzx+dwTzy+dwTzz)*(1.d0-vof(i,j,k))*rho_gas(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,E_E   ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
    call mpi_allreduce(MPI_IN_PLACE,ens1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
    call mpi_allreduce(MPI_IN_PLACE,ens2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
    call mpi_allreduce(MPI_IN_PLACE,ens2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! a2
    E_E    = E_E   /vol
    ens1_v = ens1_v/vol1
    ens2_v = ens2_v/vol2
    ens2_m = ens2_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,ur_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,ur1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,ur2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,ur2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,vr2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    call mpi_allreduce(MPI_IN_PLACE,wr2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b1
    ur_t  = sqrt(ur_t/vol)
    vr_t  = sqrt(vr_t/vol)
    wr_t  = sqrt(wr_t/vol)
    ur1_v = sqrt(ur1_v/vol1)
    vr1_v = sqrt(vr1_v/vol1)
    wr1_v = sqrt(wr1_v/vol1)
    ur2_v = sqrt(ur2_v/vol2)
    vr2_v = sqrt(vr2_v/vol2)
    wr2_v = sqrt(wr2_v/vol2)
    ur2_m = sqrt(ur2_m/mass_g)
    vr2_m = sqrt(vr2_m/mass_g)
    wr2_m = sqrt(wr2_m/mass_g)
    !
    call mpi_allreduce(MPI_IN_PLACE,um_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,um1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,um2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,um2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,vm2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    call mpi_allreduce(MPI_IN_PLACE,wm2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! b2
    um_t  = um_t/vol
    vm_t  = vm_t/vol
    wm_t  = wm_t/vol
    um1_v = um1_v/vol1
    vm1_v = vm1_v/vol1
    wm1_v = wm1_v/vol1
    um2_v = um2_v/vol2
    vm2_v = vm2_v/vol2
    wm2_v = wm2_v/vol2
    um2_m = um2_m/mass_g
    vm2_m = vm2_m/mass_g
    wm2_m = wm2_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,eps_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
    call mpi_allreduce(MPI_IN_PLACE,eps1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
    call mpi_allreduce(MPI_IN_PLACE,eps2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
    call mpi_allreduce(MPI_IN_PLACE,eps2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c1
    eps_t  = eps_t /vol
    eps1_v = eps1_v/vol1
    eps2_v = eps2_v/vol2
    eps2_m = eps2_m/mass_g
    !
    !call mpi_allreduce(MPI_IN_PLACE,prd_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    !call mpi_allreduce(MPI_IN_PLACE,prd1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    !call mpi_allreduce(MPI_IN_PLACE,prd2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    !
    !eps1_v = eps1_v/vol1
    !eps2_v = eps2_v/vol2
    !eps2_m = eps2_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,prd_t ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    call mpi_allreduce(MPI_IN_PLACE,prd1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    call mpi_allreduce(MPI_IN_PLACE,prd2_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    call mpi_allreduce(MPI_IN_PLACE,prd2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c2
    prd_t  = prd_t /vol
    prd1_v = prd1_v/vol1
    prd2_v = prd2_v/vol2
    prd2_m = prd2_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,ke_t  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,ke1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,ke2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,ke2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
    call mpi_allreduce(MPI_IN_PLACE,Re_L  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c3
    ke_t  = ke_t /vol
    ke1_v = ke1_v/vol1
    ke2_v = ke2_v/vol2
    ke2_m = ke2_m/mass_g
    ! Re_L  = Re_L/((ng(1))*(ng(2))*(ng(3)))
    Re_L  = Re_L/vol
    !
    call mpi_allreduce(MPI_IN_PLACE,Psi_mf  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c4
    call mpi_allreduce(MPI_IN_PLACE,Psi_nu  ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c5
    call mpi_allreduce(MPI_IN_PLACE,Psi_gf_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c6-a
    call mpi_allreduce(MPI_IN_PLACE,Psi_gf_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c6-b
    Psi_mf   = Psi_mf/vol
    Psi_nu   = Psi_nu/vol
    Psi_gf_v = Psi_gf_v/vol2
    Psi_gf_m = Psi_gf_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,Tp1_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
    call mpi_allreduce(MPI_IN_PLACE,Tp2_v ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
    call mpi_allreduce(MPI_IN_PLACE,Tp2_m ,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c7
    Tp1_v  = Tp1_v/vol1
    Tp2_v  = Tp2_v/vol2
    Tp2_m  = Tp2_m/mass_g
    !
    call mpi_allreduce(MPI_IN_PLACE,Tnu1_v,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,Tnu2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
    call mpi_allreduce(MPI_IN_PLACE,Tnu2_m,1,mpi_real8,mpi_sum,comm_cart,ierr) ! c8
    Tnu1_v = Tnu1_v/vol1
    Tnu2_v = Tnu2_v/vol2
    Tnu2_m = Tnu2_m/mass_g
    !
    if(myid.eq.0) then
      !
      ! overall budget
      !
      !the coloumn 18th corresponds to averaged Re_lambda   
      !the coloumn 17th corresponds to Re_lambda based on averaged k-epsilon
      !
      open(94,file=trim(datapos)//'budget_overall.out',position='append')
      write(94,'(18E15.7)') 1.d0*istep,time,um_t,vm_t,wm_t,ur_t,vr_t,wr_t, &
                            ens_t,ke_t,prd_t,eps_t,psi_nu,psi_mf,psi_gf_v,psi_gf_m, &
                            sqrt((20.0/3.0)*ke_t**2/(mu2*eps_t)),Re_L 
      close(94)
      !
      ! phase 1 - vol
      !
      open(94,file=trim(datapos)//'budget_phase1_v.out',position='append')
      write(94,'(14E15.7)') 1.d0*istep,time,um1_v,vm1_v,wm1_v,ur1_v,vr1_v,wr1_v, &
                            ens1_v,ke1_v,prd1_v,eps1_v,Tp1_v,Tnu1_v 
      close(94)
      !
      ! phase 2 - vol
      !
      open(94,file=trim(datapos)//'budget_phase2_v.out',position='append')
      write(94,'(14E15.7)') 1.d0*istep,time,um2_v,vm2_v,wm2_v,ur2_v,vr2_v,wr2_v, &
                            ens2_v,ke2_v,prd2_v,eps2_v,Tp2_v,Tnu2_v 
      close(94)
      !
      ! phase 2 - mass
      !
      open(94,file=trim(datapos)//'budget_phase2_m.out',position='append')
      write(94,'(14E15.7)') 1.d0*istep,time,um2_m,vm2_m,wm2_m,ur2_m,vr2_m,wr2_m, &
                            ens2_m,ke2_m,prd2_m,eps2_m,Tp2_m,Tnu2_m 
      close(94)
      !
    endif
    !
    return
  end subroutine energy_balance
  ! 
end module mod_output
