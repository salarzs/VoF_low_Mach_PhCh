module mod_tagging
  !
  use mpi
  use mod_common_mpi, only: myid, comm_Cart,ierr, coord,status, &
                            xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf, &
                            left,right,front,back,xhalo,yhalo
  use mod_param     , only: dl,dims,datadir_ta
  !
  implicit none
  !
  real(8) :: vof_th = 0.5
  !
  contains
  !
  !subroutine droplet_tagging(n,dli,vof,u,v,w,istep)
  !subroutine droplet_tagging(n,dli,vof,u,v,w,tmp,istep)
  subroutine droplet_tagging(n,dli,vof,u,v,w,tmp,phi,istep)
  !subroutine droplet_tagging(n,dli,vof,delta,kappa,u,v,w,tmp,phi,tmpge,tmple,sca,mflux,pth,istep,time)
    !
    implicit none
    !
    integer, intent(in), dimension(3)           :: n
    real(8), intent(in), dimension(3)           :: dli
    real(8), intent(in), dimension( 0:, 0:, 0:) :: vof!,delta,kappa
    real(8), intent(in), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in), dimension(-2:,-2:,-2:) :: phi
    !real(8), intent(in), dimension( 0:, 0:, 0:) :: tmpge,tmple
    !real(8), intent(in), dimension( 0:, 0:, 0:) :: sca
    !real(8), intent(in), dimension( 0:, 0:, 0:) :: mflux
    !real(8), intent(in)                         :: pth
    integer, intent(in)                         :: istep
    !real(8), intent(in)                         :: time
    !
    integer(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vofTag
    real(8),    dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vof_f
    real(8)                                           :: rdnInit
    integer                                           :: i,j,k,b, nd, cellcount
    integer(4)                                        :: id, idst
    integer(4), allocatable,    dimension(:,:)        :: procShare, share,idShare, orderField
    integer(4), allocatable,    dimension(:)          :: myid_out, drid
    integer(8), allocatable,    dimension(:)          :: idloc
    real(8),    allocatable,    dimension(:)          :: xd,yd,zd, ud, vd, wd, vold
    character(len=3) ::proidch
    ! 
    call random_number(rdnInit)
    ! idst = floor(rdnInit*100000+100000*myid) ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    idst = 0 ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    id = idst+1
    vof_f = vof
    vofTag = 0
    nd = 0
    ! main loop for tagging droplets 
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if (vof_f(i,j,k).ge.vof_th) then
            vof_f(i,j,k)  = 0.0
            vofTag(i,j,k) = id
            !call recursiveTagging(vof_f,vofTag,id,i,j,k,n)
            call recursiveTagging(n,i,j,k,id,vof_f,vofTag)
            id = id+1
            nd = nd+1
          endif
        enddo
      enddo
    enddo

    ! Update halo
    call updthalo((/n(1),n(2)/),1,vofTag)
    call updthalo((/n(1),n(2)/),2,vofTag)

    ! write(proidch,'(i3.3)') myid
    ! open(90+myid,file='proc'//proidch)
    ! do i=0,n(1)+1
    !   do j=0,n(2)+1
    !     do k=0,n(3)+1
    !       write(90+myid,'(4I5.4)') i,j,k,vofTag(i,j,k)
    !     enddo
    !   enddo
    ! enddo
    ! close(92)

    ! print*,'max vof', maxval(vof_f), myid

    ! allocate droplet properties
    allocate(xd(nd),yd(nd),zd(nd),ud(nd),vd(nd),wd(nd),vold(nd), & 
             idloc(nd),share(nd,6),myid_out(nd),procShare(nd,6), &
             idShare(nd,6), orderField(nd,6),drid(nd))
    xd = 0.
    yd = 0.
    zd = 0.
    ud = 0.
    vd = 0.
    wd = 0.
    vold = 0.
    id = 0
    cellcount = 0
    idShare = 0
    procShare = -1
    share = 0

    ! Loop for properties characterization
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          if (vofTag(i,j,k).gt.0) then
            id       = vofTag(i,j,k)-idst
            drid(id) = id
            idloc(id)= vofTag(i,j,k) 
            xd(id)   = xd(id) + vof(i,j,k)*(i+coord(1)*n(1)-0.5)*dl(1)
            yd(id)   = yd(id) + vof(i,j,k)*(j+coord(2)*n(2)-0.5)*dl(2)
            zd(id)   = zd(id) + vof(i,j,k)*(k              -0.5)*dl(3)
            ud(id)   = ud(id) + u(i,j,k)*vof(i,j,k)
            vd(id)   = vd(id) + v(i,j,k)*vof(i,j,k)
            wd(id)   = wd(id) + w(i,j,k)*vof(i,j,k)
            vold(id) = vold(id) + vof(i,j,k)
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if ((i.eq.1).and. share(id,1).eq.0 .and. vofTag(0,j,k).ne.0) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = vofTag(0,j,k)
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            else if ((i.eq.n(1)).and. share(id,2).eq.0 .and.vofTag(n(1)+1,j,k).ne.0) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = vofTag(n(1)+1,j,k)
            else if ((j.eq.1) .and. share(id,3).eq.0 .and. vofTag(i,0,k).ne.0) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = vofTag(i,0,k)
              ! print*, 'j1', idShare(id,:)
            else if ((j.eq.n(2)) .and. share(id,4).eq.0 .and. vofTag(i,n(2)+1,k).ne.0) then
              share(id,4)        = 1
              procShare(id,4) = back
              idShare(id,4)   = vofTag(i,n(2)+1,k)
              ! print*, 'n2', idShare(id,:)
            else if ((k.eq.1) .and. share(id,5).eq.0) then
              share(id,5)        = 1
              procShare(id,5) = myid
              idShare(id,5)   = vofTag(i,j,0)
            else if ((k.eq.n(3)).and. share(id,6).eq.0) then
              share(id,6)        = 1
              procShare(id,6) = myid
              idShare(id,6)   = vofTag(i,j,n(3)+1)
            endif
          endif
        enddo
      enddo
    enddo 

    do i=1,nd
      do j=1,6
        orderField(i,j) = myid*100+i*10+J
        ! print*,'proc', procShare(i,:), myid
        ! print*, nd, myid
        ! print*,'id  ', idShare(i,:), myid
      enddo

      ! print*,'id  ', idShare(i,:), myid
      ! print*,'proc', procShare(i,:), myid
    enddo

    myid_out = 0
    do id=1,nd
      myid_out(id) = myid
      xd(id)   = xd(id)/vold(id)    
      yd(id)   = yd(id)/vold(id)    
      zd(id)   = zd(id)/vold(id)    
      ud(id)   = ud(id)/vold(id)    
      vd(id)   = vd(id)/vold(id)    
      wd(id)   = wd(id)/vold(id)   
      ! print*, xd(id),yd(id),zd(id),vold(id),id,myid 
      vold(id) = vold(id)*dl(1)*dl(2)*dl(3)
    enddo 

    ! Writing results
    call write_output_r8(istep,'xpos',xd  , nd)
    call write_output_r8(istep,'ypos',yd  , nd)
    call write_output_r8(istep,'zpos',zd  , nd)
    call write_output_r8(istep,'uvel',ud  , nd)
    call write_output_r8(istep,'vvel',vd  , nd)
    call write_output_r8(istep,'wvel',wd  , nd)
    call write_output_r8(istep,'vold',vold, nd)
    call write_output_i1(istep,'drid',drid, nd)
    call write_output_i1(istep,'idmy',myid_out, nd)
    call write_output_i6(istep,'proc',transpose(procShare), nd)
    call write_output_i6(istep,'orde',transpose(orderField), nd)
    call write_output_i6(istep,'idsh',transpose(idShare), nd)
    !
    deallocate(xd,yd,zd,ud,vd,wd,vold, &
               procShare,idShare)
    !
    return
  end subroutine droplet_tagging
  !
  recursive subroutine recursiveTagging(n,i,j,k,id,vof_f,vofTag)
    !
    implicit none
    !
    integer ,   intent(in   ), dimension(3)        :: n
    integer ,   intent(in   )                      :: i,j,k
    integer(4), intent(inout)                      :: id
    real(8)   , intent(inout), dimension(0:,0:,0:) :: vof_f
    integer(8), intent(inout), dimension(0:,0:,0:) :: vofTag
    !
    integer :: ii,jj,kk
    ! 
    do kk=k-1,k+1
      do jj=j-1,j+1
        do ii=i-1,i+1
          !
          if((vof_f(ii,jj,kk).gt.vof_th).and.(ii.ge.1   ).and.(jj.ge.1   ).and.(kk.ge.1   ).and.  &
                                             (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3))) then
            vof_f(ii,jj,kk)  = 0.d0
            vofTag(ii,jj,kk) = id
            ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
            call recursiveTagging(n,ii,jj,kk,id,vof_f,vofTag)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine recursiveTagging
  !
  subroutine write_output_r8(istep,varname,var,nd)
    !
    implicit none
    !
    integer,         intent(in)                :: istep
    character(len=*),intent(in)                :: varname
    real(8),         intent(in), dimension(1:) :: var
    integer,         intent(in)                :: nd
    !
    integer                                    :: offset,fh
    character(7)                               :: fldnum
    integer(kind=MPI_OFFSET_KIND)              :: filesize,disp,nreals,nreals_myid
    character(len=1024)                        :: fname
    !
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    !
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_r8
  !
  subroutine write_output_i6(istep,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                   :: istep, nd
    character(len=4), intent(in)                   :: varname
    integer(4)      , intent(in), dimension(1:,1:) :: var
    !
    integer                                        :: offset, fh, ndtot,nelem
    integer, save                                  :: subarray
    character(7)                                   :: fldnum
    integer(kind=MPI_OFFSET_KIND)                  :: filesize,nreals,nreals_myid, displacement,disp
    character(len=1024)                            :: fname
    ! 
    nelem = 6
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)    
    call MPI_FILE_SET_VIEW(fh,disp*nelem*4, MPI_INTEGER, MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i6
  !
  subroutine write_output_i1(istep,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep, nd
    character(len=4), intent(in)                :: varname
    integer(4),       intent(in), dimension(1:) :: var
    !
    integer                                     :: offset, fh, nelem
    character(7)                                :: fldnum
    integer(kind=MPI_OFFSET_KIND)               :: filesize,disp,nreals,nreals_myid
    character(len=1024)                         :: fname
    ! 
    nelem = 1
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*4*nelem, MPI_INTEGER,MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i1
  !
subroutine updthalo(n,idir,p)
  implicit none
  integer , dimension(2), intent(in) :: n
  integer , intent(in) :: idir
  integer(8), dimension(0:,0:,0:), intent(inout) :: p
#ifdef USE_CUDA
  attributes(managed) :: p
  integer :: istat
#endif
  integer :: i,j,k,n_1,n_2
  !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
  !
  !  this subroutine updates the halos that store info
  !  from the neighboring computational sub-domain
  !
  select case(idir)
  case(1) ! x direction
    !if( .false. ) then
    if(dims(1) .eq.  1) then
#ifdef USE_CUDA
      n_1=n(1)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1 ,j,k) = p(  1,j,k)
          p(0     ,j,k) = p(n_1,j,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(n(1)+1,:,:) = p(   1,:,:)
      p(0     ,:,:) = p(n(1),:,:)
      !$OMP END WORKSHARE
#endif
    else
      n_1=n(1)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          xsl_buf(j,k) = p(  1,j,k)
          xsr_buf(j,k) = p(n_1,j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(xsl_buf(0,0), size( xsl_buf ),MPI_REAL8,left ,0, &
        xrr_buf(0,0), size( xrr_buf ),MPI_REAL8,right,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(xsr_buf(0,0), size( xsr_buf ),MPI_REAL8,right,0, &
        xrl_buf(0,0), size( xrl_buf ),MPI_REAL8,left ,0, &
        comm_cart,status,ierr)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1,j,k) = xrr_buf(j,k)
          p(0    ,j,k) = xrl_buf(j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      !call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
      !                  p(n(1)+1,0,0),1,xhalo,right,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
      !                  p(0   ,0,0),1,xhalo,left ,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  case(2) ! y direction
    !if( .false. ) then
    if( dims(2) .eq.  1 ) then
#ifdef USE_CUDA
      n_2=n(2)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = p(i,  1,k)
          p(i,    0,k) = p(i,n_2,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(:,n(2)+1,:)  = p(:, 1  ,:)
      p(:, 0     ,:) = p(:,n(2),:)
      !$OMP END WORKSHARE
#endif
    else
      n_2=n(2)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          ysl_buf(i,k) = p(i,  1,k)
          ysr_buf(i,k) = p(i,n_2,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(ysl_buf(0,0), size( ysl_buf ),MPI_REAL8,front,0, &
        yrr_buf(0,0), size( yrr_buf ),MPI_REAL8,back ,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(ysr_buf(0,0), size( ysr_buf ),MPI_REAL8,back ,0, &
        yrl_buf(0,0), size( yrl_buf ),MPI_REAL8,front,0, &
        comm_cart,status,ierr)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = yrr_buf(i,k)
          p(i,    0,k) = yrl_buf(i,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
      !                  p(0,n(2)+1,0),1,yhalo,back ,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
      !                  p(0,0   ,0),1,yhalo,front,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  end select

  return
end subroutine updthalo

end module
