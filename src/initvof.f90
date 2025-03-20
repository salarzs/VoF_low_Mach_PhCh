module mod_initvof
  !
  use mod_param     , only: inivof,lx,ly,lz,cbcvof,xcc,ycc,zcc,rd,n_dp!,small_tg
  use decomp_2d
  use mod_common_mpi, only: myid,ierr,coord
  !
  implicit none
  !
  real(8), parameter, public :: b_th = 2.d0
  !
  private
  public initvof_ma!,initvof_other
  !
  contains
  !
  !
  !subroutine initvof(n,dli,vof,vof_index)
  !  !
  !  ! computes initial conditions for the volume of fluid field
  !  !
  !  implicit none
  !  !
  !  integer, intent(in ), dimension(3)        :: n
  !  real(8), intent(in ), dimension(3)        :: dli
  !  real(8), intent(out), dimension(0:,0:,0:) :: vof
  !  real(8), intent(out), dimension(0:,0:,0:) :: vof_index
  !  !
  !  integer :: i,j,k,q,ii,jj,kk,p
  !  real(8) :: x,y,z,xl,yl,zl,xx,yy,zz
  !  real(8) :: sdist,sdistmin
  !  real(8) :: eps,epsbox
  !  real(8), dimension(3) :: dl,dlbox
  !  integer, dimension(3) :: iperiod
  !  real(8) :: grid_vol_ratio
  !  !
  !  integer, parameter :: nbox = 50
  !  !
  !  dl(:) = dli(:)**(-1)
  !  dlbox(:) = dl(:)/(1.d0*nbox)
  !  !
  !  eps    = 0.5d0*sqrt(dl(1)**2    + dl(2)**2    + dl(3)**2   )
  !  epsbox = 0.5d0*sqrt(dlbox(1)**2 + dlbox(2)**2 + dlbox(3)**2)
  !  !    
  !  grid_vol_ratio = product(dlbox(:))/product(dl(:))
  !  iperiod(:) = 1
  !  !
  !  vof(:,:,:) = 0.d0
  !  !
  !  do k=1,n(3)
  !    z = (k-0.5d0)*dl(3)
  !    do j=1,n(2)
  !      y = (j+coords(2)*n(2)-0.5d0)*dl(2)
  !      do i=1,n(1)
  !        x = (i+coords(1)*n(1)-0.5d0)*dl(1)
  !        !
  !        sdistmin = max(lx,ly,lz)*2.d0
  !        !
  !        do p=1,n_dp ! loop over all the dispersed phases
  !          !
  !          do kk = -1,1
  !            do jj = -1,1
  !              do ii = -1,1
  !                !
  !                sdist = sqrt( (x+ii*iperiod(1)*lx-xc(p))**2 + &
  !                              (y+jj*iperiod(2)*ly-yc(p))**2 + &
  !                              (z+kk*iperiod(3)*lz-zc(p))**2 ) - rd(p)
  !                if(abs(sdist).lt.sdistmin) sdistmin = sdist
  !                !
  !              enddo
  !            enddo
  !          enddo
  !          !
  !          ! update the tag for the inside bulk region
  !          !
  !          if(sdist.lt.0.d0) then
  !            vof_index(i,j,k) = 1.d0*p
  !          endif
  !          !
  !          sdist = sdistmin
  !          !
  !          if(     sdist.lt.-eps ) then
  !            vof(i,j,k) = 1.d0
  !          elseif( sdist.gt. eps ) then
  !            vof(i,j,k) = 0.d0
  !          else
  !            zl = z-0.5d0*dl(3)
  !            yl = y-0.5d0*dl(2)
  !            xl = x-0.5d0*dl(1)
  !            !
  !            do kk=1,nbox
  !              zz = zl + kk*dlbox(3)
  !              do jj=1,nbox
  !                yy = yl + jj*dlbox(2)
  !                do ii=1,nbox
  !                  xx = xl + ii*dlbox(1)
  !                  !
  !                  sdist = sqrt((xx-xc(p))**2 + (yy-yc(p))**2 + (zz-zc(p))**2) - rd(p)
  !                  !
  !                  if(sdist.lt.-epsbox) then
  !                    vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
  !                    vof_index(i,j,k) = 1.d0*p
  !                  endif
  !                  !
  !                enddo
  !              enddo
  !            enddo
  !            !
  !          endif
  !          !
  !        enddo
  !        !
  !      enddo
  !    enddo
  !  enddo
  !  !
  !  return
  !end subroutine initvof
  !
  subroutine initvof_ma(n,dli,n_dp,xc,yc,zc,rc,vof)
    !
    ! computes initial conditions for the volume of fluid field
    !
    implicit none
    !
    integer, intent(in ), dimension(3)        :: n
    real(8), intent(in ), dimension(3)        :: dli
    integer, intent(in )                      :: n_dp
    real(8), intent(in ), dimension(n_dp)     :: xc,yc,zc,rc
    real(8), intent(out), dimension(0:,0:,0:) :: vof
    !
    real(8), dimension(3) :: dl,dlbox
    integer, dimension(3) :: iperiod
    integer, dimension(2) :: nx_b,ny_b,nz_b
    integer :: i,j,k,q,ii,jj,kk,p
    real(8) :: x,y,z,xl,yl,zl,xx,yy,zz
    real(8) :: sdist,sdistmin
    real(8) :: eps,epsbox
    real(8) :: grid_vol_ratio
    !
    integer, parameter :: nbox = 50
    !
    dl(:) = dli(:)**(-1)
    dlbox(:) = dl(:)/(1.d0*nbox)
    !
    eps    = 0.5d0*sqrt(dl(1)**2    + dl(2)**2    + dl(3)**2   )
    epsbox = 0.5d0*sqrt(dlbox(1)**2 + dlbox(2)**2 + dlbox(3)**2)
    !    
    grid_vol_ratio = product(dlbox(:))/product(dl(:))
    iperiod(:) = 1
    !
    vof(:,:,:) = 0.d0
    !
    do k=1,n(3)
      z = (k-0.5d0)*dl(3)
      do j=1,n(2)
        y = (j+coord(2)*n(2)-0.5d0)*dl(2)
        do i=1,n(1)
          x = (i+coord(1)*n(1)-0.5d0)*dl(1)
          !
          sdistmin = max(lx,ly,lz)*2.d0
          !
          do p=1,n_dp ! loop over all the dispersed phases
            !
            !print*, n_dp,p
            call bub_lim(n(1),coord(1),xc(p),rc(p),dl(1),nx_b)
            call bub_lim(n(2),coord(2),yc(p),rc(p),dl(2),ny_b)
            call bub_lim(n(3),0       ,zc(p),rc(p),dl(3),nz_b)
            !
            do kk = -1,1
              do jj = -1,1
                do ii = -1,1
                  !
                  sdist = sqrt( (x+ii*iperiod(1)*lx-xc(p))**2 + &
                                (y+jj*iperiod(2)*ly-yc(p))**2 + &
                                (z+kk*iperiod(3)*lz-zc(p))**2 ) - rc(p)
                  if(abs(sdist).lt.sdistmin) sdistmin = sdist
                  !
                enddo
              enddo
            enddo
            !
            ! update the tag for the inside bulk region
            !
            !if(sdist.lt.0.d0) then
            !  vof_index(i,j,k) = 1.d0*p
            !endif
            !
            sdist = sdistmin
            !
            if(     sdist.lt.-eps ) then
              vof(i,j,k) = 1.d0
            elseif( sdist.gt. eps ) then
              vof(i,j,k) = 0.d0
            else
              zl = z-0.5d0*dl(3)
              yl = y-0.5d0*dl(2)
              xl = x-0.5d0*dl(1)
              !
              do kk=1,nbox
                zz = zl + kk*dlbox(3)
                do jj=1,nbox
                  yy = yl + jj*dlbox(2)
                  do ii=1,nbox
                    xx = xl + ii*dlbox(1)
                    !
                    sdist = sqrt((xx-xc(p))**2 + (yy-yc(p))**2 + (zz-zc(p))**2) - rc(p)
                    !
                    if(sdist.lt.-epsbox) then
                      vof(i,j,k) = vof(i,j,k) + grid_vol_ratio
                      !vof_index(i,j,k) = 1.d0*p
                    endif
                    !
                  enddo
                enddo
              enddo
              !
            endif
            !
          enddo
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine initvof_ma
  !
  subroutine bub_lim(n,coord,xc,r,dx,nx_b)
    !
    implicit none
    !
    integer, intent(in )               :: n,coord
    real(8), intent(in )               :: xc,r,dx
    integer, intent(out), dimension(2) :: nx_b
    !
    integer :: lb,ub
    !
    lb = floor((xc-r)/dx) - 1  
    if(    lb.le.(n*coord+1)) then
      nx_b(1) = 1
    elseif(lb.ge.(n*coord+n)) then
      nx_b(1) = n
    else
      nx_b(1) = lb-n*coord 
    endif
    !
    ub = ceiling((xc+r)/dx) + 1
    if (ub.le.(n*coord+1)) then
      nx_b(2) = 1
    elseif (ub.ge.(n*coord+n)) then
      nx_b(2) = n
    else
      nx_b(2) = ub-n*coord 
    endif
    !
    return
  end subroutine bub_lim
  !
end module mod_initvof
