module mod_initgrid
  !
  use mod_param, only: pi
  !
  implicit none
  !
  private
  public initgrid
  !
  contains
  !
  subroutine initgrid(inivel,n,gr,lz,dzc,dzf,zc,zf,asr)
    !
    ! initializes the non-uniform grid along z 
    ! (i.e., the non-parallelized direction)
    !
    implicit none
    !
    character(len=3), intent(in )                    :: inivel
    integer         , intent(in )                    :: n
    real(8)         , intent(in )                    :: gr,lz
    real(8)         , intent(out), dimension(-2:n+3) :: dzc,dzf,zc,zf,asr
    !
    real(8) :: z0
    integer :: k
    procedure (), pointer :: gridpoint => null()
    !
    select case(inivel)
    case('zer','log','poi','cou')
      gridpoint => gridpoint_cluster_two_end
    case('hcl','hcp')
      gridpoint => gridpoint_cluster_one_end
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces zf
    !
    do k=1,n
      z0  = (k-0.d0)/(1.d0*n)
      call gridpoint(gr,z0,zf(k))
      zf(k) = zf(k)*lz
    enddo
    zf(0) = 0.d0
    !
    ! step 2) determine grid spacing between faces dzf
    !
    do k=1,n
      dzf(k) = zf(k)-zf(k-1)
    enddo
    dzf(0  ) = dzf(1)
    dzf(n+1) = dzf(n)
    !
    ! step 3) determine grid spacing between centers dzc
    !
    do k=0,n
      dzc(k) = 0.5d0*(dzf(k)+dzf(k+1))
    enddo
    dzc(n+1) = dzc(n-1)
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0)    = -dzc(0)/2.d0
    zf(0)    = 0.d0
    do k=1,n+1
      zc(k)  = zc(k-1) + dzc(k-1)
      zf(k)  = zf(k-1) + dzf(k)
    enddo
    !
    ! step 5) extension to -1,-2 and n+2,n+3 for dzf and dzc
    !
    dzf(-1 ) = dzf(2  )
    dzf(-2 ) = dzf(3  )
    dzf(n+2) = dzf(n-1)
    dzf(n+3) = dzf(n-2)
    !
    dzc(-1 ) = dzc(1  )
    dzc(-2 ) = dzc(2  )
    dzc(n+2) = dzc(n-2)
    dzc(n+3) = dzc(n-3)
    !
    ! step 6) extension to -1,-2 and n+2,n+3 for zf and zc
    !
    ! apart from zf(-2) and zc(-2), we use zc(k)=zc(k-1)+dzc(k-1) and zf(k) = zf(k-1)+dzf(k)
    !
    zf(-2 )  = zf(0  ) - (dzf( 0)+dzf(-1))
    zf(-1 )  = zf(-2 ) + dzf(-1 )
    zf(n+2)  = zf(n+1) + dzf(n+2)
    zf(n+3)  = zf(n+2) + dzf(n+3)
    !
    zc(-2 )  = zc(0  ) - (dzc(-1)+dzc(-2))
    zc(-1 )  = zc(-2 ) + dzc(-2 )
    zc(n+2)  = zc(n+1) + dzc(n+1)
    zc(n+3)  = zc(n+2) + dzc(n+2)
    !
    ! step 7) aspect-ratio
    !
    do k=-2,n+2
      asr(k) = dzf(k)/dzf(k+1)
    enddo
    asr(n+3) = asr(n+2)
    !
    return
  end subroutine initgrid
  !
  ! grid stretching functions 
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi 
  !
  subroutine gridpoint_cluster_two_end(alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    !
    real(8), intent(in ) :: alpha,z0
    real(8), intent(out) :: z
    !
    if(alpha.ne.0.d0) then
      z = 0.5d0*(1.d0+tanh((z0-0.5d0)*alpha)/tanh(alpha/2.d0))
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_two_end
  !
  subroutine gridpoint_cluster_one_end(alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    !
    real(8), intent(in ) :: alpha,z0
    real(8), intent(out) :: z
    !
    if(alpha.ne.0.d0) then
      z = 1.0d0*(1.d0+tanh((z0-1.0d0)*alpha)/tanh(alpha/1.d0))
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_one_end
  !
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    !
    real(8), intent(in ) :: alpha,z0
    real(8), intent(out) :: z
    !
    if(alpha.ne.0.d0) then
      if(    z0.le.0.5d0) then 
        z = 0.5d0*(1.d0-1.d0+tanh(2.d0*alpha*(z0-0.d0))/tanh(alpha))
      elseif(z0.gt.0.5d0) then
        z = 0.5d0*(1.d0+1.d0+tanh(2.d0*alpha*(z0-1.d0))/tanh(alpha))
      endif
    else
      z = z0
    endif
    !
    return
  end subroutine gridpoint_cluster_middle
  !
end module mod_initgrid
