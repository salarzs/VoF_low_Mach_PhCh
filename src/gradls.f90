module mod_gradls
  !
  implicit none
  !
  private 
  public  :: gradls,weno5,houc5,upwd1
  !
  contains
  !
  subroutine gradls(method,n,qmin,dli,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    character(len=*)                                      :: method
    integer, intent(in ), dimension(3)                    :: n
    integer, intent(in )                                  :: qmin
    real(8), intent(in ), dimension(3)                    :: dli
    logical, intent(in )                                  :: is_f
    real(8), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: ux,uy,uz
    real(8), intent(out), dimension(     :,     :,     :) :: dphidt
    !
    select case(method)
    case('houc5')
      call houc5(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    case('weno5')
      call weno5(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    case('upwd1')
      call upwd1(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    end select
    !
    return
  end subroutine gradls
  !
  subroutine weno5(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    real(8), parameter, dimension(3,3) :: c     = 1.d0/6.d0*reshape((/ 2.,-7.,11., &
                                                                      -1., 5., 2., &
                                                                       2., 5.,-1./),shape(c))
    real(8), parameter, dimension(3)   :: sigma = (/.1d0,.6d0,.3d0/) 
    real(8), parameter                 :: eps   = 1.e-6
    !
    integer, intent(in ), dimension(3)                    :: n
    real(8), intent(in ), dimension(3)                    :: dli
    integer, intent(in )                                  :: qmin
    logical, intent(in )                                  :: is_f
    real(8), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: ux,uy,uz
    real(8), intent(out), dimension(     :,     :,     :) :: dphidt
    !
    real(8), dimension(-2:2) :: f
    real(8), dimension(3) :: beta,we,dfdlh
    real(8) :: dphidx,dphidy,dphidz
    integer :: a,i,j,k,p
    real(8) :: uxc,uyc,uzc,sum_we
    integer :: q
    !
    if(is_f) then
      q = 1
    else
      q = 0
    endif
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
#ifdef TWOD
          dphidx = 0.d0
#else
          uxc = 0.5d0*(ux(i-q,j,k)+ux(i,j,k))
          a = nint(sign(1.d0,uxc))
          f(-2) = a*(phi(i-2*a,j,k) - phi(i-3*a,j,k))*dli(1)
          f(-1) = a*(phi(i-1*a,j,k) - phi(i-2*a,j,k))*dli(1)
          f( 0) = a*(phi(i+0*a,j,k) - phi(i-1*a,j,k))*dli(1)
          f( 1) = a*(phi(i+1*a,j,k) - phi(i+0*a,j,k))*dli(1)
          f( 2) = a*(phi(i+2*a,j,k) - phi(i+1*a,j,k))*dli(1)
          beta(1) = (13./12.)*(   f(-2) - 2.*f(-1) +    f( 0))**2 + &
                    ( 1./4. )*(   f(-2) - 4.*f(-1) + 3.*f( 0))**2
          beta(2) = (13./12.)*(   f(-1) - 2.*f( 0) +    f( 1))**2 + &
                    ( 1./4. )*(   f(-1)            -    f( 1))**2
          beta(3) = (13./12.)*(   f( 0) - 2.*f( 1) +    f( 2))**2 + &
                    ( 1./4. )*(3.*f( 0) - 4.*f( 1) +    f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:)    = sigma(:)*(1.d0+(abs(beta(1)-beta(3))/(eps+beta(:))))
          sum_we   = sum(we(:)) ! keep here
          we(:)    = we(:)/sum_we
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidx   = sum(we(:)*dfdlh(:))
#endif
          !
          uyc = 0.5d0*(uy(i,j-q,k)+uy(i,j,k))
          a = nint(sign(1.d0,uyc))
          f(-2) = a*(phi(i,j-2*a,k) - phi(i,j-3*a,k))*dli(2)
          f(-1) = a*(phi(i,j-1*a,k) - phi(i,j-2*a,k))*dli(2)
          f( 0) = a*(phi(i,j+0*a,k) - phi(i,j-1*a,k))*dli(2)
          f( 1) = a*(phi(i,j+1*a,k) - phi(i,j+0*a,k))*dli(2)
          f( 2) = a*(phi(i,j+2*a,k) - phi(i,j+1*a,k))*dli(2)
          beta(1) = (13./12.)*(   f(-2) - 2.*f(-1) +    f( 0))**2 + &
                    ( 1./4. )*(   f(-2) - 4.*f(-1) + 3.*f( 0))**2
          beta(2) = (13./12.)*(   f(-1) - 2.*f( 0) +    f( 1))**2 + &
                    ( 1./4. )*(   f(-1)            -    f( 1))**2
          beta(3) = (13./12.)*(   f( 0) - 2.*f( 1) +    f( 2))**2 + &
                    ( 1./4. )*(3.*f( 0) - 4.*f( 1) +    f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:)    = sigma(:)*(1.d0+(abs(beta(1)-beta(3))/(eps+beta(:))))
          sum_we   = sum(we(:)) ! keep here
          we(:)    = we(:)/sum_we
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidy   = sum(we(:)*dfdlh(:))
          !
          uzc = 0.5d0*(uz(i,j,k-q)+uz(i,j,k))
          a = nint(sign(1.d0,uzc))
          f(-2) = a*(phi(i,j,k-2*a) - phi(i,j,k-3*a))*dli(3)
          f(-1) = a*(phi(i,j,k-1*a) - phi(i,j,k-2*a))*dli(3)
          f( 0) = a*(phi(i,j,k+0*a) - phi(i,j,k-1*a))*dli(3)
          f( 1) = a*(phi(i,j,k+1*a) - phi(i,j,k+0*a))*dli(3)
          f( 2) = a*(phi(i,j,k+2*a) - phi(i,j,k+1*a))*dli(3)
          beta(1) = (13./12.)*(   f(-2) - 2.*f(-1) +    f( 0))**2 + &
                    ( 1./4. )*(   f(-2) - 4.*f(-1) + 3.*f( 0))**2
          beta(2) = (13./12.)*(   f(-1) - 2.*f( 0) +    f( 1))**2 + &
                    ( 1./4. )*(   f(-1)            -    f( 1))**2
          beta(3) = (13./12.)*(   f( 0) - 2.*f( 1) +    f( 2))**2 + &
                    ( 1./4. )*(3.*f( 0) - 4.*f( 1) +    f( 2))**2
          !we(:) = sigma(:)/(beta(:)+eps)**2
          we(:)    = sigma(:)*(1.d0+(abs(beta(1)-beta(3))/(eps+beta(:))))
          sum_we   = sum(we(:)) ! keep here
          we(:)    = we(:)/sum_we
          dfdlh(1) = sum(c(:,1)*f(-2:0))
          dfdlh(2) = sum(c(:,2)*f(-1:1))
          dfdlh(3) = sum(c(:,3)*f( 0:2))
          dphidz   = sum(we(:)*dfdlh(:))
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine weno5
  !
  subroutine houc5(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    real(8), parameter, dimension(-3:2) :: c = 1.d0/60.*(/-2.,15.,-60.,20.,30.,-3./)
    !
    integer, intent(in ), dimension(3)                    :: n
    real(8), intent(in ), dimension(3)                    :: dli
    integer, intent(in )                                  :: qmin
    logical, intent(in )                                  :: is_f
    real(8), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: ux,uy,uz
    real(8), intent(out), dimension(     :,     :,     :) :: dphidt
    !
    real(8) :: dphidx,dphidy,dphidz
    integer :: a,i,j,k,q,qq,p
    real(8) :: uxc,uyc,uzc
    !
    if(is_f) then
      qq = 1
    else
      qq = 0
    endif
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          uxc = 0.5d0*(ux(i-qq,j,k)+ux(i,j,k))
          a = nint(sign(1.d0,uxc))
          dphidx = 0.d0
          do q=-3,2
            dphidx = dphidx + a*(phi(i+q*a,j,k)*c(q))*dli(1)
          enddo
          !
          uyc = 0.5d0*(uy(i,j-qq,k)+uy(i,j,k))
          a = nint(sign(1.d0,uyc))
          dphidy = 0.d0
          do q=-3,2
            dphidy = dphidy + a*(phi(i,j+q*a,k)*c(q))*dli(2)
          enddo
          !
          uzc = 0.5d0*(uz(i,j,k)+uz(i,j,k-qq))
          a = nint(sign(1.d0,uzc))
          dphidz = 0.d0
          do q=-3,2
            dphidz = dphidz + a*(phi(i,j,k+q*a)*c(q))*dli(3)
          enddo
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine houc5
  !
  subroutine upwd1(n,dli,qmin,is_f,phi,ux,uy,uz,dphidt)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                    :: n
    real(8), intent(in ), dimension(3)                    :: dli
    integer, intent(in )                                  :: qmin
    logical, intent(in )                                  :: is_f
    real(8), intent(in ), dimension(   -2:,   -2:,   -2:) :: phi
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: ux,uy,uz
    real(8), intent(out), dimension(     :,     :,     :) :: dphidt
    !
    real(8) :: dphidx,dphidy,dphidz
    integer :: a,i,j,k,q,p
    real(8) :: uxc,uyc,uzc
    !
    if(is_f) then
      q = 1
    else
      q = 0
    endif
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          uxc = 0.5d0*(ux(i-q,j,k)+ux(i,j,k))
          a = nint(sign(1.d0,uxc))
          dphidx = a*( phi(i,j,k) - phi(i-a,j,k) )*dli(1)
          !
          uyc = 0.5d0*(uy(i,j-q,k)+uy(i,j,k))
          a = nint(sign(1.d0,uyc))
          dphidy = a*( phi(i,j,k) - phi(i,j-a,k) )*dli(2)
          !
          uzc = 0.5d0*(uz(i,j,k-q)+uz(i,j,k))
          a = nint(sign(1.d0,uzc))
          dphidz = a*( phi(i,j,k) - phi(i,j,k-a) )*dli(3)
          !
          dphidt(i,j,k) = - (uxc*dphidx + uyc*dphidy + uzc*dphidz)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine upwd1
  !
end module mod_gradls
