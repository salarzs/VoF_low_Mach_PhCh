module mod_vof
  !
  ! --> written by Marco Edoardo Rosti based on:
  !     Ii, Satoshi, et al. "An interface capturing method with a continuous function: 
  !     The THINC method with multi-dimensional reconstruction." 
  !     Journal of Computational Physics 231.5 (2012): 2328-2358.
  ! --> modified by Pedro Costa for compatibility
  !     note: implemented for constant grid spacing
  ! --> modified by Nicolo Scapin to improve efficiency
  !
  use mod_initvof, only: b_th
  !
  implicit none
  !
  integer , parameter :: rp    = 8
  !real(rp), parameter :: b_th  = 1.0d0, qu_th  = 1.0d0, &
  real(rp), parameter :: qu_th  = 1.0d0, &
                         limit = 1.d-8, limitf = 1.d-5
  real(rp), parameter :: b_thi = 1.d0/b_th
  real(rp), parameter :: small = 1e-16
  !
  private
  public  :: advvof,update_vof,update_property
  !
  contains
  !
  subroutine advvof(n,dli,dt,lvg,ldz,dzc,dzf,ug,vg,wg,vof,nor,cur,kappa,d_thinc)
    !
    use mod_param, only: cbcvof,bcvof
    use mod_bound, only: boundp
    !
    ! VoF advection
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)              :: n
    real(rp), intent(in   ), dimension(3)              :: dli
    real(rp), intent(in   )                            :: dt
    integer , intent(in   )                            :: lvg,ldz
    real(rp), intent(in   ), dimension(ldz:)           :: dzc,dzf
    real(rp), intent(in   ), dimension(lvg:,lvg:,lvg:) :: ug,vg,wg ! interface velocity
    real(rp), intent(inout), dimension(0:,0:,0:)       :: vof
    real(rp), intent(inout), dimension(0:,0:,0:,1:)    :: nor,cur
    real(rp), intent(inout), dimension(0:,0:,0:)       :: kappa    ! curvature
    real(rp), intent(inout), dimension(0:,0:,0:)       :: d_thinc
    !
    real(rp), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: dvof1,dvof2,dvof3!cc,dvof3
    real(rp), dimension(0:n(1)  ,0:n(2)  ,0:n(3)  ) :: flux
    real(rp), dimension(3) :: dl
    integer :: i,j,k,im,jm,km
    !
    dl(:) = dli(:)**(-1)
    !
    !do k=1,n(3)
    !  do j=1,n(2)
    !    do i=1,n(1)
    !      if(  vof(i,j,k).gt.0.5d0) then
    !        cc(i,j,k) = 1.d0
    !      else
    !        cc(i,j,k) = 0.d0
    !      endif
    !    enddo
    !  enddo
    !enddo
    !call boundp(cbcvof,n,bcvof,dl,dzc,dzf,cc)
    !
    ! flux in x
    !
#ifndef TWOD
    call cmpt_vof_flux(n,dli(1),dt,lvg,vof,nor,cur,d_thinc,1,ug,flux) 
#endif
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          im = i-1
          !
#ifdef TWOD
          dvof1(i,j,k) = vof(i,j,k)
#else
          dvof1(i,j,k) = (vof(i,j,k)-(flux(i,j,k)-flux(im,j,k))*dli(1))/(1.d0-dt*dli(1)*(ug(i,j,k)-ug(im,j,k)))
#endif
          !
        enddo
      enddo
    enddo
    !
    ! update vof
    !
    call filter_vof(n,dvof1)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,dvof1)
    call update_vof(n,dli,ldz,dzc,dzf,dvof1,nor,cur,kappa,d_thinc)
    !
    ! flux in y
    !
    call cmpt_vof_flux(n,dli(2),dt,lvg,dvof1,nor,cur,d_thinc,2,vg,flux) 
    !
    do k=1,n(3)
      do j=1,n(2)
        jm = j-1
        do i=1,n(1)
          !
          dvof2(i,j,k) = (dvof1(i,j,k)-(flux(i,j,k)-flux(i,jm,k))*dli(2))/(1.d0-dt*dli(2)*(vg(i,j,k)-vg(i,jm,k)))
          !
        enddo
      enddo
    enddo
    !
    ! update vof
    !
    call filter_vof(n,dvof2)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,dvof2)
    call update_vof(n,dli,ldz,dzc,dzf,dvof2,nor,cur,kappa,d_thinc)
    !
    ! flux in z
    !
    call cmpt_vof_flux(n,dli(3),dt,lvg,dvof2,nor,cur,d_thinc,3,wg,flux) 
    !
    do k=1,n(3)
      km = k-1
      do j=1,n(2)
        do i=1,n(1)
          !
          dvof3(i,j,k) = (dvof2(i,j,k)-(flux(i,j,k)-flux(i,j,km))*dli(3))/(1.d0-dt*dli(3)*(wg(i,j,k)-wg(i,j,km)))
          !
        enddo
      enddo
    enddo
    !
    ! update vof
    !
    call filter_vof(n,dvof3)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,dvof3)
    !call update_vof(n,dli,ldz,dzc,dzf,dvof3,nor,cur,kappa,d_thinc) ! it can be skipped
    !
    ! divergence correction step
    !
    do k=1,n(3)
      km = k-1
      do j=1,n(2)
        jm = j-1
        do i=1,n(1)
          im = i-1
          !
          vof(i,j,k) = dvof3(i,j,k) - dt*( dvof1(i,j,k)*dli(1)*(ug(i,j,k)-ug(im,j,k)) + &
                                           dvof2(i,j,k)*dli(2)*(vg(i,j,k)-vg(i,jm,k)) + &
                                           dvof3(i,j,k)*dli(3)*(wg(i,j,k)-wg(i,j,km)) )
          !
#ifdef VAP_MASS
          vof(i,j,k) = vof(i,j,k)/( &
                                    1.d0-dt*( &
                                              dli(1)*(ug(i,j,k)-ug(im,j,k)) + &
                                              dli(2)*(vg(i,j,k)-vg(i,jm,k)) + &
                                              dli(3)*(wg(i,j,k)-wg(i,j,km))   &
                                            ) &
                                  )
#endif
          !
        enddo
      enddo
    enddo
    !
    ! update vof
    !
    call filter_vof(n,vof)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,vof)
    call update_vof(n,dli,ldz,dzc,dzf,vof,nor,cur,kappa,d_thinc)
    !
    return
  end subroutine advvof
  !
  subroutine update_vof(n,dli,ldz,dzc,dzf,vof,nor,cur,kappa,d_thinc)
    !
    use mod_param, only: cbcvof,bcvof
    use mod_bound, only: boundp
    !
    implicit none
    !
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dli
    integer , intent(in )                         :: ldz
    real(rp), intent(in ), dimension(ldz:)        :: dzc,dzf
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    real(rp), intent(out), dimension(0:,0:,0:   ) :: d_thinc
    !
    real(rp), dimension(3) :: dl
    integer :: p,i,j,k
    !
    dl(:) = dli(:)**(-1)
    !
    call cmpt_nor_curv(n,dli,vof,nor,cur,kappa)
    do p=1,3
      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,nor(:,:,:,p))
    enddo
    do p=1,6
      call boundp(cbcvof,n,bcvof,dl,dzc,dzf,cur(:,:,:,p))
    enddo
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,kappa)
    !
    call cmpt_d_thinc(n,nor,cur,vof,d_thinc)
    call boundp(cbcvof,n,bcvof,dl,dzc,dzf,d_thinc)
    !
    return
  end subroutine update_vof
  !
  subroutine update_property(n,prop12,vof,prop)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(2)        :: prop12
    real(rp), intent(in ), dimension(0:,0:,0:) :: vof
    real(rp), intent(out), dimension(0:,0:,0:) :: prop
    !
    integer :: i,j,k
    ! 
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          prop(i,j,k) = vof(i,j,k)*prop12(1)+(1.d0-vof(i,j,k))*prop12(2)
        enddo
      enddo
    enddo
    !
    return
  end subroutine update_property
  !
  subroutine filter_vof(n,vof)
    !
    implicit none
    !
    integer , intent(in   ), dimension(3)        :: n
    real(rp), intent(inout), dimension(0:,0:,0:) :: vof
    !
    integer :: i,j,k
    !
    ! the next "brute force" seems needed only for 
    ! some challenging cases (e.g., mixing layer)
    !   --> please, comment in advvof if you do not need
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          vof(i,j,k) = min(max(0.d0,vof(i,j,k)),1.d0) 
        enddo
      enddo
    enddo
    !
    return
  end subroutine filter_vof
  !
  subroutine cmpt_nor_curv(n,dli,vof,nor,cur,kappa)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(3)           :: dli
    real(rp), intent(in ), dimension(0:,0:,0:   ) :: vof
    real(rp), intent(out), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(out), dimension(0:,0:,0:   ) :: kappa
    !
    real(rp), dimension(8) :: nx,ny,nz,mx,my,mz
    real(rp), dimension(3) :: dl
    real(rp) :: norm
    integer  :: i,j,k,p
    !
    dl(:) = dli(:)**(-1)
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
          do p=1,8
            norm  = sqrt(mx(p)**2+my(p)**2+mz(p)**2+small)
            nx(p) = mx(p)/norm
            ny(p) = my(p)/norm
            nz(p) = mz(p)/norm
          enddo
          !
          ! compute the normal vector
          !
          nor(i,j,k,1) = 0.125d0*sum(mx(:))
          nor(i,j,k,2) = 0.125d0*sum(my(:))
          nor(i,j,k,3) = 0.125d0*sum(mz(:))
          norm         = sqrt(nor(i,j,k,1)**2+nor(i,j,k,2)**2+nor(i,j,k,3)**2+small)
          nor(i,j,k,1) = nor(i,j,k,1)/norm
          nor(i,j,k,2) = nor(i,j,k,2)/norm
          nor(i,j,k,3) = nor(i,j,k,3)/norm
          ! 
          ! compute the curvature tensor
          !
          cur(i,j,k,1) = ((nx(1)+nx(2)+nx(3)+nx(4))-(nx(5)+nx(6)+nx(7)+nx(8)))*dl(1)*0.25d0
          cur(i,j,k,2) = ((ny(1)+ny(3)+ny(5)+ny(7))-(ny(2)+ny(4)+ny(6)+ny(8)))*dl(2)*0.25d0
          cur(i,j,k,3) = ((nz(1)+nz(2)+nz(5)+nz(6))-(nz(3)+nz(4)+nz(7)+nz(8)))*dl(3)*0.25d0
          cur(i,j,k,4) = ((ny(1)+ny(2)+ny(3)+ny(4))-(ny(5)+ny(6)+ny(7)+ny(8)))*dl(2)*0.25d0*0.5d0+&
                         ((nx(1)+nx(3)+nx(5)+nx(7))-(nx(2)+nx(4)+nx(6)+nx(8)))*dl(1)*0.25d0*0.5d0
          cur(i,j,k,5) = ((nz(1)+nz(2)+nz(3)+nz(4))-(nz(5)+nz(6)+nz(7)+nz(8)))*dl(3)*0.25d0*0.5d0+&
                         ((nx(1)+nx(2)+nx(5)+nx(6))-(nx(3)+nx(4)+nx(7)+nx(8)))*dl(1)*0.25d0*0.5d0
          cur(i,j,k,6) = ((nz(1)+nz(3)+nz(5)+nz(7))-(nz(2)+nz(4)+nz(6)+nz(8)))*dl(3)*0.25d0*0.5d0+&
                         ((ny(1)+ny(2)+ny(5)+ny(6))-(ny(3)+ny(4)+ny(7)+ny(8)))*dl(2)*0.25d0*0.5d0
          !
          kappa(i,j,k) = -(cur(i,j,k,1)*dli(1)**2+cur(i,j,k,2)*dli(2)**2+cur(i,j,k,3)*dli(3)**2) ! curvature
          cur(i,j,k,:) = cur(i,j,k,:)*qu_th
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_nor_curv
  !
  subroutine cmpt_d_thinc(n,nor,cur,vof,d_thinc)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)           :: n
    real(rp), intent(in ), dimension(0:,0:,0:,1:) :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:)    :: vof
    real(rp), intent(out), dimension(0:,0:,0:)    :: d_thinc
    !
    real(rp), dimension(3,3) :: av,bv,ev
    real(rp), dimension(3) :: nor_v,cv,dv
    real(rp), dimension(6) :: cur_v
#ifdef TWOD
    real(rp) :: fm,fp,a2,b2,c2
#else
    real(rp) :: fpp,fmm,fpm,fmp,a4,b4,c4,d4,e4
#endif
    real(rp) :: aa,qq,surf
    real(rp) :: dtemp
    integer  :: i,j,k,p
    !
    real(rp), parameter :: rmm = 0.5d0*(1.d0-1.d0/sqrt(3.d0)), &
                           rpp = 0.5d0*(1.d0+1.d0/sqrt(3.d0))
    !
    real(rp), parameter, dimension(3,3) :: xmm = reshape((/0.d0,rmm,rmm, &
                                                           rmm,0.d0,rmm, &
                                                           rmm,rmm,0.d0/),shape(xmm))
    real(rp), parameter, dimension(3,3) :: xmp = reshape((/0.d0,rmm,rpp, &
                                                           rmm,0.d0,rpp, &
                                                           rmm,rpp,0.d0/),shape(xmp))
    real(rp), parameter, dimension(3,3) :: xpm = reshape((/0.d0,rpp,rmm, &
                                                           rpp,0.d0,rmm, &
                                                           rpp,rmm,0.d0/),shape(xpm))
    real(rp), parameter, dimension(3,3) :: xpp = reshape((/0.d0,rpp,rpp, &
                                                           rpp,0.d0,rpp, &
                                                           rpp,rpp,0.d0/),shape(xpp))
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if((vof(i,j,k).le.limit).or.(vof(i,j,k).ge.1.d0-limit)) then
            !
            d_thinc(i,j,k) = -1000.d0
            !
          else
            !
            nor_v(:) = nor(i,j,k,:)
            cur_v(:) = cur(i,j,k,:)
            cv(:) = 1.d0
            cv(maxloc(abs(nor_v(:)))) = 0.d0
            !
            ! calculation of the coefficients for the surface function
            !
            ev(:,1) = (/0.d0,cur_v(2),cur_v(3)/)*0.5d0
            ev(:,2) = (/cur_v(1),0.d0,cur_v(3)/)*0.5d0
            ev(:,3) = (/cur_v(1),cur_v(2),0.d0/)*0.5d0
            av(:,1) = (/0.d0,0.d0,cur_v(6)/)
            av(:,2) = (/0.d0,cur_v(5),0.d0/)
            av(:,3) = (/cur_v(4),0.d0,0.d0/)
            bv(:,1) = (/0.d0,nor_v(2)-0.5d0*(cur_v(2)+cur_v(6)),nor_v(3)-0.5d0*(cur_v(3)+cur_v(6))/)
            bv(:,2) = (/nor_v(1)-0.5d0*(cur_v(1)+cur_v(5)),0.d0,nor_v(3)-0.5d0*(cur_v(3)+cur_v(5))/)
            bv(:,3) = (/nor_v(1)-0.5d0*(cur_v(1)+cur_v(4)),nor_v(2)-0.5d0*(cur_v(2)+cur_v(4)),0.d0/)
            dv(1)   = (nor_v(1)-0.5d0*cv(1)*(cur_v(1)+cv(2)*cur_v(4)+cv(3)*cur_v(5)))
            dv(2)   = (nor_v(2)-0.5d0*cv(2)*(cur_v(2)+cv(1)*cur_v(4)+cv(3)*cur_v(6)))
            dv(3)   = (nor_v(3)-0.5d0*cv(3)*(cur_v(3)+cv(1)*cur_v(5)+cv(2)*cur_v(6)))
            !
            ! build the polynomial equation and find its roots
            !
            aa = 0.d0
            qq = 0.d0
            !
#ifdef TWOD
            !
            ! --> 2D: quadratic equation
            !
            fm = 0.d0
            fp = 0.d0
            do p=2,3
              !
              aa = aa + (1.d0-cv(p))*exp(2.d0*b_th*dv(p))
              qq = qq + (1.d0-cv(p))*exp(2.d0*b_th*dv(p)*(2.d0*vof(i,j,k)-1.d0))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fm   = fm + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fp   = fp + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              !
            enddo
            !
            a2 = aa*fm*fp*(aa-qq)
            b2 = aa*(fm+fp)*(1.d0-qq)
            c2 = 1.d0-aa*qq
            call solve_quad(c2,b2,a2,dtemp)
            !
#else
            !
            ! --> 3D: quartic equation
            !
            fmm = 0.d0
            fmp = 0.d0
            fpm = 0.d0
            fpp = 0.d0
            do p=1,3
              !
              aa = aa + (1.d0-cv(p))*exp(2.d0*b_th*dv(p))
              qq = qq + (1.d0-cv(p))*exp(4.d0*b_th*dv(p)*(2.d0*vof(i,j,k)-1.d0))
              !
              surf = sum(ev(:,p)*xmm(:,p)**2) + sum(bv(:,p)*xmm(:,p)) + &
                     av(1,p)*xmm(1,p)*xmm(2,p)+av(2,p)*xmm(1,p)*xmm(3,p)+av(3,p)*xmm(2,p)*xmm(3,p)
              fmm  = fmm + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              surf = sum(ev(:,p)*xmp(:,p)**2) + sum(bv(:,p)*xmp(:,p)) + &
                     av(1,p)*xmp(1,p)*xmp(2,p)+av(2,p)*xmp(1,p)*xmp(3,p)+av(3,p)*xmp(2,p)*xmp(3,p)
              fmp  = fmp + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              surf = sum(ev(:,p)*xpm(:,p)**2) + sum(bv(:,p)*xpm(:,p)) + &
                     av(1,p)*xpm(1,p)*xpm(2,p)+av(2,p)*xpm(1,p)*xpm(3,p)+av(3,p)*xpm(2,p)*xpm(3,p)
              fpm  = fpm + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              surf = sum(ev(:,p)*xpp(:,p)**2) + sum(bv(:,p)*xpp(:,p)) + &
                     av(1,p)*xpp(1,p)*xpp(2,p)+av(2,p)*xpp(1,p)*xpp(3,p)+av(3,p)*xpp(2,p)*xpp(3,p)
              fpp  = fpp + (1.d0-cv(p))*exp(2.d0*b_th*surf)
              !
            enddo
            !
            a4 = aa**2*(fmm*fmp*fpm*fpp)*(aa**2-qq)
            b4 = aa**2*(fmm*fpm*fpp+fmp*fpm*fpp+fmm*fmp*fpm+fmm*fmp*fpp)*(aa-qq)
            c4 = aa**2*(fpm*fpp+fmm*fpm+fmp*fpm+fmm*fpp+fmm*fmp+fmp*fpp)*(1.d0-qq)
            d4 = aa*(fpm+fpp+fmm+fmp)*(1.d0-aa*qq)
            e4 = 1.d0-aa**2*qq
            !call solve_quar(e4,d4,c4,b4,a4,dtemp)
            call solve_quar_paper(e4,d4,c4,b4,a4,dtemp)
            !
#endif
            !
            d_thinc(i,j,k) = 0.5d0*b_thi*log(dtemp)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_d_thinc
  !
  subroutine cmpt_vof_flux(n,dli,dt,lvg,vof,nor,cur,d_thinc,dir,vel,flux)
    !
    implicit none
    !
    integer , intent(in ), dimension(3)              :: n
    real(rp), intent(in )                            :: dli
    real(rp), intent(in )                            :: dt
    integer , intent(in )                            :: lvg
    real(rp), intent(in ), dimension(0:,0:,0:   )    :: vof
    real(rp), intent(in ), dimension(0:,0:,0:,1:)    :: nor,cur
    real(rp), intent(in ), dimension(0:,0:,0:   )    :: d_thinc
    integer , intent(in )                            :: dir
    real(rp), intent(in ), dimension(lvg:,lvg:,lvg:) :: vel
    real(rp), intent(out), dimension(  0:,  0:,  0:) :: flux
    !
    real(rp), dimension(3)   :: nor_v,cv
    real(rp), dimension(6)   :: cur_v
    real(rp), dimension(6,3) :: f2_l
    real(rp), dimension(8,3) :: xv
    real(rp), dimension(6)   :: crd
    integer , dimension(3)   :: f1_l
    real(rp) :: xa,xb,ya,yb,za,zb,cf,vel_sign,dl
    real(rp) :: a,b,rm1,rp1,rm2,rp2
    real(rp) :: cxxa,cyya,czza,cxya,cyza,cxza,a100,a010,a001,sum_int_func
    real(rp) :: surf_a,surf_b
    integer  :: ind2,ii,jj,kk,i,j,k,p,q
    !
    real(rp), parameter :: gm = -1.d0/sqrt(3.d0), &
                           gp = +1.d0/sqrt(3.d0)
    !
    dl = 1.d0/dli
    !
    select case(dir)
     case(1)
      crd(:) = (/0.d0,0.d0,0.d0,1.d0,0.d0,1.d0/) ! along x
     case(2) 
      crd(:) = (/0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/) ! along y
     case(3)
      crd(:) = (/0.d0,1.d0,0.d0,1.d0,0.d0,0.d0/) ! along z
    end select
    !
    ! initialize two auxiliary arrays
    !
    f1_l(:)   = 0
    f2_l(:,:) = 0.d0
    !
    do k=0,n(3)
      do j=0,n(2)
        do i=0,n(1)
          !
          ! decide the upwind path
          !
          vel_sign = sign(1.d0,vel(i,j,k))
          !
          ! compute the indexes, ii,jj,kk
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          f1_l(dir) = nint(0.5d0*(-vel_sign+1.d0)) ! i.e., f1_l(dir) = 0 (vel>=0) or 1 (vel<0)
          ii = i+f1_l(1)
          jj = j+f1_l(2)
          kk = k+f1_l(3)
          !
          ! compute the integration extrema, xa,xb,ya,yb,za,zb
          !  note: they have to computed for each (i,j,k) since
          !        they depend on the local velocity and direction
          !
          cf = 0.5d0*(vel_sign+1.d0) ! i.e., cf = 1.0 (vel>=0) or 0.0 (vel<0)
          f2_l(2*dir-1+f1_l(dir),dir) = cf-dt*dli*vel(i,j,k) ! we touch only the integration direction
          f2_l(2*dir+0-f1_l(dir),dir) = cf
          !
          xa = crd(1)+f2_l(1,dir)
          xb = crd(2)+f2_l(2,dir)
          ya = crd(3)+f2_l(3,dir)
          yb = crd(4)+f2_l(4,dir)
          za = crd(5)+f2_l(5,dir)
          zb = crd(6)+f2_l(6,dir)
          !
          ! decide dir2
          !
          if((vof(ii,jj,kk).le.limit).or.(vof(ii,jj,kk).ge.1.d0-limit)) then
            !
            flux(i,j,k) = vof(ii,jj,kk)*(xb-xa)*(yb-ya)*(zb-za)
            !
          else
            !
            nor_v(:) = nor(ii,jj,kk,:)
            cur_v(:) = cur(ii,jj,kk,:)
            ind2 = maxloc(abs(nor_v(:)),1)
            !
            if(    ind2.eq.1) then
              !
              a    = xa
              b    = xb
              rm1  = 0.50d0*((yb-ya)*gm+(ya+yb))
              rp1  = 0.50d0*((yb-ya)*gp+(ya+yb))
              rm2  = 0.50d0*((zb-za)*gm+(za+zb))
              rp2  = 0.50d0*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (yb-ya)*(zb-za)
              !
              xv(1,:) = (/b,rm1,rm2/) 
              xv(2,:) = (/a,rm1,rm2/)
              xv(3,:) = (/b,rm1,rp2/)
              xv(4,:) = (/a,rm1,rp2/)
              xv(5,:) = (/b,rp1,rm2/)
              xv(6,:) = (/a,rp1,rm2/)
              xv(7,:) = (/b,rp1,rp2/)
              xv(8,:) = (/a,rp1,rp2/)
              !
            elseif(ind2.eq.2) then
              !
              a    = ya
              b    = yb
              rm1  = 0.50d0*((xb-xa)*gm+(xa+xb))
              rp1  = 0.50d0*((xb-xa)*gp+(xa+xb))
              rm2  = 0.50d0*((zb-za)*gm+(za+zb))
              rp2  = 0.50d0*((zb-za)*gp+(za+zb))
              flux(i,j,k) = (xb-xa)*(zb-za)
              !
              xv(1,:) = (/rm1,b,rm2/) 
              xv(2,:) = (/rm1,a,rm2/)
              xv(3,:) = (/rm1,b,rp2/)
              xv(4,:) = (/rm1,a,rp2/)
              xv(5,:) = (/rp1,b,rm2/)
              xv(6,:) = (/rp1,a,rm2/)
              xv(7,:) = (/rp1,b,rp2/)
              xv(8,:) = (/rp1,a,rp2/)
              !
            elseif(ind2.eq.3) then
              !
              a    = za
              b    = zb
              rm1  = 0.50d0*((xb-xa)*gm+(xa+xb))
              rp1  = 0.50d0*((xb-xa)*gp+(xa+xb))
              rm2  = 0.50d0*((yb-ya)*gm+(ya+yb))
              rp2  = 0.50d0*((yb-ya)*gp+(ya+yb))
              flux(i,j,k) = (xb-xa)*(yb-ya)
              !
              xv(1,:) = (/rm1,rm2,b/) 
              xv(2,:) = (/rm1,rm2,a/)
              xv(3,:) = (/rm1,rp2,b/)
              xv(4,:) = (/rm1,rp2,a/)
              xv(5,:) = (/rp1,rm2,b/)
              xv(6,:) = (/rp1,rm2,a/)
              xv(7,:) = (/rp1,rp2,b/)
              xv(8,:) = (/rp1,rp2,a/)
              !
            endif
            !
            cv(:)    = 1.d0
            cv(ind2) = 0.d0
            !
            cxxa = cv(1)*0.5d0*cv(1)*cur_v(1)
            cyya = cv(2)*0.5d0*cv(2)*cur_v(2)
            czza = cv(3)*0.5d0*cv(3)*cur_v(3)
            cxya = cv(1)*cv(2)*cur_v(4)
            cxza = cv(1)*cv(3)*cur_v(5)
            cyza = cv(2)*cv(3)*cur_v(6)
            a100 = (nor_v(1)-0.5d0*cv(1)*(cur_v(1)+cv(2)*cur_v(4)+cv(3)*cur_v(5)))
            a010 = (nor_v(2)-0.5d0*cv(2)*(cur_v(2)+cv(1)*cur_v(4)+cv(3)*cur_v(6)))
            a001 = (nor_v(3)-0.5d0*cv(3)*(cur_v(3)+cv(1)*cur_v(5)+cv(2)*cur_v(6)))
            !
            sum_int_func = 0.d0
            do p=1,7,2
              !
              q = p+1
              surf_a = cxxa*xv(p,1)*xv(p,1) + cyya*xv(p,2)*xv(p,2) + czza*xv(p,3)*xv(p,3) + &
                       cxya*xv(p,1)*xv(p,2) + cyza*xv(p,2)*xv(p,3) + cxza*xv(p,3)*xv(p,1) + &
                       a100*xv(p,1)+a010*xv(p,2)+a001*xv(p,3)
              surf_b = cxxa*xv(q,1)*xv(q,1) + cyya*xv(q,2)*xv(q,2) + czza*xv(q,3)*xv(q,3) + &
                       cxya*xv(q,1)*xv(q,2) + cyza*xv(q,2)*xv(q,3) + cxza*xv(q,3)*xv(q,1) + &
                       a100*xv(q,1)+a010*xv(q,2)+a001*xv(q,3)
              !
              sum_int_func = sum_int_func + &
                             (b-a+b_thi/(nor_v(ind2))*log( &
                             cosh(b_th*(surf_a+d_thinc(ii,jj,kk))) / &
                             cosh(b_th*(surf_b+d_thinc(ii,jj,kk)))))
              !
            enddo
            !
            flux(i,j,k) = 0.125d0*flux(i,j,k)*sum_int_func
            !
          endif
          !
          flux(i,j,k) = flux(i,j,k)*vel_sign*dl
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine cmpt_vof_flux
  !
  subroutine solve_quar_paper(a0,a1,a2,a3,a4,x)
    !
    ! subroutine to solve the quartic equation
    !
    implicit none
    !
    real(rp), intent(in ) :: a0,a1,a2,a3,a4
    real(rp), intent(out) :: x
    !
    real(rp) :: b3,b2,b1,b0,a4i
    real(rp) :: c2,c1,c0
    real(rp) :: z1,check,a,b
    real(rp) :: aa,bb,cc,dd,dd_s
    ! 
    ! First change the coefficients in the form of Eq. (B.2)
    !
    a4i = 1._rp/a4
    b3  = a3*a4i
    b2  = a2*a4i
    b1  = a1*a4i
    b0  = a0*a4i
    !
    ! Calculate the coefficients of Eq. (B.5)
    !  note: in c1 it is not b2 (as in the paper), but b3!
    !
    c2  = -b2
    c1  = b1*b3-4._rp*b0
    c0  = b0*(4._rp*b2-b3**2)-b1**2
    !
    ! Calculate z1, Eq. (B.7)
    !
    a   = -c2**2/9._rp+c1/3._rp
    b   = 2._rp*c2**3/27._rp-c1*c2/3._rp+c0
    check = b**2+4._rp*a**3
    !
    if(check.ge.0._rp) then
      z1 = sign(1._rp,(-b+sqrt(check)))*abs(0.5_rp*(-b+sqrt(check)))**(1._rp/3._rp) + &
           sign(1._rp,(-b-sqrt(check)))*abs(0.5_rp*(-b-sqrt(check)))**(1._rp/3._rp) - &
           c2/3._rp
    else
      z1 = 2._rp*sqrt(-a)*cos(atan(sqrt(-check)/(-b))/3._rp)-c2/3._rp
    endif
    !
    ! Find new coefficients, Eq. (B.9)
    !
    aa   = 0.5_rp*b3
    bb   = 0.5_rp*z1
    dd_s = bb**2-b0 ! as said in the paper, this should be always positive but still we put a check
    if(dd_s.le.0._rp) then
      dd = limit    ! to avoid singularity in case of imaginary solution (which we filter out)
    else
      dd = sqrt(dd_s)
    endif
    cc = (-0.5_rp*b1+aa*bb)/dd
    !
    ! Finally, calculate solution from Eq. (B.11)
    ! Be aware of a small error in the paper, i.e. 
    ! the (+) sign inside the square root should be (-)
    !
    x = 0.5_rp*(-(aa-cc) + sqrt((aa-cc)**2-4._rp*(bb-dd)))
    !
    return
  end subroutine solve_quar_paper
  !
  subroutine solve_quar(e,d,c,b,a,x)
    !
    ! subroutine to solve the quartic equation
    !
    implicit none
    !
    real(rp), intent(in ) :: e,d,c,b,a
    real(rp), intent(out) :: x
    !
    real(rp)   :: p,q,d0,d1
    complex(8) :: sqrt1,sqrt2,sqrt3,sqrt4,ss,qq,qt
    complex(8), dimension(4) :: xx
    ! 
    ! parameters for the quartic equation
    !
    x     = 0.d0
    p     = (8.*a*c-3.*b**2)/(8.*a**2)
    q     = (b**3-4.*a*b*c+8.*a**2*d)/(8.*a**3)
    d0    = c**2-3.*b*d+12.*a*e
    d1    = 2.*c**3-9.*b*c*d+27.*b**2*e+27.*a*d**2-72*a*c*e
    !
    ! root #1
    !
    sqrt1 = d1**2-4*d0**3
    qt    = 0.5d0*(d1+sqrt(sqrt1))
    qq    = (qt)**(1./3.)
    !
    ! root #2
    !
    sqrt2 = -2./3.*p+(qq+d0/qq)/(3.*a)
    ss    = 0.5d0*sqrt(sqrt2)
    !
    ! root #3
    !
    sqrt3 = -4.*ss**2-2.*p+q/ss
    xx(1) = -b/(4.*a)-ss+0.5*sqrt(sqrt3)
    if((abs((aimag(xx(1)))).eq.0.)) x = (real(xx(1)))
    xx(2) = -b/(4.*a)-ss-0.5*sqrt(sqrt3)
    if((abs((aimag(xx(2)))).eq.0.).and.((real(xx(2))).gt.x)) x=(real(xx(2)))
    !
    ! root #4
    !
    sqrt4 = -4.*ss**2-2.*p-q/ss
    xx(3) = -b/(4.*a)+ss+0.5*sqrt(sqrt4)
    if((abs((aimag(xx(3)))).eq.0.).and.((real(xx(3))).gt.x)) x=(real(xx(3)))
    xx(4) = -b/(4.*a)+ss-0.5*sqrt(sqrt4)
    if((abs((aimag(xx(4)))).eq.0.).and.((real(xx(4))).gt.x)) x=(real(xx(4)))
    !
    if(x.lt.0.d0) then
      !
      if(abs((aimag(xx(1)))).lt.limit) xx(1)=(real(xx(1)))
      if((abs((aimag(xx(1)))).eq.0.)) x=(real(xx(1)))
      if(abs((aimag(xx(2)))).lt.limit) xx(2)=(real(xx(2)))
      if((abs((aimag(xx(2)))).eq.0.).and.((real(xx(2))).gt.x)) x=(real(xx(2)))
      if(abs((aimag(xx(3)))).lt.limit) xx(3)=(real(xx(3)))
      if((abs((aimag(xx(3)))).eq.0.).and.((real(xx(3))).gt.x)) x=(real(xx(3)))
      if(abs((aimag(xx(4)))).lt.limit) xx(4)=(real(xx(4)))
      if((abs((aimag(xx(4)))).eq.0.).and.((real(xx(4))).gt.x)) x=(real(xx(4)))
      ! 
      if(x.lt.0.) then
        write(*,*) 'OPS THINC 3D 1',xx(1),xx(2),xx(3),xx(4)
      endif
      !
    endif
    ! 
    if(x.eq.0.d0) then
      !
      !  write(*,*) 'ZER1',xx(1),xx(2),xx(3),xx(4)
      if(abs((aimag(xx(1)))).lt.limit) xx(1)=(real(xx(1)))
      if((abs((aimag(xx(1)))).eq.0.)) x=(real(xx(1)))
      if(abs((aimag(xx(2)))).lt.limit) xx(2)=(real(xx(2)))
      if((abs((aimag(xx(2)))).eq.0.).and.((real(xx(2))).gt.x)) x=(real(xx(2)))
      if(abs((aimag(xx(3)))).lt.limit) xx(3)=(real(xx(3)))
      if((abs((aimag(xx(3)))).eq.0.).and.((real(xx(3))).gt.x)) x=(real(xx(3)))
      if(abs((aimag(xx(4)))).lt.limit) xx(4)=(real(xx(4)))
      if((abs((aimag(xx(4)))).eq.0.).and.((real(xx(4))).gt.x)) x=(real(xx(4)))
      ! 
      if(x.eq.0.d0) then
        ! write(*,*) 'OPS THINC 3D 2',xx(1),xx(2),xx(3),xx(4)
        !
        if(abs((aimag(xx(1)))).lt.limitf) xx(1)=(real(xx(1)))
        if((abs((aimag(xx(1)))).eq.0.)) x=(real(xx(1)))
        if(abs((aimag(xx(2)))).lt.limitf) xx(2)=(real(xx(2)))
        if((abs((aimag(xx(2)))).eq.0.).and.((real(xx(2))).gt.x)) x=(real(xx(2)))
        if(abs((aimag(xx(3)))).lt.limitf) xx(3)=(real(xx(3)))
        if((abs((aimag(xx(3)))).eq.0.).and.((real(xx(3))).gt.x)) x=(real(xx(3)))
        if(abs((aimag(xx(4)))).lt.limitf) xx(4)=(real(xx(4)))
        if((abs((aimag(xx(4)))).eq.0.).and.((real(xx(4))).gt.x)) x=(real(xx(4)))
        !
        if(x.eq.0.d0) then
          write(*,*) 'OPS THINC 3D 2',xx(1),xx(2),xx(3),xx(4)
          x = -1000.d0
        endif
        !
      endif
      !
    endif
    !
    return
  end subroutine solve_quar
  !
  subroutine solve_quad(a0,a1,a2,x)
    !
    ! subroutine to solve the quadratic equation
    !
    implicit none
    !
    real(rp), intent(in ) :: a0,a1,a2
    real(rp), intent(out) :: x
    !
    real(rp)              :: x1,x2
    ! 
    x1 = 0.5d0*(-a1+sqrt(a1**2-4.d0*a2*a0))/a2
    x2 = 0.5d0*(-a1-sqrt(a1**2-4.d0*a2*a0))/a2
    x  = max(x1,x2)
    !
    return
  end subroutine solve_quad
  !
end module mod_vof
