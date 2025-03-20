module mod_initflow
  !
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp
  use mod_common_mpi, only: ierr,coord,myid
  use mod_param     , only: dims,pi,dx,dy,dz,lx,ly,lz,tmpl0,tmpg0,sinit,uref, &
                            pi,k0_freq,abc,amp_a,amp_b,amp_n
  use mod_output    , only: out0d
  use mod_rks       , only: ab2_scal1
  use mod_thermo    , only: mass_fraction
  !
  implicit none
  !
  private
  public  :: initvel,initsca,inittmp,add_noise
  !
  contains
  !
  subroutine initvel(inivel,n,zclzi,dzclzi,dzflzi,norm,norm_u,time,vof,u,v,w,p)
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    !
    character(len=3), intent(in)                          :: inivel
    integer         , intent(in ), dimension(3)           :: n
    real(8)         , intent(in ), dimension(-2:)         :: zclzi,dzclzi,dzflzi
    real(8)         , intent(in )                         :: norm,norm_u,time
    real(8)         , intent(in ), dimension( 0:, 0:, 0:) :: vof
    real(8)         , intent(out), dimension(-2:,-2:,-2:) :: u,v,w
    real(8)         , intent(out), dimension( 0:, 0:, 0:) :: p
    !
    real(8), allocatable, dimension(:) :: u1d
    !real(8), allocatable, dimension(:,:) :: u2d
    integer :: i,j,k
    real(8) :: q,amp,t_per
    logical :: is_noise,is_mean,is_pair
    real(8) :: xc,yc,zc,xf,yf,zf
    !
    t_per = lx/(2.d0*pi*abc(1))
    !
    allocate(u1d(n(3)))
    is_noise = .false.
    !is_noise = .true.
    is_mean  = .false.
    !is_pair  = .true.
    is_pair  = .false.
    q = .5d0
    select case(inivel)
    case('cou')
      call couette(   q,n(3),zclzi,norm,u1d)
    case('poi')
      call poiseuille(q,n(3),zclzi,norm,u1d)
      is_mean=.true.
    case('zer')
      u1d(:) = 0.
!    case('log')
!      call log_profile(q,n(3),zclzi,visc,u1d)
!      is_noise = .true.
!      is_mean = .true.
!    case('hcl')
!      deallocate(u1d)
!      allocate(u1d(2*n(3)))
!      call log_profile(q,2*n(3),zclzi,visc,u1d)
!      is_noise = .true.
!      is_mean=.true.
    case('hcp')
      deallocate(u1d)
      allocate(u1d(2*n(3)))
      call poiseuille(q,2*n(3),zclzi,norm,u1d)
      is_mean = .true.
    case('tgv')
      do k=1,n(3)
        zc = zclzi(k)*2.d0*pi
        do j=1,n(2)
          yc = (j+coord(2)*n(2)-.5d0)*dy/ly*2.d0*pi
          yf = (j+coord(2)*n(2)-.0d0)*dy/ly*2.d0*pi
          do i=1,n(1)
            xc = (i+coord(1)*n(1)-.5d0)*dx/lx*2.d0*pi
            xf = (i+coord(1)*n(1)-.0d0)*dx/lx*2.d0*pi
            u(i,j,k) =  sin(xf)*cos(yc)*cos(zc)
            v(i,j,k) = -cos(xc)*sin(yf)*cos(zc)
            w(i,j,k) = 0.d0
            p(i,j,k) = 0.d0!(cos(2.d0*xc)+cos(2.d0*yc))*(cos(2.d0*zc)+2.d0)/16.d0
          enddo
        enddo
      enddo
    case('khi')
      do k=1,n(3)
        zc = 1.*k*dz-0.5*dz - 0.5*lz
        zf = 1.*k*dz-0.0*dz - 0.5*lz
        do j=1,n(2)
          yc = 1.*(j+coord(2)*n(2))*dy-0.5*dy -0.5*ly
          yf = 1.*(j+coord(2)*n(2))*dy-0.0*dy -0.5*ly
          do i=1,n(1)
            !
            u(i,j,k) = 0.d0
            v(i,j,k) = tanh(2.*zc*28.0)+&
                       0.001*(exp(-zc**2*14.0**2)*(cos(4.0*pi*yf)+0.1*sin(10.0*pi*yf))*(-2.0*zc*28.0**2))
            w(i,j,k) = -0.001*(-exp(-zf**2*28.0**2)*(sin(4.0*pi*yc)*4.0*pi-0.1*cos(10.0*pi*yc)*10.0*pi))
            p(i,j,k) = 0.d0
            !
          enddo
        enddo
      enddo
    case('mer')
      do k=1,n(3)
        zc = (k-.5d0)*dz 
        zf = (k-.0d0)*dz 
        do j=1,n(2)
          yc = (j+coord(2)*n(2)-0.5d0)*dy 
          yf = (j+coord(2)*n(2)-0.0d0)*dy 
          do i=1,n(1)
            xc = (i+coord(1)*n(1)-0.5d0)*dx 
            xf = (i+coord(1)*n(1)-0.0d0)*dx 
             if(yf.lt.0.5d0*ly) then
               u(i,j,k) = 0.d0
               v(i,j,k) = +uref*0.5d0*(vof(i,j,k)+vof(i,j+1,k))
               w(i,j,k) = 0.d0
             else
               u(i,j,k) = 0.d0
               !v(i,j,k) = 0.01d0*uref*0.5d0*(vof(i,j,k)+vof(i,j+1,k))
               v(i,j,k) = -uref*0.5d0*(vof(i,j,k)+vof(i,j+1,k))
               w(i,j,k) = 0.d0
             endif
          enddo
        enddo
      enddo
    case('abc')
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            zc = (k              -0.5)*dz/lz*2.*pi
            yc = (j+coord(2)*n(2)-0.5)*dy/ly*2.*pi
            yf = (j+coord(2)*n(2)-0.0)*dy/ly*2.*pi
            xc = (i+coord(1)*n(1)-0.5)*dx/lx*2.*pi
            xf = (i+coord(1)*n(1)-0.0)*dx/lx*2.*pi
            !
            amp = 1.0!(amp_a+(amp_b/2.d0)*sin((2.d0*pi/(amp_n*t_per))*time))
            !amp = (amp_a+(amp_b/2.d0)*cos((2.d0*pi/(amp_n*t_per))*time))
            u(i,j,k) = (abc(1)*sin(1.d0*k0_freq*zc) + abc(3)*cos(1.d0*k0_freq*yc))*amp
            v(i,j,k) = (abc(2)*sin(1.d0*k0_freq*xc) + abc(1)*cos(1.d0*k0_freq*zc))*amp
            w(i,j,k) = (abc(3)*sin(1.d0*k0_freq*yc) + abc(2)*cos(1.d0*k0_freq*xc))*amp
            p(i,j,k) = 0.0!(cos(2._rp*xc)+cos(2._rp*yc))*(cos(2._rp*zc)+2._rp)/16._rp
            !
          enddo
        enddo
      enddo
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial velocity field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to ierrs in the case file ***'
      if(myid.eq.0) print*, '    check setup.h90'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    if(inivel.ne.'tgv'.and.inivel.ne.'mer'.and.inivel.ne.'khi'.and.inivel.ne.'abc') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0.d0
            w(i,j,k) = 0.d0
            p(i,j,k) = 0.d0
          enddo
        enddo
      enddo
    endif
    !
    if(is_mean) then
      call set_mean(n,1.d0,dzflzi,u(1:n(1),1:n(2),1:n(3)))
    endif
    if(is_pair) then
      !
      ! initialize a streamwise vortex pair for a fast transition
      ! to turbulence in a pressure-driven channel:
      !        psi(x,y,z)  = f(z)*g(x,y), with
      !        f(z)        = (1-z**2)**2, and
      !        g(x,y)      = y*exp[-(16x**2-4y**2)]
      ! (x,y,z) --> (streamwise, spanwise, wall-normal) directions
      !
      ! see Henningson and Kim, JFM 1991
      !
      !do k=1,n(3)
      !  zc = 2.d0*zclzi(k) - 1.d0 ! z rescaled to be between -1 and +1
      !  zf = 2.d0*(zclzi(k) + .5d0*dzflzi(k)) - 1.d0
      !  do j=1,n(2)
      !    yc = ((coord(2)*n(2)+j-0.5)*dy-.5d0*ly)*2.d0/lz
      !    yf = ((coord(2)*n(2)+j-0.0)*dy-.5d0*ly)*2.d0/lz
      !    do i=1,n(1)
      !      xc = ((coord(1)*n(1)+i-0.5)*dx-.5d0*lx)*2.d0/lz
      !      xf = ((coord(1)*n(1)+i-0.0)*dx-.5d0*lx)*2.d0/lz
      !      u(i,j,k) = u1d(k)
      !      v(i,j,k) =  1.d0 * fz(zc)*dgxy(yf,xc)*uref*1.5d0
      !      w(i,j,k) = -1.d0 * gxy(yc,xc)*dfz(zf)*uref*1.5d0
      !      p(i,j,k) = 0.d0
      !    enddo
      !  enddo
      !enddo
      !
      ! Taylor-Green vortices
      !
      !do k=1,n(3)
      !  zc = (k-0.5d0)*dz/lz*2.*pi
      !  zf = (k-0.0d0)*dz/lz*2.*pi
      !  do j=1,n(2)
      !    yc = (j+coord(2)-0.50)*dy/ly*2.*pi
      !    yf = (j+coord(2)-0.00)*dy/ly*2.*pi
      !    do i=1,n(1)
      !      xc = (i+coord(1)-0.5d0)*dx/lx*2.*pi
      !      xf = (i+coord(1)-0.0d0)*dx/lx*2.*pi
      !      u(i,j,k) = 0.d0
      !      v(i,j,k) =  sin(xc)*cos(yf)*cos(zc)*uref
      !      w(i,j,k) = -cos(xc)*sin(yc)*cos(zf)*uref
      !      p(i,j,k) = 0.!(cos(2.*xc)+cos(2.*yc))*(cos(2.*zc)+2.)/16.
      !    enddo
      !  enddo
      !enddo
      !
      ! Andreas trick to promote transition in RB
      !
      do k=1,n(3)/2
        zc = 4.d0*zclzi(k) - 1.d0 ! z rescaled to be between -1 and +1
        zf = 4.d0*(zclzi(k) + .5d0*dzflzi(k)) - 1.d0
        do j=1,n(2)
          yc = ((coord(2)*n(2)+j-0.5)*dy-.5d0*ly)*2.d0/lz
          yf = ((coord(2)*n(2)+j-0.0)*dy-.5d0*ly)*2.d0/lz
          do i=1,n(1)
            xc = ((coord(1)*n(1)+1-0.5)*dx-.5d0*lx)*2.d0/lz
            xf = ((coord(1)*n(1)+1-0.0)*dx-.5d0*lx)*2.d0/lz
            !
            u(i,j,k) = u1d(k)
            v(i,j,k) =  norm_u * fz(zc)*dgxy(yf,xc)*uref*1.5d0
            w(i,j,k) = -norm_u * gxy(yc,xc)*dfz(zf)*uref*1.5d0
            p(i,j,k) = 0.d0
            !
          enddo
        enddo
      enddo
      !
      do k=n(3)/2,n(3)
        zc = 4.d0*zclzi(k) - 3.d0 ! z rescaled to be between -1 and +1
        zf = 4.d0*(zclzi(k) + .5d0*dzflzi(k)) - 3.d0
        do j=1,n(2)
          yc = ((coord(2)*n(2)+j-0.5)*dy-.5d0*ly)*2.d0/lz
          yf = ((coord(2)*n(2)+j-0.0)*dy-.5d0*ly)*2.d0/lz
          do i=1,n(1)
            xc = ((coord(1)*n(1)+1-0.5)*dx-.5d0*lx)*2.d0/lz
            xf = ((coord(1)*n(1)+1-0.0)*dx-.5d0*lx)*2.d0/lz
            !
            u(i,j,k) = u1d(k)
            v(i,j,k) =  norm_u * fz(zc)*dgxy(yf,xc)*uref*1.5
            w(i,j,k) = -norm_u * gxy(yc,xc)*dfz(zf)*uref*1.5
            p(i,j,k) = 0.d0
            !
          enddo
        enddo
      enddo
      !
    endif
    !
    ! Add noise to the final velocity to further trigger
    ! transition
    !
    if(is_noise) then
#ifdef TWOD
      u(:,:,:) = 0.d0
#else
      call add_noise(n,123,norm_u*uref,u(1:n(1),1:n(2),1:n(3)))
#endif
      call add_noise(n,456,norm_u*uref,v(1:n(1),1:n(2),1:n(3)))
      call add_noise(n,789,norm_u*uref,w(1:n(1),1:n(2),1:n(3)))
    endif
    !
    deallocate(u1d)
    return
  end subroutine initvel
  !
  subroutine inittmp(initmp,n,dli,phi,vof,tmp)
    !
    ! computes initial conditions for the temperature field
    !
    implicit none
    !
    character(len=3), intent(in )                         :: initmp
    integer         , intent(in ), dimension(3)           :: n
    real(8)         , intent(in ), dimension(3)           :: dli
    real(8)         , intent(in ), dimension(-2:,-2:,-2:) :: phi
    real(8)         , intent(in ), dimension( 0:, 0:, 0:) :: vof
    real(8)         , intent(out), dimension(-2:,-2:,-2:) :: tmp
    !
    real(8) :: zc
    integer :: i,j,k
    !
    select case(initmp)
    case('uni')
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            tmp(i,j,k) = tmpl0
          enddo
        enddo
      enddo
    case('sin')
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            tmp(i,j,k) = tmpl0*vof(i,j,k)+tmpg0*(1.d0-vof(i,j,k))
            !
          enddo
        enddo
      enddo
    case('lin')
      !do k=1,n(3)
      !  zc = (k-0.5d0)*dz
      !  do j=1,n(2)
      !    do i=1,n(1)
      !      !
      !      tmp(i,j,k) = tb_b + ((tb_t-tb_b)/lz)*(zc-0.d0)
      !      !
      !    enddo
      !  enddo
      !enddo
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial temperature field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to ierrs in the case file ***'
      if(myid.eq.0) print*, '    check setup.h90'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    !if(noise) then
    ! call add_noise(n,123,norm_t,tmp(1:n(1),1:n(2),1:n(3)))
    !endif
    !
#ifdef TWOD
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          tmp(i,j,k) = tmp(n(1)/2,j,k) ! this ensure (u,v,w)=(0,v,w)
        enddo
      enddo
    enddo
#endif
    !
    return
  end subroutine inittmp
  !
  subroutine initsca(inisca,n,dli,dzci,dzfi,cbc,bc,rhog,d_m,pth,u,v,w,phi,tmpge,tmple,sca)
    !
    ! computes initial conditions for the scalar field
    !
    implicit none
    !
    character(len=3), intent(in )                         :: inisca
    integer         , intent(in ), dimension(3)           :: n
    real(8)         , intent(in ), dimension(3)           :: dli
    real(8)         , intent(in ), dimension(-2:)         :: dzci,dzfi
    character(len=1), intent(in ), dimension(0:1,3)       :: cbc
    real(8)         , intent(in ), dimension(0:1,3)       :: bc
    real(8)         , intent(in ), dimension( 0:, 0:, 0:) :: rhog,d_m
    real(8)         , intent(in )                         :: pth
    real(8)         , intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8)         , intent(in ), dimension(-2:,-2:,-2:) :: phi
    real(8)         , intent(in ), dimension( 0:, 0:, 0:) :: tmpge,tmple
    real(8)         , intent(out), dimension( 0:, 0:, 0:) :: sca
    !
    integer :: i,j,k,istep
#ifdef IMPDIFF
    integer, parameter :: nstep_sca = 20
#else
    integer, parameter :: nstep_sca = 50000
#endif
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: scae
    real(8) :: sca_av,dt,d_mMax
    real(8), dimension(3) :: dl
    real(8), dimension(-2:n(3)+3) :: dzc,dzf
    !
    dl  = dli**(-1)
    dzc = dzci**(-1)
    dzf = dzfi**(-1)
    !
    select case(inisca)
    case('std')
      !
      ! a. preinitialization
      !
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            if(phi(i,j,k).gt.0.d0) then
              sca(i,j,k)  = mass_fraction(pth,tmpl0)
              scae(i,j,k) = mass_fraction(pth,tmpl0)
            else
              sca(i,j,k)  = sinit
              scae(i,j,k) = sinit
            endif
          enddo
        enddo
      enddo
      call boundp(cbc,n,bc,dl,dzc,dzf,sca )
      call boundp(cbc,n,bc,dl,dzc,dzf,scae)
      !
      ! b. calculation of the steady-state solution
      !
      if(myid.eq.0) print*, 'calculation of the steady state solution'
#ifdef IMPDIFF
      dt = 100.d0
#else
      d_mMax = maxval(d_m(1:n(1),1:n(2),1:n(3)))
      call mpi_allreduce(MPI_IN_PLACE,d_mMax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      dt = 0.010d0/(min(dli(1)**2,dli(2)**2,dli(3)**2)*d_mMax)
#endif
      do istep=1,nstep_sca
        !
        if(myid.eq.0) print*, 'Step = ', istep
        !
        call ab2_scal1(n,dli,dzci,dzfi,cbc,bc,dt,rhog,d_m,pth,u,v,w,phi,tmpge,tmple,scae,sca)
        call boundp(cbc,n,bc,dl,dzc,dzf,sca)
        sca_av = (1.d0/(1.d0*n(1)*dims(1)*n(2)*dims(2)*n(3)))*sum(sca(1:n(1),1:n(2),1:n(3)))
        call mpi_allreduce(MPI_IN_PLACE,sca_av,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        if(myid.eq.0) then
          call out0d('data/sca_av.out',2,(/1.d0*istep,sca_av/))
        endif
      !
      enddo
    case default
      if(myid.eq.0) print*, 'ERROR: invalid name for initial scalar field'
      if(myid.eq.0) print*, ''
      if(myid.eq.0) print*, '*** Simulation abortited due to ierrs in the case file ***'
      if(myid.eq.0) print*, '    check setup.h90'
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)
      call exit
    end select
    !
    return
  end subroutine initsca
  !
  subroutine add_noise(n,iseed,norm,p)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)              :: n
    integer, intent(in   )                            :: iseed
    real(8), intent(in   )                            :: norm
    real(8), intent(inout), dimension(n(1),n(2),n(3)) :: p
    !
    integer(4), allocatable, dimension(:) :: seed
    real(8) :: rn
    integer, dimension(3) :: ng
    integer :: i,j,k,ii,jj
    !
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    ng(:) = n(:)
    ng(1:2) = ng(1:2)*dims(1:2)
    do k=1,ng(3)
      do j=1,ng(2)
        jj = j-coord(2)*n(2)
        do i=1,ng(1)
          ii = i-coord(1)*n(1)
          call random_number(rn)
          if(ii.ge.1.and.ii.le.n(1) .and. &
             jj.ge.1.and.jj.le.n(2) ) then
             p(ii,jj,k) = p(ii,jj,k) + 2.d0*(rn-.5d0)*norm
          endif
        enddo
      enddo
    enddo
    !
    return
  end subroutine add_noise
  !
  subroutine set_mean(n,mean,dzlzi,p)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)               :: n
    real(8), intent(in   ), dimension(-2:)             :: dzlzi 
    real(8), intent(in   )                             :: mean
    real(8), intent(inout), dimension(n(1),n(2),n(3))  :: p
    !
    real(8) :: meanold
    integer :: i,j,k
    meanold = 0.d0
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:meanold)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          meanold = meanold + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,meanold,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    meanold = meanold/(1.d0*n(1)*dims(1)*n(2)*dims(2))
    !
    if(meanold.ne.0.d0) then
    !$OMP WORKSHARE
    p(:,:,:) = p(:,:,:)/meanold*mean
    !$OMP END WORKSHARE
    endif
    !
    return
  end subroutine set_mean
  !
  subroutine couette(q,n,zc,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    !
    real(8), intent(in )                 :: q
    integer, intent(in )                 :: n
    real(8), intent(in ), dimension(-2:) :: zc
    real(8), intent(in )                 :: norm
    real(8), intent(out), dimension(n)   :: p
    !
    integer :: k
    real(8) :: z
    do k=1,n
      z    = zc(k)!1.d0*((k-1)+q)/(1.d0*n)
      p(k) = .5d0*(1.d0-2.d0*z)/norm
    enddo
    !
    return
  end subroutine couette
  !
  subroutine poiseuille(q,n,zc,norm,p)
    !
    implicit none
    !
    real(8), intent(in )                 :: q
    integer, intent(in )                 :: n
    real(8), intent(in ), dimension(-2:) :: zc
    real(8), intent(in )                 :: norm
    real(8), intent(out), dimension(n)   :: p
    !
    integer :: k
    real(8) :: z
    !
    ! plane poiseuille profile normalized by the bulk velocity
    !
    do k=1,n
      z    = zc(k)!1.d0*((k-1)+q)/(1.d0*n)
      p(k) = 6.d0*z*(1.d0-z)/norm
    enddo
    !
    return
  end subroutine poiseuille
  !
  !subroutine log_profile(q,n,zc,visc,p)
  !  implicit none
  !  real(8), intent(in)   :: q
  !  integer, intent(in)   :: n
  !  real(8), intent(in), dimension(0:) :: zc
  !  real(8), intent(in)   :: visc
  !  real(8), intent(out), dimension(n) :: p
  !  integer :: k
  !  real(8) :: z,reb,retau ! z/lz and bulk Reynolds number
  !  reb = rey
  !  retau = 0.09*reb**(0.88) ! from Pope's book
  !  do k=1,n/2
  !    z    = zc(k)*2.*retau!1.d0*((k-1)+q)/(1.d0*n)*2.*retau
  !    p(k) = 2.5d0*log(z) + 5.5d0
  !    if (z.le.11.6d0) p(k)=z
  !    p(n+1-k) = p(k)
  !  enddo
  !  return
  !end subroutine log_profile
  !
  ! functions to initialize the streamwise vortex pair
  ! (explained above)
  !
  function fz(zc)
  real(8), intent(in) :: zc
  real(8) :: fz
    fz = ((1.d0-zc**2)**2)
  end function
  !
  function dfz(zc)
  real(8), intent(in) :: zc
  real(8) :: dfz
    dfz = -4.d0*zc*((1.d0-zc**2)**2)
  end function
  !
  function gxy(xc,yc)
  real(8), intent(in) :: xc,yc
  real(8) :: gxy
    gxy = yc*exp(-4.d0*(4.d0*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  real(8), intent(in) :: xc,yc
  real(8) :: dgxy
    dgxy = exp(-4.d0*(4.d0*xc**2+yc**2))*(1.d0-8.d0*yc**2)
  end function
  !
end module mod_initflow
