module mod_rk
  !
  use mod_common_mpi
  use mod_mom
  use mod_param, only: gacc,cbcpre,lx,ly,lz,turb_type,abc,k0_wave,nu2,k0_freq,turb_type,c_or_t, &
                       amp_a,amp_b,amp_n
  !
  implicit none
  !
  public  :: ab2_mom
  private 
  !
  contains
  !
  subroutine ab2_mom(istep,is_first,n,dli,dzci,dzfi,dt,dto,u,v,w,tmp,rho_gas,pold,pp,curv,vof,mu,rho,rho0, &
                                              dudtrko,dvdtrko,dwdtrko,up,vp,wp,time)
    !
    ! second order Adams-Bashforth scheme
    ! for time integration of the momentum equation
    !
    implicit none
    ! 
    integer, intent(in   )                         :: istep
    logical, intent(inout)                         :: is_first
    integer, intent(in   ), dimension(3)           :: n
    real(8), intent(in   )                         :: dt,dto
    real(8), intent(in   ), dimension(3)           :: dli
    real(8), intent(in   ), dimension(-2:)         :: dzci,dzfi
    real(8), intent(in   ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in   ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: rho_gas,pold,pp
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: curv,vof,mu,rho
    real(8), intent(in   )                         :: rho0
    real(8), intent(inout), dimension( 0:, 0:, 0:) :: dudtrko,dvdtrko,dwdtrko
    real(8), intent(out  ), dimension( 0:, 0:, 0:) :: up,vp,wp
    real(8), intent(out  )                         :: time
    !
    real(8), dimension(n(1),n(2),n(3)) :: dudtrk,dvdtrk,dwdtrk 
    real(8) :: factor1,factor2,factor12
    real(8) :: rho_av
    integer :: i,j,k
    !
    factor1 = dt*(1.d0+0.5d0*(dt/dto))
    factor2 = dt*(    -0.5d0*(dt/dto))
    if(is_first) then
      factor1 = dt*(1.d0)
      factor2 = dt*(0.d0)
      is_first = .false.
    endif
    factor12 = factor1+factor2
    !
    ! calculate the average density (to subtract the net gravitational force per unit mass).
    ! For a physical and a numerical explanation see the 
    ! lecture notes (DNS-Solver.pdf) of prof. Gretar Tryggvason
    ! eq. 4.15, pag. 44 at the following links:
    ! https://www3.nd.edu/~gtryggva/MultiphaseDNS/index.html
    ! https://www3.nd.edu/~gtryggva/MultiphaseDNS/DNS-Solver.pdf
    !
    rho_av = 0.d0
    if(((cbcpre(0,1)//cbcpre(1,1).eq.'PP').and.gacc(1).ne.0.d0).or. &
       ((cbcpre(0,2)//cbcpre(1,2).eq.'PP').and.gacc(2).ne.0.d0).or. &
       ((cbcpre(0,3)//cbcpre(1,3).eq.'PP').and.gacc(3).ne.0.d0)     ) then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            !
            rho_av = rho_av + rho(i,j,k)/(dli(1)*dli(2)*dzfi(k))
            !
          enddo
        enddo
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE,rho_av,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      rho_av = rho_av/(lx*ly*lz)
    endif
    !
#ifdef TWOD
    dudtrk(1:n(1),1:n(2),1:n(3)) = 0.d0
#else
    call momxad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,u,v,w,mu,rho,dudtrk)
#endif
    call momyad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,u,v,w,mu,rho,dvdtrk)
    call momzad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,u,v,w,mu,rho,dwdtrk)
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          up(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k) 
          wp(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          !    
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
    !
#ifdef TWOD
    dudtrk(1:n(1),1:n(2),1:n(3)) = 0.d0
#else
    call momxp(n(1),n(2),n(3),dli(1),pold,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dudtrk)
#endif
    call momyp(n(1),n(2),n(3),dli(2),pold,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dvdtrk)
    call momzp(n(1),n(2),n(3),dzci  ,pold,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dwdtrk)
    !
#ifdef TURB_FORCING
    !call turb_forc_src(turb_type,c_or_t,n(1),n(2),n(3),1.d0/dli(1),1.d0/dli(2),1.d0/dli(3),lx,ly,lz, &
    !                   abc(1),abc(2),abc(3),f0_t,k0_freq,rho,dudtrk,dvdtrk,dwdtrk)
    !if(istep.lt.3000000) then 
    !if(mod(istep,icheck).eq.0) then   
    call turb_forc_src(turb_type,c_or_t,n(1),n(2),n(3),1.d0/dli(1),1.d0/dli(2),1.d0/dli(3),lx,ly,lz, &
                       abc(1),abc(2),abc(3),amp_a,amp_b,amp_n, &
                       time,nu2,k0_wave,k0_freq,rho,dudtrk,dvdtrk,dwdtrk)
    !endif
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          up(i,j,k) = up(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor12*dwdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine ab2_mom
  !
end module mod_rk
