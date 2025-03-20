module mod_rks
  !
  use mod_moms      , only: momtad,rhs_pth,momsad_flm
  use mod_common_mpi, only: ierr,comm_cart
  use mpi
  !
  implicit none
  !
  private
  public  :: cmpt_pth,ab2_tmp,ab2_scal1
  !
  contains
  !
  subroutine cmpt_pth(n,dli,dt,vof,kappa,mflux,rho_gas,d_lg,rho1,tmp,tmpge,sca,scae,mgas,mfxt,pth,dpthdt_n,dpthdt)
    !
    ! compute the thermodynamic pressure
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)           :: n
    real(8), intent(in   ), dimension(3)           :: dli
    real(8), intent(in   )                         :: dt
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: vof,kappa,mflux
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: rho_gas,d_lg
    real(8), intent(in   )                         :: rho1
    real(8), intent(in   ), dimension(-2:,-2:,-2:) :: tmp 
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: tmpge
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: sca,scae
    real(8), intent(inout)                         :: mgas,mfxt
    real(8), intent(inout)                         :: pth,dpthdt_n
    real(8), intent(inout)                         :: dpthdt
    !
    real(8) :: mfxto,dpthdt_no,integral
    integer :: i,j,k
    !
    ! a. calculation of quantities at previous time step
    !
    mfxto     = mfxt
    dpthdt_no = dpthdt_n
    !
    ! b. calculation of the normalized time derivative of pth
    !
    call rhs_pth(n,dli,vof,mflux,kappa,rho_gas,d_lg,rho1,tmp,tmpge,sca,scae,pth,dpthdt_n,mfxt,integral)
    !
    ! c. calculation of pth (note: we use a trapezoidal rule to compute the integral in time of mflux)
    !
    !pth = pth*exp(dt*0.5d0*(dpthdt_n+dpthdt_no)) ! method 1
    !
    mgas = mgas + (dt*0.5d0*(mfxt+mfxto))         
    pth  = mgas/integral                          ! method 2 
    !
    ! d. calculation of dpthdt
    !
    dpthdt = pth*dpthdt_n
    !
    return
  end subroutine cmpt_pth
  !
  subroutine ab2_scal1(n,dli,dzci,dzfi,cbc,bc,dt,rhog,d_m,pth,u,v,w,phi,tmpge,tmple,scae,s)
    !
    ! first order Euler scheme
    ! for time integration of the scalar field.
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)           :: n
    real(8)         , intent(in   ), dimension(3)           :: dli
    real(8)         , intent(in   ), dimension(-2:)         :: dzci,dzfi
    character(len=1), intent(in   ), dimension(0:1,3)       :: cbc
    real(8)         , intent(in   ), dimension(0:1,3)       :: bc
    real(8)         , intent(in   )                         :: dt
    real(8)         , intent(in   ), dimension( 0:, 0:, 0:) :: rhog,d_m
    real(8)         , intent(in   )                         :: pth
    real(8)         , intent(in   ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8)         , intent(in   ), dimension(-2:,-2:,-2:) :: phi
    real(8)         , intent(in   ), dimension( 0:, 0:, 0:) :: tmpge,tmple
    real(8)         , intent(in   ), dimension( 0:, 0:, 0:) :: scae
    real(8)         , intent(inout), dimension( 0:, 0:, 0:) :: s
    !
    real(8), dimension(n(1),n(2),n(3)) :: dsdtrk
    real(8) :: factor1
    integer :: i,j,k
    !
    factor1 = 1.0d0*dt
    call momsad_flm(n,dli,+0,pth,rhog,d_m,u,v,w,s,phi,tmpge,tmple,dsdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(phi(i,j,k).gt.0.d0) then
            !s(i,j,k)      = 1.d0
            s(i,j,k)      = scae(i,j,k)
          else
            s(i,j,k)      = s(i,j,k) + dsdtrk(i,j,k)*factor1
#ifdef IMPDIFF
            dsdtrk(i,j,k) = s(i,j,k)
#endif
          endif
          !
        enddo
      enddo
    enddo
    !
#ifdef IMPDIFF
    call helmholtz_mg(n,dli,dzci,cbc,bc,-factor1,pth,phi,tmp,scae,rhog,d_m,dsdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(phi(i,j,k).gt.0.d0) then
            !s(i,j,k) = 1.d0
            s(i,j,k) = scae(i,j,k)
          else
            s(i,j,k) = dsdtrk(i,j,k)
          endif
          !
        enddo
      enddo
    enddo
#endif
    !
    return
  end subroutine ab2_scal1
  !
  subroutine ab2_tmp(is_first,n,dli,dt,dto,vof,mflux,rhog,d_lg,kappa,u,v,w,sca,tmp,dtmpdtrkold,dpthdto)
    !
    use mod_param, only: rho1,cp1,cp2,m1,m2,lheat,delta_cp
    !
    ! Adams-Bashforth scheme for temperature equation 
    !
    implicit none
    !
    logical, intent(inout)                         :: is_first
    integer, intent(in   ), dimension(3)           :: n
    real(8), intent(in   ), dimension(3)           :: dli
    real(8), intent(in   )                         :: dt,dto
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: vof,mflux
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: rhog,d_lg,kappa
    real(8), intent(in   ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in   ), dimension( 0:, 0:, 0:) :: sca
    real(8), intent(inout), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(inout), dimension( 0:, 0:, 0:) :: dtmpdtrkold
    real(8), intent(in   )                         :: dpthdto
    ! 
    real(8), dimension(n(1),n(2),n(3)) :: dtmpdtrk
    real(8), dimension(8) :: mx,my,mz
    real(8) :: factor1,factor2,factor12
    real(8) :: dvofdx,dvofdy,dvofdz,rhocpci,cpci,scr_te1,scr_te2,scr_te3
    real(8) :: dscadxp,dscadxm,dscadyp,dscadym,dscadzp,dscadzm, &
               dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm
    integer :: i,j,k,ip,jp,kp,im,jm,km
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
    call momtad(n,dli,u,v,w,rhog,vof,kappa,tmp,dtmpdtrk)
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
          dscadxp  = (sca(ip,j,k)-sca(i ,j,k))*dli(1)
          dscadxm  = (sca(i ,j,k)-sca(im,j,k))*dli(2)
          dscadyp  = (sca(i,jp,k)-sca(i,j ,k))*dli(2)
          dscadym  = (sca(i,j ,k)-sca(i,jm,k))*dli(2)
          dscadzp  = (sca(i,j,kp)-sca(i,j,k ))*dli(3)
          dscadzm  = (sca(i,j,k )-sca(i,j,km))*dli(3)
          !
          dtmpdxp  = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
          dtmpdxm  = (tmp(i ,j,k)-tmp(im,j,k))*dli(1)
          dtmpdyp  = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
          dtmpdym  = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
          dtmpdzp  = (tmp(i,j,kp)-tmp(i,j,k ))*dli(3)
          dtmpdzm  = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
          !
          !rhocpci = vof(i,j,k)/(rho1*cp1) + ( 1.d0-vof(i,j,k) )/(rhog(i,j,k)*cp2)
          !cpci    = vof(i,j,k)/(cp1     ) + ( 1.d0-vof(i,j,k) )/(cp2            )
          rhocpci = (vof(i,j,k)*(rho1*cp1) + ( 1.d0-vof(i,j,k) )*(rhog(i,j,k)*cp2))**(-1.d0)
          cpci    = (vof(i,j,k)*(cp1     ) + ( 1.d0-vof(i,j,k) )*(cp2            ))**(-1.d0)
          !
#ifdef TWOD
          mx(1:8) = 0.d0
#else
          !i+1/2 j+1/2 k+1/2
          mx(1)=((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mx(2)=((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k+1)+vof(i+1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1)))*dli(1)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mx(3)=((vof(i+1,j  ,k  )+vof(i+1,j+1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1)))*dli(1)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mx(4)=((vof(i+1,j  ,k  )+vof(i+1,j-1,k  )+vof(i+1,j  ,k-1)+vof(i+1,j-1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mx(5)=((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j+1,k+1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j+1,k+1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mx(6)=((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k+1)+vof(i  ,j-1,k+1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mx(7)=((vof(i  ,j  ,k  )+vof(i  ,j+1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j+1,k-1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j+1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j+1,k-1)))*dli(1)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mx(8)=((vof(i  ,j  ,k  )+vof(i  ,j-1,k  )+vof(i  ,j  ,k-1)+vof(i  ,j-1,k-1))-&
                 (vof(i-1,j  ,k  )+vof(i-1,j-1,k  )+vof(i-1,j  ,k-1)+vof(i-1,j-1,k-1)))*dli(1)*0.25d0
#endif
          !
          !i+1/2 j+1/2 k+1/2
          my(1)=((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k+1/2
          my(2)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1))-&
                 (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1)))*dli(2)*0.25d0
          !i+1/2 j+1/2 k-1/2
          my(3)=((vof(i  ,j+1,k  )+vof(i+1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)))*dli(2)*0.25d0
          !i+1/2 j-1/2 k-1/2
          my(4)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1))-&
                 (vof(i  ,j-1,k  )+vof(i+1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(2)*0.25d0
          !i-1/2 j+1/2 k+1/2
          my(5)=((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k+1/2
          my(6)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1))-&
                 (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1)))*dli(1)*0.25d0
          !i-1/2 j+1/2 k-1/2
          my(7)=((vof(i  ,j+1,k  )+vof(i-1,j+1,k  )+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)))*dli(2)*0.25d0
          !i-1/2 j-1/2 k-1/2
          my(8)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1))-&
                 (vof(i  ,j-1,k  )+vof(i-1,j-1,k  )+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(2)*0.25d0
          !
          !i+1/2 j+1/2 k+1/2
          mz(1)=((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i+1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  )))*dli(3)*0.25d0
          !i+1/2 j-1/2 k+1/2
          mz(2)=((vof(i  ,j  ,k+1)+vof(i+1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i+1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  )))*dli(3)*0.25d0
          !i+1/2 j+1/2 k-1/2
          mz(3)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j+1,k  )+vof(i+1,j+1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i+1,j+1,k-1)))*dli(3)*0.25d0
          !i+1/2 j-1/2 k-1/2
          mz(4)=((vof(i  ,j  ,k  )+vof(i+1,j  ,k  )+vof(i  ,j-1,k  )+vof(i+1,j-1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i+1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i+1,j-1,k-1)))*dli(3)*0.25d0
          !i-1/2 j+1/2 k+1/2
          mz(5)=((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j+1,k+1)+vof(i-1,j+1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  )))*dli(3)*0.25d0
          !i-1/2 j-1/2 k+1/2
          mz(6)=((vof(i  ,j  ,k+1)+vof(i-1,j  ,k+1)+vof(i  ,j-1,k+1)+vof(i-1,j-1,k+1))-&
                 (vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  )))*dli(3)*0.25d0
          !i-1/2 j+1/2 k-1/2
          mz(7)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j+1,k  )+vof(i-1,j+1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j+1,k-1)+vof(i-1,j+1,k-1)))*dli(3)*0.25d0
          !i-1/2 j-1/2 k-1/2
          mz(8)=((vof(i  ,j  ,k  )+vof(i-1,j  ,k  )+vof(i  ,j-1,k  )+vof(i-1,j-1,k  ))-&
                 (vof(i  ,j  ,k-1)+vof(i-1,j  ,k-1)+vof(i  ,j-1,k-1)+vof(i-1,j-1,k-1)))*dli(3)*0.25d0
          !
          dvofdx = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dvofdy = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dvofdz = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          !
          ! --> scr_te1 = source term due to phase change (only at the interface)
          ! --> scr_te2 = source term due to change in th. pre (only in the gas phase)
          ! --> scr_te3 = source term due to change in composition, "enthalpy diffusion" (only in the gas phase)
          !
          scr_te1 = ( lheat )*mflux(i,j,k)*sqrt(dvofdx**2+dvofdy**2+dvofdz**2)*rhocpci
          scr_te2 = dpthdto*(1.d0-vof(i,j,k))*rhocpci
          scr_te3 = d_lg(i,j,k)*((dscadxp+dscadxm)*(dtmpdxp+dtmpdxm) + &
                                 (dscadyp+dscadym)*(dtmpdyp+dtmpdym) + &
                                 (dscadzp+dscadzm)*(dtmpdzp+dtmpdzm))*0.25d0*(1.d0-vof(i,j,k))*delta_cp*cpci
          !
          tmp(i,j,k)         = tmp(i,j,k) + factor1*dtmpdtrk(i,j,k) + factor2*dtmpdtrkold(i,j,k) + &
                               factor12*( - scr_te1 + scr_te2 + scr_te3 )
          dtmpdtrkold(i,j,k) = dtmpdtrk(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine ab2_tmp
  !
end module mod_rks
