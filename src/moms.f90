module mod_moms
  !
  use mod_gradls    , only: weno5
  use mod_common_mpi, only: coord,comm_cart,ierr
  use mod_param     , only: ru,cp2,m1,m2,rho1,cp1,delta_cp
  use mod_thermo    , only: thermo_rhog,thermo_d_lg,mass_fraction
  use mpi
  !
  implicit none
  !
  real(8), parameter :: theta_thr = 0.25d0
  !
  private
  public  :: momtad,rhs_pth,momsad_flm
  !
  contains
  !
  subroutine momsad_flm(n,dli,q,pth,rhog,d_lg,u,v,w,s,phi,tmpge,tmple,dsdt)
    !
    ! first-order scheme to compute the gradient
    ! on an irregular domain (see Sato and Niceno, JCP2013 for the reference)
    ! (modified to be more general as it checks for each i,j,k both the 
    ! (im,ip; jm,jp; km,kp) cells so that the upwind/central scheme is always 
    ! possible to compute)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in )                         :: pth
    integer, intent(in )                         :: q
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rhog,d_lg
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: s
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: phi
    real(8), intent(in ), dimension(-q:,-q:,-q:) :: tmpge,tmple
    real(8), intent(out), dimension(  :,  :,  :) :: dsdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(8) :: rdl_xp ,rdl_xm ,rdl_yp ,rdl_ym ,rdl_zp ,rdl_zm
    real(8) :: udsdx,vdsdy,wdsdz,uc,vc,wc
    real(8) :: xx,yy,zz!,dlp,dlm
    real(8) :: theta,tmpl_ip,tmpl_im,tmpg_ip,tmpg_im,ss_ip,ss_im, &
               dlg_im,dlg_ip,rhog_im,rhog_ip
    real(8), dimension(3) :: dl
    !
    dl = dli**(-1.d0)
    !
    do k=1,n(3)
      kp = k+1
      km = k-1
      zz = (k-0.5d0)*dl(3)
      do j=1,n(2)
        jp = j+1
        jm = j-1
        yy = (j+coord(2)*n(2)-0.5d0)*dl(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          xx = (i+coord(1)*n(1)-0.5d0)*dl(1)
          !
          if(phi(i,j,k).gt.0.d0) then
            !
            dsdt(i,j,k) = 0.d0
            !
          else
            !
            ! Advection contribution
            !
            uc = 0.5d0*( u(im,j,k)+u(i,j,k) )
            vc = 0.5d0*( v(i,jm,k)+v(i,j,k) )
            wc = 0.5d0*( w(i,j,km)+w(i,j,k) )
            !
            dsdxp = (s(ip,j,k)-s(i ,j,k))*dli(1)
            dsdxm = (s(i ,j,k)-s(im,j,k))*dli(1)
            dsdyp = (s(i,jp,k)-s(i,j ,k))*dli(2)
            dsdym = (s(i,j ,k)-s(i,jm,k))*dli(2)
            dsdzp = (s(i,j,kp)-s(i,j,k ))*dli(3)
            dsdzm = (s(i,j,k )-s(i,j,km))*dli(3)
            !
            ! Diffusion contribution
            !
            rdl_xp = 0.50d0*(rhog(ip,j,k)*d_lg(ip,j,k)+rhog(i ,j,k)*d_lg(i ,j,k))
            rdl_xm = 0.50d0*(rhog(i ,j,k)*d_lg(i ,j,k)+rhog(im,j,k)*d_lg(im,j,k))
            rdl_yp = 0.50d0*(rhog(i,jp,k)*d_lg(i,jp,k)+rhog(i,j ,k)*d_lg(i,j ,k))
            rdl_ym = 0.50d0*(rhog(i,j ,k)*d_lg(i,j ,k)+rhog(i,jm,k)*d_lg(i,jm,k))
            rdl_zp = 0.50d0*(rhog(i,j,kp)*d_lg(i,j,kp)+rhog(i,j,k )*d_lg(i,j,k ))
            rdl_zm = 0.50d0*(rhog(i,j,k )*d_lg(i,j,k )+rhog(i,j,km)*d_lg(i,j,km))
            !
            ! now we change the stencil discretization in those cells
            ! cut by the inteface, following a dimension by dimension approach
            ! 
            ! along x
            !
            if(     phi(ip,j,k)*phi(i,j,k).lt.0.d0.and.phi(im,j,k)*phi(i,j,k).gt.0.d0 ) then ! im,i same side, ip not
              theta   = abs(phi(i,j,k))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmpl_ip = (tmple(i,j,k)*abs(phi(ip,j,k))+tmple(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmpg_ip = (tmpge(i,j,k)*abs(phi(ip,j,k))+tmpge(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              rhog_ip = thermo_rhog(pth,tmpg_ip,ss_ip)
              dlg_ip  = thermo_d_lg(pth,tmpg_ip,ss_ip)
              !rdl_xp  = 0.50d0*(rhog_ip*dlg_ip+rhog(i ,j,k)*d_lg(i ,j,k))
              if(theta.ge.theta_thr) then
                rdl_xp = (1.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(i ,j,k)*d_lg(i ,j,k))/theta
                rdl_xp = 0.5d0*(rdl_xp + rhog(i ,j,k)*d_lg(i ,j,k))
                dsdxp  = - ( s(i ,j,k)-ss_ip )*dli(1)/theta
              else
                rdl_xp = (2.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(im,j,k)*d_lg(im,j,k))/(1.d0+theta)
                rdl_xp = 0.5d0*(rdl_xp + rhog(i ,j,k)*d_lg(i ,j,k))
                dsdxp  = - ( s(im,j,k)-ss_ip )*dli(1)/(1.d0+theta)
              endif
            elseif( phi(im,j,k)*phi(i,j,k).lt.0.d0.and.phi(ip,j,k)*phi(i,j,k).gt.0.d0 ) then ! ip,i same side, im not
              theta   = abs(phi(i,j,k))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              tmpl_im = (tmple(i,j,k)*abs(phi(im,j,k))+tmple(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              tmpg_im = (tmpge(i,j,k)*abs(phi(im,j,k))+tmpge(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              rhog_im = thermo_rhog(pth,tmpg_im,ss_im)
              dlg_im  = thermo_d_lg(pth,tmpg_im,ss_im)
              !rdl_xm  = 0.50d0*(rhog_im*dlg_im+rhog(i ,j,k)*d_lg(i ,j,k))
              if(theta.ge.theta_thr) then
                rdl_xm = (1.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(i ,j,k)*d_lg(i ,j,k))/theta
                rdl_xm = 0.5d0*(rdl_xm + rhog(i ,j,k)*d_lg(i ,j,k))
                dsdxm  = + ( s(i ,j,k)-ss_im )*dli(1)/theta
              else
                rdl_xm = (2.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(ip,j,k)*d_lg(ip,j,k))/(1.d0+theta)
                rdl_xm = 0.5d0*(rdl_xm + rhog(i ,j,k)*d_lg(i ,j,k))
                dsdxm  = + ( s(ip,j,k)-ss_im )*dli(1)/(1.d0+theta)
              endif
            elseif( phi(im,j,k)*phi(i,j,k).lt.0.d0.and.phi(ip,j,k)*phi(i,j,k).lt.0.d0 ) then ! ip,im same side, i not
              rdl_xp = 0.d0
              rdl_xm = 0.d0
              dsdxp  = 0.d0
              dsdxm  = 0.d0
            endif
            !
            ! along y
            !
            if(     phi(i,jp,k)*phi(i,j,k).lt.0.d0.and.phi(i,jm,k)*phi(i,j,k).gt.0.d0 ) then ! jm,j same side, jp not
              theta   = abs(phi(i,j,k))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmpl_ip = (tmple(i,j,k)*abs(phi(i,jp,k))+tmple(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmpg_ip = (tmpge(i,j,k)*abs(phi(i,jp,k))+tmpge(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              rhog_ip = thermo_rhog(pth,tmpg_ip,ss_ip)
              dlg_ip  = thermo_d_lg(pth,tmpg_ip,ss_ip)
              !rdl_yp  = 0.50d0*(rhog_ip*dlg_ip+rhog(i,j ,k)*d_lg(i,j ,k))
              if(theta.ge.theta_thr) then
                rdl_yp = (1.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(i,j ,k)*d_lg(i,j ,k))/theta
                rdl_yp = 0.5d0*(rdl_yp + rhog(i,j ,k)*d_lg(i,j ,k))
                dsdyp  = - ( s(i,j ,k)-ss_ip )*dli(2)/theta
              else
                rdl_yp = (2.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(i,jm,k)*d_lg(i,jm,k))/(1.d0+theta)
                rdl_yp = 0.5d0*(rdl_yp + rhog(i,j ,k)*d_lg(i,j ,k))
                dsdyp  = - ( s(i,jm,k)-ss_ip )*dli(2)/(1.d0+theta)
              endif
            elseif( phi(i,jm,k)*phi(i,j,k).lt.0.d0.and.phi(i,jp,k)*phi(i,j,k).gt.0.d0 ) then ! jp,j same side, jm not
              theta   = abs(phi(i,j,k))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              tmpl_im = (tmple(i,j,k)*abs(phi(i,jm,k))+tmple(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              tmpg_im = (tmpge(i,j,k)*abs(phi(i,jm,k))+tmpge(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              rhog_im = thermo_rhog(pth,tmpg_im,ss_im)
              dlg_im  = thermo_d_lg(pth,tmpg_im,ss_im)
              !rdl_ym  = 0.50d0*(rhog_im*dlg_im+rhog(i,j ,k)*d_lg(i,j ,k))
              if(theta.ge.theta_thr) then
                rdl_ym = (1.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(i,j ,k)*d_lg(i,j ,k))/theta
                rdl_ym = 0.5d0*(rdl_ym + rhog(i,j ,k)*d_lg(i,j ,k))
                dsdym  = + ( s(i,j ,k)-ss_im )*dli(2)/theta
              else
                rdl_ym = (2.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(i,jp,k)*d_lg(i,jp,k))/(1.d0+theta)
                rdl_ym = 0.5d0*(rdl_ym + rhog(i,j ,k)*d_lg(i,j ,k))
                dsdym  = + ( s(i,jp,k)-ss_im )*dli(2)/(1.d0+theta)
              endif
            elseif( phi(i,jp,k)*phi(i,j,k).lt.0.d0.and.phi(i,jp,k)*phi(i,j,k).lt.0.d0 ) then ! jp,jm same side, j not
              rdl_yp = 0.d0
              rdl_ym = 0.d0
              dsdyp  = 0.d0
              dsdym  = 0.d0
            endif
            !
            ! along z
            !
            if(     phi(i,j,kp)*phi(i,j,k).lt.0.d0.and.phi(i,j,km)*phi(i,j,k).gt.0.d0 ) then ! km,k same side, kp not
              theta   = abs(phi(i,j,k))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmpl_ip = (tmple(i,j,k)*abs(phi(i,j,kp))+tmple(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmpg_ip = (tmpge(i,j,k)*abs(phi(i,j,kp))+tmpge(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              ss_ip   = mass_fraction(pth,tmpl_ip)
              rhog_ip = thermo_rhog(pth,tmpg_ip,ss_ip)
              dlg_ip  = thermo_d_lg(pth,tmpg_ip,ss_ip)
              !rdl_zp  = 0.50d0*(rhog_ip*dlg_ip+rhog(i,j,k )*d_lg(i,j,k ))
              if(theta.gt.theta_thr) then
                rdl_zp = (1.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(i,j,k )*d_lg(i,j,k ))/theta
                rdl_zp = 0.5d0*(rdl_zp + rhog(i,j,k )*d_lg(i,j,k ))
                dsdzp  = - ( s(i,j,k )-ss_ip )*dli(3)/theta
              else
                rdl_zp = (2.d0*rhog_ip*dlg_ip+(theta-1.d0)*rhog(i,j,km)*d_lg(i,j,km))/(1.d0+theta)
                rdl_zp = 0.5d0*(rdl_zp + rhog(i,j,k )*d_lg(i,j,k ))
                dsdzp  = - ( s(i,j,km)-ss_ip )*dli(3)/(1.d0+theta)
              endif
            elseif( phi(i,j,km)*phi(i,j,k).lt.0.d0.and.phi(i,j,kp)*phi(i,j,k).gt.0.d0 ) then ! kp,k same side, km not
              theta   = abs(phi(i,j,k))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              tmpl_im = (tmple(i,j,k)*abs(phi(i,j,km))+tmple(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              tmpg_im = (tmpge(i,j,k)*abs(phi(i,j,km))+tmpge(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              ss_im   = mass_fraction(pth,tmpl_im)
              rhog_im = thermo_rhog(pth,tmpg_im,ss_im)
              dlg_im  = thermo_d_lg(pth,tmpg_im,ss_im)
              !rdl_zm  = 0.50d0*(rhog_im*dlg_im+rhog(i,j,k )*d_lg(i,j,k ))
              if(theta.gt.theta_thr) then
                rdl_zm = (1.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(i,j,k )*d_lg(i,j,k ))/theta
                rdl_zm = 0.5d0*(rdl_zm + rhog(i,j,k )*d_lg(i,j,k ))
                dsdzm  = + ( s(i,j,k )-ss_im )*dli(3)/theta
              else
                rdl_zm = (2.d0*rhog_im*dlg_im+(theta-1.d0)*rhog(i,j,kp)*d_lg(i,j,kp))/(1.d0+theta)
                rdl_zm = 0.5d0*(rdl_zm + rhog(i,j,k )*d_lg(i,j,k ))
                dsdzm  = + ( s(i,j,kp)-ss_im )*dli(3)/(1.d0+theta)
              endif
            elseif( phi(i,j,km)*phi(i,j,k).lt.0.d0.and.phi(i,j,kp)*phi(i,j,k).lt.0.d0 ) then ! kp,km same side, k not
              rdl_zp = 0.d0
              rdl_zm = 0.d0
              dsdzp  = 0.d0
              dsdzm  = 0.d0
            endif
            !
#ifndef TWOD
            udsdx = 0.5d0*(uc+abs(uc))*dsdxm + 0.5d0*(uc-abs(uc))*dsdxp
#endif
            vdsdy = 0.5d0*(vc+abs(vc))*dsdym + 0.5d0*(vc-abs(vc))*dsdyp
            wdsdz = 0.5d0*(wc+abs(wc))*dsdzm + 0.5d0*(wc-abs(wc))*dsdzp
            !
            dsdt(i,j,k) = - ( udsdx + vdsdy + wdsdz ) + &
                            ( &
#ifndef TWOD
                            (rdl_xp*dsdxp-rdl_xm*dsdxm)*dli(1) + &
#endif
                            (rdl_yp*dsdyp-rdl_ym*dsdym)*dli(2) + &
                            (rdl_zp*dsdzp-rdl_zm*dsdzm)*dli(3) )/rhog(i,j,k)
            ! 
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine momsad_flm
  !
  subroutine momtad(n,dli,u,v,w,rhog,vof,kappa,s,dsdt)
    !
    ! energy equation solved using the whole domain formulation
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rhog,vof,kappa
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: s
    real(8), intent(out), dimension(  :,  :,  :) :: dsdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    real(8) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    real(8) :: rhocpci
    integer :: qmin 
    !
    qmin = abs(lbound(u,1))
    !
    call weno5(n,dli,qmin,.true.,s,u,v,w,dsdt)
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
          kappaxp = 0.5d0*(kappa(ip,j,k)+kappa(i,j,k))
          kappaxm = 0.5d0*(kappa(im,j,k)+kappa(i,j,k))
          kappayp = 0.5d0*(kappa(i,jp,k)+kappa(i,j,k))
          kappaym = 0.5d0*(kappa(i,jm,k)+kappa(i,j,k))
          kappazp = 0.5d0*(kappa(i,j,kp)+kappa(i,j,k))
          kappazm = 0.5d0*(kappa(i,j,km)+kappa(i,j,k))
          !
          !rhocpci = vof(i,j,k)/(rho1*cp1) + ( 1.d0-vof(i,j,k) )/(rhog(i,j,k)*cp2)
          rhocpci = (vof(i,j,k)*(rho1*cp1) + ( 1.d0-vof(i,j,k) )*(rhog(i,j,k)*cp2))**(-1.d0)
          !
          dsdxp  = (s(ip,j,k)-s(i ,j,k))*dli(1)
          dsdxm  = (s(i ,j,k)-s(im,j,k))*dli(2)
          dsdyp  = (s(i,jp,k)-s(i,j ,k))*dli(2)
          dsdym  = (s(i,j ,k)-s(i,jm,k))*dli(2)
          dsdzp  = (s(i,j,kp)-s(i,j,k ))*dli(3)
          dsdzm  = (s(i,j,k )-s(i,j,km))*dli(3)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + &
                        rhocpci*( &
#ifndef TWOD
                                  (kappaxp*dsdxp-kappaxm*dsdxm)*dli(1) + &
#endif
                                  (kappayp*dsdyp-kappaym*dsdym)*dli(2) + &
                                  (kappazp*dsdzp-kappazm*dsdzm)*dli(3)   &
                                ) 
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine momtad
  !
  subroutine rhs_pth(n,dli,vof,mflux,kappa,rho_gas,d_lg,rho1,tmp,tmpge,sca,scae,pth,dpthdt_n,mfxt,integral)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: vof,mflux,kappa
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: rho_gas,d_lg
    real(8), intent(in )                         :: rho1
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: tmpge
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: sca,scae
    real(8), intent(in )                         :: pth
    real(8), intent(out)                         :: dpthdt_n,mfxt,integral
    !
    real(8), dimension(8) :: mx,my,mz
    real(8), dimension(3) :: dl
    integer :: i,j,k,im,jm,km,ip,jp,kp
    real(8) :: dvofdx,dvofdy,dvofdz
    real(8) :: diff_t,diff_s
    real(8) :: s_in,t_in,rhog_in,rhol_in
    real(8) :: m_avg,term_den,term_pch
    real(8) :: dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm
    real(8) :: dscadxp,dscadxm,dscadyp,dscadym,dscadzp,dscadzm
    real(8) :: rhdlgxp,rhdlgxm,rhdlgyp,rhdlgym,rhdlgzp,rhdlgzm
    real(8) :: vofxp,vofxm,vofyp,vofym,vofzp,vofzm
    real(8) :: kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm
    !
    dl(:) = dli(:)**(-1.d0)
    !
    dpthdt_n = 0.0d0
    term_den = 0.0d0
    mfxt     = 0.0d0
    integral = 0.0d0
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
          ! a. calculation of the interfacial area
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
          dvofdx = 0.125d0*(mx(1)+mx(2)+mx(3)+mx(4)+mx(5)+mx(6)+mx(7)+mx(8))
          dvofdy = 0.125d0*(my(1)+my(2)+my(3)+my(4)+my(5)+my(6)+my(7)+my(8))
          dvofdz = 0.125d0*(mz(1)+mz(2)+mz(3)+mz(4)+mz(5)+mz(6)+mz(7)+mz(8))
          !
          ! b. vapour mass effects
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
          m_avg   = (sca(i,j,k)/m1+(1.d0-sca(i,j,k))/m2)**(-1.d0)
          diff_s  = ((dscadxp*rhdlgxp-dscadxm*rhdlgxm)*dli(1) + &
                     (dscadyp*rhdlgyp-dscadym*rhdlgym)*dli(2) + &
                     (dscadzp*rhdlgzp-dscadzm*rhdlgzm)*dli(3))*(1.d0/rho_gas(i,j,k))*m_avg*(1.d0/m1-1.d0/m2)
          !
          ! c. thermal effects
          ! 
          dtmpdxm = (tmp(i ,j,k)-tmp(im,j,k))*dli(1)
          dtmpdxp = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
          dtmpdym = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
          dtmpdyp = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
          dtmpdzm = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
          dtmpdzp = (tmp(i,j,kp)-tmp(i,j ,k))*dli(3)
          !
          kappaxm = 0.5d0*(kappa(i,j,k)+kappa(im,j,k))
          kappaxp = 0.5d0*(kappa(i,j,k)+kappa(ip,j,k))
          kappaym = 0.5d0*(kappa(i,j,k)+kappa(i,jm,k))
          kappayp = 0.5d0*(kappa(i,j,k)+kappa(i,jp,k))
          kappazm = 0.5d0*(kappa(i,j,k)+kappa(i,j,km))
          kappazp = 0.5d0*(kappa(i,j,k)+kappa(i,j,kp))
          ! 
          diff_t  = ( (kappaxp*dtmpdxp-kappaxm*dtmpdxm)*dli(1) + &
                      (kappayp*dtmpdyp-kappaym*dtmpdym)*dli(2) + &
                      (kappazp*dtmpdzp-kappazm*dtmpdzm)*dli(3) )*(ru/(cp2*m_avg)) + &
                    delta_cp*d_lg(i,j,k)*( &
                             (dscadxp+dscadxm)*(dtmpdxp+dtmpdxm) + &
                             (dscadyp+dscadym)*(dtmpdyp+dtmpdym) + &
                             (dscadzp+dscadzm)*(dtmpdzp+dtmpdzm)   &
                                         )*0.25d0*rho_gas(i,j,k)*(ru/(cp2*m_avg)) 
          ! 
          ! d. compute an estimation of the t_in and s_in at the interface, but distributed in a band
          !
          s_in    = scae(i,j,k)
          t_in    = tmpge(i,j,k) 
          rhog_in = thermo_rhog(pth,t_in,s_in) ! we should provide an estimation of the interfacial value
          rhol_in = rho1
          !
          term_pch = mflux(i,j,k)*(1.d0/rhog_in-1.d0/rhol_in)*sqrt(dvofdx**2.0+dvofdy**2.0+dvofdz**2.0)
          !
          ! e. compute the denominator
          !
          term_den = term_den + (1.d0-vof(i,j,k))*(1.d0-(ru/(cp2*m_avg)))
          !
          ! f. put all together
          !
          dpthdt_n = dpthdt_n + term_pch + (pth*diff_s + diff_t)*(1.d0-vof(i,j,k))*(1.d0/pth)
          !
          ! h. compute two quantities for pth
          !
          mfxt = mfxt + mflux(i,j,k)*sqrt(dvofdx**2+dvofdy**2+dvofdz**2)*dl(1)*dl(2)*dl(3)
          !
          integral = integral + (1.d0-vof(i,j,k))*(m_avg/ru)*dl(1)*dl(2)*dl(3)/tmp(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    call mpi_allreduce(MPI_IN_PLACE,term_den,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,dpthdt_n,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mfxt    ,1,mpi_real8,mpi_sum,comm_cart,ierr)
    call mpi_allreduce(MPI_IN_PLACE,integral,1,mpi_real8,mpi_sum,comm_cart,ierr)
    !
    dpthdt_n = dpthdt_n/term_den
    !
    return
  end subroutine rhs_pth
  !
end module mod_moms
