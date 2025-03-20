module mod_mom
  !
  use mod_common_mpi, only: coord
  use mod_funcs     , only: interp_g
  use mod_param     , only: sigmaca,gacc,pi
#ifndef LOW_MACH
  use mod_param     , only: beta_g_th,rho2_0,rho1,beta_l_th,tmp0
#endif 
  !
  implicit none
  !
  public  :: momxad,momyad,momzad,momxp,momyp,momzp,turb_forc_src
  private
  !
  contains
  !
  subroutine momxad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,u,v,w,mu,rho,dudt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in )                         :: dxi,dyi,dzi
    real(8), intent(in ), dimension(-2)          :: dzci,dzfi
    !real(8), intent(in ), dimension( 0:, 0:, 0:) :: u,v,w
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: mu,rho
    !real(8), intent(out), dimension( 0:, 0:, 0:) :: dudt
    real(8), intent(out), dimension(  :,  :,  :) :: dudt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k,q
    real(8) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(8) :: ug_u,vg_u,wg_u,uc,vc,wc
    real(8) :: dudxp,dudxm,dvdxp,dvdxm,dudyp,dudym, &
               dudzp,dudzm,dwdxp,dwdxm
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divxm,divxp
    real(8) :: rhox
    real(8), dimension(-2:2) :: vec
    !
    !$omp parallel default(shared) &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uuip,uuim,uvjp,uvjm,uwkp,uwkm,dudxp,dudxm,dudyp,dudym,dudzp,dudzm)
    !$omp do
    do k=1,nz
      kp = k+1
      km = k-1
      do j=1,ny
        jp = j+1
        jm = j-1
        do i=1,nx
          ip = i+1
          im = i-1
          ! 
          ! Advection
          ! 
          !uuip  = 0.25d0*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          !uuim  = 0.25d0*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          !uvjp  = 0.25d0*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          !uvjm  = 0.25d0*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          !uwkp  = 0.25d0*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          !uwkm  = 0.25d0*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          !
          uc = u(i,j,k)
          do q=-2,2
            vec(q) = u(i+q,j,k)
          enddo
          ug_u = interp_g(vec,uc,dxi)
          !
          vc = 0.25d0*( v(i,j,k)+v(ip,j,k)+v(i,jm,k)+v(ip,jm,k) )
          do q=-2,2
            vec(q) = u(i,j+q,k)
          enddo
          vg_u = interp_g(vec,vc,dyi)
          !
          wc = 0.25d0*( w(i,j,k)+w(ip,j,k)+w(i,j,km)+w(ip,j,km) )
          do q=-2,2
            vec(q) = u(i,j,k+q)
          enddo
          wg_u = interp_g(vec,wc,dzi)
          !
          ! Diffusion
          ! 
          dudxp = (u(ip,j ,k)-u(i ,j ,k))*dxi
          dudxm = (u(i ,j ,k)-u(im,j ,k))*dxi
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(ip,jm,k)-v(i ,jm,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(i ,j ,k)-u(i ,jm,k))*dyi
          dudzp = (u(i ,j,kp)-u(i ,j,k ))*dzi
          dudzm = (u(i ,j,k )-u(i ,j,km))*dzi
          dwdxp = (w(ip,j,k )-w(i ,j,k ))*dxi
          dwdxm = (w(ip,j,km)-w(i ,j,km))*dxi
          !
#ifdef LOW_MACH
          divxp = (w(ip,j,k)-w(ip,j ,km))*dzi + &
                  (v(ip,j,k)-v(ip,jm,k ))*dyi + &
                  (u(ip,j,k)-u(i ,j ,k ))*dxi
          !
          divxm = (w(i ,j,k)-w(i ,j ,km))*dzi + &
                  (v(i ,j,k)-v(i ,jm,k ))*dyi + &
                  (u(i ,j,k)-u(im,j ,k ))*dxi
#else
          divxp = 0.d0
          divxm = 0.d0
#endif
          !
          muxp = mu(ip,j,k)
          muxm = mu(i ,j,k)
          muyp = 0.25d0*(mu(i,j,k)+mu(i,jp,k)+mu(ip,jp,k)+mu(ip,j,k))
          muym = 0.25d0*(mu(i,j,k)+mu(i,jm,k)+mu(ip,jm,k)+mu(ip,j,k))
          muzp = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j,kp)+mu(ip,j,k))
          muzm = 0.25d0*(mu(i,j,k)+mu(i,j,km)+mu(ip,j,km)+mu(ip,j,k))
          !
          ! Momentum balance
          !
          rhox = 0.5d0*(rho(  ip,j,k)+rho(  i,j,k))
          !
          dudt(i,j,k) = &
                        !
                        !dxi*( -uuip + uuim ) + &
                        !dyi*( -uvjp + uvjm ) + &
                        !dzi*( -uwkp + uwkm ) + &
                        !u(i,j,k)*(divxp+divxm) + &
                        !
                        - (ug_u + vg_u + wg_u) + &
                        !
                        ( &
                        dxi*((dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm - (2.d0/3.d0)*(muxp*divxp-muxm*divxm)) + &
                        dyi*((dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym                                      ) + &
                        dzi*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm                                      )   &
                        )/rhox
          !
        enddo
      enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momxad
  !
  subroutine momyad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,u,v,w,mu,rho,dvdt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in )                         :: dxi,dyi,dzi
    real(8), intent(in ), dimension(-2)          :: dzci,dzfi
    !real(8), intent(in ), dimension( 0:, 0:, 0:) :: u,v,w
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: mu,rho
    !real(8), intent(out), dimension( 0:, 0:, 0:) :: dvdt
    real(8), intent(out), dimension(  :,  :,  :) :: dvdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k,q
    real(8) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(8) :: ug_v,vg_v,wg_v,uc,vc,wc
    real(8) :: dvdxp,dvdxm,dudyp,dudym,dvdyp,dvdym, &
               dvdzp,dvdzm,dwdyp,dwdym
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divym,divyp
    real(8) :: rhoy
    real(8), dimension(-2:2) :: vec
    !
    !$omp parallel default(shared) &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uvip,uvim,vvjp,vvjm,wvkp,wvkm,dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm)
    !$omp do
    do k=1,nz
      kp = k+1
      km = k-1
      do j=1,ny
        jp = j+1
        jm = j-1
        do i=1,nx
          ip = i+1
          im = i-1
          !
          ! Advection
          !
          !uvip  = 0.25d0*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          !uvim  = 0.25d0*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          !vvjp  = 0.25d0*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          !vvjm  = 0.25d0*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          !wvkp  = 0.25d0*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          !wvkm  = 0.25d0*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          !
          uc = 0.25d0*( u(i,j,k)+u(i,jp,k)+u(im,j,k)+u(im,jp,k) )
          do q=-2,2
            vec(q) = v(i+q,j,k)
          enddo
          ug_v = interp_g(vec,uc,dxi)
          !
          vc = v(i,j,k) 
          do q=-2,2
            vec(q) = v(i,j+q,k)
          enddo
          vg_v = interp_g(vec,vc,dyi)
          !
          wc = 0.25d0*( w(i,j,k)+w(i,jp,k)+w(i,j,km)+w(i,jp,km) )
          do q=-2,2
            vec(q) = v(i,j,k+q)
          enddo
          wg_v = interp_g(vec,wc,dzi)
          !
          ! Diffusion
          !
          dvdxp = (v(ip,j ,k)-v(i ,j ,k))*dxi
          dvdxm = (v(i ,j ,k)-v(im,j ,k))*dxi
          dudyp = (u(i ,jp,k)-u(i ,j ,k))*dyi
          dudym = (u(im,jp,k)-u(im,j ,k))*dyi
          dvdyp = (v(i,jp ,k)-v(i ,j ,k))*dyi
          dvdym = (v(i,j  ,k)-v(i ,jm,k))*dyi
          dvdzp = (v(i,j,kp) -v(i,j,k  ))*dzi
          dvdzm = (v(i,j,k ) -v(i,j,km ))*dzi
          dwdyp = (w(i,jp,k )-w(i,j,k  ))*dyi
          dwdym = (w(i,jp,km)-w(i,j,km ))*dyi
          !
#ifdef LOW_MACH
          divyp = (w(i,jp,k)-w(i ,jp,km))*dzi + &
                  (v(i,jp,k)-v(i ,j ,k ))*dyi + &
                  (u(i,jp,k)-u(im,jp,k ))*dxi
          !
          divym = (w(i,j ,k)-w(i ,j ,km))*dzi + &
                  (v(i,j ,k)-v(i ,jm,k ))*dyi + &
                  (u(i,j ,k)-u(im,j ,k ))*dxi
#else
          divyp = 0.d0
          divym = 0.d0
#endif
          !
          muxp = 0.25d0*(mu(i,j,k)+mu(ip,j,k)+mu(ip,jp,k)+mu(i,jp,k))
          muxm = 0.25d0*(mu(i,j,k)+mu(im,j,k)+mu(im,jp,k)+mu(i,jp,k))
          muyp = mu(i,jp,k)
          muym = mu(i,j ,k)
          muzp = 0.25d0*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,kp)+mu(i,j,kp))
          muzm = 0.25d0*(mu(i,j,k)+mu(i,jp,k)+mu(i,jp,km)+mu(i,j,km))
          !
          ! Momentum balance
          !
          rhoy = 0.5d0*(rho(  i,jp,k)+rho  (i,j,k))
          !
          dvdt(i,j,k) = & 
                        !
                        !dxi*( -uvip + uvim ) + &
                        !dyi*( -vvjp + vvjm ) + & 
                        !dzi*( -wvkp + wvkm ) + &
                        !v(i,j,k)*(divyp+divym) + &
                        !
                        - ( ug_v + vg_v + wg_v) + &
                        !
                        ( &
                        dxi*((dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm                                      ) + &
                        dyi*((dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym - (2.d0/3.d0)*(muyp*divyp-muym*divym)) + &
                        dzi*((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm                                      )   &
                        )/rhoy
          !
        enddo
      enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momyad
  !
  subroutine momzad(nx,ny,nz,dxi,dyi,dzi,dzci,dzfi,u,v,w,mu,rho,dwdt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in )                         :: dxi,dyi,dzi
    real(8), intent(in ), dimension(-2)          :: dzci,dzfi
    !real(8), intent(in ), dimension( 0:, 0:, 0:) :: u,v,w
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: u,v,w
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: mu,rho
    !real(8), intent(out), dimension( 0:, 0:, 0:) :: dwdt
    real(8), intent(out), dimension(  :,  :,  :) :: dwdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k,q
    real(8) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(8) :: dwdxp,dwdxm,dudzp,dudzm,dwdyp,dwdym, &
               dvdzp,dvdzm,dwdzp,dwdzm
    real(8) :: ug_w,vg_w,wg_w,uc,vc,wc
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divzm,divzp
    real(8) :: rhoz
    real(8), dimension(-2:2) :: vec
    !
    !$omp parallel default(shared) &
    !$omp&private(i,j,k,ip,jp,kp,im,jm,km) &
    !$omp&private(uwip,uwim,vwjp,vwjm,wwkp,wwkm,dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm)
    !$omp do
    do k=1,nz
      kp = k+1
      km = k-1
      do j=1,ny
        jp = j+1
        jm = j-1
        do i=1,nx
          ip = i+1
          im = i-1
          !
          ! Advection
          ! 
          !uwip  = 0.25d0*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          !uwim  = 0.25d0*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          !vwjp  = 0.25d0*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          !vwjm  = 0.25d0*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          !wwkp  = 0.25d0*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          !wwkm  = 0.25d0*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          !   
          uc = 0.25d0*( u(i,j,k)+u(i,j,kp)+u(im,j,k)+u(im,j,kp) )
          do q=-2,2
            vec(q) = w(i+q,j,k)
          enddo
          ug_w = interp_g(vec,uc,dxi)
          !
          vc = 0.25d0*( v(i,j,k)+v(i,j,kp)+v(i,jm,k)+v(i,jm,kp) )
          do q=-2,2
            vec(q) = w(i,j+q,k)
          enddo
          vg_w = interp_g(vec,vc,dyi)
          !
          wc = w(i,j,k)
          do q=-2,2
            vec(q) = w(i,j,k+q)
          enddo
          wg_w = interp_g(vec,wc,dzi)
          !
          ! Diffusion
          ! 
          dwdxp = (w(ip,j,k )-w(i ,j,k))*dxi
          dwdxm = (w(i ,j,k )-w(im,j,k))*dxi
          dudzp = (u(i ,j,kp)-u(i ,j,k))*dzi
          dudzm = (u(im,j,kp)-u(im,j,k))*dzi
          dwdyp = (w(i,jp,k )-w(i,j ,k))*dyi
          dwdym = (w(i,j ,k )-w(i,jm,k))*dyi
          dvdzp = (v(i,j ,kp)-v(i,j ,k))*dzi
          dvdzm = (v(i,jm,kp)-v(i,jm,k))*dzi
          dwdzp = (w(i,j,kp )-w(i,j,k ))*dzi
          dwdzm = (w(i,j,k  )-w(i,j,km))*dzi
          !
#ifdef LOW_MACH
          divzp = (w(i,j,kp)-w(i,j ,k ))*dzi + &
                  (v(i,j,kp)-v(i,jm,kp))*dyi + &
                  (u(i,j,kp)-u(im,j,kp))*dxi
          
          divzm = (w(i,j,k)-w(i ,j ,km))*dzi + &
                  (v(i,j,k)-v(i ,jm,k ))*dyi + &
                  (u(i,j,k)-u(im,j ,k ))*dxi
#else
          divzp = 0.d0
          divzm = 0.d0
#endif
          !
          muxp = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(ip,j ,kp)+mu(ip,j ,k))
          muxm = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(im,j ,kp)+mu(im,j ,k))
          muyp = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jp,kp)+mu(i ,jp,k))
          muym = 0.25d0*(mu(i,j,k)+mu(i,j,kp)+mu(i ,jm,kp)+mu(i ,jm,k))
          muzp = mu(i,j,kp)
          muzm = mu(i,j,k )
          !
          ! Momentum balance
          !
          rhoz = 0.5d0*(rho(  i,j,kp)+rho(  i,j,k))
          !
          dwdt(i,j,k) = &
                        !
                        !dxi*( -uwip + uwim ) + &
                        !dyi*( -vwjp + vwjm ) + &
                        !dzi*( -wwkp + wwkm ) + &
                        !w(i,j,k)*(divzp+divzm) + &
                        !
                        - ( ug_w + vg_w + wg_w) + &
                        !
                        ( &
                        dxi*((dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm                                      ) + &
                        dyi*((dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym                                      ) + &
                        dzi*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm - (2.d0/3.d0)*(muzp*divzp-muzm*divzm)) & 
                        )/rhoz
          !
        enddo
      enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momzad
  !
  subroutine momxp(nx,ny,nz,dxi,p,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dudt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in )                         :: dxi 
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: p,pp,curv,vof,rho,rho_gas
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in )                         :: rho0,rho_av 
    real(8), intent(out), dimension(  :,  :,  :) :: dudt
    !
    real(8) :: rhox,rho0i,curvx
#ifndef LOW_MACH
    real(8) :: rhox_ob,voff,tmpf,rhogf
#endif
    integer :: i,j,k
    integer :: ip
    !
    rho0i = 1.d0/rho0
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,ip) &
    !$OMP PRIVATE(rhox,curvx) &
    !$OMP SHARED(nx,ny,nz,dxi,p,rho,curv,dudt)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          ip = i + 1
          !
          rhox    = 0.5d0*(rho( ip,j,k)+rho( i,j,k))
          curvx   = 0.5d0*(curv(ip,j,k)+curv(i,j,k))
#ifdef LOW_MACH
          !rhogf   = 0.5d0*(rho_gas(ip,j,k)+rho_gas(i,j,k))
#else
          voff    = 0.5d0*(vof( ip,j,k)+vof( i,j,k))
          tmpf    = 0.5d0*(tmp( ip,j,k)+tmp( i,j,k))
          rhogf   = rho2_0*(1.d0-beta_g_th*(tmpf-tmp0))
#endif
          !rhox_ob = voff*rho1*(1.d0-beta_l_th*(tmpf-tmp0)) + (1.d0-voff)*rhogf
          !
          dudt(i,j,k) = - ( p(ip,j,k)-p(i,j,k) )*dxi !- dpdl(1)
          dudt(i,j,k) = ( gacc(1)*(rhox-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
          !dudt(i,j,k) = ( gacc(1)*(rhox_ob-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
                          dxi*sigmaca*curvx*(vof(ip,j,k)-vof(i,j,k)) &
                        )/rhox + &
#ifdef FFT_DIRECT
                        dudt(i,j,k)*rho0i  &
                        - (1.d0/rhox-rho0i)*( pp(ip,j,k)-pp(i,j,k) )*dxi 
#else
                        dudt(i,j,k)/rhox
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momxp
  !
  subroutine momyp(nx,ny,nz,dyi,p,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dvdt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in )                         :: dyi 
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: p,pp,curv,vof,rho,rho_gas
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in )                         :: rho0,rho_av 
    real(8), intent(out), dimension(  :,  :,  :) :: dvdt
    !
    real(8) :: rhoy,rho0i,curvy
#ifndef LOW_MACH
    real(8) :: rhoy_ob,voff,tmpf,rhogf
#endif
    integer :: i,j,k
    integer :: jp
    !
    rho0i = 1.d0/rho0
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,jp) &
    !$OMP PRIVATE(rhoy,curvy) &
    !$OMP SHARED(nx,ny,nz,dyi,p,rho,curv,dvdt)
    do k=1,nz
      do j=1,ny
        jp = j + 1
        do i=1,nx
          !
          rhoy    = 0.5d0*(rho( i,jp,k)+rho( i,j,k))
          curvy   = 0.5d0*(curv(i,jp,k)+curv(i,j,k))
#ifdef LOW_MACH
          !rhogf   = 0.5d0*(rho_gas(i,jp,k)+rho_gas(i,j,k))
#else
          voff    = 0.5d0*(vof( i,jp,k)+vof( i,j,k))
          tmpf    = 0.5d0*(tmp( i,jp,k)+tmp( i,j,k))
          rhogf   = rho2_0*(1.d0-beta_g_th*(tmpf-tmp0))
#endif
          !rhoy_ob = voff*rho1*(1.d0-beta_l_th*(tmpf-tmp0)) + (1.d0-voff)*rhogf
          !
          dvdt(i,j,k) = - ( p(i,jp,k)-p(i,j,k) )*dyi !- dpdl(2)
          dvdt(i,j,k) = ( gacc(2)*(rhoy-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
          !dvdt(i,j,k) = ( gacc(2)*(rhoy_ob-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
                          dyi*sigmaca*curvy*(vof(i,jp,k)-vof(i,j,k)) &
                        )/rhoy + &
#ifdef FFT_DIRECT
                        dvdt(i,j,k)*rho0i &
                        - (1.d0/rhoy-rho0i)*( pp(i,jp,k)-pp(i,j,k) )*dyi 
#else
                        dvdt(i,j,k)/rhoy
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momyp
  !
  subroutine momzp(nx,ny,nz,dzci,p,pp,curv,vof,rho,rho_gas,tmp,rho0,rho_av,dwdt)
    !
    implicit none
    !
    integer, intent(in )                         :: nx,ny,nz
    real(8), intent(in ), dimension(-2:)         :: dzci
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: p,pp,curv,vof,rho,rho_gas
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: tmp
    real(8), intent(in )                         :: rho0,rho_av
    real(8), intent(out), dimension(  :,  :,  :) :: dwdt
    !
    real(8) :: rhoz,rho0i,curvz
#ifndef LOW_MACH
    real(8) :: rhoz_ob,voff,tmpf,rhogf
#endif
    integer :: kp
    integer :: i,j,k
    !
    rho0i = 1.d0/rho0
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k,kp) &
    !$OMP PRIVATE(rhoz,curvz) &
    !$OMP SHARED(nx,ny,nz,p,rho,curv,dwdt,dzci)
    do k=1,nz
      kp = k + 1
      do j=1,ny
        do i=1,nx
          !
          rhoz    = 0.5d0*(rho( i,j,kp)+rho( i,j,k))
          curvz   = 0.5d0*(curv(i,j,kp)+curv(i,j,k))
          !
#ifdef LOW_MACH
          !rhogf   = 0.5d0*(rho_gas(i,j,kp)+rho_gas(i,j,k))
#else
          voff    = 0.5d0*(vof( i,j,kp)+vof( i,j,k))
          tmpf    = 0.5d0*(tmp( i,j,kp)+tmp( i,j,k))
          rhogf   = rho2_0*(1.d0-beta_g_th*(tmpf-tmp0))
#endif
          !rhoz_ob = voff*rho1*(1.d0-beta_l_th*(tmpf-tmp0)) + (1.d0-voff)*rhogf
          !
          dwdt(i,j,k) = - ( p(i,j,kp)-p(i,j,k) )*dzci(k) !- dpdl(3)
          dwdt(i,j,k) = ( gacc(3)*(rhoz-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
          !dwdt(i,j,k) = ( gacc(3)*(rhoz_ob-rho_av) + & ! in case it is periodic, subtract the mean gravitational force per unit mass
                          dzci(k)*sigmaca*curvz*(vof(i,j,kp)-vof(i,j,k)) &
                        )/rhoz + &
#ifdef FFT_DIRECT
                        dwdt(i,j,k)*rho0i &
                        - (1.d0/rhoz-rho0i)*( pp(i,j,kp)-pp(i,j,k) )*dzci(k) 
#else
                        dwdt(i,j,k)/rhoz
#endif
          !
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    return
  end subroutine momzp
  !

  subroutine turb_forc_src(turb_type,c_or_t,nx,ny,nz,dx,dy,dz,lx,ly,lz,abc1,abc2,abc3,amp_a,amp_b,amp_n, &
                           time,nu2,k0_wave,k0_freq,rho,dudt,dvdt,dwdt)
    !
    use mod_param, only: pi
    !
    implicit none
    !
    character(len=3), intent(in   )                      :: turb_type  
    character(len=3), intent(in   )                      :: c_or_t 
    integer         , intent(in   )                      :: nx,ny,nz
    real(8)         , intent(in   )                      :: dx,dy,dz
    real(8)         , intent(in   )                      :: lx,ly,lz
    real(8)         , intent(in   )                      :: abc1,abc2,abc3
    real(8)         , intent(in   )                      :: amp_a,amp_b,amp_n
    real(8)         , intent(in   )                      :: time,nu2,k0_wave
    integer         , intent(in   )                      :: k0_freq ! 1, 2, 3, ...
    real(8)         , intent(in   ), dimension(0:,0:,0:) :: rho
    real(8)         , intent(inout), dimension(1:,1:,1:) :: dudt,dvdt,dwdt
    !
    real(8) :: rhox,rhoy,rhoz
    real(8) :: f0_t,tper
    real(8) :: xc,yc,zc,xf,yf
    integer :: i,j,k,ip,jp,kp
    !
    f0_t = nu2*k0_wave**2
    tper = 1.d0
    select case(c_or_t)
    case('uni')
      f0_t = f0_t
      tper = tper
    case('tdp')
      tper = (lx/(2.d0*pi*abc1))
      f0_t = f0_t*( amp_a+0.5d0*amp_b*sin(time*((2.d0*pi)/(amp_n*tper))) )
    end select
    !
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          zc = (k            -0.5d0)*dz/lz*2.d0*pi
          yc = (j+coord(2)*ny-0.5d0)*dy/ly*2.d0*pi
          yf = (j+coord(2)*ny-0.0d0)*dy/ly*2.d0*pi
          xc = (i+coord(1)*nx-0.5d0)*dx/lx*2.d0*pi
          xf = (i+coord(1)*nx-0.0d0)*dx/lx*2.d0*pi
          !
          rhox = 0.5d0*(rho(i+1,j,k)+rho(i,j,k))
          rhoy = 0.5d0*(rho(i,j+1,k)+rho(i,j,k))
          rhoz = 0.5d0*(rho(i,j,k+1)+rho(i,j,k))
          !
          if(    turb_type.eq.'tgv') then
            dudt(i,j,k) = dudt(i,j,k) + (f0_t*sin(k0_freq*xf)*cos(k0_freq*yc)*cos(k0_freq*zc))
            dvdt(i,j,k) = dvdt(i,j,k) - (f0_t*cos(k0_freq*xc)*sin(k0_freq*yf)*cos(k0_freq*zc))
          elseif(turb_type.eq.'abc') then
            dudt(i,j,k) = dudt(i,j,k) + f0_t*(abc1*sin(k0_freq*zc) + abc3*cos(k0_freq*yc)) 
            dvdt(i,j,k) = dvdt(i,j,k) + f0_t*(abc2*sin(k0_freq*xc) + abc1*cos(k0_freq*zc)) 
            dwdt(i,j,k) = dwdt(i,j,k) + f0_t*(abc3*sin(k0_freq*yc) + abc2*cos(k0_freq*xc)) 
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine turb_forc_src 
  !
end module mod_mom
