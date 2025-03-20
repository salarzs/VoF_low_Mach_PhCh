module mod_mom
  !
  implicit none
  !
  private
  public momxad,momyad,momzad
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
    real(8), intent(out), dimension( 0:, 0:, 0:) :: dudt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    real(8) :: dudxp,dudxm,dvdxp,dvdxm,dudyp,dudym, &
               dudzp,dudzm,dwdxp,dwdxm
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divxm,divxp
    real(8) :: rhox
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
          uuip  = 0.25d0*( u(ip,j,k)+u(i,j,k) )*( u(ip,j ,k )+u(i,j ,k ) )
          uuim  = 0.25d0*( u(im,j,k)+u(i,j,k) )*( u(im,j ,k )+u(i,j ,k ) )
          uvjp  = 0.25d0*( u(i,jp,k)+u(i,j,k) )*( v(ip,j ,k )+v(i,j ,k ) )
          uvjm  = 0.25d0*( u(i,jm,k)+u(i,j,k) )*( v(ip,jm,k )+v(i,jm,k ) )
          uwkp  = 0.25d0*( u(i,j,kp)+u(i,j,k) )*( w(ip,j ,k )+w(i,j ,k ) )
          uwkm  = 0.25d0*( u(i,j,km)+u(i,j,k) )*( w(ip,j ,km)+w(i,j ,km) )
          !
          !divxp = (W(i+1,j,k)-W(i+1,j,k-1))*dzi + &
          !        (V(i+1,j,k)-V(i+1,j-1,k))*dyi + &
          !        (U(i+1,j,k)-U(i+1-1,j,k))*dxi
          !divxm = (W(i,j,k)-W(i,j,k-1))*dzi + &
          !        (V(i,j,k)-V(i,j-1,k))*dyi + &
          !        (U(i,j,k)-U(i-1,j,k))*dxi
          !
          divxp = (w(ip,j,k)-w(ip,j ,km))*dzi + &
                  (v(ip,j,k)-v(ip,jm,k ))*dyi + &
                  (u(ip,j,k)-u(i ,j ,k ))*dxi
          !
          divxm = (w(i ,j,k)-w(i ,j ,km))*dzi + &
                  (v(i ,j,k)-v(i ,jm,k ))*dyi + &
                  (u(i ,j,k)-u(im,j ,k ))*dxi
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
                        dxi*( -uuip + uuim ) + &
                        dyi*( -uvjp + uvjm ) + &
                        dzi*( -uwkp + uwkm ) + &
                        u(i,j,k)*(divxp+divxm) + &
                        !
                        ( &
                        dxi*((dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm - (2.d0/3.d0)*(muxp*divxp-muxm*divxp)) + &
                        dyi*((dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym                                      ) + &
                        dzi*((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm                                      )   &
                        )/rhox
                        !
          !
          ! Momentum balance
          !
          !advnU(i,j,k) = dxi*( -uuip + uuim ) + &
          !               dyi*( -uvjp + uvjm ) + &
          !               dzi*( -uwkp + uwkm ) + &
          !               U(i,j,k)*(divxp+divxm)*(vof(i,j,k)+vof(i+1,j,k))*0.25
          !
          !dfunU(i,j,k) = 0.5*(vi(i,j,k)+vi(i+1,j,k))*((dudxp-dudxm)*dxi + &
          !                                            (dudyp-dudym)*dyi + &
          !                                            (dudzp-dudzm)*dzi + &
          !                                   (-2./3.)*(divxp-divxm)*dxi*(vof(i,j,k)+vof(i+1,j,k))*0.0)
          ! 
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
    real(8), intent(out), dimension( 0:, 0:, 0:) :: dvdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    real(8) :: dvdxp,dvdxm,dudyp,dudym,dvdyp,dvdym, &
               dvdzp,dvdzm,dwdyp,dwdym
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divym,divyp
    real(8) :: rhoy
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
          uvip  = 0.25d0*( u(i ,j,k)+u(i ,jp,k) )*( v(i,j,k )+v(ip,j ,k) )
          uvim  = 0.25d0*( u(im,j,k)+u(im,jp,k) )*( v(i,j,k )+v(im,j ,k) )
          vvjp  = 0.25d0*( v(i,j,k )+v(i,jp,k)  )*( v(i,j,k )+v(i ,jp,k) )
          vvjm  = 0.25d0*( v(i,j,k )+v(i,jm,k)  )*( v(i,j,k )+v(i ,jm,k) )
          wvkp  = 0.25d0*( w(i,j,k )+w(i,jp,k)  )*( v(i,j,kp)+v(i ,j ,k) )
          wvkm  = 0.25d0*( w(i,j,km)+w(i,jp,km) )*( v(i,j,km)+v(i ,j ,k) )
          !
          divyp = (w(i,jp,k)-w(i ,jp,km))*dzi + &
                  (v(i,jp,k)-v(i ,j ,k ))*dyi + &
                  (u(i,jp,k)-u(im,jp,k ))*dxi
          !
          divym = (w(i,j ,k)-w(i ,j ,km))*dzi + &
                  (v(i,j ,k)-v(i ,jm,k ))*dyi + &
                  (u(i,j ,k)-u(im,j ,k ))*dxi
          !
          !divyp = (W(i,j+1,k)-W(i,j+1,k-1))*dzi + &
          !        (V(i,j+1,k)-V(i,j+1-1,k))*dyi + &
          !        (U(i,j+1,k)-U(i-1,j+1,k))*dxi
          !divym = (W(i,j,k)-W(i,j,k-1))*dzi + &
          !        (V(i,j,k)-V(i,j-1,k))*dyi + &
          !        (U(i,j,k)-U(i-1,j,k))*dxi
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
                        dxi*( -uvip + uvim ) + &
                        dyi*( -vvjp + vvjm ) + & 
                        dzi*( -wvkp + wvkm ) + &
                        v(i,j,k)*(divyp+divym) + &
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
    real(8), intent(out), dimension( 0:, 0:, 0:) :: dwdt
    !
    integer :: im,ip,jm,jp,km,kp,i,j,k
    real(8) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    real(8) :: dwdxp,dwdxm,dudzp,dudzm,dwdyp,dwdym, &
               dvdzp,dvdzm,dwdzp,dwdzm
    real(8) :: muxm,muxp,muym,muyp,muzm,muzp
    real(8) :: divzm,divzp
    real(8) :: rhoz
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
          uwip  = 0.25d0*( w(i,j,k)+w(ip,j,k) )*( u(i ,j ,k)+u(i ,j ,kp) )
          uwim  = 0.25d0*( w(i,j,k)+w(im,j,k) )*( u(im,j ,k)+u(im,j ,kp) )
          vwjp  = 0.25d0*( w(i,j,k)+w(i,jp,k) )*( v(i ,j ,k)+v(i ,j ,kp) )
          vwjm  = 0.25d0*( w(i,j,k)+w(i,jm,k) )*( v(i ,jm,k)+v(i ,jm,kp) )
          wwkp  = 0.25d0*( w(i,j,k)+w(i,j,kp) )*( w(i ,j ,k)+w(i ,j ,kp) )
          wwkm  = 0.25d0*( w(i,j,k)+w(i,j,km) )*( w(i ,j ,k)+w(i ,j ,km) )
          !    
          !divzp = (W(i,j,k+1)-W(i,j,k-1+1))*dzi + &
          !        (V(i,j,k+1)-V(i,j-1,k+1))*dyi + &
          !        (U(i,j,k+1)-U(i-1,j,k+1))*dxi
          !divzm = (W(i,j,k)-W(i,j,k-1))*dzi + &
          !        (V(i,j,k)-V(i,j-1,k))*dyi + &
          !        (U(i,j,k)-U(i-1,j,k))*dxi
          !
          divzp = (w(i,j,kp)-w(i,j ,k ))*dzi + &
                  (v(i,j,kp)-v(i,jm,kp))*dyi + &
                  (u(i,j,kp)-u(im,j,kp))*dxi
          !
          divzm = (w(i,j,k)-w(i ,j ,km))*dzi + &
                  (v(i,j,k)-v(i ,jm,k ))*dyi + &
                  (u(i,j,k)-u(im,j ,k ))*dxi
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
                        dxi*( -uwip + uwim ) + &
                        dyi*( -vwjp + vwjm ) + &
                        dzi*( -wwkp + wwkm ) + &
                        W(i,j,k)*(divzp+divzm) + &
                        !
                        ( &
                        dxi*((dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm                                      ) + &
                        dyi*((dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym                                      ) + &
                        dzi*((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm - (2.d0/3.d0)*(muzp*divzp-muzm*divzm)) & 
                        )/rhoz
                        !
         ! advnW(i,j,k) = dxi*( -uwip + uwim ) + &
         !                dyi*( -vwjp + vwjm ) + &
         !                dzi*( -wwkp + wwkm ) + &
         !                W(i,j,k)*(divzp+divzm)*(vof(i,j,k)+vof(i,j,k+1))*0.25
    
         ! dfunW(i,j,k) = 0.5*(vi(i,j,k)+vi(i,j,k+1))*((dwdxp-dwdxm)*dxi + &
         !                                             (dwdyp-dwdym)*dyi + &
         !                                             (dwdzp-dwdzm)*dzi + &
         !                                    (-2./3.)*(divzp-divzm)*dzi*(vof(i,j,k)+vof(i,j,k+1))*0.0)
    
        enddo
      enddo
    enddo
    !$omp end parallel
    !
    return
  end subroutine momzad
  !
end module mod_mom
