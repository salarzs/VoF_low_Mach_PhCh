module mod_vof_to_ls
  !
  ! module written by Francesco De Vita 
  ! modified by Nicolo Scapin for compatibility
  !
  ! This module allows to reconstruct a level-set distance field from a given VoF field. 
  ! The procedure is based on two papers:
  !     1. Russo and Smereka (JCP2000) for the re-distancing algorithm;
  !     2. Albadawi et. al (IJMF2013) for the initial value of the distance function
  !        and the time-step restriction.
  !
  ! Note: ensure that dl(1)=dl(2)=dl(3) to employ this module.
  !
  implicit none
  !
  private 
  public  :: vof_to_ls
  !
  contains
  !
  subroutine vof_to_ls(n,dli,dzci,dzfi,cbcphi,bcphi,vof,ls)
    !
    use mod_bound, only: boundsb
    !
    implicit none
    !
    integer         , intent(in ), dimension(3)           :: n
    real(8)         , intent(in ), dimension(3)           :: dli
    character(len=1), intent(in ), dimension(0:1,3)       :: cbcphi
    real(8)         , intent(in ), dimension(0:1,3)       :: bcphi
    real(8)         , intent(in ), dimension(-2:)         :: dzci,dzfi
    real(8)         , intent(in ), dimension( 0:, 0:, 0:) :: vof
    real(8)         , intent(out), dimension(-2:,-2:,-2:) :: ls
    !
    real(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: hs
    real(8), dimension(  n(1)  ,  n(2)  ,  n(3)  ) :: lsSgn,G,Dd
    real(8), dimension(-2:n(3)+3) :: dzc,dzf
    real(8), dimension(3) :: dl
    real(8) :: dlmin,dlmini
    real(8) :: d1,d2,d3,d4,d5,d6,d7,d8,d9,dt_ps
    integer :: i,j,k,ip,jp,kp,im,jm,km,s
    !
    integer, parameter :: niter = 60 ! sixty iteration 0.1*60 = 6 cells
    real(8), parameter :: small = 1.0e-12
    logical, parameter :: first_order = .false., second_order = .true.
    !
    dl  = dli**(-1)
    dzc = dzci**(-1)
    dzf = dzfi**(-1)
    !dt_ps  = 0.1d0*dl(1)
    !dt_ps  = 0.1d0*dl(2)
    dlmin  = min(minval(dl(:)),minval(dzf(:)))
    dlmini = 1.d0/dlmin
    dt_ps  = 0.1d0*dlmin
    !
    ! a. initial value of the distance function taken from VoF function
    !     see Albadawi et. al (IJMF2013)
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          !
          !hs(i,j,k) = (2.0d0*vof(i,j,k)-1.0d0)*dl(1)*0.75d0
          !hs(i,j,k) = (2.0d0*vof(i,j,k)-1.0d0)*dl(2)*0.75d0
          hs(i,j,k) = (2.0d0*vof(i,j,k)-1.0d0)*dlmin*0.75d0
          !
        end do
      end do
    end do
    !
    ! b. sgn function based on the initial value hs
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
          ls(i,j,k) = hs(i,j,k) ! initial value (Note: probably, lsSign is not needed to be stored)
          if    (ls(i,j,k).gt.0.d0) then 
            lsSgn(i,j,k) = +1.d0
          elseif(ls(i,j,k).lt.0.d0) then
            lsSgn(i,j,k) = -1.d0
          else
            lsSgn(i,j,k) = 0.d0
          endif
          !
          d1 = 0.5d0*abs(hs(ip,j,k)-hs(im,j,k)) 
          d2 = 0.5d0*abs(hs(i,jp,k)-hs(i,jm,k)) 
          d3 = 0.5d0*abs(hs(i,j,kp)-hs(i,j,km)) 
          d4 = abs(hs(ip,j,k)-hs(i ,j,k))
          d5 = abs(hs(i,jp,k)-hs(i,j ,k))
          d6 = abs(hs(i,j,kp)-hs(i,j,k ))
          d7 = abs(hs(i ,j,k)-hs(im,j,k))
          d8 = abs(hs(i,j ,k)-hs(i,jm,k))
          d9 = abs(hs(i,j,k )-hs(i,j,km))
          !
          !Dd(i,j,k) = dl(1)*hs(i,j,k)/max(sqrt(d1**2+d2**2+d3**3),sqrt(d4**2+d5**2+d6**2),sqrt(d7**2+d8**2+d9**2),small)
          !Dd(i,j,k) = dl(2)*hs(i,j,k)/max(sqrt(d1**2+d2**2+d3**3),sqrt(d4**2+d5**2+d6**2),sqrt(d7**2+d8**2+d9**2),small)
          Dd(i,j,k) = dlmin*hs(i,j,k)/max(sqrt(d1**2+d2**2+d3**3),sqrt(d4**2+d5**2+d6**2),sqrt(d7**2+d8**2+d9**2),small)
          !
        end do
      end do
    end do
    call boundsb(cbcphi,n,bcphi,dl,dzc,dzf,ls)
    !
    ! c: apply the redistancing algorithm
    ! 
    do s = 1,niter
      !  
      ! c1. choose the order of accuracy for the gradient
      !
      if(first_order ) call first_order_gradient( n,dli,hs,ls,G)
      if(second_order) call second_order_gradient(n,dli,hs,ls,G)
      !
      ! c2. advance in time
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
            if( hs(i,j,k)*hs(im,j,k).lt.0.d0.or.hs(i,j,k)*hs(ip,j,k).lt.0.d0.or. &
                hs(i,j,k)*hs(i,jm,k).lt.0.d0.or.hs(i,j,k)*hs(i,jp,k).lt.0.d0.or. &
                hs(i,j,k)*hs(i,j,km).lt.0.d0.or.hs(i,j,k)*hs(i,j,kp).lt.0.d0 ) then
              !
              !ls(i,j,k) = ls(i,j,k)-(dt_ps*dli(1))*(lsSgn(i,j,k)*abs(ls(i,j,k))-Dd(i,j,k))
              !ls(i,j,k) = ls(i,j,k)-(dt_ps*dli(2))*(lsSgn(i,j,k)*abs(ls(i,j,k))-Dd(i,j,k))
              ls(i,j,k) = ls(i,j,k)-(dt_ps*dlmini)*(lsSgn(i,j,k)*abs(ls(i,j,k))-Dd(i,j,k))
              !
            else
              !
              ls(i,j,k) = ls(i,j,k)-dt_ps*lsSgn(i,j,k)*G(i,j,k)
              !
            endif
            !
          end do
        end do
      end do
      call boundsb(cbcphi,n,bcphi,dl,dzc,dzf,ls)
      !
    end do
    !
    return
  end subroutine vof_to_ls
  ! 
  subroutine first_order_gradient(n,dli,hs,ls,G)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: hs
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: ls
    real(8), intent(out), dimension(  :,  :,  :) :: G
    !
    integer :: i,j,k!,c_t
    integer :: ip,jp,kp,im,jm,km
    real(8) :: a,b,c,d,e,f,ap,am,bp,bm,cp,cm,dp,dm,ep,em,fp,fm
    real(8) :: c_t 
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
          ! compute 
          !  -> a ,b ,c ,d ,e ,f
          !  -> ap,bp,cp,dp,ep,fp
          !  -> am,bm,cm,dm,em,fm
          !
          a  = dli(1)*(ls(i,j,k )-ls(im,j,k))
          b  = dli(1)*(ls(ip,j,k)-ls(i,j,k ))
          c  = dli(2)*(ls(i,j,k )-ls(i,jm,k))
          d  = dli(2)*(ls(i,jp,k)-ls(i,j,k ))
          e  = dli(3)*(ls(i,j,k )-ls(i,j,km))
          f  = dli(3)*(ls(i,j,kp)-ls(i,j,k ))
          !
          ap = max(a,0.d0)
          am = min(a,0.d0)
          bp = max(b,0.d0)
          bm = min(b,0.d0)
          cp = max(c,0.d0)
          cm = min(c,0.d0)
          dp = max(d,0.d0)
          dm = min(d,0.d0)
          ep = max(e,0.d0)
          em = min(e,0.d0)
          fp = max(f,0.d0)
          fm = min(f,0.d0)
          !
          ! compute the gradient, G
          !
          c_t      = sign(1.d0,hs(i,j,k))
          G(i,j,k) = 0.5d0*(c_t+abs(c_t))*(sqrt(max(ap**2,bm**2)+max(cp**2,dm**2)+max(ep**2,fm**2))-1.d0) - &
                     0.5d0*(c_t-abs(c_t))*(sqrt(max(am**2,bp**2)+max(cm**2,dp**2)+max(em**2,fp**2))-1.d0)
          !if    (hs(i,j,k).gt.0.d0) then
          !  G(i,j,k) = sqrt(max(ap**2,bm**2)+max(cp**2,dm**2)+max(ep**2,fm**2))-1.d0
          !elseif(hs(i,j,k).lt.0.d0) then
          !  G(i,j,k) = sqrt(max(am**2,bp**2)+max(cm**2,dp**2)+max(em**2,fp**2))-1.d0
          !endif
          !
        end do
      end do
    end do
    !
    return
  end subroutine first_order_gradient
  !
  subroutine second_order_gradient(n,dli,hs,ls,G)
    !
    implicit none
    !
    integer, intent(in ), dimension(3)           :: n
    real(8), intent(in ), dimension(3)           :: dli
    real(8), intent(in ), dimension( 0:, 0:, 0:) :: hs
    real(8), intent(in ), dimension(-2:,-2:,-2:) :: ls
    real(8), intent(out), dimension(  :,  :,  :) :: G
    !
    integer :: i,j,k,p
    real(8) :: c_t
    real(8) :: a,b,c,d,e,f,ap,am,bp,bm,cp,cm,dp,dm,ep,em,fp,fm,ccp,ccm
    real(8), dimension(-2:2)      :: xp,yp,zp,fp_i,fp_j,fp_k
    real(8), dimension(-2:2,-2:2) :: Phi_i,Phi_j,Phi_k
    real(8), dimension(3) :: dl
    !
    dl = dli**(-1) 
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          ! x-direction
          !
          do p = -2,2
            xp(p) = p*dl(1)
            fp_i(p) = ls(i+p,j,k)
          end do
          !
          ! y-direction
          !
          do p = -2,2
            yp(p) = p*dl(2)
            fp_j(p) = ls(i,j+p,k)
          end do
          !
          ! z-direction
          !
          do p = -2,2
            zp(p) = p*dl(3)
            fp_k(p) = ls(i,j,k+p)
          end do
          !
          ! compute the table of divided differences
          !
          do p = -2,1
            Phi_i(p,p+1) = (fp_i(p+1)-fp_i(p))/(xp(p+1)-xp(p))
            Phi_j(p,p+1) = (fp_j(p+1)-fp_j(p))/(yp(p+1)-yp(p))
            Phi_k(p,p+1) = (fp_k(p+1)-fp_k(p))/(zp(p+1)-zp(p))
          end do
          do p = -2,0
            Phi_i(p,p+2) = (Phi_i(p+1,p+2)-Phi_i(p,p+1))/(xp(p+2)-xp(p))
            Phi_j(p,p+2) = (Phi_j(p+1,p+2)-Phi_j(p,p+1))/(yp(p+2)-yp(p))
            Phi_k(p,p+2) = (Phi_k(p+1,p+2)-Phi_k(p,p+1))/(zp(p+2)-zp(p))
          end do
          !
          ! compute 
          !  -> a ,b ,c ,d ,e ,f
          !  -> ap,bp,cp,dp,ep,fp
          !  -> am,bm,cm,dm,em,fm
          !
          ccm = minmod(Phi_i(-2,0),Phi_i(-1,1))
          ccp = minmod(Phi_i(-1,1),Phi_i( 0,2))
          a = Phi_i(-1,0) + ccm*(xp(0)-xp(-1))
          b = Phi_i( 0,1) + ccp*(xp(0)-xp( 1))
          !
          ccm = minmod(Phi_j(-2,0),Phi_j(-1,1))
          ccp = minmod(Phi_j(-1,1),Phi_j( 0,2))
          c = Phi_j(-1,0) + ccm*(yp(0)-yp(-1))
          d = Phi_j( 0,1) + ccp*(yp(0)-yp( 1))
          !
          ccm = minmod(Phi_k(-2,0),Phi_k(-1,1))
          ccp = minmod(Phi_k(-1,1),Phi_k( 0,2))
          e = Phi_k(-1,0) + ccm*(zp(0)-zp(-1))
          f = Phi_k( 0,1) + ccp*(zp(0)-zp( 1))
          !
          ap = max(a,0.d0)
          am = min(a,0.d0)
          bp = max(b,0.d0)
          bm = min(b,0.d0)
          cp = max(c,0.d0)
          cm = min(c,0.d0)
          dp = max(d,0.d0)
          dm = min(d,0.d0)
          ep = max(e,0.d0)
          em = min(e,0.d0)
          fp = max(f,0.d0)
          fm = min(f,0.d0)
          !
          ! compute the gradient, G
          !
          c_t      = sign(1.d0,hs(i,j,k))
          G(i,j,k) = 0.5d0*(c_t+abs(c_t))*(sqrt(max(ap**2,bm**2)+max(cp**2,dm**2)+max(ep**2,fm**2))-1.d0) - &
                     0.5d0*(c_t-abs(c_t))*(sqrt(max(am**2,bp**2)+max(cm**2,dp**2)+max(em**2,fp**2))-1.d0)
          !if    (hs(i,j,k).gt.0.d0) then
          !  G(i,j,k) = sqrt(max(ap**2,bm**2)+max(cp**2,dm**2)+max(ep**2,fm**2))-1.d0
          !elseif(hs(i,j,k).lt.0.d0) then
          !  G(i,j,k) = sqrt(max(am**2,bp**2)+max(cm**2,dp**2)+max(em**2,fp**2))-1.d0
          !endif
          !
        end do
      end do
    end do
    !
    return
  end subroutine second_order_gradient
  !
  function minmod(alpha,beta) result(c)
    !
    implicit none
    !
    real(8), intent(in) :: alpha, beta
    real(8)             :: c
    !
    if(alpha*beta.gt.0.d0) then
      if(abs(alpha).le.abs(beta)) then
        c = alpha
      else
        c = beta
      endif
    else
      c = 0.d0
    endif
    !
    return
  end function minmod
  !
end module mod_vof_to_ls
