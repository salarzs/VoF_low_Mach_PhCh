module mod_tagging
  !
  use mpi
  use mod_common_mpi, only: myid,comm_cart,ierr,coord,status, &
                            left,right,front,back,xhalo,yhalo, &
                            xsl_buf, xrl_buf, ysl_buf, yrl_buf, xsr_buf, xrr_buf, ysr_buf, yrr_buf
  use mod_param     , only: dl,dims,lheat,datadir_ta,tmpg0
  use mod_thermo    , only: mass_fraction,thermo_kag,thermo_d_lg,thermo_rhog
  !
  !implicit none
  !
  !real(8), parameter :: vof_th = 0.5d0
  !real(8), parameter :: vof_th = 0.10d0
  !integer, parameter :: num    = 3
  !integer, parameter :: rdir   = +1
  real(8) :: vof_th = 0.10d0
  integer :: num    = 3
  integer :: rdir   = +1
  !
  !private
  !public  :: droplet_tagging
  !
  contains
  !
  subroutine droplet_tagging(n,e,h,dli,vof,delta,kappa,u,v,w,phi,tmp,tmpge,tmple,sca,mflux,pth,istep,time)
    !
    implicit none
    !
    integer, intent(in), dimension(3)        :: n
    integer, intent(in)                      :: e,h
    real(8), intent(in), dimension(3)        :: dli
    real(8), intent(in), dimension(0:,0:,0:) :: vof,delta,kappa
    real(8), intent(in), dimension(e:,e:,e:) :: u,v,w
    real(8), intent(in), dimension(h:,h:,h:) :: tmp
    real(8), intent(in), dimension(h:,h:,h:) :: phi
    real(8), intent(in), dimension(0:,0:,0:) :: tmpge,tmple
    real(8), intent(in), dimension(0:,0:,0:) :: sca
    real(8), intent(in), dimension(0:,0:,0:) :: mflux
    real(8), intent(in)                      :: pth
    integer, intent(in)                      :: istep
    real(8), intent(in)                      :: time
    !
    integer(8), dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vofTag
    real(8)   , dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1) :: vof_f
    real(8)   , dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1,num) :: phi_v
    real(8)   , dimension(8) :: mx,my,mz 
    real(8) :: norm
    real(8) :: dtmpdxp,dtmpdxm,dtmpdyp,dtmpdym,dtmpdzp,dtmpdzm, &
               lap_kap,kappaxp,kappaxm,kappayp,kappaym,kappazp,kappazm, &
               dvofdx,dvofdy,dvofdz,tmpg_i,tmpl_i,s_in,d,yg_i
    integer :: i,j,k,b,nd,cellcount,ip,jp,kp,im,jm,km,p
    integer(4) :: id, idst
    integer(4), allocatable, dimension(:,:) :: procShare, share,idShare, orderField
    integer(4), allocatable, dimension(:)   :: myid_out, drid
    integer(8), allocatable, dimension(:)   :: idloc
    real(8),    allocatable, dimension(:)   :: xd,yd,zd,ud,vd,wd,vold,tmpd, & 
                                               area,smid,smad,smed, &
                                               rhgd,dlgd,kagd, &
                                               mfxd,qkfd,qksd,tmpi
    real(8), allocatable, dimension(:,:) :: tinf_v,yinf_v,rcoun
    integer, allocatable, dimension(:,:) :: counter_p2
    integer, allocatable, dimension(:)   :: counter_p1
    character(len=7) :: p_ch
    !
    !call random_number(rdnInit)
    ! idst = floor(rdnInit*100000+100000*myid) ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    idst   = 0 ! avoid possible ID duplication between processors. May fail for n_d >>10000???
    id     = idst+1
    vof_f  = vof
    vofTag = 0
    nd     = 0
    !
    ! main loop for tagging droplets --> vof_tag 
    !
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          if(vof_f(i,j,k).ge.vof_th) then
            !
            vof_f(i,j,k)  = 0.d0
            vofTag(i,j,k) = id
            call recursiveTagging(n,i,j,k,id,vof_f,vofTag)
            id = id+1
            nd = nd+1
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    ! Update halo
    !
    call updthalo((/n(1),n(2)/),1,vofTag)
    call updthalo((/n(1),n(2)/),2,vofTag)
    !vofTag_r = 1.d0*vofTag
    !call updthalo((/n(1),n(2)/),1,vofTag_r)
    !call updthalo((/n(1),n(2)/),2,vofTag_r)
    !vofTag   = ceiling(vofTag_r,8)
    !
    !call boundper_int(vofTag,time)
    !
    ! write(proidch,'(i3.3)') myid
    ! open(90+myid,file='proc'//proidch)
    ! do i=0,n(1)+1
    !   do j=0,n(2)+1
    !     do k=0,n(3)+1
    !       write(90+myid,'(4I5.4)') i,j,k,vofTag(i,j,k)
    !     enddo
    !   enddo
    ! enddo
    ! close(92)
    ! print*,'max vof', maxval(vof_f), myid
    ! allocate droplet properties
    allocate(xd(nd),yd(nd),zd(nd),ud(nd),vd(nd),wd(nd),vold(nd), & 
             area(nd),tmpd(nd),qkfd(nd),qksd(nd),smid(nd),smad(nd),smed(nd), &
             mfxd(nd),rhgd(nd),dlgd(nd),kagd(nd),tmpi(nd), &
             tinf_v(nd,num),yinf_v(nd,num), &
             counter_p1(nd),counter_p2(nd,num), &
             idloc(nd),share(nd,6),myid_out(nd),procShare(nd,6), &
             idShare(nd,6),orderField(nd,6),drid(nd))
    !
    ! volumetric quantities
    !
    xd   = 0.d0
    yd   = 0.d0
    zd   = 0.d0
    ud   = 0.d0
    vd   = 0.d0
    wd   = 0.d0
    vold = 0.d0
    tmpd = 0.d0
    !
    ! surface quantities
    !
    area = 0.d0
    tmpi = 0.d0
    qkfd = 0.d0
    qksd = 0.d0
    mfxd = 0.d0
    smid = 1.1d0 ! random number
    smad = 0.d0
    smed = 0.d0
    rhgd = 0.d0
    dlgd = 0.d0
    kagd = 0.d0
    tinf_v = 0.d0
    yinf_v = 0.d0
    counter_p1 = 0
    counter_p2 = 0
    !
    ! processors stuff
    ! 
    id        = 0
    cellcount = 0
    idShare   = 0
    procShare = -1
    share     = 0
    !
    ! 1. first compute the iso-surfaces away from the interface
    !    in the gas phase
    !
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)+1
          !
          do p=1,num
            phi_v(i,j,k,p) = phi(i,j,k) + p*rdir*dl(2)
          enddo
          !
        enddo
      enddo
    enddo
    !
    ! Loop for properties characterization
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
          if(vofTag(i,j,k).gt.0) then
            !
            id        = vofTag(i,j,k)-idst
            drid(id)  = id
            idloc(id) = vofTag(i,j,k) 
            !
            ! compute the quantities at inf
            !
            if(phi(i,j,k).lt.0.d0) then
              !
              do p=1,num
                !
                ! along x
                !
                if( phi_v(ip,j,k,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(ip,j,k,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(ip,j,k,p))+tmp(ip,j,k)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(ip,j,k,p))+sca(ip,j,k)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                if( phi_v(im,j,k,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(im,j,k,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(im,j,k,p))+tmp(im,j,k)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(im,j,k,p))+sca(im,j,k)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                !
                ! along y
                !
                if( phi_v(i,jp,k,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(i,jp,k,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,jp,k,p))+tmp(i,jp,k)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(i,jp,k,p))+sca(i,jp,k)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                if( phi_v(i,jm,k,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(i,jm,k,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,jm,k,p))+tmp(i,jm,k)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(i,jm,k,p))+sca(i,jm,k)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                !
                ! along z
                !
                if( phi_v(i,j,kp,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(i,j,kp,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,j,kp,p))+tmp(i,j,kp)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(i,j,kp,p))+sca(i,j,kp)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                if( phi_v(i,j,km,p)*phi_v(i,j,k,p).lt.0. ) then
                  d         = (abs(phi_v(i,j,km,p))+abs(phi_v(i,j,k,p)))
                  tmpg_i    = (tmp(i,j,k)*abs(phi_v(i,j,km,p))+tmp(i,j,km)*abs(phi_v(i,j,k,p)))/d
                  yg_i      = (sca(i,j,k)*abs(phi_v(i,j,km,p))+sca(i,j,km)*abs(phi_v(i,j,k,p)))/d
                  tinf_v(id,p) = tinf_v(id,p) + tmpg_i
                  yinf_v(id,p) = yinf_v(id,p) + yg_i
                  counter_p2(id,p) = counter_p2(id,p) + 1
                endif
                !
              enddo
              !
            endif
            ! 
            norm = delta(i,j,k)
            !
            ! maximum, minimum, mean interfacial values
            !
            if( phi(ip,j,k)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(ip,j,k))+tmple(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(ip,j,k))+tmpge(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            if( phi(im,j,k)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(im,j,k))+tmple(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(im,j,k))+tmpge(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            !
            if( phi(i,jp,k)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(i,jp,k))+tmple(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(i,jp,k))+tmpge(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            if( phi(i,jm,k)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(i,jm,k))+tmple(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(i,jm,k))+tmpge(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            !
            if( phi(i,j,kp)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(i,j,kp))+tmple(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(i,j,kp))+tmpge(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            if( phi(i,j,km)*phi(i,j,k).lt.0. ) then
              tmpl_i = (tmple(i,j,k)*abs(phi(i,j,km))+tmple(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              tmpg_i = (tmpge(i,j,k)*abs(phi(i,j,km))+tmpge(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              s_in   = mass_fraction(pth,tmpl_i)
              smad(id)    = max(s_in,smad(id))
              smid(id)    = min(s_in,smid(id))
              smed(id)    = smed(id) + s_in
              counter_p1(id) = counter_p1(id) + 1
              kagd(id)    = kagd(id) + thermo_kag(tmpg_i)
              tmpi(id)    = tmpi(id) + tmpl_i
              rhgd(id)    = rhgd(id) + thermo_rhog(pth,tmpg_i,s_in)
              dlgd(id)    = dlgd(id) + thermo_d_lg(pth,tmpg_i,s_in)
            endif
            !
            ! Here we compute the quantities for each droplet
            !
            ! 1. position
            !
            xd(id) = xd(id) + vof(i,j,k)*(i+coord(1)*n(1)-0.5d0)*dl(1)
            yd(id) = yd(id) + vof(i,j,k)*(j+coord(2)*n(2)-0.5d0)*dl(2)
            zd(id) = zd(id) + vof(i,j,k)*(k              -0.5d0)*dl(3)
            !
            ! 2. velocity
            !
            ud(id) = ud(id) + u(i,j,k)*vof(i,j,k)
            vd(id) = vd(id) + v(i,j,k)*vof(i,j,k)
            wd(id) = wd(id) + w(i,j,k)*vof(i,j,k)
            !
            ! 3. temperature
            !
            tmpd(id) = tmpd(id) + tmp(i,j,k)*vof(i,j,k)
            !
            ! 4. volume and area
            ! 
            vold(id) = vold(id) + vof(i,j,k)
            area(id) = area(id) + norm*dl(1)*dl(2)*dl(3)
            !
            ! 5. mflux
            !
            mfxd(id) = mfxd(id) + mflux(i,j,k)*norm*dl(1)*dl(2)*dl(3)
            !
            ! 6. latent/sensible heat
            !
            kappaxp = 0.5d0*(kappa(ip,j,k)+kappa(i,j,k))
            kappaxm = 0.5d0*(kappa(im,j,k)+kappa(i,j,k))
            kappayp = 0.5d0*(kappa(i,jp,k)+kappa(i,j,k))
            kappaym = 0.5d0*(kappa(i,jm,k)+kappa(i,j,k))
            kappazp = 0.5d0*(kappa(i,j,kp)+kappa(i,j,k))
            kappazm = 0.5d0*(kappa(i,j,km)+kappa(i,j,k))
            !
            dtmpdxp = (tmp(ip,j,k)-tmp(i ,j,k))*dli(1)
            dtmpdxm = (tmp(i ,j,k)-tmp(im,j,k))*dli(2)
            dtmpdyp = (tmp(i,jp,k)-tmp(i,j ,k))*dli(2)
            dtmpdym = (tmp(i,j ,k)-tmp(i,jm,k))*dli(2)
            dtmpdzp = (tmp(i,j,kp)-tmp(i,j,k ))*dli(3)
            dtmpdzm = (tmp(i,j,k )-tmp(i,j,km))*dli(3)
            !
            lap_kap = (kappaxp*dtmpdxp-kappaxm*dtmpdxm)*dli(1) + &
                      (kappayp*dtmpdyp-kappaym*dtmpdym)*dli(2) + &
                      (kappazp*dtmpdzp-kappazm*dtmpdzm)*dli(3)
            !
            qkfd(id) = qkfd(id) + lap_kap*vof(i,j,k)*dl(1)*dl(2)*dl(3)      ! integ. over the liq. vol.: (1/m)*(W/m*K)*(K/m)*m^3 = W
            qksd(id) = qksd(id) + mflux(i,j,k)*lheat*norm*dl(1)*dl(2)*dl(3) ! integ. over the liq. vol.: (kg/m^2*s)*(J/kg)*(1/m)*m^3 = J/s = W
            !
            ! print*, vof(i,j,k), i,j,k
            cellcount = 1
            if((i.eq.1).and.share(id,1).eq.0.and.vofTag(0,j,k).ne.0) then
              ! print*, 'i0i0i0i0i'
              share(id,1)     = 1
              procShare(id,1) = left
              idShare(id,1)   = vofTag(0,j,k)
              ! if (vofTag(0,j,k).ne.0)      idShare(id,1)   = vofTag(0,j,k)
            elseif((i.eq.n(1)).and.share(id,2).eq.0.and.vofTag(n(1)+1,j,k).ne.0) then
              ! print*, 'n1n1n1n1n1'
              share(id,2)     = 1
              procShare(id,2) = right
              idShare(id,2)   = vofTag(n(1)+1,j,k)
            elseif((j.eq.1).and.share(id,3).eq.0.and.vofTag(i,0,k).ne.0) then
              share(id,3)     = 1
              procShare(id,3) = front
              idShare(id,3)   = vofTag(i,0,k)
              ! print*, 'j1', idShare(id,:)
            elseif((j.eq.n(2)).and.share(id,4).eq.0.and.vofTag(i,n(2)+1,k).ne.0) then
              share(id,4)     = 1
              procShare(id,4) = back
              idShare(id,4)   = vofTag(i,n(2)+1,k)
              ! print*, 'n2', idShare(id,:)
            elseif((k.eq.1).and.share(id,5).eq.0) then
              share(id,5)     = 1
              procShare(id,5) = myid
              idShare(id,5)   = vofTag(i,j,0)
            elseif((k.eq.n(3)).and.share(id,6).eq.0) then
              share(id,6)     = 1
              procShare(id,6) = myid
              idShare(id,6)   = vofTag(i,j,n(3)+1)
            endif
            !
          endif
          !
        enddo
      enddo
    enddo 
    !
    print*, "ciao"
    do i=1,nd
      do j=1,6
        orderField(i,j) = myid*100+i*10+J
        ! print*,'proc', procShare(i,:), myid
        ! print*, nd, myid
        ! print*,'id  ', idShare(i,:), myid
      enddo
      ! print*,'id  ', idShare(i,:), myid
      ! print*,'proc', procShare(i,:), myid
    enddo
    !
    myid_out = 0
    do id=1,nd
      !
      myid_out(id) = myid
      xd(id)   = xd(id)/vold(id)    
      yd(id)   = yd(id)/vold(id)    
      zd(id)   = zd(id)/vold(id)    
      ud(id)   = ud(id)/vold(id)    
      vd(id)   = vd(id)/vold(id)    
      wd(id)   = wd(id)/vold(id)   
      tmpd(id) = tmpd(id)/vold(id)
      vold(id) = vold(id)*dl(1)*dl(2)*dl(3)
      qkfd(id) = qkfd(id)!/vold(id) ! we multiply by dl(1)*dl(2)*dl(3) when we compute qkfd
      area(id) = area(id)
      qksd(id) = qksd(id)!/area(id)
      mfxd(id) = mfxd(id)!/area(id)
      !
      tmpi(id) = tmpi(id)/(1.d0*counter_p1(id))
      smid(id) = smid(id)
      smed(id) = smed(id)/(1.d0*counter_p1(id))
      smad(id) = smad(id)
      rhgd(id) = rhgd(id)/(1.d0*counter_p1(id))
      dlgd(id) = dlgd(id)/(1.d0*counter_p1(id))
      kagd(id) = kagd(id)/(1.d0*counter_p1(id))
      !
      do p=1,num
        !
        if(counter_p2(id,p).eq.0) then
          tinf_v(id,p) = tmpg0
          yinf_v(id,p) = 0.d0
        else
          tinf_v(id,p) = tinf_v(id,p)/(1.d0*counter_p2(id,p))
          yinf_v(id,p) = yinf_v(id,p)/(1.d0*counter_p2(id,p))
        endif
        !
      enddo
      !
    enddo 
    !
    ! Writing results
    !
    call write_output_r8(istep,'xpos',xd  ,nd,'sum') ! position-x
    call write_output_r8(istep,'ypos',yd  ,nd,'sum') ! position-y
    call write_output_r8(istep,'zpos',zd  ,nd,'sum') ! position-z
    call write_output_r8(istep,'uvel',ud  ,nd,'sum') ! velocity-x
    call write_output_r8(istep,'vvel',vd  ,nd,'sum') ! velocity-y
    call write_output_r8(istep,'wvel',wd  ,nd,'sum') ! velocity-z
    call write_output_r8(istep,'tmpd',tmpd,nd,'sum') ! temperature of the liquid
    call write_output_r8(istep,'vold',vold,nd,'sum') ! volume
    !
    call write_output_r8(istep,'area',area,nd,'sum') ! area
    call write_output_r8(istep,'tmpi',tmpi,nd,'sum') ! surface temperature of the liquid
    call write_output_r8(istep,'mfxd',mfxd,nd,'sum') ! mflux
    call write_output_r8(istep,'qkfd',qkfd,nd,'sum') ! sensible heat
    call write_output_r8(istep,'qksd',qksd,nd,'sum') ! latent heat
    call write_output_r8(istep,'smid',smid,nd,'min') ! minimum Y_G
    call write_output_r8(istep,'smed',smed,nd,'sum') ! mean Y_G 
    call write_output_r8(istep,'smad',smad,nd,'max') ! maximum Y_G
    call write_output_r8(istep,'rhgd',rhgd,nd,'sum') ! gas density at Gamma
    call write_output_r8(istep,'dlgd',dlgd,nd,'sum') ! liquid-gas d_lg at Gamma
    call write_output_r8(istep,'kagd',kagd,nd,'sum') ! kag at Gamma
    !
    do p=1,num
      write(p_ch,'(i7.7)') p
      call write_output_r8(istep,'tinf_'//p_ch//'_v',tinf_v(:,p),nd,'sum') ! t at inf
      call write_output_r8(istep,'yinf_'//p_ch//'_v',yinf_v(:,p),nd,'sum') ! y at inf
    enddo
    !
    call write_output_i1(istep,'drid',drid,nd)
    call write_output_i1(istep,'idmy',myid_out,nd)
    !
    call write_output_i6(istep,'proc',transpose(procShare ),nd)
    call write_output_i6(istep,'orde',transpose(orderField),nd)
    call write_output_i6(istep,'idsh',transpose(idShare   ),nd)
    !
    !deallocate(xd,yd,zd,ud,vd,wd,vold,procShare,idShare)
    deallocate(xd,yd,zd,ud,vd,wd,vold,area, & 
               tmpd,qkfd,qksd,smid,smad,smed, &
               tmpi,mfxd,rhgd,dlgd,kagd, &
               tinf_v,yinf_v, &
               counter_p1,counter_p2, &
               idloc,share,myid_out,procShare, &
               idShare,orderField,drid)
    !
    return
  end subroutine droplet_tagging
  !
  recursive subroutine recursiveTagging(n,i,j,k,id,vof_f,vofTag)
    !
    implicit none
    !
    integer ,   intent(in   ), dimension(3)        :: n
    integer ,   intent(in   )                      :: i,j,k
    integer(4), intent(inout)                      :: id
    real(8)   , intent(inout), dimension(0:,0:,0:) :: vof_f
    integer(8), intent(inout), dimension(0:,0:,0:) :: vofTag
    !
    integer :: ii,jj,kk
    ! 
    do kk=k-1,k+1
      do jj=j-1,j+1
        do ii=i-1,i+1
          !
          if((vof_f(ii,jj,kk).gt.vof_th).and.(ii.ge.1   ).and.(jj.ge.1   ).and.(kk.ge.1   ).and.  &
                                             (ii.le.n(1)).and.(jj.le.n(2)).and.(kk.le.n(3))) then
            vof_f(ii,jj,kk)  = 0.d0
            vofTag(ii,jj,kk) = id
            ! print*, ii,jj,kk, vof_f(ii,jj,kk), vofTag(ii,jj,kk)
            call recursiveTagging(n,ii,jj,kk,id,vof_f,vofTag)
            !
          endif
          !
        enddo
      enddo
    enddo
    !
    return
  end subroutine recursiveTagging
  !
  subroutine write_output_r8(istep,varname,var,nd,mpi_type)
    !
    implicit none
    !
    integer,         intent(in)                :: istep
    character(len=*),intent(in)                :: varname
    real(8),         intent(in), dimension(1:) :: var
    integer,         intent(in)                :: nd
    character(len=3),intent(in)                :: mpi_type
    !
    integer                                    :: offset,fh
    character(7)                               :: fldnum
    integer(kind=MPI_OFFSET_KIND)              :: filesize,disp,nreals,nreals_myid
    character(len=1024)                        :: fname
    !
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    !
    !if(    mpi_type.eq.'max') then
    !  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    !elseif(mpi_type.eq.'sum') then
    !  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    !elseif(mpi_type.eq.'min') then
    !  call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr) ! why max does not work?
    !else
    !  if(myid.eq.0) print*, "ierr in tagging"
    !  call exit
    !endif
    !
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*8, MPI_DOUBLE_PRECISION,MPI_DOUBLE_PRECISION, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_r8
  !
  subroutine write_output_i6(istep,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                   :: istep, nd
    character(len=4), intent(in)                   :: varname
    integer(4)      , intent(in), dimension(1:,1:) :: var
    !
    integer                                        :: offset, fh, ndtot,nelem
    integer, save                                  :: subarray
    character(7)                                   :: fldnum
    integer(kind=MPI_OFFSET_KIND)                  :: filesize,nreals,nreals_myid, displacement,disp
    character(len=1024)                            :: fname
    ! 
    nelem = 6
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)    
    call MPI_FILE_SET_VIEW(fh,disp*nelem*4, MPI_INTEGER, MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE_all(fh,var,nelem*nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i6
  !
  subroutine write_output_i1(istep,varname,var,nd)
    !
    implicit none
    !
    integer,          intent(in)                :: istep, nd
    character(len=4), intent(in)                :: varname
    integer(4),       intent(in), dimension(1:) :: var
    !
    integer                                     :: offset, fh, nelem
    character(7)                                :: fldnum
    integer(kind=MPI_OFFSET_KIND)               :: filesize,disp,nreals,nreals_myid
    character(len=1024)                         :: fname
    ! 
    nelem = 1
    write(fldnum,'(i7.7)') istep
    disp = 0
    !fname = 'data/post/tagging/'//varname//'fld'//fldnum//'.bin'
    fname = trim(datadir_ta)//varname//'fld'//fldnum//'.bin'
    call MPI_SCAN(nd,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    disp = disp-nd
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*4*nelem, MPI_INTEGER,MPI_INTEGER, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,var,nd,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    return
  end subroutine write_output_i1
  !
  !subroutine updthalo(n,idir,p)
  !  !
  !  implicit none
  !  !
  !  integer   , intent(in   ), dimension(2)        :: n
  !  integer   , intent(in   )                      :: idir
  !  integer(8), intent(inout), dimension(0:,0:,0:) :: p
  !  !
  !  !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
  !  !
  !  !  this subroutine updates the halos that store info
  !  !  from the neighboring computational sub-domain
  !  !
  !  select case(idir)
  !  case(1) ! x direction
  !    call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
  !                      p(n(1)+1,0,0),1,xhalo,right,0, &
  !                      comm_cart,status,ierr)
  !    call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
  !                      p(0   ,0,0),1,xhalo,left ,0, &
  !                      comm_cart,status,ierr)
  !       !call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
  !       !               comm_cart,requests(2),error)
  !       !call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
  !       !               comm_cart,requests(1),error)
  !       !call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
  !       !               comm_cart,requests(4),error)
  !       !call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
  !       !               comm_cart,requests(3),error)
  !       !call MPI_WAITALL(4, requests, statuses, error)
  !  case(2) ! y direction
  !    call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
  !                      p(0,n(2)+1,0),1,yhalo,back ,0, &
  !                      comm_cart,status,ierr)
  !    call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
  !                      p(0,0   ,0),1,yhalo,front,0, &
  !                      comm_cart,status,ierr)
  !       !call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
  !       !               comm_cart,requests(1),error)
  !       !call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
  !       !               comm_cart,requests(2),error)
  !       !call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
  !       !               comm_cart,requests(3),error)
  !       !call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
  !       !               comm_cart,requests(4),error)
  !       !call MPI_WAITALL(4, requests, statuses, error)
  !  end select
  !  return
  !end subroutine updthalo
  !

subroutine updthalo(n,idir,p)
  implicit none
  integer , dimension(2), intent(in) :: n
  integer , intent(in) :: idir
  integer(8), dimension(0:,0:,0:), intent(inout) :: p
#ifdef USE_CUDA
  attributes(managed) :: p
  integer :: istat
#endif
  integer :: i,j,k,n_1,n_2
  !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
  !
  !  this subroutine updates the halos that store info
  !  from the neighboring computational sub-domain
  !
  select case(idir)
  case(1) ! x direction
    !if( .false. ) then
    if(dims(1) .eq.  1) then
#ifdef USE_CUDA
      n_1=n(1)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1 ,j,k) = p(  1,j,k)
          p(0     ,j,k) = p(n_1,j,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(n(1)+1,:,:) = p(   1,:,:)
      p(0     ,:,:) = p(n(1),:,:)
      !$OMP END WORKSHARE
#endif
    else
      n_1=n(1)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          xsl_buf(j,k) = p(  1,j,k)
          xsr_buf(j,k) = p(n_1,j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(xsl_buf(0,0), size( xsl_buf ),MPI_REAL8,left ,0, &
        xrr_buf(0,0), size( xrr_buf ),MPI_REAL8,right,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(xsr_buf(0,0), size( xsr_buf ),MPI_REAL8,right,0, &
        xrl_buf(0,0), size( xrl_buf ),MPI_REAL8,left ,0, &
        comm_cart,status,ierr)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do j=lbound(p,2),ubound(p,2)
          p(n_1+1,j,k) = xrr_buf(j,k)
          p(0    ,j,k) = xrl_buf(j,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      !call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
      !                  p(n(1)+1,0,0),1,xhalo,right,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
      !                  p(0   ,0,0),1,xhalo,left ,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  case(2) ! y direction
    !if( .false. ) then
    if( dims(2) .eq.  1 ) then
#ifdef USE_CUDA
      n_2=n(2)
      !$cuf kernel do(2) <<<*,*>>>
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = p(i,  1,k)
          p(i,    0,k) = p(i,n_2,k)
        enddo
      enddo
#else
      !$OMP WORKSHARE
      p(:,n(2)+1,:)  = p(:, 1  ,:)
      p(:, 0     ,:) = p(:,n(2),:)
      !$OMP END WORKSHARE
#endif
    else
      n_2=n(2)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          ysl_buf(i,k) = p(i,  1,k)
          ysr_buf(i,k) = p(i,n_2,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !
      call MPI_SENDRECV(ysl_buf(0,0), size( ysl_buf ),MPI_REAL8,front,0, &
        yrr_buf(0,0), size( yrr_buf ),MPI_REAL8,back ,0, &
        comm_cart,status,ierr)
      call MPI_SENDRECV(ysr_buf(0,0), size( ysr_buf ),MPI_REAL8,back ,0, &
        yrl_buf(0,0), size( yrl_buf ),MPI_REAL8,front,0, &
        comm_cart,status,ierr)
#ifdef USE_CUDA
      !$cuf kernel do(2) <<<*,*>>>
#endif
      do k=lbound(p,3),ubound(p,3)
        do i=lbound(p,1),ubound(p,1)
          p(i,n_2+1,k) = yrr_buf(i,k)
          p(i,    0,k) = yrl_buf(i,k)
        enddo
      enddo
      !@cuf istat = cudaDeviceSynchronize()
      !call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
      !                  p(0,n(2)+1,0),1,yhalo,back ,0, &
      !                  comm_cart,status,ierr)
      !call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
      !                  p(0,0   ,0),1,yhalo,front,0, &
      !                  comm_cart,status,ierr)
      !!call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
      !!               comm_cart,requests(1),error)
      !!call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
      !!               comm_cart,requests(2),error)
      !!call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
      !!               comm_cart,requests(3),error)
      !!call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
      !!               comm_cart,requests(4),error)
      !!call MPI_WAITALL(4, requests, statuses, error)
    endif
  end select

  return
end subroutine updthalo
end module mod_tagging
