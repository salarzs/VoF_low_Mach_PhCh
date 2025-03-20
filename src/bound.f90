module mod_bound
  !
  use mpi
  use mod_common_mpi, only: ierr,status,comm_cart,left,right,front,back,xhalo,yhalo,xhalo_3p,yhalo_3p,coord     
  !
  implicit none
  !
  private
  public bounduvw,bounduvw_b,boundp,boundsb,updt_rhs_b
  !
  contains
  !
  subroutine bounduvw(cbc,n,bc,isoutflow,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions on smaller stencil
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3,3)     :: cbc
    integer         , intent(in   ), dimension(3)           :: n 
    real(8)         , intent(in   ), dimension(0:1,3,3)     :: bc
    logical         , intent(in   ), dimension(0:1,3)       :: isoutflow
    real(8)         , intent(in   ), dimension(3)           :: dl
    real(8)         , intent(in   ), dimension(-2:)         :: dzc,dzf
    real(8)         , intent(inout), dimension( 0:, 0:, 0:) :: u,v,w
    !
    integer :: q,idir,sgn,ioutflowdir,qmin
    !
    qmin = abs(lbound(u,1)) ! absolute module of the minimum starting index of u,v,w
    !
    call updthalo((/n(1),n(2)/),1,u)
    call updthalo((/n(1),n(2)/),2,u)
    call updthalo((/n(1),n(2)/),1,v)
    call updthalo((/n(1),n(2)/),2,v)
    call updthalo((/n(1),n(2)/),1,w)
    call updthalo((/n(1),n(2)/),2,w)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1,1),0,n(1),1,.false.,bc(0,1,1),dl(1),u)
      call set_bc(cbc(0,1,2),0,n(1),1,.true. ,bc(0,1,2),dl(1),v)
      call set_bc(cbc(0,1,3),0,n(1),1,.true. ,bc(0,1,3),dl(1),w)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1,1),1,n(1),1,.false.,bc(1,1,1),dl(1),u)
      call set_bc(cbc(1,1,2),1,n(1),1,.true. ,bc(1,1,2),dl(1),v)
      call set_bc(cbc(1,1,3),1,n(1),1,.true. ,bc(1,1,3),dl(1),w)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2,1),0,n(2),2,.true. ,bc(0,2,1),dl(2),u)
      call set_bc(cbc(0,2,2),0,n(2),2,.false.,bc(0,2,2),dl(2),v)
      call set_bc(cbc(0,2,3),0,n(2),2,.true. ,bc(0,2,3),dl(2),w)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2,1),1,n(2),2,.true. ,bc(1,2,1),dl(2),u)
      call set_bc(cbc(1,2,2),1,n(2),2,.false.,bc(1,2,2),dl(2),v)
      call set_bc(cbc(1,2,3),1,n(2),2,.true. ,bc(1,2,3),dl(2),w)
    endif
    call set_bc(cbc(0,3,1),0,n(3),3,.true. ,bc(0,3,1),dzc(0)   ,u)
    call set_bc(cbc(0,3,2),0,n(3),3,.true. ,bc(0,3,2),dzc(0)   ,v)
    call set_bc(cbc(0,3,3),0,n(3),3,.false.,bc(0,3,3),dzf(0)   ,w)
    call set_bc(cbc(1,3,1),1,n(3),3,.true. ,bc(1,3,1),dzc(n(3)),u)
    call set_bc(cbc(1,3,2),1,n(3),3,.true. ,bc(1,3,2),dzc(n(3)),v)
    call set_bc(cbc(1,3,3),1,n(3),3,.false.,bc(1,3,3),dzf(n(3)),w)
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n,ioutflowdir,qmin,dl,dzf,u,v,w)
        endif
      enddo
    enddo
    !
    return
  end subroutine bounduvw
  !
  subroutine bounduvw_b(cbc,n,bc,isoutflow,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions on a larger stencil
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3,3)     :: cbc
    integer         , intent(in   ), dimension(3)           :: n 
    real(8)         , intent(in   ), dimension(0:1,3,3)     :: bc
    logical         , intent(in   ), dimension(0:1,3)       :: isoutflow
    real(8)         , intent(in   ), dimension(3)           :: dl
    real(8)         , intent(in   ), dimension(-2:)         :: dzc,dzf
    real(8)         , intent(inout), dimension(-2:,-2:,-2:) :: u,v,w
    !
    integer :: q,idir,sgn,ioutflowdir,qmin
    !
    qmin = abs(lbound(u,1)) ! absolute module of the minium starting index of u,v,w
    !
    call updthalo_3p((/n(1),n(2)/),1,u)
    call updthalo_3p((/n(1),n(2)/),2,u)
    call updthalo_3p((/n(1),n(2)/),1,v)
    call updthalo_3p((/n(1),n(2)/),2,v)
    call updthalo_3p((/n(1),n(2)/),1,w)
    call updthalo_3p((/n(1),n(2)/),2,w)
    !
    do q=0,qmin
    if(left.eq.MPI_PROC_NULL) then
      select case(cbc(0,1,1)) ! u in x,0
      case('D')
        u(0-q     ,:,:) = 0.d0
      case('N')
        u(0-q     ,:,:) = u(1+q   ,:,:) 
      end select
    endif
    if(right.eq.MPI_PROC_NULL) then
      select case(cbc(1,1,1)) ! u in x,nx+1
      case('D')
        !u(n(1)+1,:,:) = u(n(1)-1,:,:) ! not used
        u(n(1)+q  ,:,:) = 0.d0
      case('N')
        u(n(1)+1+q,:,:) = u(n(1)-q,:,:)
      end select
    endif
    if(front.eq.MPI_PROC_NULL) then
      select case(cbc(0,2,1)) ! u in y,0
      case('N')
        u(:,0-q     ,:) =  u(:,1+q   ,:)
      case('D')
        u(:,0-q     ,:) = -u(:,1+q   ,:)
      end select
    endif
    if(back .eq.MPI_PROC_NULL) then
      select case(cbc(1,2,1)) ! u in y,ny+1
      case('N')
        u(:,n(2)+1+q,:) =  u(:,n(2)-q,:)
      case('D')
        u(:,n(2)+1+q,:) = -u(:,n(2)-q,:)
      end select
    endif
    select case(cbc(0,3,1)) ! u in z,0
    case('P')
      u(:,:,0-q     ) =  u(:,:,n(3)-q)
    case('N')
      u(:,:,0-q     ) =  u(:,:,1+q   )
    case('D')
      u(:,:,0-q     ) = -u(:,:,1+q   )
    end select
    select case(cbc(1,3,1)) ! u in z,nz+1
    case('P')
      u(:,:,n(3)+1+q) =  u(:,:,1+q   )
    case('N')
      u(:,:,n(3)+1+q) =  u(:,:,n(3)-q)
    case('D')
      u(:,:,n(3)+1+q) = -u(:,:,n(3)-q)
    end select
    if(left.eq.MPI_PROC_NULL) then
      select case(cbc(0,1,2)) ! v in x,0
      case('N')
        v(0-q     ,:,:) =  v(1+q   ,:,:)
      case('D')
        v(0-q     ,:,:) = -v(1+q   ,:,:)
      end select
    endif
    if(right.eq.MPI_PROC_NULL) then
      select case(cbc(1,1,2)) ! v in x,nx+1
      case('N')
        v(n(1)+1+q,:,:) =  v(n(1)-q,:,:)
      case('D')
        v(n(1)+1+q,:,:) = -v(n(1)-q,:,:)
      end select
    endif
    if(front.eq.MPI_PROC_NULL) then
      select case(cbc(0,2,2)) ! v in y,0
      case('D')
        v(:,0-q     ,:) = 0.d0
      case('N')
        v(:,0-q     ,:) = v(:,1+q   ,:)
        !v(:,0-q     ,:) = v(:,0+q   ,:)
        !if(q.eq.0) v(:,0,:)=v(:,1,:)
      end select
    endif
    if(back .eq.MPI_PROC_NULL) then
      select case(cbc(1,2,2)) ! v in y,ny+1
      case('D')
        v(:,n(2)+1,:) = v(:,n(2)-1,:) ! not used
        !v(:,n(2)+q  ,:) = 0.d0
      case('N')
        v(:,n(2)+1+q,:) = v(:,n(2)-q,:)
        !v(:,n(2)+q,:) = v(:,n(2)-q,:)
        !if(q.eq.0) v(:,n(2),:)=v(:,n(2)-1,:)
      end select
    endif
    select case(cbc(0,3,2)) ! v in z,0
    case('P')
      v(:,:,0-q     ) =  v(:,:,n(3)-q)
    case('N')
      v(:,:,0-q     ) =  v(:,:,1+q   )
    case('D')
      v(:,:,0-q     ) = -v(:,:,1+q   )
    end select
    select case(cbc(1,3,2)) ! v in z,nz+1
    case('P')
      v(:,:,n(3)+1+q) =  v(:,:,1+q   )
    case('N')
      v(:,:,n(3)+1+q) =  v(:,:,n(3)-q)
    case('D')
      v(:,:,n(3)+1+q) = -v(:,:,n(3)-q)
    end select
    if(left.eq.MPI_PROC_NULL) then
      select case(cbc(0,1,3)) ! w in x,0
      case('N')
        w(0-q     ,:,:) =  w(1+q   ,:,:)
      case('D')
        w(0-q     ,:,:) = -w(1+q   ,:,:)
      end select
    endif
    if(right.eq.MPI_PROC_NULL) then
      select case(cbc(1,1,3)) ! w in x,nx+1
      case('N')
        w(n(1)+1+q,:,:) =  w(n(1)-q,:,:)
      case('D')
        w(n(1)+1+q,:,:) = -w(n(1)-q,:,:)
      end select
    endif
    if(front.eq.MPI_PROC_NULL) then
      select case(cbc(0,2,3)) ! w in y,0
      case('N')
        w(:,0-q     ,:) =  w(:,1+q   ,:)
      case('D')
        w(:,0-q     ,:) = -w(:,1+q   ,:)
      end select
    endif
    if(back .eq.MPI_PROC_NULL) then
      select case(cbc(1,2,3)) ! w in y,ny+1
      case('N')
        w(:,n(2)+1+q,:) =  w(:,n(2)-q,:)
      case('D')
        w(:,n(2)+1+q,:) = -w(:,n(2)-q,:)
      end select
    endif
    select case(cbc(0,3,3)) ! w in z,0
    case('P')
      w(:,:,0-q     ) =  w(:,:,n(3)-q)
    case('D')
      w(:,:,0-q     ) = 0.d0
    case('N')
      w(:,:,0-q     ) = w(:,:,1+q   )
      !w(:,:,0-q     ) =  w(:,:,0+q   )
      !if(q.eq.0) w(:,:,0)=w(:,:,1)
    end select
    select case(cbc(1,3,3)) ! w in z,nz+1
    case('P')
      w(:,:,n(3)+1+q) = w(:,:,1+q   )
    case('D')
      !w(:,:,n(3)+1) = w(:,:,n(3)-1) ! not used
      w(:,:,n(3)+q  ) = 0.d0
    case('N')
      w(:,:,n(3)+1+q) = w(:,:,n(3)-q)
      !w(:,:,n(3)+q) =  w(:,:,n(3)-q)
      !if(q.eq.0) w(:,:,n(3))=w(:,:,n(3)-1)
    end select
    enddo
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n,ioutflowdir,qmin,dl,dzf,u,v,w)
        endif
      enddo
    enddo
    !
    return
  end subroutine bounduvw_b
  !
  subroutine boundp(cbc,n,bc,dl,dzc,dzf,p)
    !
    ! imposes pressure and scalar (on smaller stencil) boundary conditions
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)    :: cbc
    integer         , intent(in   ), dimension(3)        :: n 
    real(8)         , intent(in   ), dimension(0:1,3)    :: bc
    real(8)         , intent(in   ), dimension(3)        :: dl
    real(8)         , intent(in   ), dimension(-2:)      :: dzc,dzf
    real(8)         , intent(inout), dimension(0:,0:,0:) :: p
    !
    call updthalo((/n(1),n(2)/),1,p)
    call updthalo((/n(1),n(2)/),2,p)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1),0,n(1),1,.true.,bc(0,1),dl(1),p)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1),1,n(1),1,.true.,bc(1,1),dl(1),p)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2),0,n(2),2,.true.,bc(0,2),dl(2),p)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2),1,n(2),2,.true.,bc(1,2),dl(2),p)
    endif
    call set_bc(cbc(0,3),0,n(3),3,.true.,bc(0,3),dzc(0)   ,p)
    call set_bc(cbc(1,3),1,n(3),3,.true.,bc(1,3),dzc(n(3)),p)
    !
    return
  end subroutine boundp
  !
  subroutine boundsb(cbc,n,bc,dl,dzc,dzf,s)
    !
    ! imposes scalar boundary conditions (defined on a larger stencil)
    !
    implicit none
    !
    character(len=1), intent(in   ), dimension(0:1,3)       :: cbc
    integer         , intent(in   ), dimension(3)           :: n 
    real(8)         , intent(in   ), dimension(0:1,3)       :: bc
    real(8)         , intent(in   ), dimension(3)           :: dl
    real(8)         , intent(in   ), dimension(-2:)         :: dzc,dzf
    real(8)         , intent(inout), dimension(-2:,-2:,-2:) :: s
    !
    integer :: q
    !
    call updthalo_3p((/n(1),n(2)/),1,s)
    call updthalo_3p((/n(1),n(2)/),2,s)
    !
    do q=0,2
    if(left.eq.MPI_PROC_NULL) then
      select case(cbc(0,1)) ! in x,0
      case('N')
        s(0-q   ,:,:) =  s(1   +q,:,:)
      case('D')
        s(0-q   ,:,:) = -s(1   +q,:,:) + 2.d0*bc(0,1)
      case('L')
        s(0-q   ,:,:) = 2.d0*s(1-q  ,:,:) - s(2-q   ,:,:)
      end select
    endif
    if(right.eq.MPI_PROC_NULL) then
      select case(cbc(1,1)) ! in x,nx
      case('N')
        s(n(1)+1+q,:,:) =  s(n(1)-q,:,:)
      case('D')
        s(n(1)+1+q,:,:) = -s(n(1)-q,:,:) + 2.d0*bc(1,1)
      case('L')
        s(n(1)+1+q,:,:) = 2.d0*s(n(1)+q,:,:) - s(n(1)-1+q,:,:)
      end select
    endif
    if(front.eq.MPI_PROC_NULL) then
      select case(cbc(0,2)) ! in y,0
      case('N')
        s(:,0-q   ,:) =  s(:,1+q ,:)
      case('D')
        s(:,0-q   ,:) = -s(:,1+q ,:) + 2.d0*bc(0,2)
      case('L')
        s(:,0-q   ,:) = 2.d0*s(:,1-q ,:) - s(:,2-q   ,:)
      end select
    endif
    if(back .eq.MPI_PROC_NULL) then
      select case(cbc(1,2)) ! in y,ny
      case('N')
        s(:,n(2)+1+q,:) =  s(:,n(2)-q,:)
      case('D')
        s(:,n(2)+1+q,:) = -s(:,n(2)-q,:) + 2.d0*bc(1,2)
      case('L')
        s(:,n(2)+1+q,:) = 2.d0*s(:,n(2)+q,:) - s(:,n(2)-1+q,:)
      end select
    endif
    select case(cbc(0,3)) ! in z,0
    case('P')
      s(:,:,0-q   ) =  s(:,:,n(3)-q)
    case('N')
      s(:,:,0-q   ) =  s(:,:,1+q )
    case('D')
      s(:,:,0-q   ) = -s(:,:,1+q ) + 2.d0*bc(0,3)
    case('L')
      s(:,:,0-q   ) = 2.d0*s(:,:,1-q ) - s(:,:,2-q   )
    end select
    select case(cbc(1,3)) ! in z,nz
    case('P')
      s(:,:,n(3)+1+q) =  s(:,:,1+q )
    case('N')
      s(:,:,n(3)+1+q) =  s(:,:,n(3)-q)
    case('D')
      s(:,:,n(3)+1+q) = -s(:,:,n(3)-q) + 2.d0*bc(1,3)
    case('L')
      s(:,:,n(3)+1+q) = 2.d0*s(:,:,n(3)+q) - s(:,:,n(3)-1+q)
    end select
    enddo
    !
    return
  end subroutine boundsb
  !
  subroutine set_bc(ctype,ibound,n,idir,centered,rvalue,dr,p)
    !
    implicit none
    !
    character(len=1), intent(in   )                      :: ctype
    integer         , intent(in   )                      :: ibound,n,idir
    logical         , intent(in   )                      :: centered
    real(8)         , intent(in   )                      :: rvalue,dr
    real(8)         , intent(inout), dimension(0:,0:,0:) :: p
    !
    real(8) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and.centered) then
      factor = 2.d0*factor
      sgn    = -1.d0
    endif
    if(ctype.eq.'N'.and.centered) then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1.d0
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !p(0  ,:,:) = p(n,:,:)
        !p(n+1,:,:) = p(1,:,:)
      case(2)
        !p(:,0  ,:) = p(:,n,:)
        !p(:,n+1,:) = p(:,1,:)
      case(3)
        !$OMP WORKSHARE
        p(:,:,0  ) = p(:,:,n)
        p(:,:,n+1) = p(:,:,1)
        !$OMP END WORKSHARE
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0  ,:,:) = factor+sgn*p(1,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n+1,:,:) = factor+sgn*p(n,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0  ,:) = factor+sgn*p(:,1,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n+1,:) = factor+sgn*p(:,n,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0  ) = factor+sgn*p(:,:,1)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n+1) = factor+sgn*p(:,:,n)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0,:,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n  ,:,:) = factor
            p(n+1,:,:) = p(n-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n  ,:) = factor
            p(:,n+1,:) = p(:,n-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n  ) = factor
            p(:,:,n+1) = p(:,:,n-1)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0,  :,:) = 1.d0*factor + p(1  ,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n  ,:,:) = 1.d0*factor + p(n-1,:,:)
            p(n+1,:,:) = 2.d0*factor + p(n-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0  ,:) = 1.d0*factor + p(:,1  ,:) 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n  ,:) = 1.d0*factor + p(:,n-1,:)
            p(:,n+1,:) = 2.d0*factor + p(:,n-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0  ) = 1.d0*factor + p(:,:,1  )
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n  ) = 1.d0*factor + p(:,:,n-1)
            p(:,:,n+1) = 2.d0*factor + p(:,:,n-1)
            !$OMP END WORKSHARE
          endif
        end select
      endif
    end select
    !
    return
  end subroutine set_bc
  !
  subroutine outflow(n,idir,qmin,dl,dzf,u,v,w)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)                    :: n
    integer, intent(in   )                                  :: idir
    integer, intent(in   )                                  :: qmin
    real(8), intent(in   ), dimension(3)                    :: dl
    real(8), intent(in   ), dimension(-2:)                  :: dzf
    real(8), intent(inout), dimension(-qmin:,-qmin:,-qmin:) :: u,v,w
    !
    real(8) :: dx,dy,dxi,dyi
    real(8), dimension(-2:n(3)+3) :: dzfi
    integer :: i,j,k,q
    !
    dx   = dl(1)     
    dxi  = dl(1)**(-1)
    dy   = dl(2)     
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    ! Note: here we assume that in outflow boundaries, 
    !       there is not phase change or any other volume sources.
    !
    select case(idir)
    case(1) ! x direction, right
      if(right.eq.MPI_PROC_NULL) then
        !i = n(1) + 1
        i = n(1) + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(n,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,n(3)
          do j=1,n(2)
            do q=0,qmin+1 ! note that the value at u(i+qmin+1,j,k) is not used
              u(i+q,j,k) = u(i-1-q,j,k) - dx*((v(i+q,j,k)-v(i+q,j-1,k))*dyi+(w(i+q,j,k)-w(i+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    case(2) ! y direction, back
      if(back.eq.MPI_PROC_NULL) then
        !j = n(2) + 1
        j = n(2) + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(n,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,n(3)
          do i=1,n(1)
            !do q=0,qmin
            do q=0,qmin+1 ! note that the value at v(i,j+qmin+1,k) is not used
              v(i,j+q,k) = v(i,j-1-q,k) - dy*((u(i,j+q,k)-u(i-1,j+q,k))*dxi+(w(i,j+q,k)-w(i,j+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(3) ! z direction, top
      !k = n(3) + 1
      k = n(3) + 0
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,q) &
      !$OMP SHARED(n,k,qmin,u,v,w,dzf,dxi,dyi)
      do j=1,n(2)
        do i=1,n(1)
          do q=0,qmin+1 ! note that the value at w(i,j,k+qmin+1) is not used
            w(i,j,k+q) = w(i,j,k-1-q) - dzf(k+q)*((u(i,j,k+q)-u(i-1,j,k+q))*dxi+(v(i,j,k+q)-v(i,j-1,k+q))*dyi)
          enddo
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        i = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k,q) &
        !$OMP SHARED(n,i,qmin,u,v,w,dx,dyi,dzfi)
        do k=1,n(3)
          do j=1,n(2)
            do q=0,qmin
              u(i-q,j,k) = u(i+1+q,j,k) + dx*((v(i+1+q,j,k)-v(i+1+q,j-1,k))*dyi+(w(i+1+q,j,k)-w(i+1+q,j,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-2) ! y direction, front
      if(front.eq.MPI_PROC_NULL) then
        j = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k,q) &
        !$OMP SHARED(n,j,qmin,u,v,w,dy,dxi,dzfi)
        do k=1,n(3)
          do i=1,n(1)
            do q=0,qmin
              v(i,j-q,k) = v(i,j+1+q,k) + dy*((u(i,j+1+q,k)-u(i-1,j+1+q,k))*dxi+(w(i,j+1+q,k)-w(i,j+1+q,k-1))*dzfi(k))
            enddo
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-3) ! z direction, bottom
      k = 0
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j,q) &
      !$OMP SHARED(n,k,qmin,u,v,w,dzf,dxi,dyi)
      do j=1,n(2)
        do i=1,n(1)
          do q=0,qmin
            w(i,j,k-q) = w(i,j,k+1+q) + dzf(k-q)*((u(i,j,k+1+q)-u(i-1,j,k+1+q))*dxi+(v(i,j,k+1+q)-v(i,j-1,k+1+q))*dyi)
          enddo
        enddo
      enddo 
      !$OMP END PARALLEL DO
    end select
    !
    return
  end subroutine outflow
  !
  subroutine inflow(n,idir,qmin,dl,dzf,vel2d,u,v,w)
    !
    implicit none
    !
    integer, intent(in   ), dimension(3)                 :: n
    integer, intent(in   )                               :: idir
    integer, intent(in   )                               :: qmin
    real(8), intent(in   ), dimension(3)                 :: dl
    real(8), intent(in   ), dimension(-2:)               :: dzf
    real(8), intent(in   ), dimension(0:,0:)             :: vel2d
    real(8), intent(inout), dimension(-qmin,-qmin,-qmin) :: u,v,w
    !
    real(8) :: dx,dy,dxi,dyi
    real(8), dimension(-2:n(3)+3) :: dzfi
    integer :: i,j,k,q
    !
    dx   = dl(1)
    dxi  = dl(1)**(-1)
    dy   = dl(2)
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    !
    select case(idir)
      case(1) ! x direction
        if(left.eq.MPI_PROC_NULL) then
          i = 0
          do k=1,n(3)
            do j=1,n(2)
              do q=0,qmin
                u(i-q,j,k) = vel2d(j,k)
              enddo
            enddo
          enddo 
        endif
      case(2) ! y direction
        j = 0
        if(front.eq.MPI_PROC_NULL) then
          do k=1,n(3)
            do i=1,n(1)
              do q=0,qmin
                v(i,j-q,k) = vel2d(i,k)
              enddo
            enddo
          enddo 
        endif
      case(3) ! z direction
        k = 0
        do j=1,n(2)
          do i=1,n(1)
            do q=0,qmin
              w(i,j,k-q) = vel2d(i,j)
            enddo
          enddo
        enddo 
    end select
    !
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,rhsbx,rhsby,rhsbz,p)
    !
    implicit none
    !
    character       , intent(in   ), dimension(3)        :: c_or_f
    character(len=1), intent(in   ), dimension(0:1,3)    :: cbc
    integer         , intent(in   ), dimension(3)        :: n
    real(8)         , intent(in   ), dimension(:,:,0:)   :: rhsbx,rhsby,rhsbz
    real(8)         , intent(inout), dimension(0:,0:,0:) :: p
    !
    integer, dimension(3) :: q
    integer :: idir
    q(:) = 0
    !
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(left.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1        ,1:n(2),1:n(3)) = p(1        ,1:n(2),1:n(3)) + rhsbx(:,:,0)
      !$OMP END WORKSHARE
    endif  
    if(right.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(n(1)-q(1),1:n(2),1:n(3)) = p(n(1)-q(1),1:n(2),1:n(3)) + rhsbx(:,:,1)
      !$OMP END WORKSHARE
    endif
    if(front.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),1        ,1:n(3)) = p(1:n(1),1        ,1:n(3)) + rhsby(:,:,0)
      !$OMP END WORKSHARE
    endif
    if(back.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),n(2)-q(2),1:n(3)) = p(1:n(1),n(2)-q(2),1:n(3)) + rhsby(:,:,1)
      !$OMP END WORKSHARE
    endif
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1        ) = p(1:n(1),1:n(2),1        ) + rhsbz(:,:,0)
    p(1:n(1),1:n(2),n(3)-q(3)) = p(1:n(1),1:n(2),n(3)-q(3)) + rhsbz(:,:,1)
    !$OMP END WORKSHARE
    !
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(n,idir,p)
    !
    implicit none
    !
    integer, intent(in   ), dimension(2)        :: n
    integer, intent(in   )                      :: idir
    real(8), intent(inout), dimension(0:,0:,0:) :: p
    !
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
                        p(n(1)+1,0,0),1,xhalo,right,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
                        p(0   ,0,0),1,xhalo,left ,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
         !               comm_cart,requests(2),error)
         !call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
         !               comm_cart,requests(1),error)
         !call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
         !               comm_cart,requests(4),error)
         !call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
         !               comm_cart,requests(3),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
                        p(0,n(2)+1,0),1,yhalo,back ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
                        p(0,0   ,0),1,yhalo,front,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
         !               comm_cart,requests(1),error)
         !call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
         !               comm_cart,requests(2),error)
         !call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
         !               comm_cart,requests(3),error)
         !call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
         !               comm_cart,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    end select
    return
  end subroutine updthalo
  !
  subroutine updthalo_3p(n,idir,p)
    !
    implicit none
    !
    integer, intent(in   ), dimension(2)           :: n
    integer, intent(in   )                         :: idir
    real(8), intent(inout), dimension(-2:,-2:,-2:) :: p
    !
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  This subroutine updates the 3-point halo regions 
    !  that store info from the neighboring computational 
    !  sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(1     ,-2,-2),1,xhalo_3p,left ,0, &
                        p(n(1)+1,-2,-2),1,xhalo_3p,right,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(n(1)-2,-2,-2),1,xhalo_3p,right,0, &
                        p(-2    ,-2,-2),1,xhalo_3p,left ,0, &
                        comm_cart,status,ierr)
    case(2) ! y direction
      call MPI_SENDRECV(p(-2,1     ,-2),1,yhalo_3p,front,0, &
                        p(-2,n(2)+1,-2),1,yhalo_3p,back ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(-2,n(2)-2,-2),1,yhalo_3p,back ,0, &
                        p(-2,-2    ,-2),1,yhalo_3p,front,0, &
                        comm_cart,status,ierr)
    end select
    !
    return
  end subroutine updthalo_3p
  !
end module mod_bound
