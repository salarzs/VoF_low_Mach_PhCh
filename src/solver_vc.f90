module mod_solver_vc
  ! 
  use mod_param     , only: dims!,toll_pslm_it,toll_pslm_mg
  use mod_common_mpi, only: coord,myid,ierr,comm_cart
  use mod_thermo    , only: mass_fraction,thermo_d_lg,thermo_rhog
  use decomp_2d
  !
  implicit none
  !
  integer, parameter :: HYPRESolverSMG      = 1, &
                        HYPRESolverPFMG     = 2, &
                        HYPRESolverGMRES    = 3, &
                        HYPRESolverBiCGSTAB = 4
  integer, parameter :: maxiter = 50
  real(8), parameter :: tol = 1e-08, maxError = tol
  !integer(8) ::  HYPRESolverType = HYPRESolverSMG  ! less efficient, but more robust
  integer(8) ::  HYPRESolverType = HYPRESolverPFMG ! more efficient, but less robust
  integer(8) ::  grid,stencil,solver,mat,vec,x,precond,precond_id
  !
  private
  public init_solver_mg,solver_mg,helmholtz_mg,destroy_solver_mg
  !
  contains
  !
  subroutine init_solver_mg(n,dli,dzci,cbc,bc)
    !
    ! initializiation of an Helmholtz/Poisson solver 
    ! for a Poisson equation with variable coefficients
    !
    implicit none
    !
    integer         , intent(in), dimension(3)     :: n
    real(8)         , intent(in), dimension(3)     :: dli
    real(8)         , intent(in), dimension(-2:)   :: dzci
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(8)         , intent(in), dimension(0:1,3) :: bc
    !
    ! hypre solver parameters/variables
    !
    real(8) :: final_res_norm,tol
    integer :: extra,local_size
    integer :: mpi_comm
    !
    integer, parameter :: ndims = 3, nstencil = 7
    integer, dimension(ndims) :: ijlower, ijupper
    integer, dimension(ndims) :: periods
    integer, dimension(ndims,nstencil) :: offsets
    integer, dimension(nstencil) :: stencil_indices
    integer :: nvalues,num_iterations
    integer :: q,qq
    integer, dimension(3) :: ng
    !
    ng(1:2) = n(1:2)*dims(1:2)
    ng(3)   = n(3)
    periods(:)  = 0
    do q=1,3
      do qq=0,1
        select case(cbc(qq,q))
        case('N')
        case('D')
        case('P')
          periods(q) = ng(q)
        end select
      enddo
    enddo
    mpi_comm = comm_cart
    !
    ! initialize matrix
    !
    ! 1 - setup the grid
    !
    !   1.1 - create empty 2D grid object
    !
    call HYPRE_StructGridCreate(mpi_comm,ndims,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    !
    !   1.2 - set grid extents and assemble
    !
    ijlower(1:2) = (coord(:)  )*n(1:2)+1
    ijupper(1:2) = (coord(:)+1)*n(1:2)
    ijlower(3  ) = 1
    ijupper(3  ) = n(3)
    call HYPRE_StructGridSetExtents(grid,ijlower,ijupper,ierr)
    call HYPRE_StructGridAssemble(  grid,ierr)
    !
    ! 2 - setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(ndims,nstencil,stencil,ierr)
    offsets = reshape((/ 0, 0, 0, &
                        -1, 0, 0, &
                         1, 0, 0, &
                         0,-1, 0, &
                         0, 1, 0, &
                         0, 0,-1, &
                         0, 0, 1 /),shape(offsets))
    do q = 1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q),ierr)
    enddo
    !
    ! 3 - create coefficient matrix and rhs vector
    !
    call HYPRE_StructMatrixCreate(mpi_comm,grid,stencil,mat,ierr)
    call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)                ! set the symmetric property of the matrix 
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(mpi_comm,grid,vec,ierr)
    call HYPRE_StructVectorInitialize(vec,ierr)
    call HYPRE_StructVectorCreate(mpi_comm,grid,x,ierr)
    call HYPRE_StructVectorInitialize(x,ierr)
    !
    ! setup solver, and solve
    ! note: this part was taken from the Paris Simulator code
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      call HYPRE_StructSMGCreate(mpi_comm, solver, ierr)
      call HYPRE_StructSMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructSMGSetTol(solver, MaxError, ierr)
      call hypre_structSMGsetLogging(solver, 1, ierr)
      call HYPRE_StructSMGSetPrintLevel(solver,1,ierr) 
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGCreate(mpi_comm, solver, ierr)
      call HYPRE_StructPFMGSetMaxIter(solver, maxiter, ierr)
      call HYPRE_StructPFMGSetTol(solver, MaxError, ierr)
      call HYPRE_structPFMGsetLogging(solver, 1, ierr)
      call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr) 
      call HYPRE_StructPFMGSetRelChange(solver, 1, ierr) 
      ! Relaxiation Method: 2 is the fastest if symm matrix 
      ! 0: Jacobi
      ! 1: Weighted Jacobi (default)
      ! 2: Red/Black Gauss-Seidel (symmetric: RB pre- and post-relaxation)
      ! 3: Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
      call HYPRE_StructPFMGSetRelaxType(solver,1, ierr) 
      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverGMRES .or. & 
             HYPRESolverType .eq. HYPRESolverBiCGSTAB   ) then
      if (HYPRESolverType .eq. HYPRESolverGMRES) then 
        call HYPRE_StructGMRESCreate(mpi_comm, solver,ierr)
        call HYPRE_StructGMRESSetMaxIter(solver, maxiter,ierr)
        call HYPRE_StructGMRESSetTol(solver, MaxError, ierr)
        !call HYPRE_StructGMRESSetLogging(solver, 1 ,ierr)
      elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
        call HYPRE_StructBiCGSTABCreate(mpi_comm, solver,ierr)
        call HYPRE_StructBiCGSTABSetMaxIter(solver, maxiter,ierr)
        call HYPRE_StructBiCGSTABSetTol(solver, MaxError, ierr)
      endif
      ! Use PFMG as preconditioner
      call HYPRE_StructPFMGCreate(mpi_comm, precond, ierr)
      call HYPRE_StructPFMGSetMaxIter(precond,10, ierr)
      call HYPRE_StructPFMGSetTol(precond,0.0,ierr)
      call HYPRE_StructPFMGSetZeroGuess(precond,ierr)
      call HYPRE_StructPFMGSetRelChange(precond,1,ierr) 
      call HYPRE_StructPFMGSetRelaxType(precond,2,ierr) 
      precond_id = 1   ! Set PFMG as preconditioner
      if (HYPRESolverType .eq. HYPRESolverGMRES) then 
        call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
      elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
        call HYPRE_StructBiCGSTABSetPrecond(solver,precond_id,precond,ierr)
      endif
    endif
    !
    return 
  end subroutine init_solver_mg
  !
  subroutine solver_mg(n,dli,dzci,cbc,bc,alpha,pg,p)
    !
    ! Poisson solver for a Poisson equation with variable
    ! coeffcients
    ! Note: to have an Helmholtz solver just add +1 to cc and 
    !       multiply by c1 = dt
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)        :: n
    real(8)         , intent(in   ), dimension(3)        :: dli
    real(8)         , intent(in   ), dimension(-2:)      :: dzci
    character(len=1), intent(in   ), dimension(0:1,3)    :: cbc
    real(8)         , intent(in   ), dimension(0:1,3)    :: bc
    real(8)         , intent(in   ), dimension(0:,0:,0:) :: alpha,pg
    real(8)         , intent(inout), dimension(0:,0:,0:) :: p
    !
    integer, dimension(3) :: ng
    real(8) :: cxp,cxm,cyp,cym,czp,czm,cc,rhs
    real(8) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
    integer :: i,j,k,ip,jp,kp,im,jm,km,ii,jj,kk,q,qq
    real(8), dimension(0:1,3) :: factor,sgn
    !
    ! hypre solver parameters/variables
    !
    real(8) :: final_res_norm,tol
    integer :: extra,local_size
    integer :: mpi_comm
    !
    integer, parameter :: ndims = 3, nstencil = 7
    integer, dimension(ndims) :: ijlower, ijupper
    integer, dimension(ndims) :: periods
    integer, dimension(ndims,nstencil) :: offsets
    integer, dimension(nstencil) :: stencil_indices
    real(8), allocatable, dimension(:) :: matvalues,vecvalues, &
                                          guessvalues
    integer :: nvalues,num_iterations
    !
    ng(1:2) = n(1:2)*dims(1:2)
    ng(3)   = n(3)
    factor(:,:) = 0.d0
    sgn(   :,:) = 0.d0
    periods(:)  = 0
    !
    ! set of the boundary conditions
    !
    do q=1,3
      do qq=0,1
        !
        select case(cbc(qq,q))
        case('N')
          factor(qq,q) = 1.d0/dli(q)*bc(qq,q)
          sgn(   qq,q) = 1.d0
        case('D')
          factor(qq,q) = -2.d0*bc(qq,q)
          sgn(   qq,q) = -1.d0
        case('P')
          factor(qq,q) = 0.d0
          sgn(   qq,q) = 0.d0
          periods(q)   = ng(q)
        end select
        !
      enddo
    enddo
    mpi_comm = comm_cart
    !
    ijlower(1:2) = (coord(:)  )*n(1:2)+1
    ijupper(1:2) = (coord(:)+1)*n(1:2)
    ijlower(3  ) = 1
    ijupper(3  ) = n(3)
    !
    ! preliminaries for setting up the coefficient matrix
    !
    stencil_indices = (/0,1,2,3,4,5,6/)
    nvalues = product(n(:))*nstencil ! number of grid points*number of stencil entries
    allocate(matvalues(nvalues))
    matvalues(:) = 0.d0
    nvalues = product(n(:))
    allocate(vecvalues(nvalues))
    vecvalues(:) = 0.d0
    allocate(guessvalues(nvalues))
    guessvalues(:) = 0.d0
    !
    ! compute stencil coefficients and rhs
    !
    q = 0
    do k=1,n(3)
      kp = k+1
      km = k-1
      kk = k
      do j=1,n(2)
        jp = j+1
        jm = j-1
        jj = j+coord(2)*n(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          ii = i+coord(1)*n(1)
          !
          q = q+1
          !
#ifdef TWOD
          alphaxp = 0.d0
          alphaxm = 0.d0
#else
          alphaxp = (0.5d0*(alpha(ip,j,k)+alpha(i,j,k)))**(-1)
          alphaxm = (0.5d0*(alpha(im,j,k)+alpha(i,j,k)))**(-1)
#endif
          alphayp = (0.5d0*(alpha(i,jp,k)+alpha(i,j,k)))**(-1)
          alphaym = (0.5d0*(alpha(i,jm,k)+alpha(i,j,k)))**(-1)
          alphazp = (0.5d0*(alpha(i,j,kp)+alpha(i,j,k)))**(-1)
          alphazm = (0.5d0*(alpha(i,j,km)+alpha(i,j,k)))**(-1)
          !
          cxm = alphaxm*dli(1)**2
          cxp = alphaxp*dli(1)**2
          cym = alphaym*dli(2)**2
          cyp = alphayp*dli(2)**2
          czm = alphazm*dli(3)**2 ! change to dzci
          czp = alphazp*dli(3)**2 ! change to dzci
          cc  = -(cxm+cxp+cym+cyp+czm+czp)
          rhs = p(i,j,k)
          !
#ifndef TWOD
          if(periods(1).eq.0) then
            if(    ii.eq.  1) then
              rhs = rhs + cxm*factor(0,1)
              cc  = cc  + sgn(0,1)*cxm
              cxm = 0.d0
            elseif(ii.eq.ng(1)) then
              rhs = rhs + cxp*factor(1,1)
              cc  = cc  + sgn(1,1)*cxp
              cxp = 0.d0
            endif
          endif
#endif
          if(periods(2).eq.0) then
            if(    jj.eq.  1 ) then
              rhs = rhs + cym*factor(0,2)
              cc  = cc  + sgn(0,2)*cym
              cym = 0.d0
            elseif(jj.eq.ng(2)) then
              rhs = rhs + cyp*factor(1,2)
              cc  = cc  + sgn(1,2)*cyp
              cyp = 0.d0
            endif
          endif
          if(periods(3).eq.0) then
            if(    kk.eq.  1) then
              rhs = rhs + czm*factor(0,3)
              cc  = cc  + sgn(0,3)*czm
              czm = 0.d0
            elseif(kk.eq.ng(3)) then
              rhs = rhs + czp*factor(1,3)
              cc  = cc  + sgn(1,3)*czp
              czp = 0.d0
            endif
          endif
          !
          matvalues((q-1)*nstencil+1) = cc
          matvalues((q-1)*nstencil+2) = cxm
          matvalues((q-1)*nstencil+3) = cxp
          matvalues((q-1)*nstencil+4) = cym
          matvalues((q-1)*nstencil+5) = cyp
          matvalues((q-1)*nstencil+6) = czm
          matvalues((q-1)*nstencil+7) = czp
          vecvalues(q               ) = rhs
          guessvalues(q             ) = pg(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    call HYPRE_StructMatrixSetBoxValues(mat,ijlower,ijupper,nstencil, &
                                        stencil_indices,matvalues, &
                                        ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    call HYPRE_StructVectorSetBoxValues(vec,ijlower,ijupper, &
                                        vecvalues,ierr)
    call HYPRE_StructVectorAssemble(vec,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(x,ijlower,ijupper, &
                                        guessvalues,ierr)
    call HYPRE_StructVectorAssemble(x,ierr)
    deallocate(guessvalues)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructSMGSetTol(solver,tol,ierr)
      call HYPRE_StructSMGSetRelChange(solver, 0,ierr)
      call HYPRE_StructSMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructSMGSetNumPostRelax(solver,1,ierr)
      call HYPRE_StructSMGSetLogging(solver, 1,ierr)
      call HYPRE_StructSMGSetup(solver, mat, vec, x,ierr)
      call HYPRE_StructSMGSolve(solver, mat, vec, x,ierr)
      call HYPRE_StructSMGDestroy(solver,ierr)
      !call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
      !call HYPRE_StructSMGSetup(solver, mat, vec, x, ierr)
      !call HYPRE_StructSMGSolve(solver, mat, vec, x, ierr)
      !call HYPRE_StructSMGDestroy(solver,ierr)
      !call HYPRE_StructSMGGetNumIterations(solver, num_iterations,ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructPFMGSetTol(solver,tol,ierr)
      call HYPRE_StructPFMGSetRelChange(solver, 0,ierr)
      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetLogging(solver, 1,ierr)
      call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGDestroy(solver,ierr)
      !call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
      !call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
      !call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
      !call HYPRE_StructPFMGDestroy(solver,ierr)
      !call HYPRE_StructPFMGGetNumIterations(solver, num_iterations,ierr)
    elseif (HYPRESolverType .eq. HYPRESolverGMRES) then 
      call HYPRE_StructGMRESCreate(mpi_comm, solver,ierr)
      call HYPRE_StructGMRESSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESDestroy(solver,ierr)
      !call HYPRE_StructGMRESGetNumIterations(solver, num_iterations,ierr)
    elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
      call HYPRE_StructBiCGSTABCreate(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABDestroy(solver, mat, vec, x, ierr)
      !call HYPRE_StructBiCGSTABGetNumIterations(solver, num_iterations,ierr)
    endif ! HYPRESolverType
    !
    ! end of part based on the Paris Simulator code
    !
    ! fecth results
    !
    call HYPRE_StructVectorGetBoxValues(x,ijlower,ijupper,vecvalues,ierr)
    q = 0
    p(:,:,:) = 0.d0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q = q + 1
          p(i,j,k) = vecvalues(q)
        enddo
      enddo
    enddo
    deallocate(vecvalues,matvalues)
    !
    return 
  end subroutine solver_mg
  !
  subroutine helmholtz_mg(n,dli,dzci,cbc,bc,c1,pth,phi,tmp,sg,alpha1,alpha2,dsdt)
    !
    ! Helmholtz solver for the species mass fraction equation with prescribed Dirichlet BC 
    ! at the vapour/liquid interface Y = Yi, where Yi is obtained from thermodynamic
    ! considerations.
    ! it solves the equation in the domain sign(phi)*abs(phi)<0
    ! see Gibou et al. JCP 2002
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)           :: n
    real(8)         , intent(in   ), dimension(3)           :: dli
    real(8)         , intent(in   ), dimension(-2)          :: dzci
    character(len=1), intent(in   ), dimension(0:1,3)       :: cbc
    real(8)         , intent(in   ), dimension(0:1,3)       :: bc
    real(8)         , intent(in   )                         :: c1
    real(8)         , intent(in   )                         :: pth
    real(8)         , intent(in   ), dimension(-2:,-2:,-2:) :: phi,tmp
    real(8)         , intent(in   ), dimension( 0:, 0:, 0:) :: sg
    real(8)         , intent(in   ), dimension( 0:, 0:, 0:) :: alpha1,alpha2
    real(8)         , intent(inout), dimension(  :,  :,  :) :: dsdt
    !
    real(8) :: theta,si,tmp_i,xx,yy,zz,xi,yi,zi
    integer, dimension(3) :: ng
    real(8), dimension(3) :: dl
    real(8) :: cxp,cxm,cyp,cym,czp,czm,cc,rhs
    real(8) :: alpha1_i,alpha2_i,alpha1_alpha2_g
    real(8) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
    integer :: i,j,k,ip,im,jp,jm,kp,km,ii,jj,kk,q,qq
    real(8), dimension(0:1,3) :: factor,sgn
    integer :: solverid 
    !
    ! hypre solver parameters/variables
    !
    real(8) :: final_res_norm,tol
    integer :: extra,local_size
    integer :: mpi_comm
    !
    integer, parameter :: ndims = 3, nstencil = 7
    integer, dimension(ndims) :: ijlower, ijupper
    integer, dimension(ndims) :: periods
    integer, dimension(ndims,nstencil) :: offsets
    integer, dimension(nstencil) :: stencil_indices
    real(8), allocatable, dimension(:) :: matvalues,vecvalues,guessvalues
    integer :: nvalues,num_iterations
    !
    dl = dli**(-1)
    ng(1:2) = n(1:2)*dims(1:2)
    ng(3)   = n(3)
    factor(:,:) = 0.d0
    sgn(   :,:) = 0.d0
    periods(:)  = 0
    !
    ! set of the boundary conditions
    !
    do q=1,3
      do qq=0,1
        !
        select case(cbc(qq,q))
        case('N')
          factor(qq,q) = 1.d0/dli(q)*bc(qq,q)
          sgn(   qq,q) = 1.d0
        case('D')
          factor(qq,q) = -2.d0*bc(qq,q)
          sgn(   qq,q) = -1.d0
        case('P')
          factor(qq,q) = 0.d0
          sgn(   qq,q) = 0.d0
          periods(q)   = ng(q)
        end select
        !
      enddo
    enddo
    mpi_comm = comm_cart
    !
    ijlower(1:2) = (coord(:)  )*n(1:2)+1
    ijupper(1:2) = (coord(:)+1)*n(1:2)
    ijlower(3  ) = 1
    ijupper(3  ) = n(3)
    !
    ! preliminaries for setting up the coefficient matrix
    !
    stencil_indices = (/0,1,2,3,4,5,6/)
    nvalues = product(n(:))*nstencil ! number of grid points*number of stencil entries
    allocate(matvalues(nvalues))
    matvalues(:)   = 0.d0
    nvalues = product(n(:))
    allocate(vecvalues(nvalues))
    vecvalues(:)   = 0.d0
    allocate(guessvalues(nvalues))
    guessvalues(:) = 0.d0
    !
    ! compute stencil coefficients and rhs
    !
    q = 0
    do k=1,n(3)
      kp = k+1
      km = k-1
      kk = k
      zz = (kk-0.5d0)*dl(3)
      do j=1,n(2)
        jp = j+1
        jm = j-1
        jj = j+coord(2)*n(2)
        yy = (jj-0.5d0)*dl(2)
        do i=1,n(1)
          ip = i+1
          im = i-1
          ii = i+coord(1)*n(1)
          xx = (ii-0.5d0)*dl(1)
          !
          q = q + 1
          if(phi(i,j,k).lt.0.d0) then
            !
            alphaxp = 0.5d0*(alpha1(ip,j,k)*alpha2(ip,j,k)+alpha1(i,j,k)*alpha2(i,j,k))
            alphaxm = 0.5d0*(alpha1(im,j,k)*alpha2(im,j,k)+alpha1(i,j,k)*alpha2(i,j,k))
            alphayp = 0.5d0*(alpha1(i,jp,k)*alpha2(i,jp,k)+alpha1(i,j,k)*alpha2(i,j,k))
            alphaym = 0.5d0*(alpha1(i,jm,k)*alpha2(i,jm,k)+alpha1(i,j,k)*alpha2(i,j,k))
            alphazp = 0.5d0*(alpha1(i,j,kp)*alpha2(i,j,kp)+alpha1(i,j,k)*alpha2(i,j,k))
            alphazm = 0.5d0*(alpha1(i,j,km)*alpha2(i,j,km)+alpha1(i,j,k)*alpha2(i,j,k))
            !
#ifdef TWOD
            cxm = 0.d0
            cxp = 0.d0
#else
            cxm = alphaxm*dli(1)**2*c1/alpha1(i,j,k)
            cxp = alphaxp*dli(1)**2*c1/alpha1(i,j,k)
#endif
            cym = alphaym*dli(2)**2*c1/alpha1(i,j,k)
            cyp = alphayp*dli(2)**2*c1/alpha1(i,j,k)
            czm = alphazm*dli(3)**2*c1/alpha1(i,j,k)
            czp = alphazp*dli(3)**2*c1/alpha1(i,j,k)
            cc  = -(cxm+cxp+cym+cyp+czm+czp) + 1.d0
            rhs = dsdt(i,j,k)
            !
#ifndef TWOD
            if(periods(1).eq.0) then
              if(    ii.eq.  1) then
                rhs = rhs + cxm*factor(0,1)
                cc  = cc  + sgn(0,1)*cxm
                cxm = 0.d0
              elseif(ii.eq.ng(1)) then
                rhs = rhs + cxp*factor(1,1)
                cc  = cc  + sgn(1,1)*cxp
                cxp = 0.d0
              endif
            endif
#endif
            if(periods(2).eq.0) then
              if(    jj.eq.  1 ) then
                rhs = rhs + cym*factor(0,2)
                cc  = cc  + sgn(0,2)*cym
                cym = 0.d0
              elseif(jj.eq.ng(2)) then
                rhs = rhs + cyp*factor(1,2)
                cc  = cc  + sgn(1,2)*cyp
                cyp = 0.d0
              endif
            endif
            if(periods(3).eq.0) then
              if(    kk.eq.  1) then
                rhs = rhs + czm*factor(0,3)
                cc  = cc  + sgn(0,3)*czm
                czm = 0.d0
              elseif(kk.eq.ng(3)) then
                rhs = rhs + czp*factor(1,3)
                cc  = cc  + sgn(1,3)*czp
                czp = 0.d0
              endif
            endif
            !
            ! now we change the stencil discretization in those cells
            ! cut by the inteface, following a dimension by dimension approach
            ! 
#ifndef TWOD
            !
            ! along x
            !
            !if(phi(im,j,k).gt.0.d0.and.(ii.ne.1.and.periods(1).eq.0)) then
            if(phi(im,j,k).gt.0.d0.and.(ii.ne.1.or.periods(1).eq.0)) then
            !if(phi(im,j,k).gt.0.d0.and.(ii.ne.1)) then
            !if(phi(im,j,k).gt.0.d0) then
              xi       = (xx*abs(phi(im,j,k))+(xx-dl(1))*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(im,j,k))+tmp(im,j,k)*abs(phi(i,j,k)))/(abs(phi(im,j,k))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(xx-xi)*dli(1)
              rhs      = rhs
              cc       = cc+cxm ! we remove the coefficient from the diagonal
              if(theta.gt.dl(1)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxm = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(1)**2*c1/alpha1(i,j,k)
                rhs = rhs-cxm*si/theta
                cc  = cc-(cxm/theta)
              endif
              !if(theta.gt.dl(1)) then
              !  alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
              !  cxm = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(1)**2*c1/alpha1(i,j,k)
              !  rhs = rhs - cxm*si/theta
              !  cc  = cc+cxm-(cxm/theta)
              !else 
              !  rhs = rhs
              !  cc  = cc+cxm
              !endif
              cxm = 0.d0
            endif
            !if(phi(ip,j,k).gt.0.d0.and.(ii.ne.ng(1).and.periods(1).eq.0)) then
            if(phi(ip,j,k).gt.0.d0.and.(ii.ne.ng(1).or.periods(1).eq.0)) then
            !if(phi(ip,j,k).gt.0.d0.and.(ii.ne.ng(1))) then
            !if(phi(ip,j,k).gt.0.d0) then
              xi       = (xx*abs(phi(ip,j,k))+(xx+dl(1))*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(ip,j,k))+tmp(ip,j,k)*abs(phi(i,j,k)))/(abs(phi(ip,j,k))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(xx-xi)*dli(1)
              rhs      = rhs
              cc       = cc+cxp ! we remove the coefficient from the diagonal
              if(theta.gt.dl(1)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxp = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(1)**2*c1/alpha1(i,j,k)
                rhs = rhs-cxp*si/theta
                cc  = cc-(cxp/theta)
              endif
              !if(theta.gt.dl(1)) then
              !  rhs = rhs - cxp*si/theta
              !  cc  = cc+cxp-(cxp/theta)
              !else
              !  rhs = rhs
              !  cc  = cc+cxp
              !endif
              cxp = 0.d0
            endif
#endif
            !
            ! along y
            !
            !if(phi(i,jm,k).gt.0.d0.and.(jj.ne.1.and.periods(2).eq.0)) then
            if(phi(i,jm,k).gt.0.d0.and.(jj.ne.1.or.periods(2).eq.0)) then
            !if(phi(i,jm,k).gt.0.d0.and.(jj.ne.1)) then
            !if(phi(i,jm,k).gt.0.d0) then
              yi       = (yy*abs(phi(i,jm,k))+(yy-dl(2))*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(i,jm,k))+tmp(i,jm,k)*abs(phi(i,j,k)))/(abs(phi(i,jm,k))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(yy-yi)*dli(2)
              rhs      = rhs
              cc       = cc+cym ! we remove the coefficient from the diagonal
              if(theta.gt.dl(2)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxp = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(2)**2*c1/alpha1(i,j,k)
                rhs = rhs-cym*si/theta
                cc  = cc-(cym/theta)
              endif
              !if(theta.gt.dl(2)) then
              !  rhs = rhs - cym*si/theta
              !  cc  = cc+cym-(cym/theta)
              !else
              !  rhs = rhs
              !  cc  = cc+cym
              !endif
              cym = 0.d0
            endif
            !if(phi(i,jp,k).gt.0.d0.and.(jj.ne.ng(2).and.periods(2).eq.0)) then
            if(phi(i,jp,k).gt.0.d0.and.(jj.ne.ng(2).or.periods(2).eq.0)) then
            !if(phi(i,jp,k).gt.0.d0.and.(jj.ne.ng(2))) then
            !if(phi(i,jp,k).gt.0.d0) then
              yi       = (yy*abs(phi(i,jp,k))+(yy+dl(2))*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(i,jp,k))+tmp(i,jp,k)*abs(phi(i,j,k)))/(abs(phi(i,jp,k))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(yy-yi)*dli(2)
              rhs      = rhs
              cc       = cc+cyp ! we remove the coefficient from the diagonal
              if(theta.gt.dl(2)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxp = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(2)**2*c1/alpha1(i,j,k)
                rhs = rhs-cyp*si/theta
                cc  = cc-(cyp/theta)
              endif
              !if(theta.gt.dl(2)) then
              !  rhs = rhs - cyp*si/theta
              !  cc  = cc+cyp-(cyp/theta)
              !else
              !  rhs = rhs
              !  cc  = cc+cyp
              !endif
              cyp = 0.d0
            endif
            !
            ! along z
            !
            !if(phi(i,j,km).gt.0.d0.and.(kk.ne.1.and.periods(3).eq.0)) then
            if(phi(i,j,km).gt.0.d0.and.(kk.ne.1.or.periods(3).eq.0)) then
            !if(phi(i,j,km).gt.0.d0.and.(kk.ne.1)) then
            !if(phi(i,j,km).gt.0.d0) then
              zi       = (zz*abs(phi(i,j,km))+(zz-dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(i,j,km))+tmp(i,j,km)*abs(phi(i,j,k)))/(abs(phi(i,j,km))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(zz-zi)*dli(3)
              rhs      = rhs
              cc       = cc+czm ! we remove the coefficient from the diagonal
              if(theta.gt.dl(3)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxp = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(3)**2*c1/alpha1(i,j,k)
                rhs = rhs-czm*si/theta
                cc  = cc-(czm/theta)
              endif
              !if(theta.gt.dl(3)) then
              !  rhs = rhs - czm*si/theta
              !  cc  = cc+czm-(czm/theta)
              !else
              !  rhs = rhs
              !  cc  = cc+czm
              !endif
              czm = 0.d0
            endif
            !if(phi(i,j,kp).gt.0.d0.and.(kk.ne.ng(3).and.periods(3).eq.0)) then
            if(phi(i,j,kp).gt.0.d0.and.(kk.ne.ng(3).or.periods(3).eq.0)) then
            !if(phi(i,j,kp).gt.0.d0.and.(kk.ne.ng(3))) then
            !if(phi(i,j,kp).gt.0.d0) then
              zi       = (zz*abs(phi(i,j,kp))+(zz+dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              tmp_i    = (tmp(i,j,k)*abs(phi(i,j,kp))+tmp(i,j,kp)*abs(phi(i,j,k)))/(abs(phi(i,j,kp))+abs(phi(i,j,k)))
              si       = mass_fraction(pth,tmp_i)
              alpha1_i = thermo_rhog(pth,tmp_i,si)
              alpha2_i = thermo_d_lg(pth,tmp_i,si)
              theta    = abs(zz-zi)*dli(3)
              rhs      = rhs
              cc       = cc+czp ! we remove the coefficient from the diagonal
              if(theta.gt.dl(3)) then
                alpha1_alpha2_g = (alpha1_i*alpha2_i + (theta-1.d0)*alpha1(i,j,k)*alpha2(i,j,k))/theta
                cxp = 0.5d0*(alpha1(i,j,k)*alpha2(i,j,k)+alpha1_alpha2_g)*dli(3)**2*c1/alpha1(i,j,k)
                rhs = rhs-czp*si/theta
                cc  = cc-(czp/theta)
              endif
              !if(theta.gt.dl(3)) then
              !  rhs = rhs - czp*si/theta
              !  cc  = cc+czp-(czp/theta)
              !else
              !  rhs = rhs
              !  cc  = cc+czp
              !endif
              czp = 0.d0
            endif
            !
          else
            !
            cc  = 1.d0
            cxm = 0.d0
            cym = 0.d0
            czm = 0.d0
            cxp = 0.d0
            cyp = 0.d0
            czp = 0.d0
            rhs = dsdt(i,j,k)
            !
          endif
          !
          matvalues( (q-1)*nstencil+1) = cc
          matvalues( (q-1)*nstencil+2) = cxm
          matvalues( (q-1)*nstencil+3) = cxp
          matvalues( (q-1)*nstencil+4) = cym
          matvalues( (q-1)*nstencil+5) = cyp
          matvalues( (q-1)*nstencil+6) = czm
          matvalues( (q-1)*nstencil+7) = czp
          vecvalues  (q              ) = rhs
          guessvalues(q              ) = sg(i,j,k)
          !
        enddo
      enddo
    enddo
    !
    call HYPRE_StructMatrixSetBoxValues(mat,ijlower,ijupper,nstencil, &
                                        stencil_indices,matvalues, &
                                        ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    call HYPRE_StructVectorSetBoxValues(vec,ijlower,ijupper, &
                                        vecvalues,ierr)
    call HYPRE_StructVectorAssemble(vec,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(x,ijlower,ijupper, &
                                        guessvalues,ierr)
    call HYPRE_StructVectorAssemble(x,ierr)
    deallocate(guessvalues)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructSMGSetTol(solver,tol,ierr)
      call HYPRE_StructSMGSetRelChange(solver, 0,ierr)
      call HYPRE_StructSMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructSMGSetNumPostRelax(solver,1,ierr)
      call HYPRE_StructSMGSetLogging(solver, 1,ierr)
      call HYPRE_StructSMGSetup(solver, mat, vec, x,ierr)
      call HYPRE_StructSMGSolve(solver, mat, vec, x,ierr)
      call HYPRE_StructSMGDestroy(solver,ierr)
      !call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
      !call HYPRE_StructSMGSetup(solver, mat, vec, x, ierr)
      !call HYPRE_StructSMGSolve(solver, mat, vec, x, ierr)
      !call HYPRE_StructSMGDestroy(solver,ierr)
      !call HYPRE_StructSMGGetNumIterations(solver, num_iterations,ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructPFMGSetTol(solver,tol,ierr)
      call HYPRE_StructPFMGSetRelChange(solver, 0,ierr)
      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetLogging(solver, 1,ierr)
      call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
      call HYPRE_StructPFMGDestroy(solver,ierr)
      !call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
      !call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
      !call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
      !call HYPRE_StructPFMGDestroy(solver,ierr)
      !call HYPRE_StructPFMGGetNumIterations(solver, num_iterations,ierr)
    elseif (HYPRESolverType .eq. HYPRESolverGMRES) then 
      call HYPRE_StructGMRESCreate(mpi_comm, solver,ierr)
      call HYPRE_StructGMRESSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructGMRESDestroy(solver,ierr)
      !call HYPRE_StructGMRESGetNumIterations(solver, num_iterations,ierr)
    elseif (HYPRESolverType .eq. HYPRESolverBiCGSTAB) then 
      call HYPRE_StructBiCGSTABCreate(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSetup(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABSolve(solver, mat, vec, x, ierr)
      call HYPRE_StructBiCGSTABDestroy(solver, mat, vec, x, ierr)
      !call HYPRE_StructBiCGSTABGetNumIterations(solver, num_iterations,ierr)
    endif ! HYPRESolverType
    !
    ! end of part based on the Paris Simulator code
    !
    ! fecth results
    !
    call HYPRE_StructVectorGetBoxValues(x,ijlower,ijupper,vecvalues,ierr)
    q = 0
    dsdt(:,:,:) = 0.d0
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q = q + 1
          if(phi(i,j,k).lt.0.d0) then
            dsdt(i,j,k) = vecvalues(q)
          endif
        enddo
      enddo
    enddo
    deallocate(vecvalues,matvalues)
    !
    return 
  end subroutine helmholtz_mg
  !
  subroutine destroy_solver_mg
    !
    implicit none
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( HYPRESolverType .eq. HYPRESolverSMG ) then 
      !call HYPRE_StructSMGDestroy(solver, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverPFMG ) then  
      !call HYPRE_StructPFMGDestroy(solver, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverGMRES ) then  
      !call HYPRE_StructGMRESDestroy(solver, ierr)
      call HYPRE_StructPFMGDestroy(precond, ierr)
    elseif ( HYPRESolverType .eq. HYPRESolverBiCGSTAB ) then  
      !call HYPRE_StructBiCGSTABDestroy(solver, ierr)
      call HYPRE_StructPFMGDestroy(precond, ierr)
    endif ! HYPRESolverType
    !
    ! end of part based on the Paris Simulator code
    !
    call HYPRE_StructGridDestroy(grid,ierr)
    call HYPRE_StructStencilDestroy(stencil,ierr)
    call HYPRE_StructMatrixDestroy(mat,ierr)
    call HYPRE_StructVectorDestroy(vec,ierr)
    call HYPRE_StructVectorDestroy(x,ierr)
    !
    return 
  end subroutine destroy_solver_mg
  !
!  subroutine helmholtz(n,dli,cbc,bc,c1,phi,t,alpha,dsdt)
!    !
!    ! Helmholtz solver for the species mass fraction equation with prescribed Dirichlet BC 
!    ! at the vapour/liquid interface Y = Yi, where Yi is obtained from thermodynamic
!    ! considerations.
!    ! it solves the equation in the domain phi<0
!    ! see Gibou et al. JCP 2002
!    !
!    ! note: old implementation which includes creation and destruction of the matrix,rhs vector and solution vector
!    !       at each time-step (not very efficient, but it works)
!    !
!    implicit none
!    integer, intent(in), dimension(3) :: n
!    real(8), intent(in), dimension(3) :: dli
!    character(len=1), intent(in), dimension(0:1,3) :: cbc
!    real(8), intent(in), dimension(0:1,3) :: bc
!    real(8), intent(in) :: c1
!    real(8), intent(in), dimension(-2:,-2:,-2:) :: phi,t
!    real(8), intent(in) :: alpha
!    real(8), intent(inout), dimension(:,:,:) :: dsdt
!    integer, dimension(3) :: ng
!    real(8), dimension(3) :: dl
!    real(8) :: cxp,cxm,cyp,cym,czp,czm,cc,rhs
!    real(8) :: xx,yy,zz,xi,yi,zi
!    real(8) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
!    integer :: i,j,k,ii,jj,kk,q,qq
!    real(8) :: theta
!    real(8), dimension(0:1,3) :: factor,sgn
!    !
!    ! hypre solver parameters/variables
!    !
!    real(8) :: final_res_norm,tol
!    integer :: mpi_comm
!    !
!    integer(8) grid
!    integer(8) stencil
!    integer, dimension(3) :: ijlower, ijupper
!    integer, parameter :: ndims = 3, nstencil = 7
!    integer, dimension(ndims) :: periods
!    integer, dimension(3,nstencil) :: offsets
!    integer(8) mat,vec,x
!    integer, dimension(nstencil) :: stencil_indices
!    real(8), allocatable, dimension(:) :: matvalues,vecvalues
!    integer(8) solver
!    integer :: nvalues,num_iterations,solverid
!    real(8) :: ti,si,coeff_t
!    !
!#ifdef DEBUG
!    if(myid.eq.0) print*,'solving Helmholtz equation...'
!#endif
!    !
!    ng(1:2) = n(1:2)*dims(1:2)
!    ng(3)   = n(3)
!    factor(:,:) = 0.d0
!    sgn(   :,:) = 0.d0
!    periods(:)  = 0
!    do q=1,3
!      do qq=0,1
!        select case(cbc(qq,q))
!        case('N')
!          factor(qq,q) = 1.d0/dli(q)*bc(qq,q)
!          sgn(   qq,q) = 1.d0
!        case('D')
!          factor(qq,q) = -2.d0*bc(qq,q)
!          sgn(   qq,q) = -1.d0
!        case('P')
!          factor(qq,q) = 0.d0
!          sgn(   qq,q) = 0.d0
!          periods(q) = ng(q)
!        end select
!      enddo
!    enddo
!    mpi_comm = comm_cart
!    dl(:) = 1.d0/dli(:)
!    !
!    ! initialize matrix
!    !
!    ! 1 - setup the grid
!    !
!    !   1.1 - create empty 2D grid object
!    !
!    call HYPRE_StructGridCreate(mpi_comm,ndims,grid,ierr)
!    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
!    !
!    !   1.2 - set grid extents and assemble
!    !
!    ijlower(1:2) = (coord(:)  )*n(1:2)+1
!    ijupper(1:2) = (coord(:)+1)*n(1:2)
!    ijlower(3  ) = 1
!    ijupper(3  ) = n(3)
!    call HYPRE_StructGridSetExtents(grid,ijlower,ijupper,ierr)
!    call HYPRE_StructGridAssemble(  grid,ierr)
!    !
!    ! 2 - setup the finite-difference stencil
!    !
!    call HYPRE_StructStencilCreate(ndims,nstencil,stencil,ierr)
!    offsets = reshape((/ 0, 0, 0, &
!                        -1, 0, 0, &
!                         1, 0, 0, &
!                         0,-1, 0, &
!                         0, 1, 0, &
!                         0, 0,-1, &
!                         0, 0, 1 /),shape(offsets))
!    do q = 1,nstencil
!     call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q), &
!                                        ierr)
!    enddo
!    !
!    ! 3 - create coefficient matrix and rhs vector
!    !
!    call HYPRE_StructMatrixCreate(mpi_comm,grid,stencil,mat,ierr)
!    call HYPRE_StructMatrixInitialize(mat,ierr)
!    call HYPRE_StructVectorCreate(mpi_comm,grid,vec,ierr)
!    call HYPRE_StructVectorInitialize(vec,ierr)
!    !
!    ! preliminaries for setting up the coefficient matrix
!    !
!    stencil_indices = (/0,1,2,3,4,5,6/)
!    nvalues = product(n(:))*nstencil ! number of grid points*number of stencil entries
!    allocate(matvalues(nvalues))
!    matvalues(:) = 0.d0
!    nvalues = product(n(:))
!    allocate(vecvalues(nvalues))
!    vecvalues(:) = 0.d0
!    !
!    ! compute stencil coefficients and rhs
!    !
!    q = 0
!    do k=1,n(3)
!      kk = k
!      zz = (kk-0.5)*dl(3)
!      do j=1,n(2)
!        jj = j+coord(2)*n(2)
!        yy = (jj-0.5)*dl(2)
!        do i=1,n(1)
!          ii = i+coord(1)*n(1)
!          xx = (ii-0.5)*dl(1)
!          q = q + 1
!          if(phi(i,j,k).lt.0.d0) then
!            !
!            alphaxp = 0.5d0*(alpha+alpha)
!            alphaxm = 0.5d0*(alpha+alpha)
!            alphayp = 0.5d0*(alpha+alpha)
!            alphaym = 0.5d0*(alpha+alpha)
!            alphazp = 0.5d0*(alpha+alpha)
!            alphazm = 0.5d0*(alpha+alpha)
!            !
!#ifdef TWOD
!            cxm = 0.d0
!            cxp = 0.d0
!#else
!            cxm = alphaxm*dli(1)**2*c1
!            cxp = alphaxp*dli(1)**2*c1
!#endif
!            cym = alphaym*dli(2)**2*c1
!            cyp = alphayp*dli(2)**2*c1
!            czm = alphazm*dli(3)**2*c1
!            czp = alphazp*dli(3)**2*c1
!            cc  = -(cxm+cxp+cym+cyp+czm+czp) + 1.d0
!            rhs = dsdt(i,j,k)
!            !
!#ifndef TWOD
!            if(periods(1).eq.0) then
!              if(    ii.eq.  1) then
!                rhs = rhs + cxm*factor(0,1)
!                cc = cc + sgn(0,1)*cxm
!                cxm = 0.d0
!              elseif(ii.eq.ng(1)) then
!                rhs = rhs + cxp*factor(1,1)
!                cc = cc + sgn(1,1)*cxp
!                cxp = 0.d0
!              endif
!            endif
!#endif
!            if(periods(2).eq.0) then
!              if(    jj.eq.  1 ) then
!                rhs = rhs + cym*factor(0,2)
!                cc = cc + sgn(0,2)*cym
!                cym = 0.d0
!              elseif(jj.eq.ng(2)) then
!                rhs = rhs + cyp*factor(1,2)
!                cc = cc + sgn(1,2)*cyp
!                cyp = 0.d0
!              endif
!            endif
!            if(periods(3).eq.0) then
!              if(    kk.eq.  1) then
!                rhs = rhs + czm*factor(0,3)
!                cc = cc + sgn(0,3)*czm
!                czm = 0.d0
!              elseif(kk.eq.ng(3)) then
!                rhs = rhs + czp*factor(1,3)
!                cc = cc + sgn(1,3)*czp
!                czp = 0.d0
!              endif
!            endif
!            !
!#ifndef TWOD
!            if(phi(i-1,j,k).gt.0.d0.and.(ii.ne.1.and.periods(1).eq.0)) then
!              xi = (xx*abs(phi(i-1,j,k))+(xx-dl(1))*abs(phi(i,j,k)))/(abs(phi(i-1,j,k))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i-1,j,k))+t(i-1,j,k)*abs(phi(i,j,k)))/(abs(phi(i-1,j,k))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(xx-xi)*dli(1)
!              if(theta.gt.dl(1)) then
!                rhs = rhs - cxm*si/theta
!                cc  = cc+cxm-(cxm/theta)
!              else 
!                rhs = rhs
!                cc  = cc+cxm
!              endif
!              cxm = 0.d0
!            endif
!            if(phi(i+1,j,k).gt.0.d0.and.(ii.ne.ng(1).and.periods(1).eq.0)) then
!              xi = (xx*abs(phi(i+1,j,k))+(xx+dl(1))*abs(phi(i,j,k)))/(abs(phi(i+1,j,k))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i+1,j,k))+t(i+1,j,k)*abs(phi(i,j,k)))/(abs(phi(i+1,j,k))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(xx-xi)*dli(1)
!              if(theta.gt.dl(1)) then
!                rhs = rhs - cxp*si/theta
!                cc  = cc+cxp-(cxp/theta)
!              else
!                rhs = rhs
!                cc  = cc+cxp
!              endif
!              cxp = 0.d0
!            endif
!#endif
!            if(phi(i,j-1,k).gt.0.d0.and.(jj.ne.1.and.periods(2).eq.0)) then
!              yi = (yy*abs(phi(i,j-1,k))+(yy-dl(2))*abs(phi(i,j,k)))/(abs(phi(i,j-1,k))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i,j-1,k))+t(i,j-1,k)*abs(phi(i,j,k)))/(abs(phi(i,j-1,k))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(yy-yi)*dli(2)
!              if(theta.gt.dl(2)) then
!                rhs = rhs - cym*si/theta
!                cc  = cc+cym-(cym/theta)
!              else
!                rhs = rhs
!                cc  = cc+cym
!              endif
!              cym = 0.d0
!            endif
!            if(phi(i,j+1,k).gt.0.d0.and.(jj.ne.ng(2).and.periods(2).eq.0)) then
!              yi = (yy*abs(phi(i,j+1,k))+(yy+dl(2))*abs(phi(i,j,k)))/(abs(phi(i,j+1,k))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i,j+1,k))+t(i,j+1,k)*abs(phi(i,j,k)))/(abs(phi(i,j+1,k))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(yy-yi)*dli(1)
!              if(theta.gt.dl(2)) then
!                rhs = rhs - cyp*si/theta
!                cc  = cc+cyp-(cyp/theta)
!              else
!                rhs = rhs
!                cc  = cc+cyp
!              endif
!              cyp = 0.d0
!            endif
!            !
!            if(phi(i,j,k-1).gt.0.d0.and.(kk.ne.1.and.periods(3).eq.0)) then
!              zi = (zz*abs(phi(i,j,k-1))+(zz-dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,k-1))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i,j,k-1))+t(i,j,k-1)*abs(phi(i,j,k)))/(abs(phi(i,j,k-1))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(zz-zi)*dli(3)
!              if(theta.gt.dl(3)) then
!                rhs = rhs - czm*si/theta
!                cc  = cc+czm-(czm/theta)
!              else
!                rhs = rhs
!                cc  = cc+czm
!              endif
!              czm = 0.d0
!            endif
!            if(phi(i,j,k+1).gt.0.d0.and.(kk.ne.ng(3).and.periods(3).eq.0)) then
!              zi = (zz*abs(phi(i,j,k+1))+(zz+dl(3))*abs(phi(i,j,k)))/(abs(phi(i,j,k+1))+abs(phi(i,j,k)))
!              ti = (t(i,j,k)*abs(phi(i,j,k+1))+t(i,j,k+1)*abs(phi(i,j,k)))/(abs(phi(i,j,k+1))+abs(phi(i,j,k)))
!              si = mass_fraction(ti)
!              theta = abs(zz-zi)*dli(3)
!              if(theta.gt.dl(3)) then
!                rhs = rhs - czp*si/theta
!                cc  = cc+czp-(czp/theta)
!              else
!                rhs = rhs
!                cc  = cc+czp
!              endif
!              czp = 0.d0
!            endif
!            !
!          else
!            cc  = 1.d0
!            cxm = 0.d0
!            cym = 0.d0
!            czm = 0.d0
!            cxp = 0.d0
!            cyp = 0.d0
!            czp = 0.d0
!            rhs = dsdt(i,j,k)
!          endif
!          !
!          matvalues((q-1)*nstencil+1) = cc
!          matvalues((q-1)*nstencil+2) = cxm
!          matvalues((q-1)*nstencil+3) = cxp
!          matvalues((q-1)*nstencil+4) = cym
!          matvalues((q-1)*nstencil+5) = cyp
!          matvalues((q-1)*nstencil+6) = czm
!          matvalues((q-1)*nstencil+7) = czp
!          vecvalues(q               ) = rhs
!        enddo
!      enddo
!    enddo
!    call HYPRE_StructMatrixSetBoxValues(mat,ijlower,ijupper,nstencil, &
!                                        stencil_indices,matvalues, &
!                                        ierr)
!    call HYPRE_StructMatrixAssemble(mat,ierr)
!    call HYPRE_StructVectorSetBoxValues(vec,ijlower,ijupper, &
!                                        vecvalues,ierr)
!    call HYPRE_StructVectorAssemble(vec,ierr)
!    !
!    ! create soluction vector
!    !
!    call HYPRE_StructVectorCreate(mpi_comm,grid,x,ierr)
!    call HYPRE_StructVectorInitialize(x,ierr)
!    vecvalues(:) = 0.d0
!    call HYPRE_StructVectorSetBoxValues(x,ijlower,ijupper, &
!                                        vecvalues,ierr)
!    call HYPRE_StructVectorAssemble(x,ierr)
!    !
!    ! setup solver, and solve
!    !
!    solverid = 2
!    if(solverid.eq.0) then
!      call HYPRE_StructPCGCreate(mpi_comm,solver,ierr)
!      !call HYPRE_StructPFMGSetMaxIter(solver,50,ierr)
!      call HYPRE_StructPCGSetTol(solver,tol,ierr)
!      call HYPRE_StructPCGSetPrintLevel(solver,2,ierr)
!      call HYPRE_StructPCGSetup(solver,mat,vec,x,ierr)
!      call HYPRE_StructPCGSolve(solver,mat,vec,x,ierr)
!      call HYPRE_StructPCGDestroy(solver,ierr)
!    elseif(solverid.eq.1) then
!      call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
!      call HYPRE_StructSMGSetMemoryUse(solver, 0,ierr)
!      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
!      call HYPRE_StructSMGSetTol(solver,tol,ierr)
!      call HYPRE_StructSMGSetRelChange(solver, 0,ierr)
!      call HYPRE_StructSMGSetNumPreRelax(solver, 1,ierr)
!      call HYPRE_StructSMGSetNumPostRelax(solver, 1,ierr)
!      call HYPRE_StructSMGSetLogging(solver, 1,ierr)
!      call HYPRE_StructSMGSetup(solver, mat, vec, x,ierr)
!      call HYPRE_StructSMGSolve(solver, mat, vec, x,ierr)
!      call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
!!      call HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &
!!           final_res_norm,ierr)
!      if(myid.eq.0) print*, 'Number of iterations: ', num_iterations
!      call HYPRE_StructSMGDestroy(solver,ierr)
!    elseif(solverid.eq.2) then
!      call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
!      !call HYPRE_StructPFMGSetMemoryUse(solver, 0,ierr)
!      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
!      call HYPRE_StructPFMGSetTol(solver,tol,ierr)
!      call HYPRE_StructPFMGSetRelChange(solver, 0,ierr)
!      call HYPRE_StructPFMGSetNumPreRelax(solver, 1,ierr)
!      call HYPRE_StructPFMGSetNumPostRelax(solver, 1,ierr)
!      call HYPRE_StructPFMGSetLogging(solver, 1,ierr)
!      call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
!      call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
!      call HYPRE_StructPFMGDestroy(solver,ierr)
!    endif
!    !
!    ! fetch results
!    !
!    call HYPRE_StructVectorGetBoxValues(x,ijlower,ijupper,vecvalues, &
!                                        ierr)
!    q = 0
!    dsdt(:,:,:) = 0.d0
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          q = q + 1
!          if(phi(i,j,k).lt.0.d0) then
!            dsdt(i,j,k) = vecvalues(q)
!          endif
!        enddo
!      enddo
!    enddo
!    call HYPRE_StructGridDestroy(grid,ierr)
!    call HYPRE_StructStencilDestroy(stencil,ierr)
!    call HYPRE_StructMatrixDestroy(mat,ierr)
!    call HYPRE_StructVectorDestroy(vec,ierr)
!    call HYPRE_StructVectorDestroy(x,ierr)
!    deallocate(vecvalues,matvalues)
!    !
!#ifdef DEBUG
!    if(myid.eq.0) print*,'Helmholtz equation solved!'
!#endif
!    !
!    return 
!  end subroutine helmholtz
!  !
!  subroutine solver_mg1(n,dli,dzci,cbc,bc,alpha,pg,p)
!    !
!    ! Helmholtz/poisson solver for a Poisson equation with variable
!    ! coefficients
!    !
!    ! note: old implementation which includes creation and destruction of the matrix,rhs vector and solution vector
!    !       at each time-step (not very efficient, but it works)
!    !
!    implicit none
!    integer, intent(in), dimension(3) :: n
!    real(8), intent(in), dimension(3) :: dli
!    real(8), intent(in), dimension(-2:) :: dzci
!    character(len=1), intent(in), dimension(0:1,3) :: cbc
!    real(8)         , intent(in), dimension(0:1,3) :: bc
!    real(8), intent(in   ), dimension(0:,0:,0:) :: alpha,pg
!    real(8), intent(inout), dimension(0:,0:,0:) :: p
!    integer, dimension(3) :: ng
!    real(8) :: cxp,cxm,cyp,cym,czp,czm,cc,rhs
!    real(8) :: alphaxp,alphaxm,alphayp,alphaym,alphazp,alphazm
!    integer :: i,j,k,ii,jj,kk,q,qq
!    real(8), dimension(0:1,3) :: factor,sgn
!    !
!    ! hypre solver parameters/variables
!    !
!    real(8) :: final_res_norm,tol
!    integer :: mpi_comm
!    !
!    integer(8) grid
!    integer(8) stencil
!    integer, dimension(3) :: ijlower, ijupper
!    integer, parameter :: ndims = 3, nstencil = 7
!    integer, dimension(ndims) :: periods
!    integer, dimension(3,nstencil) :: offsets
!    integer(8) mat,vec,x
!    integer, dimension(nstencil) :: stencil_indices
!    real(8), allocatable, dimension(:) :: matvalues,vecvalues, &
!                                          guessvalues
!    integer(8) solver
!    integer :: nvalues,num_iterations,solverid
!    !
!#ifdef DEBUG
!    if(myid.eq.0) print*,'solving Helmholtz equation...'
!#endif
!    !
!    ng(1:2) = n(1:2)*dims(1:2)
!    ng(3)   = n(3)
!    factor(:,:) = 0.d0
!    sgn(   :,:) = 0.d0
!    periods(:)  = 0
!    do q=1,3
!      do qq=0,1
!        select case(cbc(qq,q))
!        case('N')
!          factor(qq,q) = 1.d0/dli(q)*bc(qq,q)
!          sgn(   qq,q) = 1.d0
!        case('D')
!          factor(qq,q) = -2.d0*bc(qq,q)
!          sgn(   qq,q) = -1.d0
!        case('P')
!          factor(qq,q) = 0.d0
!          sgn(   qq,q) = 0.d0
!          periods(q) = ng(q)
!        end select
!      enddo
!    enddo
!    mpi_comm = comm_cart
!    !
!    ! initialize matrix
!    !
!    ! 1 - setup the grid
!    !
!    !   1.1 - create empty 2D grid object
!    !
!    call HYPRE_StructGridCreate(mpi_comm,ndims,grid,ierr)
!    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
!    !
!    !   1.2 - set grid extents and assemble
!    !
!    ijlower(1:2) = (coord(:)  )*n(1:2)+1
!    ijupper(1:2) = (coord(:)+1)*n(1:2)
!    ijlower(3  ) = 1
!    ijupper(3  ) = n(3)
!    call HYPRE_StructGridSetExtents(grid,ijlower,ijupper,ierr)
!    call HYPRE_StructGridAssemble(  grid,ierr)
!    !
!    ! 2 - setup the finite-difference stencil
!    !
!    call HYPRE_StructStencilCreate(ndims,nstencil,stencil,ierr)
!    offsets = reshape((/ 0, 0, 0, &
!                        -1, 0, 0, &
!                         1, 0, 0, &
!                         0,-1, 0, &
!                         0, 1, 0, &
!                         0, 0,-1, &
!                         0, 0, 1 /),shape(offsets))
!    do q = 1,nstencil
!     call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q), &
!                                        ierr)
!    enddo
!    !
!    ! 3 - create coefficient matrix and rhs vector
!    !
!    call HYPRE_StructMatrixCreate(mpi_comm,grid,stencil,mat,ierr)
!    call HYPRE_StructMatrixInitialize(mat,ierr)
!    call HYPRE_StructVectorCreate(mpi_comm,grid,vec,ierr)
!    call HYPRE_StructVectorInitialize(vec,ierr)
!    call HYPRE_StructVectorCreate(mpi_comm,grid,x,ierr)
!    call HYPRE_StructVectorInitialize(x,ierr)
!    !
!    ! preliminaries for setting up the coefficient matrix
!    !
!    stencil_indices = (/0,1,2,3,4,5,6/)
!    nvalues = product(n(:))*nstencil ! number of grid points*number of stencil entries
!    allocate(matvalues(nvalues))
!    matvalues(:) = 0.d0
!    nvalues = product(n(:))
!    allocate(vecvalues(nvalues))
!    vecvalues(:) = 0.d0
!    allocate(guessvalues(nvalues))
!    guessvalues(:) = 0.d0
!    !
!    ! compute stencil coefficients and rhs
!    !
!    q = 0
!    do k=1,n(3)
!      kk = k
!      do j=1,n(2)
!        jj = j+coord(2)*n(2)
!        do i=1,n(1)
!          ii = i+coord(1)*n(1)
!          q = q + 1
!          !
!#ifdef TWOD
!          alphaxp = 0.d0
!          alphaxm = 0.d0
!#else
!          alphaxp = (0.5d0*(alpha(i+1,j,k)+alpha(i,j,k)))**(-1)
!          alphaxm = (0.5d0*(alpha(i-1,j,k)+alpha(i,j,k)))**(-1)
!#endif
!          alphayp = (0.5d0*(alpha(i,j+1,k)+alpha(i,j,k)))**(-1)
!          alphaym = (0.5d0*(alpha(i,j-1,k)+alpha(i,j,k)))**(-1)
!          alphazp = (0.5d0*(alpha(i,j,k+1)+alpha(i,j,k)))**(-1)
!          alphazm = (0.5d0*(alpha(i,j,k-1)+alpha(i,j,k)))**(-1)
!          cxm = alphaxm*dli(1)**2
!          cxp = alphaxp*dli(1)**2
!          cym = alphaym*dli(2)**2
!          cyp = alphayp*dli(2)**2
!          czm = alphazm*dli(3)**2 ! change to dzci
!          czp = alphazp*dli(3)**2 ! change to dzci
!          cc  = -(cxm+cxp+cym+cyp+czm+czp)
!          rhs = p(i,j,k)
!          !
!#ifndef TWOD
!          if(periods(1).eq.0) then
!            if(    ii.eq.  1) then
!              rhs = rhs + cxm*factor(0,1)
!              cc = cc + sgn(0,1)*cxm
!              cxm = 0.d0
!            elseif(ii.eq.ng(1)) then
!              rhs = rhs + cxp*factor(1,1)
!              cc = cc + sgn(1,1)*cxp
!              cxp = 0.d0
!            endif
!          endif
!#endif
!          if(periods(2).eq.0) then
!            if(    jj.eq.  1 ) then
!              rhs = rhs + cym*factor(0,2)
!              cc = cc + sgn(0,2)*cym
!              cym = 0.d0
!            elseif(jj.eq.ng(2)) then
!              rhs = rhs + cyp*factor(1,2)
!              cc = cc + sgn(1,2)*cyp
!              cyp = 0.d0
!            endif
!          endif
!          if(periods(3).eq.0) then
!            if(    kk.eq.  1) then
!              rhs = rhs + czm*factor(0,3)
!              cc = cc + sgn(0,3)*czm
!              czm = 0.d0
!            elseif(kk.eq.ng(3)) then
!              rhs = rhs + czp*factor(1,3)
!              cc = cc + sgn(1,3)*czp
!              czp = 0.d0
!            endif
!          endif
!          matvalues((q-1)*nstencil+1) = cc
!          matvalues((q-1)*nstencil+2) = cxm
!          matvalues((q-1)*nstencil+3) = cxp
!          matvalues((q-1)*nstencil+4) = cym
!          matvalues((q-1)*nstencil+5) = cyp
!          matvalues((q-1)*nstencil+6) = czm
!          matvalues((q-1)*nstencil+7) = czp
!          vecvalues(q               ) = rhs
!          guessvalues(q             ) = pg(i,j,k)
!        enddo
!      enddo
!    enddo
!    !
!    call HYPRE_StructMatrixSetBoxValues(mat,ijlower,ijupper,nstencil, &
!                                        stencil_indices,matvalues, &
!                                        ierr)
!    call HYPRE_StructMatrixAssemble(mat,ierr)
!    call HYPRE_StructVectorSetBoxValues(vec,ijlower,ijupper, &
!                                        vecvalues,ierr)
!    call HYPRE_StructVectorAssemble(vec,ierr)
!    !
!    ! create soluction vector
!    !
!    !vecvalues(:) = 0.d0
!    call HYPRE_StructVectorSetBoxValues(x,ijlower,ijupper, &
!                                        guessvalues,ierr)
!    call HYPRE_StructVectorAssemble(x,ierr)
!    !
!    ! setup solver, and solve
!    !
!    solverid = 2
!    if(solverid.eq.0) then
!      call HYPRE_StructPCGCreate(mpi_comm,solver,ierr)
!      !call HYPRE_StructPFMGSetMaxIter(solver,50,ierr)
!      call HYPRE_StructPCGSetTol(solver,tol,ierr)
!      call HYPRE_StructPCGSetPrintLevel(solver,2,ierr)
!      call HYPRE_StructPCGSetup(solver,mat,vec,x,ierr)
!      call HYPRE_StructPCGSolve(solver,mat,vec,x,ierr)
!      call HYPRE_StructPCGDestroy(solver,ierr)
!    elseif(solverid.eq.1) then
!      call HYPRE_StructSMGCreate(mpi_comm, solver,ierr)
!      call HYPRE_StructSMGSetMemoryUse(solver, 0,ierr)
!      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
!      call HYPRE_StructSMGSetTol(solver,tol,ierr)
!      call HYPRE_StructSMGSetRelChange(solver, 0,ierr)
!      call HYPRE_StructSMGSetNumPreRelax(solver, 1,ierr)
!      call HYPRE_StructSMGSetNumPostRelax(solver, 1,ierr)
!      call HYPRE_StructSMGSetLogging(solver, 1,ierr)
!      call HYPRE_StructSMGSetup(solver, mat, vec, x,ierr)
!      call HYPRE_StructSMGSolve(solver, mat, vec, x,ierr)
!      call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
!!      call HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &
!!           final_res_norm,ierr)
!      if(myid.eq.0) print*, 'Number of iterations: ', num_iterations
!      call HYPRE_StructSMGDestroy(solver,ierr)
!    elseif(solverid.eq.2) then
!      call HYPRE_StructPFMGCreate(mpi_comm, solver,ierr)
!      !call HYPRE_StructPFMGSetMemoryUse(solver, 0,ierr)
!      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
!      call HYPRE_StructPFMGSetTol(solver,tol,ierr)
!      call HYPRE_StructPFMGSetRelChange(solver, 0,ierr)
!      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
!      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
!      call HYPRE_StructPFMGSetLogging(solver, 1,ierr)
!      call HYPRE_StructPFMGSetup(solver, mat, vec, x,ierr)
!      call HYPRE_StructPFMGSolve(solver, mat, vec, x,ierr)
!      call HYPRE_StructPFMGDestroy(solver,ierr)
!    endif
!    !
!    ! fetch results
!    !
!    call HYPRE_StructVectorGetBoxValues(x,ijlower,ijupper,vecvalues, &
!                                        ierr)
!    q = 0
!    p(:,:,:) = 0.d0
!    do k=1,n(3)
!      do j=1,n(2)
!        do i=1,n(1)
!          q = q + 1
!          p(i,j,k) = vecvalues(q)
!        enddo
!      enddo
!    enddo
!    call HYPRE_StructGridDestroy(grid,ierr)
!    call HYPRE_StructStencilDestroy(stencil,ierr)
!    call HYPRE_StructMatrixDestroy(mat,ierr)
!    call HYPRE_StructVectorDestroy(vec,ierr)
!    call HYPRE_StructVectorDestroy(x,ierr)
!    deallocate(vecvalues,matvalues,guessvalues)
!    !
!    return 
!  end subroutine solver_mg1
   !
end module mod_solver_vc
