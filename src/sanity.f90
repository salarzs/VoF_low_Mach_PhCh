module mod_sanity
  !
  use iso_c_binding , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr
  use mod_correc    , only: correc
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_initsolver, only: initsolver
  use mod_param     , only: small
  use mod_solver    , only: solver
  !
  implicit none
  !
  private
  public test_sanity
  !
  contains
  !
  subroutine test_sanity(ng,n,dims,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                         dli,dzci,dzfi)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    !
    integer         , intent(in), dimension(3)         :: ng,n
    integer         , intent(in), dimension(2)         :: dims
    character(len=1), intent(in), dimension(0:1,3,3)   :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)     :: cbcpre
    real(8)         , intent(in), dimension(0:1,3,3)   :: bcvel
    real(8)         , intent(in), dimension(0:1,3)     :: bcpre
    logical         , intent(in), dimension(0:1,3)     :: is_outflow
    logical         , intent(in), dimension(3)         :: is_forced
    real(8)         , intent(in), dimension(3)         :: dli
    real(8)         , intent(in), dimension(-2:n(3)+3) :: dzci,dzfi
    !
    logical :: passed
    !
    call chk_dims(ng,dims,passed);                 if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed); if(.not.passed) call abortit
    call chk_outflow(cbcpre,is_outflow,passed);    if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced  ,passed);   if(.not.passed) call abortit 
    !call chk_solvers(n,dli,dzci,dzfi,cbcvel,cbcpre,bcvel,bcpre,is_outflow,passed)
    !if(.not.passed) call abortit
    !
    return
  end subroutine test_sanity
  !
  subroutine chk_dims(ng,dims,passed)
    !
    implicit none
    !
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(2) :: dims
    logical, intent(out) :: passed
    logical :: passed_loc
    passed = .true.
    passed_loc = all(mod(ng(:),2).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot, jtot and ktot should be even.'
    passed = passed.and.passed_loc
    passed_loc = all(mod(ng(1:2),dims(1:2)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: itot and jtot should be divisable by dims(1) and dims(2), respectively.'
    passed = passed.and.passed_loc
    passed_loc = (mod(ng(2),dims(1)).eq.0).and.(mod(ng(3),dims(2)).eq.0)
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: jtot should be divisable by both dims(1) and dims(2), and &
                     &ktot should be divisable by dims(2)'
    passed = passed.and.passed_loc
    !
    return
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,bcvel,bcpre,passed)
    !
    implicit none
    !
    character(len=1), intent(in ), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in ), dimension(0:1,3  ) :: cbcpre
    real(8)         , intent(in ), dimension(0:1,3,3) :: bcvel
    real(8)         , intent(in ), dimension(0:1,3  ) :: bcpre
    logical         , intent(out)                     :: passed
    character(len=2) :: bc01v,bc01p                   
    integer :: ivel,idir                              
    logical :: passed_loc                             
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v.eq.'PP').or. &
                                      (bc01v.eq.'ND').or. &
                                      (bc01v.eq.'DN').or. &
                                      (bc01v.eq.'NN').or. &
                                      (bc01v.eq.'DD') )
      enddo
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p.eq.'PP').or. &
                                    (bc01p.eq.'ND').or. &
                                    (bc01p.eq.'DN').or. &
                                    (bc01p.eq.'NN').or. &
                                    (bc01p.eq.'DD') )
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.' 
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v.eq.'PP'.and.bc01p.eq.'PP').or. &
                                    (bc01v.eq.'ND'.and.bc01p.eq.'DN').or. &
                                    (bc01v.eq.'DN'.and.bc01p.eq.'ND').or. &
                                    (bc01v.eq.'DD'.and.bc01p.eq.'NN').or. &
                                    (bc01v.eq.'NN'.and.bc01p.eq.'DD') )
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir).eq.0.d0).and.(bcpre(1,idir).eq.0.d0))
    enddo
    if(myid.eq.0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.d0).'
    passed = passed.and.passed_loc
    !
    return 
  end subroutine chk_bc
  !
  subroutine chk_outflow(cbcpre,is_outflow,passed)
    !
    implicit none
    !
    logical         , intent(in), dimension(0:1,3  ) :: is_outflow
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    logical         , intent(out) :: passed
    integer :: idir,ibound
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and outflow BC
    !
    do idir=1,3
      do ibound = 0,1
        passed = passed.and. &
                 (cbcpre(ibound,idir).eq.'D'.and.(is_outflow(ibound,idir))) .or. &
                 (.not.is_outflow(ibound,idir))
      enddo
    enddo
    if(myid.eq.0.and.(.not.passed)) &
      print*, 'ERROR: Dirichlet pressure BC should be an outflow direction; check the BC or is_outflow in bc.h90.'
    !
    return 
  end subroutine chk_outflow
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    !
    implicit none
    !
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(in), dimension(3) :: is_forced
    logical         , intent(out) :: passed
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and forcing BC
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir).eq.'PP')
      endif
    enddo
    if(myid.eq.0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in bc.h90.'
    !
    return 
  end subroutine chk_forcing
  !
  subroutine abortit
    !
    implicit none
    !
    if(myid.eq.0) print*, ''
    if(myid.eq.0) print*, '*** Simulation aborted due to ierrs in the case file ***'
    if(myid.eq.0) print*, '    check bc.h90 and setup.h90'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    call exit
    !
    return
  end subroutine abortit
  !
end module mod_sanity
