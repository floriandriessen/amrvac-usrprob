!==============================================================================
! Make an initial condition using a 1D CAK wind model by interpolating it on
! the simulation grid. This is a modification of AMRVAC mod_oneblock module.
! Works with 1D, 2D, 3D (M)HD radiation-driven wind models.
!
! NOTE: if input grid < simulation grid the outer wind will be held constant at
!       values of outermost density and radial velocity of input wind.
!
! To make AMRVAC aware of this module include a 'local.make' file with rules.
!
! Coded up by Florian Driessen: 2018-2022 (PhD).
!==============================================================================
module mod_init_cakwind

  use mod_kind_parameter, only: dp

  implicit none

  ! Number of radial cells
  integer :: nrcells

  ! Arrays holding radial grid and hydrodynamic variables
  real(dp), allocatable :: xoneblock(:), woneblock(:,:)

contains

!==============================================================================
! Read in a relaxed 1D CAK profile stored in a .blk file that is produced
! with the 1D CAK wind code. Assumed input format: rho, vr, ...
!==============================================================================
  subroutine read_oned_cakwind(filename)

    use mpi
    use mod_global_parameters, only: domain_nx1, mype, npe, icomm, ierrmpi
    use mod_comm_lib,          only: mpistop

    ! Subroutine argument
    character(len=*), intent(in) :: filename

    ! Local variables
    integer  :: ir, unit = 94
    logical  :: alive
    !--------------------------------------------------------------------------

    ! Master does the reading
    if (mype == 0) then
      inquire(file=trim(filename), exist=alive)

      if (alive) then
        open(unit, file=trim(filename), status='unknown')
      else
        call mpistop('Input file you want to use cannot be found!')
      endif

      print*,'Wind initialisation with input file: ', trim(filename)
      print*

      ! The header information:
      read(unit,*) ! skip variables
      read(unit,*) nrcells
      read(unit,*) ! skip time information

      if (nrcells < domain_nx1) then
        print*,'Input grid contains less radial cells than simulation grid'
        print*,'NOTE: constant extrapolation rho, vr done in outer wind'
      endif

      ! Allocate and read the grid and variables
      allocate(xoneblock(nrcells))
      allocate(woneblock(nrcells,1:2))

      do ir = 1,nrcells
        read(unit,*) xoneblock(ir), woneblock(ir,:)
      enddo

      close(unit)
    endif

    call MPI_BARRIER(icomm,ierrmpi)

    ! Broadcast what mype=0 read
    if (npe > 1) then
      call MPI_BCAST(nrcells,1,MPI_INTEGER,0,icomm,ierrmpi)

      if (mype /= 0) then
        allocate(xoneblock(nrcells))
        allocate(woneblock(nrcells,1:2))
      endif

      call MPI_BCAST(xoneblock,nrcells,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
      call MPI_BCAST(woneblock,nrcells*2,MPI_DOUBLE_PRECISION,0,icomm,ierrmpi)
    endif

  end subroutine read_oned_cakwind

!==============================================================================
! Linear interpolation of 1D CAK profile from .blk file onto simulation grid.
!==============================================================================
  subroutine interp_oned_cakwind_on_grid(local_r,interp_rho,interp_vr)

    ! Subroutine arguments
    real(dp), intent(in)  :: local_r
    real(dp), intent(out) :: interp_rho, interp_vr

    ! Local arguments
    integer  :: ic1, ic11, ic21
    real(dp) :: xd1, rp
    !--------------------------------------------------------------------------

    rp = local_r

    ! Find correct block index to get correct value
    ic1 = minloc(abs(rp - xoneblock(:)), dim=1, mask=.true.)

    if (xoneblock(ic1) <= rp) then
      ic11 = ic1
    else
      ic11 = ic1 - 1
    endif
    ic21 = ic11 + 1

    ! Constant extrapolation when outside the range of input radial grid
    if (ic11 <= 1) then
      ic11 = 1
      ic21 = ic11 + 1
      rp   = xoneblock(ic1)
    endif

    if (ic21 >= nrcells) then
      ic21 = ic21
      ic11 = ic21 - 1
      rp   = xoneblock(ic11)
    endif

    ! Interpolate
    xd1        = (rp - xoneblock(ic11)) / (xoneblock(ic21) - xoneblock(ic11))
    interp_rho = woneblock(ic11,1) * (1.0_dp - xd1) + woneblock(ic21,1) * xd1
    interp_vr  = woneblock(ic11,2) * (1.0_dp - xd1) + woneblock(ic21,2) * xd1

  end subroutine interp_oned_cakwind_on_grid

end module mod_init_cakwind
