!==============================================================================
! 2.5D model for simulating an Alfven wave-driven wind driven at the stellar
! surface by torsional magnetic field perturbations.
! Setup following the descriptions as given in the works of Ofman & Davila
! (1997), ApJ, 476 [solar] and Airapetian+ (2000), ApJ, 965 [red supergiant].
!
! Coded up by Florian Driessen (2025).
!==============================================================================
module mod_usr

  use mod_mhd
  use mod_constants,      only: const_c, const_G, const_LSun, const_MSun, &
       const_RSun, mp_cgs, kB_cgs
  use mod_kind_parameter, only: dp

  implicit none

  ! User input parameters
  real(dp) :: mstar_sol, rstar_sol, twind_cgs, b0_cgs, rho0_cgs, vdrive_cgs

  ! Dimensionless variables for computations
  real(dp) :: mstar, rstar, b0, rho0, pth0, twind, alpha, valf0, talf, csound
  real(dp) :: vdrive, omegadrive, rcrit, v0_parker, vinf_parker, const_Fr

contains

!==============================================================================
! This routine should set user methods and activate the physics module.
!==============================================================================
  subroutine usr_init

    call usr_params_read(par_files)

    unit_length   = rstar_sol * const_RSun
    unit_density  = rho0_cgs
    unit_velocity = b0_cgs / sqrt(4.0_dp*dpi*rho0_cgs)

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_bound
    usr_gravity        => stellar_gravity
    usr_aux_output     => set_extravar_output
    usr_add_aux_names  => set_extravarnames_output

    call set_coordinate_system("spherical_2.5D")
    call mhd_activate()

    if (mhd_energy) call mpistop("ERROR: No support for energy equation!")

  end subroutine usr_init

!==============================================================================
! Read in the usr.par file with the problem specific list.
!==============================================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n
    !--------------------------------------------------------------------------

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, b0_cgs, rho0_cgs, &
         vdrive_cgs

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

  end subroutine usr_params_read

!==============================================================================
! Compute some quantities of interest (in CGS) before making unitless.
!==============================================================================
  subroutine initglobaldata_usr

    character(len=*), parameter :: label_fmt = '(A30, " = ", ES12.6)'
    real(dp) :: mstar_cgs, rstar_cgs, valf0_cgs, talf_cgs, logg_cgs, mumol
    real(dp) :: csound_cgs, hscale_cgs, vesc_cgs, omegadrive_cgs
    real(dp) :: v0_parker_cgs, rcrit_cgs, vinf_parker_cgs
    real(dp) :: const_Lu, const_Eu, vesc, Ggrav
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar wind quantities
    valf0_cgs  = b0_cgs / sqrt(4*dpi*rho0_cgs)
    talf_cgs   = rstar_cgs / valf0_cgs
    logg_cgs   = log10(const_G * mstar_cgs/rstar_cgs**2.0_dp)
    mumol      = (1.0_dp + 4.0_dp*He_abundance) / (2.0_dp + 3.0_dp*He_abundance)
    csound_cgs = sqrt(twind_cgs * kB_cgs/(mumol * mp_cgs))
    hscale_cgs = csound_cgs**2.0_dp / 10.0_dp**logg_cgs
    vesc_cgs   = sqrt(2.0_dp * const_G * mstar_cgs / rstar_cgs)
    alpha      = rstar_cgs / hscale_cgs
    omegadrive_cgs = 12.0_dp / talf_cgs

    ! Parker wind quantities
    rcrit_cgs = 0.5_dp * const_G * mstar_cgs / csound_cgs**2.0_dp
    v0_parker_cgs = csound_cgs * (0.5_dp * vesc_cgs/csound_cgs)**2.0_dp &
         * exp(-0.5_dp*vesc_cgs**2.0_dp / csound_cgs**2.0_dp + 1.5_dp)
    vinf_parker_cgs = 2.0_dp*csound_cgs * sqrt(log(xprobmax1))

    ! Dimensionless fluid numbers
    ! Choice unit normalisations link Lundquist number directly to resistivity
    const_Fr = valf0_cgs**2.0_dp * rstar_cgs / (const_G * mstar_cgs)
    const_Eu = csound_cgs**2.0_dp / (mhd_gamma * valf0_cgs**2.0_dp)
    const_Lu = 1.0_dp / mhd_eta

    if (mype == 0 .and..not.convert) then
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '   Unity quantities   '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'unit length', unit_length
      write(*, label_fmt) 'unit density', unit_density
      write(*, label_fmt) 'unit velocity', unit_velocity
      write(*, label_fmt) 'unit numberdensity', unit_numberdensity
      write(*, label_fmt) 'unit pressure', unit_pressure
      write(*, label_fmt) 'unit temperature', unit_temperature
      write(*, label_fmt) 'unit magneticfield', unit_magneticfield
      write(*, label_fmt) 'unit time', unit_time
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '   Stellar and wind parameters in CGS units    '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'M/Msun', mstar_cgs/const_MSun
      write(*, label_fmt) 'R/Rsun', rstar_cgs/const_RSun
      write(*, label_fmt) 'Twind', twind_cgs
      write(*, label_fmt) 'Magnetic field strength', b0_cgs
      write(*, label_fmt) 'Mean molecular weight', mumol
      write(*, label_fmt) 'Helium abundance', He_abundance
      write(*, label_fmt) 'Alfven speed at surface', valf0_cgs
      write(*, label_fmt) 'Alfven travel time', talf_cgs
      write(*, label_fmt) 'log(g)', logg_cgs
      write(*, label_fmt) 'h/Rstar', hscale_cgs/rstar_cgs
      write(*, label_fmt) 'alpha=Rstar/h', alpha
      write(*, label_fmt) 'isothermal csound', csound_cgs
      write(*, label_fmt) 'vesc', vesc_cgs
      write(*, label_fmt) 'Omega driving', omegadrive_cgs
      write(*, label_fmt) 'Period driving (days)', &
           (2.0_dp * dpi) / (24.0_dp * 3600.0_dp * omegadrive_cgs)
      write(*, label_fmt) 'Parker rcrit/Rstar', rcrit_cgs/rstar_cgs
      write(*, label_fmt) 'Parker v0', v0_parker_cgs
      write(*, label_fmt) 'Parker vinf', vinf_parker_cgs
    endif

    ! Make some variables needed for computations dimensionless
    mstar       = mstar_cgs / (unit_density * unit_length**3.0_dp)
    rstar       = rstar_cgs / unit_length
    b0          = b0_cgs / unit_magneticfield
    csound      = csound_cgs / unit_velocity
    twind       = twind_cgs / unit_temperature
    rho0        = rho0_cgs / unit_density
    pth0        = rho0 * twind
    valf0       = valf0_cgs / unit_velocity
    talf        = talf_cgs / unit_time
    vdrive      = vdrive_cgs / unit_velocity
    omegadrive  = omegadrive_cgs * unit_time
    Ggrav       = const_G * unit_density * unit_time**2.0_dp
    v0_parker   = v0_parker_cgs / unit_velocity
    rcrit       = rcrit_cgs / unit_length
    vinf_parker = vinf_parker_cgs / unit_velocity
    vesc        = vesc_cgs / unit_velocity

    mhd_adiab = pth0 / rho0**mhd_gamma

    if (mype == 0 .and. .not.convert) then
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '  Dimensionless computation quantities  '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'Gravitational constant', Ggrav
      write(*, label_fmt) 'Mstar', mstar
      write(*, label_fmt) 'Rstar', rstar
      write(*, label_fmt) 'Surface density', rho0
      write(*, label_fmt) 'Surface pressure', pth0
      write(*, label_fmt) 'Wind temperature', twind
      write(*, label_fmt) 'Surface magnetic field', b0
      write(*, label_fmt) 'Surface Alfven speed', valf0
      write(*, label_fmt) 'Alfven timescale', talf
      write(*, label_fmt) 'Surface driving speed', vdrive
      write(*, label_fmt) 'Driving frequency', omegadrive
      write(*, label_fmt) 'Isothermal sound speed', csound
      write(*, label_fmt) 'Escape speed', vesc
      write(*,*)
      write(*, label_fmt) 'Parker surface speed', v0_parker
      write(*, label_fmt) 'Parker critical speed', csound
      write(*, label_fmt) 'Parker terminal wind speed', vinf_parker
      write(*,*)  
      write(*, label_fmt) 'adiabatic constant', mhd_adiab
      write(*, label_fmt) 'adiabatic gamma', mhd_gamma
      write(*, label_fmt) 'MHD resistivity', mhd_eta
      write(*, label_fmt) 'Lundquist no.', const_Lu
      write(*, label_fmt) 'Froude no.', const_Fr
      write(*, label_fmt) 'Euler no.', const_Eu 
    endif

    ! For cgs output in convert stage
    if (convert .and. saveprim) then
      w_convert_factor(rho_)   = unit_density
      w_convert_factor(mom(:)) = unit_velocity
      w_convert_factor(mag(:)) = unit_magneticfield
    endif

  end subroutine initglobaldata_usr

!==============================================================================
! Initial conditions start from spherically symmetric 1-D isothermal atmosphere
! in hydrostatic equilibrium with a purely radial magnetic field.
!==============================================================================
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    ! Local arguments
    integer  :: ir, itheta
    real(dp) :: vin, machparker, local_r, rho_static

    real(dp), parameter :: rho_r = 1.0_dp/9.0_dp, theta_s = 0.15_dp
    !--------------------------------------------------------------------------

    select case (iprob)
    case(1) ! Coronal hole of Ofman & Davila (1997), eqs. 5-6
      do ir = ixOmin1,ixOmax1
        local_r = x(ir,ixOmin2,1)

        do itheta = ixOmin2,ixOmax2
          rho_static = rho0 * exp(alpha * (1.0_dp / x(ir,itheta,1) - 1.0_dp))
          w(ir,itheta,rho_) = rho_static / rho_r * (1.0_dp - (1.0_dp - rho_r) &
               * exp(-( (x(ir,itheta,2) - 0.5_dp*dpi)/theta_s )**4.0_dp) )
        enddo
      enddo

    case(2) ! Parker wind solution
      do ir = ixOmin1,ixOmax1
        local_r = x(ir,ixOmin2,1)

        machparker = compute_parker_velocity(local_r, rcrit)

        w(ir^%1ixO^S,mom(1)) = machparker * csound
        w(ir^%1ixO^S,rho_) = rho0 * exp(-0.5_dp*(machparker*csound)**2.0_dp) &
             * exp(-alpha * (1.0_dp - 1.0_dp / x(ir^%1ixO^S,1)))
      enddo

    case default
      call mpistop('initial_conditions: choose iprob = 1,2')
    end select

    ! Mass density from hydrostatic isothermal atmosphere
    w(ixO^S,mom(1)) = 0.0_dp

    ! Velocity components
    w(ixO^S,mom(2)) = 0.0_dp
    w(ixO^S,mom(3)) = 0.0_dp

    ! Magnetic field components
    w(ixO^S,mag(1)) = b0 / x(ixO^S,1)**2.0_dp
    w(ixO^S,mag(2)) = 0.0_dp
    w(ixO^S,mag(3)) = 0.0_dp

    ! If using Dedner+(2002) divergence cleaning
    if (mhd_glm) w(ixO^S,psi_) = 0.0_dp

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

!==============================================================================
! Determine the outflowing transonic solution of the Parker wind.
!==============================================================================
  function compute_parker_velocity(r, rcritp) result(mach)

    ! Function arguments
    real(dp), intent(in) :: r, rcritp

    ! Local variables
    integer,  parameter :: maxiter = 20
    real(dp), parameter :: eps = 1.0e-7_dp
    integer  :: i
    real(dp) :: mach, delta, fx, fxp, x0, xout
    !--------------------------------------------------------------------------

    if (r <= rcritp) then
      x0 = v0_parker / csound
    else
      x0 = 2.0_dp
    endif

    xout = x0
    i = 0
    do while (i < maxiter)

      fx = xout**2.0_dp - 2.0_dp*log(xout) - 4.0_dp*log(r/rcritp) &
          - 4.0_dp*rcritp/r + 3.0_dp
      fxp = 2.0_dp*xout - 2.0_dp / xout
      delta = -fx / fxp
      xout = xout + delta

      if (abs(delta) < eps) exit
      i = i + 1
    enddo

    if (i > maxiter) then
      write(*,*) "compute_parker_velocity: no convergence of Newton method"
      call mpistop('QUITTING')
    endif

    mach = xout

  end function compute_parker_velocity

!==============================================================================
! Special user boundary conditions at inner radial boundary 
!   rho, vr, vphi (extrapolated); vtheta, Br, Btheta, Bphi (fixed)
!==============================================================================
  subroutine special_bound(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: ixI^L, ixB^L, iB
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir
    !--------------------------------------------------------------------------

    select case (iB)
    case(1)

      call mhd_to_primitive(ixI^L,ixI^L,w,x)

      ! Inner boundary is subsonic, sub-Alfvenic such that there are 3 outgoing
      ! characteristics and 4 incoming characteristics

      ! Note: convention here is from flow-based perspective
      ! 3 outgoing characteristics: first order (linear) extrapolation done
      do ir = ixEmax1,ixEmin1,-1
        w(ir^%1ixB^S,rho_)   = 2.0_dp*w(ir+1^%1ixB^S,rho_) &
             - w(ir+2^%1ixB^S,rho_)
        w(ir^%1ixB^S,mom(1)) = 2.0_dp*w(ir+1^%1ixB^S,mom(1)) &
             - w(ir+2^%1ixB^S,mom(1))
        w(ir^%1ixB^S,mom(3)) = 2.0_dp*w(ir+1^%1ixB^S,mom(3)) &
             - w(ir+2^%1ixB^S,mom(3))
      enddo

      w(ixB^S,mom(1)) = max(0.0_dp, w(ixB^S,mom(1)))

      ! 4 incoming characteristics: physical values assigned
      ! Polar velocity component
      w(ixB^S,mom(2)) = 0.0_dp

      ! Magnetic field components
      w(ixB^S,mag(1)) = b0
      w(ixB^S,mag(2)) = 0.0_dp
      w(ixB^S,mag(3)) = vdrive / valf0 * cos(omegadrive * qt)

      if (mhd_glm) w(ixB^S,psi_) = 0.0_dp

      call mhd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==============================================================================
! Compute stellar gravity.
!==============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer,  intent(in)  :: ixI^L, ixO^L
    real(dp), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(dp), intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixO^S,:) = 0.0_dp

    ! Only in radial direction
    gravity_field(ixO^S,1) = -1.0_dp / (const_Fr * x(ixO^S,1)**2.0_dp)

  end subroutine stellar_gravity

!==============================================================================
! Computes and stores additional variables of interest in convert stage.
!==============================================================================
  subroutine set_extravar_output(ixI^L, ixO^L, w, x, normconv)

    ! Subroutine arguments
    integer,  intent(in) :: ixI^L, ixO^L
    real(dp), intent(in) :: x(ixI^S,1:ndim)
    real(dp)             :: w(ixI^S,nw+nwauxio)
    real(dp)             :: normconv(0:nw+nwauxio)

    ! Local variable
    integer  :: idirmin
    real(dp) :: divbboy(ixI^S), gradvalf(ixI^S), current(ixI^S,7-2*ndir:3)
    !--------------------------------------------------------------------------

    ! Output the Alfven speed by summing squared Bfield in each direction
    if (B0field) then
      w(ixO^S,nw+1) = sqrt( sum((w(ixO^S,mag(:)) &
           + block%B0(ixO^S,:,0))**2, dim=ndim+1) / w(ixO^S,rho_) )
    else
      w(ixO^S,nw+1) = sqrt(sum(w(ixO^S,mag(:))**2, dim=ndim+1)/w(ixO^S,rho_))
    endif

    ! Output divB, for Tanaka splitting this will be divB1
    call get_divb(w,ixI^L,ixO^L,divbboy)
    w(ixO^S,nw+2) = divbboy(ixO^S)

    ! Ohmic heating rate density
    call get_current(w,ixI^L,ixO^L,idirmin,current)
    w(ixO^S,nw+3) = sum(current(ixO^S,:)**2, dim=ndim+1) / mhd_eta

    ! Gradient of Alfven speed
    !call gradient(w(ixO^S,nw+1),ixI^L,ixO^L,1,gradvalf)
    !w(ixO^S,nw+4) = gradvalf(ixO^S)

  end subroutine set_extravar_output

!==============================================================================
! Additional auxiliary io variables: 
!   Alfven velocity, divergence B, Ohmic heating rate density
!==============================================================================
  subroutine set_extravarnames_output(varnames)

    character(len=*) :: varnames
    !--------------------------------------------------------------------------

    varnames = 'vAlf divB Hohm'

  end subroutine set_extravarnames_output

end module mod_usr
