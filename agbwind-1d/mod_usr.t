!==============================================================================
! 1D model for simulating a pulsation-driven wind from a Mira-like star.
! Setup mainly following the work of Bowen (1988), ApJ 329.
! 
! Pulsation only happens temporally not spatially. Additionally, the velocity
! perturbation is applied at the photosphere (= stellar radius), not in a point
! below the photosphere.
!
! Coded up by Florian Driessen (2025).
!==============================================================================
module mod_usr

  use mod_hd
  use mod_constants,      only: const_LSun, const_MSun, const_RSun, const_c, &
       const_years, mH_cgs, const_G, const_sigma
  use mod_kind_parameter, only: dp

  implicit none

  ! User input parameters
  real(dp) :: mstar_sol, v_amp_cgs, puls_period_day, teff_cgs, tcond1_cgs
  real(dp) :: rhostar_cgs, beta, grad_to_ggrav
  logical  :: use_lya_cooling

  ! Dimensionless computation variables
  real(dp) :: lstar, rstar, rhostar, v_amp, puls_period, tcond1, clight
  real(dp) :: teff, delta_bowen, gmstar, kmax_bowen, v_amp0
  real(dp) :: coeff_bowen_cool, coeff_lya_cool

  ! Mean molecular weight
  real(dp), parameter :: mu_bar = 1.26_dp

  ! Initial boundary velocity amplitude of 1 cm/s; Bowen (1988), p. 302
  real(dp), parameter :: v_amp0_cgs = 1.0_dp

  ! ISM values for floor temperature and density (1 H particle/cm^3)
  real(dp), parameter :: floor_temperature_cgs = 10.0_dp
  real(dp), parameter :: floor_density_cgs = mH_cgs

  ! Exponent for density power-law stratification
  real(dp), parameter :: alpha = 10.0_dp

  ! Cooling 'time constant' (g s cm^-3), width condensation temperature range;
  ! Bowen (1988), table 1
  real(dp), parameter :: cool_bowen_cgs = 1.0e-5_dp
  real(dp), parameter :: delta_bowen_cgs = 60.0_dp

  ! Extra variable index
  integer :: itemp_, iqbow_, iqlya_, grad_

contains

!==============================================================================
! This routine should set user methods and activate the physics module.
!==============================================================================
  subroutine usr_init()

    call params_read(par_files)

    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid  => initial_conditions
    usr_special_bc     => special_boundary
    usr_gravity        => stellar_gravity
    usr_source         => cooling_sources
    usr_process_grid   => update_temperature
    usr_get_dt         => special_dt
    usr_add_aux_names  => set_extravarnames_output
    usr_aux_output     => set_extravar_output

    call set_coordinate_system("spherical")
    call hd_activate()
    
    ! Add user auxiliary variable
    itemp_ = var_set_auxvar("Tgas", "Tgas")

    ! Extra variables for output
    iqbow_ = var_set_extravar("Q_bowen", "Q_bowen")
    iqlya_ = var_set_extravar("Q_lya", "Q_lya")
    grad_  = var_set_extravar("grad", "grad")

    if ( (hd_radiative_cooling .and. rc_fl%coolcurve == 'SPEX') .and. &
         use_lya_cooling ) &
         call mpistop('Turn off Lya cooling. Included in SPEX cooling curve.')

  end subroutine usr_init
  
!==============================================================================
! Read in the usr.par file with a user problem specific list.
!==============================================================================
  subroutine params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n
    !--------------------------------------------------------------------------

    namelist /usr_list/  mstar_sol, v_amp_cgs, puls_period_day, teff_cgs, &
         beta, rhostar_cgs, tcond1_cgs, grad_to_ggrav, &
         use_lya_cooling

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, usr_list, end=111)
      111 close(unitpar)
    end do

  end subroutine params_read

!==============================================================================
! Initialise global variables to be used during run. This subroutine is read
! at the first iteration and after a restart.
!==============================================================================
  subroutine initglobaldata_usr

    use mod_hd_phys, only: rc_fl
    use mod_radiative_cooling

    ! Local variables
    real(dp) :: unit_ggrav, unit_lum, unit_opacity, unit_energy
    real(dp) :: lstar_cgs, mstar_cgs, log_rstar_sol, rstar_cgs, puls_period_cgs
    real(dp) :: log_g_cgs, hscale_cgs, csound_iso_cgs, kmax_bowen_cgs
    real(dp) :: mstar, csound_iso

    real(dp), parameter :: const_day_cgs = 8.64e4_dp
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun

    ! Ostlie & Cox (1986), ApJ 311, sec. IVa
    log_rstar_sol = ( log10(puls_period_day) + 1.92_dp &
         + 0.73_dp*log10(mstar_sol) ) / 1.86_dp
    rstar_cgs     = 10.0_dp**log_rstar_sol * const_RSun

    ! Stellar and wind properties
    lstar_cgs       = 4.0_dp*dpi * const_sigma * rstar_cgs**2 * teff_cgs**4
    puls_period_cgs = puls_period_day * const_day_cgs
    csound_iso_cgs  = sqrt(teff_cgs * kB_cgs / (mu_bar * mH_cgs))
    log_g_cgs       = log10(const_G * mstar_cgs / rstar_cgs**2.0_dp)
    hscale_cgs      = csound_iso_cgs**2.0_dp / 10.0_dp**log_g_cgs

    ! Ratio of g_rad/g_grav; Bowen (1988), eq. 6b
    kmax_bowen_cgs = grad_to_ggrav * 4.0_dp*dpi * const_G * mstar_cgs &
         * const_c / lstar_cgs

    ! Code units
    unit_length        = rstar_cgs
    unit_density       = rhostar_cgs
    unit_time          = puls_period_cgs
    unit_velocity      = unit_length / unit_time
    unit_numberdensity = unit_density / (mu_bar * mH_cgs)
    unit_mass          = unit_density * unit_length**3.0_dp
    unit_pressure      = unit_density * unit_velocity**2.0_dp
    unit_temperature   = mu_bar * mH_cgs * unit_velocity**2.0_dp / kB_cgs
    unit_opacity       = unit_length**2.0_dp / unit_mass

    unit_ggrav   = unit_density * unit_time**2.0_dp
    unit_lum     = unit_density * unit_length**5.0_dp / unit_time**3.0_dp
    unit_energy  = unit_length**2.0_dp * unit_mass / unit_time**2.0_dp

    lstar       = lstar_cgs / unit_lum
    mstar       = mstar_cgs / unit_mass
    rstar       = rstar_cgs / unit_length
    rhostar     = rhostar_cgs / unit_density
    teff        = teff_cgs / unit_temperature
    v_amp       = v_amp_cgs / unit_velocity
    v_amp0      = v_amp0_cgs / unit_velocity
    puls_period = puls_period_cgs / unit_time
    tcond1      = tcond1_cgs / unit_temperature
    gmstar      = const_G * unit_ggrav * mstar
    clight      = const_c / unit_velocity
    csound_iso  = csound_iso_cgs / unit_velocity
 
    kmax_bowen       = kmax_bowen_cgs / unit_opacity
    delta_bowen      = delta_bowen_cgs / unit_temperature
    coeff_bowen_cool = cool_bowen_cgs / (unit_density * unit_time)

    ! Coefficient of Lyman alpha cooling (erg cm3/s); Spitzer (1978)
    coeff_lya_cool = 7.5e-19_dp &
         / (unit_energy / (unit_time * unit_numberdensity))

    ! Reset some AMRVAC values
    small_density     = floor_density_cgs / unit_density
    small_temperature = floor_temperature_cgs / unit_temperature
    small_pressure    = small_density * small_temperature

    if (hd_radiative_cooling) rc_fl%tlow = small_temperature

    if (mype == 0 .and. .not.convert) then
      print*, '============================'
      print*, '   Unity quantities (cgs)   '
      print*, '============================'
      print*, 'unit length        = ', unit_length
      print*, 'unit density       = ', unit_density
      print*, 'unit velocity      = ', unit_velocity
      print*, 'unit numberdensity = ', unit_numberdensity
      print*, 'unit pressure      = ', unit_pressure
      print*, 'unit temperature   = ', unit_temperature
      print*, 'unit time          = ', unit_time
      print*, 'unit mass          = ', unit_mass
      print*, 'unit opacity       = ', unit_opacity
      print*
      print*, '======================================'
      print*, '   Problem parameters in cgs units    '
      print*, '======================================'
      print*, 'Lstar/Lsun         = ', lstar_cgs / const_LSun
      print*, 'Mstar/Msun         = ', mstar_cgs / const_MSun
      print*, 'Rstar/Rsun         = ', rstar_cgs / const_RSun
      print*, 'Tstar              = ', teff_cgs
      print*, 'log(g)             = ', log_g_cgs
      print*, 'H/Rstar            = ', hscale_cgs / rstar_cgs
      print*, 'Mean mol. weight   = ', mu_bar
      print*, 'P pulsation        = ', puls_period_cgs
      print*, 'Bowen k_max        = ', kmax_bowen_cgs
      print*, 'Bowen delta        = ', delta_bowen_cgs
      print*, 'iso .csound        = ', csound_iso_cgs
      print*, 'adiabatic gamma    = ', hd_gamma
      print*, 'Unit time in years = ', unit_time / const_years
      print*
      print*, 'small density     = ', small_density * unit_density
      print*, 'small pressure    = ', small_pressure * unit_pressure
      print*, 'small temperature = ', small_temperature * unit_temperature 
      print*
      print*, '========================================'
      print*, '    Dimensionless AMRVAC quantities     '
      print*, '========================================'
      print*, 'Extra computed unit quantities:'
      print*, '   unit luminosity = ', unit_lum
      print*, '   unit cgravity   = ', unit_ggrav
      print*, 'Lstar        = ', lstar
      print*, 'Mstar        = ', mstar
      print*, 'Rstar        = ', rstar
      print*, 'Tstar        = ', teff
      print*, 'Tcond1       = ', tcond1
      print*, 'Tfloor       = ', small_temperature
      print*, 'rhosurface   = ', rhostar
      print*, 'G*Mstar      = ', gmstar
      print*, 'iso. csound  = ', csound_iso
      print*, 'clight       = ', clight
      print*
      print*, 'Vsurf amp.   = ', v_amp
      print*, 'Vsurf amp.0  = ', v_amp0
      print*, 'P pulsation  = ', puls_period
      print*, 'Bowen k_max  = ', kmax_bowen
      print*, 'Bowen delta  = ', delta_bowen
      print*, 'Bowen cool coeff. = ', coeff_bowen_cool
      print*, 'Lya cool coeff.   = ', coeff_lya_cool
      print*
      print*, 'small density     = ', small_density
      print*, 'small pressure    = ', small_pressure
      print*, 'small temperature = ', small_temperature
      print*
    end if

  end subroutine initglobaldata_usr

!==============================================================================
! Initial conditions from 1D spherically symmetric, hydrostatic equilibrium.
!==============================================================================
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: ixI^L, ixO^L
    real(dp), intent(in)    :: x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !--------------------------------------------------------------------------

    ! Radial velocity
    w(ixI^S,mom(1)) = 0.0_dp

    ! Mass density; Boulangier+ (2019), eq. 10
    w(ixI^S,rho_) = rhostar * (rstar / x(ixI^S,1))**alpha

    ! Temperature; Boulangier+ (2019), eq. 11
    w(ixI^S,itemp_) = teff * (rstar / x(ixI^S,1))**beta

    ! Pressure from ideal gas
    w(ixI^S,p_) = w(ixI^S,rho_) * w(ixI^S,itemp_)

    where (w(ixI^S,rho_) < small_density) w(ixI^S,rho_) = small_density
    where (w(ixI^S,p_) < small_pressure) w(ixI^S,p_) = small_pressure
    where (w(ixI^S,itemp_) < small_temperature) &
         w(ixI^S,itemp_) = small_temperature
          
    call hd_to_conserved(ixI^L, ixI^L, w, x)

  end subroutine initial_conditions

!==============================================================================
! Special user boundary conditions at inner radial boundary:
!   rho, vr, T, p (fixed)
!==============================================================================
  subroutine special_boundary(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: ixI^L, ixB^L, iB
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    real(dp) :: v_amp_bc
    !--------------------------------------------------------------------------

    ! Slowly build up pulsation over 20 periods to avoid big shock wave
    if (qt < 20.0_dp*puls_period) then
      v_amp_bc = v_amp0 + qt * (v_amp - v_amp0) / (20.0_dp*puls_period)
    else
      v_amp_bc = v_amp
    endif

    select case (iB)
    case(1)

      w(ixB^S,rho_)   = rhostar
      w(ixB^S,itemp_) = teff
      w(ixB^S,mom(1)) = v_amp_bc * cos(2.0_dp*dpi * qt / puls_period)
      w(ixB^S,p_)     = rhostar * teff

      call hd_to_conserved(ixI^L, ixB^L, w, x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_boundary

!==============================================================================
! Extra cooling processes for the wind from Bowen, Spitzer, and Hartquist.
!==============================================================================
  subroutine cooling_sources(qdt, ixI^L, ixO^L, iw^LIM, qtC, wCT, qt, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: ixI^L, ixO^L, iw^LIM
    real(dp), intent(in)    :: qdt, qtC, qt
    real(dp), intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    integer  :: ix^D
    real(dp) :: tempeq, q_bowen, q_lya, n_H, n_e, kappa_dust
    real(dp) :: eth(ixO^S), eth_min(ixO^S), ptherm(ixI^S), grad(ixO^S)
    !--------------------------------------------------------------------------

    ! Radiative acceleration due to dust
    {do ix^DB = ixO^LIM^DB\}
      ! Dust temperature is radiative equilibrium temperature
      tempeq = teff * (rstar / x(ix^D,1))**beta

      kappa_dust = kmax_bowen &
            / (1.0_dp + exp((tempeq - tcond1) / delta_bowen))

      grad(ix^D) = kappa_dust * lstar / (4.0_dp*dpi*clight * x(ix^D,1)**2.0_dp)

      w(ix^D,grad_) = grad(ix^D)

      w(ix^D,mom(1)) = w(ix^D,mom(1)) + qdt * grad(ix^D) * wCT(ix^D,rho_)
      w(ix^D,e_)     = w(ix^D,e_) + qdt * grad(ix^D) * wCT(ix^D,mom(1))
    {enddo^D&\}

    ! Minimum thermal energy based on starting state of time integrator
    call hd_get_pthermal(w,x,ixI^L,ixO^L,ptherm)
    eth_min(ixO^S) = w(ixO^S,rho_) * small_temperature / (hd_gamma - 1.0_dp)
    eth(ixO^S) = ptherm(ixO^S) / (hd_gamma - 1.0_dp) - eth_min(ixO^S)

    {do ix^DB = ixO^LIM^DB\}
      ! Thermal cooling of gas to radiative equilibrium; Bowen (1988), eq. 3b
      ! (multiply with rho for units of energy density per unit time)
      tempeq = teff * (rstar / x(ix^D,1))**beta

      q_bowen = (1.5_dp*(wCT(ix^D,itemp_) - tempeq) * wCT(ix^D,rho_) &
           / coeff_bowen_cool) * wCT(ix^D,rho_)

      ! Restrict maximum cooling
      q_bowen = min( q_bowen, eth(ix^D) / qdt )

      w(ix^D,iqbow_) = q_bowen
      w(ix^D,e_) = w(ix^D,e_) - qdt * q_bowen
    {enddo^D&\}

    ! Lyman alpha cooling of Spitzer (1978)
    if (use_lya_cooling) then
      call hd_get_pthermal(w,x,ixI^L,ixO^L,ptherm)
      eth_min(ixO^S) = w(ixO^S,rho_) * small_temperature / (hd_gamma - 1.0_dp)
      eth(ixO^S) = ptherm(ixO^S) / (hd_gamma - 1.0_dp) - eth_min(ixO^S)
      
      {do ix^DB = ixO^LIM^DB\}
        if (wCT(ix^D,itemp_)*unit_temperature < 6000.0_dp) cycle

        ! Assume n_e = 0.001*n_H; Bowen (1988), IIb, p. 301
        ! 20% of the radiation is lost; Bowen (1988), IIc(iii), p. 304
        n_H = wCT(ix^D,rho_)
        n_e = 1.0e-3_dp*n_H

        q_lya = 0.8_dp*coeff_lya_cool * n_e * n_H &
              * exp(-118400.0_dp / (wCT(ix^D,itemp_)*unit_temperature))

        ! Restrict maximum cooling
        q_lya = min( q_lya, eth(ix^D) / qdt )

        w(ix^D,iqlya_) = q_lya
        w(ix^D,e_) = w(ix^D,e_) - qdt * q_lya
      {enddo^D&\}
    endif

  end subroutine cooling_sources

!==============================================================================
! Ensure before new advection step that temperature is from hydrodynamic state
! of last iteration.
!==============================================================================
  subroutine update_temperature(igrid, level, ixI^L, ixO^L, qt, w, x)

    ! Subroutine arguments
    integer,  intent(in)    :: igrid, level, ixI^L, ixO^L
    real(dp), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(dp), intent(inout) :: w(ixI^S,1:nw)
    !--------------------------------------------------------------------------

    call hd_to_primitive(ixI^L,ixO^L,w,x)

    w(ixO^S,itemp_) = w(ixO^S,p_) / w(ixO^S,rho_)

    where (w(ixO^S,itemp_) < small_temperature) &
         w(ixO^S,itemp_) = small_temperature

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine update_temperature

!==============================================================================
! Gravity field of the central star.
!==============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer,          intent(in)  :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(in)  :: wCT(ixI^S,1:nw)
    double precision, intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixI^S,1:ndim) = 0.0_dp

    gravity_field(ixO^S,1) = -gmstar / x(ixO^S,1)**2.0_dp

  end subroutine stellar_gravity

!==============================================================================
! Adjust the hydrodynamic timestep for pulsations.
!==============================================================================
  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

    ! Subroutine arguments
    integer,          intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    ! Local variable
    real(dp) :: max_grad
    !--------------------------------------------------------------------------

    dtnew = bigdouble

    ! Radiative force due to dust
    max_grad = maxval(abs(w(ixO^S,grad_)) / block%ds(ixO^S,1))
    max_grad = max(max_grad, smalldouble)
    dtnew    = min(dtnew, 1.0_dp / sqrt(max_grad))

    ! Stellar pulsation; Bowen (1988), IIb, p. 300
    dtnew = min(dtnew, puls_period * 2.0e-2_dp)

  end subroutine special_dt

!==============================================================================
! Computes and stores additional variables of interest in convert stage.
!==============================================================================
  subroutine set_extravar_output(ixI^L, ixO^L, w, x, normconv)

    ! Subroutine arguments
    integer,  intent(in) :: ixI^L,ixO^L
    real(dp), intent(in) :: x(ixI^S,1:ndim)
    real(dp)             :: w(ixI^S,1:nw+nwauxio)
    real(dp)             :: normconv(0:nw+nwauxio)
    !--------------------------------------------------------------------------

    ! Mass-loss rate (Msun/yr)
    w(ixO^S,nw+1) = 4.0_dp*dpi * x(ixO^S,1)**2.0_dp * w(ixO^S,mom(1)) &
           * unit_mass / unit_time * const_years / const_MSun
     
  end subroutine set_extravar_output

!==============================================================================
! Additional auxiliary io variable: mass-loss rate.
!==============================================================================
  subroutine set_extravarnames_output(varnames)

    ! Subroutine argument
    character(len=*) :: varnames
    !--------------------------------------------------------------------------

    varnames = 'Mdot'

  end subroutine set_extravarnames_output
  
end module mod_usr
