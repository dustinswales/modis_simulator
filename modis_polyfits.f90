!ifort -L/usr/local/ifort/lib -lnetcdff -I/usr/local/ifort/include modis_polyfits.f90 -o modis_polyfits.o
program modis_polyfits
  use netcdf
  implicit none
  ! Algorithmic parameters
  real,parameter :: &
       re_water_min = 4.,  & ! Minimum effective radius (liquid)
       re_water_max = 30., & ! Maximum effective radius (liquid)
       re_ice_min   = 5.,  & ! Minimum effective radius (ice)
       re_ice_max   = 90.    ! Minimum effective radius (ice)
  integer,parameter :: &
       num_trial_res = 15, & ! Increase to make the linear pseudo-retrieval of size more accurate
       phaseIsLiquid = 1,  & !
       phaseIsIce    = 2     !

  ! Testing parameters
  real, dimension(2), parameter :: &
       re_liq_range = (/1.,35./),  & ! Liquid particle-size range (microns)
       re_ice_range = (/1.,95./),  & ! Ice particle-size range (microns)
       tau_range    = (/1.,20./)     ! Optical-depth range (1)

  ! Local varaibles
  integer :: iSize, nRe_liq, nRe_ice, iTau, nTau, ncid, status, dimID1, dimID2, &
       dimID3, varID(20)
  real,dimension(:), allocatable :: &
       g_liq_old, g_liq_new, ssa_liq_old, ssa_liq_new, re_liq, &
       g_ice_old, g_ice_new, ssa_ice_old, ssa_ice_new, re_ice, &
       tau
  real, dimension(:,:), allocatable :: &
       refl_ice_old, refl_ice_new, refl_liq_old, refl_liq_new, &
       sizeICE_old,  sizeICE_new,  sizeLIQ_old,  sizeLIQ_new
  real, dimension(num_trial_res) :: &
       predicted_Refl_ice_old, predicted_Refl_ice_new, &
       predicted_Refl_liq_old, predicted_Refl_liq_new
  
  ! ##################################################################################
  ! Initialization (from cosp_modis_interface.F90_init())
  ! ##################################################################################
  integer :: i
  real,dimension(num_trial_res) :: &
       trial_re_w, & ! Near-IR optical params vs size for retrieval scheme (liquid)
       trial_re_i    ! Near-IR optical params vs size for retrieval scheme (ice)
  real,dimension(num_trial_res) :: &
       g_w_new,        & ! Assymettry parameter for size retrieval (liquid)
       g_i_new,        & ! Assymettry parameter for size retrieval (ice)
       w0_w_new,       & ! Single-scattering albedo for size retrieval (liquid)
       w0_i_new,       & ! Single-scattering albedo for size retrieval (ice)
       g_w_old,        & ! Assymettry parameter for size retrieval (liquid)
       g_i_old,        & ! Assymettry parameter for size retrieval (ice)
       w0_w_old,       & ! Single-scattering albedo for size retrieval (liquid)
       w0_i_old          ! Single-scattering albedo for size retrieval (ice)
  
  ! Precompute near-IR optical params vs size for retrieval scheme    
  trial_re_w(1:num_trial_res) = re_water_min + (re_water_max - re_water_min) /         &
       (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
  trial_re_i(1:num_trial_res) = re_ice_min   + (re_ice_max -   re_ice_min) /           &
       (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
  
  ! Initialize estimates for size retrieval
  g_w_new(1:num_trial_res)  = get_g_nir_new(  phaseIsLiquid, trial_re_w(1:num_trial_res))
  w0_w_new(1:num_trial_res) = get_ssa_nir_new(phaseIsLiquid, trial_re_w(1:num_trial_res))
  g_i_new(1:num_trial_res)  = get_g_nir_new(  phaseIsIce,    trial_re_i(1:num_trial_res))
  w0_i_new(1:num_trial_res) = get_ssa_nir_new(phaseIsIce,    trial_re_i(1:num_trial_res))
  g_w_old(1:num_trial_res)  = get_g_nir_old(  phaseIsLiquid, trial_re_w(1:num_trial_res))
  w0_w_old(1:num_trial_res) = get_ssa_nir_old(phaseIsLiquid, trial_re_w(1:num_trial_res))
  g_i_old(1:num_trial_res)  = get_g_nir_old(  phaseIsIce,    trial_re_i(1:num_trial_res))
  w0_i_old(1:num_trial_res) = get_ssa_nir_old(phaseIsIce,    trial_re_i(1:num_trial_res))
  ! ##################################################################################
  ! ##################################################################################

  ! Problem size
  nRe_liq = re_liq_range(2) - re_liq_range(1) + 1
  nRe_ice = re_ice_range(2) - re_ice_range(1) + 1
  nTau    = tau_range(2)    - tau_range(1)    + 1

  ! Allocate space
  allocate(tau(nTau), re_liq(nRe_liq), re_ice(nRe_ice),            &
           refl_ice_old(nTau,nRe_ice), refl_ice_new(nTau,nRe_ice), &
           refl_liq_old(nTau,nRe_ice), refl_liq_new(nTau,nRe_ice), &
           g_liq_old(nRe_liq),         g_liq_new(nRe_liq),         &
           g_ice_old(nRe_ice),         g_ice_new(nRe_ice),         &
           ssa_liq_old(nRe_liq),       ssa_liq_new(nRe_liq),       &
           ssa_ice_old(nRe_ice),       ssa_ice_new(nRe_ice),       &
           sizeICE_old(nTau,nRe_ice),  sizeICE_new(nTau,nRe_ice),  &
           sizeLIQ_old(nTau,nRe_liq),  sizeLIQ_new(nTau,nRe_liq))
  
  ! Compute optical properties
  ! Liquid
  do iSize=re_liq_range(1),re_liq_range(2)
     re_liq(iSize)      = iSize
     g_liq_old(iSize)   = get_g_nir_old(  phaseIsLiquid,real(iSize))
     g_liq_new(iSize)   = get_g_nir_new(  phaseIsLiquid,real(iSize))
     ssa_liq_old(iSize) = get_ssa_nir_old(phaseIsLiquid,real(iSize))
     ssa_liq_new(iSize) = get_ssa_nir_new(phaseIsLiquid,real(iSize))
  enddo
  ! Ice
  do iSize=re_ice_range(1),re_ice_range(2)
     re_ice(iSize)      = iSize
     g_ice_old(iSize)   = get_g_nir_old(  phaseIsIce,real(iSize))
     g_ice_new(iSize)   = get_g_nir_new(  phaseIsIce,real(iSize))
     ssa_ice_old(iSize) = get_ssa_nir_old(phaseIsIce,real(iSize))
     ssa_ice_new(iSize) = get_ssa_nir_new(phaseIsIce,real(iSize))
  enddo
  ! Optical depth
  tau(1:5) = [0.3,0.5,0.9,1.3,3.0]
  do iTau=tau_range(1)+5,tau_range(2)
     tau(iTau) = real(iTau)
  enddo
  
  ! Run forward calculation(s), perform size retrieval...
  do iTau=tau_range(1),tau_range(2)
     ! Ice clouds
     do iSize=re_ice_range(1),re_ice_range(2)
        ! Old fits
        refl_ice_old(iTau,iSize)  = compute_toa_reflectace(1, tau(iTau), g_ice_old(iSize), ssa_ice_old(iSize))
        predicted_refl_ice_old(:) = two_stream_reflectance(tau(iTau), g_i_old, w0_i_old)
        sizeICE_old(iTau,iSize)   = interpolate_to_min(trial_re_i, predicted_Refl_ice_old, refl_ice_old(iTau,iSize))
        ! New fits
        refl_ice_new(iTau,iSize)  = compute_toa_reflectace(1, tau(iTau), g_ice_new(iSize), ssa_ice_new(iSize))
        predicted_refl_ice_new(:) = two_stream_reflectance(tau(iTau), g_i_new, w0_i_new)
        sizeICE_new(iTau,iSize)   = interpolate_to_min(trial_re_i, predicted_Refl_ice_new, refl_ice_new(iTau,iSize))
     enddo
     ! Liquid clouds
     do iSize=re_liq_range(1),re_liq_range(2)
        ! Old fits
        refl_liq_old(iTau,iSize)  = compute_toa_reflectace(1, tau(iTau), g_liq_old(iSize), ssa_liq_old(iSize))
        predicted_refl_liq_old(:) = two_stream_reflectance(tau(iTau), g_w_old, w0_w_old)
        sizeLIQ_old(iTau,iSize)   = interpolate_to_min(trial_re_w, predicted_Refl_liq_old, refl_liq_old(iTau,iSize))
        ! New fits
        refl_liq_new(iTau,iSize)  = compute_toa_reflectace(1, tau(iTau), g_liq_new(iSize), ssa_liq_new(iSize))
        predicted_refl_liq_new(:) = two_stream_reflectance(tau(iTau), g_w_new, w0_w_new)
        sizeLIQ_new(iTau,iSize)   = interpolate_to_min(trial_re_w, predicted_Refl_liq_new, refl_liq_new(iTau,iSize))
     enddo
  enddo
  
  ! Write output to netcdf.
  status = nf90_create("modis_polyfits.nc",nf90_clobber,ncid)
  status = nf90_def_dim(ncid,'nre_liq',nRe_liq,dimID1)
  status = nf90_def_dim(ncid,'nre_ice',nRe_ice,dimID2)
  status = nf90_def_dim(ncid,'ntau',   nTau,   dimID3)
  status = nf90_def_var(ncid,'re_liq',     nf90_float,(/dimID1/),varID(1))
  status = nf90_def_var(ncid,'g_liq_new',  nf90_float,(/dimID1/),varID(2))
  status = nf90_def_var(ncid,'g_liq_old',  nf90_float,(/dimID1/),varID(3))
  status = nf90_def_var(ncid,'ssa_liq_new',nf90_float,(/dimID1/),varID(4))
  status = nf90_def_var(ncid,'ssa_liq_old',nf90_float,(/dimID1/),varID(5))
  status = nf90_def_var(ncid,'re_ice',     nf90_float,(/dimID2/),varID(6))
  status = nf90_def_var(ncid,'g_ice_new',  nf90_float,(/dimID2/),varID(7))
  status = nf90_def_var(ncid,'g_ice_old',  nf90_float,(/dimID2/),varID(8))
  status = nf90_def_var(ncid,'ssa_ice_new',nf90_float,(/dimID2/),varID(9))
  status = nf90_def_var(ncid,'ssa_ice_old',nf90_float,(/dimID2/),varID(10))
  status = nf90_def_var(ncid,'tau',        nf90_float,(/dimID3/),varID(11))
  status = nf90_def_var(ncid,'sizeLIQ_old',nf90_float,(/dimID3,dimID1/),varID(12))
  status = nf90_def_var(ncid,'sizeLIQ_new',nf90_float,(/dimID3,dimID1/),varID(13))
  status = nf90_def_var(ncid,'sizeICE_old',nf90_float,(/dimID3,dimID2/),varID(14))
  status = nf90_def_var(ncid,'sizeICE_new',nf90_float,(/dimID3,dimID2/),varID(15))
  status = nf90_enddef(ncid)
  status = nf90_put_var(ncid,varID(1),re_liq)
  status = nf90_put_var(ncid,varID(2),g_liq_new)
  status = nf90_put_var(ncid,varID(3),g_liq_old)
  status = nf90_put_var(ncid,varID(4),ssa_liq_new)
  status = nf90_put_var(ncid,varID(5),ssa_liq_old)
  status = nf90_put_var(ncid,varID(6),re_ice)
  status = nf90_put_var(ncid,varID(7),g_ice_new)
  status = nf90_put_var(ncid,varID(8),g_ice_old)
  status = nf90_put_var(ncid,varID(9),ssa_ice_new)
  status = nf90_put_var(ncid,varID(10),ssa_ice_old)
  status = nf90_put_var(ncid,varID(11),tau)
  status = nf90_put_var(ncid,varID(12),sizeLIQ_old)
  status = nf90_put_var(ncid,varID(13),sizeLIQ_new)
  status = nf90_put_var(ncid,varID(14),sizeICE_old)
  status = nf90_put_var(ncid,varID(15),sizeICE_new)
  status = nf90_close(ncid)
  
contains

  ! ########################################################################################
  pure function interpolate_to_min(x, y, yobs)
    ! INPUTS
    real,intent(in),dimension(num_trial_res) :: x, y 
    real,intent(in)                          :: yobs
    ! OUTPUTS
    real                                     :: interpolate_to_min
    ! LOCAL VARIABLES
    real, dimension(num_trial_res)           :: diff
    integer                                      :: nPoints, minDiffLoc, lowerBound, upperBound
    
    ! Given a set of values of y as y(x), find the value of x that minimizes abs(y - yobs)
    !   y must be monotonic in x
    
    nPoints = size(y)
    diff(1:num_trial_res) = y(1:num_trial_res) - yobs
    minDiffLoc = minloc(abs(diff), dim = 1) 
    
    if(minDiffLoc == 1) then 
       lowerBound = minDiffLoc
      upperBound = minDiffLoc + 1
    else if(minDiffLoc == nPoints) then
      lowerBound = minDiffLoc - 1
      upperBound = minDiffLoc
    else
      if(diff(minDiffLoc-1) * diff(minDiffLoc) < 0) then
        lowerBound = minDiffLoc-1
        upperBound = minDiffLoc
      else 
        lowerBound = minDiffLoc
        upperBound = minDiffLoc + 1
      end if 
    end if 
    
    if(diff(lowerBound) * diff(upperBound) < 0) then     
      !
      ! Interpolate the root position linearly if we bracket the root
      !
      interpolate_to_min = x(upperBound) - & 
                           diff(upperBound) * (x(upperBound) - x(lowerBound)) / (diff(upperBound) - diff(lowerBound))
    else 
      interpolate_to_min = -999.
    end if 
    

  end function interpolate_to_min  
  
  ! ########################################################################################
  elemental function get_g_nir_old (phase, re)
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    ! INPUTS
    integer, intent(in) :: phase
    real,intent(in) :: re
    ! OUTPUTS
    real            :: get_g_nir_old 
    ! LOCAL VARIABLES(parameters)
    real, dimension(3), parameter :: &
         ice_coefficients         = (/ 0.7432,  4.5563e-3, -2.8697e-5 /), & 
         small_water_coefficients = (/ 0.8027, -1.0496e-2,  1.7071e-3 /), & 
         big_water_coefficients   = (/ 0.7931,  5.3087e-3, -7.4995e-5 /) 
    
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
       if(re < 8.) then 
          get_g_nir_old = fit_to_quadratic(re, small_water_coefficients)
          if(re < re_water_min) get_g_nir_old = fit_to_quadratic(re_water_min, small_water_coefficients)
       else
          get_g_nir_old = fit_to_quadratic(re,   big_water_coefficients)
          if(re > re_water_max) get_g_nir_old = fit_to_quadratic(re_water_max, big_water_coefficients)
       end if
    else
       get_g_nir_old = fit_to_quadratic(re, ice_coefficients)
       if(re < re_ice_min) get_g_nir_old = fit_to_quadratic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_g_nir_old = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if
  end function get_g_nir_old
  
  ! ########################################################################################
  elemental function get_ssa_nir_old (phase, re)
    ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    
    ! INPUTS
    integer, intent(in) :: phase
    real,intent(in) :: re
    ! OUTPUTS
    real            :: get_ssa_nir_old
    ! LOCAL VARIABLES (parameters)
    real, dimension(4), parameter :: ice_coefficients   = (/ 0.9994, -4.5199e-3, 3.9370e-5, -1.5235e-7 /)
    real, dimension(3), parameter :: water_coefficients = (/ 1.0008, -2.5626e-3, 1.6024e-5 /) 
    
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
       get_ssa_nir_old = fit_to_quadratic(re, water_coefficients)
       if(re < re_water_min) get_ssa_nir_old = fit_to_quadratic(re_water_min, water_coefficients)
       if(re > re_water_max) get_ssa_nir_old = fit_to_quadratic(re_water_max, water_coefficients)
    else
       get_ssa_nir_old = fit_to_cubic(re, ice_coefficients)
       if(re < re_ice_min) get_ssa_nir_old = fit_to_cubic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_ssa_nir_old = fit_to_cubic(re_ice_max, ice_coefficients)
    end if
    
  end function get_ssa_nir_old
  
  ! ########################################################################################
  elemental function get_g_nir_new (phase, re)
    !
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    !
    integer, intent(in) :: phase
    real,    intent(in) :: re
    real :: get_g_nir_new 
    real, dimension(3), parameter :: ice_coefficients         = (/ 0.7490, 6.5153e-3, -5.4136e-5 /), &
                                         small_water_coefficients = (/ 1.0364, -8.8800e-2, 7.0000e-3 /)
    real, dimension(4), parameter :: big_water_coefficients   = (/ 0.6035, 2.8993e-2, -1.1051e-3, 1.5134e-5 /)
    ! approx. fits from MODIS Collection 6 LUT scattering calculations for 3.7 µm channel size retrievals
    if(phase == phaseIsLiquid) then 
       if(re < 7.) then
          get_g_nir_new = fit_to_quadratic(re, small_water_coefficients)
          if(re < re_water_min) get_g_nir_new = fit_to_quadratic(re_water_min, small_water_coefficients)
       else
          get_g_nir_new = fit_to_cubic(re, big_water_coefficients)
          if(re > re_water_max) get_g_nir_new = fit_to_cubic(re_water_max, big_water_coefficients)
       end if
    else
       get_g_nir_new = fit_to_quadratic(re, ice_coefficients)
       if(re < re_ice_min) get_g_nir_new = fit_to_quadratic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_g_nir_new = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if
  end function get_g_nir_new
  
  ! ########################################################################################
  elemental function get_ssa_nir_new (phase, re)
    integer, intent(in) :: phase
    real,    intent(in) :: re
    real                :: get_ssa_nir_new
    !
    ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    !
    real, dimension(4), parameter :: ice_coefficients   = (/ 0.9625, -1.8069e-2, 3.3281e-4,-2.2865e-6/)
    real, dimension(3), parameter :: water_coefficients = (/ 1.0044, -1.1397e-2, 1.3300e-4 /)
    
    ! approx. fits from MODIS Collection 6 LUT scattering calculations
    if(phase == phaseIsLiquid) then
       get_ssa_nir_new = fit_to_quadratic(re, water_coefficients)
       if(re < re_water_min) get_ssa_nir_new = fit_to_quadratic(re_water_min, water_coefficients)
       if(re > re_water_max) get_ssa_nir_new = fit_to_quadratic(re_water_max, water_coefficients)
    else
       get_ssa_nir_new = fit_to_cubic(re, ice_coefficients)
       if(re < re_ice_min) get_ssa_nir_new = fit_to_cubic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_ssa_nir_new = fit_to_cubic(re_ice_max, ice_coefficients)
    end if
  end function get_ssa_nir_new
  
  ! ########################################################################################
  pure function fit_to_cubic(x, coefficients) 
    ! INPUTS
    real,               intent(in) :: x
    real, dimension(4), intent(in) :: coefficients
    ! OUTPUTS
    real                           :: fit_to_cubic  
    
    fit_to_cubic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3) + x * coefficients(4)))
  end function fit_to_cubic
    
  ! ########################################################################################
  pure function fit_to_quadratic(x, coefficients) 
    ! INPUTS
    real,               intent(in) :: x
    real, dimension(3), intent(in) :: coefficients
    ! OUTPUTS
    real                           :: fit_to_quadratic
    
    fit_to_quadratic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3)))
  end function fit_to_quadratic

  ! ########################################################################################
  ! Radiative transfer
  ! ########################################################################################
  pure function compute_toa_reflectace(nLevels,tau, g, w0)
    ! This wrapper reports reflectance only and strips out non-cloudy elements from the 
    ! calculation
    
    ! INPUTS
    integer,intent(in)                     :: nLevels
    real,intent(in),dimension(nLevels) :: tau, g, w0
    ! OUTPUTS
    real                               :: compute_toa_reflectace
    ! LOCAL VARIABLES
    logical, dimension(nLevels)                   :: cloudMask
    integer, dimension(count(tau(1:nLevels) > 0)) :: cloudIndicies
    real,dimension(count(tau(1:nLevels) > 0)) :: Refl,Trans
    real                                      :: Refl_tot, Trans_tot
    integer                                       :: i
    
    cloudMask(1:nLevels) = tau(1:nLevels) > 0. 
    cloudIndicies = pack((/ (i, i = 1, nLevels) /), mask = cloudMask) 
    do i = 1, size(cloudIndicies)
       call two_stream(tau(cloudIndicies(i)), g(cloudIndicies(i)), w0(cloudIndicies(i)), Refl(i), Trans(i))
    end do
    
    call adding_doubling(count(tau(1:nLevels) > 0),Refl(:), Trans(:), Refl_tot, Trans_tot)  
    
    compute_toa_reflectace = Refl_tot
    
  end function compute_toa_reflectace
  
  ! ########################################################################################
  pure subroutine two_stream(tauint, gint, w0int, ref, tra) 
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    ! INPUTS
    real, intent(in)  :: tauint, gint, w0int
    ! OUTPUTS
    real, intent(out) :: ref, tra
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real,parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den, th
    
    ! Compute reflectance and transmittance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    f   = gint**2
    tau = (1. - w0int * f) * tauint
    w0  = (1. - f) * w0int / (1. - w0int * f)
    g   = (gint - f) / (1. - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7. - w0* (4. + 3. * g)) / 4.
    gamma2 = -(1. - w0* (4. - 3. * g)) / 4.
    gamma3 =  (2. - 3.*g*xmu) / 4.
    gamma4 =   1. - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
  
          ref = rh / (1. + gamma1 * tau)
          tra = 1. - ref       
      else if(beam == 2) then
          ref = gamma1*tau/(1. + gamma1*tau)
          tra = 1. - ref
      endif
    else
      ! Non-conservative scattering
      a1 = gamma1 * gamma4 + gamma2 * gamma3
      a2 = gamma1 * gamma3 + gamma2 * gamma4

      rk = sqrt(gamma1**2 - gamma2**2)
      
      r1 = (1. - rk * xmu) * (a2 + rk * gamma3)
      r2 = (1. + rk * xmu) * (a2 - rk * gamma3)
      r3 = 2. * rk *(gamma3 - a2 * xmu)
      r4 = (1. - (rk * xmu)**2) * (rk + gamma1)
      r5 = (1. - (rk * xmu)**2) * (rk - gamma1)
      
      t1 = (1. + rk * xmu) * (a1 + rk * gamma4)
      t2 = (1. - rk * xmu) * (a1 - rk * gamma4)
      t3 = 2. * rk * (gamma4 + a1 * xmu)
      t4 = r4
      t5 = r5

      beta = -r5 / r4         
  
      e1 = min(rk * tau, 500.) 
      e2 = min(tau / xmu, 500.) 
      
      if (beam == 1) then
         den = r4 * exp(e1) + r5 * exp(-e1)
         ref  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         den = t4 * exp(e1) + t5 * exp(-e1)
         th  = exp(-e2)
         tra = th-th*w0*(t1*exp(e1)-t2*exp(-e1)-t3*exp(e2))/den
      elseif (beam == 2) then
         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         ref = (gamma2*(1.-ef2))/((rk+gamma1)*(1.-beta*ef2))
         tra = (2.*rk*ef1)/((rk+gamma1)*(1.-beta*ef2))
      endif
    end if
  end subroutine two_stream

  ! ########################################################################################
  elemental function two_stream_reflectance(tauint, gint, w0int)
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    
    ! INPUTS
    real, intent(in) :: tauint, gint, w0int
    ! OUTPUTS
    real             :: two_stream_reflectance
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real,parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den

    f   = gint**2
    tau = (1. - w0int * f) * tauint
    w0  = (1. - f) * w0int / (1. - w0int * f)
    g   = (gint - f) / (1. - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7. - w0* (4. + 3. * g)) / 4.
    gamma2 = -(1. - w0* (4. - 3. * g)) / 4.
    gamma3 =  (2. - 3.*g*xmu) / 4.
    gamma4 =   1. - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
          two_stream_reflectance = rh / (1. + gamma1 * tau)
      elseif (beam == 2) then
          two_stream_reflectance = gamma1*tau/(1. + gamma1*tau)
      endif
        
    else    !

        ! Non-conservative scattering
         a1 = gamma1 * gamma4 + gamma2 * gamma3
         a2 = gamma1 * gamma3 + gamma2 * gamma4

         rk = sqrt(gamma1**2 - gamma2**2)
         
         r1 = (1. - rk * xmu) * (a2 + rk * gamma3)
         r2 = (1. + rk * xmu) * (a2 - rk * gamma3)
         r3 = 2. * rk *(gamma3 - a2 * xmu)
         r4 = (1. - (rk * xmu)**2) * (rk + gamma1)
         r5 = (1. - (rk * xmu)**2) * (rk - gamma1)
         
         t1 = (1. + rk * xmu) * (a1 + rk * gamma4)
         t2 = (1. - rk * xmu) * (a1 - rk * gamma4)
         t3 = 2. * rk * (gamma4 + a1 * xmu)
         t4 = r4
         t5 = r5

         beta = -r5 / r4         
         
         e1 = min(rk * tau, 500.) 
         e2 = min(tau / xmu, 500.) 
         
         if (beam == 1) then
           den = r4 * exp(e1) + r5 * exp(-e1)
           two_stream_reflectance  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         elseif (beam == 2) then
           ef1 = exp(-e1)
           ef2 = exp(-2*e1)
           two_stream_reflectance = (gamma2*(1.-ef2))/((rk+gamma1)*(1.-beta*ef2))
         endif
           
      end if
  end function two_stream_reflectance 

  ! ########################################################################################
  pure subroutine adding_doubling (npts,Refl, Tran, Refl_tot, Tran_tot)      
    ! Use adding/doubling formulas to compute total reflectance and transmittance from 
    ! layer values
    
    ! INPUTS
    integer,intent(in)                  :: npts
    real,intent(in),dimension(npts) :: Refl,Tran
    ! OUTPUTS
    real,intent(out)                :: Refl_tot, Tran_tot
    ! LOCAL VARIABLES
    integer :: i
    real, dimension(npts) :: Refl_cumulative, Tran_cumulative
    
    Refl_cumulative(1) = Refl(1)
    Tran_cumulative(1) = Tran(1)    
    
    do i=2, npts
       ! place (add) previous combined layer(s) reflectance on top of layer i, w/black surface (or ignoring surface):
       Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i)*(Tran_cumulative(i-1)**2)/(1. - Refl_cumulative(i-1) * Refl(i))
       Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) / (1. - Refl_cumulative(i-1) * Refl(i))
    end do
    
    Refl_tot = Refl_cumulative(size(Refl))
    Tran_tot = Tran_cumulative(size(Refl))
    
  end subroutine adding_doubling

end program modis_polyfits

