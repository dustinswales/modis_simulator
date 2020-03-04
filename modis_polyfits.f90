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
       phaseIsLiquid = 1,  & !
       phaseIsIce    = 2     !

  ! Testing parameters
  real, dimension(2), parameter :: &
       re_liq_range = (/1.,40./), & ! Range (microns) to compute optical properties for
       re_ice_range = (/1.,100./)
  
  ! Local varaibles
  integer :: ij, nRe_liq, nRe_ice, ncid, status, dimID1, dimID2, varID(10)
  real,dimension(re_liq_range(2)-re_liq_range(1)+1) :: &
       g_liq_old, g_liq_new, ssa_liq_old, ssa_liq_new, re_liq
  real,dimension(re_ice_range(2)-re_ice_range(1)+1) :: &
       g_ice_old, g_ice_new, ssa_ice_old, ssa_ice_new, re_ice

  nRe_liq = re_liq_range(2) - re_liq_range(1)+1
  nRe_ice = re_ice_range(2) - re_ice_range(1)+1
  
  ! Liquid
  do ij=re_liq_range(1),re_liq_range(2)
     re_liq(ij)      = ij
     g_liq_old(ij)   = get_g_nir_old(  phaseIsLiquid,real(ij))
     g_liq_new(ij)   = get_g_nir_new(  phaseIsLiquid,real(ij))
     ssa_liq_old(ij) = get_ssa_nir_old(phaseIsLiquid,real(ij))
     ssa_liq_new(ij) = get_ssa_nir_new(phaseIsLiquid,real(ij))
  enddo
  ! Ice
  do ij=re_ice_range(1),re_ice_range(2)
     re_ice(ij)      = ij
     g_ice_old(ij)   = get_g_nir_old(  phaseIsIce,real(ij))
     g_ice_new(ij)   = get_g_nir_new(  phaseIsIce,real(ij))
     ssa_ice_old(ij) = get_ssa_nir_old(phaseIsIce,real(ij))
     ssa_ice_new(ij) = get_ssa_nir_new(phaseIsIce,real(ij))
  enddo

  ! Write to outputs
  status = nf90_create("modis_polyfits.nc",nf90_clobber,ncid)
  status = nf90_def_dim(ncid,'nre_liq',nRe_liq,dimID1)
  status = nf90_def_dim(ncid,'nre_ice',nRe_ice,dimID2)
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
  status = nf90_close(ncid)
  
contains
  
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


end program modis_polyfits

