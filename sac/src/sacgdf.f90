module sacgdf
  use hdf5



contains
  
  subroutine sac_gdf_read_eqpar(file_id, dimensionality)
    ! Convert simulation parameters to the eqpar array
    use common_variables
    use hdf5, only: HID_T, h5gopen_f, h5gclose_f
    use gdf_helpers, only: read_attribute

    implicit none

    integer(HID_T), intent(in) :: file_id
    integer, intent(in) :: dimensionality
    

    call h5gopen_f(file_id, 'simulation_parameters', g_id, error)
    call read_attribute(g_id, 'gamma', eqpar(gamma_))
    call read_attribute(g_id, 'eta', eqpar(eta_))
    call read_attribute(g_id, 'gravity0', eqpar(grav0_))
    call read_attribute(g_id, 'gravity1', eqpar(grav1_))
    if (dimensionality .GE. 2) then
       call read_attribute(g_id, 'gravity2', eqpar(grav2_))
    end if
    if (dimensionality .EQ. 3) then
       call read_attribute(g_id, 'gravity3', eqpar(grav3_))
    end if
    call read_attribute(g_id, 'nu', eqpar(nu_))
    call h5close_f(g_id, error)


  end subroutine simulation_params_eqpar_3D
  
end module sacgdf
