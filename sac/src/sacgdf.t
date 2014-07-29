module sacgdf
  use hdf5



contains

  subroutine sacgdf_write_eqpar(file_id, dimensionality)
    ! Convert simulation parameters to the eqpar array
    use common_variables
    use hdf5, only: HID_T, h5gopen_f, h5gclose_f
    use helpers_hdf5, only: create_attribute

    implicit none

    integer(HID_T), intent(in) :: file_id
    integer, intent(in) :: dimensionality
    
    integer(HID_T) :: g_id
    integer :: error

    call h5gopen_f(file_id, 'simulation_parameters', g_id, error)
    call create_attribute(g_id, 'gamma', eqpar(gamma_))
    call create_attribute(g_id, 'eta', eqpar(eta_))
    call create_attribute(g_id, 'gravity0', eqpar(grav0_))
    call create_attribute(g_id, 'gravity1', eqpar(grav1_))
    ! Read the extra parameters only if we are 2D or 3D
    {^IFTWOD
    call create_attribute(g_id, 'gravity2', eqpar(grav2_))
    }
    {^IFTHREED
    call create_attribute(g_id, 'gravity3', eqpar(grav3_))
    }
    call create_attribute(g_id, 'nu', eqpar(nu_))
    call h5close_f(g_id, error)

  end subroutine sacgdf_write_eqpar
  
  subroutine sacgdf_read_eqpar(file_id, dimensionality)
    ! Convert simulation parameters to the eqpar array
    use common_variables
    use hdf5, only: HID_T, h5gopen_f, h5gclose_f
    use helpers_hdf5, only: read_attribute

    implicit none

    integer(HID_T), intent(in) :: file_id
    integer, intent(in) :: dimensionality
    
    integer(HID_T) :: g_id
    integer :: error
    

    call h5gopen_f(file_id, 'simulation_parameters', g_id, error)
    call read_attribute(g_id, 'gamma', eqpar(gamma_))
    call read_attribute(g_id, 'eta', eqpar(eta_))
    call read_attribute(g_id, 'gravity0', eqpar(grav0_))
    call read_attribute(g_id, 'gravity1', eqpar(grav1_))
    ! Read the extra parameters only if we are 2D or 3D
    {^IFTWOD
    call read_attribute(g_id, 'gravity2', eqpar(grav2_))
    }
    {^IFTHREED
    call read_attribute(g_id, 'gravity3', eqpar(grav3_))
    }
    call read_attribute(g_id, 'nu', eqpar(nu_))
    call h5close_f(g_id, error)

  end subroutine sacgdf_read_eqpar
  
end module sacgdf
