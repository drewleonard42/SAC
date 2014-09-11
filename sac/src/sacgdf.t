module sacgdf
  use hdf5

contains
  
  subroutine sacgdf_write_file(file_id, gdf_rd, gdf_sp, field_types)
    use hdf5, only: HID_T
    use gdf, only: gdf_write_file, gdf_root_datasets_T, gdf_parameters_T, gdf_field_type_T
    implicit none
    
    integer(HID_T), intent(inout) :: file_id
    type(gdf_root_datasets_T), intent(inout) :: gdf_rd
    type(gdf_parameters_T), intent(inout) :: gdf_sp
    type(gdf_field_type_T), dimension(:), intent(inout) :: field_types

    call gdf_write_file(file_id, "Sheffield Advanced Code", "GDF Testing version", &
         gdf_rd, gdf_sp, field_types)
    
    call sacgdf_write_eqpar(file_id)
    
  end subroutine sacgdf_write_file


  subroutine sacgdf_read_file(file_id, software_name, software_version, &
       gdf_rd, gdf_sp, field_types)
    use hdf5, only: HID_T
    use gdf, only: gdf_read_file, gdf_root_datasets_T, gdf_parameters_T, gdf_field_type_T
    implicit none
    
    integer(HID_T), intent(inout) :: file_id
    type(gdf_root_datasets_T), intent(inout) :: gdf_rd
    type(gdf_parameters_T), intent(inout) :: gdf_sp
    type(gdf_field_type_T), dimension(:), allocatable, intent(inout) :: field_types
    character(len=*), intent(out) :: software_name, software_version
    
    call gdf_read_file(file_id, software_name, software_version, &
         gdf_rd, gdf_sp, field_types)

    call sacgdf_read_eqpar(file_id)
    
  end subroutine sacgdf_read_file





  subroutine sacgdf_write_eqpar(file_id)
    ! Convert simulation parameters to the eqpar array
    use common_variables
    use hdf5, only: HID_T, h5gopen_f, h5gclose_f
    use helpers_hdf5, only: create_attribute

    implicit none

    integer(HID_T), intent(in) :: file_id
    
    integer(HID_T) :: g_id
    integer :: error
    real(kind=8), dimension(:), pointer :: r_ptr
    integer(kind=4), dimension(:), pointer :: i4_ptr
    integer(kind=4), dimension(1), target :: it_arr

    call h5gopen_f(file_id, 'simulation_parameters', g_id, error)

    r_ptr => eqpar(gamma_:gamma_)
    call create_attribute(g_id, 'gamma', r_ptr)

    r_ptr => eqpar(eta_:eta_)
    call create_attribute(g_id, 'eta', r_ptr)

    r_ptr => eqpar(grav0_:grav0_)
    call create_attribute(g_id, 'gravity0', r_ptr)

    r_ptr => eqpar(grav1_:grav1_)
    call create_attribute(g_id, 'gravity1', r_ptr)

    ! Read the extra parameters only if we are 2D or 3D
    {^IFTWOD
    r_ptr => eqpar(grav2_:grav2_)
    call create_attribute(g_id, 'gravity2', r_ptr)

    }
    {^IFTHREED
    r_ptr => eqpar(grav2_:grav2_)
    call create_attribute(g_id, 'gravity2', r_ptr
    r_ptr => eqpar(grav3_:grav3_)
    call create_attribute(g_id, 'gravity3', 
    }

    r_ptr => eqpar(nu_:nu_)
    call create_attribute(g_id, 'nu', r_ptr)

    it_arr(1) = it
    i4_ptr => it_arr
    call create_attribute(g_id, 'current_iteration', i4_ptr)

    call h5gclose_f(g_id, error)

  end subroutine sacgdf_write_eqpar

  subroutine sacgdf_make_field_types(field_types)
    
    use common_variables, only: nw
    use gdf, only: gdf_field_type_T
    implicit none

    type(gdf_field_type_T), dimension(nw), intent(inout) :: field_types

    select case (nw)
    case (7)
       {^IFONED call sacgdf_make_field_types_1D(field_types) }
    case (10)
       {^IFTWOD call sacgdf_make_field_types_2D(field_types) }
    case (13)
       {^IFTHREED call sacgdf_make_field_types_3D(field_types) }
    case default
       call die("Wrong nw")
    end select
    
  end subroutine sacgdf_make_field_types
       

 subroutine sacgdf_make_field_types_1D(field_types)

    use gdf, only: gdf_field_type_T
    implicit none

    type(gdf_field_type_T), dimension(7), intent(inout) :: field_types
    
    ! Write Field Type information
    call field_types(1)%init()
    field_types(1)%field_to_cgs = 0.001
    field_types(1)%staggering = 0
    field_types(1)%field_units = "kg m^{-3}"
    field_types(1)%variable_name = "density_pert"
    field_types(1)%field_name = "Pertubation Density"

    call field_types(2)%init()
    field_types(2)%field_to_cgs = 0.001
    field_types(2)%staggering = 0
    field_types(2)%field_units = "kg m^{-3}"
    field_types(2)%variable_name = "density_bg"
    field_types(2)%field_name = "Background Density"

    call field_types(3)%init()
    field_types(3)%field_to_cgs = 100
    field_types(3)%staggering = 0
    field_types(3)%field_units = "ms^{-1}"
    field_types(3)%variable_name = "velocity_x"
    field_types(3)%field_name = "x Component of Velocity"

    call field_types(4)%init()
    field_types(4)%field_to_cgs = 10
    field_types(4)%staggering = 0
    field_types(4)%field_units = "Pa"
    field_types(4)%variable_name = "internal_energy_pert"
    field_types(4)%field_name = "Pertubation Internal Energy"

    call field_types(5)%init()
    field_types(5)%field_to_cgs = 10
    field_types(5)%staggering = 0
    field_types(5)%field_units = "Pa"
    field_types(5)%variable_name = "internal_energy_bg"
    field_types(5)%field_name = "Background Internal Energy"

    call field_types(6)%init()
    field_types(6)%field_to_cgs = 10000.0
    field_types(6)%staggering = 0
    field_types(6)%field_units = "T"
    field_types(6)%variable_name = "mag_field_x_bg"
    field_types(6)%field_name = "x Component of Background Magnetic Field"

    call field_types(7)%init()
    field_types(7)%field_to_cgs = 10000.0
    field_types(7)%staggering = 0
    field_types(7)%field_units = "T"
    field_types(7)%variable_name = "mag_field_x_pert"
    field_types(7)%field_name = "x Component of Pertubation Magnetic Field"

  end subroutine sacgdf_make_field_types_1D

  subroutine sacgdf_make_field_types_2D(field_types)

    use gdf, only: gdf_field_type_T
    use common_variables, only: nw
    implicit none

    type(gdf_field_type_T), dimension(10), intent(inout) :: field_types

    call sacgdf_make_field_types_1D(field_types(1:7))

    call field_types(8)%init()
    field_types(8)%field_to_cgs = 100
    field_types(8)%staggering = 0
    field_types(8)%field_units = "ms^{-1}"
    field_types(8)%variable_name = "velocity_y"
    field_types(8)%field_name = "y Component of Velocity"

    call field_types(9)%init()
    field_types(9)%field_to_cgs = 10000.0
    field_types(9)%staggering = 0
    field_types(9)%field_units = "T"
    field_types(9)%variable_name = "mag_field_y_bg"
    field_types(9)%field_name = "y Component of Background Magnetic Field"

    call field_types(10)%init()
    field_types(10)%field_to_cgs = 10000.0
    field_types(10)%staggering = 0
    field_types(10)%field_units = "T"
    field_types(10)%variable_name = "mag_field_y_pert"
    field_types(10)%field_name = "y Component of Pertubation Magnetic Field"

  end subroutine sacgdf_make_field_types_2D

  subroutine sacgdf_make_field_types_3D(field_types)

    use gdf, only: gdf_field_type_T
    use common_variables, only: nw
    implicit none

    type(gdf_field_type_T), dimension(13), intent(inout) :: field_types

    call sacgdf_make_field_types_2D(field_types(1:10))

    call field_types(11)%init()
    field_types(11)%field_to_cgs = 10000.0
    field_types(11)%staggering = 0
    field_types(11)%field_units = "T"
    field_types(11)%variable_name = "mag_field_z_pert"
    field_types(11)%field_name = "z Component of Pertubation Magnetic Field"

    call field_types(12)%init()
    field_types(12)%field_to_cgs = 10000.0
    field_types(12)%staggering = 0
    field_types(12)%field_units = "T"
    field_types(12)%variable_name = "mag_field_z_bg"
    field_types(12)%field_name = "z Component of Background Magnetic Field"

    call field_types(13)%init()
    field_types(13)%field_to_cgs = 100
    field_types(13)%staggering = 0
    field_types(13)%field_units = "ms^{-1}"
    field_types(13)%variable_name = "velocity_z"
    field_types(13)%field_name = "z Component of Velocity"

  end subroutine sacgdf_make_field_types_3D
  
  subroutine sacgdf_read_eqpar(file_id)
    ! Convert simulation parameters to the eqpar array
    use common_variables
    use hdf5, only: HID_T, h5gopen_f, h5gclose_f
    use helpers_hdf5, only: read_attribute

    implicit none

    integer(HID_T), intent(in) :: file_id
    
    integer(HID_T) :: g_id
    integer :: error

    real(kind=8), dimension(:), pointer :: r_ptr
    integer(kind=4), dimension(:), pointer :: i4_ptr
    integer(kind=4), dimension(1), target :: it_arr
    

    call h5gopen_f(file_id, 'simulation_parameters', g_id, error)
    r_ptr => eqpar(gamma_:gamma_)
    call read_attribute(g_id, 'gamma', r_ptr)

    r_ptr => eqpar(eta_:eta_)
    call read_attribute(g_id, 'eta', r_ptr)

    r_ptr => eqpar(grav0_:grav0_)
    call read_attribute(g_id, 'gravity0', r_ptr)

    r_ptr => eqpar(grav1_:grav1_)
    call read_attribute(g_id, 'gravity1', r_ptr)

    ! Read the extra parameters only if we are 2D or 3D
    {^IFTWOD
    r_ptr => eqpar(grav2_:grav2_)
    call read_attribute(g_id, 'gravity2', r_ptr)
    }
    {^IFTHREED
    r_ptr => eqpar(grav2_:grav2_)
    call read_attribute(g_id, 'gravity2', r_ptr)
    r_ptr => eqpar(grav3_:grav3_)
    call read_attribute(g_id, 'gravity3', r_ptr)
    }

    r_ptr => eqpar(nu_:nu_)
    call read_attribute(g_id, 'nu', r_ptr)

    i4_ptr => it_arr
    call read_attribute(g_id, 'current_iteration', i4_ptr)
    it = it_arr(1)

    call h5gclose_f(g_id, error)

  end subroutine sacgdf_read_eqpar

  subroutine build_x_array(ix^L, nx, left_edge, right_edge, x)
    ! Construct the x array (which is cell centred) from the left_edge and right_edge
    ! values, which are also assumed to be cell centred.
    implicit none

    integer, intent(in) :: ix^L
    integer, dimension(^ND), intent(in) :: nx
    real(kind=8), dimension(^ND), intent(in) :: left_edge, right_edge
    real(kind=8), dimension(:^D&,:), intent(inout) :: x

    integer :: ix^D
    real(kind=8) :: dx^D
    
!!$    ! Set ixmin = 1
!!$    ixmin^D=1;
!!$    ! Set ixmax to ixmin+nx
!!$    ixmax^D=ixmin^D+nx(^D)-1;

    dx^D=(right_edge(^D)-left_edge(^D))/nx(^D);

    {
    forall(ix^DD=ixmin^DD:ixmax^DD)
       x(ix^DD,^D)= ((ix^D-ixmin^D)*right_edge(^D) + (ixmax^D-ix^D)*left_edge(^D)) / (ixmax^D-ixmin^D) 
    end forall
    \}

  end subroutine build_x_array

  subroutine sacgdf_read_datasets(place, plist_id, w, ix^L)
    use hdf5, only: HID_T
    use gdf_datasets
    use common_variables, only: ixGhi^D, ixGlo^D, nw

    implicit none

    integer(HID_T), intent(inout) :: place
    integer(HID_T), intent(inout) :: plist_id         !< Property list identifier
    integer, intent(in) :: ix^L
    real(kind=8), intent(inout):: w(ixG^T,nw)


    call sacgdf_read_datasets_1D(place, plist_id, w, ix^L)
    {^IFTWOD   call sacgdf_read_datasets_2D(place, plist_id, w, ix^L) }
    {^IFTHREED call sacgdf_read_datasets_3D(place, plist_id, w, ix^L) }

  end subroutine sacgdf_read_datasets

  subroutine sacgdf_read_datasets_1D(place, plist_id, w, ix^L)
    use hdf5, only: HID_T
    use gdf_datasets, only: read_real8_dataset
    use common_variables, only: rho_, rhob_, e_, eb_, m1_, b1_, bg1_
    use common_variables, only: ixGhi^D, ixGlo^D, nw, nx

    implicit none

    integer(HID_T), intent(inout) :: place
    integer(HID_T), intent(inout) :: plist_id         !< Property list identifier
    integer, intent(in) :: ix^L
    real(kind=8), intent(inout) :: w(ixG^T,nw)

    real(kind=8), dimension(:, :, :), allocatable :: wdata3D

    ! Density pert
    call read_real8_dataset(place, 'density_pert', wdata3D)
    w(ix^S, rho_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Denisty bg
    call read_real8_dataset(place, 'density_bg', wdata3D)
    w(ix^S, rhob_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Velocity
    call read_real8_dataset(place, 'velocity_x', wdata3D)
    w(ix^S, m1_) = reshape(wdata3D, nx) * (w(ix^S, rho_) + w(ix^S, rhob_))
    deallocate(wdata3D)

    ! internal energy pert
    call read_real8_dataset(place, 'internal_energy_pert', wdata3D)
    w(ix^S, e_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! internal energy bg
    call read_real8_dataset(place, 'internal_energy_bg', wdata3D)
    w(ix^S, eb_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Mag field pert
    call read_real8_dataset(place, 'mag_field_x_pert', wdata3D)
    w(ix^S, b1_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Mag field bg
    call read_real8_dataset(place, 'mag_field_x_bg', wdata3D)
    w(ix^S, bg1_) = reshape(wdata3D, nx)
    deallocate(wdata3D)
    
  end subroutine sacgdf_read_datasets_1D

{^IFTWOD
  subroutine sacgdf_read_datasets_2D(place, plist_id, w, ix^L)
    use hdf5, only: HID_T
    use gdf_datasets, only: read_real8_dataset
    use common_variables, only: rho_, rhob_, m2_, b2_, bg2_
    use common_variables, only: ixGhi^D, ixGlo^D, nw, nx

    implicit none

    integer(HID_T), intent(inout) :: place
    integer(HID_T), intent(inout) :: plist_id         !< Property list identifier
    integer, intent(in) :: ix^L
    real(kind=8), intent(inout) :: w(ixG^T,nw)

    real(kind=8), dimension(:, :, :), allocatable :: wdata3D

    ! Velocity
    call read_real8_dataset(place, 'velocity_y', wdata3D)
    w(ix^S, m2_) = reshape(wdata3D, nx) * (w(ix^S, rho_) + w(ix^S, rhob_))
    deallocate(wdata3D)

     ! Mag field pert
    call read_real8_dataset(place, 'mag_field_y_pert', wdata3D)
    w(ix^S, b2_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Mag field bg
    call read_real8_dataset(place, 'mag_field_y_bg', wdata3D)
    w(ix^S, bg2_) = reshape(wdata3D, nx)
    deallocate(wdata3D)
    
  end subroutine sacgdf_read_datasets_2D
}
{^IFTHREED
  subroutine sacgdf_read_datasets_3D(place, plist_id, w, ix^L)
    use hdf5, only: HID_T
    use gdf_datasets, only: read_real8_dataset
    use common_variables, only: rho_, rhob_, m3_, b3_, bg3_
    use common_variables, only: ixGhi^D, ixGlo^D, nw, nx

    implicit none

    integer(HID_T), intent(inout) :: place
    integer(HID_T), intent(inout) :: plist_id         !< Property list identifier
    integer, intent(in) :: ix^L
    real(kind=8), intent(inout) :: w(ixG^T,nw)

    real(kind=8), dimension(:, :, :), allocatable :: wdata3D

    ! Velocity
    call read_real8_dataset(place, 'velocity_y', wdata3D)
    w(ix^S, m2_) = reshape(wdata3D, nx) * (w(ix^S, rho_) + w(ix^S, rhob_))
    deallocate(wdata3D)

     ! Mag field pert
    call read_real8_dataset(place, 'mag_field_y_pert', wdata3D)
    w(ix^S, b2_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Mag field bg
    call read_real8_dataset(place, 'mag_field_y_bg', wdata3D)
    w(ix^S, bg2_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Velocity
    call read_real8_dataset(place, 'velocity_z', wdata3D)
    w(ix^S, m3_) = reshape(wdata3D, nx) * (w(ix^S, rho_) + w(ix^S, rhob_))
    deallocate(wdata3D)

     ! Mag field pert
    call read_real8_dataset(place, 'mag_field_z_pert', wdata3D)
    w(ix^S, b3_) = reshape(wdata3D, nx)
    deallocate(wdata3D)

    ! Mag field bg
    call read_real8_dataset(place, 'mag_field_z_bg', wdata3D)
    w(ix^S, bg3_) = reshape(wdata3D, nx)
    deallocate(wdata3D)
    
  end subroutine sacgdf_read_datasets_3D
}
end module sacgdf
