PROGRAM testread
  USE MPI
  USE HDF5
  IMPLICIT NONE
  !! SAC File Params 
  INTEGER :: unitini = 10
  INTEGER, PARAMETER :: ndim = 3, nw = 13, ixmin1=1, ixmin2=1, ixmin3=1 !not par
  INTEGER, PARAMETER :: ixmax1=64, ixmax2=64, ixmax3=32 !not par
  !! Read in varibles
  CHARACTER*79 :: filehead,varnames
  INTEGER :: it, neqpar, ios, nx(ndim)
  INTEGER :: ieqpar, neqparin, eqparextra, i, cnt
  REAL(8) :: eqpar(7), w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,nw)
  REAL(8) :: x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,ndim)
  REAL(8) :: t

  CHARACTER(len=4) :: fproc,ext
  CHARACTER(len=400) :: fname = "/fastdata/smq11sjm/3D_data/3D_tube128_Slog_p240_C_A10_B5_np020204.out"
  CHARACTER(len=400) :: outfname = "/fastdata/smq11sjm/3D_data/3D_tube128_Slog_p240_C_A10_B5_600s_np020204.h5"
  CHARACTER(len=400) :: rank_fname
  CHARACTER(len=3*2+2) :: np_str
  CHARACTER(len=10) :: wname

  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: dset_id, hslab_id      ! Dataset identifier 
  INTEGER(HID_T) :: xf_space, wf_space     ! Dataspace identifier in file 
  INTEGER(HID_T) :: xm_space, wm_space      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id,splist_id,gplist_id, aplist_id    ! Property list identifier 

  INTEGER(HSIZE_T),  DIMENSION(4) :: xcount, wcount  
  INTEGER(HSSIZE_T), DIMENSION(4) :: offset 
  INTEGER(HSIZE_T),  DIMENSION(4) :: stride
  INTEGER(HSIZE_T),  DIMENSION(4) :: xblock, wblock

  INTEGER(HID_T) :: file, sac_group, time_group, xspace, xset, wspace, wset
  INTEGER(HSIZE_T), DIMENSION(4) :: xdims = (/128,128,128,ndim/), np_xdims = (/ixmax1,ixmax2,ixmax3,3/)
  INTEGER(HSIZE_T), DIMENSION(4) :: wdims = (/128,128,128,nw/), np_wdims = (/ixmax1,ixmax2,ixmax3,nw/)

  !All attr varibles
  INTEGER(HID_T) :: attr_id       ! Current Attribute identifier
  INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
  ! String attr varibles
  INTEGER(HID_T) :: aspace_str, astr_type
  INTEGER(HSIZE_T), DIMENSION(1) :: astr_dims = (/1/)
  INTEGER     ::   astr_rank = 1
  INTEGER(SIZE_T) :: astr_len = 80
  CHARACTER(LEN=80), DIMENSION(1) ::  astr_data
  !Set up attribute dataspace for string integer and float
  INTEGER(HSIZE_T), DIMENSION(1) :: aint_dims = (/1/)
  INTEGER     ::   aint_rank = 1
  INTEGER(HID_T) ::  aspace_int
  !Set up attribute dataspace for string integer and float
  INTEGER(HSIZE_T), DIMENSION(1) :: adbl_dims = (/1/)
  INTEGER     ::   adbl_rank = 1
  REAL(8), DIMENSION(1) ::  adbl_data
  INTEGER(HID_T) :: aspace_dbl
  !Set up attribute dataspace for nx
  INTEGER(HSIZE_T), DIMENSION(1) :: anx_dims = (/ndim/)
  INTEGER(HID_T) ::  aspace_nx
  !Set up attribute dataspace for eqpar
  INTEGER(HSIZE_T), DIMENSION(1) :: aeq_dims = (/7/)
  INTEGER(HID_T) :: aspace_eqpar

  INTEGER :: ierr, rank, size,np1,np2,np3,name_len,error,fierr,comm,info,coords(3), rank2

  !Initilise MPI communications
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(comm, rank, ierr)
  CALL MPI_COMM_SIZE(comm, size, ierr)
  !Give each process a unique file id
  unitini = unitini + rank
  !Extract the file extension
  name_len = LEN_TRIM(fname)
  ext = fname(name_len-3:name_len)
  !Extract the grid size from the filenames
  np_str = fname(name_len-3-3*2:name_len-4)
  READ (np_str, "(I2,I2,I2)"), np1,np2,np3
  WRITE(fproc, "(A1,I0.3)") "_", rank
  rank_fname = fname(:name_len-4)//fproc//ext
  PRINT*, rank, "opens", rank_fname
  !Open the binary files and read in
  OPEN(unitini, file=rank_fname, status='old', form='unformatted')
  !Read first step
  CALL read_sac_step(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,unitini, ndim, nw,&
       filehead,varnames,it, t, ndim, neqpar, nw, nx,ieqpar, neqparin, eqparextra,eqpar, w, x, ios)




  !Determine the processor location in the grid and calculate offsets etc.
  stride(:) = 1 
  xcount(:) =  (/ixmax1,ixmax2,ixmax3,ndim/)
  xblock = (/np1,np2,np3,1/)  

  coords(:) = 0
  IF ( rank .LT. np1) THEN
     coords(1) = rank
  ELSE IF (rank .EQ. np1) THEN
     coords(2) = 1
  ELSE IF (rank .LT. np1 + np2) THEN
     coords(1) = rank - np1
     coords(2) = rank - np2
  ELSE
     coords(3) = rank / (np1 * np2)
     rank2 = rank - (coords(3) * np3)
     IF ( rank2 .LT. np1) THEN
        coords(1) = rank2
     ELSE IF (rank2 .EQ. np1) THEN
        coords(2) = 1
     ELSE IF (rank2 .LT. np1 + np2) THEN
        coords(1) = rank2 - np1
        coords(2) = rank2 - np2
     END IF
  END IF

  offset(1:3) = coords * nx
  offset(4) = 0


  !
  ! Initialize HDF5 library and Fortran interfaces.
  !
  CALL h5open_f(error) 

  ! 
  ! Setup file access property list with parallel I/O access.
  !
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

  CALL h5pcreate_f(H5P_GROUP_ACCESS_F, gplist_id, error)
  CALL h5pcreate_f(H5P_ATTRIBUTE_CREATE_F, aplist_id, error)

  !
  ! Create the file collectively.
  ! 
  CALL h5fcreate_f(TRIM(outfname), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! DEFINE Datatypes and spaces for eveything !!!!!!!!!!!
!!!!! Single (str, int, dbl), and x and w arrays !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Create dataspace and 80len character datatype
  CALL h5screate_simple_f(astr_rank, astr_dims, aspace_str, error)
  CALL h5tcopy_f(H5T_NATIVE_CHARACTER, astr_type, error)
  CALL h5tset_size_f(astr_type, astr_len, error)

  !Create Dataspace for int
  CALL h5screate_simple_f(aint_rank, aint_dims, aspace_int, error)
  !Create Dataspace for dbl
  CALL h5screate_simple_f(adbl_rank, adbl_dims, aspace_dbl, error)
  !Create Dataspace for nx
  CALL h5screate_simple_f(aint_rank, anx_dims, aspace_nx, error)
  !Create Dataspace for dbl
  CALL h5screate_simple_f(adbl_rank, aeq_dims, aspace_eqpar, error)
  !
  ! Create the data space for the  dataset. 
  !
  CALL h5screate_simple_f(ndim+1, xdims, xf_space, error)
  CALL h5screate_simple_f(ndim+1, np_xdims, xm_space, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Write the filehead to a file attribute
  CALL h5acreate_f(file_id, "filehead", astr_type, aspace_str, attr_id, error, acpl_id = aplist_id)
  astr_data(1) = filehead !Setdata
  CALL h5awrite_f(attr_id, astr_type, astr_data, data_dims, error)
  CALL h5aclose_f(attr_id, error)

  !Write A description of the file to a file attribute
  CALL h5acreate_f(file_id, "filedesc", astr_type, aspace_str, attr_id, fierr)
  astr_data(1) = "This is the initial conditions for a self-similar flux tube"
  CALL h5awrite_f(attr_id, astr_type, astr_data, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Create a group!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL h5gcreate_f(file_id, "/SACdata", sac_group, error, gapl_id = gplist_id)
  !Write attributes to main group
  !neqpar, eqpar, nt, its, ndim
  !neqpar
  CALL h5acreate_f(sac_group, "neqpar", H5T_NATIVE_INTEGER, aspace_int, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, neqpar, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)
  !eqpar
  CALL h5acreate_f(sac_group, "eqpar", H5T_NATIVE_DOUBLE, aspace_eqpar, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, eqpar, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)

  !ndim
  CALL h5acreate_f(sac_group, "ndim", H5T_NATIVE_INTEGER, aspace_int, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, ndim, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)
  !nx
  CALL h5acreate_f(sac_group, "nx", H5T_NATIVE_INTEGER, aspace_nx, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, nx, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)
  !
  ! Create chunked dataset.
  !
  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  CALL h5pset_chunk_f(plist_id, ndim+1, np_xdims, error)
  CALL h5dcreate_f(sac_group, "x", H5T_NATIVE_DOUBLE, xf_space, dset_id, error, plist_id)
  CALL h5sclose_f(xf_space, error)
  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !

  ! 
  ! Select hyperslab in the file.
  !
  CALL h5dget_space_f(dset_id, xf_space, error)
  CALL h5sselect_hyperslab_f (xf_space, H5S_SELECT_SET_F, offset, xcount, error)

  !
  ! Create property list for collective dataset write
  !
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  !
  ! Write X dataset collectively. 
  !
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, np_xdims, error, &
       file_space_id = xf_space, mem_space_id = xm_space, xfer_prp = plist_id)
  !
  ! Close dataspaces.
  !
  CALL h5sclose_f(xf_space, error)
  CALL h5sclose_f(xm_space, error)
  !
  ! Close the dataset.
  !
  CALL h5dclose_f(dset_id, error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DO W bit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Create a group!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL h5gcreate_f(sac_group, "wseries", time_group, error, gapl_id = gplist_id)

  !Write varnames and nw to time series group
  CALL h5acreate_f(time_group, "varnames", astr_type, aspace_str, attr_id, fierr, acpl_id = aplist_id)
  astr_data(1) = varnames
  CALL h5awrite_f(attr_id, astr_type, astr_data, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)

  CALL h5acreate_f(time_group, "nw", H5T_NATIVE_INTEGER, aspace_int, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, nw, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)


  !Create a dataspace for w on this processor and for one varible
  !
  ! Create the data space for the  dataset. 
  !
  CALL h5screate_simple_f(ndim+1, wdims, wf_space, error)
  CALL h5screate_simple_f(ndim+1, np_wdims, wm_space, error)

  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  wcount(:) =  (/64,64,32,nw/)
  wblock = (/np1,np2,np3,1/)  

  CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  CALL h5pset_chunk_f(plist_id, ndim+1, np_wdims, error)

  
  !
  ! Create property list for collective dataset write
  !
  CALL h5pcreate_f(H5P_DATASET_XFER_F, splist_id, error) 
  CALL h5pset_dxpl_mpio_f(splist_id, H5FD_MPIO_COLLECTIVE_F, error)
  
  cnt = 1
  DO i=1,1000
     IF (rank .EQ. 0)  PRINT*, i, cnt

     !Stop if EOF
     IF (ios<0) EXIT

     !Generate a name for the dataset
     WRITE (wname, "(A1,I0.5)") "w", cnt

     !
     ! Create chunked dataset.
     !

     CALL h5dcreate_f(time_group, wname, H5T_NATIVE_DOUBLE, wf_space, dset_id, error, plist_id)
     CALL h5sclose_f(wf_space, error)

     ! 
     ! Select hyperslab in the file.
     !
     CALL h5dget_space_f(dset_id, wf_space, error)
     CALL h5sselect_hyperslab_f (wf_space, H5S_SELECT_SET_F, offset, wcount, error)
     
     !
     ! Write W dataset collectively. 
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, w, np_wdims, error, &
          file_space_id = wf_space, mem_space_id = wm_space, xfer_prp = splist_id)

     !Write it and t to dataset
     !t
     CALL h5acreate_f(dset_id, "t", H5T_NATIVE_DOUBLE, aspace_dbl, attr_id, fierr, acpl_id = aplist_id)
     CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, t, data_dims, fierr)
     CALL h5aclose_f(attr_id, error)
     !iter
     CALL h5acreate_f(dset_id, "it", H5T_NATIVE_INTEGER, aspace_int, attr_id, fierr, acpl_id = aplist_id)
     CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, it, data_dims, fierr)
     CALL h5aclose_f(attr_id, error)



     cnt = cnt + 1
     CALL read_sac_step(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,unitini, ndim, nw,&
          filehead,varnames,it, t, ndim, neqpar, nw, nx,ieqpar, neqparin, eqparextra,eqpar, w, x, ios)
     !
     ! Close the dataset.
     !
     CALL h5dclose_f(dset_id, error)
  END DO

  !Write final time and number of iterations in simulation to top group.
  !nt
  CALL h5acreate_f(sac_group, "final t", H5T_NATIVE_DOUBLE, aspace_dbl, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, t, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)
  !nits
  CALL h5acreate_f(sac_group, "nt", H5T_NATIVE_INTEGER, aspace_int, attr_id, fierr, acpl_id = aplist_id)
  CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, it, data_dims, fierr)
  CALL h5aclose_f(attr_id, error)


  !
  ! Close dataspaces.
  !
  CALL h5sclose_f(wf_space, error)
  CALL h5sclose_f(wm_space, error)
  !
  ! Close the dataset.
  !
  !CALL h5dclose_f(dset_id, error)

  CLOSE(unitini)
  !
  ! Close the property lists.
  !
  CALL h5pclose_f(plist_id, error)
  CALL h5pclose_f(gplist_id, error)
  CALL h5pclose_f(aplist_id, error)
  !
  !Close the groups
  !
  CALL h5gclose_f(sac_group, error)
  CALL h5gclose_f(time_group,error)
  !
  ! Close the attribute dataspaces
  !
  !  CALL h5aclose_f(attr_id, error)
  CALL h5sclose_f(aspace_str, error)
  CALL h5sclose_f(aspace_int, error)
  CALL h5sclose_f(aspace_dbl, error)
  CALL h5sclose_f(aspace_nx, error)
  CALL h5sclose_f(aspace_eqpar, error)
  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id, error)

  !
  ! Close FORTRAN interfaces and HDF5 library.
  !
  CALL h5close_f(error)



  CALL MPI_FINALIZE(ierr)
  PRINT*, "Terminated:", rank
END PROGRAM Testread

SUBROUTINE read_sac_step(ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,unitini, ndim, nw,&
fileheadini,varnamesini,it, t, ndimini, neqparini, nwini, nx,ieqpar, neqparin, eqparextra,eqpar, w, x, ios)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
  INTEGER, INTENT(IN) :: unitini, ndim, nw
  
  CHARACTER*79, INTENT(OUT) :: fileheadini,varnamesini
  INTEGER, INTENT(OUT) :: it, ndimini, neqparini, nwini, nx(ndim)
  INTEGER, INTENT(OUT) :: ieqpar, neqparin, eqparextra, ios
  REAL(8), INTENT(OUT) :: eqpar(7), w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,nw)
  REAL(8), INTENT(OUT) :: x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,ndim), t
  INTEGER :: snapshot
  INTEGER :: idim, iw
  
  READ(unitini,iostat=ios) fileheadini
  IF(ios<0) RETURN                ! Cycle until the last recorded state
  READ(unitini,iostat=ios) it, t, ndimini, neqparini, nwini

  READ(unitini,iostat=ios) nx
!!$    READ(unitini,iostat=ios) (eqpar(ieqpar), ieqpar=1, neqparin),&
!!$       (eqparextra, ieqpar = neqparin+1, neqparini)
  READ(unitini,iostat=ios) eqpar

  READ(unitini,iostat=ios) varnamesini

  READ(unitini,iostat=ios) (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       idim), idim = 1, ndim)

  ! To conform savefileout_bin we use loop for iw
  DO iw=1,nw
     READ(unitini,iostat=ios) w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,iw)
  END DO
  IF(ios/=0)THEN
     PRINT*, "ERROR in binary file read:", ios, unitini
  END IF


END SUBROUTINE read_sac_step
