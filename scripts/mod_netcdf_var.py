def mod_netcdf_var(file_name,var_name,new_values):

  import h5py
  import os
  import netCDF4

  # os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
  #
  # f1 = h5py.File(file_name, 'r+')        # open the file
  # data = f1[var_name]                    # load the data
  # data[...] = data[...] + new_values     # assign new values to data
  # f1.close()

  dset = netCDF4.Dataset(file_name, 'r+')
  dset[var_name][:] = new_values

  dset.close()

  return True
