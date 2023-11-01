#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_file_lists(pattern: str):
  """
  Creates two file lists from a given pattern.

  Creates list of file names from the pattern. A second list replaces
  'torque' in filenames of the first list by 'magfie_param'.

  Example:
    [files, files2] = get_file_lists(pattern = r'driftorbit([0-9]+)_torque\.out')

  input:
  ------
  pattern: (raw) string, the file pattern used to create the file lists.
    Expected to contain 'torque'.

  output:
  -------
  files: list of file names matching the pattern.
  files2: list same as files, but with 'torque' replaced by
    'magfie_param'. This means the lists should have the same number of
    elements.

  sideeffects:
  none
  """

  import os
  import re

  pattern = re.compile(pattern)

  files = os.listdir()
  files = [f for f in files if pattern.match(f)]
  files2 = [f.replace('torque', 'magfie_param') for f in files]
  files3 = [f.replace('_torque', '') for f in files]

  return [files, files2, files3]


def read_magfie_param(infilename: str, data_dictionary: dict):
  """
  Read a magfie_param file and append data to dictionary.

  Assumptions on magfie_param file:
  Might combine lines with and without data.
  Lines without data should have at most 1 whitespace.
  Data lines must have at least three fields seperated by whitespaces.
  The first field is ignored, the second field is assumed to be the name
  of the entry for the dictionary, and the last field is assumed to be
  the value.
  Values are converted to floating point values if possible.

  Examples of (data) lines (without the whitespace at the left):
   test_magfie: T [eV]    =    2129.6338999999939
   test_magfie: m0        =    0.0000000000000000
   -------------------------
   test_magfie: pertfile  =  T
  Here the third line whould be ignored, as is the first field
  ("test_magfie:"). 'T', 'm0' and 'pertfile' would be used as names for
  the fields, while ~2129.63, 0 and 'T' would be appended to the values.

  input:
  ------
  infilename: string, name (and path) of the file to read.
  data_dictionary: dictionary, where the value is a list, to which the
    data is appended.

  output:
  -------
  the modified data_dictionary.

  sideeffects:
  ------------
  none
  """
  with open(infilename) as f:
    for line in f:
      if len(line.split()) >= 3:
        parts = line.split()
        try:
          d = float(parts[-1])
        except:
          d = parts[-1]
        if parts[1] in data_dictionary.keys():
          data_dictionary[parts[1]].append(d)
        else:
          data_dictionary[parts[1]] = [d]

  return data_dictionary


def load_files_from_lists(files: list, files2: list):
  """
  Load data from two lists of files, with assumptions about content.

  Given two list of files, load data from the text files.
  The files in the first list are assumed to have numerical data in
  matrix form, i.e. every line has same number of entries/numbers, not
  necessarily characters. The first column is used to sort the data from
  booth files. As file listing probably uses lexicographic ordering this
  is probably necessary to have entries in suitable order.
  For the files in the second list it is assumed that there
  are at least thirteen lines in the file, with only the last value in
  the thirteenth line used as floating point number.

  input:
  ------
  files: list, first list of file names from which to read data.
  files2: list, second list of file names from which to read data.
    Note, that it is assumed that non-empty file number n in list files
    corresponds to entry number n in list files2 and that list files2
    has at least as much entries as list files.

  output:
  -------
  data: list of lists, the data read from list files.
  q: list, data read from list files2.
  """
  import numpy as np

  data = []
  q = []
  kf = -1
  data_magfie_param = {}

  for f in files:
    kf = kf+1
    dat = np.loadtxt(f)
    if dat.size == 0:
      continue

    data.append(dat)

    with open(files2[kf]) as f2:
      read_magfie_param(files2[kf], data_magfie_param)

  data = np.array(data)
  q = np.array(data_magfie_param["q"])
  # Sort datapoints, file listing uses different order (e.g. 10 before
  # 2), thus probably required.
  order = np.argsort(data[:,0])
  data = data[order,:]
  q = q[order]
  for k in data_magfie_param.keys():
    data_magfie_param[k] = np.array(data_magfie_param[k])[order]

  return [data, q]


def load_transport_coefficient_files_from_list(files: list):
  """Load data froma files, with minor assumptions about content.

  Given a list of filenames, data from these files is loaded.
  Only assumption about the content, is that first column should be used
  for sorting.

  input:
  ------
  """
  import numpy as np

  data = []

  any_found = False

  for f in files:
    try:
      dat = np.loadtxt(f)
      any_found = True
    except OSError as e:
      continue

    if dat.size == 0:
      continue

    data.append(dat)

  data = np.array(data)

  if any_found:
    # Sort datapoints, file listing uses different order (e.g. 10 before
    # 2), thus probably required.
    order = np.argsort(data[:,0])
    data = data[order,:]

  return data


def collect_torque_data(infilepattern: str, outfilename: str):
  """
  Collect torque data from radial scan of neo-rt in hdf5 file.

  Loads data from all torque and magfie_param files matching the
  infilepattern.
  Some additional quantities are calculated.
  The quantities are then written to an hdf5 file.

  Example:
    collect_torque_data(infilepattern=r'driftorbit([0-9]+)_torque\.out',
                        outfilename='neo-rt_out.h5')

  input:
  ------
  infilepattern: (raw) string, (used to create) the pattern which
    determines the input files.
    Note that 'torque' is expected to be a part of the name, and that
    also files with 'magfie_param' instead, exist.
  outfilename: string, name(+path) of file to which data is written.

  output:
  -------
  none

  sideeffects:
  ------------
  Creates new file.
  """
  import numpy as np
  import scipy.integrate as spint
  import scipy.interpolate as spi
  import scipy.optimize as spo

  import hdf5tools

  [files, files2, files3] = get_file_lists(pattern = infilepattern)

  data = []
  q = []
  kf = -1

  smin = 1.0e-4
  smax = 1.0e+0

  [data, q] = load_files_from_lists(files, files2)
  data_transport = load_transport_coefficient_files_from_list(files3)

  iotaspl = spi.splrep(data[:,0],1.0/q)
  iotaint = spint.quad(lambda x: spi.splev(x, iotaspl), 0.0, 1.0)[0]

  condi = (data[:,0]<smax)*(data[:,0]>smin)
  data = data[condi]
  q = q[condi]

  s = data[:,0]
  rhopol = np.sqrt([spi.splint(0.0, sk, iotaspl)/iotaint for sk in s])

  dVds = data[:,1]*1e-6
  sbox = np.append(smin,np.append(s[:-1]+np.diff(s)/2.0,smax))
  ds = np.diff(sbox)
  Mt = data[:,2]
  dTphi_int_ds = np.sum(data[:,3:6],1)*1e-7
  Tphi = dTphi_int_ds/dVds

  spli = spi.InterpolatedUnivariateSpline(rhopol**2, np.cumsum(dTphi_int_ds*ds))
  dTphi_int_dspol = spli.derivative()

  qeval = lambda x: 1.0/spi.splev(x, iotaspl)

  with hdf5tools.get_hdf5file_new(outfilename) as h5f:
    # Torque data (and q)
    h5f.create_dataset('s', data=s)
    h5f['s'].attrs['unit'] = '1'
    h5f.create_dataset('spol', data=rhopol**2)
    h5f['spol'].attrs['unit'] = '1'
    h5f.create_dataset('q', data=qeval(s))
    h5f['q'].attrs['unit'] = '1'
    h5f['q'].attrs['comment'] = 'safety factor'
    h5f.create_dataset('dVds', data=dVds)
    h5f['dVds'].attrs['unit'] = 'm^3/1'
    h5f.create_dataset('Mt', data=Mt)
    h5f['Mt'].attrs['unit'] = '1'
    h5f['Mt'].attrs['comment'] = 'toroidal mach number'
    h5f.create_dataset('dTphi_int_ds', data=dTphi_int_ds)
    h5f['dTphi_int_ds'].attrs['unit'] = 'Nm/?'
    h5f['dTphi_int_ds'].attrs['comment'] = 'derivative of toroidal torque density'
    h5f.create_dataset('dTphi_int_dspol', data=dTphi_int_dspol(rhopol**2))
    h5f['dTphi_int_dspol'].attrs['unit'] = 'Nm/?'
    h5f['dTphi_int_dspol'].attrs['comment'] = 'derivative of toroidal torque density'
    h5f.create_dataset('Tphi', data=Tphi)
    h5f['Tphi'].attrs['unit'] = 'Nm/m^3'
    h5f['Tphi'].attrs['comment'] = 'toroidal torque density'

    if data_transport.shape != (0,):
      # Transport coefficients, minus one due to zero based index.
      h5f.create_dataset('D11_copassing', data=data_transport[:,2-1])
      h5f['D11_copassing'].attrs['unit'] = '?'
      h5f.create_dataset('D11_counterpassing', data=data_transport[:,3-1])
      h5f['D11_counterpassing'].attrs['unit'] = '?'
      h5f.create_dataset('D11_trapped', data=data_transport[:,4-1])
      h5f['D11_trapped'].attrs['unit'] = '?'
      h5f.create_dataset('D11', data=data_transport[:,5-1])
      h5f['D11'].attrs['unit'] = '?'
      h5f['D11'].attrs['comment'] = 'sum of co+counter+trapped'
      h5f.create_dataset('D12_copassing', data=data_transport[:,6-1])
      h5f['D12_copassing'].attrs['unit'] = '?'
      h5f.create_dataset('D12_counterpassing', data=data_transport[:,7-1])
      h5f['D12_counterpassing'].attrs['unit'] = '?'
      h5f.create_dataset('D12_trapped', data=data_transport[:,8-1])
      h5f['D12_trapped'].attrs['unit'] = '?'
      h5f.create_dataset('D12', data=data_transport[:,9-1])
      h5f['D12'].attrs['unit'] = '?'
      h5f['D12'].attrs['comment'] = 'sum of co+counter+trapped'


if __name__ == "__main__":
  """
  Script to collect torque data from radial scan of neo-rt in hdf5 file.

  Loads data from all torque and magfie_param files matching the
  pattern.
  Some additional quantities are calculated.
  The quantities are then written to an hdf5 file.
  """

  collect_torque_data(infilepattern=r'driftorbit([0-9]+)_torque\.out',
                        outfilename='neo-rt_out.h5')
