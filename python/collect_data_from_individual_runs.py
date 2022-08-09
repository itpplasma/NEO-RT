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
  files2 = [f.replace('torque','magfie_param') for f in files]

  return [files, files2]


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

  for f in files:
    kf = kf+1
    dat = np.loadtxt(f)
    if dat.size == 0:
      continue

    data.append(dat)

    with open(files2[kf]) as f2:
      for kl in range(12):
        f2.readline()
      dat2 = f2.readline()
      q.append(float(dat2.split()[-1]))

  data = np.array(data)
  q = np.array(q)
  # Sort datapoints, file listing uses different order (e.g. 10 before
  # 2), thus probably required.
  order = np.argsort(data[:,0])
  data = data[order,:]
  q = q[order]

  return [data, q]


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

  [files, files2] = get_file_lists(pattern = infilepattern)

  data = []
  q = []
  kf = -1

  smin = 1.0e-4
  smax = 1.0e+0

  [data, q] = load_files_from_lists(files, files2)

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

  with hdf5tools.get_hdf5file_new('neo-rt_out.h5') as h5f:
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
