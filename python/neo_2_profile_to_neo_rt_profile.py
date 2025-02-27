#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

sys.path.append("/var/tmp/ert/code/NEO-2/python/src/neo2_util/")

def neo_2_profile_to_neo_rt_profile(input_neo2_profile_name: str, input_neo2_outputfile_name: str,
    s_min: float, s_max: float, number_s_points: int):
  """
  input:
  ------
  input_neo2_profile_name: string name of the profile file of neo-2 to
    be used as input (and which was input to neo-2).
  input_neo2_outputfile_name: string name of the output file of neo-2,
    collected over the flux surfaces. Used here as input to get toroidal
    mach number.
  s_min, s_max: floating point numbers, minimum and maximum values for
    output s-grid to use. Note that these values will be included in the
    grid.
  number_s_points: integer, number of flux surfaces the output should
    have.

  output:
  -------
  None

  sideeffects:
  ------------
  Files profile.in and plasma.in are written. If they already exist,
  they are overwritten.
  """
  import numpy as np
  from scipy.interpolate import InterpolatedUnivariateSpline as Spline

  import hdf5tools

  outputfilename_plasma = 'plasma.in'
  outputfilename_profile = 'profile.in'

  dalton_in_kg = 1.66053906660e-27
  conversion_gram_to_kg = 1.0e-3
  conversion_erg_to_ev = 6.2415e11
  conversion_erg_to_joule = 1.0e-7
  conversion_m_per_s_to_cm_per_s = 1.0e2

  index_electron = 0
  index_ion1 = 1
  index_ion2 = 2

  # read in data from neo-2 profile file
  with hdf5tools.get_hdf5file(input_neo2_profile_name) as input_profile:
    number_input_file_points = input_profile['num_radial_pts'][0]
    s_input = input_profile['boozer_s'][:]

    charge_ion1 = input_profile['species_def'][0, index_ion1, 0]
    # ~ density_electron = input_profile['n_prof'][0, :]
    density_input_ion1 = input_profile['n_prof'][index_ion1, :]
    mass_ion1 = input_profile['species_def'][1, index_ion1, 0]
    temperature_input_electron = input_profile['T_prof'][index_electron, :]
    temperature_input_ion1 = input_profile['T_prof'][index_ion1, :]
    velocity_input = input_profile['Vphi'][:]
    if (input_profile['num_species'][0] > 2):
      charge_ion1 = input_profile['species_def'][0, index_ion2, 0]
      density_input_ion2 = input_profile['n_prof'][index_ion2, :]
      mass_ion1 = input_profile['species_def'][1, index_ion2, 0]
      temperature_input_ion2 = input_profile['T_prof'][index_ion2, :]
    else:
      charge_ion2 = charge_ion1
      density_input_ion2 = 0.0*density_input_ion1 # Scaling to have same size/shape, but zero value
      mass_ion2 = mass_ion1
      temperature_input_ion2 = temperature_input_ion1 / temperature_input_ion1 # Scaling to have same size/shape, but value of one.

  with hdf5tools.get_hdf5file(input_neo2_outputfile_name) as input_neo2_output:
    mtoverr = input_neo2_output['MtOvR'][:, index_ion1]
    R0 = input_neo2_output['R0'][:]
    s_neo2_output = input_neo2_output['boozer_s'][:]

  # checks
  if (s_min < s_input[0]):
    print('WARNING: s_min ({}) smaller than lower bound in input profile ({})!'.format(s_min, s_input[0]))
  if (s_max > s_input[-1]):
    print('WARNING: s_max ({}) larger than upper bound in input profile ({})!'.format(s_max, s_input[-1]))
  if (number_s_points > number_input_file_points):
    print('WARNING: requested number of points ([}) larger than number of point in input ({})!'.format(number_s_points, number_input_file_points))

  # calculate derived quantities and make conversions
  thermal_velocity = np.sqrt(2*temperature_input_ion1*conversion_erg_to_joule/(mass_ion1*conversion_gram_to_kg))*conversion_m_per_s_to_cm_per_s
  # ~ mach_number_input = velocity_input / thermal_velocity
  mass_ion1 = mass_ion1*conversion_gram_to_kg / dalton_in_kg
  mass_ion2 = mass_ion2*conversion_gram_to_kg / dalton_in_kg

  temperature_input_ion1 = conversion_erg_to_ev*temperature_input_ion1
  temperature_input_ion2 = conversion_erg_to_ev*temperature_input_ion2
  temperature_input_electron = conversion_erg_to_ev*temperature_input_electron

  mach_number_input = mtoverr*R0

  # interpolation
  s = np.linspace(s_min, s_max, number_s_points)
  # ~ mach_number = Spline(s_input, mach_number_input)(s)
  mach_number = Spline(s_neo2_output, mach_number_input)(s)
  density_ion1 = Spline(s_input, density_input_ion1)(s)
  density_ion2 = Spline(s_input, density_input_ion2)(s)
  temperature_ion1 = Spline(s_input, temperature_input_ion1)(s)
  temperature_ion2 = Spline(s_input, temperature_input_ion2)(s)
  temperature_electron = Spline(s_input, temperature_input_electron)(s)

  # write files
  with open('profile.in', 'w') as profile, open('plasma.in', 'w') as plasma:
    # plasma.in gets a header:
    plasma.write('\% N am1 am2 Z1 Z2\n')
    plasma.write('{} {} {} {} {}\n'.format(number_s_points, mass_ion1, mass_ion2, charge_ion1, charge_ion2))
    plasma.write('\% s ni_1[cm^-3] ni_2[cm^-3] Ti_1[eV] Ti_2[eV] Te[eV]\n')

    for k in range(number_s_points):
      plasma.write('{: 17.10e} {: 17.10e} {: 17.10e} {: 17.10e} {: 17.10e} {: 17.10e}\n'.format(
          s[k], density_ion1[k], density_ion2[k], temperature_ion1[k], temperature_ion2[k], temperature_electron[k]))
      profile.write('{: 17.10e} {: 17.10e} {: 17.10e}\n'.format(
          s[k], mach_number[k], thermal_velocity[k]))


if __name__ == "__main__":
  neo_2_profile_to_neo_rt_profile(
    '/temp/grassl_g/GPEC_NEO2_AUG_30835/RUN_stor_lag6/aug_30835.in',
    '/temp/grassl_g/GPEC_NEO2_AUG_30835/RUN_stor_lag6/aug_30835.out',
      0.05, 0.9, 100)
