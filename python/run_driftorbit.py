#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def get_profile_data_from_file_data(profile_data, flux_surface_number: int):
  """
  input:
  ------
  profile_data: 2d array-like structure, i.e. one that can be indexed
    with two indices. First index is for flux surface, second determines
    the quantity.
  flux_surface_number: flux surface index requested.

  output:
  -------
  List with the requested values (s, M_t, vth, epsm).

  sideeffects:
  ------------
  None
  """
  s       = profile_data[flux_surface_number,0]
  M_t     = profile_data[flux_surface_number,1]
  vth     = profile_data[flux_surface_number,2]
  #epsm    = profile_data[flux_surface_number,3]
  epsm = 1.0e+0

  return [s, M_t, vth, epsm]


def get_profile_data_for_flux_surface(profile_file_name: str,
    flux_surface_number: int):
  """Get profile data for flux surface with given number.

  Return profile data for a given flux surface number from a profile
  file.

  input:
  ------
  profile_file_name: string, contains the complete name (including file
    extension) from which to get the data.
  flux_surface_number: integer.

  output:
  -------
  List with two elements. First is list with actual data requested,
  second entry is the complete data of the file (as numpy array?). This
  is done so it can be reused (via other function) if data for more than
  one flux surface is required.

  sideeffects:
  ------------
  None
  """

  import numpy as np

  profile_data = np.loadtxt(profile_file_name)

  return [get_profile_data_from_file_data(profile_data, flux_surface_number), profile_data]


def run_single_flux_surface(executable_name: str, template_file_name: str, runname: str,
    s: float, M_t: float, vth: float, epsm: float):
  """Run code for given template parameter values.

  This function will run the executable for a given set of template
  parameters. This is done by replacing the templates with the actual
  values and then use this as input.

  input:
  ------
  executable_name: string, name of the executable to run. Is assumed
    to be found in the current directory.
  template_file_name: string, name of the template file.
  runname: string, name (without ending) to use for the input file of
    the executable, and which is passed to the executable (and thus also
    used as prefix for the output files).
  s, M_t, vth, epsm: floating point numbers, the values that should be
    used for the respective template parameters.

  output:
  -------
  none

  sideeffects:
  ------------
  Creates an input file for the executable.
  Uses subprocess to run executable, and thus files are produced,
  including those for standard and error output.
  """

  from subprocess import run
  import sys
  import re

  dic = {
      '<S_TOKEN>': s,
      '<M_T_TOKEN>': M_t,
      '<VTH_TOKEN>': vth,
      '<EPSM_TOKEN>': epsm
      }

  pattern = re.compile('|'.join(dic.keys()))

  with open(template_file_name,'r') as tempf, open(runname+'.in','w') as outf:
    for line in tempf:
      result = pattern.sub(lambda x: str(dic[x.group()]), line)
      outf.write(result)

  with open(runname+'.log', 'w') as log, open(runname+'.err', 'w') as err:
    retcode = run('./{} {}'.format(executable_name, runname),
        shell=True, stdout=log, stderr=err)


def run_multiple_flux_surfaces(executable_name: str,
    profile_file_name: str, template_file_name: str, base_output_file_name: str,
    index_low: int, index_high: int):
  """Make simulations for range of flux surfaces.

  Run simulations for a range of flux surfaces from a given profile.
  For this a template input file is changed according to the
  corresponding parameters from the profile (file).

  input:
  ------
  executable_name: string, name of the executable to run. Is assumed
    to be found in the current directory.
  profile_file_name: string, contains the complete name (including file
    extension) from which to get the data.
  template_file_name: string, name of the template file.
  base_output_file_name: string, name (without ending) to use for the input file of
    the executable, and which is passed to the executable (and thus also
    used as prefix for the output files).
    Eeach flux surface gets its own number appended to this.
  index_low, index_high: integers, range of indices to consider. As
    usual for python lower end is included, higher end excluded.
    Note: actual valid indices are required, i.e. -1 is not allowed.

  output:
  -------
  None

  sideeffects:
  ------------
  There is one call to run_single_flux_surface per flux surfaces, so
  same as for those function.
  """
  # The profile data is assumed to be constant, and thus only needed
  # once. We ignore the other return values from this call.
  [[s, M_t, vth, epsm], profile_data] = get_profile_data_for_flux_surface(profile_file_name, index_low)

  for k in range(index_low, index_high):
    print('Processing: index {} of {}'.format(k+1, index_high))
    [s, M_t, vth, epsm] = get_profile_data_from_file_data(profile_data, k)

    runname = '{}{}'.format(base_output_file_name, k)
    run_single_flux_surface(executable_name, template_file_name, runname, s, M_t, vth, epsm)


if __name__ == "__main__":
  import argparse
  import os

  parser = argparse.ArgumentParser(
      description="Populate drift-orbit inputs and run neo_rt.x over one or more flux surfaces.")
  parser.add_argument(
      "lower_index",
      type=int,
      help="First flux-surface index to process (zero-based)")
  parser.add_argument(
      "upper_index",
      type=int,
      nargs="?",
      help="Optional upper bound (exclusive). If omitted only lower_index is run.")
  parser.add_argument(
      "--exe",
      default=os.environ.get("NEORT_EXECUTABLE", "neo_rt.x"),
      help="Path to the neo_rt executable (default: %(default)s or $NEORT_EXECUTABLE)")
  parser.add_argument(
      "--template",
      default="driftorbit.in.template",
      help="Template input file to populate (default: %(default)s)")
  parser.add_argument(
      "--profile",
      default="profile.in",
      help="Profile file providing s, M_t, and vth columns (default: %(default)s)")
  parser.add_argument(
      "--prefix",
      default="driftorbit",
      help="Prefix for generated input/output files (default: %(default)s)")

  args = parser.parse_args()

  template_file_name = args.template
  profile_file_name = args.profile
  executable_name = args.exe
  base_output_file_name = args.prefix

  if args.upper_index is None:
    fsnum = args.lower_index
    runname = f"{base_output_file_name}{fsnum}"
    [[s, M_t, vth, epsm], profile_data] = get_profile_data_for_flux_surface(profile_file_name, fsnum)
    run_single_flux_surface(executable_name, template_file_name, runname, s, M_t, vth, epsm)
  else:
    index_low = args.lower_index
    index_high = args.upper_index
    run_multiple_flux_surfaces(
        executable_name,
        profile_file_name,
        template_file_name,
        base_output_file_name,
        index_low,
        index_high)
