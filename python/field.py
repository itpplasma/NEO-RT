# Reads files of the type
#
# CC Boozer-coordinate data file Version 0.1 Author J.Geiger Created: 26.07.2010
# CC based on calculation: tok_examples/boozer.tok02
# CC <beta> = 0.006432
# CC Configurationname = tok02, aspect ratio ca 3.8, q_0 > 1, q_a < 3
#  m0b  n0b nsurf nper flux/[Tm^2]     a/[m]     R/[m]
#  18    0    63   1   -1.330000e+00   0.46000   1.64377
#        s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)
#                             [A]            [A]   dp/ds,[Pa] (dV/ds)/nper
#   2.3438E-02  8.9997E-01 -1.7754E+07 -2.3937E+04 -1.7231E+04 -7.7240E+00
#     m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     bmn/[T]
#      0    0  1.81092169e+00  0.00000000e+00  0.00000000e+00  1.96335464e+00
#      1    0  6.96313408e-02  7.05114198e-02  5.75628561e-06 -7.54078850e-02
#      2    0  1.81691876e-03  1.83233922e-03  1.41838683e-06 -5.24507105e-04
#      3    0  5.26925565e-05  5.19496929e-05  7.22883916e-08 -8.14119880e-06
#      4    0  4.88665890e-06  4.87741134e-06  1.50368389e-09 -3.73627634e-06
#      5    0  5.16747065e-07  5.21371270e-07  1.37963561e-10 -3.73512433e-07
#      6    0  2.81138594e-08  2.83759676e-08  2.60830060e-11 -1.33593560e-08
#      7    0  1.13594171e-09  1.14171750e-09  2.18323901e-12  1.55087102e-08
#      8    0  9.25116789e-11  9.37984697e-11  1.30458799e-13  2.45353472e-09
#      9    0  7.90964516e-12  8.10274666e-12  1.13545431e-14  2.23212394e-10
#     10    0  2.62928568e-13  2.78831157e-13  1.16740303e-15  2.22889855e-11
#     11    0 -9.04137876e-15 -6.13203065e-15  8.74728807e-17  2.77964898e-12
#     12    0 -2.07819872e-15 -5.77446077e-16  5.23684075e-18  2.89611276e-13
#     13    0  2.79290480e-16 -3.37728977e-17  3.77584789e-19  2.15289736e-14
#     14    0 -2.78726695e-15  5.87637577e-17  3.47084985e-20 -8.15397559e-16
#     15    0 -4.39752401e-15 -2.75929453e-17  1.07202607e-21 -1.34522897e-15
#     16    0  6.97358837e-16 -7.95804395e-17 -6.45862622e-21 -3.03974546e-15
#     17    0 -2.21350716e-15  1.73689188e-16  1.30694537e-20 -1.57651395e-15
#     18    0 -3.35495520e-15  5.49690501e-17  5.02264068e-21 -5.63108126e-15
#       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)
#                            [A]            [A]   dp/ds,[Pa] (dV/ds)/nper
#  3.9063E-02  8.9969E-01 -1.7740E+07 -4.0657E+04 -1.8846E+04 -7.7199E+00
#    m    n        r/[m]           z/[m] (phib-phi)*nper/twopi     bmn/[T]
#  ...

import numpy as np

class Field:
    """Base class holding electromagnetic field data"""

    def init(self):
        # Initialize empty arrays
        self.A = np.zeros(3)  # Covariant     components of A
        self.B = np.zeros(3)  # Contravariant components of B
        self.h = np.zeros(3)  # Covariant components of h = B/Bmod
        self.Bmod = 0.0
        self.sqrtg = 0.0

    def evaluate(self, x):
        """Evaluates field at position x and stores result in class members"""
        pass


# TODO: implement based on NEO-RT
class BoozerField(Field):
    def __init__(self, filename='in_file'):
        super(BoozerField, self).__init__()
        # TODO: 
        # 1) Read Boozer file of type "inp_swi=8" see do_magfie_standalone.f90, boozer_read
        # 2) Store constant quantities like psi_pr based on flux as class members
        # 3) Construct splines for 1D quantities e.g. with scipy.interpolate.InterpolatedUnivariateSpline
        #       s         iota  curr_pol/nper    curr_tor    pprime   sqrt g(0,0)
        self.psi_pr = 0.0  # TODO

    def evaluate(self, x):
        """Evaluates field at position x = (s, \vartheta, \varphi) in Boozer coordinates
        and stores result in class members"""
        
        # In any flux coordinates, A_s = 0, A_th = A_th(s), A_ph = A_ph(s) are
        # independent of flux angles th and ph.

        self.A[0] = 0.0                # A_s = 0, since flux coordinates
        self.A[1] = -self.psi_pr*x[0]  # A_\vartheta = -psi_tor, see d'haeseleer
        # TODO: implement A_\varphi via spline antiderivative to integrate and
        # \psi_{pol} = \int d\psi_{pol}/ds ds 
        #            = \int d\psi_{tor}/ds * d\psi_{pol}/d\psi_{tor} ds.
        #            = \int psi_pr*iota ds
        # \psi_{pol}(s=0) should vanish.
        # Acov[2] = ...

        # TODO: implement the rest with call to libneo_rt_mc.magfie for now
        
