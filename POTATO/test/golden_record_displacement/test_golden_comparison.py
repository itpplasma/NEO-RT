#!/usr/bin/env python3
"""Small numerical oracles for the displacement golden comparator."""
import unittest

import numpy as np

from test_potato_golden import RTOL, compare_rows


class GoldenComparisonTests(unittest.TestCase):
    def test_equal_opposite_corruption(self):
        reference = np.array([[0., 1.], [1., 2.]])
        actual = np.array([[0., 101.], [1., -98.]])
        np.testing.assert_array_equal(actual.sum(axis=0), reference.sum(axis=0))
        self.assertIsNotNone(compare_rows(actual, reference, [0]))

    def test_values_remain_associated_with_coordinates(self):
        reference = np.array([[0., 1.], [1., 2.]])
        actual = np.array([[0., 2.], [1., 1.]])
        # Even signed sums and absolute sums are identical after this swap.
        np.testing.assert_array_equal(actual.sum(axis=0), reference.sum(axis=0))
        np.testing.assert_array_equal(abs(actual).sum(axis=0), abs(reference).sum(axis=0))
        self.assertIsNotNone(compare_rows(actual, reference, [0]))

    def test_reordering_and_portable_rounding(self):
        reference = np.array([[0., 1.], [1., 2.]])
        actual = reference[::-1].copy()
        actual[:, 1] *= 1 + RTOL / 2
        self.assertIsNone(compare_rows(actual, reference, [0]))
        actual[:, 1] = reference[::-1, 1] * (1 + 2 * RTOL)
        self.assertIsNotNone(compare_rows(actual, reference, [0]))

    def test_small_modes_use_same_row_torque_scale(self):
        reference = np.array([[0., 10., 1e-12], [1., 1., 0.]])
        actual = reference.copy()
        actual[0, 2] += 1e-10
        self.assertIsNone(compare_rows(actual, reference, [0], [[1, 2]]))
        actual[1, 2] = 2 * RTOL
        self.assertIsNotNone(compare_rows(actual, reference, [0], [[1, 2]]))

    def test_shape_and_nonfinite_values_fail(self):
        reference = np.array([[0., 1.], [1., 2.]])
        self.assertIsNotNone(compare_rows(reference[:1], reference, [0]))
        for value in (np.nan, np.inf, -np.inf):
            actual = reference.copy()
            actual[0, 1] = value
            self.assertIsNotNone(compare_rows(actual, reference, [0]))


if __name__ == "__main__":
    unittest.main()
