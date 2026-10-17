"""Independent readback of actual native matrices; no NEO-RT solver imports."""
import h5py
import numpy as np


def check_receipt(baseline, diagnostic):
    with h5py.File(baseline) as old, h5py.File(diagnostic) as new:
        def unchanged(name, obj):
            if isinstance(obj, h5py.Dataset) and name != "created":
                np.testing.assert_array_equal(obj[()], new[name][()], err_msg=name)
        old.visititems(unchanged)
        group = new["mixture/conserving_matrix_receipt"]
        def matrix(name):
            return (group[name + "_real"][()] +
                    1j * group[name + "_imag"][()]).T
        mass = group["mass"][()].T
        collision = matrix("collision")
        streaming = group["streaming"][()].T
        precession = group["precession"][()].T
        source, hamiltonian = matrix("source"), matrix("hamiltonian")
        coefficient = matrix("coefficient").ravel()
        cross = matrix("source_coefficient")
        operator = collision + streaming + 1j * precession
        np.testing.assert_array_equal(group["dof_offsets"][()], [1, 18, 35, 52])
        pairs = []
        for k, indices in enumerate(((1, 1), (1, 2), (2, 2),
                                     (1, 3), (2, 3), (3, 3)), 1):
            pair = group[f"pair_{k}"]
            np.testing.assert_array_equal(pair["species_indices"][()], indices)
            block = (pair["collision_real"][()] +
                     1j * pair["collision_imag"][()]).T
            np.testing.assert_array_equal(block, block.conj().T)
            excluded = [j for j in range(51) if j // 17 + 1 not in indices]
            np.testing.assert_array_equal(block[excluded, :], 0)
            np.testing.assert_array_equal(block[:, excluded], 0)
            assert np.linalg.eigvalsh(block).min() >= -2e-13 * np.linalg.norm(block)
            pairs.append(block)
        scale = np.linalg.norm(collision)
        assert np.linalg.norm(sum(pairs) - collision) < 2e-14 * scale
        assert np.linalg.norm(operator @ coefficient - source.sum(axis=1)) < (
            2e-12 * max(np.linalg.norm(source.sum(axis=1)),
                        np.linalg.norm(operator) * np.linalg.norm(coefficient)))
        np.testing.assert_allclose(cross.sum(axis=1), coefficient, rtol=2e-10,
                                   atol=2e-13 * np.max(np.abs(coefficient)))
        values, vectors = np.linalg.eigh(mass)
        inverse_sqrt = (vectors / np.sqrt(values)) @ vectors.T
        singular = np.linalg.svd(inverse_sqrt @ operator @ inverse_sqrt,
                                 compute_uv=False)
        np.testing.assert_allclose(group["singular_values"][()], singular,
                                   rtol=2e-10)
        for species, name in enumerate(("D", "T", "e")):
            data = new[f"species/{name}/transport"]
            columns = hamiltonian[:, 3*species:3*species+3]
            work = (data["conserving_work_real"][()] +
                    1j * data["conserving_work_imag"][()]).ravel()
            np.testing.assert_allclose(columns.conj().T @ coefficient, work,
                                       rtol=2e-12)
            source_work = (data["conserving_cross_work_real"][()] +
                           1j * data["conserving_cross_work_imag"][()]).T
            np.testing.assert_allclose(columns.conj().T @ cross, source_work,
                                       rtol=2e-12)
