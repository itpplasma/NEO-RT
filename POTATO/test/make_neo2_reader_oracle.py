import sys

import h5py
import numpy as np


with h5py.File(sys.argv[1], "w", libver="latest") as handle:
    vectors = {
        "boozer_s": [0.0, 0.5, 1.0],
        "aiota": [0.5, 0.5, 0.5],
        "Er": [-2.0, 0.0, 2.0],
    }
    for name, values in vectors.items():
        dataset = handle.create_dataset(name, data=np.array(values, dtype=np.float64))
        dataset.attrs["lbounds"] = np.array([1], dtype=np.int32)
        dataset.attrs["ubounds"] = np.array([3], dtype=np.int32)

    matrices = {
        "n_spec": np.arange(1.0, 7.0, dtype=np.float64),
        "T_spec": np.arange(11.0, 17.0, dtype=np.float64),
        "MtOvR": np.arange(21.0, 27.0, dtype=np.float64),
    }
    for name, values in matrices.items():
        dataset = handle.create_dataset(name, data=values.reshape(2, 3))
        dataset.attrs["lbounds"] = np.array([1, 1], dtype=np.int32)
        dataset.attrs["ubounds"] = np.array([3, 2], dtype=np.int32)
