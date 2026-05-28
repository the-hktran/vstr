import importlib.util
import pathlib
import sys
import tempfile
import types
import unittest

import h5py
import numpy as np


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
MOL_PATH = REPO_ROOT / "nmode" / "mol.py"


def load_mol_module():
    if "test_nmode_mol" in sys.modules:
        return sys.modules["test_nmode_mol"]

    vstr_pkg = types.ModuleType("vstr")
    vstr_pkg.__path__ = []
    sys.modules.setdefault("vstr", vstr_pkg)

    utils_pkg = types.ModuleType("vstr.utils")
    utils_pkg.__path__ = []
    init_funcs = types.ModuleType("vstr.utils.init_funcs")
    constants = types.ModuleType("vstr.utils.constants")
    constants.AU_TO_INVCM = 1.0
    constants.ANGSTROM_TO_AU = 1.0
    utils_pkg.init_funcs = init_funcs
    utils_pkg.constants = constants
    sys.modules["vstr.utils"] = utils_pkg
    sys.modules["vstr.utils.init_funcs"] = init_funcs
    sys.modules["vstr.utils.constants"] = constants

    ff_pkg = types.ModuleType("vstr.ff")
    ff_pkg.__path__ = []
    force_field = types.ModuleType("vstr.ff.force_field")
    force_field.ScaleFC_me = lambda *args, **kwargs: None
    sys.modules["vstr.ff"] = ff_pkg
    sys.modules["vstr.ff.force_field"] = force_field

    cpp_pkg = types.ModuleType("vstr.cpp_wrappers")
    cpp_pkg.__path__ = []
    vhci_jf_pkg = types.ModuleType("vstr.cpp_wrappers.vhci_jf")
    vhci_jf_pkg.__path__ = []
    vhci_mod = types.ModuleType("vstr.cpp_wrappers.vhci_jf.vhci_jf_functions")
    vhci_mod.VCISparseHamNMode = lambda *args, **kwargs: np.zeros((1, 1))
    sys.modules["vstr.cpp_wrappers"] = cpp_pkg
    sys.modules["vstr.cpp_wrappers.vhci_jf"] = vhci_jf_pkg
    sys.modules["vstr.cpp_wrappers.vhci_jf.vhci_jf_functions"] = vhci_mod

    spectra_pkg = types.ModuleType("vstr.spectra")
    spectra_pkg.__path__ = []
    dipole_mod = types.ModuleType("vstr.spectra.dipole")
    dipole_mod.GetDipole = lambda *args, **kwargs: 0.0
    sys.modules["vstr.spectra"] = spectra_pkg
    sys.modules["vstr.spectra.dipole"] = dipole_mod

    perf_utils = types.ModuleType("vstr.utils.perf_utils")

    class DummyTimer:
        def __init__(self, *args, **kwargs):
            return None

        def start(self, *args, **kwargs):
            return None

        def stop(self, *args, **kwargs):
            return None

    perf_utils.TIMER = DummyTimer
    sys.modules["vstr.utils.perf_utils"] = perf_utils

    sys.modules.setdefault("scipy", types.ModuleType("scipy"))
    sys.modules.setdefault("numdifftools", types.ModuleType("numdifftools"))
    sys.modules.setdefault("xfacpy", types.ModuleType("xfacpy"))

    pyscf_pkg = types.ModuleType("pyscf")
    pyscf_pkg.__path__ = []
    sys.modules["pyscf"] = pyscf_pkg
    sys.modules["pyscf.gto"] = types.ModuleType("pyscf.gto")
    sys.modules["pyscf.scf"] = types.ModuleType("pyscf.scf")
    sys.modules["pyscf.cc"] = types.ModuleType("pyscf.cc")

    spec = importlib.util.spec_from_file_location("test_nmode_mol", MOL_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["test_nmode_mol"] = module
    spec.loader.exec_module(module)
    return module


def hermitian_tensor(order, ngridpts, seed):
    nstates = ngridpts ** order
    upper = np.triu_indices(nstates)
    matrix = np.zeros((nstates, nstates), dtype=float)
    values = np.arange(upper[0].size, dtype=float) + seed
    matrix[upper] = values
    matrix[(upper[1], upper[0])] = values
    return matrix.reshape((ngridpts,) * (2 * order))


class DummyNM:
    def __init__(self):
        self.x0 = np.zeros(3)
        self.V0 = 0.0
        self.mu0 = 0.0
        self.nm_coeff = np.eye(3)


class TestCompactNModeIntegralStorage(unittest.TestCase):
    def setUp(self):
        self.mod = load_mol_module()
        self.tmpdir = tempfile.TemporaryDirectory()
        self.ints_file = pathlib.Path(self.tmpdir.name) / "ints.h5"

    def tearDown(self):
        self.tmpdir.cleanup()

    def make_molecule(self):
        mol = self.mod.Molecule.__new__(self.mod.Molecule)
        mol.IntsFile = str(self.ints_file)
        mol.Order = 3
        mol.OrderPlus = None
        mol.Nm = 3
        mol.ngridpts = 2
        mol.use_onemode_states = True
        mol.doTCIResidual = False
        mol.nm = DummyNM()
        mol.Frequencies = np.array([1.0, 2.0, 3.0])
        mol.onemode_coeff = [np.eye(2) for _ in range(mol.Nm)]
        mol.onemode_eig = [np.array([0.0, 1.0]) for _ in range(mol.Nm)]
        mol.ints = [np.empty(0, dtype=object) for _ in range(5)]
        mol.ints[0] = np.stack([hermitian_tensor(1, mol.ngridpts, 10 + i) for i in range(mol.Nm)])
        mol.ints[1] = np.empty((mol.Nm, mol.Nm), dtype=object)
        mol.ints[2] = np.empty((mol.Nm, mol.Nm, mol.Nm), dtype=object)
        for i in range(mol.Nm):
            for j in range(i, mol.Nm):
                mol.ints[1][i, j] = hermitian_tensor(2, mol.ngridpts, 100 + 10 * i + j)
        for i in range(mol.Nm):
            for j in range(i, mol.Nm):
                for k in range(j, mol.Nm):
                    mol.ints[2][i, j, k] = hermitian_tensor(3, mol.ngridpts, 1000 + 100 * i + 10 * j + k)
        return mol

    def test_save_and_read_use_compact_canonical_storage(self):
        mol = self.make_molecule()
        mol.SaveIntegrals()

        with h5py.File(self.ints_file, "r") as handle:
            self.assertEqual(handle["ints"].attrs["storage"], "canonical_hermitian_packed")
            self.assertEqual(sorted(handle["ints/2"].keys()), ["1_1", "1_2", "1_3", "2_2", "2_3", "3_3"])
            self.assertEqual(sorted(handle["ints/3"].keys()), ["1_1_1", "1_1_2", "1_1_3", "1_2_2", "1_2_3", "1_3_3", "2_2_2", "2_2_3", "2_3_3", "3_3_3"])
            self.assertEqual(handle["ints/2/1_3"].ndim, 1)
            self.assertEqual(handle["ints/2/1_3"].shape[0], self.mod._nmode_packed_size(2, mol.ngridpts))

        read_mol = self.make_molecule()
        read_mol.ints = [np.asarray([]) for _ in range(5)]
        read_mol.ReadIntegralsAsArrays()

        two_mode_idx = self.mod._nmode_mode_index((0, 2))
        three_mode_idx = self.mod._nmode_mode_index((0, 1, 2))
        np.testing.assert_allclose(
            read_mol.ints[1][two_mode_idx],
            self.mod._pack_nmode_integral(mol.ints[1][0, 2]),
        )
        np.testing.assert_allclose(
            read_mol.ints[2][three_mode_idx],
            self.mod._pack_nmode_integral(mol.ints[2][0, 1, 2]),
        )

    def test_read_integrals_restores_permuted_dense_tensors(self):
        mol = self.make_molecule()
        mol.SaveIntegrals()

        read_mol = self.make_molecule()
        read_mol.ints = [np.asarray([]) for _ in range(5)]
        read_mol.ReadIntegrals()

        expected_two_mode = np.transpose(mol.ints[1][0, 2], (1, 0, 3, 2))
        expected_three_mode = np.transpose(mol.ints[2][0, 1, 2], (2, 0, 1, 5, 3, 4))
        np.testing.assert_allclose(read_mol.ints[1][2, 0], expected_two_mode)
        np.testing.assert_allclose(read_mol.ints[2][2, 0, 1], expected_three_mode)

    def test_get_ints_uses_canonical_pairs_for_two_mode_storage(self):
        nm = types.SimpleNamespace(
            nmodes=3,
            freqs=np.ones(3),
            mol=types.SimpleNamespace(doSaveIntsOTF=False),
        )
        nmode = self.mod.NModePotential(nm)

        call_count = {"value": 0}

        def potential_2mode(i, j, qi, qj):
            call_count["value"] += 1
            return np.full((1, 1), i + j + qi + qj)

        nmode.nm.potential_2mode = potential_2mode
        nmode.get_heg = lambda ngridpts, optimized=False, ngridpts0=None: ([np.array([0.0])] * 3, [np.array([[1.0]])] * 3)

        ints = nmode.get_ints(2, ngridpts=1, onemode_coeff=[np.array([[1.0]]) for _ in range(3)])

        self.assertEqual(call_count["value"], 6)
        np.testing.assert_allclose(ints[1, 0], ints[0, 1].transpose(1, 0, 3, 2))
        np.testing.assert_allclose(ints[2, 1], ints[1, 2].transpose(1, 0, 3, 2))


if __name__ == "__main__":
    unittest.main()
