"""Check FORD's actual source selection when documentation tooling is installed."""

from contextlib import chdir
from importlib.util import find_spec
from pathlib import Path
import unittest


@unittest.skipUnless(find_spec("ford"), "documentation tooling is not installed")
class DocsSourceSelectionTest(unittest.TestCase):
    def test_sources_are_owned_and_use_one_field_backend(self):
        from ford import load_settings
        from ford.fortran_project import find_all_files

        root = Path(__file__).resolve().parents[1]
        _, settings = load_settings((root / "doc.md").read_text(), root)
        settings.normalise_paths(root)
        with chdir(root):
            files = find_all_files(settings)

        self.assertTrue(files)
        self.assertTrue(all(path.is_relative_to(root / "src") for path in files))
        self.assertIn(root / "src/diag/neort_diag.f90", files)
        self.assertIn(root / "src/do_magfie_standalone.f90", files)
        self.assertNotIn(root / "src/do_magfie_neo.f90", files)


if __name__ == "__main__":
    unittest.main()
