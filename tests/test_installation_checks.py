"""Installation checks must not load models unless explicitly requested."""

from contextlib import redirect_stdout
import importlib.util
from io import StringIO
from pathlib import Path
import os
import sys
import unittest
from unittest.mock import patch


APP_DIR = Path(__file__).resolve().parents[1] / "gradio_app"
spec = importlib.util.spec_from_file_location("installation", APP_DIR / "test_installation.py")
installation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installation)


class InstallationTest(unittest.TestCase):
    def test_default_checks_do_not_load_a_model(self):
        original_cwd = Path.cwd()
        try:
            with (
                patch.object(sys, "argv", ["test_installation.py"]),
                patch.object(installation, "check_dependencies", return_value=True),
                patch.object(installation, "check_data_files", return_value=True),
                patch.object(installation, "test_sequence_access", return_value=True),
                patch.object(installation, "test_model_loading", return_value=True) as load,
                redirect_stdout(StringIO()),
            ):
                self.assertEqual(installation.main(), 0)
                load.assert_not_called()
                self.assertEqual(Path.cwd(), APP_DIR)
        finally:
            os.chdir(original_cwd)

    def test_model_loading_is_available_when_requested(self):
        original_cwd = Path.cwd()
        try:
            with (
                patch.object(sys, "argv", ["test_installation.py", "--load-model"]),
                patch.object(installation, "check_dependencies", return_value=True),
                patch.object(installation, "check_data_files", return_value=True),
                patch.object(installation, "test_sequence_access", return_value=True),
                patch.object(installation, "test_model_loading", return_value=False) as load,
                redirect_stdout(StringIO()),
            ):
                self.assertEqual(installation.main(), 1)
                load.assert_called_once_with()
        finally:
            os.chdir(original_cwd)


if __name__ == "__main__":
    unittest.main()
