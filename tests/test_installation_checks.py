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
    def run_main(self, argv, model_loads):
        original_cwd = Path.cwd()
        try:
            with (
                patch.object(sys, "argv", ["test_installation.py", *argv]),
                patch.object(installation, "check_dependencies", return_value=True),
                patch.object(installation, "check_data_files", return_value=True),
                patch.object(installation, "test_sequence_access", return_value=True),
                patch.object(installation, "test_model_loading", return_value=model_loads) as load,
                redirect_stdout(StringIO()),
            ):
                return installation.main(), load.call_count, Path.cwd()
        finally:
            os.chdir(original_cwd)

    def test_default_checks_do_not_load_a_model(self):
        self.assertEqual(self.run_main([], model_loads=True), (0, 0, APP_DIR))

    def test_model_loading_is_available_when_requested(self):
        rc, calls, _ = self.run_main(["--load-model"], model_loads=False)
        self.assertEqual((rc, calls), (1, 1))


if __name__ == "__main__":
    unittest.main()
