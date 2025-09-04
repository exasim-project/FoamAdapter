

import subprocess
import sys
import shutil
import pytest
import os

# Fixture to run Allclean before and Allrun after each test
@pytest.fixture(autouse=True)
def run_allrun_allclean(request):
    test_dir = str(request.fspath.dirname)
    allclean = os.path.join(test_dir, "Allclean")
    allrun = os.path.join(test_dir, "Allrun")
    if os.path.isfile(allrun):
        subprocess.run([allrun], cwd=test_dir, check=False)
    yield
    if os.path.isfile(allclean):
        subprocess.run([allclean], cwd=test_dir, check=False)

def test_icofoam_serial():
    # Run the solver in serial mode
    result = subprocess.run([
        sys.executable, '-m', 'foamadapter.cli.app', 'solver', 'icofoam'
    ], capture_output=False, text=True)
    assert result.returncode == 0
@pytest.mark.skipif(not shutil.which("mpirun"), reason="mpirun not available")
def test_icofoam_parallel():
    # Run the solver in parallel mode (should handle -parallel argument)
    result = subprocess.run([
        'mpirun', '-np', '4', sys.executable, '-m', 'foamadapter.cli.app', 'solver', 'icofoam', '-parallel'
    ], capture_output=True, text=True)
    assert result.returncode == 0
