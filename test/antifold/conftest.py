#!/usr/bin/env python3

"""
Pytest configuration file for AntiFold tests.
"""

import os
import shutil
import tempfile
from datetime import datetime

import pytest
import torch

from rfantibody.config import PathConfig

# Option is defined in the root conftest.py

# Get test paths for this module
_test_paths = PathConfig.get_test_paths('antifold')


@pytest.fixture(scope="session", autouse=True)
def check_gpu():
    """Check if CUDA is available.

    AntiFold tests validate interfaces and inference correctness, not bit-exact
    reference outputs, so any CUDA-capable GPU is acceptable. GPU-specific
    reference directories are still honored by the ``ref_dir`` fixture when an
    A4000 or H100 is detected.
    """
    if not torch.cuda.is_available():
        pytest.skip("No GPU found, AntiFold tests require a CUDA-capable GPU")

    gpu_info = torch.cuda.get_device_properties(0)
    print(f"Running tests on {gpu_info.name} GPU")
    if 'A4000' in gpu_info.name:
        print("Using A4000-specific reference outputs")
    elif 'H100' in gpu_info.name:
        print("Using H100-specific reference outputs")
    else:
        print("Using generic reference outputs (no GPU-specific set for this device)")


@pytest.fixture(scope="session")
def output_dir(request):
    """
    Create and provide a temporary directory for test results.

    By default, uses a system temporary directory that will be automatically
    cleaned up. If --keep-outputs is specified, saves to a timestamped directory.
    """
    keep_outputs = request.config.getoption("--keep-outputs", default=False)

    if keep_outputs:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = _test_paths['outputs'].parent / f"example_outputs_{timestamp}"
        os.makedirs(output_path, exist_ok=True)
        print(f"Saving test outputs to: {output_path}")
        return str(output_path)
    else:
        temp_dir = tempfile.TemporaryDirectory(prefix="rfantibody_antifold_test_")
        request.config._rfantibody_antifold_temp_dir = temp_dir
        return temp_dir.name


@pytest.fixture(scope="session")
def ref_dir():
    """
    Provide the reference directory path based on GPU type.

    Uses GPU-specific references when running on a supported GPU (A4000 or H100).
    """
    base_ref_dir = _test_paths['references']

    if torch.cuda.is_available():
        gpu_info = torch.cuda.get_device_properties(0)
        if 'A4000' in gpu_info.name:
            return str(base_ref_dir / "A4000_references")
        elif 'H100' in gpu_info.name:
            return str(base_ref_dir / "H100_references")

    return str(base_ref_dir)


@pytest.fixture(scope="session")
def clean_output_dir(output_dir):
    """Clean the output directory before tests"""
    if os.path.exists(output_dir):
        for f in os.listdir(output_dir):
            file_path = os.path.join(output_dir, f)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Error cleaning {file_path}: {e}")

    yield
