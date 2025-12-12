"""Tests for dependencies module."""
import pytest
from unittest.mock import patch, MagicMock
from snippy_ng.dependencies import Dependency, PythonDependency
from snippy_ng.exceptions import (
    MissingDependencyError,
    InvalidDependencyError,
    InvalidDependencyVersionError,
)


def test_dependency_format_version_requirements_no_constraints():
    """Test formatting version requirements with no constraints."""
    dep = Dependency("test_tool")
    assert dep.format_version_requirements() == "test_tool"


def test_dependency_format_version_requirements_exact_version():
    """Test formatting version requirements with exact version."""
    dep = Dependency("test_tool", version="1.2.3")
    assert dep.format_version_requirements() == "test_tool =1.2.3"


def test_dependency_format_version_requirements_min_version():
    """Test formatting version requirements with minimum version."""
    dep = Dependency("test_tool", min_version="1.0.0")
    assert dep.format_version_requirements() == "test_tool >=1.0.0"


def test_dependency_format_version_requirements_max_version():
    """Test formatting version requirements with maximum version."""
    dep = Dependency("test_tool", max_version="2.0.0")
    assert dep.format_version_requirements() == "test_tool <=2.0.0"


def test_dependency_format_version_requirements_less_than():
    """Test formatting version requirements with less than constraint."""
    dep = Dependency("test_tool", less_then="2.0.0")
    assert dep.format_version_requirements() == "test_tool <2.0.0"


def test_dependency_format_version_requirements_multiple():
    """Test formatting version requirements with multiple constraints."""
    dep = Dependency("test_tool", min_version="1.0.0", max_version="2.0.0")
    result = dep.format_version_requirements()
    assert "test_tool" in result
    assert ">=1.0.0" in result
    assert "<=2.0.0" in result


@patch('snippy_ng.dependencies.which')
def test_dependency_check_missing(mock_which):
    """Test checking a missing dependency."""
    mock_which.return_value = None
    dep = Dependency("missing_tool")
    
    with pytest.raises(MissingDependencyError, match="Could not find dependency missing_tool"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_check_valid_version(mock_run, mock_which):
    """Test checking a dependency with valid version."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 1.2.3\n")
    
    dep = Dependency("test_tool", min_version="1.0.0")
    version = dep.check()
    assert str(version) == "1.2.3"


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_check_version_too_old(mock_run, mock_which):
    """Test checking a dependency with version too old."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 0.9.0\n")
    
    dep = Dependency("test_tool", min_version="1.0.0")
    with pytest.raises(InvalidDependencyError, match="minimum version allowed is 1.0.0"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_check_version_too_new(mock_run, mock_which):
    """Test checking a dependency with version too new."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 2.5.0\n")
    
    dep = Dependency("test_tool", max_version="2.0.0")
    with pytest.raises(InvalidDependencyError, match="maximum version allowed is 2.0.0"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_check_version_less_than(mock_run, mock_which):
    """Test checking a dependency with less than constraint."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 2.0.0\n")
    
    dep = Dependency("test_tool", less_then="2.0.0")
    with pytest.raises(InvalidDependencyError, match="version must be less than 2.0.0"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_check_exact_version_mismatch(mock_run, mock_which):
    """Test checking a dependency with exact version mismatch."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 1.2.4\n")
    
    dep = Dependency("test_tool", version="1.2.3")
    with pytest.raises(InvalidDependencyError, match="version must be 1.2.3"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_get_version_no_match(mock_run, mock_which):
    """Test getting version when pattern doesn't match."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="weird output\n")
    
    dep = Dependency("test_tool")
    with pytest.raises(InvalidDependencyVersionError, match="Could not extract version"):
        dep.check()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_custom_version_pattern(mock_run, mock_which):
    """Test dependency with custom version pattern."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="v2.17-r123\n")
    
    dep = Dependency("test_tool", version_pattern=r"v(\d+\.\d+)", min_version="2.0")
    version = dep.check()
    assert str(version) == "2.17"


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_dependency_no_version_arg(mock_run, mock_which):
    """Test dependency with no version argument."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="Version: 1.5.0\n")
    
    dep = Dependency("test_tool", version_arg=None)
    version = dep.check()
    assert str(version) == "1.5.0"


@patch('importlib.metadata.version')
def test_python_dependency_found(mock_version):
    """Test checking a Python dependency that exists."""
    mock_version.return_value = "1.2.3"
    
    dep = PythonDependency("test_package")
    version = dep.check()
    assert str(version) == "1.2.3"


@patch('importlib.metadata.version')
def test_python_dependency_not_found(mock_version):
    """Test checking a Python dependency that doesn't exist."""
    from importlib.metadata import PackageNotFoundError
    mock_version.side_effect = PackageNotFoundError()
    
    dep = PythonDependency("missing_package")
    with pytest.raises(MissingDependencyError, match="Could not find dependency missing_package"):
        dep.check()


@patch('importlib.metadata.version')
def test_python_dependency_version_constraints(mock_version):
    """Test Python dependency with version constraints."""
    mock_version.return_value = "1.0.0"
    
    dep = PythonDependency("test_package", min_version="1.5.0")
    with pytest.raises(InvalidDependencyError, match="minimum version allowed is 1.5.0"):
        dep.check()


def test_dependency_invalid_version_parse():
    """Test dependency with unparseable version."""
    dep = Dependency("test_tool")
    
    with pytest.raises(InvalidDependencyVersionError, match="Could not parse version"):
        dep._base_validator("not-a-version")
