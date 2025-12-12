"""Tests for the main Snippy class."""
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from snippy_ng.snippy import Snippy
from snippy_ng.stages.base import BaseStage, BaseOutput
from snippy_ng.dependencies import Dependency
from snippy_ng.exceptions import DependencyError, MissingOutputError
from pydantic import Field


class MockOutput(BaseOutput):
    test_file: Path


class MockStage(BaseStage):
    test_param: str = Field("test", description="Test parameter")
    _dependencies = [Dependency("test_tool", min_version="1.0.0")]

    @property
    def output(self):
        return MockOutput(test_file=Path(f"{self.prefix}.test"))

    @property
    def commands(self):
        return []


def test_snippy_init_empty():
    """Test Snippy initialization with no stages."""
    snippy = Snippy()
    assert snippy.stages == []


def test_snippy_init_with_stages():
    """Test Snippy initialization with stages."""
    stage = MockStage()
    snippy = Snippy(stages=[stage])
    assert len(snippy.stages) == 1
    assert snippy.stages[0] == stage


def test_snippy_add_stage():
    """Test adding a stage to Snippy."""
    snippy = Snippy()
    stage = MockStage()
    snippy.add_stage(stage)
    assert len(snippy.stages) == 1
    assert snippy.stages[0] == stage


def test_snippy_dependencies():
    """Test getting dependencies from stages."""
    stage1 = MockStage()
    stage2 = MockStage()
    snippy = Snippy(stages=[stage1, stage2])
    deps = snippy.dependencies
    assert len(deps) == 2  # Two stages, each with one dependency


def test_snippy_citations():
    """Test getting citations from dependencies."""
    dep1 = Dependency("tool1", citation="Citation 1")
    dep2 = Dependency("tool2", citation="Citation 2")
    dep3 = Dependency("tool3", citation="Citation 1")  # Duplicate
    
    stage1 = MockStage()
    stage1._dependencies = [dep1, dep2]
    stage2 = MockStage()
    stage2._dependencies = [dep3]
    
    snippy = Snippy(stages=[stage1, stage2])
    citations = snippy.citations
    
    assert len(citations) == 2  # Duplicates removed
    assert "Citation 1" in citations
    assert "Citation 2" in citations


def test_snippy_citations_excludes_none():
    """Test that citations excludes dependencies with no citation."""
    dep1 = Dependency("tool1", citation="Citation 1")
    dep2 = Dependency("tool2", citation=None)
    dep3 = Dependency("tool3", citation="")
    
    stage = MockStage()
    stage._dependencies = [dep1, dep2, dep3]
    
    snippy = Snippy(stages=[stage])
    citations = snippy.citations
    
    assert len(citations) == 1
    assert "Citation 1" in citations


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_snippy_validate_dependencies_success(mock_run, mock_which):
    """Test successful dependency validation."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 1.5.0\n")
    
    stage = MockStage()
    snippy = Snippy(stages=[stage])
    
    # Should not raise an exception
    snippy.validate_dependencies()


@patch('snippy_ng.dependencies.which')
def test_snippy_validate_dependencies_missing(mock_which):
    """Test dependency validation with missing dependency."""
    mock_which.return_value = None
    
    stage = MockStage()
    snippy = Snippy(stages=[stage])
    
    with pytest.raises(DependencyError):
        snippy.validate_dependencies()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_snippy_validate_dependencies_invalid_version(mock_run, mock_which):
    """Test dependency validation with invalid version."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 0.5.0\n")
    
    stage = MockStage()
    snippy = Snippy(stages=[stage])
    
    with pytest.raises(DependencyError):
        snippy.validate_dependencies()


@patch('snippy_ng.dependencies.which')
@patch('snippy_ng.dependencies.subprocess.run')
def test_snippy_validate_dependencies_skip_duplicate(mock_run, mock_which):
    """Test that duplicate dependencies are only checked once."""
    mock_which.return_value = "/usr/bin/test_tool"
    mock_run.return_value = MagicMock(stdout="test_tool version 1.5.0\n")
    
    stage1 = MockStage()
    stage2 = MockStage()
    snippy = Snippy(stages=[stage1, stage2])
    
    snippy.validate_dependencies()
    
    # Should only be called once for the duplicate dependency
    assert mock_run.call_count == 1


def test_snippy_run_no_stages():
    """Test running Snippy with no stages."""
    snippy = Snippy()
    
    with pytest.raises(ValueError, match="No stages to run"):
        snippy.run()


def test_snippy_run_missing_output(tmp_path):
    """Test that missing outputs raise an error."""
    stage = MockStage()
    snippy = Snippy(stages=[stage])
    
    with pytest.raises(MissingOutputError, match="Expected output 'test_file'"):
        snippy.run()


def test_snippy_set_working_directory(tmp_path):
    """Test setting working directory."""
    import os
    original_dir = os.getcwd()
    
    snippy = Snippy()
    snippy.set_working_directory(str(tmp_path))
    
    assert os.getcwd() == str(tmp_path)
    
    # Reset
    os.chdir(original_dir)


def test_snippy_welcome():
    """Test welcome message."""
    stage1 = MockStage()
    stage2 = MockStage()
    snippy = Snippy(stages=[stage1, stage2])
    
    # Should not raise an exception
    snippy.welcome()


def test_snippy_goodbye():
    """Test goodbye message."""
    dep = Dependency("tool1", citation="Citation 1")
    stage = MockStage()
    stage._dependencies = [dep]
    
    snippy = Snippy(stages=[stage])
    snippy.start_time = 0
    snippy.end_time = 10
    
    # Should not raise an exception
    snippy.goodbye()



