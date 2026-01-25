from pathlib import Path
from typing import Optional
import json


class Metadata:
    """Class to hold metadata information for a Snippy-ng run."""

    def __init__(self, path: Optional[Path] = None, **metadata_overrides):
        """Initialize Metadata with optional path to metadata file."""
        self.path = path
        self._static_data = metadata_overrides

    @property
    def reference(self) -> str:
        """Get the reference file name from metadata."""
        return self.get("reference")

    @property
    def format(self) -> str:
        """Get the reference format from metadata."""
        return self.get("format")

    @property
    def num_sequences(self) -> int:
        """Get the number of sequences in the reference from metadata."""
        return self.get("num_sequences")

    @property
    def total_length(self) -> int:
        """Get the total length of the reference from metadata."""
        return self.get("total_length")

    @property
    def num_features(self) -> int:
        """Get the number of features in the reference from metadata."""
        return self.get("num_features")

    @property
    def prefix(self) -> str:
        """Get the prefix used for the reference from metadata."""
        return self.get("prefix")

    @property
    def datetime(self) -> str:
        """Get the datetime of the reference preparation from metadata."""
        return self.get("datetime")

    def get(self, key: str):
        """Get a metadata value by key."""
        if key in self._static_data:
            return self._static_data[key]

        with open(self.path, "r") as f:
            data = json.load(f)
        return data.get(key)
