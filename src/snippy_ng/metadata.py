from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Dict
import json


@dataclass
class ReferenceMetadata:
    path: Optional[Path] = None

    _overrides: Dict[str, Any] = field(default_factory=dict, init=False, repr=False)
    _data: Optional[Dict[str, Any]] = field(default=None, init=False, repr=False)
    _loaded: bool = field(default=False, init=False, repr=False)

    FIELDS = (
        "reference",
        "format",
        "num_sequences",
        "total_length",
        "num_features",
        "prefix",
        "datetime",
        "version",
    )

    def __init__(self, path: Optional[Path] = None, **kwargs: Any):
        self.path = path
        self._overrides = {}
        self._data = None
        self._loaded = False

        # Treat constructor kwargs as runtime overrides (like your old _reference fields)
        for k, v in kwargs.items():
            if k not in self.FIELDS:
                raise TypeError(f"Unknown metadata field: {k}")
            self._overrides[k] = v

    def _ensure_loaded(self) -> None:
        if self._loaded:
            return
        self._loaded = True

        if self.path is None:
            self._data = {}
            return

        try:
            with self.path.open("r", encoding="utf-8") as f:
                obj = json.load(f)
            if not isinstance(obj, dict):
                raise ValueError("Metadata JSON must be an object/dict.")
            self._data = obj
        except FileNotFoundError:
            self._data = {}
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in metadata file: {self.path}") from e

    def get(self, key: str, default: Any = None) -> Any:
        if key not in self.FIELDS:
            raise KeyError(f"Unknown metadata field: {key}")

        if key in self._overrides:
            return self._overrides[key]

        self._ensure_loaded()
        assert self._data is not None
        return self._data.get(key, default)

    def set(self, key: str, value: Any, *, persist: bool = False) -> None:
        if key not in self.FIELDS:
            raise KeyError(f"Unknown metadata field: {key}")
        self._overrides[key] = value
        if persist:
            self.write()

    def write(self) -> None:
        """
        Write metadata to disk.
        This merges cached file data with overrides, but gives overrides precedence.
        """
        if self.path is None:
            raise ValueError("Metadata path is not set.")

        self._ensure_loaded()
        assert self._data is not None

        merged = dict(self._data)
        merged.update(self._overrides)

        # Only write known fields (prevents accidental junk keys)
        out = {k: merged.get(k) for k in self.FIELDS}

        tmp = self.path.with_suffix(self.path.suffix + ".tmp")
        tmp.parent.mkdir(parents=True, exist_ok=True)

        with tmp.open("w", encoding="utf-8") as f:
            json.dump(out, f, indent=2, sort_keys=True)
            f.write("\n")

        tmp.replace(self.path)

        # sync internal state
        self._data = out
        self._overrides.clear()

    def validate(self) -> None:
        missing = [k for k in self.FIELDS if self.get(k) is None]
        if missing:
            raise ValueError(f"Missing required metadata field(s): {', '.join(missing)}")

    # convenience properties
    @property
    def reference(self) -> Optional[str]:
        return self.get("reference")

    @property
    def format(self) -> Optional[str]:
        return self.get("format")

    @property
    def num_sequences(self) -> Optional[int]:
        return self.get("num_sequences")

    @property
    def total_length(self) -> Optional[int]:
        return self.get("total_length")

    @property
    def num_features(self) -> Optional[int]:
        return self.get("num_features")

    @property
    def prefix(self) -> Optional[str]:
        return self.get("prefix")

    @property
    def datetime(self) -> Optional[str]:
        return self.get("datetime")

    @property
    def version(self) -> Optional[str]:
        return self.get("version")
