from __future__ import annotations

import os
from typing import Any, Callable, Generic, Mapping, TypeVar

from pydantic.fields import FieldInfo
from pydantic_core import PydanticUndefined


TRUE_VALUES = {"1", "true", "t", "yes", "y", "on"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "off", ""}
ENV_PREFIX = "SNIPPY_NG_"

T = TypeVar("T")


def normalize_envvar_name(name: str) -> str:
    normalized = name.strip().upper()
    if normalized.startswith(ENV_PREFIX):
        return normalized
    return f"{ENV_PREFIX}{normalized}"


def candidate_envvar_names(name: str) -> tuple[str, ...]:
    return (normalize_envvar_name(name),)


def _get_env_value(env_var: str, *, environ: Mapping[str, str] | None = None) -> str | None:
    source = os.environ if environ is None else environ
    for candidate in candidate_envvar_names(env_var):
        if candidate in source:
            return source[candidate]
    return None


def _parse_bool_raw(env_var: str, raw_value: str) -> bool:
    normalized = raw_value.strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False

    allowed = ", ".join(sorted(TRUE_VALUES | FALSE_VALUES))
    raise ValueError(f"Invalid value for {normalize_envvar_name(env_var)}: {raw_value!r}. Expected one of: {allowed}")


def _parse_int_raw(env_var: str, raw_value: str) -> int:
    try:
        return int(raw_value.strip())
    except ValueError as exc:
        raise ValueError(f"Invalid integer value for {normalize_envvar_name(env_var)}: {raw_value!r}") from exc


def _parse_float_raw(env_var: str, raw_value: str) -> float:
    try:
        return float(raw_value.strip())
    except ValueError as exc:
        raise ValueError(f"Invalid float value for {normalize_envvar_name(env_var)}: {raw_value!r}") from exc


def parse_bool_envvar(
    env_var: str,
    *,
    default: bool = False,
    environ: Mapping[str, str] | None = None,
) -> bool:
    raw_value = _get_env_value(env_var, environ=environ)
    if raw_value is None:
        return default
    return _parse_bool_raw(env_var, raw_value)


def parse_int_envvar(
    env_var: str,
    *,
    default: int = 0,
    environ: Mapping[str, str] | None = None,
) -> int:
    raw_value = _get_env_value(env_var, environ=environ)
    if raw_value is None:
        return default
    return _parse_int_raw(env_var, raw_value)


def parse_float_envvar(
    env_var: str,
    *,
    default: float = 0.0,
    environ: Mapping[str, str] | None = None,
) -> float:
    raw_value = _get_env_value(env_var, environ=environ)
    if raw_value is None:
        return default
    return _parse_float_raw(env_var, raw_value)


def parse_str_envvar(
    env_var: str,
    *,
    default: str = "",
    environ: Mapping[str, str] | None = None,
) -> str:
    raw_value = _get_env_value(env_var, environ=environ)
    if raw_value is None:
        return default
    return raw_value


def _default_description(env_var: str) -> str:
    return f"Value loaded from environment variable {normalize_envvar_name(env_var)}"


def _infer_parser(env_var: str, default: Any) -> Callable[[str], Any]:
    if isinstance(default, bool):
        return lambda raw_value: _parse_bool_raw(env_var, raw_value)
    if isinstance(default, int):
        return lambda raw_value: _parse_int_raw(env_var, raw_value)
    if isinstance(default, float):
        return lambda raw_value: _parse_float_raw(env_var, raw_value)
    if isinstance(default, str):
        return lambda raw_value: raw_value
    raise TypeError(
        f"Cannot infer parser for {normalize_envvar_name(env_var)} from default value {default!r}. "
        "Provide a supported default type or use a typed helper like bool_envvar()."
    )


class EnvVar(FieldInfo, Generic[T]):
    env_var: str
    fallback_default: T | Any
    fallback_default_factory: Callable[[], T] | None
    parser: Callable[[str], T]

    def __init__(
        self,
        env_var: str,
        default: T | Any = PydanticUndefined,
        *,
        default_factory: Callable[[], T] | None = None,
        description: str | None = None,
        parser: Callable[[str], T] | None = None,
        **kwargs,
    ) -> None:
        normalized_env_var = normalize_envvar_name(env_var)
        inferred_default = default_factory() if default is PydanticUndefined and default_factory is not None else default

        self.env_var = normalized_env_var
        self.fallback_default = default
        self.fallback_default_factory = default_factory
        self.parser = parser or _infer_parser(env_var, inferred_default)

        super().__init__(
            default_factory=self.get,
            description=description or _default_description(env_var),
            **kwargs,
        )

    def get(self, *, environ: Mapping[str, str] | None = None) -> T:
        raw_value = _get_env_value(self.env_var, environ=environ)
        if raw_value is not None:
            return self.parser(raw_value)
        if self.fallback_default_factory is not None:
            return self.fallback_default_factory()
        if self.fallback_default is not PydanticUndefined:
            return self.fallback_default
        raise ValueError(f"Environment variable {self.env_var} is required")

    @property
    def value(self) -> T:
        return self.get()


class BoolEnvVar(EnvVar[bool]):

    @property
    def enabled(self) -> bool:
        return self.get()

    def __bool__(self) -> bool:
        return self.get()


def envvar(
    env_var: str,
    *,
    default: T,
    description: str | None = None,
    parser: Callable[[str], T],
) -> EnvVar[T]:
    return EnvVar(env_var=env_var, default=default, description=description, parser=parser)


def bool_envvar(env_var: str, *, default: bool = False, description: str | None = None) -> BoolEnvVar:
    return BoolEnvVar(
        env_var=env_var,
        default=default,
        description=description,
        parser=lambda raw_value: _parse_bool_raw(env_var, raw_value),
    )


def int_envvar(env_var: str, *, default: int = 0, description: str | None = None) -> EnvVar[int]:
    return envvar(
        env_var,
        default=default,
        description=description,
        parser=lambda raw_value: _parse_int_raw(env_var, raw_value),
    )


def float_envvar(env_var: str, *, default: float = 0.0, description: str | None = None) -> EnvVar[float]:
    return envvar(
        env_var,
        default=default,
        description=description,
        parser=lambda raw_value: _parse_float_raw(env_var, raw_value),
    )


def str_envvar(env_var: str, *, default: str = "", description: str | None = None) -> EnvVar[str]:
    return envvar(
        env_var,
        default=default,
        description=description,
        parser=lambda raw_value: raw_value,
    )