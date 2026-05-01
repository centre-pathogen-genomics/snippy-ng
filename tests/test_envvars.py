from pydantic import BaseModel

from snippy_ng.envvars import (
    EnvVar,
    bool_envvar,
    candidate_envvar_names,
    float_envvar,
    int_envvar,
    normalize_envvar_name,
    parse_bool_envvar,
    parse_float_envvar,
    parse_int_envvar,
    parse_str_envvar,
    str_envvar,
)


def test_parse_bool_envvar_uses_default_when_missing():
    assert parse_bool_envvar("SNIPPY_NG_MISSING_FLAG", default=True, environ={}) is True
    assert parse_bool_envvar("SNIPPY_NG_MISSING_FLAG", default=False, environ={}) is False


def test_parse_bool_envvar_accepts_truthy_and_falsy_values():
    assert parse_bool_envvar("FLAG", environ={"SNIPPY_NG_FLAG": "true"}) is True
    assert parse_bool_envvar("FLAG", environ={"SNIPPY_NG_FLAG": "ON"}) is True
    assert parse_bool_envvar("FLAG", environ={"SNIPPY_NG_FLAG": "0"}) is False
    assert parse_bool_envvar("FLAG", environ={"SNIPPY_NG_FLAG": " no "}) is False


def test_parse_bool_envvar_rejects_invalid_values():
    try:
        parse_bool_envvar("FLAG", environ={"SNIPPY_NG_FLAG": "sometimes"})
    except ValueError as exc:
        assert "Invalid value for SNIPPY_NG_FLAG" in str(exc)
    else:
        raise AssertionError("Expected ValueError for invalid boolean envvar value")


def test_parse_int_envvar_supports_default_and_conversion():
    assert parse_int_envvar("SNIPPY_NG_THREADS", default=8, environ={}) == 8
    assert parse_int_envvar("THREADS", environ={"SNIPPY_NG_THREADS": " 12 "}) == 12


def test_parse_float_envvar_supports_default_and_conversion():
    assert parse_float_envvar("SNIPPY_NG_RATIO", default=1.5, environ={}) == 1.5
    assert parse_float_envvar("RATIO", environ={"SNIPPY_NG_RATIO": " 2.75 "}) == 2.75


def test_parse_str_envvar_supports_default_and_raw_strings():
    assert parse_str_envvar("SNIPPY_NG_LABEL", default="default", environ={}) == "default"
    assert parse_str_envvar("LABEL", environ={"SNIPPY_NG_LABEL": "  keep spacing  "}) == "  keep spacing  "


def test_normalize_envvar_name_auto_prefixes_and_uppercases():
    assert normalize_envvar_name("debug") == "SNIPPY_NG_DEBUG"
    assert normalize_envvar_name("SNIPPY_NG_local_bcftools_csq") == "SNIPPY_NG_LOCAL_BCFTOOLS_CSQ"


def test_candidate_envvar_names_include_prefixed_and_unprefixed_forms():
    assert candidate_envvar_names("LOCAL_BCFTOOLS_CSQ") == ("SNIPPY_NG_LOCAL_BCFTOOLS_CSQ",)


def test_envvar_defaults_description_when_missing():
    test_bool = bool_envvar("TEST_BOOL")

    assert test_bool.description == "Value loaded from environment variable SNIPPY_NG_TEST_BOOL"


def test_define_envvars_read_current_environment(monkeypatch):
    test_bool = bool_envvar("TEST_BOOL", description="Toggle a test boolean")
    test_int = int_envvar("TEST_INT", default=3, description="Configure integer test value")
    test_float = float_envvar("TEST_FLOAT", default=1.0, description="Configure float test value")
    test_str = str_envvar("TEST_STR", default="fallback", description="Configure string test value")

    monkeypatch.setenv("SNIPPY_NG_TEST_BOOL", "1")
    monkeypatch.setenv("SNIPPY_NG_TEST_INT", "42")
    monkeypatch.setenv("SNIPPY_NG_TEST_FLOAT", "0.125")
    monkeypatch.setenv("SNIPPY_NG_TEST_STR", "configured")

    assert test_bool.env_var == "SNIPPY_NG_TEST_BOOL"
    assert test_int.env_var == "SNIPPY_NG_TEST_INT"
    assert test_bool.description == "Toggle a test boolean"
    assert test_bool.enabled is True
    assert test_int.value == 42
    assert test_float.value == 0.125
    assert test_str.value == "configured"
def test_envvar_field_reads_current_environment_on_model_creation(monkeypatch):
    class ExampleModel(BaseModel):
        enabled: bool = EnvVar("EXAMPLE_ENABLED", default=True)
        retries: int = EnvVar("EXAMPLE_RETRIES", default=3)

    monkeypatch.setenv("SNIPPY_NG_EXAMPLE_ENABLED", "0")
    monkeypatch.setenv("SNIPPY_NG_EXAMPLE_RETRIES", "8")

    model = ExampleModel()

    assert model.enabled is False
    assert model.retries == 8