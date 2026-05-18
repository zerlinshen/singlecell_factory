"""Unit tests for ModuleContractError (US-W5-FOUNDATION)."""
from workflow.modular._contract_violation import ModuleContractError


def test_import():
    assert ModuleContractError is not None


def test_instantiation():
    err = ModuleContractError("test")
    assert str(err) == "test"


def test_is_exception_subclass():
    assert issubclass(ModuleContractError, Exception)


def test_can_be_raised_and_caught():
    try:
        raise ModuleContractError("contract violated")
    except ModuleContractError as exc:
        assert "contract violated" in str(exc)
    else:
        raise AssertionError("ModuleContractError was not raised")
