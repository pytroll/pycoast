"""The conftest file."""

from pytest import hookimpl


@hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Add test status in the report for fixtures to use."""
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()
    # store test results for each phase of a call, which can
    # be "setup", "call", "teardown"
    setattr(item, "rep_" + rep.when, rep)
