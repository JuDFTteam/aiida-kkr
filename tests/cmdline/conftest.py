# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Fixtures for the command line interface.
Most of these fixtures are taken from the aiida-quantum espresso package since they
are needed to mock commands running aiida process
"""
import pytest


def mock_launch_process(*_, **__):
    """Mock the :meth:`~aiida_kkr.cmdline.util.utils.launch_process` to be a no-op."""
    return


@pytest.fixture
def run_cli_command():
    """Run a `click` command with the given options.

    The call will raise if the command triggered an exception or the exit code returned is non-zero.
    """

    def _run_cli_command(command, options=None, raises=None):
        """Run the command and check the result.

        :param command: the command to invoke
        :param options: the list of command line options to pass to the command invocation
        :param raises: optionally an exception class that is expected to be raised
        """
        import traceback
        from click.testing import CliRunner

        runner = CliRunner()
        result = runner.invoke(command, options or [])

        if raises is not None:
            assert result.exception is not None, result.output
            assert result.exit_code != 0
        else:
            assert result.exception is None, ''.join(traceback.format_exception(*result.exc_info))
            assert result.exit_code == 0, result.output

        result.output_lines = [line.strip() for line in result.output.split('\n') if line.strip()]

        return result

    return _run_cli_command


@pytest.fixture
def run_cli_process_launch_command(run_cli_command, monkeypatch):
    """Run a process launch command with the given options.

    The call will raise if the command triggered an exception or the exit code returned is non-zero.

    :param command: the command to invoke
    :param options: the list of command line options to pass to the command invocation
    :param raises: optionally an exception class that is expected to be raised
    """

    def _inner(command, options=None, raises=None):
        """Run the command and check the result."""
        from aiida_kkr.cmdline.util import utils
        monkeypatch.setattr(utils, 'launch_process', mock_launch_process)
        return run_cli_command(command, options, raises)

    return _inner


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return a `Code` instance configured to run calculations of given entry point on localhost `Computer`."""

    def _fixture_code(entry_point_name):
        from aiida.orm import Code
        return Code(input_plugin_name=entry_point_name, remote_computer_exec=[fixture_localhost, '/bin/ls'])

    return _fixture_code
