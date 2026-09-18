#!/usr/bin/env python
"""Upload path built by kkr_flex_wc.move_kkrflex_files (no database needed).

Regression guard for issue #180: on a computer whose work directory is a template such as
`/scratch/{username}/aiida_run`, the Green's function was uploaded to a directory named
literally `{username}`. The path is now expanded from the transport that the step opens
anyway. See also KkrimpCalculation.get_remote_symlink, which builds the same path on the
reading side.
"""

import os
from types import SimpleNamespace

from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
from aiida_kkr.calculations.kkrimp import KkrimpCalculation

REMOTE_PATH = '/scratch/someuser/aiida_run/abc/123'
UUID_RETRIEVED = '20ac6bc9-9fed-4c1e-beb0-a4c7282e1cd4'
WHOAMI = 'iff003user'


class _RecordingConnection:
    """Transport stand-in that records the directory it was asked to create."""

    def __init__(self):
        self.made_dirs = []

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False

    def whoami(self):
        return WHOAMI

    def isdir(self, path):
        return False

    def makedirs(self, path, ignore_existing=False):
        self.made_dirs.append(path)

    def symlink(self, source, destination):
        pass

    def copyfile(self, source, destination):
        pass

    def remove(self, path):
        pass


def _workchain_stub(workdir, connection):
    """Minimal stand-in for a kkr_flex_wc instance whose KKRFLEX calculation finished."""
    computer = SimpleNamespace(
        label='testcomputer',
        get_workdir=lambda: workdir,
        get_transport=lambda: connection,
    )
    flexrun = SimpleNamespace(
        is_finished_ok=True,
        computer=computer,
        outputs=SimpleNamespace(
            remote_folder=SimpleNamespace(get_remote_path=lambda: REMOTE_PATH),
            retrieved=SimpleNamespace(uuid=UUID_RETRIEVED),
        ),
    )
    return SimpleNamespace(
        ctx=SimpleNamespace(retrieve_kkrflex=False, flexrun=flexrun),
        report=lambda message: None,
    )


def test_upload_path_expands_username_template():
    """A templated work directory is expanded with the remote login name, not left literal."""
    connection = _RecordingConnection()
    kkr_flex_wc.move_kkrflex_files(_workchain_stub('/scratch/{username}/aiida_run', connection))

    expected = os.path.join('/scratch', WHOAMI, 'aiida_run', KkrimpCalculation._DIRNAME_GF_UPLOAD, UUID_RETRIEVED)
    assert connection.made_dirs == [expected]
    assert '{username}' not in connection.made_dirs[0]


def test_upload_path_unchanged_for_literal_workdir():
    """A work directory without a placeholder is passed through untouched.

    This is the backwards-compatibility assertion: every setup that works today has a literal
    work directory, and str.format is the identity function on a string with no placeholder,
    so the path such a setup gets must be byte-identical to the one it got before the fix.
    """
    connection = _RecordingConnection()
    kkr_flex_wc.move_kkrflex_files(_workchain_stub('/work/jara0191/jw782093/aiida', connection))

    expected = os.path.join('/work/jara0191/jw782093/aiida', KkrimpCalculation._DIRNAME_GF_UPLOAD, UUID_RETRIEVED)
    assert connection.made_dirs == [expected]


def test_no_upload_when_files_are_retrieved():
    """With retrieve_kkrflex set, the step returns before opening any connection."""
    connection = _RecordingConnection()
    stub = _workchain_stub('/scratch/{username}/aiida_run', connection)
    stub.ctx.retrieve_kkrflex = True
    kkr_flex_wc.move_kkrflex_files(stub)

    assert connection.made_dirs == []


# run test manually
if __name__ == '__main__':
    test_upload_path_expands_username_template()
    test_upload_path_unchanged_for_literal_workdir()
    test_no_upload_when_files_are_retrieved()
    print('ok')
