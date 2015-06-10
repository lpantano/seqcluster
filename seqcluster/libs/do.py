"""Centralize running of external commands, providing logging and tracking. Integrated from bcbio package with some changes.
"""
import collections
import os
import subprocess
import logging

logger = logging.getLogger("run")


def run(cmd, data=None, checks=None, region=None, log_error=True,
        log_stdout=False):
    """Run the provided command, logging details and checking for errors.
    """
    try:
        logger.debug(" ".join(str(x) for x in cmd) if not isinstance(cmd, basestring) else cmd)
        _do_run(cmd, checks, log_stdout)
    except:
        if log_error:
            logger.info("error at command")
        raise

def find_bash():
    for test_bash in [find_cmd("bash"), "/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")

def find_cmd(cmd):
    try:
        return subprocess.check_output(["which", cmd]).strip()
    except subprocess.CalledProcessError:
        return None

def _normalize_cmd_args(cmd):
    """Normalize subprocess arguments to handle list commands, string and pipes.
    Piped commands set pipefail and require use of bash to help with debugging
    intermediate errors.
    """
    if isinstance(cmd, basestring):
        # check for standard or anonymous named pipes
        if cmd.find(" | ") > 0 or cmd.find(">(") or cmd.find("<("):
            return "set -o pipefail; " + cmd, True, find_bash()
        else:
            return cmd, True, None
    else:
        return [str(x) for x in cmd], False, None

def _do_run(cmd, checks, log_stdout=False):
    """Perform running and check results, raising errors for issues.
    """
    cmd, shell_arg, executable_arg = _normalize_cmd_args(cmd)
    s = subprocess.Popen(cmd, shell=shell_arg, executable=executable_arg,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, close_fds=True)
    debug_stdout = collections.deque(maxlen=100)
    while 1:
        line = s.stdout.readline()
        if line:
            debug_stdout.append(line)
            if log_stdout:
                logger.debug(line.rstrip())
            else:
                logger.debug(line.rstrip())
        exitcode = s.poll()
        if exitcode is not None:
            for line in s.stdout:
                debug_stdout.append(line)
            if exitcode is not None and exitcode != 0:
                error_msg = " ".join(cmd) if not isinstance(cmd, basestring) else cmd
                error_msg += "\n"
                error_msg += "".join(debug_stdout)
                s.communicate()
                s.stdout.close()
                raise subprocess.CalledProcessError(exitcode, error_msg)
            else:
                break
    s.communicate()
    s.stdout.close()
    # Check for problems not identified by shell return codes
    if checks:
        for check in checks:
            if not check():
                raise IOError("External command failed")
