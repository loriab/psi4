#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

# /usr/bin/env python
# vim:ft=python

# reap-sow helper script
# executed after the sow step
# AJ
import os
import re
import sys
import subprocess
import time

root_directory = os.path.dirname(os.path.realpath(__file__))


def sowList(first_out):
    """read the output from the 'sow' step to the list of files
    to run before the 'reap' step

    """
    find_cmd = re.compile("^(#|\s)\s+ psi4 -i (?P<infile>(?P<tag>[a-zA-Z]+)-[a-z0-9]+(-[a-z0-9]+)?\.in)\s+-o (?P<outfile>[a-zA-Z]+-[a-z0-9]+(-[a-z0-9]+)?\.out)")  # noqa: E501
    the_list = []
    with open(first_out, 'r') as sow_out:
        for line in sow_out:
            matchobj = find_cmd.match(line)
            if (matchobj):
                the_in = matchobj.group('infile')
                the_out = matchobj.group('outfile')
                the_tag = matchobj.group('tag')
                the_list.append((the_in, the_out))
    if the_list:
        master_in = the_tag + "-master.in"
        master_out = the_tag + "-master.out"
        return the_list, master_in, master_out
    else:
        return None


def storeTests(testDir):
    """Store the tests in string"""
    testFile = os.path.join(testDir, "tests")
    with open(testFile) as F:
        retTests = F.read()
    return retTests


def addTests(theMaster, theTests):
    """ Append the tests to the master file, it must be created first"""
    with open(theMaster, 'a') as F:
        F.write(theTests)


def runFiles(psi4, theFileList, runningDir, psi4datadir):
    """ run the input files generated by the first input """
    for infile, outfile in theFileList:
        thisInFile = os.path.join(runningDir, infile.strip())
        thisOutFile = os.path.join(runningDir, outfile.strip())
        cmd = [psi4, "-i", thisInFile, "-o", thisOutFile, "-l", psi4datadir]
        retcode = subprocess.call(cmd)
        if retcode is not None:
            sys.stdout.write('Multi-invocation %s exited with status %i\n' % (cmd, retcode))
            if retcode != 0:
                sys.exit(1)


def runMaster(psi4, inMasterFile, outMasterFile, logfile, psi4datadir, append=False):
    """ run the master file and record the output in the logfile,
    the same way that runtest.py does it.

    """
    cmd = [psi4, "-i", inMasterFile, "-o", outMasterFile, "-l", psi4datadir]
    if append:
        cmd.append("-a")
    try:
        loghandle = open(logfile, 'a')
    except IOError as e:
        print("""%s can't write to %s: %s""" % (__name__, logfile, e))
        sys.exit(1)
    try:
        retcode = subprocess.Popen(cmd, bufsize=0, stdout=subprocess.PIPE, universal_newlines=True)
    except OSError as e:
        sys.stderr.write('Command %s execution failed: %s\n' % cmd, e.strerror)
        sys.exit(1)
    p4out = ''
    while True:
        data = retcode.stdout.readline()
        if not data:
            break
        sys.stdout.write(data)  # screen
        loghandle.write(data)  # file
        loghandle.flush()
        p4out += data  # string
    loghandle.close()
    while True:
        retcode.poll()
        exstat = retcode.returncode
        if exstat is not None:
            if exstat != 0:
                sys.exit(exstat)
            return exstat
        time.sleep(0.1)


def main(first_input, first_output, logfile, psi4, psi4datadir):
    # check the psi4 path is there and is executable
    input_dir = os.path.dirname(first_input)
    output_dir = os.path.dirname(first_output)
    if os.path.isfile(psi4) and os.access(psi4, os.X_OK):
        # get the list of intermediate input files
        # files list and master file are not abolute paths just names
        result = sowList(first_output)
        if result is None:
            # OPT mode
            # get the commands from tests
            tests = storeTests(input_dir)
            oiter = 0
            while oiter < 20:
                result = sowList('OPT-master.in')
                if result is None:
                    break
                else:
                    files_list, master_in, master_out = result
                    if files_list[-1] == ('OPT-master.in', 'OPT-master.out'):
                        files_list.pop()
                    # set the "reapmode" master file full path
                    reap_master_in = os.path.join(output_dir, master_in.strip())
                    reap_master_out = os.path.join(output_dir, master_out.strip())
                    # append the tests
                    addTests(reap_master_in, tests)
                    # run intermediates
                    runFiles(psi4, files_list, output_dir, psi4datadir)
                    # run the master
                    runMaster(psi4, reap_master_in, reap_master_out, logfile,
                              psi4datadir, append=False if oiter == 0 else True)
                oiter += 1
            # success
            sys.exit(0)
        else:
            files_list, master_in, master_out = result
            # get the commands from tests
            tests = storeTests(input_dir)
            # set the "reapmode" master file full path
            reap_master_in = os.path.join(output_dir, master_in.strip())
            reap_master_out = os.path.join(output_dir, master_out.strip())
            # append the tests
            addTests(reap_master_in, tests)
            # run intermediates
            runFiles(psi4, files_list, output_dir, psi4datadir)
            # run the master
            runMaster(psi4, reap_master_in, reap_master_out, logfile, psi4datadir)
            # success
            sys.exit(0)
    else:
        print("""check psi4 path  %s is not executable, or was not found""" % (psi4))
        sys.exit(1)


if __name__ == '__main__':
    main(*sys.argv[1:])
