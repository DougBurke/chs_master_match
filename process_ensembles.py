#!/usr/bin/env python

"""
Usage:

  ./process_ensembles.py
      svdqafile centroidfile ensemblestk stackmap outdir
      --mrgsrc3dir <dirname>
      --verbose

Aim:

Run through the initial CHS master-match stage - identification
of master hulls and a proposed shape - for the ensembles in the
ensembles stack. The stackmap file provides a mapping between
ensemble and stack ids (it has columns ensemble and stack),
and outdir is the base directory for the output. The output
directory will be created if need be.

The svdqafile argument points to a file which lists all the stacks
that went to SVD QA, one stack per line. The first column is used.
This file is used to fill out the STKSVDQA column.

The centroidfile is a list of stack,component,use_centroid lines
which indicate which stack-level hulls to use in the centroid
calculation. Any stack-level hull not in this file is assumed
to be included.

The outputs are
   <outdir>/<ensemble>
   <outdir>/log.<ensemble>

If the output directory exists then the ensemble will be skipped.
The log file will be over-written if it exists but the output
directory does not.

The --verbose flag just displays the command being run

"""

import os
from subprocess import CalledProcessError, check_output, STDOUT
import sys

import six

import stk


def doit(ensemblestk, stackmapfile, outdir,
         svdqafile, centroidfile, mrgsrc3dir="",
         verbose=False):
    """Process the ensembles to find master hulls.

    If the output directory for an ensemble exists it is skipped.

    Parameters
    ----------
    ensemblestk : str
        The stack of ensembles - e.g. "ens00000500_001,..."
        or "@ensmap.dat".
    stackmapfile : str
        The file must have columns ensemble and stack, with
        one stack per line.
    outdir : str
        The name of the output directory; it will be created
        if it does not exist.
    svdqafile : str
        The name of the file containing the stack ids
        that went to SVD QA. The first column of this file is used
        as the stack id. The full path is written to the SVDQAFIL
        header keyword.
    centroidfile : str or None, optional
        If given theh the name of the file containing the stack,cpt,
        include_centroid information (a partial list). Used to set
        up the USE_CEN column.
    mrgsrc3dir : str, optional
        Passed through as the --mrgsrc3dir value if not empty.
    verbose : bool, optional
        If set then displays the command being run.
    """

    # Assume the directory name does not have leading or trailing
    # white space.
    #
    mrgsrc3dir = mrgsrc3dir.strip()
    if mrgsrc3dir != "" and not os.path.isdir(mrgsrc3dir):
        raise IOError("Unable to find " +
                      "mrgsrc3dir={}".format(mrgsrc3dir))

    # Look for the tool in the same directory as this script
    dirname = os.path.dirname(__file__)
    toolname = os.path.join(dirname, "chs_create_initial_masters.py")
    if not os.path.exists(toolname):
        raise IOError("Unable to find '{}'".format(toolname))

    ensnames = stk.build(ensemblestk)
    ntot = len(ensnames)
    if ntot == 0:
        print("No ensembles found. That's surprising.")
        return

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise IOError("outdir '{}' ".format(outdir) +
                          "is not a directory!")

    else:
        print("Creating: {}".format(outdir))
        os.mkdir(outdir)

    failed = []
    for i, ensname in enumerate(ensnames):

        sys.stdout.write("[{}/{}] {}\n".format(i + 1, ntot, ensname))
        sys.stdout.flush()

        dirname = os.path.join(outdir, ensname)
        if os.path.exists(dirname):
            print("  - skipping as {} exists".format(dirname))
            if not os.path.isdir(dirname):
                print("    WARNING: not a directory!")

            continue

        logfile = os.path.join(outdir, 'log.' + ensname)
        if os.path.exists(logfile):
            os.remove(logfile)

        args = ['python', toolname, svdqafile, centroidfile,
                stackmapfile, ensname, dirname]
        if mrgsrc3dir != "":
            args.extend(["--mrgsrc3dir", mrgsrc3dir])

        if verbose:
            # assume no protection/quoting needed
            print(">> {}".format(" ".join(args)))

        try:
            out = check_output(args, stderr=STDOUT)
            if not six.PY2:
                out = out.decode('ascii')

        except CalledProcessError as exc:
            out = "ERROR: ensemble={}\n{}\n".format(ensname, exc) + \
                "\n"
            if six.PY2:
                out += exc.output
            else:
                out += exc.output.decode('ascii')

            sys.stdout.write("    FAILED\n")
            sys.stdout.flush()
            failed.append(ensname)

        with open(logfile, 'w') as fh:
            fh.write(out)

    nfail = len(failed)
    if nfail == 0:
        print("All ran successfully.")
        return

    print("")

    if nfail == ntot:
        print("*** They ALL failed!\n")
    elif nfail == 1:
        print("*** There was one failure:")
    else:
        print("*** There were {} failures:".format(nfail))

    if nfail != ntot:
        for i, ensname in enumerate(failed):
            print("  {}/{}  {}".format(i + 1, nfail, ensname))

    sys.exit(1)


help_str = "Run chs_create_initial_masters on a set of ensembles."

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("svdqafile",
                        help="The list of stacks that went to SVD QA")
    parser.add_argument("centroidfile",
                        help="stack,cpt,include_centroid information")
    parser.add_argument("ensemblestk",
                        help="The ensembles to process, as a stack")
    parser.add_argument("stackmap", type=str,
                        help="The mapping betweenensemble and stack ids")
    parser.add_argument("outdir", type=str,
                        help="The output directory (must not exist)")

    parser.add_argument("--mrgsrc3dir",
                        default="",
                        help="The mrgsrc3 directory: default %(default)s")

    parser.add_argument("--verbose", action="store_true",
                        help="Displays the command being run")

    args = parser.parse_args(sys.argv[1:])

    doit(args.ensemblestk, args.stackmap, args.outdir,
         svdqafile=args.svdqafile,
         centroidfile=args.centroidfile,
         mrgsrc3dir=args.mrgsrc3dir,
         verbose=args.verbose)
