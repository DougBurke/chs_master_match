#!/usr/bin/env python

"""
Usage:

  ./process_images.py indir

Aim:

Run chs_create_review_images_mpl.py on each indir/<ensemble>/. The
output is logged to indir/<ensemble>/log.plot.


"""

import glob
import os
from subprocess import CalledProcessError, check_output, STDOUT
import sys

import six


def doit(indir, mrgsrc3dir="", stkevt3dir="", stkfov3dir="",
         xmdat3dir="", verbose=False):
    """Process the ensembles to create the review images.

    Parameters
    ----------
    indir : str
        Process indir/<ensemble>/master*fits.
    mrgsrc3dir : str, optional
        Passed through as the --mrgsrc3dir value if not empty.
    stkevt3dir : str, optional
        Passed through as the --stkevt3dir value if not empty.
    stkfov3dir : str, optional
        Passed through as the --stkfov3dir value if not empty.
    xmdat3dir : str, optional
        Passed through as the --xmdat3dir value if not empty.
    verbose : bool, optional
        If set then displays the command being run.

    """

    # Look for the tool in the same directory as this script
    dirname = os.path.dirname(__file__)
    toolname = os.path.join(dirname, "chs_create_review_images_mpl.py")
    if not os.path.exists(toolname):
        raise IOError("Unable to find '{}'".format(toolname))

    infiles = glob.glob(os.path.join(indir, "*", "master*fits"))
    ntot = len(infiles)
    if ntot == 0:
        print("No ensemble master-hull files found. That's surprising.")
        return

    infiles = sorted(infiles)

    failed = []
    for i, infile in enumerate(infiles):

        outdir, _ = os.path.split(infile)
        _, ensemble = os.path.split(outdir)

        sys.stdout.write("[{}/{}] {}\n".format(i + 1, ntot, ensemble))
        sys.stdout.flush()

        logfile = os.path.join(outdir, 'log.plot')
        if os.path.exists(logfile):
            os.remove(logfile)

        args = ['python', toolname, infile, outdir]

        for pval, pname in [(mrgsrc3dir, 'mrgsrc3dir'),
                            (stkevt3dir, 'stkevt3dir'),
                            (stkfov3dir, 'stkfov3dir'),
                            (xmdat3dir, 'xmdat3dir'),
                            ]:
            if pval != "":
                args.extend(["--{}".format(pname), pval])

        if verbose:
            # assume no protection/quoting needed
            print(">> {}".format(" ".join(args)))
        try:
            out = check_output(args, stderr=STDOUT)
            if not six.PY2:
                out = out.decode('ascii')

        except CalledProcessError as exc:
            out = "ERROR: ensemble={}\n{}\n".format(ensemble, exc) + \
                "\n"
            if six.PY2:
                out += exc.output
            else:
                out += exc.output.decode('ascii')

            sys.stdout.write("    FAILED\n")
            sys.stdout.flush()
            failed.append(ensemble)

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
        for i, ensemble in enumerate(failed):
            print("  {}/{}  {}".format(i + 1, nfail, ensemble))

    sys.exit(1)


help_str = "Run chs_create_review_images_mpl on a set of ensembles."

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description=help_str,
                                     prog=sys.argv[0])

    parser.add_argument("indir",
                        help="The input directory containing the ensembles")

    parser.add_argument("--mrgsrc3dir",
                        default="",
                        help="The mrgsrc3 directory: default %(default)s")
    parser.add_argument("--stkevt3dir",
                        default="",
                        help="The stkevt3 directory: default %(default)s")
    parser.add_argument("--stkfov3dir",
                        default="",
                        help="The stkfov3 directory: default %(default)s")
    parser.add_argument("--xmdat3dir",
                        default="",
                        help="The xmdat3 directory: default %(default)s")

    parser.add_argument("--verbose", action="store_true",
                        help="Displays the command being run")

    args = parser.parse_args(sys.argv[1:])

    doit(args.indir,
         mrgsrc3dir=args.mrgsrc3dir,
         stkevt3dir=args.stkevt3dir,
         stkfov3dir=args.stkfov3dir,
         xmdat3dir=args.xmdat3dir,
         verbose=args.verbose)
