#!/usr/bin/env python

import os
import sys
import argparse

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()


def sys_call(cmd):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    output = ""
    while not proc.poll():
        line = process_readline(proc.stdout.readline())
        if line:
            if log:
                log.info(line)
            else:
                output += line + "\n"
        if proc.returncode is not None:
            break

    for line in proc.stdout.readlines():
        line = process_readline(line)
        if line:
            if log:
                log.info(line)
            else:
                output += line + "\n"

    if proc.returncode:
        log.error("system call for: \"%s\" finished abnormally, OS return value: %d" % (cmd, proc.returncode))
    return output


def parse_args():
    parser = argparse.ArgumentParser(description='Monomer Inference Problem: complement monomers set')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('monomers', help='fasta-file with monomers')
    parser.add_argument('-o', '--out-dir', dest="outdir", help='output directory (default=\'.\')', default=".", required=False)
    parser.add_argument('--len', help='the monomer length (default=171)', type=int, default=171, required=False)
    parser.add_argument('--lenDiv', '--max-length-divergence',
                        help='the maximum differ from length (default= 0.02*len)',
                        type=int, default=-1, required=False)
    parser.add_argument('--resDiv', '--max-resolved-divergence', help='max divergence in identity for resolve block (default=5%)',
                        type=float, default=5, required=False)
    parser.add_argument('--maxDiv','--max-divergence', help='max divergence in identity for monomeric-block (default=25%)',
                        type=float, default=25, required=False)
    return parser.parse_args()


def main():
    log.log("Start Monomer Inference")
    args = parse_args()
    if args.lenDiv == -1:
        args.lenDiv = int(args.len * 0.02)
    args.sequences = os.path.abspath(args.sequences)
    args.monomers = os.path.abspath(args.monomers)

    #get path to current script
    current_script_path = os.path.abspath(__file__)

    #path to the run_decomposer script
    sd_script_path = os.path.join(current_script_path, "..", "run_decomposer.py")
    log.log("Path to run_decomposer: " + sd_script_path)

    #create output_dir if not exists
    args.outdir = os.path.abspath(args.outdir)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    #switch working directory to output dir
    os.chdir(args.outdir)

    monomer_set_complete = False
    iter_id = 0
    while (not monomer_set_complete):
        log.log("Start iteration " + str(iter_id))
        #create for current iteration
        iter_outdir = os.path.join(args.outdir, "iter_" + str(iter_id))
        if not os.path.exists(iter_outdir):
            os.makedirs(iter_outdir)
        os.chdir(iter_outdir)

        # run string decomposer
        sys_call(["python3", sd_script_path, args.sequences, args.monomers])

        # parse output csv file
        monomer_set_complete = True
        iter_id += 1


if __name__ == "__main__":
    main()