#!/usr/bin/python

import sys
import os
import yaml
import shutil
import argparse
import subprocess
from traceback import print_exc

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)
        sys.stdout.flush()

    def warn(self, s):
        msg = "WARNING: " + s + "\n"
        self.text += msg
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

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('info', metavar='CONFIG_FILE', type=str,  help='a path to .yaml file with test description')
    args = parser.parse_args()
    return args


def load_info(info_filename):
    info = yaml.load(open(info_filename, 'r'))
    return info


def compile():
    log.log("Current working dir: " + str(os.getcwd()))

    os.chdir(os.path.join(os.getcwd(), ".."))
    err_code = os.system('make')
    os.chdir(os.path.join(os.getcwd(), "sd"))
    return err_code


def run_script(dataset_info):
    log.log("Running script " + dataset_info["script"] + " with args: " + str(dataset_info["args"]))
    err_code = subprocess.run([sys.executable, dataset_info["script"]] + dataset_info["args"]).returncode
    return err_code

try:
    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args()
    dataset_info = load_info(args.info)
    log.log(str(dataset_info))

    # compile
    ecode = compile()
    if ecode != 0:
        log.err("Compilation finished abnormally with exit code " + str(ecode))
        sys.exit(3)

    ecode = run_script(dataset_info)
    if ecode != 0:
        log.err("Running script finished abnormally with exit code " + str(ecode))
        sys.exit(2)

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(1)