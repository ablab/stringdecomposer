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

try:
    sys.stderr = sys.stdout
    exit_code = 0
    args = parse_args()
    dataset_info = load_info(args.info)
    log.log(str(dataset_info))
    working_dir = os.getcwd()

except:
    log.err("The following unexpected error occured during the run:")
    print_exc()
    sys.exit(1)