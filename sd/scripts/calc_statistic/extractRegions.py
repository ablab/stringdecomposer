r'''
Extract blocks for specific monomers
'''

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--SD_tsv", dest="SD_tsv")
    parser.add_argument("--seq", dest="seq")

    return parser.parse_args()


def main():
    args = parse_args()


if __name__ =="__main__":
    main()