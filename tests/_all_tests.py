# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

import logging

import os
import sys
import unittest

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))


from sd.standard_logger import get_logger
from sd.utils.git import get_git_revision_short_hash


def main():
    logger = get_logger('test.log',
                        logger_name='SD',
                        level=logging.DEBUG,
                        filemode='w',
                        stdout=False)
    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    unittest.main()


if __name__ == "__main__":
   main()
