# (c) 2020 by Authors
# This file is a part of SD program.
# see LICENSE file

from bisect import bisect_left
from itertools import islice

import numpy as np


def list2str(lst, sep=' '):
    return sep.join(str(e) for e in lst)


def listEls2str(lst):
    return [str(e) for e in lst]


def fst_iterable(iterable):
    return next(iter(iterable))
