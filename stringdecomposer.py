#!/usr/bin/env python3

import os
import sys

sd_root = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, sd_root)

from stringdecomposer.main import main
sys.exit(main())
