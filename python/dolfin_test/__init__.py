# -*- coding: utf-8 -*-
"""Main module for DOLFIN"""


import sys

# Store dl open flags to restore them after import
stored_dlopen_flags = sys.getdlopenflags()


# Fix dlopen flags (may need reorganising)
import sys
if "linux" in sys.platform:
    # FIXME: What with other platforms?
    try:
        from ctypes import RTLD_NOW, RTLD_GLOBAL
    except ImportError:
        RTLD_NOW = 2
        RTLD_GLOBAL = 256
    sys.setdlopenflags(RTLD_NOW | RTLD_GLOBAL)
del sys

#import dolfin_test.cpp

# Reset dl open flags
#import sys
#sys.setdlopenflags(stored_dlopen_flags)
#del sys
