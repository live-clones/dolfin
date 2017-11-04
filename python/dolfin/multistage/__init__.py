# -*- coding: utf-8 -*-
"""The multistage module of dolfin"""

from dolfin.multistage import multistagescheme
from dolfin.multistage import multistagesolvers
from dolfin.multistage import rushlarsenschemes

from dolfin.multistage.multistagescheme import *  # noqa
from dolfin.multistage.multistagesolvers import *  # noqa
from dolfin.multistage.rushlarsenschemes import *  # noqa

# NOTE: The automatic documentation system in DOLFIN requires to _not_
# define classes or functions within this file. Use separate modules
# for that purpose.

__all__ = multistagescheme.__all__ + multistagesolvers.__all__ + \
    rushlarsenschemes.__all__
