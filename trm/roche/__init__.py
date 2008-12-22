"""
Roche geometry routines

This module provides some easy-to-use routines for computations
to do with Roche geometry. All routines work in coordinates scaled
by the binary separation, and time units such that the binary angular
frequency = 1.

Functions
=========

ingress_egress -- calculate the ingress/egress phase of a given point.
lobe1          -- the primary star's Roche lobe
lobe2          -- the secondary star's Roche lobe
stream         -- the gas stream
strmnx         -- position and velocity of n-th turning point of gas stream
vlobe1         -- the primary star's Roche lobe, velocity space
vlobe2         -- the secondary star's Roche lobe, velocity space
vstream        -- gas stream in velocity coordinates
xl1            -- L1 position
xl2            -- L2 position
xl3            -- L3 position
"""

import sys
sys.path.append('.')
from _roche import *
