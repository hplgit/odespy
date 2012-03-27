#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
This file contains the constants used in Assimulo.
"""

#Return flags
"""
DEF ID_OK       = 0
DEF ID_DISCARD  = 1
DEF ID_EVENT    = 2
DEF ID_COMPLETE = 3
DEF ID_FAIL     = -1
"""

#Verbosity levels
DEF DEBUG    = 10
DEF INFO     = 20
DEF WARNING  = 30
DEF ERROR    = 40
DEF CRITICAL = 50

#Backward compability
QUIET    = CRITICAL
WHISPER  = ERROR
NORMAL   = WARNING
LOUD     = INFO
SCREAM   = DEBUG


DEF ID_OK       = 0
DEF ID_DISCARD  = 1
DEF ID_EVENT    = 2
DEF ID_COMPLETE = 3
DEF ID_FAIL     = -1

ID_PY_OK       = 0
ID_PY_DISCARD  = 1
ID_PY_EVENT    = 2
ID_PY_COMPLETE = 3
ID_PY_FAIL     = -1

