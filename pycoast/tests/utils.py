#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2022 PyCoast Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Utilities for working with tests."""

import os
from contextlib import contextmanager
from pathlib import Path
from typing import Union


@contextmanager
def set_directory(path: Union[str, Path]):
    """Change to new directory within the context.

    Args:
        path: The path to the new cwd

    Yields:
        None
    """
    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)
