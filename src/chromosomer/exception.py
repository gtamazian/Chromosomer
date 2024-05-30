#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


class Error(Exception):
    """
    The class describes a basic error that may occur in any of the
    Chromosomer-related routines.
    """

    pass


class MapError(Error):
    """
    The class describes an error that may occur while working with a
    fragment __map object.
    """

    pass


class AlignmentToMapError(Error):
    """
    The class describes an error that may occur while creating a
    fragment __map from alignments.
    """

    pass
