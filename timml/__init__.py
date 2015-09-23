"""
Copyright (C), 2015, Mark Bakker.
Mark Bakker, Delft University of Technology
mark dot bakker at tudelft dot nl

TimML is a computer program for the simulation of steady-state
multiaquifer flow with analytic elements and consists of a
library of Python scripts and FORTRAN extensions.
"""
# Import all classes and functions
from ml import Model, Model3D
from mlaquifer import Aquifer, PolygonInhom, CircleInhomData, EllipseInhomData
from mlconstant import Constant
from mluflow import Uflow, Triple
from mlwell import Well
from mllinesink import LineSink, HeadLineSink, ResLineSink, LineSinkDitch
from mlcircareasink import CircAreaSink
from mllinedoubletimp import LineDoubletImp
from mlpolyareasink import PolyAreaSink
from mlinhom import MakeInhomPolySide, MakeInhomogeneity
from mlcircinhom import CircleInhom
from mltrace import traceline
from mlutil import timcontour, timlayout, timtracelines, timvertcontour, \
    capturezone

# Import version number
from version import __version__


__author__ = 'Mark Bakker'
