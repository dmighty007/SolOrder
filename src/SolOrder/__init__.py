"""Structural order parameters for water simulations."""

from importlib.metadata import version

from .dihedral import OrderParameter
from .entropy import eFingerprint
from .structure import StructuralOrder, HydrogenBonds
from .extended import AngularDistribution, VoronoiGeometry, IceClassification, RingStatistics
from .descriptors import GraphDescriptor, TemplateMatch

__version__ = version("SolOrder")
__all__ = [
    "OrderParameter", "eFingerprint", "StructuralOrder", "HydrogenBonds",
    "AngularDistribution", "VoronoiGeometry", "IceClassification", "RingStatistics",
    "GraphDescriptor", "TemplateMatch",
]
