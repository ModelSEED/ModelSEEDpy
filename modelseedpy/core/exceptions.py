# -*- coding: utf-8 -*-

# Adding a few exception classes to handle different types of errors in a central file
class ModelSEEDError(Exception):
    """Error in ModelSEED execution logic"""

    pass


class FeasibilityError(Exception):
    """Error in FBA formulation"""

    def __init__(self, message):
        super(FeasibilityError, self).__init__(message)


class PackageError(Exception):
    """Error in package manager"""

    pass


class GapfillingError(Exception):
    """Error in model gapfilling"""
    pass

class ObjectError(Exception):
    """Error in the construction of a base KBase object"""
    pass 

class ParameterError(Exception):
    """Error in a parameterization"""
    pass 

class ObjectAlreadyDefinedError(Exception):
    pass

class NoFluxError(Exception):
    """Error for FBA solutions"""
    pass

class ObjectiveError(Exception):
    """Erroneous assignment of a secondary objective via a constraint"""
    pass

class ModelError(Exception):
    """Errors in a model that corrupt the simulation"""
    pass
