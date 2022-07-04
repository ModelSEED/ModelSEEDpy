
#Adding a few exception classes to handle different types of errors in a central file
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