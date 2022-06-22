
#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation"""
    def __init__(self, message):
        super(FeasibilityError, self).__init__(message)