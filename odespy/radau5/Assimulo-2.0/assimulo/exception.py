
class AssimuloException(Exception):
    pass 

class TerminateSimulation(AssimuloException):
    pass

class DiscardValue(AssimuloException):
    pass

class Explicit_ODE_Exception(AssimuloException):
    pass
    
class ODE_Exception(AssimuloException):
    pass

class Implicit_ODE_Exception(AssimuloException):
    pass    
    
