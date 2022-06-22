# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class RevBinPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"reversible binary",{"revbin":"reaction"},{"revbinF":"reaction","revbinR":"reaction"})

    def build_package(self,filter = None):
        for reaction in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if filter == None or reaction.id in filter:
                self.build_variable(reaction)
                self.build_constraint(reaction)
    
    def build_variable(self,object):
        return BaseFBAPkg.build_variable(self,"revbin",0,1,"binary",object)
    
    def build_constraint(self,object):
        #-1000 * revbin(i) + forv(i) <= 0
        BaseFBAPkg.build_constraint(self,"revbinF",None,0,{self.variables["revbin"][object.id]:-1000,object.forward_variable:1},object)
        #1000 * revbin(i) + revv(i) <= 1000
        return BaseFBAPkg.build_constraint(self,"revbinR",None,1000,{self.variables["revbin"][object.id]:1000,object.reverse_variable:1},object)
