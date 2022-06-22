# -*- coding: utf-8 -*-

from __future__ import absolute_import
from modelseedpy.core.exceptions import PackageError

import sys
import inspect

class MSPackageManager:
    """
    Class for organizing FBA package objects
    """
    pkgmgrs = {}
    
    @staticmethod
    def get_pkg_mgr(model,create_if_missing = True):
        if model in MSPackageManager.pkgmgrs:
            return  MSPackageManager.pkgmgrs[model]
        elif create_if_missing:
            MSPackageManager.pkgmgrs[model] = MSPackageManager(model)
            return MSPackageManager.pkgmgrs[model]
        else:
            return None
    
    def __init__(self, model):
        self.model = model
        self.packages = {}
        self.available_packages = {}
        for name, obj in inspect.getmembers(sys.modules["modelseedpy"]):
            if name != "BaseFBAPkg" and name[-3:] == "Pkg":
                self.available_packages[name] = obj
        
    def list_available_packages(self):
        return list(self.available_packages.keys())
    
    def list_active_packages(self):
        return list(self.packages.keys())
    
    def addpkgobj(self,object):
        classname = type(object).__name__
        if classname not in self.packages:
            self.packages[classname] = object
        elif self.packages[classname] != object:
            raise PackageError("Package with name "+classname+" already in model!")
        return self.packages[classname]
    
    def addpkgs(self,packages):
        for package in packages:
            if package not in self.available_packages:
                raise PackageError("Package "+package+" does not exist!")
            if package not in self.packages:
                self.packages[package] = self.available_packages[package](self.model)
    
    def getpkg(self,package,create_if_missing = True):
        if package not in self.packages:
            if create_if_missing:
                self.addpkgs([package])
            else:
                return None
        return self.packages[package]