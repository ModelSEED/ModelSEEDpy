# -*- coding: utf-8 -*-
import logging
from cobra.core.dictlist import DictList
from modelseedpy.core.msmodelutl import MSModelUtil

logger = logging.getLogger(__name__)


class MediaCompound:
    def __init__(self, compound_id, lower_bound, upper_bound, concentration=None):
        self.id = compound_id
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.concentration = concentration

    @property
    def maxFlux(self):
        # TODO: will be removed later just for old methods
        return -self.lower_bound

    @property
    def minFlux(self):
        # TODO: will be removed later just for old methods
        return -self.upper_bound

    def get_mdl_exchange_hash(self, model_or_mdlutl):
        modelutl = model_or_mdlutl
        if not isinstance(model_or_mdlutl, MSModelUtil):
            modelutl = MSModelUtil.get(model_or_mdlutl)
        mets = modelutl.find_met(self.id)
        output = {}
        exchange_hash = modelutl.exchange_hash()
        for met in mets:
            if met in exchange_hash:
                output[met] = exchange_hash[met]
        return output


class MSMedia:
    def __init__(self, media_id, name=""):
        self.id = media_id
        self.name = name
        self.mediacompounds = DictList()

    @staticmethod
    def from_dict(media_dict):
        """
        Either dict with exchange bounds (example: {'cpd00027': (-10, 1000)}) or
        just absolute value of uptake (example: {''cpd00027': 10})
        :param media_dict:
        :return:
        """
        media = MSMedia("media")
        media_compounds = []
        for cpd_id, v in media_dict.items():
            if isinstance(v, tuple):
                media_compounds.append(MediaCompound(cpd_id, v[0], v[1]))
            else:
                media_compounds.append(MediaCompound(cpd_id, -v, 1000))
        media.mediacompounds += media_compounds
        return media

    def get_media_constraints(self, cmp="e0"):
        """
        Parameters:
            cmp (str): compound suffix (model compartment)
        Returns:
            dict(str) -> (float,float): compound_ids mapped to lower/upper bound
        """
        media = {}
        for compound in self.mediacompounds:
            met_id = compound.id
            if cmp is not None:
                met_id += "_" + cmp
            media[met_id] = (compound.lower_bound, compound.upper_bound)
        return media

    def find_mediacpd(self, cpd_id):
        for cpd in self.mediacompounds:
            if cpd.id == cpd_id:
                return cpd
        return None

    def merge(self, media, overwrite_overlap=False):
        new_cpds = []
        for cpd in media.mediacompounds:
            newcpd = MediaCompound(
                cpd.id, -cpd.maxFlux, -cpd.minFlux, cpd.concentration
            )
            if newcpd.id in self.mediacompounds:
                if overwrite_overlap:
                    self.mediacompounds[newcpd.id] = newcpd
            else:
                new_cpds.append(newcpd)
        self.mediacompounds += new_cpds
