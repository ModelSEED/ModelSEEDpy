# -*- coding: utf-8 -*-
import logging
from cobra.core.dictlist import DictList
from modelseedpy.core.msmodelutl import MSModelUtil

logger = logging.getLogger(__name__)


class MediaCompound:
    def __init__(self, compound_id, lower_bound, upper_bound, concentration=None, name=None):
        self.id = compound_id
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.concentration = concentration
        self.name = name

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
        self.media_ref = None  # Reference to the media object in the model
        self.mediacompounds = DictList()

    @staticmethod
    def from_dict(media_dict):
        """
        Create MSMedia from a dictionary in various formats.

        Supported formats:
        - Minimal: {'cpd00027': 10} - compound_id mapped to uptake rate (upper_bound)
        - Bounds: {'cpd00027': (-10, 1000)} - compound_id mapped to (lower_bound, upper_bound)
        - Complete: {'cpd00027': {'id': 'cpd00027', 'lower_bound': -10, 'upper_bound': 1000,
                                   'concentration': 5.0, 'name': 'Glucose'}}

        Parameters:
            media_dict (dict): Dictionary in one of the supported formats

        Returns:
            MSMedia: A new MSMedia instance
        """
        media = MSMedia("media")
        media_compounds = []

        for cpd_id, v in media_dict.items():
            if isinstance(v, dict):
                # Complete format - dictionary with all fields
                media_compounds.append(MediaCompound(
                    v.get('id', cpd_id),
                    v.get('lower_bound', -10),
                    v.get('upper_bound', 1000),
                    concentration=v.get('concentration'),
                    name=v.get('name')
                ))
            elif isinstance(v, tuple):
                # Bounds format - tuple of (lower_bound, upper_bound)
                media_compounds.append(MediaCompound(cpd_id, v[0], v[1]))
            else:
                # Minimal format - just the uptake value (upper_bound)
                media_compounds.append(MediaCompound(cpd_id, -v, 1000))

        media.mediacompounds += media_compounds
        return media

    @staticmethod
    def from_kbase_object(media_object):
        """
        Create MSMedia from KBase media object.
        :param media_object: KBase media object
        :return: MSMedia instance
        """
        media_id = media_object.id
        media_ref = None
        media_name = media_object.name if media_object.name else media_id
        if media_object.info is not None:
            media_id = media_object.info.id
            media_name = media_object.info.id
            media_ref = media_object.info.reference
        output = MSMedia(media_id, name=media_name)
        output.media_ref = media_ref
        media_compounds = []
        for mediacpd in media_object.mediacompounds:
            newmediacpd = MediaCompound(mediacpd.id, -1*mediacpd.maxFlux, -1*mediacpd.minFlux, concentration=mediacpd.concentration)
            media_compounds.append(newmediacpd)
        output.mediacompounds += media_compounds
        return output

    def to_dict(self, output_type="minimal"):
        """
        Convert MSMedia to a dictionary in various formats.

        Parameters:
            output_type (str): Type of output format. Options are:
                - "minimal": Dict of compound_id -> upper_bound (default)
                - "bounds": Dict of compound_id -> (lower_bound, upper_bound)
                - "complete": Dict of compound_id -> {all MediaCompound fields}

        Returns:
            dict: Dictionary representation of the media in the specified format

        Examples:
            Minimal: {'cpd00027': 1000}
            Bounds: {'cpd00027': (-10, 1000)}
            Complete: {'cpd00027': {'id': 'cpd00027', 'lower_bound': -10,
                                    'upper_bound': 1000, 'concentration': 5.0,
                                    'name': 'Glucose'}}
        """
        output = {}

        if output_type == "minimal":
            for compound in self.mediacompounds:
                output[compound.id] = compound.upper_bound

        elif output_type == "bounds":
            for compound in self.mediacompounds:
                output[compound.id] = (compound.lower_bound, compound.upper_bound)

        elif output_type == "complete":
            for compound in self.mediacompounds:
                output[compound.id] = {
                    'id': compound.id,
                    'lower_bound': compound.lower_bound,
                    'upper_bound': compound.upper_bound,
                    'concentration': compound.concentration,
                    'name': compound.name
                }

        else:
            raise ValueError(f"Invalid output_type '{output_type}'. Must be 'minimal', 'bounds', or 'complete'.")

        return output

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

    def add_compound(self, compound_id, lower_bound, upper_bound, concentration=None, name=None):
        """
        Add a compound to the media.

        Parameters:
            compound_id (str): The ID of the compound to add
            lower_bound (float): Lower bound for the compound flux
            upper_bound (float): Upper bound for the compound flux
            concentration (float, optional): Concentration of the compound. Defaults to None.
            name (str, optional): Name of the compound. Defaults to None.

        Returns:
            MediaCompound: The newly created MediaCompound object
        """
        new_compound = MediaCompound(
            compound_id,
            lower_bound,
            upper_bound,
            concentration=concentration,
            name=name
        )
        self.mediacompounds.append(new_compound)
        return new_compound

    def remove_compounds(self, compound_ids):
        """
        Remove compounds from the media by their IDs.

        Parameters:
            compound_ids (list): List of compound IDs to remove from the media

        Returns:
            list: List of removed MediaCompound objects
        """
        removed_compounds = []
        compounds_to_keep = []

        for compound in self.mediacompounds:
            if compound.id in compound_ids:
                removed_compounds.append(compound)
            else:
                compounds_to_keep.append(compound)

        # Replace the mediacompounds list with only the compounds to keep
        self.mediacompounds = DictList(compounds_to_keep)

        return removed_compounds

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

    def copy(self):
        """
        Create a deep copy of the MSMedia object.

        Returns:
            MSMedia: A new MSMedia instance with copied media compounds
        """
        new_media = MSMedia(self.id, name=self.name)
        new_media.media_ref = self.media_ref

        # Copy all media compounds
        copied_compounds = []
        for compound in self.mediacompounds:
            copied_compound = MediaCompound(
                compound.id,
                compound.lower_bound,
                compound.upper_bound,
                concentration=compound.concentration
            )
            copied_compounds.append(copied_compound)

        new_media.mediacompounds += copied_compounds
        return new_media
