from cobra.core.dictlist import DictList


class MediaCompound:

    def __init__(self, compound_id, lower_bound, upper_bound, concentration=None):
        self.id = compound_id
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.concentration = concentration

    @property
    def maxFlux(self):
        # TODO: will be removed later just for old methods
        return -1 * self.lower_bound

    @property
    def minFlux(self):
        # TODO: will be removed later just for old methods
        return -1 * self.upper_bound


class MSMedia:

    def __init__(self, media_id):
        self.id = media_id
        self.mediacompounds = DictList()

    @staticmethod
    def from_dict(d):
        """
        Either dict with exchange bounds (example: {'cpd00027': (-10, 1000)}) or
        just absolute value of update (example: {''cpd00027': 10})
        :param d:
        :return:
        """
        media = MSMedia('media')
        media_compounds = []
        for cpd_id in d:
            v = d[cpd_id]
            if type(v) is tuple:
                media_compounds.append(MediaCompound(cpd_id, v[0], v[1]))
            else:
                media_compounds.append(MediaCompound(cpd_id, -1 * v, 1000))

        media.mediacompounds += media_compounds
        return media

    def get_media_constraints(self, cmp='e0'):
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
                met_id += '_' + cmp
            media[met_id] = (compound.lower_bound, compound.upper_bound)

        return media
