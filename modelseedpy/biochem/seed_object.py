class ModelSEEDObject:

    def __init__(self, data, api=None):
        self.data = data
        self.api = api

    @property
    def id(self):
        return self.data['id']

    @property
    def name(self):
        return self.data['name']
