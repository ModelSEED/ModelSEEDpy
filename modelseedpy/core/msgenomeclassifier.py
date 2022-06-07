from modelseedpy.ml.predict_phenotype import create_indicator_matrix


class MSGenomeClassifier:
    def __init__(self, model, model_features):
        self.features = model_features
        self.model = model

    @staticmethod
    def extract_features_from_genome(genome, ontology_term):
        """
        :param genome: ModelSEED Genome to classify
        :param ontology_term: Ontology Term to classify (Example: RAST)
        :return:
        """
        features = set()
        for feat in genome.features:
            if ontology_term in feat.ontology_terms:
                features.update(feat.ontology_terms[ontology_term])
        return {'genome': list(features)}

    def classify(self, genome, ontology_term='RAST'):
        roles = self.extract_features_from_genome(genome, ontology_term)
        indicator_df, master_role_list = create_indicator_matrix(roles, self.features)
        predictions_numerical = self.model.predict(indicator_df[master_role_list].values)
        return predictions_numerical[0]


def load_classifier_from_folder(directory, filename):
    """
    TEMPORARY SOLUTION TO LOAD AN EXISTING CLASSIFIER
    :param directory:
    :param filename:
    :return:
    """
    import pickle
    import json
    with open(f'{directory}/{filename}.pickle', 'rb') as fh:
        model_filter = pickle.load(fh)
    with open(f'{directory}/{filename}_features.json', 'r') as fh:
        features = json.load(fh)

    return MSGenomeClassifier(model_filter, features)
