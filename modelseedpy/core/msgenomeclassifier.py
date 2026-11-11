# -*- coding: utf-8 -*-
import numpy as np
from modelseedpy.core.msgenome import MSGenome


class MSLinearModel:
    """
    Linear classifier stored as plain arrays (no pickle, no scikit-learn needed to predict).

    predict(X) = classes[argmax(X @ coef.T + intercept)], the decision rule of
    sklearn LinearSVC / LogisticRegression for more than two classes.
    """

    def __init__(self, coef, intercept, classes):
        self.coef = np.asarray(coef, dtype=np.float64)
        self.intercept = np.asarray(intercept, dtype=np.float64)
        self.classes_ = np.asarray(classes)
        if self.coef.ndim != 2 or self.coef.shape[0] != len(self.classes_):
            raise ValueError(
                f"coef shape {self.coef.shape} does not match {len(self.classes_)} classes"
            )

    def decision_function(self, X):
        return np.asarray(X, dtype=np.float64) @ self.coef.T + self.intercept

    def predict(self, X):
        return self.classes_[np.argmax(self.decision_function(X), axis=1)]


class MSGenomeClassifier:
    def __init__(self, model, model_features):
        self.features = list(model_features)
        self.model = model
        self._feature_index = {f: i for i, f in enumerate(self.features)}

    @staticmethod
    def extract_features_from_genome(genome, ontology_term):
        """

        :param genome: ModelSEED Genome to classify
        :param ontology_term: Ontology Term to classify (Example: RAST)
        :return:
        """
        features = set()
        for feature in genome.features:
            if ontology_term in feature.ontology_terms:
                features.update(feature.ontology_terms[ontology_term])
        return {"genome": list(features)}

    def indicator_vector(self, roles):
        """0/1 row vector of the model features present in roles"""
        x = np.zeros((1, len(self.features)), dtype=np.int64)
        idx = [self._feature_index[r] for r in set(roles) if r in self._feature_index]
        x[0, idx] = 1
        return x

    def classify(self, genome_or_roles, ontology_term="RAST"):
        """
        :param genome_or_roles: MSGenome, or dict {genome_id: [roles]} (the first genome is classified)
        :param ontology_term: ontology of the genome features holding the roles (RAST)
        :return: predicted class label
        """
        if isinstance(genome_or_roles, MSGenome):
            genome_or_roles = self.extract_features_from_genome(
                genome_or_roles, ontology_term
            )
        roles = next(iter(genome_or_roles.values()))
        return self.model.predict(self.indicator_vector(roles))[0]


def load_classifier_from_npz(path):
    """
    Load a classifier saved as .npz with arrays: model_type ('linear'), coef, intercept,
    classes and features (role order of the coef columns).
    """
    with np.load(path, allow_pickle=False) as data:
        model_type = str(data["model_type"])
        if model_type != "linear":
            raise ValueError(f"unsupported classifier model_type: {model_type}")
        model = MSLinearModel(
            data["coef"], data["intercept"], data["classes"].astype(str)
        )
        features = data["features"].astype(str).tolist()
    return MSGenomeClassifier(model, features)


def load_classifier_from_folder(path, filename):
    """
    TEMPORARY SOLUTION TO LOAD AN EXISTING CLASSIFIER
    :param path:
    :param filename: {filename}.npz, or {filename}.pickle + {filename}_features.json
    :return:
    """
    import os
    import pickle
    import json

    if os.path.exists(f"{path}/{filename}.npz"):
        return load_classifier_from_npz(f"{path}/{filename}.npz")
    with open(f"{path}/{filename}.pickle", "rb") as fh:
        model_filter = pickle.load(fh)
    with open(f"{path}/{filename}_features.json", "r") as fh:
        features = json.load(fh)

    return MSGenomeClassifier(model_filter, features)
