# -*- coding: utf-8 -*-
import logging
import os
import requests
import pickle
import json
from configparser import ConfigParser


logger = logging.getLogger(__name__)
project_dir = os.path.abspath(os.path.dirname(__file__)) + "/"
config = ConfigParser()
config.read(project_dir + "/config.cfg")


def get_or_download_file(filename, k, value, config):
    folder_path = f"{project_dir}/" + config.get(k, value)
    file_path = f"{folder_path}/{filename}"
    if not os.path.exists(folder_path):
        logger.warning("mkdir: %s", folder_path)
        os.makedirs(folder_path)
    if not os.path.exists(file_path):
        logger.warning("downloading data file to: %s", file_path)
        url = "https://bioseed.mcs.anl.gov/~fliu/modelseedpy/" + filename
        r = requests.get(url, allow_redirects=True)
        with open(file_path, "wb") as fh:
            fh.write(r.content)
    return file_path


def get_file(filename, k, value):
    return get_or_download_file(filename, k, value, config)


def get_classifier(classifier_id):
    """
    Load a genome classifier from the classifier folder (data/ml).

    Classifiers shipped with the package are {classifier_id}.npz (numpy arrays, no pickle).
    Legacy classifiers are {classifier_id}.pickle + {classifier_id}_features.json, downloaded
    on first use; those pickles need the scikit-learn version that wrote them
    (scikit-learn 1.2.0 for the *_01_17_2023 models).
    """

    npz = f"{project_dir}/{config.get('data', 'classifier_folder')}/{classifier_id}.npz"
    if os.path.exists(npz):
        from modelseedpy.core.msgenomeclassifier import load_classifier_from_npz

        return load_classifier_from_npz(npz)

    logger.warning(
        f"{npz} not found. fallback to old classifier. (requires scikit-learn 1.2.0)"
    )

    cls_pickle = get_file(f"{classifier_id}.pickle", "data", "classifier_folder")
    cls_features = get_file(
        f"{classifier_id}_features.json", "data", "classifier_folder"
    )
    with open(cls_pickle, "rb") as fh:
        model_filter = pickle.load(fh)
    with open(cls_features, "r") as fh:
        features = json.load(fh)

    from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier

    return MSGenomeClassifier(model_filter, features)


def get_template(template_id):
    # we need a mstemplate object!
    template_file = get_file(f"{template_id}.json", "data", "template_folder")
    with open(template_file, "r") as fh:
        return json.load(fh)
