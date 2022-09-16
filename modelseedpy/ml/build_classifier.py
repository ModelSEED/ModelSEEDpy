# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from sklearn.utils import class_weight
from sklearn.metrics import confusion_matrix


def unload_training_set(training_set_object):
    """
    Take a Training Set Object and extracts "Genome Name", "Genome Reference", "Phenotype",
    and "Phenotype Enumeration" and places into a DataFrame

    class_enumeration will be something like: {'N': 0, 'P': 1} ie. phenotype --> number

    Parameter
    ---------
    current_ws : str
        current_ws
    training_set_name : str
        training set to use for app
    """

    phenotype = training_set_object[0]["data"]["classification_type"]
    classes_sorted = training_set_object[0]["data"]["classes"]

    class_enumeration = {}  # {'N': 0, 'P': 1}
    for index, _class in enumerate(classes_sorted):
        class_enumeration[_class] = index

    training_set_object_data = training_set_object[0]["data"]["classification_data"]
    training_set_object_reference = training_set_object[0]["path"][0]

    _names = []
    _references = []
    _phenotypes = []
    _phenotype_enumeration = []

    for genome in training_set_object_data:
        _names.append(genome["genome_name"])
        _references.append(genome["genome_ref"])
        _phenotypes.append(genome["genome_classification"])

        _enumeration = class_enumeration[genome["genome_classification"]]
        _phenotype_enumeration.append(_enumeration)

    uploaded_df = pd.DataFrame(
        data={
            "Genome Name": _names,
            "Genome Reference": _references,
            "Phenotype": _phenotypes,
            "Phenotype Enumeration": _phenotype_enumeration,
        }
    )

    return phenotype, class_enumeration, uploaded_df, training_set_object_reference


def save_classifier(scratch, folder_name, file_name, dfu_utils, classifier, x, y):
    import os
    import pickle

    with open(os.path.join(scratch, folder_name, "data", file_name), "wb") as fh:
        main_clf = classifier.fit(x, y)
        pickle.dump(main_clf, fh, protocol=2)

    shock_id, handle_id = dfu_utils._upload_to_shock(
        os.path.join(scratch, folder_name, "data", file_name)
    )
    return shock_id, handle_id


def execute_classifier(
    current_ws, common_classifier_information, current_classifier_object, folder_name
):
    """
    Creates k=splits number of classifiers and then generates a confusion matrix that averages
    over the predicted results for all of the classifiers.

    Generates statistics for each classifier (saved in classification_report_dict)

    Saves each classifier object as a pickle file and then uploads that object to shock, saves the object's
    shock handle into the KBASE Categorizer Object (https://narrative.kbase.us/#spec/type/KBaseClassifier.GenomeCategorizer)

    Calls function to create png of confusion matrix

    Saves information for callback of build_classifier in individual_classifier_info

    Parameter
    ---------
    current_ws : str
        current_ws
    common_classifier_information : dict
        information that is common to the current classifier that is going to be built
    current_classifier_object: dict
        information that is specific to the current classifier that is going to be built
    folder_name:
        folder name gives the location to save classifier and images
    """

    individual_classifier_info = {}
    matrix_size = len(common_classifier_information["class_list_mapping"])
    cnf_matrix_proportion = np.zeros(shape=(matrix_size, matrix_size))

    classifier = current_classifier_object["classifier_to_execute"]

    for c in range(common_classifier_information["splits"]):
        X_train = common_classifier_information["whole_X"][
            common_classifier_information["list_train_index"][c]
        ]
        y_train = common_classifier_information["whole_Y"][
            common_classifier_information["list_train_index"][c]
        ]
        X_test = common_classifier_information["whole_X"][
            common_classifier_information["list_test_index"][c]
        ]
        y_test = common_classifier_information["whole_Y"][
            common_classifier_information["list_test_index"][c]
        ]

        # do class reweighting specifically for GaussianNB¶
        if current_classifier_object["classifier_type"] == "gaussian_nb":
            # https://datascience.stackexchange.com/questions/13490/how-to-set-class-weights-for-imbalanced-classes-in-keras
            unique_classes = np.unique(y_train)
            class_weights = class_weight.compute_class_weight(
                "balanced", unique_classes, y_train
            )

            dict_class_to_weight = {
                curr_class: curr_weight
                for curr_class, curr_weight in zip(unique_classes, class_weights)
            }
            sample_weight = [dict_class_to_weight[curr_class] for curr_class in y_train]

            classifier.fit(X_train, y_train, sample_weight=sample_weight)
        else:
            classifier.fit(X_train, y_train)

        y_pred = classifier.predict(X_test)

        cnf = confusion_matrix(
            y_test,
            y_pred,
            labels=list(common_classifier_information["class_list_mapping"].values()),
        )
        cnf_f = cnf.astype("float") / cnf.sum(axis=1)[:, np.newaxis]
        for i in range(len(cnf)):
            for j in range(len(cnf)):
                cnf_matrix_proportion[i][j] += cnf_f[i][j]

    # get statistics for the last case made
    # diagonal entries of cm are the accuracies of each class
    target_names = list(common_classifier_information["class_list_mapping"].keys())
    classification_report_dict = classification_report(
        y_test, y_pred, target_names=target_names, output_dict=True
    )

    # save down classifier object in pickle format
    # no more saving!
    """
    shock_id, handle_id = save_classifier(self.scratch, folder_name,
                                          current_classifier_object["classifier_name"] + ".pickle",
                                          classifier,
                                          common_classifier_information["whole_X"],
                                          common_classifier_information["whole_Y"])
    """

    classifier_object = {
        "classifier_id": "",
        "classifier_type": current_classifier_object["classifier_type"],
        "classifier_name": current_classifier_object["classifier_name"],
        "classifier_data": "",  # saved in shock
        # 'classifier_handle_ref': handle_id,
        "classifier_description": common_classifier_information["description"],
        "lib_name": "sklearn",
        "attribute_type": common_classifier_information["attribute_type"],
        "number_of_attributes": len(
            common_classifier_information["attribute_data"]
        ),  # size of master_role_list
        "attribute_data": common_classifier_information["attribute_data"],
        "class_list_mapping": common_classifier_information["class_list_mapping"],
        "number_of_genomes": len(common_classifier_information["whole_Y"]),
        "training_set_ref": common_classifier_information["training_set_ref"],
    }

    """
    obj_save_ref = self.ws_client.save_objects({'workspace': current_ws,
                                                'objects': [{
                                                    'type': 'KBaseClassifier.GenomeCategorizer',
                                                    'data': classifier_object,
                                                    'name': current_classifier_object["classifier_name"],
                                                    'provenance': self.ctx['provenance']
                                                }]
                                                })[0]
    """

    # information for call back
    """
    individual_classifier_info = {"classifier_name": current_classifier_object["classifier_name"],
                                  "classifier_ref": obj_save_ref,
                                  "accuracy": classification_report_dict["accuracy"]}
    """
    cm = np.round(
        cnf_matrix_proportion / common_classifier_information["splits"] * 100.0, 1
    )
    title = "CM: " + current_classifier_object["classifier_type"]
    # self.plot_confusion_matrix(cm, title, current_classifier_object["classifier_name"],
    #                           list(common_classifier_information["class_list_mapping"].keys()), folder_name)

    return classification_report_dict, individual_classifier_info, classifier
