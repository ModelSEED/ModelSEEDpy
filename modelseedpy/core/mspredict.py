# -*- coding: utf-8 -*-
import logging
from typing import Tuple
from enum import Enum
from modelseedpy.helpers import get_template, get_classifier
from modelseedpy.core.mstemplate import MSTemplateBuilder, MSTemplate
from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)


class MSGenomeClass(Enum):
    P = "Gram Positive"
    N = "Gram Negative"
    C = "Cyano"
    A = "Archaea"


DEFAULT_CLASSIFIER = "linear_svc_ACNP_RAST_full_11_11_2026"


class MSPredict:
    def __init__(self, classifier=DEFAULT_CLASSIFIER):
        """
        :param classifier: classifier id in the data/ml folder. The default is a linear SVM on
            all RAST roles (33,382), stored as numpy arrays; scikit-learn is not needed.
            The previous default, "knn_ACNP_RAST_filter_01_17_2023", is still loadable
            (downloaded on first use; its pickle requires scikit-learn 1.2.0).
        """
        logger.debug(f"Initializing MSPredict with classifier: {classifier}")
        self.genome_classifier = get_classifier(classifier)

    def predict(self, genome: MSGenome) -> MSGenomeClass:
        """
        :param genome: MSGenome - The genome to predict the class for.
        :return: MSGenomeClass - The genome class for the given genome.
        Predicts the genome class for a given genome.
        """
        logger.debug(f"Predicting genome class:")

        genome_class_str = self.genome_classifier.classify(genome)

        # Map the classifier result to the appropriate enum
        if genome_class_str == "P":
            return MSGenomeClass.P
        elif genome_class_str == "N":
            return MSGenomeClass.N
        elif genome_class_str == "C":
            return MSGenomeClass.C
        elif genome_class_str == "A":
            return MSGenomeClass.A
        else:
            raise ValueError(f"Unknown genome class: {genome_class_str}")

    def auto_select_template(
        self, genome_class: MSGenomeClass
    ) -> Tuple[MSTemplate, MSTemplate]:
        """
        :param genome_class: MSGenomeClass - The genome class to select the template for.
        :return: Tuple[MSTemplate, MSTemplate] - The template for the genome class.
        Auto-selects the template for a given genome class.
        """
        template_genome_scale_map = {
            MSGenomeClass.N: "template_gram_neg",
            MSGenomeClass.C: "template_gram_neg",
            MSGenomeClass.A: "template_gram_neg",
            MSGenomeClass.P: "template_gram_pos",
        }
        template_core_map = {
            MSGenomeClass.A: "template_core",
            MSGenomeClass.C: "template_core",
            MSGenomeClass.N: "template_core",
            MSGenomeClass.P: "template_core",
        }

        if (
            genome_class in template_genome_scale_map
            and genome_class in template_core_map
        ):
            d_template_core = get_template(template_core_map[genome_class])
            d_template_genome_scale = get_template(
                template_genome_scale_map[genome_class]
            )
            template_core = MSTemplateBuilder.from_dict(d_template_core).build()
            template_genome_scale = MSTemplateBuilder.from_dict(
                d_template_genome_scale
            ).build()
            return template_core, template_genome_scale
        else:
            raise ValueError(f"Unknown genome class: {genome_class}")
