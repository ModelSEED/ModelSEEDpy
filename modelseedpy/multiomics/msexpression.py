# -*- coding: utf-8 -*-
import logging
from typing import Optional, Union, TYPE_CHECKING

#from numpy._core.numeric import True_
import pandas as pd
import numpy as np
import re
import copy
from cobra.core.dictlist import DictList
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr, GPR
from cobra import Solution
from ast import And, BitAnd, BitOr, BoolOp, Expression, Name, NodeTransformer, Or

from sympy.logic import false
from modelseedpy.core.msgenome import MSGenome, MSFeature
from modelseedpy.core.msmodelutl import MSModelUtil

logger = logging.getLogger(__name__)

def compute_gene_score(expr, values, default, datatype):
    # Handle tuple return from parse_gpr() in newer COBRApy versions
    if isinstance(expr, tuple) and len(expr) == 2:
        expr = expr[0]  # Extract GPR object from (GPR, frozenset) tuple

    if isinstance(expr, (Expression, GPR)):
        return compute_gene_score(expr.body, values, default, datatype)
    elif isinstance(expr, Name):
        if expr.id in values:
            return values[expr.id]
        else:
            return default
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            best = None
            total = None
            for subexpr in expr.values:
                value = compute_gene_score(subexpr, values, default, datatype)
                if value != None:
                    if datatype == "NormalizedRatios":
                        diff = abs(value - 1)
                        if best == None or diff > best:
                            best = diff
                            total = value
                    elif datatype == "RelativeAbundance" or datatype == "FPKM" or datatype == "TPM" or datatype == "AbsoluteAbundance":
                        if total == None:
                            total = 0
                        total += value
            return total
        elif isinstance(op, And):
            best = None
            best_value = None
            for subexpr in expr.values:
                value = compute_gene_score(subexpr, values, default, datatype)
                if value != None:
                    if datatype == "NormalizedRatios":
                        diff = abs(value - 1)
                        if best == None or diff > best:
                            best = diff
                            best_value = value
                    elif datatype == "RelativeAbundance" or datatype == "FPKM" or datatype == "TPM" or datatype == "AbsoluteAbundance":
                        if best == None or value < best:
                            best = value
                            best_value = value
            return best_value
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    elif expr is None:
        return default
    else:
        raise TypeError("unsupported operation  " + repr(expr))

class MSCondition:
    def __init__(self, id, parent):
        self.id = id
        self.parent = parent

    def value_at_zscore(self, zscore: float) -> Optional[float]:
        """Calculate the value at a given z-score for this condition.

        Args:
            zscore: The z-score threshold

        Returns:
            The value at the specified z-score, or None if no data
        """
        if self.id not in self.parent._data.columns:
            return None

        values = self.parent._data[self.id].dropna()
        if len(values) == 0:
            return None

        mean = values.mean()
        std_dev = values.std()
        return mean + (zscore * std_dev)

    def lowest_value(self) -> Optional[float]:
        """Get the minimum value for this condition.

        Returns:
            The minimum value, or None if no data
        """
        if self.id not in self.parent._data.columns:
            return None

        values = self.parent._data[self.id].dropna()
        if len(values) == 0:
            return None

        return float(values.min())

    def highest_value(self) -> Optional[float]:
        """Get the maximum value for this condition.

        Returns:
            The maximum value, or None if no data
        """
        if self.id not in self.parent._data.columns:
            return None

        values = self.parent._data[self.id].dropna()
        if len(values) == 0:
            return None

        return float(values.max())

    def average_value(self) -> Optional[float]:
        """Get the mean value for this condition.

        Returns:
            The mean value, or None if no data
        """
        if self.id not in self.parent._data.columns:
            return None

        values = self.parent._data[self.id].dropna()
        if len(values) == 0:
            return None

        return float(values.mean())

    def sum_value(self) -> float:
        """Get the sum of all values for this condition.

        Returns:
            The sum of all values (0.0 if no data)
        """
        if self.id not in self.parent._data.columns:
            return 0.0

        return float(self.parent._data[self.id].sum(skipna=True))

class MSExpressionFeature:
    def __init__(self, feature, parent):
        self.id = feature.id
        self.feature = feature
        self.parent = parent

    def add_value(self, condition: Union[str, 'MSCondition'], value: float) -> None:
        """Add a value for a specific condition.

        Args:
            condition: MSCondition object or condition ID string
            value: The expression value to store
        """
        # Resolve condition to condition_id
        if isinstance(condition, str):
            condition_id = condition
        else:
            condition_id = condition.id

        # Ensure feature row exists in parent DataFrame
        if self.id not in self.parent._data.index:
            self.parent._data.loc[self.id, :] = np.nan

        # Ensure condition column exists in parent DataFrame
        if condition_id not in self.parent._data.columns:
            self.parent._data[condition_id] = np.nan

        # Set the value
        self.parent._data.loc[self.id, condition_id] = value

    def get_value(self, condition: Union[str, 'MSCondition'], convert_to_relative_abundance: bool = False) -> Optional[float]:
        """Get the expression value for a specific condition.

        Args:
            condition: MSCondition object or condition ID string
            convert_to_relative_abundance: If True, convert value to relative abundance

        Returns:
            The expression value, or None if not found
        """
        # Resolve condition to condition_id and condition object
        if isinstance(condition, str):
            if condition not in self.parent.conditions:
                logger.warning(
                    "Condition " + condition + " not found in expression object!"
                )
                return None
            condition_id = condition
            condition_obj = self.parent.conditions.get_by_id(condition)
        else:
            condition_id = condition.id
            condition_obj = condition

        # Check if value exists in DataFrame
        if self.id not in self.parent._data.index or condition_id not in self.parent._data.columns:
            logger.info(
                "Condition " + condition_id + " has no value in " + self.feature.id
            )
            return None

        # Get value from DataFrame and convert NaN to None
        value = self.parent._data.loc[self.id, condition_id]

        # Handle duplicate indices - if loc returns a Series, take the last value
        if isinstance(value, pd.Series):
            value = value.iloc[-1]

        if pd.isna(value):
            logger.info(
                "Condition " + condition_id + " has no value in " + self.feature.id
            )
            return None

        # Apply relative abundance conversion if requested
        if convert_to_relative_abundance:
            if self.parent.type == "AbsoluteAbundance":
                value = value / condition_obj.column_sum
            elif self.parent.type == "FPKM":
                value = value / condition_obj.column_sum
            elif self.parent.type == "TPM":
                value = value / condition_obj.column_sum
            elif self.parent.type == "Log2":
                value = 2 ** (value - condition_obj.lowest_value()) / 2 ** (condition_obj.column_sum - self.parent.features.len() * condition_obj.lowest_value())

        return value

class MSExpression:
    def __init__(self, type):
        self.type = type  # RelativeAbundance, AbsoluteAbundance, FPKM, TPM, Log2, NormalizedRatios
        self.object = None
        self.features = DictList()
        self.conditions = DictList()
        self._data = pd.DataFrame()
        self._data.index.name = 'feature_id'

    @staticmethod
    def from_msexpression(msexpression: 'MSExpression') -> 'MSExpression':
        """Create a copy of an existing MSExpression object.

        Args:
            msexpression: The MSExpression object to copy

        Returns:
            A new MSExpression object with the same data
        """
        new_expression = MSExpression(msexpression.type)
        new_expression.object = msexpression.object
        # Copy features
        for feature in msexpression.features:
            new_expression.features.append(MSExpressionFeature(feature.feature, new_expression))
        # Copy conditions
        for condition in msexpression.conditions:
            new_expression.conditions.append(MSCondition(condition.id, new_expression))
        # Copy data DataFrame
        new_expression._data = msexpression._data.copy()
        return new_expression
    
    def from_dataframe(
        df: pd.DataFrame,
        genome_or_model: Union['MSGenome', 'Model'],
        create_missing_features: bool = False,
        ignore_columns: list = None,
        description_column: Optional[str] = None,
        id_column: Optional[str] = None,
        id_translation: Optional[dict] = None,
        type: str = "RelativeAbundance"
    ) -> 'MSExpression':
        """Create an MSExpression object from a pandas DataFrame.

        Args:
            df: DataFrame with feature IDs and condition values
            genome: MSGenome object (optional)
            create_missing_features: If True, create features not in genome
            ignore_columns: List of column names to ignore
            description_column: Name of column containing descriptions
            id_column: Name of column containing feature IDs (default: first column)
            type: Expression data type (RelativeAbundance, AbsoluteAbundance, FPKM, TPM, Log2)

        Returns:
            MSExpression object with data loaded from DataFrame
        """
        if ignore_columns is None:
            ignore_columns = []

        expression = MSExpression(type)
        if genome_or_model is None:
            expression.object = MSGenome()
            create_missing_features = True
        else:
            expression.object = genome_or_model

        # Identify columns
        headers = list(df.columns)
        if id_column is None:
            id_column = headers[0]

        # Identify condition columns
        conditions = []
        description_present = description_column is not None and description_column in headers

        for header in headers:
            if header == id_column:
                continue
            elif header == description_column:
                continue
            elif header not in ignore_columns:
                conditions.append(header)
                if header not in expression.conditions:
                    expression.conditions.append(MSCondition(header, expression))
                # Initialize metadata attributes
                expression.conditions.get_by_id(header).column_sum = 0
                expression.conditions.get_by_id(header).feature_count = 0

        # Add features to the expression object
        valid_feature_ids = []
        for index, row in df.iterrows():
            gene_id = row[id_column]
            if id_translation is not None and gene_id in id_translation:
                gene_id = id_translation[gene_id]
            description = None
            if description_present:
                description = row[description_column]
            protfeature = expression.add_feature(
                gene_id, create_missing_features, description=description
            )
            if protfeature is not None:
                valid_feature_ids.append(protfeature.id)

        # Bulk load data into DataFrame
        if len(valid_feature_ids) > 0 and len(conditions) > 0:
            # Apply ID translation to the dataframe's ID column if provided
            if id_translation is not None:
                # Create a translated version of the ID column for filtering and indexing
                df_translated = df.copy()
                df_translated['_translated_id'] = df_translated[id_column].map(
                    lambda x: id_translation.get(x, x)
                )
                # Filter using translated IDs
                data_df = df_translated[df_translated['_translated_id'].isin(valid_feature_ids)].copy()
                # Set index using translated IDs
                data_df = data_df.set_index('_translated_id')
            else:
                # Extract numeric data columns without translation
                data_df = df[df[id_column].isin(valid_feature_ids)].copy()
                data_df = data_df.set_index(id_column)
            data_df = data_df[conditions]

            # Convert to numeric, coercing errors to NaN
            for col in conditions:
                data_df[col] = pd.to_numeric(data_df[col], errors='coerce')

            # Assign to expression._data
            expression._data = data_df
            expression._data.index.name = 'feature_id'

        return expression
    
    @staticmethod
    def from_spreadsheet(
        filename: str,
        sheet_name: Union[str, int] = 0,
        skiprows: int = 0,
        genome_or_model: Union['MSGenome', 'Model'] = None,
        create_missing_features: bool = True,
        id_translation: Optional[dict] = None,
        ignore_columns: list = None,
        description_column: Optional[str] = None,
        id_column: Optional[str] = None,
        type: str = "RelativeAbundance"
    ) -> 'MSExpression':
        """Create an MSExpression object from an Excel spreadsheet.

        Args:
            filename: Path to Excel file
            sheet_name: Sheet name or index (default: 0)
            skiprows: Number of rows to skip at start
            genome: MSGenome object (optional)
            create_missing_features: If True, create features not in genome
            ignore_columns: List of column names to ignore
            description_column: Name of column containing descriptions
            id_column: Name of column containing feature IDs (default: first column)
            type: Expression data type (RelativeAbundance, AbsoluteAbundance, FPKM, TPM, Log2)

        Returns:
            MSExpression object with data loaded from spreadsheet
        """
        df = pd.read_excel(filename, sheet_name=sheet_name, skiprows=skiprows)
        return MSExpression.from_dataframe(
            df,
            genome_or_model=genome_or_model,
            create_missing_features=create_missing_features,
            ignore_columns=ignore_columns,
            description_column=description_column,
            id_column=id_column,
            id_translation=id_translation,
            type=type
        )

    @staticmethod
    def from_gene_feature_file(
        filename: str,
        genome: Optional['MSGenome'] = None,
        create_missing_features: bool = False,
        ignore_columns: list = None,
        description_column: Optional[str] = None,
        sep: str = "\t",
        id_column: Optional[str] = None,
        type: str = "RelativeAbundance"
    ) -> 'MSExpression':
        """Create an MSExpression object from a delimited text file.

        Args:
            filename: Path to delimited file
            genome: MSGenome object (optional)
            create_missing_features: If True, create features not in genome
            ignore_columns: List of column names to ignore
            description_column: Name of column containing descriptions
            sep: Field delimiter (default: tab)
            id_column: Name of column containing feature IDs (default: first column)
            type: Expression data type (RelativeAbundance, AbsoluteAbundance, FPKM, TPM, Log2)

        Returns:
            MSExpression object with data loaded from file
        """
        df = pd.read_csv(filename, sep=sep)
        return MSExpression.from_dataframe(
            df,
            genome=genome,
            create_missing_features=create_missing_features,
            ignore_columns=ignore_columns,
            description_column=description_column,
            id_column=id_column,
            type=type
        )

    @staticmethod
    def load_from_dict(
        data_dict: dict,
        genome_or_model: Union['MSGenome', 'Model'],
        value_type: str,
        create_missing_features: bool = False
    ) -> 'MSExpression':
        """Create an MSExpression object from a dictionary.

        The dictionary should have feature IDs as keys and nested dictionaries
        mapping condition IDs to expression values as values.

        Example:
            data = {
                "gene1": {
                    "condition1": 10.5,
                    "condition2": 20.3
                },
                "gene2": {
                    "condition1": 8.2,
                    "condition2": 15.7
                }
            }
            expr = MSExpression.load_from_dict(
                data,
                genome_or_model=genome,
                value_type="Log2"
            )

        Args:
            data_dict: Dictionary with feature IDs as keys and condition-value
                dictionaries as values
            genome_or_model: MSGenome object (for gene expression) or Model
                object (for reaction expression). Required.
            value_type: Expression data type (RelativeAbundance, AbsoluteAbundance,
                FPKM, TPM, Log2, NormalizedRatios). Required.
            create_missing_features: If True, create features not in genome/model

        Returns:
            MSExpression object with data loaded from dictionary
        """
        # Create expression object
        expression = MSExpression(value_type)
        expression.object = genome_or_model

        # Convert dictionary to DataFrame
        # The dictionary format is: {feature_id: {condition_id: value, ...}, ...}
        if not data_dict:
            return expression

        # Convert to DataFrame with feature IDs as index and conditions as columns
        data_df = pd.DataFrame.from_dict(data_dict, orient='index').T

        if 'Feature ID' not in data_df.columns:
            # If 'Feature ID' is the index, reset it
            data_df = data_df.reset_index()
            if 'index' in data_df.columns:
                data_df = data_df.rename(columns={'index': 'Feature ID'})

        return MSExpression.from_dataframe(
            df=data_df,
            genome_or_model=genome_or_model,
            create_missing_features=create_missing_features,
            type=value_type,
            id_column='Feature ID'
        )

    def add_feature(
        self,
        id: str,
        create_gene_if_missing: bool = False,
        description: Optional[str] = None
    ) -> Optional['MSExpressionFeature']:
        """Add a feature to the expression object.

        Args:
            id: Feature ID
            create_gene_if_missing: If True, create the gene in genome if missing
            description: Optional feature description

        Returns:
            MSExpressionFeature object, or None if feature not found
        """
        if id in self.features:
            return self.features.get_by_id(id)
        feature = None
        # Check if object is MSGenome (gene expression) or Model (reaction expression)
        if isinstance(self.object, MSGenome):
            if self.object.search_for_gene(id) is None:
                if create_gene_if_missing:
                    self.object.features.append(MSFeature(id, ""))
            feature = self.object.search_for_gene(id)
        else:
            # Assume it's a COBRApy Model with reactions
            if hasattr(self.object, 'reactions') and id in self.object.reactions:
                feature = self.object.reactions.get_by_id(id)
            elif hasattr(self.object.model, 'reactions') and id in self.object.model.reactions:
                feature = self.object.model.reactions.get_by_id(id)
        if feature is None:
            logger.warning(
                "Feature referred by expression " + str(id)+ " not found in genome object!"
            )
            return None
        if feature.id in self.features:
            return self.features.get_by_id(feature.id)
        protfeature = MSExpressionFeature(feature, self)
        self.features.append(protfeature)
        return protfeature

    def get_value(
        self,
        feature: Union[str, 'MSExpressionFeature'],
        condition: Union[str, 'MSCondition']
    ) -> Optional[float]:
        """Get expression value for a feature and condition.

        Args:
            feature: MSExpressionFeature object or feature ID string
            condition: MSCondition object or condition ID string

        Returns:
            The expression value, or None if not found
        """
        if isinstance(feature, str):
            if feature not in self.features:
                logger.warning(
                    "Feature " + feature + " not found in expression object!"
                )
                return None
            feature = self.features.get_by_id(feature)
        return feature.get_value(condition)

    def build_reaction_expression(self, model, default: Optional[float] = None) -> 'MSExpression':
        """Build reaction-level expression from gene-level expression using GPR rules.

        Args:
            model: COBRApy Model object
            default: Default value for missing genes

        Returns:
            MSExpression object with reaction-level expression data
        """
        # Creating the expression and features
        rxnexpression = MSExpression(self.type)
        rxnexpression.object = model
        for rxn in model.reactions:
            if len(rxn.genes) > 0:
                rxnexpression.add_feature(rxn.id)
        for condition in self.conditions:
            newcondition = MSCondition(condition.id, rxnexpression)
            rxnexpression.conditions.append(newcondition)

        # Pulling the gene values from the current expression using DataFrame
        values = {}
        for gene in model.genes:
            # First, try to find the gene directly in expression features
            # This handles cases where expression was loaded without a genome
            if gene.id in self.features:
                feature = self.features.get_by_id(gene.id)
            else:
                # Fallback: search through the genome object (supports aliases)
                feature = self.object.search_for_gene(gene.id)
                if feature is None:
                    logger.debug(
                        "Model gene " + gene.id + " not found in genome or expression"
                    )
                    continue
                if feature.id not in self.features:
                    logger.debug(
                        "Model gene " + gene.id + " in genome but not in expression"
                    )
                    continue
                feature = self.features.get_by_id(feature.id)

            # Extract expression values for this feature
            for condition in self.conditions:
                if condition.id not in values:
                    values[condition.id] = {}
                # Get value from DataFrame instead of feature.values dictionary
                if feature.id in self._data.index and condition.id in self._data.columns:
                    value = self._data.loc[feature.id, condition.id]
                    if not pd.isna(value):
                        values[condition.id][gene.id] = value

        # Computing the reaction level values
        for condition in rxnexpression.conditions:
            for feature in rxnexpression.features:
                tree = GPR().from_string(str(feature.feature.gene_reaction_rule))
                feature.add_value(
                    condition, compute_gene_score(tree, values[condition.id], default, self.type)
                )
        return rxnexpression
    
    def get_dataframe(self, reset_index: bool = False) -> pd.DataFrame:
        """Get a DataFrame with expression data.

        Args:
            reset_index: If True, move feature_id from index to column (default: False)

        Returns:
            DataFrame with feature IDs as index (or column if reset_index=True)
            and conditions as columns
        """
        if reset_index:
            return self._data.reset_index()
        else:
            return self._data.copy()

    def add_row_to_data(
        self,
        feature_id: str,
        values: Optional[dict] = None,
        default_value: float = 0.0,
        create_feature_if_missing: bool = True,
        overwrite: bool = False
    ) -> bool:
        """Add a row to the expression data for a feature.

        This method adds both a row to the underlying DataFrame and creates
        the corresponding MSExpressionFeature object in the features list.

        Args:
            feature_id: The ID of the feature to add
            values: Optional dictionary mapping condition IDs to values.
                    If None, all conditions will use default_value.
                    If a condition is missing from the dict, it will use default_value.
            default_value: Default value to use for conditions not in values dict (default: 0.0)
            create_feature_if_missing: If True, create the feature in the genome/model
                                       if it doesn't exist (default: True)
            overwrite: If True, overwrite existing data if the feature already exists
                       (default: False)

        Returns:
            True if the row was added/updated, False if the feature already exists
            and overwrite is False

        Example:
            # Add a gene with 0.0 for all conditions
            expression.add_row_to_data("gene123")

            # Add a gene with specific values
            expression.add_row_to_data("gene456", {"condition1": 0.5, "condition2": 0.8})

            # Add a gene with default value of 1.0
            expression.add_row_to_data("gene789", default_value=1.0)

            # Overwrite existing data for a gene
            expression.add_row_to_data("gene123", {"condition1": 0.9}, overwrite=True)
        """
        # Check if feature already exists in the data
        exists_in_data = feature_id in self._data.index
        exists_in_features = feature_id in self.features

        if exists_in_data or exists_in_features:
            if not overwrite:
                logger.warning(f"Feature {feature_id} already exists in expression data")
                return False
            else:
                logger.info(f"Overwriting existing data for feature {feature_id}")

        # Build the row values
        row_values = {}
        for condition in self.conditions:
            if values is not None and condition.id in values:
                row_values[condition.id] = values[condition.id]
            else:
                row_values[condition.id] = default_value

        # Add/update row in DataFrame
        self._data.loc[feature_id] = row_values

        # If feature already exists in features list, we're done (just updated data)
        if exists_in_features:
            logger.info(f"Updated data for feature {feature_id} with {len(row_values)} condition values")
            return True

        # Create and add the MSExpressionFeature
        # First, try to get/create the underlying feature object
        feature = None
        if isinstance(self.object, MSGenome):
            feature = self.object.search_for_gene(feature_id)
            if feature is None and create_feature_if_missing:
                self.object.features.append(MSFeature(feature_id, ""))
                feature = self.object.search_for_gene(feature_id)
        elif hasattr(self.object, 'reactions') and feature_id in self.object.reactions:
            feature = self.object.reactions.get_by_id(feature_id)
        elif hasattr(self.object, 'model') and hasattr(self.object.model, 'reactions'):
            if feature_id in self.object.model.reactions:
                feature = self.object.model.reactions.get_by_id(feature_id)

        # If we still don't have a feature, create a simple placeholder
        if feature is None:
            if create_feature_if_missing:
                # Create a simple feature-like object
                class SimpleFeature:
                    def __init__(self, fid):
                        self.id = fid
                feature = SimpleFeature(feature_id)
            else:
                logger.warning(f"Could not find or create feature {feature_id}")
                # Still keep the data row, just don't add to features list
                return True

        # Add to features list
        expr_feature = MSExpressionFeature(feature, self)
        self.features.append(expr_feature)

        logger.info(f"Added feature {feature_id} to expression data with {len(row_values)} condition values")
        return True
    
    def translate_data(self, target_type: str) -> 'MSExpression':
        """Translate expression data to a different type.

        Args:
            target_type: Target expression type (RelativeAbundance, AbsoluteAbundance, FPKM, TPM, Log2)

        Returns:
            New MSExpression object with translated data
        """
        # Create a copy of the current expression object
        new_expression = MSExpression.from_msexpression(self)
        new_expression.type = target_type
        # Perform translation based on source and target types
        for condition in self.conditions:
            denominator = None
            for feature in self.features:
                value = feature.get_value(condition)
                if value is not None:
                    if self.type == "AbsoluteAbundance" or self.type == "FPKM" or self.type == "TPM":
                        if target_type == "RelativeAbundance":
                            if condition.sum_value() > 0.01:
                                transformed_value = value / condition.sum_value()
                            else:
                                transformed_value = 0
                        elif target_type == "NormalizedRatios":
                            if condition.highest_value() > 0.01:
                                transformed_value = value / condition.highest_value()
                            else:
                                transformed_value = 0
                        else:
                            raise ValueError(
                                f"Translation from {self.type} to {target_type} not supported"
                            )
                        new_expression._data.loc[feature.id, condition.id] = transformed_value
                    elif self.type == "Log2":
                        if target_type == "RelativeAbundance":
                            if denominator is None:
                                denominator = 0
                                for ftr in self.features:
                                    denominator += 2 ** ftr.get_value(condition)
                            numerator = 2 ** value
                            if denominator > 0.01:
                                transformed_value = numerator / denominator
                            else:
                                transformed_value = 0
                        else:
                            raise ValueError(
                                f"Translation from {self.type} to {target_type} not supported"
                            )                            
                        new_expression._data.loc[feature.id, condition.id] = transformed_value
                    else:
                        raise ValueError(
                            f"Translation from {self.type} to {target_type} not supported"
                        )
                    
        return new_expression

    def average_expression_replicates(self, strain_list: list) -> 'MSExpression':
        """Average expression replicates for each strain.

        Takes an MSExpression object with replicate columns (e.g., ACN2586_1, ACN2586_2, ...)
        and averages them to create single columns per strain (e.g., ACN2586).

        Args:
            strain_list: List of strain names (e.g., ["ACN2586", "ACN2821", ...])

        Returns:
            New MSExpression object with averaged data per strain

        Raises:
            ValueError: If no data found for any strain in the list
        """
        try:
            # Access the underlying DataFrame
            expression_df = self._data.copy()

            # Create new DataFrame for averaged data
            averaged_data = {}

            # Keep the index (gene/protein IDs)
            averaged_data['index'] = expression_df.index

            # For each strain, find and average its replicates
            for strain in strain_list:
                # Find columns that match this strain pattern (e.g., ACN2586_1, ACN2586_2, ...)
                replicate_cols = [col for col in expression_df.columns if col.startswith(f"{strain}_")]

                if replicate_cols:
                    # Average the replicates
                    averaged_data[strain] = expression_df[replicate_cols].mean(axis=1)
                    logger.info(f"Averaged {len(replicate_cols)} replicates for strain {strain}")
                else:
                    # No replicates found - check if strain column exists as-is
                    if strain in expression_df.columns:
                        averaged_data[strain] = expression_df[strain]
                        logger.info(f"No replicates found for {strain}, using existing column")
                    else:
                        logger.warning(f"No data found for strain {strain}")

            # Create new DataFrame from averaged data
            averaged_df = pd.DataFrame(averaged_data)
            averaged_df.set_index('index', inplace=True)

            # Create a deep copy of the expression object
            averaged_expression = copy.deepcopy(self)

            # Replace the data with averaged data
            averaged_expression._data = averaged_df

            # Update conditions list to match new columns
            # Clear and rebuild conditions using proper MSCondition class
            averaged_expression.conditions = DictList()
            for strain in strain_list:
                if strain in averaged_df.columns:
                    condition = MSCondition(strain, averaged_expression)
                    averaged_expression.conditions.append(condition)

            logger.info(f"Created averaged expression data with {len(averaged_expression.conditions)} conditions")

            return averaged_expression

        except Exception as e:
            logger.error(f"Error averaging expression replicates: {str(e)}")
            raise

    def fit_model_expression_to_data(
        self,
        model: 'MSModelUtil',
        condition: str,
        default_coef: float = 0.00001,
        activation_threshold: float = None,  # cshenry 10/16/2026: Changed default for activation to None so activation will be off by default
        deactivation_threshold: float = 0.000001,
        on_coef_override: float = None,
        off_coef_override: float = None,
        analyze_solution_essentiality: bool = False
    ) -> Solution:
        """Fit metabolic model fluxes to expression data using threshold-based constraints.

        This function integrates gene or reaction expression data with a metabolic model
        to predict fluxes that are consistent with observed expression patterns. Reactions
        with high expression are encouraged (activated), while reactions with low expression
        are discouraged (deactivated).

        The function automatically handles:
        - Conversion from gene-level to reaction-level expression (if needed)
        - Transformation of expression data types to RelativeAbundance
        - Construction of activation/deactivation dictionaries
        - Integration with ExpressionActivationPkg for constrained optimization

        Parameters
        ----------
        model : MSModelUtil
            The metabolic model to fit to expression data. Must contain the same reactions
            as the expression model (if expression is already at reaction level).
        condition : str
            The condition ID whose expression values will be used for fitting. Must exist
            in the expression data conditions.
        activation_threshold : float, optional
            Minimum relative abundance value for a reaction to be considered "active".
            Reactions with expression above this threshold will be encouraged in the
            optimization. Default: 0.002 (0.2% relative abundance).
        deactivation_threshold : float, optional
            Maximum relative abundance value for a reaction to be considered "inactive".
            Reactions with expression below this threshold will be discouraged in the
            optimization. Default: 0.000001 (0.0001% relative abundance).

        Returns
        -------
        cobra.Solution
            The optimized FBA solution with expression-based flux constraints. Contains:
            - objective_value: The optimized objective value
            - fluxes: pandas Series of reaction fluxes
            - status: Optimization status ('optimal', 'infeasible', etc.)

        Raises
        ------
        ValueError
            If the expression model does not match the input model (when expression is
            already at reaction level).
        ValueError
            If the specified condition does not exist in the expression data.
        ValueError
            If activation_threshold is less than or equal to deactivation_threshold.
        ValueError
            If no reactions have expression values above the activation threshold.
        RuntimeError
            If the optimization fails or produces an infeasible solution.

        Notes
        -----
        - This function does NOT modify the original expression data or model
        - All model modifications occur within a context manager and are reverted
        - Reactions without expression data or reactions between the specified thresholds are deactivated base on the default_coef argument
        - The function uses ExpressionActivationPkg internally for constraint building

        See Also
        --------
        MSExpression.build_reaction_expression : Convert gene to reaction expression
        ExpressionActivationPkg.build_package : Build expression-based constraints
        MSModelUtil.test_single_condition : Run FBA on a single condition

        References
        ----------
        .. [1] Becker, S. A., & Palsson, B. O. (2008). Context-specific metabolic networks
           are consistent with experiments. PLoS computational biology, 4(5), e1000082.
        """
        # cshenry 10/16/2026: Checking that model is MSModelUtil and converting if not
        if not isinstance(model, MSModelUtil):
            model = MSModelUtil(model)
        
        # Task 1.6: Initial logging
        logger.info(f"Fitting model flux to expression or fitness data for condition: {condition}")

        # Task 1.4: Threshold validation
        if activation_threshold is not None and activation_threshold <= deactivation_threshold:
            raise ValueError(
                f"activation_threshold ({activation_threshold}) must be greater than "
                f"deactivation_threshold ({deactivation_threshold})"
            )

        # Task 1.5: Condition validation
        if condition not in self.conditions:
            available_conditions = [c.id for c in self.conditions]
            raise ValueError(
                f"Condition '{condition}' not found in expression data. "
                f"Available conditions: {available_conditions}"
            )

        # Task 2.1-2.4: Genome vs model detection and conversion
        if isinstance(self.object, MSGenome):
            # Task 2.2: Genome-level expression - convert to reaction-level
            rxn_expression = self.build_reaction_expression(model.model)
        else:
            # Task 2.3: Model-level expression - validate model match
            expr_rxns = set(self.object.model.reactions.list_attr('id'))
            model_rxns = set(model.model.reactions.list_attr('id'))

            if expr_rxns != model_rxns:
                # Task 2.4: Model mismatch error
                missing_in_expr = list(model_rxns - expr_rxns)[:10]
                missing_in_model = list(expr_rxns - model_rxns)[:10]
                raise ValueError(
                    f"Models must match when fitting flux to data. "
                    f"Expression model has {len(expr_rxns)} reactions, "
                    f"input model has {len(model_rxns)} reactions. "
                    f"Missing in expression: {missing_in_expr}. "
                    f"Missing in model: {missing_in_model}"
                )

            rxn_expression = self
        
        # Task 2.5-2.13: Expression type transformation
        if rxn_expression.type != "RelativeAbundance" and rxn_expression.type != "NormalizedRatios":
            raise ValueError(
                f"Reaction expression must be in terms of relative abundance or normalized ratios"
            )

        # Initialize empty dictionaries
        on_hash = {}
        off_hash = {}

        # Iterate through reactions and build dictionaries
        for rxn in model.model.reactions:
            expr_value = rxn_expression.get_value(rxn.id, condition)
            if expr_value is None:
                continue
            if rxn_expression.type == "NormalizedRatios":
                if activation_threshold is not None and abs(1 - expr_value) >= activation_threshold:
                    on_hash[rxn.id] = 10 * (1 - expr_value)
                elif abs(1-expr_value) <= deactivation_threshold:
                    off_hash[rxn.id] = 10*(1 - expr_value)
            else:
                if activation_threshold is not None and expr_value > activation_threshold:
                    if activation_threshold != 0:
                        on_hash[rxn.id] = (expr_value - activation_threshold)/activation_threshold
                    else:
                        on_hash[rxn.id] = expr_value+1
                elif expr_value < deactivation_threshold:
                    if deactivation_threshold != 0:
                        off_hash[rxn.id] = (deactivation_threshold - expr_value)/deactivation_threshold
                    else:
                        off_hash[rxn.id] = expr_value+1

        print("On:", on_hash)
        print("Off:", off_hash)

        # Log dictionary sizes
        logger.info(f"Identified {len(on_hash)} reactions for activation (above threshold {activation_threshold})")
        logger.info(f"Identified {len(off_hash)} reactions for deactivation (below threshold {deactivation_threshold})")

        # Access package manager
        pkgmgr = model.pkgmgr

        # Get ExpressionActivationPkg
        expr_pkg = pkgmgr.getpkg("ExpressionActivationPkg")

        # Use context manager for transient modifications
        output = {"on_on":[],"on_off":[], "off_on":[], "off_off":[],"none_on":[],"none_off":[],"on_on_reduced":[],"off_on_reduced":[],"none_on_reduced":[],"solution":None}
        original_objective = model.model.objective
        with model.model:
            # Task 4.4: Build package with dictionaries
            expr_pkg.build_package(on_hash, off_hash, other_coef=default_coef, on_coeff=on_coef_override, off_coeff=off_coef_override)
            # Task 4.5: Execute optimization
            output["objective"] = model.model.objective
            output["solution"] = model.model.optimize()
            for rxn in model.model.reactions:
                if rxn.id in on_hash:
                    if abs(output["solution"].fluxes[rxn.id]) > 1e-6:
                        output["on_on"].append(rxn.id)
                    else:
                        output["on_off"].append(rxn.id)
                elif rxn.id in off_hash:
                    if abs(output["solution"].fluxes[rxn.id]) > 1e-6:
                        output["off_on"].append(rxn.id)
                    else:
                        output["off_off"].append(rxn.id)
                else:
                    if abs(output["solution"].fluxes[rxn.id]) > 1e-6:
                        output["none_on"].append(rxn.id)
                    else:
                        output["none_off"].append(rxn.id)

            # Task 4.6: Validate solution status
            if output["solution"].status != "optimal":
                raise RuntimeError(
                    f"Optimization failed with status: {output['solution'].status}. "
                    f"The model may be infeasible with the given expression constraints."
                )

            # Task 4.7: Log optimization result
            logger.info(f"Optimization completed with objective value: {output['solution'].objective_value}")

        # Categorize reactions by flux
        zero_flux_rxns = []
        active_rxns = []
        
        for rxn_id, flux in output["solution"].fluxes.items():
            if rxn_id not in [r.id for r in model.model.reactions]:
                continue
            if abs(flux) <= 1e-9:
                zero_flux_rxns.append(rxn_id)
            else:
                active_rxns.append((rxn_id, flux))
        
        print(f"  Zero-flux reactions: {len(zero_flux_rxns)}")
        print(f"  Active reactions: {len(active_rxns)}")
        
        with model.model:
            #model.model.objective = original_objective
            # Set zero-flux reactions to have zero bounds
            for rxn_id in zero_flux_rxns:
                rxn = model.model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            
            # Get baseline growth with constrained model
            output["baseline_growth"] = model.model.optimize().objective_value
            
            # Test each active reaction knockout
            essentiality_results = {}
            essential_count = 0
            reduced_count = 0
            
            for rxn_id, original_flux in active_rxns:
                rxn = model.model.reactions.get_by_id(rxn_id)
                
                # Save original bounds
                orig_lb = rxn.lower_bound
                orig_ub = rxn.upper_bound
                
                # Knock out the reaction
                rxn.lower_bound = 0
                rxn.upper_bound = 0
                
                # Optimize
                ko_solution = model.model.optimize()
                
                if ko_solution.status == 'optimal':
                    ko_growth = ko_solution.objective_value
                    growth_ratio = ko_growth / baseline_growth if baseline_growth > 0 else 0
                else:
                    ko_growth = 0
                    growth_ratio = 0
                
                # Categorize impact
                if growth_ratio < 0.01:
                    impact = "essential"
                    essential_count += 1
                elif growth_ratio < 0.95:
                    impact = "reduced"
                    reduced_count += 1
                else:
                    impact = "dispensable"
                
                essentiality_results[rxn_id] = {
                    "expression_data_status":"none",
                    "original_flux": original_flux,
                    "ko_growth": ko_growth,
                    "growth_ratio": growth_ratio,
                    "impact": impact
                }
                if rxn_id in on_hash:
                    essentiality_results[rxn_id]["expression_data_status"] = "on"
                    if growth_ratio < 0.95:
                        output["on_on_reduced"].append(rxn_id)
                elif rxn_id in off_hash:
                    essentiality_results[rxn_id]["expression_data_status"] = "off"
                    if growth_ratio < 0.95:
                        output["off_on_reduced"].append(rxn_id)
                else:
                    essentiality_results[rxn_id]["expression_data_status"] = "none"
                    if growth_ratio < 0.95:
                        output["none_on_reduced"].append(rxn_id)

                # Restore original bounds
                rxn.lower_bound = orig_lb
                rxn.upper_bound = orig_ub
        
        print(f"  Essential reactions: {essential_count}")
        print(f"  Reduced growth reactions: {reduced_count}")
        print(f"  Dispensable reactions: {len(active_rxns) - essential_count - reduced_count}")
        
        output["baseline_growth"] = baseline_growth
        output["zero_flux_count"] = len(zero_flux_rxns)
        output["active_count"] = len(active_rxns)
        output["essential_count"] = essential_count
        output["reduced_count"] = reduced_count
        output["reactions"] = essentiality_results

        # Task 4.8: Return solution
        return output


    def fit_flux_to_mutant_growth_rate_data(
        self,
        model: 'MSModelUtil',
        condition: str,
        default_coef: float = 0.00001,
        activation_threshold: float = 0.90,
        deactivation_threshold: float = 0.95,
        on_coef_override: float = None,
        off_coef_override: float = None,
        use_activation_constraints: bool = False
    ) -> Solution:
        """Fit metabolic model fluxes to mutant growth rate data using threshold-based constraints
        """
        output = {"on_on":[],"on_off":[], "off_on":[], "off_off":[],"none_on":[],"none_off":[],
                  "on_genes":[],"off_genes":[],"genes_on_on":{},"genes_on_off":{},"genes_off_on":{},"genes_off_off":{},"genes_on_norxn":[],"genes_off_norxn":[],"solution":None,}
        if not isinstance(model, MSModelUtil):
            model = MSModelUtil(model)

        logger.info(f"Fitting model flux to mutant growth rate data for condition: {condition}")

        if activation_threshold >= deactivation_threshold:
            raise ValueError(
                f"activation_threshold ({activation_threshold}) must be less than "
                f"deactivation_threshold ({deactivation_threshold})"
            )

        if condition not in self.conditions:
            available_conditions = [c.id for c in self.conditions]
            raise ValueError(
                f"Condition '{condition}' not found in expression data. "
                f"Available conditions: {available_conditions}"
            )

        if self.type != "NormalizedRatios":
            raise ValueError(
                f"Expression must be in terms of normalized ratios"
            )

        # Initialize empty dictionaries
        on_hash = {}
        off_hash = {}

        # Track gene-level on/off status
        on_genes_set = set()
        off_genes_set = set()

        # First pass: categorize all genes by their expression values
        for feature in self.features:
            gene_id = feature.id if hasattr(feature, 'id') else str(feature)
            expr_value = self.get_value(gene_id, condition)
            if expr_value is not None:
                if expr_value <= activation_threshold:
                    on_genes_set.add(gene_id)
                elif expr_value >= deactivation_threshold:
                    off_genes_set.add(gene_id)

        # Iterate through reactions and build dictionaries
        gene_rxns = {}
        for rxn in model.model.reactions:
            # Check if the reaction ID is in the expression data
            if rxn.id in self.features:
                expr_value = self.get_value(rxn.id, condition) # No genes if reaction has direct expression data
            else:
                lowest_value = None
                inducing_gene = None
                for gene in rxn.genes:
                    gene_id = gene.id  # Extract string ID from Gene object
                    if gene_id in self.features:
                        expr_value = self.get_value(gene_id, condition)
                        if expr_value is not None:
                            if lowest_value is None or expr_value < lowest_value:
                                lowest_value = expr_value
                                inducing_gene = gene_id
                if lowest_value is None:
                    expr_value = None
                else:
                    expr_value = lowest_value
            if expr_value is None:
                continue
            if expr_value <= activation_threshold:
                on_hash[rxn.id] = 1 + 10 * (activation_threshold - expr_value)
            elif expr_value >= deactivation_threshold:
                off_hash[rxn.id] = 1 + 10 * (expr_value-deactivation_threshold)
            for gene in rxn.genes:
                if gene.id not in gene_rxns:
                    gene_rxns[gene.id] = {}
                gene_rxns[gene.id][rxn.id] = 1

        # Store gene lists in output
        output["on_genes"] = list(on_genes_set)
        output["off_genes"] = list(off_genes_set)

        print("On:", on_hash)
        print("Off:", off_hash)

        # Log dictionary sizes
        logger.info(f"Identified {len(on_hash)} reactions for activation (above threshold {activation_threshold})")
        logger.info(f"Identified {len(off_hash)} reactions for deactivation (below threshold {deactivation_threshold})")

        # Get ExpressionActivationPkg
        expr_pkg = model.pkgmgr.getpkg("ExpressionActivationPkg")

        # Use context manager for transient modifications
        original_objective = model.model.objective
        with model.model:
            expr_pkg.build_package(on_hash, off_hash, other_coef=default_coef, on_coeff=on_coef_override, off_coeff=off_coef_override,use_activation_constraints=use_activation_constraints)
            output["solution"] = model.model.optimize()
            if output["solution"].status != "optimal":
                raise RuntimeError(
                    f"Optimization failed with status: {output['solution'].status}. "
                    f"The model may be infeasible with the given expression constraints."
                )
            for rxn in model.model.reactions:
                if abs(output["solution"].fluxes[rxn.id]) > 1e-6:
                    if rxn.id in on_hash:
                        output["on_on"].append(rxn.id)
                    else:
                        output["on_off"].append(rxn.id)
                else:
                    if rxn.id in on_hash:
                        output["on_off"].append(rxn.id)
                    else:
                        output["off_off"].append(rxn.id)
            for gene in on_genes_set:
                if gene not in gene_rxns:
                    output["genes_on_norxn"].append(gene)
                    continue
                on_rxn_found = false
                for rxn_id  in gene_rxns[gene]:
                    if rxn_id in on_hash:
                        on_rxn_found = True
                        break
                if on_rxn_found:
                    if gene not in output["genes_on_on"]:
                        output["genes_on_on"][gene] = {"on_rxns":[],"off_rxns":[]}
                    for rxn_id in gene_rxns[gene]:
                        if rxn_id in on_hash:
                            output["genes_on_on"][gene]["on_rxns"].append(rxn_id)
                        elif rxn_id in off_hash:
                            output["genes_on_on"][gene]["off_rxns"].append(rxn_id)
                else:
                    if gene not in output["genes_on_off"]:
                        output["genes_on_off"][gene] = {"on_rxns":[],"off_rxns":[]}
                    for rxn_id in gene_rxns[gene]:
                        if rxn_id in on_hash:
                            output["genes_on_off"][gene]["on_rxns"].append(rxn_id)
                        elif rxn_id in off_hash:
                            output["genes_on_off"][gene]["off_rxns"].append(rxn_id)
            for gene in off_genes_set:
                if gene not in gene_rxns:
                    output["genes_off_norxn"].append(gene)
                    continue
                on_rxn_found = false
                for rxn_id in gene_rxns[gene]:
                    if rxn_id in on_hash:
                        on_rxn_found = True
                        break
                if on_rxn_found:
                    if gene not in output["genes_off_on"]:
                        output["genes_off_on"][gene] = {"on_rxns":[],"off_rxns":[]}
                    for rxn_id in gene_rxns[gene]:
                        if rxn_id in on_hash:
                            output["genes_off_on"][gene]["on_rxns"].append(rxn_id)
                        elif rxn_id in off_hash:
                            output["genes_off_on"][gene]["off_rxns"].append(rxn_id)
                else:
                    if gene not in output["genes_off_off"]:
                        output["genes_off_off"][gene] = {"on_rxns":[],"off_rxns":[]}
                    for rxn_id in gene_rxns[gene]:
                        if rxn_id in on_hash:
                            output["genes_off_off"][gene]["on_rxns"].append(rxn_id)
                        elif rxn_id in off_hash:
                            output["genes_off_off"][gene]["off_rxns"].append(rxn_id)

            logger.info(f"Optimization completed with objective value: {output['solution'].objective_value}")

        # Task 4.8: Return solution
        return output

    def fit_flux_to_proteomics_fold_change_data(
        self,
        model: 'MSModelUtil',
        reference_condition: str,
        target_condition: str,
        reference_flux: dict,
        zero_flux: float = 0.001,
        least_squares: bool = True,
        fold_change_thresholds: list = [1.0, 2.0, 3.0],
        quadratic_formulation: bool = True 
    ) -> dict:
        """Fit metabolic model fluxes to proteomics fold change data.

        This function uses proteomics log2 data to compute fold changes between a
        reference and target condition, then fits model fluxes to match those fold
        changes relative to a reference flux distribution.

        The function assumes the expression data is in log2 form, so fold change
        is computed as: 2^(target_value - reference_value)

        Parameters
        ----------
        model : MSModelUtil
            The metabolic model to fit to proteomics data.
        reference_condition : str
            The condition ID for the reference state (e.g., "wildtype", "control").
        target_condition : str
            The condition ID for the target state (e.g., "mutant", "treatment").
        reference_flux : dict
            Dictionary of {reaction_id: flux_value} representing the baseline flux
            distribution to scale by fold changes. Typically from pFBA on reference.
        zero_flux : float, optional
            Default flux value to use for reactions with zero reference flux but
            non-zero predicted flux. Default: 0.001
        least_squares : bool, optional
            If True, use least squares fitting (minimize sum of squared deviations).
            If False, use linear fitting. Default: True
        fold_change_thresholds : list, optional
            List of fold change thresholds (in log2 units) for categorizing
            significant changes. Default: [1.0, 2.0, 3.0] (2x, 4x, 8x fold changes)

        Returns
        -------
        dict
            Results dictionary containing:
            - 'solution': cobra.Solution object with fitted fluxes
            - 'target_flux': dict of {rxn_id: target_flux} used for fitting
            - 'fold_changes': dict of {rxn_id: fold_change} computed from proteomics
            - 'flux_changes': dict categorizing reactions by flux change magnitude:
                - 'increased_1std': reactions with >1 stddev increase
                - 'increased_2std': reactions with >2 stddev increase
                - 'increased_3std': reactions with >3 stddev increase
                - 'decreased_1std': reactions with >1 stddev decrease
                - 'decreased_2std': reactions with >2 stddev decrease
                - 'decreased_3std': reactions with >3 stddev decrease
            - 'protein_changes': dict categorizing proteins:
                - 'significant_implemented': proteins with significant fold change
                  whose reactions could match the change
                - 'significant_not_implemented': proteins with significant fold change
                  whose reactions could NOT match the change
            - 'reaction_matching': dict categorizing reactions:
                - 'matched': reactions that matched their gene fold change demands
                - 'not_matched': reactions that could not match gene demands
            - 'statistics': dict with summary statistics

        Raises
        ------
        ValueError
            If reference_condition or target_condition not in expression data.
        ValueError
            If expression data is not in log2 form.
        RuntimeError
            If optimization fails.

        Notes
        -----
        The function uses FluxFittingPkg internally to set up the optimization
        problem. Target fluxes are computed as: reference_flux * fold_change
        where fold_change = 2^(log2_target - log2_reference)
        """
        import numpy as np
        from modelseedpy.core.msmodelutl import MSModelUtil

        if not isinstance(model, MSModelUtil):
            model = MSModelUtil(model)

        logger.info(f"Fitting flux to proteomics fold change: {reference_condition} -> {target_condition}")

        # Validate conditions
        if reference_condition not in self.conditions:
            available = [c.id for c in self.conditions]
            raise ValueError(f"Reference condition '{reference_condition}' not found. Available: {available}")

        if target_condition not in self.conditions:
            available = [c.id for c in self.conditions]
            raise ValueError(f"Target condition '{target_condition}' not found. Available: {available}")

        # Compute fold changes for each reaction
        # Since data is log2, fold_change = 2^(target - reference)
        fold_changes = {}
        gene_fold_changes = {}

        # First compute gene-level fold changes
        for feature in self.features:
            gene_id = feature.id if hasattr(feature, 'id') else str(feature)
            ref_value = self.get_value(gene_id, reference_condition)
            target_value = self.get_value(gene_id, target_condition)

            if ref_value is not None and target_value is not None:
                # Log2 difference gives log2(fold_change)
                log2_fc = target_value - ref_value
                gene_fold_changes[gene_id] = {
                    'log2_fc': log2_fc,
                    'fold_change': 2 ** log2_fc,
                    'ref_value': ref_value,
                    'target_value': target_value
                }

        # Map gene fold changes to reactions (use minimum fold change like GPR AND logic)
        for rxn in model.model.reactions:
            if len(rxn.genes) == 0:
                continue

            rxn_fold_change = None
            contributing_genes = []

            for gene in rxn.genes:
                gene_id = gene.id
                if gene_id in gene_fold_changes:
                    fc = gene_fold_changes[gene_id]['fold_change']
                    contributing_genes.append(gene_id)
                    if rxn_fold_change is None or fc < rxn_fold_change:
                        rxn_fold_change = fc

            if rxn_fold_change is not None:
                fold_changes[rxn.id] = {
                    'fold_change': rxn_fold_change,
                    'log2_fc': np.log2(rxn_fold_change),
                    'contributing_genes': contributing_genes
                }

        # Compute target fluxes based on reference flux * fold change
        target_flux = {}
        for rxn_id, fc_data in fold_changes.items():
            if rxn_id in reference_flux:
                ref_flux = reference_flux[rxn_id]
                fold_change = fc_data['fold_change']
                if fold_change > 3:
                    fold_change = 3
                elif fold_change < 0.33:
                    fold_change = 0.33
                if abs(ref_flux) > 1e-9:
                    target_flux[rxn_id] = ref_flux * fc_data['fold_change']
                else:
                    # Use zero_flux as baseline if reference is zero but fold change suggests activity
                    if fc_data['fold_change'] > 1.5:  # If significantly upregulated
                        target_flux[rxn_id] = zero_flux * fc_data['fold_change']
            else:
                # Reaction not in reference flux - use zero_flux if upregulated
                if fc_data['fold_change'] > 1.5:
                    target_flux[rxn_id] = zero_flux * fc_data['fold_change']

        logger.info(f"Computed target fluxes for {len(target_flux)} reactions")

        # Set up flux fitting optimization
        output = {
            'solution': None,
            'target_flux': target_flux,
            'fold_changes': fold_changes,
            'gene_fold_changes': gene_fold_changes,
            'flux_changes': {
                'increased_1std': [],
                'increased_2std': [],
                'increased_3std': [],
                'decreased_1std': [],
                'decreased_2std': [],
                'decreased_3std': []
            },
            'protein_changes': {
                'significant_implemented': {},
                'significant_not_implemented': {}
            },
            'reaction_matching': {
                'matched': [],
                'not_matched': []
            },
            'statistics': {}
        }

        # Use FluxFittingPkg for optimization
        flux_fit_pkg = model.pkgmgr.getpkg("FluxFittingPkg")

        with model.model:
            flux_fit_pkg.build_package({
                'target_flux': target_flux,
                'set_objective': 1,
                'totalflux': 0,
                'quadratic_formulation': quadratic_formulation
            })

            # Optimize
            solution = model.model.optimize()

            if solution.status != 'optimal':
                raise RuntimeError(f"Optimization failed with status: {solution.status}")

            output['solution'] = solution

            # Analyze flux changes
            fitted_fluxes = solution.fluxes.to_dict()
            flux_deviations = []

            for rxn_id in target_flux:
                if rxn_id in fitted_fluxes:
                    target = target_flux[rxn_id]
                    fitted = fitted_fluxes[rxn_id]
                    if abs(target) > 1e-9:
                        deviation = (fitted - target) / abs(target)
                        flux_deviations.append(deviation)

            # Compute statistics
            if flux_deviations:
                flux_std = np.std(flux_deviations)
                flux_mean = np.mean(flux_deviations)
                output['statistics']['flux_deviation_mean'] = flux_mean
                output['statistics']['flux_deviation_std'] = flux_std
                output['statistics']['flux_deviation_count'] = len(flux_deviations)
            else:
                flux_std = 1.0
                output['statistics']['flux_deviation_std'] = 0

            # Categorize flux changes by standard deviations
            for rxn_id in reference_flux:
                if rxn_id not in fitted_fluxes:
                    continue

                ref_flux = reference_flux[rxn_id]
                fitted_flux = fitted_fluxes[rxn_id]

                if abs(ref_flux) < 1e-9:
                    ref_flux = zero_flux

                flux_ratio = fitted_flux / ref_flux if abs(ref_flux) > 1e-9 else 0
                log2_change = np.log2(abs(flux_ratio)) if flux_ratio != 0 else 0

                # Categorize by fold change thresholds
                for i, threshold in enumerate(fold_change_thresholds):
                    if log2_change > threshold:
                        output['flux_changes'][f'increased_{i+1}std'].append({
                            'rxn_id': rxn_id,
                            'ref_flux': ref_flux,
                            'fitted_flux': fitted_flux,
                            'log2_change': log2_change
                        })
                    elif log2_change < -threshold:
                        output['flux_changes'][f'decreased_{i+1}std'].append({
                            'rxn_id': rxn_id,
                            'ref_flux': ref_flux,
                            'fitted_flux': fitted_flux,
                            'log2_change': log2_change
                        })

            # Analyze protein/gene changes and implementation
            for gene_id, gc_data in gene_fold_changes.items():
                log2_fc = gc_data['log2_fc']

                # Check if this gene's fold change is significant (>1 threshold)
                if abs(log2_fc) < fold_change_thresholds[0]:
                    continue

                # Find reactions associated with this gene
                gene_rxns = []
                for rxn in model.model.reactions:
                    if any(g.id == gene_id for g in rxn.genes):
                        gene_rxns.append(rxn.id)

                if not gene_rxns:
                    output['protein_changes']['significant_not_implemented'][gene_id] = {
                        'log2_fc': log2_fc,
                        'fold_change': gc_data['fold_change'],
                        'reason': 'no_reactions_in_model'
                    }
                    continue

                # Check if any of the gene's reactions showed matching flux change
                implemented = False
                for rxn_id in gene_rxns:
                    if rxn_id in fitted_fluxes and rxn_id in reference_flux:
                        ref_flux = reference_flux[rxn_id]
                        fitted_flux = fitted_fluxes[rxn_id]

                        if abs(ref_flux) > 1e-9:
                            flux_log2_change = np.log2(abs(fitted_flux / ref_flux)) if fitted_flux != 0 else -10

                            # Check if flux change direction matches protein change
                            if (log2_fc > 0 and flux_log2_change > 0.5) or \
                               (log2_fc < 0 and flux_log2_change < -0.5):
                                implemented = True
                                break

                if implemented:
                    output['protein_changes']['significant_implemented'][gene_id] = {
                        'log2_fc': log2_fc,
                        'fold_change': gc_data['fold_change'],
                        'reactions': gene_rxns
                    }
                else:
                    output['protein_changes']['significant_not_implemented'][gene_id] = {
                        'log2_fc': log2_fc,
                        'fold_change': gc_data['fold_change'],
                        'reactions': gene_rxns,
                        'reason': 'flux_could_not_match'
                    }

            # Analyze reaction matching
            for rxn_id, fc_data in fold_changes.items():
                expected_fc = fc_data['fold_change']

                if rxn_id not in fitted_fluxes or rxn_id not in reference_flux:
                    output['reaction_matching']['not_matched'].append({
                        'rxn_id': rxn_id,
                        'expected_fc': expected_fc,
                        'reason': 'missing_flux_data'
                    })
                    continue

                ref_flux = reference_flux[rxn_id]
                fitted_flux = fitted_fluxes[rxn_id]

                if abs(ref_flux) < 1e-9:
                    actual_fc = abs(fitted_flux / zero_flux) if abs(fitted_flux) > 1e-9 else 1.0
                else:
                    actual_fc = abs(fitted_flux / ref_flux) if abs(fitted_flux) > 1e-9 else 0.0

                # Check if actual fold change is within 50% of expected
                fc_ratio = actual_fc / expected_fc if expected_fc > 0 else 0
                if 0.5 <= fc_ratio <= 2.0:
                    output['reaction_matching']['matched'].append({
                        'rxn_id': rxn_id,
                        'expected_fc': expected_fc,
                        'actual_fc': actual_fc,
                        'contributing_genes': fc_data['contributing_genes']
                    })
                else:
                    output['reaction_matching']['not_matched'].append({
                        'rxn_id': rxn_id,
                        'expected_fc': expected_fc,
                        'actual_fc': actual_fc,
                        'contributing_genes': fc_data['contributing_genes'],
                        'reason': 'fold_change_mismatch'
                    })

        # Summary statistics
        output['statistics']['total_reactions_with_fold_change'] = len(fold_changes)
        output['statistics']['total_genes_with_fold_change'] = len(gene_fold_changes)
        output['statistics']['reactions_matched'] = len(output['reaction_matching']['matched'])
        output['statistics']['reactions_not_matched'] = len(output['reaction_matching']['not_matched'])
        output['statistics']['proteins_implemented'] = len(output['protein_changes']['significant_implemented'])
        output['statistics']['proteins_not_implemented'] = len(output['protein_changes']['significant_not_implemented'])

        logger.info(f"Flux fitting complete. Matched {output['statistics']['reactions_matched']}/{output['statistics']['total_reactions_with_fold_change']} reactions")
        logger.info(f"Proteins implemented: {output['statistics']['proteins_implemented']}, not implemented: {output['statistics']['proteins_not_implemented']}")

        return output

    def plot_expression_distributions(
        self,
        conditions: list = None,
        figsize: tuple = None,
        bins: int = 50,
        show_kde: bool = True,
        show_median: bool = True,
        show_wildtype_estimate: bool = True,
        title_suffix: str = "",
        color: str = 'skyblue',
        kde_color: str = 'red',
        save_path: str = None,
        return_stats: bool = False
    ) -> dict:
        """Plot distribution of expression values for each condition.

        Creates histogram plots with optional KDE curves showing the distribution
        of expression values across features for each condition. Useful for
        understanding data quality, identifying outliers, and estimating
        wildtype/baseline expression levels.

        Parameters
        ----------
        conditions : list, optional
            List of condition IDs to plot. If None, plots all conditions.
        figsize : tuple, optional
            Figure size as (width, height). If None, automatically calculated
            based on number of conditions.
        bins : int, optional
            Number of histogram bins. Default: 50
        show_kde : bool, optional
            If True, overlay a kernel density estimate curve. Default: True
        show_median : bool, optional
            If True, show vertical line at median value. Default: True
        show_wildtype_estimate : bool, optional
            If True, estimate and show wildtype value as the peak density
            (mode of the KDE). Default: True
        title_suffix : str, optional
            Suffix to add to subplot titles (e.g., "(Gene-level)"). Default: ""
        color : str, optional
            Histogram bar color. Default: 'skyblue'
        kde_color : str, optional
            KDE line color. Default: 'red'
        save_path : str, optional
            If provided, save the figure to this path. Default: None
        return_stats : bool, optional
            If True, return statistics dictionary along with figure. Default: False

        Returns
        -------
        dict
            Dictionary containing:
            - 'figure': matplotlib Figure object
            - 'axes': array of matplotlib Axes objects
            - 'statistics': dict of per-condition statistics (if return_stats=True)
              Each condition has: count, mean, median, std, min, max, q25, q75,
              q95, wildtype_estimate (peak density value)

        Examples
        --------
        >>> # Basic usage - plot all conditions
        >>> result = expression.plot_expression_distributions()
        >>> result['figure'].savefig('distributions.png')

        >>> # Plot specific conditions with statistics
        >>> result = expression.plot_expression_distributions(
        ...     conditions=['glucose', 'acetate'],
        ...     return_stats=True,
        ...     save_path='my_distributions.png'
        ... )
        >>> print(result['statistics']['glucose']['wildtype_estimate'])

        >>> # Customize appearance
        >>> result = expression.plot_expression_distributions(
        ...     bins=30,
        ...     color='lightcoral',
        ...     kde_color='blue',
        ...     title_suffix='(Reaction-level)'
        ... )

        Notes
        -----
        The wildtype estimate is computed as the x-value where the KDE reaches
        its maximum density. This represents the most common expression value
        and is often a good estimate of the baseline/wildtype level when most
        genes are unaffected.

        The function requires matplotlib and scipy to be installed.
        """
        try:
            import matplotlib.pyplot as plt
            from scipy import stats as scipy_stats
        except ImportError as e:
            raise ImportError(
                "matplotlib and scipy are required for plotting. "
                "Install with: pip install matplotlib scipy"
            ) from e

        # Determine conditions to plot
        if conditions is None:
            conditions = [c.id for c in self.conditions]

        # Validate conditions exist
        valid_conditions = []
        for cond in conditions:
            if cond in [c.id for c in self.conditions]:
                valid_conditions.append(cond)
            else:
                logger.warning(f"Condition '{cond}' not found in expression data, skipping")

        if not valid_conditions:
            raise ValueError("No valid conditions to plot")

        n_conditions = len(valid_conditions)

        # Calculate grid dimensions
        if n_conditions <= 3:
            n_rows, n_cols = 1, n_conditions
        elif n_conditions <= 6:
            n_rows, n_cols = 2, 3
        elif n_conditions <= 9:
            n_rows, n_cols = 3, 3
        else:
            n_cols = 4
            n_rows = (n_conditions + n_cols - 1) // n_cols

        # Set figure size if not provided
        if figsize is None:
            figsize = (5 * n_cols, 4 * n_rows)

        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)

        # Flatten axes for easy iteration
        if n_conditions == 1:
            axes = np.array([axes])
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

        # Store statistics for each condition
        statistics = {}

        for idx, condition in enumerate(valid_conditions):
            ax = axes[idx]

            # Get data for this condition
            if condition not in self._data.columns:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12)
                ax.set_title(f'{condition.capitalize()} {title_suffix}'.strip())
                continue

            data = self._data[condition].dropna()

            if len(data) < 2:
                ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center',
                       transform=ax.transAxes, fontsize=12)
                ax.set_title(f'{condition.capitalize()} {title_suffix}'.strip())
                continue

            # Calculate statistics
            stats = {
                'count': len(data),
                'mean': float(data.mean()),
                'median': float(data.median()),
                'std': float(data.std()),
                'min': float(data.min()),
                'max': float(data.max()),
                'q25': float(data.quantile(0.25)),
                'q75': float(data.quantile(0.75)),
                'q95': float(data.quantile(0.95))
            }

            # Create histogram
            ax.hist(data, bins=bins, alpha=0.6, color=color,
                   edgecolor='black', density=True)

            # Calculate KDE and find peak density (wildtype estimate)
            kde = scipy_stats.gaussian_kde(data)
            x_range = np.linspace(data.min(), data.max(), 1000)
            kde_values = kde(x_range)
            peak_density_idx = np.argmax(kde_values)
            wildtype_estimate = float(x_range[peak_density_idx])
            stats['wildtype_estimate'] = wildtype_estimate

            # Plot KDE curve
            if show_kde:
                ax.plot(x_range, kde_values, color=kde_color, linewidth=2, label='KDE')

            # Add vertical lines for statistics
            if show_median:
                ax.axvline(stats['median'], color='green', linestyle='--',
                          linewidth=2, label=f"Median: {stats['median']:.3f}")

            if show_wildtype_estimate:
                ax.axvline(wildtype_estimate, color='orange', linestyle='--',
                          linewidth=2, label=f"WT est: {wildtype_estimate:.3f}")

            # Formatting
            ax.set_xlabel('Expression Value', fontsize=10)
            ax.set_ylabel('Density', fontsize=10)
            ax.set_title(f'{condition.capitalize()} {title_suffix}'.strip(),
                        fontsize=12, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

            statistics[condition] = stats

        # Hide unused subplots
        for idx in range(n_conditions, len(axes)):
            axes[idx].axis('off')

        plt.tight_layout()

        # Save if path provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Distribution plot saved to {save_path}")

        # Build return dictionary
        result = {
            'figure': fig,
            'axes': axes[:n_conditions]
        }

        if return_stats:
            result['statistics'] = statistics

        return result

    def print_distribution_statistics(self, conditions: list = None) -> pd.DataFrame:
        """Print and return distribution statistics for expression data.

        Computes and displays summary statistics for expression values
        across all features for each condition.

        Parameters
        ----------
        conditions : list, optional
            List of condition IDs to analyze. If None, analyzes all conditions.

        Returns
        -------
        pd.DataFrame
            DataFrame with statistics for each condition including:
            count, mean, median, std, min, max, q25, q75, q95

        Examples
        --------
        >>> stats_df = expression.print_distribution_statistics()
        >>> print(stats_df)

        >>> # Analyze specific conditions
        >>> stats_df = expression.print_distribution_statistics(['glucose', 'acetate'])
        """
        if conditions is None:
            conditions = [c.id for c in self.conditions]

        records = []

        print("Expression Distribution Statistics:")
        print("=" * 80)

        for condition in conditions:
            if condition not in self._data.columns:
                logger.warning(f"Condition '{condition}' not found, skipping")
                continue

            data = self._data[condition].dropna()

            if len(data) == 0:
                logger.warning(f"No data for condition '{condition}', skipping")
                continue

            stats = {
                'condition': condition,
                'count': len(data),
                'mean': data.mean(),
                'median': data.median(),
                'std': data.std(),
                'min': data.min(),
                'max': data.max(),
                'q25': data.quantile(0.25),
                'q75': data.quantile(0.75),
                'q95': data.quantile(0.95)
            }
            records.append(stats)

            print(f"\n{condition.upper()}:")
            print(f"  Count:      {stats['count']}")
            print(f"  Mean:       {stats['mean']:.4f}")
            print(f"  Median:     {stats['median']:.4f}")
            print(f"  Std Dev:    {stats['std']:.4f}")
            print(f"  Min:        {stats['min']:.4f}")
            print(f"  Max:        {stats['max']:.4f}")
            print(f"  25th %ile:  {stats['q25']:.4f}")
            print(f"  75th %ile:  {stats['q75']:.4f}")
            print(f"  95th %ile:  {stats['q95']:.4f}")

        print("\n" + "=" * 80)

        return pd.DataFrame.from_records(records)
