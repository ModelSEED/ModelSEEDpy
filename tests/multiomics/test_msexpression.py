# -*- coding: utf-8 -*-
"""
Comprehensive test suite for MSExpression refactoring.

This test module covers the DataFrame-based refactoring of the MSExpression class,
ensuring all functionality works correctly with the new pandas-based data storage.
"""

import unittest
import pandas as pd
import numpy as np
import tempfile
import os
from modelseedpy.multiomics.msexpression import (
    MSExpression,
    MSExpressionFeature,
    MSCondition,
    compute_gene_score
)
from modelseedpy.core.msgenome import MSGenome, MSFeature
from cobra import Model, Reaction, Gene


class TestDataFrameInitialization(unittest.TestCase):
    """TC-1: Test DataFrame initialization and structure."""

    def test_msexpression_init_creates_empty_dataframe(self):
        """Test that MSExpression.__init__() creates empty DataFrame with named index."""
        expr = MSExpression("RelativeAbundance")
        self.assertIsInstance(expr._data, pd.DataFrame)
        self.assertEqual(expr._data.index.name, 'feature_id')
        self.assertEqual(len(expr._data), 0)
        self.assertEqual(len(expr._data.columns), 0)

    def test_dataframe_index_name_is_feature_id(self):
        """Test that DataFrame index is named 'feature_id'."""
        expr = MSExpression("FPKM")
        self.assertEqual(expr._data.index.name, 'feature_id')

    def test_dataframe_is_private_attribute(self):
        """Test that _data is a private attribute (single underscore)."""
        expr = MSExpression("TPM")
        self.assertTrue(hasattr(expr, '_data'))
        self.assertIsInstance(expr._data, pd.DataFrame)


class TestDataLoadingFromDataFrame(unittest.TestCase):
    """TC-2: Test data loading from DataFrame."""

    def setUp(self):
        """Create test data for loading tests."""
        self.test_df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [10.5, 8.2, 15.3],
            'condition2': [20.3, 15.7, 12.1],
            'condition3': [5.1, 3.4, 7.8]
        })

        # Create a genome with features
        self.genome = MSGenome()
        self.genome.features.append(MSFeature('gene1', ''))
        self.genome.features.append(MSFeature('gene2', ''))
        self.genome.features.append(MSFeature('gene3', ''))

    def test_from_dataframe_basic_loading(self):
        """Test basic loading from DataFrame."""
        expr = MSExpression.from_dataframe(
            self.test_df,
            genome=self.genome,
            create_missing_features=False,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        self.assertEqual(len(expr.features), 3)
        self.assertEqual(len(expr.conditions), 3)
        self.assertEqual(expr._data.shape, (3, 3))

    def test_from_dataframe_bulk_loading(self):
        """Test that data is loaded in bulk, not row-by-row."""
        expr = MSExpression.from_dataframe(
            self.test_df,
            genome=self.genome,
            id_column='gene_id'
        )

        # Verify all values are loaded correctly
        self.assertAlmostEqual(expr.get_value('gene1', 'condition1'), 10.5)
        self.assertAlmostEqual(expr.get_value('gene2', 'condition2'), 15.7)
        self.assertAlmostEqual(expr.get_value('gene3', 'condition3'), 7.8)

    def test_from_dataframe_with_ignore_columns(self):
        """Test loading with ignore_columns parameter."""
        df = self.test_df.copy()
        df['metadata'] = ['A', 'B', 'C']

        expr = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            ignore_columns=['metadata'],
            id_column='gene_id'
        )

        self.assertEqual(len(expr.conditions), 3)
        self.assertNotIn('metadata', [c.id for c in expr.conditions])

    def test_from_dataframe_with_description_column(self):
        """Test loading with description_column parameter."""
        df = self.test_df.copy()
        df['description'] = ['Desc1', 'Desc2', 'Desc3']

        expr = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            description_column='description',
            id_column='gene_id'
        )

        self.assertEqual(len(expr.conditions), 3)
        self.assertNotIn('description', [c.id for c in expr.conditions])

    def test_from_dataframe_default_id_column(self):
        """Test that first column is used as ID if id_column not specified."""
        expr = MSExpression.from_dataframe(
            self.test_df,
            genome=self.genome
        )

        self.assertEqual(len(expr.features), 3)
        self.assertIn('gene1', expr.features)

    def test_from_dataframe_creates_missing_features(self):
        """Test create_missing_features parameter."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene_new'],
            'condition1': [10.5, 8.2, 5.0]
        })

        expr = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            create_missing_features=True,
            id_column='gene_id'
        )

        # gene_new should be created
        self.assertEqual(len(self.genome.features), 4)

    def test_from_dataframe_handles_nan_values(self):
        """Test that NaN values are handled correctly."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],
            'condition1': [10.5, np.nan],
            'condition2': [np.nan, 15.7]
        })

        expr = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id'
        )

        # NaN should be converted to None in API
        self.assertIsNone(expr.get_value('gene2', 'condition1'))
        self.assertIsNone(expr.get_value('gene1', 'condition2'))


class TestFeatureAndConditionManagement(unittest.TestCase):
    """TC-3: Test feature and condition management."""

    def setUp(self):
        """Create test expression object."""
        self.genome = MSGenome()
        self.genome.features.append(MSFeature('gene1', ''))
        self.genome.features.append(MSFeature('gene2', ''))

        self.expr = MSExpression('RelativeAbundance')
        self.expr.object = self.genome

    def test_add_feature_creates_msexpressionfeature(self):
        """Test that add_feature creates MSExpressionFeature without values dict."""
        feature = self.expr.add_feature('gene1')

        self.assertIsInstance(feature, MSExpressionFeature)
        self.assertEqual(feature.id, 'gene1')
        self.assertFalse(hasattr(feature, 'values'))

    def test_add_value_updates_dataframe(self):
        """Test that MSExpressionFeature.add_value updates parent DataFrame."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        feature.add_value(condition, 10.5)

        self.assertEqual(self.expr._data.loc['gene1', 'cond1'], 10.5)

    def test_add_value_creates_row_if_missing(self):
        """Test that add_value creates feature row if missing."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        feature.add_value(condition, 10.5)

        self.assertIn('gene1', self.expr._data.index)

    def test_add_value_creates_column_if_missing(self):
        """Test that add_value creates condition column if missing."""
        feature = self.expr.add_feature('gene1')
        condition = MSCondition('cond1', self.expr)

        feature.add_value(condition, 10.5)

        self.assertIn('cond1', self.expr._data.columns)

    def test_get_value_from_dataframe(self):
        """Test that get_value retrieves from DataFrame."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        feature.add_value(condition, 10.5)

        value = feature.get_value(condition)
        self.assertAlmostEqual(value, 10.5)

    def test_get_value_converts_nan_to_none(self):
        """Test that get_value converts NaN to None."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        self.expr._data.loc['gene1', 'cond1'] = np.nan

        value = feature.get_value(condition)
        self.assertIsNone(value)

    def test_get_value_with_string_condition(self):
        """Test get_value with condition ID string."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        feature.add_value(condition, 10.5)

        value = feature.get_value('cond1')
        self.assertAlmostEqual(value, 10.5)

    def test_msexpression_get_value(self):
        """Test MSExpression.get_value method."""
        condition = MSCondition('cond1', self.expr)
        self.expr.conditions.append(condition)

        feature = self.expr.add_feature('gene1')
        feature.add_value(condition, 10.5)

        value = self.expr.get_value('gene1', 'cond1')
        self.assertAlmostEqual(value, 10.5)


class TestStatisticalMethods(unittest.TestCase):
    """TC-4: Test statistical methods using pandas operations."""

    def setUp(self):
        """Create test expression with data."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3', 'gene4'],
            'condition1': [10.0, 20.0, 30.0, 40.0],
            'condition2': [5.0, 10.0, 15.0, 20.0]
        })

        genome = MSGenome()
        for i in range(1, 5):
            genome.features.append(MSFeature(f'gene{i}', ''))

        self.expr = MSExpression.from_dataframe(
            df,
            genome=genome,
            id_column='gene_id'
        )

    def test_condition_average_value(self):
        """Test MSCondition.average_value using pandas."""
        cond1 = self.expr.conditions.get_by_id('condition1')
        avg = cond1.average_value()
        self.assertAlmostEqual(avg, 25.0)

    def test_condition_lowest_value(self):
        """Test MSCondition.lowest_value using pandas."""
        cond1 = self.expr.conditions.get_by_id('condition1')
        lowest = cond1.lowest_value()
        self.assertAlmostEqual(lowest, 10.0)

    def test_condition_highest_value(self):
        """Test MSCondition.highest_value using pandas."""
        cond1 = self.expr.conditions.get_by_id('condition1')
        highest = cond1.highest_value()
        self.assertAlmostEqual(highest, 40.0)

    def test_condition_sum_value(self):
        """Test MSCondition.sum_value using pandas."""
        cond1 = self.expr.conditions.get_by_id('condition1')
        total = cond1.sum_value()
        self.assertAlmostEqual(total, 100.0)

    def test_condition_value_at_zscore(self):
        """Test MSCondition.value_at_zscore calculation."""
        cond1 = self.expr.conditions.get_by_id('condition1')
        # Mean = 25, StdDev = 11.18 (approximately)
        value_at_z1 = cond1.value_at_zscore(1.0)
        self.assertGreater(value_at_z1, 25.0)
        self.assertLess(value_at_z1, 40.0)

    def test_statistical_methods_with_nan_values(self):
        """Test that statistical methods handle NaN correctly."""
        self.expr._data.loc['gene1', 'condition1'] = np.nan

        cond1 = self.expr.conditions.get_by_id('condition1')
        avg = cond1.average_value()

        # Should compute average of remaining 3 values: (20 + 30 + 40) / 3 = 30
        self.assertAlmostEqual(avg, 30.0)

    def test_statistical_methods_empty_condition(self):
        """Test statistical methods with empty condition."""
        # Add empty condition
        self.expr._data['empty_cond'] = np.nan
        self.expr.conditions.append(MSCondition('empty_cond', self.expr))

        cond = self.expr.conditions.get_by_id('empty_cond')
        self.assertIsNone(cond.average_value())
        self.assertIsNone(cond.lowest_value())
        self.assertIsNone(cond.highest_value())
        self.assertEqual(cond.sum_value(), 0.0)


class TestValueRetrievalAndManipulation(unittest.TestCase):
    """TC-5: Test value retrieval and manipulation."""

    def setUp(self):
        """Create test expression."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],
            'condition1': [10.0, 20.0],
            'condition2': [5.0, 10.0]
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))
        genome.features.append(MSFeature('gene2', ''))

        self.expr = MSExpression.from_dataframe(
            df,
            genome=genome,
            id_column='gene_id'
        )

    def test_get_value_with_feature_string(self):
        """Test get_value with feature ID string."""
        value = self.expr.get_value('gene1', 'condition1')
        self.assertAlmostEqual(value, 10.0)

    def test_get_value_with_condition_string(self):
        """Test get_value with condition ID string."""
        feature = self.expr.features.get_by_id('gene1')
        value = feature.get_value('condition1')
        self.assertAlmostEqual(value, 10.0)

    def test_get_value_missing_feature(self):
        """Test get_value with missing feature."""
        value = self.expr.get_value('nonexistent', 'condition1')
        self.assertIsNone(value)

    def test_get_value_missing_condition(self):
        """Test get_value with missing condition."""
        feature = self.expr.features.get_by_id('gene1')
        value = feature.get_value('nonexistent')
        self.assertIsNone(value)

    def test_add_value_multiple_times(self):
        """Test that add_value can update existing values."""
        feature = self.expr.features.get_by_id('gene1')
        condition = self.expr.conditions.get_by_id('condition1')

        feature.add_value(condition, 100.0)
        value = feature.get_value(condition)
        self.assertAlmostEqual(value, 100.0)


class TestGPRIntegration(unittest.TestCase):
    """TC-6: Test GPR integration and reaction expression building."""

    def setUp(self):
        """Create test model and expression."""
        # Create a simple model
        self.model = Model('test_model')

        # Add genes
        gene1 = Gene('gene1')
        gene2 = Gene('gene2')
        gene3 = Gene('gene3')
        self.model.genes.extend([gene1, gene2, gene3])

        # Add reactions with GPR rules
        rxn1 = Reaction('rxn1')
        rxn1.gene_reaction_rule = 'gene1'

        rxn2 = Reaction('rxn2')
        rxn2.gene_reaction_rule = 'gene1 and gene2'

        rxn3 = Reaction('rxn3')
        rxn3.gene_reaction_rule = 'gene1 or gene2'

        self.model.add_reactions([rxn1, rxn2, rxn3])

        # Create genome
        self.genome = MSGenome()
        self.genome.features.append(MSFeature('gene1', ''))
        self.genome.features.append(MSFeature('gene2', ''))
        self.genome.features.append(MSFeature('gene3', ''))

        # Create gene expression
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [10.0, 20.0, 30.0]
        })

        self.gene_expr = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id'
        )

    def test_build_reaction_expression_basic(self):
        """Test building reaction expression from gene expression."""
        rxn_expr = self.gene_expr.build_reaction_expression(self.model, default=0.0)

        self.assertEqual(len(rxn_expr.features), 3)
        self.assertEqual(len(rxn_expr.conditions), 1)

    def test_build_reaction_expression_or_rule(self):
        """Test OR rule in GPR (sum of gene values)."""
        rxn_expr = self.gene_expr.build_reaction_expression(self.model, default=0.0)

        # rxn3 has 'gene1 or gene2', should be 10 + 20 = 30
        value = rxn_expr.get_value('rxn3', 'condition1')
        self.assertAlmostEqual(value, 30.0)

    def test_build_reaction_expression_and_rule(self):
        """Test AND rule in GPR (min of gene values)."""
        rxn_expr = self.gene_expr.build_reaction_expression(self.model, default=0.0)

        # rxn2 has 'gene1 and gene2', should be min(10, 20) = 10
        value = rxn_expr.get_value('rxn2', 'condition1')
        self.assertAlmostEqual(value, 10.0)

    def test_build_reaction_expression_single_gene(self):
        """Test single gene GPR rule."""
        rxn_expr = self.gene_expr.build_reaction_expression(self.model, default=0.0)

        # rxn1 has 'gene1', should be 10
        value = rxn_expr.get_value('rxn1', 'condition1')
        self.assertAlmostEqual(value, 10.0)

    def test_build_reaction_expression_uses_dataframe(self):
        """Test that build_reaction_expression accesses _data DataFrame."""
        rxn_expr = self.gene_expr.build_reaction_expression(self.model, default=0.0)

        # Verify the reaction expression has a DataFrame
        self.assertIsInstance(rxn_expr._data, pd.DataFrame)
        self.assertGreater(len(rxn_expr._data), 0)


class TestDataExport(unittest.TestCase):
    """TC-7: Test data export methods."""

    def setUp(self):
        """Create test expression."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [10.0, 20.0, 30.0],
            'condition2': [5.0, 10.0, 15.0]
        })

        genome = MSGenome()
        for i in range(1, 4):
            genome.features.append(MSFeature(f'gene{i}', ''))

        self.expr = MSExpression.from_dataframe(
            df,
            genome=genome,
            id_column='gene_id'
        )

    def test_get_dataframe_default_format(self):
        """Test get_dataframe returns index format by default."""
        df = self.expr.get_dataframe()

        self.assertEqual(df.index.name, 'feature_id')
        self.assertIn('gene1', df.index)
        self.assertIn('condition1', df.columns)

    def test_get_dataframe_reset_index(self):
        """Test get_dataframe with reset_index=True."""
        df = self.expr.get_dataframe(reset_index=True)

        self.assertIn('feature_id', df.columns)
        self.assertIn('condition1', df.columns)

    def test_get_dataframe_returns_copy(self):
        """Test that get_dataframe returns a copy, not reference."""
        df = self.expr.get_dataframe()
        df.loc['gene1', 'condition1'] = 999.0

        # Original should be unchanged
        original_value = self.expr.get_value('gene1', 'condition1')
        self.assertAlmostEqual(original_value, 10.0)

    def test_get_dataframe_preserves_column_order(self):
        """Test that column order is preserved."""
        df = self.expr.get_dataframe()
        columns = list(df.columns)

        self.assertEqual(columns, ['condition1', 'condition2'])


class TestEdgeCasesAndErrorHandling(unittest.TestCase):
    """TC-8: Test edge cases and error handling."""

    def test_empty_dataframe_loading(self):
        """Test loading empty DataFrame."""
        df = pd.DataFrame(columns=['gene_id', 'condition1'])
        genome = MSGenome()

        expr = MSExpression.from_dataframe(df, genome=genome, id_column='gene_id')

        self.assertEqual(len(expr.features), 0)
        self.assertEqual(len(expr._data), 0)

    def test_single_feature_single_condition(self):
        """Test with single feature and single condition."""
        df = pd.DataFrame({
            'gene_id': ['gene1'],
            'condition1': [10.0]
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))

        expr = MSExpression.from_dataframe(df, genome=genome, id_column='gene_id')

        self.assertEqual(len(expr.features), 1)
        self.assertEqual(len(expr.conditions), 1)
        self.assertAlmostEqual(expr.get_value('gene1', 'condition1'), 10.0)

    def test_all_nan_condition(self):
        """Test condition with all NaN values."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],
            'condition1': [np.nan, np.nan]
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))
        genome.features.append(MSFeature('gene2', ''))

        expr = MSExpression.from_dataframe(df, genome=genome, id_column='gene_id')

        cond = expr.conditions.get_by_id('condition1')
        self.assertIsNone(cond.average_value())

    def test_mixed_numeric_non_numeric_values(self):
        """Test that non-numeric values are converted to NaN."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],
            'condition1': ['10.0', 'invalid']
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))
        genome.features.append(MSFeature('gene2', ''))

        expr = MSExpression.from_dataframe(df, genome=genome, id_column='gene_id')

        value1 = expr.get_value('gene1', 'condition1')
        value2 = expr.get_value('gene2', 'condition1')

        self.assertAlmostEqual(value1, 10.0)
        self.assertIsNone(value2)

    def test_duplicate_feature_ids(self):
        """Test handling of duplicate feature IDs."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene1', 'gene2'],
            'condition1': [10.0, 20.0, 30.0]
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))
        genome.features.append(MSFeature('gene2', ''))

        expr = MSExpression.from_dataframe(df, genome=genome, id_column='gene_id')

        # Should handle duplicates gracefully (last value wins in pandas)
        self.assertIsNotNone(expr.get_value('gene1', 'condition1'))

    def test_nonexistent_id_column(self):
        """Test behavior with nonexistent id_column."""
        df = pd.DataFrame({
            'gene_id': ['gene1'],
            'condition1': [10.0]
        })

        genome = MSGenome()

        with self.assertRaises(KeyError):
            MSExpression.from_dataframe(
                df,
                genome=genome,
                id_column='nonexistent'
            )

    def test_feature_not_in_genome(self):
        """Test behavior when feature not found in genome."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene_missing'],
            'condition1': [10.0, 20.0]
        })

        genome = MSGenome()
        genome.features.append(MSFeature('gene1', ''))

        expr = MSExpression.from_dataframe(
            df,
            genome=genome,
            create_missing_features=False,
            id_column='gene_id'
        )

        # Only gene1 should be loaded
        self.assertEqual(len(expr.features), 1)
        self.assertNotIn('gene_missing', expr.features)


class TestTypeAnnotations(unittest.TestCase):
    """TC-9: Test type annotations and signatures."""

    def test_from_dataframe_type_annotations(self):
        """Test that from_dataframe has proper type annotations."""
        import inspect
        sig = inspect.signature(MSExpression.from_dataframe)

        # Check key parameters have annotations
        self.assertIn('df', sig.parameters)
        self.assertIn('genome', sig.parameters)
        # Check return annotation exists
        self.assertIsNot(sig.return_annotation, inspect.Signature.empty)

    def test_get_value_type_annotations(self):
        """Test that get_value has proper type annotations."""
        import inspect
        sig = inspect.signature(MSExpression.get_value)

        self.assertIn('feature', sig.parameters)
        self.assertIn('condition', sig.parameters)

    def test_statistical_methods_type_annotations(self):
        """Test that statistical methods have type annotations."""
        import inspect

        methods = [
            MSCondition.average_value,
            MSCondition.lowest_value,
            MSCondition.highest_value,
            MSCondition.sum_value
        ]

        for method in methods:
            sig = inspect.signature(method)
            # Should have return annotation
            self.assertIsNotNone(sig.return_annotation)


class TestComputeGeneScore(unittest.TestCase):
    """TC-10: Test compute_gene_score helper function."""

    def test_compute_gene_score_single_gene(self):
        """Test compute_gene_score with single gene."""
        from cobra.core.gene import parse_gpr

        gpr = parse_gpr('gene1')
        values = {'gene1': 10.0}

        score = compute_gene_score(gpr, values, default=0.0)
        self.assertAlmostEqual(score, 10.0)

    def test_compute_gene_score_or_operation(self):
        """Test compute_gene_score with OR operation."""
        from cobra.core.gene import parse_gpr

        gpr = parse_gpr('gene1 or gene2')
        values = {'gene1': 10.0, 'gene2': 20.0}

        score = compute_gene_score(gpr, values, default=0.0)
        self.assertAlmostEqual(score, 30.0)  # Sum

    def test_compute_gene_score_and_operation(self):
        """Test compute_gene_score with AND operation."""
        from cobra.core.gene import parse_gpr

        gpr = parse_gpr('gene1 and gene2')
        values = {'gene1': 10.0, 'gene2': 20.0}

        score = compute_gene_score(gpr, values, default=0.0)
        self.assertAlmostEqual(score, 10.0)  # Min

    def test_compute_gene_score_missing_gene_uses_default(self):
        """Test compute_gene_score with missing gene uses default."""
        from cobra.core.gene import parse_gpr

        gpr = parse_gpr('gene1')
        values = {}

        score = compute_gene_score(gpr, values, default=5.0)
        self.assertAlmostEqual(score, 5.0)


class TestFitModelFluxToData(unittest.TestCase):
    """Test fit_model_flux_to_data function."""

    def setUp(self):
        """Create test models and expression data for testing."""
        # Create a simple model with reactions and genes
        self.model = Model('test_model')

        # Add genes
        gene1 = Gene('gene1')
        gene2 = Gene('gene2')
        gene3 = Gene('gene3')
        self.model.genes.extend([gene1, gene2, gene3])

        # Add reactions with GPR rules
        from cobra import Metabolite
        met_a = Metabolite('met_a', compartment='c')
        met_b = Metabolite('met_b', compartment='c')
        met_c = Metabolite('met_c', compartment='c')

        rxn1 = Reaction('rxn1')
        rxn1.gene_reaction_rule = 'gene1'
        rxn1.add_metabolites({met_a: -1, met_b: 1})
        rxn1.bounds = (-10, 10)

        rxn2 = Reaction('rxn2')
        rxn2.gene_reaction_rule = 'gene2'
        rxn2.add_metabolites({met_b: -1, met_c: 1})
        rxn2.bounds = (-10, 10)

        rxn3 = Reaction('rxn3')
        rxn3.gene_reaction_rule = 'gene3'
        rxn3.add_metabolites({met_c: -1, met_a: 1})
        rxn3.bounds = (-10, 10)

        # Add biomass reaction
        biomass = Reaction('biomass')
        biomass.add_metabolites({met_b: -1})
        biomass.bounds = (0, 10)

        self.model.add_reactions([rxn1, rxn2, rxn3, biomass])
        self.model.objective = 'biomass'

        # Create genome for gene expression
        self.genome = MSGenome()
        self.genome.features.append(MSFeature('gene1', ''))
        self.genome.features.append(MSFeature('gene2', ''))
        self.genome.features.append(MSFeature('gene3', ''))

    def test_fit_model_flux_basic(self):
        """Task 5.3: Test basic functionality with default thresholds."""
        # Create gene expression data with high expression for all genes
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2]  # All above default activation threshold
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        # Create MSModelUtil (mock basic functionality)
        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Call fit_model_flux_to_data
        solution = expression.fit_model_flux_to_data(
            model=model_util,
            condition='condition1'
        )

        # Verify solution is optimal
        self.assertEqual(solution.status, 'optimal')
        self.assertIsNotNone(solution.objective_value)

    def test_fit_model_flux_invalid_thresholds(self):
        """Task 5.7: Test invalid threshold validation."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Test with activation <= deactivation
        with self.assertRaises(ValueError) as context:
            expression.fit_model_flux_to_data(
                model=model_util,
                condition='condition1',
                activation_threshold=0.001,
                deactivation_threshold=0.002
            )

        self.assertIn('must be greater than', str(context.exception))

    def test_fit_model_flux_missing_condition_error(self):
        """Task 5.6: Test missing condition error."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Test with non-existent condition
        with self.assertRaises(ValueError) as context:
            expression.fit_model_flux_to_data(
                model=model_util,
                condition='nonexistent_condition'
            )

        self.assertIn('not found in expression data', str(context.exception))
        self.assertIn('Available conditions', str(context.exception))

    def test_fit_model_flux_genome_to_model_conversion(self):
        """Task 5.4: Test genome to model conversion."""
        # Create gene expression (genome-level)
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        # Verify it's genome-level
        self.assertIsInstance(expression.object, MSGenome)

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Call fit_model_flux_to_data - should convert to reaction level
        solution = expression.fit_model_flux_to_data(
            model=model_util,
            condition='condition1'
        )

        # Verify optimization succeeded
        self.assertEqual(solution.status, 'optimal')

    def test_fit_model_flux_model_mismatch_error(self):
        """Task 5.5: Test model mismatch error."""
        # Create reaction-level expression for model A
        rxn_expr = MSExpression('RelativeAbundance')
        rxn_expr.object = self.model
        for rxn in self.model.reactions:
            rxn_expr.add_feature(rxn.id)

        cond = MSCondition('condition1', rxn_expr)
        rxn_expr.conditions.append(cond)

        # Add some values
        for rxn in self.model.reactions:
            rxn_expr.features.get_by_id(rxn.id).add_value(cond, 0.5)

        # Create a different model (Model B)
        model_b = Model('test_model_b')
        gene_b = Gene('gene_b')
        model_b.genes.append(gene_b)

        from cobra import Metabolite
        met_b = Metabolite('met_b', compartment='c')
        rxn_b = Reaction('rxn_b')
        rxn_b.gene_reaction_rule = 'gene_b'
        rxn_b.add_metabolites({met_b: -1})
        model_b.add_reactions([rxn_b])

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util_b = MSModelUtil.get(model_b)

        # Try to fit with mismatched model
        with self.assertRaises(ValueError) as context:
            rxn_expr.fit_model_flux_to_data(
                model=model_util_b,
                condition='condition1'
            )

        self.assertIn('Models must match', str(context.exception))

    def test_fit_model_flux_type_transformations(self):
        """Task 5.8: Test type transformations to RelativeAbundance."""
        # Test FPKM transformation
        df_fpkm = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [100.0, 200.0, 300.0]  # FPKM values
        })

        expression_fpkm = MSExpression.from_dataframe(
            df_fpkm,
            genome=self.genome,
            id_column='gene_id',
            type='FPKM'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Should transform FPKM to RelativeAbundance
        solution = expression_fpkm.fit_model_flux_to_data(
            model=model_util,
            condition='condition1',
            activation_threshold=0.1,  # Lower threshold for transformed data
            deactivation_threshold=0.01
        )

        self.assertEqual(solution.status, 'optimal')

        # Verify original data is unchanged
        self.assertEqual(expression_fpkm.type, 'FPKM')
        self.assertAlmostEqual(expression_fpkm.get_value('gene1', 'condition1'), 100.0)

    def test_fit_model_flux_custom_thresholds(self):
        """Task 5.9: Test custom activation/deactivation thresholds."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.05, 0.001]  # High, medium, low expression
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Use custom thresholds
        solution = expression.fit_model_flux_to_data(
            model=model_util,
            condition='condition1',
            activation_threshold=0.1,  # gene1 should be activated
            deactivation_threshold=0.01  # gene3 should be deactivated
        )

        self.assertEqual(solution.status, 'optimal')

    def test_fit_model_flux_no_activated_reactions(self):
        """Task 5.10: Test error when no reactions above threshold."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.0001, 0.0001, 0.0001]  # All very low
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Use very high threshold - no reactions should be activated
        with self.assertRaises(ValueError) as context:
            expression.fit_model_flux_to_data(
                model=model_util,
                condition='condition1',
                activation_threshold=10.0
            )

        self.assertIn('No reactions have expression values above', str(context.exception))

    def test_fit_model_flux_all_reactions_activated(self):
        """Task 5.11: Test with very low threshold - all reactions activated."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Use very low thresholds
        solution = expression.fit_model_flux_to_data(
            model=model_util,
            condition='condition1',
            activation_threshold=0.00001,
            deactivation_threshold=0.000001
        )

        # Should still succeed with all reactions activated
        self.assertEqual(solution.status, 'optimal')

    def test_fit_model_flux_missing_expression_data(self):
        """Task 5.12: Test that reactions without expression data are ignored."""
        # Only provide expression for gene1 and gene2, not gene3
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2'],  # gene3 missing
            'condition1': [0.5, 0.3]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance',
            create_missing_features=False
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Should handle missing data gracefully
        solution = expression.fit_model_flux_to_data(
            model=model_util,
            condition='condition1'
        )

        # Optimization should still work
        self.assertEqual(solution.status, 'optimal')

    def test_fit_model_flux_multiple_conditions(self):
        """Task 5.14: Test multiple conditions - verify each returns different fluxes."""
        df = pd.DataFrame({
            'gene_id': ['gene1', 'gene2', 'gene3'],
            'condition1': [0.5, 0.3, 0.2],
            'condition2': [0.2, 0.5, 0.3],
            'condition3': [0.3, 0.2, 0.5]
        })

        expression = MSExpression.from_dataframe(
            df,
            genome=self.genome,
            id_column='gene_id',
            type='RelativeAbundance'
        )

        from modelseedpy.core.msmodelutl import MSModelUtil
        model_util = MSModelUtil.get(self.model)

        # Run for multiple conditions
        solutions = {}
        for cond in ['condition1', 'condition2', 'condition3']:
            solutions[cond] = expression.fit_model_flux_to_data(
                model=model_util,
                condition=cond
            )

        # Verify all succeeded
        for cond, sol in solutions.items():
            self.assertEqual(sol.status, 'optimal')

        # Verify original data unchanged
        self.assertEqual(expression.type, 'RelativeAbundance')
        self.assertAlmostEqual(expression.get_value('gene1', 'condition1'), 0.5)


if __name__ == '__main__':
    unittest.main()
