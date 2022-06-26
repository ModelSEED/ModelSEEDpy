import cobra
from tests.test_data.mock_data import mock_model_ecoli_core
from modelseedpy.fbapkg.revbinpkg import RevBinPkg


def test_init():
    rev_bin = RevBinPkg(model=mock_model_ecoli_core(True))
    assert type(rev_bin.model) is cobra.core.Model
    assert type(rev_bin.name) is str
    assert type(rev_bin.variable_types) is dict
    assert type(rev_bin.constraint_types) is dict
    assert 'revbin' in list(rev_bin.variables.keys())
    assert 'revbinR' in list(rev_bin.constraints.keys())
    assert 'revbinF' in list(rev_bin.constraints.keys())


def test_build_package():
    # build_variable parameters
    lower_bound = 0
    upper_bound = 1
    var_type = "binary"
    constraint_types = {"F": [None, 0],'R': [None, 1000]}

    # execute the function
    rev_bin = RevBinPkg(model=mock_model_ecoli_core(True))
    rev_bin.build_package()

    # execute the function and assert results of the function
    for reaction in rev_bin.model.reactions:
        rev_bin_var = rev_bin.variables["revbin"][reaction.id]
        assert rev_bin_var
        assert rev_bin_var.ub == upper_bound
        assert rev_bin_var.lb == lower_bound
        assert rev_bin_var.type == var_type

        for constraint_type in constraint_types:
            rev_bin_cons = rev_bin.constraints['revbin{}'.format(constraint_type)][reaction.id]
            assert rev_bin_cons
            assert rev_bin_cons.lb == constraint_types[constraint_type][0]
            assert rev_bin_cons.ub == constraint_types[constraint_type][1]
            assert rev_bin_cons.name == '{}_revbin{}'.format(reaction.id, constraint_type)
