# import cobrakbase
import os
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'KCYWCM5Q3ZFAGQKDG3ESMJXGHQPK7UMN'
kbase = cobrakbase.KBaseAPI(token)

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
from modelseedpy.fbapkg.revbinpkg import RevBinPkg
revbin = RevBinPkg(model = model)    

# ------------------------ test the RevBin Package ---------------------------------------
        
def test_init():
    
    # assert results of the model 
    assert type(revbin.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(revbin.name) is str
    assert type(revbin.variable_types) is dict
    assert type(revbin.constraint_types) is dict
    assert 'revbin' in list(revbin.variables.keys())
    assert 'revbinR' in list(revbin.constraints.keys())
    assert 'revbinF' in list(revbin.constraints.keys())
    

def test_build_package():
    # build_variable parameters
    lower_bound = 0
    upper_bound = 1
    var_type = "binary"
    constraint_types = {"F": [None, 0],
                        'R': [None, 1000]}
    
    # execute the function
    revbin.build_package()

    # execute the function and assert results of the function
    for reaction in revbin.model.reactions:
        revbin_var = revbin.variables["revbin"][reaction.id]
        assert revbin_var
        assert revbin_var.ub == upper_bound
        assert revbin_var.lb == lower_bound
        assert revbin_var.type == 'binary'
        
        for type in constraint_types:
            revbin_cons = revbin.constraints['revbin{}'.format(type)][reaction.id]
            assert revbin_cons
            assert revbin_cons.lb == constraint_types[type][0]
            assert revbin_cons.ub == constraint_types[type][1]    
            assert revbin_cons.name == '{}_revbin{}'.format(reaction.id, type)