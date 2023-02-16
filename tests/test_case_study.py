import sys

import cobra
from rdkit import RDLogger

from boimmgpy import ROOT_DIR
from boimmgpy import RedundantCaseSolver

DATABASE_CONF = ROOT_DIR + "/configs/database_settings.conf"

def test_case_study_rrc(uri,user,password):

    with open(DATABASE_CONF,"w") as file:
        file.write("uri="+str(uri)+"\n")
        file.write("user="+str(user)+"\n")
        file.write("password="+str(password))

    try:
        RDLogger.DisableLog('rdApp.*')

        model = cobra.io.read_sbml_model(ROOT_DIR + "/case_studies/redundant_representation_case/iJR904_mapped.xml")
        components = ["cpd00214", "cpd03847", "cpd05274", "cpd25615", "cpd05237"]

        solver = RedundantCaseSolver(model, "BiGG")

        solver.swap_from_generic(["cpd22513", "cpd15649"], components, True)
        # solver.generateISAreactions()

        os.remove(DATABASE_CONF)

    except:
        os.remove(DATABASE_CONF)
        raise Exception("Not well run")





if __name__=="__main__":
    import os
    os.getcwd()
    uri=sys.argv[1]
    user = sys.argv[2]
    password = sys.argv[3]
    test_case_study_rrc(uri, user, password)
