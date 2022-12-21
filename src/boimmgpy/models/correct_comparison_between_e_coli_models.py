import cobra

from src.boimmgpy import definitions


def map_e_coli_granulated_model():
    #################### pe ##################33
    # model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iAF1260b.xml")
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iML1515.xml")
    # model = cobra.io.read_sbml_model("C:/Users/Joao/Desktop/thesis_jcapela/BOIMMGpy/models/iAF1260b.xml")

    for c in ["c", "p"]:
        pe181 = model.metabolites.get_by_id("pe181_" + c)
        pe181.annotation["seed.compound"] = "cpd25142"
        pe181.annotation["boimmg.compound"] = "C_BOIMMG_8715"

        pe161 = model.metabolites.get_by_id("pe161_" + c)
        pe161.annotation["seed.compound"] = "cpd31682"
        pe161.annotation["boimmg.compound"] = "C_BOIMMG_12583"

        pe16 = model.metabolites.get_by_id("pe160_" + c)
        pe16.annotation["boimmg.compound"] = "C_BOIMMG_12595"

        pe14 = model.metabolites.get_by_id("pe140_" + c)
        pe14.annotation["seed.compound"] = "cpd33856"
        pe14.annotation["boimmg.compound"] = "C_BOIMMG_12604"

        pe141 = model.metabolites.get_by_id("pe141_" + c)
        pe141.annotation["boimmg.compound"] = "C_BOIMMG_12585"

        ############### pg ################

        pe141 = model.metabolites.get_by_id("pg141_" + c)
        pe141.annotation["boimmg.compound"] = "C_BOIMMG_7456"

        pe181 = model.metabolites.get_by_id("pg181_" + c)
        pe181.annotation["boimmg.compound"] = "C_BOIMMG_7408"

        pe141 = model.metabolites.get_by_id("pg140_" + c)
        pe141.annotation["seed.compound"] = "cpd35089"
        pe141.annotation["boimmg.compound"] = "C_BOIMMG_7474"

        pe160 = model.metabolites.get_by_id("pg160_" + c)
        pe160.annotation["boimmg.compound"] = "C_BOIMMG_427"

        pe161 = model.metabolites.get_by_id("pg161_" + c)
        pe161.annotation["boimmg.compound"] = "C_BOIMMG_7454"

        ############### ps ################

        if c == "c":


            ps160 = model.metabolites.get_by_id("ps160_" + c)
            ps160.annotation["boimmg.compound"] = "C_BOIMMG_729826"

            ps161 = model.metabolites.get_by_id("ps161_" + c)
            ps161.annotation["boimmg.compound"] = "C_BOIMMG_729871"

            ps140 = model.metabolites.get_by_id("ps140_" + c)
            ps140.annotation["boimmg.compound"] = "C_BOIMMG_729819"

            ps141 = model.metabolites.get_by_id("ps141_" + c)
            ps141.annotation["boimmg.compound"] = "C_BOIMMG_729870"

            ps180 = model.metabolites.get_by_id("ps180_" + c)
            ps180.annotation["boimmg.compound"] = "C_BOIMMG_729833"

            ps181 = model.metabolites.get_by_id("ps181_" + c)
            ps181.annotation["boimmg.compound"] = "C_BOIMMG_729912"

        ############## pgp ##############

        pgp161 = model.metabolites.get_by_id("pgp161_" + c)
        pgp161.annotation["boimmg.compound"] = "C_BOIMMG_296512"

        pgp140 = model.metabolites.get_by_id("pgp140_" + c)
        pgp140.annotation["boimmg.compound"] = "C_BOIMMG_296532"

        pgp141 = model.metabolites.get_by_id("pgp141_" + c)
        pgp141.annotation["boimmg.compound"] = "C_BOIMMG_296513"

        pgp160 = model.metabolites.get_by_id("pgp160_" + c)
        pgp160.annotation["boimmg.compound"] = "C_BOIMMG_296523"

        pgp181 = model.metabolites.get_by_id("pgp181_" + c)
        pgp181.annotation["boimmg.compound"] = "C_BOIMMG_293736"

        ################### cdpdag ######################3

        if c == "c":
            cdp141 = model.metabolites.get_by_id("cdpdtdec7eg_" + c)
            cdp141.annotation["boimmg.compound"] = "C_BOIMMG_574506"

            cdp181 = model.metabolites.get_by_id("cdpdodec11eg_" + c)
            cdp181.annotation["boimmg.compound"] = "C_BOIMMG_540741"

            # cdp180 = model.metabolites.get_by_id("cdpdodecg_" + c)
            # cdp180.annotation["boimmg.compound"] = "C_BOIMMG_540741"

            cdp160 = model.metabolites.get_by_id("cdpdhdecg_" + c)
            cdp160.annotation["seed.compound"] = "cpd23593"
            cdp160.annotation["boimmg.compound"] = "C_BOIMMG_239898"

            cdp140 = model.metabolites.get_by_id("cdpdtdecg_" + c)
            cdp140.annotation["seed.compound"] = "cpd32277"
            cdp140.annotation["boimmg.compound"] = "C_BOIMMG_574421"

            cdp161 = model.metabolites.get_by_id("cdpdhdec9eg_" + c)
            cdp161.annotation["boimmg.compound"] = "C_BOIMMG_524486"

        ############# PA ###############

        pa141 = model.metabolites.get_by_id("pa140_" + c)
        pa141.annotation["boimmg.compound"] = "C_BOIMMG_38435"

        pa141 = model.metabolites.get_by_id("pa141_" + c)
        pa141.annotation["boimmg.compound"] = "C_BOIMMG_38417"

        pa160 = model.metabolites.get_by_id("pa160_" + c)
        pa160.annotation["boimmg.compound"] = "C_BOIMMG_38416"

        pa161 = model.metabolites.get_by_id("pa161_" + c)
        pa161.annotation["boimmg.compound"] = "C_BOIMMG_423"

        pa180 = model.metabolites.get_by_id("pa181_" + c)
        pa180.annotation["seed.compound"] = "cpd33760"
        pa180.annotation["boimmg.compound"] = "C_BOIMMG_38371"

        ##############3 12dgr #################3

        dgr140 = model.metabolites.get_by_id("12dgr140_" + c)
        dgr140.annotation["seed.compound"] = "cpd16468"
        dgr140.annotation["boimmg.compound"] = "C_BOIMMG_53297"

        dgr141 = model.metabolites.get_by_id("12dgr141_" + c)
        dgr141.annotation["boimmg.compound"] = "C_BOIMMG_53244"

        dgr160 = model.metabolites.get_by_id("12dgr160_" + c)
        dgr160.annotation["seed.compound"] = "cpd26605"
        dgr160.annotation["boimmg.compound"] = "C_BOIMMG_61409"

        dgr161 = model.metabolites.get_by_id("12dgr161_" + c)
        dgr161.annotation["boimmg.compound"] = "C_BOIMMG_61238"

        dgr180 = model.metabolites.get_by_id("12dgr180_" + c)
        dgr180.annotation["boimmg.compound"] = "C_BOIMMG_70311"

        dgr181 = model.metabolites.get_by_id("12dgr181_" + c)
        dgr181.annotation["boimmg.compound"] = "C_BOIMMG_70100"

        ############## agpg ##############

        if c == "p":
            clpn161 = model.metabolites.get_by_id("clpn161_" + c)
            clpn161.annotation["boimmg.compound"] = "C_BOIMMG_341201"

            clpn160 = model.metabolites.get_by_id("clpn160_" + c)
            clpn160.annotation["boimmg.compound"] = "C_BOIMMG_322789"

            clpn181 = model.metabolites.get_by_id("clpn181_" + c)
            clpn181.annotation["boimmg.compound"] = "C_BOIMMG_347940"

            agpg140 = model.metabolites.get_by_id("1agpg140_" + c)
            agpg140.annotation["boimmg.compound"] = "C_BOIMMG_452"

            agpg141 = model.metabolites.get_by_id("1agpg141_" + c)
            agpg141.annotation["boimmg.compound"] = "C_BOIMMG_16364"

            agpg160 = model.metabolites.get_by_id("1agpg160_" + c)
            agpg160.annotation["boimmg.compound"] = "C_BOIMMG_314"

            agpg161 = model.metabolites.get_by_id("1agpg161_" + c)
            agpg161.annotation["boimmg.compound"] = "C_BOIMMG_16363"

            agpg181 = model.metabolites.get_by_id("1agpg181_" + c)
            agpg181.annotation["boimmg.compound"] = "C_BOIMMG_16317"



        ############## apg ##############

        if c == "c":

            apg140 = model.metabolites.get_by_id("apg140_" + c)
            apg140.annotation["boimmg.compound"] = "C_BOIMMG_765781"

            apg141 = model.metabolites.get_by_id("apg141_" + c)
            apg141.annotation["boimmg.compound"] = "C_BOIMMG_765700"

            apg161 = model.metabolites.get_by_id("apg161_" + c)
            apg161.annotation["boimmg.compound"] = "C_BOIMMG_765716"

            apg160 = model.metabolites.get_by_id("apg160_" + c)
            apg160.annotation["boimmg.compound"] = "C_BOIMMG_765924"

            apg181 = model.metabolites.get_by_id("apg181_" + c)
            apg181.annotation["boimmg.compound"] = "C_BOIMMG_765860"

        ############# agpe #######################

        if c == "p":
            agpe140 = model.metabolites.get_by_id("1agpe140_" + c)
            agpe140.annotation["boimmg.compound"] = "C_BOIMMG_27943"

            agpe141 = model.metabolites.get_by_id("1agpe141_" + c)
            agpe141.annotation["boimmg.compound"] = "C_BOIMMG_27925"

            agpe160 = model.metabolites.get_by_id("1agpe160_" + c)
            agpe160.annotation["boimmg.compound"] = "C_BOIMMG_291"

            agpe161 = model.metabolites.get_by_id("1agpe161_" + c)
            agpe161.annotation["seed.compound"] = "cpd25196"
            agpe161.annotation["boimmg.compound"] = "C_BOIMMG_27923"

            agpe181 = model.metabolites.get_by_id("1agpe181_" + c)
            agpe181.annotation["seed.compound"] = "cpd25197"
            agpe181.annotation["boimmg.compound"] = "C_BOIMMG_27877"

        ###################3 agpc ###############3


        # cobra.io.write_sbml_model(model,"iAF1260b_mapped.xml")
        cobra.io.write_sbml_model(model, "../case_studies/redundant_representation_case/iML1515_mapped.xml")

if __name__ == "__main__":
    map_e_coli_granulated_model()






