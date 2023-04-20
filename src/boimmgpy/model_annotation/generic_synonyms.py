
from boimmgpy.database.accessors.database_access_manager import DatabaseAccessManager


def set_compound_synonyms():
    """Insert synonyms in the genericall compounds to allow the proper work of the annotator class
    """
    driver = DatabaseAccessManager(conf_file_path="my_database.conf").connect()
    synonym_list=["triacylglycerol","1,2-diacyl-sn-glycerol","1-phosphatidyl-1d-myo-inositol3-phosphate","1-phosphatidyl-1d-myo-inositol4-phosphate","1,2-diacyl-3-o-galactosyl-sn-glycerol","digalactosyl-diacylglycerol","monogalactosyldiacylglycerol","digalactosyldiacylglycerol"]
    compound_list=["triacyl-sn-glycerol","1,2-diacyl-sn-glycerol","1,2-diacyl-sn-glycero-3-phospho-1D-myo-inositol-3-phosphate","1,2-diacyl-sn-glycero-3-phospho-1D-myo-inositol-4-phosphate","1,2-diacyl-3-beta-D-galactosyl-sn-glycerol","1,2-diacyl-3-O-[alpha-D-galactosyl-(1->6)-beta-D-galactosyl]-sn-glycerol","1,2-diacyl-3-beta-D-galactosyl-sn-glycerol","1,2-diacyl-3-O-[alpha-D-galactosyl-(1->6)-beta-D-galactosyl]-sn-glycerol"]
    with driver.session() as session1:
        for i in range(len(synonym_list)):
            lipid_compound = compound_list[i]
            synonym = synonym_list[i]
            session1.run('MERGE (s: Synonym {synonym:"%s"})' % str(synonym))
            session1.run("match (l:Compound),(s:Synonym) where l.name=$sl_id and "
                         "s.synonym=$synonym merge (s)-[:is_synonym_of]->(l)", synonym=synonym,
                         sl_id=lipid_compound)
            


if __name__ == '__main__':
    set_compound_synonyms()