# Biochemical cOmplex data Integration in Metabolic Models at Genome-scale (BOIMMG)

Several ontologies were designed to represent the structural relationships between chemical species. However, they fail in capturing the functional relationships between compounds and their biosynthetic precursors. This is highly relevant in the context of GSM modelling, as it eases the process of converting generic biosynthetic pathways into specific ones and vice-versa. Moreover, to the best of our knowledge, there is not any computational approach to integrate such information.

BOIMMG is a novel and modular approach aiming at tackling several issues in the representation of lipids in GSM models.  

![alt text](boimmg_pipeline.PNG)
* Figure 1 - BOIMMG's pipeline - 1 - Integration of several databases; 2 -Semi-automated knowledge expansion 3 - integration of those ontologies into GSM models

In this frame, you are going to find the 2nd and 3rd modules only. The 2nd module is present in the **ontologies_generators** folder inside **boimmgpy** folder, whereas the latter is in a package format in **boimmgpy/service/**.

If one wants to test for the case studies, the **case_studies** folder is provided. However, one would need either the database with relevant information or the access to it, which is unfortunately not provided yet.

However, a web-service is available [here](https://boimmg.bio.di.uminho.pt/).

## Redundant representation case

###*Escherichia coli* electron-transfer quinones' generalization:
1. Quinones in model: cpd15560
1. Quinones to introduce: cpd11669
1. Click in the button **Add more quinones to swap**
1. Quinones in model: cpd15500
1. Quinones to introduce: cpd11451
1. Metabolite format: BiGG
1. Upload the file **boimmgpy/models/iML1515.xml**
    
### *Escherichia coli*'s previously generalized model granulation to ubiquinone-5
1. Quinones in model: cpd11669
1. Quinones to introduce: C_BOIMMG_749196
1. Metabolite format: BiGG
1. Upload the file **case_studies/simple_representation_case/iML1515_generalized.xml**

### *Saccharomyces cerevisiae*'s generalization
1. Quinones in model: cpd15290
1. Quinones to introduce: cpd11669
1. Metabolite format: BiGG
1. Upload the file **boimmgpy/models/iMM904.xml**

### *Saccharomyces cerevisiae*'s previously generalized model granulation to ubiquinone-6
1. Quinones in model: cpd11669
1. Quinones to introduce: cpd15290
1. Metabolite format: BiGG
1. Upload the file **case_studies/simple_representation_case/iMM904_generalized.xml**

## Redundant representation case

### *Escherichia coli*'s lipid granulation

Firstly, the iJR904 *Escherichia coli* model was mapped to contain database cross-references complaint with BOIMMG's.
Such operation can be found in the following python script: **case_studies/redundant_representation_case/map_ecoli_model**

1. Generic lipid in model: cpd22513 
1. Fatty acids: cpd00214,cpd03847,cpd05274,cpd25615,cpd05237
1. Click in the button **Add more lipids to granulate**
1. Generic lipid in model: cpd15649 
1. Fatty acids: cpd00214,cpd03847,cpd05274,cpd25615,cpd05237
1. Mixed components: same components
1. Metabolite format: BiGG
1. Upload the file **case_studies/redundant_representation_case/iJR904_mapped.xml**

