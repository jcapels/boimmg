import math

import cobra
import re
from sympy.solvers import solve
from collections import defaultdict

from boimmgpy import definitions


def map_iJR904():
    model = cobra.io.read_sbml_model(definitions.ROOT_DIR + "/models/iJR904.xml")

    cardiolipin = model.metabolites.get_by_id("clpn_EC_c")

    phosphatidylserine = model.metabolites.get_by_id("ps_EC_c")

    phosphatidylethanolamine = model.metabolites.get_by_id("pe_EC_c")

    phosphatidylglycerol =  model.metabolites.get_by_id("pg_EC_c")

    phosphatidylglycerophosphate =  model.metabolites.get_by_id("pgp_EC_c")

    CDPdiacylglycerol = model.metabolites.get_by_id("cdpdag_EC_c")

    phosphatidate = model.metabolites.get_by_id("pa_EC_c")

    diacylglycerol = model.metabolites.get_by_id("12dgr_EC_c")

    acyl_glycerophosphoglycerol = model.metabolites.get_by_id("agpg_EC_c")

    acyl_phosphatidylglycerol = model.metabolites.get_by_id("apg_EC_c")

    acyl_glycerophosphoethanolamine = model.metabolites.get_by_id("agpe_EC_c")

    acyl_glycerophosphocholine = model.metabolites.get_by_id("agpc_EC_c")






    ################## mapping ###################

    acyl_glycerophosphocholine.annotation = {}
    acyl_glycerophosphocholine.annotation["seed.compound"] = "cpd21879"

    acyl_phosphatidylglycerol.annotation = {}
    acyl_phosphatidylglycerol.annotation["seed.compound"] = "cpd15649"

    acyl_glycerophosphoglycerol.annotation = {}
    acyl_glycerophosphoglycerol.annotation["seed.compound"] = "cpd19396"

    acyl_glycerophosphoethanolamine.annotation = {}
    acyl_glycerophosphoethanolamine.annotation["seed.compound"] = "cpd12547"

    cardiolipin.annotation = {}
    cardiolipin.annotation["seed.compound"] = "cpd22513"

    phosphatidylserine.annotation = {}
    phosphatidylserine.annotation["seed.compound"] = "cpd27342"

    phosphatidylethanolamine.annotation = {}
    phosphatidylethanolamine.annotation["seed.compound"] = "cpd27339"

    phosphatidylglycerol.annotation = {}
    phosphatidylglycerol.annotation["seed.compound"] = "cpd27340"

    phosphatidylglycerophosphate.annotation = {}
    phosphatidylglycerophosphate.annotation["seed.compound"] = "cpd27341"

    CDPdiacylglycerol.annotation = {}
    CDPdiacylglycerol.annotation["seed.compound"] = "cpd22517"

    phosphatidate.annotation = {}
    phosphatidate.annotation["seed.compound"] = "cpd27368"

    diacylglycerol.annotation = {}
    diacylglycerol.annotation["seed.compound"] = "cpd26855"

    # return model
    cobra.io.write_sbml_model(model, "iJR904_mapped.xml")

def balance_reaction(eq):
    Ls=list('abcdefghijklmnopqrstuvwxyz')

    Ss,Os,Es,a,i=defaultdict(list),Ls[:],[],1,1

    for p in eq.split('->'):
        for k in p.split('+'):
            c = [Ls.pop(0), 1]
            for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',k):
                m=1 if m=='' else int(m)
                a*=m
                d=[c[0],c[1]*m*i]
                Ss[e][:0],Es[:0]=[d],[[e,d]]
        i=-1

    Ys=dict((s,eval('Symbol("'+s+'")')) for s in Os if s not in Ls)

    Qs=[eval('+'.join('%d*%s'%(c[1],c[0]) for c in Ss[s]),{},Ys) for s in Ss]+[Ys['a']-a]

    k=solve(Qs,*Ys)

    if k:
        N=[k[Ys[s]] for s in sorted(Ys)]

        g=N[0]

        for a1, a2 in zip(N[0::2],N[1::2]):
            g=math.gcd(g,a2)

        N=[i/g for i in N]

        pM=lambda c: str(c) if c!=1 else ''

        print('->'.join('+'.join(pM(N.pop(0))+str(t) for t in p.split('+')) for p in eq.split('->')))
    else: print('Nope!')


if __name__=="__main__":
    # balance_reaction("H2O+C38H73O13P2+H->C38H75O10P+HO4P")
    map_iJR904()





