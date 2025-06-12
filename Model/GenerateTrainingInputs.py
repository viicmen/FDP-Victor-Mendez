#
# FUNCTIONS FOR GENERATING THE RANDOM INPUTS FOR TRAINIG THE MODEL
#

import json
import os
import math
import random
import glob
import pandas as pd

def PosOrNeg(output_tomtom, json_file):
'''
Parses output of "tomtom.sh" and classifies each query-target match as positive (p-value < 0.01) or negative (p-value > 0.5). 
It stores the classification in a JSON dictionary.
'''
    dic = {}
    with open(output_tomtom, 'r') as f:
        f.readline()
        line = f.readline()
        while len(line)>3:

            line = line.strip().split()
            query_id = line[0]
            target_id = line[1]
            pvalue = float(line[3])
            if pvalue > 0.5: pvalue = 0
            elif pvalue < 0.01: pvalue = 1
            else:
                line = f.readline()
                continue
            if query_id in dic:
                dic[query_id][pvalue].append(target_id)
            else:
                dic[query_id] = {1: [], 0: []}
                dic[query_id][pvalue].append(target_id)

            line = f.readline()

    with open(json_file, 'w') as f:
        json.dump(dic,f)


def GetIDP(output_search, json_file):
'''
Reads pairwise sequence identity percentages from "mmseqs.sh" output and stores them in a nested dictionary (query → target → identity), saved as a JSON file.
'''
    dic={}

    with open(output_search, 'r') as f:

        f.readline()
        for line in f:
            line = line.strip().split()
            query_id = line[0]
            target_id = line[1]
            idp = float(line[2])

            if query_id in dic:
                dic[query_id][target_id] = idp
            else:
                dic[query_id] = {target_id: idp}

    with open(json_file, 'w') as f:
        json.dump(dic,f)


def ItersPerFam(domanins_file, json_file):
'''
Counts how many sequences belong to each DBD family and assigns a number of iterations for sampling based on family size. 
Outputs the result to a JSON file.
'''
	families = ['AFT', 'AP2', 'ARID', 'AT_hook', 'B3', 'BEN', 'BES1_N',
            'BrkDBD', 'bZIP', 'CBFB_NFYA', 'CBFD_NFYB_HMF', 'CENPB',
            'CG-1', 'COE1_DBD', 'Copper-fist', 'CP2', 'CSD', 'CUT',
            'CXC', 'DM', 'DUF260', 'DUF573', 'E2F_TDP', 'EIN3', 'Ets',
            'FAR1', 'FLYWCH', 'Forkhead', 'GAGA_bind', 'GATA', 'GCM',
            'GCR1_C', 'GRAS', 'GTF2I', 'HLH', 'HMG', 'Homeodomain',
            'HPD', 'HSF_DNA-bind', 'HTH_psq', 'IBD', 'IRF', 'KilA-N',
            'LAG1-DNAbind', 'LOB', 'MADF_DNA_bdg', 'MBD', 'MH1', 'Myb',
            'NAM', 'NDT80_PhoG', 'Nrf1_DNA-bind', 'P53', 'PAX', 'PCC',
            'Pou', 'Prox1', 'PurA', 'RFX_DNA_binding', 'RHD_DNA_bind',
            'RRM_1', 'Runt', 'SAM_LFY', 'SAND', 'SBP', 'SPT2', 'SRF-TF',
            'STAT_bind', 'Stb3', 'STE', 'T-box', 'TBP', 'TCP', 'TCR',
            'TEA', 'THAP', 'TIG', 'UNKNOWN', 'Vhr1', 'WRC', 'WRKY',
            'zf-BED', 'zf-C2H2', 'zf-C2HC', 'zf-C4', 'zf-CCCH',
            'zf-CXXC', 'ZF-HD_dimer', 'Zn_clus', 'Zn_ribbon_Dof']
            
	dic={}
	for fam in families:
	    dic[fam] = 0

	with open(domains_file) as f:

	    f.readline()
	    for line in f.readlines():
	        dic[line.split()[1]] += 1

	for fam in families:

	    if dic[fam] >= 300: dic[fam] = 10
	    elif dic[fam] >= 100: dic[fam] = 20
	    elif dic[fam] >= 50: dic[fam] = 30
	    elif dic[fam] >= 25: dic[fam] = 40
	    elif dic[fam] > 0: dic[fam] = math.ceil(1000/dic[fam])

	with open(json_file, 'w') as f:
	    json.dump(dic, f)


def GetInputsNN(dbd, id_per, classifier, db):
'''
Selects five positive and five negative PWM examples for a given DBD from NN predictions, retrieves identity percentages, and returns paths and IDP values.
'''
    positive = []
    negative = []
    pos_idp = []
    neg_idp = []

    for _ in range(5):

        pos = random.choice(classifier["1"])
        neg = random.choice(classifier["0"])

        pos_idp.append(float(id_per[dbd][pos]) if pos in id_per[dbd] else 0.1)
        neg_idp.append(float(id_per[dbd][neg]) if neg in id_per[dbd] else 0.1)
        positive.append(db+"NN_r/pwms/"+pos)
        negative.append(db+"NN_r/pwms/"+neg)

    return positive, pos_idp, negative, neg_idp

def GetInputsModcre(dbd, idp_modcre, classifier, modcre_pwms, db):
'''
Similar to GetInputsNN, but uses PWMs from a Modcre dataset. 
It selects five positives and five negatives, retrieves their IDPs, and returns the paths and values.
'''
    positive = []
    negative = []
    pos_idp = []
    neg_idp = []

    motifs = [os.path.basename(f) for f in glob.glob(modcre_pwms+dbd+"*meme")]

    for _ in range(5):

        pos = random.choice(motifs)
        pos_seq = "_".join(pos.split(":")[-1].split("_")[1:3])
        pos = pos[:-5]
          
        if dbd not in idp_modcre:
            pos_idp.append(0)

        elif pos_seq in idp_modcre[dbd]:
            pos_idp.append(idp_modcre[dbd][pos_seq])

        else:
            pos_idp.append(0.1)

        positive.append(db+"Modcre/pwms/"+pos)

    for _ in range(5):

        neg = []

        neg = random.choice(classifier["0"])
        neg = [os.path.basename(f) for f in glob.glob(modcre_pwms+neg+"*meme")]
        if neg:
            neg = random.choice(neg)
            neg_seq = "_".join(neg.split(":")[-1].split("_")[1:3])
            neg = neg[:-5]

        if dbd not in idp_modcre:
            neg_idp.append(0)

        elif neg_seq in idp_modcre[dbd]:
            neg_idp.append(idp_modcre[dbd][neg_seq])

        else:
            neg_idp.append(0.1)

        negative.append(db+"Modcre/pwms/"+neg)

    return positive, pos_idp, negative, neg_idp


def CreateRandomInputs(classifier, idp_nn, idp_modcre, output_dir, modcre_pwms, domains_file, iters, db, section):
'''
Creates training data for the CNN by sampling positive/negative examples from both NN and Modcre predictions for each DBD, using the number of iterations per family. 
The data is saved in JSONL format for each DBD.
'''
    f = open(classifier, "r")
    classifier = json.load(f)
    f.close()

    f = open(idp_nn, "r")
    idp_nn = json.load(f)
    f.close()

    f = open(idp_modcre, "r")
    idp_modcre = json.load(f)
    f.close()

    domains = pd.read_csv(domains_file, sep="\t")

    f = open(iters, "r")
    iters = json.load(f)
    f.close()

    families = list(iters.keys())

    for dbd in list(classifier.keys()):
        dom = domains[domains["seq_id"] == dbd]["domains"].values[0]

        f = open(output_dir+dbd+".jsonl", "w")

        dic = {"1": {"nn": {"ids":[], "ipd":[]}, "modcre": {"ids":[], "ipd":[]}}, "0": {"nn": {"ids":[], "ipd":[]}, "modcre": {"ids":[], "ipd":[]}}, "family": families.index(dom)}

        for i in range(iters[dom]):

            dic["1"]["nn"]["ids"], dic["1"]["nn"]["ipd"], dic["0"]["nn"]["ids"], dic["0"]["nn"]["ipd"] = GetInputsNN(dbd, idp_nn, classifier[dbd], db)
            dic["1"]["modcre"]["ids"], dic["1"]["modcre"]["ipd"], dic["0"]["modcre"]["ids"], dic["0"]["modcre"]["ipd"] = GetInputsModcre(dbd, idp_modcre, classifier[dbd], modcre_pwms, db)

            f.write(json.dumps(dic))
            f.write("\n")

        f.close()
