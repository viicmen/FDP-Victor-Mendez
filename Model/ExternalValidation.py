#
# FUNCTIONS FOR GENERATING THE INPUTS FOR TESTING
#

import torch
import random
import json
from GenerateTraininginputs import CreateLogo, TensorToMeme
def CreateIDPRanges(file_path, out_json):
    id_ranges = [(100, 80), (80, 60), (60, 40), (40, 20), (20, 0)]
    result = {}
    def get_range(pident):
        for high, low in id_ranges:
            if high >= pident > low:
                return f"{high}-{low}"
    with open(file_path, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            query_id = parts[0]
            target_id = parts[1]
            print(query_id, target_id)
            if query_id == target_id: continue
            pident = float(parts[2])*100
            print(pident)
            range_key = get_range(pident)
            if query_id not in result:
                result[query_id] = {f"{h}-{l}": [] for h, l in id_ranges}
            result[query_id][range_key].append(target_id)
    
    with open(out_json, 'w') as f:
        json.dump(result)
        
    return result

def FindModcrePreds(pwm_dir, pident_dic, out_json):
    
    pwms = os.listdir(pwm_dir)
    new = {}
    for m, ranges in pident_dic.items():
        new[m] = {}
        motifs = [x for x in pwms if x.startswith(m)]
        for r, ids in ranges.items():
            new[m][r]  = []
            for i in ids:
                for j in motifs:
                    if i in j:
                        new[m][r].append(j)
                        
    with open(out_json, 'w') as f:
        json.dump(new)
        
    return new

def GetPMW_Test(pwm_id):

    pwm_id:

    with open(pwm_id + '.meme', 'r') as pwm:

        lines = [l for l in pwm.readlines() if l != '\n']

    i = 0
    while lines[i].split()[0][0] != 'l':
        i+=1

    if lines[-1][0] == 'U': lines = lines[:-1]
    probs = [[float(p) for p in l.strip().split()] for l in lines[i+1:]]

    right = (30-len(probs))//2
    left = 30-len(probs)-right

    right = [[0.25 for _ in range(4)] for _ in range(right)]
    left = [[0.25 for _ in range(4)] for _ in range(left)]

    probs = left + probs + right

def CreateInputs_Test(pwms, idp):

    pwms = [GetPMW(random.choice(pwms)) for _ in range(5)]
    input_tensor = torch.tensor([[pwms[i][j]+[idp[i]] for j in range(40)] for i in range(50)])
    input_tensor = input_tensor.permute(1, 0, 2)
    return input_tensor

def CreateOneHot(family):

    one_hot = torch.zeros(1, 7)
    if family == 'AP2': one_hot[0][0] = 1
    elif family == 'bZIP': one_hot[0][1] = 1
    elif family == 'HLH': one_hot[0][2] = 1
    elif family == 'Homeodomain': one_hot[0][3] = 1
    elif family == 'Myb': one_hot[0][4] = 1
    elif family == 'zf-C2H2': one_hot[0][5] = 1
    elif family == 'Zn_clus': one_hot[0][6] = 1
    return one_hot

def PredictPWM(model_path, motif_id, nn_idp, nn_json, modcre_idp, modcre_json, family):
    
    with open('nn_json') as f, open('modcre_json') as g:
        nn_dic = json.load(f)
        modcre_dic = json.load(g)
        
    model = torch.load(model_path)

    ipds = {"100_80":90 , "80_60":70, "60-40":50, "40_20":30, "20_0":10}
    input_tensor_nn = CreateInput(nn_dic[motif_id][nn_idp], ipds[nn_idp])
    input_tensor_modcre = CreateInput(modcre_dic[motif_id][modcre_idp], ipds[modcre_idp])
    onehot = CreateOneHot(family)

    pred = model(input_tensor_nn, input_tensor_modcre, onehot)
    CreateLogo(pred, f'test_logos/{motif_id}_pred.png')
    TensorToMeme(pred, motif_id, f'test_memes/{motif_id}_pred.png')
    
    
    
    
