import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import os
import json
import optuna
from optuna.visualization.matplotlib import plot_optimization_history
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from torch.utils.data import Dataset
from sklearn.model_selection import StratifiedShuffleSplit


def GetPvalue(motif_id, real_path, pred_path):
'''
Executes TomTom for comparing the Real motif against the Prediction and retrieves the p-value from the output.
'''
    os.system(f'tomtom -thresh 1 {real_path}/{motif_id}.meme {pred_path}/{motif_id}.meme')

    f = open('tomtom_out/tomtom.tsv', 'r')
    f.readline()
    line = f.readline().split()
    pvalue = line[5]

    return pvalue


def TensorToMeme(pwm, pwm_id, out_file):
'''
converts a predicted Position Weight Matrix (PWM), given as a PyTorch tensor, into the MEME motif format. 
'''
    matrix = pwm.detach().numpy()
    f = open(out_file, 'w')
    f.write('MEME version 4\n\n')
    f.write('ALPHABET= ACGT\n\n')
    f.write('Background letter frequencies\n')
    f.write('A 0.25 C 0.25 G 0.25 T 0.25\n\n')
    f.write(f'MOTIF {pwm_id}\n\n')
    f.write(f'letter-probability matrix: w= {str(len(matrix))}\n')

    for pos in matrix:
        for n in pos:
            f.write(str(pos)+' ')
        f.write('\n')
    f.write('URL')
    f.close()

def CreateLogo(pwm_pred, name):
'''
Generates a sequence logo image from a predicted PWM tensor using the logomaker library. 
It first converts the tensor to a NumPy array and into a DataFrame with columns for A, C, G, and T. 
It then transforms the PWM into an information content matrix and uses this to draw the logo. 
The resulting plot is saved as a PNG image.
'''
    pwm_np = pwm_pred.detach().numpy()
    df = pd.DataFrame(pwm_np, columns=['A', 'C', 'G', 'T'])
    info_df = logomaker.transform_matrix(df, from_type='probability', to_type='information')
    plt.figure(figsize=(10,4))
    logo = logomaker.Logo(info_df)
    plt.title("Sequence Logo de PWM Predicha")
    plt.xlabel("Posición")
    plt.ylabel("Información")
    plt.tight_layout()
    plt.savefig(name+'.png')
    plt.close()


class CustomDataset(Dataset):
'''
class is a PyTorch dataset implementation designed to store and manage TF input data for model training and evaluation. 
The class stores the input tensors, the corresponding labels and subgroup identifiers, which represent different TF families or source databases. 
The AddInput method allows new data instances to be added, including their real PWM, model inputs, and metadata. 
The dataset supports indexing and length querying, enabling it to be used with PyTorch’s DataLoader.
''' 
    def __init__(self):

        self.data = []
        self.labels = []
        self.subgroups = []

    def AddInput(self, input_id, data, pwm_real, label, db):

        input_tensor_nn, input_tensor_modcre, one_hot = DataLoader(data, label)
        fam = f"{data['family']}_{db}_{label}"
        self.data.append([input_id, input_tensor_nn, input_tensor_modcre, one_hot, pwm_real})
        self.labels.append(label)
        self.subgroups.append(fam)

    def __getitem__(self, idx):

        x = self.data[idx]
        y = self.labels[idx]

        return x, y

    def __len__(self):

        return len(self.data)

def LoadCustomDataset():
'''
Loads input data and real PWMs. It reads motif identifiers from directory listings, retrieves their real PWMs from .meme files, and loads associated input features from the generated JSONL files. 
Each example is added to an instance of CustomDataset using the AddInput method. 
'''
    dataset = CustomDataset()

    for db in ['JASPAR', 'cis-bp', 'hocomoco']:
        for label in ['1','0']:
            for db in dbs:
                inputs = [f[:-5] for f in os.listdir(db+'/NN_nr/pwms/')]

                for input in inputs:
                    pwm_real = torch.tensor(GetPMW(f'{db}/NN_nr/pwms/{input}', False))

                    with open(f'{db}/inputs/{input}.jsonl', encoding='utf-8') as f:
                        for line in f:
                            data = json.loads(line)

                            dataset.AddInput(input, data, pwm_real, label, db)

    return dataset


def GetPMW(pwm_id, pad=True):
'''
Loads a PWM from a MEME-format file, returning the PWM as a nested list of nucleotide probabilities. 
If the PWM is shorter than 30 positions, it can be padded symmetrically on both sides with uniform probabilities to create a standardized input length. 
This padding allows all PWMs to be input into a neural network with a fixed input shape. If the PWM is missing, a default 30x4 matrix of zeros is returned.
'''
    if pwm_id:

        with open(pwm_id + '.meme', 'r') as pwm:

            lines = [l for l in pwm.readlines() if l != '\n']

        i = 0
        while lines[i].split()[0][0] != 'l':
            i+=1

        if lines[-1][0] == 'U': lines = lines[:-1]
        probs = [[float(p) for p in l.strip().split()] for l in lines[i+1:]]

        if pad:
            right = (30-len(probs))//2
            left = 30-len(probs)-right

            right = [[0.25 for _ in range(4)] for _ in range(right)]
            left = [[0.25 for _ in range(4)] for _ in range(left)]

            probs = left + probs + right

        return probs

    return [[0.0 for _ in range(4)] for _ in range(40)]

def DataLoader(data, label):
'''
Constructs model inputs from a data dictionary and a label key ('1' or '0'). 
It generates two input tensors, one from NN another from ModCRE motifs, and a one-hot encoded vector representing the TF family.
These are the required input componenets for the neural network model to make predictions.
'''
    input_tensor_nn = CreateInput(data[label]['nn']['ids'], data[label]['nn']['ipd'])
    input_tensor_modcre = CreateInput(data[label]['modcre']['ids'], data[label]['modcre']['ipd'])
    one_hot = CreateOneHot(data['family'])

    return input_tensor_nn, input_tensor_modcre, one_hot

def CreateInput(pwms, idp):
'''
Takes a list of PWM file paths and a list of identity values and constructs a tensor suitable for input into a neural network. 
Each PWM is loaded, and the identity value is appended to each position. 
The result is a 30x5x5 tensor representing the 5 pwms at every positions plus their identity.
'''
    pwms = [GetPMW(pwm_id) for pwm_id in pwms]
    input_tensor = torch.tensor([[pwms[i][j]+[idp[i]] for j in range(30)] for i in range(5)])
    input_tensor = input_tensor.permute(1, 0, 2)
    return input_tensor

def CreateOneHot(family):
'''
Converts a TF family identifier into a one-hot encoded tensor. 
This is a fixed-size vector in which only the position corresponding to the given family is set to 1
'''
    one_hot = torch.zeros(1, 7)
    if family == 1: one_hot[0][0] = 1       # AP2
    elif family == 8: one_hot[0][1] = 1     # bZIP
    elif family == 34: one_hot[0][2] = 1    # HLH
    elif family == 36: one_hot[0][3] = 1    # Homeodomain
    elif family == 48: one_hot[0][4] = 1    # Myb
    elif family == 82: one_hot[0][5] = 1    # zf-C2H2
    elif family == 88: one_hot[0][6] = 1    # Zn_clus
    return one_hot



class FC_PWM(nn.Module):
'''
The Fully Connected Neural Network arquitecture as a Pytorch Module Class
'''
    def __init__(self, fc_size, hl_size, drop_out):
        super(FC_PWM, self).__init__()

        flatten_dim = 30 * 5 * 5

        self.fc_nn_1 = nn.Linear(flatten_dim, fc_size)
        self.fc_nn_2 = nn.Linear(fc_size, fc_size//2)

        self.fc_mod_1 = nn.Linear(flatten_dim, fc_size)
        self.fc_mod_2 = nn.Linear(fc_size, fc_size//2)

        self.fc_combined_1 = nn.Linear(fc_size//2 * 5 + 7, hl_size)
        self.fc_combined_2 = nn.Linear(hl_size, 30 * 4)

        self.drop = nn.Dropout(drop_out)

    def forward(self, x_nn, x_mod, one_hot):
        x_nn = x_nn.reshape(1, -1)
        x_mod = x_mod.reshape(1, -1)

        x_nn = F.relu(self.fc_nn_1(x_nn))
        x_nn = self.drop(x_nn)
        x_nn = F.relu(self.fc_nn_2(x_nn))
        x_nn = self.drop(x_nn)

        x_mod = F.relu(self.fc_mod_1(x_mod))
        x_mod = self.drop(x_mod)
        x_mod = F.relu(self.fc_mod_2(x_mod))
        x_mod = self.drop(x_mod)

        x_comb = torch.cat((x_nn, x_mod, one_hot), dim=1)
        x_comb = F.relu(self.fc_combined_1(x_comb))
        x_out = F.relu(self.fc_combined_2(x_comb))

        out = x_out.view(-1, 30, 4)
        out = F.softmax(out, dim=2)
        return out


class CNN_PWM(nn.Module):
'''
The Convolutional Neural Network arquitecture as a Pytorch Module Class
'''
    def __init__(self, hl_size, drop_out):
        super(CNN_PWM, self).__init__()

        self.conv1_nn = nn.Conv2d(in_channels=30 out_channels=32, kernel_size=(3,3), padding=1)
        self.pool1_nn = nn.MaxPool2d(kernel_size=(2,1))
        self.conv2_nn = nn.Conv2d(in_channels=32, out_channels=16, kernel_size=(5,3), padding=(2,1))
        self.pool2_nn = nn.MaxPool2d(kernel_size=(2,1))

        self.conv1_mod = nn.Conv2d(in_channels=30, out_channels=32, kernel_size=(3,3), padding=1)
        self.pool1_mod = nn.MaxPool2d(kernel_size=(2,1))
        self.conv2_mod = nn.Conv2d(in_channels=32, out_channels=16, kernel_size=(5,3), padding=(2,1))
        self.pool2_mod = nn.MaxPool2d(kernel_size=(2,1))

        self.fc1 = nn.Linear(16 * 2 * 5 + 7, hl_size)
        self.fc2 = nn.Linear(hl_size, 30 * 4)

        self.drop = nn.Dropout(drop_out)

    def forward(self, x_nn, x_mod, one_hot):

        x_nn = F.relu(self.conv1_nn(x_nn))
        x_nn = self.pool1_nn(x_nn)
        x_nn = self.drop(x_nn)
        x_nn = F.relu(self.conv2_nn(x_nn))
        x_nn = self.pool2_nn(x_nn)
        x_nn = self.drop(x_nn)


        x_mod = F.relu(self.conv1_mod(x_mod))
        x_mod = self.pool1_mod(x_mod)
        x_mod = self.drop(x_mod)
        x_mod = F.relu(self.conv2_mod(x_mod))
        x_mod = self.pool2_mod(x_mod)
        x_mod = self.drop(x_mod)

        x_cat = torch.cat((x_nn, x_mod), dim=1)
        x_flat = x_cat.view(1,-1)
        x_comb = torch.cat((x_flat, one_hot), dim=1)

        x_fc1 = F.relu(self.fc1(x_comb))
        x_fc2 = F.relu(self.fc2(x_fc1))

        out = x_fc2.view(-1, 30, 4)

        out = F.softmax(out, dim=2)
        return out

def JSDScore(pwm1, pwm2):
'''
Calculates the maximum similarity between the real and the predicted PWMs using Jensen-Shannon Divergence. 
It aligns the real across the predicted one, computes JSD at each position, and returns the highest similarity score as 1 - JSD / log(2), ensuring values range from 0 (different) to 1 (identical).
'''
    pwm1 = pwm1 + 1e-10
    pwm2 = pwm2 + 1e-10
    L1, L2 = pwm1.size(0), pwm2.size(0)
    max_score = torch.tensor(-float('inf'))

    for offset in range(0, L1 - L2 + 1):
        seg1 = pwm1[offset : offset + L2]

        m = 0.5 * (seg1 + pwm2)
        kl_p = torch.sum(seg1 * (torch.log(seg1 / m)), dim=1)
        kl_q = torch.sum(pwm2 * (torch.log(pwm2 / m)), dim=1)
        jsd = 0.5 * (kl_p + kl_q).mean()

        score = 1 - jsd/torch.log(torch.tensor(2))
        max_score = torch.max(max_score, score)

    return max_score

def EntropyPenalty(pwm, penalty_strength):
'''
Penalises PWMs with high entropy.
'''
    entropy = -torch.sum(pwm * torch.log(pwm + 1e-10), dim=-1)
    max_entropy = torch.log(torch.tensor(4.0))
    penalty = torch.mean(entropy / max_entropy)
    return penalty_strength * penalty

def KLPenalty(pwm, penalty_strength):
'''
Penalises PWMs with peaks or close to a uniform distribution.
'''
    uniform = torch.ones_like(pwm) / 4.0
    kl_div = torch.sum(pwm * torch.log(pwm / uniform), dim=-1)
    return penalty_strength * torch.mean(kl_div)

def NoPenalty(pwm, penalty_strength):
    return 0

opt_map = {
    'Adam': torch.optim.Adam,
    'SGD': torch.optim.SGD,
    'ASGD': torch.optim.ASGD,
    'RMSprop': torch.optim.RMSprop,
    'Adagrad': torch.optim.Adagrad
}
crit_map = {
    'BCELoss': nn.BCELoss(),
    'MSELoss': nn.MSELoss(),
    'L1Loss': nn.L1Loss(),
}
penalty_map  ={
    'entropy': EntropyPenalty,
    'kl': KLPenalty,
    'none': NoPenalty
}

def objective(trial):
'''
Defines the search space for Optuna, exploring various hyperparameters. 
For each trial, it trains the model on a training split and computes the average loss based on the JSD between the predicted and real PWMs. 
'''
    optimizer_name = trial.suggest_categorical('optimizer', list(opt_map.keys()))
    criterion_name = trial.suggest_categorical('criterion', list(crit_map.keys()))
    lr = trial.suggest_categorical('lr', [1e-4, 1e-3, 1e-2, 1e-1])
    drop_out = trial.suggest_categorical('drop_out', [0.1, 0.2, 0.3, 0.4])
    fc_size = trial.suggest_categorical('fc_size', [1024, 512, 256])
    hl_size = trial.suggest_categorical('hl_size', [128, 256, 512])
    penalty_name = trial.suggest_categorical('penalty',['entropy', 'kl', 'none'])
    penalty_strength = trial.suggest_categorical('penalty_strength', [0.05, 0.1, 0.2])
    model_ = trial.suggest_categorical('model', ['FC', 'CNN'])

    if model_ == 'FC':
        model = FC_PWM(fc_size=fc_size, hl_size=hl_size, drop_out=drop_out)
    else:
        model = CNN_PWM(hl_size=hl_size, drop_out=drop_out)

    optimizer = opt_map[optimizer_name](model.parameters(), lr=lr)
    criterion = crit_map[criterion_name]
    penalty = penalty_map[penalty_name]

    model.train()
    total_loss = 0.0
    total_samples = 0

    splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.25, random_state=99)

    for train_idx, test_idx in splitter.split(dataset, subgroups.families):

        for idx in train_idx:

            data, label = dataset[idx]

            input_id, input_tensor_nn, input_tensor_modcre, one_hot, pwm_real = data

            pwm_pred = model(input_tensor_nn, input_tensor_modcre, one_hot).squeeze(0)
            score = JSDScore(pwm_pred, pwm_real)
            target = torch.tensor(float(label))
            loss = criterion(score, target) + penalty(pwm_pred, penalty_strength)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            total_samples += 1

        avg_loss = total_loss / total_samples

    return avg_loss


def TrainEvalBestParams(pwms_dir):
'''
Trains the model using the optimal set of hyperparameters. 
It evaluates performance on a test split, logging scores, losses, and sequence logos for each TF. 
It saves these outputs and stores the final trained model to disk.
'''
    optimizer_name = 'Adam'
    criterion_name = 'MSELoss'
    lr = 0.001
    drop_out = 0.1
    fc_size = 1024
    hl_size = 512
    penalty_name = 'none'
    penalty_strength = 0.05
    model_ = 'CNN'

    if model_ == 'FC':
        model = FC_PWM(fc_size=fc_size, hl_size=hl_size, drop_out=drop_out)
    else:
        model = CNN_PWM(hl_size=hl_size, drop_out=drop_out)

    optimizer = opt_map[optimizer_name](model.parameters(), lr=lr)
    criterion = crit_map[criterion_name]
    penalty = penalty_map[penalty_name]
    model.train()

    splitter = StratifiedShuffleSplit(n_splits=1, test_size=0.25, random_state=99)

    for train_idx, test_idx in splitter.split(dataset, dataset.subgroups):

        train_loss = 0.0
        train_samples = 0
        for idx in train_idx:

            data, label = dataset[idx]

            input_id, input_tensor_nn, input_tensor_modcre, one_hot, pwm_real = data

            pwm_pred = model(input_tensor_nn, input_tensor_modcre, one_hot).squeeze(0)
            score = jsd_alignment(pwm_pred, pwm_real)
            target = torch.tensor(float(label))
            loss = criterion(score, target) + penalty(pwm_pred, penalty_strength)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if label == '1':
                train_loss += loss.item()
                train_samples += 1

        avg_train_loss = train_loss / train_samples

        model.eval()
        output = open('output_eval.txt', 'w')
        output.write('Motif_id\tScore\tLoss\tFamily\tDatabase\n')
        eval_loss = 0
        eval_samples = 0
        score_list = []
        id_list = []

        with torch.no_grad():


            for idx in test_idx:

                data, label = dataset[idx]

                if label == '1':
                    input_id, input_tensor_nn, input_tensor_modcre, one_hot, pwm_real = data

                    pwm_pred = model(input_tensor_nn, input_tensor_modcre, one_hot).squeeze(0)
                    score = jsd_alignment(pwm_pred, pwm_real)
                    score_list.append(score.item())
                    target = torch.tensor(float(label))
                    loss = criterion(score, target) + penalty(pwm_pred, penalty_strength)

                    fam, db, _ = dataset.subgroups[idx].split('_')
                    output.write(f'{input_id}\t{score.item()}\t{loss.item()}\t{fam}\t{db}\n')
                    eval_loss += loss.item()
                    eval_samples += 1
                    if input_id not in set(id_list):
                        CreateLogo(pwm_pred, input_id)
                        TensorToMeme(pwm_pred, input_id, 'iv_memes')
                        GetPvalue(pwms_dir, f'{pwms_dir}/{input_id}', f'iv_memes/{input_id}')
                    id_list.append(input_id)

        avg_eval_loss = eval_loss / eval_samples
        output.close()
    torch.save(model, 'final_model.pth')

    return avg_train_loss, avg_eval_loss, auc
    
def OptunaOptimization(trials):
'''
Launches the hyperparameter search using Optuna for a specified number of trials. 
It then displays the best parameter set and lowest loss achieved, and saves a plot of the optimization history.
'''
    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=trials)

    print('Best value:', study.best_value)
    print('Best params:', study.best_params)

    fig=plot_optimization_history(study).figure
    fig.set_size_inches(12, 6)
    plt.tight_layout()
    fig.savefig("history_optimization_optuna.png",dpi=300)
    
