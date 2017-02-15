import numpy as np
import time
from sklearn.linear_model import LogisticRegression as lr


def read_similarities():  # all keys are order sensitive
    chemical_sim = open('data/chemicalsimilarity.csv')
    chemical_sim_dict = {}
    line = chemical_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        chemical_sim_dict[db_id1, db_id2] = float(similarity)
        line = chemical_sim.readline()
    chemical_sim.close()

    atc_sim = open('data/atcsimilarity.csv')
    atc_sim_dict = {}
    line = atc_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        atc_sim_dict[db_id1, db_id2] = float(similarity)
        line = atc_sim.readline()
    atc_sim.close()

    go_sim = open('data/gosimilarity.csv')
    go_sim_dict = {}
    line = go_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        go_sim_dict[db_id1, db_id2] = float(similarity)
        line = go_sim.readline()
    go_sim.close()

    dist_sim = open('data/distsimilarity.csv')
    dist_sim_dict = {}
    line = dist_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        dist_sim_dict[db_id1, db_id2] = float(similarity)
        line = dist_sim.readline()
    dist_sim.close()

    ligand_sim = open('data/ligandsimilarity.csv')
    ligand_sim_dict = {}
    line = ligand_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        ligand_sim_dict[db_id1, db_id2] = float(similarity)
        line = ligand_sim.readline()
    ligand_sim.close()

    seq_sim = open('data/seqsimilarity.csv')
    seq_sim_dict = {}
    line = seq_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        seq_sim_dict[db_id1, db_id2] = float(similarity)
        line = seq_sim.readline()
    seq_sim.close()

    sideeffect_sim = open('data/sideeffectsimilarity.csv')
    sideeffect_sim_dict = {}
    line = sideeffect_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        sideeffect_sim_dict[db_id1, db_id2] = float(similarity)
        line = sideeffect_sim.readline()
    sideeffect_sim.close()

    common_key = (chemical_sim_dict.keys() & dist_sim_dict.keys() & go_sim_dict.keys() & seq_sim_dict.keys())
    # print("num of common key: ", common_key.__len__())
    chemical_sim_dict = {key: chemical_sim_dict[key] for key in common_key}
    dist_sim_dict = {key: dist_sim_dict[key] for key in common_key}
    go_sim_dict = {key: go_sim_dict[key] for key in common_key}
    seq_sim_dict = {key: seq_sim_dict[key] for key in common_key}

    return atc_sim_dict, chemical_sim_dict, dist_sim_dict, go_sim_dict, ligand_sim_dict, seq_sim_dict, sideeffect_sim_dict


def read_interacts():
    # read interaction data
    interaction_file = open('data/interacts.csv')
    interact_dict = {}
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        if int(interact_level) == 1:
            interact_dict[db_id1, db_id2] = int(interact_level)  # use multiple keys
        line = interaction_file.readline()
    interaction_file.close()
    return interact_dict


def select_pairs(interacts, common_key):
    # select interacted/non-interacted pairs as train set

    # first remove interacts without similarities

    for key in interacts.copy():  # use copy to avoid RuntimeError: dictionary changed size during iteration
        if key not in common_key:
            interacts.pop(key)

    #
    num = 0
    num_non_inter = 0
    train_pair = {}
    for key in interacts:
        train_pair[key] = 1
        num += 1

    # the number of non-interacted pairs should be same as interacted pairs
    sim_key = (chemical_sim.keys() & dist_sim.keys() & go_sim.keys() & seq_sim.keys())
    for key in sim_key:
        if num_non_inter < num:
                train_pair[key] = 0
                num_non_inter += 1
        else:
            break

    # num_atc, num_chem, num_dist, num_go, num_lig, num_seq, num_side = [0, 0, 0, 0, 0, 0, 0]
    # for key in train_pair:
    #     if train_pair[key] == 1:
    #         if key not in atc_sim.keys():
    #             num_atc += 1  # out
    #         if key not in chemical_sim.keys():
    #             num_chem += 1
    #         if key not in dist_sim.keys():
    #             num_dist += 1
    #         if key not in go_sim.keys():
    #             num_go += 1
    #         if key not in ligand_sim.keys():
    #             num_lig += 1  # out
    #         if key not in seq_sim.keys():
    #             num_seq += 1
    #         if key not in sideeffect_sim.keys():
    #             num_side += 1  # out
    # print(num_atc, num_chem, num_dist, num_go, num_lig, num_seq, num_side)
    return train_pair


def aggregate_similarities(chemical_sim, dist_sim, go_sim, seq_sim, weight):
    weighted_sim = {}
    for key in chemical_sim:
        weighted_sim[key] = weight[0]*chemical_sim[key] + weight[1]*dist_sim[key] + \
                            weight[2]*go_sim[key] + weight[3]*seq_sim[key]
    return weighted_sim


def log_reg(weighted_sim, samples):
    num = 0
    train_set = []
    for key in samples:

        id1 = key[0]
        candidate = []  # a list of candidate pairs ab1~abn

        for simkey in weighted_sim:
            if id1 in simkey and weighted_sim[simkey] > 0.2:
                candidate.append(simkey)

        max_sim = 0
        # max_sim_key = 0
        for interkey in samples:
            if interkey != key:
                if id1 in interkey and weighted_sim[interkey]>max_sim:
                    # max_sim_key = interkey
                    max_sim = weighted_sim[interkey]
        # print(max_sim)
        prob = 0
        num_cand = 0
        for item in candidate:
            if item in samples:
                prob += weighted_sim[item]
                num_cand += 1
        # if max_sim is not 0:
        #     print(max_sim)
        if num_cand is not 0:
            prob = float(prob) / num_cand
        train_set.append([prob, max_sim, samples[key]])
        # print(num)
        num += 1
    train = np.array(train_set)
    target = train.transpose()[2]
    train = train.transpose()[0:2]
    train = train.transpose()
    logistic = lr()
    logistic.fit(train, target)
    re = logistic.predict(train)

    correct = 0
    for index in range(0, re.__len__()):
        if re[index] == target[index]:
            correct += 1
    return correct

## test functions ##
start = time.time()
atc_sim, chemical_sim, dist_sim, go_sim, ligand_sim, seq_sim, sideeffect_sim = read_similarities()
interacts = read_interacts()
common_key = (chemical_sim.keys() & dist_sim.keys() & go_sim.keys() & seq_sim.keys() & interacts.keys())
samples = select_pairs(interacts, common_key)

start = time.time()
correct = 0
num = 1
step = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
results = []
for w1 in step:
    w4 = 0.5 - w1
    for w2 in step:
        w3 = 0.5 - w2
        weight = [w1, w2, w3, w4]
        weighted_sim = aggregate_similarities(chemical_sim, dist_sim, go_sim, seq_sim, weight)
        tmp = log_reg(weighted_sim, samples)
        if tmp > correct:
            correct = tmp
            results.append([weight, tmp])
            optimized_w = weight
        print('loop ', num, '....')
        num += 1


print(correct)
end = time.time()
print(end-start)

result1 = [[0.4, 0.1, 0.2, 0.3], 444] #12,34
result2 = [[0.45, 0.2, 0.05, 0.3], 447] #13,24
result3 = [[0.15, 0.05, 0.45, 0.35], 456] #14, 23

# baseline = [0.25, 0.25, 0.25, 0.25], fold1:0.4696/0.9159; fold2:0.4705/0.9047; fold3:0.5243/0.9325
