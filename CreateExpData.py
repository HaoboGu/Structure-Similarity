import numpy as np
import time
import os
import csv
import random

def read_interacts(filename, pos_interacts, neg_interacts):
    # read interaction data
    interaction_file = open(filename)
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        if int(interact_level) == 1:
            pos_interacts[db_id1, db_id2] = int(interact_level)  # use multiple keys
        elif int(interact_level) == 0:
            neg_interacts[db_id1, db_id2] = int(interact_level)  # use multiple keys
        line = interaction_file.readline()
    interaction_file.close()
    return pos_interacts, neg_interacts


def create_data(pos_interacts, neg_interacts):
    pos = []
    neg = []
    all = []
    ids = []
    for key in pos_interacts:
        line = str(key[0]) +'\t'+ str(key[1]) +'\t'+ str(pos_interacts[key])+'\n'
        pos.append(line)
        all.append(line)
        line = str(key[0]) +'\t'+ str(key[1]) +'\n'
        ids.append(line)

    for key in neg_interacts:
        line = str(key[0]) + '\t' + str(key[1]) + '\t' + str(neg_interacts[key]) + '\n'
        neg.append(line)
        all.append(line)
        line = str(key[0]) + '\t' + str(key[1]) + '\n'
        ids.append(line)
    # random.shuffle(all)
    # for item in all:
    #     tup = (item[0], item[1])
    #     ids.append(tup)
    return pos, neg, all, ids

def divide_data(pos_interacts, neg_interacts):
    # pos: 3624 = 362*9(3258)+366; neg: 85394=8539*9(76851)+8543
    pos_keys = list(pos_interacts.keys())
    random.shuffle(pos_keys)
    pos_dicts = []
    for i in range(0, 9):
        p = {pos_key: pos_interacts[pos_key] for pos_key in pos_keys[362*i:362*(i+1)]}
        pos_dicts.append(p)
    p = {pos_key: pos_interacts[pos_key] for pos_key in pos_keys[3258:]}
    pos_dicts.append(p)
    
    neg_keys = list(neg_interacts.keys())
    random.shuffle(neg_keys)
    neg_dicts = []
    for i in range(0, 9):
        p = {neg_key: neg_interacts[neg_key] for neg_key in neg_keys[8539*i:8539*(i+1)]}
        neg_dicts.append(p)
    p = {neg_key: neg_interacts[neg_key] for neg_key in neg_keys[76851:]}
    neg_dicts.append(p)

    return pos_dicts, neg_dicts

def write_interacts(pos_interacts, neg_interacts, path):
    if not os.path.exists(path):
        os.makedirs(path)

    pos, neg, all, ids = create_data(pos_interacts, neg_interacts)
    filename = path + '/interacts.csv'
    with open(filename, 'w') as f:
        f.writelines(all)
    filename = path + '/interacts_negatives.csv'
    with open(filename, 'w') as f:
        f.writelines(neg)
    filename = path + '/interacts_positives.csv'
    with open(filename, 'w') as f:
        f.writelines(pos)
    filename = path + '/interactsids.csv'
    with open(filename, 'w') as f:
        f.writelines(ids)





##########################  main  #####################################
pos_interacts = {}
neg_interacts = {}
for i in range(1, 11):
    filename = 'data/all/' + str(i) + '/interacts.csv'
    pos_interacts, neg_interacts = read_interacts(filename, pos_interacts, neg_interacts)


# divide


for j in range(0, 30):
    pos_dict, neg_dict = divide_data(pos_interacts, neg_interacts)
    for i in range(1, 11):
        # ten folds
        path = 'data/traningdata/all' + str(j) + '/' + str(i)
        write_interacts(pos_dict[i-1], neg_dict[i-1], path)



