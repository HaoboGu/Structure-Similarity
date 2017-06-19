import os
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
def read_all_interacts(filename, interacts, interact_lst):
    interaction_file = open(filename)
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        interacts[db_id1, db_id2] = int(interact_level)
        interact_lst.append((db_id1, db_id2, interact_level))
        line = interaction_file.readline()
    interaction_file.close()
    return interacts, interact_lst


def divide_data(pos_interacts, neg_interacts):  # divide pos_interacts and neg_interacts to 10 folds
    # pos: 3978 = 397*9(3573)+405; neg: 94932 = 9493*9(85437)+9495

    # divide positive interacts
    pos_keys = list(pos_interacts.keys())
    random.shuffle(pos_keys)  # shuffle positive interacts
    pos_dicts = []
    # save first 9 folds
    for i in range(0, 9):
        p = {pos_key: pos_interacts[pos_key] for pos_key in pos_keys[397*i:397*(i+1)]}
        pos_dicts.append(p)
    # save the rest to 10-th fold
    p = {pos_key: pos_interacts[pos_key] for pos_key in pos_keys[3573:]}
    pos_dicts.append(p)

    # divide negative interacts
    neg_keys = list(neg_interacts.keys())
    random.shuffle(neg_keys)
    neg_dicts = []
    for i in range(0, 9):
        p = {neg_key: neg_interacts[neg_key] for neg_key in neg_keys[9493*i:9493*(i+1)]}
        neg_dicts.append(p)
    p = {neg_key: neg_interacts[neg_key] for neg_key in neg_keys[85437:]}
    neg_dicts.append(p)

    return pos_dicts, neg_dicts
def divide_data2(interacts):
    # divide all interacts to 10 folds
    keys = list(interacts.keys())
    random.shuffle(keys)  # shuffle positive interacts
    interact_dicts = []
    # 5 folds with 9892 drugs
    for i in range(0, 5):
        p = {key: interacts[key] for key in keys[9892 * i:9892 * (i + 1)]}
        interact_dicts.append(p)
    # 5 folds with 9890 drugs
    for i in range(0, 5):
        p = {key: interacts[key] for key in keys[(49460+9890 * i):(49460+9890 * (i + 1))]}
        interact_dicts.append(p)
    s= 0
    for item in interact_dicts:
        s = s + item.__len__()
    print(s)
    return interact_dicts


def create_data(pos_interacts, neg_interacts):
    pos = []
    neg = []
    all = []
    ids = []
    # read all interacts
    for key in pos_interacts:
        line = str(key[0]) +'\t'+ str(key[1]) +'\t'+ str(pos_interacts[key])+'\n'
        all.append(line)

    for key in neg_interacts:
        line = str(key[0]) + '\t' + str(key[1]) + '\t' + str(neg_interacts[key]) + '\n'
        all.append(line)

    random.shuffle(all)  # shuffle interacts

    # save pos/neg interacts
    for index in range(0, all.__len__()):
        id1, id2, interacts = all[index].split('\t')
        line = id1 + '\t' + id2 + '\t' + interacts
        if interacts == '1\n':
            pos.append(line)
        elif interacts == '0\n':
            neg.append(line)
        line = id1 + '\t' + id2 + '\n'
        ids.append(line)
    return pos, neg, all, ids
def create_data2(interact_dict):
    pos = []
    neg = []
    all = []
    ids = []
    # read all interacts
    for key in interact_dict:
        line = str(key[0]) + '\t' + str(key[1]) + '\t' + str(interact_dict[key]) + '\n'
        all.append(line)
        if interact_dict[key] == 1:
            pos.append(line)
        elif interact_dict[key] == 0:
            neg.append(line)
        line = str(key[0]) + '\t' + str(key[1]) + '\n'
        ids.append(line)
    return pos, neg, all, ids


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
def write_interacts2(interact_dicts, path):
    if not os.path.exists(path):
        os.makedirs(path)
    pos, neg, all, ids = create_data2(interact_dicts)
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


# #### 先生成一个数据集
# interacts = {}
# interact_lst = []
# l = 0
# last = 0
# for i in range(1, 11):
#     filename = 'data/all/' + str(i) + '/interacts.csv'
#     # pos_interacts, neg_interacts = read_interacts(filename, pos_interacts, neg_interacts)
#     print(filename)
#     interacts, interact_lst = read_all_interacts(filename, interacts, interact_lst)
# # divide data
# lst = []
# for key in interacts:
#     lst.append((key[0],key[1],interacts[key]))
# # random.shuffle(lst)  #
# samples_each_fold = 9891
# interact_dicts = []
# for i in range(0, 10):
#     fold = lst[i*samples_each_fold:samples_each_fold*(i+1)]
#     interact_dicts.append(fold)
#
# # interact_lists = divide_data2(interacts)
# j=2
# for i in range(1, 11):
#     interact_dict = interact_dicts[i-1]
#     # ten folds
#     path = 'data/all_dataset' + str(j) + '/all/' + str(i)
#     # write_interacts(pos_dict[i-1], neg_dict[i-1], path)
#     if not os.path.exists(path):
#         os.makedirs(path)
#
#     pos = []
#     neg = []
#     all = []
#     ids = []
#     # read all interacts
#     for item in interact_dict:
#         line = item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\n'
#         all.append(line)
#         if item[2] == 1:
#             pos.append(line)
#         elif item[2] == 0:
#             neg.append(line)
#         line = item[0] + '\t' + item[1] + '\n'
#         ids.append(line)
#
#     filename = path + '/interacts.csv'
#     with open(filename, 'w') as f:
#         f.writelines(all)
#     filename = path + '/interacts_negatives.csv'
#     with open(filename, 'w') as f:
#         f.writelines(neg)
#     filename = path + '/interacts_positives.csv'
#     with open(filename, 'w') as f:
#         f.writelines(pos)
#     filename = path + '/interactsids.csv'
#     with open(filename, 'w') as f:
#         f.writelines(ids)
    ####



#########################  main  #####################################
pos_interacts = {}
neg_interacts = {}
interacts = {}
l=0
last = 0
for i in range(1, 11):
    filename = 'data/all/' + str(i) + '/interacts.csv'
    # pos_interacts, neg_interacts = read_interacts(filename, pos_interacts, neg_interacts)
    print(filename)
    interacts = read_all_interacts(filename, interacts)


# divide


# 30 datasets
for j in range(1, 31):
    # pos_dict, neg_dict = divide_data(pos_interacts, neg_interacts)
    # j=1
    interact_lists = divide_data2(interacts)
    for i in range(1, 11):
        # ten folds
        path = 'data/all_dataset' + str(j) + '/all/' + str(i)
        # write_interacts(pos_dict[i-1], neg_dict[i-1], path)
        write_interacts2(interact_lists[i-1], path)



