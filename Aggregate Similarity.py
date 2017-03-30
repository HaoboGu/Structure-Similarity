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

def read_interacts(filename):
    # read interaction data
    interaction_file = open(filename)
    interact_dict = {}
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        if int(interact_level) == 1:
            interact_dict[db_id1, db_id2] = int(interact_level)  # use multiple keys
        line = interaction_file.readline()
    interaction_file.close()
    return interact_dict

def select_pairs(interacts, sim_keys):
    # select interacted/non-interacted pairs as train set

    num = 0
    num_non_inter = 0
    train_pair = {}
    for key in interacts:
        train_pair[key] = 1
        num += 1

    # the number of non-interacted pairs should be same as interacted pairs
    for key in sim_keys:
        if key not in train_pair.keys():
            if num_non_inter < num:
                    train_pair[key] = 0
                    num_non_inter += 1
            else:
                break
    return train_pair

def create_avg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    true_value = []
    train_set = []
    for id1, id2 in samples.keys():  # test drug a and c
        sa, sc, sd, sg, ss = [0, 0, 0, 0, 0]
        na, nc, nd, ng, ns = [0, 0, 0, 0, 0]
        for key in sim_keys:  # key: (a, b)
            if (id1 in key and id1 != key[1]) and (id2, key[1]) in interacts:
                # "id1 in key and id1 != key[1]" means key[0] is a(id1), key[1] is b
                # "(id2,key[1]) in interacts" means b-c interact
                # 此处如果添加一个sim值的限制的话，会导致每个sim参与平均的n的不同，否则都相同
                # n=0意味着没有满足条件的b，使得ab有sim，bc反应
                thres = 0.4
                if atc_sim[key] > thres:
                    sa += atc_sim[key]
                    na += 1
                if chemical_sim[key] > thres:
                    sc += chemical_sim[key]
                    nc += 1
                if dist_sim[key] > thres:
                    sd += dist_sim[key]
                    nd += 1
                if go_sim[key] > thres:
                    sg += go_sim[key]
                    ng += 1
                if seq_sim[key] > thres:
                    ss += seq_sim[key]
                    ns += 1
        if na != 0:
            fa = sa / na
        else:
            fa = 0
        if nc != 0:
            fc = sc / nc
        else:
            fc = 0
        if nd != 0:
            fd = sd / nd
        else:
            fd = 0
        if ng != 0:
            fg = sg / ng
        else:
            fg = 0
        if ns != 0:
            fs = ss / ns
        else:
            fs = 0
        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    # for i in range(0, feature_atc.__len__()):
    #     print(feature_atc[i], feature_chem[i], feature_dist[i], feature_go[i], feature_seq[i], true_value[i])
    return train_set

def create_max_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    true_value = []
    train_set = []
    for id1, id2 in samples.keys():  # test drug a and c
        ma, mc, md, mg, ms = [0, 0, 0, 0, 0]
        for key in sim_keys:  # key: (a, b)
            if (id1 in key and id1 != key[1]) and (id2, key[1]) in interacts:
                # "id1 in key and id1 != key[1]" means key[0] is a(id1), key[1] is b
                # "(id2,key[1]) in interacts" means b-c interact
                thres = 0.4
                if atc_sim[key] > ma and atc_sim[key] > thres:
                    ma = atc_sim[key]
                if chemical_sim[key] > mc and chemical_sim[key] > thres:
                    mc = chemical_sim[key]
                if dist_sim[key] > md and dist_sim[key] > thres:
                    md = dist_sim[key]
                if go_sim[key] > mg and go_sim[key] > thres:
                    mg = go_sim[key]
                if seq_sim[key] > ms and seq_sim[key] > thres:
                    ms = seq_sim[key]
        feature_atc.append(ma)
        feature_chem.append(mc)
        feature_dist.append(md)
        feature_go.append(mg)
        feature_seq.append(ms)
        true_value.append(samples[id1, id2])
        train_set.append([ma, mc, md, mg, ms, samples[id1, id2]])

    return train_set

## test functions ##   atc, chemical, dist, go, seq
start = time.time()

for i in range(1,11):
    filename = 'data/interacts_r' + str(i) + '.csv'
    atc_sim, chemical_sim, dist_sim, go_sim, ligand_sim, seq_sim, sideeffect_sim = read_similarities()
    interacts = read_interacts(filename)
    sim_keys = (chemical_sim.keys() & atc_sim.keys())  # keys of chem, dist, go, seq are same
    samples = select_pairs(interacts, sim_keys)
    # print(samples)
    train_set = create_max_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim)
    train = np.array(train_set)
    target = train.transpose()[5]
    train = train.transpose()[0:5]
    train = train.transpose()
    # print(train)
    logistic = lr(random_state=1)
    logistic.fit(train, target)
    print(logistic.coef_)
    end = time.time()
    print(end-start)

