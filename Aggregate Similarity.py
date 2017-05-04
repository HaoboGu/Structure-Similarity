import numpy as np
import time
import math
from sklearn.linear_model import LogisticRegression as lr
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score as aupr_score
from sklearn.preprocessing import normalize
from sklearn.linear_model import LogisticRegressionCV as lrcv
import sklearn.feature_selection as fs
import matplotlib.pyplot as plt

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

def create_avg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    true_value = []
    train_set = []
    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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
        dict_atc[(id1, id2)] = fa
        dict_chem[(id1, id2)] = fc
        dict_dist[(id1, id2)] = fd
        dict_go[(id1, id2)] = fg
        dict_seq[(id1, id2)] = fs

        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    # for i in range(0, feature_atc.__len__()):
    #     print(feature_atc[i], feature_chem[i], feature_dist[i], feature_go[i], feature_seq[i], true_value[i])
    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_avg_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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
                    # sa += atc_sim[key]
                    sa += 2*atc_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    na += 1
                if chemical_sim[key] > thres:
                    sc += 2*chemical_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    nc += 1
                if dist_sim[key] > thres:
                    sd += 2*dist_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    nd += 1
                if go_sim[key] > thres:
                    sg += 2*go_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    ng += 1
                if seq_sim[key] > thres:
                    ss += 2*seq_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
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

        dict_atc[(id1, id2)] = fa
        dict_chem[(id1, id2)] = fc
        dict_dist[(id1, id2)] = fd
        dict_go[(id1, id2)] = fg
        dict_seq[(id1, id2)] = fs

        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    # for i in range(0, feature_atc.__len__()):
    #     print(feature_atc[i], feature_chem[i], feature_dist[i], feature_go[i], feature_seq[i], true_value[i])
    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_sum_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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

        dict_atc[(id1, id2)] = sa
        dict_chem[(id1, id2)] = sc
        dict_dist[(id1, id2)] = sd
        dict_go[(id1, id2)] = sg
        dict_seq[(id1, id2)] = ss

        feature_atc.append(sa)
        feature_chem.append(sc)
        feature_dist.append(sd)
        feature_go.append(sg)
        feature_seq.append(ss)
        true_value.append(samples[id1, id2])
        train_set.append([sa,sc,sd,sg,ss,samples[id1,id2]])
    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_sum_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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
                    sa += 2*atc_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    na += 1
                if chemical_sim[key] > thres:
                    sc += 2*chemical_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    nc += 1
                if dist_sim[key] > thres:
                    sd += 2*dist_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    nd += 1
                if go_sim[key] > thres:
                    sg += 2*go_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    ng += 1
                if seq_sim[key] > thres:
                    ss += 2*seq_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    ns += 1

        dict_atc[(id1, id2)] = sa
        dict_chem[(id1, id2)] = sc
        dict_dist[(id1, id2)] = sd
        dict_go[(id1, id2)] = sg
        dict_seq[(id1, id2)] = ss

        feature_atc.append(sa)
        feature_chem.append(sc)
        feature_dist.append(sd)
        feature_go.append(sg)
        feature_seq.append(ss)
        true_value.append(samples[id1, id2])
        train_set.append([sa,sc,sd,sg,ss,samples[id1,id2]])
    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_max_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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

        dict_atc[(id1, id2)] = ma
        dict_chem[(id1, id2)] = mc
        dict_dist[(id1, id2)] = md
        dict_go[(id1, id2)] = mg
        dict_seq[(id1, id2)] = ms

        feature_atc.append(ma)
        feature_chem.append(mc)
        feature_dist.append(md)
        feature_go.append(mg)
        feature_seq.append(ms)
        true_value.append(samples[id1, id2])
        train_set.append([ma, mc, md, mg, ms, samples[id1, id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_max_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}

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
                    ma = 2*atc_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    # print('sim:', atc_sim[key], ',degree:', degree[key[1]], ', avg_deg:', avg_deg, ', feature:', ma)
                if chemical_sim[key] > mc and chemical_sim[key] > thres:
                    mc = 2*chemical_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if dist_sim[key] > md and dist_sim[key] > thres:
                    md = 2*dist_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if go_sim[key] > mg and go_sim[key] > thres:
                    mg = 2*go_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if seq_sim[key] > ms and seq_sim[key] > thres:
                    ms = 2*seq_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))

        dict_atc[(id1, id2)] = ma
        dict_chem[(id1, id2)] = mc
        dict_dist[(id1, id2)] = md
        dict_go[(id1, id2)] = mg
        dict_seq[(id1, id2)] = ms

        feature_atc.append(ma)
        feature_chem.append(mc)
        feature_dist.append(md)
        feature_go.append(mg)
        feature_seq.append(ms)
        true_value.append(samples[id1, id2])
        train_set.append([ma, mc, md, mg, ms, samples[id1, id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_min_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}
    true_value = []
    train_set = []
    for id1, id2 in samples.keys():  # test drug a and c
        ma, mc, md, mg, ms = [1, 1, 1, 1, 1]
        for key in sim_keys:  # key: (a, b)
            if (id1 in key and id1 != key[1]) and (id2, key[1]) in interacts:
                # "id1 in key and id1 != key[1]" means key[0] is a(id1), key[1] is b
                # "(id2,key[1]) in interacts" means b-c interact
                thres = 0.4
                if atc_sim[key] < ma and atc_sim[key] > thres:
                    ma = atc_sim[key]
                if chemical_sim[key] < mc and chemical_sim[key] > thres:
                    mc = chemical_sim[key]
                if dist_sim[key] < md and dist_sim[key] > thres:
                    md = dist_sim[key]
                if go_sim[key] < mg and go_sim[key] > thres:
                    mg = go_sim[key]
                if seq_sim[key] < ms and seq_sim[key] > thres:
                    ms = seq_sim[key]

        dict_atc[(id1, id2)] = ma
        dict_chem[(id1, id2)] = mc
        dict_dist[(id1, id2)] = md
        dict_go[(id1, id2)] = mg
        dict_seq[(id1, id2)] = ms
        feature_atc.append(ma)
        feature_chem.append(mc)
        feature_dist.append(md)
        feature_go.append(mg)
        feature_seq.append(ms)
        true_value.append(samples[id1, id2])
        train_set.append([ma, mc, md, mg, ms, samples[id1, id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_min_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}
    true_value = []
    train_set = []
    for id1, id2 in samples.keys():  # test drug a and c
        ma, mc, md, mg, ms = [0, 0, 0, 0, 0]
        for key in sim_keys:  # key: (a, b)
            if (id1 in key and id1 != key[1]) and (id2, key[1]) in interacts:
                # "id1 in key and id1 != key[1]" means key[0] is a(id1), key[1] is b
                # "(id2,key[1]) in interacts" means b-c interact
                thres = 0.4
                if atc_sim[key] < ma and atc_sim[key] > thres:
                    ma = 2*atc_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                    # print('sim:', atc_sim[key], ',degree:', degree[key[1]], ', avg_deg:', avg_deg, ', feature:', ma)
                if chemical_sim[key] < mc and chemical_sim[key] > thres:
                    mc = 2*chemical_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if dist_sim[key] < md and dist_sim[key] > thres:
                    md = 2*dist_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if go_sim[key] < mg and go_sim[key] > thres:
                    mg = 2*go_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))
                if seq_sim[key] < ms and seq_sim[key] > thres:
                    ms = 2*seq_sim[key]/(1+math.exp(-((degree[key[1]]-avg_deg)/avg_deg)))

        dict_atc[(id1, id2)] = ma
        dict_chem[(id1, id2)] = mc
        dict_dist[(id1, id2)] = md
        dict_go[(id1, id2)] = mg
        dict_seq[(id1, id2)] = ms
        feature_atc.append(ma)
        feature_chem.append(mc)
        feature_dist.append(md)
        feature_go.append(mg)
        feature_seq.append(ms)
        true_value.append(samples[id1, id2])
        train_set.append([ma, mc, md, mg, ms, samples[id1, id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_avdeg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []
    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}
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
                    # sa += atc_sim[key]
                    sa += degree[key[1]]
                    na += 1
                if chemical_sim[key] > thres:
                    sc += degree[key[1]]
                    nc += 1
                if dist_sim[key] > thres:
                    sd += degree[key[1]]
                    nd += 1
                if go_sim[key] > thres:
                    sg += degree[key[1]]
                    ng += 1
                if seq_sim[key] > thres:
                    ss += degree[key[1]]
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
        dict_atc[(id1, id2)] = fa
        dict_chem[(id1, id2)] = fc
        dict_dist[(id1, id2)] = fd
        dict_go[(id1, id2)] = fg
        dict_seq[(id1, id2)] = fs
        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    # for i in range(0, feature_atc.__len__()):
    #     print(feature_atc[i], feature_chem[i], feature_dist[i], feature_go[i], feature_seq[i], true_value[i])
    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}
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
                    # sa += atc_sim[key]
                    na += 1
                if chemical_sim[key] > thres:
                    nc += 1
                if dist_sim[key] > thres:
                    nd += 1
                if go_sim[key] > thres:
                    ng += 1
                if seq_sim[key] > thres:
                    ns += 1
        if na != 0:
            fa = degree[key[0]] / avg_deg
        else:
            fa = 0
        if nc != 0:
            fc = degree[key[0]] / avg_deg
        else:
            fc = 0
        if nd != 0:
            fd = degree[key[0]] / avg_deg
        else:
            fd = 0
        if ng != 0:
            fg = degree[key[0]] / avg_deg
        else:
            fg = 0
        if ns != 0:
            fs = degree[key[0]] / avg_deg
        else:
            fs = 0

        dict_atc[(id1, id2)] = fa
        dict_chem[(id1, id2)] = fc
        dict_dist[(id1, id2)] = fd
        dict_go[(id1, id2)] = fg
        dict_seq[(id1, id2)] = fs
        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def create_ratio_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg):
    feature_atc = []
    feature_chem = []
    feature_dist = []
    feature_seq = []
    feature_go = []

    dict_atc = {}
    dict_chem = {}
    dict_dist = {}
    dict_go = {}
    dict_seq = {}
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
                    # sa += atc_sim[key]
                    na += 1
                if chemical_sim[key] > thres:
                    nc += 1
                if dist_sim[key] > thres:
                    nd += 1
                if go_sim[key] > thres:
                    ng += 1
                if seq_sim[key] > thres:
                    ns += 1
        if na != 0:
            fa = degree[key[0]] / na
        else:
            fa = 0
        if nc != 0:
            fc = degree[key[0]] / nc
        else:
            fc = 0
        if nd != 0:
            fd = degree[key[0]] / nd
        else:
            fd = 0
        if ng != 0:
            fg = degree[key[0]] / ng
        else:
            fg = 0
        if ns != 0:
            fs = degree[key[0]] / ns
        else:
            fs = 0

        dict_atc[(id1, id2)] = fa
        dict_chem[(id1, id2)] = fc
        dict_dist[(id1, id2)] = fd
        dict_go[(id1, id2)] = fg
        dict_seq[(id1, id2)] = fs
        feature_atc.append(fa)
        feature_chem.append(fc)
        feature_dist.append(fd)
        feature_go.append(fg)
        feature_seq.append(fs)
        true_value.append(samples[id1, id2])
        train_set.append([fa,fc,fd,fg,fs,samples[id1,id2]])

    return train_set, dict_atc, dict_chem, dict_dist, dict_go, dict_seq

def compute_avg_deg(interacts):
    degree = {}
    for id1, id2 in interacts:
        if id1 not in degree.keys():
            degree[id1] = 1
        else:
            degree[id1] += 1

        if id2 not in degree.keys():
            degree[id2] = 1
        else:
            degree[id2] += 1
    for key in degree:
        degree[key] = int(degree[key]/2)
    l = list(degree.values())
    avg_deg = int(sum(l)/len(l))
    # print(avg_deg)
    return degree, avg_deg

def norm(dict):
    max_val = max(dict.values())
    # print('max: ', max_val)
    for key in dict:
        dict[key] = dict[key]/max_val
    return dict

def save_data(dict, filename):
    file = open(filename, 'w')
    norm(dict)
    for key in dict:
        line = str(key[0]) + '\t' + str(key[1]) + '\t' + str(dict[key]) + '\n'
        file.write(line)
    file.close()

start = time.time()
filename = 'data/interacts_all.csv'
atc_sim, chemical_sim, dist_sim, go_sim, ligand_sim, seq_sim, sideeffect_sim = read_similarities()
interacts = read_interacts(filename)
sim_keys = (chemical_sim.keys() & atc_sim.keys())  # keys of chem, dist, go, seq are same
samples = select_pairs(interacts, sim_keys)
degree, avg_deg = compute_avg_deg(samples)
train_set, avg_atc, avg_chem, avg_dist, avg_go, avg_seq = create_avg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts)
train_set2, max_atc, max_chem, max_dist, max_go, max_seq = create_max_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts)
train_set3, avdeg_atc, avdeg_chem, avdeg_dist, avdeg_go, avdeg_seq = create_avdeg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)
train_set4, max_deg_atc, max_deg_chem, max_deg_dist, max_deg_go, max_deg_seq = create_max_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)
train_set5, avg_deg_atc, avg_deg_chem, avg_deg_dist, avg_deg_go, avg_deg_seq = create_avg_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)
train_set6, sum_atc, sum_chem, sum_dist, sum_go, sum_seq = create_sum_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts)
train_set7, sum_deg_atc, sum_deg_chem, sum_deg_dist, sum_deg_go, sum_deg_seq = create_sum_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)
train_set8, min_atc, min_chem, min_dist, min_go, min_seq = create_min_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts)
train_set9, ratio_atc, ratio_chem, ratio_dist, ratio_go, ratio_seq = create_ratio_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)
train_set10, deg_atc, deg_chem, deg_dist, deg_go, deg_seq = create_deg_features(sim_keys, samples, atc_sim, chemical_sim, dist_sim, go_sim, seq_sim, interacts, degree, avg_deg)


train1 = np.array(train_set)
target = train1.transpose()[5]
# train1 = train1.transpose()[0:5]  # train1: avg_sim
train1 = train1.transpose()[0:1]
train1 = normalize(train1, 'max')

train2 = np.array(train_set2)
# train2 = train2.transpose()[0:5]  # train2: max_sim
train2 = train2.transpose()[0:2]
train2 = normalize(train2, 'max')

train3 = np.array(train_set3)
# train3 = train3.transpose()[0:5]  # train3: avg_degree
train31 = train3.transpose()[0:2]
train32 = train3.transpose()[4:5]
train3 = np.concatenate((train31, train32), 0)
train3 = normalize(train3, 'max')

train4 = np.array(train_set4)
train4 = train4.transpose()[0:5]  # train4: max_sim*deg
train4 = normalize(train4, 'max')

train5 = np.array(train_set5)
# train5 = train5.transpose()[0:5]  # train5: avg_sim*deg
train5 = train5.transpose()[0:1]
train5 = normalize(train5, 'max')

train6 = np.array(train_set6)
# train6 = train6.transpose()[0:5]  # train6: sum_sim
train6 = train6.transpose()[0:1]
train6 = normalize(train6, 'max')

train7 = np.array(train_set7)
# train7 = train7.transpose()[0:5]  # train7: sum_(sim*deg)
train71 = train7.transpose()[0:1]
train72 = train7.transpose()[2:4]
train7 = np.concatenate((train71, train72), 0)
train7 = normalize(train7, 'max')

train8 = np.array(train_set8)
# train8 = train8.transpose()[0:5]  # train8: min_sim
train81 = train8.transpose()[0:2]
train82 = train8.transpose()[4:5]
train8 = np.concatenate((train81, train82), 0)
train8 = normalize(train8, 'max')

train9 = np.array(train_set9)
# train9 = train9.transpose()[0:5]  # train9: ratio(degree_a/num_b)
train9 = train9.transpose()[4:5]
train9 = normalize(train9, 'max')

train10 = np.array(train_set10)
# train10 = train10.transpose()[0:5]  # train10: degree
train10 = train10.transpose()[2:5]  # train10: degree
train10 = normalize(train10, 'max')

t = np.concatenate((train1, train2, train3, train5, train6, train7, train8, train9, train10),0)
t = t.transpose()
train = t

# maxCs = []
# for kk in range(1,51):
#     train_new = fs.SelectKBest(fs.f_classif, k=kk).fit_transform(t, target)
#     train = train_new
C = [10, 100]
l = lrcv(Cs=C, cv=10, scoring='average_precision', penalty='l2', random_state=1, solver='liblinear', n_jobs=-1)
l.fit(train, target)
re = l.predict(train)
a = l.scores_[1].transpose()
print('average scores for Cs of 10 folds: ', [a[i].mean() for i in range(0, 2)])
# print('scores: ', l.scores_)
print('coefs for Cs: ', l.coef_)
# for i in range(0, 50):
#     print(float(l.coef_[0][i]))
# print('k=', kk, ' with best score ', max(a[i].mean() for i in range(0, 8)))
    # maxCs.append(max(a[i].mean() for i in range(0, 8)))

end = time.time()
print(end-start)

########################### Save data ###############################
save_data(avg_atc, 'avg_atc.csv')
save_data(avg_go, 'avg_go.csv')
save_data(max_atc, 'max_atc.csv')
save_data(max_chem, 'max_chem.csv')
save_data(max_dist, 'max_dist.csv')
save_data(avdeg_atc, 'avdeg_atc.csv')
save_data(avdeg_chem, 'avdeg_chem.csv')
save_data(avdeg_go, 'avdeg_go.csv')
save_data(avdeg_seq, 'avdeg_seq.csv')
save_data(avg_deg_atc, 'avg_deg_atc.csv')
save_data(sum_atc, 'sum_atc.csv')
save_data(sum_deg_atc, 'sum_deg_atc.csv')
save_data(sum_deg_dist, 'sum_deg_dist.csv')
save_data(sum_deg_go, 'sum_deg_go.csv')
save_data(min_atc, 'min_atc.csv')
save_data(min_chem, 'min_chem.csv')
save_data(min_dist, 'min_dist.csv')
save_data(min_seq, 'min_seq.csv')
save_data(ratio_atc, 'ratio_atc.csv')
save_data(ratio_chem, 'ratio_chem.csv')
save_data(ratio_seq, 'ratio_seq.csv')
save_data(deg_dist, 'deg_dist.csv')
save_data(deg_go, 'deg_go.csv')
save_data(deg_dist, 'deg_dist.csv')
print('over')
#####################################################################

# plt.interactive(False)
# ax = range(1,51)
# plt.plot(ax, maxCs, 'r.')
# plt.show()

# logistic = lr(penalty='l1', C=1.0, random_state=1, solver='liblinear', n_jobs=-1)
# logistic.fit(train, target)
# re = logistic.predict(train)
# auroc = roc_auc_score(target, re)
# aupr = aupr_score(target, re)
#
# print('roc: ', auroc)
# print('aupr: ', aupr)
# print(logistic.coef_)


############################# save integrated similarity######################################
integrated_sim = {}
for key in avg_atc:
    integrated_sim[key] = 1.6*avg_atc[key] + 1.91*avg_chem[key] + 0.23*avg_dist[key] + 0.68*avg_go[key] + 0.726*avg_seq[key]