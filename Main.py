# -*- coding:utf-8 -*-

import random
import time
import xlsxwriter
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
import csv

class Similarity:
    def __init__(self, med_molregno1=0, med_molregno2=0, maccs=0, fcfp4=0, ecfp4=0, topo=0, weighted_sim=0):
        self.med_molregno1 = med_molregno1
        self.med_molregno2 = med_molregno2
        self.maccs = maccs
        self.ecfp4 = ecfp4
        self.fcfp4 = fcfp4
        self.topo = topo
        self.weighted_sim = weighted_sim

    def get_simtable(self):
        return [self.med_molregno1, self.med_molregno2, self.maccs, self.ecfp4, self.fcfp4, self.topo, self.weighted_sim]

    def from_simtable(self, table):
        self.med_molregno1 = table[0]
        self.med_molregno2 = table[1]
        self.maccs = float(table[2])
        self.ecfp4 = float(table[3])
        self.fcfp4 = float(table[4])
        self.topo = float(table[5])
        self.weighted_sim = float(table[6])

    @staticmethod
    def read_similarities():
        similarities = []
        sim_file = open('result.txt')
        while 1:
            s = Similarity()
            line = sim_file.readline()
            if not line:
                break
            s.from_simtable(line.split())
            # s.print()
            similarities.append(s)
        sim_file.close()
        return similarities

    @staticmethod
    def read_sims_to_dict():
        sim_file = open('result.txt')
        maccs_dict = {}
        ecfp4_dict = {}
        fcfp4_dict = {}
        topo_dict = {}
        while 1:
            line = sim_file.readline()
            if not line:
                break
            table = line.split()
            key = table[0] + ' ' + table[1]
            maccs_dict[key] = float(table[2])
            ecfp4_dict[key] = float(table[3])
            fcfp4_dict[key] = float(table[4])
            topo_dict[key] = float(table[5])
        sim_file.close()
        return [maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict]

    @staticmethod
    def read_pairs():
        pairs = []
        pairs_file = open('pairs.txt')
        while 1:
            line = pairs_file.readline()
            if not line:
                break
            item = line.split()  # item[0]: molregno1; item[1]: most similar mec's molregno; item[2]: similarity
            pairs.append(item)
        return pairs


class ChnMed:  # Chinese medicine class

    def __init__(self, lst):
        self.id = lst[0]
        self.chn_name = lst[1]
        self.chn_word_id = lst[2]
        self.component = lst[3]
        self.description = lst[4]
        self.chn_description = lst[5]

    # read chinese medicine data
    @staticmethod
    def read_chn_med():
        chn_med_file = open('CMedc.txt')
        chn_med = []
        while 1:
            line = chn_med_file.readline()
            if not line:
                break
            row = line.split()
            med = ChnMed(row)
            chn_med.append(med)
        chn_med_file.close()
        return chn_med

    def chn_str(self):
        return str(self.id) + ' ' + str(self.chn_name) + ' ' + str(self.chn_word_id) + ' ' +\
               str(self.component) + ' ' + str(self.description) + ' ' + str(self.chn_description)

    @staticmethod
    def write_chn_med(chn_med):
        file = open('CMedc1.txt', 'w')
        for item in chn_med:
            line = item.chn_str()
            file.write(line + '\n')
        file.close()


class WstMed:  # Western medicine class

    def __init__(self, lst):
        self.id = lst[0]  # drugs.com id
        self.name = lst[1]
        self.component = lst[2]
        self.molregno = lst[3]  # CHEMBL id
        self.smiles = lst[4]  # store medicine's SMILES notation, rather than mol object

    def wst_str(self):
        return str(self.id) + ' ' + str(self.name) + ' ' + str(self.component) + ' ' +\
               str(self.molregno) + ' ' + str(self.smiles)

    @staticmethod
    # read western medicine data
    def read_wstmed_to_dict():
        wst_med_file = open('WMedc.txt')
        wstmed_molregno_dict = {}
        wstmed_id_dict = {}
        while 1:
            line = wst_med_file.readline()
            if not line:
                break
            row = line.split()
            med = WstMed(row)
            wstmed_molregno_dict[med.molregno] = med
            wstmed_id_dict[med.id] = med
        wst_med_file.close()
        return [wstmed_molregno_dict, wstmed_id_dict]

    @staticmethod
    # read western medicine data
    def read_wstmed():
        wst_med_file = open('WMedc.txt')
        wst_med = []
        while 1:
            line = wst_med_file.readline()
            if not line:
                break
            row = line.split()
            med = WstMed(row)
            wst_med.append(med)
        wst_med_file.close()
        return wst_med


class Interaction:  # interaction between western medicines

    def __init__(self, lst):

        self.id = lst[0]
        self.medicine_id1 = lst[1]
        self.medicine_name1 = lst[2]
        self.medicine_id2 = lst[3]
        self.medicine_name2 = lst[4]
        self.interaction_level = lst[5]

    # read interaction data
    @staticmethod
    def read_interactions():
        interaction_file = open('interactions.txt')
        interactions = []
        while 1:
            line = interaction_file.readline()
            if not line:
                break
            row = line.split()
            inter = Interaction(row)
            interactions.append(inter)
        interaction_file.close()
        return interactions

    def read_interactions_to_dict(wstmed_id):
        interaction_file = open('interactions.txt')
        interactions_dict = {}
        while 1:
            line = interaction_file.readline()
            if not line:
                break
            row = line.split()
            if row[1] in wstmed_id.keys() and row[3] in wstmed_id.keys():  # consider interactions between 1366 drugs
                key = row[1] + ' ' + row[3]
                interactions_dict[key] = row[5]  # Key = drugs.com id strings
        interaction_file.close()
        return interactions_dict

    @staticmethod
    def write_interactions(interactions):
        interaction_file = open('interactions.txt', 'w')
        # interactions = []
        for item in interactions:
            line = item.interaction_str()
            interaction_file.write(line + '\n')
        interaction_file.close()

    def interaction_str(self):
        return self.id + ' ' + self.medicine_id1 + ' ' + self.medicine_name1 + ' ' + self.medicine_id2 + ' ' +\
               self.medicine_name2 + ' ' + self.interaction_level




class Validation:

    def __init__(self, wst_med, similarities, interaction):
        self.wst_med = wst_med
        self.sim = similarities  # key: molregno
        self.interaction = interaction  # key: drugs.com id
        self.train_set = []  # 90%
        self.validation_set = []  # 10%
        self.train_inters = {}  # molregno1 + molregno2: interaction level, all interactions between drugs in train set
        self.maccs_pair_mol = {}  # drug molregno: similar drug's molregno
        self.ecfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.fcfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.topo_pair_mol = {}  # drug molregno: similar drug's molregno
        self.maccs_pair_id = {}  # drug molregno: similar drug's id
        self.topo_pair_id = {}  # drug molregno: similar drug's id
        self.ecfp4_pair_id = {}  # drug molregno: similar drug's id
        self.fcfp4_pair_id = {}  # drug molregno: similar drug's id
        self.maccs = {}
        self.ecfp4 = {}
        self.fcfp4 = {}
        self.topo = {}
        self.mol_by_id = {}  # find molregno by drugs' id
        self.id_by_mol = {}  # find id by molregno
        self.index_array = np.zeros(1366)
        self.train_interactions = {}
        self.validation_interactions = {}
        self.inter_matrix = np.zeros(shape=(1366,1366))

    def input_sims(self, maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict):
        self.maccs = maccs_dict
        self.ecfp4 = ecfp4_dict
        self.fcfp4 = fcfp4_dict
        self.topo = topo_dict

    def sim_by_mol(self, mol1, mol2, sim_type=0):  # sim_type: 0-maccs, 1-ecfp4, 2-fcfp4, 3-topo
        key = mol1 + ' ' + mol2
        key2 = mol2 + ' ' + mol1
        if sim_type == 0:
            if key in self.maccs.keys():
                return self.maccs[key]
            elif key2 in self.maccs.keys():
                return self.maccs[key2]
            else:
                print("maccs_sim_by_mol error: no key ", key)
        elif sim_type == 1:
            if key in self.ecfp4.keys():
                return self.ecfp4[key]
            elif key2 in self.ecfp4.keys():
                return self.ecfp4[key2]
            else:
                print("ecfp4_sim_by_mol error: no key ", key)
        elif sim_type == 2:
            if key in self.fcfp4.keys():
                return self.fcfp4[key]
            elif key2 in self.fcfp4.keys():
                return self.fcfp4[key2]
            else:
                print("fcfp4_sim_by_mol error: no key ", key)
        elif sim_type == 3:
            if key in self.topo.keys():
                return self.topo[key]
            elif key2 in self.topo.keys():
                return self.topo[key2]
            else:
                print("topo_sim_by_mol error: no key ", key)
        else:
            print("similarity type error!!!!!")

    def interaction_by_id(self, id1, id2):
        key1 = id1 + ' ' + id2
        key2 = id2 + ' ' + id1
        if key1 in self.interaction.keys():
            # return int(self.interaction[key1])
            return 1
        elif key2 in self.interaction.keys():
            return 1
            # return int(self.interaction[key2])
        else:
            return 0

    def divide_data(self):
        self.train_set = []
        self.validation_set = []
        index = random.sample(range(0, 1366), 136)  # randomly select 1/10 data as test_set
        flag = 0
        for i in self.wst_med:
            if flag in index:
                self.validation_set.append(self.wst_med[i])
            else:
                self.train_set.append(self.wst_med[i])
            flag += 1

    def create_pairs_for_data_set(self):  # for training process
        for train1_index in self.wst_med:
            maxmaccs = 0
            maxecfp = 0
            maxfcfp = 0
            maxtopo = 0
            train1 = self.wst_med[train1_index]
            for train2 in self.train_set:
                if train1 != train2:
                    maccs = self.sim_by_mol(train1.molregno, train2.molregno, 0)
                    ecfp = self.sim_by_mol(train1.molregno, train2.molregno, 1)
                    fcfp = self.sim_by_mol(train1.molregno, train2.molregno, 2)
                    topo = self.sim_by_mol(train1.molregno, train2.molregno, 3)
                    if maccs >= maxmaccs:
                        maxmaccs = maccs
                        self.maccs_pair_mol[train1.molregno] = train2.molregno
                        self.maccs_pair_id[train1.molregno] = train2.id
                    if ecfp >= maxecfp:
                        maxecfp = ecfp
                        self.ecfp4_pair_mol[train1.molregno] = train2.molregno
                        self.ecfp4_pair_id[train1.molregno] = train2.id
                    if fcfp >= maxfcfp:
                        maxfcfp = fcfp
                        self.fcfp4_pair_mol[train1.molregno] = train2.molregno
                        self.fcfp4_pair_id[train1.molregno] = train2.id
                    if topo >= maxtopo:
                        maxtopo = topo
                        self.topo_pair_mol[train1.molregno] = train2.molregno
                        self.topo_pair_id[train1.molregno] = train2.id

    def create_interactions_train_set(self):  # find all interactions between train set
        for d1 in self.train_set:
            for d2 in self.train_set:
                if d1 != d2:
                    key = d1.id + ' ' + d2.id
                    if key in self.interaction.keys():
                        self.train_inters[key] = self.interaction[key]

    def create_id_mol_dict(self):
        for key in self.wst_med:
            self.mol_by_id[self.wst_med[key].id] = self.wst_med[key].molregno

    def create_mol_id_dict(self):
        for key in self.wst_med:
            self.id_by_mol[self.wst_med[key].molregno] = self.wst_med[key].id

    def link_sim(self, d1, d2):  # create training array of drug d1 and d2
        # find interaction lvl between d1, d2
        inter = self.interaction_by_id(d1.id, d2.id)
        if 1:
            # calculate sim feature using (sim(s1,d1) + sim(s2,d2))/2 * interaction lvl(s1,s2)
            s1_mol = self.maccs_pair_mol[d1.molregno]
            s2_mol = self.maccs_pair_mol[d2.molregno]
            s1_id = self.maccs_pair_id[d1.molregno]
            s2_id = self.maccs_pair_id[d2.molregno]
            maccs1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=0)
            maccs2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=0)
            feature1 = (float(maccs1) + float(maccs2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            # feature1 = (float(maccs1) + float(maccs2)) / 2

            s1_mol = self.ecfp4_pair_mol[d1.molregno]
            s2_mol = self.ecfp4_pair_mol[d2.molregno]
            s1_id = self.ecfp4_pair_id[d1.molregno]
            s2_id = self.ecfp4_pair_id[d2.molregno]
            ecfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=1)
            ecfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=1)
            feature2 = (float(ecfp41) + float(ecfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            # feature2 = (float(ecfp41) + float(ecfp42)) / 2

            s1_mol = self.fcfp4_pair_mol[d1.molregno]
            s2_mol = self.fcfp4_pair_mol[d2.molregno]
            s1_id = self.fcfp4_pair_id[d1.molregno]
            s2_id = self.fcfp4_pair_id[d2.molregno]
            fcfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=2)
            fcfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=2)
            feature3 = (float(fcfp41) + float(fcfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            # feature3 = (float(fcfp41) + float(fcfp42)) / 2

            s1_mol = self.topo_pair_mol[d1.molregno]
            s2_mol = self.topo_pair_mol[d2.molregno]
            s1_id = self.topo_pair_id[d1.molregno]
            s2_id = self.topo_pair_id[d2.molregno]
            topo1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=3)
            topo2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=3)
            feature4 = (float(topo1) + float(topo2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            # feature4 = (float(topo1) + float(topo2)) / 2
            return [feature1, feature2, feature3, feature4, inter]
        else:
            return [0, 0, 0, 0, 0]

    def create_train_array(self, portion):
        ar1 = []
        ar2 = []
        ar3 = []
        ar4 = []
        inters = []
        for d1 in v.train_set:
            for d2 in v.train_set:
                if d1 != d2:
                    f1, f2, f3, f4, inter = v.link_sim(d1, d2)
                    # if f1!=0 and f2!=0 and f3!=0 and f4!=0:
                    if inter != 0:
                        ar1.append(f1)
                        ar2.append(f2)
                        ar3.append(f3)
                        ar4.append(f4)
                        inters.append(inter)
                    else:
                        index = random.sample(range(0, 100), 1)  # randomly select 1/10 data as test_set
                        if index[0] > portion:  # 25% zeros in training set
                            ar1.append(f1)
                            ar2.append(f2)
                            ar3.append(f3)
                            ar4.append(f4)
                            inters.append(inter)

        tr = [ar1, ar2, ar3, ar4]
        tr = [list(x) for x in zip(*tr)]  # transpose
        return [tr, inters]

    def create_val_array(self):
        ar1 = []
        ar2 = []
        ar3 = []
        ar4 = []
        inters = []
        for d1 in v.validation_set:
            # for d2 in v.validation_set:
            for d2 in v.train_set:
                if d1 != d2:
                    f1, f2, f3, f4, inter = v.link_sim(d1, d2)
                    if 1:
                        ar1.append(f1)
                        ar2.append(f2)
                        ar3.append(f3)
                        ar4.append(f4)
                        inters.append(inter)
        val = [ar1, ar2, ar3, ar4]
        val = [list(x) for x in zip(*val)]  # transpose
        return [val, inters]

    def logistic_regression(self, portion):
        # find interactions in train set
        # self.create_interactions_train_set()

        # for pairs in validation set, find most similar pairs
        # self.create_pairs_for_train_set()
        # self.create_pairs_for_validation_set()
        self.create_pairs_for_data_set()

        # create training array
        tr, inters = self.create_train_array(portion)
        self.tr = tr

        # train logistic regression model
        lr = LogisticRegression(solver='sag')
        lr.fit(tr, inters)
        # svm = LinearSVC()
        # print('start fitting')
        # svm.fit(tr, inters)
        # print('fit completed')

        # create validation array
        val, inters = self.create_val_array()
        self.val = val
        self.result = lr.predict(val)
        prob_re = lr.predict_proba(val)
        # self.result = svm.predict(val)
        # prob_re= prob_re.transpose()
        # print(prob_re.__len__(), inters.__len__())
        # fpr_grd_lm, tpr_grd_lm, _ = roc_curve(inters, prob_re[0])
        # plt.plot(fpr_grd_lm, tpr_grd_lm, label='GBT + LR')
        # validation
        same = 0
        unsame = 0
        num = 0
        for i in range(0, inters.__len__()):
            num += 1
            if int(self.result[i]) == inters[i]:
                same += 1
            else:
                unsame += 1

        TP = 0  # predict 1, actual 1
        FP = 0  # predict 1, actual 0
        TN = 0  # predict 0, actual 0
        FN = 0  # predict 0, actual 1
        for i in range(0, inters.__len__()):
            if int(self.result[i]) != 0 or inters[i] != 0:
                # print(self.result[i], inters[i])
                if int(self.result[i]) == int(inters[i]):
                    TP += 1
                elif int(self.result[i]) != 0 and inters[i] == 0:
                    FP += 1
                elif inters[i] != 0 and int(self.result[i])==0:
                    FN += 1
            elif int(self.result[i]) == 0 and inters[i] == 0:
                TN += 1
        print('TP:', TP)
        print('FP:', FP)
        print('TN:', TN)
        print('FN:', FN)
        precision = TP/(TP+FP)
        recall = TP/(TP+FN)

        print('precision:', precision)
        print('recall:', recall)
        print('f-score: ', 2*precision*recall/(precision + recall))
        print(same, unsame, num)
        print(same / num)
        return 0

    def find_most_similar_link(self, d1, d2):
        max_link_maccs = 0
        max_link_ecfp4 = 0
        max_link_fcfp4 = 0
        max_link_topo = 0
        summaccs = 0
        sumecfp = 0
        sumfcfp = 0
        sumtopo = 0
        maccs = []
        ecfp = []
        fcfp = []
        topo = []
        i = 0
        for link_key in self.train_inters:
            id1, id2 = link_key.split()
            if id1 != d1.id and id1 != d2.id:
                if id2 != d1.id and id2 != d2.id:
                    i = i + 1
                    link_maccs = (self.sim_by_mol(self.mol_by_id[id1], d1.molregno, 0) +
                                  self.sim_by_mol(self.mol_by_id[id2], d2.molregno, 0)) / 2.0
                    link_ecfp = (self.sim_by_mol(self.mol_by_id[id1], d1.molregno, 1) +
                                 self.sim_by_mol(self.mol_by_id[id2], d2.molregno, 1)) / 2.0
                    link_fcfp = (self.sim_by_mol(self.mol_by_id[id1], d1.molregno, 2) +
                                 self.sim_by_mol(self.mol_by_id[id2], d2.molregno, 2)) / 2.0
                    link_topo = (self.sim_by_mol(self.mol_by_id[id1], d1.molregno, 3) +
                                 self.sim_by_mol(self.mol_by_id[id2], d2.molregno, 3)) / 2.0
                    maccs.append(link_maccs)
                    ecfp.append(link_ecfp)
                    fcfp.append(link_fcfp)
                    topo.append(link_topo)
                    # if link_maccs >= max_link_maccs:
                    #     max_link_maccs = link_maccs
                    #     # result_id1 = id1
                    #     # result_id2 = id2
                    # if link_ecfp >= max_link_ecfp4:
                    #     max_link_ecfp4 = link_ecfp
                    # if link_fcfp >= max_link_fcfp4:
                    #     max_link_fcfp4 = link_fcfp
                    # if link_topo >= max_link_topo:
                    #     max_link_topo = link_topo
        maccsar = np.array(maccs)
        ecfpar = np.array(ecfp)
        fcfpar = np.array(fcfp)
        topoar = np.array(topo)
        max_link_maccs = maccsar.max()
        max_link_ecfp4 = ecfpar.max()
        max_link_fcfp4 = fcfpar.max()
        max_link_topo = topoar.max()
        max_link_maccs = (max_link_maccs - maccsar.mean())/maccsar.std()
        max_link_ecfp4 = (max_link_ecfp4 - ecfpar.mean())/ecfpar.std()
        max_link_fcfp4 = (max_link_fcfp4 - fcfpar.mean())/fcfpar.std()
        max_link_topo = (max_link_topo - topoar.mean())/topoar.std()
        # print(result_id1, result_id2, max_link_sim)
        return [max_link_maccs, max_link_ecfp4, max_link_fcfp4, max_link_topo]

    def fail_attempt(self):
        train_set = v.train_set[0:50]
        i = 0
        inters = []
        tr = []
        numof1 = 0
        for d1 in train_set:
            for d2 in train_set:
                i += 1
                print(i, '-th process....')
                inter = v.interaction_by_id(d1.id, d2.id)
                if inter != 0:
                    numof1 += 1
                inters.append(inter)
                start = time.time()
                feature = v.find_most_similar_link(d1, d2)
                print(feature, inter)
                end = time.time()
                print('time: ', end - start)
                tr.append(feature)
        print('1 in 100 interactions:', numof1)

        # tr = [list(x) for x in zip(*tr)]  # transpose
        lr = LogisticRegression(solver='sag', max_iter=10000)
        lr.fit(tr, inters)

        i = 0
        test_set = v.validation_set[20:40]
        val = []
        numof1 = 0
        val_inters = []
        for d1 in test_set:
            for d2 in test_set:
                i += 1
                print(i, '-th process....')
                inter = v.interaction_by_id(d1.id, d2.id)
                if inter != 0:
                    numof1 += 1
                val_inters.append(inter)
                featureval = v.find_most_similar_link(d1, d2)
                # print('val:', featureval)
                val.append(featureval)
        # val = [list(x) for x in zip(*val)]  # transpose
        for i in val:
            if np.isnan(i[1]):
                val.remove(i)

        v.result = lr.predict(val)
        print('val_inters', val_inters)
        print('predicted result', v.result)

    def create_index_array(self):
        # create index_array
        self.index_array = np.zeros(1366)
        i = 0
        for key in v.wst_med:
            self.index_array[i] = key
            i += 1

    def divide_interactions(self):
        self.train_interactions = {}
        self.validation_interactions = {}
        num = self.interaction.__len__()//10
        index = random.sample(range(0, self.interaction.__len__()), num)
        flag = 0
        for key in self.interaction:
            if flag in index:
                self.validation_interactions[key] = float(self.interaction[key])
            else:
                self.train_interactions[key] = float(self.interaction[key])
            flag += 1

    def get_inter_matrix(self):  # train interactions matrix
        # create index_array
        # self.create_index_array()
        self.inter_matrix = np.zeros(shape=(1366, 1366))
        for key in self.train_interactions:
            id1, id2 = key.split()
            row = np.where(self.index_array == float(id1))[0][0]
            col = np.where(self.index_array == float(id2))[0][0]
            self.inter_matrix[row][col] = float(self.train_interactions[key])
        for i in range(0, 1366):
            self.inter_matrix[i][i] = 0

    def sim_matrix(self):
        # maccs, ecfp, fcfp, topo matrix
        self.maccs_matrix = np.zeros(shape=(1366, 1366))
        self.ecfp_matrix = np.zeros(shape=(1366, 1366))
        self.fcfp_matrix = np.zeros(shape=(1366, 1366))
        self.topo_matrix = np.zeros(shape=(1366, 1366))
        self.create_mol_id_dict()
        # self.create_index_array()
        for key in self.maccs:
            mol1, mol2 = key.split()
            id1 = self.id_by_mol[mol1]
            id2 = self.id_by_mol[mol2]
            row = np.where(self.index_array == float(id1))[0][0]
            col = np.where(self.index_array == float(id2))[0][0]
            self.maccs_matrix[row][col] = float(self.maccs[key])
            self.ecfp_matrix[row][col] = float(self.ecfp4[key])
            self.fcfp_matrix[row][col] = float(self.fcfp4[key])
            self.topo_matrix[row][col] = float(self.topo[key])
        for index in range(0, 1366):
            self.maccs_matrix[index][index] = 0
            self.ecfp_matrix[index][index] = 0
            self.fcfp_matrix[index][index] = 0
            self.topo_matrix[index][index] = 0

    def create_predict_matrix(self, inter_matrix, sim_matrix):
        # m12 = inter_matrix.dot(sim_matrix)  # * is element-wise multiply
        m12 = np.zeros(shape=(1366, 1366))
        st = sim_matrix.transpose()
        for row in range(0,1366):
            for col in range(0,1366):
                m12[row][col] = (inter_matrix[row]*st[col]).max()
        m12t = m12.transpose()
        pos_bigger = m12t > m12
        return m12 - np.multiply(m12, pos_bigger) + np.multiply(m12t, pos_bigger)

    def matrix_approach(self):
        for key in v.interaction:
            v.interaction[key] = 1
        v.create_index_array()
        v.divide_interactions()
        v.get_inter_matrix()
        v.sim_matrix()
        re = v.create_predict_matrix(v.inter_matrix, v.maccs_matrix)

        # transform re matrix to re_interactions dict
        re_list = []
        re_interactions = {}
        for row in range(0, 1366):
            for col in range(0, 1366):
                id1 = str(int(v.index_array[row]))
                id2 = str(int(v.index_array[col]))
                key = id1 + ' ' + id2
                re_list.append([id1,id2,re[row][col]])
                re_interactions[key] = re[row][col]
        re_list = np.array(re_list)
        re_list = re_list[re_list.transpose()[2].argsort(), :]
        re_list = re_list[::-1]

        # count TP TN FN FP and precision, recall
        TP = 0  # predict 1, actual 1
        TN = 0  # predict 0, actual 0
        FN = 0  # predict 0, actual 1
        FP = 0  # predict 1, actual 0
        for item in re_list:
            key = item[0] + ' ' + item[1]
            if key in v.validation_interactions.keys():
                # print(v.validation_interactions[key], item[2])
                if float(item[2]) > 0.8:
                    TP += 1
                else:
                    FN += 1
            elif key not in v.train_interactions.keys():
                if float(item[2]) > 0.8:
                    FP += 1
        # for key in v.validation_interactions:
        #     # id1, id2 = key.split()
        #     # row = np.where(v.index_array == float(id1))[0][0]
        #     # col = np.where(v.index_array == float(id2))[0][0]
        #     # print(row, col)
        #     # if key not in re_interactions.keys():
        #     #     id1, id2 = key.split()
        #     #     key = id2 + ' ' + id1
        #     # print(v.validation_interactions[key], re_interactions[key])
        #     if v.validation_interactions[key] != 0 and re_interactions[key] >0.8:
        #         TP += 1
        #     elif v.validation_interactions[key] != 0 and re_interactions[key] <0.8:
        #         FN += 1
        # for key in re_interactions:
        #     if re_interactions[key] != 0:
        #         if key not in v.validation_interactions.keys():
        #             if key not in v.train_interactions.keys():
        #             # if v.validation_interactions[key] == 0:
        #                 FP += 1
            #
        print('TP:', TP)
        print('FP:', FP)
        print('TN:', TN)
        print('FN:', FN)
        precision = TP/(TP+FP)
        recall = TP/(TP+FN)

        print('precision:', precision)
        print('recall:', recall)
        print('f-score: ', 2*precision*recall/(precision + recall))

    def hehe_approach(self):
        for key in self.interaction:
            self.interaction[key] = 1
        self.create_pairs_for_data_set()  # create the most similar drug's mol/id according four similarities
        allre = {}
        result = {}
        for item in self.validation_set:
            # item = self.validation_set[key]
            pair_id = self.maccs_pair_id[item.molregno]
            pair_mol = self.maccs_pair_mol[item.molregno]
            # print(pair_mol)

            for key in self.wst_med:
                target = self.wst_med[key]
                if target.id != item.id and pair_id != target.id:
                    key = target.id + ' ' + item.id
                    pair_target_key = pair_id + ' ' + target.id
                    pair_target_key2 = target.id + ' ' + pair_id
                    if pair_target_key in self.interaction.keys():
                        result[key] = self.interaction[pair_target_key] * self.sim_by_mol(pair_mol, target.molregno, 0)
                    elif pair_target_key2 in self.interaction.keys():  # try another key
                        result[key] = self.interaction[pair_target_key2] * self.sim_by_mol(pair_mol, target.molregno, 0)
                    else:
                        result[key] = 0
                else:
                    key = target.id + ' ' + item.id
                    result[key] = 0
            allre = result
            # print('get result')

        TP = 0  # predict 1, actual 1
        TN = 0  # predict 0, actual 0
        FN = 0  # predict 0, actual 1
        FP = 0  # predict 1, actual 0
        for item in self.validation_set:
            for key in self.wst_med:
                target = self.wst_med[key]
                if target.id != item.id:
                    key = target.id + ' ' + item.id
                    key2 = item.id + ' ' + target.id
                    if key in self.interaction.keys():
                        if allre[key] >0:#== 1:  # key in self.interaction.keys() -> interaction[key] must be 1
                            TP += 1
                        else:
                            FN += 1
                    elif key2 in self.interaction.keys():
                        if allre[key] >0:#== 1:
                            TP += 1
                        else:
                            FN += 1
                    elif allre[key] >0:#== 1:
                        FP += 1
                    elif allre[key] <=0:#== 0:
                        TN += 1
        print('TP:', TP)
        print('FP:', FP)
        print('TN:', TN)
        print('FN:', FN)
        precision = TP / (TP + FP)
        recall = TP / (TP + FN)

        print('precision:', precision)
        print('recall:', recall)
        print('f-score: ', 2 * precision * recall / (precision + recall))
        print((TP+TN)/(TP+TN+FP+FN))
        return allre






start = time.time()

maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict = Similarity.read_sims_to_dict()  # 1864590
wstmed_molregno, wstmed_id = WstMed.read_wstmed_to_dict()  # search western medicine via molregno
interaction_dict = Interaction.read_interactions_to_dict(wstmed_id)  # 128535 inters from totally 853272 of 4696 drugs

v = Validation(wstmed_id, maccs_dict, interaction_dict)
v.divide_data()
v.input_sims(maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict)

count = {}
for interkey in v.interaction:
    id1, id2 = interkey.split()
    if id1 not in count.keys():
        count[id1] = 1
    else:
        count[id1] = count[id1] + 1

    if id2 not in count.keys():
        count[id2] = 1
    else:
        count[id2] = count[id2] + 1
# a = v.hehe_approach()
# v.sim = maccs_dict
# v.create_pairs_for_validation_set()
# v.create_pairs_for_train_set()
# v.create_interactions_train_set()
portions = [40, 50, 60, 70, 80]
for item in portions:
    v.logistic_regression(item)
# v.logistic_regression(75)
# v.create_interactions_train_set()
# v.create_id_mol_dict()
# start = time.time()
# v.hehe_approach()
# re = 0

# train_interactions = 0
# validation_interactions = 0
#
# re = np.array(list(count.values()))
# re.tofile('a.csv', ',')

end = time.time()
print('time: ', end - start)
