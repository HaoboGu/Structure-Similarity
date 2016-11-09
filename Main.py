# -*- coding:utf-8 -*-

import time
import random
from concurrent.futures._base import _AcquireFutures
from sklearn import svm

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
        self.sim = similarities
        self.interaction = interaction
        self.pairs = []
        self.train_set = []
        self.validation_set = []
        self.train_inters = {}
        self.val_train_sim = {}
        self.train_train_sim = {}

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

    def __get_similarity(self, med1_id, med2_id):
        for item in self.sim:
            if item.drug1_id == med1_id and item.drug2_id == med2_id:
                return item.weighted_sim
        return -1

    def logistic_regression(self):
        # find interactions in train set
        for d1 in self.train_set:
            for d2 in self.train_set:
                if d1 != d2:
                    key = d1.id + ' ' + d2.id
                    if key in self.interaction.keys():
                        self.train_inters[key] = self.interaction[key]
        # for pairs in validation set, find most similar pairs
        for d1 in self.validation_set:
            for d2 in self.validation_set:
                if d1 != d2:
                    #TODO
                    continue  #

        return 0

    def create_pair_for_validation_set(self):
        for val in self.validation_set:
            maxsim = 0
            for train in self.train_set:
                key = val.molregno + ' ' + train.molregno
                key2 = train.molregno + ' ' + val.molregno
                if key in self.sim.keys():
                    if self.sim[key] > maxsim:
                        maxsim = self.sim[key]
                        self.val_train_sim[val.molregno] = train.molregno
                elif key2 in self.sim.keys():
                    if self.sim[key2] > maxsim:
                        maxsim = self.sim[key2]
                        self.val_train_sim[val.molregno] = train.molregno
                else:
                    print('error: cannot find similarity between', val.molregno, train.molregno)

    def create_pair_for_train_set(self):  # for training process
        for t1 in self.train_set:
            maxsim = 0
            for t2 in self.train_set:
                if t1 != t2:
                    key = t1.molregno + ' ' + t2.molregno
                    key2 = t2.molregno + ' ' + t1.molregno
                    if key in self.sim.keys():
                        if self.sim[key] > maxsim:
                            maxsim = self.sim[key]
                            self.train_train_sim[t1.molregno] = t2.molregno
                    elif key2 in self.sim.keys():
                        if self.sim[key2] > maxsim:
                            maxsim = self.sim[key2]
                            self.train_train_sim[t1.molregno] = t2.molregno
                    else:
                        print('error: cannot find similarity between', t1.molregno, t2.molregno)




start = time.time()

maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict = Similarity.read_sims_to_dict()  # 1864590
wstmed_molregno, wstmed_id = WstMed.read_wstmed_to_dict()  # search western medicine via molregno
interaction_dict = Interaction.read_interactions_to_dict(wstmed_id)  # 128535 inters from totally 853272 of 4696 drugs

v = Validation(wstmed_id, maccs_dict, interaction_dict)
v.divide_data()
v.create_pair_for_validation_set()
v.create_pair_for_train_set()
# v.logistic_regression()
end = time.time()
print('time: ', end - start)

for i in v.train_train_sim:
    key1 = v.train_train_sim[i] + ' ' + i
    key2 = i + ' ' + v.train_train_sim[i]
    if key1 in v.sim.keys():
        print(v.sim[key1])
    elif key2 in v.sim.keys():
        print(v.sim[key2])
    else:
        print('eeeeeeeeeeeeeeeeeeeee')