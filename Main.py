# -*- coding:utf-8 -*-

import time
import random
from concurrent.futures._base import _AcquireFutures
from sklearn.linear_model import LogisticRegression

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
        self.pairs = []
        self.train_set = []  # 90%
        self.validation_set = []  # 10%
        self.train_inters = {}  # molregno1 + molregno2: interaction level, all interactions between drugs in train set
        self.val_maccs_pair_mol = {}  # drug molregno: similar drug's molregno
        self.val_topo_pair_mol = {}  # drug molregno: similar drug's molregno
        self.val_ecfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.val_fcfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.train_maccs_pair_mol = {}  # drug molregno: similar drug's molregno
        self.train_ecfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.train_fcfp4_pair_mol = {}  # drug molregno: similar drug's molregno
        self.train_topo_pair_mol = {}  # drug molregno: similar drug's molregno
        self.val_maccs_pair_id = {}  # drug molregno: similar drug's id
        self.val_topo_pair_id = {}  # drug molregno: similar drug's id
        self.val_ecfp4_pair_id = {}  # drug molregno: similar drug's id
        self.val_fcfp4_pair_id = {}  # drug molregno: similar drug's id
        self.train_maccs_pair_id = {}  # drug molregno: similar drug's id
        self.train_ecfp4_pair_id = {}  # drug molregno: similar drug's id
        self.train_fcfp4_pair_id = {}  # drug molregno: similar drug's id
        self.train_topo_pair_id = {}  # drug molregno: similar drug's id
        self.maccs = {}
        self.ecfp4 = {}
        self.fcfp4 = {}
        self.topo = {}

    def input_sims(self, maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict):
        self.maccs = maccs_dict
        self.ecfp4 = ecfp4_dict
        self.fcfp4 = fcfp4_dict
        self.topo = topo_dict

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

    def logistic_regression(self):
        # find interactions in train set
        self.create_interactions_train_set()
        # for pairs in validation set, find most similar pairs
        self.create_pairs_for_train_set()
        self.create_pairs_for_validation_set()
        # create training array
        for d1 in self.train_set:
            for d2 in self.train_set:
                feature = self.link_sim(d1, d2)
                print(feature)
        # Logistic Regression
        # lr = LogisticRegression(multi_class='multinomial')
        # lr.fit()
        # lr.predict()
        return 0

    def create_pairs_for_validation_set(self):
        for val in self.validation_set:
            maxmaccs = 0
            maxecfp = 0
            maxfcfp = 0
            maxtopo = 0
            for train in self.train_set:
                maccs = self.sim_by_mol(val.molregno, train.molregno, 0)
                ecfp = self.sim_by_mol(val.molregno, train.molregno, 1)
                fcfp = self.sim_by_mol(val.molregno, train.molregno, 2)
                topo = self.sim_by_mol(val.molregno, train.molregno, 3)
                # print('maccs,ecfp,fcfp,topo: ', maccs, ecfp, fcfp, topo)
                if maccs >= maxmaccs:
                    maxmaccs = maccs
                    self.val_maccs_pair_mol[val.molregno] = train.molregno
                    self.val_maccs_pair_id[val.molregno] = train.id
                if ecfp >= maxecfp:
                    maxecfp = ecfp
                    self.val_ecfp4_pair_mol[val.molregno] = train.molregno
                    self.val_ecfp4_pair_id[val.molregno] = train.id
                if fcfp >= maxfcfp:
                    maxfcfp = fcfp
                    self.val_fcfp4_pair_mol[val.molregno] = train.molregno
                    self.val_fcfp4_pair_id[val.molregno] = train.id
                if topo >= maxtopo:
                    maxtopo = topo
                    self.val_topo_pair_mol[val.molregno] = train.molregno
                    self.val_topo_pair_id[val.molregno] = train.id

    def create_pairs_for_train_set(self):  # for training process
        for train1 in self.train_set:
            maxmaccs = 0
            maxecfp = 0
            maxfcfp = 0
            maxtopo = 0
            for train2 in self.train_set:
                if train1 != train2:
                    maccs = self.sim_by_mol(train1.molregno, train2.molregno, 0)
                    ecfp = self.sim_by_mol(train1.molregno, train2.molregno, 1)
                    fcfp = self.sim_by_mol(train1.molregno, train2.molregno, 2)
                    topo = self.sim_by_mol(train1.molregno, train2.molregno, 3)
                    if maccs >= maxmaccs:
                        maxmaccs = maccs
                        self.train_maccs_pair_mol[train1.molregno] = train2.molregno
                        self.train_maccs_pair_id[train1.molregno] = train2.id
                    if ecfp >= maxecfp:
                        maxecfp = ecfp
                        self.train_ecfp4_pair_mol[train1.molregno] = train2.molregno
                        self.train_ecfp4_pair_id[train1.molregno] = train2.id
                    if fcfp >= maxfcfp:
                        maxfcfp = fcfp
                        self.train_fcfp4_pair_mol[train1.molregno] = train2.molregno
                        self.train_fcfp4_pair_id[train1.molregno] = train2.id
                    if topo >= maxtopo:
                        maxtopo = topo
                        self.train_topo_pair_mol[train1.molregno] = train2.molregno
                        self.train_topo_pair_id[train1.molregno] = train2.id

    def create_interactions_train_set(self):  # find all interactions between train set
        for d1 in self.train_set:
            for d2 in self.train_set:
                if d1 != d2:
                    key = d1.id + ' ' + d2.id
                    if key in self.interaction.keys():
                        self.train_inters[key] = self.interaction[key]

    def link_sim(self, d1, d2):  # create training array of drug d1 and d2
        # find interaction lvl between d1, d2
        inter = self.interaction_by_id(d1.id, d2.id)
        if 1:
            # calculate sim feature using (sim(s1,d1) + sim(s2,d2))/2 * interaction lvl(s1,s2)
            s1_mol = self.train_maccs_pair_mol[d1.molregno]
            s2_mol = self.train_maccs_pair_mol[d2.molregno]
            s1_id = self.train_maccs_pair_id[d1.molregno]
            s2_id = self.train_maccs_pair_id[d2.molregno]
            maccs1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=0)
            maccs2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=0)
            # if float(self.interaction_by_id(s1_id, s2_id))!=0:
            #     print(float(maccs1), float(maccs2), float(self.interaction_by_id(s1_id, s2_id)))
            feature1 = (float(maccs1) + float(maccs2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.train_ecfp4_pair_mol[d1.molregno]
            s2_mol = self.train_ecfp4_pair_mol[d2.molregno]
            s1_id = self.train_ecfp4_pair_id[d1.molregno]
            s2_id = self.train_ecfp4_pair_id[d2.molregno]
            ecfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=1)
            ecfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=1)
            feature2 = (float(ecfp41) + float(ecfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.train_fcfp4_pair_mol[d1.molregno]
            s2_mol = self.train_fcfp4_pair_mol[d2.molregno]
            s1_id = self.train_fcfp4_pair_id[d1.molregno]
            s2_id = self.train_fcfp4_pair_id[d2.molregno]
            fcfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=2)
            fcfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=2)
            feature3 = (float(fcfp41) + float(fcfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.train_topo_pair_mol[d1.molregno]
            s2_mol = self.train_topo_pair_mol[d2.molregno]
            s1_id = self.train_topo_pair_id[d1.molregno]
            s2_id = self.train_topo_pair_id[d2.molregno]
            topo1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=3)
            topo2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=3)
            feature4 = (float(topo1) + float(topo2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            return [feature1, feature2, feature3, feature4, inter]
        else:
            return [0, 0, 0, 0, 0]

    def link_sim_val(self, d1, d2):  # create similar array of drug d1 and d2
        # find interaction lvl between d1, d2
        inter = self.interaction_by_id(d1.id, d2.id)
        if 1:
            # calculate sim feature using (sim(s1,d1) + sim(s2,d2))/2 * interaction lvl(s1,s2)
            s1_mol = self.val_maccs_pair_mol[d1.molregno]
            s2_mol = self.val_maccs_pair_mol[d2.molregno]
            s1_id = self.val_maccs_pair_id[d1.molregno]
            s2_id = self.val_maccs_pair_id[d2.molregno]
            maccs1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=0)
            maccs2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=0)
            # if float(self.interaction_by_id(s1_id, s2_id))!=0:
            #     print(float(maccs1), float(maccs2), float(self.interaction_by_id(s1_id, s2_id)))
            feature1 = (float(maccs1) + float(maccs2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.val_ecfp4_pair_mol[d1.molregno]
            s2_mol = self.val_ecfp4_pair_mol[d2.molregno]
            s1_id = self.val_ecfp4_pair_id[d1.molregno]
            s2_id = self.val_ecfp4_pair_id[d2.molregno]
            ecfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=1)
            ecfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=1)
            feature2 = (float(ecfp41) + float(ecfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.val_fcfp4_pair_mol[d1.molregno]
            s2_mol = self.val_fcfp4_pair_mol[d2.molregno]
            s1_id = self.val_fcfp4_pair_id[d1.molregno]
            s2_id = self.val_fcfp4_pair_id[d2.molregno]
            fcfp41 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=2)
            fcfp42 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=2)
            feature3 = (float(fcfp41) + float(fcfp42)) * float(self.interaction_by_id(s1_id, s2_id)) / 2

            s1_mol = self.val_topo_pair_mol[d1.molregno]
            s2_mol = self.val_topo_pair_mol[d2.molregno]
            s1_id = self.val_topo_pair_id[d1.molregno]
            s2_id = self.val_topo_pair_id[d2.molregno]
            topo1 = self.sim_by_mol(s1_mol, d1.molregno, sim_type=3)
            topo2 = self.sim_by_mol(s2_mol, d2.molregno, sim_type=3)
            feature4 = (float(topo1) + float(topo2)) * float(self.interaction_by_id(s1_id, s2_id)) / 2
            return [feature1, feature2, feature3, feature4, inter]
        else:
            return [0, 0, 0, 0, 0]

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
            return self.interaction[key1]
        elif key2 in self.interaction.keys():
            return self.interaction[key2]
        else:
            return 0


start = time.time()

maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict = Similarity.read_sims_to_dict()  # 1864590
wstmed_molregno, wstmed_id = WstMed.read_wstmed_to_dict()  # search western medicine via molregno
interaction_dict = Interaction.read_interactions_to_dict(wstmed_id)  # 128535 inters from totally 853272 of 4696 drugs

v = Validation(wstmed_id, maccs_dict, interaction_dict)
v.divide_data()
v.input_sims(maccs_dict, ecfp4_dict, fcfp4_dict, topo_dict)
# v.sim = maccs_dict
v.create_pairs_for_validation_set()
v.create_pairs_for_train_set()
v.create_interactions_train_set()
# v.logistic_regression()


train_array = []
num = 0
ar1 = []
ar2 = []
ar3 = []
ar4 = []
inters = []
for d1 in v.train_set:
    for d2 in v.train_set:
        if d1 != d2:
            # if v.interaction_by_id(d1.id,d2.id) != 0:
            #     print(d1.id, d2.id, v.interaction_by_id(d1.id,d2.id))
            f1, f2, f3, f4, inter = v.link_sim(d1, d2)
            # if f1!=0 and f2!=0 and f3!=0 and f4!=0:
            if 1:
                ar1.append(f1)
                ar2.append(f2)
                ar3.append(f3)
                ar4.append(f4)
                inters.append(inter)
                # print(f1, f2, f3, f4, inter)
            # if feature != 0:
            #     train_array.append(feature)
            #     print(feature)
                num += 1
#

train_array = [ar1, ar2, ar3, ar4, inters]
tr = [ar1, ar2, ar3, ar4]
a = [list(x) for x in zip(*tr)]

lr = LogisticRegression(solver='sag', multi_class='multinomial')
lr.fit(a, inters)
re = lr.predict(a)
same = 0
unsame = 0
num = 0
for i in range(0, inters.__len__()):

    if int(re[i]) != 0:
        num += 1
        if int(re[i]) == inters[i]:
            same += 1
        else:
            print(re[i], inters[i])
            unsame += 1
print(same, unsame, num)
print(same/num)

end = time.time()
print('time: ', end - start)

train_array = []
num = 0
ar1 = []
ar2 = []
ar3 = []
ar4 = []
inters = []
for d1 in v.validation_set:
    for d2 in v.validation_set:
        if d1 != d2:
            # if v.interaction_by_id(d1.id,d2.id) != 0:
            #     print(d1.id, d2.id, v.interaction_by_id(d1.id,d2.id))
            f1, f2, f3, f4, inter = v.link_sim_val(d1, d2)
            # if f1!=0 and f2!=0 and f3!=0 and f4!=0:
            if 1:
                ar1.append(f1)
                ar2.append(f2)
                ar3.append(f3)
                ar4.append(f4)
                inters.append(inter)
                # TODO: all inters are 0
                # print(f1, f2, f3, f4, inter)
            # if feature != 0:
            #     train_array.append(feature)
            #     print(feature)
                num += 1
#

train_array = [ar1, ar2, ar3, ar4, inters]
tr = [ar1, ar2, ar3, ar4]
a = [list(x) for x in zip(*tr)]

re = lr.predict(a)
same = 0
unsame = 0
num = 0
for i in range(0, inters.__len__()):
    num+=1
    if re[i] == inters[i]:
        same += 1
    else:
        unsame += 1
print(same, unsame, num)
print(same/num)