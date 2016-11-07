# -*- coding:utf-8 -*-

import time
import random
from concurrent.futures._base import _AcquireFutures
from sklearn import svm

class Similarity:
    def __init__(self, med_id1=0, med_id2=0, maccs=0, fcfp4=0, ecfp4=0, topo=0, weighted_sim=0):
        self.med_molregno1 = med_id1
        self.med_molregno2 = med_id2
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
        self.test_set = []
        self.validation_set = []

    def divide_data(self):
        self.test_set = []
        self.validation_set = []
        index = random.sample(range(0, 1366), 136)  # randomly select 1/10 data as test_set
        for i in range(0, self.wst_med.__len__()):
            if i not in index:
                self.validation_set.append(self.wst_med[i])
            else:
                self.test_set.append(self.wst_med[i])

    def __get_similarity(self, med1_id, med2_id):
        for item in self.sim:
            if item.drug1_id == med1_id and item.drug2_id == med2_id:
                return item.weighted_sim
        return -1

    def get_interaction_level(self, med_id2, inter_levels):
        for item in inter_levels:
            if med_id2 == item[1]:
                return float(item[2])
        return 0

    def generate_interaction_levels(self, med_id):
        interaction_levels = []  # [med_id, inter-med_id, inter_level]
        for item in self.interaction:
            if item.medicine_id1 == med_id:
                interaction_levels.append([med_id, item.medicine_id2, item.interaction_level])
            elif item.medicine_id2 == med_id:
                interaction_levels.append([med_id, item.medicine_id1, item.interaction_level])
        return interaction_levels

    def save_pair2(self, filename):
        self.generate_pairs()
        # write pair result to file
        file = open(filename, 'w')
        for item in self.pairs:
            line = str(item[0]) + ' ' + str(item[1]) + ' ' + str(item[2])
            file.write(line + '\n')
        file.close()

    def generate_pairs(self):
        # only find pair in test_set
        pairs = []
        for med in self.test_set:
            max_sim = 0
            # paired = 0
            for sim_item in self.sim:
                if sim_item.med_molregno1 == med.molregno:
                    if sim_item.weighted_sim > max_sim:
                        for vali in self.validation_set:
                            if sim_item.med_molregno2 == vali.molregno:
                                max_sim = sim_item.weighted_sim
                                paired = vali
            pairs.append([med, paired, max_sim])  # structure: [med_obj, paired_med_obj, similarity]

        # another approach, too slow
        # pairs = []
        # for med in self.test_set:
        #     max_sim = 0
        #     for vali in self.validation_set:
        #         for sim_item in self.sim:
        #             if sim_item.med_molregno1 == med.molregno and sim_item.med_molregno2 == vali.molregno:
        #                 if sim_item.weighted_sim > max_sim:
        #                     max_sim = sim_item.weighted_sim
        #                     paired = vali
        #     pairs.append([med, paired, max_sim])  # structure: [med_obj, paired_med_obj, similarity]

        self.pairs = pairs

    def find_pair2(self, med_molregno):  # for validation set
        for item in self.pairs:
            if item[0].molregno == med_molregno:
                return [item[1], float(item[2])]

    def save_pair(self, filename):
        # find the most similar medcs and save
        pairs = []
        for med in self.wst_med:
            max_sim = 0
            paired_molregno = 0
            for sim_item in self.sim:
                if sim_item.med_molregno1 == med.molregno:
                    if sim_item.weighted_sim > max_sim:
                        max_sim = sim_item.weighted_sim
                        paired_molregno = sim_item.med_molregno2
            pairs.append([med.molregno, paired_molregno, max_sim])

        # write pair result to file
        file = open(filename, 'w')
        for item in pairs:
            line = str(item[0]) + ' ' + str(item[1]) + ' ' + str(item[2])
            file.write(line + '\n')
        file.close()

    def compute_conflict_index(self):
        result = []
        same = 0
        total = 0
        zero_level = 0
        tt_level = 0
        zero_real_level = 0
        for test_med in self.test_set:
            pair_med, pair_sim = self.find_pair2(test_med.molregno)  # find most similar drug in validation set
            interaction_levels = self.generate_interaction_levels(pair_med.id)  # get pair_med's all interactions
            inter_level_validate = self.generate_interaction_levels(test_med.id)  # get test_med's all interactions
            print('pair:', test_med.id, pair_med.id, pair_sim)
            for known_med in self.validation_set:
                inter_level = self.get_interaction_level(known_med.id, interaction_levels)
                real_level = self.get_interaction_level(known_med.id, inter_level_validate)
                if inter_level == 0:
                    zero_level += 1
                tt_level += 1
                if real_level == 0:
                    zero_real_level += 1
                if inter_level != 0:
                    c_index = inter_level * pair_sim
                    if inter_level == real_level:
                        same += 1
                    if real_level != 0:
                        total += 1
                    print(pair_med.id, known_med.id, inter_level, c_index, real_level)
                    # store test_med, pair_med, knowm_med and conflict index
                    result.append([test_med, pair_med, known_med, c_index, real_level])
        print(same, total, float(same)/float(total))
        print('zero_lvl:', zero_level)
        print('zero_real_lvl:', zero_real_level)
        print('tt_num_lvl:', tt_level)
        return result

    def compute_validate_para(self, conflict_index):
        num_conflict = 0
        ratio = 0.6
        ratio_recall = 0
        return ratio_recall

    def validate(self):
        self.generate_pairs()
        conflict_index = self.compute_conflict_index()
        ratio = self.compute_validate_para(conflict_index)
        print(ratio)

start = time.time()

similarities = Similarity.read_similarities()
interactions = Interaction.read_interactions()
wst_med = WstMed.read_wstmed()
# pairs = Similarity.read_pairs()

# chn_med = ChnMed.read_chn_med()
v = Validation(wst_med, similarities, interactions)
v.divide_data()
v.generate_pairs()
conflict_index = v.compute_conflict_index()
# ratio = v.compute_validate_para(conflict_index)
# v.validate()
# v.find_and_save_pair('pairs.txt')
# result = v.validate()

################# Generate Similarities #################
# sims = SimOperation.sim_table()
# SimOperation.write_similarities(sims)
# num = 1
# for x in range(1, chn_med.__len__()):
#     if chn_med[x].chn_name != chn_med[x-1].chn_name:
#         num += 1

end = time.time()
print('time:', end - start, 's')
