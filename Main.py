# -*- coding:utf-8 -*-

import time
import random


class Similarity:
    def __init__(self, drug1_id=0, drug2_id=0, maccs=0, fcfp4=0, ecfp4=0, topo=0, weighted_sim=0):
        self.drug1_id = drug1_id
        self.drug2_id = drug2_id
        self.maccs = maccs
        self.ecfp4 = ecfp4
        self.fcfp4 = fcfp4
        self.topo = topo
        self.weighted_sim = weighted_sim

    def get_simtable(self):
        return [self.drug1_id, self.drug2_id, self.maccs, self.ecfp4, self.fcfp4, self.topo, self.weighted_sim]

    def from_simtable(self, table):
        self.drug1_id = table[0]
        self.drug2_id = table[1]
        self.maccs = table[2]
        self.ecfp4 = table[3]
        self.fcfp4 = table[4]
        self.topo = table[5]
        self.weighted_sim = table[6]

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

    def __get_interaction_level(self, med1_id, med2_id):
        for item in self.interaction:
            if item.medicine_id1 == med1_id and item.medicine_id2 == med2_id:
                return item.interaction_level
        return -1

    def conflict_index(self, med1_id, med2_id): 
        # Just compute the conflict index between 2 medicines
        conflict = []
        sim = self.__get_similarity(med1_id, med2_id)
        # print 'sim:', sim
        interaction_level = self.__get_interaction_level(med1_id, med2_id)
        # print 'inter_level', interaction_level
        conflict_index = sim * interaction_level
        # Meaning of doing this? Without Chinese medicine, should we
        # find the most similar medicine then predict whether it conflict with
        # a known medicine? 
        return conflict_index


start = time.time()

similarities = Similarity.read_similarities()
interactions = Interaction.read_interactions()
wst_med = WstMed.read_wstmed()
chn_med = ChnMed.read_chn_med()
# v = Validation(medicine, similarities, interactions)
# v.divide_data()
# c_index = v.conflict_index(v.test_set[0].ID, v.validation_set[0].ID)
# # print c_index


################# Generate Similarities #################
# sims = SimOperation.sim_table()
# SimOperation.write_similarities(sims)

end = time.time()
print(end - start)
