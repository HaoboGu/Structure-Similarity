# -*- coding:utf-8 -*-

from rdkit import Chem
from rdkit.DataStructs import cDataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
import os
# import time
import random


def read_molfiles():
    path = 'medicine_structure'  # folder name of drugs
    folders = os.listdir(path)
    num = folders.__len__()
    data_set = []
    for i in range(0, num):
        filename = path + '/' + folders[i] + '/molfile.txt'
        data_set.append(Medicine(folders[i], Chem.MolFromMolFile(filename)))
    return data_set


class SimOperation:  # compute similarities, read/write sims, create sim_table
    p = {'maccs': 0.25, 'ecfp4': 0.25, 'fcfp4': 0.25, 'topo': 0.25, 'target': 0, 'mechanism': 0}

    @classmethod
    def get_similarity(cls, m):  # get weighted similarity of 2 compounds, with weights p
        target = 0
        mechanism = 0
        maccs = cls.similarity_maccs(m)
        ecfp4 = cls.similarity_ecfp4(m)
        fcfp4 = cls.similarity_fcfp4(m)
        topo = cls.similarity_topo(m)
        return cls.p['maccs']*maccs + cls.p['ecfp4']*ecfp4 + cls.p['fcfp4']*fcfp4 + cls.p['topo']*topo + \
            cls.p['target']*target + cls.p['mechanism']*mechanism

    @staticmethod
    def similarity_maccs(m):  # compute similarity using MACCS fingerprint
        fp1 = MACCSkeys.GenMACCSKeys(m[0])
        fp2 = MACCSkeys.GenMACCSKeys(m[1])
        return cDataStructs.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_ecfp4(m):  # ECFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2)  # 2 is the radius of the atom environments considered
        fp2 = AllChem.GetMorganFingerprint(m[1], 2)
        return cDataStructs.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_fcfp4(m):  # FCFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2, useFeatures=True)
        fp2 = AllChem.GetMorganFingerprint(m[1], 2, useFeatures=True)
        return cDataStructs.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_topo(m):  # Topological Fingerprints
        fp1 = FingerprintMols.FingerprintMol(m[0])
        fp2 = FingerprintMols.FingerprintMol(m[1])
        return cDataStructs.TanimotoSimilarity(fp1, fp2)

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

    @classmethod
    def sim_table(cls):  # create similarity table
        data = read_molfiles()
        result = []
        # compute similarity
        for d1 in data:
            for d2 in data:
                if d1 != d2:
                    result.append([d1.ID, d2.ID,
                                   cls.similarity_maccs([d1.mol, d2.mol]),
                                   cls.similarity_ecfp4([d1.mol, d2.mol]),
                                   cls.similarity_fcfp4([d1.mol, d2.mol]),
                                   cls.similarity_topo([d1.mol, d2.mol]),
                                   cls.get_similarity([d1.mol, d2.mol])])
                # print d2.ID
        return result

    @staticmethod
    def write_similarities(sim):  # sim: a list of Similarity obj
        # Write the result to .txt file
        f = open('result.txt', 'w')
        for item in sim:
            # sim_table = item.get_simtable()
            tmp = ''.join([str(i) + ' ' for i in item])
            f.write(tmp + '\n')
        f.close()


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

    # def printsim(self):
    #     print self.drug1_id, self.drug2_id, self.maccs, self.ecfp4, self.fcfp4, self.topo, self.weighted_sim


class Medicine:  # Medicine used to read mol files

    def __init__(self, med_id, mol):
        self.ID = med_id
        self.mol = mol


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
        print 'sim:', sim
        interaction_level = self.__get_interaction_level(med1_id, med2_id)
        print 'inter_level', interaction_level
        conflict_index = sim * interaction_level
        # TODO: The meaning of doing this? Without Chinese medicine, should we
        # find the most similar medicine then predict whether it conflict with
        # a known medicine? 
        return conflict_index
    
    def compute_CI(self, drug1, drug2):  # CI: Conflict Index
        return


class ChnMed:  # Chinese medicine class

    def __init__(self, lst):
        self.id = lst[0]
        self.chn_name = lst[1]
        self.chn_word_id = lst[2]
        self.component = lst[3]
        self.description = lst[4]
        self.chn_description = lst[5]


class WstMed:  # Western medicine class

    def __init__(self, lst):
        self.id = lst[0]
        self.name = lst[1]
        self.component = lst[2]
        self.molregno = lst[3]


class Interaction:  # interaction between western medicines

    def __init__(self, lst):
        self.id = lst[0]
        self.medicine_id1 = lst[1]
        self.medicine_name1 = lst[2]
        self.medicine_id2 = lst[3]
        self.medicine_name2 = lst[4]
        self.interaction_level = lst[5]


# read chinese medicine data
def read_chnmed():
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

# read interaction data
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

similarities = SimOperation.read_similarities()
interactions = read_interactions()
medicine = read_molfiles()
v = Validation(medicine, similarities, interactions)
v.divide_data()
c_index = v.conflict_index(v.test_set[0].ID, v.validation_set[0].ID)
print c_index


################# Generate Similarities #################
# sims = SimOperation.sim_table()
# SimOperation.write_similarities(sims)
