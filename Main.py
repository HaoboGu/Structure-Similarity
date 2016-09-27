# -*- coding:utf-8 -*-

from rdkit import Chem
from rdkit.DataStructs import cDataStructs as ds
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
        data_set.append(Drug(folders[i], Chem.MolFromMolFile(filename)))
    return data_set


class SimOperation:
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
        return ds.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_ecfp4(m):  # ECFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2)  # 2 is the radius of the atom environments considered
        fp2 = AllChem.GetMorganFingerprint(m[1], 2)
        return ds.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_fcfp4(m):  # FCFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2, useFeatures=True)
        fp2 = AllChem.GetMorganFingerprint(m[1], 2, useFeatures=True)
        return ds.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def similarity_topo(m):  # Topological Fingerprints
        fp1 = FingerprintMols.FingerprintMol(m[0])
        fp2 = FingerprintMols.FingerprintMol(m[1])
        return ds.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def read_similarities():
        similarities = []
        file = open('result.txt')
        while 1:
            s = Similarity()
            line = file.readline()
            if not line:
                break
            s.from_simtable(line.split())
            # s.print()
            similarities.append(s)
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
                print(d2.ID)
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

    def print(self):
        print(self.drug1_id, self.drug2_id, self.maccs, self.ecfp4, self.fcfp4, self.topo, self.weighted_sim)


class Drug:

    def __init__(self, drug_id, mol):
        self.ID = drug_id
        self.mol = mol


class Validation:

    def __init__(self, drug, similarities):
        self.drug = drug
        self.sim = similarities
        self.test_set = []
        self.validation_set = []

    def divide_data(self):
        self.test_set = []
        self.validation_set = []
        index = random.sample(range(0, 1366), 136)  # randomly select 1/10 data as test_set
        for i in range(0, self.drug.__len__()):
            if i not in index:
                self.validation_set.append(self.data[i])
            else:
                self.test_set.append(self.data[i])

    def find_conflict(self, drug_id):  # Find similar drugs of a drug
        conflict = []
        for item in self.sim:
            if item.drug1_id == drug_id and item.weighted_sim > 0.5:  # temporary , a more concise criterion is needed
                item.print()
                conflict.append(item)

    def compute_CI(self, drug1, drug2):  # CI: Conflict Index
        return


# drugs = read_molfiles()
# sims = SimOperation.read_similarities()
# v = Validation(drugs, sims)


################# Generate Similarities #################

sims = SimOperation.sim_table()
SimOperation.write_similarities(sims)

