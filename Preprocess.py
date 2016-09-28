# -*- coding:utf-8 -*-

from rdkit import Chem
from rdkit.DataStructs import cDataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
import os


class SimOperation:  # compute similarities, create sim_table, write sims
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

    @classmethod
    def sim_table(cls):  # create similarity table from mol files
        data = MedicineStructure.read_molfiles()
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
    def write_similarities(sim):  # sim: a list of similarities
        # Write the result to .txt file
        f = open('result.txt', 'w')
        for item in sim:
            tmp = ''.join([str(i) + ' ' for i in item])
            f.write(tmp + '\n')
        f.close()


class MedicineStructure:  # Medicine used to read mol files

    def __init__(self, med_id, mol):
        self.ID = med_id
        self.mol = mol

    @staticmethod
    def read_molfiles():
        path = 'medicine_structure'  # folder name of medicines
        folders = os.listdir(path)
        num = folders.__len__()
        data_set = []
        for i in range(0, num):
            filename = path + '/' + folders[i] + '/molfile.txt'
            data_set.append(MedicineStructure(folders[i], Chem.MolFromMolFile(filename)))
        return data_set
