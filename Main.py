from rdkit import Chem
from rdkit.DataStructs import cDataStructs as ds
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
import os
import time
# import xlrd
import xlsxwriter
import random


class Similarity:
    def __init__(self):
        # p: weights of similarities  # p is a dict, represents weights of similarities
        self.p = {'MACCS': 0.25, 'ECFP4': 0.25, 'FCFP4': 0.25, 'Topo': 0.25, 'target': 0, 'mechanism': 0}

    def SimilarityMACCS(self, m):  # compute similarity using MACCS fingerprint
        fp1 = MACCSkeys.GenMACCSKeys(m[0])
        fp2 = MACCSkeys.GenMACCSKeys(m[1])
        return ds.TanimotoSimilarity(fp1, fp2)

    def SimilarityECFP4(self, m):  # ECFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2)  # 2 is the radius of the atom environments considered
        fp2 = AllChem.GetMorganFingerprint(m[1], 2)
        return ds.TanimotoSimilarity(fp1, fp2)

    def SimilarityFCFP4(self, m):  # FCFP4 fingerprint
        fp1 = AllChem.GetMorganFingerprint(m[0], 2, useFeatures=True)
        fp2 = AllChem.GetMorganFingerprint(m[1], 2, useFeatures=True)
        return ds.TanimotoSimilarity(fp1, fp2)

    def SimilarityTopo(self, m):  # Topological Fingerprints
        fp1 = FingerprintMols.FingerprintMol(m[0])
        fp2 = FingerprintMols.FingerprintMol(m[1])
        return ds.TanimotoSimilarity(fp1, fp2)

    def getSimilarity(self, m):  # get weighted similarity of 2 compounds, with weights p
        target = 0
        mechanism = 0
        MACCS = self.SimilarityMACCS(m)
        ECFP4 = self.SimilarityECFP4(m)
        FCFP4 = self.SimilarityFCFP4(m)
        Topo = self.SimilarityTopo(m)
        return self.p['MACCS']*MACCS + self.p['ECFP4']*ECFP4 + self.p['FCFP4']*FCFP4 + self.p['Topo']*Topo + self.p['target']*target + self.p['mechanism']*mechanism


class Drug:

    def __init__(self, ID, mol):
        self.ID = ID
        self.mol = mol


class DataPreprocessing:
    dir = 'medicine_structure'  # folder name of drugs

    def __init__(self):
        self.folders = os.listdir(self.dir)
        self.num = self.folders.__len__()
        self.validation_set = []
        self.test_set = []

    def __ReadMolFiles(self):  # Read molfile, for comparing similarity, private
        for i in range(0, self.num):
            filename = self.dir + '/' + self.folders[i] + '/molfile.txt'
            self.validation_set.append(Drug(self.folders[i], Chem.MolFromMolFile(filename)))

    def __DivideData(self):  # divide data into 90% known_set and 10% test_set
        index = random.sample(range(0, 1366), 136)  # randomly select 136 nums in range 0~1365
        for i in index:
            self.test_set.append(self.validation_set[i])  # add to testset
        for item in self.test_set:
            self.validation_set.remove(item)  # remove this item from mols dataset

    def Processing(self):  # Process original data, test_set and validation_set are obtained
        self.__ReadMolFiles()
        self.__DivideData()

    def ReadData(self):
        self.__ReadMolFiles()
        return self.validation_set


class Validation:

    def __init__(self):
        dp = DataPreprocessing()
        dp.Processing()
        self.validation_set = dp.validation_set
        self.test_set = dp.test_set
        self.pairs = {}
        self.similarities = {}
        self.data = []

    def WriteSimilarites(self):
        dp = DataPreprocessing()
        data = dp.ReadData()
        S = Similarity()
        result = []
        # compute similarity
        for d1 in data:
            for d2 in data:
                if d1 != d2:
                    result.append([d1.ID, d2.ID, S.SimilarityTopo([d1.mol, d2.mol]),
                                   S.SimilarityECFP4([d1.mol, d2.mol]), S.SimilarityFCFP4([d1.mol, d2.mol]),
                                   S.SimilarityMACCS([d1.mol, d2.mol]), S.getSimilarity([d1.mol, d2.mol])])

        print(result)  # result is the table of similaries

        # Write the result to .txt file
        f = open('result.txt', 'w')
        for item in result:
            tmp = ''.join([str(i) + ' ' for i in item])
            f.write(tmp + '\n')
        f.close()

    def ReadSimilarities(self):
        self.data = []
        file = open('result')
        while 1:
            line = file.readline()
            if not line:
                break
            self.data.append(line.split())
        return self.data

    def FindConflict(self, drugID):  # Find similar drugs of a drug
        conflict = []
        for item in self.data:
            if item[0] == drugID and item[6] > 0.5:  # temporary 0.5, a more concise criterion is needed
                print(item)
                conflict.append(item)


################## Generate Similarities #################


# Write the result to .xlsx file
# wb = xlsxwriter.Workbook('similarites.xlsx')
# sheet1 = wb.add_worksheet('Data')
# sheet1.write(0, 0, 'Drug1 ID')  # table structure
# sheet1.write(0, 1, 'Drug2 ID')
# sheet1.write(0, 2, 'Similarity1')
# sheet1.write(0, 3, 'Similarity2')
# sheet1.write(0, 4, 'Similarity3')
# sheet1.write(0, 5, 'Similarity4')
# sheet1.write(0, 6, 'Weighted Similarity')
# row = 1
# for item in result:
#     sheet1.write(row, 0, item[0])
#     sheet1.write(row, 1, item[1])
#     sheet1.write(row, 2, item[2])
#     sheet1.write(row, 3, item[3])
#     sheet1.write(row, 4, item[4])
#     sheet1.write(row, 5, item[5])
#     sheet1.write(row, 6, item[6])
#     row += 1
# wb.close()
