from rdkit import Chem
from rdkit.DataStructs import cDataStructs as ds
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
import os
import time
import xlrd
import xlsxwriter
import random


class Similarity:
    p = {'MACCS': 0.25, 'ECFP4': 0.25, 'FCFP4': 0.25, 'Topo': 0.25, 'target': 0, 'mechanism': 0}  # weights of similarities

    def __init__(self, p):
        self.p = p  # p is a dict, represents weights of similarities

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

    def RealData(self):
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

    def FindPairs(self):  # Find the most similar drugs of test_set in validation_set
        p = {'MACCS': 0.25, 'ECFP4': 0.25, 'FCFP4': 0.25, 'Topo': 0.25, 'target': 0, 'mechanism': 0}
        S = Similarity(p)
        for tested_drug in self.test_set:
            sim = 0
            for known_drug in self.validation_set:
                tmp = S.getSimilarity([tested_drug.mol, known_drug.mol])
                if tmp > sim:
                    sim = tmp
                    self.pairs[tested_drug.ID] = known_drug.ID
            self.similarities[tested_drug.ID] = sim

        # self.pairs = dict(zip(drug1, drug2))

###################################### Main script #############################################





start = time.time()
p = {'MACCS': 0.25, 'ECFP4': 0.25, 'FCFP4': 0.25, 'Topo': 0.25, 'target': 0, 'mechanism': 0}  # weights of similarities
#data = xlrd.open_workbook('interactions.xlsx')
#sheet = data.sheet_by_index(0)
# V = Validation()
# V.FindPairs()

# write excel
# wb = xlwt.Workbook(encoding='ascii')
# sheet1 = wb.add_sheet('Data')
# sheet1.write(0, 0, 'Tested drug')
# sheet1.write(0, 1, 'Most similar drug')
# sheet1.write(0, 2, 'Similarity')
# for row, caption in enumerate(V.pairs):
#     sheet1.write(row+1, 0, caption)
#     sheet1.write(row+1, 1, V.pairs[caption])
#     sheet1.write(row+1, 2, V.similarities[caption])
# wb.save('pairs.xls')

dp = DataPreprocessing()
# dp.Processing()
data = dp.RealData()
S = Similarity(p)


result = []
for d1 in data:
    for d2 in data:
        if d1 != d2:
            result.append([d1.ID, d2.ID, S.getSimilarity([d1.mol, d2.mol])])
# for test in dp.test_set
#     for vali in dp.validation_set:
#         result.append([test.ID, vali.ID, S.getSimilarity([test.mol, vali.mol])])

wb = xlsxwriter.Workbook('sim.xlsx')
sheet1 = wb.add_worksheet('Data')
sheet1.write(0, 0, 'Drug1 ID')
sheet1.write(0, 1, 'Drug2 ID')
sheet1.write(0, 2, 'Similarity')
row = 1
for item in result:
    sheet1.write(row, 0, item[0])
    sheet1.write(row, 1, item[1])
    sheet1.write(row, 2, item[2])
    row += 1
#wb.save('sim.xls')
wb.close()

end = time.time()
print(end - start)
