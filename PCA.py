from sklearn.decomposition import PCA
import numpy as np

def read_sim_matrix():
    sim = open('data/chemicalsimilarity.csv')
    index_dict = {}
    line = sim.readline()
    sim_dict = {}
    num_of_drug = 0
    while line:
        id1, id2, similarity = line[0:-1].split('\t')
        if id1 not in index_dict.keys():
            index_dict[id1] = num_of_drug  # store index of drug
            num_of_drug += 1
        if id2 not in index_dict.keys():
            index_dict[id2] = num_of_drug  # store index of drug
            num_of_drug += 1
        sim_dict[id1, id2] = float(similarity)
        line = sim.readline()
    sim.close()
    matrix = np.zeros([num_of_drug, num_of_drug])
    for key in sim_dict:
        index1 = index_dict[key[0]]
        index2 = index_dict[key[1]]
        matrix[index1, index2] = sim_dict[key]
        matrix[index2, index1] = sim_dict[key]
    return matrix

def read_similarities():
    chemical_sim = open('data/chemicalsimilarity.csv')
    chemical_sim_dict = {}
    line = chemical_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        chemical_sim_dict[db_id1, db_id2] = float(similarity)
        line = chemical_sim.readline()
    chemical_sim.close()

    atc_sim = open('data/atcsimilarity.csv')
    atc_sim_dict = {}
    line = atc_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        atc_sim_dict[db_id1, db_id2] = float(similarity)
        line = atc_sim.readline()
    atc_sim.close()

    go_sim = open('data/gosimilarity.csv')
    go_sim_dict = {}
    line = go_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        go_sim_dict[db_id1, db_id2] = float(similarity)
        line = go_sim.readline()
    go_sim.close()

    dist_sim = open('data/distsimilarity.csv')
    dist_sim_dict = {}
    line = dist_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        dist_sim_dict[db_id1, db_id2] = float(similarity)
        line = dist_sim.readline()
    dist_sim.close()
    
    ligand_sim = open('data/ligandsimilarity.csv')
    ligand_sim_dict = {}
    line = ligand_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        ligand_sim_dict[db_id1, db_id2] = float(similarity)
        line = ligand_sim.readline()
    ligand_sim.close()

    seq_sim = open('data/seqsimilarity.csv')
    seq_sim_dict = {}
    line = seq_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        seq_sim_dict[db_id1, db_id2] = float(similarity)
        line = seq_sim.readline()
    seq_sim.close()
    
    sideeffect_sim = open('data/sideeffectsimilarity.csv')
    sideeffect_sim_dict = {}
    line = sideeffect_sim.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        sideeffect_sim_dict[db_id1, db_id2] = float(similarity)
        line = sideeffect_sim.readline()
    sideeffect_sim.close()
    return atc_sim_dict, chemical_sim_dict, dist_sim_dict, go_sim_dict, ligand_sim_dict, seq_sim_dict, sideeffect_sim_dict

atc_sim, chemical_sim, dist_sim, go_sim, ligand_sim, seq_sim, sideeffect_sim = read_similarities()

pca = PCA(n_components=7)


print(atc_sim.__len__(), chemical_sim.__len__(),dist_sim.__len__(),go_sim.__len__(),ligand_sim.__len__(),seq_sim.__len__(),sideeffect_sim.__len__())

num = 0
for key in atc_sim:
    if key not in dist_sim.keys():
        print(key)
        num += 1
print(num)

intersection = (atc_sim.keys() & chemical_sim.keys() & dist_sim.keys() & go_sim.keys() & ligand_sim.keys() & seq_sim.keys() & sideeffect_sim.keys())
atc_sim = {key: atc_sim[key] for key in intersection}
chemical_sim = {key: chemical_sim[key] for key in intersection}
dist_sim = {key: dist_sim[key] for key in intersection}
go_sim = {key: go_sim[key] for key in intersection}
ligand_sim = {key: ligand_sim[key] for key in intersection}
seq_sim = {key: seq_sim[key] for key in intersection}
sideeffect_sim = {key: sideeffect_sim[key] for key in intersection}

all_sims = []
for key in atc_sim:
    sims = [atc_sim[key], chemical_sim[key], dist_sim[key], go_sim[key], ligand_sim[key], seq_sim[key], sideeffect_sim[key]]
    all_sims.append(sims)

np.array(all_sims)
X = np.array(all_sims)
pca = PCA(n_components=7)
pca.fit(X)
print(pca.components_)
print(pca.explained_variance_ratio_)





