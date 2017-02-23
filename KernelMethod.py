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

    # read interacts
    interaction_file = open('data/interacts.csv')
    interact_dict = {}
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        if int(interact_level) == 1:
            interact_dict[db_id1, db_id2] = int(interact_level)  # use multiple keys
        line = interaction_file.readline()
    interaction_file.close()
    inter_matrix = np.zeros([num_of_drug, num_of_drug])
    for key in interact_dict:
        index1 = index_dict[key[0]]
        index2 = index_dict[key[1]]
        inter_matrix[index1, index2] = interact_dict[key]
        inter_matrix[index2, index1] = interact_dict[key]
    return matrix, inter_matrix


def read_result():
    re = open('kernel_re.txt')
    line = re.readline()

    while line:
        print(line)
        re_m = np.array()
        string = line[0:-1].split(',')
        string = list(map(float, string))
        a = np.array(string)
        re_m = np.concatenate((re_m, a))

        line =  re.readline()

m, im = read_sim_matrix()
np.savetxt('im.txt', im, '%f')