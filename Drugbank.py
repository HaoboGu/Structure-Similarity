import random
import numpy as np
import time
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

def read_drugbank_data():
    # read interaction data
    interaction_file = open('data/interacts.csv')
    interact_dict = {}
    line = interaction_file.readline()
    while line:
        db_id1, db_id2, interact_level = line[0:-1].split('\t')
        interact_dict[db_id1, db_id2] = int(interact_level)  # use multiple keys
        line = interaction_file.readline()
    interaction_file.close()

    # read similarity data
    similarity_file = open('data/chemicalsimilarity.csv')
    similarity_dict = {}
    line = similarity_file.readline()
    while line:
        db_id1, db_id2, similarity = line[0:-1].split('\t')
        similarity_dict[db_id1, db_id2] = float(similarity)
        line = similarity_file.readline()
    similarity_file.close()
    return interact_dict, similarity_dict


class Validation:

    def __init__(self, interact_dict, similarity_dict):
        self.interaction = interact_dict
        self.similarity = similarity_dict
        self.train_set = {}
        self.validation_set = {}
        self.sim_link = {}
        self.positive_train = {}
        self.max_sim_with_positive_link = {}
        self.max_sim_with_positive_link_for_val = {}

    def divide_data(self):
        self.train_set = {}
        self.validation_set = {}
        index = random.sample(range(0, 9892), 989)  # randomly select 1/10 interactions as test_set
        flag = 0
        for i in self.interaction:
            if flag in index:
                self.validation_set[i] = self.interaction[i]
            else:
                self.train_set[i] = self.interaction[i]
            flag += 1
        # create known ddi dict:
        for key in self.train_set:
            if self.train_set[key] == 1:
                self.positive_train[key] = 1


    def compute_link_sim(self, key1, key2):
        link_sim1 = (self.similarity[key1[0], key2[0]] + self.similarity[key1[1], key2[1]]) / 2.0
        link_sim2 = (self.similarity[key1[0], key2[1]] + self.similarity[key1[1], key2[0]]) / 2.0
        return max(link_sim1, link_sim2)

    def create_simlink(self):
        self.sim_link = {}
        # num = 1
        for inter_key in self.train_set:
            max_link_sim = 0
            for inter_key2 in self.positive_train:
                if inter_key[0] in inter_key2 and inter_key[1] in inter_key2:
                    continue
                else:
                    link_sim = self.compute_link_sim(inter_key, inter_key2)
                    if link_sim > max_link_sim:
                        max_link_sim = link_sim
                        self.sim_link[inter_key] = inter_key2
                        self.max_sim_with_positive_link[inter_key] = max_link_sim
            # print('iter', num)
            # num += 1

    def create_simlink_for_val(self):
        self.sim_link = {}
        # num = 1
        for inter_key in self.validation_set:
            max_link_sim = 0
            for inter_key2 in self.positive_train:
                if inter_key[0] in inter_key2 and inter_key[1] in inter_key2:
                    continue
                else:
                    link_sim = self.compute_link_sim(inter_key, inter_key2)
                    if link_sim > max_link_sim:
                        max_link_sim = link_sim
                        # self.sim_link[inter_key] = inter_key2
                        self.max_sim_with_positive_link_for_val[inter_key] = max_link_sim
        sim_list = []
        inter_list = []
        for inter_key in self.validation_set:
            feature = self.max_sim_with_positive_link_for_val[inter_key]
            sim_list.append(feature)
            inter_list.append(self.validation_set[inter_key])
        return sim_list, inter_list

    def create_train_array(self):
        sim_list = []
        inter_list = []
        num = 0
        for inter_key in self.train_set:
            if self.train_set[inter_key] == 1:
                feature = self.max_sim_with_positive_link[inter_key]
                sim_list.append(feature)
                inter_list.append(self.train_set[inter_key])
                num += 1
        print('num of positive samples in train set: ', num)
        num = num * 3
        for inter_key in self.train_set:
            if self.train_set[inter_key] == 0:
                feature = self.max_sim_with_positive_link[inter_key]
                sim_list.append(feature)
                inter_list.append(self.train_set[inter_key])
                num = num - 1
            if num == 0:
                break
        return sim_list, inter_list


    def lr(self, sim_list, inter_list):
        lr = LogisticRegression(solver='sag')
        sim_list = np.array(sim_list)
        sim_list = sim_list.reshape(sim_list.shape[0], 1)
        inter_list = np.array(inter_list)
        inter_list = inter_list.reshape(inter_list.shape[0], 1)
        lr.fit(sim_list, inter_list)
        val_sim, val_inter = self.create_simlink_for_val()
        val_sim = np.array(val_sim)
        val_sim = val_sim.reshape(val_sim.shape[0], 1)
        val_inter = np.array(val_inter).reshape(val_inter.__len__(), 1)
        result = lr.predict(val_sim)
        prob_re = lr.predict_proba(val_sim)
        prob_re = prob_re.transpose()
        auroc = roc_auc_score(val_inter, prob_re[1])
        print('roc score:', auroc)
        return result, prob_re, val_inter

start = time.time()
interact_dict, sim_dict = read_drugbank_data()
v = Validation(interact_dict, sim_dict)
v.divide_data()
v.create_simlink()
sim_list, inter_list = v.create_train_array()

result, prob_re, val_inter = v.lr(sim_list, inter_list)

TP = 0  # predict 1, actual 1
FP = 0  # predict 1, actual 0
TN = 0  # predict 0, actual 0
FN = 0  # predict 0, actual 1
for i in range(0, 989):
    if result[i] == 0 and result[i] == 0:
        TN += 1
    elif result[i] == 0 and val_inter[i] == 1:
        FN += 1
    elif result[i] == 1 and val_inter[i] == 0:
        FP += 1
    elif result[i] == 1 and val_inter[i] == 1:
        TP += 1
print('tp:', TP, ' tn:', TN, ' fp:', FP, ' fn:', FN)
precision = TP / (TP + FP)
recall = TP / (TP + FN)

print('precision:', precision)
print('recall:', recall)
print('f-score: ', 2 * precision * recall / (precision + recall))
end = time.time()
print(end-start)