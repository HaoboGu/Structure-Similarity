# import os
# import random
#
#
# def read_all_interacts(filename, interact_lst):  # read interactions from a file
#     interaction_file = open(filename)
#     line = interaction_file.readline()  # read lines from original interaction file
#     while line:
#         db_id1, db_id2, interact_level = line[0:-1].split('\t')  # get two drugs' ids and interaction level
#         interact_lst.append((db_id1, db_id2, interact_level))  # create a list of tuples (id1, id2, interaction)
#         line = interaction_file.readline()  # read next line
#     interaction_file.close()
#     return interact_lst
#
# def two_folds():
#     interact_lst = []
#     for i in range(1, 11):
#         filename = 'data/all/' + str(i) + '/interacts.csv'
#         lst = read_all_interacts(filename, interact_lst)
#     random.shuffle(lst)
#     samples_each_fold = 49455  # number of samples in each fold：9891*10 = 98910 = 315*314
#     interact_lists = []  # a list of 10 interaction lists
#     for i in range(0, 2):  # 10 folds
#         fold = lst[i * samples_each_fold:samples_each_fold * (i + 1)]  # interactions in each fold
#         interact_lists.append(fold)
#
#     # write to dataset
#     j = 33
#     for i in range(1, 3):
#         interact_list = interact_lists[i - 1]  # interaction list of fold i
#         path = 'data/all_dataset' + str(j) + '/all/' + str(i)  # path: data/all_datasetj/all/(1~10)
#         if not os.path.exists(path):
#             os.makedirs(path)
#
#         pos = []  # list of positive interactions
#         neg = []  # list of negative interactions
#         all = []  # list of all interactions
#         ids = []  # list of ids
#         # read all interacts
#         for item in interact_list:
#             line = item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\n'  # create line to write in interaction file
#             if item[2] == '1':  # positive interactions
#                 pos.append(line)
#                 all.append(line)
#                 line = item[0] + '\t' + item[1] + '\n'  # create line to write in ids file
#                 ids.append(line)
#             elif item[2] == '0':  # negative interactions
#                 neg.append(line)
#                 all.append(line)
#                 line = item[0] + '\t' + item[1] + '\n'
#                 ids.append(line)
#
#         filename = path + '/interacts.csv'  # file of all interactions
#         with open(filename, 'w') as f:
#             f.writelines(all)
#         filename = path + '/interacts_negatives.csv'  # file of negative interactions
#         with open(filename, 'w') as f:
#             f.writelines(neg)
#         filename = path + '/interacts_positives.csv'  # file of positive interactions
#         with open(filename, 'w') as f:
#             f.writelines(pos)
#         filename = path + '/interactsids.csv'  # file of ids
#         with open(filename, 'w') as f:
#             f.writelines(ids)
#
# def exchange():
#     # 读取两个文件夹，从中间切开交换
#     f1 = 'data/all_dataset33/all/1/interacts.csv'
#     f2 = 'data/all_dataset33/all/2/interacts.csv'
#     interacts_1 = []
#     interacts_1 = read_all_interacts(f1, interacts_1)
#     interacts_2 = []
#     interacts_2 = read_all_interacts(f2, interacts_2)
#     new_1 = interacts_1[0:24727]
#     new_1.extend(interacts_2[24727:49455])
#     new_2 = interacts_2[0:24727]
#     new_2.extend(interacts_1[24727:49455])
#
#     interact_lists = []
#     interact_lists.append(new_1)
#     interact_lists.append(new_2)
#     # write to dataset
#     j = 34
#     for i in range(1, 3):
#         interact_list = interact_lists[i - 1]  # interaction list of fold i
#         path = 'data/all_dataset' + str(j) + '/all/' + str(i)  # path: data/all_datasetj/all/(1~10)
#         if not os.path.exists(path):
#             os.makedirs(path)
#
#         pos = []  # list of positive interactions
#         neg = []  # list of negative interactions
#         all = []  # list of all interactions
#         ids = []  # list of ids
#         # read all interacts
#         for item in interact_list:
#             line = item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\n'  # create line to write in interaction file
#             if item[2] == '1':  # positive interactions
#                 pos.append(line)
#                 all.append(line)
#                 line = item[0] + '\t' + item[1] + '\n'  # create line to write in ids file
#                 ids.append(line)
#             elif item[2] == '0':  # negative interactions
#                 neg.append(line)
#                 all.append(line)
#                 line = item[0] + '\t' + item[1] + '\n'
#                 ids.append(line)
#
#         filename = path + '/interacts.csv'  # file of all interactions
#         with open(filename, 'w') as f:
#             f.writelines(all)
#         filename = path + '/interacts_negatives.csv'  # file of negative interactions
#         with open(filename, 'w') as f:
#             f.writelines(neg)
#         filename = path + '/interacts_positives.csv'  # file of positive interactions
#         with open(filename, 'w') as f:
#             f.writelines(pos)
#         filename = path + '/interactsids.csv'  # file of ids
#         with open(filename, 'w') as f:
#             f.writelines(ids)
# # Generate folds
# interact_lst = []
# for i in range(1, 11):  # from original dataset data/all/(1~10)/interacts.csv read all interaction data
#     filename = 'data/all/' + str(i) + '/interacts.csv'
#     lst = read_all_interacts(filename, interact_lst)
#
#
#
#
# # shuffle the list
# # If this line is commented out, the order of interactions will be same as the original dataset
# # If keep this line, the result will be wrong
# random.shuffle(lst)
#
# # divide list
# samples_each_fold = 9891  # number of samples in each fold：9891*10 = 98910 = 315*314
# interact_lists = []  # a list of 10 interaction lists
# for i in range(0, 10):  # 10 folds
#     fold = lst[i*samples_each_fold:samples_each_fold*(i+1)]  # interactions in each fold
#     interact_lists.append(fold)
#
# # write to dataset
# j = 3
# for i in range(1, 11):
#     interact_list = interact_lists[i - 1]  # interaction list of fold i
#     path = 'data/all_dataset'+ str(j) +'/all/' + str(i)  # path: data/all_datasetj/all/(1~10)
#     if not os.path.exists(path):
#         os.makedirs(path)
#
#     pos = []  # list of positive interactions
#     neg = []  # list of negative interactions
#     all = []  # list of all interactions
#     ids = []  # list of ids
#     # read all interacts
#     for item in interact_list:
#         line = item[0] + '\t' + item[1] + '\t' + str(item[2]) + '\n'  # create line to write in interaction file
#         if item[2] == '1':  # positive interactions
#             pos.append(line)
#             all.append(line)
#             line = item[0] + '\t' + item[1] + '\n'  # create line to write in ids file
#             ids.append(line)
#         elif item[2] == '0':  # negative interactions
#             neg.append(line)
#             all.append(line)
#             line = item[0] + '\t' + item[1] + '\n'
#             ids.append(line)
#
#     filename = path + '/interacts.csv'  # file of all interactions
#     with open(filename, 'w') as f:
#         f.writelines(all)
#     filename = path + '/interacts_negatives.csv'  # file of negative interactions
#     with open(filename, 'w') as f:
#         f.writelines(neg)
#     filename = path + '/interacts_positives.csv'  # file of positive interactions
#     with open(filename, 'w') as f:
#         f.writelines(pos)
#     filename = path + '/interactsids.csv'  # file of ids
#     with open(filename, 'w') as f:
#         f.writelines(ids)


# end

import os
import sys
import itertools as it
import pandas as pd
import numpy as np
import csv



def create(dataset_name):
	#Change this line accordingly based on where you place this script
	data_dir = os.path.join('data', 'all')

	interactions = pd.DataFrame()

	for i in range(1, 11):
		interacts_file = os.path.join(data_dir, str(i), 'interacts.csv')

		interacts = set()
		with open(interacts_file, 'r') as f:
			reader = csv.reader(f, delimiter='\t')
			for row in reader:
				if (row[1], row[0], row[2]) not in interacts:
					interacts.add(tuple(row))

		interacts = [list(t) for t in interacts]
		ints_df = pd.DataFrame(interacts, columns=['d1', 'd2', 'val'])
		ints_df['fold'] = i
		interactions = pd.concat([interactions, ints_df])

	interactions['fold'] = np.random.permutation(interactions.fold)


	for i in range(1, 11):
		out = os.path.join(data_dir, '..', '..', dataset_name, 'all', str(i))

		if not os.path.exists(out):
			os.makedirs(out)

		valid_interacts = interactions.loc[interactions.fold == i, ['d1', 'd2', 'val']]
		symmetric = valid_interacts[['d2', 'd1', 'val']]
		symmetric = symmetric.rename(columns={'d2':'d1', 'd1':'d2'})

		ints = pd.concat([valid_interacts, symmetric])
		ids = ints[['d1', 'd2']]

		ints.to_csv(os.path.join(out, 'interacts.csv'), sep='\t', header=False, index=False)
		ids.to_csv(os.path.join(out, 'interactsids.csv'), sep='\t', header=False, index=False)


for i in range(1,31):
	dataset_name = 'all_dataset' + str(i)
	create(dataset_name)
