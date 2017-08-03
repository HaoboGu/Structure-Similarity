
def read(filename):
    f = open(filename)
    line = f.readline()
    print(line)
    folds_res = line.split(',')
    f.close()
    return folds_res[1:11]


for i in range(1, 13):
    filename = 'output/psl/auc_all_all_dataset' + str(i) + '.txt'
    folds_res = read(filename)

    lines = []
    write_line = ''
    for item in folds_res:
        write_line = write_line + item + '\t'
    write_line = write_line + '\n'
    lines.append(write_line)

    write_file = 'output.txt'
    with open(write_file, 'a') as f:
        f.writelines(lines)

