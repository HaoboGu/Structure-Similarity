# -*- coding:utf-8 -*-

import xlrd


def ProcessCMedc():
    wb = xlrd.open_workbook('data/CMedc.xls')
    sheet = wb.sheet_by_index(0)
    chn_med = []
    for index in range(0, sheet.nrows):
        row = sheet.row_values(index)
        for j in range(3, 6):
            row[j] = row[j].replace(' ', '')  # remove spaces
        chn_med.append(row)

    file = open('CMedc.txt', 'w')
    for item in chn_med:
        line = ''.join([str(s) + ' ' for s in item])
        file.write(line + '\n')
    file.close()



def ProcessInteraction():
    wb = xlrd.open_workbook('data/interactions.xlsx')
    sheet = wb.sheet_by_index(0)
    interactions = []
    for index in range(0, sheet.nrows):
        row = sheet.row_values(index)
        row[2] = row[2].replace(' ', '')  # remove spaces in drug1's name
        row[4] = row[4].replace(' ', '')  # remove spaces in drug2's name
        interactions.append(row)
    file = open('Interactions.txt', 'w')
    for item in interactions:
        line = ''.join([str(s) + ' ' for s in item])
        file.write(line + '\n')
    file.close()

def ProcessWstMedc():
    wb = xlrd.open_workbook('data/WMedc.xls')
    sheet = wb.sheet_by_index(0)
    wst_med = []
    for index in range(0, sheet.nrows):
        row = sheet.row_values(index)
        row[1] = row[1].replace(' ', '')  # remove spaces in western drug's name
        wst_med.append(row)
    file = open('WMedc.txt', 'w')
    for item in wst_med:
        line = ''.join([str(s) + ' ' for s in item])
        file.write(line + '\n')
    file.close()


# ProcessCMedc()
# ProcessInteraction()
# ProcessWstMedc()
