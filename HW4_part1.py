#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:15:58 2019

@author: mateogomez
"""



import numpy as np
import sys
import csv

from scipy.stats import kstest

# open csv file
def read_data(table):
	input_file = sys.argv[1]

    
    
    with open(input_file, 'r') as f:
        fd = csv.reader(f,delimiter=",")
        fd = list(fd)
    #table = np.array(fd[:,2:len(fd[0])-1], dtype=float)
    for c in range(len(fd[0])):
        tmp = []
        if c == 1:
            continue
        for r in range(1,len(fd)):
            if fd[r][c] == '':
                tmp.append(0)
            else:
                tmp.append(float(fd[r][c]))
        table.append(np.array(tmp))

def MovingP(pricing, alpha, MovingMat):
# get moving average 
    for row in pricing:
        line = []
        for t in range(len(row)):
            if t == 0:
                line.append(row[t])
            else:
                line.append(alpha*line[t-1] + (1-alpha)*row[t])
        MovingMat.append(line)

def mean_std_kstest(MovingMat, means, stds, kstests):
    groups = [0,94,189,284,378,472,567,662,756]

    for asset in MovingMat:
        mean_each_asset = []
        std_each_asset = []
        kstest_each_asset = []
        for i in range(len(groups)-1):
            three_month_moving_average = asset[groups[i]: groups[i+1]]
            mean = np.mean(three_month_moving_average)
            std = np.std(three_month_moving_average)
            mean_each_asset.append(mean)
            std_each_asset.append(std)
            if std != 0:
                three_month_moving_average = (three_month_moving_average-mean)/std
            kstest_each_asset.append(kstest(three_month_moving_average, 'norm'))
        kstests.append(kstest_each_asset)
        means.append(mean_each_asset)
        stds.append(std_each_asset)

def main():
	price = []
	read_data(price)
	alpha = 0.5
	MovingMat = []
	MovingP(price, alpha, MovingMat)
	MovingMat = np.array(MovingMat)

	means = []
	stds = []
	kstests = []
	
	mean_std_kstest(MovingMat, means, stds, kstests)

	valid_count = 0
	asset_count = 0
	for i in kstests:
		for j in i:
			asset_count += 1
			if j[0] <= 0.05:
				valid_count += 1
	print("among ", asset_count, " assets, we have: ", valid_count, " assets can be assumed as normal distribution.")

if __name__ == "__main__":
	main()

