import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import re

# Import data (Pandas is nicer than with open...)
data = pd.read_csv('./data/Spinal cord C Test.txt', sep=',')
ranges = np.asanyarray(data[['Latest Ref Rng & Units']].T)[0, 1:]
data = data.drop('Latest Ref Rng & Units', axis=1).replace(' <0.2', 0.1)
lowerBound = []
upperBound = []
for i in range(0, len(ranges)):
    ranges[i] = ranges[i].replace('<=', '0-')
    ranges[i] = re.sub(r"[^\d.-]+", "", ranges[i])
    if ranges[i] == '': ranges[i] = '0'
    lower = float(ranges[i].partition('-')[0])
    lowerBound.append(lower)
    upper = ranges[i].partition('-')[2]
    if upper == '': upper = lower
    upper = float(upper)
    upperBound.append(upper)

# Fill in missings with neighbor observations:
data = data.fillna(method='ffill', axis=1)

dataArr = np.asanyarray(data)
for i in range(1,len(dataArr)):
    for j in range(1,len(dataArr[1])):
        dataArr[i][j] = float(dataArr[i][j])

def convert24(time):
    time = time.strip(' ')
    if 'AM' in time:
        stripped = time.strip(' AM')
        part = stripped.partition(':')
        if int(part[0]) == 12:
            part[0] = 0
        return int(part[0])*60 + int(part[2])
    elif time[:2] == '12':
        stripped = time.strip(' PM')
        part = stripped.partition(':')
        return int(part[0])*60 + int(part[2])
    elif 'PM' in time:
        part = time.partition(':')
        return (int(part[0]) + 12)*60 + int(part[2].strip(' PM'))

minutes = dataArr[0].copy()[1:]
for i in range(0, len(minutes)):
    minutes[i] = convert24(minutes[i])

times = list(data.columns)[1:]

# Get earliest date (last entry in this application)
times[len(times)-1] = datetime.strptime(times[len(times)-1].strip(' '), "%m/%d/%Y")

# Loop on rest of dates and compare to the first one, also includes hour & minute differences on those days.
for i in range(0, len(times)-1):
    times[i] = abs(times[len(times)-1] - datetime.strptime(times[i].strip(' '), "%m/%d/%Y")).days * 24 * 60
    times[i] += (minutes[i] - minutes[len(times)-1])

# Set earliest entry as 0 for the starting point
times[len(times)-1] = 0
times = np.asarray(times).astype(float)

dataArr = np.delete(dataArr, (0), axis=0)
labels = dataArr[:,:1].reshape(1, len(dataArr))[0]
dataArr = dataArr[:,1:].astype(float)

means = {}
# Compute row means
for i in range(0, len(dataArr)):
    means[labels[i]] = np.average(dataArr[i])

devs = {}

# Compute row standard deviations
for i in range(0, len(dataArr)):
    devs[labels[i]] = np.std(dataArr[i])

def remove_outliers(dataArr):
    ret = []
    test = np.copy(dataArr)
    for i in range(0,len(dataArr)):
        mean = np.average(dataArr[i])
        std = np.std(dataArr[i])
        vals = []
        for j in range(len(dataArr[i])):
            if dataArr[i,j] > mean + 3 * std:
                vals.append(float(mean + 3 * std)) 
            elif dataArr[i,j] < mean - 3 * std:
                vals.append(float(mean - 3 * std)) 
            else:
                vals.append(float(dataArr[i,j]))
        ret.append(vals)
    return np.asarray(ret, dtype='float64')

# Remove outlier points
dataArr = remove_outliers(dataArr)

    
#Sliding window
# This function computes the sliding window average. 
# Expanded the work for our new data size of 33
def moving_average(data):
    noStr = data #define noStr to be the sodium 
    window_size = 3
    moving_averages = [] #array to keep the new "long wave" values
    length = len(noStr)

    ##Checking to see if the variables are defined properly
    #print(len(noStr))
    #print(len(np.flip(times)))

    # Loop through the array t o
    #consider every window of size 3
    i = 0
    while  i < len(noStr)-window_size+1:

        # Calculate the average of current window
        window_average = round(np.sum(noStr[i:i+window_size]) / window_size, 2)

        # Store the average of current
        # window in moving average list
        moving_averages.append(window_average)

        # Shift window to right by one position
        i += 1

    #print(noStr[23:26]) # pulling from the orginal data to identify the 2 values we need to solve for < t_n -> t_n-2 >   
    moving_averages.append(noStr[length-2]) # First assign f(t_24) from no
    moving_averages.append(noStr[length-1])
    moving_averages[length-2] = round((moving_averages[length-3]+moving_averages[length-1])/2,2)
    return(moving_averages) 

# This function adds data points between indices 23 and 25
def interpData(data):
    data = np.flip(data)
    new_length = 11 
    interp= np.linspace(data[23], data[25], new_length)
    newData = np.delete(data, (23,24,25))
    newData = np.append(newData, interp)
    return newData

# Interpolating the data for each biomarker
interpDataArr = []
for i in range(0,len(dataArr)):
    interpD = interpData(dataArr[i])
    interpDataArr.append(interpD)



newTime = interpData(times)

slidingData = []
for i in range(0, len(interpDataArr)):
    #noStr1 = np.flip(interpDataArr[i])
    move_avg = moving_average(interpDataArr[i])
    slidingData.append(move_avg)


#slidingData = []
#for i in range(0, len(dataArr)):
    #noStr1 = dataArr[i]
    #move_avg = moving_average(noStr1)
    #slidingData.append(move_avg)
    


cor = np.corrcoef(dataArr)
for i in range(0,19):
    for j in range(0,19):
        if (i != j and abs(cor[i,j] > 0.8)):
            print(labels[i], ',', labels[j])
            print(cor[i,j])

#sodium = dataArr[0]
# noStr = dataArr[0]
# plt.plot(times, noStr, 'bo-')
# plt.axhline(y = lowerBound[0], color = 'r', linestyle = 'dotted',label = 'Min')
# plt.axhline(y = upperBound[0], color = 'b', linestyle = 'dotted',label = 'Min')
# plt.title(labels[0])

# plt.figure(figsize=(40,40))
# for i in range(1,19):
#     noStr = dataArr[i]
#     plt.subplot(6,3,i)
#     plt.axhline(y = lowerBound[i], color = 'r', linestyle = 'dotted',label = 'Min')
#     plt.axhline(y = upperBound[i], color = 'b', linestyle = 'dotted',label = 'Min')
#     plt.plot(times, noStr,'bo-')
#     plt.title(labels[i])