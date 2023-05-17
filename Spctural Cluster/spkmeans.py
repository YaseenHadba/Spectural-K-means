import pandas as pd
import numpy as np
import kmeansspk as kmp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("k", type=int)
parser.add_argument("goal", type=str)
parser.add_argument("input", type=str)
args = parser.parse_args()


def help2(vector1, vector2):#calulating the (diffecrence between two vectors)^2
    
    arr = [0 for i in range(len(vector2))]
    for i in range(len(vector2)):
        arr[i] = vector1[i] - vector2[i]
    sum2 = 0
    for i in range(len(vector2)):
        sum2 += pow(arr[i], 2)
    return sum2

def divide(vector, num):  # dividing all elements in a vector
    
    for i in range(len(vector)):
        vector[i] = vector[i] / num
    return vector

def starting_centroids(k, vectors): #intialaizing random centroids
    
    N = len(vectors)
    i=0
    indexes = [0 for i in range(k)]
    np.random.seed(0)
    indexes[0] = np.random.choice(N)
    mu0 = vectors[indexes[0]]
    mu = np.zeros((k, len(vectors[0])))
   
    mu[0] = mu0
    while i < k:
        i+=1
        if(i==k) :break
        D = [0 for s in range(N)]
        for l in range(0,N):
            minimum = help2(vectors[l], mu[0])
            for j in range(0,i):
                calculation = help2(vectors[l], mu[j])
                if calculation < minimum:
                    minimum = calculation
            D[l] = minimum
        sum3 = 0
        for m in range(0,N):
            sum3 = sum3 + D[m]
        p = [0 for i in range(N)]
        for l in range(0,N):
            p[l] = D[l] / sum3
        indexes[i] = np.random.choice(N, p=p)
        mu[i] = vectors[indexes[i]]
    print(','.join(map(str, indexes)))
    return mu, indexes


def printMatrix(matrix):
    
    for i in range(len(matrix)):
        print(','.join([format(matrix[i][j], ".4f") for j in range(len(matrix[i]))]))

def centroids_result(k, max_iter, vectors, initial_centroids): #calculating final centroids and printing them
    
    initial_centroids = initial_centroids.tolist()
    vectors = vectors.tolist()
    centroids = kmp.calcspk(k,max_iter,k, vectors, initial_centroids)
    printMatrix(centroids)


def execute( k,goal, max_iter, vectors):
    if goal == "wam":
        Matrix = kmp.fit(vectors, "wam", len(vectors), len(vectors[0]))
    elif goal == "ddg":
        Matrix = kmp.fit(vectors, "ddg", len(vectors), len(vectors[0]))
    elif goal == "lnorm":
        Matrix = kmp.fit(vectors,"lnorm", len(vectors), len(vectors[0]))
    elif goal == "jacobi":
        Matrix = kmp.fit(vectors, "jacobi", len(vectors), len(vectors[0]))
        for i in range(len(Matrix[0])):
            if(Matrix[0][i]<0 and Matrix[0][i]> -0.00005): Matrix[0][i]=0
    elif goal == "spk":
        TMatrix = kmp.fitgetT(vectors, k, len(vectors), len(vectors[0]))
        TMatrix = np.array(TMatrix)
        k = len(TMatrix[0])
        centroids = starting_centroids(k, TMatrix)[0]    
        centroids_result(k, max_iter,TMatrix, centroids)
    else:
        print("invalid input")
        return
    if goal != 'spk': printMatrix(Matrix)


def valid(input):
    if not (str(input[0]).isdigit()):
        return 0
    if(input[0]<0): return 0
    return 1

#reading files and checking inputs
input_args = [args.k ,args.goal, args.input]
file = open(args.input,"r")
vectors = []
curr_line = 0
for line in file:
    str_line = line.split(",")
    int_line = []
    for i in range(len(str_line)):
        int_line.append(float(str_line[i]))
    vectors.append(int_line)  
    if line != "\n":
        curr_line += 1
        
if(len(vectors) < args.k) :
   print("Invalid Input!")
if(not valid(input_args)):
    print("Invalid Input!")

execute(args.k,args.goal, 300, vectors)
