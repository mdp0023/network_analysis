# this is where I can mess with sparse matrices to see how they are manipulated

from numpy import multiply
from scipy.sparse import csr_matrix
import numpy as np
from scipy import sparse
import time
import sys
from scipy.sparse import random
import scipy.sparse._sparsetools as _sparsetools
from scipy.sparse._sputils import (get_index_dtype, upcast)
from multiprocessing import Pool


data1 = np.array([[0,0,0,0,5],
                  [0,0,0,3,0],
                  [0,0,7,0,0],
                  [4,5,6,7,0],
                  [0,8,0,2,0]])
data2 = np.array([[600, 0, 0, 0, 10],
                  [0, 0, 0, 6, 0],
                  [0, 0, 14, 0, 0],
                  [8, 10, 12, 14, 0],
                  [0, 16, 0, 4, 0]])
data3 = np.array([[17, 0, 0, 0, 5],
                  [0, 0, 0, 3, 0],
                  [0, 0, 7, 0, 0],
                  [4, 12, 19, 7, 0],
                  [0, 8, 0, 15, 0]])



sparse1 = sparse.lil_array(data1)
for i, row in enumerate(sparse1.rows):
    for j in row:
        print(i,j)


sparse1 = sparse.csr_array(data1,dtype='float32')
sparse2 = sparse.csr_array(data2,dtype='float32')
sparse3 = sparse.csr_array(data3, dtype='float32')
addition = sparse.csr_array(
    ([1 for i in sparse2.data], sparse2.indices, sparse2.indptr))
alpha_array=0.15
beta_array=4





weight_array_iter = sparse3.multiply(sparse1.multiply(
    sparse2.power(-1)).power(beta_array).multiply(alpha_array)+addition)
# print(weight_array_iter)


weight_array_iter = sparse3 * \
    (1+alpha_array*(sparse1/sparse2)**beta_array)
# print(weight_array_iter)


# individual math test
sparse1 = sparse.rand(15000, 15000, density=0.0002, format='lil')
sparse2 = sparse.rand(15000, 15000, density=0.0002, format='lil')
sparse3 = sparse.rand(15000, 15000, density=0.0002, format='lil')
addition = sparse.rand(15000, 15000, density=0.0002, format='lil')

# n = [[int(OD_matrix[n][0]), OD_matrix[n][1], int(OD_matrix[n][2]), 'path', 'weight_array'] for n in range(0, len(OD_matrix))]
n = [[0,1,2,'string','string',sparse1,sparse1] for n in range(0,10086)]

# unique_links=[[400,600],[500,70],[90,4],[2345,9824],[10234,19],[9543,25],[0,0],[1,2],[12345,7644],[2345,9876]]
# for unique_link in unique_links:
#     sparse1[unique_link[0], unique_link[1]]=10
#     sparse2[unique_link[0], unique_link[1]]=10
#     sparse3[unique_link[0], unique_link[1]]=10
#     addition[unique_link[0], unique_link[1]] = 10

# addition=sparse.csr_array(addition)


# time1 = time.time()
# sparse1=sparse.csr_array(sparse1)
# sparse2=sparse.csr_array(sparse2)
# sparse3 = sparse.csr_array(sparse3)
# weight_array_iter = sparse3.multiply(sparse1.multiply(sparse2.power(-1)).power(beta_array).multiply(alpha_array)+addition)
# time2=time.time()

# sparse1=sparse.lil_array(sparse1)
# sparse2=sparse.lil_array(sparse2)
# sparse3 = sparse.lil_array(sparse3)

# time3=time.time()
# weight_array_iter = sparse.lil_array(sparse3)
# for unique_link in unique_links:
#     weight_array_iter[unique_link[0], unique_link[1]] = sparse3[
#         unique_link[0], unique_link[1]] * (1+alpha_array*(sparse1[unique_link[0], unique_link[1]]/sparse2[unique_link[0], unique_link[1]])**beta_array)
# time4=time.time()

# print(time2-time1)
# print(time4-time3)



# # parallelization attempt
# # Create two sparse matrices
# sparse1 = sparse.rand(15000, 15000, density=0.0002, format='csr')
# sparse2 = sparse.rand(15000, 15000, density=0.0002, format='csr')
# sparse3 = sparse.rand(15000, 15000, density=0.0002, format='csr')
# addition = sparse.rand(15000, 15000, density=0.0002, format='csr')

# # sparse1 = sparse.lil_array(data1,dtype='float32')
# # sparse2 = sparse.lil_array(data2,dtype='float32')
# # sparse3 = sparse.lil_array(data3, dtype='float32')
# # addition = addition.tolil()


# # Define the number of processes to use
# num_processes = 2

# def multiply_matrices(args):
#     sparse1a, sparse2a, sparse3a, additiona = args
#     c = sparse1a.multiply(sparse2a.power(-1)).power(beta_array).multiply(alpha_array)
#     d=c+additiona
#     weight_array_iter = sparse3a.multiply(d)
#     return weight_array_iter

# time1=time.time()
# # Create a pool of processes
# pool = Pool(processes=num_processes)
# # Split the matrices into chunks based on the number of processes
# chunks = [(sparse1[i::num_processes], 
#            sparse2[i::num_processes],
#            sparse3[i::num_processes],
#            addition[i::num_processes]) for i in range(num_processes)]
# time1a=time.time()
# # Use the pool to parallelize the matrix multiplication
# results = pool.map(multiply_matrices, chunks)
# time1b=time.time()
# # Combine the results into a single sparse matrix
# result = sparse.vstack(results)
# # Close the pool of processes
# pool.close()
# pool.join()

# time2=time.time()

# c = sparse1.multiply(sparse2.power(-1)).power(beta_array).multiply(alpha_array)
# d = c+addition
# weight_array_iter = sparse3.multiply(d)
# time3=time.time()


# print(time1a-time1)
# print(time1b-time1a)
# print(time2-time1b)
# print(time2-time1)
# print(time3-time2)



