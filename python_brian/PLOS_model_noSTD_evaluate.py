# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:48:05 2015

Compare the transfer functions of individual neurons across multiple values of
q2.  Tests whether or not individual neurons would be a good gauge of the 
stability of the network.

This is a redone FFT for a constant input rather than a Schroeder multisine.

@author: grant
"""
import brian as bn
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

num_files = 10
q_base = 0.3
q = np.array([0.26,0.27,0.28,0.29,0.3,0.5])
#q = np.array([0.26])
num_cells = 3200

# Load appropriate files
hold =  hold=scipy.io.loadmat('nostd_fft_p20_qee' + str(q[0]) + '_' + str(1) + '.mat')
T_step = hold['T_step']

# Determine frequencies for the FFT and integral 
freqs = np.squeeze(np.fft.fftfreq(hold['Pe_output'][:,4000:].shape[1],T_step))
fo = 0.5
fend = 5.5
inds_int = (freqs >= fo) & (freqs <= fend)

# Compute FFT and take the integral under the FFT
x_fft_total = []
t_func_total = []
t_func_cell_ave = np.zeros((num_cells,q.shape[0]));
t_func_int = np.zeros((num_files*num_cells,q.shape[0]))
t_func_other = np.zeros((num_files*num_cells,q.shape[0]))
for j,q_hold in enumerate(q):
  t_func_mean = []
  x_fft_all = []
  for i in range(1,num_files+1):
    hold =  hold=scipy.io.loadmat('./nostd_fft_p20_qee' + str(q_hold) + '_' + str(i) + '.mat')
    x = hold['Pe_glut'][:,4000:]
    x_fft = np.fft.fft(x)
    
    t_func = np.abs(x_fft)
    
    t_func_hold = t_func[:,inds_int]
    t_func_int[(i-1)*num_cells:i*num_cells,j] = np.trapz(t_func_hold,dx = freqs[1]-freqs[0],axis=1)
    t_func_cell_ave[:,j] = t_func_cell_ave[:,j] +  t_func_int[(i-1)*num_cells:i*num_cells,j]
    
    x_fft_all.append(np.mean(np.abs(x_fft),axis=0))
    t_func_mean.append(np.mean(t_func,axis=0))
    
    

  x_fft_total.append(sum(x_fft_all)/float(num_files))
  t_func_total.append(sum(t_func_mean)/float(num_files))

t_func_cell_ave = t_func_cell_ave / (num_files*1.0)
hold = t_func_cell_ave - np.tile(t_func_cell_ave[:,4],(6,1)).T 
t_func_int_mean = np.mean(hold,axis=0)
t_func_int_std = np.std(hold,axis=0)
inds = (freqs >= 1.0) & (freqs <= 15.0)



holder = {'t_func_total':t_func_total,'output_fft_total':x_fft_total,'t_func_cell_ave':t_func_cell_ave,
          't_func_int_mean':t_func_int_mean,'t_func_int_std':t_func_int_std,
          't_func_int':t_func_int,'freqs':freqs,'inds':inds}
scipy.io.savemat('nostd_tfunc_p20_test', mdict=holder)





plt.subplot(231)
for i in range(q.shape[0]):
  plt.plot(freqs[inds],t_func_total[i][inds])
plt.title('Transfer Function Magnitude')
plt.xlabel('Magnitude (abs)')
plt.ylabel('Frequency (Hz)')

plt.subplot(232)
for i in range(q.shape[0]):
  plt.plot(freqs[inds],x_fft_total[i][inds])
plt.title('Output FFT')
plt.xlabel('Magnitude (abs)')
plt.ylabel('Frequency (Hz)')


plt.subplot(233)
plt.errorbar(q,t_func_int_mean,yerr=t_func_int_std)
plt.title('Mean Distance from q=0.30 Transfer Function')
plt.xlabel('qee')
plt.ylabel('Mean Distance')

#plt.subplot(234)
#plt.hist((t_func_cell_ave[:,0]-t_func_cell_ave[:,4])/10.0,bins=20)
#plt.title('FFT as a function of q')

plt.show()








