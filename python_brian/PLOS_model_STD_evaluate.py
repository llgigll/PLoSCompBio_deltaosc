# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 15:48:05 2015

Compare the transfer functions of individual neurons across multiple values of
q2.  Tests whether or not individual neurons would be a good gauge of the 
stability of the network.

@author: grant
"""
import brian as bn
import time
import numpy as np
import numpy.matlib 
import matplotlib.pyplot as plt
import pickle
import scipy.io

num_files = 10
delta_u_base = 0
delta_u = [-0.075,-0.07,-0.06,-0.05,-0.025,0.0,0.15]
num_cells = 3200

hold=scipy.io.loadmat('./std_fft_p20_delta_u' + str(delta_u[0]) + '_' + str(1) + '.mat')
T_step = hold['T_step']
freqs = np.squeeze(np.fft.fftfreq(hold['Pe_output'][:,2000:].shape[1],T_step))
fo = 3.0
fend = 8.0
inds_int = (freqs >= fo) & (freqs <= fend)


x_fft_total = []
t_func_total = []
t_func_cell_ave = np.zeros((num_cells,np.size(delta_u)));
t_func_int = np.zeros((num_files*num_cells,np.size(delta_u)))
t_func_other = np.zeros((num_files*num_cells,np.size(delta_u)))
#t_func_intbase = np.zeros((num_files*num_cells,1))
for j,delta_u_hold in enumerate(delta_u):
  t_func_mean = []
  x_fft_all = []
  for i in range(1,num_files+1):
    hold =  hold=scipy.io.loadmat('./std_fft_p20_delta_u' + str(delta_u_hold) + '_' + str(i) + '.mat')
    x = hold['Pe_glut'][:,2000:]
    x_fft = np.fft.fft(x)
    
    t_func = np.abs(x_fft)
    
    t_func_hold = t_func[:,inds_int]
    t_func_int[(i-1)*num_cells:i*num_cells,j] = np.trapz(t_func_hold,dx = freqs[1]-freqs[0],axis=1)
    t_func_cell_ave[:,j] = t_func_cell_ave[:,j] +  t_func_int[(i-1)*num_cells:i*num_cells,j]

    
    x_fft_all.append(np.mean(np.abs(x_fft),axis=0))
    t_func_mean.append(np.mean(t_func,axis=0))
    
  x_fft_total.append(sum(x_fft_all)/float(num_files))
  t_func_total.append(sum(t_func_mean)/float(num_files))
  if delta_u_hold == delta_u_base:
    val = np.trapz(t_func_total[j][inds_int],dx = freqs[1]-freqs[0])
#print val

t_func_cell_ave = t_func_cell_ave / (num_files*1.0)
hold = t_func_cell_ave - np.tile(t_func_cell_ave[:,5],(len(delta_u),1)).T 
t_func_int_mean = np.mean(hold,axis=0)
t_func_int_std = np.std(hold,axis=0)



inds = (freqs >= 0.2)# & (freqs <= 10.0)

holder = {'t_func_total':t_func_total,'output_fft_total':x_fft_total,'t_func_cell_ave':t_func_cell_ave,
          't_func_int_mean':t_func_int_mean,'t_func_int_std':t_func_int_std,
          't_func_int':t_func_int,'freqs':freqs,'inds':inds}
scipy.io.savemat('PNAS_goldman_std_tfunc_p20', mdict=holder)

plt.subplot(231)
for i in range(np.size(delta_u)):
  plt.plot(freqs[inds],t_func_total[i][inds])
plt.title('Transfer Function Magnitude')
plt.xlabel('Magnitude (abs)')
plt.ylabel('Frequency (Hz)')

plt.subplot(232)
for i in range(np.size(delta_u)):
  plt.plot(freqs[inds],x_fft_total[i][inds])
plt.title('Output FFT')
plt.xlabel('Magnitude (abs)')
plt.ylabel('Frequency (Hz)')


plt.subplot(233)
plt.errorbar(delta_u,t_func_int_mean,yerr=t_func_int_std)
plt.title('Mean Distance from q=0.30 Transfer Function')
plt.xlabel('qee')
plt.ylabel('Mean Distance')

plt.subplot(234)
plt.hist((t_func_cell_ave[:,0]-t_func_cell_ave[:,4])/10.0)
plt.title('FFT as a function of q')

plt.show()









