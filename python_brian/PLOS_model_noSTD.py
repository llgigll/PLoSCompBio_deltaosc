# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:35:55 2014

Compute the FFT of the activity in the Pe population for multiple randomly
connected networks.

@author: Grant
"""

#from Spike_Stats import *
import brian as bn
import time
import numpy as np
import scipy.io
import matplotlib.pyplot as plt

def fft_nostd(qee,run_num,new_connectivity,osc,rep):
  
  #bn.seed(int(time.time()))
#  bn.seed(1412958308+2)
  bn.reinit_default_clock()
  bn.defaultclock.dt = 0.5*bn.ms
  
#==============================================================================
# Define constants for the model.
#==============================================================================
  fft_file = './nostd_fft_p20_'
  rate_file = './nostd_rate_p20_'
  
  if osc:
    T = 8.0 * bn.second
  else:
    T = 3.5 * bn.second
  n_tsteps = T / bn.defaultclock.dt
  fft_start = 3.0 * bn.second / bn.defaultclock.dt  # Time window for the FFT computation
  #run_num = 10
  ro = 1.2*bn.Hz
#==============================================================================
#   Need to do all others besides 0.2 and 0.5
#==============================================================================
  print qee
  print run_num
  print new_connectivity
  print rep
  qie = 0.3 # Fraction of NMDA receptors for e to i connections
  
  k = 0.65
  Jeo_const = 1.0#*bn.mV # Base strength of o (external) to e connections
  
  Ne = 3200 # number of excitatory neurons
  Ni = 800 # number of inhibitory neurons
  No = 2000 # number of external neurons
  N = Ne + Ni
  
  pcon = 0.2 # probability of connection
  
  Jee = 5.0/(Ne*pcon)
  Jie = 5.0/(Ne*pcon)
  Jii = k*5.0/(Ni*pcon)
  Jei = k*5.0/(Ni*pcon)
  Jeo = 1.0
  
  
  El = -60.0*bn.mV # leak reversal potential
  Vreset = -52.0*bn.mV # reversal potential
  Vthresh = -40.0*bn.mV # spiking threshold
  
  tref = 2.0*bn.ms # refractory period
  te = 20.0*bn.ms # membrane time constant of excitatory neurons
  ti = 10.0*bn.ms # membrane time constant of inhibitory neruons
  tee_ampa = 10.0*bn.ms # time const of ampa currents at excitatory neurons
  tee_nmda = 100.0*bn.ms # time const of nmda currents at excitatory neurons
  tie_ampa = 10.0*bn.ms  # time const of ampa currents at inhibitory neurons
  tie_nmda = 100.0*bn.ms # time const of nmda currents at inhibitory neurons
  tii_gaba = 10.0*bn.ms # time const of GABA currents at inhibitory neurons
  tei_gaba = 10.0*bn.ms # time const of GABA currents at excitatory neurons
  teo_input = 100.0*bn.ms
  
#==============================================================================
# Define model structure
#==============================================================================
  
  model = '''
  dV/dt = (-(V-El)+J_ampa*I_ampa+J_nmda*I_nmda-J_gaba*I_gaba+J_input*I_input+eta+eta_corr)/tm : bn.volt
  dI_ampa/dt = -I_ampa/t_ampa : bn.volt
  dI_nmda/dt = -I_nmda/t_nmda : bn.volt
  dI_gaba/dt = -I_gaba/t_gaba : bn.volt
  dI_input/dt = (-I_input+mu)/t_input : bn.volt
  J_ampa : 1
  J_nmda : 1
  J_gaba : 1
  J_input : 1
  mu : bn.volt
  eta : bn.volt
  eta_corr : bn.volt
  tm : bn.second
  t_ampa : bn.second
  t_nmda : bn.second
  t_gaba : bn.second
  t_input : bn.second
  '''
  
  P_reset = "V=-52*bn.mV"
  
  Se_model = '''
  we_ampa : bn.volt
  we_nmda : bn.volt
  '''
  
  Se_pre = ('I_ampa += we_ampa','I_nmda += we_nmda')
  
  
  Si_model = '''
  wi_gaba : bn.volt
  '''
  
  Si_pre = 'I_gaba += wi_gaba'  
  
  So_model = '''
  wo_input : bn.volt
  '''
  
  So_pre = 'I_input += wo_input'
  
#==============================================================================
# Define populations
#==============================================================================
  
  P = bn.NeuronGroup(N, model, threshold = Vthresh, reset = P_reset, refractory = tref)
  
  Pe = P[0:Ne]
  Pe.tm = te
  Pe.t_ampa = tee_ampa
  Pe.t_nmda = tee_nmda
  Pe.t_gaba = tei_gaba
  Pe.t_input = teo_input
  Pe.I_ampa = 0*bn.mV
  Pe.I_nmda = 0*bn.mV
  Pe.I_gaba = 0*bn.mV
  Pe.I_input = 0*bn.mV
  Pe.V = (np.random.rand(Pe.V.size)*12-52)*bn.mV
  
  
  Pi = P[Ne:(Ne+Ni)]
  Pi.tm = ti
  Pi.t_ampa = tie_ampa
  Pi.t_nmda = tie_nmda
  Pi.t_gaba = tii_gaba
  Pi.t_input = teo_input
  Pi.I_ampa = 0*bn.mV
  Pi.I_nmda = 0*bn.mV
  Pi.I_gaba = 0*bn.mV
  Pi.I_input = 0*bn.mV
  Pi.V = (np.random.rand(Pi.V.size)*12-52)*bn.mV
  
  Pe.J_ampa = Jee*(1-qee)#*SEE1
  Pe.J_nmda = Jee*qee#*SEE1
  
  Pi.J_ampa = Jie*(1-qie)#*SEE1
  Pi.J_nmda = Jie*qie#*SEE1
  
  Pe.J_gaba = Jei
  Pi.J_gaba = Jii
  
  Pe.J_input = Jeo
  Pi.J_input = Jeo
  
#==============================================================================
# Define inputs
#==============================================================================
  if osc:
    Pe.mu = 2.0*bn.mV
    holder = np.zeros((n_tsteps,))
    t_freq =   np.linspace(0,10,n_tsteps)
    
    fo = 0.2 # Smallest frequency in the signal
    fe = 10.0 # Largest frequency in the signal
    F = int(fe/0.2)
    for m in range(1,F+1):
      holder = holder + np.cos(2*np.pi*m*fo*t_freq-m*(m-1)*np.pi/F)  
    holder = holder / np.max(holder)
    Pe.eta = bn.TimedArray(0.0*bn.mV * holder)#, dt=0.5*bn.ms)
    Pe.eta_corr = 0*bn.mV
    
    Background_eo = bn.PoissonInput(Pe, N=1000, rate=1.0*bn.Hz, weight=0.2*bn.mV, state='I_input')
    Background_io = bn.PoissonInput(Pi, N=1000, rate=1.05*bn.Hz, weight=0.2*bn.mV, state='I_input')
    
    Pi.mu = 0*bn.mV
    Pi.eta = 0*bn.mV#, dt=0.5*bn.ms)
    Pi.eta_corr = 0 * bn.mV
    
    Po = bn.PoissonGroup(No,rates=0*bn.Hz)
  
  else:
    Background_eo = bn.PoissonInput(Pe, N=1000, rate=1.0*bn.Hz, weight=0.2*bn.mV, state='I_input')
    Background_io = bn.PoissonInput(Pi, N=1000, rate=1.05*bn.Hz, weight=0.2*bn.mV, state='I_input')
    
    holder_pe = np.zeros((n_tsteps,))
    time_steps = np.linspace(0,T/bn.second,n_tsteps)  
    holder_pe[time_steps < 0.5] = 0.0 *bn.mV
    holder_pe[time_steps >= 0.5] = 3.0*bn.mV # 35.0/Jeo *bn.mV #25
    Pe.mu = bn.TimedArray(holder_pe)
    
    def firing_function(t,ro):
      if t>0.5*bn.second and t<3.5*bn.second:
        return 0.0 * bn.Hz
      else:
        return 0.0*bn.Hz
    
#    Pe.mu = 0*bn.mV
    Pe.eta = 0*bn.mV#, dt=0.5*bn.ms)  
    Pi.mu =  0.0*bn.mV
    Pi.eta = 0*bn.mV#, dt=0.5*bn.ms)
    
    Po = bn.PoissonGroup(No,rates=lambda t:firing_function(t,ro))  
  


#==============================================================================
# Define synapses  
#==============================================================================
  
  See = bn.Synapses(Pe, Pe, model = Se_model, pre = Se_pre)
  Sie = bn.Synapses(Pe, Pi, model = Se_model, pre = Se_pre)
  
  Sei = bn.Synapses(Pi, Pe, model = Si_model, pre = Si_pre)
  Sii = bn.Synapses(Pi, Pi, model = Si_model, pre = Si_pre)
  
  Seo = bn.Synapses(Po, Pe, model = So_model, pre = So_pre)  
      
#==============================================================================
#  Define monitors
#==============================================================================
  
  Pe_mon_V = bn.StateMonitor(Pe,'V',timestep = 1,record=True)
  Pe_mon_eta = bn.StateMonitor(Pe,'eta',timestep = 1,record=True)
  Pe_mon_ampa = bn.StateMonitor(Pe,'I_ampa',timestep = 1,record=True)
  Pe_mon_nmda = bn.StateMonitor(Pe,'I_nmda',timestep = 1,record=True)
  Pe_mon_gaba = bn.StateMonitor(Pe,'I_gaba',timestep = 1,record=True)
  Pe_ratemon = bn.PopulationRateMonitor(Pe,bin=10.0*bn.ms)
  

#==============================================================================
# Define random connections
#==============================================================================
  
  
  if new_connectivity:
    See.connect_random(Pe,Pe,sparseness=pcon)
    Sie.connect_random(Pe,Pi,sparseness=pcon)
    Sii.connect_random(Pi,Pi,sparseness=pcon)
    Sei.connect_random(Pi,Pe,sparseness=pcon)
    Seo.connect_random(Po,Pe,sparseness=pcon)
  
    print 'Saving'
    See.save_connectivity('./See_connections_nostd_saver_p20'+str(run_num))
    Sie.save_connectivity('./Sie_connections_nostd_saver_p20'+str(run_num))
    Sii.save_connectivity('./Sii_connections_nostd_saver_p20'+str(run_num))
    Sei.save_connectivity('./Sei_connections_nostd_saver_p20'+str(run_num))
    Seo.save_connectivity('./Seo_connections_nostd_saver_p20'+str(run_num))
  
  else:
    print 'Loading'
    See.load_connectivity('./See_connections_nostd_saver_p20'+str(run_num))
    Sie.load_connectivity('./Sie_connections_nostd_saver_p20'+str(run_num))
    Sii.load_connectivity('./Sii_connections_nostd_saver_p20'+str(run_num))
    Sei.load_connectivity('./Sei_connections_nostd_saver_p20'+str(run_num))
    Seo.load_connectivity('./Seo_connections_nostd_saver_p20'+str(run_num))
  
  
  See.we_ampa = 1.0*bn.mV/tee_ampa
  See.we_nmda = 1.0*bn.mV/tee_nmda
  
  
  Sie.we_ampa = 1.0*bn.mV/tie_ampa
  Sie.we_nmda = 1.0*bn.mV/tie_nmda
  
  Sei.wi_gaba = 1.0*bn.mV/tei_gaba
  Sii.wi_gaba = 1.0*bn.mV/tii_gaba
  
  Seo.wo_input = 1.0*bn.mV/teo_input
  
  
#==============================================================================
# Run model
#==============================================================================
  
  
  timer = 0*bn.second
  t_start = time.time()
  bn.run(T, report='graphical')
  timer = timer + T
  print '-------------------------------------------------------'
  print 'Time is ' + str(timer)+' seconds'
  t_end = time.time()
  print 'Time to compute last ' +str(T)+' seconds is: ' + \
        str(t_end - t_start) + ' seconds'
  print '-------------------------------------------------------\n'
  
 
#==============================================================================
# Save into a Matlab file 
#==============================================================================
  
  if osc:
    Pe_output = Pe.J_ampa[0]*Pe_mon_ampa.values+Pe.J_nmda[0]*Pe_mon_nmda.values-Pe.J_gaba[0]*Pe_mon_gaba.values
    Pe_output = Pe_output[:,fft_start:,]
    Pe_glut = Pe.J_ampa[0]*Pe_mon_ampa.values+Pe.J_nmda[0]*Pe_mon_nmda.values
    Pe_glut = Pe_glut[:,fft_start:,]
    Pe_gaba = Pe.J_gaba[0]*Pe_mon_gaba.values[:,fft_start:,]
     
    Pe_V = Pe_mon_V.values[:,fft_start:,]
    Pe_input = Pe_mon_eta[:,fft_start:,]
    T_step = bn.defaultclock.dt
    
    holder = {'Pe_output':Pe_output,'Pe_input':Pe_input,'Pe_V':Pe_V,'Pe_glut':Pe_glut,'Pe_gaba':Pe_gaba,'T_step':T_step}
    scipy.io.savemat(fft_file+'qee'+str(qee)+'_'+str(rep), mdict=holder)
  
  else:
    holder = {'Pe_rate':Pe_ratemon.rate,'Pe_time':Pe_ratemon.times}
    scipy.io.savemat(rate_file+'qee_'+str(qee)+'_'+str(run_num)+'rep'+str(rep), mdict=holder)


  
#==============================================================================
#  Run the goldman model and save output
#==============================================================================

if __name__ == '__main__':
    
  q_hold = [0.27,0.28,0.29,0.30,0.5]
#  q_hold = [0.5]
  runs = [1]
  num_reps = [1,2,3,4,5,6,7,8,9,10]
  new_connectivity = False
  osc = True
  
  for i in q_hold:
    for j in runs:
      for k in num_reps:
        fft_nostd(i,j,new_connectivity,osc,k)





































