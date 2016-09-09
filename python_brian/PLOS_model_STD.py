# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:35:55 2014

This script shifts the synaptic depression into two equations for each neuron
rather than computing it for every synapse. This is possible since every
excitatory neuron only has two discrete possibilities for STD.

@author: Grant
"""
#from Spike_Stats import *
import brian as bn
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.io

#def goldman_model(T):
'''
Based on the Lim and Goldman model.
'''

def fft_std(delta_u,run_num,new_connectivity,osc,rep):
  #bn.seed(int(time.time()))
  bn.reinit_default_clock()
  #bn.seed(1412958308+2)
  bn.defaultclock.dt = 0.5*bn.ms
  
  #==============================================================================
  # Define constants for the model.
  #==============================================================================
  fft_file = './std_fft_p20_'
  rate_file = './std_rate_p20_'
  print delta_u
  print run_num
  print new_connectivity
  print rep
  
  if osc:
    T = 5.5 * bn.second
  else:
    T = 2.5 * bn.second
  n_tsteps = T / bn.defaultclock.dt
  fft_start = 0.5 * bn.second / bn.defaultclock.dt  # Time window for the FFT computation
  ro = 1.2*bn.Hz
  
  SEE1 = 1.0
  SEE2 = 1.0
  qee1 = 1.00 # Fraction of NMDA receptors for e to e connections
  qee2 = 0.00
  qie1 = 1.00 # Fraction of NMDA receptors for e to i connections
  qie2 = 0.00
  
  uee1 = 0.2-delta_u
  uee2 = 0.2+delta_u
  uie1 = 0.2
  uie2 = 0.2
  trec1 = 1000.0*bn.ms
  trec2 = 1000.0*bn.ms
  
  k = 0.65
  #Jeo_const = 1.0#*bn.mV # Base strength of o (external) to e connections
  
  Ne = 3200 # number of excitatory neurons
  Ni = 800 # number of inhibitory neurons
  No = 20000 # number of external neurons
  N = Ne + Ni
  
  pcon = 0.2 # probability of connection
  
  Jee = 10.0/(Ne*pcon)
  Jie = 10.0/(Ne*pcon)
  Jii = k*10.0/(Ni*pcon)
  Jei = k*10.0/(Ni*pcon)
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
  dV/dt = (-(V-El)+J_ampa1*I_ampa1+J_nmda1*I_nmda1+J_ampa2*I_ampa2+J_nmda2*I_nmda2-J_gaba*I_gaba+J_input*I_input+eta)/tm : bn.volt
  dI_ampa1/dt = -I_ampa1/t_ampa : bn.volt
  dI_nmda1/dt = -I_nmda1/t_nmda : bn.volt
  dI_ampa2/dt = -I_ampa2/t_ampa : bn.volt
  dI_nmda2/dt = -I_nmda2/t_nmda : bn.volt
  dI_gaba/dt = -I_gaba/t_gaba : bn.volt
  dI_input/dt = (-I_input+mu)/t_input : bn.volt
  dx1/dt = (1-x1)/t1_rec : 1
  dx2/dt = (1-x2)/t2_rec : 1
  u1 : 1
  t1_rec : bn.second
  u2 : 1
  t2_rec : bn.second
  mu : bn.volt
  eta : bn.volt
  J_ampa1 : 1
  J_nmda1 : 1
  J_ampa2 : 1
  J_nmda2 : 1
  J_gaba : 1
  J_input : 1
  tm : bn.second
  t_ampa : bn.second
  t_nmda : bn.second
  t_gaba : bn.second
  t_input : bn.second
  '''
  
  P_reset = "V=-52*bn.mV;x1+=-u1*x1;x2+=-u2*x2"
  
  Se_model = '''
  we_ampa1 : bn.volt
  we_nmda1 : bn.volt
  we_ampa2 : bn.volt
  we_nmda2 : bn.volt
  '''
  
  Se_pre = ('I_ampa1 += x1_pre*we_ampa1','I_nmda1 += x1_pre*we_nmda1',
            'I_ampa2 += x2_pre*we_ampa2','I_nmda2 += x2_pre*we_nmda2')
  
  
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
  Pe.I_ampa1 = 0*bn.mV
  Pe.I_nmda1 = 0*bn.mV
  Pe.I_ampa2 = 0*bn.mV
  Pe.I_nmda2 = 0*bn.mV
  Pe.I_gaba = 0*bn.mV
  Pe.I_input = 0*bn.mV
  Pe.V = (np.random.rand(Pe.V.size)*12-52)*bn.mV
  
  Pe.x1 = 1.0
  Pe.x2 = 1.0
  Pe.u1 = uee1
  Pe.u2 = uee2
  Pe.t1_rec = trec1
  Pe.t2_rec = trec2
  
  
  
  Pi = P[Ne:(Ne+Ni)]
  Pi.tm = ti
  Pi.t_ampa = tie_ampa
  Pi.t_nmda = tie_nmda
  Pi.t_gaba = tii_gaba
  Pi.t_input = teo_input
  Pi.I_ampa1 = 0*bn.mV
  Pi.I_nmda1 = 0*bn.mV
  Pi.I_ampa2 = 0*bn.mV
  Pi.I_nmda2 = 0*bn.mV
  Pi.I_gaba = 0*bn.mV
  Pi.I_input = 0*bn.mV
  Pi.V = (np.random.rand(Pi.V.size)*12-52)*bn.mV
  
  Pi.x1 = 1.0
  Pi.x2 = 1.0
  Pi.u1 = 0.0
  Pi.u2 = 0.0
  Pi.t1_rec = 1.0
  Pi.t2_rec = 1.0
  
  
  Pe.J_ampa1 = Jee*(1-qee1)#*SEE1
  Pe.J_nmda1 = Jee*qee1#*SEE1
  Pe.J_ampa2 = Jee*(1-qee2)#*SEE2
  Pe.J_nmda2 = Jee*qee2#*SEE2
  
  Pi.J_ampa1 = Jie*(1-qie2)#*SEE2
  Pi.J_nmda1 = Jie*qie2#*SEE2
  Pi.J_ampa2 = Jie*(1-qie1)#*SEE1
  Pi.J_nmda2 = Jie*qie1#*SEE1
  
  Pe.J_gaba = Jei
  Pi.J_gaba = Jii
  
  Pe.J_input = Jeo
  Pi.J_input = Jeo
  
  #==============================================================================
  # Define inputs
  #==============================================================================
  
  if osc:
    Pe.mu = 12.0*bn.mV
    holder = np.zeros((n_tsteps,))
    t_freq =   np.linspace(0,10,n_tsteps)
    
    fo = 0.2 # Smallest frequency in the signal
    fe = 10.0 # Largest frequency in the signal
    F = int(fe/0.2)
    for m in range(1,F+1):
      holder = holder + np.cos(2*np.pi*m*fo*t_freq-m*(m-1)*np.pi/F)  
    holder = holder / np.max(holder)
    Pe.eta = bn.TimedArray(0.0*bn.mV * holder)#, dt=0.5*bn.ms)
    
    Background_eo = bn.PoissonInput(Pe, N=1000, rate=1.05*bn.Hz, weight=0.2*bn.mV, state='I_input')
    Background_io = bn.PoissonInput(Pi, N=1000, rate=1.0*bn.Hz, weight=0.2*bn.mV, state='I_input')
    
    Pi.mu = 0*bn.mV
    Pi.eta = 0*bn.mV#, dt=0.5*bn.ms)
    
    Po = bn.PoissonGroup(No,rates=0*bn.Hz)
  else:

    Background_eo = bn.PoissonInput(Pe, N=1000, rate=1.05*bn.Hz, weight=0.2*bn.mV, state='I_input')
    Background_io = bn.PoissonInput(Pi, N=1000, rate=1.0*bn.Hz, weight=0.2*bn.mV, state='I_input')    
    holder_pe = np.zeros((n_tsteps,))
    time_steps = np.linspace(0,T/bn.second,n_tsteps)  
    holder_pe[time_steps < 0.5] = 0.0 *bn.mV
    holder_pe[time_steps >= 0.5] = 6.0 *bn.mV #25
    holder_pe[time_steps > 1.5] = 0.0 *bn.mV #25
    Pe.mu = bn.TimedArray(holder_pe)

    def firing_function(t,ro):
      if t>0.5*bn.second and t<3.5*bn.second:
        return 0.0 * bn.Hz
      else:
        return 0.0*bn.Hz
    
    Pe.eta = 0*bn.mV#, dt=0.5*bn.ms)  
    Pi.mu =  0.0*bn.mV
    Pi.eta = 0*bn.mV#, dt=0.5*bn.ms)
    
    Po = bn.PoissonGroup(No,rates=lambda t:firing_function(t,ro))  
  
  
  
  #==============================================================================
  # Define synapses  
  #==============================================================================
  
  See1 = bn.Synapses(Pe, Pe, model = Se_model, pre = Se_pre)
  See2 = bn.Synapses(Pe, Pe, model = Se_model, pre = Se_pre)
  Sie1 = bn.Synapses(Pe, Pi, model = Se_model, pre = Se_pre)
  Sie2 = bn.Synapses(Pe, Pi, model = Se_model, pre = Se_pre)
  
  Sei = bn.Synapses(Pi, Pe, model = Si_model, pre = Si_pre)
  Sii = bn.Synapses(Pi, Pi, model = Si_model, pre = Si_pre)
  
  Seo = bn.Synapses(Po, Pe, model = So_model, pre = So_pre)
  
  #==============================================================================
  # Define random connections
  #==============================================================================
  
  
  if new_connectivity:
    See1.connect_random(Pe,Pe,sparseness=pcon/2.0)
    See2.connect_random(Pe,Pe,sparseness=pcon/2.0)
    Sie1.connect_random(Pe,Pi,sparseness=pcon/2.0)  
    Sie2.connect_random(Pe,Pi,sparseness=pcon/2.0)
    Sii.connect_random(Pi,Pi,sparseness=pcon)
    Sei.connect_random(Pi,Pe,sparseness=pcon)
    Seo.connect_random(Po,Pe,sparseness=pcon)

  
    print 'Saving'
    See1.save_connectivity('./See1_connections_std_saver_p20_'+str(run_num))
    See2.save_connectivity('./See2_connections_std_saver_p20_'+str(run_num))
    Sie1.save_connectivity('./Sie1_connections_std_saver_p20_'+str(run_num))
    Sie2.save_connectivity('./Sie2_connections_std_saver_p20_'+str(run_num))
    Sii.save_connectivity('./Sii_connections_std_saver_p20_'+str(run_num))
    Sei.save_connectivity('./Sei_connections_std_saver_p20_'+str(run_num))
    Seo.save_connectivity('./Seo_connections_std_saver_p20_'+str(run_num))
  else:
    print 'Loading'
    See1.load_connectivity('./See1_connections_std_saver_p20_'+str(run_num))
    See2.load_connectivity('./See2_connections_std_saver_p20_'+str(run_num))
    Sie1.load_connectivity('./Sie1_connections_std_saver_p20_'+str(run_num))
    Sie2.load_connectivity('./Sie2_connections_std_saver_p20_'+str(run_num))
    Sii.load_connectivity('./Sii_connections_std_saver_p20_'+str(run_num))
    Sei.load_connectivity('./Sei_connections_std_saver_p20_'+str(run_num))
    Seo.load_connectivity('./Seo_connections_std_saver_p20_'+str(run_num))
  
  
  See1.we_ampa1 = SEE1*1.0*bn.mV/tee_ampa
  See1.we_nmda1 = SEE1*1.0*bn.mV/tee_nmda
  See1.we_ampa2 = 0.0*bn.mV/tee_ampa
  See1.we_nmda2 = 0.0*bn.mV/tee_nmda
  
  
  See2.we_ampa1 = 0.0*bn.mV/tee_ampa
  See2.we_nmda1 = 0.0*bn.mV/tee_nmda
  See2.we_ampa2 = SEE2*1.0*bn.mV/tee_ampa
  See2.we_nmda2 = SEE2*1.0*bn.mV/tee_nmda
  
  
  Sie1.we_ampa1 = 0.0*bn.mV/tie_ampa
  Sie1.we_nmda1 = 0.0*bn.mV/tie_nmda
  Sie1.we_ampa2 = SEE1*1.0*bn.mV/tie_ampa
  Sie1.we_nmda2 = SEE1*1.0*bn.mV/tie_nmda
  
  
  Sie2.we_ampa1 = SEE2*1.0*bn.mV/tie_ampa
  Sie2.we_nmda1 = SEE2*1.0*bn.mV/tie_nmda
  Sie2.we_ampa2 = 0.0*bn.mV/tie_ampa
  Sie2.we_nmda2 = 0.0*bn.mV/tie_nmda
  
  
  Sei.wi_gaba = 1.0*bn.mV/tei_gaba
  Sii.wi_gaba = 1.0*bn.mV/tii_gaba
  
  Seo.wo_input = 1.0*bn.mV/teo_input

  
  #==============================================================================
  #  Define monitors
  #==============================================================================
  
  
  Pe_mon_V = bn.StateMonitor(Pe,'V',timestep = 10,record=True)
  Pe_mon_eta = bn.StateMonitor(Pe,'eta',timestep = 1,record=True)
  Pe_mon_ampa1 = bn.StateMonitor(Pe,'I_ampa1',timestep = 1,record=True)
  Pe_mon_nmda1 = bn.StateMonitor(Pe,'I_nmda1',timestep = 1,record=True)
  Pe_mon_ampa2 = bn.StateMonitor(Pe,'I_ampa2',timestep = 1,record=True)
  Pe_mon_nmda2 = bn.StateMonitor(Pe,'I_nmda2',timestep = 1,record=True)
  Pe_mon_gaba = bn.StateMonitor(Pe,'I_gaba',timestep = 1,record=True)
  Pe_mon_input = bn.StateMonitor(Pe,'I_input',timestep = 10,record=True)
  See1_mon_x = bn.StateMonitor(Pe,'x1',timestep = 10,record=True)
  See2_mon_x = bn.StateMonitor(Pe,'x2',timestep = 10,record=True)
    
  Pe_ratemon = bn.PopulationRateMonitor(Pe,bin=10.0*bn.ms)
  Pi_ratemon = bn.PopulationRateMonitor(Pi,bin=10.0*bn.ms)
  
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
  
  
  
  
  Pe_mon_ampa1_vals = Pe.J_ampa1[0]*np.mean(Pe_mon_ampa1.values.T,axis=1)
  Pe_mon_nmda1_vals = Pe.J_nmda1[0]*np.mean(Pe_mon_nmda1.values.T,axis=1)
  Pe_mon_ampa2_vals = Pe.J_ampa2[0]*np.mean(Pe_mon_ampa2.values.T,axis=1)
  Pe_mon_nmda2_vals = Pe.J_nmda2[0]*np.mean(Pe_mon_nmda2.values.T,axis=1)
  Pe_mon_ampa_vals = Pe_mon_ampa1_vals + Pe_mon_ampa2_vals
  Pe_mon_nmda_vals = Pe_mon_nmda1_vals + Pe_mon_nmda2_vals
  
  Pe_mon_gaba_vals = Pe.J_gaba[0]*np.mean(Pe_mon_gaba.values.T,axis=1)
  Pe_mon_input_vals = Pe.J_input[0]*np.mean(Pe_mon_input.values.T,axis=1)
  Pe_mon_V_vals = np.mean(Pe_mon_V.values.T,axis=1)
  
  Pe_mon_all_vals = Pe_mon_ampa_vals+Pe_mon_nmda_vals-Pe_mon_gaba_vals
  
  See1_mon_x_vals = np.mean(See1_mon_x.values.T,axis=1)
  See2_mon_x_vals = np.mean(See2_mon_x.values.T,axis=1)
  

  
  #==============================================================================
  # Save into a Matlab file 
  #==============================================================================
  
  
  if osc:

    Pe_output = Pe.J_ampa1[0]*Pe_mon_ampa1.values+Pe.J_nmda1[0]*Pe_mon_nmda1.values + \
    Pe.J_ampa2[0]*Pe_mon_ampa2.values+Pe.J_nmda2[0]*Pe_mon_nmda2.values-Pe.J_gaba[0]*Pe_mon_gaba.values
    Pe_output = Pe_output[:,fft_start:,]  
    Pe_V = Pe_mon_V.values[:,fft_start:,]
    Pe_glut = Pe.J_ampa1[0]*Pe_mon_ampa1.values+Pe.J_nmda1[0]*Pe_mon_nmda1.values + \
    Pe.J_ampa2[0]*Pe_mon_ampa2.values+Pe.J_nmda2[0]*Pe_mon_nmda2.values
    Pe_glut = Pe_glut[:,fft_start:,]
    Pe_gaba = Pe.J_gaba[0]*Pe_mon_gaba.values
    Pe_gaba = Pe_gaba[:,fft_start:,]
      
    Pe_input = Pe_mon_eta[:,fft_start:,]
    T_step = bn.defaultclock.dt
    
    holder = {'Pe_output':Pe_output,'Pe_input':Pe_input,'Pe_V':Pe_V,'Pe_glut':Pe_glut,'Pe_gaba':Pe_gaba,'T_step':T_step}
    scipy.io.savemat(fft_file+'delta_u'+str(delta_u)+'_'+str(rep), mdict=holder)
  else:
    holder = {'Pe_rate':Pe_ratemon.rate,'Pe_time':Pe_ratemon.times,'uee1':uee1,
              'uee2':uee2,'uie1':uie1,'uie2':uie2}
    scipy.io.savemat(rate_file+'delta_q_'+str(delta_u)+'_'+str(run_num)+'rep'+str(rep), mdict=holder)
  bn.clear(erase=True, all=True)
  
  #==============================================================================
  #  Run model and save output
  #==============================================================================
  
if __name__ == '__main__':
#  delta_q_hold = [-0.08,-0.06,-0.04,-0.02,0,0.15]
  delta_q_hold = [-0.075,-0.07,-0.06,-0.05,-0.025,0.0,0.15]
  runs = [1]
  num_reps = [6,7,8,9,10]
  new_connectivity = False
  osc = True
  for k in num_reps:
    for j in runs:
      for i in delta_q_hold:
        fft_std(i,j,new_connectivity,osc,k)
#        plt.show()
  
  
  
  