"""
This is a Python 2.7 program.

This program uses PyIMSL, which requires Python 2.6 but seems to 
work with Python 2.7.

usage:
module load pyimsl
[jlynch@krathys python]$ python netsim.py
"""
import heapq
import itertools
import operator
import random
import sys

import numpy


class Neuron:
  def __init__(self, V_t):
    self.V_t = V_t
    self.V_t_spikes = 0.0
    self.refractoryTimeRemaining = int(0)
    self.outgoingConnectionList = []

class ExcitatoryNeuron(Neuron):
  Jext = 0.0
  decayCoefficient = 0.0
  globalActivity = 0.0
  def __init__(self, V_t):
    Neuron.__init__(self, V_t)

  def getJIndex(self):
    return 0
 
  def getDecayCoefficient(self):
    return ExcitatoryNeuron.decayCoefficient
 
  def getExternalSpikeInput(self, spikeCount):
    return ExcitatoryNeuron.Jext * spikeCount
 
  def accumulateGlobalActivity(self):
    ExcitatoryNeuron.globalActivity += self.V_t

class InhibitoryNeuron(Neuron):
  Jext = 0.0
  decayCoefficient = 0.0
  globalActivity = 0.0
  def __init__(self, V_t):
    Neuron.__init__(self, V_t)

  def getJIndex(self):
    return 1
 
  def getDecayCoefficient(self):
    return InhibitoryNeuron.decayCoefficient
 
  def getExternalSpikeInput(self, spikeCount):
    return InhibitoryNeuron.Jext * spikeCount
 
  def accumulateGlobalActivity(self):
    InhibitoryNeuron.globalActivity += self.V_t

class Spike:
  def __init__(self, targetNeuron, I):
    self.targetNeuron = targetNeuron
    self.arrivalTime = int(-1)
    self.I = I

class Connection:
  def __init__(self, targetNeuron, sourceNeuron, efficacy, delay):
    self.targetNeuron = targetNeuron
    self.sourceNeuron = sourceNeuron
    self.efficacy = efficacy
    self.delay = delay
    self.spike = Spike(targetNeuron, I = self.efficacy)
    self.sourceNeuron.outgoingConnectionList.append(self)

class NeuralNetwork:
  def __init__(self, neuronCount):
    self.dt = 0.01 # ms per time step

    self.spikeThreshold = 20.0 # mV
    self.neuronRefractoryTime = int(2.0 / self.dt)
    self.neuronRestingPotential = 10.0 # mV

    self.Delta = 0.1
    self.J_EE = 0.21
    self.J_EE_width = self.Delta * self.J_EE
    self.J_EI = 3.0 * self.J_EE
    self.J_EI_width = self.Delta * self.J_EI
    self.J_IE = 2 * self.J_EE
    self.J_IE_width = self.Delta * self.J_IE
    self.J_II = 3 * self.J_IE
    self.J_II_width = self.Delta * self.J_II

    self.tau_d = 1.0 / self.dt
    self.s_d = 0.5 / self.dt

    self.tau_E_ms = 10.0
    self.tau_E_ts = self.tau_E_ms / self.dt
    self.tau_I_ms = 5.0
    self.tau_I_ts = self.tau_I_ms / self.dt

    self.N = neuronCount
    self.c = 0.2
    self.excitatoryFraction = 0.8
    self.N_E = int(self.N * self.excitatoryFraction)
    self.N_I = self.N - self.N_E

    self.nu_E_ms = 13.0/1000.0
    self.externalInputLambda = self.c * self.N_E * self.nu_E_ms * self.dt

    self.neuronList = []
    self.neuronList.extend(
      [ExcitatoryNeuron(numpy.random.normal(self.neuronRestingPotential, 5.0)) 
        for i in xrange(self.N_E)]
    )
    self.neuronList.extend(
      [InhibitoryNeuron(numpy.random.normal(self.neuronRestingPotential, 5.0))
        for i in xrange(self.N_I)]
    )
    # seed the numpy rng
    self.J = [[(self.J_EE, self.J_EE_width), (self.J_IE, self.J_IE_width)],[(-1.0*self.J_EI, self.J_EI_width), (-1.0*self.J_II, self.J_II_width)]]
    self.connectionList = [self.buildConnection(targetNeuron, sourceNeuron) 
      for (targetNeuron, sourceNeuron) in itertools.permutations(self.neuronList, 2) if numpy.random.uniform() < self.c]

    self.spikeList = []
    self.spikeHeap = []

    ExcitatoryNeuron.decayCoefficient = 1.0 - self.dt/self.tau_E_ms
    ExcitatoryNeuron.Jext = self.J_EE
    InhibitoryNeuron.decayCoefficient = 1.0 - self.dt/self.tau_I_ms
    InhibitoryNeuron.Jext = self.J_IE


  def printNetworkParameters(self):
    print('simulation time step dt        : {}'.format(self.dt))
    print('external input lambda          : {}'.format(self.externalInputLambda))
    print('excitatory time constant tau_E : {} ms'.format(self.tau_E_ms))
    print('excitatory time constant tau_E : {} time steps'.format(self.tau_E_ts))
    print('inhibitory time constant tau_I : {} ms'.format(self.tau_I_ms))
    print('inhibitory time constant tau_I : {} time steps'.format(self.tau_I_ts))
    print('J_EE                           : {}'.format(self.J_EE))
    print('J_IE                           : {}'.format(self.J_IE))
    print('J_EI                           : {}'.format(self.J_EI))
    print('J_II                           : {}'.format(self.J_II))
    print('tau_d                          : {} time steps'.format(self.tau_d))
    print('s_d                            : {} time steps'.format(self.s_d))
    print('neuron refractory time         : {} time steps'.format(self.neuronRefractoryTime))
    print('excitatory neuron decay coefficient: {}'.format(ExcitatoryNeuron.decayCoefficient))
    print('inhibitory neuron decay coefficient: {}'.format(InhibitoryNeuron.decayCoefficient))


  def buildConnection(self, targetNeuron, sourceNeuron):
    (J, J_width) = self.J[targetNeuron.getJIndex()][sourceNeuron.getJIndex()]
    #print('J: {} J_width: {}'.format(J, J_width))
    efficacy = numpy.random.normal(J, J_width)
    delay = int(numpy.random.uniform(self.tau_d - self.s_d, self.tau_d + self.s_d))
    return Connection(targetNeuron, sourceNeuron, efficacy = efficacy, delay = delay)


  def simulate(self, timeSteps):
    with open('network_simulation_data.txt', 'w') as outputFile:
      for t in xrange(timeSteps):
        ExcitatoryNeuron.globalActivity = 0.0
        InhibitoryNeuron.globalActivity = 0.0
        # accumulate local spikes
        #deliveredSpikeList = []
        #for spike in self.spikeList:
        #  if spike.arrivalTime == t:
        #    spike.targetNeuron.V_t_spikes += spike.I
        #    spike.arrivalTime = int(-1)
        #    deliveredSpikeList.append(spike)
        #  else:
        #    # the spikes are in sorted order so we can stop accumulating as soon
        #    # as we see one that arrives later than now
        #    break
        
        while len(self.spikeHeap) > 0 and self.spikeHeap[0][0] == t:
          (arrivalTime, spike) = heapq.heappop(self.spikeHeap)
          spike.targetNeuron.V_t_spikes += spike.I
          spike.arrivalTime = int(-1)
          #deliveredSpikeList.append(spike)
          
        #print('  delivered {} spikes'.format(len(deliveredSpikeList)))
        #totalInternalSpikeArrivalCount = len(deliveredSpikeList)
        #for spike in deliveredSpikeList:
        #  self.spikeList.remove(spike)
            
        # apply leakage
        # deliver external spikes
        totalExternalSpikeCount = int(0)
        totalInternalSpikingNeuronCount = int(0)
        globalActivity = 0.0
        newSpikeCount = int(0)
        for neuron in self.neuronList:
          if neuron.refractoryTimeRemaining > 0:
            neuron.refractoryTimeRemaining -= 1;
          else:
            neuron.V_t *= neuron.getDecayCoefficient()
            externalSpikeCount = numpy.random.poisson(self.externalInputLambda)
            totalExternalSpikeCount += externalSpikeCount
            neuron.V_t += neuron.getExternalSpikeInput(externalSpikeCount)
            neuron.V_t += neuron.V_t_spikes
            neuron.V_t_spikes = 0.0
            if neuron.V_t > self.spikeThreshold:
              totalInternalSpikingNeuronCount += 1
              neuron.refractoryTimeRemaining = self.neuronRefractoryTime
              neuron.V_t = self.neuronRestingPotential
              for connection in neuron.outgoingConnectionList:
                connection.spike.arrivalTime = t + connection.delay
                #self.spikeList.append(connection.spike)
                heapq.heappush(self.spikeHeap, (connection.spike.arrivalTime, connection.spike))
                newSpikeCount += 1
          globalActivity += neuron.V_t
          neuron.accumulateGlobalActivity()
        # sort the spike list by arrival time
        #if newSpikeCount > 0:
        #  self.spikeList.sort(key=operator.attrgetter('arrivalTime'))

        outputFile.write('{} {} {} {}\n'.format(t, totalInternalSpikingNeuronCount, ExcitatoryNeuron.globalActivity, InhibitoryNeuron.globalActivity))

        # print some numbers
        if t % 100 == 0:
          print('time step {}'.format(t))
          print('  total external spike count         : {}'.format(totalExternalSpikeCount))
          print('  total internal spiking neurons     : {}'.format(totalInternalSpikingNeuronCount))
          #print('  total internal spike arrivals      : {}'.format(totalInternalSpikeArrivalCount))
          print('  global activity                    : {} mV'.format(globalActivity))
          print('  excitatory global activity         : {} mV'.format(ExcitatoryNeuron.globalActivity))
          print('  inhibitory global activity         : {} mV'.format(InhibitoryNeuron.globalActivity))


def netsim():
  # get number of neurons and time steps from command line
  neuronCount = int(sys.argv[1])
  timeStepCount = int(sys.argv[2])
  nn = NeuralNetwork(neuronCount)
  nn.printNetworkParameters()
  print('neural network has')
  print('  {} excitatory neurons'.format(nn.N_E))
  print('  {} inhibitory neurons'.format(nn.N_I))
  print('  {} connections'.format(len(nn.connectionList)))

  print('simulating {} time steps'.format(timeStepCount))
  nn.simulate(timeStepCount)

if __name__ == "__main__":
  netsim()
