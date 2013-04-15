#ifndef NEURALNETWORK_H
#define NEURALNETWORK_H

#include <iostream>
#include <string>
#include <vector>

#include <boost/mpi.hpp>

#include <boost/math/distributions.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
//#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "Neuron.h"
#include "Connection.h"

class NeuralNetwork {
  private:
    static const double connectionProbability  = 0.2;
    static const double excitatoryFraction     = 0.8;

    static const double dt                     = 0.01;           // ms per time step

    static const int    neuronRefractoryTime   = (int) (2.0/0.01); // (2.0/dt) in time steps
    static const double spikeThreshold         = 20.0;           // mV
    static const double neuronRestingPotential = 10.0;           // mv

    const double J_EE;
    const double J_IE;
    const double J_EI;
    const double J_II;
    double Delta;

    static const double tau_E_ms               = 10.0;
    static const double tau_I_ms               = 5.0;

    static const double tau_d                  = 1.0 / 0.01; // 1.0 / dt
    static const double s_d                    = 0.5 / 0.01; // 0.5 / dt

    static const double nu_E_ms                = 13.0/1000.0;
    double              externalInputDistMean;

    int         excitatoryNeuronCount;
    int         inhibitoryNeuronCount;
    int         connectionCount;

    SpikeList localSpikeList;
    SpikeList remoteSpikeList;
    const int maxSpikeDelay;

    typedef std::vector<Neuron*> NeuronPointerList;
    NeuronPointerList neuronList;
    typedef std::vector<LocalNeuron*> LocalNeuronPointerList;
    LocalNeuronPointerList localNeuronList;

    boost::random::mt19937                     rng;
    //boost::random::uniform_01<boost::random::mt19937&, double>                uniformDist;
    boost::random::uniform_real_distribution<> uniformRealDist;
    boost::random::uniform_int_distribution<>  spikeDelayDistribution;

    typedef boost::random::normal_distribution<> RandomNormalDistribution;
    RandomNormalDistribution      distributionJ_EE;
    RandomNormalDistribution      distributionJ_IE;
    RandomNormalDistribution      distributionJ_EI;
    RandomNormalDistribution      distributionJ_II;

    boost::mpi::communicator world;

    int totalSpikingNeuronCount;

  public:
    NeuralNetwork(NeuronPointerList::size_type n, int seed, double j_ee, double j_ie, double j_ei, double j_ii);

    void printNetworkParameters(std::ostream& out);
    void simulate(int timeSteps);
    void writeOutputFiles(std::string connectivityFilePath, std::string activityFilePath);

    int getSpikingNeuronCount() { return totalSpikingNeuronCount; }
};

#endif
