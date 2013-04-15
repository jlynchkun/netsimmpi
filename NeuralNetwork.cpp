#include <fstream>
#include <sstream>

#include <boost/mpi/collectives.hpp>

#include "NeuralNetwork.h"

NeuralNetwork::NeuralNetwork(const NeuronPointerList::size_type neuronCount, int seed, double j_ee, double j_ie, double j_ei, double j_ii) :
  J_EE(j_ee),
  J_IE(j_ie),
  J_EI(j_ei),
  J_II(j_ii),
  //rng(world.rank()), this compiles but does it work? does not seem to work
  rng(),
  //uniformDist(rng),
  uniformRealDist(0.0, 1.0),
  spikeDelayDistribution(tau_d - s_d, tau_d + s_d),
  distributionJ_EE(J_EE, 0.1*J_EE),
  distributionJ_IE(J_IE, 0.1*J_IE),
  distributionJ_EI(-J_EI, 0.1*J_EI),
  distributionJ_II(-J_II, 0.1*J_II),
  neuronList(),
  localNeuronList(),
  maxSpikeDelay(int(tau_d + s_d)),
  localSpikeList(),
  remoteSpikeList(),
  totalSpikingNeuronCount(0) {

    std::cout << "process " << world.rank() << " using seed " << seed << std::endl;
    rng.seed(seed); // this works

    std::cout << "maximum spike delay: " << maxSpikeDelay << std::endl;
    std::cout << "creating local spike list" << std::endl;
    for ( int i = 0; i < maxSpikeDelay; i++ ) {
      localSpikeList.push_back(new std::vector<double>());
    }

    // only need one remote spike list
    std::cout << "creating remote spike list" << std::endl;
    for ( int r = 0; r < 1; r++ ) {
      remoteSpikeList.push_back(new std::vector<double>());
    }

    externalInputDistMean = connectionProbability * neuronCount * excitatoryFraction * nu_E_ms * dt;

    ExcitatoryNeuron::Jexternal = J_EE;
    InhibitoryNeuron::Jexternal = J_IE;
    ExcitatoryNeuron::decayCoefficient = 1.0 - dt/tau_E_ms;
    InhibitoryNeuron::decayCoefficient = 1.0 - dt/tau_I_ms;

    std::cout << "creating " << neuronCount << " neurons" << std::endl;    
    excitatoryNeuronCount = 0;
    inhibitoryNeuronCount = 0;
    for ( int i = 0; i < neuronCount; i++ ) {
      // use the mpi process size and rank to determine which neurons are local
      // use i % mpiProcessSize
      //  i   i % mpiProcessSize=1  i % mpiProcessSize=2  i % mpiProcessSize=3
      //  0     0                     0                     0
      //  1     0                     1                     1
      //  2     0                     0                     2
      //  9     0                     1                     0
      // 10     0                     0                     1
      // if i % mpiProcessSize == (mpiProcessRank-1) then neuron i is local
      // the process rank of neuron i is easy to predict
      //   processRankNeuronI = (i % mpiProcessSize)
      //   if processRankNeuronI == mpiProcessRank then neuron i is local

      int mpiProcessRankOfNeuronI = (i % world.size());
      //if ( i < 10 ) {
      //  std::cout << "neuron " << i << " has mpiProcessRank " << mpiProcessRankOfNeuronI << std::endl;
      //}
      if ( mpiProcessRankOfNeuronI == world.rank() ) {
        // neuron i is local
        //if ( i < 10 ) {
          //std::cout << "neuron " << i << " in process with rank " << world.rank() << " is local" << std::endl;
        //}
        if ( uniformRealDist(rng) < excitatoryFraction ) {
          ExcitatoryNeuron* en = new ExcitatoryNeuron(i, mpiProcessRankOfNeuronI, &localSpikeList);
          neuronList.push_back(en);
          localNeuronList.push_back(en);
          excitatoryNeuronCount++;
        }
        else {
          InhibitoryNeuron* in = new InhibitoryNeuron(i, mpiProcessRankOfNeuronI, &localSpikeList);
          neuronList.push_back(in);
          localNeuronList.push_back(in);
          inhibitoryNeuronCount++;
        }
      }
      else {
        //if (i < 10) {
        //  std::cout << "neuron " << i << " in process with rank " << world.rank() << " is remote" << std::endl;
        //}
        // neuron i is remote
        RemoteNeuron* rn = new RemoteNeuron(i, mpiProcessRankOfNeuronI, &remoteSpikeList);
        neuronList.push_back(rn);
      }
    }

    // send the type id of all local neurons to the other processes
    std::vector<int> neuronTypeIdBuffer;
    for ( int i = 0; i < localNeuronList.size(); i++ ) {
      neuronTypeIdBuffer.push_back(localNeuronList.at(i)->getId());
      neuronTypeIdBuffer.push_back(localNeuronList.at(i)->getTypeId());
    }
    std::vector< std::vector<int> > remoteNeuronTypeIdBufferList;
    boost::mpi::all_gather(world, neuronTypeIdBuffer, remoteNeuronTypeIdBufferList);
    std::vector< std::vector<int> >::const_iterator x = remoteNeuronTypeIdBufferList.begin();
    for ( int allgather_rank = 0; allgather_rank < world.size(); allgather_rank++ ) {
      std::cout << "process " << world.rank()
                << " reading remote neuron type id buffer from mpi process with rank " 
                << allgather_rank << std::endl;
      // do not do anything with this process's neurons
      if ( allgather_rank == world.rank() ) {
        std::cout << "skipping local neurons in process " << world.rank() << std::endl;
      }
      else {
        std::vector<int> remoteNeuronTypeIdBuffer = remoteNeuronTypeIdBufferList[allgather_rank];
        std::vector<int>::const_iterator y = remoteNeuronTypeIdBuffer.begin();
        while ( y != remoteNeuronTypeIdBuffer.end() ) {
          int remoteNeuronId = *y;
          y++;
          int remoteNeuronTypeId = *y;
          y++; 
          //std::cout << "process with rank " << world.rank()
          //          << " setting remote neuron " << remoteNeuronId 
          //          << " type to " << remoteNeuronTypeId
          //          << std::endl;
          RemoteNeuron* remoteNeuron = (RemoteNeuron*) neuronList[remoteNeuronId];
          neuronList[remoteNeuronId]->setTypeId(remoteNeuronTypeId);
        }
      }
    }

    std::cout << "creating connections" << std::endl;
    std::cout << localNeuronList.size() << " local neurons" << std::endl;
    std::cout << neuronList.size() << " total neurons" << std::endl;
    RandomNormalDistribution* distributionJ_xx[2][2];
    distributionJ_xx[0][0] = &distributionJ_EE;
    distributionJ_xx[1][0] = &distributionJ_IE;
    distributionJ_xx[0][1] = &distributionJ_EI;
    distributionJ_xx[1][1] = &distributionJ_II;
    int connectionCounts[2][2];
    connectionCounts[0][0] = 0;
    connectionCounts[1][0] = 0;
    connectionCounts[0][1] = 0;
    connectionCounts[1][1] = 0;

    connectionCount = 0;
    int localConnectionCount = 0;
    int remoteConnectionCount = 0;
    for (int i = 0; i < localNeuronList.size(); i++) {
      LocalNeuron* sourceNeuron = localNeuronList[i];
      for (int j = 0; j < neuronList.size(); j++) {
        Neuron* targetNeuron = neuronList[j];
        if (sourceNeuron->getId() == targetNeuron->getId()) {
          // do not connect a neuron to itself
        }
        else if ( uniformRealDist(rng) < connectionProbability ) {
          connectionCount++;

          int delay = spikeDelayDistribution(rng);
          int targetNeuronTypeId = targetNeuron->getTypeId();
          int sourceNeuronTypeId = sourceNeuron->getTypeId();
          connectionCounts[targetNeuronTypeId][sourceNeuronTypeId]++;
          RandomNormalDistribution* efficacyDistribution = distributionJ_xx[targetNeuron->getTypeId()][sourceNeuron->getTypeId()];
          double efficacy = (*efficacyDistribution)(rng);
          sourceNeuron->addConnection(
            new Connection(sourceNeuron, targetNeuron, efficacy, delay)
          );
          if ( sourceNeuron->getMpiRank() == targetNeuron->getMpiRank() ) {
            localConnectionCount++;
          }
          else {
            remoteConnectionCount++;
          }
        }
        else {
         // no connection between these neurons
        }
      }
    }

    std::cout << "process " << world.rank() << " : " << localConnectionCount << " local connections" << std::endl;
    std::cout << "process " << world.rank() << " : " << remoteConnectionCount << " remote connections" << std::endl;

    std::cout << "process " << world.rank() << " : " << connectionCounts[0][0] << " ee connections" << std::endl;
    std::cout << "process " << world.rank() << " : " << connectionCounts[1][0] << " ie connections" << std::endl;
    std::cout << "process " << world.rank() << " : " << connectionCounts[0][1] << " ei connections" << std::endl;
    std::cout << "process " << world.rank() << " : " << connectionCounts[1][1] << " ii connections" << std::endl;
    std::cout << "process " << world.rank() << " NeuralNetwork initialized" << std::endl;

}

void NeuralNetwork::printNetworkParameters(std::ostream& out) {
  out << "           neuron count: " << neuronList.size() << std::endl;
  out << "        neuron capacity: " << neuronList.capacity() << std::endl;
  out << "excitatory neuron count: " << excitatoryNeuronCount << std::endl;
  out << "inhibitory neuron count: " << inhibitoryNeuronCount << std::endl;
  out << "       connection count: " << connectionCount << std::endl;
}

void NeuralNetwork::simulate(int timeSteps) {
  std::stringstream debugFilePath;
  debugFilePath << "netsimmpi-debug-process-" << world.rank() << ".txt";
  std::ofstream debugFile(debugFilePath.str().c_str());

  // debug
  for ( int i = 0; i < 10; i++) {
    Neuron* neuron = localNeuronList[i];
    debugFile << neuron->getTypeId() << ' ';
  }
  debugFile << std::endl;
  // debug

  boost::random::poisson_distribution<> externalSpikeDistribution(externalInputDistMean);

  int remoteSpikeCount = 0;

  int verbose = 0;

  for (int t = 1; t <= timeSteps; t++) {
    if ( t % 10000 == 0 ) {
    //if (t % 1 == 0) {
      std::cout << "process " << world.rank() << " t: " << t * dt << " ms total spiking neuron count: " << totalSpikingNeuronCount << std::endl;
    }
    // handle local spikes
    // deliver spikes arriving now
    // spikes from remote neurons have already been added to the appropriate list
    int spikingNeuronCount = 0;
    int receivedSpikeCount = 0;
    int spikeListIndex = t % maxSpikeDelay;
    if (verbose) std::cout << "spike list index " << spikeListIndex << std::endl;
    std::vector<double>* currentSpikeList = localSpikeList[spikeListIndex];
    if ( ! currentSpikeList->empty() ) {
      int c = 0;
      while (c < currentSpikeList->size()) {
        //std::cout << "  connection " << c << std::endl;
        // original (local) implementation
        //Connection* connection = currentSpikeList->at(c);
        //connection->target->receiveSpike(connection->efficacy);
        int neuronId = currentSpikeList->at(c);
        c++;
        double efficacy = currentSpikeList->at(c);
        c++;
        neuronList[neuronId]->receiveSpike(efficacy);
        receivedSpikeCount++;
        if (verbose) std::cout << "t " << t << " rank " << world.rank() << ": local neuron " << neuronId << " in world.rank " << world.rank() << " received spike with efficacy " << efficacy << std::endl;
      }
      currentSpikeList->clear();
    }
    if (verbose) std::cout << "t " << t << " rank " << world.rank() << ": " << receivedSpikeCount << " spike(s) received " << std::endl;
    //std::cout << "updating neuron activity" << std::endl;
    for ( int n = 0; n < localNeuronList.size(); n++ ) {
      //std::cout << "  n " << n << std::endl;
      LocalNeuron* neuron = localNeuronList.at(n);
      if ( neuron->refractoryTimeRemaining > 0 ) {
        neuron->refractoryTimeRemaining--;
        neuron->setActivity(neuronRestingPotential);
      }
      else {
        neuron->updateActivity(externalSpikeDistribution(rng));
        if ( neuron->getActivity() > spikeThreshold ) {
          // count number of neurons that spike this time step
          spikingNeuronCount++;
          totalSpikingNeuronCount++;
          if (verbose) std::cout << "t " << t << " rank " << world.rank() << ": local neuron " << neuron->getId() << " spiked" << std::endl;
          neuron->spikeNow(t, maxSpikeDelay);
          neuron->refractoryTimeRemaining = neuronRefractoryTime;
        }
      }
      // debug
      if ( n < 10 ) {
        debugFile << neuron->getActivity() << ' ';
      }
    }
    debugFile << std::endl;
    if (verbose) std::cout << "t " << t << " rank " << world.rank() << ": " << spikingNeuronCount << " neuron(s) spiked" << std::endl;

    //std::cout << "all_gather" << std::endl;
    // send spikes to remote neurons
    std::vector< std::vector<double> > remoteSpikeProcessList;
    if (verbose) std::cout << "t " << t << " rank " << world.rank() << " before all_gather remoteSpikeProcessList has size " << remoteSpikeProcessList.size() << std::endl;
    if (verbose) std::cout << "t " << t << " rank " << world.rank() << " sending remoteSpikeList[0] with size " << remoteSpikeList[0]->size() << std::endl;
    boost::mpi::all_gather(world, *remoteSpikeList[0], remoteSpikeProcessList);
    if (verbose) std::cout << "t " << t << " rank " << world.rank() << " after all_gather remoteSpikeProcessList has size " << remoteSpikeProcessList.size() << std::endl;

    for ( int allgather_rank = 0; allgather_rank < world.size(); allgather_rank++ ) {
      //std::cout << "process " << world.rank()
      //          << " reading remote spike list from mpi process with rank " 
      //          << allgather_rank << std::endl;
      // do not do anything with this process's neurons
      if ( allgather_rank == world.rank() ) {
        if (verbose) std::cout << "t " << t << " rank " << world.rank() << " skipping local spikes from process with rank " << allgather_rank << std::endl;
      }
      else {
        std::vector<double> _remoteSpikeList = remoteSpikeProcessList[allgather_rank];
        if (verbose) std::cout << "t " << t << " rank " << world.rank() << " remote spike list from process with rank " << allgather_rank << " has " << _remoteSpikeList.size() / 4 << " spikes (list size " << _remoteSpikeList.size() << ")" << std::endl;
        int s = 0;
        while ( s < _remoteSpikeList.size() ) {
          int sourceNeuronId = _remoteSpikeList[s];
          s++;
          int targetNeuronId = _remoteSpikeList[s];
          s++;
          double efficacy = _remoteSpikeList[s];
          s++;
          int arrivalTimeStep = _remoteSpikeList[s];
          s++;
          if ( arrivalTimeStep <= t ) {
            std::cout << "arrival time step " << arrivalTimeStep << " is less than current time step " << t << std::endl;
          }
          if (verbose) std::cout << "t " << t << " rank " << world.rank() << ": remote neuron " << sourceNeuronId << " sent a spike to neuron " << targetNeuronId << " with efficacy " << efficacy << " and arrival time " << arrivalTimeStep << std::endl;
          if ( (targetNeuronId % world.size()) == world.rank() ) {
            // the target neuron is local
            if ( (sourceNeuronId % world.size() == world.rank()) ) {
              std::cout << "source neuron is local !!!" << std::endl;
            }
            else {
              if (verbose) std::cout << "  neuron " << targetNeuronId << " is local to process with rank " << world.rank() << std::endl;
              LocalNeuron* localNeuron = (LocalNeuron*) neuronList[targetNeuronId];
              localNeuron->scheduleSpikeArrival(sourceNeuronId, arrivalTimeStep, maxSpikeDelay, efficacy);
              remoteSpikeCount++;
            }
          }
          else {
            // this neuron is not local
            if (verbose) std::cout << "  neuron " << targetNeuronId << "is NOT local to process with rank " << world.rank() << std::endl;            
          }
        }
      }
    }
    remoteSpikeList[0]->clear();
  }
  std::cout << "process " << world.rank() << " has remote spike count " << remoteSpikeCount << std::endl;
  std::cout << "process " << world.rank() << " has total spiking neuron count " << totalSpikingNeuronCount << std::endl;
}

void NeuralNetwork::writeOutputFiles(std::string connectivityFilePath, std::string activityFilePath) {
  std::ofstream connectivityFile(connectivityFilePath.c_str());
  std::ofstream activityFile(activityFilePath.c_str());

  for (int n = 0; n < localNeuronList.size(); n++) {
    LocalNeuron* neuron = localNeuronList.at(n);
    connectivityFile << neuron->getId() << ' ' << neuron->getTypeId();
    for (int c = 0; c < neuron->getOutgoingConnections().size(); c++) {
      connectivityFile << ' ' << neuron->getOutgoingConnections().at(c)->target->getId();
    }
    connectivityFile << std::endl;

    activityFile << neuron->getTypeId();
    for (int s = 0; s < neuron->getSpikeTimeList().size(); s++) {
      activityFile << ' ' << neuron->getSpikeTimeList().at(s);
    }
    activityFile << std::endl;
  }

  connectivityFile.close();
  activityFile.close();
}
