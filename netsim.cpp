#include <iostream>
#include <string>
#include <sstream>

#include <functional>

#include <boost/mpi.hpp>

#include "NeuralNetwork.h"

// do not load the intel compiler on compute nodes
// load module gcc-4.5.2 and module boost
// mpic++ -I/share/apps/boost/include netsim.cpp NeuralNetwork.cpp Neuron.cpp -L/share/apps/boost/lib -lboost_mpi -lboost_serialization
// mpirun -np 2 a.out

int main(int argc, char* argv[]) {
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  std::cout << "argc " << argc << std::endl;
  std::cout << "argv[0] " << argv[0] << std::endl;
  std::cout << "argv[1] " << argv[1] << std::endl; // simulation time in seconds
  std::cout << "argv[2] " << argv[2] << std::endl; // seed multiplier
  std::cout << "argv[3] " << argv[3] << std::endl; // output directory

  std::stringstream secondsString(argv[1]);
  int seconds;
  secondsString >> seconds;
  std::cout << "simulation time " << seconds << " (s)" << std::endl;

  std::stringstream seedString(argv[2]);
  int seed;
  seedString >> seed;
  seed = seed + world.rank();
  std::cout << "seed " << seed << std::endl;

  char* outputDirPath = argv[3];
  std::cout << "output directory path "  << outputDirPath << std::endl;

  std::cout << "process count : " << world.size() << std::endl;
  std::cout << "process rank  : " << world.rank() << std::endl;  

  // do several 2 second simulations with different J
  float Jmultiplier = 1.75;  // making this double causes problems with mpi ???
  //for ( Jmultiplier = 1.7; Jmultiplier < 2.0; Jmultiplier += 0.1 ) {
    std::stringstream activityFileName;
    activityFileName << outputDirPath << "/nn-activity-" << seconds << "s-process-" << world.rank() << "of" << world.size() << ".txt";
    std::stringstream connectivityFileName;
    connectivityFileName << outputDirPath << "/nn-connectivity-" << seconds << "s-process-" << world.rank() << "of" << world.size() << ".txt";

    //double m = 1.80; // 1.72 for one process works !!!! this broke the program
    //double Jee = 0.21;
    //double Jie = Jmultiplier * 0.21;
    //double Jei = 3.0 * 0.21;
    //double Jii = 3.0 * Jmultiplier * 0.21;
    //int seedMultiplier = 13;
    NeuralNetwork nn = NeuralNetwork(7500, seed, 0.21, Jmultiplier * 0.21, 3.0 * 0.21, 3.0 * Jmultiplier * 0.21);

    nn.printNetworkParameters(std::cout);

    // testing -- 500 steps to get the first spikes
    // testing -- 1000 steps to look for divergence between 1 and 2 processes
    //nn.simulate(100*10);
    // production:
    nn.simulate(100*1000 * seconds);

    nn.writeOutputFiles(connectivityFileName.str(), activityFileName.str());
  //}

  int rootProcessRank = 0;
  if ( world.rank() == rootProcessRank ) {
    int spikingNeuronCountAllProcesses;
    boost::mpi::reduce(world, nn.getSpikingNeuronCount(), spikingNeuronCountAllProcesses, std::plus<int>(), rootProcessRank);
    std::cout << "all processes total spiking neuron count: " << spikingNeuronCountAllProcesses << std::endl;
  }
  else {
    boost::mpi::reduce(world, nn.getSpikingNeuronCount(), std::plus<int>(), rootProcessRank);
  }
  return 0;
}
