# build the netsimmpi program
mpic++ -I/share/apps/boost/include netsim.cpp NeuralNetwork.cpp Neuron.cpp -L/share/apps/boost/lib -lboost_mpi -lboost_serialization -o netsimmpi.out
