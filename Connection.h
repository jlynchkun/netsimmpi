#ifndef CONNECTION_H
#define CONNECTION_H

#include <vector>

class Neuron;

class Connection {
  public:
    Neuron* source;
    Neuron* target;
    const double efficacy;
    const int delay;

    Connection(Neuron* s, Neuron* t, double e, int d): source(s), target(t), efficacy(e), delay(d) {}
};

typedef std::vector<Connection*> ConnectionPointerList;
//typedef std::vector< ConnectionPointerList* > SpikeList;
typedef std::vector< std::vector<double>* > SpikeList;

#endif
