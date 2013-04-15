#ifndef NEURON_H
#define NEURON_H

#include <stdexcept>
#include <iostream>
#include <vector>

#include "Connection.h"

class Neuron {
  private:
    const int i;
    const int myMpiProcessRank;
    //std::vector<int> spikeTimeList;
    SpikeList* spikeList;

  public:
    Neuron(const int _i, const int _myMpiProcessRank, SpikeList* _spikeList) : 
      i(_i), myMpiProcessRank(_myMpiProcessRank), spikeList(_spikeList) {}

    int getId() { return i; }

    int getMyMpiProcessRank() { return myMpiProcessRank; }

    virtual std::vector<int>& getSpikeTimeList() { throw std::runtime_error("Neuron::getSpikeTimeList() should never be called");  }
    virtual void addSpikeTime(int t) { throw std::runtime_error("Neuron::addSpikeTime(int) should never be called");  }
    // change to virtual and define in local and remote neurons
    //virtual std::vector<int>& getSpikeTimeList() = 0; //{ return spikeTimeList; }
    // change to virtual and define in local and remote neurons
    //virtual void addSpikeTime(int t) = 0; //{ spikeTimeList.push_back(t); }

    SpikeList* getSpikeList() { return spikeList; }

    int getMpiRank() { return myMpiProcessRank; }

    virtual void scheduleSpikeArrival(int sourceNeuronId, int arrivalTimeStep, int maxSpikeDelay, double efficacy) = 0;
    virtual void receiveSpike(double efficacy) { throw std::runtime_error("Neuron::receiveSpike(double) should never be called"); }
    virtual int getTypeId() = 0;
    virtual void setTypeId(int typeId) {throw std::runtime_error("Neuron::setTypeId(int) should never be called"); }
};

class LocalNeuron : public Neuron {
  private:
    int id;
    ConnectionPointerList connectionList;

    std::vector<int> spikeTimeList;
 
    double V_t;
    double V_spike;
    double V_ext_spike;

  public:
    int refractoryTimeRemaining;

    LocalNeuron(const int i, const int myMpiProcessRank, SpikeList* spikeList) : 
      Neuron(i, myMpiProcessRank, spikeList), V_t(0.0), V_spike(0.0), refractoryTimeRemaining(0) {}
    virtual int getTypeId() = 0;
    virtual double getExternalJ() = 0;
    virtual double getDecayCoefficient() = 0;

    std::vector<int>& getSpikeTimeList() { return spikeTimeList; }
    void addSpikeTime(int t) { spikeTimeList.push_back(t); }

    // local neurons put the spike in the appropriate spike list
    void scheduleSpikeArrival(int sourceNeuronId, int arrivalTimeStep, int maxSpikeDelay, double efficacy) {
      int futureSpikeListIndex = arrivalTimeStep % maxSpikeDelay;
      //std::cout << " local neuron " << getId() << " scheduled a spike from neuron " << sourceNeuronId <<  " with arrival time " << arrivalTimeStep << " (" << futureSpikeListIndex << ") and efficacy " << efficacy << std::endl;
      getSpikeList()->at(futureSpikeListIndex)->push_back(getId());
      getSpikeList()->at(futureSpikeListIndex)->push_back(efficacy);
    }
    void receiveSpike(double efficacy) { addSpike(efficacy); }
    void addSpike(double I) { V_spike += I; }
    double getActivity() { return V_t; }
    void setActivity(double V) { V_t = V; }
    void updateActivity(int externalSpikeCount) {
      V_ext_spike = getExternalJ() * externalSpikeCount;
      V_t *= getDecayCoefficient();
      V_t += V_ext_spike + V_spike;
      V_spike = 0.0;
    }
    void addConnection(Connection* c) { connectionList.push_back(c); }

    std::vector<Connection*>& getOutgoingConnections() { return connectionList; }

    void spikeNow(int t, int maxSpikeDelay) {
      setActivity(50.0);
      addSpikeTime(t);

      if ( ! connectionList.empty() ) {
        for ( int c = 0; c < connectionList.size(); c++ ) {
          Connection* connection = connectionList.at(c);
          //int futureSpikeListIndex = (t + connection->delay) % maxSpikeDelay;
          // here is the original (local) implementation:
          ////spikeList->at(futureSpikeListIndex)->push_back(connection);
          // rather than using the source neuron spike list use the target
          //   neuron spike list
          // rather than connection perhaps store
          //   target neuron global id
          //   efficacy
          // rather than insert spike information directly let the target neuron
          //   decide how to handle it -- local and remote neurons handle things
          //   differently
          connection->target->scheduleSpikeArrival(getId(), t+connection->delay, maxSpikeDelay, connection->efficacy);
        }
      }
    }

};

class ExcitatoryNeuron : public LocalNeuron {
  public:
    static const int typeId = 0;
    static double Jexternal;
    static double decayCoefficient;

    ExcitatoryNeuron(int i, int myMpiProcessRank, SpikeList* spikeList) : 
      LocalNeuron(i, myMpiProcessRank, spikeList) {}
    int getTypeId() { return typeId; }
    double getExternalJ() { return Jexternal; }
    double getDecayCoefficient() { return decayCoefficient; }
};

class InhibitoryNeuron : public LocalNeuron {
  public:
    static const int typeId = 1;
    static double Jexternal;
    static double decayCoefficient;

    InhibitoryNeuron(int i, int myMpiProcessRank, SpikeList* spikeList) : 
      LocalNeuron(i, myMpiProcessRank, spikeList) {}
    int getTypeId() { return typeId; }
    double getExternalJ() { return Jexternal; }
    double getDecayCoefficient() { return decayCoefficient; }
};

class RemoteNeuron : public Neuron {
  private:
    int typeId;
  public:
    RemoteNeuron(int i, int myMpiProcessRank, SpikeList* spikeList) :
      Neuron(i, myMpiProcessRank, spikeList), typeId(-1) {}
    int getTypeId() { return typeId; }
    void setTypeId(int _typeId) { typeId = _typeId; }
    // remote neurons put the spike in one big spike list
    // always use the first list in the spike list
    void scheduleSpikeArrival(int sourceNeuronId, int arrivalTimeStep, int maxSpikeDelay, double efficacy) {
      //std::cout << " remote neuron " << getId() << " scheduled a spike from neuron " << sourceNeuronId <<  " with arrival time " << arrivalTimeStep << " and efficacy " << efficacy << std::endl;
      getSpikeList()->at(0)->push_back(sourceNeuronId);
      getSpikeList()->at(0)->push_back(getId());
      getSpikeList()->at(0)->push_back(efficacy);
      getSpikeList()->at(0)->push_back(arrivalTimeStep);
    }

};

#endif
