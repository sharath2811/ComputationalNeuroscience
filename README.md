# ComputationalNeuroscience
Build a model of a Neural Network using Izekewich Neurons to show to process of learning in the brain.

The overarching question for my model was “What is the most efficient network configuration for a Neural Network?”.

The hypothesis was “The size of a network depends on the strength of the connections between them.”

I used a 1000 neuron Izekevich Network with 800 excitatory and 200 inhibitory neurons. 

I tested this hypothesis using the concept behind Hebbian Learning. My input for the network was a normally distributed random thalamic current. This input represents a stimuli for a particular task. 

The input randomly causes firing of some neurons in the network.  Everytime a neuron fires, all the connections to that neuron increases by a small step. After the simulation is run, a threshold was set for the strength of the connections and all the neurons falling below this threshold was removed from the network. 

The output for the network is the connection strengths of the neurons. The resulting network is a network that represents learning for that particular stimulus.
