# Bachelor's thesis
### Comparison of Recurrent Neural Networks Trained with Standard and Biologically Plausible Learning Algorithms

This repository contains the code for training and analysis developed for my Bachelor’s degree thesis in Physics, focused on the study of recurrent neural networks (RNNs) as dynamical systems and their training through biologically plausible learning rules.

The thesis addresses a central question at the intersection of computational neuroscience, statistical physics, and machine learning:
to what extent can recurrent neural networks be trained using learning algorithms that respect biological constraints while still producing effective and stable behavior?

Standard training methods for RNNs, such as backpropagation through time, are computationally efficient but rely on mechanisms—global error signals, symmetric connections, and continuous supervision—that are unlikely to be implemented in biological neural circuits. In contrast, this work investigates reward-modulated Hebbian learning rules, which rely only on information locally available at the synapse and on delayed, sparse reward signals.

In particular, this project implements and studies a biologically plausible learning algorithm proposed by Miconi, which combines Hebbian plasticity, stochastic neural perturbations, and delayed reward signals. A fully connected recurrent neural network operating in a near-chaotic regime is trained to perform a delayed nonmatch-to-sample decision-making task, a minimal cognitive task requiring context-dependent memory, temporal integration, and nonlinear mixed selectivity.

The analysis focuses on:
* Learning performance and convergence
* Changes in synaptic weight distributions
* Network dynamics before and after training
* Emergence of nonlinear mixed selectivity in neural activity

This repository includes:

C++ code for numerical implementation of the RNN, the biologically plausible learning algorithm, the delayed nonmatch-to-sample task

Python analysis code for weights, dynamics, and selectivity

Overall, this work supports the idea that biologically realistic learning rules can successfully train recurrent neural networks, reproducing neural dynamics observed in cognitive tasks while maintaining physical and biological plausibility.