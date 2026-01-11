# Wigner-Crystal-in-a-Langevin-Heat-Bath-Perturbation-Behavior

A Wigner crystal is a lattice structure that occurs naturally in low density electron gas and was first discovered in 2019. 
Wigner crystals occur naturally in the low electron density quantum dots and modeling their 'thermal melting' behavior is relevant to the understanding of quantum phase transitions in 2D semiconductors, which is the focus of research at the Condensed Matter Group under Dr. Ilya Esterlis of the University of Wisconsin-Madison.

Simulation results at a temperature 0.003C match perturbation behavior observed in literature for the melting of a 26 electron Wigner crystal [See Ref]

Coloumb iteractions are calculated using modular arithmetic to reduce computation time and double-counting errors. 
A modified Verlet algorithm known as 'Leapfrogging' is used with a Gaussian-distributed noise term to incorporate thermal effects of the Langevin heat bath.



[Ref]: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.2667
