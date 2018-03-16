# Performance-Estimation-Problems-For-Newton

This code can be used to reproduce the results from the work (on [arXiv](https://arxiv.org/abs/1709.05191)):

> [1] de Klerk, Etienne, Francois Glineur, and Adrien Taylor. "Worst-case convergence analysis of gradient and Newton methods through semidefinite programming performance estimation." arXiv:1709.05191 (2017).

## Getting started

To use the code, download the repository and execute the files on a one-by-one basis. The code makes use of the Symbolic Computation Toolbox of Matlab.


## List of files

- [`ExactLS_distance_validations`](ExactLS_distance_validations.m) Symbolic verification of the worst-case bound (in terms of distance to the optimal solution) obtained when using exact linesearch.
- [`ExactLS_gradientnorm_validations`](ExactLS_gradientnorm_validations.m) Symbolic verification of the worst-case bound (in terms of gradient norm) obtained when using exact linesearch.

- [`Fixedstep_distance_validations_small`](Fixedstep_distance_validations_small.m) Symbolic verification of the worst-case bound (in terms of distance to the optimal solution) obtained when using fixed stepsizes.
- [`Fixedstep_gradientnorm_validations_small`](Fixedstep_gradientnorm_validations_small.m) Symbolic verification of the worst-case bound (in terms of distance to the optimal solution) obtained when using fixed stepsizes.
- [`Fixedstep_funcvalues_validations`](Fixedstep_funcvalues_validations.m) Symbolic verification of the worst-case bound (in terms of objective function accuracy) obtained when using fixed stepsizes.

## Authors
- [**Etienne de Klerk**](https://sites.google.com/site/homepageetiennedeklerk/)
- [**François Glineur**](https://perso.uclouvain.be/francois.glineur/)
- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
