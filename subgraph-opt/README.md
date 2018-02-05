# subgraph-opt
fast optimization for subgraph problems

aistats-code/ contains Matlab code for experiments in AISTATS paper. *auc*.m are main scripts, starOpt.m contains the main subgraph algorithm.

TODO:
  * Optimize starOpt.m to work with large scale problems in terms of memory and run-time
    - Only keep low rank sketch for Y using expv
    - Optimize output of paramQ_star to use sparse matrices, and figure out C similarly to obtain an optimized QX
    - Give function handle to eigs that does multiplication by input matrix
  * Run new experiments on larger simulated datasets, e.g. geometric graphs with shape anomaly to evaluate scaling better
    - Look at real graph datasets from Leskovec's SNAP page
  * Clarify proof of Theorem 4.1 and provide sketch for statements after the theorem, regarding nearly-linear iteration time
  * Maybe extend to anchor-less case by using a different constraint. We currently have
      L_star = L_{K_S} + u u^T
    for some vector u that denotes current input to root/outputs to other nodes. Instead at each iteration select a u randomly, then round to solution at the end by projecting (as it may include multi clusters).
