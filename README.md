# parameter_estimation_hidden_gaussian

These files represent the numerical analysis of the two examples in the paper paper "On parameter estimation of the hidden Gaussian process in perturbed SDE.".
The paper is available via https://arxiv.org/pdf/1904.09750.pdf or https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-15/issue-1/On-parameter-estimation-of-the-hidden-Gaussian-process-in-perturbed/10.1214/20-EJS1788.full

We have implemented both examples in the programming language R. 
Additionally, we have implemented the Two-step MLE-process which is discribed in the discussion section (compare also the paper "On the multi-step MLE-process for ergodic diffusion.")
and the well-known MLE approach for parameter estimation (for Example 2 only; in Example 1, this construction works analogously).

The implementation is far from being perfect or optimal in some sense. For instance, one can substitue all for-loops by one over each file. 
However, these files are working in a reasonable time (except for the Two-step MLE-process, which has heavy computational disadvantages by construction).

Questions, remarks or suggestions for improvement are highly appreciated! Please send me an e-mail via jonasewers@yahoo.de or jonas.ewers@uni-muenster.de
