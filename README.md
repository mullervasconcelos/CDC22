# Jammer  without channel sensing

**Matlab code for the jammer without channel sensing in the paper "Robust remote estimation over the collision channel in the presence of an intelligent jammer"**

## functions

* **OptimalJammingProbability.m**: Optimal jamming probability for the  jammer without channel sensing beta*. Here, c=1, d=1, and X ~ N(0,1).
* **OptimalJammingProbability_VS_d.m**: Optimal jamming probability for the  jammer without channel sensing beta* as a function of d. Here, c=1, and X ~ N(0,1).
* **OptimalJammingProbability_VS_var.m**: Optimal jamming probabilities beta* as a function of \sigma^2. Here, c=1, d=1 and X ~ N(0,sigma^2).
* **OptimalJammingProbability_VS_cd.m**: Optimal jamming probability for the  jammer without channel sensing beta*as a function of  $c$ and $d$ for **Fig. 2**. Here, X ~ N(0,1). 

# Reactive Jammer

**Matlab code for the reactive jammer in the paper "Robust remote estimation over the collision channel in the presence of an intelligent jammer"**

## functions

* **Main_approximate_FNE.m**: PGA-CCP algorithm for **Table I**.
* **Optimalalphabeta_VS_d.m**: Optimal jamming probabilities alpha* and beta* as a function of $d$. Here, c=1, and X ~ N(0,1).
* **Optimalalphabeta_VS_var**: Optimal jamming probabilities alpha* and beta* as a function of $\sigma^2$ for **Fig. 4**. Here, c=1, d=1 and X ~ N(0,sigma^2).
* **PGA_CCP_Algorithm**: Convergence of PGA-CCP for **Fig. 6**
* **GDA_Algorithm.m**: Convergence of GDA for **Fig. 6**

## auxiliary functions

* **FirstNashEqulibiumChecker.m**: Check whether approximate First Nash Equilibrium is satisfied
* **grad_PGA.m** Gradients for PGA
* **grad_CCP.m**: Gradients for CCP
* **grad_GD.m**: Gradients for GD



