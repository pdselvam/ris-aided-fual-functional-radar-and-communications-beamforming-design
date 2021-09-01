# RIS-aided Dual-functional Radar and Communication Beamforming Design

This is the project submitted in fulfillment of requirements for the degree of MSc Communications and Signal Processing at Imperial College London.

## Content

### Techniques

- Communication and radar spectrum sharing
- Reconfigurable intelligent surface
- MIMO radar system
- MISO multi-user wireless communication system
- Weighted Minimum Mean Square Error (WMMSE) optimization framework
- Fractional programming
- Semidefinite Relaxation (SDR)

### Contributions

- We propose a joint beamforming design at BS and RIS that maximizes WSR and the probing power at target for both separated and shared deployment. To our knowledge, this is the first study of WSR maximising for a RIS-aided DFRC system. As in [[1](https://ieeexplore.ieee.org/abstract/document/9200993)], we also consider probing power as radar metric to make a clear tradeoff comparison. However, the logarithm form and additional radar signal term of WSR, the quadratic form of probing power, and the quadratic equality power constraint make the optimization of joint beamforming rather non-convex and elusive.
- This is also the first work that investigates the novel group or fully connected RIS model in RIS-aided DFRC. This novel RIS model is capable of enhancing the SNR performance especially in the Rayleigh fading channel [[2](https://ieeexplore.ieee.org/abstract/document/9514409)], but its potential benefits for enlarging achievable region of WSR and probing power in DFRC are studied in this project.
- We propose an AO algorithm for the proposed non-convex design in separated deployment. As was proposed in [[1](https://ieeexplore.ieee.org/abstract/document/9200993)], the active beamforming problem is reformulated to a semidefinite programming (SDP) problem using WMMSE framework. The passive beamforming problem is reformulated to an unconstrained problem based on Lagrangian dual transform [[3](https://ieeexplore.ieee.org/abstract/document/8310563)], quadratic transform [[4](https://ieeexplore.ieee.org/abstract/document/8314727)], and scattering-reactance relationship [[2](https://ieeexplore.ieee.org/abstract/document/9514409)]. Compared with the passive beamforming optimization method in [[5](https://ieeexplore.ieee.org/abstract/document/8982186)], this method is an extended version that considers an additional radar signal term and generalized RIS model.
- We propose another AO algorithm to solve the non-convex design for shared deployment. The optimal active beamforming can be obtained by WMMSE framework and semidefinite relaxation (SDR) [[6](https://ieeexplore.ieee.org/abstract/document/5447068)]. In contrast to the majorization-minimization (MM) method used in [[1](https://ieeexplore.ieee.org/abstract/document/9200993)], which requires several steps to reach the optimal values, this SDR method only needs one step and therefore has lower complexity. The passive beamforming is optimized using the same method as separated deployment.

## Running the simulations

### Prerequisites

- [MATLAB](https://uk.mathworks.com/products/matlab.html)
- [CVX](http://cvxr.com/cvx/)

### Launch

##### Separated Deployment

```
>> cd simulation_separated\
```

For the `tradoff and beampattern in Rayleigh channel` plot,  run

```
>> comparison_rician_0
```

For the `tradoff and beampattern in LOS channel` plot, run

```
>> comparison_rician_1000
```

For the `effect of the number of reflecting elements ` plot, run

```
>> comparison_elements
```

To evaluate the `convergence of algorithm`, run

```
>> ris_aided_convergence
```



##### Shared Deployment

```
cd simulation_shared\
```

For the `comparison of eigenvalue decomposition and Gaussian randomization` plot, run

```
>> eigenvalue_gaussian_comparison
```

 Others are the same as those in separated deployment.

