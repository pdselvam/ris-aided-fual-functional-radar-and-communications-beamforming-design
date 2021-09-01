# RIS-aided Dual-functional Radar and Communication Beamforming Design

This is the project submitted in fulfillment of requirements for the degree of MSc Communications and Signal Processing at Imperial College London.

## Content

### Techniques

- Communication and radar spectrum sharing
- Reconfigurable intelligent surface
- MIMO radar system
- MISO multi-user wireless communication system
- Weighted Minimum Mean Square Error (WMMSE) optimization framework [[1](https://ieeexplore.ieee.org/abstract/document/4712693)]
- Fractional programming [[2](https://ieeexplore.ieee.org/abstract/document/8314727), [3](https://ieeexplore.ieee.org/abstract/document/8310563)]
- Semidefinite Relaxation (SDR) [[4](https://ieeexplore.ieee.org/abstract/document/5447068)]

### Contributions

-  To our knowledge, this is the first study of WSR maximising for a RIS-aided DFRC system. 
- This is also the first work that investigates the novel group or fully connected RIS model [[5](https://ieeexplore.ieee.org/abstract/document/9514409)] in RIS-aided DFRC.
- Two alternating optimization algorithms for the proposed non-convex beamforming design.

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

