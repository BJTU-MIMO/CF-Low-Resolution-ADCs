# CF-Low-Resolution-ADCs

# Cell-Free Massive MIMO With Low-Resolution ADCs Over Spatially Correlated Channels

This is a code package is related to the following scientific article:

Jiayi Zhang, Jing Zhang, and Bo Ai, "[Cell-Free Massive MIMO With Low-Resolution ADCs Over Spatially Correlated Channels](https://ieeexplore.ieee.org/document/9148882)," *ICC 2020 - 2020 IEEE International Conference on Communications (ICC)*.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. *We encourage you to also perform reproducible research!*


## Abstract of Article

Cell-free massive multiple-input multiple-output (MIMO) is a promising technology for future wireless networks. One main challenge of realizing practical cell-free massive MIMO is the high power consumption and huge hardware cost for employing high-resolution analog-to-digital converters (ADCs). To tackle this issue, a promising solution is to use low-resolution ADCs. In this paper, we investigate the cell-free massive MIMO system with low-resolution ADCs over spatially correlated channels. We generalize three levels of receiver cooperation and derive a tight closed-form expression of the spectral efficiency (SE) for a centralized receiver cooperation with large-scale fading decoding as a function of the resolution of ADCs. We also investigate the impact of spatial correlation magnitude on the sum SE. Moreover, we proposed a low-complexity power control method for maximizing the sum SE. Numerical results show that the centralized receiver cooperation needs more quantization bits to achieve the ideal performance and the proposed power control is efficient for improving the system performance.

## Content of Code Package

The package generates the simulation SE results which are used in Figure 1, Figure 2, and Figure 3. To be specific:

- `main`: Main function;
- `Z_functionChannelEstimates`: Perform MMSE channel estimation;
- `Z_Func_LSFD_CorrelatedSMMSE`: Compute the SE by the analytical expression using LSFD;
- `Z_Func_OptLSFD_CorrelatedSMMSE`: Compute the SE by the analytical expression using LSFD with power control optimization;
- `Z_Q_ceshi1_functionComputeSE_AP_uplink`: Compute the SE using low ADC;

See each file for further documentation.


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
