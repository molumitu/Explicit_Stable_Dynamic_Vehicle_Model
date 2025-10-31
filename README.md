# [Vehicle System Dynamics] An Explicit Discrete-Time Dynamic Vehicle Model with Assured Numerical Stability

## Method
Let vehicle dynamical state is $[X, Y, U, V, \varphi, \omega]$, discrete time step is $T_s$, The formulation of our dynamic vehicle model is
$$
\left[\begin{array}{c}
X_{k+1} \\
Y_{k+1} \\
U_{k+1} \\
V_{k+1} \\
\varphi_{k+1} \\
\omega_{k+1}
\end{array}\right]=\left[\begin{array}{c}
X_{k}+T_{\text{s}}\left(U_{k} \cos \varphi_{k}-V_{k} \sin \varphi_{k}\right) \\
Y_{k}+T_{\text{s}}\left(V_{k} \cos \varphi_{k}+U_{k} \sin \varphi_{k}\right) \\
U_{k}+T_{\text{s}} a_{k} \\
\dfrac{m U_{k} V_{k}+T_{\text{s}}\left(l_{\text{f}}k_{\text{f}}-l_{\text{r}} k_{\text{r}}\right) \omega_{k}-T_{\text{s}} k_{\text{f}} \delta_{k} U_{k}-T_{\text{s}} m U_{k}^{2} \omega_{k}}{m U_{k}-T_{\text{s}}\left(k_{\text{f}}+k_{\text{r}}\right)} \\
\varphi_{k}+T_{\text{s}} \omega_{k} \\
\dfrac{I_{\text{z}} U_{k} \omega_{k}+T_{\text{s}}\left(l_{\text{f}} k_{\text{f}}-l_{\text{r}} k_{\text{r}}\right) V_{k}-T_{\text{s}} l_{\text{f}} k_{\text{f}} \delta_{k} U_{k}}{I_{\text{z}} U_{k}-T_{\text{s}}\left(l_{\text{f}}^{2} k_{\text{f}}+l_{\text{r}}^{2} k_{\text{r}}\right)}
\end{array}\right]
$$

## Simulation
We conduct an open-loop test (`open_loop_test.py`) 
- **Control input**: A two-stage step steering.
- **Baseline**: The existing dynamic vehicle model, discretized via the Forward Euler method, serves as the baseline for comparison. 
- **Ground truth**: To provide a reliable reference without relying on complex software like Carsim, we utilized a Simulink model employing the ODE45 solver to generate the ground truth trajectory.
- **Results**: The simulation results with varying discrete time steps (1ms, 50ms, and 100ms) are presented below.

##### Time step = 1 ms
![Time step = 1ms](results/ts_1ms.png)
##### Time step = 50 ms
![Time step = 50ms](results/ts_50ms.png)
##### Time step = 100 ms
![Time step = 100ms](results/ts_100ms.png)

## Bibtex citation
If you find our work useful, please cite our paper:
```text
@article{zhan2024explicit,
  title={An explicit discrete-time dynamic vehicle model with assured numerical stability},
  author={Zhan, Guojian and Ge, Qiang and Gao, Haoyu and Yin, Yuming and Zhao, Bin and Eben Li, Shengbo},
  journal={Vehicle System Dynamics},
  pages={1--24},
  year={2024},
  publisher={Taylor \& Francis}
}
```