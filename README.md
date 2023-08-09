# Transformed Model Reduction Experiments
  Experiments in `Transformed model reduction for partial differential equations with sharp inner layers` submitted to SISC 2023.
- [Transformed Model Reduction Experiments](#transformed-model-reduction-experiments)
  - [Projects Description](#projects-description)
    - [File description](#file-description)
  - [Usage Instructions](#usage-instructions)
    - [ac2d](#ac2d)
    - [burgers1d](#burgers1d)
  - [Authors and Contact](#authors-and-contact)
  - [Copyright and License](#copyright-and-license)

## Projects Description  
  `ac2d` contains experiments code for 2D Allen-Cahn equation, and `burgers1d` contains code for 1D Burgers equation. 
### File description
- `./ac2d`
  - `Fp_bad.m`: a bad DEIM approximation function.
  - `Fp_bd1.m`: speed-up DEIM approximation of nonlinear-function in slow variable equation with the boundary condition $\frac{\partial^2 v}{\partial \mathbf{n}^2}=0$. 
  - `Fp_bd2.m`: speed-up DEIM approximation of nonlinear-function in slow variable equation with the boundary condition $\frac{\partial v}{\partial \mathbf{n}}=1$.
  - `initial.m`: set up spatial region, initial function, parameters like $\varepsilon$, etc.
  - `main.m`: main program of 2D Allen-Cahn equation.
  - `MakeFigures.m`: show the results.
  - `myOffline.m`: offline procedure in POD-qDEIM model reduction.
  - `myOnline.m`: online procedure in POD-qDEIM model reduction.
  - `POD.m`: POD function.
  - `SetuFOM.m`: set the full order model for origional variable $u$.
  - `SetuROM.m`: using the reduce basis from `myOffline` to set the reduced order model for origional variable $u$.
  - `SetvFOM.m`: set the full order model for slow variable $v$.
  - `SetvROM.m`: using the reduce basis from `myOffline` to set the reduced order model for slow variable $v$.
  - `TTY_RK.m`: a function of explicit RK4. 

- `./burgers1d`
  - `DEIM.m`: q-DEIM procedure.
  - `DFWCNS.m`: using WCNS(weighted compact nonlinear scheme) to generate the nonlinear term in Burgers equation.
  - `initialElements.m`: initial all the parameters, initial functions, etc. for both origional variable $u$ and slow variable $v$.
  - `main.m`: main program of 1D Burgers equation.
  - `makeFigures.m`: show the results.
  - `offline.m`: offline procedure in POD-qDEIM model reduction.
  - `online.m`: online procedure in POD-qDEIM model reduction.
  - `POD.m`: POD function.
  - `v2u.m`: a function receive transform equation and slow variable $v$ (approximate) solution then generate ($U^v_{appr}$) $U^v$.
## Usage Instructions  

### ac2d
  User can modify the initial function and spatial region and other parameters like $\varepsilon$ in `initial.m` file. This project contains 2 boundary conditions which can be switched by changing `bdCase` in `initial.m`. The POD and DEIM base number can be modified in `main.m` file. Run the `main.m` file to get the (approximate) solution of $u$ and $v$.

  If user wants to use other boundary conditions, the DEIM approximate nonlinear function should be rewritten, or using the rough DEIM approximation `Fp_bad.m`. If the other discretization schemes are using, `SetuFOM.m` and `SetvFOM.m` should be rewritten. 
  
### burgers1d
  The initial functions and parameters like $\varepsilon$, POD and DEIM base number are modified in `main.m`, run the `main.m` file to get the (approximate) solution of $u$ and $v$. If other schemes and new boundary conditions are using, the `initialElements.m` file need to be modified. 
  
## Authors and Contact
  - Tianyou Tang, LSEC, ICMSEC, NCMIS,  Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China. Contact: tangtianyou@ lsec.cc.ac.cn. 
  - Xianmin Xu. LSEC, ICMSEC, NCMIS,  Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China. Corresponding author. Contact: xmxu@ lsec.cc.ac.cn. 
## Copyright and License
This project is distributed with the MIT license which translates roughly that you can use it however you want and for whatever reason you want. All the information regarding support, copyright and the license can be found in the LICENSE file.
