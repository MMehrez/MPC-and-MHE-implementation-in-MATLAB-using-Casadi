## Description
This is a Python implementation for model predictive control (MPC) ported from the MATLAB code provided.


The implementation uses the Casadi Python Package for numerical optimization and {numpy + matplotlib} for visualization. 
A holonomic mobile robot (Mecanum wheeled omnnidirectional mobile robot (MWOMR)) is used as a system for the implementation. 
A helpful paper explaining the transfer function: [link from researchgate.net](https://www.researchgate.net/publication/334319114_Model_Predictive_Control_for_a_Mecanum-wheeled_robot_in_Dynamical_Environments)


Casadi can be downloaded here https://web.casadi.org/


## Content
1. `mpc_code.py` &rightarrow; The main Python script for MPC
2. `simulation_code.py` &rightarrow; A helper file implementing the visualization used in MPC code

## Requirements
Before running the codes, make sure that you have Python3.5+ installed on your computer alongside the following packages:
- [CasADi](https://web.casadi.org/get/)
- [NumPy](https://numpy.org/install/)
- [Matplotlib](https://matplotlib.org/users/installing.html)

## Contributors
- Dr. Mohamed Mehrez (Origianl author/owner) </br>
  m.mehrez.said@gmail.com  </br>
  [github.com/MMehrez](github.com/MMehrez) </br>
- Mohammad Abdussalam </br>
  mohammadahmad01722@gmail.com </br>
  [github.com/Mohammad1722](github.com/Mohammad1722)
- Yassmin Hesham </br>
  yasmin98hesham@gmail.com </br>
  [github.com/Yasmin-Hesham](github.com/Yasmin-Hesham)
