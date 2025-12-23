# Attitude control systems

A repository of algorithms exploring the use of **control systems applied to attitude control** in spacecrafts.

The following algorithms are available:
- [Attitude control PID](./attitude_control_pid)


## Attitude control PID 

An **object-oriented programming (OOP) code for attitude control of an spacecraft** in low orbit arround a celestial body. It offers a model for orbital and attitude dynamics, attitude control and simulation.

### ⚙️ Features

- **Orbital dynamics**: simplified orbit propagation (2-body problem) and reference vectors generation.
- **Attitude dynamics**: rotational dynamics with Runge-Kutta 4th order (RK4) numerical integration.
- **Attitude control**: digital proportional-integral-derivative (PID) controller.
- **Animation and analysis**: 3D animation of orbital and attitude moviment. 
- **Software structure**: object-oriented architeture with independent classes and configuration file.

### 🎞️ Animation

- **Orbit and attitude animation**:
![Orbit and Attitude](animation_example.gif)

- **Attitude errors**:
  <p align="center">
  <img width="733" height="395" alt="attitude_errors" src="https://github.com/user-attachments/assets/ae198956-3ff5-4c10-b2cb-45d6ffe91146" />
  </p>

- **References and output angles**:
  <p align="center">
  <img width="733" height="395" alt="references" src="https://github.com/user-attachments/assets/33efd5c7-ea06-49cf-83fd-3c6d09b54e39" />
  </p>

This algorithm offers a simplified and modular version of the problem, heavily inspired by **Spacecraft-Attitude-Control-System** in

https://github.com/brunopinto900/Spacecraft-Attitude-Control-System
