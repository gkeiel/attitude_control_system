import os
from acs_classes import Orbit, AttitudeDynamics, AttitudeController, Animator
from acs_config import *
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# instance objects and give arguments from acs_config.py
orbit    = Orbit(r0, v0, mu_earth)
att_dyn  = AttitudeDynamics(J)
att_ctrl = AttitudeController(k_p, k_i, k_d, p=5)

# initialize data
r_data, DCM_ref_data, DCM_data, e_data = [], [], [], []

for k in range(N):
    
    # orbit evolution
    r, v = orbit.step(dt)
    
    # attitude control signal    
    DCM_ref, e, tau_u = att_ctrl.compute(r, v, att_dyn.DCM, dt)
    
    # disturbance
    tau_g = AttitudeDynamics.gravity_tau(r, att_dyn.DCM, J, mu_earth)
    
    # real tau
    tau   = tau_u +tau_g
    
    # attitude dynamics
    w, DCM = att_dyn.step(tau, dt)
       
    # store data    
    r_data.append(r.copy())
    DCM_data.append(DCM.copy())
    DCM_ref_data.append(DCM_ref.copy())
    e_data.append(e.copy())
    
# 3D animation
animation = Animator(r_data, DCM_data, DCM_ref_data, e_data, skip=5)
animation.animate()
animation.plotting(dt)