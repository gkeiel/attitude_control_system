import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from acs_config import *


# =====================================================
#  Orbit
# =====================================================
class Orbit:
    def __init__(self, r, v, mu):
        self.r  = r.copy()
        self.v  = v.copy()
        self.mu = mu
        
    def step(self, dt):
        """
        simple orbit propagation:
            r, v: position and velocity
            mu:   gravitational parameter
            dt:   time step
        """
        a      = -self.mu*self.r/np.linalg.norm(self.r)**3  # acceleration
        self.v = self.v +a*dt                               # numerical integration v
        self.r = self.r +self.v*dt                          # numerical integration r
        return self.r, self.v
    
    @staticmethod
    def reference_vectors(r, v):
        """
        generate reference vectors:
            nadir
            along
            cross
        """
        nadir = -r/np.linalg.norm(r)
        along =  v/np.linalg.norm(v)
        cross = np.cross(r, v)/np.linalg.norm(np.cross(r, v))
        return nadir, along, cross
       

# =====================================================
#  AttitudeDynamics
# =====================================================
class AttitudeDynamics:
    def __init__(self, J):
        self.J   = J
        self.w   = np.zeros(3)
        self.DCM = np.eye(3)
        
    @staticmethod
    def skew(w):
        w_x, w_y, w_z = w      # angular velocity vector
        return np.array([
            [ 0,  -w_z,  w_y],
            [ w_z,   0, -w_x],
            [-w_y, w_x,  0]
        ])
        
    @staticmethod
    def gravity_tau(r, DCM, J, mu):
        r_hat_eci  = r/np.linalg.norm(r)
        r_hat_body = DCM@r_hat_eci
        return 3*mu/np.linalg.norm(r)**3*np.cross(r_hat_body, J@r_hat_body)

    @staticmethod
    def dcm_to_euler(DCM):
        theta = -np.arcsin(DCM[0,2])                                         # pitch
        phi   =  np.arctan2(DCM[1,2]/np.cos(theta), DCM[2,2]/np.cos(theta))  # roll
        psi   =  np.arctan2(DCM[0,1]/np.cos(theta), DCM[0,0]/np.cos(theta))  # yaw
        return phi, theta, psi
    
    @staticmethod
    def euler_history(DCM_data):
        roll_data, pitch_data, yaw_data = [], [], []
        for DCM in DCM_data:
            phi, theta, psi = AttitudeDynamics.dcm_to_euler(DCM)
            roll_data.append(phi)
            pitch_data.append(theta)
            yaw_data.append(psi)
        return np.array(roll_data), np.array(pitch_data), np.array(yaw_data)

    def omega_dot(self, torque):
        return np.linalg.inv(self.J) @ (torque -np.cross(self.w, self.J@self.w))

    def step(self, torque, dt):
        def f(w, DCM):
            w_dot   = np.linalg.inv(self.J) @ (torque - np.cross(w, self.J@w))
            DCM_dot = DCM@self.skew(w)
            return w_dot, DCM_dot

        w, DCM = self.w, self.DCM
        
        k1_w, k1_D = f(w, DCM)
        k2_w, k2_D = f(w +0.5*dt*k1_w, DCM +0.5*dt*k1_D)
        k3_w, k3_D = f(w +0.5*dt*k2_w, DCM +0.5*dt*k2_D)
        k4_w, k4_D = f(w +dt*k3_w,     DCM +dt*k3_D)
        
        self.w     = w   +(k1_w +2*k2_w +2*k3_w +k4_w)*dt/6
        self.DCM   = DCM +(k1_D +2*k2_D +2*k3_D +k4_D)*dt/6
        
        # re-orthonormalization
        u, _, vt = np.linalg.svd(self.DCM)
        self.DCM = u@vt

        return self.w, self.DCM


# =====================================================
#  AttitudeController
# =====================================================
class AttitudeController:
    def __init__(self, k_p, k_i, k_d, p=5):
        self.k_p  = k_p
        self.k_i  = k_i
        self.k_d  = k_d
        self.p    = p
        self.e_1  = np.zeros(3)
        self.u_i1 = np.zeros(3)
        self.u_d1 = np.zeros(3)
    
    
    def compute(self, r, v, DCM, T):
        """
        attitude controller
        arguments:
            r, v : inertial position and velocity
            DCM  : attitude
            T    : sample time
        return:
            DCM_ref: DCM reference
            u_sat:   saturated control signal
        """
        # DCM reference
        z_ref   = -r/np.linalg.norm(r)          # nadir
        x_ref   =  v/np.linalg.norm(v)          # along-track
        y_ref   = np.cross(z_ref, x_ref)
        y_ref  /= np.linalg.norm(y_ref)
        x_ref   = np.cross(y_ref, z_ref)
        DCM_ref = np.column_stack((x_ref, y_ref, z_ref))

        # tracking error
        D = DCM_ref.T@DCM
        e = 0.5*np.array([
            D[2,1] -D[1,2],
            D[0,2] -D[2,0],
            D[1,0] -D[0,1]
        ])
                
        # control law
        u_p   = self.k_p*e
        u_i   = self.u_i1 +self.k_i*T*(e +self.e_1)/2
        a     = (2 -self.p*T)/(2 +self.p*T)
        b     = 2*self.p*self.k_d/(2 +self.p*T)
        u_d   = a*self.u_d1 +b*(e -self.e_1)
        u     = -u_p -u_i -u_d

        # real control signal
        u_sat = np.clip(u, -100, 100)
        if np.any(np.abs(u) >= 100): u_i = self.u_i1
        
        # update past values
        self.e_1  = e
        self.u_i1 = u_i
        self.u_d1 = u_d
        return DCM_ref, e, u_sat


# =====================================================
#  Animator
# =====================================================
class Animator:
    def __init__(self, r_data, DCM_data, DCM_ref_data, e_data, skip, lim=8000):
        self.r_data       = r_data[::skip]
        self.DCM_data     = DCM_data[::skip]
        self.DCM_ref_data = DCM_ref_data[::skip]
        self.e_data       = e_data[::skip]
        
        self.fig      = plt.figure()
        self.ax       = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlim([-lim, lim])
        self.ax.set_ylim([-lim, lim])
        self.ax.set_zlim([-lim, lim])
        
        self.ax.set_title("Orbit and Attitude Animation")
        self.ax.set_xlabel("x [km]")
        self.ax.set_ylabel("y [km]")
        self.ax.set_zlabel("z [km]")
        self.sat_point, = self.ax.plot([], [], [], 'ko')
        self.orbit_line, = self.ax.plot([], [], [], 'k-', linewidth=1, alpha=0.7)
        self.orbit_line.set_alpha(0.5)
        self.orbit_line.set_color("red")

        # Earth
        u = np.linspace(0, 2*np.pi, 60)
        v = np.linspace(0, np.pi, 30)
        x_earth = r_earth*np.outer(np.cos(u), np.sin(v))
        y_earth = r_earth*np.outer(np.sin(u), np.sin(v))
        z_earth = r_earth*np.outer(np.ones_like(u), np.cos(v))
        self.ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.5, linewidth=0)
    
    def update(self, i):
        for c in self.ax.collections[1:]: c.remove()
    
        r   = self.r_data[i]/1e3                          # [km]
        DCM = self.DCM_data[i]
        
        # sattelite
        self.sat_point.set_data([r[0]], [r[1]])
        self.sat_point.set_3d_properties([r[2]])

        # orbit up to current time
        r_hist = np.array(self.r_data[:i+1]) / 1e3
        self.orbit_line.set_data(r_hist[:,0], r_hist[:,1])
        self.orbit_line.set_3d_properties(r_hist[:,2])
    
        scale = 1000
        self.ax.quiver(*r, *(DCM[:,0]*scale), color='r')  # x_body
        self.ax.quiver(*r, *(DCM[:,1]*scale), color='g')  # y_body
        self.ax.quiver(*r, *(DCM[:,2]*scale), color='b')  # z_body
        return self.sat_point,
    
    def animate(self, interval=20):
        self.ani = FuncAnimation(self.fig, self.update, frames=len(self.r_data), interval=interval)
        self.ani.save("data/animation.gif", writer="pillow", fps=20)
        plt.show()
    
    def plotting(self, dt):        
        # attitude error
        e_data = np.array(self.e_data)
        t      = np.arange(len(e_data))*dt
        plt.figure(figsize=(10,6))
        plt.plot(t, e_data[:,0], label=r"$e_x$")
        plt.plot(t, e_data[:,1], label=r"$e_y$")
        plt.plot(t, e_data[:,2], label=r"$e_z$")
        plt.xlabel("Time [s]")
        plt.ylabel("Attitude error")
        plt.title("Attitude error vector")
        plt.legend()
        plt.grid(True)
        plt.savefig("data/attitude_errors.png", dpi=300, bbox_inches="tight")
        plt.close()
        
        # Euler angles
        phi, theta, psi = [], [], []
        phi_r, theta_r, psi_r = [], [], []

        for D, Dref in zip(self.DCM_data, self.DCM_ref_data):
            p, t, y    = AttitudeDynamics.dcm_to_euler(D)
            pr, tr, yr = AttitudeDynamics.dcm_to_euler(Dref)

            phi.append(p);   theta.append(t);   psi.append(y)
            phi_r.append(pr);theta_r.append(tr);psi_r.append(yr)
            t = np.arange(len(phi)) * dt

        phi     = np.array(phi)
        theta  = np.array(theta)
        psi    = np.array(psi)
        phi_r  = np.array(phi_r)
        theta_r= np.array(theta_r)
        psi_r  = np.array(psi_r)

        t = np.arange(len(phi)) * dt
        plt.figure(figsize=(10,8))
        plt.subplot(3,1,1)
        plt.plot(t, np.degrees(phi), label="Output")
        plt.plot(t, np.degrees(phi_r), "--", label="Reference")
        plt.ylabel("[deg]")
        plt.title("Roll")
        plt.grid(True)
        plt.legend()

        plt.subplot(3,1,2)
        plt.plot(t, np.degrees(theta), label="Output")
        plt.plot(t, np.degrees(theta_r), "--", label="Reference")
        plt.ylabel("[deg]")
        plt.title("Pitch")
        plt.grid(True)
        plt.legend()

        plt.subplot(3,1,3)
        plt.plot(t, np.degrees(psi), label="Output")
        plt.plot(t, np.degrees(psi_r), "--", label="Reference")
        plt.ylabel("[deg]")
        plt.title("Yaw")
        plt.xlabel("Time [s]")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig("data/references.png", dpi=300, bbox_inches="tight")
        plt.close()
