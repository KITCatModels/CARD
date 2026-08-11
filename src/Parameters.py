import numpy as np


class Parameters:
    # constructor
    def __init__(self):
        self.L = 1                              # [m]
        self.D = 20e-3                          # [m]
        self.d_cat = 3e-3                       # [m]
        self.R = 8.314                          # [J/Kmol]
        self.Mw_i = np.array([0.002, 0.001])    # [kg/mol]    # 0: A, 1: B 
        self.k_0 = 10                           # [m^3/(s kg_cat)]
        self.E_A = 20e3                         # [J/mol]
        self.delta_H = -700                     # [J/mol]
        self.nu = np.array([-1, 2])             # [-]
        self.eps = 0.5                          # [-]
        self.U = 50                             # [W/m^2K]
        self.T_wall = 320                       # [K]
        self.cp_f = 1000                        # [J/kgK]
        self.cp_cat = 700                       # [J/kg/K]
        self.rho_cat = 800                      # [kg/m^3]
        self.D_A = 5e-3                         # [m^2/s]
        self.mu = 1e-5                          # [Pa s]

        # spatial discretization (axial)
        self.n_axial = 200

    # method
    def f_delta_z(self):
        delta_z = 1 / self.n_axial
        return delta_z

    def f_z_c(self):
        z_c = np.linspace(
            0 + self.f_delta_z() / 2, 1 - self.f_delta_z() / 2, self.n_axial
        )
        return z_c
