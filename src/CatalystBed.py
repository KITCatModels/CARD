import casadi as casADi

from .Properties import Properties


class CatalystBed(Properties):
    def __init__(self, w_i_in, T_in, u_in, p_in):
        super().__init__()
        # boundary conditions (inlet)
        self.w_i_in = w_i_in  # [-]   # A, B
        self.T_in = T_in  # [K]
        self.u_in = u_in  # [m/s]
        self.p_in = p_in  # [Pa]

    def create_dae(self):
        # define differential variables
        w_i = casADi.SX.sym("w_i", (self.n_axial, 2))  # [mol/m^3]      # A, B
        T = casADi.SX.sym("T", self.n_axial)  # [K]

        # define algebraic variables
        u = casADi.SX.sym("u", self.n_axial)  # [m/s]
        p = casADi.SX.sym("p", self.n_axial)  # [Pa]

        # define mass conservation as algebraic equation
        ae_conti = casADi.SX.sym("ae_conti", self.n_axial)
        for z in range(self.n_axial):
            ae_conti[z] = u[z] - self.u_in * self.f_rho(
                self.w_i_in, self.T_in, self.p_in
            ) / self.f_rho(w_i[z, :].T, T[z], p[z])

        # define pressure drop equation as algebraic equation
        ae_p = casADi.SX.sym("ae_p", self.n_axial)
        for z in range(1, self.n_axial):
            ae_p[z] = (
                1 / (self.L * self.f_delta_z()) * (p[z] - p[z - 1])
                + (1 - self.eps) ** 2
                / self.eps**3
                * (150 * self.mu)
                / self.d_cat**2
                * u[z]
                + (1 - self.eps)
                / self.eps**3
                * 1.75
                / self.d_cat
                * self.f_rho(w_i[z, :].T, T[z], p[z])
                * u[z] ** 2
            )

        # boundary condition (inlet, z=0)
        z = 0
        ae_p[z] = (
            1 / (self.L * self.f_delta_z()) * (p[z] - self.p_in)
            + (1 - self.eps) ** 2 / self.eps**3 * (150 * self.mu) / self.d_cat**2 * u[z]
            + (1 - self.eps)
            / self.eps**3
            * 1.75
            / self.d_cat
            * self.f_rho(w_i[z, :].T, T[z], p[z])
            * u[z] ** 2
        )

        # define species conservation as ode
        ode_wi = casADi.SX.sym("ode_wi", (self.n_axial, 2))
        for comp in range(len(self.nu)):
            for z in range(1, self.n_axial):
                ode_wi[z, comp] = -(w_i[z, comp] - w_i[z - 1, comp]) * u[z] / (
                    self.L * self.f_delta_z() * self.eps
                ) + self.Mw_i[comp] / (
                    self.f_rho(w_i[z, :].T, T[z], p[z]) * self.eps
                ) * (
                    1 - self.eps
                ) * self.rho_cat * self.nu[
                    comp
                ] * self.f_eta(
                    w_i[z, :].T, T[z], p[z]
                ) * self.f_r(
                    w_i[z, :].T, T[z], p[z]
                )

        # boundary condition (inlet, z=0)
        z = 0
        for comp in range(len(self.nu)):
            ode_wi[z, comp] = -(w_i[z, comp] - self.w_i_in[comp]) * u[z] / (
                self.L * self.f_delta_z() * self.eps
            ) + self.Mw_i[comp] / (self.f_rho(w_i[z, :].T, T[z], p[z]) * self.eps) * (
                1 - self.eps
            ) * self.rho_cat * self.nu[
                comp
            ] * self.f_eta(
                w_i[z, :].T, T[z], p[z]
            ) * self.f_r(
                w_i[z, :].T, T[z], p[z]
            )

        # define energy conservation as ode
        ode_T = casADi.SX.sym("ode_T", self.n_axial)
        for z in range(1, self.n_axial):
            ode_T[z] = (
                1
                / (self.f_rho_cp_eff(w_i[z, :].T, T[z], p[z]))
                * (
                    -(T[z] - T[z - 1])
                    * (u[z] * self.f_rho(w_i[z, :].T, T[z], p[z]) * self.cp_f)
                    / (self.L * self.f_delta_z())
                    - 4 * self.U / self.D * (T[z] - self.T_wall)
                    + self.f_eta(w_i[z, :].T, T[z], p[z])
                    * self.f_r(w_i[z, :].T, T[z], p[z])
                    * (-self.delta_H)
                    * (1 - self.eps)
                    * self.rho_cat
                )
            )

        # boundary condition (inlet, z=0)
        z = 0
        ode_T[z] = (
            1
            / (self.f_rho_cp_eff(w_i[z, :].T, T[z], p[z]))
            * (
                -(T[z] - self.T_in)
                * (u[z] * self.f_rho(w_i[z, :].T, T[z], p[z]) * self.cp_f)
                / (self.L * self.f_delta_z())
                - 4 * self.U / self.D * (T[z] - self.T_wall)
                + self.f_eta(w_i[z, :].T, T[z], p[z])
                * self.f_r(w_i[z, :].T, T[z], p[z])
                * (-self.delta_H)
                * (1 - self.eps)
                * self.rho_cat
            )
        )

        # reshape differential variable and ode into 1D vector
        w_i_reshaped = casADi.reshape(w_i, self.n_axial * len(self.nu), 1)
        ode_wi_reshaped = casADi.reshape(ode_wi, self.n_axial * len(self.nu), 1)
        x = casADi.vertcat(w_i_reshaped, T)
        ode = casADi.vertcat(ode_wi_reshaped, ode_T)
        z = casADi.vertcat(u, p)
        ae = casADi.vertcat(ae_conti, ae_p)

        # DAE struct
        dae = {"x": x, "z": z, "ode": ode, "alg": ae}
        return dae
