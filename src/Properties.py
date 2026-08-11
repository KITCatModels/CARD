import casadi as casADi

from .Parameters import Parameters


class Properties(Parameters):
    def __init__(self):
        super().__init__()

    # define implicit variables
    def f_Mw(self, w_i):
        unity = casADi.DM.ones(len(self.Mw_i))
        Mw = 1 / (casADi.dot(w_i / casADi.SX(self.Mw_i), unity))
        return Mw

    def f_cf(self, T, p):
        cf = p / (self.R * T)
        return cf

    def f_rho(self, w_i, T, p):
        rho = self.f_cf(T, p) * self.f_Mw(w_i)
        return rho

    def f_ci(self, w_i, T, p):
        c_i = w_i * self.f_Mw(w_i) / self.Mw_i * self.f_cf(T, p)
        return c_i

    def f_r(self, w_i, T, p):
        r = self.f_ci(w_i, T, p)[0] * self.k_0 * casADi.exp(-self.E_A / (self.R * T))
        return r

    def f_rho_cp_eff(self, w_i, T, p):
        rho_cp_eff = (
            self.eps * self.f_rho(w_i, T, p) * self.cp_f
            + (1 - self.eps) * self.rho_cat * self.cp_cat
        )
        return rho_cp_eff

    def f_thiele(self, w_i, T, p):
        thiele = (
            self.d_cat
            / 2
            * (
                (self.f_r(w_i, T, p) * self.rho_cat)
                / (self.D_A * self.f_ci(w_i, T, p)[0])
            )
            ** 0.5
        )
        return thiele

    def f_eta(self, w_i, T, p):
        eta = (
            3
            / self.f_thiele(w_i, T, p)
            * (1 / casADi.tanh(self.f_thiele(w_i, T, p)) - 1 / self.f_thiele(w_i, T, p))
        )
        return eta
