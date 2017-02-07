import numpy
import scipy
import scipy.constants as const
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fun_read_xvg import read_xvg_energy
import pycse.orgmode as org

charge_per_atom = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30]
charge_per_atom += [-1, -2, -3, -4, -5, -6, -8, -10, -12, -15, -20, -30]
charge_per_atom.sort()

c_atom_to_sigma = lambda x: x*2/(2.465e-8**2*scipy.sin(scipy.pi/3))


# Convert the adhesion energy from

A_c = 15.1e-18                  # area of the whole plane in m^2

f_base = "../data/MD/E_int_%d.xvg"


vdW_tot = []
vdW_err = []
coulomb_tot = []
coulomb_err = []
potential_tot = []
potential_err = []
coul_LR = []

f_0 = f_base % 0
data = read_xvg_energy(f_0)
vdw0 = data["LJ(SR)"]["Average"] + data["Disper.corr."]["Average"]
coul0 = data["Coulomb(SR)"]["Average"] + data["Coul.recip."]["Average"]
potential0 = data["Potential"]["Average"]
coul_LR_0 = data["Coul.recip."]["Average"]

for e in charge_per_atom:
    f_n = f_base % e
    # print(f_n)
    data = read_xvg_energy(f_n)
    vdw = data["LJ(SR)"]["Average"] + data["Disper.corr."]["Average"]
    vdw_err = data["LJ(SR)"]["RMSD"] + data["Disper.corr."]["RMSD"]
    # coul = data["Coulomb(SR)"]["Average"]
    coul = data["Coulomb(SR)"]["Average"] + data["Coul.recip."]["Average"]
    coul_err = data["Coulomb(SR)"]["RMSD"] + data["Coul.recip."]["RMSD"]
    _coul_LR = data["Coul.recip."]["Average"]
    potential = data["Potential"]["Average"]
    potential_err_ = data["Potential"]["RMSD"]
    # print(vdw, coul)
    vdW_tot.append(vdw-vdw0)
    coulomb_tot.append(coul-coul0)
    vdW_err.append(vdw_err)
    coulomb_err.append(coul_err)
    # potential_tot.append(potential-potential0-_coul_LR)
    potential_tot.append(potential-potential0)
    potential_err.append(potential_err_)
    # coul_LR.append(_coul_LR)

charge_per_atom = numpy.array(charge_per_atom)*0.001
# sigma = c_atom_to_sigma(charge_per_atom)
n_2D = c_atom_to_sigma(charge_per_atom)/10**13
vdW_tot = numpy.array(vdW_tot)/A_c/const.N_A*10**6
vdW_err = numpy.array(vdW_err)/A_c/const.N_A*10**6
coulomb_tot = numpy.array(coulomb_tot)/A_c/const.N_A*10**6
coulomb_err = numpy.array(coulomb_err)/A_c/const.N_A*10**6
potential_tot = numpy.array(potential_tot)/A_c/const.N_A*10**6
potential_err = numpy.array(potential_err)/A_c/const.N_A*10**6
# nn = numpy.linspace(-5, 5, 100)
# params = numpy.polyfit(n_2D, vdW_tot, 2)
# f = numpy.poly1d(params)
# vv = f(nn)

def plot_Phi_charge(fig, error=False):
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()           # For the charge
    # ax3 = ax1.twinx()           # For the surface tension
    l_tot = ax1.plot(n_2D, potential_tot, 's-',
             label=r"$\Delta \Phi_{\mathrm{2D}}^{or}$")
    l_vdw = ax1.plot(n_2D, vdW_tot, 's-',
             label=r"$\Delta \Phi_{\mathrm{LJ}}$")
    l_cl = ax1.plot(n_2D, coulomb_tot, 's-',
             label=r"$\Delta \Phi_{\mathrm{CL}}$")
    if error is True:
        ax1.fill_between(sigma/10**13,
                     vdW_tot-vdW_err, vdW_tot+vdW_err,
                     alpha=0.2, facecolor="blue")
        ax1.fill_between(sigma/10**13,
                     coulomb_tot-coulomb_err, coulomb_tot+coulomb_err,
                     alpha=0.2, facecolor="orange")
        ax1.fill_between(sigma/10**13,
                     potential_tot-potential_err, potential_tot+potential_err,
                     alpha=0.2, facecolor="green")
    # ax1.plot(nn, vv, color=l_vdw[0].get_color(), alpha=0.6)
    ax1.set_xlabel(r"$\sigma_{\mathrm{2D}}$ ($10^{13}$ $e\cdot$cm$^{-2}$)")
    ax1.set_ylabel(r"$\Delta \Phi$ (mJ$\cdot$m$^{-2}$)")
    ax1.legend(loc=1, bbox_to_anchor=(0.95, 1.0))
    ax1.set_xlim(-5, 5)
    ax1.set_ylim(-10, 15)
    # Change the second x axis

    ax2_ticks = numpy.linspace(-0.012, 0.012, 7)
    ax2.set_xticks(c_atom_to_sigma(ax2_ticks)/10**13)
    ax2.set_xticklabels(list(map(lambda s: "%.0f" % s, ax2_ticks*1000)))
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlabel("$\sigma_{\mathrm{2D}}$ (10$^{-3}$ $e$/atom)")
    fig.tight_layout(pad=0)

# ax1.set_xlim(-20, 20)



# ax2_ticks = numpy.linspace(-0.03, 0.03, 7)
# ax2.set_xlim(ax1.get_xlim())
# ax2.set_xticks(c_atom_to_sigma(ax2_ticks)/10**13)
# ax2.set_xticklabels(list(map(str, ax2_ticks)))
# # ax2.plot(charge_per_atom, potential_tot, alpha=0)
# ax2.set_xlabel("Unit charge per atom", labelpad=10)


# # print(ax1.get_ylim())
# # print(ax1.get_yticks())
# ax3.set_yticks(ax1.get_yticks())
# ax3.set_ylim(ax1.get_ylim())
# ax3_yticks = ax1.get_yticks()/A_c/const.N_A*10**6
# ax3.set_yticklabels(list(map(lambda a: "%.1f"%a, ax3_yticks)))
# # ax3.plot(sigma/10**13, potential_tot/A_c/const.N_A*1000, alpha=0.0)
# ax3.set_ylabel(r"$\Delta\gamma_{\mathrm{WG}}$ [mJ$\cdot$m$^{-2}$]", labelpad=-2)


# org.figure(plt.savefig("../img/e-vdw.png"))
fig = plt.figure()

if __name__ == "__main__":
    plot_Phi_charge(fig)
    org.figure(plt.savefig("../img/e-vdw-2.pdf"))
