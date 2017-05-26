import numpy, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.constants as const
import scipy
import pycse.orgmode as org

charge_per_atom = [-12, 0, 12]
name = ["neg0.012", "0.000", "0.012"]

c_atom_to_sigma = lambda x: x*2/(2.465e-8**2*scipy.sin(scipy.pi/3))
z_gr = 2.177

f_charge_base = "charge_int_{}_large2.xvg"
f_charge_water = "charge_int_water-surf.xvg"

f_dens_base = "density_int_{}_large2.xvg"
f_dens_water = "density_int_water-surf.xvg"

charge_per_atom.sort()

c_water = numpy.genfromtxt(f_charge_water, delimiter=(12, 17), skip_header=19)
d_water = numpy.genfromtxt(f_dens_water, delimiter=(12, 17), skip_header=19)

# ax1.plot(c_water[:, 0] - z_gr, c_water[:, 1], label="Water Only")

def plot_den(fig, what="mass"):
    ax = fig.add_subplot(111)
    if what is "mass":
        for index, c in enumerate(charge_per_atom):
            d_sys = numpy.genfromtxt(f_dens_base.format(name[index]),
                                     delimiter=(12, 17), skip_header=19)
            ax.plot(d_sys[:, 0] - z_gr,
                    d_sys[:, 1], label=r"%d$\times10^{-3}$ $e$/atom" % (c))
        ax.set_ylabel(r"$\rho_{\mathrm{w}}$ (kg$\cdot$m$^{-3}$)")
        ax.set_xlabel(r"$z$ (nm)")
        ax.set_xlim(0, 1.5)
        ax.legend(loc=0)
    elif what is "charge":
        for index, c in enumerate(charge_per_atom):
            c_sys = numpy.genfromtxt(f_charge_base.format(name[index]),
                                     delimiter=(12, 17), skip_header=19)
            ax.plot(c_sys[:, 0] - z_gr, c_sys[:, 1],
                    label=r"%d$\times10^{-3}$ $e$/atom" % (c) )
        ax.set_ylabel(r"$\delta_{\mathrm{w}}$ ($e\cdot$nm$^{-3}$)")
        ax.set_xlabel(r"$z$ (nm)")
        ax.set_xlim(0, 1.5)
        ax.legend(loc=0)

    fig.tight_layout(pad=0)

if __name__ == "__main__":
    fig = plt.figure()
    plot_den(fig, what="mass")
    org.figure(plt.savefig("img/density_m.pdf"))
    plt.cla()
    fig = plt.figure()
    plot_den(fig, what="charge")
    org.figure(plt.savefig("img/density_c.pdf"))




