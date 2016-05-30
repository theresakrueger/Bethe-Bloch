"""
Calculation of Bethe-Bloch energy loss with:
- density correction
- cut off energy

for electrons as target particles

MPV calculation with Landau theory
"""

from math import sqrt, log
sqr = lambda x: x*x
import matplotlib.pyplot as plt
import numpy as np

""" Constants and formulas to be used for pion """
K = 0.307075 # MeV g^-1 cm^2
z = 1 # pion charge
Z = 14
A = 28.09
Z_over_A = np.float(Z/A) # Silicon
M = 0.511 # MeV is a charged pion mass
m_e = 0.511 # MeV is the electron mass
rho = 2.336 # g cm^-3 is silicon density
gamma = lambda p: sqrt(1 + sqr(p / M))
beta = lambda p: sqrt(1 - 1 / sqr(gamma(p)))
beta_gamma = lambda p: p / m_e
T_max = lambda p: 2 * m_e * sqr(beta_gamma(p)) / (1 + 2 * gamma(p) * m_e / M + sqr(m_e / M)) # MeV
#h_omega_p = 31.05 # eV
h_omega_p = 0.00003105 # MeV

T_kin = lambda p: sqrt( sqr(p) + sqr(m_e) ) - m_e

T_cut = 2.0 # MeV
I = 0.000173 # MeV
#I = 173 # eV

""" Constants of density correction """
C = 4.44
a = 0.1492
m = 3.25
X1 = 2.87
X0 = 0.2014
delta0 = 0.14
X = lambda p: np.log10(beta_gamma(p))
#f1 = lambda x: delta0 * 10**(2*(x-X0)) # conductors pdg
f1 = lambda x: 0 # non conductors pdg
f2 = lambda x: 2 * x * np.log(10) - C + (a * ((X1 - x)**m))
f3 = lambda x: 2 * x * np.log(10) + C
delta_full = lambda x: np.piecewise(x, [x < X0 , x >= X0], [f1, f2])

""" Thickness of absorber """
thickness_um = 50
x = thickness_um/10000.0  # x in cm

epsilon = (K * rho * Z * x) / (2 * A) # Has to have units MeV

E_mean_delta_corr = lambda p:  epsilon * (log(2 * m_e * T_max(p) * sqr(beta_gamma(p)) / sqr(I)) - 2 * sqr(beta(p)) - delta_full( X(p) ) ) / sqr(beta(p)) # 
E_mean = lambda p:  epsilon * (log(2 * m_e * T_max(p) * sqr(beta_gamma(p)) / sqr(I)) - 2 * sqr(beta(p)) ) / sqr(beta(p)) # 
E_mean_delta_corr_Tcut = lambda p:  epsilon * (log(2 * m_e * T_cut * sqr(beta_gamma(p)) / sqr(I)) - 2 * sqr(beta(p)) * (1 + T_cut/T_max(p)) - delta_full( X(p) ) ) / sqr(beta(p)) # 
E_MPV_LVB = lambda p:  epsilon * ( log( ( 2 * m_e * sqr( beta_gamma(p) ) ) / I ) + log( epsilon / ( I * sqr(beta(p)) ) ) +  0.2 - sqr( beta(p) ) - delta_full( X(p) ) )  / sqr(beta(p)) # pdg
E_MPV_limit = lambda p:  epsilon * ( log( ( 2 * m_e * epsilon ) / sqr(h_omega_p) )  +  0.2 )


f = lambda x: ( ( 2 * x - 1 ) / sqr(x) ) * log( 2 ) - ( sqr( (x - 1) / x ) / 8 )

E_mean_electrons = lambda p:  epsilon * ( log( ( m_e * sqr( beta_gamma(p) ) * T_max(p) ) / sqr(I) ) - sqr( beta(p) ) - delta_full( X(p) ) + f( gamma(p) ) )  / sqr(beta(p)) # pdg

I_eV = I*1000000.0
E_low_energy_approx = lambda x:  7.8 * 1E7 * Z * log( ( 1.166 * x / I_eV ) )  / ( A * x ) # pdg


ps = np.logspace(-2, 2, 1000) # MeV

dE_over_dX_pions = np.zeros(len(ps))
dE_over_dX_Tcut_pions = np.zeros(len(ps))

beta_gammas = np.zeros(len(ps))
Xs = np.zeros(len(ps))
density_correction = np.zeros(len(ps))
betas = np.zeros(len(ps))
ps_GeV = np.zeros(len(ps))

E_means = np.zeros(len(ps))
E_mean_density_corr = np.zeros(len(ps))
E_mean_density_corr_Tcut = np.zeros(len(ps))

E_MPVs_Landau = np.zeros(len(ps))
E_MPVs_Landau_limit = np.zeros(len(ps))

E_means_electrons = np.zeros(len(ps))
T_kins = np.zeros(len(ps))

E_low_energy_approxs = np.zeros(len(ps))
T_kins_eV = np.zeros(len(ps))

for pi in range(0, len(ps)):
	# dE_over_dX_pions[pi] = dE_over_dX(ps[pi])
	# dE_over_dX_Tcut_pions[pi] = dE_over_dX_Tcut(ps[pi])

	beta_gammas[pi] = beta_gamma(ps[pi])
	density_correction[pi] = delta_full(ps[pi])
	betas[pi] = beta(ps[pi])
	
	Xs[pi] = X( ps[pi] )	
	E_mean_density_corr[pi] = E_mean_delta_corr( ps[pi] ) * 1000
	E_mean_density_corr_Tcut[pi] = E_mean_delta_corr_Tcut( ps[pi] ) * 1000
	E_means[pi] = E_mean( ps[pi] ) * 1000
	E_MPVs_Landau[pi] = E_MPV_LVB( ps[pi] ) * 1000
	E_MPVs_Landau_limit[pi] = E_MPV_limit( ps[pi] ) * 1000
	ps_GeV[pi] = ps[pi]/1000.0

	T_kins[pi] = T_kin( ps[pi] )
	
	E_means_electrons[pi] = E_mean_electrons( ps[pi] ) * 1000
	
	T_kins_eV[pi] = T_kin( ps[pi] ) * 1000000.0
	E_low_energy_approxs[pi] = E_low_energy_approx( T_kins_eV[pi] ) / 1000.0
	
print ' betas = ', betas
dE_pions_450um = np.multiply(dE_over_dX_pions, 0.045*1000)
dE_pions_Tcut_450um = np.multiply(dE_over_dX_Tcut_pions, 0.045*1000)


fig = plt.figure('Density correction control plot of X=log10(beta(p)*gamma(p))', figsize=(10, 5))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,2,1)
my_font = 16
plt.tick_params(axis='both', which='major', labelsize=my_font)
plt.plot(ps, Xs, label='X', marker='o', color='Gray', linestyle='', linewidth=5, alpha=0.5)
plt.grid(True, which="both", ls="-", color='0.65')
plt.axhline(X0, label = 'X0', color='Black', linewidth=2)
plt.axhline(X1, label = 'X1', color='Red', linewidth=2)
plt.ylabel(r'X = log($\beta\gamma$ (p))', fontsize = my_font)
plt.xlabel('Momentum p [MeV]', fontsize = my_font)
plt.legend(loc='center right')

ax2 = fig.add_subplot(1,2,2)
my_font = 16
plt.tick_params(axis='both', which='major', labelsize=my_font)
plt.plot(ps, T_kins, label='', marker='o', color='Gray', linestyle='', linewidth=5, alpha=0.5)
plt.grid(True, which="both", ls="-", color='0.65')
plt.ylabel('Kinetic energy electron [MeV]', fontsize = my_font)
plt.xlabel('Momentum p [MeV]', fontsize = my_font)
#plt.legend(loc='center right')
ax2.set_yscale('log')
ax2.set_xscale('log')
plt.tight_layout()

plt.savefig('Bethe-control-plot-electrons.png')
plt.savefig('Bethe-control-plot-electrons.svg')
#plt.show()
plt.clf()
plt.close()

factor_thickness = thickness_um/10000.0 # 1/5 * 10^(-2+6-1)   #(1cm/50um)
factor = (rho*1000.0)*factor_thickness
print 'Factor = ', factor
a = np.loadtxt("NIST-data-for-electrons.txt")
kin = a[:,0]
coll = np.multiply(a[:,1], factor )
rad = np.multiply(a[:,2], factor )
tot = np.multiply(a[:,3], factor )
dens = np.multiply(a[:,4], factor )


fig, ax1 = plt.subplots()
fig.patch.set_facecolor('white')
fig.set_size_inches(8,6)

plt.tick_params(axis='both', which='major', labelsize=my_font)
ax1.plot(T_kins, E_means_electrons, label='Bethe', marker='', color='Blue', linestyle=':', linewidth=2)

plt.plot(kin, coll, label='NIST Collision')
plt.plot(kin, rad, label='NIST Radiative')
plt.plot(kin, tot, label='NIST Total')

E_low_energy_approxs_per_50 = np.multiply(E_low_energy_approxs, factor )
print 'E_low_energy_approxs = ', E_low_energy_approxs
plt.plot(T_kins, E_low_energy_approxs_per_50, label='Low energy approx', linewidth='2')

plt.grid(True, which="both", ls="-", color='0.65')
ax1.set_xscale('log')
ax1.set_yscale('log')
plt.ylabel('Energy loss of electrons in '+str(thickness_um)+' um Si in keV', fontsize = my_font)
plt.xlabel('Kinetic energy [MeV]', fontsize = my_font)
plt.legend(loc='upper center')

# ax2 = ax1.twinx()
# energy_to_charge = 3.6
# s2 = np.divide(E_MPVs_Landau, energy_to_charge)
# ax2.plot(ps_GeV, s2, marker='', color='Blue', linestyle='-', linewidth=2)
# ax2.set_ylim(10/energy_to_charge, 45/energy_to_charge)
# ax2.set_ylabel('Signal in 1000 e/h pairs', fontsize = my_font)

plt.ylim(10, 1000)
plt.xlim(0.001, 200)
plt.tick_params(axis='both', which='major', labelsize=my_font)
plt.tight_layout()
plt.savefig('Bethe-for-thickness-'+str(thickness_um)+'um-electrons.png')
plt.savefig('Bethe-for-thickness-'+str(thickness_um)+'um-electrons.svg')
plt.show()

plt.clf()
plt.close()
