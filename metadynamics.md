# Introduction to Metadynamics

Metadynamics is an enhanced sampling method that uses a set of Collective
Variables (CVs) to overcome insufficient sampling due to high barriers.
Additionally, it offers a way to calculate the Free Energy Surface (FES) w.r.t.
the CVs. 

This documentation will cover the theory of metadynamics and well-tempered
metadynamics (WTMetaD), the choice of parameters and instructions on how to
perform metadynamics with GROMACS patched with Plumed. Especially the
theoretical parts are kept short, for more information refer to the papers
mentioned in [Sources](./metadynamics.md#sources). For additional topics
like [general analysis of a metadynamics simulation](./metad_analysis.md),
[convergence analysis](./block_averaging.md) or
[reweighting](./metad_reweighting.md), refer to their documentations.

## Sources:

- [Original Paper Introducing
  Metadynamics](https://www.nature.com/articles/s42254-020-0153-0)
- [Paper Introducing Well-Tempered
  Metadynamics](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.020603)
- [Review on Theory of
  Practice](https://medien.umbreitkatalog.de/pdfzentrale/978/111/840/Leseprobe_l_9781118407776.pdf)

# Theoretical Background

The quantity of interest for MD simulations is usually the Helmholtz free energy
$F$ (NVT ensemble) or the Gibbs free energy $G$ (NPT ensemble). The free energy
$F$, in the following used as a general term for both quantities, can be
calculated from a trajectory. Given the CV $s$, a histogram $N(s)$ of the
visited conformations in CV space can be created and plugged into the following
equation along with the temperature $T$ and the Boltzmann constant $k_B$:

$$ F(s) = -k_B T \ln N(s) $$

This calculation can be done for any trajectory produced by an MD simulation,
however, we have to guarantee, that the phase space was sampled sufficiently.
This is unfortunately often not the case, as critical conformational changes
happen at a time scale that cannot be reached by free MD simulations (time scale
problem), as there are high barriers between the conformations. Metadynamics
aims to overcome this limitation by adding a bias potential $V$ (blue in the
following figure) to the potential energy of the system in the space spanned by
a set of CVs. $V$ is a history-dependent bias potential built from a sum of
Gaussians, that are often referred to as hills. This bias potential pushes the
system out of a local minimum so that additional configurations are sampled.
After all metastable states in the CV space are filled with Gaussians, the
system can move freely between all CV states (bottom right) and the negative of
the deposited bias potentials can be used to estimate the FES.

![illustration_metad](./metad_fes.png)

A commonly used variant of standard metadynamics is **well-tempered
metadynamics** (WTMetaD). Here, the height of the Gaussians is scaled down
during the simulation, leading to a smoother FES. It needs the additional
parameter $\Delta T$ that is often reflected by the bias factor $\gamma =
\frac{T + \Delta T}{T}$. This leads to the following expression for the free
energy (neglecting an additive constant):

$$ F(s) = - \frac{\gamma}{\gamma -1} V(s) $$

# Guidelines for Choice of Parameters

For standard metadynamics, there are three parameters that need to be chosen:
the width $\sigma$ and height $w_G$ of the Gaussians (hills) as well as the time
interval $\tau_G$ between the deposition. The **width** can be estimated by
creating a histogram of the CV in an unbiased simulation. It is advised to use
one-half or one-third of its standard deviation for the hill width. We have
found, that too wide Gaussians lead to artificial minima at upper walls. Too
small widths lead to a slower filling time and maybe stability problems.

The **hill height** $w_G$ should be a fraction of the thermal energy $k_B T$ and
can be chosen as $0.1$ kcal mol$^{-1}$ for simulations at room temperature.

The **deposition time** $\tau_G$ should be estimated from the autocorrelation
time of the CVs extracted from a short unbiased simulation. First, the
autocorrelation function $R(\tau)$ is calculated as:

$$ R(\tau) = \frac{\langle (X_t - \langle X \rangle) (X_{t+\tau} - \langle X
\rangle) \rangle}{\langle (\delta X)^2 \rangle} $$

where $X_t$ is a time series, for example values of the CVs, with the
expectation value $\langle X \rangle$ and the variance $\langle (\delta X)^2
\rangle$. The value of $R$ at $\tau = 0$ is 1 and it decays to zero as the
values increase. It therefore measures the correlation between the values $\tau$
time steps apart, which is why $\tau$ is also called lag time. The
autocorrelation function should rapidly decay proportional to $\exp(-t/\tau)$
and the so-called integrated autocorrelation time can be calculated as:

$$ \tau_\text{int} = \frac{1}{2} + \sum_{t=1}^{t_{\text{max}}} R(t) $$

It is advised to use a slightly higher time than the autocorrelation time, as
the system should be able to relax after the deposition of a Gaussian, to avoid
positioning the following on top of the first one. However, using the "default"
of $1$ ps seems to be okay for most simulations.

This "default" value probably comes from the integrated autocorrelation time of
alanine dipeptide. Running a short free MD simulation and tracking the backbone
torsion $\phi$, yields $\tau_\text{int} = 0.75$ ps, so a slightly higher value
for the deposition time of 1 ps seems reasonable here. The autocorrelation
function and time can be easily calculated with a Python package:

```python
# Import a module for the calculation of the autocorrelation function (acf)
from statsmodels.tsa.stattools import acf

# Read in the values as an array, here called phi

# Calculate the acf and integrated autocorrelation time (act_int)
lags = 1000  # Amount of lag times
acf_values = acf(phi, nlags=lags)
act_int = 1 + 2 * np.sum(acf_values)
```

Finally, the parameter $\Delta T$ needs to be chosen for WTMetaD. The sum $T +
\Delta T$ can be interpreted as an **effective temperature** of the CVs.
Therefore, the sum should be in the range of the barrier heights. If the barrier
height is expected to be around 8 kcal mol$^{-1}$, this corresponds to 33.5 kJ
mol^{-1} and the bias factor can be calculated for a temperature $T$ of 300 K
as:

$$ \gamma = \frac{T + \Delta T}{T} = \frac{33 500\,\mathrm{J\,mol^{-1}} / 8.314 \mathrm{J\,mol^{-1}
K^{-1}}}{300 \mathrm{K}} = 13.4 $$

So a bias factor of 14-15 should be sufficient to pass the energy barrier. This
ensures, that the hill height is only decreased, after the barrier can be
overcome and all relevant states are sampled.

# Instructions for Metadynamics with GROMACS and Plumed

This guide builds on the previous one for setting up a QM/MM simulation with
GROMACS and DFTB+. For the metadynamics part, Plumed is needed. If the
metadynamics is performed without a QM region, so only with classical MD, the
instructions need to be modified accordingly.

Secondly, please note that the instructions are written for a system with a
disulfide exchange reaction and the distances between the sulfur atoms as CVs.

## Set Up the Working Directory and Submit the Metadynamics Simulation

Copy the following files to the working directory:

- `make_tpr.sh`
- `editPlumed4walker.sh`
- `submit_Mwalker.sh`
- `plumed_metad.dat`
- `r02` (use with velocity generation so the walkers "go" into different directions)


Copy the following files from the directory where the QM/MM simulation was set
up:

- `dftb_in.hsd`
- coordinate file (`.gro`)
- topology file (`topol.top`)
- index file (`index.ndx`)


Edit the `plumed_metad.dat` file as follows:

- Change the atom numbers under 'generalization' (instructions in the Plumed-file).
- Change the arguments under 'Restraints for metadynamic CV's' to the distance
  between the anionic sulfur and the sulfur atom it attacks and the disulfide
  atoms. So if 'cym' attacks 'ss1', use d1,d3 and otherwise d2,d3.
- Change the argument under 'Restraint angles between S' to the angle of the
  attacked sulfur atom. So if 'cym' attacks 'ss1', use angle1 and otherwise
  angle2. If the respective angle of the start structure in the .gro-file is
  lower than 130 degrees, change the value after 'AT=' to the angle in radiant.
- Change under 'Metadynamics' the number of walkers you are using (`WALKERS_N`)
  and the arguments to the distances you chose at 'Restraints for metadynamics
  CV's'.
- Change the name `plumed_metad.dat` to `plumed.dat`.

Copy everything into a backup directory.

Edit the number of walkers in the file `make_tpr.sh` (begins with 0, so 24
walkers would be `0..23`) and create the .tpr-files with the bash script.

Make directories for the walkers and move the .tpr-files into them with the bash
script `editPlumed4walker.sh`. Edit the number of walkers in this file as well.
In the directories there should now be the files `dftb_in.hsd`, `metad_$i.tpr`
and `plumed_$i.dat` where `$i` stands for the respective walker.

Edit the `submit_Mwalker.sh` as follows:

- Adjust the number of walkers after 'nproc' and give the simulation a name after '-N'.
- Adjust the number of walkers also after `n_walker`.
- Comment / uncomment export options.

Submit the metadynamics simulation as follows:

```	
qsub submit_Mwalker.sh
```

## Analyze the Metadynamics Simulation

### Check If Some Walker Crashed and Restart It:

- find out on which node the simulation is running (qstat -> queue), for example
  cn19

- switch directly to that node with `ssh cn19`

- print the number how many simulations there are in total (`pgreb gmx_d | wc
  -l`)

- to check if some walker doesn't work properly enter the following command and
  check if the numbers are similar:

```
cd /scratch/user1/GMX_[job-ID]/
WALC 24
```

- to check which walkers are still running, enter `ps -xo cmd`.

- to check which process ID belongs to which walker, use `ps -xo pid,cmd`.

- If one walker is missing, go to `/scratch/user1/GMX_[job-ID]`. You can find
  the job-ID in qstat and $i stands for the respective walker number. Then use
  `kick_walker $i` to restart the walker.

### Plot the Free Energy Surface (FES)

One of the outputs of every metadynamics simulation in Plumed is one (or more in
the case of multiple walkers) `HILLS` file. It provides information on the
deposited Gaussians every time one is deposited. The first column is the time in
ps, at which the hill was set. The frequency of the deposition was set in the
Plumed input file with the keyword `PACE` (see
[deposition time](#guidelines-for-choice-of-parameters)). The next column
is the value of the CV, where the hill is centered. This is followed by
the width along the CV and one column declaring the height of the hills. In the
case of more than one CV, the center and width along all CVs are documented in
additional columns. Finally, the bias factor can be found in the last column.

The first couple lines of a `HILLS` file for a standard metadynamics of alanine
dipeptide with $\phi$ as a CV looks like this:

```
#! FIELDS time phi sigma_phi height biasf
#! SET multivariate false
#! SET kerneltype stretched-gaussian
#! SET min_phi -pi
#! SET max_phi pi
0.5     -1.325708292819164                    0.3                    1.2                     -1
  1     -2.812889543092375                    0.3                    1.2                     -1
```

The hills are deposited every 0.5 ps, the first one at a CV value of $\phi
\approx 1.3$ rad. The width is 0.3 rad and the height 1.2 kJ/mol. The bias
factor for a standard metadynamics is always -1.

The `HILLS` file for a WTMetaD looks exactly the same, but with another bias
factor. Additionally, the hill height is re-scaled with the bias factor. How
exactly this takes place is not clear until now, but will have an effect later
on.

Plumed has an extra tool called `sum_hills`, that can be used to process these
`HILLS` files to free energy surfaces (FES). For a one-dimensional metadynamics,
it can be used as:

```
plumed sum_hills --hills HILLS
```

There are numerous additional flags, that can be looked up in the
[manual](https://www.plumed.org/doc-v2.9/user-doc/html/sum_hills.html).
For multiple walker metadynamics, there is a script linked in the `~/.bashrc`
called `SUM_HILLS_2D_maag`, that makes it easy to process all `HILLS` files at
once by adding the number of walkers after it. So for 24 walkers, you can run

```
SUM_HILLS_2D_maag 24
```

The resulting file `fes.dat` can be plotted with one of the scripts provided in
the repository called `Python_scripts` and gives the free energy surface.

#### Working Mechanism of `sum_hills`

The information in the Plumed manual regarding `sum_hills` its working mechanism
is unfortunately quite poor. With a bit of reverse engineering, this is what I
have found out so far:

`sum_hills` reads the `HILLS` file (or files in the case of multiple walkers)
and calculates the Gaussians according to:

$$ G(s) = w_G \cdot \exp \left( - \frac{(s - s_0)^2}{2 \sigma ^2} \right) $$

Therefore, it first decides either internally or with the `--bin` flag on the
amount of bins for the CV (or CVs) $s$ and calculates with the equation above
and the information from the `HILLS` file regarding the hill height $w_G$, width
$\sigma$ and center $s_0$, the value of the hill in each bin. The negative sum
over all Gaussians leads directly to the free energy:

$$ F(s) = - \sum_{s_0} G(s) $$

The behavior of `sum_hills` in the example of alanine dipeptide and the CV
$\phi$ and a bin size of 60 can thus be mimicked by the following Python script:

```python
import numpy as np

HILLS = "./HILLS"

phi = np.loadtxt(HILLS, skiprows=5, usecols=1)
width = np.loadtxt(HILLS, skiprows=5, usecols=2)
height = np.loadtxt(HILLS, skiprows=5, usecols=3)

min_phi = np.min(phi)
max_phi = np.max(phi)

bins_phi = np.linspace(min_phi, max_phi, 61)
bin_centers = (bins_phi[:-1] + bins_phi[1:]) / 2
free_energy = np.zeros(60)

for i, value in enumerate(phi):
    gauss = height[i] * np.exp(-((bin_centers - phi[i]) ** 2) / (2 * width[i] ** 2))
    free_energy -= gauss
```

Most interestingly, this script can be applied to both standard and WTMetaD,
meaning that the hill height in the `HILLS` file has been re-scaled in such a
way that the multiplication with $\gamma / (\gamma - 1)$ to get the free energy
is not necessary.

It is also possible to project a higher-dimensional FES on a lower dimension by
integrating biased CVs out. So if we have performed a metadynamics simulation
with two CVs but want to look at the FES w.r.t. only one of them, `sum_hills`
provides this with the `--idw` flag. If the CVs are called `d1` and `d2` and we
want to project the FES on the CV `d1`, this can be done by:

```bash
sum_hills --hills HILLS --idw d1 --outfile fes_d1.dat
```

I have not been able to reproduce exactly the projected FES, but I am pretty
sure that PLUMED uses here the concept of the marginal free energy, which allows
to integrate one (or more) CV $s$ out, projecting it to the other CVs $s'$:

$$ F_\text{proj}(s') = -k_B T \ln \int \exp[-F(s)/ k_B T] \text{d} s $$

Technically, this should yield the same result as
[reweighting](./metad_reweighting.md). I would interpret it as a sign of
non-convergence, if these two methods give deviating free energies. Another
disadvantage of this method is, that both $s$ and $s'$ (so the CV that is
integrated out and the one the FES is projected onto) have to be biased during
the metadynamics simulation, whereas you can choose an arbitrary CV with the
reweighting methods.
