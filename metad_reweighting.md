# Reweighting Protocol for Well-Tempered Metadynamics Simulations in Plumed Using $c(t)$

Sources:

- [Plumed Metadynamics](https://www.plumed.org/doc-v2.9/user-doc/html/_m_e_t_a_d.html)
- [Plumed Masterclass](https://www.plumed.org/doc-v2.9/user-doc/html/master-_i_s_d_d-2.html#master-ISDD-2-ex-3)
- [Plumed Masterclass Video](https://www.youtube.com/watch?v=UgTIGQxJtOc)
- [Tiwary et al. 2015](https://pubs.acs.org/doi/epdf/10.1021/jp504920s?ref=article_openPDF)
- [Valsson et al. 2016](https://pubmed.ncbi.nlm.nih.gov/26980304/)

## Introduction

Metadynamics is a powerful enhanced sampling technique that can be used to
explore the free energy landscape using collective variables (CVs). The free
energy surface (FES) can be described as the natural logarithm of the
equilibrium distribution $P(s)$ of the CVs:

$$ F(s) = - \frac{1}{\beta} \ln P(s) $$

where $\beta = 1/k_BT$ is the inverse temperature. However, the FES obtained
like that may not sample all states of interest. To overcome this limitation, a
bias potential $V(s)$ is introduced to enhance the sampling of the CVs. The
biased distribution $P_V(s)$ is then given by:

$$ P_V(s) = \frac{\exp[-\beta(F(s)+V(s)]}{\int \exp[-\beta(F(s)+V(s)] ds} $$

In metadynamics, the bias potential is constructed as a sum of Gaussians
centered at the current CV value. In **well-tempered metadynamics** (WTMetaD),
the Gaussians are iteratively added to the bias potential with a decreasing
height. This results in a time-dependent bias potential:

$$ V(s,t) = -(1-\frac{1}{\gamma}) F(s) + c(t) $$

where $\gamma$ is the bias factor ($\gamma \rightarrow 1$ yields the unbiased
ensemble and for $\gamma \rightarrow \infty$ the bias potential equals the
negative free energy) and the term $c(t)$ that is **independent of the CVs**.
This term is also called reweighting factor as it is essential for the
reweighting of WTMetaD. To obtain the value of an **arbitrary function** $O(R)$,
that depends on the atom coordinates $R$ and was not biased during the
metadynamics, we need to reweigh the biased ensemble. The biased probability
distribution $P_V(R,t)$ of this arbitrary function can be written as

$$ P_V(R,t) = P(R) \cdot \exp[-\beta V(s,t) - c(t)] $$

where $P(R)$ is the unbiased probability distribution of the function $O(R)$.
The expression for the $P(R)$ can be obtained by dividing with the exponential
function, which removes the negative sign in its exponent:

$$ P(R) = P_V(R,t) \cdot \exp[\beta V(s,t) - c(t)] $$

Therefore, the average value of the function $O(R)$ can be calculated as (the
average is performed over the metadynamics run):

$$ \langle O(R) \rangle = \langle O(R) \cdot \exp(\beta[V(s,t)-c(t)]) \rangle_V $$

Another use of reweighting with $c(t)$ is, that differing free energy surfaces
obtained directly with the bias potential $V(s,t)$ or with the expression
$V(s,t)-c(t)$ are a clear sign of the simulation not converging.

> This results in the following **workflow**:
> 
> 1. Calculate the value of $c(t)$ **during the WTMetaD** for every frame (or the
>    expression $V(s,t)-c(t)$).
> 2. Calculate the weight $w_k = \exp[\beta V(s,t) - c(t)]$ for each frame $k$.
> 3. Create a histogram of the unbiased variable $f$ and reweight each bin that
>    corresponds to the frame $k$ by its weight $w_k$.
> 4. Normalize the histogram (yielding $P(f)$) and calculate the free energy with
>    $F(f) = -k \ln P(f)$.

## PLUMED Specifics

In PLUMED, the expression $V(s,t) - c(t)$ is calculated with the keyword
`CALC_RCT` of the `METAD` module and can be accessed with `rbias`. The
reweighting factor on its own can be printed with `rct`. The keyword `bias` on
the other hand only prints the instantaneous value of $V(s,t)$.

### Plumed Template

In summary, a minimal template for a Plumed input file to make a WTMetaD that
can be reweighted for a different CV could look like this:

```plumed
biased_cv:                      # insert your cv here
unbiased_cv:                    # keep track of your other cv of interest
metad: METAD ARG=biased_cv ...
   PACE=                        # insert frequency for hill additioN
   HEIGHT=                      # insert initial hill height
   BIASFACTOR=                  # insert bias factor for wtmetad
   TEMP=                        # insert temperature, required for bias factor
   SIGMA=                       # insert gaussian width
   FILE=HILLS                   # filename for hills files
   GRID_MIN=                    # minimum value for cvs
   GRID_MAX=                    # maximum value for cvs
   CALC_RCT                     # calculate the reweighting factor
...
# Print all relevant data to a COLVAR file
PRINT ARG=biased_cv,unbiased_cv,metad.rbias FILE=COLVAR_REWEIGHT STRIDE=1
```

### Important Information

The whole scheme is also valid to reconstruct the FES w.r.t. the biased CV only,
according to [Tiwary et al.
2015](https://pubs.acs.org/doi/epdf/10.1021/jp504920s?ref=article_openPDF) this
even counts as a first estimate of the convergence, if reweighting with $c(t)$
yields the same FES as using `sum_hills` (so $F=-V$).

Depending on the complexity of the system, it can make a big difference, if the
initial part of the simulation is discarded or not. It should not make a huge
difference, if the simulation has converged for a long time compared to the time
it took for the convergence, but it should be tested.

If the `rbias` is calculated during the WTMetaD simulation, the reweighting and
even the estimation of convergence as explained in the introduction is valid, as
$c(t)$ then is time-dependent. However, if the `rbias` is determined with a
rerun and Plumed `driver`, $c(t)$ becomes a constant and therefore does not
contribute to the FES (a constant value cancels out as soon as the probability
is calculated). Reweighting is still possible, but only with the [scheme
illustrated later](#reweighting-assuming-a-constant-bias).

## Calculate the Reweighted Free Energy Surface

For the calculation of the reweighted FES, I have created different Python
scripts for different very specific cases. It is planned to structure them
differently in the future in a more general way. But for now, I will present
some key snippets. The available Python script with a short explanation on the
use case can be found in the [Python_scripts repo](https://github.com/maike713/Python_scripts).

### Calculation of the Weights

The first step is to calculate the weights for the reweighting. The weights $w$
are given by:

$$ w = \exp(\beta[V(s,t)-c(t)]) $$

In accordance to the equations above, the unbiased probability distribution can
be obtained by multiplying the biased probability distribution with the weights.
If the instantaneous value of $V(s,t)-c(t)$ was printed by the PLUMED keyword
`rbias`, this can be imported and the weights for each frame can be calculated:

```python
# import data    
data = np.loadtxt('COLVAR_REWEIGHT', comments=['#', '!'])    
time = data[:, 0]    
unbiased_cv = data[:, 2]    
rbias = data[:, 3]    
    
# Constants    
kBT = 2.5  # for 300 K in kJ/mol, so 8.314*300/1000
    
# Weights    
weights = np.exp(rbias / kBT)
```

### Calculation of the Reweighted Free Energy Surface

After the weights have been calculated, the reweighted FES can be calculated
with the `numpy.histogram` function. The weights can be passed directly to the
function and then the histogram can be normalized to obtain a probability
distribution. For two CVs, the function `numpy.histogram2d` can be used in the
same way:

```python
# Reweighted 1D histogram    
cv_sorted = np.sort(unbiased_cv)
cv_bins = np.linspace(cv_sorted[0], cv_sorted[-1], 100)
hist, cv_edges = np.histogram(unbiased_cv, bins=cv_bins, weights=weights)
# Calculate the centers of each bin for plotting
cv_centers = (cv_edges[:-1] + cv_edges[1:]) / 2

# Normalize the histogram
probabilities = hist / np.sum(hist)
```

Using these probabilites, we can now calculate the free energy with the formula
from above. As the FES can only be calculated with a constant additive and we
are only interested in relative energies, we set the minimum value and any
infinities to zero. At the end, we can optionally convert the free energy from
kJ/mol to kcal/mol:

```python
# Calculate the free energy
fes = -kBT * np.log(probabilities)

# Set the minimum free energy to zero
fes -= np.min(fes)

# Handle infinities
fes[np.isinf(fes)] = 0

# Convert the free energy from kJ/mol to kcal/mol
fes = fes * 0.239006
```

### Dealing with Multiple Walkers

When the metadynamics was run with multiple walkers, we need to calculate
reweighting factors for each walker separately. Assuming we are interested in
the FES with respect to the CV `d1` and `d3`, we can accumulate the values of
these CVs and their biases for each walker in an array:

```python
# Initialize  arrays to accumulate data from all walkers    
unbiased_cv = []
bias_values = []
                      
# Loop through each walker to read CV data and bias potentials
for i in range(0, N_WALKERS):    
    cv_data = np.loadtxt(f'WALKER_{i}/COLVAR_REWEIGHT', comments=['!', '#'])    
    time = cv_data[:, 0]    
    cv1 = cv_data[:, 2]    
    rbias = cv_data[:, 3]

    # Append the data to the arrays
    unbiased_cv.extend(cv1)
    bias_values.extend(rbias)
```

Afterwards, the weights and FES can be calculated as shown before.

# Reweighting Assuming a Constant Bias

Sources:
- [Plumed Metadynamics](https://www.plumed.org/doc-v2.9/user-doc/html/_m_e_t_a_d.html)
- [Plumed Masterclass](https://www.plumed.org/doc-v2.9/user-doc/html/master-_i_s_d_d-2.html#master-ISDD-2-ex-3)
- [Plumed Masterclass Video](https://www.youtube.com/watch?v=UgTIGQxJtOc)
- [Bonomi et al.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21305)

## Introduction

All of the instructions for reweighting with $c(t)$ can be transferred to a
slightly less complicated reweighting scheme, that only uses the potential
obtained at the end of the simulation $V(s)$. Here, we assume that the bias is
constant throughout the simulation, which can be a good simplification
neglecting (or discarding) the beginning of the metadynamics simulation, where
the bias potential for a value $s$ of the CV still changes. In a nutshell, this
reweighting technique results, if $c(t)$ is neglected or a constant, but this
will be clear in the following explanations.

Please note that this reweighting method is strictly only valid for WTMetaD and
not for standard metadynamics. Additionally, this method has some disadvantages:

- The reweighting yields only a **joint distribution** $P(s,f)$ between the
  biased CV $s$ and the unbiased variable $f$.
- Therefore, only the joint histogram $N_t(s,f)$ of $s$ and $f$ is available.
- If only the value of the unbiased variable is of interest, this has to be
  projected on this variable only **after** calculating the joint histogram.

As specified in the [Bonomi](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21305),
the probability distribution $P(f)$ of the unbiased variable $f$ is only
available **after** calculating the joint histogram $N_t(s,f)$:

$$ P(f) = \frac{\sum_s \exp[\beta V(s,t) N_t (s,f)]}{\sum_{s,f} \exp[\beta
V(s,t) N_t (s,f]} $$

> This affects the **workflow** as follows:
> 
> 1. Set up a joint histogram $N_t(s,f)$ of the biased CV $s$ and the unbiased
>    variable $f$.
> 2. For each frame reweight the bin that corresponds to $(s_k,f_k)$ by the
>    weight $w = \exp[\beta V(s,t_k)]$.
> 3. Normalize the histogram and calculate the FES according to $F(s,f) = -k \ln
>    P(s,f)$.
> 4. If we are interested how some variable $f$ changes with the biased variable
>    $s$, we are finished here. If we are only interested in the unbiased variable
>    $f$, we need to project the FES onto the variable $f$ (see equation for
>    $P(f)$).


## PLUMED Specifics

As we neglect the time-dependence of the bias potential (and therefore $c(t)$)
completely, we only need $V(s)$ for each value $s$ at the end of the simulation,
which can be accessed by `metad.bias`:

```
# The content of the PLUMED input file is the same, only the PRINT statement
needs to be slightly modified and CALC_RCT is now not necessary
PRINT ARG=biased_cv,unbiased_cv,metad.bias FILE=COLVAR_REWEIGHT STRIDE=1
```

The final bias potential can be calculated during the metadynamics run or in the
post-processing with the PLUMED `driver` tool.

### Calculate the Reweighted Free Energy Surface

The expression for the weights is again similar to the one before but without
$c(t)$:

$$ w = \exp [\beta V(s)] $$

The Python snippets from above can be applied in the same way as for $c(t)$,
with the difference that instead of the `metad.rbias`, the `metad.bias` needs to
be imported and used.

Finally, it should be mentioned that the concept of the weights was introduced
by Tiwary et al. in the [Paper with the reweighting method for
$c(t)$](https://pubs.acs.org/doi/epdf/10.1021/jp504920s?ref=article_openPDF),
which was released after the one by [Bonomi et
al.](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21305) and shortly introduces the
reweighting with a constant bias. The equation for the reweighted free
energy $\widetilde{F}$ as found in the [Bonomi
paper](https://onlinelibrary.wiley.com/doi/10.1002/jcc.21305) is as follows:

$$ \widetilde{F}_N (s,t) = -T \ln N(s,t) - V(s,t) + T \log \int ds' N(s',t) $$

In this equation, it is assumed that $k_B = 1$ and hence is missing. The last sum
with the integral can be neglected in our case, as it is an additive constant
and vanishes by setting the minimum value of the free energy to zero. Including
$k_B$ and neglecting the constant, the equation an be written as:

$$ \widetilde{F}_N(s,t) = -\frac{1}{\beta} \ln N(s,t) - V(s,t) $$

Replacing the histograms $N$ with the probability distributions $P$ and
inserting the relation $F = -\frac{1}{\beta} \ln P$ leads to:

$$ -\frac{1}{\beta} \ln P = -\frac{1}{\beta} \ln P_V - V(s,t) $$

Multiplying by $-\beta$ and applying the exponential function to both sides, we
get:

$$ P = P_V \cdot \exp[\beta V(s,t)] $$

which is exactly the same equation as proposed by Tiwary et al. but without
$c(t)$. So weights can be applied here as well, although they were technically
introduced later.
