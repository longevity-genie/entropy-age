# The models, in equations

This page states the two simulation models of Meyer & Schumacher (2024) exactly, adds the closed
forms the paper does not spell out, and explains why the reimplementation can be much faster than
the original code while sampling from the same distributions.

## 1. Bulk model: bounded random walk (Fig. 1)

A *ground state* $g \in [0,1]^F$ (F features) is the youngest sample. A sample of simulated age $a$ is

$$x^{(0)} = g + \varepsilon_0,\quad \varepsilon_0 \sim \mathcal N(0, 0.01^2),\qquad
x^{(t)} = B\!\left(x^{(t-1)} + \sigma\,\xi_t\right),\quad \xi_t \sim \mathcal N(0, I),\quad t = 1..a,$$

with one of three boundary treatments $B$:

| `Bound` | what happens | result |
|---|---|---|
| `NONE` | nothing, an unbounded random walk | $\mathbb E[x^{(a)}] = g$ for every $a$, so **no clock** (paper: Pearson −0.05) |
| `CLIP` | clip to $[0,1]$ after every step | features near a bound can only move inwards: **regression to the mean**, a clock with r ≈ 0.99 |
| `LOGIT` | add noise on the logit scale, map back with expit at the end | same idea without hard edges (paper Fig. 1, r ≈ 0.95) |

The upstream function `random_epi_logit_everystep` (logit → add → expit on every step) is **the
same model as `LOGIT`**, because $\operatorname{logit}(\operatorname{expit}(z)) = z$. Our test suite
checks that both give identical samples from the same random stream.

Age information exists only because the boundary makes the *expected* trajectory move: features
starting near 0 drift up and near 1 drift down. The elastic net exploits exactly this, so its
coefficients are positive for low ground values and negative for high ones (Fig. 1f in our figure).

`sigma` can be a vector: **per-feature stochasticity** is a first-class input.

## 2. Single-cell methylation: a two-state Markov chain per CpG (Figs. 3-6)

Each CpG $j$ in each cell is methylated (1) or not (0). Per time step:

* a methylated site stays methylated with probability $E_m$ (maintenance efficiency);
* an unmethylated site becomes methylated with probability $E_d$ (de novo efficiency).

For one cell the probability $p_t$ of being methylated obeys $p_{t+1} = E_m p_t + E_d (1 - p_t)$, so
with $r = E_m - E_d$:

$$p_t = z + (p_0 - z)\,r^{\,t},\qquad z = \frac{E_d}{1 - E_m + E_d}.$$

* $z$ is the site's **equilibrium** (Pfeifer et al. 1990); with the paper's universal setting
  $E_d = 1 - E_m$ every site relaxes to $z = 0.5$.
* $\lambda = -\ln r \approx (1 - E_m) + E_d$ is the **relaxation rate**: the per-step rate at which a
  site forgets its ground state. We call it the site's *entropy rate*. Inverting:
  $E_d = z(1 - e^{-\lambda})$, $1 - E_m = (1 - z)(1 - e^{-\lambda})$
  (`MaintenanceRates.from_rate_and_equilibrium`).

### Bulk level and its noise

A sample has $N$ cells (1,000 in the paper). The bulk value is the methylated fraction. Starting
from $k_0 = \lfloor N g_j \rfloor$ methylated cells, cells evolve independently, so

$$\mathbb E[m_t] = z + (k_0/N - z)\,r^{\,t},\qquad
\operatorname{Var}[m_t] = \frac{k_0\,p^{(1)}_t(1-p^{(1)}_t) + (N - k_0)\,p^{(0)}_t(1-p^{(0)}_t)}{N^2},$$

with $p^{(1)}_t = z + (1 - z) r^t$ and $p^{(0)}_t = z - z r^t$ (`binomial_variance`). The age signal
a clock reads is the deterministic relaxation toward $z$; the finite number of cells adds binomial
noise around it. That is why prediction is near perfect at $E_m = 99.9\,\%$ (slow, steady relaxation)
and degrades both for low maintenance (sites hit equilibrium early, the curve plateaus) and for very
high maintenance (the signal per step drowns in binomial noise).

### Why the new simulator is exact and fast

The cells at a site are exchangeable and independent, so the count after one step is exactly

$$k_{t+1} = \operatorname{Binomial}(k_t, E_m) + \operatorname{Binomial}(N - k_t, E_d).$$

Simulating counts instead of flipping each of $N$ cells costs $O(F)$ instead of $O(F N)$ per step.
A second saving is in scheduling: the original re-simulates every age from scratch (work
$\propto \sum_a a = O(a_{\max}^2)$ per replicate set); we process samples in age order and apply
step $t$ only to samples older than $t$, which still gives independent trajectories. The explicit
per-cell simulation is kept (`method="cells"`) and the tests check that both match the closed-form
mean and variance.

### Empirical maintenance rates (the paper's heuristic)

The paper takes the oldest sample as the equilibrium $z_j$, fixes $E_m$ at its lower limit and solves
$E_d = (1 - E_m) z/(1 - z)$; where that $E_d$ exceeds its upper limit, a random $E_d$ below the limit
is drawn and $E_m$ solved instead (`MaintenanceRates.from_equilibrium`). The best limits for Horvath's
clock were $E_m > 97\,\%$ and $E_d < 5\,\%$. We clip $z$ away from exactly 0 and 1, where the
upstream code divides by zero.

From cross-sectional data we can do better than one old sample: `fit_relaxation` fits
$m(t) = z + (m_0 - z)e^{-\lambda t}$ to every feature (Section 4).

## 3. Clocks and age transforms

Clocks are elastic nets (`train_clock`). The paper's settings are in `ClockConfig.paper()`
(CV over $\ell_1$ ratios 0.1..0.9) and `ClockConfig.fixed_alpha()` (real-data clocks: $\alpha = 1$).

| transform | used for | forward | inverse |
|---|---|---|---|
| `LinearRescale(1, 73)` | human stochastic clock | min-max to $[-120, 280]$ | back to steps |
| `LogLogRelativeAge(67)` | pan-mammalian clocks | $-\ln(-\ln(a/a_{\max}))$ | $e^{-e^{-y}}$ (relative age) |
| `HorvathAge()` | Horvath 2013 | log below age 20, linear above | the published anti-transform |

## 4. Extension: selection among cells

`simulate_with_selection` keeps explicit cells because selection couples the sites within a cell.
Per step:

1. every site of every cell epimutates with its own $(E_m, E_d)$;
2. a fraction `turnover` of cells dies at random and is replaced by daughters of cells drawn with
   probability $\propto \exp\!\big(-s \sum_j w_j\,\mathbb 1[x_{ij} \ne x^{\text{ref}}_{ij}]\big)$, where
   $w_j > 0$ marks constrained sites and the reference is the founder state of the cell's lineage
   (selection against *change*) or an explicit functional state;
3. daughters optionally pick up replication errors (one extra epimutation pass scaled by `division_error`).

For a constrained site under strong enough selection, the deviating fraction settles at a
**mutation–selection balance** of order $\mu / (s \cdot \text{turnover})$ instead of relaxing to the
Markov equilibrium $z$. A protected site therefore reaches a *nearby* plateau quickly. Its fitted λ can
be *higher* than that of a neutral site, while the amount it changes is much smaller. Measure
protection with the change amplitude (`window_change`, entropy slope), not with λ alone.

## 5. Per-feature stochasticity estimators

For a cross-sectional cohort (one sample per person, known ages) `stochasticity.py` provides:

* `drift_metrics`: OLS slope and Spearman ρ with age; the age slope of binary entropy
  $H(m) = -m\log_2 m - (1-m)\log_2(1-m)$, which is positive when a site moves toward 0.5; the age slope of
  between-individual dispersion, $|{\rm residual}| / \sqrt{\bar m (1 - \bar m)}$, scaled so that
  bounded features are comparable; and the log ratio of residual variance (old / young).
* `fit_relaxation`: for every candidate λ on a log grid, $z$ and $m_0$ are linear least squares
  (design $[1 - e^{-\lambda t}, e^{-\lambda t}]$), solved for all features at once; λ minimises the
  error. It flags fits at the grid edges and equilibria outside $[0,1]$, and reports `window_change`.
* `compare_groups`: gene set vs. all other features, matched within quantile bins of a baseline
  (difference weighted by in-group size; labels permuted within bins).
* `paired_comparison`: the same features in two cell types (Wilcoxon signed-rank).
* `regress_out`: the paper's cell-composition correction.
