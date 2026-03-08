# Open-Economy 3-Equation Model — Carlin & Soskice (2023)

R implementation of the open-economy three-equation macroeconomic model
presented in Chapter 11 of *Macroeconomics: Institutions, Instability,
and the Financial System* by Wendy Carlin and David Soskice (OUP, 2023).

Developed as teaching material for advanced undergraduate economics courses.

## Model overview

The model extends the closed-economy IS-PC-MR framework to a small open
economy with a flexible exchange rate. It consists of five equations:

| Equation | Expression | Description |
|---|---|---|
| Phillips Curve (PC) | π_t = π_{t-1} + α(y_t − y_e) | Inflation dynamics |
| Monetary Rule (MR) | y_t − y_e = −αβ(π_t − π^T) | Central bank loss minimization |
| Open-economy IS | y_t = A − ar_{t-1} + bq_{t-1} | Aggregate demand with net exports |
| Real UIP | r_t − r* = q^E_{t+1} − q_t | Uncovered interest parity condition |
| RX curve | y_t − y_e = −(a + b/(1−λ))(r_{t-1} − r*) | Policy implementation (replaces IS) |

The RX curve is flatter than the IS curve because the exchange rate channel
amplifies the effect of interest rate changes on aggregate demand.

## Repository structure
```
├── modelo_funciones.R       # Core model: parameters, equations, simulation
│                            # engine, and ggplot2 visualization functions
└── modelo_simulaciones.R    # Shock simulations and guided experiments
```

## Shocks covered

- **Inflation shock** — temporary upward shift in the Phillips curve
- **Permanent demand shock** — shift in autonomous demand (positive or negative)
- **Supply shock** — permanent shift in equilibrium output y_e (ERU curve)
- **World interest rate shock** — change in r* (relevant for small open economies)

Each simulation produces:
- Impulse response functions (output gap, inflation, interest rate, real exchange rate)
- PC–MR diagram (inflation–output space)
- IS–RX diagram (interest rate–output space)
- AD–ERU diagram (real exchange rate–output space)

## Key results illustrated

1. The central bank raises the interest rate **less** in the open economy than
   in the closed economy — the exchange rate channel does part of the
   stabilization work (RX is flatter than IS).
2. The real exchange rate **overshoots** its new medium-run equilibrium on
   impact, then gradually returns — a result of rational expectations in the
   forex market combined with sticky prices (Dornbusch, 1976).
3. Permanent demand and supply shocks change the **equilibrium real exchange
   rate** q̄, while the real interest rate remains anchored at r*.

## Requirements
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork"))
```

Tested on R ≥ 4.2.

## Usage
```r
source("modelo_funciones.R")   # Load model core
source("modelo_simulaciones.R") # Run all simulations
```

All parameters are set in a single `param` list at the top of
`modelo_simulaciones.R`. Edit values there and re-run to explore
how the dynamics change.

## References

- Carlin, W. & Soskice, D. (2023). *Macroeconomics: Institutions,
  Instability, and the Financial System*. Oxford University Press. Chapter 11.
- Dornbusch, R. (1976). Expectations and Exchange Rate Dynamics.
  *Journal of Political Economy*, 84(6), 1161–1176.

