# Method: why every condition is measured twice

## The problem

Muscle force is driven through a sigmoid on the rhythm generator's population
firing rate:

```
df = 1 / (1 + exp(-K_F * (r_rgf / 100 - X0_F)))
```

with `K_F = 30` and `X0_F = 1.42`, i.e. a steep threshold at 142 Hz.

RG-F's firing band is not fixed. It moves with the network:

| network | RG-F median |
|---|---|
| m = 0.40 (480 neurons) | 141 Hz |
| m = 0.30 (360 neurons) | 123 Hz |
| m = 0.20 (240 neurons) | 83 Hz |
| 100 neurons total | 40 Hz |

A threshold fixed at 142 Hz therefore sits above the band of every reduced
network and reads a flexor that is not there. The error is **systematic, not
random** — it grows with how far the band has moved — so it penalises small
networks specifically, which is exactly the shape of bias that invents a
"neurons matter more than synapses" result.

Measured size of the error, same runs, two gates:

| condition | fixed gate | own gate | error |
|---|---|---|---|
| k=0.25 m=0.30 | -0.579 | -0.896 | **0.317** |
| k=0.15 m=0.30 | -0.599 | -0.796 | 0.197 |
| k=0.25 m=0.40 | -0.921 | -0.946 | 0.025 |
| k=0.20 m=0.40 | -0.943 | -0.943 | 0.000 |

## The protocol

**Pass A** runs every condition at the stock gate. This measures where the band
sits. Its `corr_F` values are not evidence about the circuit.

**Plan** takes the median prescription across each condition's seeds:

```
X0F_rec = (p10 + p90) / 2 / 100
```

One gate per condition, not per run, so seed-to-seed variability stays in the
measurement instead of being absorbed by a per-run readout.

**Pass B** re-runs every condition at its own gate. This is the result.

The gate is not purely a readout — force drives Ia afferents, which feed back
into the network — but one pass converges: pass B's prescription reproduces its
own input to within 0.01.

## Metrics

| metric | meaning |
|---|---|
| `corr_F` | corr(force_E, force_F). The headline number. Depends on the gate. |
| `corr_RG` | corr(rate_RG-E, rate_RG-F). **Gate-independent** — valid in both passes. |
| `ampF` | flexor force p99 − p1 |
| `troughF` | flexor force p1. How completely the off-phase switches off. |
| `rgf_span` | RG-F p90 − p10, the operating band |
| `settling_state` | drift check by thirds of the run |

**`corr_RG` is the diagnostic.** It is measured on the rates before any readout,
so it says whether the half-centre alternates at all, independent of how the
gate is set. In a correct pass B, `corr_F` tracks `corr_RG` to within 0.02. A
large gap means the readout is losing a rhythm the circuit is producing.

`troughF` binds before `corr_F` does. Configurations that hold baseline sit at
3.7–3.9; degraded ones climb to 5.5–7.5 while `corr_F` still looks acceptable.

## Rules that came out of getting this wrong

- **8 seeds minimum.** Two of three originally recommended operating points were
  single lucky draws from distributions with SD 0.05–0.06.
- **Compare against the real k=1.00 baseline** (-0.953), not against the
  reduced model's own reference cell.
- **Report `corr_RG` alongside `corr_F`.** If they disagree, the readout is the
  problem, not the circuit.
- **One thread.** NEST seeds per virtual process; thread count changes results.
