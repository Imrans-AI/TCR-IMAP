# Table S9. Training-size-matched control: estimands, diagnostics and memberships

**Table S9.** The control specified in `results/review/SIZE_MATCHED_CONTROL_PLAN.md`, frozen before any model was
trained, and reported in `SIZE_MATCHED_CONTROL_RESULT.md`. For each threshold and fold the number of training
complexes retained by the similarity filter was drawn instead at random from the unrestricted pool, with no
similarity filter, on the same test complexes; five draws with seeds fixed in advance (1001–1005), three training
seeds each, 225 fits in total.

## S9a. Estimands

The **principal** estimate is the mean over draws of the per-complex difference between the model and the pair-frequency prior, bootstrapped over the 69 clu70 clusters (20,000 resamples, seed 42), and is **conditional on the five realized training draws**.
The **draw range** is a five-point spread over a shared test set and is **not a confidence interval**.
The **paired contrast** is computed per complex on the same 174 test complexes under the same cluster resamples.
The **ensemble** quantity averages predictions across draws before scoring; because average precision is nonlinear it is a different estimand and is reported as **secondary**.
All p values fall below the 1/20,000 = 5 × 10⁻⁵ floor and are reported as p < 0.001.

| threshold | principal Δ_rand | lo | hi | draw min | draw max | draw range | Δ_sim | lo  | hi  | paired contrast | lo   | hi   | ensemble Δ_rand (secondary) | lo    | hi    |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.95 | 0.2400 | 0.1827 | 0.2832 | 0.2302 | 0.2472 | 0.0170 | 0.0777 | 0.0336 | 0.1166 | −0.1623 | −0.2090 | −0.1163 | 0.2517 | 0.1906 | 0.2979 |
| 0.85 | 0.2364 | 0.1748 | 0.2822 | 0.2290 | 0.2419 | 0.0129 | 0.0428 | 0.0101 | 0.0747 | −0.1936 | −0.2434 | −0.1364 | 0.2524 | 0.1825 | 0.3038 |
| 0.7 | 0.2296 | 0.1687 | 0.2747 | 0.2134 | 0.2550 | 0.0416 | 0.0274 | 0.0011 | 0.0542 | −0.2021 | −0.2455 | −0.1471 | 0.2430 | 0.1794 | 0.2910 |

## S9b. Similarity and training-size diagnostics

Three distinct quantities, kept separate and never conflated.
The **maximum over all evaluated test complexes** of a test complex's similarity to its fold's training set.
The **mean of the five fold maxima**, which is smaller.
The **mean over test complexes of each complex's nearest-training similarity**, smaller again.
Training sizes are exact per fold; any mean is labelled as a mean.

| regime | threshold | draw seed | training per fold | training mean | max over all test complexes | mean of fold maxima | mean per-test nearest-training |
|---|---|---|---|---|---|---|---|
| similarity-filtered | 0.95 | - | 91,94,96,95,94 | 94.0 | 0.9483 | 0.9330 | 0.7190 |
| size-matched random | 0.95 | 1001 | 91,94,96,95,94 | 94.0 | 1.0000 | 1.0000 | 0.8531 |
| size-matched random | 0.95 | 1002 | 91,94,96,95,94 | 94.0 | 1.0000 | 1.0000 | 0.8654 |
| size-matched random | 0.95 | 1003 | 91,94,96,95,94 | 94.0 | 1.0000 | 1.0000 | 0.8610 |
| size-matched random | 0.95 | 1004 | 91,94,96,95,94 | 94.0 | 1.0000 | 1.0000 | 0.8640 |
| size-matched random | 0.95 | 1005 | 91,94,96,95,94 | 94.0 | 1.0000 | 1.0000 | 0.8548 |
| similarity-filtered | 0.85 | - | 86,86,88,93,86 | 87.8 | 0.8413 | 0.8376 | 0.6823 |
| size-matched random | 0.85 | 1001 | 86,86,88,93,86 | 87.8 | 1.0000 | 1.0000 | 0.8568 |
| size-matched random | 0.85 | 1002 | 86,86,88,93,86 | 87.8 | 1.0000 | 1.0000 | 0.8480 |
| size-matched random | 0.85 | 1003 | 86,86,88,93,86 | 87.8 | 1.0000 | 1.0000 | 0.8595 |
| size-matched random | 0.85 | 1004 | 86,86,88,93,86 | 87.8 | 1.0000 | 1.0000 | 0.8624 |
| size-matched random | 0.85 | 1005 | 86,86,88,93,86 | 87.8 | 1.0000 | 1.0000 | 0.8535 |
| similarity-filtered | 0.7 | - | 72,69,76,78,73 | 73.6 | 0.6992 | 0.6964 | 0.6495 |
| size-matched random | 0.7 | 1001 | 72,69,76,78,73 | 73.6 | 1.0000 | 1.0000 | 0.8435 |
| size-matched random | 0.7 | 1002 | 72,69,76,78,73 | 73.6 | 1.0000 | 1.0000 | 0.8483 |
| size-matched random | 0.7 | 1003 | 72,69,76,78,73 | 73.6 | 1.0000 | 1.0000 | 0.8563 |
| size-matched random | 0.7 | 1004 | 72,69,76,78,73 | 73.6 | 1.0000 | 1.0000 | 0.8307 |
| size-matched random | 0.7 | 1005 | 72,69,76,78,73 | 73.6 | 1.0000 | 1.0000 | 0.8297 |
| unrestricted | 1.01 | - | 139,139,139,139,140 | 139.2 | 1.0000 | 1.0000 | 0.8897 |

## S9c. Per-draw results

| threshold | draw seed | n | model | prior | Δ | lo | hi | p |
|---|---|---|---|---|---|---|---|---|
| 0.95 | 1001 | 174 | 0.5780 | 0.3478 | 0.2302 | 0.1693 | 0.2760 | 0.0001 |
| 0.95 | 1002 | 174 | 0.5961 | 0.3489 | 0.2472 | 0.1880 | 0.2936 | 0.0001 |
| 0.95 | 1003 | 174 | 0.5950 | 0.3528 | 0.2422 | 0.1825 | 0.2883 | 0.0001 |
| 0.95 | 1004 | 174 | 0.5964 | 0.3507 | 0.2457 | 0.1863 | 0.2900 | 0.0001 |
| 0.95 | 1005 | 174 | 0.5903 | 0.3554 | 0.2349 | 0.1758 | 0.2791 | 0.0001 |
| 0.95 | mean_over_draws | 174 | 0.6087 | 0.3570 | 0.2517 | 0.1906 | 0.2979 | 0.0001 |
| 0.85 | 1001 | 174 | 0.5832 | 0.3518 | 0.2314 | 0.1719 | 0.2761 | 0.0001 |
| 0.85 | 1002 | 174 | 0.5920 | 0.3500 | 0.2419 | 0.1730 | 0.2932 | 0.0001 |
| 0.85 | 1003 | 174 | 0.5935 | 0.3525 | 0.2410 | 0.1782 | 0.2880 | 0.0001 |
| 0.85 | 1004 | 174 | 0.5896 | 0.3511 | 0.2385 | 0.1745 | 0.2874 | 0.0001 |
| 0.85 | 1005 | 174 | 0.5827 | 0.3536 | 0.2290 | 0.1669 | 0.2761 | 0.0001 |
| 0.85 | mean_over_draws | 174 | 0.6099 | 0.3575 | 0.2524 | 0.1825 | 0.3038 | 0.0001 |
| 0.7 | 1001 | 174 | 0.5761 | 0.3468 | 0.2293 | 0.1643 | 0.2778 | 0.0001 |
| 0.7 | 1002 | 174 | 0.5798 | 0.3499 | 0.2299 | 0.1658 | 0.2782 | 0.0001 |
| 0.7 | 1003 | 174 | 0.6073 | 0.3524 | 0.2550 | 0.1926 | 0.3025 | 0.0001 |
| 0.7 | 1004 | 174 | 0.5590 | 0.3456 | 0.2134 | 0.1545 | 0.2595 | 0.0001 |
| 0.7 | 1005 | 174 | 0.5646 | 0.3444 | 0.2202 | 0.1528 | 0.2693 | 0.0001 |
| 0.7 | mean_over_draws | 174 | 0.5995 | 0.3565 | 0.2430 | 0.1794 | 0.2910 | 0.0001 |

Full subset memberships — the training PDB identifiers drawn for every (threshold, fold, draw) — are in
`results/review/size_matched_control/memberships.tsv` (75 subsets) and in the deposit.
