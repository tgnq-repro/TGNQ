# weibo.rda

This file contains the **Sina Weibo network panel dataset** used in the empirical analysis of the paper:

*Two-way Homogeneity Pursuit for Quantile Network Vector Autoregression*.

It provides the response panel, network structure, and covariates needed to fit the two-way grouped network quantile (TGNQ) model.

---

## 1. Main objects

After loading `weibo.rda`, you will find the following key objects:

### 1.1 Ymat

- Type: numeric matrix  
- Dimension: N × T  
- Content: user activity panel over time.  
- Each element Ymat[i, t] is the log-transformed posting frequency of user *i* on day *t*, constructed as  
  Yᵢₜ = log(1 + Zᵢₜ), where Zᵢₜ is the raw number of posts.

In the analysis of the paper, N = 355 users and T = 78 days (from 2014‑01‑01 to 2014‑03‑19).

---

### 1.2 Amat

- Type: numeric matrix  
- Dimension: N × N  
- Content: adjacency matrix of the directed follower network.

Amat[i, j] = 1 if user *i* follows user *j*, and 0 otherwise.  
This matrix encodes the observed network structure and is later row-normalized to construct the weight matrix used in the TGNQ model.

---

### 1.3 Xi

- Type: numeric matrix  
- Dimension: N × 8  
- Content: user-level static covariates, including an intercept.  
- Column order:

  1. Intercept  
     - Constant 1 for all users; used to estimate the overall intercept in the regression.

  2. gender  
     - Binary indicator of the user’s gender (e.g. 1 = male, 0 = female).

  3. loc_beijing  
     - Binary indicator: 1 if the user’s registered location is Beijing, 0 otherwise.

  4. loc_shanghai  
     - Binary indicator: 1 if the user’s registered location is Shanghai, 0 otherwise.

  5. desc_len  
     - (Log) length or standardized length of the user’s profile description.

  6. Nweibo  
     - (Log) total number of Weibo posts by the user (cumulative post count).

  7. public_account  
     - Binary indicator for whether the account is a verified/public account (e.g. 1 = verified/public, 0 = ordinary).

  8. An additional user-level feature used in the empirical analysis (see the column name in `Xi` for the exact definition).

These covariates are time-invariant and represent baseline user characteristics.

---

### 1.4 Xt

- Type: numeric matrix  
- Dimension: T × 2  
- Content: day-level covariates (common to all users), capturing calendar effects.  
- Columns:

  1. Weekend  
     - Indicator for Saturday or Sunday (1 if the day is a weekend, 0 otherwise).

  2. Holiday  
     - Indicator for national holidays (1 if the day is a public holiday, 0 otherwise).

These variables are used to model systematic temporal patterns in user activity, such as lower activity on weekends and holidays.

---

## 2. Summary

In total, `weibo.rda` provides:

- A log-transformed activity panel `Ymat` (N × T),
- A directed follower network `Amat` (N × N),
- User-level static covariates `Xi` (N × 8, with the first column being an intercept and the remaining columns encoding gender, location, description length, total posts, public account indicator, etc.),
- Day-level covariates `Xt` (T × 2, Weekend and Holiday indicators).

These objects together form the complete dataset used in the real-data application of the TGNQ model.