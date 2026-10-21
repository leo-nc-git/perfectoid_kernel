
# Approximations of elements of a perfectoid ideal

  

This project is a SageMath script for searching approximations of elements in $(t-1)_\text{pfd}$, that are congruent to elements of the form $t^{1/p^k}-1$. It was written with the help of GitHub copilot and GPT Sol 5.5.

  

The main script is [kernel_elements.sage](kernel_elements.sage). It builds finite-level algebras, solves linear systems coming from the evaluation map, and allows the user to explore compatible root towers and Frobenius jumps interactively.

  

## Overview

  

At a finite level $(M, L)$, the program works with the quotient rings

  

- $B_L = \mathbb{F}_p[a]/(a^{p^L})$

- $A _{L,M}= B_L[Y]/(Y^{p^M})$

  

These should be read as finite truncations of the integral-level objects

  
- $B_L'=\mathbb{Z}[\zeta_{p^{L+1}}] / (\zeta_p-1)$

- $A_{L,M}'=B_L'[t^{1/p^{M}}]/(t-1)$


in the following sense: 
  
- $B_L = B_L'/(\zeta_p-1)$

- $A_{L,M} = A_{L,M}'/(\zeta_p-1)$

  
with the identification

- $a = \zeta_p^{L+1}-1$

- $Y = t^{1/p^M}-1$


In the displayed solutions, $e$ is used to denote $\zeta_p-1$, and more precisely $e^{1/p^k}$ is the class of $\zeta_p^{1/p^{k+1}}-1$.


The program defines an evaluation map $E : A_{L,M} \to B_{L}^N$, where $N = p^M$, by evaluating a polynomial in $Y$ at the points

$(1+u)^i - 1, \quad i = 0, \ldots, N-1,$

where $u = a^{p^{L-M+1}}$ is the residue class of $\zeta_p^{L+1}-1$ in $B$.

  

When taking the colimit of $A_{L,M}'$ in $M$ and $L$ and then $p$-completing at the integral level, this evaluation map is not injective, and its kernel is the perfectoidization of the ideal $(t-1)$. Going modulo $p$ (more precisely modulo $\zeta_p-1$ for computational convenience), we can reconstruct approximations of kernel elements. The process is to take a $p$-th root of $t-1$ that lies in the kernel, find a compatible sequence of $p$-th roots by elements in the kernel, and then use the sharp map to lift the element.

  

Therefore the code searches for finite-level approximations of elements in $(t-1)_\mathrm{pfd}$ that are congruent to expressions of the form $t^{1/p^k}-1$.

  

At a chosen finite level and initial root index $n \ge  1$, the program searches for an element $x_n$ such that

  

- $E(x_n) = 0$

- $x_n = Y^{p^{M-n}} \bmod a^{p^{L-n}}$

The congruence implies $x_n^{p^n}=0$, since the ring has characteristic $p$ and $Y^{p^M}=a^{p^L}=0$. No separate nilpotence equations are added to the initial search. In fact, the representative $f_n\in B_L[Y]$ of degree less than $p^M$ satisfies $f_n^{p^n}=Y^{p^M}$ before quotienting by $Y^{p^M}$.

If the search is successful, the user can inspect the resulting affine family, choose a displayed solution, and then continue by taking Frobenius jumps to larger indices. The user may also increase the finite level $(M, L)$ as needed if the current approximation is too small.

  

The script is designed for interactive exploration: it keeps exact-level computations cached, supports user-controlled enlargement of the finite level, and does not behave like a checkpoint/resume system.

  

---

  

## Mathematical intent

  

In rough terms, the code looks for an element $x_n$ such that:

  

- $E(x_n) = 0$

- $x_n = Y^{p^{M-n}} \bmod a^{p^{L-n}}$ 

  

Once such an element is found, it can be used as a starting point for new Frobenius jumps:

 

- from $x_n$, search for a deeper root $x_{n+s}$ satisfying

- $E(x_{n+s}) = 0$

- $x_{n+s}^{p^s} = x_n$

- $x_{n+s} = Y^{p^{M-(n+s)}} \bmod a^{p^{L-(n+s)}}$

  

The script prints these roots in the readable $(e,t)$ notation, where

  

- $a = e^{1/p^L}$

- $Y = t^{1/p^M} - 1$

  

---

  

## Exact integral powers

Whenever the current solution $`x_n`$ is displayed, the program also automatically
computes and displays $`\widetilde{x}_n^{p^n}`$ in characteristic zero, including
when displaying a cached tower or browsing other solutions. Earlier roots in
the tower retain their existing $`(e,t)`$ display; the integral power is shown
for the current (deepest) root.

The integral notation is:

- $`T=t^{1/p^M}`$, so $`T^{p^M}=t`$.
- $`u=\zeta_{p^{L+1}}`$, a primitive $`p^{L+1}`$-st root of unity.

For the unique reduced representative

```math
x_n = \sum_{\substack{0 \le i < p^L \\ 0 \le j < p^M}} c_{ij}a^iY^j,
```

choose integer coefficients $`\widetilde{c}_{ij}\in\{0,\ldots,p-1\}`$ and form:

```math
\widetilde{x}_n = \sum_{i,j}\widetilde{c}_{ij}(u-1)^i(T-1)^j
\quad\text{in}\quad
\mathbb{Z}[u][T].
```

Here $`\mathbb{Z}[u]=\mathbb{Z}[U]/(\Phi_{p^{L+1}}(U))`$. The coefficients are
reduced using the cyclotomic relation. The variable $`T`$ remains free,
so terms such as $`T^{p^M}-1=t-1`$ remain visible.

The computation uses the stored $`(a,Y)`$ representative, before any
Frobenius simplifications in characteristic $`p`$ used in the $`(e,t)`$ display.
All binomial cross terms and all powers of $`T`$ are retained in the
integral output. The finite search still uses its original quotient rings.

This is a specified integral lift and a finite-stage sharp-map approximation.
Changing the lift convention can change the exact finite-stage power.
The degree of the power is $`p^n\deg(\widetilde{x}_n)`$ for a nonzero lift,
so the result can have many terms and large integer coefficients at higher
levels. The integral ring is constructed only when needed, and the last
power is kept in memory for repeated display; existing on-disk caches
continue to store the computations in characteristic $`p`$.

---

## Running the script

  

From the project directory:

  

```bash

sage  kernel_elements.sage

```

  

The script is interactive and asks for:

  

1. prime $p$

2. polynomial depth $M$

3. ramification depth $L$

4. initial root index $n \ge  1$

  

Example:

  

```text

Enter a prime p: 2

Enter polynomial depth M: 2

Enter ramification depth L: 2

Enter initial root index n >= 1: 1

```

  

---

  

## Main workflow

  

### 1. Level setup

  

The script builds a `FiniteLevel` object for the chosen $(M, L)$. This creates the finite-level quotient rings and the matrices used in the search.

  

It also initializes the evaluation matrix $E$ and the Frobenius matrix $F$.

  

### 2. Cache check

  

Before doing a fresh search, the script checks whether a previously computed result is already available in the cache directory.

  

If so, the user sees a choice like:

  

-  `[D]` display the cached solution

-  `[F]` rebuild the full affine family

-  `[Q]` quit

  

The cache stores one concrete root tower for the current level, not the full affine family of all solutions at that level.

  

### 3. Direct root search

  

If no cached tower is used, the program attempts to find a direct solution $x_n$ at the chosen finite level. If the requested root depth is incompatible with the current $(M, L)$, it asks the user to enlarge the level.

  

For example, if $n > M$ or $n > L$, the script prints a message like:

  

```text

The search for x_n requires M >= n and L >= n.

```

  

and asks whether to enlarge $M$, $L$, or both.

  

### 4. Affine family exploration

  

When a direct family is found, it builds an `AffineSolutionFamily` object. This stores:

  

- one particular solution

- the kernel directions, computed only when requested

- a way to generate other solutions by index

  

The user then sees a menu such as:

  

-  `[U]` use the displayed solution

-  `[N]` display the next solution

-  `[R]` display a random non-particular solution

-  `[C]` compute the affine dimension and count the solutions

-  `[Q]` stop without selecting a solution

The first solution can be displayed and accepted without computing the kernel basis. The `[N]`, `[R]`, and `[C]` actions compute that basis when needed. Each displayed root tower is cached; `[Q]` stops exploration while preserving the last displayed tower for a future run.

### 5. Frobenius jumps

  

After a root tower is selected, the program enters the main interactive loop and asks for a jump depth:

  

```text

Enter a Frobenius jump depth n >= 1 from x_k ([Q] to stop):

```

  

The program then searches for a deeper root compatible with the current tower and the needed Frobenius condition.

  

---

  

## Cache behavior

  

The script uses a cache directory named `.sage-cache` by default.

Set `KERNEL_ELEMENTS_CACHE_DIR` to use a different local cache directory. Cache filenames and Sage object formats remain compatible with earlier runs.

It stores:

  

- cached matrices for $E$ and $F$

- cached Frobenius powers such as $F^r$

- cached root towers at exact finite levels

  

Examples of stored objects include files such as:

  

-  `E_p2_M2_L2_v1.sobj`

-  `F_p2_M2_L2_v1.sobj`

-  `Fpower_p2_M2_L2_d2_v1.sobj`

-  `roots_p2_M2_L2.sobj`

  

The cache is exact-level aware: it is keyed by $(p, M, L)$ and by the relevant matrix or Frobenius power.

  

It is a reuse mechanism for exact finite-level computations.

Matrix entries are checked for the correct type, coefficient field, and dimensions when loaded. Cached root towers are checked for consecutive indices starting at 1, the correct finite-level ring, evaluation-kernel membership, congruences, and Frobenius compatibility. Invalid or unreadable entries are ignored and can be recomputed.

Cache writes use a temporary file in the same directory and atomically replace the previous entry after the write completes. Frobenius iterates are constructed directly from their action on monomials rather than by matrix exponentiation.

The repository also includes a read-only [public-cache](public-cache/README.md) of large evaluation matrices, totaling about 583 MiB. Matrix loading tries the local cache first, then the public cache, and finally computes a missing matrix into the local cache. Root towers remain local. Public objects larger than 32 MiB are stored in chunks, assembled in temporary space when needed, and checked against the sizes, checksums, and matrix metadata in the catalog. The public directory is resolved relative to the script; `KERNEL_ELEMENTS_PUBLIC_CACHE_DIR` overrides its location.

To validate and export large matrices already computed in `.sage-cache/`, run `sage tools/export_public_cache.py`. This is an explicit maintenance step; interactive searches never write to the public cache. See the public-cache README for selection rules, options, and validation limits.

---



## Level enlargement and retry logic

  

If the current finite level is too small for the requested root, the script offers to enlarge it.

  

Available actions are:

  

-  `[M]` increase $M$

-  `[L]` increase $L$

-  `[B]` increase both $M$ and $L$

-  `[Q]` stop

  

After enlargement, the previous tower is lifted to the new level and verified before continuing the search.

  
The user can dynamically increase the finite approximation until a root becomes visible.

The initial search keeps the requested root index after enlargement and retries it at the new level.

  

---

  

## Important notes

  
### The root tower is a coherent sequence

  

When a solution is found, the program reconstructs the earlier roots by applying the Frobenius.

  

### Available choices during execution

  

Even when a cached tower exists, the script still asks whether to:

  

- display it,

- rebuild the full affine family,

- or quit.

  

This keeps the user in charge of whether they want the cheap cached view or a full fresh family reconstruction.

  

---

  

## Typical session flow

  

A typical run looks like this:

  

1. Enter $(p, M, L, n)$

2. Check cache for exact-level roots

3. If cache exists, choose display or rebuild

4. Otherwise search directly for $x_n$

5. If level is too small, enlarge $(M, L)$ and retry

6. Pick a family solution or explore other solutions

7. Save the root tower to cache

8. Ask for a Frobenius jump depth and continue

  

---

  

## Files in this project

  

- [kernel_elements.sage](kernel_elements.sage) — main Sage script

- [tests/test_kernel_elements.py](tests/test_kernel_elements.py) — mathematical and interactive regression checks

- [public-cache/](public-cache/README.md) — shared large evaluation matrices and their checksummed catalog

- [tools/export_public_cache.py](tools/export_public_cache.py) — validate and export local matrix precomputations

-  `.sage-cache/` — cache directory for exact-level matrices and root towers

  

---

  

## Practical advice

  

- Start with small values of $p$, $M$, and $L$.
- The situation is exponential. With usual computers memory will crash almost certainly once the computation are too large. 
- Keep the cache directory in mind when experimenting: reusing exact-level matrices can save substantial time.
  

---

## Regression checks

Run the tests with SageMath:

```bash
sage tests/test_kernel_elements.py
```

The tests use temporary caches. They check matrix actions against ring operations, direct Frobenius iterates against matrix powers, exact polynomial lifts, successful compatible jumps, level transport, local and public cache validation, chunk assembly and checksums, interrupted writes, the export tool, and the interactive startup and solution-selection flows.


