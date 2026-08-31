
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

  

The program defines an evaluation map $E : A \to B^N$, where $N = p^M$, by evaluating a polynomial in $Y$ at the points

$(1+u)^i - 1, \quad i = 0, \ldots, N-1,$

where $u = a^{p^{L-M+1}}$ is the residue class of $\zeta_p^{L+1}-1$ in $B$.

  

When taking the colimit of $A_{L,M}'$ in $M$ and $L$ and then $p$-completing at the integral level, this evaluation map is not injective, and its kernel is the perfectoidization of the ideal $(t-1)$. Going modulo $p$ (more precisely modulo $\zeta_p-1$ for computational convenience), we can reconstruct approximations of kernel elements. The process is to take a $p$-th root of $t-1$ that lies in the kernel, find a compatible sequence of $p$-th roots by elements in the kernel, and then use the sharp map to lift the element.

  

Therefore the code searches for finite-level approximations of elements in $(t-1)_\mathrm{pfd}$ that are congruent to expressions of the form $t^{1/p^k}-1$.

  

At a chosen finite level and initial root index $n \ge  1$, the program searches for an element $x_n$ such that

  

- $E(x_n) = 0$

- $x_n^{p^n} = 0$

- $x_n = Y^{p^{M-n}} \bmod a^{p^{L-n}}$

  

If the search is successful, the user can inspect the resulting affine family, choose a displayed solution, and then continue by taking Frobenius jumps to larger indices. The user may also increase the finite level $(M, L)$ as needed if the current approximation is too small.

  

The script is designed for interactive exploration: it keeps exact-level computations cached, supports user-controlled enlargement of the finite level, and does not behave like a checkpoint/resume system.

  

---

  

## Mathematical intent

  

In rough terms, the code looks for an element $x_n$ such that:

  

- $E(x_n) = 0$

- $x_n^{p^n} = 0$ 

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

- the kernel directions

- a way to generate other solutions by index

  

The user then sees a menu such as:

  

-  `[U]` use the displayed solution

-  `[N]` display the next solution

-  `[R]` display a random non-particular solution

-  `[Q]` stop without selecting a solution

  

### 5. Frobenius jumps

  

After a root tower is selected, the program enters the main interactive loop and asks for a jump depth:

  

```text

Enter a Frobenius jump depth n >= 1 from x_k ([Q] to stop):

```

  

The program then searches for a deeper root compatible with the current tower and the needed Frobenius condition.

  

---

  

## Cache behavior

  

The script uses a cache directory named `.sage-cache` by default.

  

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

  

---

  

## Important notes

  

### Cache is not a checkpoint

  

The script does not save a full session state and resume later from it. Instead, it saves exact-level computations and root towers so that repeated work can be avoided.
  

### The root tower is a coherent sequence

  

When a solution is found, the program reconstructs the earlier roots by applying the Frobenius.

  

### The user remains in control

  

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

-  `.sage-cache/` — cache directory for exact-level matrices and root towers

  

---

  

## Practical advice

  

- Start with small values of $p$, $M$, and $L$.
- The situation is exponential. With usual computers memory will crash almost certainly once the computation are too large. 
- Keep the cache directory in mind when experimenting: reusing exact-level matrices can save substantial time.
  

---


