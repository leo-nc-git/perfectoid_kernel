"""
Search for kernel elements at finite levels.

At level (M, L):

B = F_p[a]/(a^(p^L)).
This can be seen as Z[zeta_p^{L+1}]/(zeta_p-1), where a is the residue class of zeta_p^{L+1}-1. L is the "ramification depth".

In displayed solutions, e is used to denote zeta_p-1. More precisely, e^{1/p^k} is the class of zeta_p^{1/p^{k+1}}-1.

A = B[Y]/(Y^(p^M)).
This can be seen as Z[zeta_p^{L+1}][t^{1/p^M}]/(t-1, zeta_p-1), where Y is the residue class of t^{1/p^M}-1. M is the "polynomial depth".

The program defines an evaluation map E: A -> B^N, where N = p^M, by evaluating a polynomial in Y at the points
    (1 + u)^i - 1, for i = 0, ..., N-1,
where u = a^(p^(L-M+1)) is the residue class of zeta_p^{L+1}-1 in B.

When taking the colimit in M and L and then p-completing at the integral level, this evaluation map is not injective, and its kernel is the perfectoidization of the ideal (t-1).
Going modulo p (actually modulo zeta_p-1 for computational convenience), we can in principle reconstruct approximations of the kernel elements.

The process to reach an exact element of the kernel which is t-1 modulo p goes as follows:
- take a p-th root of t-1 that lies in the kernel;
- find a compatible sequence of p-th roots of this element by elements in the kernel;
- use the sharp map to lift the element. The sharp map sends a compatible Frobenius sequence (a_nbar) to lim_n a_n^{p^n}, where the limit is taken in the p-adic topology.

Therefore the code can find approximations of such elements.

The user can increase the finite level (M, L) as needed to find a solution.

When the user inputs n >= 1, the search looks for an element x_n such that

    E(x_n) = 0, (In the kernel)
    x_n = Y^(p^(M-n)) mod (a^(p^(L-n))), (The correct root level.)

The congruence implies x_n^(p^n) = 0 in characteristic p, since
Y^(p^M) = a^(p^L) = 0, so nilpotence needs no separate equation.

If the search is successful, the user sees a solution and may display other ones if desired. From this solution, one can then "reuse" it.

Inputting a new jump depth s allows the search for a new x_(s+n) from a previous x_n by imposing the following conditions on x_(s+n):

    E(x_(s+n)) = 0,
    x_(s+n)^(p^n) = x_s,
    x_(s+n) = Y^(p^(M-s-n)) mod (a^(p^(L-s-n))).
"""

import os
from random import randrange


class AffineSolutionFamily:
    """A particular solution plus the right kernel of a linear system."""

    def __init__(self, level, system_matrix, particular_vector, validator):
        self.level = level
        self.system_matrix = system_matrix
        self.particular_vector = particular_vector
        self.validator = validator
        self._kernel_basis = None
        self._last_displayed = None
        self._last_displayed_index = None

    def kernel_basis(self):
        """Compute the homogeneous solution space only when it is requested."""
        if self._kernel_basis is None:
            print("Computing the affine solution directions...")
            kernel = self.system_matrix.right_kernel()
            self._kernel_basis = list(kernel.basis())
        return self._kernel_basis

    def dimension(self):
        return Integer(len(self.kernel_basis()))

    def number_of_solutions(self):
        return self.level.p ** self.dimension()

    def element(self, index=0):
        """
        Return solution number index without constructing all solutions.

        The base-p digits of index are the coefficients of the kernel basis.
        Index 0 is the particular solution returned by solve_right().
        """
        index = Integer(index)
        if index < 0:
            raise ValueError("The solution index must be nonnegative.")

        solution_vector = vector(self.level.F, self.particular_vector)

        if index != 0:
            basis = self.kernel_basis()
            total = self.level.p ** len(basis)
            if index >= total:
                raise IndexError("Solution index is outside the affine family.")

            remaining = index
            for direction in basis:
                coefficient = self.level.F(remaining % self.level.p)
                if coefficient != 0:
                    solution_vector += coefficient * direction
                remaining //= self.level.p

        x = self.level.vector_to_A(solution_vector)
        self.validator(x)
        return x

    def random_nonparticular_index(self):
        total = self.number_of_solutions()
        if total == 1:
            return None
        return Integer(randrange(1, int(total)))


class FiniteLevel:
    """All objects and matrices depending on one pair (M, L)."""

    matrix_cache_version = 1

    def __init__(self, p, M, L):
        self.p = Integer(p)
        self.M = Integer(M)
        self.L = Integer(L)

        if not self.p.is_prime():
            raise ValueError("p must be prime.")
        if self.M < 1:
            raise ValueError("Need M >= 1.")
        if self.L < 0:
            raise ValueError("Need L >= 0.")
        if self.L < self.M - 1:
            raise ValueError("Need L >= M-1 so that the evaluation parameter is visible.")

        self.F = GF(self.p)

        # B = F_p[a]/(a^(p^L))
        self.R = PolynomialRing(self.F, "a_polynomial")
        a_polynomial = self.R.gen()
        self.B = self.R.quotient(a_polynomial ** (self.p ** self.L))
        self.a = self.B.gen()

        # A = B[Y]/(Y^(p^M))
        self.S = PolynomialRing(self.B, "Y_polynomial")
        Y_polynomial = self.S.gen()
        self.N = self.p ** self.M
        self.A = self.S.quotient(Y_polynomial ** self.N)
        self.Y = self.A.gen()

        self.dim_B = self.p ** self.L
        self.dim_A = self.dim_B * self.N

        # u = epsilon^(1/p^(M-1))
        eval_exp = self.p ** (self.L - self.M + 1)
        self.u = self.a ** eval_exp
        self.evaluation_points = [
            (1 + self.u) ** i - 1 for i in range(self.N)
        ]

        self.basis_A = [
            self.A(self.a ** i) * self.Y ** j
            for j in range(self.N)
            for i in range(self.dim_B)
        ]

        print("\nBuilding finite level (M, L) =", (self.M, self.L))
        print("Dimension over F_p:", self.dim_A)
        self._build_matrices()

    def E(self, x):
        polynomial = x.lift()
        return [polynomial(point) for point in self.evaluation_points]

    def B_to_vector(self, b):
        polynomial = b.lift()
        return vector(self.F, [polynomial[i] for i in range(self.dim_B)])

    def A_to_vector(self, x):
        polynomial = x.lift()
        entries = []

        # Basis order: first Y-degree j, then a-degree i.
        for j in range(self.N):
            entries.extend(self.B_to_vector(polynomial[j]))

        return vector(self.F, entries)

    def vector_to_A(self, v):
        x = self.A(0)
        index = 0

        for j in range(self.N):
            for i in range(self.dim_B):
                coefficient = v[index]
                if coefficient != 0:
                    x += self.A(coefficient * self.a ** i) * self.Y ** j
                index += 1

        return x

    def E_to_vector(self, x):
        entries = []
        for value in self.E(x):
            entries.extend(self.B_to_vector(value))
        return vector(self.F, entries)

    def _matrix_cache_path(self, name):
        cache_directory = os.environ.get("KERNEL_ELEMENTS_CACHE_DIR", ".sage-cache")
        if not os.path.isdir(cache_directory):
            os.makedirs(cache_directory)
        filename = "%s_p%s_M%s_L%s_v%s.sobj" % (name, self.p, self.M, self.L, self.matrix_cache_version)
        return os.path.join(cache_directory, filename)

    def _load_matrix(self, name):
        path = self._matrix_cache_path(name)
        if not os.path.exists(path):
            return None

        try:
            cached_matrix = load(path)
        except (OSError, EOFError, ValueError):
            print("Ignoring unreadable cached matrix %s." % path)
            return None

        if cached_matrix.nrows() != self.dim_A or cached_matrix.ncols() != self.dim_A:
            print("Ignoring cached matrix with the wrong dimensions: %s" % path)
            return None
        return cached_matrix

    def _save_matrix(self, name, matrix_to_save):
        save(matrix_to_save, self._matrix_cache_path(name))

    def _iterated_frobenius_cache_path(self, depth):
        cache_directory = os.environ.get("KERNEL_ELEMENTS_CACHE_DIR", ".sage-cache")
        if not os.path.isdir(cache_directory):
            os.makedirs(cache_directory)
        filename = "Fpower_p%s_M%s_L%s_d%s_v%s.sobj" % (self.p, self.M, self.L, depth, self.matrix_cache_version)
        return os.path.join(cache_directory, filename)

    def _load_frobenius_power(self, depth):
        path = self._iterated_frobenius_cache_path(depth)
        if not os.path.exists(path):
            return None

        try:
            cached_power = load(path)
        except (OSError, EOFError, ValueError):
            print("Ignoring unreadable cached Frobenius iterate at %s." % path)
            return None

        if cached_power.nrows() != self.dim_A or cached_power.ncols() != self.dim_A:
            print("Ignoring cached Frobenius iterate with the wrong dimensions: %s" % path)
            return None
        return cached_power

    def _save_frobenius_power(self, depth, matrix_to_save):
        save(matrix_to_save, self._iterated_frobenius_cache_path(depth))

    def _build_E_matrix(self):
        """Build E from coefficient vectors, avoiding quotient-ring evaluations."""
        entries = {}

        for point_index, point in enumerate(self.evaluation_points):
            point_power = self.B(1)
            for Y_degree in range(self.N):
                point_coefficients = point_power.lift()

                for a_degree in range(self.dim_B):
                    column = Y_degree * self.dim_B + a_degree
                    for output_degree in range(a_degree, self.dim_B):
                        coefficient = point_coefficients[output_degree - a_degree]
                        if coefficient != 0:
                            row = point_index * self.dim_B + output_degree
                            entries[(row, column)] = coefficient

                point_power *= point

        return matrix(self.F, self.N * self.dim_B, self.N * self.dim_B, entries, sparse=True)

    def _build_F_matrix(self):
        """Build Frobenius directly from its action on monomial basis elements."""
        entries = {}

        for Y_degree in range(self.N):
            for a_degree in range(self.dim_B):
                column = Y_degree * self.dim_B + a_degree
                frobenius_a_degree = a_degree * self.p
                frobenius_Y_degree = Y_degree * self.p

                if frobenius_a_degree < self.dim_B and frobenius_Y_degree < self.N:
                    row = frobenius_Y_degree * self.dim_B + frobenius_a_degree
                    entries[(row, column)] = 1

        return matrix(self.F, self.dim_A, self.dim_A, entries, sparse=True)

    def _build_matrices(self):
        self.E_matrix = self._load_matrix("E")
        if self.E_matrix is None:
            print("Building matrix of E...", flush=True)
            self.E_matrix = self._build_E_matrix()
            self._save_matrix("E", self.E_matrix)
            print("Finished building matrix of E.", flush=True)
        else:
            print("Loaded cached matrix of E.", flush=True)

        self.F_matrix = self._load_matrix("F")
        if self.F_matrix is None:
            print("Building Frobenius matrix...", flush=True)
            self.F_matrix = self._build_F_matrix()
            self._save_matrix("F", self.F_matrix)
            print("Finished building Frobenius matrix.", flush=True)
        else:
            print("Loaded cached Frobenius matrix.", flush=True)

        self._reduction_matrices = {}
        self._frobenius_powers = {1: self.F_matrix}

    def iterated_frobenius_matrix(self, depth):
        """Matrix of x |-> x^(p^depth), cached by depth."""
        depth = Integer(depth)
        if depth < 0:
            raise ValueError("Frobenius depth must be nonnegative.")
        if depth == 0:
            return identity_matrix(self.F, self.dim_A, sparse=True)

        key = int(depth)
        if key in self._frobenius_powers:
            return self._frobenius_powers[key]

        cached_power = self._load_frobenius_power(depth)
        if cached_power is not None:
            print("Loaded cached Frobenius iterate F^%s." % depth)
            self._frobenius_powers[key] = cached_power
            return cached_power

        print("Building Frobenius iterate F^%s..." % depth)
        self._frobenius_powers[key] = self.F_matrix ** depth
        self._save_frobenius_power(depth, self._frobenius_powers[key])
        return self._frobenius_powers[key]

    def reduction_power(self, r):
        """Return q_r = p^(L-r), so a^q_r = e^(1/p^r)."""
        r = Integer(r)
        if r < 1 or r > self.M:
            raise ValueError("Need 1 <= r <= M.")
        if self.L < r:
            raise ValueError("Need L >= r for the congruence modulus.")
        return self.p ** (self.L - r)

    def reduction_matrix(self, r):
        """
        Matrix of
            A -> A/(a^q_r),  q_r = p^(L-r).

        For every Y-degree j, it keeps the coefficients of
            1, a, ..., a^(q_r-1).
        """
        r = Integer(r)
        key = int(r)
        if key in self._reduction_matrices:
            return self._reduction_matrices[key]

        q = self.reduction_power(r)
        Rmat = matrix(self.F, self.N * q, self.dim_A, sparse=True)

        for j in range(self.N):
            for i in range(q):
                row = j * q + i
                column = j * self.dim_B + i
                Rmat[row, column] = 1

        self._reduction_matrices[key] = Rmat
        return Rmat

    def reduction_target(self, r):
        r = Integer(r)
        q = self.reduction_power(r)
        target = vector(self.F, self.N * q)

        # The coefficient of a^0 * Y^(p^(M-r)) is 1.
        Y_degree = self.p ** (self.M - r)
        target[Y_degree * q] = 1
        return target

    def zero_root_family(self, depth):
        """
        Return the affine family of x_depth satisfying

            E(x_depth) = 0,
            x_depth = Y^(p^(M-depth)) mod a^(p^(L-depth)).

        The congruence implies x_depth^(p^depth) = 0; the validator
        checks this consequence without adding redundant linear equations.
        """
        depth = Integer(depth)
        if depth < 1:
            raise ValueError("The root index must be at least 1.")
        if self.M < depth or self.L < depth:
            return None

        reduction_matrix = self.reduction_matrix(depth)
        rhs_E = vector(self.F, self.E_matrix.nrows())
        rhs_congruence = self.reduction_target(depth)
        rhs = vector(self.F, list(rhs_E) + list(rhs_congruence))

        system_matrix = self.E_matrix.stack(reduction_matrix)

        try:
            solution = system_matrix.solve_right(rhs)
        except ValueError:
            return None

        def validator(x):
            assert all(value == 0 for value in self.E(x))
            assert x ** (self.p ** depth) == 0
            assert reduction_matrix * self.A_to_vector(x) == rhs_congruence

        return AffineSolutionFamily(self, system_matrix, solution, validator)

    def initialization_family(self):
        """Backward-compatible alias for the x_1 startup case."""
        return self.zero_root_family(1)

    def direct_root_family(self, x_current, current_index, depth):
        """
        Return the affine family for a direct jump from x_s to x_(s+n):

            E(x_(s+n)) = 0,
            x_(s+n)^(p^n) = x_s,
            x_(s+n) = Y^(p^(M-s-n)) mod a^(p^(L-s-n)).
        """
        current_index = Integer(current_index)
        depth = Integer(depth)
        if current_index < 1:
            raise ValueError("The current root index must be at least 1.")
        if depth < 1:
            raise ValueError("Direct Frobenius depth must be at least 1.")

        target_index = current_index + depth
        if target_index > self.M or target_index > self.L:
            return None

        iterated_F = self.iterated_frobenius_matrix(depth)
        reduction_matrix = self.reduction_matrix(target_index)
        rhs_E = vector(self.F, self.E_matrix.nrows())
        rhs_F = self.A_to_vector(x_current)
        rhs_congruence = self.reduction_target(target_index)
        rhs = vector(self.F, list(rhs_E) + list(rhs_F) + list(rhs_congruence))

        system_matrix = self.E_matrix.stack(iterated_F)
        system_matrix = system_matrix.stack(reduction_matrix)

        try:
            solution = system_matrix.solve_right(rhs)
        except ValueError:
            return None

        def validator(x):
            assert all(value == 0 for value in self.E(x))
            assert x ** (self.p ** depth) == x_current
            assert reduction_matrix * self.A_to_vector(x) == rhs_congruence

        return AffineSolutionFamily(self, system_matrix, solution, validator)


def lift_element(x, old_level, new_level):
    """
    Transport x from (M,L) to a larger (M',L') by
        a_old |-> a_new^(p^(L'-L)),
        Y_old |-> Y_new^(p^(M'-M)).
    """
    if old_level.p != new_level.p:
        raise ValueError("Cannot change p while lifting elements.")
    if new_level.M < old_level.M or new_level.L < old_level.L:
        raise ValueError("The new values M and L must not decrease.")

    scale_a = old_level.p ** (new_level.L - old_level.L)
    scale_Y = old_level.p ** (new_level.M - old_level.M)
    old_Y_polynomial = x.lift()
    result = new_level.A(0)

    for j in range(old_level.N):
        old_a_polynomial = old_Y_polynomial[j].lift()
        for i in range(old_level.dim_B):
            coefficient = old_a_polynomial[i]
            if coefficient != 0:
                new_coefficient = new_level.F(coefficient)
                result += new_level.A(new_coefficient * new_level.a ** (i * scale_a)) * new_level.Y ** (j * scale_Y)

    return result


def verify_tower(level, roots):
    """Verify every condition again, especially after changing level."""
    for r in sorted(roots):
        x = roots[r]
        assert all(value == 0 for value in level.E(x))
        reduction_matrix = level.reduction_matrix(r)
        assert reduction_matrix * level.A_to_vector(x) == level.reduction_target(r)
        if r == 1:
            assert x ** level.p == 0
        else:
            assert x ** level.p == roots[r - 1]


def element_in_t_e(level, x):
    """
    Return a readable translation of x from the variables (a,Y) to (e,t):

        a = e^(1/p^L),
        Y = t^(1/p^M)-1.

    Powers of Y are simplified using Frobenius.  If j = p^v*m with
    p not dividing m, then

        (t^(1/p^M)-1)^j = (t^(1/p^(M-v))-1)^m.

    The result is a display string. Computations still take place in A.
    """
    Y_polynomial = x.lift()
    terms = []

    for j in range(level.N):
        a_polynomial = Y_polynomial[j].lift()

        for i in range(level.dim_B):
            coefficient = a_polynomial[i]
            if coefficient == 0:
                continue

            factors = []

            # Omit a coefficient equal to 1 unless it is the constant term.
            if coefficient != 1 or (i == 0 and j == 0):
                factors.append(str(coefficient))

            if i != 0:
                e_exponent = QQ(i) / (level.p ** level.L)
                factors.append("e^(%s)" % e_exponent)

            if j != 0:
                # Extract the largest Frobenius power p^v dividing j.
                remaining_power = Integer(j)
                frobenius_valuation = Integer(0)
                while remaining_power % level.p == 0:
                    remaining_power //= level.p
                    frobenius_valuation += 1

                simplified_M = level.M - frobenius_valuation
                t_denominator = level.p ** simplified_M
                t_expression = "t^(1/%s)-1" % t_denominator

                if remaining_power == 1:
                    # Parentheses are needed only when another factor precedes it.
                    if factors:
                        factors.append("(" + t_expression + ")")
                    else:
                        factors.append(t_expression)
                else:
                    factors.append("(" + t_expression + ")^" + str(remaining_power))

            terms.append("*".join(factors))

    if not terms:
        return "0"

    return " + ".join(terms)


def print_root(level, r, x):
    """Print a root and, if relevant, all previous roots obtained by Frobenius descent."""
    if r > 1:
        print("x_%s in (e,t) =" % r, element_in_t_e(level, x))
        print("Previous roots obtained by Frobenius descent:")
        for lower_index in range(r - 1, 0, -1):
            lower_root = x ** (level.p ** (r - lower_index))
            print("x_%s in (e,t) =" % lower_index, element_in_t_e(level, lower_root))
    else:
        print("x_%s in (e,t) =" % r, element_in_t_e(level, x))


def save_displayed_root_tower(level, family, target_index, x_target):
    """Persist the solution family currently displayed to the user, even if they quit after inspection."""
    roots = {target_index: x_target}
    for r in range(target_index, 1, -1):
        roots[r - 1] = roots[r] ** level.p

    verify_tower(level, roots)
    save_level_roots(level, roots)
    print("Saved displayed root tower for (p, M, L) = (%s, %s, %s)." % (level.p, level.M, level.L))
    return roots


def choose_solution(level, family, r):
    """Interactively inspect and choose an element of an affine family."""
    current_index = Integer(0)
    next_index = Integer(1)
    current = family.element(current_index)

    total = family.number_of_solutions()
    print("\nA solution was found at ramification depth %s and polynomial depth %s." % (level.L, level.M))
    print("Number of affine solutions at this level:", total)
    print("Solution index:", current_index)
    print_root(level, r, current)

    roots = {r: current}
    for lower_index in range(r, 1, -1):
        roots[lower_index - 1] = roots[lower_index] ** level.p
    verify_tower(level, roots)
    save_level_roots(level, roots)

    while True:
        print("\nChoose what to do with this affine solution family:")
        print("  [U] use the displayed solution")
        print("  [N] display the next solution")
        print("  [R] display a random non-particular solution")
        print("  [Q] stop without selecting a solution")
        choice = input("> ").strip().upper()

        family._last_displayed = current
        family._last_displayed_index = current_index

        if choice == "U":
            return current

        if choice == "Q":
            return None

        if choice not in ("N", "R"):
            print("Unknown choice.")
            continue

        dimension = family.dimension()
        total = family.number_of_solutions()
        print("Affine solution dimension:", dimension)
        print("Number of solutions:", total)

        if total == 1:
            print("The displayed solution is unique at this finite level.")
            continue

        if choice == "N":
            if next_index >= total:
                print("Every solution has already been traversed.")
                continue
            current_index = next_index
            next_index += 1
        else:
            current_index = family.random_nonparticular_index()

        current = family.element(current_index)
        print("\nSolution index:", current_index)
        print_root(level, r, current)


def extend_tower_from_direct_root(level, roots, x, current_index, target_index):
    """Add a direct root and recover every skipped intermediate root."""
    current_index = Integer(current_index)
    target_index = Integer(target_index)
    extended_roots = dict(roots)
    extended_roots[target_index] = x

    for r in range(target_index, current_index, -1):
        extended_roots[r - 1] = extended_roots[r] ** level.p

    assert extended_roots[current_index] == roots[current_index]
    verify_tower(level, extended_roots)
    return extended_roots


def read_increment(name):
    answer = input("Increase " + name + " by [default 1]: ").strip()
    increment = Integer(answer) if answer else Integer(1)
    if increment < 1:
        raise ValueError("The increment must be positive.")
    return increment


def choose_larger_level(level, r):
    """Ask how to enlarge the finite level. Return (M,L), or None."""
    while True:
        print("\nChoose an action:")
        print("  [M] increase M")
        print("  [L] increase L")
        print("  [B] increase both M and L")
        print("  [Q] stop")
        choice = input("> ").strip().upper()

        new_M = level.M
        new_L = level.L

        try:
            if choice == "M":
                new_M += read_increment("M")
            elif choice == "L":
                if r > level.M:
                    print("Increasing L alone cannot make x_%s visible." % r)
                    continue
                new_L += read_increment("L")
            elif choice == "B":
                new_M += read_increment("M")
                new_L += read_increment("L")
            elif choice == "Q":
                return None
            else:
                print("Unknown choice.")
                continue
        except ValueError as error:
            print(error)
            continue

        if new_L < new_M - 1:
            print("This choice violates L >= M-1; increase L as well.")
            continue

        return (new_M, new_L)


def print_tower(level, roots):
    print("\nCompatible kernel tower at (M, L) =", (level.M, level.L))
    print("Here a = e^(1/%s) and Y = t^(1/%s)-1." % (level.p ** level.L, level.p ** level.M))
    for r in sorted(roots, reverse=True):
        print()
        if r == max(roots):
            print_root(level, r, roots[r])
        else:
            print("x_%s in (e,t) =" % r, element_in_t_e(level, roots[r]))


def rebuild_larger_level(level, roots, required_r):
    """Increase (M,L), transport the current tower, and verify it."""
    new_parameters = choose_larger_level(level, required_r)
    if new_parameters is None:
        return None

    new_M, new_L = new_parameters
    old_level = level
    new_level = FiniteLevel(level.p, new_M, new_L)
    new_roots = {
        index: lift_element(value, old_level, new_level)
        for index, value in roots.items()
    }
    verify_tower(new_level, new_roots)
    print("Existing roots transported and verified.")
    return (new_level, new_roots)


def initialization_search(level):
    """Find and select x_1 before any direct root searches."""
    roots = {}

    while True:
        if level.L < 1:
            print("\nInitialization requires L >= 1.")
            family = None
        else:
            q = level.reduction_power(1)
            print("\nSearching for x_1...")
            print("Imposing E(x_1) = 0 and x_1^p = 0.")
            print("Imposing x_1 = Y^(p^(M-1)) mod a^%s = e^(1/p)." % q)
            family = level.initialization_family()

        if family is not None:
            x_1 = choose_solution(level, family, 1)
            if x_1 is None:
                if family._last_displayed is not None:
                    roots = {1: family._last_displayed}
                    verify_tower(level, roots)
                    save_level_roots(level, roots)
                    return (level, roots, True)
                return (level, roots, False)

            roots = {1: x_1}
            verify_tower(level, roots)
            return (level, roots, True)

        print("No initialization solution at the current finite level.")
        rebuilt = rebuild_larger_level(level, roots, 1)
        if rebuilt is None:
            return (level, roots, False)
        level, roots = rebuilt
        # Retry initialization at the new level.


def read_jump_depth(current_index):
    """Read n for a direct jump from x_s to x_(s+n), or return None."""
    while True:
        answer = input("\nEnter a Frobenius jump depth n >= 1 from x_%s ([Q] to stop): " % current_index).strip()
        if answer.upper() == "Q":
            return None

        try:
            depth = Integer(answer)
        except (TypeError, ValueError):
            print("Please enter an integer.")
            continue

        if depth < 1:
            print("Need n >= 1.")
            continue
        return depth


def direct_jump_search(level, roots, current_index, depth):
    """Search directly from x_s to x_(s+n) at a user-selected depth n."""
    current_index = Integer(current_index)
    depth = Integer(depth)
    target_index = current_index + depth

    while True:
        if target_index > level.M or target_index > level.L:
            print("\nThe jump from x_%s to x_%s requires M >= %s and L >= %s." % (current_index, target_index, target_index, target_index))
            family = None
        else:
            q = level.reduction_power(target_index)
            print("\nSearching directly for x_%s with x_%s^(p^%s) = x_%s..." % (target_index, target_index, depth, current_index))
            print("Imposing x_%s = Y^(p^(M-%s)) mod a^%s = e^(1/p^%s)." % (target_index, target_index, q, target_index))
            family = level.direct_root_family(roots[current_index], current_index, depth)

        if family is not None:
            x = choose_solution(level, family, target_index)
            if x is None:
                if family._last_displayed is not None:
                    roots = extend_tower_from_direct_root(level, roots, family._last_displayed, current_index, target_index)
                    print("\nThe displayed root was saved for the exact level cache.")
                    save_level_roots(level, roots)
                    return (level, roots, True)
                return (level, roots, False)

            roots = extend_tower_from_direct_root(level, roots, x, current_index, target_index)
            print("\nThe skipped intermediate roots were recovered by Frobenius.")
            return (level, roots, True)

        print("No direct solution at the current finite level.")
        rebuilt = rebuild_larger_level(level, roots, target_index)
        if rebuilt is None:
            return (level, roots, False)
        level, roots = rebuilt
        # Retry the same jump at the new level.


def startup_search(level, requested_n=None):
    """Search for an initial root x_n directly, with n >= 1. If n == 1, this agrees with the old initialization."""
    cached_roots = load_level_roots(level)

    while True:
        if requested_n is None:
            raw = input("\nEnter initial root index n >= 1 ([Q] to stop): ").strip()
            if raw.upper() == "Q":
                return (level, {}, False)

            try:
                n = Integer(raw)
            except (TypeError, ValueError):
                print("Please enter a positive integer.")
                continue
        else:
            n = requested_n
            requested_n = None

        if n < 1:
            print("Need n >= 1.")
            continue

        if cached_roots is not None and n > 1 and n in cached_roots:
            roots = {idx: cached_roots[idx] for idx in sorted(cached_roots) if idx <= n}
            while True:
                choice = input("\nA cached solution exists for x_%s. [D] display it, [F] rebuild the full affine family, or [Q] quit: " % n).strip().upper()
                if choice == "D":
                    verify_tower(level, roots)
                    return (level, roots, True)
                if choice == "F":
                    break
                if choice == "Q":
                    return (level, {}, False)
                print("Please choose [D], [F], or [Q].")

        if n == 1:
            return initialization_search(level)

        if n > level.M or n > level.L:
            print("\nThe search for x_%s requires M >= %s and L >= %s." % (n, n, n))
            rebuilt = rebuild_larger_level(level, {}, n)
            if rebuilt is None:
                return (level, {}, False)
            level, _ = rebuilt
            cached_roots = load_level_roots(level)
            continue

        family = level.zero_root_family(n)
        if family is None:
            print("No direct x_%s solution at the current finite level." % n)
            rebuilt = rebuild_larger_level(level, {}, n)
            if rebuilt is None:
                return (level, {}, False)
            level, _ = rebuilt
            cached_roots = load_level_roots(level)
            continue

        x_n = choose_solution(level, family, n)
        if x_n is None:
            if family._last_displayed is not None:
                roots = {idx: family._last_displayed for idx in [n]}
                for r in range(n, 1, -1):
                    roots[r - 1] = roots[r] ** level.p
                verify_tower(level, roots)
                save_level_roots(level, roots)
                return (level, roots, True)
            return (level, {}, False)

        roots = {n: x_n}
        for r in range(n, 1, -1):
            roots[r - 1] = roots[r] ** level.p
        verify_tower(level, roots)
        return (level, roots, True)


def level_matrix_path(level, name):
    cache_directory = os.environ.get("KERNEL_ELEMENTS_CACHE_DIR", ".sage-cache")
    if not os.path.isdir(cache_directory):
        os.makedirs(cache_directory)
    return os.path.join(cache_directory, "%s_p%s_M%s_L%s_v%s.sobj" % (name, level.p, level.M, level.L, level.matrix_cache_version))


def level_roots_path(level):
    cache_directory = os.environ.get("KERNEL_ELEMENTS_CACHE_DIR", ".sage-cache")
    if not os.path.isdir(cache_directory):
        os.makedirs(cache_directory)
    return os.path.join(cache_directory, "roots_p%s_M%s_L%s.sobj" % (level.p, level.M, level.L))


def save_level_roots(level, roots):
    payload = {Integer(index): root for index, root in roots.items()}
    save(payload, level_roots_path(level))
    print("Saved cached roots for (p, M, L) = (%s, %s, %s) at %s" % (level.p, level.M, level.L, level_roots_path(level)))


def load_level_roots(level):
    path = level_roots_path(level)
    if not os.path.exists(path):
        return None

    try:
        payload = load(path)
    except (OSError, EOFError, ValueError):
        print("Ignoring unreadable cached roots at %s." % path)
        return None

    if not isinstance(payload, dict):
        print("Ignoring malformed cached roots at %s." % path)
        return None

    roots = {Integer(index): root for index, root in payload.items()}
    if not roots:
        return None
    return roots


def interactive_search():
    p = Integer(input("Enter a prime p: "))
    M = Integer(input("Enter polynomial depth M: "))
    L = Integer(input("Enter ramification depth L: "))
    initial_n = Integer(input("Enter initial root index n >= 1: "))

    level = FiniteLevel(p, M, L)
    cached_roots = load_level_roots(level)

    if cached_roots is not None and initial_n in cached_roots:
        print("\nA solution is already in the cache for this level.")
        while True:
            choice = input("Do you want to [D]isplay the cached solution, [F] rebuild the full affine family, or [Q] quit? ").strip().upper()
            if choice == "D":
                roots = {idx: cached_roots[idx] for idx in sorted(cached_roots) if idx <= initial_n}
                current_index = max(roots)
                break
            if choice == "F":
                level, roots, initialized = startup_search(level, initial_n)
                if not initialized:
                    print("\nStartup search stopped before choosing an initial root.")
                    return
                current_index = max(roots)
                save_level_roots(level, roots)
                break
            if choice == "Q":
                return
            print("Please choose [D], [F], or [Q].")
    else:
        level, roots, initialized = startup_search(level, initial_n)
        if not initialized:
            print("\nStartup search stopped before choosing an initial root.")
            return
        current_index = max(roots)
        save_level_roots(level, roots)

    while True:
        print_tower(level, roots)
        depth = read_jump_depth(current_index)
        if depth is None:
            break

        level, roots, success = direct_jump_search(level, roots, current_index, depth)
        if not success:
            break

        current_index += depth
        save_level_roots(level, roots)

    print_tower(level, roots)


interactive_search()
