"""Run with sage tests/test_kernel_elements.py, or Sage's Python and unittest."""

import builtins
import contextlib
import io
import hashlib
import json
import os
import sys
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import sage.all as sage


# Sage's Python CLI leaves __file__ pointing at sage.all rather than this file.
TEST_FILE = Path(sys.argv[0] if __name__ == "sage.all" else __file__).resolve()
SOURCE = TEST_FILE.parents[1] / "kernel_elements.sage"
NAMESPACE = dict(vars(sage), __name__="kernel_elements", __file__=str(SOURCE))
exec(compile(SOURCE.read_text(), str(SOURCE), "exec"), NAMESPACE)
FiniteLevel = NAMESPACE["FiniteLevel"]


class KernelElementsTests(unittest.TestCase):
    def setUp(self):
        self.directory = self.enterContext(tempfile.TemporaryDirectory())
        self.public_directory = Path(self.directory) / "public"
        self.enterContext(patch.dict(os.environ, KERNEL_ELEMENTS_CACHE_DIR=self.directory,
                                    KERNEL_ELEMENTS_PUBLIC_CACHE_DIR=str(self.public_directory)))
        self.output = self.enterContext(contextlib.redirect_stdout(io.StringIO()))

    def tower(self):
        # This level admits x_2, so compatibility tests include successful jumps.
        level = FiniteLevel(2, 4, 3)
        x2 = level.zero_root_family(2).element()
        return level, {1: x2 ** 2, 2: x2}

    def publish(self, level, kind="E", depth=None, chunked=False):
        path = Path(level._matrix_cache_path(kind) if depth is None
                    else level._iterated_frobenius_cache_path(depth))
        data = path.read_bytes()
        self.public_directory.mkdir(exist_ok=True)
        blocks = [data[:len(data) // 2], data[len(data) // 2:]] if chunked else [data]
        parts = []
        for index, block in enumerate(blocks):
            filename = "%s.part%04d" % (path.name, index) if chunked else path.name
            (self.public_directory / filename).write_bytes(block)
            parts.append({"filename": filename, "size_bytes": len(block),
                          "sha256": hashlib.sha256(block).hexdigest()})
        entry = dict(level.matrix_cache_metadata(kind, depth), parts=parts,
                     size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        manifest = {"format_version": NAMESPACE["PUBLIC_CACHE_FORMAT_VERSION"],
                    "basis_order": NAMESPACE["MATRIX_BASIS_ORDER"],
                    "matrices": {path.name: entry}}
        manifest_path = self.public_directory / "manifest.json"
        manifest_path.write_text(json.dumps(manifest))
        return path, manifest_path, manifest

    def test_matrices_match_ring_operations_and_old_frobenius_powers(self):
        for p, M, L in [(2, 1, 0), (2, 1, 1), (2, 2, 2), (2, 3, 3), (3, 2, 2)]:
            level = FiniteLevel(p, M, L)
            self.assertFalse(hasattr(level, "basis_A"))
            for j in range(level.N):
                for i in range(level.dim_B):
                    x = level.A(level.a ** i) * level.Y ** j
                    v = level.A_to_vector(x)
                    self.assertEqual(level.vector_to_A(v), x)
                    self.assertEqual(level.E_matrix * v, level.E_to_vector(x))
                    self.assertEqual(level.F_matrix * v, level.A_to_vector(x ** p))
            for depth in range(max(M, L) + 2):
                self.assertEqual(level.iterated_frobenius_matrix(depth), level.F_matrix ** depth)
            self.assertEqual(level._build_F_matrix(1000000), level.F_matrix ** max(M, L))
            with self.assertRaises(ValueError):
                level.iterated_frobenius_matrix(-1)

    def test_initial_families_and_canonical_polynomial_lifts(self):
        for p, M, L in [(2, 2, 2), (2, 4, 3), (3, 2, 2)]:
            level = FiniteLevel(p, M, L)
            for n in range(1, min(M, L) + 1):
                R = level.reduction_matrix(n)
                Fn = level.iterated_frobenius_matrix(n)
                self.assertEqual(Fn, Fn * R.transpose() * R)
                self.assertEqual(Fn * R.transpose() * level.reduction_target(n), 0)
                family = level.zero_root_family(n)
                if family is None:
                    continue
                for v in [family.particular_vector] + [
                    family.particular_vector + d for d in family.kernel_basis()
                ]:
                    x = level.vector_to_A(v)
                    family.validator(x)
                    self.assertEqual(x.lift() ** (p ** n), level.S.gen() ** (p ** M))

    @staticmethod
    def integral_coefficient_list(coefficient):
        polynomial = coefficient.lift()
        if hasattr(polynomial, "univariate_polynomial"):
            polynomial = polynomial.univariate_polynomial()
        return [sage.ZZ(c) for c in polynomial.list()]

    def test_integral_power_retains_binomial_terms_and_free_t_variable(self):
        level = FiniteLevel(2, 2, 2)
        ring = NAMESPACE["integral_lift_ring"](level, 1)
        T = ring.gen()
        u = ring.base_ring().gen()
        self.assertEqual(ring.characteristic(), 4)
        self.assertEqual(u ** 4 + 1, 0)
        self.assertNotEqual(T ** 4, 1)
        # Keep T^4=t, reducing only the coefficients modulo 2^(1+1).
        power = NAMESPACE["integral_power"](level, 1, level.Y ** 2)
        self.assertEqual(power, T ** 4 + 2 * T ** 2 + 1)
        higher_power = NAMESPACE["integral_power"](level, 2, level.Y ** 3)
        higher_T = higher_power.parent().gen()
        self.assertEqual(higher_power.parent().characteristic(), 8)
        self.assertEqual(higher_power.degree(), 12)
        self.assertEqual(higher_power, (higher_T - 1) ** 12)

    def test_integral_lift_uses_nonnegative_coefficient_representatives(self):
        level = FiniteLevel(3, 2, 2)
        ring = NAMESPACE["integral_lift_ring"](level)
        T = ring.gen()
        u = ring.base_ring().gen()
        x = level.A(2 * level.a) * level.Y + level.A(level.a ** 2) * level.Y ** 2
        A = ring(u - 1) * (T - 1)
        B = ring((u - 1) ** 2) * (T - 1) ** 2
        self.assertEqual(NAMESPACE["integral_lift"](level, x), 2 * A + B)
        modular_ring = NAMESPACE["integral_lift_ring"](level, 1)
        modular_T = modular_ring.gen()
        modular_u = modular_ring.base_ring().gen()
        A = modular_ring(modular_u - 1) * (modular_T - 1)
        B = modular_ring((modular_u - 1) ** 2) * (modular_T - 1) ** 2
        self.assertEqual(NAMESPACE["integral_power"](level, 1, x),
                         8 * A ** 3 + 12 * A ** 2 * B + 6 * A * B ** 2 + B ** 3)

    def test_integral_solution_powers_reduce_and_match_exact_powers_modulo_ideal(self):
        for p, M, L, n in [(2, 2, 2, 1), (2, 4, 3, 2), (3, 2, 2, 1),
                           (5, 1, 1, 1), (3, 2, 2, 2)]:
            level = FiniteLevel(p, M, L)
            family = level.zero_root_family(n)
            if family is None:
                # Also test higher odd-prime precision on a general element.
                roots = [level.Y + level.A(2 * level.a) * level.Y ** 2]
            else:
                roots = [family.element(0), family.element(1)]
            K = sage.CyclotomicField(p ** (L + 1), "v")
            v = K.gen()
            I = K.ideal((v ** (p ** L) - 1) ** (n + 1))
            KP = sage.PolynomialRing(K, "V")
            V = KP.gen()
            BP = sage.PolynomialRing(level.B, "S")
            S = BP.gen()
            for index, x in enumerate(roots):
                with self.subTest(p=p, M=M, L=L, n=n, index=index):
                    lifted = NAMESPACE["integral_lift"](level, x, n=n)
                    power = NAMESPACE["integral_power"](level, n, x)
                    # Reducing further modulo xi_p-1 recovers the original
                    # polynomial representative, with T still free.
                    for integral, expected in [
                            (lifted, x.lift()(S - 1)),
                            (power, x.lift()(S - 1) ** (p ** n))]:
                        reduced = BP([
                            sum((level.B(c) * (1 + level.a) ** i
                                 for i, c in enumerate(self.integral_coefficient_list(b))),
                                level.B(0))
                            for b in integral.list()])
                        self.assertEqual(reduced, expected)
                    # Compare full polynomials against an exact number-field
                    # computation, coefficient by coefficient in the true ideal.
                    expected_lift = sum(
                        (sum((sage.ZZ(c) * (v - 1) ** i
                              for i, c in enumerate(b.lift().list())), K(0))
                         * (V - 1) ** j
                         for j, b in enumerate(x.lift().list())), KP(0))
                    actual_power = KP([
                        sum((c * v ** i
                             for i, c in enumerate(self.integral_coefficient_list(b))),
                            K(0))
                        for b in power.list()])
                    difference = actual_power - expected_lift ** (p ** n)
                    self.assertTrue(all(c in I for c in difference.list()))

    def test_odd_prime_modulus_is_the_exact_cyclotomic_ideal(self):
        for p, M, L, n in [(3, 2, 2, 2), (5, 1, 1, 1)]:
            level = FiniteLevel(p, M, L)
            ring = NAMESPACE["integral_lift_ring"](level, n)
            coefficients = ring.base_ring()
            u = coefficients.gen()
            pi = u ** (p ** L) - 1
            self.assertEqual(pi ** (n + 1), 0)
            self.assertNotEqual(pi ** n, 0)
            self.assertNotEqual(ring.gen() ** (p ** M), 1)
            if p == 3:
                self.assertEqual(coefficients(9), 0)
                self.assertNotEqual(coefficients(3), 0)
            else:
                self.assertEqual(coefficients(5), 0)

    def test_modular_power_does_not_construct_the_exact_lift_ring(self):
        factory = NAMESPACE["integral_lift_ring"]
        for p in [2, 3]:
            level = FiniteLevel(p, 2, 2)
            def modular_only(current_level, n=None):
                self.assertIsNotNone(n, "The automatic power must use modular arithmetic.")
                return factory(current_level, n)
            with patch.dict(NAMESPACE, integral_lift_ring=modular_only):
                NAMESPACE["integral_power"](level, 1, level.Y)
            self.assertFalse(hasattr(level, "_integral_lift_ring"))

    def test_supplied_x2_power_is_modulo_eight_with_t_free(self):
        level, _ = self.tower()
        a, Y = level.a, level.Y
        x = (level.A(a ** 4 + a ** 5) * Y
             + level.A(a ** 2 + a ** 3) * Y ** 2
             + level.A(a ** 3) * Y ** 3 + Y ** 4)
        power = NAMESPACE["integral_power"](level, 2, x)
        T = power.parent().gen()
        coefficients = power.parent().base_ring()
        self.assertEqual(coefficients(8), 0)
        self.assertNotEqual(coefficients(4), 0)
        self.assertEqual(power.degree(), 16)
        self.assertEqual(power.leading_coefficient(), 1)
        self.assertEqual(power[0][0], 7)
        remainder = power - (T ** 16 - 1)
        self.assertTrue(all(c % 2 == 0 for b in remainder.list()
                            for c in self.integral_coefficient_list(b)))
        self.assertTrue(all(0 <= c < 8 for b in power.list()
                            for c in self.integral_coefficient_list(b)))

    def test_integral_display_uses_current_root_index_and_defines_modulus(self):
        level, roots = self.tower()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            NAMESPACE["print_root"](level, 2, roots[2])
        displayed = output.getvalue()
        self.assertIn("T = t^(1/16)", displayed)
        self.assertIn("u is a primitive 16-th root of unity", displayed)
        self.assertIn("Integral power modulo (xi_2-1)^3 in Z[u][T]:", displayed)
        self.assertIn("Coefficient modulus: 2^3 = 8.", displayed)
        self.assertIn("xi_2 = u^8", displayed)
        self.assertIn("T^16 = t", displayed)
        self.assertIn("T is a free polynomial variable", displayed)
        self.assertIn("tilde(x_2)^(2^2) =", displayed)
        self.assertNotIn("tilde(x_1)", displayed)
        self.assertIn(str(NAMESPACE["integral_power"](level, 2, roots[2])), displayed)
        odd_level = FiniteLevel(3, 2, 2)
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            NAMESPACE["print_integral_power"](odd_level, 1, odd_level.Y ** 3)
        self.assertIn("modulo (xi_3-1)^2", output.getvalue())
        self.assertIn("xi_3 = u^9", output.getvalue())

    def test_successful_jump_and_level_transport(self):
        level, roots = self.tower()
        NAMESPACE["verify_tower"](level, roots)
        family = level.direct_root_family(roots[1], 1, 1)
        self.assertIsNotNone(family)
        for v in [family.particular_vector] + [
            family.particular_vector + d for d in family.kernel_basis()
        ]:
            family.validator(level.vector_to_A(v))
        bigger = FiniteLevel(2, 4, 4)
        transported = {r: NAMESPACE["lift_element"](x, level, bigger) for r, x in roots.items()}
        NAMESPACE["verify_tower"](bigger, transported)

    def test_cache_round_trips(self):
        level, roots = self.tower()
        NAMESPACE["save_level_roots"](level, roots)
        self.assertEqual(NAMESPACE["load_level_roots"](level), roots)
        level.iterated_frobenius_matrix(2)
        reloaded = FiniteLevel(2, 4, 3)
        self.assertEqual(reloaded.E_matrix, level.E_matrix)
        self.assertEqual(reloaded.iterated_frobenius_matrix(2), level.F_matrix ** 2)

    def test_invalid_root_caches_are_rejected(self):
        level, roots = self.tower()
        other = FiniteLevel(3, 1, 1)
        incompatible = dict(roots)
        incompatible[1] += level.A(level.a ** 7) * level.Y
        payloads = [[], {}, {"1": roots[1]}, {0: roots[1]}, {2: roots[2]},
                    {1: roots[1], 3: roots[2]}, {1: other.Y}, {1: level.A(0)},
                    {1: level.Y}, {1: sage.ZZ(0)}, incompatible]
        path = NAMESPACE["level_roots_path"](level)
        for payload in payloads:
            with self.subTest(payload=str(payload)):
                sage.save(payload, path)
                self.assertIsNone(NAMESPACE["load_level_roots"](level))
        with self.assertRaisesRegex(ValueError, "incompatible"):
            NAMESPACE["verify_tower"](level, incompatible)

    def test_bad_matrix_caches_and_unreadable_files_are_misses(self):
        level = FiniteLevel(2, 2, 2)
        path = level._matrix_cache_path("E")
        for payload in [{}, sage.identity_matrix(sage.GF(3), level.dim_A),
                        sage.identity_matrix(level.F, 1)]:
            sage.save(payload, path)
            self.assertIsNone(level._load_matrix("E"))
        for data in [b"", b"not a Sage object"]:
            Path(path).write_bytes(data)
            self.assertIsNone(level._load_matrix("E"))
        reloaded = FiniteLevel(2, 2, 2)
        self.assertEqual(reloaded.E_matrix, reloaded._build_E_matrix())

    def test_failed_atomic_write_preserves_previous_cache(self):
        path = str(Path(self.directory) / "entry.sobj")
        NAMESPACE["atomic_save"]({"old": 1}, path)
        original = Path(path).read_bytes()
        with patch.object(os, "replace", side_effect=OSError("interrupted")):
            with self.assertRaises(OSError):
                NAMESPACE["atomic_save"]({"new": 2}, path)
        self.assertEqual(Path(path).read_bytes(), original)
        self.assertEqual(list(Path(self.directory).glob(".kernel-*.sobj")), [])

        def interrupted_save(obj, temporary_path):
            Path(temporary_path).write_bytes(b"partial")
            raise OSError("interrupted serialization")

        with patch.dict(NAMESPACE, save=interrupted_save):
            with self.assertRaises(OSError):
                NAMESPACE["atomic_save"]({"new": 2}, path)
        self.assertEqual(Path(path).read_bytes(), original)
        self.assertEqual(list(Path(self.directory).glob(".kernel-*.sobj")), [])
        NAMESPACE["atomic_save"]({"new": 2}, path)
        self.assertEqual(sage.load(path), {"new": 2})

    def test_accept_or_quit_first_solution_does_not_compute_kernel(self):
        level = FiniteLevel(2, 2, 2)
        for choice in ["U", "Q"]:
            family = level.zero_root_family(1)
            with patch.object(family, "kernel_basis", side_effect=AssertionError("eager kernel")):
                with patch.object(builtins, "input", return_value=choice):
                    chosen = NAMESPACE["choose_solution"](level, family, 1)
            self.assertEqual(chosen is None, choice == "Q")
            self.assertIsNotNone(NAMESPACE["load_level_roots"](level))

    def test_count_and_next_solution_compute_directions_on_request(self):
        level = FiniteLevel(2, 2, 2)
        family = level.zero_root_family(1)
        with patch.object(builtins, "input", side_effect=["C", "N", "U"]):
            selected = NAMESPACE["choose_solution"](level, family, 1)
        self.assertIsNotNone(family._kernel_basis)
        self.assertEqual(selected, family.element(1))
        self.assertEqual(NAMESPACE["load_level_roots"](level)[1], selected)

    def test_startup_has_one_cache_prompt_for_display_or_rebuild(self):
        level = FiniteLevel(2, 2, 2)
        root = level.zero_root_family(1).element()
        NAMESPACE["save_level_roots"](level, {1: root})
        with patch.object(builtins, "input", side_effect=["D"]) as inputs:
            _, roots, selected = NAMESPACE["startup_search"](level, 1)
        self.assertTrue(selected)
        self.assertEqual(roots, {1: root})
        self.assertEqual(inputs.call_count, 1)
        with patch.object(builtins, "input", side_effect=["F", "U"]) as inputs:
            _, _, selected = NAMESPACE["startup_search"](level, 1)
        self.assertTrue(selected)
        self.assertEqual(inputs.call_count, 2)

    def test_startup_retries_same_index_after_enlargement(self):
        larger = FiniteLevel(2, 4, 3)
        for M, L in [(1, 1), (2, 2)]:
            smaller = FiniteLevel(2, M, L)
            Path(NAMESPACE["level_roots_path"](larger)).unlink(missing_ok=True)
            with patch.dict(NAMESPACE, rebuild_larger_level=lambda *args: (larger, {})):
                with patch.object(builtins, "input", side_effect=["U"]) as inputs:
                    new_level, roots, selected = NAMESPACE["startup_search"](smaller, 2)
            self.assertTrue(selected)
            self.assertIs(new_level, larger)
            self.assertEqual(max(roots), 2)
            self.assertEqual(inputs.call_count, 1)

    def test_quit_after_inspection_preserves_tower_without_selection(self):
        level = FiniteLevel(2, 2, 2)
        with patch.object(builtins, "input", side_effect=["N", "Q"]):
            _, roots, selected = NAMESPACE["startup_search"](level, 1)
        self.assertFalse(selected)
        self.assertEqual(NAMESPACE["load_level_roots"](level), roots)

    def test_quit_after_jump_preserves_compatible_displayed_tower(self):
        level, roots = self.tower()
        with patch.object(builtins, "input", side_effect=["Q"]):
            _, displayed, selected = NAMESPACE["direct_jump_search"](level, {1: roots[1]}, 1, 1)
        self.assertFalse(selected)
        self.assertEqual(displayed[1], roots[1])
        self.assertEqual(max(displayed), 2)
        self.assertEqual(NAMESPACE["load_level_roots"](level), displayed)

    def test_interactive_flow(self):
        for choices in [["2", "2", "2", "1", "U", "Q"],
                        ["2", "2", "2", "1", "D", "Q"]]:
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                with patch.object(builtins, "input", side_effect=choices):
                    NAMESPACE["interactive_search"]()
            self.assertIn("tilde(x_1)^(2^1) =", output.getvalue())
            self.assertIn("T = t^(1/4)", output.getvalue())
            self.assertIn("u is a primitive 8-th root of unity", output.getvalue())

    def test_public_cache_single_file_and_chunks_are_read_only_fallbacks(self):
        level = FiniteLevel(2, 2, 2)
        for chunked in [False, True]:
            with self.subTest(chunked=chunked):
                level._save_matrix("E", level.E_matrix)
                local_path, _, _ = self.publish(level, chunked=chunked)
                local_path.unlink()
                before = {p.name: p.read_bytes() for p in self.public_directory.iterdir()}
                with patch.object(FiniteLevel, "_build_E_matrix", side_effect=AssertionError("cache miss")):
                    reloaded = FiniteLevel(2, 2, 2)
                self.assertEqual(reloaded.E_matrix, level.E_matrix)
                self.assertFalse(local_path.exists())
                self.assertEqual(before, {p.name: p.read_bytes() for p in self.public_directory.iterdir()})

    def test_local_matrix_precedes_public_cache(self):
        level = FiniteLevel(2, 2, 2)
        self.publish(level)
        with patch.dict(NAMESPACE, load_public_cached_object=lambda *args: self.fail("public lookup")):
            self.assertEqual(level._load_matrix("E"), level.E_matrix)

    def test_corrupt_local_matrix_falls_back_to_public(self):
        level = FiniteLevel(2, 2, 2)
        local_path, _, _ = self.publish(level)
        local_path.write_bytes(b"broken")
        self.assertEqual(level._load_matrix("E"), level.E_matrix)

    def test_public_frobenius_power_fallback(self):
        level = FiniteLevel(2, 2, 2)
        expected = level.iterated_frobenius_matrix(2)
        local_path, _, _ = self.publish(level, "Fpower", 2, chunked=True)
        local_path.unlink()
        level._frobenius_powers.pop(2)
        with patch.object(level, "_build_F_matrix", side_effect=AssertionError("cache miss")):
            self.assertEqual(level.iterated_frobenius_matrix(2), expected)

    def test_public_cache_rejects_bad_metadata_and_checksums(self):
        level = FiniteLevel(2, 2, 2)
        for mutation in ["field", "basis", "version", "part_checksum", "whole_checksum", "missing", "path", "type"]:
            with self.subTest(mutation=mutation):
                level._save_matrix("E", level.E_matrix)
                local_path, manifest_path, manifest = self.publish(level, chunked=True)
                local_path.unlink()
                entry = manifest["matrices"][local_path.name]
                if mutation == "field":
                    entry["coefficient_field"] = "GF(3)"
                elif mutation == "basis":
                    manifest["basis_order"] = "other order"
                elif mutation == "version":
                    manifest["format_version"] = 99
                elif mutation == "part_checksum":
                    (self.public_directory / entry["parts"][0]["filename"]).write_bytes(b"broken")
                elif mutation == "whole_checksum":
                    entry["sha256"] = "0" * 64
                elif mutation == "missing":
                    (self.public_directory / entry["parts"][0]["filename"]).unlink()
                elif mutation == "path":
                    entry["parts"][0]["filename"] = "../outside.sobj"
                else:
                    manifest["matrices"] = []
                manifest_path.write_text(json.dumps(manifest))
                self.assertIsNone(level._load_matrix("E"))
        manifest_path.write_text("not JSON")
        self.assertIsNone(level._load_matrix("E"))

    def test_export_tool_validates_and_exports_a_small_matrix(self):
        level = FiniteLevel(2, 2, 2)
        filename = Path(level._matrix_cache_path("E")).name
        result = subprocess.run([sys.executable, str(SOURCE.parent / "tools" / "export_public_cache.py"),
                                 "--source", self.directory, "--output", str(self.public_directory), filename],
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        Path(level._matrix_cache_path("E")).unlink()
        self.assertEqual(level._load_matrix("E"), level.E_matrix)
        manifest = json.loads((self.public_directory / "manifest.json").read_text())
        self.assertEqual(manifest["matrices"][filename]["validation"]["method"],
                         "sampled_columns_and_sparse_polynomial")

    def test_export_tool_rejects_a_mathematically_wrong_matrix(self):
        level = FiniteLevel(2, 2, 2)
        path = Path(level._matrix_cache_path("E"))
        sage.save(sage.zero_matrix(level.F, level.dim_A, sparse=True), str(path))
        result = subprocess.run([sys.executable, str(SOURCE.parent / "tools" / "export_public_cache.py"),
                                 "--source", self.directory, "--output", str(self.public_directory), path.name],
                                text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("fails the ring check", result.stderr)
        self.assertFalse((self.public_directory / "manifest.json").exists())


if __name__ in ("__main__", "sage.all"):
    unittest.main(module=__name__)
