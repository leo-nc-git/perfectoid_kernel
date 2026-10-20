"""Export validated local matrices: sage tools/export_public_cache.py [options]."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import sys
import tempfile

import sage.all as sage
from sage.env import SAGE_VERSION


TOOL_FILE = Path(sys.argv[0] if __name__ == "sage.all" else __file__).resolve()
ROOT = TOOL_FILE.parents[1]
SOURCE = ROOT / "kernel_elements.sage"
NAMESPACE = dict(vars(sage), __name__="kernel_elements", __file__=str(SOURCE))
exec(compile(SOURCE.read_text(), str(SOURCE), "exec"), NAMESPACE)
MATRIX_FILENAME = re.compile(r"(E|F|Fpower)_p(\d+)_M(\d+)_L(\d+)(?:_d(\d+))?_v(\d+)\.sobj")


class ValidationLevel(NAMESPACE["FiniteLevel"]):
    """Set up the rings and evaluation points without constructing matrices."""

    def _build_matrices(self):
        pass


def validate_matrix(path, samples):
    match = MATRIX_FILENAME.fullmatch(path.name)
    if match is None:
        raise ValueError("Not a recognized matrix filename: %s" % path.name)
    kind, p, M, L, depth, version = match.groups()
    if (kind == "Fpower") != (depth is not None):
        raise ValueError("Only Fpower filenames should include a depth")
    level = ValidationLevel(int(p), int(M), int(L))
    if int(version) != level.matrix_cache_version:
        raise ValueError("Unsupported matrix cache version")
    cached = level._load_cached_matrix(str(path))
    if cached is None:
        raise ValueError("Unreadable or invalid matrix: %s" % path)
    metadata = level.matrix_cache_metadata(kind, depth)
    if kind != "E":
        if cached != level._build_F_matrix(1 if depth is None else int(depth)):
            raise ValueError("Frobenius matrix does not match the direct monomial formula")
        validation = {"method": "exact_frobenius_matrix_comparison"}
    else:
        rng = random.Random(path.name)
        columns = sorted(set([0, int(level.dim_A) - 1, int(level.dim_B) - 1]
                             + rng.sample(range(int(level.dim_A)), min(samples, int(level.dim_A)))))
        combination = level.A(0)
        for index in columns:
            j, i = divmod(index, int(level.dim_B))
            x = level.A(level.a ** i) * level.Y ** j
            if cached.column(index) != level.E_to_vector(x):
                raise ValueError("Evaluation matrix fails the ring check at column %s" % index)
            combination += level.F(rng.randrange(1, int(level.p))) * x
        if cached * level.A_to_vector(combination) != level.E_to_vector(combination):
            raise ValueError("Evaluation matrix fails the sparse polynomial check")
        validation = {"method": "sampled_columns_and_sparse_polynomial",
                      "columns_checked": len(columns), "seed": path.name}
    # Historic .sobj files do not record which Sage version originally saved them.
    metadata["validation"] = dict(validation, sage_version=SAGE_VERSION)
    return metadata


def atomic_write_bytes(path, data):
    descriptor, temporary_path = tempfile.mkstemp(prefix=".export-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_path, path)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)


def export_matrix(path, output, chunk_bytes, metadata):
    size = path.stat().st_size
    parts = []
    digest = hashlib.sha256()
    with path.open("rb") as source:
        index = 0
        while block := source.read(chunk_bytes):
            filename = path.name if size <= chunk_bytes else "%s.part%04d" % (path.name, index)
            atomic_write_bytes(output / filename, block)
            digest.update(block)
            parts.append({"filename": filename, "size_bytes": len(block),
                          "sha256": hashlib.sha256(block).hexdigest()})
            index += 1
    return dict(metadata, size_bytes=size, sha256=digest.hexdigest(), parts=parts)


def main():
    parser = argparse.ArgumentParser(description="Validate and export large cached evaluation matrices into Git-sized chunks.")
    parser.add_argument("matrices", nargs="*", help="Explicit matrix filenames in the source directory; bypass the size threshold")
    parser.add_argument("--source", type=Path, default=ROOT / ".sage-cache")
    parser.add_argument("--output", type=Path, default=ROOT / "public-cache")
    parser.add_argument("--min-mib", type=float, default=1, help="Minimum evaluation-matrix size to select (default: 1 MiB)")
    parser.add_argument("--chunk-mib", type=int, default=32, help="Maximum part size (default: 32 MiB)")
    parser.add_argument("--samples", type=int, default=16, help="Evaluation columns to sample, in addition to boundary columns")
    args = parser.parse_args()
    if args.min_mib < 0 or args.chunk_mib < 1 or args.samples < 1:
        parser.error("sizes and sample count must be positive (min-mib may be zero)")
    if args.source.resolve() == args.output.resolve():
        parser.error("source and output directories must differ")
    if any(Path(name).name != name for name in args.matrices):
        parser.error("explicit matrices must be filenames, without a directory")
    paths = ([args.source / name for name in args.matrices] if args.matrices else
             sorted(path for path in args.source.glob("E_*.sobj")
                    if path.stat().st_size > 0 and path.stat().st_size >= args.min_mib * 1024 ** 2))
    if not paths:
        parser.error("no matching nonempty matrix files")
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "manifest.json"
    manifest = {"format_version": NAMESPACE["PUBLIC_CACHE_FORMAT_VERSION"],
                "basis_order": NAMESPACE["MATRIX_BASIS_ORDER"], "matrices": {}}
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if any(previous.get(key) != manifest[key] for key in ("format_version", "basis_order")):
            raise ValueError("Existing manifest has an incompatible format")
        manifest = previous
    for path in paths:
        print("Validating %s (%.2f MiB)..." % (path.name, path.stat().st_size / 1024 ** 2), flush=True)
        metadata = validate_matrix(path, args.samples)
        entry = export_matrix(path, args.output, args.chunk_mib * 1024 ** 2, metadata)
        manifest["matrices"][path.name] = entry
        atomic_write_bytes(manifest_path, (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
        print("Exported %s in %s part(s)." % (path.name, len(entry["parts"])), flush=True)


if __name__ in ("__main__", "sage.all"):
    main()
