# Public matrix cache

This directory contains the 11 nonempty evaluation matrices from the initial local cache whose saved size is at least 1 MiB. Their combined size is about 583 MiB. These are the expensive computations most useful to share; small matrices and the inexpensive direct Frobenius matrices can be computed locally. Root towers remain in the private local cache.

The matrices are ordinary Git files. Files above 32 MiB are split into ordered `.sobj.part0000`, `.sobj.part0001`, … files, each at most 32 MiB. Keep every part listed in [manifest.json](manifest.json). Part boundaries split the serialized bytes; they do not change the matrix or its basis. The complete collection adds about 583 MiB to a full Git checkout; changing these binary assets also grows Git history.

The program searches the writable local cache first, then this directory. It reconstructs chunked Sage objects in temporary space, deletes the reconstructed file after loading, and never modifies this directory during a search. For the largest matrix, loading requires approximately 227 MiB of temporary disk space in addition to the memory needed for the matrix. A cache hit saves matrix construction, not the subsequent linear solve or kernel-basis computation.

The manifest records the exact parameters, matrix type, coefficient field, dimensions, matrix-cache version, basis order, whole-file checksum, and each part's size and SHA-256 checksum. The loader checks these before deserialization and also checks the resulting matrix's field and dimensions. Invalid or incomplete public entries are ignored, allowing recomputation into the local cache.

Export or extend the collection from the repository root:

```bash
sage tools/export_public_cache.py
```

The default source is `.sage-cache/`, the minimum saved size is 1 MiB, and the maximum part size is 32 MiB. Run the tool with the Python interpreter in your Sage environment for options, for example:

```bash
python tools/export_public_cache.py --help
python tools/export_public_cache.py --min-mib 10 --chunk-mib 32
python tools/export_public_cache.py E_p2_M6_L6_v1.sobj
```

Explicit filenames bypass the minimum-size threshold. Existing catalog entries are retained when adding matrices. Exports do not happen automatically during interactive searches.

Before exporting an evaluation matrix, the tool checks its field and dimensions, compares deterministic sampled columns (including boundary columns) against direct ring evaluation, and checks a sparse linear combination of those monomials. This detects many errors but is not an exhaustive check of every entry. Frobenius matrices, if explicitly exported, are checked exactly against the direct monomial formula. Checksums detect changes after export. The recorded Sage version is the validation version; the original generation version of historic `.sobj` files is unknown. Sage object portability across other versions is not guaranteed.

Set `KERNEL_ELEMENTS_PUBLIC_CACHE_DIR` to read a different public-cache directory. The default is this directory beside the original script, independent of the working directory. Set `KERNEL_ELEMENTS_CACHE_DIR` separately to choose the writable local cache.
