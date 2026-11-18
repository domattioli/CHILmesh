# Fair campaign environments

Execution requires one CPython 3.11 base executable shared by the parent and all
four isolated environments. Record its realpath, SHA-256, full build identity,
architecture and platform. Force OMP_NUM_THREADS, OPENBLAS_NUM_THREADS,
MKL_NUM_THREADS and NUMEXPR_NUM_THREADS to 1 before starting the parent.

The four `requirements-*-fair.lock` files are deliberately blocking placeholders.
No installation or compatibility test was authorized during offline development.
They contain no invented package pins. `validate_lock_files` rejects them.

After operator authorization on the target host, create each venv with the same
base executable, install and import-test its tool, and exercise the adapter's
actual output contract. Freeze every direct and transitive distribution as exact
`distribution==version` lines. Record native runtime/library versions separately,
CHILmesh source revision, converter script/version/options, import compatibility
results and each physical lock's SHA-256 in the environment manifest. Metadata
versions come from importlib.metadata, never module __version__ attributes.

Do not change the base interpreter to work around a compatibility failure. Missing
locks, hashes, native inventories or compatibility proof block execution. Freeze
a new campaign manifest after any environment change; retain original attempts
and disclose replacement series rationale. Synthetic locks prove validation only.
