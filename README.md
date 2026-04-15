# htsgen

Python bindings for htsgen, a synthetic data generation library for high-throughput sequencing.

## Developer build

### Prerequisites

- CMake >= 3.22
- A C++20 compiler
- Python >= 3.12
- htslib >= 1.14 (findable via pkg-config, or pass paths explicitly to cmake)

### Local Build

```sh
mkdir build && cd build
python3 -m venv <venv>  # optional
<venv>/bin/pip install pybind11-stubgen  # optional
cmake .. -DMAKE_PY=ON
cmake --build .
```

The build invokes `pybind11-stubgen` as a post-build step to generate Python type stubs, if it is available in the active Python environment. Create a venv and install pybind11-stubgen before running cmake if you want the stubs. If the cmake cache exists from a previous run with a different Python, rm `CMakeCache.txt` first.

Outputs:
- `build/htsgen.cpython-*.so` — compiled extension module
- `build/stubs/htsgen.pyi` — type stubs for LSP support (if enabled).

To use the module without installing, add the build directory to `PYTHONPATH`:

```sh
PYTHONPATH=build python3 -c "import htsgen"
```

### Pip Install

To install into a virtualenv instead:

```sh
python3 -m venv pyenv
source pyenv/bin/activate
pip install <path-to-repo>  # or use git+url
# on macOS, you might need to select
# Apple clang over Homebrew clang.
# CC=/usr/bin/clang CXX=/usr/bin/clang++ pip install <...>
```

All dependcies should be automatically resolved. The bindings are available via `import htsgen` in Python. To inspect the API, use `dir(htsgen)` or read `python/bindings.cpp`.
