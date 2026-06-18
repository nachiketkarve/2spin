# Semi-Classical 2Spin System

This repository contains C++ and Python code for analysis of a two-spin system.

## Repository structure

- `src/` - C++ source code
- `include/` - C++ header files
- `python/` - Python analysis and plotting scripts
- `data/` - dataset and parameter files
- `build/` - compiled object files and binaries (ignored by git)
- `figures/` - plots and generated visuals
- `CITATION.cff` - citation metadata

## Build

Run:

```sh
make
```

Binaries will be created under `build/bin/`.

## Clean

```sh
make clean
```

## Notes

- External dependencies such as Eigen and nlohmann/json are expected under `../Libraries/`.
- Add Python dependencies to `python/requirements.txt` if required.
