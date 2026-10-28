# Contributing to NEO-RT

## Tests

```bash
make test      # CTest suite, golden-record check and pytest
make ctest     # Fortran tests only
make pytest    # Python tests only
```

Activate the Python environment from `./setup-venv.sh` first on systems with an
externally managed Python (see [doc/building.md](doc/building.md)).

## Documentation

GitHub Actions builds the FORD documentation on pull requests and merge groups
and publishes it to GitHub Pages on pushes to `main`. Build it locally with
`make doc`.

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.
