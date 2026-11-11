# Contributing

Contributions and bug reports are welcome after the project is published.

## Development setup

```bash
python -m pip install -e '.[test]'
python -m pytest -q
```

Use a new branch for changes and submit focused pull requests. Add tests for
mathematical changes and include short reproducible examples for reported bugs.
The probability vector is always in lexicographic `argsort` pattern order; do
not change its convention without a major-version review. Every numerical
change must preserve the exact spectral identities and Markov certificate checks.

The `ordpy` adapter is optional and only asserts agreement on tie-free inputs.
No exact-byte cross-platform plot identity is promised: layout and coordinates
are fixed, while font rendering may vary.
