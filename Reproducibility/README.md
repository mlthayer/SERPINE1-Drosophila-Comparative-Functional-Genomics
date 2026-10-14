# Processed workflow reproduction

Run from the repository root with the R packages documented in README.md:

```sh
bash Reproducibility/reproduce_processed_workflow.sh
```

This route uses retained processed inputs and regenerates outputs. It does not
download large source datasets or independently repeat every source extraction.
The full route passed in an isolated source copy with the installed R runtime;
processed_result_comparison.csv records agreement of thirteen result tables.
Later integration filename cleanup was rerun with identical numeric contents.
Session/package records distinguish installed core packages from optional
inspection dependencies. See ../DATA_AVAILABILITY.md for verification limits.
