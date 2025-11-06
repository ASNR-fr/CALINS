# Troubleshooting

## Common Issues

### Import Errors

**Problem**: `ModuleNotFoundError: No module named 'calins'`

**Solution**: Ensure CALINS is properly installed:
```bash
pip install .
# Or for development:
pip install -e .
```

### Missing Dependencies

**Problem**: Import errors for numpy, pandas, plotly, etc.

**Solution**: Install all dependencies:
```bash
pip install -r requirements.txt
```

### Sensitivity File Parsing Issues

**Problem**: Error when loading SDF files

**Solutions**:
- Ensure the file path is correct and the file exists
- Check that the SDF file format is compatible (SCALE/TSUNAMI format)
- Try using `mcnp=True` flag for MCNP-generated files with negative reaction IDs
- Verify the `occurrences_rule` parameter if you have duplicate iso-reac pairs

### Covariance Matrix Issues

**Problem**: Missing covariance data for specific isotopes

**Solution**: CALINS automatically handles some cases:
- Bounded isotopes (e.g., h-poly) use data from base isotopes
- Natural isotopes (e.g., c-0) are matched to specific isotopes when available
- Check the warnings in logs for missing data

### Matrix Inversion Errors

**Problem**: Cholesky decomposition fails during assimilation

**Solutions**:
- Reduce the number of benchmark cases
- Check for highly correlated benchmarks
- Verify that experimental uncertainties are realistic
- Consider using chi-squared or C<sub>k</sub> filtering

## Getting Help

If you encounter issues not covered here:
- Check the logs in `~/.CALINS/` for detailed error messages
- Review the `USER_EXAMPLES.ipynb` notebook for working examples
- Open an issue on GitHub with:
  - Description of the problem
  - Minimal code to reproduce the issue
  - Error messages and logs
  - CALINS version and Python version