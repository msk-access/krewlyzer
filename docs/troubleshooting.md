# Troubleshooting

## Common Issues

### File Not Found Error
**Error**: `FileNotFoundError: [Errno 2] No such file or directory`
**Solution**: Ensure all input files (BAM, FASTA, BED) exist and paths are correct. Use absolute paths to avoid ambiguity.

### Permission Error
**Error**: `PermissionError: [Errno 13] Permission denied`
**Solution**: Check that you have write permissions for the output directory.

### Missing Dependencies
**Error**: `ModuleNotFoundError: No module named '...'`
**Solution**: Ensure Krewlyzer is installed in your current environment.
```bash
uv pip install krewlyzer
```
Or use the Docker image.

### Reference Mismatch
**Issue**: Results look wrong or empty.
**Solution**: Ensure your BAM files and the reference FASTA are from the same genome build (e.g., both hg19 or both hg38). Krewlyzer defaults to hg19 for provided data files.

### Memory Errors
**Issue**: Process crashes on large BAM files.
**Solution**: Increase available RAM (≥16GB recommended). Reduce the number of threads if running in parallel.
