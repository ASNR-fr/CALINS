# Developer Manual - errors.py

The implemented custom errors are:
- *EmptyParsingError*: raised in *methods* when data acquisition and formatting result in an empty vector/matrix (formatting/construction error);
- *DimError*: raised when the dimensions of multiple vectors/matrices are inconsistent for matrix multiplications;
- *SensInputError*: raised when sensitivity data are not in the correct format for acquisition and formatting.