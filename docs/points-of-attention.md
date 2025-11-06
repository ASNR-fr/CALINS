# Points of Attention

⚠️ **Important Considerations** (read this section along with the [Theory](theory/introduction.md) chapter and examples in `USER_EXAMPLES.ipynb`)

## Special Isotope Handling

- **Bounded isotopes**: Some isotopes in sensitivity files are not present in variance-covariance matrices. TSURFER (SCALE) uses data from the closest isotope, and CALINS detects these cases when the hundreds digit of the ID is ≥ 3:
  - `be-9 bounded` (ID=4309 or 4509) → associated with covariances of `be-9` (ID=4009)
  - `h-1 bounded` (ID=1901 or 1801) → associated with covariances of `h-1` (ID=1001)

- **Natural isotopes**: Isotopes with 'natural' evaluated forms:
  - `c-12 bounded` (ID=6312) → associated with covariances of `c-0` (ID=6000)

## Reaction ID Considerations

- **CAPTURE vs N,GAMMA**: Some covariance matrices have data for CAPTURE reaction (ID=101), others for N,GAMMA reaction (ID=102), which are very close. Be careful when associating these reactions.

- **MCNP negative IDs**: SDF files from MCNP calculations can contain reactions with negative IDs:
  - CALINS associates reaction `-2` to reaction `101` if the `mcnp=True` flag is enabled when creating the Case object
  - Other negative reactions trigger a warning

## Sensitivity Profile Occurrences

When an SDF file has multiple sensitivity profiles for the same isotope-reaction pair (same IDs), TSURFER (SCALE) parses profiles differently. CALINS provides three rules:
- `occurrences_rule='first'`: Use the first occurrence as the sensitivity profile
- `occurrences_rule='sum'` (default): Sum all occurrences for the profile
- `occurrences_rule='last'`: Use the last occurrence