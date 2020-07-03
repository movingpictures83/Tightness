# Tightness
# Language: R
# Input: TXT (keyword, value pairs)
# Output: Prefix
# Tested with: PluMA 1.1, R 4.0.0

PluMA plugin to compute various tightness measures on clusters within a network.

The plugin accepts as input a TXT file with (keyword, value) pairs.  Three keywords are accepted:

unthresholded: Name of the CSV file for a correlation network, assumed to not be p-value thresholded
thresholded: Name of the CSV file for a correlation network, assumed to be p-value thresholded
clusters: Name of the cluster CSV file

Tightness is computed for both networks (unthresholded and thresholded).
