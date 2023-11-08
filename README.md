# di2hap
A tool for converting diploid genotypes to haploid.

## Build and Install
```
cmake -P cmake/get-dependencies.cmake cget
mkdir build; cd build
cmake -DCMAKE_PREFIX_PATH="$(pwd)/../cget" -DCMAKE_CXX_FLAGS="-I../cget/include" -DCMAKE_BUILD_TYPE=Release ..
make
make install
```

## Usage
```
# --haploid-code is the string used in the --sex-map file to denote male samples.
di2hap input.bcf --sex-map sample_sex_map.tsv --haploid-code 1 -O bcf -o output.bcf
```
