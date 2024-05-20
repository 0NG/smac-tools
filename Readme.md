# Tools, demo, and test vectors for SMAC.

## Compile

```
git clone https://github.com/jarro2783/cxxopts.git ./src/3rd/cxxopts

git clone https://github.com/google/or-tools.git ./src/3rd/or-tools

mkdir build
cd build
CC=clang CXX=clang++ cmake ..
make
```

## Run

The help manual of `./build/bin/mac_universal` shows everything you need.

To run the demo of SMAC, use `./build/bin/smac_demo`.

## Test vectors

The whole list of test vectors is in `test-vectors.txt`.
