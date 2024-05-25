# Tools, demo, and test vectors for SMAC.

## Compile

Make sure that you have CMake, clang++, and git available in your environment. Maybe a new version of g++ also works, but I didn't use it. Remember to read the [Issues](https://github.com/0NG/smac-tools?tab=readme-ov-file#issues) when you start to compile.
```
git clone https://github.com/jarro2783/cxxopts.git ./src/3rd/cxxopts

git clone https://github.com/google/or-tools.git ./src/3rd/or-tools

mkdir build
cd build
CC=clang CXX=clang++ cmake ..
make
```

## Run

The help manual of `./build/bin/mac_universal` shows everything you need. Examples are as follows.
**Please remember that these examples are just examples!**
The runtime depends on the hardware and we did not use only 16 threads.
A stronger sigma usually means a longer runtime.
For some cases, we even used 48 threads (and about 48G RAM) on AMD Zen to run for a few days.
However, sometimes we also found that 12 threads was better than 32 threads.
It could be because of the strategy used by OR-Tools.
So, play with your hardware and try more different parameters~

To test t+3 for all sigmas in a file called `dataset.txt`, run the following command that uses 16 threads, outputs a detailed log, and stops after 57600 seconds or when a trail with 19 active sboxes is found.
```
mac_universal --file ./dataset.txt --round 3 --thread 16 --threshold 19 -t 57600 -l 1
```
If you don't specify the threshold and timeout, it will run until the optimal trail is confirmed.
```
mac_universal --file ./dataset.txt --round 3 --thread 16 -l 1
```
**Note that** a special format of the file is **required**! Otherwise you will get nothing or random stuff. The content of the file should look like the following one where each sigma starts with 'sg={'. This is due to some legacy format at the beginning of our research.
```
sg={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
sg={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
```

To test t+3 for only one sigma, simply replace the `--file` option with `-s` like
```
mac_universal -s 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 --round 3 --thread 16 --threshold 19 -t 57600 -l 1
```

To run the demo of SMAC, use `./build/bin/smac_demo`.

## Test vectors

The whole list of test vectors is in `test-vectors.txt`.

## Issues

* It's normal that some errors appear when you clone OR-Tools. No worries, just try more times.
* It's normal that you spend a little bit long time on compiling the program at the first time. Take a coffee and have a break.
* It's **random** if your complier complains about `abseil`. I don't know why and how to solve. I only had this once and a re-compilation solved. However, some people couldn't solve it for a long time.
* I use **C++20**, so please use a modern C++ compiler. Otherwise, unknown errors could raise. FYI, mine is clang++-15.
* As you can see in the [CMakeLists.txt](https://github.com/0NG/smac-tools/blob/main/CMakeLists.txt), CMake version older than 3.27 is not supported. This is mostly due to OR-Tools. Also, I use new version of CMake. To make sure that everything works correctly, please follow this requirement and update your CMake via Kitware's repo if you're using Ubuntu.
* I only tested `mac_universal` on my Intel i7-11xxx and the AMD Zen in the data center. Not sure if it works on other platform. Feel free to report any problem not listed above.
