# Modern Robotics CPP

**Author**: Jingkun Liu

A Modern C++ Package that provides the functions accompanying the book [*Modern Robotics:
Mechanics, Planning, and Control*](http://modernrobotics.org) (Kevin Lynch
and Frank Park, Cambridge University Press 2017)

Please refer the document for detailed usage of this package: [docs](https://www.allen-liu.net/modern-robotics-cpp/)

## Target Platform

![build](https://github.com/nu-jliu/modern-robotics-cpp/actions/workflows/modern-robotics.yml/badge.svg)

| Environment | Value |
| ----------- | ----- |
| OS | Ubuntu 22.04.5 LTS |
| Kernel | 6.12.10-76061203-generic |
| CMake | 3.22.1 |
| g++ | 11.4.0 |
| Armadillo | 10.8.2 |

## Covered Topics

1. Rigid-Body Motion: includes all functions that computes the rigid body transformations, that conducts all transformations on rotation and homogenous transformation matrix.

## Installation

Building from source:

```bash
git clone https://github.com/nu-jliu/modern-robotics-cpp.git && cd moder-robotics-cpp
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
sudo make install
```

## Example Usage

Example C++ program for using this library

```cpp
#include <iostream>
#include <armadillo>

#include <modern_robotics/rigid_body_motions.hpp>

int main(int argc, char * argv[]) {
    arma::mat33 R{
        {1, 2, 3},
        {4, 5, 6},
        {4, 5, 6}
    };

    arma::mat33 R_inv = mr::RotInv(R);
    std::cout << R_inv << std::endl;
}

```

To use with CMake package

```cmake
find_package(modern_robotics REQUIRED)

add_executable(example src/example.cpp)
include_directory(example PRIVATE include)
target_link_library(example PUBLIC modern_robotics::modern_robotics)
```
