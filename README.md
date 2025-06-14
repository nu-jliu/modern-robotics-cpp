# Modern Robotics C++

[![Build Status](https://github.com/nu-jliu/modern-robotics-cpp/actions/workflows/modern-robotics.yml/badge.svg)](https://github.com/nu-jliu/modern-robotics-cpp/actions)

A comprehensive Modern C++ implementation of robotics algorithms from the textbook [*Modern Robotics: Mechanics, Planning, and Control*](http://modernrobotics.org) by Kevin Lynch and Frank Park (Cambridge University Press 2017).

**Author**: Jingkun Liu

## 📖 Documentation

- **API Documentation**: [docs.allen-liu.net/modern-robotics-cpp](https://www.allen-liu.net/modern-robotics-cpp/)
- **Source Code**: [github.com/nu-jliu/modern-robotics-cpp](https://github.com/nu-jliu/modern-robotics-cpp)
- **Textbook**: [modernrobotics.org](http://modernrobotics.org)

## 🛠️ Requirements

| Component | Minimum Version | Tested Version |
|-----------|----------------|----------------|
| **OS** | Ubuntu 20.04+ | Ubuntu 22.04.5 LTS |
| **CMake** | 3.16+ | 3.22.1 |
| **Compiler** | GCC 9+ or Clang 10+ | GCC 11.4.0 |
| **Armadillo** | 9.900+ | 10.8.2 |

## 🚀 Features

This library provides a complete implementation of fundamental robotics algorithms organized by textbook chapters:

### Core Modules

- **🔄 Rigid-Body Motions** (Chapter 3)
  - Rotation matrices and homogeneous transformations
  - Screw theory and exponential coordinates
  - Adjoint representations

- **🎯 Forward Kinematics** (Chapter 4)
  - Product of exponentials formula
  - Open-chain manipulator kinematics
  - Body and space frame representations

- **⚡ Velocity Kinematics & Statics** (Chapter 5)
  - Jacobian computation and analysis
  - Velocity relationships and singularities
  - Static force analysis

- **📈 Inverse Kinematics** (Chapter 6)
  - Newton-Raphson iterative algorithms
  - Numerical inverse kinematics solutions

- **🔧 Dynamics of Open Chains** (Chapter 8)
  - Forward and inverse dynamics
  - Mass matrix computation
  - Coriolis and gravitational effects

- **📊 Trajectory Generation** (Chapter 9)
  - Point-to-point trajectory planning
  - Polynomial and trapezoidal profiles

## 📦 Installation

### Prerequisites

Install required dependencies on Ubuntu:

```bash
# Install build tools and dependencies
sudo apt update
sudo apt install build-essential cmake git
sudo apt install libarmadillo-dev
```

### Build from Source

```bash
# Clone the repository
git clone https://github.com/nu-jliu/modern-robotics-cpp.git
cd modern-robotics-cpp

# Configure and build
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTS=ON -DBUILD_DOCS=ON ..
make -j$(nproc)

# Install system-wide (optional)
sudo make install
```

### Build Options

- `-DBUILD_TESTS=ON` - Build unit tests
- `-DBUILD_DOCS=ON` - Generate API documentation  
- `-DCMAKE_BUILD_TYPE=Release` - Optimized release build

## 💡 Quick Start

### Basic Usage Example

```cpp
#include <iostream>
#include <modern_robotics/all.hpp>  // Include entire library

int main() {
    // Create a rotation matrix
    arma::mat33 R = {{1, 0, 0},
                     {0, 0, -1},
                     {0, 1, 0}};
    
    // Convert to axis-angle representation
    arma::vec3 omg;
    double theta;
    mr::AxisAng3(R, omg, theta);
    
    std::cout << "Rotation axis: " << omg.t();
    std::cout << "Rotation angle: " << theta << " rad" << std::endl;
    
    return 0;
}
```

### CMake Integration

```cmake
find_package(modern_robotics REQUIRED)

add_executable(your_robot_app src/main.cpp)
target_link_libraries(your_robot_app PUBLIC modern_robotics::modern_robotics)
```

### Testing Your Installation

```bash
# Run all tests
cd build && ctest

# Run specific module tests
./build/tests/test_rigid_body_motions
./build/tests/test_forward_kinematics
```

## 🏗️ Project Structure

```
modern-robotics-cpp/
├── include/modern_robotics/    # Public headers
│   ├── all.hpp                # Convenience header (includes everything)
│   ├── rigid_body_motions.hpp
│   ├── forward_kinematics.hpp
│   └── ...
├── src/                       # Implementation files
├── tests/                     # Unit tests
├── docs/                      # Documentation
└── examples/                  # Usage examples
```

## 🧪 Development

### Running Tests

```bash
cd build
ctest --verbose              # Run all tests with output
ctest -R rigid_body         # Run specific test group
```

### Generating Documentation

```bash
cd build
make doc                    # Generates docs in build/docs/html/
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## 📚 References

- Lynch, K. M., & Park, F. C. (2017). *Modern Robotics: Mechanics, Planning, and Control*. Cambridge University Press.
- Official textbook website: [modernrobotics.org](http://modernrobotics.org)

## 👨‍💻 Author

**Jingkun Liu**  
- GitHub: [@nu-jliu](https://github.com/nu-jliu)
- Documentation: [allen-liu.net](https://www.allen-liu.net)
