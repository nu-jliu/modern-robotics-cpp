# Modern Robotics C++

[![Ubuntu22](https://github.com/nu-jliu/modern-robotics-cpp/actions/workflows/ubuntu22.yml/badge.svg)](https://github.com/nu-jliu/modern-robotics-cpp/actions)
[![Ubuntu24](https://github.com/nu-jliu/modern-robotics-cpp/actions/workflows/ubuntu24.yml/badge.svg)](https://github.com/nu-jliu/modern-robotics-cpp/actions)

A comprehensive Modern C++ implementation of robotics algorithms from the textbook [*Modern Robotics: Mechanics, Planning, and Control*](http://modernrobotics.org) by Kevin Lynch and Frank Park (Cambridge University Press 2017).

**Author**: Jingkun Liu

## 📖 Documentation

- **API Documentation**: [docs.allen-liu.net/modern-robotics-cpp](https://www.allen-liu.net/modern-robotics-cpp/)
- **Source Code**: [github.com/nu-jliu/modern-robotics-cpp](https://github.com/nu-jliu/modern-robotics-cpp)
- **Textbook**: [modernrobotics.org](http://modernrobotics.org)

## 🛠️ Requirements

| Component | Minimum Version | Tested Version |
|-----------|----------------|----------------|
| **OS** | Ubuntu 22.04+ | Ubuntu 22.04.5 LTS |
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
  - Cubic and quintic polynomial time scaling
  - Joint space, screw motion, and Cartesian trajectories

- **🎮 Robot Control** (Chapter 11)
  - Computed torque control implementation
  - PID feedback control with feedforward
  - Control simulation and trajectory tracking

- **🔧 Utilities**
  - Mathematical constants and tolerance settings
  - Vector normalization and near-zero checking
  - Common robotics utility functions

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

### Basic Usage Examples

#### Rigid Body Motions
```cpp
#include <iostream>
#include <modern_robotics/all.hpp>

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

#### Trajectory Generation
```cpp
#include <modern_robotics/trajectory_generation.hpp>

// Generate a joint space trajectory
arma::vec thetastart = {0, 0, 0};
arma::vec thetaend = {1.57, 0.5, -0.3};
double Tf = 3.0;  // 3 seconds
size_t N = 100;   // 100 points

auto trajectory = mr::JointTrajectory(thetastart, thetaend, Tf, N, 
                                     mr::Method::Quintic);
```

#### Robot Control
```cpp
#include <modern_robotics/robot_control.hpp>

// Compute control torques using PID + feedforward
arma::vec tau = mr::ComputeTorque(thetalist, dthetalist, eint, g, 
                                 Mlist, Glist, Slist,
                                 thetalistd, dthetalistd, ddthetalistd,
                                 kp, ki, kd);
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
├── include/modern_robotics/              # Public headers
│   ├── all.hpp                          # Convenience header (includes everything)
│   ├── rigid_body_motions.hpp           # Chapter 3: Rotations & transformations
│   ├── forward_kinematics.hpp           # Chapter 4: Forward kinematics
│   ├── velocity_kinematics_and_statics.hpp # Chapter 5: Jacobians & velocity
│   ├── inverse_kinematics.hpp           # Chapter 6: Inverse kinematics
│   ├── dynamics_of_open_chains.hpp      # Chapter 8: Dynamics algorithms
│   ├── trajectory_generation.hpp        # Chapter 9: Motion planning
│   ├── robot_control.hpp                # Chapter 11: Control algorithms
│   └── utils.hpp                        # Mathematical utilities
├── src/                                 # Implementation files (.cpp)
├── tests/                               # Unit tests (Catch2 framework)
│   ├── test_rigid_body_motions.cpp
│   ├── test_forward_kinematics.cpp
│   ├── test_velocity_kinematics_and_statics.cpp
│   ├── test_inverse_kinematics.cpp
│   ├── test_dynamics_of_open_chains.cpp
│   ├── test_trajectory_generation.cpp
│   ├── test_robot_control.cpp
│   └── test_utils.cpp
├── build/                               # Build directory (generated)
├── CMakeLists.txt                       # CMake configuration
├── Doxyfile.in                         # Doxygen documentation config
├── CLAUDE.md                           # Claude Code development guide
└── LICENSE                             # MIT License
```

## 🧪 Development

### Running Tests

```bash
cd build
ctest --verbose              # Run all tests with output
ctest -R rigid_body         # Run specific test group

# Run individual test executables
./tests/test_rigid_body_motions
./tests/test_forward_kinematics
./tests/test_velocity_kinematics_and_statics
./tests/test_inverse_kinematics
./tests/test_dynamics_of_open_chains
./tests/test_trajectory_generation
./tests/test_robot_control
./tests/test_utils
```

### Generating Documentation

```bash
cd build
make doc                    # Generates docs in build/docs/html/
# Open documentation: firefox build/docs/html/index.html
```

### Code Architecture & Design

This library follows modern C++ best practices and robotics conventions:

#### Type System
- **Armadillo matrices**: `arma::mat33`, `arma::mat44`, `arma::mat66` for fixed-size matrices
- **Armadillo vectors**: `arma::vec3`, `arma::vec6`, `arma::vec` for dynamic vectors
- **Const correctness**: Functions return `const` values to prevent accidental modification
- **Namespace**: All functionality contained within `mr` namespace

#### Parameter Conventions
- Parameter names match textbook notation exactly (e.g., `Slist`, `Blist`, `thetalist`)
- **Slist**: Screw axes in space frame (6×n matrix as vector of 6-vectors)
- **Blist**: Screw axes in body frame (6×n matrix as vector of 6-vectors)
- **Mlist**: Link frames relative to previous frame (vector of 4×4 matrices)
- **Glist**: Spatial inertia matrices (vector of 6×6 matrices)

#### Mathematical Constants
- Global tolerance: `mr::tolerance = 1e-6`
- Utility functions for near-zero checking and vector normalization

#### Testing Framework
- **Catch2** framework for comprehensive unit testing
- Floating-point comparisons use `Catch::Matchers::WithinAbs(expected, TOLERANCE)`
- Each module has corresponding test file with mathematical validation

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
- Portfolio: [allen-liu.net](https://www.allen-liu.net)
