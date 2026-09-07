import os
import re

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout


class MortonConan(ConanFile):
    """Conan 2.x recipe for the header-only peclet-morton library."""

    name = "morton"
    # version: read from pyproject.toml (the single version source) in set_version() below.
    license = "MIT"
    description = "Fast Morton (Z-order) codes with O(1) arithmetic"
    homepage = "https://github.com/computational-chemical-engineering/peclet-morton"
    topics = ("morton", "z-order", "space-filling-curve", "octree", "header-only")
    url = "https://github.com/computational-chemical-engineering/peclet-morton"

    def set_version(self):
        with open(os.path.join(self.recipe_folder, "pyproject.toml"), encoding="utf-8") as f:
            m = re.search(r'^version\s*=\s*"([^"]+)"', f.read(), re.M)
        if not m:
            raise RuntimeError("pyproject.toml: no `version = \"X.Y.Z\"` line")
        self.version = m.group(1)

    settings = "os", "arch", "compiler", "build_type"
    package_type = "header-library"
    exports_sources = "include/*", "cmake/*", "CMakeLists.txt", "pyproject.toml", "LICENSE", "README.md"
    no_copy_source = True

    options = {"with_bmi2": [True, False]}
    default_options = {"with_bmi2": True}

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["MORTON_BUILD_TESTS"] = False
        tc.variables["MORTON_BUILD_BENCHMARKS"] = False
        tc.variables["MORTON_BUILD_BINDINGS"] = False
        tc.variables["MORTON_ENABLE_BMI2"] = bool(self.options.with_bmi2)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        # Header-only: no libs to link, just include dir + the CMake target.
        self.cpp_info.set_property("cmake_file_name", "morton")
        self.cpp_info.set_property("cmake_target_name", "morton::morton")
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        if self.options.with_bmi2 and self.settings.compiler in ("gcc", "clang", "apple-clang"):
            self.cpp_info.cxxflags = ["-mbmi2"]

    def package_id(self):
        # Header-only package: independent of settings.
        self.info.clear()
