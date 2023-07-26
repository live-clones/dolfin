import os
import re
import sys
import platform
import subprocess
import multiprocessing

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

if sys.version_info < (3, 5):
    print("Python 3.5 or higher required, please upgrade.")
    sys.exit(1)

VERSION = "2019.2.0.dev0"
RESTRICT_REQUIREMENTS = ">=2019.2.0.dev0,<2019.3"
UFL_RESTRICT_REQUIREMENTS = ">=2021.1.0" # UFL 2019.2.0.dev0 does not exist

REQUIREMENTS = ["numpy",
                "pkgconfig",
                "pybind11",
                "mpi4py",
                "petsc4py" ,
                "fenics-ffc",
                "fenics-ufl-legacy",
                "fenics-dijitso"]

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            if "CIRCLECI" in os.environ:
                build_args += ['--', '-j2']
            else:
                num_build_threads = max(1, multiprocessing.cpu_count() - 1)
                build_args += ['--', '-j' + str(num_build_threads)]

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp, env=env)

setup(name='fenics-dolfin',
      version=VERSION,
      author='FEniCS Project',
      description='DOLFIN Python interface',
      long_description='',
      packages=["dolfin",
                "dolfin.common",
                "dolfin.function",
                "dolfin.fem",
                "dolfin.la",
                "dolfin.io",
                "dolfin.jit",
                "dolfin.mesh",
                "dolfin.multistage",
                "dolfin.parameter",
                "dolfin_utils.meshconvert",
                "dolfin_utils.test",
                "fenics"],
      ext_modules=[CMakeExtension('dolfin.cpp')],
      cmdclass=dict(build_ext=CMakeBuild),
      install_requires=REQUIREMENTS,
      zip_safe=False)
