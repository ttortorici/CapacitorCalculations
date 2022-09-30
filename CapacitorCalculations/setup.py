from setuptools import setup, Extension

sfc_module = Extension('idcappy',
                       sources=['module.cpp'],
                       include_dirs=['include', '..\\CapacitorCalculations'],
                       library_dirs=["\\root\\project"],
                       # define_macros=[("PY_MAJOR_VERSION", "3"), ("PY_MINOR_VERSION", "10")],
                       language="c++",
                       extra_compile_args=["/std:c++20"],
                       )

setup(
    scripts=[],
    ext_modules=[sfc_module]
)