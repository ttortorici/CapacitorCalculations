from setuptools import setup, Extension

sfc_module = Extension('dielepy',
                       sources=['module.cpp'],
                       include_dirs=['include', '..\\CapacitorCalculations'],
                       library_dirs=["\\root\\project"],
                       # define_macros=[("PY_MAJOR_VERSION", "3"), ("PY_MINOR_VERSION", "10")],
                       language="c++",
                       extra_compile_args=["/std:c++20"],
                       )

setup(
    name='dielepy',
    version='0.0',
    scripts=[],
    # packages=['src'],
    description='Allows for fast calculations related to dielectric spectroscopy analysis.',
    author='Teddy Tortorici',
    author_email='etc.tortorici@gmail.com',
    ext_modules=[sfc_module]
)