from setuptools import setup, find_packages, Extension
import numpy, sys
import re


# auto-updating version code stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)


def get_requires():
    reqs = []
    for line in open("requirements.txt", "r").readlines():
        reqs.append(line)
    return reqs


setup(
    name="exoSpin",
    version=get_property("__version__", "exoSpin"),
    description="exoSpin: The tool to measure obliquities of exoplanet",
    url="https://github.com/exoAtmospheres/exoSpin",
    author="Idriss Abdoulwahab & Paulina Palma-Bifani",
    author_email="",
    license="MIT",
    packages=find_packages(),
    package_data={"": ["kernels/*.cu"]},
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="Obliquities Exoplanets",
    install_requires=get_requires(),
)