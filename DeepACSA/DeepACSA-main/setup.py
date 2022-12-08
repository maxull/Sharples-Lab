"""Deep_acsauto package definition"""

from setuptools import setup, find_packages

from __init__ import __version__

# Read long description from file
with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name="DeepACSA",
    version=__version__,
    description=(
        "Anatomical cross-sectional area evalutaion in Ultrasound images."
    ),
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/PaulRitsche/ACSAuto_DeepLearning",
    author="Paul Ritsche",
    author_email="paul.ritsche@unibas.ch",
    maintainers=["Paul Ritsche", "Philipp Wirth", "Neil Cronin"],
    maintainers_email=["paul.ritsche@unibas.ch",
                       "philipp.m.wirth@gmail.com",
                       "neil.cronin@jyu.fi"],
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: ",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Physiology",
        "Topic :: Utilities",
    ],
    entry_points={
        'console_scripts': [
            'deep_acsa = deep_acsa_gui:main',
        ],
    },
    keywords=[
        'ultrasound',
        'physiology',
        'deep learning',
        'muscle',
    ],
    project_urls={
        "Github": "https://github.com/PaulRitsche/DeepACSA.git",
    },
    packages=find_packages(),
    include_package_data=True,
    setup_requires=[
        "setuptools_git == 1.2",
    ],
)
