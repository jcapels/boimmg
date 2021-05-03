import os

import setuptools
from boimmgpy import definitions

with open("boimmgpy/README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="boimmg",
    version="0.1.8",
    author="João Capela",
    author_email="jcapels96@gmail.com",
    description="Find and tag Git commits based on version numbers in commit messages.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://initialcommit.com/projects/git-tagup",
    packages=setuptools.find_packages(),
    include_package_data=True,
    # package_dir={'project_initializer': 'boimmgpy'},
    # package_data={'project_initializer':
    #     os.listdir('configs/') +
    #      os.listdir("required_files/biochemistry/Aliases/")}

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'gitpython',
        'pandas',
        'cobra',
        'cobrababel',
        'neo4j',
        'sympy',
        'biocyc',
        'biopython',
        'requests',
        'numpy',
        'python-libsbml'
    ],

    keywords='git tag boimmg tagup tag-up version autotag auto-tag commit message',
    project_urls={
        'Homepage': 'https://github.com/jcapels/boimmg',
    },
    # entry_points={
    #     'console_scripts': [
    #         'git-tagup=git_tagup.__main__:main',
    #         'gtu=git_tagup.__main__:main',
    #     ],
    # },
)