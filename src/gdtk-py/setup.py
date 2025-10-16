import setuptools

setuptools.setup(
    name="gdtk",
    version="0.1.0",
    author="The Gas Dynamic Toolkit developers",
    license="gpl-3.0",
    project_urls={
        "Source Code": "https://github.com/gdtk-uq/gdtk/tree/master/src/lib",
        "Documentation": "https://gdtk.uqcloud.net/docs/python/loadable-library/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.7",
    install_requires=["cffi"],
    setup_requires=[
            "setuptools_git",
            "setuptools_scm",
        ],
)
