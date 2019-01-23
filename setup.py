from distutils.core import setup

setup(
    name="protfarm",
    version="0.1",
    packages=[
        "protfarm",
        "protfarm.analysis",
        "protfarm.workspace"
    ],
    install_requires=[
        "pandas",
        "numpy",
        "peseq"
    ]
)
