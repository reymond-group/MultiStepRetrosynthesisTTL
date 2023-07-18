from setuptools import setup

setup(
    name='rxnmarkcenter',
    version='0.6.0',
    packages=['rxnmarkcenter'],
    install_requires=[
        'requests',
        'importlib-metadata; python_version == "3.6"',
    ],
)