from setuptools import setup

setup(
    name='poseidon-hash',
    version='0.1.0',
    py_modules=["hash"],
    package_dir={'': 'poseidon'},
    url='https://github.com/ingonyama-zk/poseidon-hash',
    author='Ingonyama',
    author_email='ekaterina@ingonyama.com',
    description='Reference implementation in Python of Poseidon and optimized Poseidon (Neptune) hash functions',
    include_package_data=True,
    long_description=open('README.md').read() + '\n\n' + open('CHANGELOG.md').read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=[
        'pytest~=7.1.2',
        'galois~=0.0.30',
        'numpy>=1.18.4',
        'setuptools>=42.0',
    ],
    keywords=['Hash Function', 'Cryptography'],
)
