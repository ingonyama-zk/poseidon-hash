# Poseidon

[![PyPI version](https://badge.fury.io/py/poseidon-hash.svg)](https://badge.fury.io/py/poseidon-hash)

This repository contains a reference implementations for the original version of Poseidon [1] and the instantiation as a hash
function tuned for Filecoin [2].
Moreover, scripts to calculate the round numbers, the round constants, and the MDS matrices are also included.


## Theoretical introduction

The set of parameters (`t`, `M`, `p`, `alpha`) fully specify a unique instance of Poseidon
All other Poseidon parameters and constants are derived from these parameters.

The optimized Poseidon hash function is instantiated in the same way as the un-optimized algorithm, however the optimized
Poseidon algorithm requires the additional pre-computation of round constants, pre-sparse matrix and sparse matrices.
We implemented the pre-computation steps according to [2].

The hash function takes as input a list of elements from a certain field `F` and outputs o single field
element.
In addition, for this implementation the `input_rate` (size of input) parameter is required, which can take any value
from `1` to `t`. A longer input size is possible and it is a work in progress.


## Usage: generate new instance with pre-generated parameters

To create a new hash instance, you can use ready-made functions in [parameters.py](poseidon/parameters.py) that have all
the necessary parameters,
or you can set your own parameters that have already been generated.

For matrices and round constants, use the hex representation.

```python
import poseidon

def main():
    poseidon_simple, t = poseidon.parameters.case_simple()

    input_vec = [x for x in range(0, t)]
    print("Input: ", input_vec)
    poseidon_digest = poseidon_simple.run_hash(input_vec)
    print("Output: ", hex(int(poseidon_digest)))

    security_level = 128
    input_rate = 3
    t_opt = 4
    full_round = 8
    partial_round = 56
    alpha = 5
    poseidon_pre_generated = poseidon.OptimizedPoseidon(poseidon.HashType.CONSTINPUTLEN, poseidon.parameters.prime_255, 
                                                    security_level, alpha, input_rate, t_opt,
                                                    full_round=full_round, partial_round=partial_round,
                                                    rc_list=poseidon.parameters.round_constants_neptune, 
                                                    mds_matrix=poseidon.parameters.matrix_neptune)

    input_vec_2 = [x for x in range(0, t_opt - 1)]
    print("Input: ", input_vec_2)
    poseidon_output = poseidon_pre_generated.run_hash(input_vec_2)
    print("Output: ", hex(int(poseidon_output)))


if __name__ == "__main__":
    main()
```

## Usage: generate new instance without pre-generated parameters

The number of rounds, roound constants and matrices are optional parameters;
if they are not specified, they will be generated automatically.
For the same required input parameters the same optional parameters will always be generated.

```python
import poseidon

def main():
    security_level = 128
    input_rate = 3
    t = 4
    alpha = 5
    poseidon_new = poseidon.Poseidon(poseidon.parameters.prime_255, security_level, alpha, input_rate, t)

    input_vec = [x for x in range(0, t)]
    print("Input: ", input_vec)
    poseidon_output = poseidon_new.run_hash(input_vec)
    print("Output: ", hex(int(poseidon_output)))


if __name__ == "__main__":
    main()
```

## Running tests

Test vectors and parameters for non-optimised implementation are taken from [3].

Parameter generation for the optimized Poseidon is checked using test vectors that are taken (in some cases generated) from
the implementation of [2].

Running all tests

```commandline
python3 -m pytest tests/ 
```

Running tests by file

```commandline
python3 -m pytest tests/test_hash.py 
```

Running a specific test. After `::` you must enter the name of the specific test in the corresponding file

```commandline
python3 -m pytest tests/test_hash.py::test_not_optimized_poseidon
```

## License

The code is released under the MIT license. See [LICENSE](LICENSE) for more information.

## Contact

Feel free to [reach out](mailto:hi@ingonyama.com)! 


[1] *Poseidon: A New Hash Function for Zero-Knowledge Proof Systems. Cryptology ePrint Archive, Report 2019/458.
https://eprint.iacr.org/2019/458. Accepted at USENIX'21.*

[2] *Neptune specification https://github.com/filecoin-project/neptune/tree/master/spec*

[3] *Reference Implementations for Various Versions of Starkad and
Poseidon https://extgit.iaik.tugraz.at/krypto/hadeshash*
 
