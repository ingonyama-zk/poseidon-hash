import pytest as pytest
import galois
import functools
import numpy as np

import poseidon.round_constants as rc
import poseidon


@pytest.mark.parametrize("rc_expected, t ,full_round, partial_round, alpha, prime, prime_bit_len", [
    # references are generated using https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/generate_parameters_grain.sage
    # sage generate_parameters_grain.sage 1 0 64 9 8 41 0xfffffffffffffeff
    (poseidon.round_constants_64, 9, 8, 41, 3, poseidon.prime_64, 64),
    # sage generate_parameters_grain.sage 1 0 255 9 8 57 73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
    (poseidon.round_constants_test, 9, 8, 57, 5, poseidon.prime_255, 255),

])
def test_calc_round_constants(rc_expected, t, full_round, partial_round, alpha, prime, prime_bit_len):
    field_p = galois.GF(prime)
    rc_actual = rc.calc_round_constants(t, full_round, partial_round, prime, field_p, alpha, prime_bit_len)
    print([hex(int(x)) for x in rc_actual])
    assert len(rc_actual) == len(rc_expected) == t * (full_round + partial_round)
    assert functools.reduce(lambda x, y: x and y, map(lambda p, q: int(p, 16) == int(q), rc_expected, rc_actual),
                            True)


@pytest.mark.parametrize("t, prime, expected_matrix", [
    (4, poseidon.prime_255, poseidon.matrix_neptune),
])
def test_mds_matrix(t, prime, expected_matrix):
    field_p = galois.GF(prime)
    matrix = rc.mds_matrix_generator(field_p, t)
    hex_matrix = rc.get_hex_matrix_from_field_matrix(matrix, 64)

    assert np.array_equal(expected_matrix, hex_matrix)


@pytest.mark.parametrize("opt_rc_expected, rc_not_opt, half_full_round, partial_round, mds_matrix, prime, t", [
    (poseidon.optimized_round_constants_neptune, poseidon.round_constants_neptune, 4, 56, poseidon.matrix_neptune,
     poseidon.prime_255, 4),
])
def test_optimized_rc(opt_rc_expected, rc_not_opt, half_full_round, partial_round, mds_matrix, prime, t):
    field_p = galois.GF(prime)

    field_matrix = rc.get_field_matrix_from_hex_matrix(field_p, mds_matrix)

    field_const = field_p([int(x, 16) for x in rc_not_opt])
    split_rc = [field_p(x.tolist()) for x in np.array_split(field_const, len(field_const) / t)]

    opt_rc_actual = rc.optimized_rc(split_rc, half_full_round, partial_round, field_matrix)
    assert functools.reduce(lambda x, y: x and y,
                            map(lambda p, q: int(p, 16) == int(q), opt_rc_expected, opt_rc_actual), True)


@pytest.mark.parametrize("opt_mds_expected_pre,opt_mds_expected_sparse, mds_not_opt, partial_round, prime, t", [
    (poseidon.pre_matrix_neptune, poseidon.sparse_matrices_neptune, poseidon.matrix_neptune, 56, poseidon.prime_255, 4),
])
def test_optimized_matrix(opt_mds_expected_pre, opt_mds_expected_sparse, mds_not_opt, partial_round, prime, t):
    field_p = galois.GF(prime)
    field_matrix = rc.get_field_matrix_from_hex_matrix(field_p, mds_not_opt)

    pre_mds_actual, spare_mds_actual = rc.optimized_matrix(field_matrix, partial_round, field_p)

    hex_pre = rc.get_hex_matrix_from_field_matrix(pre_mds_actual, 64)
    hex_spare = [rc.get_hex_matrix_from_field_matrix(x, 64) for x in spare_mds_actual]

    assert np.array_equal(hex_pre, opt_mds_expected_pre)
    assert np.array_equal(hex_spare, opt_mds_expected_sparse)
