import pytest
from math import log2

import poseidon.round_numbers as rn
import poseidon


@pytest.mark.parametrize("t , alpha, prime, security_level, full_round, partial_round", [
    # using https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/generate_parameters_grain.sage
    (5, 3, poseidon.prime_64, 128, 8, 41),
])
def test_calc_round_numbers(t, alpha, prime, security_level, full_round, partial_round):
    actual_full_round, actual_partial_round, actual_half_full_round = rn.calc_round_numbers(log2(prime), security_level,
                                                                                            t, alpha, True)

    assert actual_full_round == full_round
    assert actual_partial_round == partial_round
    assert actual_half_full_round == int(full_round / 2)
