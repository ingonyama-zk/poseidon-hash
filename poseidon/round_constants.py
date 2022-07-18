from copy import deepcopy
import numpy as np


def get_field_matrix_from_hex_matrix(field_p, mds_matrix):
    """
    Convert matrix in hex representation to matrix in field representation.

    :param field_p: A field field_p of type galois.GF(p).
    :param list mds_matrix: 2-dim array of size t*t. Consist of elements in hex.
    :return: 2-dim array of size t*t. Consist of field elements.
    """
    n = len(mds_matrix)
    mds_matrix_field = field_p.Zeros((n, n))
    for i in range(0, n):
        for j in range(0, n):
            mds_matrix_field[i, j] = field_p(int(mds_matrix[i][j], 16))
    return mds_matrix_field


def get_hex_matrix_from_field_matrix(m, size):
    """
    Help function for tests. Convert matrix in field representation to matrix in hex representation.
    :param m:
    :param size:
    :return:
    """
    n = len(m)
    hex_matrix = [['' for _ in range(n)] for _ in range(n)]
    for i in range(0, n):
        for j in range(0, n):
            hex_matrix[i][j] = '0x{0:0{1}x}'.format(int(m[i][j]), size)
    return hex_matrix


def init_state_for_grain(alpha, p, prime_bit_len, t, full_round, partial_round):
    """
    The function generates the initial state for Grain LFSR  in a self-shrinking mode.

    Initialize the state with 80 bits b0, b1, . . . , b79, where

    (a) b0, b1 describe the field,
    (b) bi for 2 ≤ i ≤ 5 describe the S-Box,
    (c) bi for 6 ≤ i ≤ 17 are the binary representation of prime_bit_len,
    (d) bi for 18 ≤ i ≤ 29 are the binary representation of t,
    (e) bi for 30 ≤ i ≤ 39 are the binary representation of RF ,
    (f) bi for 40 ≤ i ≤ 49 are the binary representation of RP , and
    (g) bi for 50 ≤ i ≤ 79 are set to 1.

    :param int alpha: The power of S-box.
    :param int p: The prime field modulus.
    :param int prime_bit_len: The number of bits of the Poseidon prime field modulus.
    :param int t: The size of Poseidon's inner state
    :param int full_round: Number of full rounds
    :param int partial_round: Number of partial rounds
    :return: Initialized state with 80 elements of type int.
    :rtype list:
    """
    init_state = []
    # Choice of encoding for alpha, consistent with filecoin documentation except else
    if alpha == 3:
        exp_flag = 0
    elif alpha == 5:
        exp_flag = 1
    elif alpha == -1:
        exp_flag = 2
    else:
        exp_flag = 3

    init_state += [int(_) for _ in (bin(p % 2)[2:].zfill(2))]
    init_state += [int(_) for _ in (bin(exp_flag)[2:].zfill(4))]
    init_state += [int(_) for _ in (bin(prime_bit_len)[2:].zfill(12))]
    init_state += [int(_) for _ in (bin(t)[2:].zfill(12))]
    init_state += [int(_) for _ in (bin(full_round)[2:].zfill(10))]
    init_state += [int(_) for _ in (bin(partial_round)[2:].zfill(10))]
    init_state += [int(1)] * 30

    return init_state


def calc_round_constants(t, full_round, partial_round, p, field_p, alpha, prime_bit_len):
    """
    This function generates constants for addition at each round.
    From the poseidon paper:
    The round constants are generated using the Grain LFSR [23] in a self-shrinking mode:

    1. Initialize the state with 80 bits b0, b1, . . . , b79 using function init_state_for_grain.
    2. Update the bits using bi+80 = bi+62 ⊕ bi+51 ⊕ bi+38 ⊕ bi+23 ⊕ bi+13 ⊕ bi.
    3. Discard the first 160 bits.
    4. Calculate next bits and state using function calc_next_bits.

    Using this method, the generation of round constants depends on the specific instance, and thus different
    round constants are used even if some of the chosen parameters (e.g., n and t) are the same.

    If a randomly sampled integer is not in Fp, we discard this value and take the next one. Note that
    cryptographically strong randomness is not needed for the round constants, and other methods can also be used.

    :param int t: The size of Poseidon's inner state
    :param int full_round: Number of full rounds
    :param int partial_round: Number of partial rounds
    :param int p: The prime field modulus.
    :param field_p: A field field_p of type galois.GF(p).
    :param int alpha: The power of S-box.
    :param int prime_bit_len: The number of bits of the Poseidon prime field modulus.
    :return: List of field elements of size t * (full_round + partial_round).
        Each t element corresponds to one round constant.
    :rtype list:
    """
    rc_number = t * (full_round + partial_round)

    state = init_state_for_grain(alpha, p, prime_bit_len, t, full_round, partial_round)
    rc_field = []
    # Discard first 160 output bits:
    for _ in range(0, 160):
        new_bit = state[62] ^ state[51] ^ state[38] ^ state[23] ^ state[13] ^ state[0]
        state.pop(0)
        state.append(new_bit)

    while len(rc_field) < rc_number:
        state, bits = calc_next_bits(state, prime_bit_len)

        rc_int = int("".join(str(i) for i in bits), 2)
        if rc_int < p:
            rc_field.append(field_p(rc_int))

    return rc_field


def calc_next_bits(state, prime_bit_len):
    """
    Function generate new LFSR state after shifts new field_size number generated

    - Update the bits using bi+80 = bi+62 ⊕ bi+51 ⊕ bi+38 ⊕ bi+23 ⊕ bi+13 ⊕ bi.
    - Evaluate bits in pairs: If the first bit is a 1, output the second bit. If it is a 0, discard the second bit.

    :param list state: Current LFSR state
    :param int prime_bit_len: The number of bits of the Poseidon prime field modulus.
    :return: New LFSR state after shifts and new field_size number generated.
    :rtype list, list:
    """
    bits = []
    while len(bits) < prime_bit_len:
        new_bit_1 = state[62] ^ state[51] ^ state[38] ^ state[23] ^ state[13] ^ state[0]
        state.pop(0)
        state.append(new_bit_1)

        new_bit_2 = state[62] ^ state[51] ^ state[38] ^ state[23] ^ state[13] ^ state[0]
        state.pop(0)
        state.append(new_bit_2)

        if new_bit_1 == 1:
            bits.append(new_bit_2)

    return state, bits


def mds_matrix_generator(field_p, t):
    """
    This function generates a maximum distance separable (MDS) matrix,
    which is used in linear layer of Poseidon hush function.

    :param field_p: A field field_p of type galois.GF(p).
    :param int t: The size of Poseidon's inner state
    :return: 2-dim array of size t*t consist of filed elements
    :rtype:
    """
    x_vec = [field_p(ele) for ele in range(0, t)]
    y_vec = [field_p(ele) for ele in range(t, 2 * t)]

    mds_matrix = field_p.Zeros((t, t))
    for i in range(t):
        for j in range(t):
            mds_matrix[i, j] = (x_vec[i] + y_vec[j]) ** (-1)

    return mds_matrix


def optimized_rc(rc, half_full_round, partial_round, mds_matrix):
    """
    Given the round constants and MDS matrix for a Poseidon instance, we are able to derive optimized round constants
    for the corresponding optimized Poseidon algorithm. A description of the optimisation can be found
    here  https://github.com/filecoin-project/neptune/tree/master/spec under 'Optimized Round Constants'.

    Each full round is associated with t field elements, while each partial round is associated with one field element.

    :param list rc: List of pre-generated round constants of size full_round + partial_round.
    :param int half_full_round: half the number of full rounds
    :param int partial_round: Number of partial rounds
    :param mds_matrix: Pre-generated MDS matrix of size t*t consist of filed elements
    :return: List of field elements of size t*full_round + partial_round.
    :rtype list:
    """
    opt_rc_field = []
    m_inv = np.linalg.inv(mds_matrix)

    # pre round constant
    opt_rc_field.extend(rc[0])
    # half_full_round - 2 constants for full rounds
    for r in range(1, half_full_round):
        buf = np.dot(rc[r], m_inv)
        opt_rc_field.extend(buf)

    partial_const = []
    final_round = half_full_round + partial_round
    acc = deepcopy(rc[final_round])
    for r in (range(0, partial_round)):
        acc_1 = acc @ m_inv
        partial_const.append(acc_1[0])
        acc_1[0] = 0
        acc = acc_1 + rc[final_round - r - 1]

    # const for r = half_full_round - 1 round
    opt_rc_field.extend(acc @ m_inv)
    # partial_round constants for partial rounds
    opt_rc_field.extend(partial_const[::-1])

    # half_full_round - 1 constants for full rounds
    start = half_full_round + partial_round
    for r in range(1, half_full_round):
        opt_rc_field.extend(rc[start + r] @ m_inv)

    return opt_rc_field


def optimized_matrix(mds_matrix, partial_round, field_p):
    """
    A description of the optimisation can be found
    here https://github.com/filecoin-project/neptune/tree/master/spec under 'Sparse MDS Matrices'.

    :param 2-dim array mds_matrix:
    :param int partial_round: Number of partial rounds
    :param field_p: A field field_p of type galois.GF(p).
    :return: 2-dim array correspond to pre-sparse matrix of size t*t
        and list of 2-dim arrays correspond to sparce matrices each of which is of the size of t*t.
    """
    sparse_matrices = []
    m = deepcopy(mds_matrix)
    for r in range(0, partial_round):
        m_1, m_2 = sparse_factorize(m, field_p)
        sparse_matrices.append(m_2)
        m = mds_matrix @ m_1
    pre_matrix = m
    sparse_matrices.reverse()

    return pre_matrix, sparse_matrices


def sparse_factorize(m, field_p):
    """
    A description of the optimisation can be found
    here https://github.com/filecoin-project/neptune/tree/master/spec under 'Sparse MDS Matrices' .

    :param 2-dim array m: Current matrix for calculating the two new
    :param field_p: A field field_p of type galois.GF(p).
    :return: Two 2-dim array of size t*t
    """
    m_1 = deepcopy(m)
    m_1[0, :] = 0
    m_1[:, 0] = 0
    m_1[0, 0] = 1

    w = m[1:, 0]
    m_cap = m[1:, 1:]
    m_inv = np.linalg.inv(m_cap)
    w_cap = m_inv @ w

    m_2 = field_p.Identity(len(m))
    m_2[0, :] = m[0, :]
    m_2[1:, 0] = w_cap

    assert np.array_equal(m, m_1 @ m_2)
    return m_1, m_2
