def transfer_index(i, j, a=6):
    return a*j+i


def get_cz_patterns(a=6, b=12):
    cz_pattern_1 = []
    cz_pattern_2 = []
    cz_pattern_3 = []
    cz_pattern_4 = []

    for i in range(a):
        for j in range(b):
            if (i%2) == (j%2):
                if j + 1 < b:
                    cz_pattern_1.append([transfer_index(i, j, a), transfer_index(i, j+1, a)])
                if i + 1 < a:
                    cz_pattern_2.append([transfer_index(i, j, a), transfer_index(i+1, j, a)])
            else:
                if j + 1 < b:
                    cz_pattern_3.append([transfer_index(i, j, a), transfer_index(i, j+1, a)])
                if i + 1 < a:
                    cz_pattern_4.append([transfer_index(i, j, a), transfer_index(i+1, j, a)])

    return [cz_pattern_1, cz_pattern_2, cz_pattern_3, cz_pattern_4]


def check_topo(program_body, cz_patterns):
    available_topology = sum(cz_patterns, [])
    for gate in program_body:
        operator, qubit, cbit, parameter = gate
        if operator == 'BARRIER':
            continue
        if operator in ['CNOT', 'CZ']:
            if list(sorted(qubit)) not in available_topology:
                raise ValueError(f'Quantum circuit do not fit the chip topology. Connection between {qubit} is not allowed')
    return True


def can_parallel(gate, gate_list, topo):
    """
    given a gate, and a layer of gates, determine whether the gate can act simultaneously with layer.

    :param gate: (list),
        List of qubits of the gate.

    :param gate_list: (List[List]),
        A list of gates, which be in one layer, and can act simultaneously.

    :param topo: (List[List[List]])
        Parallel structure of two-qubit gates in the chip. Two gates in the same sublist can act at the same time

    """
    current = []
    gate = sorted(gate)
    for g in topo:
        if gate in g:
            current = g
            break
    if gate in gate_list:
        return False
    for new_gate in gate_list:
        if new_gate not in current:
            return False
    return True


def get_layers_topo(program_body, topo=None):
    """
    Rearrange the order of quantum gates, and divide them into several layers without changing the effect,
    each layer either contains only single gates, or only contains double gates that can act simultaneously.

    :param ir: (str),
        Qpanda program body of the circuit.

    :param topo: (List[List[List]])
        Parallel structure of two-qubit gates in the chip. Two gates in the same sublist can act at the same time

    """
    gate_layers = {}
    qubit_layer_count = {}
    measure_layer = []

    for gate in program_body:
        operator, qubit, cbit, parameter = gate
        if operator not in ['BARRIER', 'MEASURE', 'U3', 'CNOT']:
            raise ValueError(f'only CNOT, U3, BARRIER and MEASURE gates are allowed. Gate {operator} is not allowed.')
        if operator == 'BARRIER':
            layer = max([qubit_layer_count.get(q, 0) for q in qubit]) + 1
            for q in qubit:
                qubit_layer_count[q] = layer + 1
        elif operator == 'MEASURE':
            measure_layer.append(gate)
            continue
        elif isinstance(qubit, int):
            if qubit not in qubit_layer_count:
                qubit_layer_count[qubit] = 0
            layer = qubit_layer_count.get(qubit, 0)
        else:
            if all(q not in qubit_layer_count for q in qubit):
                for q in qubit:
                    qubit_layer_count[q] = 0
                layer = 0
            else:
                layer = max([qubit_layer_count.get(q, 0) for q in qubit])
                while True:
                    if layer not in gate_layers:
                        break
                    layer_gates = [i[1] for i in gate_layers[layer]]
                    if can_parallel(qubit, layer_gates, topo):
                        break
                    else:
                        layer = layer + 1
                for q in qubit:
                    qubit_layer_count[q] = layer

        if layer not in gate_layers:
            gate_layers[layer] = [gate]
        else:
            gate_layers[layer].append(gate)
        if isinstance(qubit, list) and (operator != 'BARRIER'):
            for q in qubit:
                qubit_layer_count[q] += 1
    gate_layers[max(qubit_layer_count.values()) + 1] = measure_layer
    return gate_layers


def count(program_body):
    cz_patterns = get_cz_patterns(6, 12)
    available = check_topo(program_body, cz_patterns)
    if available:
        layers = get_layers_topo(program_body, cz_patterns)
        cx_layer_count = 0
        for index, layer in layers.items():
            gate_name = [gate[0] for gate in layer]
            if ('CNOT' in gate_name) and ('U3' in gate_name):
                error_info = f'''CNOT gates should be separate from single gates. Layer {index} contains {gate_name}, which is not allowed. \nTotal layers: \n'''
                for key, value in layers.items():
                    error_info += f'{key}:{value}\n'
                raise ValueError(error_info)
            if all(i == 'CNOT' for i in gate_name):
                cx_layer_count += 1
        return cx_layer_count, layers


if __name__ == '__main__':
    cz_patterns = get_cz_patterns(6, 12)

    #Correct example
    # correct_program_body = [
    #     ('U3', 7, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 8, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 14, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 13, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 15, None, (1.5707963267949, 3.14159265358979, 1.80882157041021)),
    #     ('U3', 2, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('U3', 13, None, (0,-2.87648403286798,-2.87648403286798)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('CNOT', [8, 2], None, None),
    #     ('U3', 2, None, (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [8, 2], None, None),
    #     ('CNOT', [7, 8], None, None),
    #     ('U3', 8, None, (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [7, 8],None,  None),
    #     ('CNOT', [8, 14],None,  None),
    #     ('U3', 14,None,  (0,2.97855558343452,2.97855558343452)),
    #     ('CNOT', [8, 14],None,  None),
    #     ('CNOT', [14, 13], None, None),
    #     ('U3', 13,None,  (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [14, 13],None,  None),
    #     ('BARRIER', [7, 8, 14, 13, 15, 2], None, None),
    #     ('U3',7, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 8, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 14,None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 13,None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 2, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('MEASURE', 7, 0, None),
    #     ('MEASURE', 8, 1, None),
    #     ('MEASURE', 14, 2, None),
    #     ('MEASURE', 13, 3, None),
    #     ('MEASURE', 15, 4, None),
    #     ('MEASURE', 2, 5, None),
    #     ]
    # depth, layers = count(correct_program_body)
    # print('cx layer counts: ', depth)

    # wrong answer 1: Any gate should be transferred into U3 and CNOT gates.
    # wrong_program_body_1 = [
    #     ('H', 7, None, None),  # here!
    #     ('U3', 8, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 14, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 13, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 15, None, (1.5707963267949, 3.14159265358979, 1.80882157041021)),
    #     ('U3', 2, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('U3', 13, None, (0,-2.87648403286798,-2.87648403286798)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('CNOT', [8, 2], None, None),
    #     ('U3', 2, None, (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [8, 2], None, None),
    #     ('CNOT', [7, 8], None, None),
    #     ('U3', 8, None, (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [7, 8],None,  None),
    #     ('CNOT', [8, 14],None,  None),
    #     ('U3', 14,None,  (0,2.97855558343452,2.97855558343452)),
    #     ('CNOT', [8, 14],None,  None),
    #     ('CNOT', [14, 13], None, None),
    #     ('U3', 13,None,  (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [14, 13],None,  None),
    #     ('BARRIER', [7, 8, 14, 13, 15, 2], None, None),
    #     ('U3',7, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 8, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 14,None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 13,None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 2, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('MEASURE', 7, 0, None),
    #     ('MEASURE', 8, 1, None),
    #     ('MEASURE', 14, 2, None),
    #     ('MEASURE', 13, 3, None),
    #     ('MEASURE', 15, 4, None),
    #     ('MEASURE', 2, 5, None),
    #     ]
    # depth, layers = count(wrong_program_body_1)
    # print('cx layer counts: ', depty)

    # wrong answer 2: Should be fit with chip topology
    # wrong_program_body_2 = [
    #     ('U3', 7, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 8, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 14, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 13, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('U3', 15, None, (1.5707963267949, 3.14159265358979, 1.80882157041021)),
    #     ('U3', 2, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('U3', 13, None, (0, -2.87648403286798, -2.87648403286798)),
    #     ('CNOT', [7, 13],None,  None),
    #     ('CNOT', [8, 2], None, None),
    #     ('U3', 2, None, (0, -1.27129075497598, -1.27129075497598)),
    #     ('CNOT', [8, 2], None, None),
    #     ('CNOT', [8, 15], None, None), # here!
    #     ('U3', 15, None, (0, -1.27129075497598, -1.27129075497598)),
    #     ('CNOT', [8, 15], None, None),
    #     ('CNOT', [7, 8], None, None),
    #     ('U3', 8, None, (0, -1.27129075497598, -1.27129075497598)),
    #     ('CNOT', [7, 8], None,  None),
    #     ('CNOT', [8, 14], None,  None),
    #     ('U3', 14,None,  (0,2.97855558343452,2.97855558343452)),
    #     ('CNOT', [8, 14],None,  None),
    #     ('CNOT', [14, 13], None, None),
    #     ('U3', 13,None,  (0,-1.27129075497598,-1.27129075497598)),
    #     ('CNOT', [14, 13],None,  None),
    #     ('BARRIER', [7, 8, 14, 13, 15, 2], None, None),
    #     ('U3', 7, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 8, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 14, None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 13, None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('U3', 2, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
    #     ('MEASURE', 7, 0, None),
    #     ('MEASURE', 8, 1, None),
    #     ('MEASURE', 14, 2, None),
    #     ('MEASURE', 13, 3, None),
    #     ('MEASURE', 15, 4, None),
    #     ('MEASURE', 2, 5, None),
    #     ]
    # print('cx layer counts: ', count(wrong_program_body_2))

    # wrong answer 3: CNOT gates should be separate from single gates
    wrong_program_body_3 = [
        ('U3', 7, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
        ('U3', 8, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
        ('U3', 14, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
        ('U3', 13, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
        ('U3', 15, None, (1.5707963267949, 3.14159265358979, 1.80882157041021)),
        ('U3', 2, None, (1.5707963267949, 3.14159265358979, -3.14159265358979)),
        ('CNOT', [7, 13],None,  None),
        ('U3', 13, None, (0, -2.87648403286798, -2.87648403286798)),
        ('CNOT', [7, 13],None,  None),
        ('CNOT', [8, 2], None, None),
        ('U3', 2, None, (0, -1.27129075497598, -1.27129075497598)),
        ('CNOT', [8, 2], None, None),
        ('CNOT', [7, 8], None, None),
        ('U3', 8, None, (0, -1.27129075497598, -1.27129075497598)),
        ('CNOT', [7, 8], None,  None),
        ('CNOT', [8, 14], None,  None),
        ('U3', 14,None,  (0,2.97855558343452,2.97855558343452)),
        ('CNOT', [8, 14],None,  None),
        ('CNOT', [14, 13], None, None),
        ('U3', 13,None,  (0,-1.27129075497598,-1.27129075497598)),
        ('CNOT', [14, 13],None,  None),
        ('U3', 7, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
        ('U3', 8, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
        ('U3', 14, None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
        ('U3', 13, None,  (1.33277108317959, 4.71238898038469, 1.5707963267949)),
        ('U3', 2, None, (1.33277108317959, 4.71238898038469, 1.5707963267949)),
        ('MEASURE', 7, 0, None),
        ('MEASURE', 8, 1, None),
        ('MEASURE', 14, 2, None),
        ('MEASURE', 13, 3, None),
        ('MEASURE', 15, 4, None),
        ('MEASURE', 2, 5, None),
        ]
    depth, layers = count(wrong_program_body_3)
    print('cx layer counts: ', depth)
