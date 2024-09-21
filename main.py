from pyqpanda import PauliOperator, QProg, H, BARRIER, RZ, CNOT, RX, CPUQVM, Measure, topology_match, \
    transform_to_base_qgate, convert_qprog_to_originir
import numpy as np

from OriginIRToList import parse_originir_to_list, print_parse
from get_layers import *


def create_hamiltonian(M, aM, N, aN, m, n, A, B, alpha, beta):
    """
    :param M: x类材料种类数
    :param aM: x类材料种类数mod6，向上取整
    :param N: y类材料种类数
    :param aN: y类材料种类数mod6，向上取整
    :param m: 选中的x类材料数量
    :param n: 选中的y类材料数量
    :param A: 需要达到的目标属性一得分
    :param B: 需要达到的目标属性二得分
    :param alpha: 属性一得分矩阵
    :param beta: 属性二得分矩阵
    :return: 构造的哈密顿量
    """
    # 初始化哈密顿量
    hamiltonian, s1, s2 = PauliOperator(), PauliOperator(), PauliOperator()
    lambda1 = PauliOperator('', m * n * 11)
    lambda2 = PauliOperator('', m * n * 9)
    # 构造哈密顿量的每个项
    # alpha 和 beta 是材料组合对应的得分矩阵

    # 属性一得分的偏差平方
    for i in range(M):
        for j in range(N):
            # αij * xi * yj
            s1 += PauliOperator(f"Z{i} Z{aM + j}", alpha[i, j])

    # 总得分与A的偏差
    s1 -= A
    s1 = s1 * s1

    # 属性二得分的偏差平方
    for i in range(0, M):
        for j in range(0, N):
            # βij * xi * yj
            s2 += PauliOperator(f"Z{i} Z{aM + j}", beta[i, j])

    # 总得分与B的偏差
    s2 -= B
    s2 = s2 * s2
    #  初始化材料数量约束
    cons1, cons2 = PauliOperator(), PauliOperator()
    # ((\sum_i^M {x_i}) - m)^2
    for i in range(0, M):
        cons1 += PauliOperator(f'Z{i}', 1)
    const1 = (cons1 - PauliOperator('', m)) * (cons1 - PauliOperator('', m))
    # ((\sum_j^N {y_j}) - n)^2
    for j in range(0, N):
        cons2 += PauliOperator(f'Z{aM + j}', 1)
    const2 = (cons2 - PauliOperator('', n)) * (cons2 - PauliOperator('', n))

    hamiltonian = s1 + s2 + lambda1 * const1 + lambda2 * const2
    hamiltonian = hamiltonian.toHamiltonian(True)
    hamiltonian = sort_hamiltonian(hamiltonian)

    return hamiltonian


def sort_hamiltonian(hamiltonian):

    hamiltonian.sort(key=lambda x: (min(x[0].keys())) if len(x[0].keys()) > 0 else -1)
    # hamiltonian.sort(
    #     key=lambda x: np.sum((np.array(list(x[0].keys())) - np.mean(np.array(list(x[0].keys()))))**2)
    #     if len(x[0].keys()) > 0 else -1)
    hamiltonian.sort(
        key=lambda x: np.mean((np.array(list(x[0].keys())) - np.mean(np.array(list(x[0].keys()))))**4)
        if len(x[0].keys()) > 0 else -1)

    hamiltonian.sort(
        key=lambda x: np.mean(np.array(list(x[0].keys())))
        if len(x[0].keys()) > 0 else -1)
    return hamiltonian

def getCircuit(qubits, Hamiltonian, cbate, cgamma):
    """

    :param qubits: 量子比特
    :param Hamiltonian: 哈密顿量
    :param cbate: 待优化beta
    :param cgamma: 待优化gamma
    :return: 量子程序
    """
    prog = QProg()
    qu_num = len(qubits)
    # 制备叠加态
    for i in range(qu_num):
        prog.insert(H(qubits[i]))
    prog.insert(BARRIER(qubits))
    # Cost层哈密顿量
    for i in range(0, len(Hamiltonian)):
        if i ==0:
            preItem, nowItem, nextItem = None, Hamiltonian[i], Hamiltonian[i + 1]
            preQbs, nowQbs, nextQbs = None, list(nowItem[0].keys()), list(nextItem[0].keys())
            preLength, nowLength, nextLength = None, len(nowQbs), len(nextQbs)
            preParam, nowParam, nextParam = None, nowItem[1], nextItem[1]
        elif i == len(Hamiltonian) - 1:
            preItem, nowItem, nextItem = Hamiltonian[i - 1], Hamiltonian[i], None
            preQbs, nowQbs, nextQbs = list(preItem[0].keys()), list(nowItem[0].keys()), None
            preLength, nowLength, nextLength = len(preQbs), len(nowQbs), None
            preParam, nowParam, nextParam = preItem[1], nowItem[1], None
        else:
            preItem, nowItem, nextItem = Hamiltonian[i - 1], Hamiltonian[i], Hamiltonian[i + 1]
            preQbs, nowQbs, nextQbs = list(preItem[0].keys()), list(nowItem[0].keys()), list(nextItem[0].keys())
            preLength, nowLength, nextLength = len(preQbs), len(nowQbs), len(nextQbs)
            preParam, nowParam, nextParam = preItem[1], nowItem[1], nextItem[1]
        if preQbs == nowQbs and nowQbs == nextQbs:

            if nowLength == 0:
                continue
            elif nowLength == 1 or nowLength == 2 or nowLength == 4:
                prog.insert(RZ(qubits[nowQbs[nowLength - 1]], 2 * nowParam * cgamma))

        elif preQbs == nowQbs and nowQbs!= nextQbs:
            if nowLength == 0:
                continue
            elif nowLength == 1 or nowLength == 2 or nowLength == 4:
                prog.insert(RZ(qubits[nowQbs[nowLength - 1]], 2 * nowParam * cgamma))
                prog.insert(BARRIER(qubits))
                for j in range(nowLength - 2, -1, -1):
                    if nowQbs[j] != 5:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1]]))
                    else:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1] + 1]))

        elif preQbs != nowQbs and nowQbs == nextQbs:
            if nowLength == 0:
                continue
            elif nowLength == 1 or nowLength == 2 or nowLength == 4:
                for j in range(nowLength - 1):
                    if nowQbs[j] != 5:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1]]))
                    else:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1] + 1]))
                prog.insert(BARRIER(qubits))
                prog.insert(RZ(qubits[nowQbs[nowLength - 1]], 2 * nowParam * cgamma))
        else:
            if nowLength == 0:
                continue
            elif nowLength == 1 or nowLength == 2 or nowLength == 4:
                for j in range(nowLength - 1):
                    if nowQbs[j] != 5:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1]]))
                    else:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1] + 1]))
                prog.insert(BARRIER(qubits))
                prog.insert(RZ(qubits[nowQbs[nowLength - 1]], 2 * nowParam * cgamma))
                prog.insert(BARRIER(qubits))
                for j in range(nowLength - 2, -1, -1):
                    if nowQbs[j] != 5:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1]]))
                    else:
                        prog.insert(CNOT(qubits[nowQbs[j]], qubits[nowQbs[j + 1] + 1]))

    prog.insert(BARRIER(qubits))
    # Mixer层哈密顿量
    for j in qubits:
        prog.insert(RX(j, 2 * cbate))

    return prog

def generate_circuit():
    # 示例数据
    M, N = 5, 5  # X类材料种类数
    m, n = 3, 2  # 选择X类材料数量
    maxsize = 10
    A, B = 0.6 * m * n / 2, 0.4 * m * n / 2  # 目标属性一得分

    # 随机生成得分矩阵
    alpha = np.random.randint(1, maxsize + 2, size=(M, N)) / (maxsize + 2)
    beta = np.random.randint(1, maxsize + 2, size=(M, N)) / (maxsize + 2)
    qu_num = 12

    machine = CPUQVM()
    machine.init_qvm()
    machine.set_configure(12, 12)
    qubits = machine.qAlloc_many(qu_num)
    cubits = machine.cAlloc_many(qu_num)
    # 创建哈密顿量
    actual_M, actual_N = M + 6 - M % 6, N + 6 - N % 6
    Hp = create_hamiltonian(M, actual_M, N, actual_M, m, n, A, B, alpha, beta)
    arr1, arr2 = np.random.rand(), np.random.rand()
    print_parse(Hp)
    prog = getCircuit(qubits, Hp, arr1, arr2)
    for i in range(len(qubits)):
        prog.insert(Measure(qubits[i], cubits[i]))
    newProg, newQubits = topology_match(prog, qubits, machine, 'Config.json')
    base_prog = transform_to_base_qgate(newProg, machine, ['U3'], ['CNOT'])
    originString = convert_qprog_to_originir(base_prog, machine)
    StandardList = parse_originir_to_list(originString)
    # print(StandardList)
    return StandardList


def question2():
    with open('output.txt', 'r') as f:
        lines = f.read()
        qLst = eval(lines)
    return qLst


if __name__ == '__main__':
    lst = generate_circuit()
    depth, layers = count(lst)
    print('depth:', depth)
