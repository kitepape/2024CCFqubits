
def parse_originir_to_list(ir_string):
    # 分割输入字符串，获取每一行指令
    commands = ir_string.strip().split('\n')
    parsed_commands = []

    for command in commands:
        # 去除空格和额外的字符
        command = command.strip()

        if not command:
            continue

        # 识别指令类型和相关参数
        if command.startswith('CNOT'):
            # CNOT q[1],q[8] -> ('CNOT', [1, 8], None, None)
            parts = command.split()
            qubits = [int(q.split('[')[1].split(']')[0]) for q in parts[1].split(',')]
            parsed_commands.append(('CNOT', qubits, None, None))

        elif command.startswith('BARRIER'):
            # BARRIER q[0],q[1],... -> ('BARRIER', [0, 1, ...], None, None)
            parts = command.split()
            qubits = [int(q.split('[')[1].split(']')[0]) for q in parts[1].split(',')]
            parsed_commands.append(('BARRIER', qubits, None, None))

        elif command.startswith('U3'):
            # U3 q[8],(0,-2.30112297980874,-2.30112297980874) -> ('U3', 8, None, (0, -2.3011, -2.3011))
            parts = command.split(',')
            qubit = int(parts[0].split('[')[1].split(']')[0])
            angles = eval(','.join(parts[1:]))  # Converts the tuple string to a tuple
            parsed_commands.append(('U3', qubit, None, angles))

        elif command.startswith('MEASURE'):
            # MEASURE q[0] c[0] -> ('MEASURE', 0, 0, None)
            parts = command.split(' ')
            qubit = int(parts[1].split('[')[1].split(']')[0])
            cbit = int(parts[1].split('[')[2].split(']')[0])
            parsed_commands.append(('MEASURE', qubit, cbit, None))

    return parsed_commands

def print_parse(parsed_output, f=None):
    if f==None:
        print('[',end='')
        for item in parsed_output:
            if item != parsed_output[-1]:
                print(item,',',sep='')
            else:
                print(item,end=']\n')
    else:
        print('[',end='',file=f)
        for item in parsed_output:
            if item != parsed_output[-1]:
                print(item,',',sep='',file=f)
            else:
                print(item,end=']\n',file=f)
if __name__ == '__main__':
    # 示例输入
    ir_string = """
    CNOT q[1],q[8]
    BARRIER q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11]
    U3 q[8],(0,-2.30112297980874,-2.30112297980874)
    BARRIER q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11]
    CNOT q[1],q[8]
    MEASURE q[1],c[1]
    """
    StandList = parse_originir_to_list(ir_string)
    # 调用函数并打印结果
    print_parse(StandList)
