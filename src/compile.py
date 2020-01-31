# global variables:
species_list = []
poly_list = [] #list of lists, stores number of monomers on each polymer
config_list = []
function_list = []
polymer_list = []

###################################################################################
# compiler code

class Function(object): # types of functions
    def __init__(self, name, args, body):
        self.name = name
        self.args = args
        self.body = body # list of unprocessed strings
        
    def __eq__(self, other):
        if not isinstance(other, MyClass):
            return NotImplemented
        return self.name == other.name and self.body == other.body

class Instruction(Function):
    def __init__(self, name, args, rxns):
        self.name = name
        self.args = args
        self.rxns = rxns

    def write_rxns(self, argvals, caller_idx, idx, rxn_f):
        for rxn in self.rxns:
            rxn = rxn.replace('{i}', idx)
            rxn = rxn.replace('{i+1}', Instruction.inc(idx))
            for i, arg in enumerate(self.args):
                if argvals[i] == 'i+1':
                    argvals[i] = Instruction.inc(caller_idx)
                rxn = rxn.replace('{' + arg + '}', argvals[i])
            rxn_f.write(rxn + '\n')

    @staticmethod
    def inc(idx):
        if (idx != ''):
            i_list = idx.rpartition('.')
            if i_list[1] == '':
                return str(int(i_list[2]) + 1)
            else:
                return i_list[0] + '.' + str(int(i_list[2]) + 1)
        return ''

 # individualized per-line
 # a line has a unique argument iff caller does not contain argument
 # compute true args using this method: if line has unique not in 
 # caller's args, make that true value
 # else use true value already passed down
class Line(object):
    def __init__(self, function_str, caller):
        self.function = self.get_function(function_str)
        self.caller = caller # None if main
        self.idx = self.extract_idx(function_str)
        self.argvals = self.decide_argvals(function_str)
    
    def __eq__(self, other):
        if not isinstance(other, MyClass):
            return NotImplemented
        return self.function == other.function

    def get_function(self, function_str):
        print(function_str)
        function_str = function_str.split(':')[1].strip()
        for f in function_list:
            if (f.name == extract_nameargs(function_str)[0]):
                return f

    def extract_idx(self, function_str):
        idx = function_str.split(':')[0].strip()
        if self.caller is not None:
            idx = idx.replace('i', self.caller.idx)
        return idx
    
    def decide_argvals(self, function_str):
        argvals = []
        for arg in extract_nameargs(function_str)[1]:
            argval = arg
            if self.caller is not None:
                if 'i.' not in argval:
                    for i, caller_arg in enumerate(self.caller.function.args):
                        if caller_arg == arg:
                            # replace with caller's argument value
                            argval = self.caller.argvals[i]
                        elif ('_' + caller_arg) in arg:
                            # replace T_sigma with something like T_X
                            argval = arg.replace(caller_arg, self.caller.argvals[i])
                else:
                    argval = argval.replace('i', self.caller.idx)
            argvals.append(argval)
        return argvals
    
    # compiles the body of the current line
    def compile_rxns(self, rxn_f): 
        if isinstance(self.function, Instruction):
            caller_idx = self.caller.idx if self.caller is not None else ''
            self.function.write_rxns(self.argvals, caller_idx, self.idx, rxn_f)
        else: # Not base case:
            for function_str in self.function.body:
                line = Line(function_str, self)
                line.compile_rxns(rxn_f)

def extract_nameargs(function_str):
    name_str, arg_str = function_str.strip(')').split('(')
    return name_str, arg_str.replace(' ', '').split(',')

def write_plrxns(poly_list, rxn_f):
    push_temp = '{sigma} + bot_{sigma} -> A_{sigma} + bot_{sigma}'
    pop_temp = 'A_{sigma} + bot_{sigma} -> {sigma} + bot_{sigma}'
    for poly in poly_list:
        rxn_f.write(push_temp.replace('{sigma}', poly) + '\n')
        rxn_f.write(pop_temp.replace('{sigma}', poly) + '\n')

def write_rsrxns(n, rxn_f):
    rxn_temp = 'L_{i} + A_{sigma} -> L_{k} + I_{sigma}'
    for i in range(n):
        rxn = rxn_temp.replace('{i}', str(i+1)).replace('{sigma}', polymer_list[0]).replace('{k}', 'R.1')
        rxn_f.write(rxn + '\n')
    # for handling add-to
    rxn_f.write('L_1.1 + A_' + polymer_list[0] + ' -> L_1.1 + I_' + polymer_list[0] + '\n')

def compile_rxns(main_body, restart_body, fname_str):
    rxn_f = open(fname_str + '_rxns.txt', 'w')
    write_plrxns(polymer_list, rxn_f)
    write_rsrxns(len(main_body), rxn_f)
    for function_str in main_body:
        line = Line(function_str, None)
        line.compile_rxns(rxn_f)
    for i, function_str in enumerate(restart_body):
        line = Line(function_str, None)
        line.compile_rxns(rxn_f)

# generate list of functions
def generate_functions(file):
    protocol_f = open(file, 'r')

    mode = -1
    function_body = []
    main_body = []
    restart_body = []
    lin_str = ''

    for line in protocol_f:
        if (line.startswith('polymers:')):
            mode = 0 # polymers
            lin_str = line.replace('polymers: ', '')
            for p in lin_str.split(','):
                polymer_list.append(p.strip())
        # main protocol
        elif (line == 'main:\n'):
            mode = 1 # in main
        elif (line == 'restart:\n'):
            mode = 2 # in restart
        # internal function
        elif (line.startswith('function: ')):
            mode = 3 # in internal function
            lin_str = line.replace('function: ', '')
        # base instruction
        elif (line.startswith('instruction: ')):
            mode = 4 # in base instructionk
            lin_str = line.replace('instruction: ', '')
        elif (line == '\n'):
            if (mode == 3):
                function_list.append(Function(*extract_nameargs(lin_str.strip()), function_body))
            elif (mode == 4):
                function_list.append(Instruction(*extract_nameargs(lin_str.strip()), function_body))
            mode = -1 #
            function_body = []
        else:
            if (mode == 1):
                main_body.append(line.strip())
            elif (mode == 2):
                restart_body.append(line.strip())
            elif (mode >= 3 ):
                function_body.append(line.strip())
    if (mode == 3):
        function_list.append(Function(*extract_nameargs(lin_str.strip()), function_body))
    if (mode == 4):
        function_list.append(Instruction(*extract_nameargs(lin_str.strip()), function_body))
    return main_body, restart_body

compile_rxns(*generate_functions('example.txt'), 'example')
