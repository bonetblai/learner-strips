from typing import List
import clingo
import sys
import re
from utils.names import *
#from names import *

VARTYPE_1 = ('vartype', 1)
CONSTANT_1 = ('constant', 1)
ACTION_1 = ('action', 1)
AFFECTS_2 = ('affects', 2)
TEMPLATE_2 = ('template', 2)
INSTANTIATION_3 = ('instantiation', 3)
IRREFLEXIVE_1 = ('irreflexive', 1)
USING_PRECONDITIONS_0 = ('using_preconditions', 0)
FORBID_CONSTANTS_IN_BXX_0 = ('forbid_constants_in_Bxx', 0)

LABEL_1 = ('label', 1)
LABEL_2 = ('label', 2)
LABELNAME_3 = ('labelname', 3)

OBJECT_2 = ('object', 2)
DOMAIN_3 = ('domain', 3)
SIMPLE_3 = ('simple', 3)
CONSTANT_3 = ('constant', 3)
CONSTANT_TYPE_3 = ('constant_type', 3)
ALL_IN_1 = ('all_in', 1)
ALL_OUT_1 = ('all_out', 1)
DTG_4 = ('dtg', 4)
DTG_5 = ('dtg', 5)
PREC_4 = ('prec', 4)
VAR_3 = ('var', 3)
VAL_5 = ('val', 5)

LABELMATCH_3 = ('labelmatch', 3)
APPL_4 = ('appl', 4)
NEXT_5 = ('next', 5)


def map_to_int(symbols):
    return tuple(map(lambda x: x.number, symbols))

def map_to_number(ints):
    return tuple(map(lambda x: clingo.Number(x), ints))

def str_params(it):
    return ', '.join(('x%d' % (i) if i > 0 else 'none') for i in it)

def str_objs(it):
    return ', '.join('o%d' % (i) for i in it)


def vartype_1(t):
    return clingo.Function('vartype', [clingo.Number(t)])

def action_1(a):
    return clingo.Function('action', [clingo.Number(a)])

def var_3(inst, x, t):
    return clingo.Function('var', map_to_number([inst, x, t]))

def constant_1(c):
    return clingo.Function('constant', [clingo.Number(c)])

def constant_3(inst, c ,o):
    return clingo.Function('constant', map_to_number([inst, c, o]))

def constant_type_3(inst, t ,c):
    return clingo.Function('constant_type', map_to_number([inst, t, c]))

def all_in_1(c):
    return clingo.Function('all_in', [clingo.Number(c)])

def all_out_1(c):
    return clingo.Function('all_out', [clingo.Number(c)])

def val_5(inst, t, x, v, s):
    return clingo.Function('val', map_to_number([inst, t, x, v, s]))

def dtg_45(inst, t, v1, v2, a):
    if a == 0:
        return clingo.Function('dtg', map_to_number([inst, t, v1, v2]))
    else:
        return clingo.Function('dtg', map_to_number([inst, t, v1, v2, a]))

def label_1(l):
    return clingo.Function('label', [clingo.String(l)])

def label_2(a, l):
    return clingo.Function('label', [clingo.Number(a), clingo.String(l)])

def labelmatch_3(inst, a, b):
    return clingo.Function('labelmatch', map_to_number([inst, a, b]))

def affects_2(a, t):
    return clingo.Function('affects', [clingo.Number(a), clingo.Number(t)])

def template_2(a, at):
    return clingo.Function('template', [clingo.Number(a), clingo.String(at)])

def instantiation_3(a, at, tt):
    return clingo.Function('instantiation', [clingo.Number(a), clingo.String(at), clingo.Function('', map_to_number(tt))])

def appl_4(inst, a, vt, s):
    return clingo.Function('appl', map_to_number([inst, at]) + [map_to_number(vt)] + [clingo.Number(s)])

def next_5(inst, a, oo, s1, s2):
    return clingo.Function('next', map_to_number([inst, at]) + [map_to_number(oo)] + map_to_number([s1, s2]))

def irreflexive_1(t):
    return clingo.Function('irreflexive', [clingo.Number(t)])

def using_preconditions_0():
    return clingo.Function('using_preconditions', [])

def forbid_constants_in_Bxx_0():
    return clingo.Function('forbid_constants_in_Bxx', [])

def simple_3(inst, t, v):
    return clingo.Function('simple', map_to_number([inst, t, v]))


class STRIPSSchema:
    def __init__(self):
        self._vartypes = set()
        self._constants = {}
        self._all_in = set()
        self._all_out = set()
        self._actions = {}

        self._instances = set()
        self._objects = {}
        self._domains = {}
        self._dtgs = {}
        self._precs = {}
        self._vars = {}

        self._roots = {}
        self._appl = {}
        self._next = {}

        self._irreflexive = set()
        self._using_preconditions = False
        self._forbid_constants_in_Bxx = False
        self._simple = {}

    def add_vartype(self, t):
        self._vartypes.add(t)

    def add_constant(self, c):
        if c not in self._constants:
            self._constants[c] = {}

    def add_constant_all_in(self, c):
        self.add_constant(c)
        self._all_in.add(c)

    def add_constant_all_out(self, c):
        self.add_constant(c)
        self._all_out.add(c)

    def add_action(self, action):
        if action not in self._actions:
            self._actions[action] = dict(label=None, arity=-1, type_affected_var=None, template=None, args=None, labelmatch=None)

    def set_label(self, action, label):
        #print(f'schema:label({action},{label})')
        if action not in self._actions: self.add_action(action)
        assert self._actions[action]['label'] is None
        self._actions[action]['label'] = label

    def set_type_affected_var(self, action, t):
        #print(f'schema:affected({action},{t})')
        if action not in self._actions: self.add_action(action)
        assert self._actions[action]['type_affected_var'] is None
        self._actions[action]['type_affected_var'] = t

    def set_template(self, action, at_label):
        #print(f'schema:template({action},"{at_label}")')
        if action not in self._actions: self.add_action(action)
        assert self._actions[action]['template'] is None
        self._actions[action]['template'] = at_label

    def set_instantiation(self, action, at_label, args):
        #print(f'schema:instantiation({action},"{at_label}",{args})')
        if action not in self._actions: self.add_action(action)
        if self._actions[action]['template'] is None: self.set_template(action, at_label)
        assert self._actions[action]['template'] == at_label and self._actions[action]['args'] is None
        self._actions[action]['args'] = args
        #for t in args:
        #    if t not in self._vartypes:
        #        print(f"Error: type {t} in instantiation {args} for action {action} doesn't exist")


    def add_instance(self, inst):
        self._instances.add(inst)

    def add_object(self, inst, obj):
        #print(f'schema:object({inst},{obj})')
        if obj > 0:
            self.add_instance(inst)
            if inst not in self._objects:
                self._objects[inst] = set()
            self._objects[inst].add(obj)

    def add_domain_value(self, inst, t, v):
        #print(f'schema:domain({inst},{t},{v})')
        assert t in self._vartypes
        self.add_instance(inst)
        if inst not in self._domains: self._domains[inst] = {}
        if t not in self._domains[inst]: self._domains[inst][t] = set()
        self._domains[inst][t].add(v)

    def set_irreflexive(self, t):
        #print(f'schema:irreflexive({t})')
        self._irreflexive.add(t)

    def set_using_preconditions(self):
        self._using_preconditions = True

    def set_forbid_constants_in_Bxx(self):
        self._forbid_constants_in_Bxx = True

    def set_simple(self, inst, t, v):
        #print(f'schema:simple({inst},{t},{v})')
        if inst not in self._simple: self._simple[inst] = {}
        if t not in self._simple[inst]: self._simple[inst][t] = set()
        self._simple[inst][t].add(v)

    def add_constant_denotation(self, inst, c, v):
        #print(f'schema:constant({inst},{c},{v})')
        self.add_instance(inst)
        if c not in self._constants: self.add_constant(c)
        assert inst not in self._constants[c] or self._constants[c][inst]['value'] is None
        if inst not in self._constants[c]:
            self._constants[c][inst] = dict(value=v, type=None)
        else:
            self._constants[c][inst]['value'] = v

    def add_constant_type(self, inst, t, c):
        #print(f'schema:constant_type({inst},{t},{c})')
        self.add_instance(inst)
        if c not in self._constants: self.add_constant(c)
        assert inst not in self._constants[c] or self._constants[c][inst]['type'] is None, f'inst={inst}, t={t}, c={c}, x={self._constants}'
        if inst not in self._constants[c]:
            self._constants[c][inst] = dict(value=None, type=t)
        else:
            self._constants[c][inst]['type'] = t

    def add_dtg_edge(self, inst, t, v1, v2, a):
        #print(f'schema:dtg({inst},{t},{v1},{v2},{a})')
        assert t in self._vartypes
        self.add_instance(inst)
        if inst not in self._dtgs: self._dtgs[inst] = {}
        if t not in self._dtgs[inst]: self._dtgs[inst][t] = set()
        self._dtgs[inst][t].add((v1, v2, a))

    def add_prec(self, inst, t, x, v):
        #print(f'schema:prec({inst},{t},{x},{v})')
        assert t in self._vartypes
        self.add_instance(inst)
        if inst not in self._precs: self._precs[inst] = {}
        if t not in self._precs[inst]: self._precs[inst][t] = set()
        self._precs[inst][t].add((x, v))

    def add_variable(self, inst, x, t):
        #print(f'schema:var({inst},{x},{t})')
        assert t in self._vartypes
        self.add_instance(inst)
        if inst not in self._vars: self._vars[inst] = {}
        if t not in self._vars[inst]: self._vars[inst][t] = set()
        self._vars[inst][t].add(x)

    def add_appl(self, inst, action, args, s):
        #print(f'schema:appl({inst},{action},{args},{s})')
        if inst not in self._appl: self._appl[inst] = {}
        if action not in self._appl[inst]: self._appl[inst][action] = set()
        self._appl[inst][action].add((s, args))

    def add_next(self, inst, action, args, s1, s2):
        #print(f'schema:next({inst},{action},{args},{s1},{s2})')
        if inst not in self._next: self._next[inst] = {}
        if action not in self._next[inst]: self._next[inst][action] = {}
        assert (s1, args) not in self._next[inst][action]
        self._next[inst][action][(s1, args)] = s2

    def get_root(self): #CHECK
        return self.__root

    def set_root(self, root): #CHECK
        self.__root = root

    @classmethod
    def create_from_clingo(cls, symbols : List[clingo.Symbol]):
        schema = cls()
        others = []
        vals = []
        for symbol in symbols:
            if symbol.type != clingo.SymbolType.Function:
                others.append(symbol)
                print(f'unknown=|{symbol}|')
            if symbol.match(*VARTYPE_1):
                t, = symbol.arguments
                schema.add_vartype(t.number)
            elif symbol.match(*ACTION_1):
                a, = symbol.arguments
                schema.add_action(a.number)
            elif symbol.match(*OBJECT_2):
                inst, o = symbol.arguments
                schema.add_object(inst.number, o.number)
            elif symbol.match(*VAR_3):
                inst, x, t = symbol.arguments
                schema.add_variable(inst.number, x.number, t.number)
            elif symbol.match(*DOMAIN_3):
                inst, t, v = symbol.arguments
                schema.add_domain_value(inst.number, t.number, v.number)
            elif symbol.match(*CONSTANT_1):
                c, = symbol.arguments
                schema.add_constant(c.number)
            elif symbol.match(*CONSTANT_3):
                inst, c, v = symbol.arguments
                schema.add_constant_denotation(inst.number, c.number, v.number)
            elif symbol.match(*CONSTANT_TYPE_3):
                inst, t, c = symbol.arguments
                schema.add_constant_type(inst.number, t.number, c.number)
            elif symbol.match(*ALL_IN_1):
                c, = symbol.arguments
                schema.add_constant_all_in(c.number)
            elif symbol.match(*ALL_OUT_1):
                c, = symbol.arguments
                schema.add_constant_all_out(c.number)
            elif symbol.match(*VAL_5):
                inst, t, x, v, s = symbol.arguments
                vals.append((inst.number, t.number, x.number, v.number, s.number))
                #CHECK
                #CHECK print(f'val({inst.number},{t.number},{x.number},{v.number},{s.number})')
            elif symbol.match(*DTG_4):
                inst, t, v1, v2 = symbol.arguments
                schema.add_dtg_edge(inst.number, t.number, v1.number, v2.number, 0)
            elif symbol.match(*DTG_5):
                inst, t, v1, v2, a = symbol.arguments
                schema.add_dtg_edge(inst.number, t.number, v1.number, v2.number, a.number)
            elif symbol.match(*PREC_4):
                inst, t, x, v = symbol.arguments
                schema.add_prec(inst.number, t.number, x.number, v.number)
            elif symbol.match(*LABELNAME_3):
                inst, a, l = symbol.arguments
                label = str(l)
                assert label[0] == '"' and label[-1] == '"'
                label = label[1:-1]
                #CHECK
                #CHECK print(f'labelname({inst.number},{a.number},{label})')
            elif symbol.match(*LABEL_1):
                l, = symbol.arguments
                label = str(l)
                assert label[0] == '"' and label[-1] == '"'
                label = label[1:-1]
                #CHECK
                #CHECK print(f'label({label})')
            elif symbol.match(*LABEL_2):
                a, l = symbol.arguments
                label = str(l)
                assert label[0] == '"' and label[-1] == '"'
                label = label[1:-1]
                schema.set_label(a.number, label)
            elif symbol.match(*LABELMATCH_3):
                inst, a, b = symbol.arguments
                #CHECK schema_set_labelmatch({inst.number},{a.number},{b.number})
                #CHECK print(f'labelmatch({inst.number},{a.number},{b.number})')
            elif symbol.match(*AFFECTS_2):
                a, t = symbol.arguments
                schema.set_type_affected_var(a.number, t.number)
            elif symbol.match(*TEMPLATE_2):
                a, at = symbol.arguments
                at_label = str(at)
                assert at_label[0] == '"' and at_label[-1] == '"'
                at_label = at_label[1:-1]
                schema.set_template(a.number, at_label)
            elif symbol.match(*INSTANTIATION_3):
                a, at, tt = symbol.arguments
                args = tt.arguments
                at_label = str(at)
                assert at_label[0] == '"' and at_label[-1] == '"'
                at_label = at_label[1:-1]
                schema.set_instantiation(a.number, at_label, map_to_int(args))
            elif symbol.match(*APPL_4):
                inst, a, oo, s = symbol.arguments
                args = oo.arguments
                schema.add_appl(inst.number, a.number, map_to_int(args), s.number)
            elif symbol.match(*NEXT_5):
                inst, a, oo, s1, s2 = symbol.arguments
                args = oo.arguments
                schema.add_next(inst.number, a.number, map_to_int(args), s1.number, s2.number)
            elif symbol.match(*IRREFLEXIVE_1):
                t, = symbol.arguments
                schema.set_irreflexive(t.number)
            elif symbol.match(*USING_PRECONDITIONS_0):
                schema.set_using_preconditions()
            elif symbol.match(*FORBID_CONSTANTS_IN_BXX_0):
                schema.set_forbid_constants_in_Bxx()
            elif symbol.match(*SIMPLE_3):
                inst, t, v = symbol.arguments
                schema.set_simple(inst.number, t.number, v.number)

        #root = schema.get_root()
        #print(f'root={root}')
        #for inst, (p, oo), s in vals:
        #    if s == root:
        #        schema.add_fluent_true_val(inst, p, oo)

        return schema

    def get_string(self, val=False):
        action_args = dict(A0='(x)', A1='(x)', A2='(x,v2)', A0x='(x,y)', A1x='(x,y)', A2x='(x,v2,y)', B0='(x,v2)', B0x='(x,v2,y)', C0='(x,z)', C0x='(x,z,y)')
        out = ''
        out += f"Types: {', '.join(map(lambda t: 'T' + str(t), self._vartypes))}\n"
        if len(self._irreflexive):
            out += f"Irreflexive: {', '.join(map(lambda t: 'T' + str(t), self._irreflexive))}\n"
        if len(self._constants):
            out += f"Constants:"
            for c in self._constants:
                if not out.endswith('Constants:'): out += ','
                out += f' C{c}'
                if c in self._all_in or c in self._all_out:
                    out += '('
                    if c in self._all_in:
                        out += 'all_in'
                        if c in self._all_out: out += ','
                    if c in self._all_out:
                        out += 'all_out'
                    out += ')'
            out += '\n'
        for action in sorted(self._actions.keys()):
            template = self._actions[action]['template']
            args = self._actions[action]['args']
            label = self._actions[action]['label']
            type_affected_var = self._actions[action]['type_affected_var']
            out += f'Action {action}:'
            out += f" label={label}, vartype=T{type_affected_var}, template={template}[{','.join(map(str, args))}]{action_args[template]}"
            out += '\n'
        out += f'using_preconditions={self._using_preconditions}\n'
        out += f'forbid_constants_in_Bxx={self._forbid_constants_in_Bxx}\n'
        out += '\n'

        # instances
        for inst in sorted(self._instances):
            out += f'Instance {inst}:\n'
            out += f'  #objs = {len(self._objects[inst])}\n'

            for t in self._vartypes:
                if t in self._domains[inst]:
                    out += f"  Dom[T{t}] = [{', '.join(map(str, self._domains[inst][t]))}]"
                else:
                    out += f'  Dom[T{t}] is undefined'
                out += '\n'

            for c in self._constants.keys():
                out += f"  Constant C{c}:"
                if inst in self._constants[c]:
                    value = self._constants[c][inst]['value']
                    t = self._constants[c][inst]['type']
                    out += f' value={value}, type=T{t}'
                else:
                    out += ' not defined'
                out += '\n'

            assert inst in self._dtgs
            for t in self._vartypes:
                if t in self._dtgs[inst]:
                    out += f"  DTG[T{t}] = [{', '.join(map(str, self._dtgs[inst][t]))}]"
                else:
                    out += f'  DTG[T{t}] is undefined'
                out += '\n'

            if inst in self._precs:
                for t in self._precs[inst]:
                    out += f"  PREC[T{t}] = [{', '.join(map(str, self._precs[inst][t]))}]"
                    out += '\n'

            for t in self._vartypes:
                if t in self._vars[inst]:
                    out += f"  Var[T{t}] = [{', '.join(map(str, self._vars[inst][t]))}]"
                    out += '\n'

            if inst in self._simple:
                for t in self._simple[inst]:
                    out += f"  Simple[T{t}] = [{', '.join(map(str, self._simple[inst][t]))}]"
                    out += '\n'
            out += '\n'

            out += '  Groundings:\n'
            for action in self._actions.keys():
                assert inst in self._appl and inst in self._next
                assert action not in self._appl[inst] or action in self._next[inst]
                assert action not in self._next[inst] or action in self._appl[inst]
                if action in self._appl[inst]:
                    transitions = []
                    for (s1, args) in sorted(self._appl[inst][action]):
                        assert (s1, args) in self._next[inst][action]
                        s2 = self._next[inst][action][(s1, args)]
                        transitions.append(f"(s1={s1},[{','.join(map(str, args))}],s2={s2})")
                    out += f"    Action {action} ({len(transitions)}): {', '.join(transitions)}\n"
       
        return out

    def get_schema(self):
        symbols = []
        for t in self._vartypes:
            symbols.append(vartype_1(t))
        for c in self._constants:
            symbols.append(constant_1(c))
            if c in self._all_in:
                symbols.append(all_in_1(c))
            if c in self._all_out:
                symbols.append(all_out_1(c))
        for action in self._actions:
            template = self._actions[action]['template']
            label = self._actions[action]['label']
            type_affected_var = self._actions[action]['type_affected_var']
            args = self._actions[action]['args']
            symbols.append(action_1(action))
            symbols.append(label_2(action, label))
            symbols.append(affects_2(action, type_affected_var))
            symbols.append(template_2(action, template))
            symbols.append(instantiation_3(action, template, args))
        for t in self._irreflexive:
            symbols.append(irreflexive_1(t))
        if self._using_preconditions:
            symbols.append(using_preconditions_0())
        if self._forbid_constants_in_Bxx:
            symbols.append(forbid_constants_in_Bxx_0())
        return symbols


def is_sat(string):
    if string == 'SATISFIABLE':
        return True
    if string == 'UNSATISFIABLE':
        return False
    if string == 'UNKNOWN':
        return None
    raise ValueError('Unexpected result: %s'%(string))

def parse_clingo_string(program):
    ctl = clingo.Control()
    ctl.add('base', [], program)
    ctl.ground([('base', [])])
    return list(s.symbol for s in ctl.symbolic_atoms if s.is_fact)

def parse_clingo_out(output, firstmodel=False):
    results = {}
    empty = r'\s.*?'
    answer = r'Answer: \d+?\n(.*?)\n'
    matches = re.finditer(answer, output)
    models = list(matches)
    if len(models):
        results[SAT] = True
        if firstmodel:
            text = models[0].group(1)
        else:
            text = models[-1].group(1)
        text = ''.join(w+'.' for w in text.split())
        results[SYMBOLS] = parse_clingo_string(text)
    else:
        if output.find('UNSATISFIABLE') > -1:
            results[SAT] = False
        else:
            results[SAT] = None
    if SYMBOLS not in results:
        results[SYMBOLS] = None
    models = f'Models{empty}: (.*?)\n'
    calls = f'Calls{empty}: (.*?)\n'
    time = f'Time{empty}: (.*?)s '+r'\(Solving: (.*?)s 1st Model: (.*?)s Unsat: (.*?)s\)\n'
    cputime = f'CPU Time{empty}: (.*?)s\n'
    choices = f'Choices{empty}: (\d+)'
    conflicts = f'Conflicts{empty}: (\d+)'
    rules = f'Rules{empty}: (\d+)'
    variables = f'Variables{empty}: (\d+)'
    constraints = f'Constraints{empty}: (\d+)'
    optimization = f'Models{empty}: .*?\n{empty}Optimum{empty}: (.*?)\nOptimization{empty}: (.*?)\n'
    keys_calls = [
        (models, [(MODELS, None)]),
        (calls, [(CALLS, int)]),
        (time, [(TIME, float),(SOLVING, float), (MODEL1st, float), (TIMEUNSAT, float)]),
        (cputime, [(CPUTIME, float)]),
        (choices, [(CHOICES, int)]),
        (conflicts, [(CONFLICTS, int)]),
        (rules, [(RULES, int)]),
        (variables, [(VARIABLES, int)]),
        (constraints, [(CONSTRAINTS, int)])]
    
    for regex, groups in keys_calls:
        match = re.search(regex, output)
        if match == None:
            continue
        for (key, call), value in zip(groups, match.groups()):
            results[key] = call(value) if call != None else value
    
    optmatch = re.search(optimization, output)
    if optmatch != None:
        results[OPTIMUM] = optmatch.group(1)
        results[OPTIMIZATION] = optmatch.group(2)
    else:
        results[OPTIMUM] = None
        results[OPTIMIZATION] = None
    return results


if __name__ == '__main__':
    filestr = ''.join(sys.stdin.readlines())
    clingo_result = parse_clingo_out(filestr)
    schema = STRIPSSchema.create_from_clingo(clingo_result[SYMBOLS])
    model = schema.get_schema()
    print(schema.get_string())
    for symbol in model:
        print(f'{symbol}.')

