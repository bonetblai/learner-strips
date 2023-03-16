from typing import List
import clingo
import sys
import re
from utils.names import *

PRED_1 = ('pred', 1)
STATIC_1 = ('static', 1)
ACTION_1 = ('action', 1)
P_ARITY_2 = ('p_arity', 2)
P_STATIC_1 = ('p_static', 1)
A_ARITY_2 = ('a_arity', 2)
EFF_2 = ('eff', 2)
EFF_3 = ('eff', 3)
PREC_3 = ('prec', 3)
OBJ_2 = ('object', 2)
LABEL_2 = ('label', 2)
LABELNAME_3 = ('labelname', 3)
ACTLABEL_2 = ('actlabel', 2)
ROOT_1 = ('root', 1)
#VAL_1 = ('val', 1)
VAL_2 = ('val', 2)
VAL_3 = ('val', 3)
INV_SCHEMA_3 = ('inv_schema', 3)
UNEQUAL_2 = ('unequal', 2)
FUNCTIONAL_1 = ('functional', 1)
BOOL_1 = ('bool', 1)
FREEVAR_2 = ('freevar', 2)
ARG_2 = ('arg', 2)

def map_to_int(symbols):
    return tuple(map(lambda x: x.number, symbols))

def map_to_number(ints):
    return tuple(map(lambda x: clingo.Number(x), ints))

def pred_1(p):
    return clingo.Function('pred', map_to_number([p]))

def p_arity_2(p, arity):
    return clingo.Function('p_arity', map_to_number([p,arity]))

def p_static_1(p):
    return clingo.Function('p_static', map_to_number([p]))

def action_1(a):
    return clingo.Function('action', map_to_number([a]))

def a_arity_2(a, arity):
    return clingo.Function('a_arity', map_to_number([a, arity]))

def invariant_1(inv):
    return clingo.Function('invariant', map_to_number([inv]))

def inv_schema_3(inv, p, t):
    return clingo.Function('inv_schema', map_to_number([inv, p, t]))

def eff_2(a, p, args):
    return clingo.Function('eff', [clingo.Number(a), clingo.Function('', [clingo.Number(p), map_to_number(args)])])

def eff_3(a, p, args, val):
    return clingo.Function('eff', [clingo.Number(a), clingo.Function('', [clingo.Number(p),clingo.Function('',map_to_number(args))]), clingo.Number(val)])

def prec_3(a, p, args, val):
    return clingo.Function('prec', [clingo.Number(a), clingo.Function('', [clingo.Number(p),clingo.Function('',map_to_number(args))]), clingo.Number(val)])

def unequal_2(a, args):
    return clingo.Function('unequal', [clingo.Number(a), clingo.Function('',map_to_number(args))])

def functional_1(var):
    return clingo.Function('functional', map_to_number([var]))

def bool_1(b):
    return clingo.Function('bool', map_to_number([b]))

def freevar_2(a, v):
    return clingo.Function('freevar', map_to_number([a, v]))

def label_2(a, l):
    return clingo.Function('label', [clingo.Number(a), clingo.String(l)])

def actlabel_2(a, l):
    return clingo.Function('actlabel', map_to_number([a, l]))

def arg_2(a,v):
    return clingo.Function('arg',map_to_number([a,v]))

def str_params(it):
    return ', '.join(('x%d' % (i) if i > 0 else 'none') for i in it)

def str_objs(it):
    return ', '.join('o%d' % (i) for i in it)

class InvariantExp:
    def __init__(self, atoms, fixed):
        self.__atoms = atoms
        self.__fixed = fixed

    def add_atom(self, p, args, func):
        self.__atoms.append((p, args, func))

    def __str__(self):
        l = []
        for p, args, func in self.__atoms:
            if func:
                l.append('%s(%s)=x%s' % (p, str_params(args[:-1]), args[-1]))
            else:
                l.append('%s(%s)' % (p, str_params(args)))
        res = 'exactly-1{%s}' % (', '.join(l))
        if len(self.__fixed):
            res += ' for %s' % str_params(self.__fixed)
        return res


class STRIPSSchema:
    def __init__(self, actions=None, predicates=None, is_static=None, labels=None,
                 precs=None, effs=None, objects=None, invariants=None,
                 root=0, val=None, valstatic=None):
        self._actions = actions if actions != None else {}
        self.__predicates = predicates if predicates != None else {}
        self.__functional = set() # TODO create FunctionalSchema
        self.__bool = set()
        self.__predused = set()
        self.__is_static = is_static if is_static != None else {}
        self._unequal = {}
        self._labels = labels if labels != None else {}
        self._precs = precs if precs != None else {}
        self._effs = effs if effs != None else {}
        self._effs_func = {}
        self.__objects = objects if objects != None else set()
        self.__invariants = invariants if invariants != None else {}
        self.__root = root
        self.__val = val if val != None else []
        self.__valstatic = valstatic if valstatic != None else []
        self._freevar = {}
        self._arg = {}

    def get_actions(self):
        return list(self._actions.keys())

    def add_action(self, action):
        if action not in self._actions:
            self._actions[action] = 3 # default
            self._labels[action] = None
            self._precs[action] = []
            self._effs_func[action] = []
            self._effs[action] = []
            self._unequal[action] = []
            self._freevar[action] = []
            self._arg[action] = []

    def get_action_arity(self, action):
        return self._actions[action]

    def set_action_arity(self, action, arity):
        self.add_action(action)
        self._actions[action] = arity

    def add_predicate(self, predicate):
        if predicate not in self.__predicates:
            self.__predicates[predicate] = 2 # default
            self.__is_static[predicate] = None

    def add_bool(self, b):
        self.__bool.add(b)

    def add_func(self, f):
        self.__functional.add(f)

    def is_static(self, predicate):
        return self.__is_static.get(predicate, None)

    def set_static(self, predicate):
        self.add_predicate(predicate)
        if self.is_static(predicate) == False:
            raise ValueError('Incongruency over static predicate {}'.format(predicate))
        self.__is_static[predicate] = True

    def set_fluent(self, predicate):
        self.add_predicate(predicate)
        if self.is_static(predicate) == True:
            raise ValueError('Incongruency over fluent predicate {}'.format(predicate))
        self.__is_static[predicate] = False

    def get_predicate_arity(self, predicate):
        return self.__predicates[predicate]

    def set_predicate_arity(self, predicate, arity):
        self.add_predicate(predicate)
        self.__predicates[predicate] = arity

    def get_fluents(self):
        ans = []
        for p, is_static in self.__is_static.items():
            if not is_static:
                ans.append(p)
        return ans

    def get_statics(self):
        ans = []
        for p, is_static in self.__is_static.items():
            if is_static:
                ans.append(p)
        return ans

    def get_label(self, action):
        return self._labels[action]

    def set_label(self, action, label):
        self.add_action(action)
        self._labels[action] = label

    def add_prec(self, action, predicate, args, val):
        self.add_action(action)
        self.add_predicate(predicate)
        self.__predused.add(predicate)
        self._precs[action].append((predicate, args, val))

    def add_uneq(self, action, args):
        self.add_action(action)
        self._unequal[action].append(args)

    def add_eff_func(self, action, predicate, args):
        self.add_action(action)
        self.add_predicate(predicate)
        self.add_func(predicate)
        self.__predused.add(predicate)
        self._effs_func[action].append((predicate, args))

    def add_eff(self, action, predicate, args, val):
        self.add_action(action)
        self.add_predicate(predicate)
        self.__predused.add(predicate)
        self._effs[action].append((predicate, args, val))

    def add_object(self, inst, obj): #CHECK
        #print(f'add_object: inst={inst}')
        if obj > 0:
            self.__objects.add(obj)

    def get_root(self):
        return self.__root

    def set_root(self, root):
        self.__root = root

    def add_fluent_true_val(self, inst, pred, args): #CHECK
        #print(f'add_fluent_true_val: inst={inst}')
        self.add_predicate(pred)
        self.set_fluent(pred)
        for o in args:
            self.add_object(inst, o)
        self.__val.append((pred, args))

    def add_static_true_val(self, inst, pred, args): #CHECK
        #print(f'add_static_true_val: inst={inst}')
        self.add_predicate(pred)
        self.set_static(pred)
        for o in args:
            self.add_object(inst, o)
        self.__valstatic.append((pred, args))

    def get_fluent_true_vals(self):
        vals = []
        for p, objs in self.__val:
            arity = self.get_predicate_arity(p)
            vals.append((p, objs[:arity]))
        return vals

    def get_static_true_vals(self):
        vals = []
        for p, objs in self.__valstatic:
            arity = self.get_predicate_arity(p)
            vals.append((p, objs[:arity]))
        return vals

    def add_invariant(self, inv):
        if inv not in self.__invariants:
            self.__invariants[inv] = []

    def add_pred_to_invariant(self, inv, pred, t):
        self.add_invariant(inv)
        self.__predused.add(pred)
        self.__invariants[inv].append((pred, t))

    def get_used(self):
        return self.__predused

    def add_freevar(self, a, v):
        self.add_action(a)
        self._freevar[a].append(v)

    def add_arg(self, a, v):
        self.add_action(a)
        self._arg[a].append(v)

    @classmethod
    def create_from_clingo(cls, symbols : List[clingo.Symbol]):
        others = []
        schema = cls()
        vals = []
        for symbol in symbols:
            if symbol.type != clingo.SymbolType.Function:
                others.append(symbol)
            if symbol.match(*PRED_1):
                p, = symbol.arguments
                schema.add_predicate(p.number)
            elif symbol.match(*P_STATIC_1):
                p, = symbol.arguments
                schema.set_static(p.number)
            elif symbol.match(*ACTION_1):
                a, = symbol.arguments
                schema.add_action(a.number)
            elif symbol.match(*LABEL_2):
                a, l = symbol.arguments
                label = str(l)
                assert label[0] == '"' and label[-1] == '"'
                label = label[1:-1]
                schema.set_label(a.number, label)
            elif symbol.match(*LABELNAME_3):
                inst, a, l = symbol.arguments
            elif symbol.match(*P_ARITY_2):
                p, arity = symbol.arguments
                schema.set_predicate_arity(p.number, arity.number)
            elif symbol.match(*A_ARITY_2):
                a, arity = symbol.arguments
                schema.set_action_arity(a.number, arity.number)
            elif symbol.match(*EFF_2):
                a, m = symbol.arguments
                p, args = m.arguments
                schema.add_eff_func(a.number, p.number, map_to_int(args.arguments))
            elif symbol.match(*EFF_3):
                a, m, val = symbol.arguments
                p, args = m.arguments
                schema.add_eff(a.number, p.number, map_to_int(args.arguments), val.number)
            elif symbol.match(*PREC_3):
                a, m, val = symbol.arguments
                p, args = m.arguments
                schema.add_prec(a.number, p.number, map_to_int(args.arguments), val.number)
            elif symbol.match(*UNEQUAL_2):
                a, args = symbol.arguments
                schema.add_uneq(a.number, map_to_int(args.arguments))
            elif symbol.match(*OBJ_2):
                inst, o = symbol.arguments
                schema.add_object(inst.number, o.number)
            elif symbol.match(*ROOT_1):
                s, = symbol.arguments
                schema.set_root(s.number)
            elif symbol.match(*VAL_2):
                inst, k = symbol.arguments
                p, objs = k.arguments
                schema.add_static_true_val(inst.number, p.number, map_to_int(objs.arguments))
            elif symbol.match(*VAL_3):
                inst, k, s = symbol.arguments
                p, objs = k.arguments
                vals.append((inst.number, (p.number, map_to_int(objs.arguments)), s.number))
            elif symbol.match(*INV_SCHEMA_3):
                inv, p, t = symbol.arguments
                schema.add_pred_to_invariant(inv.number, p.number, t.number)
            elif symbol.match(*BOOL_1):
                p, = symbol.arguments
                schema.add_bool(p.number)
            elif symbol.match(*FUNCTIONAL_1):
                p, = symbol.arguments
                schema.add_func(p.number)
            elif symbol.match(*FREEVAR_2):
                a, v = symbol.arguments
                schema.add_freevar(a.number, v.number)
            elif symbol.match(*ARG_2):
                a, v = symbol.arguments
                schema.add_arg(a.number, v.number)

        root = schema.get_root()
        for inst, (p, oo), s in vals:
            if s == root:
                schema.add_fluent_true_val(inst, p, oo)

        return schema

    def get_string_action(self, action):
        out = ''
        out += '\ta%s(%s) ' % (action, str_params(range(1,self.get_action_arity(action)+1)))
        label = self.get_label(action)
        if label != None:
            out += 'label %s' % (label)
        out += '\n'
        str_stat = []
        str_fluent = []
        for args in self._unequal[action]:
            str_stat.append('neq(%s)' % str_params(args))
        for p, args, val in self._precs[action]:
            sign = '' if val else '-'
            arity = self.get_predicate_arity(p)
            if p in self.__functional:
                str_fluent.append('p%d(%s)%s=%s' % (p, str_params(args[:arity-1]),'' if val else '!',str_params([args[-1]])))
            else:
                str_fluent.append(sign + 'p%s(%s)' % (p, str_params(args[:arity])))
        out += '\t\tStatic: %s\n' % (', '.join(str_stat))
        out += '\t\tPre: %s\n' % (', '.join(str_fluent))
        str_eff = []
        for p, args in self._effs_func[action]:
            arity = self.get_predicate_arity(p)
            str_eff.append('p%d(%s):=%s' % (p, str_params(args[:arity-1]),str_params([args[-1]])))
        for p, args, val in self._effs[action]:
            sign = '' if val else '-'
            arity = self.get_predicate_arity(p)
            if p in self.__functional:
                str_eff.append('p%d(%s)%s=%s' % (p, str_params(args[:arity-1]),':' if val else '!',str_params([args[-1]])))
            else:
                str_eff.append(sign + 'p%s(%s)' % (p, str_params(args[:arity])))
        out += '\t\tEff: %s\n' % (', '.join(str_eff))
        return out

    def get_string_invariant_func(self, pred):
        arity = self.get_predicate_arity(pred)
        d = {}
        for (p, t) in self.__invariants[pred]:
            if t not in d:
                d[t] = InvariantExp([('p%s'%pred,list(range(1,1+arity)),True)], [3-t] if arity==2 else [])
            if not self.is_static(pred):
                d[t].add_atom('p%s'%p, [t] if arity > 1 else [], False)
        s = ''
        for inv in d.values():
            s += '\t%s\n' % (str(inv))
        return s

    def get_string_invariant(self, inv):
        out = 'Invariant {}: '.format(inv)
        invstr = []
        for p, t in self.__invariants[inv]:
            if t == 0:
                invstr.append('unary1(p{})'.format(p))
            if t == 1:
                invstr.append('unary2(p{})'.format(p))
            if t == 2:
                invstr.append('binary1(p{})'.format(p))
            if t == 3:
                invstr.append('binary2(p{})'.format(p))
        out += ' '.join(invstr) + '\n'
        return out

    def get_string(self, val=False):
        out = ''
        if len(self.__functional) == 0:
            out += 'Fluent Predicates:\n'
            for p in sorted(self.get_fluents()):
                if p in self.__predused:
                    arity = self.get_predicate_arity(p)
                    out += '\tp%s(%s)\n' % (p, str_params(range(1,arity+1)))
            out += 'Static Predicates:\n'
            for p in sorted(self.get_statics()):
                if p in self.__predused:
                    arity = self.get_predicate_arity(p)
                    out += '\tp%s(%s)\n' % (p, str_params(range(1,arity+1)))
        else:
            out += 'Fluent Predicates:\n'
            for p in sorted(self.get_fluents()):
                if p in self.__predused:
                    arity = self.get_predicate_arity(p)
                    if p in self.__functional:
                        out += '\tp%s(%s), functional\n' % (p, str_params(range(1,arity)))
                    else:
                        out += '\tp%s(%s), bool\n' % (p, str_params(range(1,arity+1)))
            out += 'Static Predicates:\n'
            for p in sorted(self.get_statics()):
                if p in self.__predused:
                    arity = self.get_predicate_arity(p)
                    if p in self.__functional:
                        out += '\tp%s(%s), functional\n' % (p, str_params(range(1,arity)))
                    else:
                        out += '\tp%s(%s), bool\n' % (p, str_params(range(1,arity+1)))
            out += 'Invariants:\n'
            for p in sorted(self.__predicates.keys()):
                if len(self.__invariants.get(p,[])) > 0:
                    out += self.get_string_invariant_func(p)
        
        if len(self.__functional) == 0:
            for inv in sorted(self.__invariants):
                out += self.get_string_invariant(inv)
        out += 'Actions:\n'
        print(self._actions)
        for a in sorted(self._actions):
            out += self.get_string_action(a)
        out += 'Objects: %s\n' % (', '.join('o%d' % (o) for o in sorted(self.__objects)))
        if val:
            out += 'Static atoms: %s\n' % (', '.join('p%d(%s)' % (p, str_objs(objs)) for p, objs in self.get_static_true_vals()))
            root = self.get_root()
            out += 'Root %d valuation: %s' % (root, ', '.join('p%d(%s)' % (p, str_objs(objs)) for p, objs in self.get_fluent_true_vals()))
        return out

    def get_schema(self):
        symbols = []
        for action, arity in self._actions.items():
            symbols.append(label_2(action, self._labels[action]))
            symbols.append(action_1(action))
            symbols.append(a_arity_2(action, arity))
            for v in self._arg[action]:
                symbols.append(arg_2(action, v))
            for v in self._freevar[action]:
                symbols.append(freevar_2(action, v))
            for p, arg in self._effs_func[action]:
                symbols.append(eff_2(action, p, arg))
            for p, arg, val in self._effs[action]:
                symbols.append(eff_3(action, p, arg, val))
            for p, arg, val in self._precs[action]:
                symbols.append(prec_3(action, p, arg, val))
            for args in self._unequal[action]:
                symbols.append(unequal_2(action, args))
        for pred, arity in self.__predicates.items():
            symbols.append(pred_1(pred))
            if pred in self.__predused:
                symbols.append(p_arity_2(pred, arity))
                if self.is_static(pred):
                    symbols.append(p_static_1(pred))
                if pred in self.__bool:
                    symbols.append(bool_1(pred))
                elif pred in self.__functional:
                    symbols.append(functional_1(pred))

        for inv in self.__invariants:
            symbols.append(invariant_1(inv))
            for p, t in self.__invariants[inv]:
                symbols.append(inv_schema_3(inv, p, t))
        
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

class FunctionalSchema(STRIPSSchema):
    def __init__(self):
        super().__init__()
        self.__actlabel = {}

    def add_action(self, action):
        if action not in self._actions:
            self._actions[action] = 3 # default
            self.__actlabel[action] = None
            self._precs[action] = []
            self._effs_func[action] = []
            self._effs[action] = []
            self._unequal[action] = []
            self._freevar[action] = []
            self._arg[action] = []

    def set_label(self, label, name):
        self._labels[label] = name
            
    def set_act_label(self, action, label):
        self.__actlabel[action] = label

    def get_act_label(self, action):
        return self.__actlabel[action]

    def get_label(self, action):
        return super().get_label(self.__actlabel[action])

    @classmethod
    def create_from_clingo(cls, symbols : List[clingo.Symbol]):
        schema = super().create_from_clingo(symbols)
        for s in symbols:
            if s.match(*ACTLABEL_2):
                a, l = map_to_int(s.arguments)
                schema.set_act_label(a, l)
        return schema

    def get_schema(self):
        symbols = super().get_schema()
        for action in self.get_actions():
            symbols.append(actlabel_2(action, self.get_act_label(action)))
        return symbols

if __name__ == '__main__':
    filestr = ''.join(sys.stdin.readlines())
    result = parse_clingo_out(filestr)
    #print('INPUT:')
    #print(filestr)
    #print('='*80)
    print(STRIPSSchema.create_from_clingo(result[SYMBOLS]).get_string())
    #for key, value in result.items():
    #    print(key, value)
