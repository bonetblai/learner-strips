import sys, os, shutil, signal, subprocess, resource, subprocess, argparse
from termcolor import colored
from pathlib import Path
from typing import List
import logging
import re

from utils.dfa          import DFA
from utils.names        import SAT, VARIABLES, CONSTRAINTS, SYMBOLS, RULES, SOLVING, CONFLICTS, CHOICES, TIME, MODEL1st
from utils.output       import parse_clingo_out as parse_clingo_out_orig
from utils.output       import STRIPSSchema     as STRIPSSchema_orig
from utils.output_mf    import parse_clingo_out as parse_clingo_out_mf
from utils.output_mf    import STRIPSSchema     as STRIPSSchema_mf
from utils.output_vars2 import parse_clingo_out as parse_clingo_out_vars
from utils.output_vars2 import STRIPSSchema     as STRIPSSchema_vars

class Benchmark:
    def __init__(self, fields):
        index = 0

        self.samples = []
        while fields[index].endswith('.lp'):
            self.samples.append(fields[index])
            index += 1

        self.num_objs_learn = int(fields[index]); index += 1
        self.max_true_atoms_learn = int(fields[index]); index += 1
        self.options = []
        while fields[index].upper() != 'VERIFY':
            self.options.append(fields[index]); index += 1

        index += 1
        self.max_num_objs_ver = int(fields[index]); index += 1
        self.max_true_atoms_ver = int(fields[index]); index += 1
        self.verify = []
        while index < len(fields) and fields[index] != 'partial':
            self.verify.append(fields[index]); index += 1
        self.partial = None
        if index < len(fields) and fields[index] == 'partial':
            index += 1
            self.partial = int(fields[index])

    def get_options(self):
        opts = self.options
        if self.partial != None:
            opts.append('-c perc={} -c perc_nodes=1'.format(self.partial))
        return ' '.join(opts)

class Stats:
    def __init__(self, samples):
        self.data = dict(samples = '_'.join(samples))
    def __str__(self):
        str_template = '{samples} {objs} {num_rules} {num_variables} {num_constraints} {num_choices} {num_conflicts} {total_time} {solve_time} {first} {solve_memory} {satisfiable}'
        str_subtemplate = '{instance} {objs} {num_rules} {num_variables} {num_constraints} {num_choices} {num_conflicts} {total_time} {solve_time} {solve_memory} {satisfiable}'
        if 'num_rules' in self.data:
            # if only verification, self.data does not contain info about solving theory
            as_str = str_template.format(**self.data)
        else:
            as_str = '{samples}'.format(**self.data)

        success = self.data['satisfiable']
        if success:
            ver_instance = None
            for ver in self.data.get('verify', []):
                #if ver['satisfiable']:
                #    as_str += ' ' + str_subtemplate.format(**ver)
                if ver_instance != ver['instance'] and not success: break
                as_str += ' ' + str_subtemplate.format(**ver)
                success = ver['satisfiable']
                ver_instance = ver['instance']
        as_str += ' Success' if success else ' Failure'

        return as_str

# clingo scripts
g_clingo = {
    'mf'   : { 'solve'           : [ #Path('mf/base2_mf.lp'),
                                     #Path('mf/base2_mf_partition.lp'),
                                     Path('mf/base2_mf_simple.lp'),
                                     Path('mf/constraints_blai_mf.lp'),
                                     Path('mf/constraints_javier_mf.lp'),
                                     Path('mf/invariants4a_mf.lp'),
                                   ],
               'inverse_actions' : [ Path('mf/inverse_actions_mf.lp') ],
               'verify'          : [ Path('mf/base2_mf.lp'),
                                     Path('mf/invariants4a_mf.lp'),
                                   ],
               'optimize'        : [ Path('mf/optimize_mf.lp') ],
               'heuristics'      : [ Path('mf/heuristics_mf.lp') ],
               'partial'         : [ Path('partial.lp') ],
             },
    'orig' : { 'solve'           : [ Path('orig/base2.lp'),
                                     Path('orig/constraints_blai.lp'),
                                     Path('orig/constraints_javier.lp'),
                                     Path('orig/invariants4a.lp'),
                                   ],
               'verify'          : [ Path('orig/base2.lp'),
                                     Path('orig/invariants4a.lp'),
                                   ],
               'optimize'        : [ Path('orig/optimize.lp') ],
               'heuristics'      : [ Path('orig/heuristics.lp') ],
             },
}

# templates
g_templates = {
    'solve'   : '{solver} --fast-exit -Wno-global-variable {lps} {_flags_} -c opt_synthesis=1 -c opt_val={opt_val} --sat-prepro={sat_prepro} --stats=2 --time-limit={time_limit}',
    'verify'  : '{solver} --fast-exit -Wno-global-variable {lps} {samples_with_path} {_flags_} -c opt_synthesis=0 -c opt_val={opt_val} --sat-prepro={sat_prepro} --stats=2 --time-limit={time_limit}',
    'partial' : '{solver} --fast-exit -Wno-global-variable {lps} {partial_sample_path}/{sample} {_flags_} -c opt_synthesis={synthesis} -c opt_val={opt_val} --sat-prepro={sat_prepro} --stats=2 --time-limit={time_limit}',
    'flags'   : { 'fixed'   : '{options} {additional_flags}',
                  'all'     : '-c num_objects={nobj} -c max_true_atoms_per_state={max_true_atoms} {options} {additional_flags}',
                },
    '_flags_' : None,
}

# logger
def get_logger(name: str, log_file: Path, level = logging.INFO):
    logger = logging.getLogger(name)
    logger.propagate = False
    logger.setLevel(level)

    # add stdout handler
    #formatter = logging.Formatter('[%(levelname)s] %(message)s')
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] [%(funcName)s:%(lineno)d] %(message)s')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    logger.addHandler(console)

    # add file handler
    if log_file != '':
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] [%(funcName)s:%(lineno)d] %(message)s')
        file_handler = logging.FileHandler(str(log_file), 'a')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger

def close_logger(logger):
    handlers = logger.handlers
    for handler in handlers:
       logger.removeHandler(handler)
       handler.close()

# options and arguments
def get_args():
    # argument parser
    default_debug_level = 0
    parser = argparse.ArgumentParser('incremental_solver.py')
    parser.add_argument('--debug_level', type=int, default=default_debug_level, help=f'Set debug level (default={default_debug_level})')
    parser.add_argument('benchmarks', type=Path, help='Filename of file containing benchmarks')
    parser.add_argument('record', type=int, help='Record index into benchmarks file')

    # preprocessing of graphs
    graphs = parser.add_argument_group('preprocessing of input graphs')
    graphs.add_argument('--inverse_actions', action='store_true', help='Identify inverse actions in input graph')
    graphs.add_argument('--label_partitioning', action='store_true', help='Identify L1/L2 label partitioning in input graph')

    # options for solver
    default_opt_val = 1
    default_incremental = None
    default_version = 'orig'
    solver = parser.add_argument_group('options for solver')
    solver.add_argument('--heuristics', action='store_true', help='Apply heuristics when solving')
    solver.add_argument('--no_optimize', action='store_true', help="Don't do optimization when solving")
    solver.add_argument('--opt_val', type=int, default=default_opt_val, choices=[1, 2, 3], help=f"Set method for choosing state valuation, only for 'mf' version (default={default_opt_val})")
    solver.add_argument('--incremental', nargs=2, metavar=('num-samples', 'max-depth'), default=default_incremental, type=int, help=f'Set options for incremental learning (default={default_incremental})')
    solver.add_argument('--version', type=str, default=default_version, choices=['kr21', 'orig', 'mf'], help=f'Set solver version (default={default_version})')

    # options for clingo
    default_threads = 6
    default_sat_prepro = 0
    clingo = parser.add_argument_group('options for Clingo')
    clingo.add_argument('--add_lp', type=Path, action='append', default=[], help=f'Add additional .lp file for solver (possibly multiple times)')
    clingo.add_argument('--add_flag', type=str, action='append', default=[], help=f'Add additional flag for solver (possibly multiple times)')
    clingo.add_argument('--threads', type=int, default=default_threads, help=f'Set number of threads for solver (default={default_threads})')
    clingo.add_argument('--sat_prepro', type=int, default=default_sat_prepro, help=f'Set SAT preprocessing option (default={default_sat_prepro})')

    # verification
    verification = parser.add_argument_group('verification')
    verification.add_argument('--onlyver', action='store_true', help='Only do verification')
    verification.add_argument('--skipver', action='store_true', help='Skip verification step')

    # paths
    default_results_path = ''
    default_sample_path = '../samples/full'
    default_partial_sample_path = '../samples/partial'
    default_solver_path = '../clingo'
    paths = parser.add_argument_group('paths')
    paths.add_argument('--results', type=Path, default=default_results_path, help=f"Path to results folders (default='{default_results_path}')")
    paths.add_argument('--sample_path', type=Path, default=default_sample_path, help=f"Path to samples (default='{default_sample_path}')")
    paths.add_argument('--partial_sample_path', type=Path, default=default_partial_sample_path, help=f"Path to partial samples (default='{default_partial_sample_path}')")
    paths.add_argument('--solver_path', type=Path, default=default_solver_path, help=f"Path to solver files (default='{default_solver_path}')")
    paths.add_argument('--remove_dir', help='Discard existing files in results folder (if exists)', action='store_true')

    # time and mem bounds
    default_mem_bound = None
    default_time_bound = 0
    default_time_bound_ver = 0
    bounds = parser.add_argument_group('bounds')
    bounds.add_argument('--mem_bound', type=int, default=default_mem_bound, help=f'Set memory bound for solver in MBs (default={default_mem_bound})')
    bounds.add_argument('--time_bound', type=int, default=default_time_bound, help=f'Set time bound for synthesis (0 means no bound, default={default_time_bound})')
    bounds.add_argument('--time_bound_ver', type=int, default=default_time_bound_ver, help=f'Set time bound for verification (0 means no bound, default={default_time_bound_ver})')

    # parse arguments
    args = parser.parse_args()
    return args

def copy_files(filenames: List[Path], target_dir: Path, logger, prefix=None):
    for fname in filenames:
        fname_with_prefix = prefix / fname if prefix else fname
        if fname_with_prefix.exists():
            if logger: logger.info(colored(f"Copy '{fname_with_prefix.name}' to '{target_dir}'", 'green'))
            shutil.copy(fname_with_prefix, target_dir)
        else:
            if logger:
                logger.error(colored(f"File '{fname_with_prefix.name}' not found", 'red', attrs=['bold']))
            else:
                print(f"Error: file '{fname_with_prefix.name}' not found", 'red', attrs=['bold'])
            exit(0)

# create and inject instance indices for samples
def create_instances_in_destination_folder(sample_path: Path, samples: List[str], samples_catalog: dict, target_dir: Path, version: str, logger):
    assert version in [ 'orig', 'mf' ]
    actions = set()
    target_dir.mkdir(parents=True, exist_ok=True)

    # check existence of samples
    all_samples_exist = True
    for sample in samples:
        if not (sample_path / Path(sample)).exists():
            all_samples_exist = False
            print(colored(f"Error: sample file '{sample_path / Path(sample)}' doesn't exist", 'red'))
    if not all_samples_exist: exit(0)

    if version == 'orig':
        samples_with_path = [ sample_path / Path(sample) for sample in samples ]
        copy_files(samples_with_path, target_dir, logger)
        # get actions in samples
        for sample in samples:
            with (sample_path / Path(sample)).open('r') as fd:
                for line in fd.readlines():
                    if line.startswith('labelname('):
                        i, label = re.search('labelname\((\d*),"(.*)"\).', line).groups()
                        actions.add(label)
    else:
        for index, sample in enumerate(samples):
            num_nodes, num_edges = 0, 0
            if samples_catalog is not None:
                assert index + 1 not in samples_catalog
                if sample in samples_catalog:
                    logger.info('Skipping {sample} as already in catalaog as instance {samples_catalog[sample]}')
                    continue
                samples_catalog[index + 1] = sample
                samples_catalog[sample] = index + 1
            with (target_dir / Path(sample)).open('w') as wfd:
                wfd.write(f'instance({index+1}).\n')
                wfd.write(f'filename({index+1},"{sample}").\n')
                with (sample_path / Path(sample)).open('r') as rfd:
                    for line in rfd.readlines():
                        if line.startswith('node('):
                            i = line.index('(')
                            wfd.write(f'node({index+1},{line[i+1:]}')
                            num_nodes += 1
                        elif line.startswith('labelname('):
                            i, label = re.search('labelname\((\d*),"(.*)"\).', line).groups()
                            wfd.write(f'labelname({index+1},{i},"{label}").\n')
                            actions.add(label)
                        elif line.startswith('edge('):
                            i = line.index('(')
                            wfd.write(f'edge({index+1},{line[i+1:]}')
                        elif line.startswith('tlabel('):
                            i = line.index('(')
                            wfd.write(f'tlabel({index+1},{line[i+1:]}')
                            num_edges += 1
                        else:
                            wfd.write(line)
            logger.info(colored(f'Created instance {index+1} for file "{sample}" with {num_nodes} node(s) and {num_edges} edge(s)', 'green'))

    with (target_dir / Path('metadata.lp')).open('w') as fd:
        fd.write(f'num_actions({len(actions)}).\n')

    logger.info(colored(f'Created metadata file for {len(samples)} instance(s)', 'green'))
    logger.info(colored(f'{len(samples)} instance(s) created, {len(actions)} action(s) {actions}', 'green'))

# resource usage
def get_process_time_in_seconds():
    info = resource.getrusage(resource.RUSAGE_SELF)
    cinfo = resource.getrusage(resource.RUSAGE_CHILDREN)
    return info.ru_utime + info.ru_stime + cinfo.ru_utime + cinfo.ru_stime

def get_subprocess_memory():
    cinfo = resource.getrusage(resource.RUSAGE_CHILDREN)
    return cinfo.ru_maxrss

def limit_process_memory(bytes):
    resource.setrlimit(resource.RLIMIT_AS, (bytes, bytes))

# subprocesses and signals
def execute_cmd(cmd, logger=None, mem_limit=None):
    global g_logger, g_children, g_default_sigxcpu_handler
    start_time = get_process_time_in_seconds()
    stdout, stderr = '', ''
    if mem_limit != None: mem_limit *= 1e6

    g_logger = logger
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=lambda: init_process(mem_limit), encoding='utf8', bufsize=1, universal_newlines=True) as p:
        g_running_children.append(p)
        for line in p.stdout:
            stdout += line
            if logger: logger.info(line.strip('\n'))
        for line in p.stderr:
            stderr += line
            if logger: logger.info(line.strip('\n'))
    g_running_children.pop()
    g_logger = None

    #p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=lambda: init_process(mem_limit), encoding='utf8')
    #stdout = ''
    #for line in p.stdout:
    #    stdout += line
    #stderr = ''
    #for line in p.stderr:
    #    stderr += line
    #p.wait()

    elapsed_time = get_process_time_in_seconds() - start_time
    max_memory = get_subprocess_memory()
    return stdout, stderr, elapsed_time, max_memory

g_logger = None
g_running_children = []
g_default_sigxcpu_handler = signal.signal(signal.SIGXCPU, signal.SIG_IGN)

def sigterm_handler(_signo, _stack_frame):
    if g_logger:
        g_logger.warning(colored('Process INTERRUPTED by SIGTERM!', 'red'))
    else:
        print(colored('Process INTERRUPTED by SIGTERM!', 'red'))
    if g_running_children:
        if g_logger:
            g_logger.info(f'Killing {len(g_running_children)} subprocess(es)...')
        else:
            print(f'Killing {len(g_running_children)} subprocess(es)...')
        for p in g_running_children:
            p.kill()
    exit(0)

def restore_sigxcpu_handler():
    signal.signal(signal.SIGXCPU, g_default_sigxcpu_handler)

def init_process(mem_limit):
    restore_sigxcpu_handler()
    if mem_limit:
        resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))

# select records from benchmarks file
def get_records(fname: Path, record):
    benchmarks = []
    with fname.open('r') as fd:
        benchmark_index = 0
        for line in fd:
            fields = line.rstrip('\n').split()
            if len(fields) == 0 or fields[0][0] == '#': continue
            if benchmark_index == record:
                benchmarks.append(Benchmark(fields))
            benchmark_index += 1
    print(colored(f"Got {len(benchmarks)} record(s) from '{fname}' using record '{record}'", 'green'))
    return benchmarks

# parse output from clingo and store stats
def parse_clingo_output(stdout, stderr, elapsed_time, max_memory, version: str):
    result = parse_clingo_out_orig(stdout) if version == 'orig' else parse_clingo_out_mf(stdout)
    stats = dict(total_time=result.get(TIME, elapsed_time),
                 solve_memory=max_memory,
                 satisfiable=result[SAT], # True, False, or None for Sat, Unsat, Unknown
                 num_rules=result.get(RULES, -1),
                 num_variables=result.get(VARIABLES, -1),
                 num_constraints=result.get(CONSTRAINTS, -1),
                 solve_time=result.get(SOLVING, -1),
                 first=result.get(MODEL1st, -1),
                 num_choices=result.get(CHOICES, -1),
                 num_conflicts=result.get(CONFLICTS, -1))
    return result, stats

# parse clingo output
def create_schema_from_symbols(result, version: str):
    symbols = result[SYMBOLS]
    schema = STRIPSSchema_orig.create_from_clingo(symbols) if version == 'orig' else STRIPSSchema_mf.create_from_clingo(symbols)
    decoded = schema.get_string(val=False)
    model = schema.get_schema()
    return schema, decoded, model

# write output
def write_output(filename: Path, output: str, logger):
    logger.info(colored(f'Writing output to {filename}', 'green'))
    with filename.open('w') as fd:
        fd.write(f'{output}\n')

def solve_and_parse_output(task, parameters, stats, logger, extra_lps: List[Path] = None):
    version = parameters['version']
    local_parameters = dict(parameters)
    local_parameters.update(synthesis=1)
    dirpath = local_parameters['dirpath']

    # copy lps to dirpath
    assert version in g_clingo
    copy_files(g_clingo[version]['solve'], dirpath, logger, prefix=parameters['solver_path'])
    if args.inverse_actions:
        copy_files(g_clingo[version]['inverse_actions'], dirpath, logger, prefix=parameters['solver_path'])
    if not args.no_optimize:
        copy_files(g_clingo[version]['optimize'], dirpath, logger, prefix=parameters['solver_path'])
    if args.heuristics:
        copy_files(g_clingo[version]['heuristics'], dirpath, logger, prefix=parameters['solver_path'])
    if task.partial != None:
        copy_files(g_clingo[version]['partial'], dirpath, logger, prefix=parameters['solver_path'])
    if extra_lps is not None:
        copy_files(extra_lps, dirpath, logger)
    copy_files(local_parameters['add_lp'], dirpath, logger)
    local_parameters.update(lps=f"'{str(dirpath)}'/*.lp")
    local_parameters.update(_flags_=g_templates['flags']['all'].format(**local_parameters))

    logger.info(colored("Solve '{samples}' with flags '{_flags_}'".format(**local_parameters), 'magenta', attrs=['bold']))
    template = g_templates['solve'] if task.partial == None else g_templates['partial']
    logger.info('cmdline=|{}|'.format(template.format(**local_parameters)))

    stdout, stderr, elapsed_time, max_memory = execute_cmd(template.format(**local_parameters), logger=logger, mem_limit=args.mem_bound)
    result, solve_stats = parse_clingo_output(stdout, stderr, elapsed_time, max_memory, version)
    solve_stats.update(objs=local_parameters['nobj'])
    logger.info(colored(f"Elapsed time {solve_stats['total_time']} second(s)", 'green'))
    logger.info(colored(f"Memory {solve_stats['solve_memory']}", 'green'))
    stats.data.update(**solve_stats)

    # save output
    write_output(dirpath / 'solver_stdout.txt', stdout, logger)
    write_output(dirpath / 'solver_stderr.txt', stderr, logger)

    if not solve_stats['satisfiable']:
        if solve_stats['satisfiable'] == False:
            comment = 'UNSATISFIABLE ... skipping verification'.format(**local_parameters)
        else:
            comment = 'INDETERMINATE ... skipping verification'.format(**local_parameters)
        logger.info(colored(comment, 'red', attrs=['bold']))
        logger.info(colored(f'Stats: {stats}', 'green'))
        logger.info(colored(f'Stats.data: {stats.data}', 'green'))
        close_logger(logger)
        with local_parameters['stats'].open('w') as fd: fd.write(str(stats) + '\n')

    return result

def incremental_solve_and_parse_output(task, parameters, stats, logger):
    assert not task.partial
    assert parameters['incremental'] is not None

    local_parameters = dict(parameters)
    local_parameters.update(synthesis=1)
    dirpath = local_parameters['dirpath']
    num_samples, depth = parameters['incremental']

    # CHECK: following line hasn't been updated!
    logger.info(colored("Incremental solve '{samples}' with flags '{_flags_}'".format(**local_parameters), 'magenta', attrs=['bold']))
    template = g_templates['solve']

    # get DFAs
    dfa_fnames = [ local_parameters['sample_path'] / Path(sample) for sample in task.samples ]
    dfas = [ DFA(fname) for fname in dfa_fnames ]
    num_dfas = len(dfas)

    solved_tasks = [ 0 ] * num_dfas
    marked_nodes = [ set() for i in range(num_dfas) ]
    while sum(solved_tasks) < num_dfas:
        logger.info(colored(f'solved_tasks: {solved_tasks}', 'green', attrs=['bold']))

        # sample trajectories
        for i, dfa in enumerate(dfas):
            if solved_tasks[i] == 0 and len(marked_nodes[i]) < dfa.num_nodes:
                # CHECK: exploration values: number of sources, length of trajectories
                sampled_sources = dfa.sample_nodes(num_samples, avoid=marked_nodes[i])
                #sampled_sources = [0,1,3,5] #CHECK
                sampled_paths = [ dfa.sample_path(src, depth, repeat=False) for src in sampled_sources ]
                logger.info(colored(f'dfa={i}: sampled_paths={sampled_paths}', 'magenta', attrs=['bold']))
                for path in sampled_paths:
                    for node in path:
                        marked_nodes[i].add(node)
                logger.info(colored(f'dfa={i}: marked_nodes: #={len(marked_nodes[i])}, nodes={sorted(list(marked_nodes[i]))}', 'magenta', attrs=['bold']))

        # create marked.lp file with marked nodes
        marked_fname = dirpath / Path('marked.lp')
        with marked_fname.open('w') as fd:
            for i in range(num_dfas):
                if len(marked_nodes[i]) > 0:
                    instance_index = local_parameters['catalog'][dfa_fnames[i].name]
                    fd.write(f'partial({instance_index},"{dfa_fnames[i].name}").\n')
                    for node in marked_nodes[i]:
                        fd.write(f'marked({instance_index},{node}).\n')

        # remove any left over last call
        local_parameters['model'].unlink(missing_ok=True)

        # call solver and parse output
        result = solve_and_parse_output(task, local_parameters, stats, logger)
        if not result[SAT]: return result

        # parse output
        schema, decoded, model = create_schema_from_symbols(result, parameters['version'])

        # write model to file
        logger.info(colored("Save '{model}'".format(**local_parameters), 'green'))
        with local_parameters['model'].open('w') as fd: fd.write('\n'.join(str(s) + '.' for s in model))

        # verify model over complete instance
        solved_tasks = []

        for i in range(num_dfas):
            result, verify_stats = verify_instance(task.samples[i], local_parameters['nobj'], local_parameters['max_true_atoms'], task, parameters, logger=logger)
            solved_tasks.append(1 if verify_stats['satisfiable'] else 0)

        # check for failure
        failure = True
        for i in range(num_dfas):
            if len(marked_nodes[i]) < dfas[i].num_nodes:
                failure = False
                break
        if failure: break

    logger.info(colored(f'solved_tasks: {solved_tasks}', 'green', attrs=['bold']))
    logger.info(colored(f'marked_nodes: {[ len(mn) for mn in marked_nodes ]}', 'green', attrs=['bold']))

    return result

def verify_instance(instance, nobj, max_true_atoms, task, parameters, logger):
    # create and populate verify folder
    version = parameters['version']
    dirpath = parameters['dirpath']
    verify_path = dirpath / Path('verify')
    create_instances_in_destination_folder(parameters['sample_path'], [instance], None, verify_path, version, logger)

    # parameters for solver
    assert version in g_clingo
    local_parameters = dict(parameters)
    verify_solver_with_path = [parameters['solver_path'] / fname for fname in g_clingo[version]['verify']]
    lps = [str(fname) for fname in verify_solver_with_path] + ["'{model}'".format(**local_parameters)]
    local_parameters.update(lps=' '.join(lps))
    local_parameters.update(sample=Path(instance))
    local_parameters.update(samples_with_path=f"'{dirpath}'/verify/{instance}")
    local_parameters.update(nobj=nobj)
    local_parameters.update(max_true_atoms=max_true_atoms)
    local_parameters.update(options=task.get_options())
    local_parameters.update(_flags_=g_templates['flags']['all'].format(**local_parameters))
    local_parameters.update(time_limit=local_parameters['time_bound_ver'])
    local_parameters.update(synthesis=0)

    # call solver and parse output
    logger.info(colored("Verify '{sample}' with flags '{_flags_}'".format(**local_parameters), 'cyan', attrs=['bold']))
    logger.info('cmdline=|{}|'.format(g_templates['verify'].format(**local_parameters)))

    stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['verify'].format(**local_parameters), logger=logger, mem_limit=args.mem_bound)
    result, verify_stats = parse_clingo_output(stdout, stderr, elapsed_time, max_memory, version)
    verify_stats.update(instance=instance, objs=nobj)

    # save output and decoded model (if any)
    write_output((verify_path / f'solver_stdout_{Path(instance).stem}_{nobj}').with_suffix('.txt'), stdout, logger)
    write_output((verify_path / f'solver_stderr_{Path(instance).stem}_{nobj}').with_suffix('.txt'), stderr, logger)
    if verify_stats['satisfiable']:
        schema, decoded, model = create_schema_from_symbols(result, version)
        with (verify_path / f'decoded_{instance}').with_suffix('.txt').open('w') as fd:
            fd.write(decoded)
    elif verify_stats['satisfiable'] == None:
        logger.info(colored(f"Indeterminate output from solver during verification of '{instance}'", 'magenta'))

    return result, verify_stats

def verify_and_parse_output(task, parameters, stats, logger):
    logger.info(colored(f'Verification for instances {task.verify}', 'cyan', attrs=['bold']))
    if stats: stats.data['verify'] = []
    complete_verification = True
    for instance in task.verify:
        # iterate over num objects, until max, looking for a verification
        successful = False
        for nobj in range(1, 1 + task.max_num_objs_ver):
            result, verify_stats = verify_instance(instance, nobj, task.max_true_atoms_ver, task, parameters, logger)
            logger.info(colored(f"Elapsed time {verify_stats['total_time']} second(s)", 'green'))
            logger.info(colored(f"Memory {verify_stats['solve_memory']}", 'green'))

            if stats: stats.data['verify'].append(verify_stats)
            if verify_stats['satisfiable']:
                successful = True
                break

        if not successful:
            comment = f"Unsuccessful verification of '{instance}' ... skipping remaining verification"
            logger.info(colored(comment, 'red', attrs=['bold']))
            if stats:
                logger.info(colored(f'Stats: {stats}', 'green'))
                logger.info(colored(f'Stats.data: {stats.data}', 'green'))
            close_logger(logger)
            if stats:
                with parameters['stats'].open('w') as fd:
                    fd.write(str(stats) + '\n')
            complete_verification = False
            break
    return complete_verification

def main(args: dict):
    if not args.benchmarks.exists():
        print(colored(f"Error: benchmarks file '{args.benchmarks}' doesn't exist", 'red', attrs=['bold']))
        exit(0)

    benchmarks = get_records(args.benchmarks, args.record)
    for task in benchmarks:
        # create stats object
        stats = Stats(task.samples)

        # add additional options
        task.options.append(f'-t {args.threads}')
        if args.heuristics:
            task.options.append(f'--heuristic=Domain')

        # base parameters
        parameters = {
            'version'             : args.version,
            'solver'              : 'clingo',
            'sample_path'         : args.sample_path,
            'partial_sample_path' : args.partial_sample_path,
            'solver_path'         : args.solver_path,
            'nobj'                : task.num_objs_learn,
            'max_true_atoms'      : task.max_true_atoms_learn,
            'time_bound'          : args.time_bound,
            'time_bound_ver'      : args.time_bound_ver,
            'additional_flags'    : ' '.join(args.add_flag),
            'options'             : task.get_options(),
            'sat_prepro'          : args.sat_prepro,
            'opt_val'             : args.opt_val,
        }
        parameters.update(flags_fixed=g_templates['flags']['fixed'].format(**parameters).strip(' '))
        parameters.update(time_limit=parameters['time_bound'])

        # construct dirpath
        samples = ' '.join(task.samples)
        dirname = (samples + ' ' + g_templates['flags']['all'].format(**parameters).strip(' ')).replace(' ', '_')
        parameters.update(dirname=dirname)
        max_filename_length = os.pathconf('.', 'PC_NAME_MAX')
        if len(dirname) > max_filename_length:
            print(colored(f"Warning: truncating dirname to OS's max filename length of {max_filename_length}, current length is {len(dirname)}!", 'red', attrs=['bold']))
            dirname = dirname[:max_filename_length]
        dirpath = args.results / Path(dirname)
        parameters.update(dirpath=dirpath)

        # create dir (if existing, skip it)
        if not args.onlyver:
            if dirpath.exists():
                if not dirpath.is_dir():
                    print(colored(f"Skipping benchmark '{dirpath}' because non-folder file with same name exists ...", 'magenta'))
                    continue
                elif args.remove_dir:
                    print(colored(f"Folder '{dirpath}' removed", 'magenta'))
                    shutil.rmtree(dirpath)
                    print(colored(f"Folder '{dirpath}' created", 'green'))
                    dirpath.mkdir(parents=True)
                else:
                    print(colored(f"Skipping benchmark '{dirpath}' because folder with same name exists (use --remove_dir to force) ...", 'magenta'))
                    continue
            else:
                print(colored(f"Folder '{dirpath}' created", 'green'))
                dirpath.mkdir(parents=True)

            assert dirpath.exists() and dirpath.is_dir()

        # set filenames
        log_fname = 'log.onlyver.0.txt' if args.onlyver else 'log.0.txt'
        stats_fname = 'stats.onlyver.0.txt' if args.onlyver else 'stats.0.txt'
        decoded_fname = 'decoded.txt'
        model_fname = 'model.lp'
        parameters.update(log=dirpath / log_fname)
        parameters.update(decoded=dirpath / decoded_fname)
        parameters.update(model=dirpath / model_fname)
        parameters.update(stats=dirpath / stats_fname)
        assert args.onlyver or (not parameters['log'].exists() and not parameters['stats'].exists())

        # setup logger and identify call
        log_level = logging.INFO if args.debug_level == 0 else logging.DEBUG
        logger = get_logger('solve', parameters['log'], log_level)
        logger.info(colored(f"Log file '{parameters['log']}' created", 'green'))
        logger.info(f'call=|{" ".join(sys.argv)}|')

        # update samples in parameters and preprocess them in order to add instance indices
        parameters.update(samples=task.samples)
        parameters.update(metadata=f"'{parameters['dirpath']}'/metadata.lp")
        parameters.update(samples_with_path=' '.join([ f"'{parameters['dirpath']}'/{sample}" for sample in task.samples ]))
        parameters.update(catalog=dict())
        create_instances_in_destination_folder(args.sample_path, task.samples, parameters['catalog'], parameters['dirpath'], parameters['version'], logger)

        # set additional files in parameters (those passed with --add)
        parameters.update(add_lp=args.add_lp)

        # solve task?
        if args.onlyver:
            # mark theory as SAT so that stats can be printed without error
            stats.data.update(satisfiable=True)
        else:
            # call solver and parse output
            parameters.update(incremental=args.incremental)
            if args.incremental:
                result = incremental_solve_and_parse_output(task, parameters, stats, logger)
            else:
                result = solve_and_parse_output(task, parameters, stats, logger)
            if not result[SAT]: continue

            # parse output
            schema, decoded, model = create_schema_from_symbols(result, parameters['version'])

            # write model to file
            logger.info(colored("Save '{model}'".format(**parameters), 'green'))
            with parameters['model'].open('w') as fd:
                fd.write('\n'.join(str(s) + '.' for s in model))
                fd.write('\n')

            # decode
            logger.info(colored("Decode '{model}' to obtain description '{decoded}'".format(**parameters), 'green'))
            with parameters['decoded'].open('w') as fd: fd.write(decoded)

        # verification?
        if not args.skipver:
            complete_verification = verify_and_parse_output(task, parameters, stats, logger)

        logger.info(colored(f'Stats: {stats}', 'green'))
        logger.info(colored(f'Stats.data: {stats.data}', 'green'))
        if not args.skipver and complete_verification:
            logger.info(colored(f'Complete verification!', 'green', attrs=['bold']))
        close_logger(logger)
        with parameters['stats'].open('w') as fd: fd.write(str(stats) + '\n')

if __name__ == '__main__':
    # setup exec name
    exec_path = Path(sys.argv[0]).parent
    exec_name = Path(sys.argv[0]).stem

    # setup proper SIGTERM handler
    signal.signal(signal.SIGTERM, sigterm_handler)

    # read arguments and do job
    args = get_args()
    main(args)

