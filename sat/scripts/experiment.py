import sys, os, shutil, signal, subprocess, resource, subprocess, argparse
from termcolor import colored
from pathlib import Path
import logging

class Benchmark:
    def __init__(self, fields):
        index = 0
        self.marker = int(fields[index]); index += 1
        self.prefix = fields[index]; index += 1

        # actions
        self.num_actions = int(fields[index]); index += 1
        self.action_max_arities = []
        for i in range(0, self.num_actions):
            self.action_max_arities.append(int(fields[index]))
            index += 1

        # atoms
        self.num_atoms = int(fields[index]); index += 1
        self.atom_max_arities = []
        for i in range(0, self.num_atoms):
            self.atom_max_arities.append(int(fields[index]))
            index += 1

        # meta-features and static predicates
        self.num_meta_features = int(fields[index]); index += 1
        self.num_static_unary_predicates = int(fields[index]); index += 1
        self.num_static_binary_predicates = int(fields[index]); index += 1

        # dfas
        self.dfas = []
        while fields[index] != 'verify':
            num_features = int(fields[index]); index += 1
            num_objects = int(fields[index]); index += 1
            lower_bound = int(fields[index]); index += 1
            upper_bound = int(fields[index]); index += 1
            dfa = Path(fields[index]); index += 1
            if dfa.suffix != '.dfa': dfa = dfa.with_suffix('.dfa')
            self.dfas.append([num_features, num_objects, lower_bound, upper_bound, dfa])
        index += 1

        # max number objects for verification
        self.max_num_objs_ver = int(fields[index]); index += 1

        # test set
        self.verify = []
        while index < len(fields):
            dfa = Path(fields[index]); index += 1
            if dfa.suffix != '.dfa': dfa = dfa.with_suffix('.dfa')
            self.verify.append(dfa)

class Stats:
    def __init__(self):
        self.data = dict()
    def __str__(self):
        str_template = '{theory} {num_variables} {num_implications} {generation_time} {solve_time} {solve_memory} {satisfiable}'
        str_subtemplate = '{vtheory} {num_variables} {num_implications} {generation_time} {solve_time} {solve_memory} {satisfiable}'
        as_str = str_template.format(**self.data)
        if self.data['satisfiable'] == 'True':
            for ver in self.data['verify']:
                as_str += ' ' + str_subtemplate.format(**ver)
        return as_str

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
    parser = argparse.ArgumentParser('experiment.py')
    parser.add_argument('--debug_level', nargs=1, type=int, default=0, help='Set debug level')
    parser.add_argument('--dfas_path', nargs=1, type=Path, default=None, help='Path to dfas')
    parser.add_argument('--extra_flags', nargs=1, type=str, default=None, help='Extra flags for encoder')
    parser.add_argument('--extra_path', nargs=1, type=Path, default=None, help='Extra path prefix for dir names')
    parser.add_argument('--extra_prefix', nargs=1, type=str, default=None, help='Extra prefix for naming files')
    parser.add_argument('--recover', help='Try to recover incomplete experiments; overrides --remove_dir', action='store_true')
    parser.add_argument('--remove_dir', help='Remove folder contents if exists', action='store_true')
    parser.add_argument('--time_bound', type=int, default=0, help='Set time bound')
    parser.add_argument('--verbose', type=int, default=0, help='Set verbosity level (meaningful values in [0..3])')
    parser.add_argument('benchmarks', type=Path, help='Filename of file containing benchmarks')
    parser.add_argument('record', type=int, help='Record index into benchmarks file')

    # parse arguments
    args = parser.parse_args()

    # implied options
    if args.recover: args.remove_dir = False

    args.extra_flags = [] if args.extra_flags is None else args.extra_flags
    args.extra_path = [Path('.')] if args.extra_path is None else args.extra_path
    args.dfas_path = [Path('.')] if args.dfas_path is None else args.dfas_path
    args.extra_prefix = '' if args.extra_prefix is None else args.extra_prefix

    return args

# templates
g_templates = {
    'theory'  : '{strips} --output {theory} {flags} {meta_layer_parameters}   {dfas}',
    'vtheory' : '{strips} --verify-meta-layer {meta} --output {vtheory}   {flags}   {meta_layer_parameters}   {vdfa}',
    'solve'   : '{ulimit} ; {solver} {gztheory} {model}',
    'verify'  : '{ulimit} ; {solver} {gzvtheory} {vmodel}',
    'decode'  : '{strips} --decode-model --input {model} --output {decoded}   {flags}   {meta_layer_parameters}   {dfas}',
    'decode_meta_layer' : '{strips} --decode-meta-layer --input {model} --output {meta}   {flags}   {meta_layer_parameters}   {dfas}',
    'gzip'    : 'gzip {gzinput}',
}

# subprocesses and signals
def get_process_time_in_seconds():
    info = resource.getrusage(resource.RUSAGE_SELF)
    cinfo = resource.getrusage(resource.RUSAGE_CHILDREN)
    return info.ru_utime + info.ru_stime + cinfo.ru_utime + cinfo.ru_stime

def get_subprocess_memory():
    cinfo = resource.getrusage(resource.RUSAGE_CHILDREN)
    return cinfo.ru_maxrss

def execute_cmd(cmd, logger=None, mem_limit=None, verbose=True):
    global g_logger, g_children, g_default_sigxcpu_handler
    start_time = get_process_time_in_seconds()
    stdout, stderr = [], []
    if mem_limit != None: mem_limit *= 1e6

    g_logger = logger
    with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=lambda: init_process(mem_limit), encoding='utf8', bufsize=1, universal_newlines=True) as p:
        g_running_children.append(p)
        for line in p.stdout:
            line = line.strip('\n')
            stdout.append(line)
            if verbose and logger: logger.info(line)
        for line in p.stderr:
            line = line.strip('\n')
            stderr.append(line)
    g_running_children.pop()
    g_logger = None

    if stderr != []:
        lines = '\n'.join(stderr)
        if logger: logger.error(colored(f"{lines}", 'red'))

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
    return benchmarks

# determine satisfiability status from output
def determine_satisfiability(output_lines):
    for i in range(1, 1 + len(output_lines)):
        if output_lines[-i] == "s INDETERMINATE":
            return "Indet"
        elif output_lines[-i] == "s UNSATISFIABLE":
            return "False"
        elif output_lines[-i] == "s SATISFIABLE":
            return "True"
    return "Indet"

# incomplete marker
def existing_marker(dirpath: Path):
    marker = dirpath / 'incomplete.marker'
    return marker.exists()

def create_marker(dirpath: Path):
    marker = dirpath / 'incomplete.marker'
    print(colored(f"Marker '{marker}' created", 'green'))
    marker.touch()

def remove_marker(dirpath: Path, logger=None):
    assert existing_marker(dirpath)
    marker = dirpath / 'incomplete.marker'
    if logger: logger.info(colored(f"Marker '{marker}' removed", 'magenta'))
    marker.unlink()

def main(args: dict):
    benchmarks = get_records(args.benchmarks, args.record)
    for task in benchmarks:
        # create stats object
        stats = Stats()

        # basic parameters
        parameters = {
            'strips' : '../src/strips',
            'solver' : 'glucose',
            'ulimit' : 'ulimit -c 0',
        }
        if args.time_bound > 0: parameters['ulimit'] = f'ulimit -c 0 -t {args.time_bound}'

        # name
        name = f'{args.extra_prefix}_{task.prefix}' if args.extra_prefix != '' else task.prefix
        parameters.update({ 'name' : name })

        # flags
        flags = [ '--regularizer exact-arities', '--disable-colors' ]
        flags.extend(args.extra_flags)
        parameters.update({ 'flags' : ' '.join(flags) })

        # dirpath
        dirpath = args.extra_path[0] / Path(parameters['name'])
        parameters.update({ 'dirpath' : dirpath })

        # create dir (if exist, skip it)
        recovery_mode = False
        if dirpath.exists():
            if not dirpath.is_dir():
                print(colored(f"Skipping benchmark '{dirpath}' because non-folder file with same name exists ...", 'magenta'))
                continue
            elif args.remove_dir:
                print(colored(f"Folder '{dirpath}' removed", 'magenta'))
                shutil.rmtree(dirpath)
                print(colored(f"Folder '{dirpath}' created", 'green'))
                dirpath.mkdir()
                create_marker(dirpath)
            elif existing_marker(dirpath) and args.recover:
                print(colored(f"Found marker in '{dirpath}'; entering recovery mode...", 'magenta'))
                recovery_mode = True
            else:
                print(colored(f"Skipping benchmark '{dirpath}' because folder with same name exists (use either --recover or --remove_dir to force)  ...", 'magenta'))
                continue
        else:
            print(colored(f"Folder '{dirpath}' created", 'green'))
            dirpath.mkdir()
            create_marker(dirpath)

        assert dirpath.exists() and dirpath.is_dir()
        assert existing_marker(dirpath)

        # set log and stats name
        log_index = 0
        if recovery_mode:
            while (dirpath / f'log.{log_index}.txt').exists():
                log_index += 1
        log_name = dirpath / f'log.{log_index}.txt'
        stats_name = dirpath / f'stats.{log_index}.txt'
        assert not log_name.exists() and not stats_name.exists()
        parameters.update({ 'log' : log_name, 'stats' : stats_name })

        # setup logger
        log_level = logging.INFO if args.debug_level == 0 else logging.DEBUG
        logger = get_logger('solve', parameters['log'], log_level)
        logger.info(colored(f"Log file '{parameters['log']}' created", 'green'))
        if recovery_mode: logger.info('Recovery mode is TRUE')

        # parameters for dfas
        dfas = []
        for dfa in task.dfas:
            p = [ str(n) for n in dfa[0:4] ] + [ str(args.dfas_path[0] / dfa[4]) ]
            dfas.append(' '.join(p))
        parameters.update({ 'dfas' : '   '.join(dfas) })

        # filenames for theory and compressed theory
        theory = dirpath / 'theory.cnf'
        gztheory = theory.with_suffix(theory.suffix + '.gz')
        parameters.update({ 'theory' : theory, 'gztheory' : gztheory })

        # parameters for meta-layer
        p = [ task.num_actions ] + task.action_max_arities + [ task.num_atoms ] + task.atom_max_arities + [ task.num_meta_features, task.num_static_unary_predicates, task.num_static_binary_predicates ]
        parameters.update({ 'meta_layer_parameters' : ' '.join([ str(n) for n in p ]) })

        # construct and compress theory
        if recovery_mode and gztheory.exists():
            logger.info(colored(f"Compressed theory '{gztheory}' found; assuming it is OK and continuing with next step...", 'magenta'))
            stats.data.update({ 'num_variables' : -1, 'num_implications' : -1 })
            stats.data.update({ 'theory' : parameters['theory'], 'generation_time' : -1 })
        else:
            # remove any left over
            assert not gztheory.exists()
            theory.unlink(missing_ok=True)

            # construct theory
            logger.info(colored(f"Generate theory '{theory}'", 'green'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['theory'].format(**parameters), logger, verbose=args.verbose > 0)
            if args.verbose > 0:
                logger.info(f'Elapsed time {elapsed_time} second(s)')
                logger.info(f'Memory {max_memory}')

            # if something wrong (i.e. no theory), continue
            if not theory.exists():
                logger.info(colored(f"No theory '{theory}' in disk ... skipping rest", 'red'))
                close_logger(logger)
                continue

            # find number of variables and implications
            stats.data.update({ 'num_variables' : -1, 'num_implications' : -1 })
            for line in stdout:
                if line[:11] == '#variables=':
                    num_variables = int(line.split('=')[1])
                    stats.data.update({ 'num_variables' : num_variables })
                elif line[:14] == '#implications=':
                    num_implications = int(line.split('=')[1])
                    stats.data.update({ 'num_implications' : num_implications })

            # compress theory
            logger.info(colored(f"Compress theory '{theory}'", 'green'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['gzip'].format(gzinput = parameters['theory']), logger)
            stats.data.update({ 'theory' : parameters['theory'], 'generation_time' : elapsed_time })

        # filenames for solving theory
        model = dirpath / 'model.cnf'
        decoded = dirpath / 'decoded.txt'
        meta = dirpath / 'meta.txt'
        parameters.update({ 'model' : model, 'decoded' : decoded, 'meta' : meta })

        # solve theory, decode it, and compress model
        if recovery_mode and meta.exists():
            logger.info(colored(f"Meta-layer '{meta}' found; assuming it is OK and continuing with next step...", 'magenta'))
            stats.data.update({ 'solve_time' : -1, 'solve_memory' : -1, 'satisfiable' : 'True' })
        else:
            # remove any left over
            assert not model.with_suffix(model.suffix + '.gz').exists()
            assert not decoded.exists()
            assert not meta.exists()
            model.unlink(missing_ok=True)

            # solve compressed theory
            logger.info(colored(f"Solve compressed theory '{gztheory}' with '{parameters['solver']}' ...", 'magenta'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['solve'].format(**parameters), logger, verbose=args.verbose > 0)
            #logger.info(f'Elapsed time {elapsed_time} second(s)')
            #logger.info(f'Memory {max_memory}')
            stats.data.update({ 'solve_time' : elapsed_time })
            stats.data.update({ 'solve_memory' : max_memory })

            # determine SAT status
            satisfiable = determine_satisfiability(stdout)
            stats.data.update({ 'satisfiable' : satisfiable })

            # if UNSAT, jump to next benchmark
            if satisfiable == 'True':
                logger.info(colored(f"Successful solve of '{gztheory}'", 'cyan', attrs=['bold']))
            else:
                logger.info(f"UNSATISFIABLE/INDETERMINATE '{theory}' ... skipping verification")
                logger.info(f'stats: {str(stats)}')
                with stats_name.open('w') as fd: fd.write(str(stats) + '\n')
                if satisfiable == 'False': remove_marker(dirpath, logger)
                close_logger(logger)
                continue

            # decode
            logger.info(colored(f"Decode '{model}' to obtain description '{decoded}' ...", 'green'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['decode'].format(**parameters), logger, verbose=args.verbose > 2)

            # decode meta-layer
            logger.info(colored(f"Decode '{model}' to obtain meta-layer '{meta}' ...", 'green'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['decode_meta_layer'].format(**parameters), logger, verbose=args.verbose > 2)

            # compress model
            logger.info(colored(f"Compress '{model}' ...", 'green'))
            stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['gzip'].format(gzinput = parameters['model']), logger)

        # verify
        complete_verification = True
        stats.data.update({ 'verify' : [] })
        for dfa in task.verify:
            complete_verification = False
            dfa_with_path = args.dfas_path[0] / dfa
            stats.data['verify'].append(dict())

            # iterate over num objects, until max, looking for a verification
            successful = False
            for n in range(1, 1 + task.max_num_objs_ver):
                # parameters for verification
                vtheory = dirpath / f'vtheory_ON_{n}_{dfa}.cnf'
                gzvtheory = vtheory.with_suffix(vtheory.suffix + '.gz')
                marker = dirpath / f'marker_{name}_ON_{n}_{dfa}'
                parameters.update({ 'n' : n, 'dfa' : dfa })
                parameters.update({ 'vdfa' : f'0 {n} 0 0 {dfa_with_path}' })
                parameters.update({ 'vtheory' : vtheory, 'gzvtheory' : gzvtheory, 'marker' : marker })
                # generate verification theory
                if recovery_mode and gzvtheory.exists():
                    logger.info(colored(f"Compressed verification theory '{gzvtheory}' found; assuming it is OK and continuing with next step...", 'magenta'))
                    stats.data['verify'][-1].update({ 'num_variables' : -1, 'num_implications' : -1 })
                    stats.data['verify'][-1].update({ 'vtheory' : g_parameters['vtheory'], 'generation_time' : -1 })
                elif recovery_mode and (marker.with_suffix('.YES').exists() or marker.with_suffix('.NO').exists()):
                    logger.info(colored(f"Marker for compressed '{gzvtheory}' found; assuming it is OK and continuing with next step...", 'magenta'))
                    successful = marker.with_suffix('.YES').exists()
                    stats.data['verify'][-1].update({ 'num_variables' : -1, 'num_implications' : -1 })
                    stats.data['verify'][-1].update({ 'vtheory' : g_parameters['vtheory'], 'generation_time' : -1 })
                    stats.data['verify'][-1].update({ 'solve_time' : -1, 'solve_memory' : -1 })
                    if successful:
                        stats.data['verify'][-1].update({ 'satisfiable' : 'True' })
                        break
                    else:
                        stats.data['verify'][-1].update({ 'satisfiable' : 'False' })
                else:
                    # remove any left over
                    assert not gzvtheory.exists()
                    vtheory.unlink(missing_ok=True)

                    # generate verification theory
                    logger.info(colored(f"Generate verification theory '{vtheory}'", 'green'))
                    stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['vtheory'].format(**parameters), logger, verbose=args.verbose > 1)
                    #logger.info(f'Elapsed time {elapsed_time} second(s)')
                    #logger.info(f'Memory {max_memory}')

                    # if something wrong (i.e. no theory), continue
                    if not vtheory.exists():
                        logger.info(colored(f"No verification theory '{vtheory}' in disk ... skipping rest", 'red'))
                        close_logger(logger)
                        break

                    # find number of variables and implications
                    stats.data['verify'][-1].update({ 'num_variables' : -1, 'num_implications' : -1 })
                    for line in stdout:
                        if line[:11] == '#variables=':
                            num_variables = int(line.split('=')[1])
                            stats.data['verify'][-1].update({ 'num_variables' : num_variables })
                        elif line[:14] == '#implications=':
                            num_implications = int(line.split('=')[1])
                            stats.data['verify'][-1].update({ 'num_implications' : num_implications })

                    # compress verification theory
                    logger.info(colored(f"Compress verification theory '{vtheory}'", 'green'))
                    stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['gzip'].format(gzinput = parameters['vtheory']), logger)
                    stats.data['verify'][-1].update({ 'vtheory' : parameters['vtheory'], 'generation_time' : elapsed_time })

                # parameters for verifier
                vmodel = dirpath / f'vmodel_ON_{n}_{dfa}.cnf'
                parameters.update({ 'vmodel' : vmodel })

                # verification
                if recovery_mode and (marker.with_suffix('.YES').exists() or marker.with_suffix('.NO').exists()):
                    if gzvtheory.exists():
                        gzvtheory.unlink()
                        logger.info(colored(f"Marker for compressed '{gzvtheory}' found; assuming it is OK and continuing with next step...", 'magenta'))
                    successful = marker.with_suffix('.YES').exists()
                    stats.data['verify'][-1].update({ 'solve_time' : -1, 'solve_memory' : -1 })
                    if successful:
                        stats.data['verify'][-1].update({ 'satisfiable' : 'True' })
                        break
                    else:
                        stats.data['verify'][-1].update({ 'satisfiable' : 'False' })
                else:
                    # remove any left over
                    vmodel.with_suffix(vmodel.suffix + '.gz').unlink(missing_ok=True)

                    # call verifier
                    logger.info(colored(f"Verify compressed theory '{gzvtheory}' with '{parameters['solver']}' ...", 'magenta'))
                    stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['verify'].format(**parameters), logger, verbose=args.verbose > 1)
                    #logger.info(f'Elapsed time {elapsed_time} second(s)')
                    #logger.info(f'Memory {max_memory}')
                    stats.data['verify'][-1].update({ 'solve_time' : elapsed_time })
                    stats.data['verify'][-1].update({ 'solve_memory' : max_memory })

                    # determine SAT status and create markers
                    satisfiable = determine_satisfiability(stdout)
                    stats.data['verify'][-1].update({ 'satisfiable' : satisfiable })
                    if satisfiable == 'True':
                        successful = True
                        marker.with_suffix('.YES').touch()
                        logger.info(colored(f"Successful verification of '{gzvtheory}'", 'cyan', attrs=['bold']))
                    elif satisfiable == 'False':
                        marker.with_suffix('.NO').touch()

                    # remove compressed theory, if SAT or UNSAT (keep it if INDET)
                    if satisfiable == 'True' or satisfiable == 'False':
                        gzvtheory.unlink()

                    # compress model if SAT, otherwise remove model
                    if successful:
                        logger.info(colored(f"Compress vmodel '{vmodel}'", 'green'))
                        stdout, stderr, elapsed_time, max_memory = execute_cmd(g_templates['gzip'].format(gzinput = parameters['vmodel']), logger)
                    else:
                        vmodel.unlink()

                    # break if satisfiable or indeterminate, continue otherwise
                    #if satisfiable == 'True' or satisfiable == 'Indet': break

                    # break if satisfiable, continue otherwise
                    if satisfiable == 'True': break

            # if not successful, stop verification and jump to next benchmark
            if not successful: 
                logger.info(colored(f"Unsuccessful verification of '{dfa}' ... skipping remaining verifications", 'red'))
                break
            else:
                complete_verification = True

        # remove marker 'incomplete.marker'
        if complete_verification:
            remove_marker(parameters['dirpath'], logger)

        # close log, and output stats
        logger.info(f'stats: {str(stats)}')
        close_logger(logger)
        with stats_name.open('w') as fd: fd.write(str(stats) + '\n')

if __name__ == '__main__':
    # setup exec name
    exec_path = Path(sys.argv[0]).parent
    exec_name = Path(sys.argv[0]).stem

    # setup proper SIGTERM handler
    signal.signal(signal.SIGTERM, sigterm_handler)

    # read arguments and do job
    args = get_args()
    main(args)

