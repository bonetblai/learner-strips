import sys, argparse
from termcolor import colored
from pathlib import Path
from typing import List
import logging

from utils.dfa import DFA

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

def get_args():
    # default values
    default_debug_level = 0
    default_log_file = 'make_lp_from_dfa.log'

    # argument parser
    parser = argparse.ArgumentParser(description='Translate .dfa file (labeled directed graph) into .lp file that can be processed by clingo.')
    parser.add_argument('--debug_level', nargs=1, type=int, default=default_debug_level, help='Set debug level')
    parser.add_argument('--log_file', nargs=1, type=Path, default=default_log_file, help=f"Set log filename (default='{default_log_file}')")
    parser.add_argument('dfa', type=Path, help='Filename or folder for .dfa containing state space graph')
    parser.add_argument('lp', type=Path, help='Filename or folder for output .lp file')
    parser.add_argument('--suppress_labels', action='store_true', help='Surpress tlabel facts in target file for blackbox edges')
    parser.add_argument('--compute_inverse', action='store_true', help='Precompute inverse relation of actions')

    # parse arguments
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # read arguments
    args = get_args()

    # setup logger
    log_level = logging.INFO if args.debug_level == 0 else logging.DEBUG
    logger = get_logger('make_lp_from_dfa', 'make_lp_from_dfa.log', log_level)

    # setup input/output filenames
    input_fnames = []
    output_fnames = []

    # if args.lp is a folder, create it
    if args.lp.suffix != '.lp':
        logger.info(f"Creating folder '{args.lp}'")
        args.lp.mkdir(parents=True, exist_ok=True)

    # if args.dfa is a folder, process all .dfa files in it
    if args.dfa.is_dir():
        assert args.lp.is_dir()
        input_fnames = [fname for fname in args.dfa.glob('*.dfa') if fname.is_file()]
        output_fnames = [(args.lp / fname.name).with_suffix('.lp') for fname in input_fnames]
    else:
        input_fnames.append(args.dfa)
        lp_fname = (args.lp / args.dfa.name).with_suffix('.lp') if args.lp.is_dir() else args.lp
        output_fnames.append(lp_fname)

    # do job
    for input_fname, output_fname in zip(input_fnames, output_fnames):
        dfa = DFA(input_fname, logger)
        with output_fname.open('w') as fd:
            logger.info(f"Writing '{output_fname}'")
            dfa.dump_as_lp(fd,args.suppress_labels,args.compute_inverse)
    close_logger(logger)

