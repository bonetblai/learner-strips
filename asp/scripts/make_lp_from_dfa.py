import sys, argparse
from termcolor import colored
from pathlib import Path
from typing import List
import logging

from dfa import DFA

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
    parser.add_argument('dfa', type=Path, help='Filename of .dfa containing state space graph')
    parser.add_argument('lp', type=Path, help='Filename or folder for output .lp file')

    # parse arguments
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # read arguments
    args = get_args()

    # setup logger
    log_level = logging.INFO if args.debug_level == 0 else logging.DEBUG
    logger = get_logger('make_lp_from_dfa', 'make_lp_from_dfa.log', log_level)

    # set lpfname taking into account that args.lp may be a folder
    lp_fname = (args.lp / args.dfa.name).with_suffix('.lp') if args.lp.is_dir() else args.lp

    # do job
    dfa = DFA(args.dfa, logger)
    with lp_fname.open('w') as fd:
        logger.info(f"Writing '{lp_fname}'")
        dfa.dump_as_lp(fd)

    close_logger(logger)

