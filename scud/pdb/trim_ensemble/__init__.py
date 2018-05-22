import sys

import trim_ensemble


def main():
    l = trim_ensemble.run(args=sys.argv[1:])
    l.close_log()
    
if __name__ == '__main__':
    main()
