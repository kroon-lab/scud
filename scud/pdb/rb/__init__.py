import sys

import rb


def main():
    l = rb.run(args=sys.argv[1:])
    l.close_log()
    
if __name__ == '__main__':
    main()
