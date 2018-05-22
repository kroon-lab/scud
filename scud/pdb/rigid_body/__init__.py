import sys

import rigid_body


def main():
    l = rigid_body.run(args=sys.argv[1:])
    l.close_log()
    
if __name__ == '__main__':
    main()
