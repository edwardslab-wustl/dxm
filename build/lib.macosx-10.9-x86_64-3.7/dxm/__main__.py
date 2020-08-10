#!usr/bin/python3
import sys

def main():
    if __package__ == '':
        import os.path
        path = os.path.dirname(os.path.dirname(__file__))
        sys.path[0:0] = [path]
    from dxm import runDXM
    sys.exit(runDXM.main())

if __name__ == '__main__':
    print('executed')
    sys.exit(main())
