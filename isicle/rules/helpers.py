import platform


def getOS():
    system = platform.system().lower()
    if system == 'darwin':
        return 'osx'
    return system


def cycles(n):
    return ['%03d' % x for x in range(1, n + 1)]


def frames(n):
    return ['%03d' % x for x in range(n)]
