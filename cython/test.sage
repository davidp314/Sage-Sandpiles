def mytest():
    S = sandlib('generic')
    c = S.add(S.max_stable(),S.max_stable()).values()
    return S, c
