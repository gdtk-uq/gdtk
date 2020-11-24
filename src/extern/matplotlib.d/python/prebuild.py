#!/usr/bin/env python
# -*- coding: utf-8 -*-


def extract_function_names(module):
    '''
    extract function names from attributes of 'module'.
    '''
    from importlib import import_module

    mod = import_module(module.__name__)
    attr_list = dir(mod)
    scope = locals()

    def iscallable(name):
        ignore_decorators = ['dedent','deprecated','silent_list', 'warn_deprecated']
        return eval('callable(mod.{})'.format(name), scope) and name not in ignore_decorators

    return filter(iscallable, attr_list)


def gen_pyplot_functions(dub_root):
    '''
    generate 'pyplot_functions.txt' for matplotlibd.pyplot.
    '''
    import matplotlib.pyplot
    from string import ascii_lowercase

    functions = filter(lambda i: i[0] != '_' or i[0] in ascii_lowercase,
                       extract_function_names(matplotlib.pyplot))

    with open(dub_root + "/views/pyplot_functions.txt", "w") as f:
        f.write("\n".join(functions))


if __name__ == '__main__':
    from sys import argv
    gen_pyplot_functions(argv[1])
