from guppylang import guppy
from guppylang.std.builtins import owned, array, exec
from guppylang.std.quantum import cx, h, measure, qubit, x, z
import numpy as np



@guppy
def foobar(q: qubit @owned):  
    return q+1

@guppy
def main():
    q = qubit()
    result = foobar(q)
    print("Result:", measure(result))

program = main.compile()
exec(program)