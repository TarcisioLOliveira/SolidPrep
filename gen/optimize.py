#!/usr/bin/python

from solid import *
import numpy as np
import sys
from pathlib import Path
import regex as re

def eval_section(match):
    brack_re = r'(?P<call>\w+\((.*?)\))\s+(?P<inside>(\w+\((.*?)\)\s*).*?);'
    brack_sub = r'\g<call>{\n\g<inside>;\n}'

    cur_str = match.group(0)
    new_str = re.sub(brack_re, brack_sub, cur_str, re.S)
    # print(new_str)
    # while cur_str != new_str:
    #     cur_str = new_str
    #     new_str = re.sub(brack_re, brack_sub, cur_str)
    # cur_str = new_str

    brack = re.match(brack_re, cur_str)
    if brack:
        call = brack.group("call")
        inside = brack.group("inside")
        cur_str = call+"{\n"+inside+";\n}"
        brack = re.match(brack_re, cur_str)
        while brack != None:
            call = brack.group("call")
            inside = brack.group("inside")
            cur_str = call+"{\n"+inside+";\n}"
            brack = re.match(brack_re, cur_str)

    print(cur_str)


    cur_str = re.sub(';',  ',', cur_str)
    cur_str = re.sub('{', '(', cur_str)
    cur_str = re.sub('}', ')', cur_str)
    cur_str = re.sub('true', 'True', cur_str)
    cur_str = re.sub('false', 'False', cur_str)

    return eval(cur_str)

def simple_parse_scad(path):
    """
    Parses a SCAD file to get the primitives and operations themselves.

    Based on solidpython's parse_scad_callables().
    """
    no_comments_re = r'(?mxs)(//.*?\n|/\*.*?\*/)'

    # Prevent eval() from causing too much havoc
    types_re = r'(?mxs)(?:union|difference|intersection|translate|rotate|color|mirror|multmatrix|scale|resize|offset|minkowski|hull)\s*\((.*?)\)\s*{(.*?)*}|(?mxs)(circle|square|polygon|sphere|cube|polyhedron|translate|rotate|color|mirror|multmatrix|scale|resize|offset|minkowski|hull)\s*\((.*?)\);'

    scad_code_str = path.read_text()

    scad_code_str = re.sub(no_comments_re, '', scad_code_str)

    calls = re.finditer(types_re, scad_code_str, re.S)
    geom = eval_section(next(calls))

    for c in calls:
        geom += eval_section(c)
        print(scad_render(geom))


simple_parse_scad(Path('square.scad'))
