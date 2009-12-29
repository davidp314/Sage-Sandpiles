from sage.misc.sageinspect import *

"""
last modified: 12/22/09

These are utilities for helping to create the Sage Sandpiles documentation.

To create the list of methods:

1. sage: m = dir(Sandpile).
2. Manually remove unwanted methods from m and alphabetize disregarding
case.
3. run create_all_refs('sandpile_methods.txt', 'Sandpile', m)
4. clean up sandpile_methods.txt by hand.
5. run create_all_targets('sandpile_targets.txt', 'Sandpile', m)
6. clean up sandpile_targets.txt by hand.
7. shorten comments for references by hand.
8. Document utility functions by hand.
9. Add default values for arguments by hand.  (Can this be coded?)
10. Repeat for Divisor and Config methods.
11. add "-divisor" to the ends of some links in order to make them unique
12. Try to make method summaries each take one line.

Will need to get rid of underscores in names for links: __neg__ should be made
neg, for instance.  To get <= or < as links, wait until the end and edit the
html or tex by hand.

Target names cannot have extra spaces.  So use, e.g., <random_tree(n,d> rather
than <random_tree(n, d)>.

13. May need to edit the resulting html file to get links with "<" rather than
"less-than" and may need to edit the tex file to get good references to figures
(that have floated), as in "Consider the graph", on page 1 or 2.  
"""

def get_short_description(my_class, m):
    r"""
    Get short description of the my_class method with name m.

    INPUT: 

    ``my_class`` - text : "Sandpile", "Config", "Divisor"
    ``m`` - text

    OUTPUT:

    text

    NOTES:

    Must first issue the command:

        sage: from sage.misc.sageinspect import *
    """
    text  = sage_getdoc(eval(my_class + '.' + m)).split('\n')
    desc = ""
    for i in text[1:]:
        if i == '':
            break
        else:
            desc += i
    return desc

def create_ref(my_class, name):
    r"""
    Create a link reference in reST.  Clicking on the reference takes you to the
    target.

    INPUT:

    ``my_class`` - text : "Sandpile", "Config", "Divisor"
    ``name`` - text

    OUTPUT:

    text
    """
    title = name+sage_getdef(eval(my_class+'.'+name))
    reference = '- :ref:`' + title
    reference += ' <' + title.replace(' ','') + '>`\n\n'
    reference += get_short_description(my_class, name) + '\n\n'
    return reference

def create_all_refs(file_name, my_class, methods):
    r"""
    Write all reST references for the my_class methods provided.

    INPUT:

    - ``file_name`` - text (for file to be written to)
    - ``my_class`` - text : "Sandpile", "Config", "Divisor"
    - ``methods`` - list of strings identifying Sandpile methods
    """
    file = open(file_name, 'w')
    for m in methods:
        ref = create_ref(my_class, m)
	file.write(ref)
    file.close()

def create_all_targets(file_name, my_class, methods):
    r"""
    Write all reST targets for the my_class methods provided.

    INPUT:

    - ``file_name`` - text (for file to be written to)
    - ``my_class`` - text : "Sandpile", "Config", "Divisor"
    - ``methods`` - list of strings identifying Sandpile methods
    """
    file = open(file_name, 'w')
    for name in methods:
        title = name + sage_getdef(eval(my_class + '.' + name))
	target = '\n---\n\n'
	target += '.. _' + title.replace(' ','') + ':' + '\n\n'
	target += '**' + title + '**' + '\n'
	target += sage_getdoc(eval(my_class + '.' + name))
	file.write(target)
    file.close()

def get_short_description(my_class, m):
    r"""
    Get short description of the my_class method with name m.

    INPUT: 

    ``my_class`` - text : "Sandpile", "Config", "Divisor"
    ``m`` - text

    OUTPUT:

    text

    NOTES:

    Must first issue the command:

        sage: from sage.misc.sageinspect import *
    """
    text  = sage_getdoc(eval(my_class + '.' + m)).split('\n')
    desc = ""
    for i in text[1:]:
        if i == '':
            break
        else:
            desc += i
    return desc

def create_ref(my_class, name):
    r"""
    Create a link reference in reST.  Clicking on the reference takes you to the
    target.

    INPUT:

    ``my_class`` - text : "Sandpile", "Config", "Divisor"
    ``name`` - text

    OUTPUT:

    text
    """
    title = name+sage_getdef(eval(my_class+'.'+name))
    reference = '- :ref:`' + title
    reference += ' <' + title.replace(' ','') + '>`\n\n'
    reference += get_short_description(my_class, name) + '\n\n'
    return reference

def create_all_refs(file_name, my_class, methods):
    r"""
    Write all reST references for the my_class methods provided.

    INPUT:

    - ``file_name`` - text (for file to be written to)
    - ``my_class`` - text : "Sandpile", "Config", "Divisor"
    - ``methods`` - list of strings identifying Sandpile methods
    """
    file = open(file_name, 'w')
    for m in methods:
        ref = create_ref(my_class, m)
	file.write(ref)
    file.close()

def create_all_targets(file_name, my_class, methods):
    r"""
    Write all reST targets for the my_class methods provided.

    INPUT:

    - ``file_name`` - text (for file to be written to)
    - ``my_class`` - text : "Sandpile", "Config", "Divisor"
    - ``methods`` - list of strings identifying Sandpile methods
    """
    file = open(file_name, 'w')
    for name in methods:
        title = name + sage_getdef(eval(my_class + '.' + name))
	target = '\n---\n\n'
	target += '.. _' + title.replace(' ','') + ':' + '\n\n'
	target += '**' + title + '**' + '\n'
	target += sage_getdoc(eval(my_class + '.' + name))
	file.write(target)
    file.close()

"""
selected Sandpile methods list, cleaned up
"""
m= ['all_k_config',
 'all_k_div',
 'betti',
 'betti_complexes',
 'burning_config',
 'burning_script',
 'canonical_divisor',
 'dict',
 'elementary_divisors',
 'groebner',
 'group_order',
 'h_vector',
 'hilbert_function',
 'ideal',
 'identity',
 'in_degree',
 'is_undirected',
 'laplacian',
 'max_stable',
 'max_stable_div',
 'max_superstables',
 'min_recurrents',
 'nonsink_vertices',
 'nonspecial_divisors',
 'num_edges',
 'out_degree',
 'points',
 'postulation',
 'recurrents',
 'reduced_laplacian',
 'reorder_vertices',
 'resolution',
 'ring',
 'sink',
 'solve',
 'superstables',
 'symmetric_recurrents',
 'unsaturated_ideal',
 'version',
 'vertices',
 'zero_config',
 'zero_div']

c = ['__add__',
 '__and__',
 '__invert__',
 '__le__',
 '__lt__',
 '__mul__',
 '__neg__',
 '__pow__',
 '__sub__',
 'add_random',
 'deg',
 'dualize',
 'equivalent_recurrent',
 'equivalent_superstable',
 'fire_script',
 'fire_unstable',
 'fire_vertex',
 'is_recurrent',
 'is_stable',
 'is_superstable',
 'is_symmetric',
 'order',
 'stabilize',
 'support',
 'unstable',
 'values']

d=['__add__',
 '__le__',
 '__lt__',
 '__neg__',
 '__sub__',
 'add_random',
 'betti',
 'Dcomplex',
 'deg',
 'dualize',
 'effective_div',
 'fire_script',
 'fire_unstable',
 'fire_vertex',
 'is_alive',
 'is_symmetric',
 'linear_system',
 'r_of_D',
 'support',
 'unstable',
 'values']
